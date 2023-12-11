// Need to use actual GNU gcc to compile this code as Apple's
// alias to gcc, clang, does not support OpenMP
// g++-13  -I ./eigen -std=c++20  test.cpp -o test

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <cmath>
#include <iostream>

// number of points of each direction; L for x, W for y, and H for z
int L;
int H;

double N; // Number of circulation cells
// Ub for constant term for u function
double Ub = 1; // cm per second
// D for diffusion term
double D = 0.185; // cm squared per second

// Used to print out a vertical cross section of the conentration
void printVertical(double **c) {
    for (int k = 1; k < H; k++) {
        for (int i = 1; i < L - 1; i++) {
            printf("%f ", c[i][k]);
        }
        printf("\n");
    }
    printf("\n");
}

double CFL(double dx, double **ux, double **uz) {

    // Need to find largest flow velocity among all the cells
    double maxU = 0;
    double U;
    for (int i = 0; i < L; i++) {
        for (int k = 0; k < H; k++) {
            U = std::sqrt(ux[i][k] * ux[i][k] + uz[i][k] * uz[i][k]);
            if (U > maxU) {
                maxU = U;
            }
        }
    }

    // Returning the smaller CFL condition between advection and diffusion
    // For small diffusion, advection CFL is normally much smaller
    return .5 * std::min(dx / maxU, (dx * dx) / (2 * D));
}

// Given indicies, returns physical x and z coordinates
std::array<double, 2> convert(int i, int k, double dx) {
    std::array<double, 2> coord;
    double x = (double(i) - double(L) * .5) * dx + (.5) * dx;
    double z = (double(k) - (double(H) - 1)) * dx + (.5) * dx;
    coord = {x, z};
    return coord;
}

// Used to check if a cell is a ghost node or not
bool isBoundary(int i, int k) {
    return (i == 0 || i == L - 1 || k == 0 || k == H - 1);
}

double initializeC(int i, int k, double dx) {
    int a = (i + k) * dx * 0;
    return 0.004 + a; // const for now, grams per liter?
}

double initializeUx(int i, int k, double dx, double LL, double HH) {
    std::array<double, 2> coord = convert(i, k, dx);
    double x = coord[0];
    double z = coord[1];
    double ux = LL * (1 / (2 * HH)) * (1 / ((2 * N) - 1)) * cos((2 * N - 1) * (M_PI) * (2 * x - LL / 2) * (1 / LL)) * cos(M_PI * (z / HH));
    ux = 0.1;
    return ux;
}

double initializeUz(int i, int k, double dx, double LL, double HH) {
    std::array<double, 2> coord = convert(i, k, dx);
    double x = coord[0];
    double z = coord[1];
    double uz = sin((2 * N - 1) * (2 * x - LL / 2) * (M_PI / LL)) * sin(M_PI * (z / HH));
    uz = 0.1;
    return uz;
}

double **allocate2DArray() {
    double **arr = new double *[L];
    for (int i = 0; i < L; i++) {
        arr[i] = new double[H];
    }
    return arr;
}

void deallocate2DArray(double **arr) {
    for (int i = 0; i < L; ++i) {
        delete[] arr[i];
    }
    delete[] arr;
}

void process2DArray(double **arr) {
    // Access and manipulate the elements of the 3D array
    for (int k = H - 1; k > -1; k--) {
        for (int i = 0; i < L; i++) {
            std::cout << arr[i][k] << " ";
        }
        printf("\n");
    }
    std::cout << std::endl;
}

void processMatrix(const Eigen::MatrixXd &mat) {
    for (int i = 0; i < mat.rows(); i++) {
        for (int j = 0; j < mat.cols(); j++) {
            std::cout << mat(i, j) << "\t";
        }
        printf("\n");
    }
    std::cout << std::endl;
}

void computeC(double **c, double dx) {
    for (int i = 0; i < L; i++)
        for (int k = 0; k < H; k++)
            c[i][k] = initializeC(i, k, dx);
}

// TODO： detect ux<0?
void computeUx(double **ux, double dx, double LL, double HH) {
    for (int i = 0; i < L; i++)
        for (int k = 0; k < H; k++)
            if (isBoundary(i, k)) {
                ux[i][k] = 0;
            } else {
                ux[i][k] = initializeUx(i, k, dx, LL, HH);
            }
}

void computeUz(double **uz, double dx, double LL, double HH) {
    for (int i = 0; i < L; i++)
        for (int k = 0; k < H; k++)
            if (isBoundary(i, k)) {
                uz[i][k] = 0;
            } else {
                uz[i][k] = initializeUz(i, k, dx, LL, HH);
            }
}

void toVector(Eigen::VectorXd &v, double **array2D) {
    for (int k = 0; k < H; k++) {
        for (int i = 0; i < L; i++) {
            v(i + k * L) = array2D[i][k];
        }
    }
}

void toMatrix(Eigen::VectorXd &v, double **array2D) {
    for (int k = 0; k < H; k++) {
        for (int i = 0; i < L; i++) {
            array2D[i][k] = v(i + k * L);
        }
    }
}

void toDx(Eigen::MatrixXd &Dx, double **ux, double dt, double dx) {
    double factor = dt / (2 * dx);
    for (int index = 0; index < L * H - 1; index++) {
        int i = index % L;
        int k = index / L;

        /*
        if (i == 0) {
            Dx(index, index) = 1; // leftmost ghost
        } else if (i == L - 1) {
            Dx(index, index) = 1; // rightmost ghost
        } else {
            Dx(index, index - 1) = 1 * ux[i - 1][k] * factor;
            Dx(index, index + 1) = -1 * ux[i + 1][k] * factor;
        }
        */

        if (isBoundary(i, k)) {
            Dx(index, index) = 0;
        } else {
            Dx(index, index - 1) = 1 * ux[i - 1][k] * factor;
            Dx(index, index + 1) = -1 * ux[i + 1][k] * factor;
        }
    }
}

void toDy(Eigen::MatrixXd &Dy, double **uy, double dt, double dx) {
    double factor = dt / (2 * dx);
    for (int index = 0; index < L * H - 1; index++) {
        int i = index % L;
        int k = index / L;

        /*
        if (k == 0) {
            Dy(index, index) = 1; // bottommost ghost
        } else if (k == H - 1) {
            Dy(index, index) = 1; // topmost ghost
        } else {
            Dy(index, index - L) = 1 * uy[i][k - 1] * factor;
            Dy(index, index + L) = -1 * uy[i][k + 1] * factor;
        }
        */

        if (isBoundary(i, k)) {
            Dy(index, index) = 0;
        } else {
            Dy(index, index - L) = 1 * uy[i][k - 1] * factor;
            Dy(index, index + L) = -1 * uy[i][k + 1] * factor;
        }
    }
}

void toDxx(Eigen::MatrixXd &Dxx, double dt, double dx) {
    double factor = D * dt / (2 * dx * dx);
    for (int index = 0; index < L * H - 1; index++) {
        int i = index % L;
        int k = index / L;
        /*
        if (i == 0) {
            Dxx(index, index) = 1; // leftmost ghost
        } else if (i == L - 1) {
            Dxx(index, index) = 1; // rightmost ghost
        } else {
            Dxx(index, index - 1) = factor;
            Dxx(index, index) = -2 * factor;
            Dxx(index, index + 1) = factor;
        }
        */

        if (isBoundary(i, k)) {
            Dxx(index, index) = 0;
        } else {
            Dxx(index, index - 1) = factor;
            Dxx(index, index) = -2 * factor;
            Dxx(index, index + 1) = factor;
        }
    }
}

void toDyy(Eigen::MatrixXd &Dyy, double dt, double dx) {
    double factor = D * dt / (2 * dx * dx);
    for (int index = 0; index < L * H - 1; index++) {
        int i = index % L;
        int k = index / L;
        /*
        if (k == 0) {
            Dyy(index, index) = 1; // leftmost ghost
        } else if (k == H - 1) {
            Dyy(index, index) = 1; // rightmost ghost
        } else {
            Dyy(index, index - L) = factor;
            Dyy(index, index) = -2 * factor;
            Dyy(index, index + L) = factor;
        }
        */
        if (isBoundary(i, k)) {
            Dyy(index, index) = 0;
        } else {
            Dyy(index, index - L) = factor;
            Dyy(index, index) = -2 * factor;
            Dyy(index, index + L) = factor;
        }
    }
}

void enforceBoundary_ADI(Eigen::VectorXd &cv) {
    for (int index = 0; index < L * H - 1; index++) {
        int i = index % L;
        int k = index / L;

        if (i == 0) {
            cv(index) = cv(index + 1);
        } else if (i == L - 1) {
            cv(index) = cv(index - 1);
        }

        if (k == 0) {
            cv(index) = cv(index + L);
        } else if (k == H - 1 || k == H - 2) {
            cv(index) = 0;
        }
    }
}

void start_ADI(double LL, double HH, double dx, double T) {
    L = int(LL / dx) + 2; // L is number of points, + 2 points for ghost node on each side
    H = int(HH / dx) + 2;

    double **c = allocate2DArray();
    computeC(c, dx);
    // process2DArray(c);
    Eigen::VectorXd cv(L * H);
    toVector(cv, c);

    double **ux = allocate2DArray();
    double **uz = allocate2DArray();
    computeUx(ux, dx, LL, HH);
    computeUz(uz, dx, LL, HH);
    std::cout << "ux and uz" << std::endl;

    double dt = CFL(dx, ux, uz);

    Eigen::MatrixXd Dx(L * H, L * H);
    Eigen::MatrixXd Dz(L * H, L * H);
    Eigen::MatrixXd Dxx(L * H, L * H);
    Eigen::MatrixXd Dzz(L * H, L * H);
    toDx(Dx, ux, dt, dx);
    std::cout << "dx" << std::endl;
    toDy(Dz, uz, dt, dx);
    std::cout << "dz" << std::endl;
    toDxx(Dxx, dt, dx);
    std::cout << "dxx" << std::endl;
    toDyy(Dzz, dt, dx);
    std::cout << "dzz" << std::endl;

    // process2DArray(ux);
    // process2DArray(uz);
    // processMatrix(Dzz);

    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(L * H, L * H);
    Eigen::MatrixXd temp1 = I + Dz - Dzz;
    Eigen::SparseMatrix<double> temp1_sparse = temp1.sparseView();
    temp1_sparse.makeCompressed();
    Eigen::MatrixXd temp2 = I + Dxx - Dx;
    Eigen::SparseMatrix<double> temp2_sparse = temp2.sparseView();
    temp2_sparse.makeCompressed();
    Eigen::MatrixXd temp3 = I + Dx - Dxx;
    Eigen::SparseMatrix<double> temp3_sparse = temp3.sparseView();
    temp3_sparse.makeCompressed();
    Eigen::MatrixXd temp4 = I + Dzz - Dz;
    Eigen::SparseMatrix<double> temp4_sparse = temp4.sparseView();
    temp4_sparse.makeCompressed();

    Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IdentityPreconditioner> solver;
    Eigen::SparseMatrix<double> II(L * H, L * H);

    solver.compute(temp1_sparse);
    II.setIdentity();
    auto A_inv = solver.solve(II);
    std::cout << "inverse 1" << std::endl;
    solver.compute(temp3_sparse);
    II.setIdentity();
    auto B_inv = solver.solve(II);
    std::cout << "inverse 2" << std::endl;

    Eigen::MatrixXd updated(L * H, L * H);
    updated = A_inv * temp2_sparse * B_inv * temp4_sparse;
    std::cout << "inverse" << std::endl;

    enforceBoundary_ADI(cv);

    double t_total = 0;
    std::cout << "We expect to have " << T / dt << " iterations" << std::endl;
    while (t_total < T) {
        cv = updated * cv;
        enforceBoundary_ADI(cv);
        t_total += dt;
    }

    int images = std::floor(T / dt) + 1;
    fprintf(stderr, "Final time: %f, dt: %f, Images: %d, Width pixels: %d, Height pixels: %d\n", T, dt, images, L - 2, H - 1);

    toMatrix(cv, c);
}

int main(int argc, char *argv[]) {

    if (argc != 7) {
        fprintf(stderr, "Usage: ./advDif2D <adams:0, ADI:1> <Container length> <Container height> <Number of circulation cells> <Grid spacing> <Final time>");
        exit(0);
    }

    double LL, HH, dx, T;

    int method = atoi(argv[1]);
    LL = atof(argv[2]);
    HH = atof(argv[3]);
    N = atof(argv[4]);
    dx = atof(argv[5]);
    T = atof(argv[6]);

    start_ADI(LL, HH, dx, T);
    return 0;
}
