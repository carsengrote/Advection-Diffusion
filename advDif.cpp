// g++ -std=c++11 -Wall -Wextra -Werror advDif.cpp -o advDif

#include "advDif.hh"
#include <iostream>
#include <cmath>
// number of points of each direction; L for x, W for y, and H for z
int L;
int W;
int H;
// Ub for constant term for u function
double Ub = 0.01;
// D for diffusion term
double D = 0.01;

int main()
{
    double LL, WW, HH, dx, T;
    std::cout << "Enter length, width, height, and dx, and total time: ";
    std::cin >> LL >> WW >> HH >> dx >> T;
    if (LL != WW)
    {
        WW = LL; // We do this because we want to keep xy-plane symmetry; if we have some other u, we can delete this
    }
    start(LL, WW, HH, dx, T);
    return 0;
}

void start(double LL, double WW, double HH, double dx, double T)
{
    double dt = CFL(dx); // How?

    L = int(LL / dx) + 1; // L is the number of points, and has both sides boundary points and no ghost nodes
    W = int(WW / dx) + 1;
    H = int(HH / dx) + 1;
    // test input:
    // std::cout << "L is: " << L << "W is: " << W << "H is: " << H << "dx is: " << dx << std::endl;

    double ***c = allocate3DArray();
    computeC(c, dx);
    // test c initialization:
    // process3DArray(c);

    double ***ux = allocate3DArray();
    double ***uy = allocate3DArray();
    double ***uz = allocate3DArray();
    computeUx(ux, dx);
    computeUy(uy, dx);
    computeUz(uz, dx);
    // test u initialization:
    // process3DArray(ux);
    // process3DArray(uy);
    // process3DArray(uz);

    double t_total = 0;
    while (t_total < T)
    {
        updateC(c, ux, uy, uz, dx, dt);
        t_total += dt;
    }
}

double CFL(double dx)
{
    return dx * dx / 10;
}

std::array<double, 3> convert(int i, int j, int k, double dx) // if I only use it once, i can pass value back, and it'll be deleted automatically
{
    std::array<double, 3> coord;
    double x = (double(i) - L / 2) * dx;
    double y = (double(j) - W / 2) * dx;
    double z = (double(k) - (H - 1)) * dx;
    coord = {x, y, z};
    // test convert: it looked good to me: I run it on a 4*4*4 with dx=1 case, and it's ideal
    // std::cout << "right: " << x << y << z << std::endl;
    return coord;
}

bool isBoundary(int i, int j, int k)
{
    return (i == 0 || i == L - 1 || j == 0 || j == W - 1 || k == 0 || k == H - 1);
}

void enforceBoundary(double ***c)
{
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < W; j++)
        {
            c[i][j][H - 1] = 0;
        }
    }
}

double initializeC(int i, int j, int k, double dx)
{
    // std::array<double, 3> coord = convert(i, j, k, dx); // for future position-dependecy initialization
    // std::cout << coord[0] << " " << coord[1] << " " << coord[2] << std::endl;

    // why do I define a? bc if you dont use ijk, there is compile error; I used the first line of this file to compile
    int a = (i + j + k) * dx * 0;
    return 0.1 + a; // const for now;
}

double initializeUx(int i, int j, int k, double dx)
{
    std::array<double, 3> coord = convert(i, j, k, dx);
    double x = coord[0];
    double y = coord[1];
    double r = std::sqrt(x * x + y * y);
    double z = coord[2];

    // fuck it, in notes, we put L as height and H is the xy plane dimension; in the code we flipped them; be extra caustious here
    double ur = Ub * sin(M_PI * r / L) * sin(M_PI * z / H);
    double ux;
    if (r < dx / 2) // handle r=0 case
    {
        ux = 0;
    }
    else
    {
        ux = ur * coord[0] / r; // x/r is the component
    }
    return ux;
}

double initializeUy(int i, int j, int k, double dx)
{
    std::array<double, 3> coord = convert(i, j, k, dx);
    double x = coord[0];
    double y = coord[1];
    double r = std::sqrt(x * x + y * y);
    double z = coord[2];

    // fuck it, in notes, we put L as height and H is the xy plane dimension; in the code we flipped them; be extra caustious here
    double ur = Ub * sin(M_PI * r / L) * sin(M_PI * z / H);
    double uy;
    if (r < dx / 2) // handle r=0 case
    {
        uy = 0;
    }
    else
    {
        uy = ur * coord[1] / r; // y/r is the component
    }
    return uy;
}

double initializeUz(int i, int j, int k, double dx)
{
    std::array<double, 3> coord = convert(i, j, k, dx);
    double r = std::sqrt(std::pow(coord[0], 2) + std::pow(coord[1], 2)); // sqrt(x^2+y^2)
    double z = coord[2];

    // fuck it, in notes, we put L as height and H is the xy plane dimension; in the code we flipped them; be extra caustious here
    double term1 = -Ub * H * cos(M_PI * r / L) * sin(M_PI * z / H);
    double term2 = 0;
    if (r < dx / 2) // handle r=0 case
    {
        term2 = -Ub * H / L * sin(M_PI * z / H); // if r=0, we have a sinc function, so it's 1; it's H/L instead of L/H because of notation definition
    }
    else
    {
        term2 = -Ub * 1 / (M_PI * r) * (H * sin(M_PI * r / L)) * sin(M_PI * z / H);
    }
    double uz = term1 + term2;
    return uz;
}

double getC(double ***c, int i, int j, int k)
{
    return c[i][j][k];
}

std::array<double, 3> getU(double ***ux, double ***uy, double ***uz, int i, int j, int k)
{
    std::array<double, 3> u;
    double x = ux[i][j][k];
    double y = uy[i][j][k];
    double z = uz[i][j][k];
    u = {x, y, z};
    return u;
}
/*here comes the big one~~~~~~~~~~~~~~~~ update the c 3d array ~~~~~~~~~~~~~~~~~~~~~~~~*/

double updateC_cell(double ***c, double ***ux, double ***uy, double ***uz, int i, int j, int k, double dx, double dt)
{
    double c_center = getC(c, i, j, k);
    std::array<double, 3> u_center = getU(ux, uy, uz, i, j, k);

    // xyz direction advection and diffusion
    double x_ad;
    double y_ad;
    double z_ad;
    double x_df;
    double y_df;
    double z_df;

    double total_ad;
    double total_df;

    // u_x_{i+1} and c_x_{i+1}
    if (i + 1 < L)
    {
        std::array<double, 3> u_west = getU(ux, uy, uz, i + 1, j, k);
        double c_west = getC(c, i + 1, j, k);
        x_ad = u_west[0] * c_west - u_center[0] * c_center; //[0] because of x-component
    }
    else
    {
        x_ad = 0; // BC requires advection to be zero
    }

    // u_y_{i+1} and c_y_{i+1}
    if (j + 1 < W)
    {
        std::array<double, 3> u_south = getU(ux, uy, uz, i, j + 1, k);
        double c_south = getC(c, i, j + 1, k);
        y_ad = u_south[1] * c_south - u_center[1] * c_center; //[1] because of y-component
    }
    else
    {
        y_ad = 0; // BC requires advection to be zero
    }

    // u_z_{i+1} and c_z_{i+1}
    if (k + 1 < H - 1) // This is we need more thinking: since top layer is just zero, probabaly we are just interested in second layer at the top?
    {
        std::array<double, 3> u_top = getU(ux, uy, uz, i, j, k + 1);
        double c_top = getC(c, i, j, k + 1);
        z_ad = u_top[2] * c_top - u_center[2] * c_center; //[2] because of z-component
    }
    else
    {
        z_ad = c_center; // This is we need more thinking: at second layer to the top, we should get rid of everything since top is just 0
    }

    total_ad = -dt / (2 * dx) * (x_ad + y_ad + z_ad);

    // c_x_{i-1} and c_x_{i} and c_x_{i+1}
    if (i == 0)
        x_df = -getC(c, i + 1, j, k) + getC(c, i, j, k);
    else if (i + 1 == L)
        x_df = -getC(c, i - 1, j, k) + getC(c, i, j, k);
    else // interior points
        x_df = -getC(c, i - 1, j, k) + 2 * getC(c, i, j, k) - getC(c, i + 1, j, k);

    // c_z_{i-1} and c_z_{i} and c_z_{i+1}
    if (j == 0)
        y_df = -getC(c, i, j + 1, k) + getC(c, i, j, k);
    else if (j + 1 == W)
        y_df = -getC(c, i, j - 1, k) + getC(c, i, j, k);
    else // interior points
        y_df = -getC(c, i, j - 1, k) + 2 * getC(c, i, j, k) - getC(c, i, j + 1, k);

    // c_y_{i-1} and c_y_{i} and c_y_{i+1}
    if (k == 0)
        z_df = -getC(c, i, j, k + 1) + getC(c, i, j, k);
    else if (k + 1 == H)
        z_df = -getC(c, i, j, k - 1) + getC(c, i, j, k);
    else // interior points
        z_df = -getC(c, i, j, k - 1) + 2 * getC(c, i, j, k) - getC(c, i, j, k + 1);

    total_df = -D * dt / (2 * dx * dx) * (x_df + y_df + z_df);

    return c[i][j][k] + total_ad + total_df;
}

void updateC(double ***c, double ***ux, double ***uy, double ***uz, double dx, double dt)
{
    for (int i = 0; i < L; i++)
    {
        for (int j = 0; j < W; j++)
        {
            for (int k = 0; k < H; k++)
            {
                c[i][j][k] = updateC_cell(c, ux, uy, uz, i, j, k, dx, dt);
            }
        }
    }
    enforceBoundary(c);
}

/****below are generated by chatgpt or some wrapper functions, just take it brief look; they're not very important****/

double ***allocate3DArray()
{
    double ***arr = new double **[L];
    for (int i = 0; i < L; i++)
    {
        arr[i] = new double *[W];
        for (int j = 0; j < W; j++)
        {
            arr[i][j] = new double[H];
        }
    }
    return arr;
}

void deallocate3DArray(double ***arr)
{
    for (int i = 0; i < L; ++i)
    {
        for (int j = 0; j < W; ++j)
        {
            delete[] arr[i][j];
        }
        delete[] arr[i];
    }
    delete[] arr;
}

void process3DArray(double ***arr)
{
    // Access and manipulate the elements of the 3D array
    for (int k = 0; k < H; k++)
    {
        for (int j = 0; j < W; j++)
        {
            for (int i = 0; i < L; i++)
            {
                std::cout << arr[i][j][k] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}

void computeC(double ***c, double dx)
{
    for (int i = 0; i < L; i++)
        for (int j = 0; j < W; j++)
            for (int k = 0; k < H; k++)
                c[i][j][k] = initializeC(i, j, k, dx);
}

void computeUx(double ***ux, double dx)
{
    for (int i = 0; i < L; i++)
        for (int j = 0; j < W; j++)
            for (int k = 0; k < H; k++)
                ux[i][j][k] = initializeUx(i, j, k, dx);
}

void computeUy(double ***uy, double dx)
{
    for (int i = 0; i < L; i++)
        for (int j = 0; j < W; j++)
            for (int k = 0; k < H; k++)
                uy[i][j][k] = initializeUy(i, j, k, dx);
}

void computeUz(double ***uz, double dx)
{
    for (int i = 0; i < L; i++)
        for (int j = 0; j < W; j++)
            for (int k = 0; k < H; k++)
                uz[i][j][k] = initializeUz(i, j, k, dx);
}





double updateC_Carson(int i,int j, int dt, double dx, double k, c,u,D){

  cXleft = c[i-1][j][k];
  uXleft = u[i-1][j][k];
  cXright = c[i+1][j][k];
  uXright = u[i+1][j][k];

  cYleft = c[i][j-1][k];
  uYleft = u[i][j-1][k];
  cYright = c[i][j+1][k];
  uYright = u[i][j+1][k];

  cZleft = c[i][j][k-1];
  uZleft = u[i][j][k-1];
  cZright = c[i][j][k+1];
  uZright = u[i][j][k+1];

  cCenter = getC(i,j,k);

  xTotalFlux  = ((dt)/(2*dx))*(uXright*cXright - uXleft*cXleft) + ((D*dt)/(dx**2))*(2*cCenter - cXleft - cXright); 
  yTotalFlux  = ((dt)/(2*dx))*(uLright*cLright - uLleft*cLleft) + ((D*dt)/(dx**2))*(2*cCenter - cLleft - cLright); 
  zTotalFlux  = ((dt)/(2*dx))*(uZright*cZright - uZleft*cZleft) + ((D*dt)/(dx**2))*(2*cCenter - cZleft - cZright); 

  return cCenter - xTotalFlux - yTotalFlux - zTotalFlux;

}
