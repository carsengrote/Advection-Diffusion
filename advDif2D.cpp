// Need to use actual GNU gcc to compile this code as Apple's
// alias to gcc, clang, does not support OpenMP
// g++-13 -fopenmp -std=c++11 -Wall -Wextra advDif2D.cpp -o advDif

#include "advDif2D.hh"
#include <iostream>
#include <cmath>
#include <omp.h>

struct ptrStruct{
  double **newC;
  double **lastCprime;
};

// number of points of each direction; L for x and H for z
int L;
int H;

double N; // Number of circulation cells
// Ub for constant term for u function
double Ub = 1; // cm per second
// D for diffusion term
double D = 0.00925; // cm squared per second

// Used to print out a vertical cross section of the conentration
// for making videos. Is a bit weirdly formatted to work with Gnuplot
// splot and is used with the videoScriptThreaded.py script
void printVertical(double **c, double dx){
    printf("%d ", L-2);
    std::array<double,2> coord = convert(1,0,dx);
    double x = coord[0];
    printf("%f ", x - .5*dx);
    x = x + dx + .5*dx;
    for (int i = 2; i < L - 2; i++){
        printf("%f ",x);
        x = x + dx;
    }
    x = x + .5*dx;
    printf("%f \n", x);
    
    coord = convert(0, 1, dx);
    double y = coord[1];
    for (int k = 1; k < H; k++){
        printf("%f ",y);
        y = y + dx;
        for (int i = 1; i < L-1; i++){
            printf("%f ", c[i][k]);
        }
        printf("\n");
   }
   printf("\n");
}

// Calculate the average concentration over 
// all of the cells
double avgC(double **c){
    double sum = 0;
    for (int i = 1; i < L-1; i++){
        for (int k = 1; k < H -1; k++){
            sum = sum + c[i][k];
        }
    }
    return sum/((L-2)*(H-2));
}

// Sum up the concentrations in every cell
double totalC(double **c){
    double sum = 0;
    for (int i = 1; i < L -1; i++){
        for (int k = 1; k < H - 1; k++){
            sum = sum + c[i][k];
        }
    }
    return sum;
}

int main(int argc,char* argv[])
{
    // Can change depending on computer this is
    // running on. Apple says I have 10 cores.
    omp_set_num_threads(10);

    // Check CLAs
    if (argc != 6){
      fprintf(stderr, "Usage: ./advDif2D <Container length> <Container height> <Number of circulation cells> <Grid spacing> <Final time>");
      exit(0);
    }
    double LL, HH, dx, T;
    // Getting CLAs 
    LL = atof(argv[1]); // Length of container
    HH = atof(argv[2]); // Height of container
    N = atof(argv[3]); // Number of circulation cells
    dx = atof(argv[4]); // Grid spacing
    T = atof(argv[5]); // Final time
    
    start(LL, HH, dx, T);
    return 0;
}

void start(double LL, double HH, double dx, double T)
{

    L = int(LL / dx) + 2; // L is number of points, + 2 points for ghost node on each side
    H = int(HH / dx) + 2;

    double **c = allocate2DArray();
    // We need multiple arrays for keeping track of the concentration because we
    // can't edit it while we're calculating the concentration at the next time
    // step, and AB-2 is a multistep method
    double **cPrime = allocate2DArray();
    double **cNext = allocate2DArray();
    double **cPrimeLast = allocate2DArray();

    computeC(c, dx);
    enforceBoundary(c);
    // test c initialization:
    //process3DArray(c);

    // initializing the velocity field
    double **ux = allocate2DArray();
    double **uz = allocate2DArray();
    computeUx(ux, dx, LL, HH);
    computeUz(uz, dx, LL, HH);
    // test u initialization:
    //process2DArray(ux);
    //process2DArray(uz);

    double dt = CFL(dx,ux,uz);
    int images = std::floor(T / dt)+1;
    fprintf(stderr,"Final time: %f, dt: %f, Images: %d, Width pixels: %d, Height pixels: %d\n",T,dt,images,L-2,H-1);
    
    double t_total = 0.0;
    struct ptrStruct ptrs;
    printVertical(c,dx); // print initial conditions 
    //printf("%f %f %f\n",t_total,avgC(c),totalC(c));
    int E = 0; // Switch for Eulers step first
    while (t_total < T)
    {
        if (t_total == 0){
            E = 1;
        }else{
            E = 0;
        }
        // Time stepping
        ptrs= updateC(c, cPrime, cPrimeLast, cNext, ux, uz, dx, dt, E);
        c = ptrs.newC;
        cPrimeLast = ptrs.lastCprime;
        t_total += dt;
        //printf("%f %f %f\n",t_total,avgC(c), totalC(c));  
        printVertical(c,dx);
    }
    
}

double CFL(double dx, double **ux, double ** uz)
{
    
    // Need to find largest flow velocity among all the cells
    double maxU = 0;
    double U;
    for (int i = 0; i < L; i++){
        for (int k = 0; k < H; k++){
            U = std::sqrt(ux[i][k]*ux[i][k] + uz[i][k]*uz[i][k]);
            if (U > maxU){
                maxU = U;
            }
        }
    }
    
    //printf("%f, %f \n", dx/maxU, (dx*dx)/D);
    // Returning the smaller CFL condition between advection and diffusion
    // For small diffusion, advection CFL is normally much smaller
    return .75*std::min(dx/maxU, (dx*dx)/(2*D));
}

// Given indicies, returns physical x and z coordinates
std::array<double, 2> convert(int i, int k, double dx)
{
    std::array<double, 2> coord;
    double x = (double(i) - double(L)*.5) * dx + (.5)*dx;
    double z = (double(k) - (double(H)- 1)) * dx + (.5)*dx;
    coord = {x, z};
    return coord;
}

// Used to check if a cell is a ghost node or not
bool isBoundary(int i, int k)
{
    return (i == 0 || i == L - 1 || k == 0 || k == H - 1);
}

// At every time step need to keep ghost nodes having same concentration.
// Keep interior nodes on the layer closest to the surface at c = 0
void enforceBoundary(double **c)
{   
    #pragma omp parallel for collapse(2)
    for(int i = 0; i < L; i++){
        for (int k = 0; k < H; k++){
            // Ghost nodes on sides of the container
            if (i == 0){ 
                c[i][k] = c[i+1][k];
            }else if (i == L-1){
                c[i][k] = c[i-1][k];
            }   
            
            // Ghost node at the top and set boundary
            // layer to zero concentration
            if (k == 0){ 
                c[i][k] = c[i][k+1];
            }else if (k == H-1 || k == H -2){
                c[i][k] = 0;
            }   
        }   
    } 
}

double initializeC(int i, int k, double dx)
{ 
    // Random number in [-1,1]
    double r = (((double) rand() / (RAND_MAX))*2) - 1;
    r = 0; // Comment out this line for random initial concentration
    
    std::array<double,2> coord = convert(i,k,dx);
    return 1  + .5*r  + 0*(-.2)*coord[1];
}

// x component of the fluid velocity at grid point (x_i,z_j)
double initializeUx(int i, int k, double dx, double LL, double HH)
{
    std::array<double,2> coord = convert(i, k, dx);
    double x = coord[0];
    double z = coord[1];
    double ux;
    // Equation (2) in the report 
    ux =  - Ub * LL * (1/(2*HH*N)) * sin((2*M_PI*N)*((x+LL/2)/LL)) * cos((M_PI*z)/HH);
    return ux;
}

// z component of the fluid velocity at grid point (x_i, z_j)
double initializeUz(int i, int k, double dx, double LL, double HH)
{
    std::array<double, 2> coord = convert(i, k, dx);
    double x = coord[0];
    double z = coord[1];
    // Equation (3) in the report
    double uz = Ub * cos((2*M_PI*N)*((x+LL/2)/LL)) * sin((z*M_PI)/HH);
    return uz;
}

double getC(double **c, int i, int k)
{
    return c[i][k];
}
// Returns a vector representing the x and z components of velocity 
// at that point in the domain
std::array<double, 2> getU(double **ux, double **uz, int i, int k)
{
    std::array<double, 2> u;
    double x = ux[i][k];
    double z = uz[i][k];
    u = {x,z};
    return u;
}

double updateC_cell(double **c, double **ux, double **uz, int i, int k, double dx)
{
    // x,z direction advection and diffusion
    double x_ad;
    double z_ad;
    double x_df;
    double z_df;

    double total_ad;
    double total_df;

    std::array<double,2> u_east = getU(ux,uz,i+1,k);
    double c_east = getC(c,i+1,k);
    std::array<double,2> u_west = getU(ux,uz,i-1,k);
    double c_west = getC(c,i-1,k);

    // Flux due to advection in the x direction
    x_ad = u_east[0]*c_east - u_west[0]*c_west;
    
    std::array<double,2> u_top = getU(ux,uz,i,k+1);
    double c_top = getC(c,i,k+1);
    std::array<double,2> u_bot = getU(ux,uz,i,k-1);
    double c_bot = getC(c,i,k-1);
    
    // Flux due to advection in the z direction
    z_ad = u_top[1]*c_top - u_bot[1]*c_bot; 

    total_ad = (-1 / (2 * dx)) * (x_ad + z_ad);
    
    // Diffusion in the x and z directions
    x_df = -getC(c, i - 1, k) + 2 * getC(c, i, k) - getC(c, i + 1, k);
    z_df = -getC(c, i, k - 1) + 2 * getC(c, i, k) - getC(c, i, k + 1);

    total_df = -D * (1 / (dx * dx)) * (x_df + z_df);

    return total_ad + total_df;
}

struct ptrStruct updateC(double **c, double **cPrime, double **cPrimeLast, double **cNext, double **ux, double **uz, double dx, double dt, int Euler)
{
    
    double d1,d2;
    // Modify method to be Euler's step if it's the first time step of the
    // simulation, otherwise will use an AB2 step
    if (Euler == 1){
        d1 = 1;
        d2 = 0;
    }else{
        d1 = 1.5;
        d2 = 0.5;
    }
    // Adam's-Bashforth two step for O(k^2) time accuracy. 
    // y_{n+2} = y_{n+1} + (3/2)dx*f(y_{n+1}) - (1/2)dx*f(y_{n}) 
    #pragma omp parallel for num_threads(10) collapse(2)
    for (int i = 1; i < L-1; i++)
    {
        for (int k = 1; k < H-1; k++)
        {            
            cPrime[i][k] = updateC_cell(c, ux, uz, i, k, dx);
            cNext[i][k] = c[i][k] + dt*(d1)*cPrime[i][k] - (dt)*(d2)*cPrimeLast[i][k];
            if (cNext[i][k] < 0){
                cNext[i][k] = 0;
            }
        }
    }
   
    // Enforce boundary conditions after each time step 
    enforceBoundary(cNext);
    
    // Organize new pointers to return
    struct ptrStruct updatedPtrs; 
    updatedPtrs.newC = cNext;
    updatedPtrs.lastCprime = cPrime;
    
    return updatedPtrs;
}

// Used to allocate a 2D array on the heap
double **allocate2DArray()
{
    double **arr = new double *[L];
    for (int i = 0; i < L; i++)
    {
        arr[i] = new double [H];
    }
    return arr;
}

// Used to print out a 2D array for debugging
void process2DArray(double **arr)
{
    for (int k = H-1; k > -1; k--)
    {
        for (int i = 0; i < L; i++)
        {
            std::cout << arr[i][k] << " ";
        }
        printf("\n");
    }
    std::cout << std::endl;
}

void computeC(double **c, double dx)
{
    for (int i = 0; i < L; i++)
        for (int k = 0; k < H; k++)
            c[i][k] = initializeC(i, k, dx);
}

void computeUx(double **ux, double dx, double LL, double HH)
{
    for (int i = 0; i < L; i++)
        for (int k = 0; k < H; k++)
            if(isBoundary(i,k)){
                ux[i][k] = 0;
            }else if(i == 1 || i == L - 2){
                ux[i][k] = 0;
            }else{
                ux[i][k] = initializeUx(i, k, dx, LL, HH);
            }
}

void computeUz(double **uz, double dx, double LL, double HH)
{
    for (int i = 0; i < L; i++)
        for (int k = 0; k < H; k++)
            if(isBoundary(i,k)){
                uz[i][k] = 0;
            }else{
                uz[i][k] = initializeUz(i, k, dx, LL, HH);
                if (k == 1 && uz[i][k] < 0){
                    uz[i][k] = 0;
                }
            }  
}
