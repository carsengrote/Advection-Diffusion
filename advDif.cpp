// g++ -std=c++11 -Wall -Wextra -Werror advDif.cpp -o advDif

#include "advDif.hh"
#include <iostream>
#include <cmath>

// number of points of each direction; L for x, W for y, and H for z
int L;
int W;
int H;
// Ub for constant term for u function
double Ub = 1; // cm per second
// D for diffusion term
double D = 0.0000185; // cm squared per second

// Used to print out a vertical cross section of the conentration
void printVertical(double *** c){
   int  i = int(L/2);
   for (int k = 1; k < H; k++){
       for (int j = 1; j < W-1; j++){
           printf("%f ", c[i][j][k]);
       }
       printf("\n");
   }
   printf("\n");

}

int main(int argc,char* argv[])
{
    if (argc != 6){
      exit(0);
    }
    double LL, WW, HH, dx, T;
    //std::cout << "Enter length, width, height, and dx, and total time: ";
    //std::cin >> LL >> WW >> HH >> dx >> T;
    LL = atof(argv[1]);
    WW = atof(argv[2]);
    HH = atof(argv[3]);
    dx = atof(argv[4]);
    T = atof(argv[5]);
    if (LL != WW)
    {
        WW = LL; // We do this because we want to keep xy-plane symmetry; if we have some other u, we can delete this
    }
    start(LL, WW, HH, dx, T);
    return 0;
}

void start(double LL, double WW, double HH, double dx, double T)
{

    L = 2*int(LL / dx) + 2; // L is number of points, + 2 points for ghost node on each side
    W = 2*int(WW / dx) + 2;
    H = int(HH / dx) + 2;
    // test input:
    // std::cout << "L is: " << L << "W is: " << W << "H is: " << H << "dx is: " << dx << std::endl;

    double ***c = allocate3DArray();
    // We need two arrays for keeping track of the concentration because we
    // can't edit it while we're calculating the concentration at the next time
    // step
    double ***cPrime = allocate3DArray();
    double ***cNext = allocate3DArray();
    double ***cNextPrime = allocate3DArray();

    computeC(c, dx);
    enforceBoundary(c);
    // test c initialization:
    //process3DArray(c);

    double ***ux = allocate3DArray();
    double ***uy = allocate3DArray();
    double ***uz = allocate3DArray();
    computeUx(ux, dx, LL, HH);
    computeUy(uy, dx, LL, HH);
    computeUz(uz, dx, LL, HH);
    
    // test u initialization:
    //process3DArray(ux);
    //process3DArray(uy);
    //process3DArray(uz);
    double dt = .5*CFL(dx,ux,uy,uz);
    int images = std::floor(T / dt)+1;
    fprintf(stderr,"Final time: %f, dt: %f, Images: %d, Width pixels: %d, Height pixels: %d\n",T,dt,images,L-2,H-1);
    double t_total = 0;

    printVertical(c); // print initial conditions 
    while (t_total < T)
    {
        c = updateC(c, cPrime, cNext, cNextPrime, ux, uy, uz, dx, dt);
        t_total += dt;

        printVertical(c);
    }
}

double CFL(double dx, double ***ux, double ***uy, double *** uz)
{
    
    // Need to find largest velocity among all the cells
    double maxU = 0;
    double U;
    for (int i = 0; i < L; i++){
        for (int j = 0; j < W; j++){
            for (int k = 0; k < H; k++){
                U = std::sqrt(ux[i][j][k]*ux[i][j][k] +  uy[i][j][k]*uy[i][j][k] + uz[i][j][k]*uz[i][j][k]);
                if (U > maxU){
                    maxU = U;
                }
            }
        }
    }

    //printf("adv CLF: %f , dif CFL: %f\n",dx/maxU,(dx*dx)/(2*D));
    return std::min(dx/maxU, (dx*dx)/(2*D));
}

std::array<double, 3> convert(int i, int j, int k, double dx) // if I only use it once, i can pass value back,
                                                              // and it'll be deleted automatically
{
    std::array<double, 3> coord;
    double x = (double(i) - double(L)*.5) * dx + (.5)*dx;
    double y = (double(j) - double(W)*.5) * dx + (.5)*dx;
    double z = (double(k) - (double(H)- 1)) * dx + (.5)*dx;
    coord = {x, y, z};
    // test convert: it looked good to me: I run it on a 4*4*4 with dx=1 case, and it's ideal
    //std::cout << "right: " << x << ","  << y << "," << z << std::endl;
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
            c[i][j][H - 1] = 0; // Surface concentration 0
            
            for (int k = 0; k < H - 1; k++)
            {
               // Zero flux at boundaries, ghost nodes are at 0 and array
               // length - 1
               if(i == 0){
                  c[i][j][k] = c[i+1][j][k];
               }
               if(i == L - 1){
                  c[i][j][k] = c[i-1][j][k];
               }
               if (j == 0){
                  c[i][j][k] = c[i][j+1][k];
               }
               if (j == W - 1){
                  c[i][j][k] = c[i][j-1][k];
               }
               if (k == 0){
                  c[i][j][k] = c[i][j][k+1];
               }
            }
        }
    }
}

double initializeC(int i, int j, int k, double dx)
{
    // std::array<double, 3> coord = convert(i, j, k, dx); // for future position-dependecy initialization
    // std::cout << coord[0] << " " << coord[1] << " " << coord[2] << std::endl;
    // why do I define a? bc if you dont use ijk, there is compile error; I used the first line of this file to compile
    int a = (i + j + k) * dx * 0;
    return 1+ a; // const for now, grams per liter? 
}

double initializeUx(int i, int j, int k, double dx, double LL, double HH)
{
    std::array<double, 3> coord = convert(i, j, k, dx);
    double x = coord[0];
    double y = coord[1];
    double r = std::sqrt(x * x + y * y);
    double z = coord[2];

    // fuck it, in notes, we put L as height and H is the xy plane dimension; 
    // in the code we flipped them; be extra caustious here
    double ur = Ub * sin(M_PI * r / LL) * cos(M_PI * z / HH);
    double ux;
    if (r < dx / 2) // handle r=0 case
    {
        ux = 0;
    }
    else
    {
        ux = ur * coord[0] / r; // x/r is the component
    }

    if (r > LL){
        return (-.1)*ux;
    }

    return ux;
}

double initializeUy(int i, int j, int k, double dx, double LL, double HH)
{

    std::array<double, 3> coord = convert(i, j, k, dx);
    double x = coord[0];
    double y = coord[1];
    double r = std::sqrt(x * x + y * y);
    double z = coord[2];

    // fuck it, in notes, we put L as height and H is the xy plane dimension;
    // in the code we flipped them; be extra caustious here
    double ur = Ub * sin(M_PI * r / LL) * sin(M_PI * z / HH);
    double uy;
    if (r < dx / 2) // handle r=0 case
    {
        uy = 0;
    }
    else
    {
        uy = ur * coord[1] / r; // y/r is the component
    }

    if (r > LL){
        return -.1*uy;
    }
    return uy;
}

double initializeUz(int i, int j, int k, double dx, double LL, double HH)
{
    std::array<double, 3> coord = convert(i, j, k, dx);
    double r = std::sqrt(std::pow(coord[0], 2) + std::pow(coord[1], 2)); // sqrt(x^2+y^2)
    double z = coord[2];

    // fuck it, in notes, we put L as height and H is the xy plane dimension; 
    // in the code we flipped them; be extra caustious here
    double term1 = -Ub * (HH/LL) * cos(M_PI * r / LL) * sin(M_PI * z / HH);
    double term2 = 0;
    if (r < dx / 2) // handle r=0 case
    {
        term2 = -Ub * HH * sin(M_PI * z / HH); // if r=0, we have a sinc function, so it's 1; 
                                               // it's H/L instead of L/H because of notation definition
    }
    else
    {
        term2 = -Ub * (1 / (M_PI * r)) * (HH * sin(M_PI * r / LL)) * sin(M_PI * z / HH);
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

double updateC_cell(double ***c, double ***ux, double ***uy, double ***uz, int i, int j, int k, double dx)
{
    // xyz direction advection and diffusion
    double x_ad;
    double y_ad;
    double z_ad;
    double x_df;
    double y_df;
    double z_df;

    double total_ad;
    double total_df;

    std::array<double,3> u_east = getU(ux,uy,uz,i+1,j,k);
    double c_east = getC(c,i+1,j,k);

    std::array<double,3> u_west = getU(ux,uy,uz,i-1,j,k);
    double c_west = getC(c,i-1,j,k);

    x_ad = u_east[0]*c_east - u_west[0]*c_west;

    std::array<double,3> u_north = getU(ux,uy,uz,i,j+1,k);
    double c_north = getC(c,i,j+1,k);
    std::array<double,3> u_south = getU(ux,uy,uz,i,j-1,k);
    double c_south = getC(c,i,j-1,k);
    
    y_ad = u_north[1]*c_north - u_south[1]*c_south;

    std::array<double,3> u_top = getU(ux,uy,uz,i,j,k+1);
    double c_top = getC(c,i,j,k+1);
    std::array<double,3> u_bot = getU(ux,uy,uz,i,j,k-1);
    double c_bot = getC(c,i,j,k-1);

    z_ad = u_top[2]*c_top - u_bot[2]*c_bot;
    
    total_ad = (-1 / (2 * dx)) * (x_ad + y_ad + z_ad);
    
    x_df = -getC(c, i - 1, j, k) + 2 * getC(c, i, j, k) - getC(c, i + 1, j, k);

    y_df = -getC(c, i, j - 1, k) + 2 * getC(c, i, j, k) - getC(c, i, j + 1, k);

    z_df = -getC(c, i, j, k - 1) + 2 * getC(c, i, j, k) - getC(c, i, j, k + 1);

    total_df = -D * (1 / (2 * dx * dx)) * (x_df + y_df + z_df);

    return total_ad + total_df;
}

double ***  updateC(double ***c, double ***cPrime, double *** cNext, double *** cNextPrime, double ***ux, double ***uy, double ***uz, double dx, double dt)
{
    for (int i = 1; i < L-1; i++)
    {
        for (int j = 1; j < W-1; j++)
        {
            for (int k = 1; k < H-1; k++)
            {            
                cPrime[i][j][k] = updateC_cell(c, ux, uy, uz, i, j, k, dx);
                cNext[i][j][k] = c[i][j][k] + dt*cPrime[i][j][k];
                if (cNext[i][j][k] < 0){
                  cNext[i][j][k] = 0;
                }
            }
        }
    }

    enforceBoundary(cNext);
    for (int i = 1; i < L-1; i++)
    {
        for (int j = 1; j < W-1; j++)
        {
            for (int k = 1; k < H-1; k++)
            {
                cNextPrime[i][j][k] = updateC_cell(cNext, ux, uy, uz, i, j, k, dx);
                cNext[i][j][k] = cNext[i][j][k] + (dt)*(1.5)*cNextPrime[i][j][k] - (dt)*(.5)*cPrime[i][j][k];
                if (cNext[i][j][k] < 0){
                  cNext[i][j][k] = 0;
                } 
            }
        }
    }

    enforceBoundary(cNext);
    return cNext;
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

void computeUx(double ***ux, double dx, double LL, double HH)
{
    for (int i = 0; i < L; i++)
        for (int j = 0; j < W; j++)
            for (int k = 0; k < H; k++)
                if(isBoundary(i,j,k)){
                    ux[i][j][k] = 0;
                }else{
                    ux[i][j][k] = initializeUx(i, j, k, dx, LL, HH);
                }
}

void computeUy(double ***uy, double dx, double LL, double HH)
{
    for (int i = 0; i < L; i++)
        for (int j = 0; j < W; j++)
            for (int k = 0; k < H; k++)
                if(isBoundary(i,j,k)){
                    uy[i][j][k] = 0;
                }else{
                    uy[i][j][k] = initializeUy(i, j, k, dx, LL, HH);
                }
}

void computeUz(double ***uz, double dx, double LL, double HH)
{
    for (int i = 0; i < L; i++)
        for (int j = 0; j < W; j++)
            for (int k = 0; k < H; k++)
                if(isBoundary(i,j,k)){
                    uz[i][j][k] = 0;
                }else{
                    uz[i][j][k] = initializeUz(i, j, k, dx, LL, HH);
                    if (k == 1 && uz[i][j][k] < 0){
                        uz[i][j][k] = 0;
                    }
                }  
}
