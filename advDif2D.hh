#ifndef ADVDIF2D_HH
#define ADVDIF2d_HH

#include <array>
// i,k are x,z direction index
// L,H are x,z direction array length

void start(double LL, double HH, double dx, double T);
std::array<double, 2> convert(int i, int k, double dx);

double initializeC(int i, int k); // return VALUE
double initializeUx(int i, int k, double LL, double HH);
double initializeUy(int i, int k, double LL, double HH);
double initializeUz(int i, int k, double LL, double HH);

double CFL(double dx, double **ux, double **uz);

bool isBoundary(int i, int k); // given index and bound
void enforceBoundary(double **c);

double getC(double **c, int i, int k);
std::array<double, 2> getU(double **ux, double **uz, int i, int k); // return ux, uy, uz at one point

struct ptrStruct updateC(double **c, double **cPrime, double **cPrimeLast, double **cNext, double **ux, double **uz, double dx, double dt, int E);

double **allocate2DArray();
void deallocate2DArray(double **arr);
void process2DArray(double **arr);
void computeC(double **c, double dx);
void computeUx(double **ux, double dx, double LL, double HH);
void computeUz(double **uz, double dx, double LL, double HH);

double totalC(double **c);
double avgC(double **c);
void printVertical(double **c, double dx);
void processU(double **ux, double **uz, double dx);

#endif
