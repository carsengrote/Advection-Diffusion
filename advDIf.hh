#ifndef ADVDIF_HH
#define ADVDIF_HH

#include <array>

// i,j,k are x,y,z direction index
// L,W,H are x,y,z direction array length

double initializeC(int i, int j, int k); // return VALUE
double initializeUx(int i, int j, int k, double dx);
double initializeUy(int i, int j, int k, double dx);
double initializeUz(int i, int j, int k, double dx);

bool isBoundary(int i, int j, int k, int L, int W, int H); // given index and bound

double getC(int i, int j, int k);
std::array<double, 3> getU(int i, int j, int k); // return ux, uy, uz at one point

double updateC(int i, int j, int k, double dx, double dt);

#endif
