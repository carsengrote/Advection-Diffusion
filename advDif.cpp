double updateC(int i,int j, int dt, double dx, double k, c,u,D){

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
