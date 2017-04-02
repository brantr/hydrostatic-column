#include <stdio.h>
#include <stdlib.h>
#include "hydrostatic_column.h"

int main(int argc, char **argv)
{
  double *rho;  //density array in z-direction
  double x;     //column location in x-direction
  double y;     //column location in y-direction
  int nz;       //number of cells along z-direction
  int ng;       //number of ghost cells
  double dz;    //cell width along z-direction

  double L_z  = 10.0; //real box is 10 kiloparsecs in z-direction

  x = 0;  // x=0
  y = 0;  // y=0

  nz = 128; //start with 128^3
  ng = 4;   //number of ghost cells

  if(argc>1)
    nz = atoi(argv[1]);

  //set cell size in z-direction
  dz = L_z / ((double) nz);

  //set the disk properties
  SetDiskProperties(&Disk);

  //get hydrostatic column
  rho = hydrostatic_column(x, y, dz, nz, ng);

  return 0;

}