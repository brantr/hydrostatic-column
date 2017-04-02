#ifndef  HYDROSTATIC_COLUMN_H
#define  HYDROSTATIC_COLUMN_H


//initial guess for disk properties
struct DiskProperties
{
  double G;       //gravitational constant in kpc^3 / Msun / ky^2
  double M_d;     //mass in Msun
  double R_g;     //gas scale length in kpc
  double Sigma_0; //central surface density in Msun/pc^2
  double Z_d;     //vertical scale of disk in kpc
  double H_g;     //guess of vertical scale height in kpc
  double T_eos;   //disk temperature at rho_eos
  double phi_0_d; //in (kpc/kyr)^2
  double B_d;     //in (kpc/kyr)^2 / cs^2
  double gamma;   //adiabatic index
  double cs;      //sound speed for gas at T_8
  double rho_eos; //density where eos is normalized (1e8 Msun/kpc^3)
  double K_eos;   //P = K_eos * rho^gamma

};

//global structure for disk properties
extern struct DiskProperties Disk;

void SetGasDensities(double *rho, double *gz, double r, double dz, int nz, int ng, struct DiskProperties *disk);

//equation of state
double P_eos(double rho, struct DiskProperties disk);

//disk radial surface density profile
double Sigma_disk(double r, struct DiskProperties disk);

//subroutine to set the initial guess for disk properties
void SetDiskProperties(struct DiskProperties *disk);

//returns an array containing a 
//column with the hydrostatic density profile
//at locations (x,y,*)
double *hydrostatic_column(double x, double y, double dz, int nz, int ng);

//returns the cell-centered vertical
//location of the cell with index k
double z_hc(int k, double dz, int nz, int ng);

#endif //HYDROSTATIC_COLUMN_H
