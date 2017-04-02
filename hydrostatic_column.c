#include "hydrostatic_column.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define sgn(v) ( ( (v) < 0 ) ? -1 : ( (v) > 0 ) )

//global structure for disk properties
struct DiskProperties Disk;

//disk radial surface density profile
double Sigma_disk(double r, struct DiskProperties disk)
{
  //return the exponential surface density
  return disk.Sigma_0 * exp(-r/disk.R_g);
}

//equation of state
double P_eos(double rho, struct DiskProperties disk)
{
  double x  = rho/disk.rho_eos; //density where temperature is normalized
  double cs = disk.cs;          //effective soundspeed at normalized density
  double gamma = disk.gamma;    //eos adiabatic index
  //return pressure
  //return cs*cs*disk.rho_eos*pow(x,gamma)/gamma;
  return disk.K_eos*pow(rho,gamma);
}

//vertical acceleration in miyamoto nagai
double gz_disk(double R, double z, struct DiskProperties disk)
{
  double a = disk.R_g;
  double b = disk.Z_d;
  double A = sqrt(b*b + z*z);
  double B = a + A;
  double C = pow(B*B + R*R, 1.5);
  double G = disk.G;
  double M = disk.M_d;

  //checked with wolfram alpha
  return -G*M*z*B/(A*C);
}

//radial acceleration in miyamoto nagai
double gr_disk(double R, double z, struct DiskProperties disk)
{
  double a = disk.R_g;
  double b = disk.Z_d;
  double A = sqrt(b*b + z*z);
  double B = a + A;
  double C = pow(B*B + R*R, 1.5);
  double G = disk.G;
  double M = disk.M_d;

  //checked with wolfram alpha
  return -G*M*R/(A*C);
}

//exponential integral
double exponential_density_profile_integral(double z_min, double z_max, struct DiskProperties disk)
{
  //int_a^b exp(-z/h) dz = h*(exp(-a/h) - exp(-b/h))
  double h = disk.H_g;
  return h*( exp(-1.*z_min/h) - exp(-1.*z_max/h) );
}

//returns an array containing a 
//column with the hydrostatic density profile
//at locations (x,y,*)
double *hydrostatic_column(double x, double y, double dz, int nz, int ng)
{
  //x is cell center in x direction
  //y is cell center in y direction
  //dz is cell width in z direction
  //nz is number of real cells
  //ng is number of ghost cells
  //total number of cells in column is nz * 2*ng

  int k;        //index along z axis
  double *rho;  //density array in column
  double *dPdz; //pressure gradient in column
  double *gz;   //vertical acceleration in column
  double *gr;   //radial acceleration in column
  double *drho; //change to density

  double z;     //cell center in z direction
  int nzt;      //total number of cells in z-direction
  double Sigma; //surface density in column
  double z_min; //bottom of the cell
  double z_max; //top of the cell
  double r;       //cylindrical radius at xy
  double Sigma_r; //surface density expected at r
  double P_a, P_b, Delta_z; //pressure gradient computation
  double gamma = Disk.gamma;
  double K = Disk.K_eos;
  double rho_floor = 1.0e-2;

  int iter = 0; //number if iterations
  int ks; //start of integrals above disk plane
  if(nz%2)
  {
    ks = ng+(nz-1)/2;
  }else{
    ks = ng + nz/2;
  }

  printf("In hydrostatic_column.\n");

  //set the cylindrical radius
  r = sqrt(x*x + y*y);

  //get the disk surface density
  Sigma_r = Sigma_disk(r, Disk);
  //printf("Surface density at r = %e is Sigma = %e\n",r,Sigma_r);

  //set the z-column size, including ghost cells
  nzt = nz + 2*ng;

  //allocate rho
  if(!(rho=(double *) calloc(nzt,sizeof(double))))
  {
    printf("Error allocating rho array of size %d.\n",nzt);
    printf("Aborting...\n");
    fflush(stdout);
    exit(-1);
  }

  //allocate pressure gradient
  if(!(dPdz=(double *) calloc(nzt,sizeof(double))))
  {
    printf("Error allocating dPdz array of size %d.\n",nzt);
    printf("Aborting...\n");
    fflush(stdout);
    exit(-1);
  }

  //allocate vertical acceleration 
  if(!(gz=(double *) calloc(nzt,sizeof(double))))
  {
    printf("Error allocating gz array of size %d.\n",nzt);
    printf("Aborting...\n");
    fflush(stdout);
    exit(-1);
  }

  //allocate radial acceleration 
  if(!(gr=(double *) calloc(nzt,sizeof(double))))
  {
    printf("Error allocating gz array of size %d.\n",nzt);
    printf("Aborting...\n");
    fflush(stdout);
    exit(-1);
  }

  //allocate density corrections 
  if(!(drho=(double *) calloc(nzt,sizeof(double))))
  {
    printf("Error allocating drho array of size %d.\n",nzt);
    printf("Aborting...\n");
    fflush(stdout);
    exit(-1);
  }


  //check z positions
  //for(k=0;k<nzt;k++)
  //{
  //  z = z_hc(k,dz,nz,ng);
  //  printf("k %d cl %e cc %e cu %e\n",k-ng,z-0.5*dz,z,z+0.5*dz);
  //}

  //compute vertical and radial
  //gravitational accelerations
  for(k=0;k<nzt;k++)
  {
    z     = z_hc(k,dz,nz,ng);
    gz[k] = gz_disk(r,z,Disk);
    gr[k] = gr_disk(r,z,Disk);
  }

  //set densities
  //SetGasDensities(rho, gz, r, dz, nz, ng, &Disk);

  //set initial guess for disk properties
  //assume the disk is an exponential vertically to start
  Sigma = 0;
  for(k=0;k<nzt;k++)
  {
    z_min  = z_hc(k,dz,nz,ng) - 0.5*dz;
    z_max  = z_hc(k,dz,nz,ng) + 0.5*dz;
    if(z_max>0)
    {
      if(z_min<0)
      {
        //in disk plane centered at z=0
        rho[k] = 2.*exponential_density_profile_integral(0, z_max, Disk);
      }else{
        //above disk plane
        rho[k] = exponential_density_profile_integral(z_min, z_max, Disk);        
      }
    }else{

      //below disk plane
      rho[k] = exponential_density_profile_integral(fabs(z_max), fabs(z_min), Disk);
    }
    Sigma += rho[k];
  }

  //renormalize density to match surface density
  for(k=0;k<nzt;k++)
  {
    rho[k] *= Sigma_r/(Sigma*dz);

    //if(k>=ks)
    //  printf("%e\t%e\n",z_hc(k,dz,nz,ng),rho[k]);
  }
  //check
  //Sigma =0;
  //for(k=0;k<nzt;k++)
  //{
  //  Sigma += rho[k]*dz;
  //}
  //printf("Sigma %e Sigma_r %e\n",Sigma,Sigma_r);


  //OK, rho is set to an exponential
  //let's adjust to make it hydrostatic

  //begin iterative process to set the density
  int flag = 1;
  double rho_new;
  double mass_loss;
  while(flag)
  {
    z = z_hc(k,dz,nz,ng);

    //adjust density (with zeros on first iteration)
    Sigma = 0;
    for(k=0;k<nzt;k++)
    {
      Sigma  += rho[k]*dz;
    }
    //printf("Sigma %e Sigma_r %e\n",Sigma,Sigma_r);
    //for(k=0;k<nzt;k++)
    //  rho[k] *= Sigma_r/Sigma;



    mass_loss = 0;
    for(k=ks;k<nzt-1;k++)
    {
      //z position
      z     = z_hc(k,dz,nz,ng);
      drho[k] = -1.*pow(rho[k],2-gamma)*(fabs(gz[k])*dz/(gamma*K));

      //printf("z %e rho %e drho %e\n",z,rho[k],drho[k]);
      if(drho[k]<-0.9*rho[k])
        drho[k] = -0.9*rho[k];
      rho_new = rho[k]+drho[k];
      mass_loss += (rho[k+1]-rho_new);
      rho[k+1] = rho_new;
      if(rho[k+1]<rho_floor)
        rho[k+1] = rho_floor;

      //compute grad P
      P_b = P_eos(rho[k+1],Disk);
      P_a = P_eos(rho[k],Disk);
      Delta_z = 1.0*dz;     

      /*
      if( (k!=0) && (k!=(nzt-1)) )
      {
        P_b = P_eos(rho[k+1],Disk);
        P_a = P_eos(rho[k-1],Disk);
        Delta_z = 2.0*dz;
      }else{
        if(k==0)
        {
          P_b = P_eos(rho[k],Disk);
          P_a = P_eos(rho[k+1],Disk);
          Delta_z = 1.0*dz;
        }else{
          P_b = P_eos(rho[k],Disk);
          P_a = P_eos(rho[k-1],Disk);
          Delta_z = 1.0*dz;
        }
      }*/
      //based on equation of state
      dPdz[k] = -1.0*sgn(z)*fabs(P_b-P_a)/Delta_z;



      //printf("dPdz %e test %e mass_loss %e\n",dPdz[k],gamma*K*pow(rho[k],gamma-1)*drho[k]/Delta_z,mass_loss);

      //if drho<0, rho*g is bigger than dPdz
      //if drho>0, rho*g is smaller than dPdz


      //if(z>0)
      //  printf("z %e\trho %e\tdPdz % e\t-rho g % e gz % e drho %e\n",z,rho[k],dPdz[k],-rho[k]*gz[k],gz[k],drho[k]);
    }
    //printf("mass_loss = %e\n",mass_loss*dz/Sigma);
    int km;
    for(k=ks;k<nzt;k++)
    {
      if(mass_loss<0)
      {
        //rho[k] *= (1 + fabs(mass_loss)/Sigma);
        rho[k] += mass_loss/((float) (nzt-1-ks+1));
      }else{
        rho[k] -= mass_loss/((float) (nzt-1-ks+1));
        //rho[k] *= (1 - fabs(mass_loss)/Sigma);
      }
      //mirror densities
      if(nz%2)
      {
        km = (ng+(nz-1)/2) - (k-ks);
      }else{
        km = ng + nz/2 - (k-ks) -1;
      }
      rho[km] = rho[k];
    }





    // dP/dz + rho g = A
    // P = K rho^gamma
    // dP/dz = gamma * K * rho^(gamma-1)  drho/dz
    // A = rho*g + gamma * K * rho^(gamma-1) drho/dz
    // A = rho( G + K * rho^(gamma-2) drho/dz)

    //printf("*****\n");
    iter++;
    //if(iter>50)
    if(fabs(mass_loss*dz)/Sigma<1.0e-3)
      flag=0;

    if(iter>50)
      printf("Error converging with iter = %d, mass_loss = %e\n",iter,mass_loss);
  }
  

  //free ancillary arrays
  free(dPdz);
  free(gz);
  free(gr);

  //return the rho array
  return rho;
}


//returns the cell-centered vertical
//location of the cell with index k
//k is indexed at 0 at the lowest ghost cell
double z_hc(int k, double dz, int nz, int ng)
{
  //checked that this works, such that the
  //if dz = L_z/nz for the real domain, then the z positions
  //are set correctly for cell centers with nz spanning
  //the real domain, and nz + 2*ng spanning the real + ghost domains
  if(!(nz%2))
  {
    //even # of cells
    return 0.5*dz + ((double) (k-ng-nz/2))*dz;
  }else{
    //odd # of cells
    return ((double) (k-ng-(nz-1)/2))*dz;
  }
}



//subroutine to set the initial guess for disk properties
void SetDiskProperties(struct DiskProperties *disk)
{
  //some constants
  double l_s = 3.086e21;                //length scale, centimeters in a kiloparsec
  double m_s = 1.99e33;                 //mass scale, g in a solar mass
  double t_s = 3.154e10;                //time scale, seconds in a kyr
  double d_s = m_s / pow(l_s,3);        //density scale, M_sun / kpc^3
  double v_s = l_s / t_s;               //velocity scale, kpc / kyr
  double p_s = d_s*v_s*v_s;             //pressure scale, M_sun / kpc kyr^2
  double G = 6.67259e-8;                //in cm^3 g^-1 s^-2
  double mp = 1.67e-24;                 //proton mass in grams
  G = G / pow(l_s,3) * m_s * t_s*t_s;   //in kpc^3 / M_sun / kyr^2
  double KB = 1.3806e-16;               //boltzmann constant in cm^2 g / s^2 K
  double M_vir = 1.0e12;                //virial mass in solar masses
  double M_d = 6.5e10;                  //disk mass in solar masses
  double M_h = M_vir - M_d;             //halo mass in solar masses
  double R_vir = 261.0;                 //MW viral radius in kpc
  double c_vir = 20.0;                  //MW concentration (seems large!)
  double R_h = R_vir / c_vir;           //halo scale radius in kpc
  double R_d = 3.5;                     //stellar disk scale length in kpc
  double Z_d = 3.5/5.0;                 //disk scale height in kpc
  double R_g = 2*R_d;                   //gas disk scale length in kpc
  double T = 1.0e5;                     //gas temperature, 10^5 K
  double v_to_kmps = l_s/t_s/100000;    //velocity scale to km/s
  double kmps_to_kpcpkyr = 1.0220122e-6;  //velocity scale to kpc/kyr
  double gamma = 5./3.;                   //adiabatic index
  
  disk->G = G;   //Grav constant in kpc^3 / M_sun / kyr^2
  disk->M_d = M_d;   //mass in Msun
  disk->R_g = R_g;   //gas scale length in kpc
  disk->Sigma_0 = 0.25*M_d/(2*M_PI*R_g*R_g); //central surface density in Msun/kpc^2
  disk->Z_d = Z_d;   //vertical scale height in kpc
  disk->T_eos = 1.0e5; //disk temperature at rho_eos
  disk->gamma = gamma;  //adiabatic index


  //isothermal sound speed at Sigma_8
  double cs = sqrt(KB*disk->T_eos/(0.6*mp))*t_s/l_s;  
  disk->cs  = cs; //gas sound speed at T_eos



  //printf("cs = %e\n",cs);

  disk->rho_eos = 1.0e7; //gas eos normalized at 1e8 Msun/kpc^3
  disk->K_eos = cs*cs*pow(disk->rho_eos,1.0-gamma)/gamma;
  //return cs*cs*disk.rho_eos*pow(x,gamma)/gamma;
  //define phi_0_d
  disk->phi_0_d = G * M_d / R_d;
  disk->B_d     = disk->phi_0_d / (cs*cs);

  //guess at an initial scale height
  //disk->H_g = (cs*cs)/(G*disk->Sigma_0);   //gas scale height in kpc
  disk->H_g = Z_d;
  printf("Disk gas scale height (guess) = %e\n",disk->H_g);
}

/*

void SetGasDensities(double *rho, double *gz, double r, double dz, int nz, int ng, struct DiskProperties *disk)
{
  //should have that, given grad P = - rho g
  //and P = K rho^gamma
  //and grad P = gamma K rho^(gamma-1) drho/dz
  //then rho^(gamma-2) drho = -gdz/(K gamma)
  //so (1/(gamma-1))*rho^(gamma-1) + C = int_z^\infty (g dz)/(K gamma)
  //and rho^(gamma-1) = (gamma-1)*(A + int_z^\infty (g dz)/(K gamma))
  //and, finally rho = [(gamma-1)*(A + int_z^\infty (g dz)/(K gamma))]^(1/(gamma-1))
  //where A is set such that \int_0^\infty rho dz = 0.5 Sigma

  //start above disk
  int k;
  int ks;
  int nzt = nz + 2*ng;
  double gmg = (disk->gamma-1)/disk->gamma;
  double gamma = disk->gamma;
  double K = disk->K_eos;
  double Sigma_r;
  double Sigma = 0;
  double A_off;
  double *dPdz;
  double z;
  double P_a, P_b, Delta_z;
  int iter = 0;
  double a_off;
  double z_ks;
  double Error;
  if(nz%2)
  {
    ks = ng+(nz-1)/2;
  }else{
    ks = ng + nz/2;
  }
  z_ks = z_hc(ks,dz,nz,ng);
  printf("check z %e\n",z_ks);
  printf("dz %e K %e gmg %e\n",dz,K,gmg);
  //exit(-1);

  //get the disk surface density
  Sigma_r = Sigma_disk(r, Disk);

  //allocate pressure gradient
  if(!(dPdz=(double *) calloc(nzt,sizeof(double))))
  {
    printf("Error allocating dPdz array of size %d.\n",nzt);
    printf("Aborting...\n");
    fflush(stdout);
    exit(-1);
  }

  A_off = 1.0;
  a_off = 0.0;
  while(1)
  {
    K = disk->K_eos;

    //above disk plane first
    //first we compute rho^(gamma-1)
    rho[nzt-1] = fabs(gz[nzt-1])*dz*gmg/K;
    printf("k %d z %e rho %e\n",nzt-1,z_hc(nzt-1,dz,nz,ng),pow(rho[nzt-1],1./(gamma-1)));
    for(k=nzt-2;k>=ks;k--)
    {
      rho[k] = rho[k+1] + fabs(gz[k+1])*dz*gmg/K; //add contribution from cell above
      
      printf("k %d z %e rho %e\n",nzt-1,z_hc(k,dz,nz,ng),pow(rho[k],1./(gamma-1)));
    }

    Sigma = 0;
    for(k=ks;k<nzt;k++)
    {
      rho[k] = pow(rho[k],1./(gamma-1)) * A_off;
      Sigma += rho[k]*dz; //add contribution to sigma
    }

    //compute pressure gradient
    Error = 0;
    for(k=ks;k<nzt;k++)
    {
      z = z_hc(k,dz,nz,ng);
      //compute grad P
      if( (k!=ks) && (k!=(nzt-1)) )
      {
        P_b = P_eos(rho[k+1],Disk);
        P_a = P_eos(rho[k-1],Disk);
        Delta_z = 2.0*dz;
      }else{
        if(k==ks)
        {
          P_b = P_eos(rho[ks],Disk);
          P_a = P_eos(rho[ks+1],Disk);
          Delta_z = 1.0*dz;
        }else{
          P_b = P_eos(rho[k],Disk);
          P_a = P_eos(rho[k-1],Disk);
          Delta_z = 1.0*dz;
        }
      }
      dPdz[k] = -1.0*sgn(z)*fabs(P_b-P_a)/Delta_z;

      if(z>disk->Z_d)
        Error += (dPdz[k] -rho[k]*gz[k])/fabs(rho[k]*gz[k]);
      printf("z %e\trho %e\tdPdz % e\t-rho g % e gz % e\n",z,rho[k],dPdz[k],-rho[k]*gz[k],gz[k]);

    }

    printf("Error = %e\n",Error);

    //need to correct normalization of rho for surface
    //density
    printf("Sigma %e Sigma_r %e\n",2*Sigma,Sigma_r);

    //try to guess A_off
    //A_off = 0.01*pow(fabs(Sigma_r-2*Sigma)/dz,(gamma-1));
    disk->K_eos *= pow(2.*Sigma/Sigma_r, gamma-1);
    //if(Sigma_r<2*Sigma)
    //{
      disk->K_eos *= 1
    //}else{
    //  disk->K_eos *= 0.9;
    //}
    //printf("A_off = %e\n",disk->K_eos);

    
    k = ks+0.1*nz;
    P_a = fabs(fabs(dPdz[k]) - fabs(rho[k]*gz[k]))/fabs(rho[k]*gz[k]);

    k = ks+0.2*nz;
    P_b = fabs(fabs(dPdz[k]) - fabs(rho[k]*gz[k]))/fabs(rho[k]*gz[k]);
    if(P_a>P_b)
    {
      a_off -= 0.05;
    }else{
      a_off += 0.05;
    }
    printf("a_off %e\n",a_off);
    

    iter++;
    if(iter>1)
      break;

  }


  exit(-1);

  free(dPdz);

}*/