#include "block_decomposition.h"
#include <stdio.h>
#include <stdlib.h>

//find the greatest prime factor of n
int greatest_prime_factor(int n)
{
  int ns = n;
  int np = 2;

  if((n==1)||(n==2))
    return n;

  while(1)
  {
    while(!(ns%np))
      ns = ns/np;
    if(ns==1)
      break;
    np+=1;
  }
  return np;
}

//find the number of blocks 
//in the x,y,z directions
int *number_of_blocks(int nproc, int nx_global, int ny_global, int nz_global)
{
  //initialize np_x, np_y, np_z
  int np_x = 1;
  int np_y = 1;
  int np_z = 1;

  //find the greatest prime factor of the number of MPI processes

  int n_gpf = greatest_prime_factor(nproc);

  //printf("nproc %d n_gpf %d\n",nproc,n_gpf);
  int n_tmp; //for re-arranging nx,ny,nz order

  //base decomposition on whether n_gpf==2
  if(n_gpf!=2)
  {  
    //if we are dealing with two dimensions,
    //we can just assign domain
    if (nz_global==1)
    {
      np_x = n_gpf;
      np_y = nproc/np_x;
      np_z = 1;
    }else{
      //we are in 3-d, so split remainder evenly
      np_x  = n_gpf;
      n_gpf = greatest_prime_factor(nproc/n_gpf);
      if(n_gpf!=2)
      {
        //the next greatest prime is odd, so just split
        np_y = n_gpf;
        np_z = nproc/(np_x*np_y);
      }else{
        //increase ny, nz round-robin
        while(np_x*np_y*np_z < nproc)
        {
          np_y*=2;
          if(np_x*np_y*np_z==nproc)
            break;
          np_z*=2;
        }
      }
    }
  }else{
    //nproc is a power of 2
    //if we are dealing with two dimensions
    //we can just assign domain
    if(nz_global==1)
    {
      np_x = n_gpf;
      np_y = nproc/np_x;
      np_z = 1;
    }else{
      //we are in 3-d, so split remainder evenly

      //increase nx, ny, nz round-robin
      while(np_x*np_y*np_z < nproc)
      {
        np_x*=2;
        if(np_x*np_y*np_z==nproc)
          break;
        np_y*=2;
        if(np_x*np_y*np_z==nproc)
          break;
        np_z*=2;
      }
    }
  }

  //reorder x, y, z

  if(np_z>np_y)
  {
    n_tmp = np_y;
    np_y  = np_z;
    np_z  = n_tmp;
  }
  
  if(np_y>np_x)
  {
    n_tmp = np_x;
    np_x  = np_y;
    np_y  = n_tmp;
  }
  
  if(np_z>np_y)
  {
    n_tmp = np_y;
    np_y  = np_z;
    np_z  = n_tmp;
  }

  //return result
  int *narr = (int *) malloc(3*sizeof(int));
  narr[0] = np_x;
  narr[1] = np_y;
  narr[2] = np_z;
  return narr;
}


//find the local sizes and offsets
//for each block
int *domain_decomposition_block(int nx_global, int ny_global, int nz_global, int nproc_x, int nproc_y, int nproc_z, int procID)
{
  int i,j,k; //indices

  int nprocs = nproc_x*nproc_y*nproc_z; //total number of processors

  int *ix = (int *) malloc(nprocs*sizeof(int));
  int *iy = (int *) malloc(nprocs*sizeof(int));
  int *iz = (int *) malloc(nprocs*sizeof(int));

  int *tiling = (int *) malloc(nproc_x*nproc_y*nproc_z*sizeof(int));

  int n = 0;

  int nx_local, ny_local, nz_local;
  int nx_local_start, ny_local_start, nz_local_start;


  for(i=0;i<nproc_x;i++)
    for(j=0;j<nproc_y;j++)
      for(k=0;k<nproc_z;k++)
      {
        ix[n] = i;
        iy[n] = j;
        iz[n] = k;

        tiling[nproc_z*nproc_y*i + nproc_z*j + k] = n;
        n+=1;
      }

  //set local x, y, z subdomain sizes 
  n = nx_global%nproc_x;

  if(!n)
  {
    //nx_global splits evenly along x procs
    nx_local = nx_global; //nproc_x
    nx_local_start = ix[procID]*nx_local;
  }else{
    nx_local = nx_global; //nproc_x
    if(ix[procID]<n)
    {
      nx_local+=1;
      nx_local_start = ix[procID]*nx_local;
    }else{
      //check nx_local_start offsets -- should n be (n-1) below?
      nx_local_start = n*(nx_local+1) + (ix[procID]-n)*nx_local;
    }
  }
    
  n = ny_global%nproc_y;
  if(!n)
  {
    //ny_global splits evenly along y procs
    ny_local = ny_global; //nproc_y;
    ny_local_start = iy[procID]*ny_local;
  }else{
    ny_local = ny_global; //nproc_y
    if(iy[procID]<n)
    {
      ny_local+=1;
      ny_local_start = iy[procID]*ny_local;
    }else{
      ny_local_start = n*(ny_local+1) + (iy[procID]-n)*ny_local;
    }
  }
    
  n = nz_global%nproc_z;
  if(!n)
  {
    //nz_global splits evenly along z procs
    nz_local = nz_global; //nproc_z
    nz_local_start = iz[procID]*nz_local;
  }else{
    nz_local = nz_global; //nproc_z
    if(iz[procID]<n)
    {
      nz_local+=1;
      nz_local_start = iz[procID]*nz_local;
    }else{
      nz_local_start = n*(nz_local+1) + (iz[procID]-n)*nz_local;
    }
  }
  free(ix);
  free(iy);
  free(iz);
  free(tiling);

  //return the answer
  int *narr = (int *) malloc(6*sizeof(int));
  narr[0] = nx_local;
  narr[1] = ny_local;
  narr[2] = nz_local;
  narr[3] = nx_local_start;
  narr[4] = ny_local_start;
  narr[5] = nz_local_start;
  return narr;
}