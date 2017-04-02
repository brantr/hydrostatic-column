#ifndef BLOCK_DECOMPOSITION_DISK_H
#define BLOCK_DECOMPOSITION_DISK_H
#include "block_decomposition.h"
#include <stdio.h>

//find the greatest prime factor of n
int greatest_prime_factor(int n);

//find the number of blocks 
//in the x,y,z directions
int *number_of_blocks(int nproc, int nx_global, int ny_global, int nz_global);

//find the local sizes and offsets
//for each block
int *domain_decomposition_block(int nx_global, int ny_global, int nz_global, int nproc_x, int nproc_y, int nproc_z, int procID);


#endif  //BLOCK_DECOMPOSITION_DISK_H
