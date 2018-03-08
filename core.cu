#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "common.h"
#include "core.h"

//#include "core_openmp.h"

extern double size;
#define cutoff 0.01  // TODO: check Piazza if this will ever change

//
// Custom classes
//

//
// get bounds in overall domain
//
void get_bounds(int nsquares_per_side, double* vertbounds, double* horzbounds)
{
    vertbounds[0] = 0.0;
    horzbounds[0] = 0.0;
    for (int i = 1; i < nsquares_per_side+1; ++i){
        vertbounds[i] = size*i/nsquares_per_side;
        horzbounds[i] = size*i/nsquares_per_side;
    }
}

double get_box_width(int nsquares_per_side) {
    return size / nsquares_per_side;
}

//
// get 9 neighbors of each box (one of those being itself). if box has less neighbors, padded with -1
//
void get_box_neighbors(int nsquares, int nsquares_per_side, int* boxneighbors)
{
    for (int i = 0; i < nsquares*9; ++i) {
        boxneighbors[i] = -1;
    }
    for (int i = 0; i < nsquares; ++i){
        int row = floor(i/nsquares_per_side);
        int col = i%nsquares_per_side;
        for (int j = -1; j <= 1; ++j){
            for (int k = -1; k <= 1; ++k){
                if (
                    ((row+j) > -1) &&
                    ((row+j) < nsquares_per_side) &&
                    ((col+k) > -1) &&
                    ((col+k) < nsquares_per_side)
                ) {
                    int neighbor_i = (row+j)*nsquares_per_side + (col+k);
                    boxneighbors[i*9 + (j+1)*3 + k+1] = neighbor_i;
                }
            }
        }
    }
}

//
// classify particles into their respective bins and store the particle numbers in bin_contents
// count how many particles are in each bin and store nparticles_per_bin
//
void bin_particles(particle_t* particledata, int* nparticles_per_bin, int* bin_contents, double* vertbounds, double* horzbounds, int nsquares_per_side, int n)
{
    for (int part = 0; part < n; ++part){
        int flag = 0;
        for (int i = 0; i < nsquares_per_side; ++i){
            for (int j = 0; j < nsquares_per_side; ++j){
                if ((particledata[part].x >= vertbounds[i]) && (particledata[part].x < vertbounds[i+1]) && (particledata[part].y >= horzbounds[j]) && (particledata[part].y < horzbounds[j+1])){
                    int bin = j*nsquares_per_side + i;
                    bin_contents[bin*n + nparticles_per_bin[bin]] = part;
                    nparticles_per_bin[bin] += 1;
                    flag = 1;
                    break;
                }
            }
        }
        if (flag == 0){
            printf("WARNING: PARTICLE %i AT x=%f y=%f WAS NOT ASSIGNED TO ANY BOX\n", part, particledata[part].x, particledata[part].y);
        }
    }
}

int get_max_nsquares_per_side() {
    return floor(size / cutoff);
}
