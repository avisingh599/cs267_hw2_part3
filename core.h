#ifndef __CORE_H_
#define __CORE_H_
#include <set>
#include <omp.h>

void get_bounds(int nsquares_per_side, double* vertbounds, double* horzbounds);
void get_box_neighbors(int nsquares, int nsquares_per_side, int* boxneighbors);
void bin_particles(particle_t* particledata, int* nparticles_per_bin, int* bin_contents, double* vertbounds, double* horzbounds, int nsquares_per_side, int n);

double get_box_width(int nsquares_per_side);
int get_max_nsquares_per_side();

void put_particle_in_box_1(
    particle_t* particles,
    int pidx,
    int* box_positions,
    double box_width,
    int nsquares_per_side
); 

void put_particle_in_box_2(
    particle_t* particles,
    int pidx,
    int* box_positions,
    int* box_iterators,
    int* box_indices,
    int* particle_indices_boxed,
    double box_width,
    int nsquares_per_side
);


int get_box_index(
    particle_t* particle,
    double box_width,
    int nsquares_per_side
);


#endif