#ifndef __CORE_H_
#define __CORE_H_

void get_bounds(int nsquares_per_side, double* vertbounds, double* horzbounds);
void get_box_neighbors(int nsquares, int nsquares_per_side, int* boxneighbors);
void bin_particles(particle_t* particledata, int* nparticles_per_bin, int* bin_contents, double* vertbounds, double* horzbounds, int nsquares_per_side, int n);

double get_box_width(int nsquares_per_side);
int get_max_nsquares_per_side();


#endif
