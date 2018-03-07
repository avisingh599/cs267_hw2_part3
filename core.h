#ifndef __CORE_H_
#define __CORE_H_
#include <set>
#include <omp.h>

void get_bounds(int nsquares_per_side, double* vertbounds, double* horzbounds);
void get_box_neighbors(int nsquares, int nsquares_per_side, int* boxneighbors);
void bin_particles(particle_t* particledata, int* nparticles_per_bin, int* bin_contents, double* vertbounds, double* horzbounds, int nsquares_per_side, int n);

//
// Custom Stuff
//
double get_box_width(int nsquares_per_side);
int get_max_nsquares_per_side();

typedef struct
{
    particle_t* p;
    int box_index;
    int next_box_index;
    int moved_boxes;
} boxed_particle_t;

class Box
{
    public:
        std::set<int> particles;
        int* lock;
};
void put_particle_in_boxes(
    boxed_particle_t* boxed_particle,
    int p_idx,
    Box** boxes,
    double box_width,
    int nsquares_per_side,
    int n
);
void update_particle_in_boxes(
    boxed_particle_t* boxed_particle,
    int p_idx,
    Box** boxes,
    double box_width,
    int nsquares_per_side,
    int n
);

int get_box_index(
    boxed_particle_t* boxed_particle,
    double box_width,
    int nsquares_per_side
);

#endif