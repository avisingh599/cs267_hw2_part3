#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <cuda.h>
#include "common.h"
#include "core.h"

#define NUM_THREADS 256
#define MAX_N_PARTS_PER_BOX 10

extern double size;
//
//  benchmarking program
//
inline int get_box_index_serial(
     particle_t* particle,
     double box_width,
     int nsquares_per_side
 ) {
     int row = floor(particle->y/box_width);
     int col = floor(particle->x/box_width);
     return col + row*nsquares_per_side;
 }

__device__ int get_box_index(
    particle_t* particle,
    double box_width,
    int nsquares_per_side
) {
    int row = floor(particle->y/box_width);
    int col = floor(particle->x/box_width);
    return col + row*nsquares_per_side;
}

__device__ void apply_force_gpu(particle_t &particle, particle_t &neighbor)
{
  double dx = neighbor.x - particle.x;
  double dy = neighbor.y - particle.y;
  double r2 = dx * dx + dy * dy;
  if( r2 > cutoff*cutoff )
      return;
  //r2 = fmax( r2, min_r*min_r );
  r2 = (r2 > min_r*min_r) ? r2 : min_r*min_r;
  double r = sqrt( r2 );

  //
  //  very simple short-range repulsive force
  //
  double coef = ( 1 - cutoff / r ) / r2 / mass;
  particle.ax += coef * dx;
  particle.ay += coef * dy;

}

__global__ void compute_forces_grid_gpu(
        particle_t * particles,
        int* box_to_particles,
        int* box_to_num_particles,
        int* boxneighbors,
        int nsquares_per_side,
        double box_width,
        int n
) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if(tid >= n) return;

    particles[tid].ax = particles[tid].ay = 0;
    int box_i = get_box_index(particles + tid, box_width, nsquares_per_side);

    for (
        int i_in_boxneighbors = box_i * 9;
        i_in_boxneighbors < (box_i+1)*9;
        i_in_boxneighbors++
    ) {
        int neighboring_box_i = boxneighbors[i_in_boxneighbors];
        if (neighboring_box_i != -1) {
            int num_neigh_parts = box_to_num_particles[neighboring_box_i];
            for (int j = 0; j < num_neigh_parts; j++) {
                int idx_2 = box_to_particles[MAX_N_PARTS_PER_BOX * neighboring_box_i + j]; 
                apply_force_gpu(
                    particles[tid],
                    particles[idx_2]
                );
            }
        }
    }
}

__global__ void rebin_particles(
        particle_t * particles,
        int* box_to_particles,
        int* box_to_num_particles,
        int* box_to_particles_next,
        int* box_to_num_particles_next,
        int* boxneighbors,
        int nsquares,
        int nsquares_per_side,
        double box_width,
        int n
) {
    int box_i = threadIdx.x + blockIdx.x * blockDim.x;
    if (box_i >= nsquares) return;
    int new_num_parts = 0;
    for (
        int i_in_boxneighbors = box_i * 9;
        i_in_boxneighbors < (box_i+1)*9;
        i_in_boxneighbors++
    ) {
        int neighboring_box_i = boxneighbors[i_in_boxneighbors];
        if (neighboring_box_i != -1) {
            int num_parts_neighbor = box_to_num_particles[neighboring_box_i];
            for (int j = 0; j < num_parts_neighbor; j++) {
                int part_i = box_to_particles[MAX_N_PARTS_PER_BOX * neighboring_box_i + j];
                int new_box_i = get_box_index(
                     particles + part_i,
                     box_width,
                     nsquares_per_side
                 );
                if (new_box_i == box_i) {
                    box_to_particles_next[MAX_N_PARTS_PER_BOX * box_i + new_num_parts] = part_i;
                    new_num_parts++;
                }
            }
        }
    }
    box_to_num_particles_next[box_i] = new_num_parts;
}


__global__ void move_gpu (particle_t * particles, int n, double size)
{

  // Get thread (particle) ID
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= n) return;

  particle_t * p = &particles[tid];
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    p->vx += p->ax * dt;
    p->vy += p->ay * dt;
    p->x  += p->vx * dt;
    p->y  += p->vy * dt;

    //
    //  bounce from walls
    //
    while( p->x < 0 || p->x > size )
    {
        p->x  = p->x < 0 ? -(p->x) : 2*size-p->x;
        p->vx = -(p->vx);
    }
    while( p->y < 0 || p->y > size )
    {
        p->y  = p->y < 0 ? -(p->y) : 2*size-p->y;
        p->vy = -(p->vy);
    }

}



int main( int argc, char **argv )
{    
    // This takes a few seconds to initialize the runtime
    cudaThreadSynchronize(); 

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );

    // GPU particle data structure
    particle_t * d_particles;
    cudaMalloc((void **) &d_particles, n * sizeof(particle_t));

    set_size( n );
    int nsquares_per_side = get_max_nsquares_per_side();
    double box_width = get_box_width(nsquares_per_side);

    init_particles( n, particles );

    int nsquares = nsquares_per_side*nsquares_per_side;
    int boxneighbors[nsquares*9];
    get_box_neighbors(nsquares, nsquares_per_side, boxneighbors);

    int *box_to_particles = (int* ) malloc(MAX_N_PARTS_PER_BOX * nsquares * sizeof(int));
    int *box_to_num_particles = (int* ) malloc(nsquares * sizeof(int));

    // Initialize the data structures serially once.
    for (int b_idx=0; b_idx < nsquares; b_idx++) {
        box_to_num_particles[b_idx] = 0;
    }
    for (int p_idx = 0; p_idx < n; p_idx++) {
        int box_i = get_box_index_serial(
            particles + p_idx,
            box_width,
            nsquares_per_side
        );
        int n_parts = box_to_num_particles[box_i];
        box_to_particles[MAX_N_PARTS_PER_BOX * box_i + n_parts] = p_idx;
        box_to_num_particles[box_i]++;
    }

    cudaThreadSynchronize();
    double copy_time = read_timer( );

    // Copy the particles to the GPU
    cudaMemcpy(d_particles, particles, n * sizeof(particle_t), cudaMemcpyHostToDevice);

    cudaThreadSynchronize();
    copy_time = read_timer( ) - copy_time;

    //
    //  simulate a number of time steps
    //
    cudaThreadSynchronize();
    double simulation_time = read_timer( );

    //copy box-y stuff to CUDA memory 
    int *d_box_to_particles_odd;
    int *d_box_to_particles_even;
    int *d_box_to_num_particles_odd;
    int *d_box_to_num_particles_even;
    int *d_boxneighbors;

    int *d_box_to_particles;
    int *d_box_to_num_particles;
    int *d_box_to_particles_next;
    int *d_box_to_num_particles_next;

    cudaMalloc(
        (void **) &d_box_to_particles_odd,
        MAX_N_PARTS_PER_BOX * nsquares * sizeof(int)
    );
    cudaMalloc((void **) &d_box_to_num_particles_odd, nsquares * sizeof(int)
    );
    cudaMalloc(
        (void **) &d_box_to_particles_even,
        MAX_N_PARTS_PER_BOX * nsquares * sizeof(int)
    );
    cudaMalloc((void **) &d_box_to_num_particles_even, nsquares * sizeof(int));
    cudaMalloc((void **) &d_boxneighbors, (nsquares*9) * sizeof(int));

    cudaMemcpy(
        d_box_to_particles_odd,
        box_to_particles,
        MAX_N_PARTS_PER_BOX * nsquares * sizeof(int),
        cudaMemcpyHostToDevice
    );
    cudaMemcpy(
        d_box_to_num_particles_odd,
        box_to_num_particles,
        nsquares * sizeof(int),
        cudaMemcpyHostToDevice
    );
    cudaMemcpy(
        d_box_to_particles_even,
        box_to_particles,
        MAX_N_PARTS_PER_BOX * nsquares * sizeof(int),
        cudaMemcpyHostToDevice
    );
    cudaMemcpy(
        d_box_to_num_particles_even,
        box_to_num_particles,
        nsquares * sizeof(int),
        cudaMemcpyHostToDevice
    );
    cudaMemcpy(d_boxneighbors, boxneighbors, (nsquares*9) * sizeof(int), cudaMemcpyHostToDevice);

    cudaThreadSynchronize();
    int n_blks_particles = (n + NUM_THREADS - 1) / NUM_THREADS;
    int n_blks_boxes = (nsquares + NUM_THREADS - 1) / NUM_THREADS;


    cudaThreadSynchronize();

    for( int step = 0; step < NSTEPS; step++ )
    {
        if (step % 2 == 0) {
            d_box_to_particles = d_box_to_particles_even;
            d_box_to_num_particles = d_box_to_num_particles_even;
            d_box_to_particles_next = d_box_to_particles_odd;
            d_box_to_num_particles_next = d_box_to_num_particles_odd;
        } else {
            d_box_to_particles = d_box_to_particles_odd;
            d_box_to_num_particles = d_box_to_num_particles_odd;
            d_box_to_particles_next = d_box_to_particles_even;
            d_box_to_num_particles_next = d_box_to_num_particles_even;
        }

        //
        //  compute forces
        //
        compute_forces_grid_gpu <<< n_blks_particles, NUM_THREADS >>> (
            d_particles,
            d_box_to_particles,
            d_box_to_num_particles,
            d_boxneighbors,
            nsquares_per_side,
            box_width,
            n
        );

        //
        //  move particles
        //
        move_gpu <<< n_blks_particles, NUM_THREADS >>> (d_particles, n, size);

        //
        //  rebin particles
        //
        rebin_particles <<< n_blks_boxes, NUM_THREADS >>> (
            d_particles,
            d_box_to_particles,
            d_box_to_num_particles,
            d_box_to_particles_next,
            d_box_to_num_particles_next,
            d_boxneighbors,
            nsquares,
            nsquares_per_side,
            box_width,
            n
        );

        //
        //  save if necessary
        //
        if( fsave && (step%SAVEFREQ) == 0 ) {
            cudaMemcpy(particles, d_particles, n * sizeof(particle_t), cudaMemcpyDeviceToHost);
            save( fsave, n, particles);
        }
    }
    cudaThreadSynchronize();
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "CPU-GPU copy time = %g seconds\n", copy_time);
    printf( "n = %d, simulation time = %g seconds\n", n, simulation_time );
    
    free( particles );
    free(box_to_particles);
    free(box_to_num_particles);
    cudaFree(d_particles);
    cudaFree(d_box_to_particles_odd);
    cudaFree(d_box_to_particles_even);
    cudaFree(d_box_to_num_particles_odd);
    cudaFree(d_box_to_num_particles_even);
    cudaFree(d_boxneighbors);
    if( fsave )
        fclose( fsave );
    
    return 0;
}
