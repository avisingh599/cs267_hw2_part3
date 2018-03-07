#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <cuda.h>
#include "common.h"
#include "core.h"

#define NUM_THREADS 256

extern double size;
//
//  benchmarking program
//

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

// __global__ void init_boxed_particles(boxed_particles_t* boxed_particles, particle_t* particles, int n)
// {

// }

__global__ void gpu_put_particle_in_boxes(boxed_particle_t* boxed_particles, 
                                          int n, 
                                          int nsquares_per_side, 
                                          double box_width, 
                                          Box** boxes)
{
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= n) return;

  boxed_particle_t* boxed_particle = boxed_particles + tid;
  put_particle_in_boxes(
      boxed_particle,
      tid,
      boxes,
      box_width,
      nsquares_per_side,
      n
  );

}

__global__ void gpu_update_particle_in_boxes(boxed_particle_t* boxed_particles, 
                                          int n, 
                                          int nsquares_per_side, 
                                          double box_width, 
                                          Box** boxes)
{
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= n) return;

  boxed_particle_t* boxed_particle = boxed_particles + tid;
  update_particle_in_boxes(
      boxed_particle,
      tid,
      boxes,
      box_width,
      nsquares_per_side,
      n
  );

}

__global__ void compute_forces_grid_gpu(particle_t * particles,
                                        boxed_particle_t* boxed_particles,
                                        Box** boxes, 
                                        int* boxneighbors,
                                        int n, 
                                        )
{
  // Get thread (particle) ID
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= n) return;

  particles[tid].ax = particles[tid].ay = 0;
  int box_idx = boxed_particles[tid].box_index;
  Box* box = boxes[box_idx];
  for (
      int i_in_boxneighbors = box_idx * 9;
      i_in_boxneighbors < (box_idx+1)*9;
      i_in_boxneighbors++
  ) {
      int neighboring_box_i = boxneighbors[i_in_boxneighbors];
      if (neighboring_box_i != -1) {
          Box* neighboring_box = boxes[neighboring_box_i];
          std::set<int>::iterator it;
          for (
              it = neighboring_box->particles.begin();
              it != neighboring_box->particles.end();
              it++
          ) {
              int neighbor_particle_i = *it;
                  apply_force_gpu(
                      particles[tid],
                      particles[neighbor_particle_i]
                      // &dmin,
                      // &davg,
                      // &navg
                  );
          }
      }
  }

}

__global__ void compute_forces_gpu(particle_t * particles, int n)
{
  // Get thread (particle) ID
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= n) return;

  particles[tid].ax = particles[tid].ay = 0;
  for(int j = 0 ; j < n ; j++)
    apply_force_gpu(particles[tid], particles[j]);

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
   
    set_size( n );
    int nsquares_per_side = get_max_nsquares_per_side();
    double box_width = get_box_width(nsquares_per_side);
    Box** boxes = new Box*[nsquares_per_side*nsquares_per_side];

    for (int x = 0; x < nsquares_per_side; x++) {
        for (int y = 0; y < nsquares_per_side; y++) {
            boxes[x+y*nsquares_per_side] = new Box();
        }
    }

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    init_particles( n, particles );

    // GPU particle data structure
    particle_t *d_particles;
    cudaMalloc((void **) &d_particles, n * sizeof(particle_t));

    cudaThreadSynchronize();
    double copy_time = read_timer( );

    // Copy the particles to the GPU
    cudaMemcpy(d_particles, particles, n * sizeof(particle_t), cudaMemcpyHostToDevice);

    cudaThreadSynchronize();
    copy_time = read_timer( ) - copy_time;

    // Initialize the boxed particles and the boxes
    boxed_particle_t *boxed_particles = (boxed_particle_t*) malloc(
        n * sizeof(boxed_particle_t)
    );
    
    for (int i = 0; i < n; i++) {
        boxed_particles[i].p = particles + i;
        boxed_particles[i].box_index = -1;
    }

    boxed_particle_t *d_boxed_particles;
    cudaMalloc((void **) &d_boxed_particles, n * sizeof(boxed_particle_t));
    cudaThreadSynchronize();
    // TODO time this thing later
    cudaMemcpy(d_boxed_particles, boxed_particles, n * sizeof(boxed_particle_t), cudaMemcpyHostToDevice);
    cudaThreadSynchronize();

    int blks = (n + NUM_THREADS - 1) / NUM_THREADS;
    gpu_put_particle_in_boxes <<< blks, NUM_THREADS >>> (d_boxed_particles, n, nsquares_per_side, box_width, boxes);

    // for (int p_idx = 0; p_idx < n; ++p_idx) {
    //     boxed_particle_t* boxed_particle = boxed_particles + p_idx;
    //     put_particle_in_boxes(
    //         boxed_particle,
    //         p_idx,
    //         boxes,
    //         box_width,
    //         nsquares_per_side,
    //         n
    //     );
    // }

    // find the neighbors of each mesh square. squares are labelled from bottom left to top right
    int nsquares = nsquares_per_side*nsquares_per_side;
    int boxneighbors[nsquares*9];
    get_box_neighbors(nsquares, nsquares_per_side, boxneighbors);


    int* d_boxneighbors; 
    cudaMalloc((void **) &d_boxneighbors, n* sizeof(int));
    cudaThreadSynchronize();
    cudaMemcpy(d_boxneighbors, boxneighbors, n * sizeof(int), cudaMemcpyHostToDevice);
    cudaThreadSynchronize();

    Box** d_boxes; 
    cudaMalloc((void **), &d_boxes, nsquares_per_side * nsquares_per_side * sizeof(Box));
    cudaThreadSynchronize();
    cudaMemcpy(d_boxes, boxes, nsquares_per_side * nsquares_per_side * sizeof(Box), cudaMemcpyHostToDevice);
    cudaThreadSynchronize();
    //
    //  simulate a number of time steps
    //
    cudaThreadSynchronize();
    double simulation_time = read_timer( );

    for( int step = 0; step < NSTEPS; step++ )
    {
        //
        //  compute forces
        //

	    int blks = (n + NUM_THREADS - 1) / NUM_THREADS;
      compute_forces_grid_gpu <<< blks, NUM_THREADS >>> (d_particles, d_boxed_particles, d_boxes, d_boxneighbors, n);
        //	compute_forces_gpu <<< blks, NUM_THREADS >>> (d_particles, n);

        //
        //  move particles
        //
	    move_gpu <<< blks, NUM_THREADS >>> (d_particles, n, size);
  
      gpu_update_particle_in_boxes <<< blks, NUM_THREADS >>> (d_boxed_particles, n, nsquares_per_side, box_width, d_boxes);

        //
        //  save if necessary
        //
        if( fsave && (step%SAVEFREQ) == 0 ) {
	    // Copy the particles back to the CPU
            cudaMemcpy(particles, d_particles, n * sizeof(particle_t), cudaMemcpyDeviceToHost);
            save( fsave, n, particles);
	}
    }
    cudaThreadSynchronize();
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "CPU-GPU copy time = %g seconds\n", copy_time);
    printf( "n = %d, simulation time = %g seconds\n", n, simulation_time );
    
    free( particles );
    cudaFree(d_particles);
    if( fsave )
        fclose( fsave );
    
    return 0;
}
