#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <cuda.h>
#include "common.h"
#include "core.cuh"

#define NUM_THREADS 256

extern double size;
//
//  benchmarking program
//

__device__ void cudaLock(int *lock) 
{                                                              
  while (atomicCAS(lock, 0, 1) != 0);                                                              
}                                                                                                  
                                                                                                   
__device__ void cudaUnlock(int *lock)
{                                                            
  atomicExch(lock, 0);                                                                              
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

__device__ void put_particle_in_box_1(
    particle_t* particles,
    int pidx,
    int* box_positions,
    double box_width,
    int nsquares_per_side
) {
    particle_t *particle = particles + pidx;
    int box_index = get_box_index(particle,box_width,nsquares_per_side);
    box_positions[box_index+1] +=1;
}


__device__ void put_particle_in_box_2(
    particle_t* particles,
    int pidx,
    int* box_positions,
    int* box_iterators,
    int* box_indices,
    int* particle_indices_boxed,
    double box_width,
    int nsquares_per_side
) {

    particle_t *particle = particles + pidx;
    int box_index = get_box_index(particle,box_width,nsquares_per_side);
    
    int store_idx = box_positions[box_index] + box_iterators[box_index]; 
    box_iterators[box_index] += 1; 

    box_indices[pidx] = box_index; 
    particle_indices_boxed[store_idx] = pidx; 
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

__global__ void set_box_positions_zero(int* box_positions, int nsquares)
{
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= (nsquares+1)) return;

  box_positions[tid] = 0;
}

__global__ void compute_box_positions(int* box_positions, int nsquares)
{
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= (nsquares +1) || tid < 1) return;

  box_positions[tid] = box_positions[tid] + box_positions[tid-1];
}

__global__ void set_box_iterators_zero(int* box_iterators, int nsquares)
{
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= nsquares) return;

  box_iterators[tid] = 0;
}

__global__ void put_particle_in_box_1_gpu(particle_t * particles,
                                        int* box_positions,
                                        double box_width,
                                        int nsquares_per_side,                                         
                                        int n)
{
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= n) return;

  put_particle_in_box_1(
            particles,
            tid,
            box_positions,
            box_width,
            nsquares_per_side
            );

}

__global__ void put_particle_in_box_2_gpu(particle_t * particles,
                                        int* box_positions,
                                        int* box_iterators,
                                        int* box_indices,
                                        int* particle_indices_boxed,
                                        double box_width, 
                                        int nsquares_per_side,                                         
                                        int n)
{
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= n) return;

  put_particle_in_box_2(
      particles,
      tid,
      box_positions,
      box_iterators,
      box_indices,
      particle_indices_boxed,
      box_width,
      nsquares_per_side
      );
}

__global__ void compute_forces_grid_gpu(particle_t * particles,
                                        int* box_positions,                                         
                                        int* boxneighbors,
                                        int* particle_indices_boxed,
                                        int* box_indices,
                                        int n
                                        )
{

  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= n) return;

  particles[tid].ax = particles[tid].ay = 0;
  int box_idx = box_indices[tid];

  for (
      int i_in_boxneighbors = box_idx * 9;
      i_in_boxneighbors < (box_idx+1)*9;
      i_in_boxneighbors++
  ) {
      int neighboring_box_i = boxneighbors[i_in_boxneighbors];
      if (neighboring_box_i != -1) {

          for (
              int it = box_positions[neighboring_box_i];
              it != box_positions[neighboring_box_i+1];
              it++
          ) {     
              int idx_2 = particle_indices_boxed[it]; 
                  apply_force_gpu(
                      particles[tid],
                      particles[idx_2]
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
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );

    // GPU particle data structure
    particle_t * d_particles;
    cudaMalloc((void **) &d_particles, n * sizeof(particle_t));

    set_size( n );
    int nsquares_per_side = get_max_nsquares_per_side();
    double box_width = get_box_width(nsquares_per_side);

    int nsquares = nsquares_per_side*nsquares_per_side;
    int boxneighbors[nsquares*9];
    get_box_neighbors(nsquares, nsquares_per_side, boxneighbors);

    init_particles( n, particles );


    int *box_indices = (int* ) malloc(n * sizeof(int)); 
    int *particle_indices_boxed = (int* ) malloc(n * sizeof(int));
    int *box_iterators = (int* ) malloc(nsquares * sizeof(int)); 
    int *box_positions = (int* ) malloc((nsquares + 1) * sizeof(int));


    // find the neighbors of each mesh square. squares are labelled from bottom left to top right

    // for (int b_idx=0; b_idx < nsquares+1; ++b_idx)
    //     box_positions[b_idx] = 0;

    // for (int b_idx=0; b_idx < nsquares; ++b_idx)
    //     box_iterators[b_idx] = 0;

    // for (int p_idx = 0; p_idx < n; ++p_idx)
    //     put_particle_in_box_1(
    //         particles,
    //         p_idx,
    //         box_positions,
    //         box_width,
    //         nsquares_per_side
    //         );

    // for (int b_idx=1; b_idx < nsquares+1; ++b_idx)
    //     box_positions[b_idx] = box_positions[b_idx] + box_positions[b_idx-1];

    // for (int p_idx = 0; p_idx < n; ++p_idx)
    //     put_particle_in_box_2(
    //         particles,
    //         p_idx,
    //         box_positions,
    //         box_iterators,
    //         box_indices,
    //         particle_indices_boxed,
    //         box_width,
    //         nsquares_per_side
    //         );

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
    int *d_box_indices;
    int *d_particle_indices_boxed;
    int *d_box_iterators;
    int *d_box_positions; 
    int* d_boxneighbors; 

    cudaMalloc((void **) &d_box_indices, n * sizeof(int));
    cudaMalloc((void **) &d_particle_indices_boxed, n * sizeof(int));
    cudaMalloc((void **) &d_box_iterators, nsquares * sizeof(int));
    cudaMalloc((void **) &d_box_positions, (nsquares+1) * sizeof(int));
    cudaMalloc((void **) &d_boxneighbors, (nsquares*9) * sizeof(int));

    cudaMemcpy(d_box_indices, box_indices, n * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_particle_indices_boxed, particle_indices_boxed, n * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_box_iterators, box_iterators, nsquares * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_box_positions, box_positions, (nsquares+1) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_boxneighbors, boxneighbors, (nsquares*9) * sizeof(int), cudaMemcpyHostToDevice);

    cudaThreadSynchronize();
    int blks_particles = (n + NUM_THREADS - 1) / NUM_THREADS;
    int blks_boxes = (nsquares + 1 + NUM_THREADS - 1) / NUM_THREADS;


    set_box_iterators_zero <<< blks_boxes, NUM_THREADS >>> (d_box_iterators, nsquares); 
    set_box_positions_zero <<< blks_boxes, NUM_THREADS >>> (d_box_positions, nsquares);
    put_particle_in_box_1_gpu <<< blks_particles, NUM_THREADS >>> (d_particles, d_box_positions, box_width, nsquares_per_side, n); 
    compute_box_positions <<< blks_boxes, NUM_THREADS >>> (d_box_positions, nsquares);
    put_particle_in_box_2_gpu <<< blks_particles, NUM_THREADS >>> (d_particles, d_box_positions, d_box_iterators, d_box_indices, d_particle_indices_boxed, box_width, nsquares_per_side, n); 
    cudaThreadSynchronize();

    for( int step = 0; step < NSTEPS; step++ )
    {
        //
        //  compute forces
        //

	    //compute_forces_gpu <<< blks_particles, NUM_THREADS >>> (d_particles, n);
      compute_forces_grid_gpu <<< blks_particles, NUM_THREADS >>> (d_particles, d_box_positions, d_boxneighbors, d_particle_indices_boxed, d_box_indices, n);

        //
        //  move particles
        //
	    move_gpu <<< blks_particles, NUM_THREADS >>> (d_particles, n, size);

      set_box_iterators_zero <<< blks_boxes, NUM_THREADS >>> (d_box_iterators, nsquares); 
      set_box_positions_zero <<< blks_boxes, NUM_THREADS >>> (d_box_positions, nsquares);
      put_particle_in_box_1_gpu <<< blks_particles, NUM_THREADS >>> (d_particles, d_box_positions, box_width, nsquares_per_side, n); 
      compute_box_positions <<< blks_boxes, NUM_THREADS >>> (d_box_positions, nsquares);
      put_particle_in_box_2_gpu <<< blks_particles, NUM_THREADS >>> (d_particles, d_box_positions, d_box_iterators, d_box_indices, d_particle_indices_boxed, box_width, nsquares_per_side, n); 

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
