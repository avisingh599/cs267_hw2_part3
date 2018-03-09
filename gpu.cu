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

  //printf("box pos set %d: %d\n", tid, box_positions[tid]);

}

__global__ void compute_box_positions(int* box_positions, int nsquares)
{ 
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= (nsquares +1) || tid < 1) return;
  ///TODO BUG!! alternate, have some global "ticker"? 
  box_positions[tid] = box_positions[tid] + box_positions[tid-1];
}

__global__ void compute_box_positions_serial(int* box_positions, int nsquares)
{ 
  //int tid = threadIdx.x + blockIdx.x * blockDim.x;
  //if(tid >= (nsquares +1) || tid < 1) return;
  ///TODO BUG!! alternate, have some global "ticker"? 
  
  for (int i=1; i<nsquares+1; i++) {
  box_positions[i] = box_positions[i] + box_positions[i-1];
  //printf("box pos %d: %d\n", i, box_positions[i]);
  }
}

__global__ void set_box_iterators_zero(int* box_iterators, int nsquares)
{  
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= nsquares) return;

  box_iterators[tid] = 0;

  //printf("box itr set %d: %d\n", tid, box_iterators[tid]);

}

__global__ void put_particle_in_box_1_gpu(particle_t * particles,
                                        int* box_positions,
                                        double box_width,
                                        int nsquares_per_side,                                         
                                        int n)
{
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= n) return;

  particle_t *particle = particles + tid;
  int box_index = get_box_index(particle,box_width,nsquares_per_side);
  //printf("tid: %d box index: %d\n", pidx, box_index);
  //box_positions[box_index+1] +=1;
  //int* add = box_positions + box_index + 1;
  //printf("tid: %d box position before: %d\n", pidx, box_positions[box_index+1]);
  //atomicAdd(add,1);
  atomicAdd(box_positions + box_index + 1,1);
  //printf("tid: %d box position after: %d\n", pidx, box_positions[box_index+1]);

}

__global__ void put_particle_in_box_2_gpu(particle_t * particles,
                                        int* box_positions,
                                        int* box_iterators,
                                        int* box_indices,
                                        int* particle_indices_boxed,
                                        int* locks,
                                        double box_width, 
                                        int nsquares_per_side,                                         
                                        int n)
{
  int tid = threadIdx.x + blockIdx.x * blockDim.x;
  if(tid >= n) return;

  particle_t *particle = particles + tid;
  int box_index = get_box_index(particle,box_width,nsquares_per_side);
  box_indices[tid] = box_index; 

  cudaLock(locks + box_index);
  int store_idx = box_positions[box_index] + box_iterators[box_index]; 
  
  //printf("store idx: %d\n", store_idx);
  //int* add = box_iterators + box_index; 
  //atomicAdd(add, 1);
  box_iterators[box_index] += 1;
  cudaUnlock(locks + box_index);
  
  particle_indices_boxed[store_idx] = tid; 


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
              it < box_positions[neighboring_box_i+1];
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
    //int boxneighbors[nsquares*9];
    int* boxneighbors = (int* ) malloc((nsquares*9) * sizeof(int)); 
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

    // cudaMemcpy(d_box_indices, box_indices, n * sizeof(int), cudaMemcpyHostToDevice);
    // cudaMemcpy(d_particle_indices_boxed, particle_indices_boxed, n * sizeof(int), cudaMemcpyHostToDevice);
    // cudaMemcpy(d_box_iterators, box_iterators, nsquares * sizeof(int), cudaMemcpyHostToDevice);
    // cudaMemcpy(d_box_positions, box_positions, (nsquares+1) * sizeof(int), cudaMemcpyHostToDevice);
    // cudaMemcpy(d_boxneighbors, boxneighbors, (nsquares*9) * sizeof(int), cudaMemcpyHostToDevice);

  
    int *d_locks; 
    cudaMalloc((void **) &d_locks, nsquares * sizeof(int));
    cudaMemset(d_locks, 0, nsquares * sizeof(int));

    cudaThreadSynchronize();
    int blks_particles = (n + NUM_THREADS - 1) / NUM_THREADS;
    int blks_boxes = (nsquares + 1 + NUM_THREADS - 1) / NUM_THREADS;


    set_box_iterators_zero <<< blks_boxes, NUM_THREADS >>> (d_box_iterators, nsquares); 
    set_box_positions_zero <<< blks_boxes, NUM_THREADS >>> (d_box_positions, nsquares);
    put_particle_in_box_1_gpu <<< blks_particles, NUM_THREADS >>> (d_particles, d_box_positions, box_width, nsquares_per_side, n); 
    
    //compute_box_positions <<< blks_boxes, NUM_THREADS >>> (d_box_positions, nsquares);
    compute_box_positions_serial <<<1,1>>> (d_box_positions, nsquares);
    
    put_particle_in_box_2_gpu <<< blks_particles, NUM_THREADS >>> (d_particles, d_box_positions, d_box_iterators, d_box_indices, d_particle_indices_boxed, d_locks, box_width, nsquares_per_side, n); 
    cudaThreadSynchronize();


    // cudaMemcpy(box_positions, d_box_positions, (nsquares +1) * sizeof(int), cudaMemcpyDeviceToHost);
    // printf("DEBUG box positions: %d %d %d \n", box_positions[0], box_positions[10], box_positions[20]);
    
    // cudaMemcpy(box_iterators, d_box_iterators, nsquares * sizeof(int), cudaMemcpyDeviceToHost);
    // printf("DEBUG box iterators: %d %d %d \n", box_iterators[0], box_iterators[10], box_iterators[20]);

    // cudaMemcpy(particle_indices_boxed, d_particle_indices_boxed, n * sizeof(int), cudaMemcpyDeviceToHost);
    // printf("DEBUG particle indices: %d %d %d \n", box_iterators[0], box_iterators[10], box_iterators[20]);
    
    for( int step = 0; step < NSTEPS; step++ )
    {
        //
        //  compute forces
        //

	    //compute_forces_gpu <<< blks_particles, NUM_THREADS >>> (d_particles, n);
      compute_forces_grid_gpu <<< blks_particles, NUM_THREADS >>> (d_particles, d_box_positions, d_boxneighbors, d_particle_indices_boxed, d_box_indices, n);

      // cudaError_t err = cudaGetLastError();
      // if (err != cudaSuccess) 
      //   printf("Error: %s\n", cudaGetErrorString(err));
      //   //
        //  move particles
        //
	    move_gpu <<< blks_particles, NUM_THREADS >>> (d_particles, n, size);

      set_box_iterators_zero <<< blks_boxes, NUM_THREADS >>> (d_box_iterators, nsquares); 
      set_box_positions_zero <<< blks_boxes, NUM_THREADS >>> (d_box_positions, nsquares);
      put_particle_in_box_1_gpu <<< blks_particles, NUM_THREADS >>> (d_particles, d_box_positions, box_width, nsquares_per_side, n); 
//    compute_box_positions <<< blks_boxes, NUM_THREADS >>> (d_box_positions, nsquares);
      compute_box_positions_serial <<<1,1>>> (d_box_positions, nsquares);
      put_particle_in_box_2_gpu <<< blks_particles, NUM_THREADS >>> (d_particles, d_box_positions, d_box_iterators, d_box_indices, d_particle_indices_boxed, d_locks, box_width, nsquares_per_side, n); 

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

    cudaMemcpy(box_positions, d_box_positions, (nsquares +1) * sizeof(int), cudaMemcpyDeviceToHost);
    printf("DEBUG box positions: %d %d %d \n", box_positions[0], box_positions[10], box_positions[20]);
    
    printf( "CPU-GPU copy time = %g seconds\n", copy_time);
    printf( "n = %d, simulation time = %g seconds\n", n, simulation_time );
    
    free( particles );
    cudaFree(d_particles);
    if( fsave )
        fclose( fsave );
    
    return 0;
}
