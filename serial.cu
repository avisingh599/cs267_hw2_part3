#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <iostream>
#include "common.h"
#include "core.h"

inline int get_box_index(
     particle_t* particle,
     double box_width,
     int nsquares_per_side
 ) {
     int row = floor(particle->y/box_width);
     int col = floor(particle->x/box_width);
     return col + row*nsquares_per_side;
 }

void put_particle_in_box_1(
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


 void put_particle_in_box_2(
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

//
//  benchmarking program
//
int main( int argc, char **argv )
{
    //int navg,nabsavg=0;
    //double davg,dmin, absmin=1.0, absavg=0.0;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }

    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );

    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    // set the grid and store mesh boundaries
    set_size( n );
    int nsquares_per_side = get_max_nsquares_per_side();
    double box_width = get_box_width(nsquares_per_side);

    int nsquares = nsquares_per_side*nsquares_per_side;
    int boxneighbors[nsquares*9];
    get_box_neighbors(nsquares, nsquares_per_side, boxneighbors);

    // initialize the particles
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    init_particles( n, particles );

    int *box_indices = (int* ) malloc(n * sizeof(int)); 
    int *particle_indices_boxed = (int* ) malloc(n * sizeof(int));
    int *box_iterators = (int* ) malloc(nsquares * sizeof(int)); 
    int *box_positions = (int* ) malloc((nsquares + 1) * sizeof(int));
    int *box_to_particles_odd = (int* ) malloc(20 * n * sizeof(int));
    int *box_to_particles_even = (int* ) malloc(20 * n * sizeof(int));
    int *box_to_num_particles_even = (int* ) malloc(20 * n * sizeof(int));
    int *box_to_num_particles_odd = (int* ) malloc(20 * n * sizeof(int));
    int *box_to_particles = NULL;
    int *box_to_num_particles = NULL;
    int *box_to_particles_next = NULL;
    int *box_to_num_particles_next = NULL;


    // find the neighbors of each mesh square. squares are labelled from bottom left to top right

    for (int b_idx=0; b_idx < nsquares+1; ++b_idx) {
        box_to_num_particles_odd[b_idx] = 0;
        box_to_num_particles_even[b_idx] = 0;
    }

    for (int p_idx = 0; p_idx < n; ++p_idx) {
        int box_i = get_box_index(
                particles + p_idx,
                box_width,
                nsquares_per_side
            );
        box_to_particles_odd[box_i] = p_idx;
        box_to_num_particles_odd[box_i]++;
        box_to_particles_even[box_i] = p_idx;
        box_to_num_particles_even[box_i]++;
    }

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );

    for( int step = 0; step < NSTEPS; step++ )
    {
        // navg = 0;
        // davg = 0.0;
        // dmin = 1.0;
        if (step % 2 == 0) {
            box_to_particles = box_to_particles_even;
            box_to_num_particles = box_to_num_particles_even;
            box_to_particles_next = box_to_particles_odd;
            box_to_num_particles_next = box_to_num_particles_odd;
        } else {
            box_to_particles = box_to_particles_odd;
            box_to_num_particles = box_to_num_particles_odd;
            box_to_particles_next = box_to_particles_even;
            box_to_num_particles_next = box_to_num_particles_even;
        }
        //
        //  compute forces
        //
        for (int i = 0; i < n; i++ )
        {
            particles[i].ax = particles[i].ay = 0;
            int box_i = box_indices[i];
            for (
                int i_in_boxneighbors = box_i * 9;
                i_in_boxneighbors < (box_i+1)*9;
                i_in_boxneighbors++
            ) {
                int neighboring_box_i = boxneighbors[i_in_boxneighbors];
                if (neighboring_box_i != -1) {
                    int num_neigh_parts = box_to_num_particles[neighboring_box_i];
                    for (int j = 0; j < num_neigh_parts; j++) {     
                        int idx_2 = box_to_particles[j]; 
                        apply_force(
                            particles[i],
                            particles[idx_2]
                        );
                    }
                }
            }
        }


        //
        //  move particles
        //
        for( int i = 0; i < n; i++ )
            move( particles[i] );


        //
        //  rebin
        //
        for (int box_i=0; box_i < nsquares; box_i++) {
            box_to_num_particles_next[box_i] = 0;
            for (
                int i_in_boxneighbors = box_i * 9;
                i_in_boxneighbors < (box_i+1)*9;
                i_in_boxneighbors++
            ) {
                if (boxneighbors[i_in_boxneighbors] != -1) {
                    int num_parts = box_to_num_particles[box_i];
                    for (int j = 0; j < num_parts; j++) {     
                        int part_i = box_to_particles[j]; 
                        int new_box_i = get_box_index(
                             particles + part_i,
                             box_width,
                             nsquares
                         );
                        if (new_box_i == box_i) {
                            box_to_particles_next[box_i] = part_i;
                            box_to_num_particles_next[box_i]++;
                        }
                    }
                }
            }
        }

        if( find_option( argc, argv, "-no" ) == -1 )
        {
          //
          // Computing statistical data
          //
          // if (navg) {
          //   absavg +=  davg/navg;
          //   nabsavg++;
          // }
          //if (dmin < absmin) absmin = dmin;

          //
          //  save if necessary
          //
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    }
    simulation_time = read_timer( ) - simulation_time;

    printf( "n = %d, simulation time = %g seconds", n, simulation_time);

    //if( find_option( argc, argv, "-no" ) == -1 )
    //{
      //if (nabsavg) absavg /= nabsavg;
    //
    //  -The minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //
    // printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
    // if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
    // if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    //}
    // printf("\n");

    //
    // Printing summary data
    //
    if( fsum)
        fprintf(fsum,"%d %g\n",n,simulation_time);

    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );
    free( particles );
    if( fsave )
        fclose( fsave );

    return 0;
}
