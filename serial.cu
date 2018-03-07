#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <iostream>
#include "common.h"
#include "core.h"


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

    //Box** boxes = new Box*[nsquares_per_side*nsquares_per_side];


    // for (int x = 0; x < nsquares_per_side; x++) {
    //     for (int y = 0; y < nsquares_per_side; y++) {
    //         boxes[x+y*nsquares_per_side] = new Box();
    //     }
    // }

    // initialize the particles
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    init_particles( n, particles );

    // Initialize the boxed particles and the boxes
    // boxed_particle_t *boxed_particles_read = (boxed_particle_t*) malloc(
    //     n * sizeof(boxed_particle_t));
    // boxed_particle_t *boxed_particles_write = (boxed_particle_t*) malloc(
    //     n * sizeof(boxed_particle_t));

    int *box_indices = (int* ) malloc(n * sizeof(int)); 
    int *particle_indices_boxed = (int* ) malloc(n * sizeof(int));
    int *box_iterators = (int* ) malloc(nsquares * sizeof(int)); 
    int *box_positions = (int* ) malloc((nsquares + 1) * sizeof(int));

    // for (int i = 0; i < n; i++) {
    //     boxed_particles_read[i].p = particles + i;
    //     boxed_particles_read[i].box_index = -1;
    // }

    // for (int i = 0; i < n; i++) {
    //     boxed_particles_copy[i].p = particles + i;
    //     boxed_particles_copy[i].box_index = -1;
    // }



    // for (int p_idx = 0; p_idx < n; ++p_idx)
    //     put_particle_in_box_1(
    //         boxed_particles_read,
    //         p_idx,
    //         box_positions,
    //         box_width,
    //         nsquares_per_side
    //         );

    // for (int b_idx=1; b_idx < nsquares+1; ++b_idx)
    //     box_positions[b_idx] = box_positions[b_idx] + box_positions[b_idx-1];

    // for (int p_idx = 0; p_idx < n; ++p_idx)
    //     put_particle_in_box_2(
    //         boxed_particles_read,
    //         boxed_particles_write,
    //         p_idx,
    //         box_positions,
    //         box_iterators,
    //         box_width,
    //         nsquares_per_side
    //         );


    for (int b_idx=0; b_idx < nsquares+1; ++b_idx)
        box_positions[b_idx] = 0;

    for (int b_idx=0; b_idx < nsquares; ++b_idx)
        box_iterators[b_idx] = 0;

    for (int p_idx = 0; p_idx < n; ++p_idx)
        put_particle_in_box_1(
            particles,
            p_idx,
            box_positions,
            box_width,
            nsquares_per_side
            );

    for (int b_idx=1; b_idx < nsquares+1; ++b_idx)
        box_positions[b_idx] = box_positions[b_idx] + box_positions[b_idx-1];

    for (int p_idx = 0; p_idx < n; ++p_idx)
        put_particle_in_box_2(
            particles,
            p_idx,
            box_positions,
            box_iterators,
            box_indices,
            particle_indices_boxed,
            box_width,
            nsquares_per_side
            );

    // for (int p_idx = 0; p_idx < n; ++p_idx) {
    //     boxed_particles_read[p_idx].p = boxed_particles_write[p_idx].p;
    //     boxed_particles_read[p_idx].box_index = boxed_particles_write[p_idx].box_index;
    // }

    // find the neighbors of each mesh square. squares are labelled from bottom left to top right

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );

    for( int step = 0; step < NSTEPS; step++ )
    {
        // navg = 0;
        // davg = 0.0;
        // dmin = 1.0;
        //
        //  compute forces
        //
        for (int i = 0; i < n; i++ )
        {
            particles[i].ax = particles[i].ay = 0;
            int box_idx = box_indices[i];
            //Box* box = boxes[box_idx];
            for (
                int i_in_boxneighbors = box_idx * 9;
                i_in_boxneighbors < (box_idx+1)*9;
                i_in_boxneighbors++
            ) {
                int neighboring_box_i = boxneighbors[i_in_boxneighbors];
                if (neighboring_box_i != -1) {
                    //Box* neighboring_box = boxes[neighboring_box_i];
                    //std::set<int>::iterator it;
                    for (
                        int it = box_positions[neighboring_box_i];
                        it != box_positions[neighboring_box_i+1];
                        it++
                    ) {     
                        int idx_2 = particle_indices_boxed[it]; 
                        // int neighbor_particle_i = *it;
                            apply_force(
                                particles[i],
                                particles[idx_2]
                            );
                        //apply_force(*boxed_particles_read[i].p, *boxed_particles_read[it].p);
                     }
                    }
                }
            }


        //
        //  move particles
        //
        for( int i = 0; i < n; i++ )
            move( particles[i] );


        for (int b_idx=0; b_idx < nsquares+1; ++b_idx)
            box_positions[b_idx] = 0;

        for (int b_idx=0; b_idx < nsquares; ++b_idx)
            box_iterators[b_idx] = 0;

        for (int p_idx = 0; p_idx < n; ++p_idx)
            put_particle_in_box_1(
                particles,
                p_idx,
                box_positions,
                box_width,
                nsquares_per_side
                );

        for (int b_idx=1; b_idx < nsquares+1; ++b_idx)
            box_positions[b_idx] = box_positions[b_idx] + box_positions[b_idx-1];

        for (int p_idx = 0; p_idx < n; ++p_idx)
            put_particle_in_box_2(
                particles,
                p_idx,
                box_positions,
                box_iterators,
                box_indices,
                particle_indices_boxed,
                box_width,
                nsquares_per_side
                );

        // Faster than doing put_particles_in_boxes since this can assume the
        // particles are already in a box.
        // for (int p_idx = 0; p_idx < n; ++p_idx){
        //     boxed_particle_t* boxed_particle = boxed_particles + p_idx;
        //     update_particle_in_boxes(
        //         boxed_particle,
        //         p_idx,
        //         boxes,
        //         box_width,
        //         nsquares_per_side,
        //         n
        //     );
        // }


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

    // for (int x = 0; x < nsquares_per_side; x++) {
    //     for (int y = 0; y < nsquares_per_side; y++) {
    //         delete boxes[x+y*nsquares_per_side];
    //     }
    // }
    // free( boxed_particles_read );
    // free( boxed_particles_write );

    //delete boxes;
    return 0;
}
