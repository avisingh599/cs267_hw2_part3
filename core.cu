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

extern double size;

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

int get_max_nsquares_per_side() {
    return floor(size / cutoff);
}
