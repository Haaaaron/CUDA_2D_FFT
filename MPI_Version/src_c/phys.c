#include "main.h"

void buffer_init(Parameters params, Init_rank_parameters info, Physical_system *system) {

    int i, j, ind;


    //allocating 2d arrays for physical system
    system->real_buffer         = malloc(params.lx*info.j_len*sizeof(double));
    system->complex_buffer.rest = malloc(params.ly*info.i_len*sizeof(fftw_complex));
    //allocating complex field zeroth case arrays
    system->complex_buffer.kx_final = malloc(params.lyhp*sizeof(fftw_complex));
    system->complex_buffer.kx_zero  = malloc(params.lyhp*sizeof(fftw_complex));

    for (j = 0; j < info.j_len; j++) {
        for ( i = 0; i < params.lx; i++) {
            ind = (j)*params.lx+i;
            system->real_buffer[ind] = 1;

        }
    }
}

