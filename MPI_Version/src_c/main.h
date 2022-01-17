#include "mpi.h"
#include <math.h>
#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <string.h> 
#include <time.h>

/*
Default parameters
------------------
*/
typedef struct {
    int lx;
    int ly;
    int lxhp; 
    int lyhp;
    int lhx;
    int lhy;
    double pi;
    double dt; 
} Parameters;

typedef struct {
    int myid;
    int numprocs;
    int lyr,lxr;
    int rcounts, scounts;
    int sumscounts, sumrcounts;
    int* rdispls;
    int* sdispls;
    int* imin,* imax;
    int* jmin,* jmax;
    int iminmyid, imaxmyid;
    int jminmyid, jmaxmyid;
    int j_len;
    int i_len;
} Init_rank_parameters;

//Structs for FFTW arrays and plans
typedef struct {
    fftw_plan* forward_plan;
    fftw_plan* backward_plan;
    fftw_complex* complex_arr;
    double* real_arr;
} fft_real_complex;

typedef struct {
    fftw_plan* forward_plan;
    fftw_plan* backward_plan;
    fftw_complex* f_arr;
    fftw_complex* b_arr;
} fft_complex_complex;

typedef struct {
    double start_time;
    double* total_time;
} Timer;

typedef struct {
    fft_real_complex psi;
    fft_real_complex zero;
    fft_real_complex lxhp;
    fft_complex_complex psic;
    double complex* send;
    double complex* recv;
    struct Time {
        Timer fft_plan;
        Timer before_MPI;
        Timer during_MPI;
        Timer after_MPI;
        Timer run_time;
        Timer without_init;
    } time;

} fft_builds;

typedef struct {
    fftw_complex* kx_zero,* kx_final,* rest;
} Complex_field;

typedef struct {
    Complex_field complex_buffer;
    Complex_field complex_nt;
    Complex_field complex_psi;
    double* real_buffer;
} Physical_system;



extern void buffer_init(Parameters params, Init_rank_parameters info, Physical_system* system);

extern void fft_init(Parameters params, Init_rank_parameters* info, fft_builds* plans);
extern void fft_forward(Parameters params, Init_rank_parameters info, fft_builds plans, double* input_arr, Complex_field output_arr);
extern void fft_backward(Parameters params, Init_rank_parameters info, fft_builds builds, Complex_field input_arr, double* output_arr);

extern void init_params(int argc, char* argv[], int* size_x, int* size_y, int* iteration_count);
extern void output_time(struct Time time, Init_rank_parameters info);
