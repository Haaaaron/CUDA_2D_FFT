#include "mpi.h"
#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <string.h>

#include "main.h"


/*
Function: fft_init
------------------

Initializes fft arrays and plans
Calculates ranks specific fft parameters
*/
void fft_init(Parameters params, Init_rank_parameters *info, fft_builds *builds) {

    MPI_Comm_size(MPI_COMM_WORLD, &info->numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &info->myid);

    int i;
    int size = info->numprocs;

    //allocating arrays for rank information
    info->sdispls   = malloc(size*sizeof(int));
    info->rdispls   = malloc(size*sizeof(int));
    info->imin      = malloc(size*sizeof(int));
    info->imax      = malloc(size*sizeof(int));
    info->jmin      = malloc(size*sizeof(int));
    info->jmax      = malloc(size*sizeof(int));

    //allocating double to store final time
    builds->time.fft_plan.total_time    = calloc(1,sizeof(double));
    builds->time.run_time.total_time    = calloc(1,sizeof(double));
    builds->time.before_MPI.total_time  = calloc(1,sizeof(double));
    builds->time.during_MPI.total_time  = calloc(1,sizeof(double));
    builds->time.after_MPI.total_time   = calloc(1,sizeof(double));
    builds->time.without_init.total_time   = calloc(1,sizeof(double));


    info->lyr=params.ly/size;
    info->lxr=params.lhx/size;
    info->rcounts=info->lxr*info->lyr;
    info->scounts=info->lxr*info->lyr;
    info->sumscounts=info->scounts*size;
    info->sumrcounts=info->rcounts*size;

    info->imin[0]=1;
    info->imax[0]=info->lxr;
    info->jmin[0]=1;
    info->jmax[0]=info->lyr;
    info->sdispls[0]=0;
    info->rdispls[0]=0;

    for (i = 1; i < size; i++) {
        info->imin[i]=info->imax[i-1]+1;
        info->imax[i]=info->imax[i-1]+info->lxr;
        info->jmin[i]=info->jmax[i-1]+1;
        info->jmax[i]=info->jmax[i-1]+info->lyr;
        info->sdispls[i]=info->sdispls[i-1]+info->scounts;
        info->rdispls[i]=info->rdispls[i-1]+info->rcounts;
    }

    //allocating fft_in_out arrays
    builds->psi.real_arr     = (double*)fftw_malloc(sizeof(double)*(params.lx));
    builds->zero.real_arr    = (double*)fftw_malloc(sizeof(double)*params.ly);
    builds->lxhp.real_arr    = (double*)fftw_malloc(sizeof(double)*params.ly);
    builds->psi.complex_arr  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*params.lxhp);
    builds->zero.complex_arr = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*params.lyhp);
    builds->lxhp.complex_arr = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*params.lyhp);
    builds->psic.f_arr       = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*params.ly);
    builds->psic.b_arr       = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*params.ly);
    builds->send             = (double complex*)fftw_malloc(sizeof(double complex)*info->sumscounts);
    builds->recv             = (double complex*)fftw_malloc(sizeof(double complex)*info->sumrcounts);


    //creating fftw builds
    builds->time.fft_plan.start_time = MPI_Wtime();
    builds->psi.forward_plan     = fftw_plan_dft_r2c_1d(params.lx, builds->psi.real_arr , builds->psi.complex_arr, FFTW_ESTIMATE);
    builds->zero.forward_plan    = fftw_plan_dft_r2c_1d(params.ly, builds->zero.real_arr, builds->zero.complex_arr, FFTW_ESTIMATE);
    builds->lxhp.forward_plan    = fftw_plan_dft_r2c_1d(params.ly, builds->lxhp.real_arr, builds->lxhp.complex_arr, FFTW_ESTIMATE);
    builds->psic.forward_plan    = fftw_plan_dft_1d(params.ly, builds->psic.f_arr, builds->psic.b_arr, FFTW_FORWARD, FFTW_ESTIMATE);
    builds->psic.backward_plan   = fftw_plan_dft_1d(params.ly, builds->psic.b_arr, builds->psic.f_arr, FFTW_BACKWARD, FFTW_ESTIMATE);
    builds->psi.backward_plan    = fftw_plan_dft_c2r_1d(params.lx, builds->psi.complex_arr , builds->psi.real_arr, FFTW_ESTIMATE);
    builds->zero.backward_plan   = fftw_plan_dft_c2r_1d(params.ly, builds->zero.complex_arr, builds->zero.real_arr, FFTW_ESTIMATE);
    builds->lxhp.backward_plan   = fftw_plan_dft_c2r_1d(params.ly, builds->lxhp.complex_arr, builds->lxhp.real_arr, FFTW_ESTIMATE);
    *(builds->time.fft_plan.total_time) += MPI_Wtime() - builds->time.fft_plan.start_time; 
    
    MPI_Barrier(MPI_COMM_WORLD);

    //more rank information
    info->jminmyid  = info->jmin[info->myid];
    info->jmaxmyid  = info->jmax[info->myid];
    info->iminmyid  = info->imin[info->myid];
    info->imaxmyid  = info->imax[info->myid];
    info->j_len = info->jmaxmyid - info->jminmyid+1;
    info->i_len = info->imaxmyid - info->iminmyid+1;

    // printf("%d %d \n",info->j_len,info->i_len);
    // exit(0);

}

/*
Function: fft_forward
------------------

Calcualtes forward fast fourier transform for complete 2d plane
Takes in a 2d input array, and returns 2d output array
*/
void fft_forward(Parameters params, Init_rank_parameters info, fft_builds builds, double *input_arr, Complex_field output_arr) {

    int i, j;
    int ind, np;
    int shift_row;
    double complex tmp_cmplx;

    builds.time.before_MPI.start_time = MPI_Wtime();
    //builds forward for each row
    for ( j = 0; j < info.j_len; j++) {
        shift_row = j*params.lx;
        // use memcpy to copy row from 2d input matrix to fftw input array
        memcpy(builds.psi.real_arr ,input_arr+shift_row ,params.lx*sizeof(double));
        fftw_execute(builds.psi.forward_plan);

        builds.psi.complex_arr[0][1]=builds.psi.complex_arr[params.lxhp-1][0];
        
        for (i = 0; i < params.lhx; i++) {
            ind=i*info.lyr+j;
            builds.send[ind] = builds.psi.complex_arr[i][0] + builds.psi.complex_arr[i][1]*I;
        }
    }
    *(builds.time.before_MPI.total_time) += MPI_Wtime() - builds.time.before_MPI.start_time;

    builds.time.during_MPI.start_time = MPI_Wtime();
    MPI_Alltoall(builds.send, info.scounts, MPI_DOUBLE_COMPLEX, builds.recv, info.rcounts, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);
    *(builds.time.during_MPI.total_time) += MPI_Wtime() - builds.time.during_MPI.start_time;

    builds.time.after_MPI.start_time = MPI_Wtime();
    switch (info.myid) {
    case 0:

        for (np = 0; np < info.numprocs; np++) {
            for (j = info.jmin[np]; j < info.jmax[np]+1; j++) {
                ind = info.rdispls[np]+j-info.jmin[np];
                tmp_cmplx = builds.recv[ind];

                builds.zero.real_arr[j-1]=creal(tmp_cmplx);
                builds.lxhp.real_arr[j-1]=cimag(tmp_cmplx);
            }
        }

        fftw_execute(builds.lxhp.forward_plan);
        fftw_execute(builds.zero.forward_plan);

        memcpy(output_arr.kx_zero , builds.zero.complex_arr, params.lyhp*sizeof(fftw_complex));
        memcpy(output_arr.kx_final, builds.lxhp.complex_arr, params.lyhp*sizeof(fftw_complex));
        
        for (i = 1; i < info.imaxmyid; i++) {
           for (np = 0; np < info.numprocs; np++) {
                for (j = info.jmin[np]; j < info.jmax[np]+1; j++) {
                    ind = info.rdispls[np]+i*info.lyr+(j-info.jmin[np]);
                    tmp_cmplx = builds.recv[ind];
                    builds.psic.f_arr[j-1][0]=creal(tmp_cmplx);
                    builds.psic.f_arr[j-1][1]=cimag(tmp_cmplx);
                }

            }
            fftw_execute(builds.psic.forward_plan);
            shift_row=i*params.ly;
            memcpy(output_arr.rest+shift_row,builds.psic.b_arr,params.ly*sizeof(fftw_complex));
        }

        break;

    default:

        for (i = 0; i < info.i_len; i++) {
            for (np = 0; np < info.numprocs; np++) {
                for (j = info.jmin[np]; j < info.jmax[np]+1; j++) {
                    ind = info.rdispls[np]+i*info.lyr+(j-info.jmin[np]);
                    tmp_cmplx = builds.recv[ind];
                    builds.psic.f_arr[j-1][0]=creal(tmp_cmplx);
                    builds.psic.f_arr[j-1][1]=cimag(tmp_cmplx);
                }
            }

            fftw_execute(builds.psic.forward_plan);

            shift_row=(i)*params.ly;
            memcpy(output_arr.rest+shift_row,builds.psic.b_arr,params.ly*sizeof(fftw_complex));
        }
        break;

    }
    *(builds.time.after_MPI.total_time) += MPI_Wtime() - builds.time.after_MPI.start_time;

}

void fft_backward(Parameters params, Init_rank_parameters info, fft_builds builds, Complex_field input_arr, double* output_arr) {

    int i, j;
    int ind, np;
    int shift_row;
    double complex tmp_cmplx;

    builds.time.before_MPI.start_time = MPI_Wtime();
    switch (info.myid) {
    case 0:
        //printf("test");
        memcpy(builds.zero.complex_arr, input_arr.kx_zero , params.lyhp*sizeof(fftw_complex));
        memcpy(builds.lxhp.complex_arr, input_arr.kx_final, params.lyhp*sizeof(fftw_complex));

        fftw_execute(builds.zero.backward_plan);
        fftw_execute(builds.lxhp.backward_plan);

        for (j = 0; j < params.ly; j++) {
            ind = j*info.lxr;
            builds.recv[ind] = builds.zero.real_arr[j] + builds.lxhp.real_arr[j]*I;
        }

        for (i = 1; i < info.imaxmyid; i++) {

            shift_row = (i)*params.ly;
            memcpy(builds.psic.b_arr,input_arr.rest+shift_row, params.ly*sizeof(fftw_complex));

            fftw_execute(builds.psic.backward_plan);

            for ( j = 0; j < params.ly; j++) {
                ind = j*info.lxr+i;
                tmp_cmplx = builds.psic.f_arr[j][0] + builds.psic.f_arr[j][1]*I;
                builds.recv[ind] = tmp_cmplx;
            }
        }
        break;

    default:

        for (i = 0; i < info.i_len; i++) {

            shift_row = (i)*params.ly;
            memcpy(builds.psic.b_arr,input_arr.rest+shift_row, params.ly*sizeof(fftw_complex));
            fftw_execute(builds.psic.backward_plan);

            for ( j = 0; j < params.ly; j++) {
                ind = j*info.lxr+(i);
                tmp_cmplx = builds.psic.f_arr[j][0] + builds.psic.f_arr[j][1]*I;
                builds.recv[ind] = tmp_cmplx;
            }
        }
        break;
    }
    *(builds.time.before_MPI.total_time) += MPI_Wtime() - builds.time.before_MPI.start_time;

    builds.time.during_MPI.start_time = MPI_Wtime();
    MPI_Alltoall(builds.recv, info.rcounts, MPI_DOUBLE_COMPLEX, builds.send, info.scounts, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);
    *(builds.time.during_MPI.total_time) += MPI_Wtime() - builds.time.during_MPI.start_time;

    builds.time.after_MPI.start_time = MPI_Wtime();
    for ( j = 0; j < info.j_len; j++) {
        for ( np = 0; np < info.numprocs; np++) {
            for ( i = info.imin[np]; i < info.imax[np]+1; i++) {
                ind = info.sdispls[np]+(j)*info.lxr+(i)-info.imin[np];
                tmp_cmplx = builds.send[ind];
                builds.psi.complex_arr[i-1][0] = creal(tmp_cmplx);
                builds.psi.complex_arr[i-1][1] = cimag(tmp_cmplx);
            }
        }

        builds.psi.complex_arr[params.lxhp-1][0] = builds.psi.complex_arr[0][1];
        builds.psi.complex_arr[params.lxhp-1][1] = 0;
        builds.psi.complex_arr[0][1] =  0;

        fftw_execute(builds.psi.backward_plan);

        shift_row = j*params.lx;


        memcpy(output_arr+shift_row,builds.psi.real_arr,params.lx*sizeof(double));
    }
    *(builds.time.after_MPI.total_time) += MPI_Wtime() - builds.time.after_MPI.start_time;

}