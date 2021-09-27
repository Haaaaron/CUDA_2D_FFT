#include "mpi.h"
#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h> 

#include "main.h"

void coeff_init(Parameters params, Init_rank_parameters info, Physical_system *system) {

    int i, j, ind;
    double x,y;
    double qq, kkx, kky, ccc;
    double qq_zero, kkx_zero, ccc_zero, psi_zero;
    double sum_psi, intermediate;

    //allocating 2d arrays for physical system
    system->psi                 = malloc(params.lx*info.j_len*sizeof(double));
    system->nt                  = malloc(params.lx*info.j_len*sizeof(double));
    system->real_buffer         = malloc(params.lx*info.j_len*sizeof(double));
    system->final_energy        = malloc(params.lx*info.j_len*sizeof(double));
    system->complex_nt.rest     = malloc(params.ly*info.i_len*sizeof(fftw_complex));
    system->complex_psi.rest    = malloc(params.ly*info.i_len*sizeof(fftw_complex));
    system->complex_buffer.rest = malloc(params.ly*info.i_len*sizeof(fftw_complex));


    //allocating coefficient arrays
    system->coeff.energy.lxhp   = malloc(params.lyhp*sizeof(double));
    system->coeff.energy.zero   = malloc(params.lyhp*sizeof(double));
    system->coeff.psi.lxhp      = malloc(params.lyhp*sizeof(double));
    system->coeff.psi.zero      = malloc(params.lyhp*sizeof(double));
    system->coeff.ntc.lxhp      = malloc(params.lyhp*sizeof(double));
    system->coeff.ntc.zero      = malloc(params.lyhp*sizeof(double));
    system->coeff.psi.rest      = malloc(params.ly*info.i_len*sizeof(double));
    system->coeff.ntc.rest      = malloc(params.ly*info.i_len*sizeof(double));
    system->coeff.energy.rest   = malloc(params.ly*info.i_len*sizeof(double));

    //allocating complex field zeroth case arrays

    system->complex_buffer.kx_final = malloc(params.lyhp*sizeof(fftw_complex));
    system->complex_buffer.kx_zero  = malloc(params.lyhp*sizeof(fftw_complex));
    system->complex_nt.kx_final     = malloc(params.lyhp*sizeof(fftw_complex));
    system->complex_nt.kx_zero      = malloc(params.lyhp*sizeof(fftw_complex));
    system->complex_psi.kx_final    = malloc(params.lyhp*sizeof(fftw_complex));
    system->complex_psi.kx_zero     = malloc(params.lyhp*sizeof(fftw_complex));



    //psi array
    for (j = 0; j < info.j_len; j++) {
        for ( i = 0; i < params.lx; i++) {
            ind = (j)*params.lx+i;
            x = (i+1)*params.dx;
            y = (j+info.jminmyid)*params.dy;
            intermediate =-0.13*(cos(sqrt(3.0)*x/2.0)*cos(y/2.0)-0.5*cos(y))+params.pm;
            sum_psi += intermediate;
            system->psi[ind] = intermediate;

        }
    }

    intermediate = sum_psi/(params.lx*info.lyr)-params.pm;
    for (j = 0; j < info.j_len; j++) {
        for ( i = 0; i < params.lx; i++) {
            ind = j*params.lx+i;

            system->psi[ind] = system->psi[ind] - intermediate;
        }
    }

    //First row has fourier coefficients that are of power 0 i.e. equating to 1
    switch (info.myid) {
    case 0:
        //zeroth case
        for (j = 0; j < params.lyhp; j++) {

            kky=(2*params.pi*j/(params.ly*params.dy));

            kkx_zero=0;
            qq_zero=-(kkx_zero*kkx_zero+kky*kky);
            ccc_zero=params.r+1+2*qq_zero+qq_zero*qq_zero;

            kkx=params.pi/params.dx;
            qq=-(kkx*kkx+kky*kky);
            ccc=params.r+1+2*qq+qq*qq;

            system->coeff.energy.zero[j]=params.VV*ccc_zero;
            system->coeff.psi.zero[j]=1.0*params.VV/(1-params.dt*qq_zero*ccc_zero);
            system->coeff.ntc.zero[j]=params.dt*qq_zero*params.VV/(1-params.dt*qq_zero*ccc_zero);

            system->coeff.energy.lxhp[j]=params.VV*ccc;
            system->coeff.psi.lxhp[j]=1.0*params.VV/(1-params.dt*qq*ccc);
            system->coeff.ntc.lxhp[j]=params.dt*qq*params.VV/(1-params.dt*qq*ccc);
        }
        
        for (i = 1; i < info.i_len; i++) {

            kkx = 2*params.pi*(i)/(params.lx*params.dx);

            for (j = 0; j < params.ly; j++) {

                ind = (i)*params.ly+j;

                if (j < params.lyhp) {
                    kky=2*params.pi*(j)/(params.ly*params.dy);

                } else {
                    kky=2*params.pi*(j-params.ly)/(params.ly*params.dy);

                }
                
                qq=-(kkx*kkx+kky*kky);
                ccc=params.r+1+2*qq+qq*qq;
                system->coeff.energy.rest[ind]=params.VV*ccc;
                system->coeff.psi.rest[ind]=1.0*params.VV/(1-params.dt*qq*ccc);
                //rdeb(system->coeff.psi.rest[ind],ind);
                system->coeff.ntc.rest[ind]=qq*params.dt*params.VV/(1-params.dt*qq*ccc);

            }
        }
        //exit(0);
        break;

    default:

        for (i = info.iminmyid-1; i < info.imaxmyid; i++) {

            kkx = 2*params.pi*(i)/(params.lx*params.dx);

            for (j = 0; j < params.ly; j++) {

                ind = (i-info.iminmyid+1)*params.ly+j;

                if (j < params.lyhp) {
                    kky=2*params.pi*(j)/(params.ly*params.dy);

                } else {
                    kky=2*params.pi*(j-params.ly)/(params.ly*params.dy);

                }
                
                qq=-(kkx*kkx+kky*kky);
                ccc=params.r+1+2*qq+qq*qq;
                system->coeff.energy.rest[ind]=params.VV*ccc;
                system->coeff.psi.rest[ind]=1.0*params.VV/(1-params.dt*qq*ccc);
                system->coeff.ntc.rest[ind]=qq*params.dt*params.VV/(1-params.dt*qq*ccc);

            }
        }
        break;
    }
}

double energy(Parameters params, Init_rank_parameters info, Physical_system system, fft_builds builds) {

    double energy;
    int i, j, ind;

    fft_forward(params, info, builds, system.psi, system.complex_buffer);

    switch (info.myid) {
    case 0:

        for (i = 0; i < params.lyhp; i++) {
            system.complex_buffer.kx_zero[i][0]  = system.complex_buffer.kx_zero[i][0]*system.coeff.energy.zero[i];
            system.complex_buffer.kx_zero[i][1]  = system.complex_buffer.kx_zero[i][1]*system.coeff.energy.zero[i];
            system.complex_buffer.kx_final[i][0] = system.complex_buffer.kx_final[i][0]*system.coeff.energy.lxhp[i];
            system.complex_buffer.kx_final[i][1] = system.complex_buffer.kx_final[i][1]*system.coeff.energy.lxhp[i];
        }

        for (i = 1; i < info.i_len; i++) {
            for ( j = 0; j < params.ly; j++) {
                ind = i*params.ly+j;
                system.complex_buffer.rest[ind][0] = system.complex_buffer.rest[ind][0]*system.coeff.energy.rest[ind];
                system.complex_buffer.rest[ind][1] = system.complex_buffer.rest[ind][1]*system.coeff.energy.rest[ind];
            }
        }
        
        break;

    default:

        for (i = 0; i < info.i_len; i++) {
            for ( j = 0; j < params.ly; j++) {
                ind = i*params.ly+j;

                system.complex_buffer.rest[ind][0] = system.complex_buffer.rest[ind][0]*system.coeff.energy.rest[ind];
                system.complex_buffer.rest[ind][1] = system.complex_buffer.rest[ind][1]*system.coeff.energy.rest[ind];
            }
        }

        break;

    }

    fft_backward(params, info, builds, system.complex_buffer, system.real_buffer);

    energy=0;
    for ( i = 0; i < params.lx*info.j_len; i++) {
        energy += (1.0/4.0)*pow(system.psi[i],4)+(1.0/2.0)*system.psi[i]*system.real_buffer[i];
    }
    energy = energy*params.VV;

    MPI_Reduce(&energy, &system.total_energy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    return system.total_energy;
}

double* update(Parameters params, Init_rank_parameters info, Physical_system system, fft_builds builds) {
    
    int i,j,ind;
    double tmp_real, tmp_imag;

    for ( i = 0; i < params.lx*info.j_len; i++)
        system.nt[i] = pow(system.psi[i],3);

    fft_forward(params, info, builds, system.nt, system.complex_nt);
    fft_forward(params, info, builds, system.psi, system.complex_psi);
    switch (info.myid) 
    {
    case 0:
        for ( i = 0; i < params.lyhp; i++) {
            system.complex_psi.kx_zero[i][0] = system.coeff.psi.zero[i]*system.complex_psi.kx_zero[i][0] 
                + system.coeff.ntc.zero[i]*system.complex_nt.kx_zero[i][0];
            system.complex_psi.kx_zero[i][1] = system.coeff.psi.zero[i]*system.complex_psi.kx_zero[i][1] 
                + system.coeff.ntc.zero[i]*system.complex_nt.kx_zero[i][1];

            system.complex_psi.kx_final[i][0] = system.coeff.psi.lxhp[i]*system.complex_psi.kx_final[i][0] 
                + system.coeff.ntc.lxhp[i]*system.complex_nt.kx_final[i][0];
            system.complex_psi.kx_final[i][1] = system.coeff.psi.lxhp[i]*system.complex_psi.kx_final[i][1] 
                + system.coeff.ntc.lxhp[i]*system.complex_nt.kx_final[i][1];
        }
        
        for ( i = 1; i < info.imaxmyid; i++) {
            for ( j = 0; j < params.ly; j++) {   
                ind = (i)*params.ly+j;
                tmp_real = system.coeff.psi.rest[ind]*system.complex_psi.rest[ind][0] 
                    + system.coeff.ntc.rest[ind]*system.complex_nt.rest[ind][0];
                tmp_imag = system.coeff.psi.rest[ind]*system.complex_psi.rest[ind][1] 
                    + system.coeff.ntc.rest[ind]*system.complex_nt.rest[ind][1];
                
                system.complex_psi.rest[ind][0] = tmp_real;
                system.complex_psi.rest[ind][1] = tmp_imag;
            }            
        }

        break;
    
    default:
        for ( i = 0; i < info.i_len; i++) {
            for ( j = 0; j < params.ly; j++) {   
                ind = i*params.ly+j;
                tmp_real = system.coeff.psi.rest[ind]*system.complex_psi.rest[ind][0] 
                    + system.coeff.ntc.rest[ind]*system.complex_nt.rest[ind][0];
                tmp_imag = system.coeff.psi.rest[ind]*system.complex_psi.rest[ind][1] 
                    + system.coeff.ntc.rest[ind]*system.complex_nt.rest[ind][1];
                
                system.complex_psi.rest[ind][0] = tmp_real;
                system.complex_psi.rest[ind][1] = tmp_imag;

            }
        }

        break;
    }    
    fft_backward(params, info, builds, system.complex_psi, system.psi);
    return system.psi;
}
