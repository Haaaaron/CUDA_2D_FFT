#include "mpi.h"
#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>

#include "main.h"

int main(int argc, char** argv) {
    
    int size_x;
    int size_y;
    int iteration_count;
    double pi=acos(-1.0), sq3=sqrt(3), at=7.255197456936871;
    //clock_t start, end, start_without_init;
    double start, end, start_without_init;

    init_params(argc, argv, &size_x, &size_y, &iteration_count);

    Init_rank_parameters    init_information;
    fft_builds              builds;
    Physical_system         system;

    Parameters params = {
        .lx     = size_x, 
        .ly     = size_y, 
        .lxhp   = size_x/2+1,
        .lyhp   = size_y/2+1, 
        .lhx    = size_x/2, 
        .lhy    = size_y/2,
        .pi     = pi, 
        .sq3    = sq3, 
        .at     = at, 
        .qqt    = 2*pi/at,
        .dt     = 0.125,
        .VV     = 1.0/(size_x*size_y),
        .dx     = at/8.0,
        .dy     = (at*sq3/2.0)/8.0,
        .r      = -0.25,
        .pm     = -0.25
    };

    //MPI Barriers so that run time is global between all ranks
    MPI_Init(&argc, &argv);
    MPI_Barrier(MPI_COMM_WORLD);
    start  = MPI_Wtime();

    //Initialization
    fft_init(params, &init_information, &builds);
    coeff_init(params, init_information, &system);

    MPI_Barrier(MPI_COMM_WORLD);
    start_without_init = MPI_Wtime();

    //Main simulation
    run_simulation(params, init_information, builds, system, iteration_count);

    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();

    if(init_information.myid==0) {
        *(builds.time.run_time.total_time) = end - start;
        *(builds.time.without_init.total_time) = end - start_without_init;
    }

    output_time(builds.time, init_information);

    MPI_Finalize();

    return 0;
}

void run_simulation(Parameters params, Init_rank_parameters info, fft_builds builds, Physical_system system, int iteration_count) {
    
    int iteration;
    system.total_energy=0;

    for ( iteration = 0; iteration < iteration_count; iteration++)
        system.psi =  update(params, info, system, builds);
    
    system.total_energy = energy(params, info, system, builds);
    if (info.myid == 0) printf("Total energy after simulation: %.16e\n", system.total_energy);

}

