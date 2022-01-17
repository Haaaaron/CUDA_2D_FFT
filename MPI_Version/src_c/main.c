#include "main.h"

int main(int argc, char** argv) {
    
    int size_x;
    int size_y;
    int iteration_count;
    double pi=acos(-1.0);
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
        .dt     = 0.125,
    };

    //MPI Barriers so that run time is global between all ranks
    MPI_Init(&argc, &argv);
    MPI_Barrier(MPI_COMM_WORLD);
    start  = MPI_Wtime();

    //Initialization
    fft_init(params, &init_information, &builds);
    buffer_init(params, init_information, &system);

    MPI_Barrier(MPI_COMM_WORLD);
    start_without_init = MPI_Wtime();

    //Main benchmark
    run_benchmark(params, init_information, builds, system, iteration_count);

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

void run_benchmark(Parameters params, Init_rank_parameters info, fft_builds builds, Physical_system system, int iteration_count) {
    
    printf("\n");
    for (int iteration = 0; iteration < iteration_count; iteration++) {
        fft_forward(params, info, builds, system.real_buffer, system.complex_buffer);
        fft_backward(params, info, builds, system.complex_buffer, system.real_buffer);
        // printf("Iteration: %d\n",iteration);
        // for (int j = 0; j < info.j_len ; j++) {
        //     for (int i = 0; i < params.lx ; i++) {
        //         printf("%.4e ", system.real_buffer[(j)*params.lx+i]);
        //     }
        //     printf("\n");
        // }
        // printf("\n");
    }
}

