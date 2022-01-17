#include "main.h"

void init_params(int argc, char* argv[], int* size_x, int* size_y, int* iteration_count) {

    switch (argc)
    {
    case 1:
        *size_x = 128;
        *size_y = 256;
        *iteration_count = 100;
        break;
    
    case 4:
        *size_x = atoi(argv[1]);
        *size_y = atoi(argv[2]);
        *iteration_count = atoi(argv[3]);
        break;

    default:
        printf("Unsupported number of command line arguments\n");
        printf("Example usage: ./TwoDMPIPFC {size_x} {size_y} {iteration_count}");
        exit(-1);

    }
}

void output_time(struct Time time, Init_rank_parameters info) {

    typedef struct {
        double avg;
        double min;
        double max;
    } final_time;

    final_time fft_plan;
    final_time before_MPI;
    final_time during_MPI;
    final_time after_MPI ;
    
    MPI_Reduce(time.fft_plan.total_time  , &(fft_plan.avg)  , 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(time.before_MPI.total_time, &(before_MPI.avg), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(time.during_MPI.total_time, &(during_MPI.avg), 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(time.after_MPI.total_time , &(after_MPI.avg) , 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    MPI_Reduce(time.fft_plan.total_time  , &(fft_plan.min)  , 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(time.before_MPI.total_time, &(before_MPI.min), 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(time.during_MPI.total_time, &(during_MPI.min), 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(time.after_MPI.total_time , &(after_MPI.min) , 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    
    MPI_Reduce(time.fft_plan.total_time  , &(fft_plan.max)  , 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(time.before_MPI.total_time, &(before_MPI.max), 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(time.during_MPI.total_time, &(during_MPI.max), 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(time.after_MPI.total_time , &(after_MPI.max) , 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (info.myid == 0) {

        fft_plan.avg   /= info.numprocs;
        before_MPI.avg /= info.numprocs; 
        during_MPI.avg /= info.numprocs;
        after_MPI.avg  /= info.numprocs;
        printf("Total run time [s]: %.16e\n", *(time.run_time.total_time));
        printf("Total run time without init [s]: %.16e\n", *(time.without_init.total_time));
        printf("Time fft plan creation, Average[s]: %lf Min[s]: %lf Max[s]: %lf\n", fft_plan.avg  , fft_plan.min  , fft_plan.max  );
        printf("Time before MPI,        Average[s]: %lf Min[s]: %lf Max[s]: %lf\n", before_MPI.avg, before_MPI.min, before_MPI.max);
        printf("Time during MPI,        Average[s]: %lf Min[s]: %lf Max[s]: %lf\n", during_MPI.avg, during_MPI.min, during_MPI.max);
        printf("Time after  MPI,        Average[s]: %lf Min[s]: %lf Max[s]: %lf\n", after_MPI .avg, after_MPI .min, after_MPI .max);
    }
}