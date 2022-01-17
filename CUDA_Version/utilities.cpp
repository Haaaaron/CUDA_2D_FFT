#include "main.h"

/*
init_arguments
--------------
Reads in command line arguments:
{dx} {dy} {iteration_count}
*/
void init_arguments(int argc, char const* argv[], int* dx, int* dy, int* iteration_count) {
    switch(argc) {
        case 1:
            *dx = 10;
            *dy = 10;
            *iteration_count = 1;
            break;
        
        case 4:
            *dx = atoi(argv[1]);
            *dy = atoi(argv[2]);
            *iteration_count = atoi(argv[3]);
            break;

        default:
            cout << "unsupported number of command line arguments \n"
                 << "Example usage: ./fft_benchmark {dx} {dy} {iteration_count} \n";
            exit(-1);
    }
}