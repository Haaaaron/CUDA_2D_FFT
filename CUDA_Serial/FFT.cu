#include <stdio.h>
#include <tuple>
#include "main.h"
#include "fft.h"

__global__ void mykernel() {
    printf("Hello from mykernel\n");
}

void forward_fft(system_2D<double>& host_system) {
    transform_system_2D device_system(host_system.get_dimensions());
    mykernel<<<1,1>>>();
    cudaDeviceSynchronize();
}