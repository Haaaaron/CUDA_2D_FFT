#include <stdio.h>
#include <tuple>
#include "main.h"
#include "fft.h"

__global__ void transpose() {
    printf("Hello from mykernel\n");
}

void forward_fft(system_2D<double>& host_system) {
    transform_system_2D device_system(host_system.get_dimensions());

    device_system.forward_transform(host_system.get_data());
    device_system.backward_transform(host_system.get_data());

    host_system.print();
    
    cudaDeviceSynchronize();
}