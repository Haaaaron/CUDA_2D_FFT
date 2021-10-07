#include <stdio.h>
#include "main.h"

__global__ void mykernel() {
    printf("Hello from mykernel\n");
}

void launch_kernel() {
    mykernel<<<1,1>>>();
    cudaDeviceSynchronize();
}