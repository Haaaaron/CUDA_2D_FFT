#include <stdio.h>
#include <tuple>
#include "main.h"
#include "fft.h"

// __global__ void transpose(cufftDoubleComplex* input_data, cufftDoubleComplex* output_data, int width, int height) {
    
//     int index = blockIdx.x*blockDim.x + threadIdx.x;
//     int tile_y = threadIdx.x/tile_dim;
//     int tile_x = threadIdx.x % tile_dim;
//     int x_index = index/width;
//     int y_index = index % width;
//     int out_index = y_index*height + x_index;
//     __shared__ cufftDoubleComplex tile[tile_dim][tile_dim+1];

//     if (index<width*height) tile[tile_y][tile_x]=input_data[index];

//     __syncthreads();

//     if(out_index < height*width) output_data[out_index]=tile[tile_x][tile_y];
// }

__global__ void transpose(cufftDoubleComplex* in, cufftDoubleComplex* out, int width, int height) {
    
    int x_index= blockIdx.x * tile_dim + threadIdx.x;
    int y_index= blockIdx.y * tile_dim + threadIdx.y;
    __shared__ cufftDoubleComplex tile[tile_dim][tile_dim+1];
    if(y_index<height && x_index<width) {
        int in_index=y_index*width + x_index; //[y][x]
        tile[threadIdx.y][threadIdx.x]=in[in_index];
        // Coalesced reads
    }
    __syncthreads();
    int new_x = blockIdx.y * tile_dim + threadIdx.x; // transpose block offset
    int new_y = blockIdx.x * tile_dim + threadIdx.y;
    if(y_index<width && x_index<height) {
        int out_index=new_y*width + new_x; //[x][y] different indices. Correct????
        // Coalesced writes
        // (Based on https://developer.nvidia.com/blog/efficient-matrix-transpose-cuda-cc/)
        out[out_index]=tile[threadIdx.x][threadIdx.y];
    }
}

void forward_fft(system_2D<double>& host_system) {
    transform_system_2D<cufftDoubleReal, cufftDoubleComplex> device_system(host_system.get_dimensions());
    device_system.forward_transform(host_system.get_data());
    device_system.inverse_transform(host_system.get_data());
    host_system.print();
    cudaDeviceSynchronize();
}