#include <stdio.h>
#include <tuple>
#include "main.h"
#include "fft.h"

__global__ void transpose(cufftDoubleComplex *idata, cufftDoubleComplex *odata, int width, int height) {
	__shared__ cufftDoubleComplex block[tile_dim][tile_dim+1];
	
	// read the matrix tile into shared memory
    // load one element per thread from device memory (idata) and store it
    // in transposed order in block[][]
	unsigned int xIndex = blockIdx.x * tile_dim + threadIdx.x;
	unsigned int yIndex = blockIdx.y * tile_dim + threadIdx.y;
	if((xIndex < width) && (yIndex < height))
	{
		unsigned int index_in = yIndex * width + xIndex;
		block[threadIdx.y][threadIdx.x] = idata[index_in];
	}

    // synchronise to ensure all writes to block[][] have completed
	__syncthreads();

	// write the transposed matrix tile to global memory (odata) in linear order
	xIndex = blockIdx.y * tile_dim + threadIdx.x;
	yIndex = blockIdx.x * tile_dim + threadIdx.y;
	if((xIndex < height) && (yIndex < width))
	{
		unsigned int index_out = yIndex * height + xIndex;
		odata[index_out] = block[threadIdx.x][threadIdx.y];
	}
}

__global__ void evolve_system_kernel(cufftDoubleComplex *out, cufftDoubleComplex *in, int size) {
	
	unsigned int ind=threadIdx.x+blockIdx.x*blockDim.x;
	cufftDoubleComplex tmp = in[ind];
	out[ind].x = tmp.x*0.25;
	out[ind].y = tmp.y*0.25;	
}

__global__ void rescale(cufftDoubleReal *out, cufftDoubleReal *in, int size) {
	
	float scale_factor = 1.0/size;
	unsigned int ind=threadIdx.x+blockIdx.x*blockDim.x;
	out[ind] = in[ind]*scale_factor;
}
