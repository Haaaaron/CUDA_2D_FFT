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

// __global__ void convolution_kernel(cufftDoubleComplex *idata, const cufftDoubleReal dx, int x, int y) {

//     unsigned int xIndex = blockIdx.x * gridDim.x + threadIdx.x;
//     int ix = xIndex/y;
//     int iy = xIndex % y;
//     double kx = 2*acos(-1)/(x*dx)*(ix > x/2 ? ix - x: ix);
//     double ky = 2*acos(-1)/(y*dy)*(iy);

//     double kx_square = 2*(1-cos(dx*kx))/(dx*dx);
//     double ky_square = 2*(1-cos(dy*ky))/(dy*dy);

//     double loc_x, loc_y;

//     if(xIndex < x*y) {
        
//         cufftDoubleComplex pppk=idata[ind];

//         double klapl=-(kx_square+ky_square);


//         pppk.x=klapl*(pppk.x);
//         pppk.y=klapl*(pppk.y);

//         idata[ind]=pppk;
//     }


// }