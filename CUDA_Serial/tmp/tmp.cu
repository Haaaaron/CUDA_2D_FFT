#include <stdio.h>
#include <stdlib.h>
#include <cufft.h>
#include <assert.h>

/********************/
/* CUDA ERROR CHECK */
/********************/
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, char *file, int line, bool abort=true)
{
    if (code != cudaSuccess) 
    {
        fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) { getchar(); exit(code); }
    }
}

/*********************/
/* CUFFT ERROR CHECK */
/*********************/
static const char *_cudaGetErrorEnum(cufftResult error)
{
    switch (error)
    {
        case CUFFT_SUCCESS:
            return "CUFFT_SUCCESS";

        case CUFFT_INVALID_PLAN:
            return "CUFFT_INVALID_PLAN";

        case CUFFT_ALLOC_FAILED:
            return "CUFFT_ALLOC_FAILED";

        case CUFFT_INVALID_TYPE:
            return "CUFFT_INVALID_TYPE";

        case CUFFT_INVALID_VALUE:
            return "CUFFT_INVALID_VALUE";

        case CUFFT_INTERNAL_ERROR:
            return "CUFFT_INTERNAL_ERROR";

        case CUFFT_EXEC_FAILED:
            return "CUFFT_EXEC_FAILED";

        case CUFFT_SETUP_FAILED:
            return "CUFFT_SETUP_FAILED";

        case CUFFT_INVALID_SIZE:
            return "CUFFT_INVALID_SIZE";

        case CUFFT_UNALIGNED_DATA:
            return "CUFFT_UNALIGNED_DATA";
    }

    return "<unknown>";
}

#define cufftSafeCall(err)      __cufftSafeCall(err, __FILE__, __LINE__)
inline void __cufftSafeCall(cufftResult err, const char *file, const int line)
{
    if( CUFFT_SUCCESS != err) {
                fprintf(stderr, "CUFFT error in file '%s', line %d\n %s\nerror %d: %s\nterminating!\n",__FILE__, __LINE__,err, \
                           _cudaGetErrorEnum(err)); \
             cudaDeviceReset(); assert(0); \
    }
}

/********/
/* MAIN */
/********/
int main() {

    cufftHandle forward_plan, inverse_plan; 

    int batch = 3;
    int rank = 1;

    int nRows = 1;
    int nCols = 10;
    int n[1] = {nCols};

    int idist = nRows*nCols;
    int odist = nRows*(nCols/2+1);

    int inembed[] = {nCols};
    int onembed[] = {nCols/2+1};

    int istride = 1;
    int ostride = 1;

    cufftSafeCall(cufftPlanMany(&forward_plan,  rank, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_D2Z, batch));

    cufftDoubleReal *h_in = (cufftDoubleReal*)malloc(sizeof(cufftDoubleReal)*nRows*nCols*batch);
    for(int i=0; i<nRows*nCols*batch; i++) h_in[i] = 1.f;

    cufftDoubleComplex* h_freq = (cufftDoubleComplex*)malloc(sizeof(cufftDoubleComplex)*nRows*(nCols/2+1)*batch);

    cufftDoubleReal* d_in;            gpuErrchk(cudaMalloc(&d_in, sizeof(cufftDoubleReal)*nRows*nCols*batch)); 
    cufftDoubleComplex* d_freq; gpuErrchk(cudaMalloc(&d_freq, sizeof(cufftDoubleComplex)*nRows*(nCols/2+1)*batch)); 

    gpuErrchk(cudaMemcpy(d_in,h_in,sizeof(cufftDoubleReal)*nRows*nCols*batch,cudaMemcpyHostToDevice));

    cufftSafeCall(cufftExecD2Z(forward_plan, d_in, d_freq));

    gpuErrchk(cudaMemcpy(h_freq,d_freq,sizeof(cufftDoubleComplex)*nRows*(nCols/2+1)*batch,cudaMemcpyDeviceToHost));

    for(int i=0; i<nRows*(nCols/2+1)*batch; i++) printf("Direct transform: %i %f %f\n",i,h_freq[i].x,h_freq[i].y); 

    cufftSafeCall(cufftPlanMany(&inverse_plan, rank, n, onembed, ostride, odist, inembed, istride, idist, CUFFT_Z2D, batch));

    cufftSafeCall(cufftExecZ2D(inverse_plan, d_freq, d_in));

    gpuErrchk(cudaMemcpy(h_in,d_in,sizeof(cufftDoubleReal)*nRows*nCols*batch,cudaMemcpyDeviceToHost));

    for(int i=0; i<nRows*nCols*batch; i++) printf("Inverse transform: %i %f \n",i,h_in[i]); 

    getchar();

}