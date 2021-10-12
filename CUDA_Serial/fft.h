#include <tuple>
#include <cuda.h>
#include <cufft.h>

class transform_system_2D {
    private:
        cufftHandle to_complex, to_real;
        cufftReal *real_buffer;
        cufftComplex *complex_buffer;
        int array_len, batch_size;

    public:
        __host__ transform_system_2D(tuple<int,int> dimensions):
        array_len(get<0>(dimensions)), batch_size(get<1>(dimensions)) {
            // Allocate device memory
            cudaMalloc((void**)&real_buffer, sizeof(cufftReal)*array_len*batch_size);
            cudaMalloc((void**)&complex_buffer, sizeof(cufftComplex)*array_len*batch_size);

            int len[1] = {array_len};
            // Create plans
            cufftPlanMany(
                &to_complex, 1, len,
                NULL, 1, array_len,
                NULL, 1, array_len,
                CUFFT_R2C, batch_size);
            
            //TODO add excpetion iif cufft plan creation fails
            //if (res != CUFFT_SUCCESS) {printf("Plan creation fail\n"); return 1}

            cufftPlanMany(
                &to_real, 1, len,
                NULL, 1, array_len,
                NULL, 1, array_len,
                CUFFT_C2R, batch_size);
        }

        __host__ void forward_transform(double* host_array) {

            cudaMemcpy(real_buffer, host_array, 
                       sizeof(cufftReal)*array_len*batch_size, 
                       cudaMemcpyHostToDevice);

            cufftExecR2C(to_complex, real_buffer, complex_buffer);
        }

        __host__ void backward_transform(double* host_array) {

            cufftExecC2R(to_real, complex_buffer, real_buffer);

            cudaMemcpy(host_array, real_buffer,
                       sizeof(cufftReal)*array_len*batch_size,
                       cudaMemcpyDeviceToHost);
        }

        
};