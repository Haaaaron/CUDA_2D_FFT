#include <tuple>
#include <cuda.h>
#include <cufft.h>

#define tile_dim 32

__global__ void transpose(cufftDoubleComplex* input_data, cufftDoubleComplex* output_data, int width, int height);

template <typename R, typename C>
class transform_system_2D {
    private:
        cufftHandle to_complex, to_real, complex_to_complex;
        R *real_buffer;
        C *complex_buffer;
        C *complex_buffer_transposed;
        int array_len, batch_size;
        int threads_per_block = tile_dim*tile_dim;
        int block_x=((array_len/2+1)+tile_dim-1)/tile_dim;
        int block_y=(batch_size+tile_dim-1)/tile_dim;

        //momentary
        C *host_arr_print_complex;
        R *host_arr_print_real;

    public:
        transform_system_2D(tuple<int,int> dimensions):
        array_len(get<0>(dimensions)), batch_size(get<1>(dimensions)) {
            // Allocate device memory
            cudaMalloc((void**)&real_buffer, sizeof(R)*array_len*batch_size);
            cudaMalloc((void**)&complex_buffer, sizeof(C)*(array_len/2+1)*batch_size);
            cudaMalloc((void**)&complex_buffer_transposed, sizeof(C)*(array_len/2+1)*batch_size);

            host_arr_print_complex = new C[(array_len/2+1)*batch_size];
            host_arr_print_real = new R[(array_len)*batch_size];

            int len_r[1] = {array_len};
            int len_c[1] = {batch_size};
            // Create plans
            cufftPlanMany(
                &to_complex, 1, len_r,
                NULL, 1, array_len,
                NULL, 1, (array_len/2+1),
                CUFFT_D2Z, batch_size);
            
            //TODO add excpetion iif cufft plan creation fails
            //if (res != CUFFT_SUCCESS) {throw exceptin()}

            cufftPlanMany(
                &to_real, 1, len_r,
                NULL, 1, (array_len/2+1),
                NULL, 1, array_len,
                CUFFT_Z2D, batch_size);
            
            cufftPlanMany(
                &complex_to_complex, 1, len_c,
                NULL, 1, batch_size,
                NULL, 1, batch_size,
                CUFFT_Z2Z, (array_len/2+1));
        }

        void forward_transform(R* host_array) {

            cudaMemcpy(real_buffer, host_array, 
                       sizeof(R)*array_len*batch_size, 
                       cudaMemcpyHostToDevice);
            cufftResult res = cufftExecD2Z(to_complex, real_buffer, complex_buffer);
            transpose<<<dim3(block_x,block_y), dim3(tile_dim, tile_dim)>>>(complex_buffer,complex_buffer_transposed, array_len/2+1, batch_size);
            cufftExecZ2Z(complex_to_complex, complex_buffer_transposed, complex_buffer_transposed, CUFFT_FORWARD);
        }

        void inverse_transform(double* host_array) {
            cufftExecZ2Z(complex_to_complex, complex_buffer_transposed, complex_buffer_transposed, CUFFT_INVERSE);
            transpose<<<dim3(block_x,block_y), dim3(tile_dim, tile_dim)>>>(complex_buffer_transposed,complex_buffer, array_len/2+1, batch_size);
            (*this).print_complex();
            cufftExecZ2D(to_real, complex_buffer, real_buffer);
            cudaMemcpy(host_array, real_buffer,
                      sizeof(R)*array_len*batch_size,
                      cudaMemcpyDeviceToHost);
        }     

        void print_complex(){        
            cudaMemcpy(host_arr_print_complex, complex_buffer,
                        sizeof(C)*(array_len/2+1)*batch_size,
                        cudaMemcpyDeviceToHost);
            cout << "Complex buffer results:" << "\n";
            for (auto i = 0; i < 10; i++) {
                for (auto j = 0; j < 10; j++) {
                    cout << host_arr_print_complex[i*10+j].x << " + I*" << host_arr_print_complex[i*6+j].y << "   ";
                }
                cout << "\n";
            }
            cout << "\n";
        }
        
        void print_real(){        
            cudaMemcpy(host_arr_print_real, real_buffer,
                        sizeof(R)*(array_len)*batch_size,
                        cudaMemcpyDeviceToHost);
            cout << "Real buffer results:" << "\n";
            for (auto i = 0; i < batch_size; i++) {
                for (auto j = 0; j < array_len; j++) {
                    cout << host_arr_print_real[i*array_len+j] << "   ";
                }
                cout << "\n";
            }
            cout << "\n";
        }
};