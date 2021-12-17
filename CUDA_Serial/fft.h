#include <tuple>
#include <cuda.h>
#include <cufft.h>
// https://enccs.github.io/openmp-gpu/gpu-architecture/
#define tile_dim 4

__global__ void transpose(cufftDoubleComplex* input_data, cufftDoubleComplex* output_data, int width, int height);

template <typename R, typename C>
class transform_system_2D {
    private:
        cufftHandle to_complex, to_real, complex_to_complex;
        R *real_buffer;
        C *complex_buffer;
        C *complex_buffer_transposed;

        int real_dx, complex_dx, real_dy, complex_dy;
        int threads_per_block = tile_dim*tile_dim;
        int block_x=(complex_dx+tile_dim-1)/tile_dim;
        int block_y=(real_dy+tile_dim-1)/tile_dim;
        
        //momentary
        C *host_arr_print_complex;
        R *host_arr_print_real;

    public:
        transform_system_2D(tuple<int,int> dimensions):
        real_dx(get<0>(dimensions)), complex_dx(get<0>(dimensions)/2+1), 
        real_dy(get<1>(dimensions)), complex_dy(get<1>(dimensions)/2+1) {
            
            // Allocate device memory
            cudaMalloc((void**)&real_buffer, sizeof(R)*real_dx*real_dy);
            cudaMalloc((void**)&complex_buffer, sizeof(C)*complex_dx*real_dy);
            cudaMalloc((void**)&complex_buffer_transposed, sizeof(C)*complex_dx*real_dy);

            host_arr_print_complex = new C[complex_dx*real_dy];
            host_arr_print_real = new R[real_dx*real_dy];

            int fft_len_real[1] = {real_dx};
            int fft_len_complex[1] = {real_dy};
            // Create plans
            cufftPlanMany(
                &to_complex, 1, fft_len_real,
                NULL, 1, real_dx,
                NULL, 1, complex_dx,
                CUFFT_D2Z, real_dy);
            
            //TODO add excpetion iif cufft plan creation fails
            //if (res != CUFFT_SUCCESS) {throw exceptin()}

            cufftPlanMany(
                &to_real, 1, fft_len_real,
                NULL, 1, complex_dx,
                NULL, 1, real_dx,
                CUFFT_Z2D, real_dy);
            
            // cufftPlanMany(
            //     &complex_to_complex, 1, len_c,
            //     NULL, 1, real_dy,
            //     NULL, 1, real_dy,
            //     CUFFT_Z2Z, complex_dx);
            
            cufftPlanMany(
                &complex_to_complex, 1, fft_len_complex,
                NULL, 1, real_dy,
                NULL, 1, real_dy,
                CUFFT_Z2Z, complex_dx/2+1);
        }

        template <typename real, typename complex>
        void forward_transform(real host_array_in, complex host_array_out) {

            cudaMemcpy(real_buffer, 
                       host_array_in, 
                       sizeof(R)*real_dx*real_dy, 
                       cudaMemcpyHostToDevice);

            cufftExecD2Z(to_complex, real_buffer, complex_buffer);
            //CUDA kernels asynchronous, timing needs to be done differently 
            transpose<<<dim3(block_x,block_y), dim3(tile_dim, tile_dim)>>>(complex_buffer,complex_buffer_transposed, complex_dx, real_dy);
            cufftExecZ2Z(complex_to_complex, complex_buffer_transposed, complex_buffer_transposed, CUFFT_FORWARD);
            cudaMemcpy(host_array_out, 
                       complex_buffer_transposed,
                       sizeof(C)*complex_dx*real_dy,
                       cudaMemcpyDeviceToHost);
        }

        template <typename complex, typename real>
        void inverse_transform(complex host_array_in, real host_array_out) {

            cudaMemcpy(complex_buffer_transposed, host_array_in,
                       sizeof(C)*complex_dx*real_dy,
                       cudaMemcpyHostToDevice);
                       
            cufftExecZ2Z(complex_to_complex, complex_buffer_transposed, complex_buffer_transposed, CUFFT_INVERSE);
            transpose<<<dim3(block_y,block_x), dim3(tile_dim, tile_dim)>>>(complex_buffer_transposed,complex_buffer, real_dy, complex_dx);
            cufftExecZ2D(to_real, complex_buffer, real_buffer);

            cudaMemcpy(host_array_out, real_buffer,
                      sizeof(R)*real_dx*real_dy,
                      cudaMemcpyDeviceToHost);
        }     

        void print_complex(){        
            cudaMemcpy(host_arr_print_complex, complex_buffer,
                        sizeof(C)*(complex_dx)*real_dy,
                        cudaMemcpyDeviceToHost);
            cout << "Complex buffer:" << "\n";
            double x, y;
            for (auto i = 0; i < real_dy; i++) {
                for (auto j = 0; j < complex_dx; j++) {
                    if (fabs(host_arr_print_complex[i*complex_dx+j].x) < 10e-6) {
                        x = 0;
                    } else {
                        x = host_arr_print_complex[i*complex_dx+j].x;
                    }
                    if (fabs(host_arr_print_complex[i*complex_dx+j].y) < 10e-6) {
                        y = 0;
                    } else {
                        y = host_arr_print_complex[i*complex_dx+j].y;
                    }
                    cout << x << " + I*" << y << "   ";
                }
                cout << "\n";
            }
            
            cout << "\n";
        }

        void print_complex_T(){        
            cudaMemcpy(host_arr_print_complex, complex_buffer_transposed,
                        sizeof(C)*(complex_dx)*real_dy,
                        cudaMemcpyDeviceToHost);
            cout << "Complex buffer T:" << "\n";
            double x, y;
            for (auto i = 0; i < complex_dx; i++) {
                for (auto j = 0; j < real_dy; j++) {
                    if (fabs(host_arr_print_complex[i*real_dy+j].x) < 10e-6) {
                        x = 0;
                    } else {
                        x = host_arr_print_complex[i*real_dy+j].x;
                    }
                    if (fabs(host_arr_print_complex[i*real_dy+j].y) < 10e-6) {
                        y = 0;
                    } else {
                        y = host_arr_print_complex[i*real_dy+j].x;
                    }
                    cout << x << " + I*" << y << "   ";
                }
                cout << "\n";
            }
            
            cout << "\n";
            
        }
        
        void print_real(){        
            cudaMemcpy(host_arr_print_real, real_buffer,
                        sizeof(R)*(real_dx)*real_dy,
                        cudaMemcpyDeviceToHost);
            cout << "Real buffer results:" << "\n";
            for (auto i = 0; i < real_dy; i++) {
                for (auto j = 0; j < real_dx; j++) {
                    cout << host_arr_print_real[i*real_dx+j] << "   ";
                }
                cout << "\n";
            }
            cout << "\n";
        }
};