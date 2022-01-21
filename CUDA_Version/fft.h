#include <iostream>
#include <stdio.h>
#include <tuple>
#include <cuda.h>
#include <cufft.h>


// https://enccs.github.io/openmp-gpu/gpu-architecture/
#define tile_dim 4

__global__ void transpose(cufftDoubleComplex* input_data, cufftDoubleComplex* output_data, int width, int height);
__global__ void evolve_system_kernel(cufftDoubleComplex *out, cufftDoubleComplex *in, int size);
__global__ void rescale(cufftDoubleReal *out, cufftDoubleReal *in, int size);

using std::cout;
using std::tuple;
using std::get;

class Timing_and_memory_utility{
private:
    float aggregate_transform;
    float aggregate_transpose;

    double free_m,total_m,used_m;
    size_t free_t,total_t;

public:
    cudaEvent_t start_transform, stop_transform;
    cudaEvent_t start_transpose, stop_transpose;

    Timing_and_memory_utility() {
        aggregate_transform = 0;
        aggregate_transpose = 0;

        cudaEventCreate(&start_transform);
        cudaEventCreate(&stop_transform);
        cudaEventCreate(&start_transpose);
        cudaEventCreate(&stop_transpose);
    }

    void aggregate(float transform, float transpose) {
        aggregate_transform += (transform-transpose);
        aggregate_transpose += transpose;
    }

    void print_time(){
        aggregate_transform = aggregate_transform*1.0e-3;
        aggregate_transpose = aggregate_transpose*1.0e-3;
        cout << "Aggregate transform time: " << aggregate_transform << "s\n" //*10e-4
             << "Aggregate transpose time: " << aggregate_transpose << "s\n";//*10e-4
    }

    void gpu_memory_usage(const char* inquiry_text) {

        cudaMemGetInfo(&free_t,&total_t);
        free_m =free_t/1048576.0 ;
        total_m=total_t/1048576.0;
        used_m=total_m-free_m;

        cout << inquiry_text << used_m << '\n';
    }
};

template <typename R, typename C>
class Transform_system_2D {
private:

    cufftHandle to_complex, to_real, complex_to_complex;
    R *real_buffer;
    C *complex_buffer;
    C *complex_buffer_transposed;
    int real_dx, complex_dx, real_dy, complex_dy;
    int threads_per_block = tile_dim*tile_dim;
    int block_x=(complex_dx+tile_dim-1)/tile_dim;
    int block_y=(real_dy+tile_dim-1)/tile_dim;
    float scale = 1/(real_dx*real_dy);

    //timing
    float transform_t;
    float transpose_t;

    //momentary
    C *host_arr_print_complex;
    R *host_arr_print_real;

public:

    Timing_and_memory_utility utility;

    //constructor
    Transform_system_2D(tuple<int,int> dimensions):
    real_dx(get<0>(dimensions)), complex_dx(get<0>(dimensions)/2+1), 
    real_dy(get<1>(dimensions)), complex_dy(get<1>(dimensions)/2+1) {

        utility.gpu_memory_usage("Used memory beofre malloc: ");

        // Allocate device memory
        cudaMalloc((void**)&real_buffer, sizeof(R)*real_dx*real_dy);
        cudaMalloc((void**)&complex_buffer, sizeof(C)*complex_dx*real_dy);
        cudaMalloc((void**)&complex_buffer_transposed, sizeof(C)*complex_dx*real_dy);

        utility.gpu_memory_usage("Used memory after malloc: ");

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

        
        cufftPlanMany(
            &to_real, 1, fft_len_real,
            NULL, 1, complex_dx,
            NULL, 1, real_dx,
            CUFFT_Z2D, real_dy);


        
        cufftPlanMany(
            &complex_to_complex, 1, fft_len_complex,
            NULL, 1, real_dy,
            NULL, 1, real_dy,
            CUFFT_Z2Z, complex_dx/2+1);

        utility.gpu_memory_usage("Used memory after plan creation: ");
                    
    }

    template <typename real>
    void copy_from(real host_array_in) {
        cudaMemcpy(real_buffer,
                    host_array_in,
                    sizeof(R)*real_dx*real_dy,
                    cudaMemcpyHostToDevice);
    }

    template <typename real>
    void copy_to(real host_array_out) {
        cudaMemcpy(host_array_out,
                    real_buffer,
                    sizeof(R)*real_dx*real_dy,
                    cudaMemcpyDeviceToHost);
    }

    template <typename complex>
    void copy_complex_from(complex host_array_in) {
        cudaMemcpy(complex_buffer_transposed, 
                    host_array_in,
                    sizeof(C)*complex_dx*real_dy,
                    cudaMemcpyHostToDevice);
    }

    template <typename complex>
    void copy_complex_to(complex host_array_out) {
        cudaMemcpy(host_array_out, 
                    complex_buffer_transposed,
                    sizeof(C)*complex_dx*real_dy,
                    cudaMemcpyDeviceToHost);
    }

    void forward_transform() {

        transform_t = 0;
        transpose_t = 0;

        //AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        cudaEventRecord(utility.start_transform, 0);
        cufftExecD2Z(to_complex, real_buffer, complex_buffer);

        cudaEventRecord(utility.start_transpose, 0);
        transpose<<<dim3(block_x,block_y), dim3(tile_dim, tile_dim)>>>(complex_buffer,complex_buffer_transposed, complex_dx, real_dy);
        cudaEventRecord(utility.stop_transpose, 0);

        cufftExecZ2Z(complex_to_complex, complex_buffer_transposed, complex_buffer_transposed, CUFFT_FORWARD);
        cudaEventRecord(utility.stop_transform,0);
        //VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV


        cudaEventSynchronize(utility.stop_transpose);
        cudaEventSynchronize(utility.stop_transform);

        cudaEventElapsedTime(&transform_t, utility.start_transform, utility.stop_transform);
        cudaEventElapsedTime(&transpose_t, utility.start_transpose, utility.stop_transpose);

        utility.aggregate(transform_t, transpose_t);
        //cudaDeviceSynchronize();
    }

    void inverse_transform() {
        transform_t = 0;
        transpose_t = 0;

        cudaEventRecord(utility.start_transform, 0);
        cufftExecZ2Z(complex_to_complex, complex_buffer_transposed, complex_buffer_transposed, CUFFT_INVERSE);
       
        cudaEventRecord(utility.start_transpose, 0);
        transpose<<<dim3(block_y,block_x), dim3(tile_dim, tile_dim)>>>(complex_buffer_transposed,complex_buffer, real_dy, complex_dx);
        cudaEventRecord(utility.stop_transpose, 0);
       
        cufftExecZ2D(to_real, complex_buffer, real_buffer);
        cudaEventRecord(utility.stop_transform,0);

        //rescale<<<(real_dx*real_dy+1024-1)/1024,1024>>>(real_buffer,real_buffer,real_dx*real_dy);

        cudaEventSynchronize(utility.stop_transpose);
        cudaEventSynchronize(utility.stop_transform);

        cudaEventElapsedTime(&transform_t, utility.start_transform, utility.stop_transform);
        cudaEventElapsedTime(&transpose_t, utility.start_transpose, utility.stop_transpose);

        utility.aggregate(transform_t, transpose_t);
        //cudaDeviceSynchronize();

    }

    // Physical system evolution
    void evolve_system() {
        evolve_system_kernel<<<dim3(block_y,block_x), dim3(tile_dim, tile_dim)>>>(complex_buffer_transposed, complex_buffer_transposed, complex_dx*real_dy);
    }

    void print_time() {
        utility.print_time();
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

        