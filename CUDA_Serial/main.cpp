#include <stdio.h>
#include <iostream>
#include <tuple>
#include <complex>
#include <math.h>
#include "main.h"
#include "fft.h"

template <typename T> std::string type_name();

using std::complex;
using std::cout;

int main(int argc, char const *argv[]) {

    const int dx=30000,dy=30000;
    const double pi = acos(-1);
    const int steps = 100;

    system_2D<double, complex<double>> host(dx,dy);
    transform_system_2D<cufftDoubleReal, cufftDoubleComplex> device(host.get_dimensions());
    
    for (auto j = 0; j < dy; j++) {
        for (auto i = 0; i < dx; i++) {
            host(i,j) = cos(2*pi*i/10);
        }
    }
    
    //host.print(1);
    device.copy_real_from_host(host.real());

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start);
    for (auto j = 0; j<steps ; j++) {
        device.forward_transform();
        device.evolve_system();
        device.inverse_transform();
    }
    cudaEventRecord(stop);

    device.copy_real_to_host(host.real());
    cout << host.real()[0] << "\n";

    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);

    cout << "Transform took "
         << milliseconds
         << " milliseconds\n";
    return 0;
}
