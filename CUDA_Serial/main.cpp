#include <stdio.h>
#include <iostream>
#include <tuple>
#include <cuda_runtime.h>
#include <complex>
#include <math.h>
#include "main.h"
#include "fft.h"

template <typename T> std::string type_name();

using std::complex;
using std::cout;

int main(int argc, char const *argv[]) {

    int dx=12,dy=30;
    system_2D<double, complex<double>> host(dx,dy);
    transform_system_2D<cufftDoubleReal, cufftDoubleComplex> device(host.get_dimensions());
    
    for (auto j = 0; j < dy; j++) {
        for (auto i = 0; i < dx; i++) {
            host(j,i) = 0;
        }
    }
    
    host(dx/2,dy/2) = 1

    host.print(1);
    device.forward_transform(host.real(),host.complex());   
    device.inverse_transform(host.complex(),host.real());
    host.print(dx*dy);

    cout << typeid(host.real()).name() << "\n";
    return 0;
}