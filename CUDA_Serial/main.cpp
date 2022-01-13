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

    const int dx=20,dy=20;
    const double pi = acos(-1);

    system_2D<double, complex<double>> host(dx,dy);
    transform_system_2D<cufftDoubleReal, cufftDoubleComplex> device(host.get_dimensions());
    
    for (auto j = 0; j < dy; j++) {
        for (auto i = 0; i < dx; i++) {
            host(i,j) = cos(2*pi*i/10);
        }
    }

    host.print(1);
    device.forward_transform(host.real(),host.complex());   
    device.inverse_transform(host.complex(),host.real());
    host.print(dx*dy);

    return 0;
}