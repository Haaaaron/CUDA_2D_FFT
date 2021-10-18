#include <stdio.h>
#include <iostream>
#include <tuple>
#include "main.h"
#include "fft.h"
using namespace std;

int main(int argc, char const *argv[]) {

    int dx=10,dy=10;
    system_2D<double> host(dx,dy);
    transform_system_2D<cufftDoubleReal, cufftDoubleComplex> device(host.get_dimensions());

    for (auto i = 0; i < dx; i++) {
        for (auto j = 0; j < dy; j++) {
            host(i,j) = 1;
        }
    }

    //for (int& n : host_system.get_data()){n = 1;};
    
    
    host.print();
    device.forward_transform(host.get_data());
    device.inverse_transform(host.get_data());
    //forward_fft(host_system);
    return 0;
}

