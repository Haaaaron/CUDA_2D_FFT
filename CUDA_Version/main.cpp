#include "main.h"
#include "fft.h"

int main(int argc, char const *argv[]) {

    const int dx=12000,dy=12000;
    const double pi = acos(-1);
    const int steps = 10;

    System_2D<double, complex<double>> host(dx,dy);
    Transform_system_2D<cufftDoubleReal, cufftDoubleComplex> device(host.get_dimensions());
    
    for (auto j = 0; j < dy; j++) {
        for (auto i = 0; i < dx; i++) {
            //host(i,j) = cos(2*pi*i/10);
            host(i,j) = 1;
        }
    }
    
    //host.print(1);

    device.copy_from(host.real());

    auto start = high_resolution_clock::now();
    for (auto j = 0; j<steps ; j++) {
        device.forward_transform();
        //device.evolve_system();
        device.inverse_transform();

    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    device.copy_to(host.real());
    
    cout << "Total time:" << duration.count()*10e-7 << " s\n";
    device.time.print_time();

    return 0;
}