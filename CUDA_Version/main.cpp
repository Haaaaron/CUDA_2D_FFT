#include "main.h"
#include "fft.h"

int main(int argc, char const *argv[]) {

    int dx; 
    int dy;
    int iteration_count;
    const double pi = acos(-1);

    init_arguments(argc, argv, &dx, &dy, &iteration_count);

    //Allocating host arrays
    System_2D<double, complex<double>> host(dx,dy);

    //Timing and initializing fft plans and allocating cuda arrays
    auto start_plans_t = high_resolution_clock::now();
    Transform_system_2D<cufftDoubleReal, cufftDoubleComplex> device(host.get_dimensions());
    auto stop_plans_t = high_resolution_clock::now();
    auto duration_plans = duration_cast<milliseconds>(stop_plans_t - start_plans_t);

    //Initializing system
    for (auto j = 0; j < dy; j++) {
        for (auto i = 0; i < dx; i++) {
            //host(i,j) = cos(2*pi*i/10.0);
            host(i,j) = 1;
        }
    }

    device.copy_from(host.real());
    host.print_R();

    auto start = high_resolution_clock::now();
    for (auto j = 0; j<iteration_count ; j++) {
        device.forward_transform();
        device.inverse_transform();
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);

    device.copy_to(host.real());
    
    cout << "Total run without init: " << duration.count()*1.0e-3 << " s\n";
    cout << "Time fft plan creation: " << duration_plans.count()*1.0e-3 << " s\n";
    device.utility.print_time();

    // double sum = 0;
    // for (auto j = 0; j < dy; j++) {
    //     for (auto i = 0; i < dx; i++) {
    //         //host(i,j) = cos(2*pi*i/10);
    //         sum += host(i,j);
    //     }
    // }

    // //host.print(1);
    // cout << sum << '\n';

    return 0;
}
