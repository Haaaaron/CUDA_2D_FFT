#include "main.h"
#include "fft.h"

int main(int argc, char const *argv[]) {

    int dx; 
    int dy;
    int iteration_count;
    const double pi = acos(-1);

    init_arguments(argc, argv, &dx, &dy, &iteration_count);

    System_2D<double, complex<double>> host(dx,dy);

    auto start_plans_t = high_resolution_clock::now();
    Transform_system_2D<cufftDoubleReal, cufftDoubleComplex> device(host.get_dimensions());
    auto stop_plans_t = high_resolution_clock::now();
    auto duration_plans = duration_cast<microseconds>(stop_plans_t - start_plans_t);

    for (auto j = 0; j < dy; j++) {
        for (auto i = 0; i < dx; i++) {
            //host(i,j) = cos(2*pi*i/10);
            host(i,j) = 1;
        }
    }
    
    //host.print(1);

    device.copy_from(host.real());

    auto start = high_resolution_clock::now();
    for (auto j = 0; j<iteration_count ; j++) {
        device.forward_transform();
        //device.evolve_system();
        device.inverse_transform();

    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);

    device.copy_to(host.real());
    
    cout << "Total run without init: " << duration.count()*10e-7 << " s\n";
    cout << "Time fft plan creation: " << duration_plans.count()*10e-7 << " s\n";
    device.time.print_time();

    return 0;
}
