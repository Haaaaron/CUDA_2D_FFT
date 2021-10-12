#include <stdio.h>
#include <iostream>
#include <tuple>
#include "main.h"

using namespace std;

int main(int argc, char const *argv[]) {
    system_2D<double> test_system(10,10);

    for (auto i = 0; i < 10; i++) {
        for (auto j = 0; j < 10; j++) {
            test_system(i,j) = 1;
        }
    }
    
    test_system.print();
    forward_fft(test_system);
    return 0;
}
