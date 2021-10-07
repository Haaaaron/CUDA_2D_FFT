#include <stdio.h>
#include <iostream>
#include <tuple>
#include "main.h"



int main(int argc, char const *argv[])
{
    system_2D<int> test_system(100,100);
    test_system(0,0) = 2;
    std::cout << test_system(0,0) << "\n";

    launch_kernel();
    return 0;
}
