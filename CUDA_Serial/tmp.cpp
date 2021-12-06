#include "fft.h"

int main(int argc, char const *argv[])
{
    transform_system_2D<cufftDoubleReal, cufftDoubleComplex> device(host.get_dimensions());

    return 0;
}
