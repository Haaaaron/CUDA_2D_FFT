PRO=false		#1 to enable profiling

CC=mpicc
CCFLAGS=-O3 -w -std=c99 -g -march=native -lfftw3 -lm -I/home/Software/fftw-3.3.10/api/#-lfftw3_omp # -I/appl/spack/v014/install-tree/gcc-9.3.0/fftw-3.3.8-kml7ii/include

EXECUTABLE = mpi_version

$(EXECUTABLE): src_c/main.c
	$(CC) -o mpi_version src_c/main.c src_c/fft.c src_c/phys.c src_c/utilities.c ${CCFLAGS}


