PRO=false		#1 to enable profiling

FC=mpif90
CC=mpicc
FCFLAGS=-O3 -fbounds-check -fopenmp -march=native -fdefault-real-8 -fdefault-double-8 -lfftw3 #-lfftw3_omp # -I/appl/spack/v014/install-tree/gcc-9.3.0/fftw-3.3.8-kml7ii/include
CCFLAGS=-O3 -w -std=c99 -g -march=native -lfftw3 -lm -I/home/Software/fftw-3.3.10/api/#-lfftw3_omp # -I/appl/spack/v014/install-tree/gcc-9.3.0/fftw-3.3.8-kml7ii/include

MPI_FFTW_Bench: src_c/main.c
	$(CC) -o MPI_FFTW_Bench src_c/main.c src_c/fft.c src_c/phys.c src_c/utilities.c ${CCFLAGS}


