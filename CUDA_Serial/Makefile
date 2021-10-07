CC = g++
CUDACC = nvcc
CUDAFLAGS = -arch=sm_86 -c
LIBDIRS=-L/usr/local/cuda-11.4/lib64

FFT.o: FFT.cu
	$(CUDACC) $(CUDAFLAGS) FFT.cu
main.o: main.cpp
	$(CC) -c main.cpp

main: main.o FFT.o
	$(CC) -o serial_cuda main.o FFT.o $(LIBDIRS) -lcudart

clean:
	rm -rf serial_cuda *.o