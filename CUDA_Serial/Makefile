objects = main.o FFT.o

all: $(objects)
	nvcc -arch=sm_86 -lcufft $(objects) -o serial_cuda

%.o: %.cpp
	nvcc -x cu -arch=sm_86 -I. -dc $< -o $@


FFT.o: FFT.cu
	nvcc -arch=sm_86 -c $< -lcufft

clean:
	rm -f *.o serial_cuda