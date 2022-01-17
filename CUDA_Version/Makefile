
SOURCES = main.cpp FFT.cu utilities.cpp
OBJECTS = main.o utilities.o FFT.o
EXECUTABLE = cuda_version

$(EXECUTABLE): $(OBJECTS)
	nvcc -arch=sm_86 -lcufft $(OBJECTS) -o $@

%.o: %.cpp
	nvcc -x cu -arch=sm_86 -I. -dc $< -o $@

%.o: %.cu
	nvcc -arch=sm_86 -c $< -lcufft

main.o: fft.h main.h
all:
	@echo $(OBJECTS)
clean:
	rm -f *.o cuda_version