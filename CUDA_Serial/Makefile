objects = main.o FFT.o

all: $(objects)
	nvcc -arch=sm_86 -lcufft $(objects) -o app

%.o: %.cpp
	nvcc -x cu -arch=sm_86 -I. -dc $< -o $@

clean:
	rm -f *.o app