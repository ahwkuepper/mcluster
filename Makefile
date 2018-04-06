# Compiler
FC = gfortran -O2 #-ffast-math
#FC = gfortran-mp-4.4 -O2 #-ffast-math
#use of -ffast-math may cause problems on some machines

INCPATH=include

CC = gcc -O2 -fopenmp -Wall #-ffast-math
#use "-D NOOMP" and remove "-fopenmp" for compilation without OpenMP

CFLAGS = -L/usr/lib/ -lgfortran
CUFLAGS= -O3 -D WITH_CUDA5 -L$(CUDA_PATH)/lib64 -lcudart -lstdc++ -lm $(CFLAGS)
CUDA_PATH = $(shell which nvcc | rev | cut -d'/' -f3- | rev)
SDK_PATH=$(CUDA_PATH)/samples

SOURCES := $(shell find . -type f -name '*.f')
INCLUDES := $(shell find . -type f -name '*.h')
OBJECTS := $(SOURCES:.f=.o)

%.o:%.f
	$(FC) -c $^ -o $@

# Default
all: mcluster

mcluster_sse: $(OBJECTS)
	$(CC) -c main.c -D SSE -lm
	$(CC) $(OBJECTS) main.o -o mcluster_sse -lm $(CFLAGS)

mcluster_gpu: $(OBJECTS) $(INCPATH)/gpupot.gpu.o
	$(CC) -c main.c -D SSE -D GPU -lm -I$(CUDA_PATH)/include
	$(CC) $(OBJECTS) main.o gpupot.gpu.o  -o mcluster_gpu $(CUFLAGS)

mcluster:
	@echo $(OBJECTS)
	$(CC) -o mcluster main.c -lm

$(INCPATH)/gpupot.gpu.o: $(INCPATH)/gpupot.gpu.cu $(INCPATH)/cuda_pointer.h
	nvcc -c $(CUFLAGS) -Xcompiler "-fPIC -O3 -Wall" -I$(SDK_PATH)/common/inc -I$(INCPATH) $(INCPATH)/gpupot.gpu.cu

clean:
	rm -f $(INCPATH)/*.o *.o mcluster_sse mcluster mcluster_gpu
