IGNORE:
FC = gfortran -O2 #-ffast-math
#FC = gfortran-mp-4.4 -O2 #-ffast-math
#use of -ffast-math may cause problems on some machines

LFLAGS = const_bse.h zdata.h 

CC = gcc -O2 -fopenmp -Wall #-ffast-math 
#use "-D NOOMP" and remove "-fopenmp" for compilation without OpenMP

#CFLAGS = -L/opt/local/lib/gcc44/ -lgfortran
CFLAGS = -L/usr/lib/ -lgfortran
CUFLAGS= -O3 -D WITH_CUDA5
CUDA_PATH = /usr/local/cuda
SDK_PATH=$(CUDA_PATH)/samples

.f.o:
	$(FC) -c $<

SOURCE = \
deltat.f evolv1.f hrdiag.f kick.f mlwind.f mrenv.f \
ran3.f star.f zcnsts.f zfuncs.f\
comenv.f corerd.f dgcore.f evolv2.f gntage.f \
instar.f mix.f rl.f

OBJECTS = $(SOURCE:.f=.o)

mcluster_sse: $(OBJECTS) $(LFLAGS)
	$(CC) -c main.c -D SSE -lm
	$(CC) $(OBJECTS) main.o -o mcluster_sse -lm $(CFLAGS) 

mcluster_ssegpu: $(OBJECTS) $(LFLAGS) gpupot.gpu.o main.c
	$(CC) -c main.c -D SSE -D GPU -lm -I$(CUDA_PATH)/include
	$(CC) $(OBJECTS) main.o gpupot.gpu.o  -o mcluster_gpu -L$(CUDA_PATH)/lib64 -lcudart -lstdc++ -lm $(CFLAGS) 

mcluster_gpu: gpupot.gpu.o main.c
	$(CC) -c main.c -D GPU -I$(CUDA_PATH)/include
	$(CC) main.o gpupot.gpu.o -o mcluster_gpu -L$(CUDA_PATH)/lib64 -lcudart -lstdc++ -lm

mcluster: 
	$(CC) -o mcluster main.c -lm

gpupot.gpu.o: gpupot.gpu.cu cuda_pointer.h
	nvcc -c $(CUFLAGS) -Xcompiler "-fPIC -O3 -Wall" -I$(SDK_PATH)/common/inc -I. gpupot.gpu.cu 

clean:
	rm -f *.o mcluster_sse mcluster mcluster_gpu
