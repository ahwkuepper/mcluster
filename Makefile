IGNORE:
FC = gfortran -O2 #-ffast-math
#FC = gfortran-mp-4.4 -O2 #-ffast-math
#use of -ffast-math may cause problems on some machines

LFLAGS = const_bse.h zdata.h 

CC = gcc -O2 -fopenmp #-ffast-math 
#use "-D NOOMP" and remove "-fopenmp" for compilation without OpenMP

#CFLAGS = -L/opt/local/lib/gcc44/ -lgfortran
CFLAGS = -L/usr/lib/ -lgfortran

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
	$(CC) $(CFLAGS) $(OBJECTS) main.o -o mcluster_sse -lm 

mcluster: 
	$(CC) -o mcluster main.c -lm

clean:
	rm *.o mcluster_sse mcluster
