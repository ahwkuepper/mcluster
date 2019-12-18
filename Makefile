OPTIMIZATION=-O3
FFLAGS=$(HEADERS:%=-I%) -fbounds-check -mcmodel=medium -fopenmp -fPIC
LIBS=-lm -lgsl -lgslcblas -lc -lstdc++
FC=gfortran $(OPTIMIZATION)

default: all

all: clean mcluster

mcluster:
	$(FC) inipar/dictionary.c inipar/iniparser.c inipar/iniparser_interface.c input.f main.c $(LIBS) $(FFLAGS) -o mcluster
#gfortran -O3 inipar/dictionary.c inipar/iniparser.c inipar/iniparser_interface.c input.f main.c -lm -lgsl -lgslcblas -fbounds-check -mcmodel=medium -fopenmp -o mcluster
clean:
	rm --f mcluster
