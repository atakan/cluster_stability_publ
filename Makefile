CC = gcc
CPP = g++
FC = gfortran
CFLAGS = -march=native
#ALL_CFLAGS = -Wall -O2 -fpeel-loops -fomit-frame-pointer -ffast-math -DNDEBUG -fopenmp $(CFLAGS)
#ALL_CFLAGS = -Wall -O0 -g -fopenmp $(CFLAGS)
ALL_CFLAGS = -Wall -O2 $(CFLAGS)

## icc: this is copied over from another project, not tested
#CC = /opt/intel/Compiler/11.1/064/bin/ia32/icc
#CFLAGS = -fno-alias -xSSE4.2 -fomit-frame-pointer -O2 -funroll-loops -openmp
#ALL_CFLAGS = -Wall -wd810,869,981,1418,1572 -DNDEBUG $(CFLAGS)
#LIBS = libnpy.a -L/opt/intel/Compiler/11.1/064/lib/ia32 -lgsl -lgslcblas -lm

all: init main

#pattern rule to compile object files from C, C++ and Fortran files
%.o : %.c Makefile Star.h
	$(CPP) $(ALL_CFLAGS) -c $< -o $@
%.o : %.cpp Makefile Star.h
	$(CPP) $(ALL_CFLAGS) -c $< -o $@
%.o : %.f Makefile Star.h
	$(FC) $(ALL_CFLAGS) -c $< -o $@

init: init.o Makefile Star.h
	$(CPP) $(ALL_CFLAGS) init.o -o init $(LIBS)

main: main.o Makefile Star.h
	$(CPP) $(ALL_CFLAGS) main.o -o main $(LIBS)

clean:
	rm -f init.o main.o Star.o init main
