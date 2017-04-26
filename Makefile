# GNU Fortran Compiler
F90    = gfortran
FFLAGS = -O0 -g -fconvert=big-endian -fno-automatic -fbacktrace \
         -ffpe-trap=invalid,zero,overflow,underflow,denormal    \
         -fPIC -fcheck=all -Wall -Wextra -fimplicit-none        \
         -Wuninitialized -pedantic -Warray-temporaries
#FFLAGS = -O3 -march=native -funroll-loops

## Intel Fortran Compiler
#F90    = ifort
#FFLAGS = -g -C -traceback -mkl

snsol.exe : main.f90
	$(F90) $(FFLAGS) -o snsol.exe main.f90

clean:
	rm -f snsol.exe *.mod *.o

