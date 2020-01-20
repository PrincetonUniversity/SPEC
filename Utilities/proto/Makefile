.PHONY: spec

all: spec

spec: spec.f90 spec_input.f90
	gfortran -ffree-line-length-300 spec_input.f90 -I/usr/include -c
	gfortran -ffree-line-length-300 spec.f90 -I/usr/include -c
	gfortran spec_input.o spec.o -lhdf5_fortran -lhdf5 -lm -lz -lpthread -o spec

test:
	./spec demo_nml.sp

clean:
	rm -f *.o *.mod spec
