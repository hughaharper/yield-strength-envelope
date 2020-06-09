#!/bin/csh	Makefile for vload

FFT     = ../../lib/libflex.a
FFLAGS= -O -Bstatic -m64
FC = gfortran

vload:	vload.o
	$(FC) $(FFLAGS) vload.o	\
		$(FFT) -L$(GMTHOME)/lib -lgmt -L$(NETCDFHOME)/lib -lnetcdf -o ../../bin/vload
	rm *.o
