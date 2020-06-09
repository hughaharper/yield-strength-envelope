#!/bin/csh	Makefile for flex2dc


NETCDF = $(NETCDFHOME)/lib
GMTLIB = $(GMTHOME)/lib

FFT     = ../../lib/libflex.a
FFLAGS= -O -Bstatic -m64
FC = gfortran

flex2dc:	flex2dc.o
	$(FC) $(FFLAGS) flex2dc.o	\
		$(FFT) -L$(GMTLIB) -lgmt -L$(NETCDF) -lnetcdf -o ../../bin/flex2dc
	rm *.o
