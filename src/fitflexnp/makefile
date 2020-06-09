FC = gfortran
FLAGS = -O2 -Bstatic 

# invert using nnls with bathymetry and gravity data, solve for 3 constants and fault slip rates
fitflex:
	$(FC) $(FLAGS) fitflexnp.f nnls.f -o ../../bin/fitflexnp -framework veclib

