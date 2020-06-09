FC = gfortran
FLAGS = -O2 -Bstatic 

# invert using nnls with gps, insar and constriant matrix and generate uncertainties, solve for 3 constants and fault slip rates
fitflex:
	$(FC) $(FLAGS) fitflexsmt.f nnls.f -o ../../bin/fitflexsmt -framework veclib

