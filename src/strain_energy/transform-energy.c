#include "litho.h"

void set_litho_defaults_(Litho *);
double yield_stress_(Litho *,)

int main (int argc, char *argv[])
{
  int i,j,nz,nx;
  double dz,dx,z,x,phi;

  unsigned int wcsw, tesw, hssw;

  Litho l;
  Litho *lptr = &l; /* pointer to litho struct */

  // check for input errors

  set_litho_defaults_(lptr);

  /* 1. Compute all the yield stresses in the transform zone */
  for(i=0;i<nz;i++){
    for(j=0;j<nx;j++){
      z =;
      x =;
      ystr[i][j] = yield_stress_(lptr,&z,&x,&wcsw,&tesw,&hssw);

    }
  }

  /* 2. Integrate the yield stress times (full) spreading rate over transform zone */
  phi = 0; // energy dissipated by transform
  for(i=0;i<nz;i++){
    for(j=0;j<nx;j++){
      phi = phi + ystr[i][j]*2*lptr->usp*dx*dz;
    }
  }

}
