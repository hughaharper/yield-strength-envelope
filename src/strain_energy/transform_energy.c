#include "litho.h"

void set_litho_defaults_(Litho *);
double yield_stress_(Litho *, double *, double *, unsigned int *, unsigned int *, unsigned int *);

int main (int argc, char *argv[])
{
  int i,j,nz,nx;
  double age_offset, spr_rate, spr_rate_ms;
  double dz,dx,z,x_a,x_b,phi;
  double ystr_a,ystr_b;
  double zmax=25000,xmax;

  double *ystr;

  unsigned int wcsw=1, tesw=0, hssw=0;

  Litho l;
  Litho *lptr = &l; /* pointer to litho struct */

  // check for input errors
  if(argc < 3){
    fprintf(stderr,"\n Usage ocean_litho_yse age (Myr) rate (m/yr) (wcsw 0/1 tesw 0/1 hssw 0/1) \n");
    exit(-1);
  }

  // inputs should be 1) age offset 2) spreading rate and 3) remaining flags
  age_offset = atof(argv[1]);
  spr_rate = atof(argv[2]);

  if(argc > 3){
    i = 3;
    while(i < argc){
      /* string match with flags */
      if(strncmp(argv[i],"wcsw",4) == 0){
        wcsw = atoi(argv[i+1]);
      }
      else if(strncmp(argv[i],"tesw",4) == 0) {
        tesw = atoi(argv[i+1]);
      }
      else if(strncmp(argv[i],"hssw",4) == 0) {
        hssw = atoi(argv[i+1]);
      }
      else {
        fprintf(stderr,"Not a valid flag\n");
        fprintf(stderr,"\n Usage ocean_litho_yse age (Myr) rate (m/yr) (wcsw 0/1 tesw 0/1 hssw 0/1) \n");
        exit(-1);
      }
      i+=2;
    }
  }

  set_litho_defaults_(lptr);
  lptr->usp = spr_rate;
  spr_rate_ms = spr_rate/365.25/24/60/60;


  dx = 100; // meters
  dz = 100;
  xmax = age_offset*1e6*spr_rate;
  nx = (int) (xmax / dx);
  nz = (int) (zmax / dz);
  fprintf(stderr,"Depth: %.1lf, Length: %.1lf \n",zmax,xmax);
  fprintf(stderr,"nz: %d, nx: %d\n",nz,nx);

  ystr = (double*)malloc(nx*nz*sizeof(double));


  /* 1. Compute all the yield stresses in the transform zone
  and ntegrate the yield stress times (full) spreading rate over transform zone
  */
  phi = 0;
  for(i=0;i<nz;i++){
    z = ((double)i+0.5)*dz;
    for(j=0;j<nx;j++){
      x_a =((double)j+0.5)*dx;
      x_b = xmax - x_a;
      // Need to compute 2 yield stresses, for either side of transform...
      ystr_a = yield_stress_(lptr,&z,&x_a,&wcsw,&tesw,&hssw);
      ystr_b = yield_stress_(lptr,&z,&x_b,&wcsw,&tesw,&hssw);
      ystr[i*nx + j] = fmax(ystr_a,ystr_b);
      phi = phi + -1*ystr[i*nx + j]*2*spr_rate_ms*dx*dz;
    }
  }

  fprintf(stdout,"%.2lf\n", phi);
}
