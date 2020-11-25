/*
 * compute temperature of oceanic lith using Morton and Sleep (1985) model
*/

#include "litho.h"

void set_litho_defaults_(Litho *);
double temp_sleep_(Litho *, double *, double *, unsigned int *);

int main (int argc, char **argv)
{
  /* variables */
  int i,j,nz,nx;
  double x,dx,z,dz,u,age;
  double T_all;
  unsigned int hssw;

  Litho l;
  Litho *lptr = &l;

  if(argc < 2){
    printf("\n Usage sleep_cooling half_spreading_rate hssw (0/1)\n");
    exit(-1);
  }

  u = atof(argv[1])/365/24/60/60; /* rate in m/yr, convert to m/s */
  hssw = atoi(argv[2]);
  set_litho_defaults_(lptr);
  lptr->usp = atof(argv[1]); // this needs to be stored as m/yr in the litho structure
  dz = 100; /* z spacing */
  dx = 100; /* x spacing */
  nz = (int) 30000/dz; /* model space only 30km by 30km, no reason to go to 125 km depth */
  nx = (int) 30000/dx;
  fprintf(stderr,"nz: %d, nx: %d\n",nz,nx);

  /* ------------------------------------------------------------------------ */

  fprintf(stdout,"0 ");
  for(j=0;j<nz;j++) {
    z=((double)j)*dz;
    fprintf(stdout,"%lf ",z);
  }

  for(j=0;j<nx;j++) {
    /* loop through x distance */
    x = ((double) j + 0.5)*dx;
    age = (x/u)/SPMYR;

    fprintf(stdout,"\n%lf ",x);

    for(i=0;i<nz;i++) {
      /* loop thru depths */
      z = ((double) i + 0.5)*dz;

      T_all = temp_sleep_(lptr,&z,&age,&hssw);

      fprintf(stdout,"%lf ",T_all);
    } /* end depth loop */
  } /* end age (distance) loop */
  return(EXIT_SUCCESS);

} /* main */
