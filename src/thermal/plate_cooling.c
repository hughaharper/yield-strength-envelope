/*
* compute plate cooling thermal model for a given plate age
*/

#include "litho.h"

void set_litho_defaults_(Litho *);
double temp_plt_(Litho *, double *, double *);

int main (int argc, char **argv)
{
  int i,j,nz=1000,nx=100;
  double age,d_age,z,dz,zp,temp,max_age;

  Litho l;
  Litho *lptr = &l;

  if(argc < 2){
    printf("\n Usage plate_cooling max_age flags\n");
    exit(-1);
  }

  max_age = atof(argv[1]);
  set_litho_defaults_(lptr);

  zp = lptr->dp;
  dz = zp/nz;
  d_age = max_age/nx;

  /* print a header line of depts */
  fprintf(stdout,"0 ");
  for(j=0;j<nz;j++) {
    z=((double)j+0.5)*dz;
    fprintf(stdout,"%lf ",z);
  }

  /* loop thru age and depth */
  for(i=0;i<nx;i++) {
    age = ((double)i+0.5)*d_age;

    fprintf(stdout,"\n%lf ",age);
    for(j=0;j<nz;j++) {

      z=((double)j+0.5)*dz;
      temp=temp_plt_(lptr,&z,&age);

      fprintf(stdout,"%lf ",temp);
    }

  }
  return(EXIT_SUCCESS);

}
