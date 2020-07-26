/*******************************************************************************
* Compute yield strength envelopes for oceanic lithosphere of variable age *
*******************************************************************************/

#include "litho.h"

void set_litho_defaults_(Litho *);
void print_litho_defaults_(Litho *);
double temp_plt_(Litho *, double *, double *);
double depth_sflr_(Litho *, double *);
double pressure_(Litho *, double *, double *, unsigned int *);
double byerlee_(Litho *, double *, double *, unsigned int *);
double ductile_(Litho *, double *, unsigned int *);

int main (int argc, char **argv)
{
  int j,nz=1000;
  double age,z,dz,dsf,temp,obp,zp,zmt,bystr,dustr;
  double ystrp,ystrm;

  /* switches */
  unsigned int wcsw = 1; /* don't include water column overburden = 0, include = 1 */
  unsigned int bysw; /* compression = 0, tension = 1 */
  unsigned int dusw; /* dorn law = 0, power law = 1 */

  Litho l;
  Litho *lptr = &l; /* pointer to litho structure */

  /* get info from command line */

  if(argc < 2){
    printf("\n Usage ocean_litho_yse age \n \n");
    exit(-1);
  }

  age = atof(argv[1]);
  set_litho_defaults_(lptr);
  /*print_litho_defaults_(lptr); */

  dsf=depth_sflr_(lptr,&age);

  zp = lptr->dp;
  dz = zp/nz;
  dusw = 0;
  zmt = zp; /* Mechanical thickness is plate thickness initially*/

  for(j=0;j<nz;j++) {
    z=((double)j+0.5)*dz;
    temp=temp_plt_(lptr,&z,&age);
    obp=pressure_(lptr,&z,&dsf,&wcsw);
    dustr=ductile_(lptr,&temp,&dusw);

    if (dustr <= 2.e8 && dusw != 1) {
      dusw = 1;
      dustr=ductile_(lptr,&temp,&dusw);
    }

    bysw = 1; /* in tension */
    bystr = byerlee_(lptr,&z,&obp,&bysw);
    ystrp = fmin(bystr,dustr);
    bysw = 0; /*in compression */
    bystr = byerlee_(lptr,&z,&obp,&bysw);
    ystrm = fmax(bystr,-dustr);

    fprintf(stdout,"%lf %lf %lf %lf %lf \n",z*1.e-3,temp,obp*1.e-6,ystrp*1.e-6,ystrm*1.e-6);

  }
  return(EXIT_SUCCESS);

}
