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

  double tc = 7.e3; /* 7 km crustal thickness */
  /* flow law for wet olivine, Karato et al. 1986 */
  /* change this to flow law for Gabbro */
  /*double wo_exp = 3.0;
  double wo_pow = 1.9e-15;
  double wo_qp = 4.2e5; */

  /* diabase */
  double wo_exp = 4.7;
  double wo_pow = 5.0e-28;
  double wo_qp = 4.82e5;

  /* switches */
  unsigned int wcsw = 1; /* don't include water column overburden = 0, include = 1 */
  unsigned int bysw; /* compression = 0, tension = 1 */
  unsigned int dusw; /* dorn law = 0, power law = 1 */
  unsigned int reset = 0; /* for switching from crust to mantle */

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
  dusw = 1;
  zmt = zp; /* Mechanical thickness is plate thickness initially*/

  /* adjust parameters for the crust */
  lptr->str_exp = wo_exp;
  lptr->str_pow = wo_pow;
  lptr->qp = wo_qp;

  for(j=0;j<nz;j++) {
    z=((double)j+0.5)*dz;

    if (z > tc && reset == 0) {
      fprintf(stderr, "reset defaults at depth: %lf \n",z);
      /* power law flow for dry olivine, Karato et al., 1986 */
      set_litho_defaults_(lptr);
      lptr->str_exp = 3.5;
      lptr->str_pow = 2.4e-16;
      lptr->qp = 5.40e5;
      dusw = 1;
      reset = 1;
    }

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
