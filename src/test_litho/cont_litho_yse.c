/*******************************************************************************
* Compute yield strength envelopes for continental lithosphere *
* Crust goes by westerly granite flow law *
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
  double z,dz,temp,obp,zp,zmt,bystr,dustr;
  double ystrp,ystrm;

  /* ductile flow law for westerly granite */
  double wg_exp = 2.4;
  double wg_pow = 1.2e-15;
  double wg_qp = 2.93e5;

  double age = 100;
  double tc = 4.e4; /* 40 km crustal thickness */
  double dsf = 0;
  /* switches */
  unsigned int wcsw = 0;
  unsigned int bysw; /* compression = 0, tension = 1 */
  unsigned int dusw; /* dorn law = 0, power law = 1 */
  unsigned int reset = 0; /* for switching from upper to lower lith */

  Litho l;
  Litho *lptr = &l; /* pointer to litho structure */

  set_litho_defaults_(lptr);

  zp = lptr->dp;
  dusw = 1;
  dz = zp/nz;
  zmt = zp; /* Mechanical thickness is plate thickness initially*/

  /* change some lith values for upper crust */
  lptr->phyd = 0.0;
  lptr->str_exp = wg_exp;
  lptr->str_pow = wg_pow;
  lptr->qp = wg_qp;

  for(j=0;j<nz;j++) {
    z=((double)j+0.5)*dz;

    if (z > tc && reset == 0) {
      fprintf(stderr, "reset defaults at depth: %lf \n",z);
      set_litho_defaults_(lptr);
      /*dusw = 0;*/
      reset = 1;
    }

    temp=temp_plt_(lptr,&z,&age);
    obp=pressure_(lptr,&z,&dsf,&wcsw);
    dustr=ductile_(lptr,&temp,&dusw);

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
