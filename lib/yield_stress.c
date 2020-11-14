#include "litho.h"

/*------------------------------------------------------------*/
/* computes yield strength according to form of Byerlee's law */
/* (values from Mueller and Philllips, 1995)                  */
/*------------------------------------------------------------*/

double depth_sflr_(Litho *, double *);
double pressure_(Litho *, double *, double *, unsigned int *);
double temp_plt_(Litho *, double *, double *);
double temp_sleep_(Litho *, double *, double *, unsigned int *);
double byerlee_(Litho *, double *, double *, unsigned int *);
double ductile_(Litho *, double *, unsigned int *);

double yield_stress_(Litho *l, double *zpt, double *xpt, unsigned int *wcsw, \
                    unsigned int *tesw, unsigned int *hssw)
{
  double age,dsf,temp,obp,dustr,bystr,ystr;

  unsigned int dusw = 1;
  unsigned int bysw = 0;
  // check for errors...

  // convert x to age in Myr
  age = *xpt / l->usp / 1e6;

  //
  dsf = depth_sflr_(l,&age);

  // choose flow law based on z > crustal thickness
  if (*zpt <= l-> dc) {
    // flow law for wet olivine
    l->str_exp = 3.0;
    l->str_pow = 1.9e-15;
    l->qp = 4.2e5;
  }
  else if (*zpt > l->dc) {
    // power law flow for dry olivine
    l->str_exp = 3.5;
    l->str_pow = 2.4e-16;
    l->qp = 5.40e5;
  }

  if (*tesw == 1) {
    temp = temp_sleep_(l,zpt,&age,hssw);
  }
  else {
    temp = temp_plt_(l,zpt,&age);
  }

  obp = pressure_(l,zpt,&dsf,wcsw);
  dustr = ductile_(l,&temp,&dusw);

  bystr = byerlee_(l,zpt,&obp,&bysw);
  ystr = fmax(bystr,-dustr);

  return ystr;
}
