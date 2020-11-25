/*******************************************************************************
* Compute yield strength envelopes for oceanic lithosphere of variable age *
*******************************************************************************/

#include "litho.h"

void set_litho_defaults_(Litho *);
void print_litho_defaults_(Litho *);
double temp_plt_(Litho *, double *, double *);
double temp_sleep_(Litho *, double *, double *, unsigned int *);
double depth_sflr_(Litho *, double *);
double pressure_(Litho *, double *, double *, unsigned int *);
double byerlee_(Litho *, double *, double *, unsigned int *);
double ductile_(Litho *, double *, unsigned int *);

int main (int argc, char *argv[])
{
  int i,j,nz=1000;
  double age,z,dz,dsf,temp,obp,zp,zmt,bystr,dustr;
  double ystrp,ystrm;

  /* switches */
  unsigned int bysw; /* compression = 0, tension = 1 */
  unsigned int dusw; /* dorn law = 0, power law = 1 */
  unsigned int wcsw = 1; /* don't include water column overburden = 0, include = 1 */
  unsigned int tesw = 1; /* switch between thermal models, 0 = plate, 1 = sleep */
  unsigned int hssw = 0; /* heat sinks, models hydrothermal circulation. 0 = no, 1 = yes */
  unsigned int flow = 1; /* crustal flow, 0 = same as mantle, 1 = wet olivine, 2 = diabase */

  unsigned int reset = 0; /* for switching from crust to mantle */

  Litho l;
  Litho *lptr = &l; /* pointer to litho structure */

  /* get info from command line */
  if(argc < 2){
    fprintf(stderr,"\n Usage ocean_litho_yse age (wcsw 0/1 tesw 0/1 hssw 0/1) \n");
    exit(-1);
  }
  /* parse command line info */
  age = atof(argv[1]);
  if(argc > 2){
    i = 2;
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
      else if(strncmp(argv[i],"flow",4) == 0) {
        flow = atoi(argv[i+1]);
      }
      else {
        fprintf(stderr,"Not a valid flag\n");
        fprintf(stderr,"\n Usage ocean_litho_yse age (wcsw 0/1 tesw 0/1 hssw 0/1 flow 0/1/2) \n");
        exit(-1);
      }
      i+=2;
    }
  }

  set_litho_defaults_(lptr);
  /*print_litho_defaults_(lptr); */

  dsf=depth_sflr_(lptr,&age);
  /* set spreading rate  in m / yr*/
  lptr->usp = 0.02;
  zp = lptr->dp;
  dz = zp/nz;
  dusw = 1;
  zmt = zp; /* Mechanical thickness is plate thickness initially*/

  /* adjust parameters for the crust */
  if (flow == 1){
    // wet olivine, Karato et al. 1986
    lptr->str_exp = 3.0;
    lptr->str_pow = 1.9e-15;
    lptr->qp = 4.2e5;
  }
  else if (flow == 2){
    // diabase, citation needed
    lptr->str_exp = 4.7;
    lptr->str_pow = 5.0e-28;
    lptr->qp = 4.82e5;
  }
  else {
    lptr->str_exp = 3.5;
    lptr->str_pow = 2.4e-16;
    lptr->qp = 5.40e5;
  }

  /*for(j=0;j<nz;j++) { */
  /* instead of going all the way to the base, just go to 25 km */
  j = 0;
  z = ((double)j+0.5)*dz;
  while (z < 25000) {
    /*z=((double)j+0.5)*dz; */

    if (z > lptr->dc && reset == 0) {
      /* power law flow for dry olivine, Karato et al., 1986 */
      set_litho_defaults_(lptr);
      lptr->str_exp = 3.5;
      lptr->str_pow = 2.4e-16;
      lptr->qp = 5.40e5;
      dusw = 1;
      reset = 1;
    }

    if (tesw == 1) {
      temp=temp_sleep_(lptr,&z,&age,&hssw);
    }
    else {
      temp=temp_plt_(lptr,&z,&age);
    }

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
    j++;
    z=((double)j+0.5)*dz;
  }
  return(EXIT_SUCCESS);

}
