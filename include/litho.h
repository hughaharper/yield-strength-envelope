/* litho.h
Structure for defining the thermal and mechanical properties of the lithosphere
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <time.h>

#define PI 3.1415926535897932
#define PI2 9.8696044010893586
#define RT 8.3144621                 /* universal gas constant */
#define SPMYR 3.15576e13              /* seconds per million year */
#define GAMA 6.67384e-11             /* gravitational constant */
#define GBAR 9.81					 /* acceleration due to gravity */

/* structure for litho */

struct LITHO {

/*  parameters for plate cooling model */
double  ts;     /*  surface temp */
double  tm;     /*  mantle temperature */
double  diff;   /*  thermal diffusivity */
double  cp;     /*  heat capacity at constant pressure */
double  dp;     /*  plate thickness */
double  usp;    /* full spreading rate [cm yr-1] */
 
/*  parameters for thermal subsidence, sediment loading, and pressure models */
double  gbar;   /*  gravitational acceleration */
double  dref;   /*  ridge crest depth */
double  alph;   /*  volume expansion coeff */
double  rw;     /*  water density */
double  rs;     /*  sediment density */
double  rc;     /*  crustal density */
double  rm;     /*  mantle density */
double  dr;     /*  ridge axis depth */
double  dw;     /*  mean seafloor depth */
double  dc;     /*  crustal thickness */
double  ds;     /*  sediment thickness */

/* parameters for brittle yield strength envelope */
double  byerlpt;  /*  tension, byerlee coeff. for pressure term, for vert. stress > 529.9 MPa */
double  byergpt;  /*  tension, byerlee coeff. for pressure term, for vert. stress < 529.9 MPa */
double  byergst;  /*  tension, byerlee value for cohesion term, for vert. stress < 529.9 MPa */
double  byergpc;  /*  compression, byerlee coeff. for pressure term, for vert. stress < 113.2 MPa */
double  byerlpc;  /*  compression, byerlee coeff. for pressure term, for vert. stress > 113.2 MPa */
double  byerlsc;  /*  compression, byerlee value for cohesion term, for vert. stress > 113.2 MPa */
double  phyd;     /*  pore pressure level */

/* elastic parameters */
double  young;  /*  young's modulus. */
double  pois;   /*  poisson's ratio. */
double  telas;  /*  temperature at base of elastic layer */

/* parameters for ductile flow law */
double  eps1;   /*  strain rate */
double  str_exp;/*  stress exponent */
double  str_pow;/*  stress amplitude factor for power law */
double  str_dor;/* stress constant for Dorn law */
double  sren_dor;/* strain rate factor for Dorn law */
double  qp;     /*  activation energy for power law */
double  qd;     /*  activation energy for Dorn */

/* depth-related parameters */
double  zmt;    /*  mechanical thickness based on ductile strength */
double  zn;     /*  depth of nodal plane */ 
double  zy;     /*  depth to top of elastic layer. */
};

typedef struct LITHO Litho;
