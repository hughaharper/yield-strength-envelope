#include "litho.h"

/*--------------------------------*/
/* compute lithospheric temperature (degree C) */
/*--------------------------------*/

double temp_sleep_(Litho *l, double *z, double *age, unsigned int *hssw)
{
  int i,k,nsum=501;
  double j,R_p,a_m,A_m,B_m,B_m1,B_m2,B_m3,B_m4;
  double T_temp, T_homo, T_full;
  double rhoc,kappa,gamma;
  double T_c,T_seg,T_m,zp,x,u;
  double T_part,b_m,g_a,g_b,g_c;
  FILE *heat_sinks;
  char cr;
  int nQ;
  double x_Q,z_Q,Q_d; /* how to make this arbitrary? */

  double z_seg=33.e3,z_crust=6.e3,latent_h=1.028e9,conduct=2.5104;
  double l_adiab=1.e-3,d_adiab=0.3e-3,melt_grad=3.e-3; /* alpha, beta */

  u = l->usp/365/24/60/60; /* m/yr to m/s */
  x = l->usp**age*1e6; /* convert m/yr to m/yr, Myr to yr */
  kappa = l->diff;
  rhoc = conduct/kappa;
  zp = l->dp;
  T_m = l->tm;
  gamma = 1 - (d_adiab*zp)/(T_m);

  T_c = T_m*(gamma - ((gamma*z_crust)/zp) + (z_crust/zp));
  T_seg = T_m*(gamma - ((gamma*z_seg)/zp) + (z_seg/zp));


  R_p = (2*kappa*PI)/(u*zp);
  T_homo = 0;
  fprintf(stderr,"Position: x - %.2f, z - %.2f\n",x,*z);
  /* ------------------------------------------------------------------------ */
  if (*hssw == 1) {
    /* read in heat sink data */
    heat_sinks = fopen("heat_sinks.xz","r");
    if (heat_sinks == NULL){
      fprintf(stderr,"Error Reading Heat Sinks File\n");
      exit(-1);
    }

    // Count Lines
    cr = getc(heat_sinks);
    while( cr != EOF ) {
      if ( cr == '\n' ) {
        nQ++;
      }
      cr = getc(heat_sinks);
    }
    rewind(heat_sinks);

  }

  /* ------------------------------------------------------------------------ */

  for(i=0;i<nsum;i++){
    j=(double)(i+1);

    a_m =(u/(2*kappa))*(1 - sqrt(1 + (R_p*R_p*j*j)));

    B_m1 = cos((j*PI*z_seg)/zp)*((1 - (z_seg/zp))*T_m*gamma - T_seg + T_m*(z_seg/zp));
    B_m2 = (1/(j*PI))*sin((j*PI*z_seg)/zp)*(T_m*gamma + melt_grad*zp - T_m);
    B_m3 = cos((j*PI*z_crust)/zp)*(T_seg - (z_seg - z_crust)*melt_grad - (latent_h/rhoc) - T_c);
    B_m4 = (1/(j*PI)) * sin((j*PI*z_crust)/zp) * (l_adiab*zp - melt_grad*zp) + (latent_h/rhoc) + T_c - l_adiab*z_crust;

    B_m = ((2*u*rhoc)/(j*PI))*(B_m1 + B_m2 + B_m3 + B_m4);

    A_m = 2/(1 + sqrt(1 + (R_p*R_p*j*j)));

    T_temp = A_m*B_m*sin((*z*j*PI)/zp)*exp(a_m*x);

    T_homo = T_homo + T_temp;
  }

  T_part = 0.;
  /* Add heat sinks if flag is there */
  if (*hssw == 1) {
    /* for k=0 to k< no. of sinks */
    for(k=0;k<nQ;k++){
      fscanf(heat_sinks,"%lf %lf %lf",&x_Q,&z_Q,&Q_d);
      /* arrays of heat sink positions */
      T_temp = 0.;
      for(i=1;i<201;i++){
        j = (double)(i);
        a_m = (u/(2*kappa))*(1 - sqrt(1 + (R_p*R_p*j*j)));
        b_m = (u/(2*kappa))*(1 + sqrt(1 + (R_p*R_p*j*j)));
        g_a = (2*Q_d)/(kappa*zp*(b_m-a_m));
        g_b = sin((j*PI*z_Q)/zp);
        if (x < x_Q) {
          g_c = (exp(b_m*x) - exp(a_m*x))*exp(-1*b_m*x_Q);
        }
        else if (x > x_Q) {
          g_c = (exp(-1*a_m*x_Q) - exp(-1*b_m*x_Q))*exp(a_m*x);
        }
        T_temp = T_temp + g_a*g_b*g_c*sin((j*PI**z)/zp);
      } /* end sum loop */
      T_part = T_part + T_temp;
    } /* end loop thru heat sinks */
    fclose(heat_sinks);
  }

  T_full = (1/(u*rhoc))*T_homo + ((*z*T_m)/zp) + T_part;
  return(T_full);
}
