#include "litho.h"

/*--------------------------------*/
/* compute lithospheric temperature (degree C) */
/*--------------------------------*/

double temp_sleep_(Litho *l, double *z, double *age)
{
  int i,nsum=50;
  double j,R_p,a_m,A_m,B_m,B_m1,B_m2,B_m3,B_m4;
  double t_ta, t_tb, temp;
  double rhoc,kappa,gamma;
  double T_c,T_seg,T_m,zp,x,u;

  double z_seg=33.e3,z_crust=7.e3,latent_h=1.028e9,conduct=2.5104;
  double l_adiab=1.e-3,d_adiab=0.3e-3,melt_grad=3.e-3; /* alpha, beta */

  u = 1e-2*l->usp/365/24/60/60; /* cm/yr to m/s */
  x = l->usp**age*1e4; /* convert cm/yr to m/yr, Myr to yr */
  kappa = l->diff;
  rhoc = conduct*kappa;
  zp = l->dp;
  T_m = l->tm;
  gamma = 1 - (d_adiab*zp)/(T_m);

  T_c = T_m*(gamma - ((gamma*z_crust)/zp) + (z_crust/zp));
  T_seg = T_m*(gamma - ((gamma*z_seg)/zp) + (z_seg/zp));


  R_p = (2*kappa*PI)/(u*zp);
  t_tb = 0;

  for(i=0;i<nsum;i++){
    j=(double)(i+1);

    a_m =(u/(2*kappa))*(1 - sqrt(1 + (R_p*R_p*j*j)));

    B_m1 = cos((j*PI*z_seg)/zp)*((1 - (z_seg/zp))*T_m*gamma - T_seg + T_m*(z_seg/zp));
    B_m2 = (1/(j*PI))*sin((j*PI*z_seg)/zp)*(T_m*gamma + melt_grad*zp - T_m);
    B_m3 = cos((j*PI*z_crust)/zp)*(T_seg - (z_seg - z_crust)*melt_grad - (latent_h/rhoc) - T_c);
    B_m4 = (1/(j*PI)) * sin((j*PI*z_crust)/zp) * (l_adiab*zp - melt_grad*zp) + (latent_h/rhoc) + T_c - l_adiab*z_crust;

    B_m = ((2*u*rhoc)/(j*PI))*(B_m1 + B_m2 + B_m3 + B_m4);

    A_m = 2/(1 + sqrt(1 + (R_p*R_p*j*j)));

    t_ta = A_m*B_m*sin((*z*j*PI)/zp)*exp(a_m*x);

    t_tb = t_tb + t_ta;
  }

  temp = (1/(u*rhoc))*t_tb + ((*z*T_m)/zp);

  return temp;
}
