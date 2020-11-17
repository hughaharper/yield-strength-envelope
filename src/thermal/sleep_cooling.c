/*
 * compute temperature of oceanic lith using Sleep (1975) model
*/

#include "litho.h"

void set_litho_defaults_(Litho *);

int main (int argc, char **argv)
{
  /* variables */
  int i,j,nz,nx,m;
  double x,dx,u,z,zp,dz,temp,t_ta,t_tb;
  double M,R_p,a_m,A_m,B_m,B_m1,B_m2,B_m3,B_m4;
  double rhoc,kappa,gamma;
  double T_c,T_seg,T_m;

  Litho l;
  Litho *lptr = &l;

  /* numerical parameters (nominal) */
  /* t_lith, Temp_mantle, spec heat */
  double z_seg=33e3,z_crust=5e3,latent_h=1.028e9,conduct=2.5104;
  double l_adiab=1e-3,d_adiab=0.3e-3,melt_grad=3e-3; /* alpha, beta */

  if(argc < 2){
    fprintf(stderr,"\n Usage sleep_cooling half_spreading_rate flags\n");
    exit(-1);
  }

  u = atof(argv[1])/365/24/60/60; /* rate in m/yr, convert to m/s */

  set_litho_defaults_(lptr);

  /*kappa = lptr->diff; */
  rhoc = 3.807e6;
  kappa = conduct/rhoc;
  zp = lptr->dp;
  T_m = lptr->tm;
  zp = 1e5;
  T_m = 1.29e3;
  gamma = 1 - (d_adiab*zp)/(T_m);

  T_c = T_m*(gamma - ((gamma*z_crust)/zp) + (z_crust/zp));
  T_seg = T_m*(gamma - ((gamma*z_seg)/zp) + (z_seg/zp));

  dz = 100;
  dx = 100;
  nz = (int) 30000/dz;
  nx = (int) 30000/dx;
  fprintf(stderr,"nz: %d, nx: %d\n",nz,nx);

  fprintf(stdout,"0 ");
  for(j=0;j<nz;j++) {
    z=((double)j)*dz;
    fprintf(stdout,"%lf ",z);
  }

  for(j=0;j<nx;j++) {
    /* loop through x distance (age) */
    x = ((double)j)*dx;

    fprintf(stdout,"\n%lf ",x);

    for(i=0;i<nz;i++) {
      /* loop thru depths */
      z = ((double)i)*dz;

      /* evaluate initial terms */
      t_tb = 0;
      R_p = (2*kappa*PI)/(u*zp);

      for(m=1;m<501;m++){
        M = (double) m;
        a_m =(u/(2*kappa))*(1 - sqrt(1 + (R_p*R_p*M*M)));

        B_m1 = cos((M*PI*z_seg)/zp)*((1 - (z_seg/zp))*T_m*gamma - T_seg + T_m*(z_seg/zp));
        B_m2 = (1/(M*PI))*sin((M*PI*z_seg)/zp)*(T_m*gamma + melt_grad*zp - T_m);
        B_m3 = cos((M*PI*z_crust)/zp)*(T_seg - (z_seg - z_crust)*melt_grad - (latent_h/rhoc) - T_c);
        B_m4 = (1/(M*PI)) * sin((M*PI*z_crust)/zp) * (l_adiab*zp - melt_grad*zp) + (latent_h/rhoc) + T_c - l_adiab*z_crust;

        B_m = ((2*u*rhoc)/(M*PI))*(B_m1 + B_m2 + B_m3 + B_m4);

        A_m = 2/(1 + sqrt(1 + (R_p*R_p*M*M)));

        t_ta = A_m*B_m*sin((M*PI*z)/zp)*exp(a_m*x);

        /* constant intrusion temp solution */
        /*t_ta = A_m*sin((M*PI*z)/zp)*exp(a_m*x)/M; */

        t_tb = t_tb + t_ta;

      }

      temp = (1/(u*rhoc))*t_tb + ((T_m*z)/zp);

      /* constant intrusion temp solution */
      /*temp = ((T_m*z)/zp) + ((2*T_m)/PI)*t_tb; */

      fprintf(stdout,"%lf ",temp);
    } /* end depth loop */
  } /* end age (distance) loop */
  return(EXIT_SUCCESS);

} /* main */
