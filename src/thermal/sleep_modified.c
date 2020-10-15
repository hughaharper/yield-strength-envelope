/*
 * compute temperature of oceanic lith using Morton and Sleep (1985) model
*/

#include "litho.h"

void set_litho_defaults_(Litho *);

int main (int argc, char **argv)
{
  /* variables */
  int i,j,nz,nx,m,k;
  double x,dx,z,dz,zp,u;
  double T_all,T_homo,T_part,T_temp;
  double M,R_p,a_m,A_m,B_m,B_m1,B_m2,B_m3,B_m4;
  double b_m,g_a,g_b,g_c;
  int nQ=17;
  double x_Q[nQ],z_Q[nQ],Q_d[nQ]; /* how to make this arbitrary? */
  double rhoc,kappa,gamma;
  double T_c,T_seg,T_m;
  FILE *heat_sinks;

  Litho l;
  Litho *lptr = &l;

  /* numerical parameters (nominal) */
  /* t_lith, Temp_mantle, spec heat */
  double z_seg=33e3,z_crust=5e3,latent_h=1.028e9,conduct=2.5104;
  double l_adiab=1e-3,d_adiab=0.3e-3,melt_grad=3e-3; /* alpha, beta */

  if(argc < 2){
    printf("\n Usage sleep_cooling half_spreading_rate flags\n");
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

  dz = 200; /* z spacing */
  dx = 200; /* x spacing */
  nz = (int) 30000/dz; /* model space only 30km by 30km, no reason to go to 125 km depth */
  nx = (int) 30000/dx;
  fprintf(stderr,"nz: %d, nx: %d\n",nz,nx);

  /* read in heat sink data */
  heat_sinks = fopen("heat_sinks.xz","r");
  if (heat_sinks == NULL){
    fprintf(stderr,"Error Reading File\n");
    exit(-1);
  }
  for(j=0;j<nQ;j++){
    fscanf(heat_sinks,"%lf %lf %lf",&x_Q[j],&z_Q[j],&Q_d[j]);
    fprintf(stderr,"x_pos: %lf z_pos: %lf Q:%lf\n",x_Q[j],z_Q[j],Q_d[j]);
  }

  fprintf(stdout,"0 ");
  for(j=0;j<nz;j++) {
    z=((double)j)*dz;
    fprintf(stdout,"%lf ",z);
  }

  for(j=0;j<nx;j++) {
    /* loop through x distance */
    x = ((double) j)*dx;

    fprintf(stdout,"\n%lf ",x);

    for(i=0;i<nz;i++) {
      /* loop thru depths */
      z = ((double) i)*dz;

      /* evaluate initial terms */
      T_homo = 0;
      R_p = (2*kappa*PI)/(u*zp);

      /* Compute the homogeneous solutions */
      for(m=1;m<501;m++){
        M = (double) m;
        a_m = (u/(2*kappa))*(1 - sqrt(1 + (R_p*R_p*M*M)));

        B_m1 = cos((M*PI*z_seg)/zp)*((1 - (z_seg/zp))*T_m*gamma - T_seg + T_m*(z_seg/zp));
        B_m2 = (1/(M*PI))*sin((M*PI*z_seg)/zp)*(T_m*gamma + melt_grad*zp - T_m);
        B_m3 = cos((M*PI*z_crust)/zp)*(T_seg - (z_seg - z_crust)*melt_grad - (latent_h/rhoc) - T_c);
        B_m4 = (1/(M*PI)) * sin((M*PI*z_crust)/zp) * (l_adiab*zp - melt_grad*zp) + (latent_h/rhoc) + T_c - l_adiab*z_crust;

        B_m = ((2*u*rhoc)/(M*PI))*(B_m1 + B_m2 + B_m3 + B_m4);

        A_m = 2/(1 + sqrt(1 + (R_p*R_p*M*M)));

        T_temp = A_m*B_m*sin((M*PI*z)/zp)*exp(a_m*x);

        T_homo = T_homo + T_temp;

      } /* end sum loop */

      /* Compute the particular solution for heat sinks */
      T_part = 0;
      /* for k=0 to k< no. of sinks */
      for(k=0;k<nQ;k++){
        /* arrays of heat sink positions */
        T_temp = 0;
        for(m=1;m<501;m++){
          M = (double) m;
          a_m = (u/(2*kappa))*(1 - sqrt(1 + (R_p*R_p*M*M)));
          b_m = (u/(2*kappa))*(1 + sqrt(1 + (R_p*R_p*M*M)));
          g_a = (2*Q_d[k])/(kappa*zp*(b_m-a_m));
          g_b = sin((M*PI*z_Q[k])/zp);
          if (x < x_Q[k]) {
            g_c = (exp(b_m*x) - exp(a_m*x))*exp(-1*b_m*x_Q[k]);
          } else if (x > x_Q[k]) {
            g_c = (exp(-1*a_m*x_Q[k]) - exp(-1*b_m*x_Q[k]))*exp(a_m*x);
          }

          T_temp = T_temp + g_a*g_b*g_c*sin((M*PI*z)/zp);
        } /* end sum loop */
        T_part = T_part + T_temp;
      } /* end loop thru heat sinks */

      /* Full solution */
      T_all = (1/(u*rhoc))*T_homo + T_part + ((T_m*z)/zp);

      fprintf(stdout,"%lf ",T_all);
    } /* end depth loop */
  } /* end age (distance) loop */
  return(EXIT_SUCCESS);

} /* main */
