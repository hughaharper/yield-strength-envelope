/***************************************************************************/
/* curv2rigid reads residual phase and computes residual topography.       */
/***************************************************************************/

/***************************************************************************
 * Modification history:                                                   *
 *                                                                         *
 ***************************************************************************/

#include"litho.h"

void set_litho_defaults_(Litho *l);
void print_litho_defaults_(Litho *l);
double temp_plt_(Litho *l, double *z, double *age);
double depth_sflr_(Litho *l, double *age);
double ductile_(Litho *l, double *temp, unsigned int *dusw);
double mechthk_(Litho *l, double *age);
double yse_moment_(Litho *l, int *zpts, double *age, double *curv, double *ival1, \
                   double *ival2, double *zmt, double *npstr, unsigned int *wcsw, \
                   unsigned int *znsw, int *rfiter, unsigned int *vsw);
void die (char *s1, char *s2);

char *USAGE = "Usage: curv2rigid age_in.grd curv_in.grd mt_in.grd \
					  rigid_out.grd bendmo_out.grd yld_out.grd \n \n"
"    age_in.grd     - name of input seafloor age file. \n"
"    curv_in.grd    - name of input curvature file \n"
"    mt_in.grd    	- name of input mechanical thickness file \n"
"    rigid_out.grd  - name of output flexural rigidity file. \n"
"    bendmo_out.grd - name of output bending moment file. \n"
"    yld_out.grd    - name of output depth of yielding file. \n\n";

int
main(int argc, char **argv)
{
     int     i;
     int	 nc, nz;	  
     double  rland, rdum;
     double  ival1, ival2, zmt, dz, dsf, zmax;
     double  age, curv, mt, dc, cbound, rigid; 
     double  npstr, bendmo;
	 
	 /* switches */
	 unsigned int verbose = 0;
	 unsigned int wcsw; /* do not include water column overburden = 0, 
                              include it = 1 */
     unsigned int znsw = 0;
     unsigned int vsw = 0;     
     
     int riter, totriter;
	 
     Litho lprops;		
     Litho *lptr = &lprops; /* pointer to litho structure */
     
     
     zmax = 8.e4;
     dz = 5.e2;
     dc = 1.e-7;
     cbound = 1e-5;
     
     nc=(int)(floor(2.0*cbound/dc));
     nz=(int)(floor(zmax/dz));
    

     if(argc < 2){
     	printf("\n Usage: mocurv age \n \n");
     	exit(-1);
     }

     age = atof(argv[1]);

     /* populate litho structure */
     set_litho_defaults_(lptr);
     if(verbose == 1) print_litho_defaults_(lptr);
        
     /* set value for in-plane stress (may be taken from user input later) */
     npstr=0.0;
        
     riter = 0;
     totriter = 0;
     
     /* test pore pressure */
     lptr->phyd = 1.0;
     /* set switch for water column overburden */
     wcsw = 1.0;
     /* move the above to set_litho_defaults */
     
     
     curv = -cbound-dc;
     
	 for (i=0; i<nc; i++){
	
        curv += dc;
        
        mt = mechthk_(lptr,&age);
        dsf=depth_sflr_(lptr,&age);              
              

		/* find depth of seafloor from thermal subsidence */
		      
        if (curv < 1.e-9 && curv >= 0.0) curv = 1.e-9;
    	if (curv > -1.e-9 && curv < 0.0) curv = -1.e-9; 
		   
              
		if(age > 5.0) {
        	ival1=dz;
            ival2=zmt+5.0e3;              
        }
        else {
         	ival1=dz;
            ival2=zmt+2.0e4;               
        }    
              
        bendmo = yse_moment_(lptr,&nz,&age,&curv,&ival1,&ival2,&mt,&npstr,&wcsw,&znsw,&riter,&vsw); 
              
        if (znsw == 0) {
        	ival2 = zmax;
            bendmo = yse_moment_(lptr,&nz,&age,&curv,&ival1,&ival2,&mt,&npstr,&wcsw,&znsw,&riter,&vsw);
              
        } 
              
        if(znsw == 0) {
        	fprintf(stderr,"Warning: root not found, age: %lf, mechanical thickness: %lf \n",age,mt);
            vsw=1;
            bendmo = yse_moment_(lptr,&nz,&age,&curv,&ival1,&ival2,&mt,&npstr,&wcsw,&znsw,&riter,&vsw);
            vsw=0;
              	
            return(EXIT_FAILURE);
        }
              
        totriter += riter;
        
        rigid = ((2.0*bendmo)/(curv*(1.0+lptr->pois)));
               
        if(rigid < 0.0) {
        	fprintf(stderr,"Warning: negative rigidity \n");
            vsw=1;
            bendmo = yse_moment_(lptr,&nz,&age,&curv,&ival1,&ival2,&zmt,&npstr,&wcsw,&znsw,&riter,&vsw);
            vsw=0;
                  
            return(EXIT_FAILURE);
        }
        	 
        	 
		fprintf(stdout,"%g %g \n",curv,bendmo);      
	 }


	 /* write the grd file */ 
	 rdum = -1.e22; 
	 rland  = -1.e22;

     return(EXIT_SUCCESS);
}

void die (char *s1, char *s2)
{
        fprintf(stderr," %s %s \n",s1,s2);
        exit(1);
}
