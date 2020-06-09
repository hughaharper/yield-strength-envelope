/*******************************************************************************
 * Compute temperature and yield strength envelopes for the oceanic lithosphere*
 *******************************************************************************/

#include "litho.h"

void set_litho_defaults_(Litho *l);
void print_litho_defaults_(Litho *l);
double temp_plt_(Litho *l, double *z, double *age);
double depth_sflr_(Litho *l, double *age);
double ductile_(Litho *l, double *temp, unsigned int *dusw);
double yse_moment_(Litho *l, double *age, double *curv, double *ival1, \
                   double *ival2, double *npstr, unsigned int *wcsw);


int main (int argc, char **argv)
{
		unsigned int verbose = 0;
        int i,j, nz = 100, nc = 1000;
        double age,bendmo,curv,dcurv,z,dz,dsf,temp,zp,zmt,dustr;
        
        double ival1, ival2, curv1, curv2, npstr;
        
        /* switches */
        unsigned int wcsw; /* do not include water column overburden = 0, 
                              include it = 1 */
        unsigned int dusw; /* dorn law = 0, power law = 1 */

        Litho l;		
        Litho *lptr = &l; /* pointer to litho structure */

	    /* get the information from the command line */
	    
        if(argc < 2){
                printf("\n Usage: mocurv age \n \n");
                exit(-1);
        }

        age = atof(argv[1]);
        
        /* insert some error handling for zero values of either curvature or age */ 

        set_litho_defaults_(lptr);
        if(verbose == 1) print_litho_defaults_(lptr);
        
        /* set value for in-plane stress (may be taken from user input later) */
        npstr=0.0;

		/* find depth of seafloor from thermal subsidence */
        dsf=depth_sflr_(lptr,&age);
        
        /* compute profile of ductile yielding to estimate "mechanical thickness" */
        zp = lptr->dp;
        dz = zp/nz;
        zmt = zp;
        
        dusw = 0;
        /* switch for (not) including water column overburden pressure */
        wcsw = 0;
        
        for(j=0;j<nz;j++){
        		 z=((double)j+0.5)*dz;
        		 
        		 temp=temp_plt_(lptr,&z,&age);
        		 dustr=ductile_(lptr,&temp,&dusw);
        		 
        		 if (dustr <= 2.e8 && dusw != 1) {
        				dusw = 1;
        				dustr=ductile_(lptr,&temp,&dusw);
        		 } 
        		 
        		 if ( dustr <= 5.e7) {
        				zmt=z;
        				break;
        		 } 
        
        }
        /* might be worth it to recode such that stresses are stored in arrays */
        
        curv1 = -1.e-5;
        curv2 = 1.e-5;
        dcurv = (curv2-curv1)/nc;
        
        curv = curv1;
        
        for(i=0;i<nc;i++) {
            curv = curv1 + ((i-1)*dcurv);
            
            /* reset initial boundaries of interval */
            ival1=1.5*dz;
            ival2=zmt-(1.5*dz);    
            bendmo = yse_moment_(lptr,&age,&curv,&ival1,&ival2,&npstr,&wcsw);
            fprintf(stdout,"%lf %lf \n",curv*1.e6,bendmo*1.e-17);
        }    
        
	return(EXIT_SUCCESS);
}
