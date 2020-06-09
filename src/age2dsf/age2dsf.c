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
double mechthk_(Litho *l, double *age);
void die (char *s1, char *s2);

char *USAGE = "Usage: age2dsf age_in.grd dsf_out.grd \n \n"
"    age_in.grd     - name of input seafloor age file. \n"
"    dsf_out.grd    - name of output seafloor depth file \n";


int
main(int argc, char **argv)
{
     int     j, k, m;
     int     xdim=2040, ydim=906; 
     float   *age, *dsf, *mt;
     double  ymin, xmin, yinc, xinc, rland, rdum;
     double  agept, dsfpt; 
     
     unsigned int verbose = 0;
     
     char    agefilename[128], dsffilename[128], title[128];
	 
     Litho lprops;		
     Litho *lptr = &lprops; /* pointer to litho structure */

	 if (argc < 3) die("\n", USAGE);
	 
	 /* prepare the output filename */
	 strcpy(agefilename, argv[1]); 
	 strcpy(dsffilename, argv[2]); 
	 
	 /* allocate the memory for the arrays */
	 age = (float *) malloc(ydim * xdim * sizeof(float));
	 dsf = (float *) malloc(ydim * xdim * sizeof(float));

	 readgrd_(age, &xdim, &ydim, &ymin, &xmin, &yinc, &xinc, &rdum, title, agefilename);

     /* populate litho structure */
     set_litho_defaults_(lptr);
     if(verbose == 1) print_litho_defaults_(lptr);
     
     fprintf(stderr,"xdim: %d ydim: %d \n",xdim,ydim);    
     
	 for (j=0; j<ydim; j++){
	     for (k=0; k<xdim; k++) { 
              m = j*xdim + k;
		      agept = (double)age[m];	
              
              dsfpt = depth_sflr_(lptr,&agept);
              
              if (dsfpt <= lptr->dref) {
              	dsf[m] = (float)lptr->dref;
              }
              else {
                dsf[m] = (float)dsfpt;
              }
              
		}  
	 fprintf(stderr,"row: %d \n",j+1);    
	 }


	 /* write the grd file */ 
	 rdum = -1.e22; 
	 rland  = -1.e22;
	 writegrd_(dsf, &xdim, &ydim, &ymin, &xmin, &yinc, &xinc, &rland, &rdum, title, dsffilename);

     return(EXIT_SUCCESS);
}

void die (char *s1, char *s2)
{
        fprintf(stderr," %s %s \n",s1,s2);
        exit(1);
}
