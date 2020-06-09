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
double mechthk_(Litho *l, double *age);
void die (char *s1, char *s2);

int
main(int argc, char **argv)
{
     int     i, na;
     double  age, mt, da, maxage; 
     
     unsigned int verbose = 0;
	 
     Litho lprops;		
     Litho *lptr = &lprops; /* pointer to litho structure */


     /* populate litho structure */
     set_litho_defaults_(lptr);
     if(verbose == 1) print_litho_defaults_(lptr);
     
     da = 0.05;
     i = 0;
     maxage = 180;
     na=(int)(floor(maxage/da));
     
	 for(i=1;i<na;i++){
	 	age = (double)(i*da);
	 	mt = mechthk_(lptr,&age);
	 	fprintf(stdout,"%lf %lf \n",age,mt);    
	 } 

     return(EXIT_SUCCESS);
}
