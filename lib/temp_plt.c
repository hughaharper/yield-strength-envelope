#include "litho.h"

/*------------------------------------------------------------------------------*/
/* computes lithospheric temperature (degree C) and ocean floor topography (m)  */
/*------------------------------------------------------------------------------*/

double temp_plt_(Litho *l, double *z, double *age)
{
        double argex,argsin,term,tsum=0.;
        double temp, tage;
        int i,nsum=50; /* number of terms in series */
        double j;
        
        tage = *age*SPMYR;
        
        for(i=0;i<nsum;i++){
            j=(double)(i+1);
            argex=l->diff*j*j*PI2*tage/(l->dp*l->dp);
            argsin=*z*j*PI/l->dp;
            term=exp(-1.0*argex)*sin(argsin)/j;
            tsum=tsum+term;
        }
        
        temp = l->ts + (l->tm-l->ts)*((*z/l->dp)+(2.0*tsum/PI));

    return temp;
}
