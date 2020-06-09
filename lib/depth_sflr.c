#include "litho.h"

/*------------------------------------------------------------------*/
/* computes seafloor topography due to thermal subsidence           */
/*------------------------------------------------------------------*/

double depth_sflr_(Litho *l, double *age)
{
        double argex,term,tsum=0.;
        double depth, tage, pref;
        int i,nsum=50; /* number of terms in series */
        double j;
        
        tage = *age*SPMYR;
        
        pref = l->rm*l->alph*l->dp*(l->tm-l->ts)/(l->rm-l->rw);
        
        for(i=0;i<nsum;i++){
            j=(double)(((2*i)+1)*((2*i)+1));
            argex=(l->diff*j*PI2*tage/(l->dp*l->dp));
            term=exp(-1.0*argex);
            tsum=tsum+term;
        }
        
        depth = l->dref + pref*(0.5-4.0*tsum/PI2);

    return depth;
}
