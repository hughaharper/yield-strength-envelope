#include "litho.h"

/*------------------------------------------------------------*/
/* computes yield strength according to form of Byerlee's law */
/* (values from Mueller and Philllips, 1995)                  */
/*------------------------------------------------------------*/

double byerlee_(Litho *l, double *z, double *obp, unsigned int *bysw)
{
        double byerstr;

        if(*bysw == 1){
			  /* tensional regime */
            if (*obp <= 5.299e8){
                 byerstr = *obp*l->byerlpt;
        	  }
        	  else {
        	     byerstr = l->byergst + *obp*l->byergpt;
        	  }
        }
        else {
            /* compressional regime  */
            if (*obp <= 1.132e8){
                 byerstr = -1.0*(*obp*l->byergpc);
        	  }
        	  else {
        	     byerstr = -1.0*(l->byerlsc + *obp*l->byerlpc);
        	  }
        }

    return byerstr;
}
