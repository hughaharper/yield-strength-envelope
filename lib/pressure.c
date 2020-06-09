#include "litho.h"

/*-----------------------------------------------------------------*/
/* computes pressure due to overlying water and rock column        */
/*-----------------------------------------------------------------*/

double pressure_(Litho *l, double *z, double *dsf, unsigned int *wcsw)
{
        double press,ps;
        /* pressure from water column overburden */
        
        if(*wcsw == 1) {
        	ps = (*dsf+l->dref)*l->rw*GBAR;
        }
        else {
        	ps = 0.0;
        }
        // press = ps + *z*(l->rc-(l->phyd*l->rw))*GBAR;
		press = ps + (*z*l->rc*GBAR)-(l->phyd*l->rw*GBAR)*(*z+*dsf+l->dref);
		
    return press;
}
