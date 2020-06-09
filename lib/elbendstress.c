#include "litho.h"

/*-----------------------------------------------------------------*/
/* computes horizontal stress in for thin elastic plate bending    */
/*-----------------------------------------------------------------*/

double elbendstress_(Litho *l, double *zloc, double *curv)
{
        double ebstress;

        ebstress=(*curv*l->young**zloc/(1-(l->pois*l->pois)));
 
    return ebstress;
}
