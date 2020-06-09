#include "litho.h"

/*-----------------------------------------------------------*/
/* computes pressure due to overlying water column           */
/*-----------------------------------------------------------*/

double ductile_(Litho *l, double *temp, unsigned int *dusw)
{
     double ductstr, tempk;
     
     tempk = *temp+273.15;
        
     switch(*dusw) {
          case 0 : /* Dorn law, diff. stress > 200 Mpa */
               ductstr=l->str_dor*(1-sqrt((tempk*RT/l->qd)*log(l->sren_dor/l->eps1)));
               break;
          case 1 :
               ductstr=pow((l->eps1/l->str_pow)*exp(l->qp/(tempk*RT)),1.0/l->str_exp);
               break;
     }  
       
     return ductstr;
}
