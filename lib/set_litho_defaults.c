#include "litho.h"

/*------------------------------------------------------*/
/* set the default parameters for the lithosphere       */
/*------------------------------------------------------*/
void set_litho_defaults_(Litho *litho)
{

litho->ts      = 0.;      /*  surface temp */
litho->tm      = 1365.;   /*  mantle temperature, from Mueller and Phillips */
litho->diff    = 8.0e-7;   /*  thermal diffusivity */
litho->cp      = 1172.;   /*  heat capacity at constant pressure */
litho->dp      = 1.25e5;  /*  plate thickness */

/*parameters for thermal subsidence, sediment loading, and pressure models */
litho->gbar    = 9.81;    /*  gravitational acceleration */
litho->dref    = 2600.;   /*  ridge crest depth */
litho->alph    = 3.1e-5;  /*  volume expansion coeff */
litho->rw      = 1025;    /*  water density */
litho->rs      = 2300;    /*  sediment density */
litho->rc      = 2800;    /*  crustal density */
litho->rm      = 3200;    /*  mantle density */
litho->dw      = 4100;    /*  mean seafloor depth */
litho->dc      = 6000;    /*  crustal thickness */
litho->ds      = 500.;    /*  sediment thickness */

/* parameters for brittle yield strength envelope
   assumes optimal fault orientations */
litho->byerlpt = 0.786;   /*  tension, byerlee coeff. for pressure term, for vert. stress > 529.9 MPa */ 
litho->byergpt = 0.679;   /*  tension, byerlee coeff. for pressure term, for vert. stress < 529.9 MPa */
litho->byergst = 5.67e7;  /*  tension, byerlee coeff. for cohesion term, for vert. stress < 529.9 MPa */
litho->byergpc = 3.68;    /*  compression, byerlee coeff. for pressure term, for vert. stress > 113.2 MPa */
litho->byerlpc = 2.12;    /*  compression, byerlee coeff. for pressure term, for vert. stress < 113.2 MPa */    
litho->byerlsc = 1.766e8; /*  compression, byerlee coeff. for cohesion term, for vert. stress < 113.2 MPa */
litho->phyd    = 1.0;     /*  pore pressure level */

/* elastic parameters */
litho->young   = 6.5e10;   /*  young's modulus. */
litho->pois    = 0.25;    /*  poisson's ratio. */
litho->telas   = 600;     /*  temperature at base of elastic layer */

/* parameters for ductile flow law */
litho->eps1    = 1.e-16;  /*  strain rate */
litho->str_exp = 3.0;     /*  stress exponent */
litho->str_pow = 7.e-14;  /*  stress amplitude factor for power law */
litho->str_dor = 8.5e9;   /* stress constant for Dorn law */
litho->sren_dor= 5.7e11;  /* strain rate factor for Dorn law */
litho->qp      = 5.23e5;  /*  activation energy for power law */
litho->qd      = 5.49e5;  /*  activation energy for Dorn */

/* parameters to be solved for when finding effective rigidity */
litho->zmt     = -1.0;
litho->zn      = -1.0; 
litho->zy      = 1000;    /*  depth to top of elastic layer. */
}

void print_litho_defaults_(Litho *litho)
{
        fprintf(stderr," \n default settings *************\n\n");
        fprintf(stderr," ts       = %g \n",litho->ts      );
        fprintf(stderr," tm       = %g \n",litho->tm      );
        fprintf(stderr," diff     = %g \n",litho->diff    );
        fprintf(stderr," cp       = %g \n",litho->cp      );
        fprintf(stderr," dp       = %g \n",litho->dp      );
        fprintf(stderr," gbar     = %g \n",litho->gbar    );
        fprintf(stderr," dref     = %g \n",litho->dref    );
        fprintf(stderr," alph     = %g \n",litho->alph    );
        fprintf(stderr," rw       = %g \n",litho->rw      );
        fprintf(stderr," rs       = %g \n",litho->rs      );
        fprintf(stderr," rc       = %g \n",litho->rc      );
        fprintf(stderr," rm       = %g \n",litho->rm      );
        fprintf(stderr," dw       = %g \n",litho->dw      );
        fprintf(stderr," dc       = %g \n",litho->dc      );
        fprintf(stderr," ds       = %g \n",litho->ds      );
        fprintf(stderr," byerlpt  = %g \n",litho->byerlpt );
        fprintf(stderr," byergpt  = %g \n",litho->byergpt );
        fprintf(stderr," byergst  = %g \n",litho->byergst );
        fprintf(stderr," byergpc  = %g \n",litho->byergpc );        
        fprintf(stderr," byerlpc  = %g \n",litho->byerlpc );
        fprintf(stderr," byerlsc  = %g \n",litho->byerlsc );
        fprintf(stderr," phyd     = %g \n",litho->phyd    );
        fprintf(stderr," young    = %g \n",litho->young   );
        fprintf(stderr," pois     = %g \n",litho->pois    );
        fprintf(stderr," telas    = %g \n",litho->telas   );
        fprintf(stderr," eps1     = %g \n",litho->eps1    );
        fprintf(stderr," str_exp  = %g \n",litho->str_exp );
        fprintf(stderr," str_pow  = %g \n",litho->str_pow );
        fprintf(stderr," str_dor  = %g \n",litho->str_dor );
        fprintf(stderr," sren_dor = %g \n",litho->sren_dor);
        fprintf(stderr," qp       = %g \n",litho->qp      );
        fprintf(stderr," qd       = %g \n",litho->qd      );
        fprintf(stderr," zmt      = %g \n",litho->zmt      );
        fprintf(stderr," zn       = %g \n",litho->zn      );
        fprintf(stderr," zy       = %g \n",litho->zy      );
}
