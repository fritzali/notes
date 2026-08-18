/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Define the components of the diagonal resistive tensor. 

  Use this function to supply the resistivity in the three directions
  \f$ \eta_{x1}\f$, \f$ \eta_{x2}\f$ and \f$ \eta_{x3}\f$.
  
  \authors T. Matsakos \n
           A. Mignone (mignone@ph.unito.it)\n
  \date    March 22, 2013
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void Resistive_eta(double *v, double x1, double x2, double x3,
                   double *J, double *eta)
/*!
 * Compute the resistive tensor components as function of the primitive
 * variables, coordinates and currents.
 *
 * \param [in]  v    array of primitive variables
 * \param [in]  x1   coordinate in the X1 direction
 * \param [in]  x2   coordinate in the X2 direction
 * \param [in]  x3   coordinate in the X3 direction
 * \param [in]  J    current components, J[IDIR], J[JDIR], J[KDIR]
 * \param [out] eta  an array containing the three components of
 *                   \f$ \tens{\eta}\f$.
 *
 *********************************************************************** */
{
 double coeff, cs, eps2, rcyl, beta, Bg[3], Bpol2, eta0;

 rcyl=x1*sin(x2);
 eps2=g_inputParam[EPS]*g_inputParam[EPS];
 coeff=2./5./eps2*(1./x1-(1.-5./2.*eps2)/rcyl);
 coeff = MAX(coeff,0.0);
 cs=eps2*coeff; /* initial sound speed*/

 #if (BACKGROUND_FIELD == YES)
       BackgroundField (x1,x2,x3,Bg);
       Bpol2=(v[BX1]+Bg[0])*(v[BX1]+Bg[0])+(v[BX2]+Bg[1])*(v[BX2]+Bg[1]);
 #else
  Bpol2=(v[BX1])*(v[BX1])+(v[BX2])*(v[BX2]);
 #endif

  beta=2.*v[PRS]/Bpol2;
      if (beta>0.5) {
       eta0=g_inputParam[ALPHAV]
       *cs*sqrt(rcyl*rcyl*rcyl)*v[TRC];
  }else {
  eta0 = 0.0;
    }	
 eta[IDIR] = eta0; 
 eta[JDIR] = eta0;
 eta[KDIR] = eta0;

}