/* /////////////////////////////////////////////////////////////////// */
/*! \file  
 *  \brief Specification of explicit first and second viscosity coefficients*/
/* /////////////////////////////////////////////////////////////////// */
#include "pluto.h"
/* ************************************************************************** */
void Visc_nu(double *v, double x1, double x2, double x3,
                        double *nu1, double *nu2)
/*! 
 *
 *  \param [in]      v  pointer to data array containing cell-centered quantities
 *  \param [in]      x1 real, coordinate value 
 *  \param [in]      x2 real, coordinate value 
 *  \param [in]      x3 real, coordinate value 
 *  \param [in, out] nu1  pointer to first viscous coefficient
 *  \param [in, out] nu2  pointer to second viscous coefficient
 *
 *  \return This function has no return value.
 * ************************************************************************** */
{
double coeff, cs, eps2, rcyl, beta, Bg[3], Bpol2;

 rcyl=x1*sin(x2);
 eps2=g_inputParam[EPS]*g_inputParam[EPS];
 coeff=2./5./eps2*(1./x1-(1.-5./2.*eps2)/rcyl);
 coeff = MAX(coeff,0.0);
 cs=eps2*coeff; /* initial sound speed*/
  
#if PHYSICS == MHD

 #if (BACKGROUND_FIELD == YES)
       BackgroundField (x1,x2,x3,Bg);
       Bpol2=(v[BX1]+Bg[0])*(v[BX1]+Bg[0])+(v[BX2]+Bg[1])*(v[BX2]+Bg[1]);
 #else
  Bpol2=(v[BX1])*(v[BX1])+(v[BX2])*(v[BX2]);
 #endif

  beta=2.*v[PRS]/Bpol2;

      if (beta>0.5) {

       *nu1=2./3.*v[RHO]*g_inputParam[BETAV]*v[TRC]*sqrt(rcyl);
  }else {
 *nu1 = 0.0;
    }
#else
       double fact=(v[RHO]>1e-3 && v[VX3]>0.05) ? 1 : 0;
       *nu1=fact*2./3.*v[RHO]*g_inputParam[BETAV]*sqrt(rcyl);
#endif    
 *nu2 = 0.0;
}
