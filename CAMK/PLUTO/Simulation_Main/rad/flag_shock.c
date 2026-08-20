/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Shock finding algorithm.

  Search and flag computational zones lying in a shock wave.
  The flagging strategy is based on two switches designed to detect
  the presence of compressive motion or shock waves in the fluid:
  \f[
     \nabla\cdot\vec{v} < 0 \qquad{\rm and}\qquad
     \Delta x\frac{|\nabla p|}{p} > \epsilon_p
  \f]
  where \f$\epsilon_p\f$ sets the shock strength.
  At the discrete level we replace the two conditions by
  \f[
    \sum_d \frac{ A_{\vec{i}+\HALF\hvec{e}_d}v_{d,\vec{i}+\HALF\hvec{e}_d}
                 -A_{\vec{i}-\HALF\hvec{e}_d}v_{d,\vec{i}-\HALF\hvec{e}_d} }
                {\Delta{\cal V}_{d,\vec{i}}}  < 0
                \qquad{\rm and}\qquad
    \sum_{d} \left|p_{\vec{i}+\hvec{e}_d} - p_{\vec{i}-\hvec{e}_d}\right|
             <
     \epsilon_p \min_d\left(p_{\vec{i}+\hvec{e}_d},
                            p_{\vec{i}-\hvec{e}_d},p_{\vec{i}}\right)
  \f]
  where \f$\hvec{i} = (i,j,k)\f$ is a vector of integer numbers
  giving the position of a computational zone, while \f$\hvec{e}_d =
  (\delta_{1d},\delta_{2d},\delta_{3d})\f$ is a unit vector in the direction
  given by \c d.
  Once a zone has been tagged as lying in a shock, different flags may be
  switched on or off to control the update strategy in these critical regions.

  This function can be called called when:
  - \c SHOCK_FLATTENING has been set to \c MULTID: in this case shocked zones
    are tagged with \c FLAG_MINMOD and \c FLAG_HLL that will later
    be used to force the reconstruction with the minmod limiter
    and the Riemann solver to HLL.
  - \c ENTROPY_SWITCH has been turned on: this flag will be checked later in
    the ConsToPrim() functions in order to recover pressure from the
    entropy density rather than from the total energy density.
    The update process is:

    - start with a vector of primitive variables  <tt> {V,s} </tt>
      where \c s is the entropy;
    - set boundary condition on \c {V};
    -  compute \c {s} from \c {V};
    - flag zones where entropy may be used (flag = 1);
    - evolve equations for one time step;
    - convert <tt> {U,S} </tt> to primitive:
      \code
        if (flag == 0) {  // Use energy
          p = p(E)
          s = s(p)
        }else{            // use entropy
          p = p(S)
          E = E(p)
        }
        \endcode

  \b Reference
     - "Maintaining Pressure Positivity in Magnetohydrodynamics Simulations"
       Balsara \& Spicer, JCP (1999) 148, 133

  \authors A. Mignone (andrea.mignone@unito.it)
  \date    Feb 28, 2023
*/
/* ///////////////////////////////////////////////////////////////////// */
#include"pluto.h"

#ifndef EPS_PSHOCK_FLATTEN
 #define EPS_PSHOCK_FLATTEN  5.0
#endif

#ifndef EPS_PSHOCK_ENTROPY
 #define EPS_PSHOCK_ENTROPY  0.05
#endif

#define NBUF   2

#if (SHOCK_FLATTENING == MULTID) || ENTROPY_SWITCH
/* *************************************************************** */
void FlagShock (const Data *d, Grid *grid)
/*!
 * \param [in,out] d     pointer to data structure
 * \param [in]     grid  pointer to grid structure
 *
 *  csk: fully rewritten to mimic P4.41
 ***************************************************************** */
{
int  i, j, k, nv;
  double eint, etot;
  double ***pt, ***rho, ***vx1, ***vx2, ***vx3, ***bx1, ***bx2, ***bx3, Bg[3];

/* -------------------------------------------------
   1. Define pointers to variables
   ------------------------------------------------- */
  rho = d->Vc[RHO];

  vx1 = d->Vc[VX1];
  vx2 = d->Vc[VX2];
  vx3 = d->Vc[VX3];

# if PHYSICS == MHD
  bx1 = d->Vc[BX1];
  bx2 = d->Vc[BX2];
  bx3 = d->Vc[BX3];
# endif  
  pt = d->Vc[PRS];

  for (k = INCLUDE_KDIR; k < NX3_TOT-INCLUDE_KDIR; k++){
  for (j = INCLUDE_JDIR; j < NX2_TOT-INCLUDE_JDIR; j++){
  for (i = INCLUDE_IDIR; i < NX1_TOT-INCLUDE_IDIR; i++){

#if (PHYSICS == MHD) && (BACKGROUND_FIELD == YES)
  BackgroundField (grid->x[IDIR][i],grid->x[JDIR][j],grid->x[KDIR][k],Bg);
#else
  Bg[0]=0;Bg[1]=0;Bg[2]=0;
#endif

       eint = pt[k][j][i]/(g_gamma-1.);
       etot = eint + 0.5*rho[k][j][i]*(vx1[k][j][i]*vx1[k][j][i]+vx2[k][j][i]*vx2[k][j][i]+vx3[k][j][i]*vx3[k][j][i]);

#if PHYSICS == MHD
       etot += 0.5*(pow(bx1[k][j][i]+Bg[0],2)+pow(bx2[k][j][i]+Bg[1],2)+pow(bx3[k][j][i]+Bg[2],2));
# endif

       if (eint < 1.e-2*etot) {
         d->flag[k][j][i]   |= FLAG_HLL;
       }
  }
  }
  }
#ifdef PARALLEL
  AL_Exchange (d->flag[0][0], SZ_uint16_t);
#endif

}
#undef  EPS_PSHOCK_ENTROPY
#endif
