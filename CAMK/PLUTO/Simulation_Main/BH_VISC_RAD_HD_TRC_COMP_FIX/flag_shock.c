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

  \b Note (this implementation): the body of FlagShock() below has been
  rewritten relative to PLUTO's standard multidimensional switch
  described above (see the "csk" note in the source). Rather than
  evaluating the discrete velocity-divergence and pressure-jump
  criteria cell-by-cell, this version flags a zone using a single,
  cheaper local diagnostic: the ratio of internal (thermal) energy
  density to total energy density (thermal + kinetic + magnetic).
  Zones where the internal energy is a very small fraction of the
  total energy (\c eint \c < \c 1.e-2*etot) are exactly the zones
  where recovering pressure from the total energy is numerically
  fragile (catastrophic cancellation between large kinetic/magnetic
  and total energy terms), which is the typical situation deep inside
  a strong shock or in a highly super-magnetosonic/high-Mach-number
  region (e.g. close to the black hole, in the funnel, or in strongly
  magnetized disk regions). Such zones are marked with \c FLAG_HLL so
  that the more diffusive, positivity-preserving HLL Riemann solver
  is used there instead of a higher-order solver that may return
  negative pressures.

  \b Reference
     - "Maintaining Pressure Positivity in Magnetohydrodynamics Simulations"
       Balsara \& Spicer, JCP (1999) 148, 133

  \authors A. Mignone (andrea.mignone@unito.it)
  \date    Feb 28, 2023
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

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
 *  \note  fully rewritten (csk) to mimic the pressure-positivity
 *         safeguard of P4.41: replaces the standard multidimensional
 *         divergence/pressure-jump shock switch with a direct
 *         internal-energy-fraction test (see file-level note above).
 ***************************************************************** */
{
  int  i, j, k, nv;
  double eint, etot;
  double ***pt, ***rho, ***vx1, ***vx2, ***vx3, ***bx1, ***bx2, ***bx3, Bg[3];

/* -------------------------------------------------
   1. Define pointers to variables.

      pt, rho, vx1-3 are the thermal pressure, density and the
      three velocity components; bx1-3 are the magnetic field
      components (MHD only). These are raw pointers into the
      Data structure's primitive-variable array, indexed as
      [k][j][i] (X3, X2, X1).
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

  /* -------------------------------------------------
     2. Background field contribution.

        In the split (BACKGROUND_FIELD) formulation, the total
        physical field is the sum of the evolved perturbation
        (bx1-3) and a fixed analytic background Bg (e.g. a
        large-scale field threading the disk/black hole
        magnetosphere). This must be added back in when computing
        the true total (magnetic) energy density used below.
     ------------------------------------------------- */

#if (PHYSICS == MHD) && (BACKGROUND_FIELD == YES)
    BackgroundField (grid->x[IDIR][i],grid->x[JDIR][j],grid->x[KDIR][k],Bg);
#else
    Bg[0]=0; Bg[1]=0; Bg[2]=0;
#endif

  /* -------------------------------------------------
     3. Internal-energy-fraction diagnostic.

        eint = p/(gamma-1)                         (thermal energy density)
        etot = eint + (1/2)*rho*v^2 [+ (1/2)*B^2]  (total energy density,
                                                     kinetic [+ magnetic])

        A zone is flagged (FLAG_HLL) whenever the thermal energy is
        less than 1% of the total energy density. In such zones the
        kinetic/magnetic energy dominates so strongly that computing
        pressure as p = (gamma-1)*(E_tot - E_kin - E_mag) from the
        evolved total energy E is prone to large cancellation errors
        and can yield negative or noisy pressures. Using HLL (a more
        diffusive but positivity-preserving Riemann solver) in these
        zones trades some accuracy for robustness, consistent with
        the Balsara & Spicer (1999) pressure-positivity strategy
        referenced above.
     ------------------------------------------------- */

     eint = pt[k][j][i]/(g_gamma-1.);
     etot = eint + 0.5*rho[k][j][i]*(vx1[k][j][i]*vx1[k][j][i]
                                     +vx2[k][j][i]*vx2[k][j][i]
                                     +vx3[k][j][i]*vx3[k][j][i]);

#if PHYSICS == MHD
     etot += 0.5*(pow(bx1[k][j][i]+Bg[0],2)
                 +pow(bx2[k][j][i]+Bg[1],2)
                 +pow(bx3[k][j][i]+Bg[2],2));
# endif

     if (eint < 1.e-2*etot) {
       d->flag[k][j][i] |= FLAG_HLL;
     }
  }
  }
  }

/* -------------------------------------------------
   4. Parallel boundary exchange of flags, so that
      neighboring MPI domains agree on which zones
      are flagged near domain boundaries.
   ------------------------------------------------- */

#ifdef PARALLEL
  AL_Exchange (d->flag[0][0], SZ_uint16_t);
#endif

}
#undef  EPS_PSHOCK_ENTROPY
#endif
