/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Specification of explicit first and second viscosity coefficients.

  This function returns the dynamical (\c nu1) and bulk (\c nu2) viscosity
  coefficients used by PLUTO's explicit viscosity module to build the
  viscous stress tensor,
  \f[
     \Pi_{ij} = -\rho\nu_1\left(\partial_i v_j+\partial_j v_i
                -\frac{2}{3}\delta_{ij}\nabla\cdot\vec{v}\right)
                -\rho\nu_2\,\delta_{ij}\,\nabla\cdot\vec{v} \,.
  \f]

  The dynamical viscosity is modeled with the standard Shakura-Sunyaev
  \f$\alpha\f$-prescription for turbulent angular-momentum transport in
  accretion disks,
  \f[
     \nu_1 = \frac{2}{3}\,\alpha\, c_s\, H \,,\qquad H \sim R_{\rm cyl}^{3/2}
  \f]
  where \f$ c_s \f$ is the local (initial) sound speed of the disk and
  \f$ R_{\rm cyl}^{3/2} \f$ plays the role of a pressure-scale-height-like
  length scale appropriate to the adopted disk equilibrium (see comments
  in the body of the function). The coefficient is switched off in regions
  where the magnetic pressure dominates over the thermal pressure
  (plasma \f$\beta \le 0.5\f$), so that \f$\alpha\f$-viscosity only acts
  in weakly magnetized / thermally-supported parts of the flow (e.g. a
  disk corona or current-free funnel should not be artificially
  viscously heated).

  Disk membership is now determined by DiskFraction() (see init.c)
  instead of the passive tracer v[TRC], since the tracer field is
  unreliable in this setup (see project notes). DiskFraction() uses a
  local density-vs-corona-profile threshold combined with a rotation
  criterion, and requires no MPI reduction.

  \authors T. Matsakos \n
           A. Mignone (mignone@ph.unito.it)
  \date    March 22, 2013
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include "modifications.h"

/* ********************************************************************* */
void Visc_nu(double *v, double x1, double x2, double x3,
             double *nu1, double *nu2)
/*!
 * Compute the first (\c nu1) and second (\c nu2) viscosity coefficients
 * as a function of the primitive variables and spatial position.
 *
 * \param [in]      v    array of primitive variables
 * \param [in]      x1   coordinate in the X1 direction (spherical radius)
 * \param [in]      x2   coordinate in the X2 direction (polar angle)
 * \param [in]      x3   coordinate in the X3 direction (azimuth)
 * \param [in,out]  nu1  pointer to the first (shear/dynamical) viscosity
 *                       coefficient
 * \param [in,out]  nu2  pointer to the second (bulk) viscosity coefficient
 *
 * \return  This function has no return value.
 *********************************************************************** */
{
  double coeff, cs, eps2, rcyl, beta, Bg[3], Bpol2;
  double disk_frac;

/* --------------------------------------------------------
   0. Geometry and background thermal profile.

      rcyl is the cylindrical radius R = r*sin(theta), i.e. the
      distance from the disk rotation axis in spherical coordinates.

      eps2 = eps^2, where eps = H/R is the (constant) disk aspect
      ratio input parameter. The combination

        coeff = (2/5/eps^2) * [ 1/r - (1 - 5/2*eps^2)/rcyl ]

      reconstructs the initial (equilibrium) enthalpy profile of a
      rotating, eps-thin torus in spherical coordinates, consistent
      with the disk setup used to initialize the simulation
      (a polytropic, radiatively-inefficient torus threaded by a
      near-Keplerian rotation law). Clipping to zero (MAX(coeff,0.0))
      prevents negative sound speeds outside the initial disk body
      (e.g. above the disk surface or inside evacuated funnel
      regions) where the analytic profile is not meant to apply.

      cs = eps^2 * coeff then gives the square of the local sound
      speed associated with that initial disk temperature profile;
      it is evaluated locally from (x1,x2) rather than read off the
      evolving flow, so it represents the background/initial thermal
      scale used to calibrate the viscosity - not the instantaneous
      sound speed of the evolved gas.
   -------------------------------------------------------- */

  rcyl = x1*sin(x2);
  eps2 = g_inputParam[EPS]*g_inputParam[EPS];
  coeff = 2./5./eps2*(1./x1 - (1. - 5./2.*eps2)/rcyl);
  coeff = MAX(coeff, 0.0);
  cs = eps2*coeff;  /* initial sound speed (squared, in code units) */

  disk_frac = DiskFraction(v, x1, x2);

#if PHYSICS == MHD

/* --------------------------------------------------------
   1. Plasma-beta gate (MHD only).

      Compute the poloidal magnetic energy density (optionally adding
      a split/background field component Bg via BackgroundField()),
      then form beta = 2*p/Bpol^2, i.e. the ratio of thermal to
      poloidal magnetic pressure.

      alpha-viscosity is only switched on where beta > 0.5, i.e.
      where the flow is not strongly magnetically dominated. This
      keeps the explicit (hydrodynamic-analog) turbulent viscosity
      confined to thermally-supported disk material, and avoids
      injecting spurious dissipation in magnetically-dominated
      regions (corona, jet/funnel) where angular-momentum transport
      is already handled by the resolved MHD stresses.
   -------------------------------------------------------- */

  #if (BACKGROUND_FIELD == YES)
    BackgroundField(x1, x2, x3, Bg);
    Bpol2 = (v[BX1]+Bg[0])*(v[BX1]+Bg[0]) + (v[BX2]+Bg[1])*(v[BX2]+Bg[1]);
  #else
    Bpol2 = (v[BX1])*(v[BX1]) + (v[BX2])*(v[BX2]);
  #endif

  beta = 2.*v[PRS]/Bpol2;

  if (beta > 0.5) {

  /* -- 1a. alpha-viscosity coefficient.

        nu1 = (2/3) * rho * alpha * cs * rcyl^(3/2) * disk_frac

        The rcyl^(3/2) factor plays the role of a local disk
        scale-height-like length (consistent with the initial
        torus equilibrium above), so that nu1 has the standard
        alpha-disk scaling nu ~ alpha * cs * H. disk_frac (valued
        0 or 1, from DiskFraction()) restricts the viscosity to
        zones identified as disk material (as opposed to the
        ambient/funnel/outflow gas), replacing the passive tracer
        v[TRC] previously used for this purpose.
     -- */

    *nu1 = 2./3.*v[RHO]*g_inputParam[ALPHAV]
           *cs*sqrt(rcyl*rcyl*rcyl)*disk_frac;
  } else {
    *nu1 = 0.0;
  }

#else /* PHYSICS != MHD : purely hydrodynamic disk */

/* --------------------------------------------------------
   1'. Pure hydro case: no magnetic field to gate the viscosity,
       so the alpha-prescription is applied unconditionally
       (still restricted to disk material via disk_frac).
   -------------------------------------------------------- */

  *nu1 = 2./3.*v[RHO]*g_inputParam[ALPHAV]
         *cs*sqrt(rcyl*rcyl*rcyl)*disk_frac;

#endif

/* --------------------------------------------------------
   2. Bulk viscosity.

      No bulk viscosity is modeled; only the shear (dynamical)
      viscosity nu1 contributes to angular-momentum transport.
   -------------------------------------------------------- */

  *nu2 = 0.0;
}
