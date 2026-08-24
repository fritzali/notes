/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Define the components of the diagonal resistive tensor.

  Use this function to supply the resistivity in the three directions
  \f$ \eta_{x1}\f$, \f$ \eta_{x2}\f$ and \f$ \eta_{x3}\f$.

  The resistivity is modeled as an isotropic, alpha-type turbulent
  diffusivity, built with the same physical prescription used for the
  companion viscosity coefficient in Visc_nu() (see visc_nu.c): a
  Shakura-Sunyaev-like scaling
  \f[
     \eta = \alpha\, c_s\, R_{\rm cyl}^{3/2}
  \f]
  evaluated from the initial (equilibrium) disk sound speed, gated by
  the local plasma \f$\beta\f$ so that magnetic diffusivity is only
  active where the flow is thermally dominated (\f$\beta>0.5\f$), and
  restricted to disk material via DiskFraction() (see init.c), which
  replaces the passive tracer v[TRC] as the disk/corona discriminator
  used previously. \c eta is returned as a diagonal tensor (isotropic
  resistivity, i.e. \f$\eta_{x1}=\eta_{x2}=\eta_{x3}\f$).

  \authors T. Matsakos \n
           A. Mignone (mignone@ph.unito.it)\n
  \date    March 22, 2013
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include "modifications.h"

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
  double disk_frac;

/* --------------------------------------------------------
   0. Background thermal profile of the initial disk.

      Same construction as in Visc_nu(): rcyl = r*sin(theta) is the
      cylindrical radius, eps = H/R is the disk aspect ratio, and
      "coeff" reconstructs the equilibrium enthalpy profile of the
      initial eps-thin rotating torus. Clipping to zero avoids
      negative/undefined sound speeds outside the initial disk body.
      cs here is the square of the local background sound speed used
      purely to set the magnitude of the diffusivity - not the
      instantaneous sound speed of the evolved gas.
   -------------------------------------------------------- */

  rcyl = x1*sin(x2);
  eps2 = g_inputParam[EPS]*g_inputParam[EPS];
  coeff = 2./5./eps2*(1./x1 - (1. - 5./2.*eps2)/rcyl);
  coeff = MAX(coeff, 0.0);
  cs = eps2*coeff;  /* initial sound speed (squared, in code units) */

  disk_frac = DiskFraction(v, x1, x2);

/* --------------------------------------------------------
   1. Plasma-beta gate.

      Poloidal magnetic energy density (optionally including a
      split/background field component via BackgroundField()) is
      used to compute beta = 2*p/Bpol^2. Resistivity is switched on
      only where beta > 0.5, mirroring the viscosity prescription:
      the effective turbulent diffusivity should vanish in strongly
      magnetized regions where a mean-field/laminar description
      (rather than a turbulent alpha-closure) is more appropriate.
   -------------------------------------------------------- */

  #if (BACKGROUND_FIELD == YES)
    BackgroundField(x1, x2, x3, Bg);
    Bpol2 = (v[BX1]+Bg[0])*(v[BX1]+Bg[0]) + (v[BX2]+Bg[1])*(v[BX2]+Bg[1]);
  #else
    Bpol2 = (v[BX1])*(v[BX1]) + (v[BX2])*(v[BX2]);
  #endif

  beta = 2.*v[PRS]/Bpol2;

  if (beta > 0.5) {

  /* -- 1a. alpha-resistivity coefficient.

        eta0 = alpha * cs * rcyl^(3/2) * disk_frac

        Same alpha-disk scaling eta ~ alpha*cs*H used for the
        viscosity, with disk_frac (from DiskFraction()) restricting
        the diffusivity to zones identified as disk material,
        replacing the passive tracer v[TRC] used previously.
     -- */

    eta0 = g_inputParam[ALPHAV]*cs*sqrt(rcyl*rcyl*rcyl)*disk_frac;
  } else {
    eta0 = 0.0;
  }

/* --------------------------------------------------------
   2. Isotropic resistive tensor.

      The same scalar diffusivity eta0 is assigned to all three
      directions, i.e. the resistive tensor is diagonal and
      isotropic: eta_x1 = eta_x2 = eta_x3 = eta0.
   -------------------------------------------------------- */

  eta[IDIR] = eta0;
  eta[JDIR] = eta0;
  eta[KDIR] = eta0;
}
