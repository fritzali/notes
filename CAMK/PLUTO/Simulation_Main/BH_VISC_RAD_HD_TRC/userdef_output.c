/* /////////////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief User-defined output routines for PLUTO simulations.

  Computes diagnostic user-defined output variables (kinematic viscosity nu,
  magnetic diffusivity num, temperature ratio Te, disk fraction diskfrac,
  and absorption/scattering opacities kappa_abs / kappa_scat) and registers
  them for file output (.dbl / VTK).

  \author A. Mignone (mignone@ph.unito.it)
  \date   Sep 2012

  \modified M. Cemeljic (miki@camk.edu.pl)
  \date   Apr 2021 / Modified 2024 (ccm)

  FURTHER MODIFIED: added kappa_abs / kappa_scat diagnostic output, via
  UserDefOpacitiesAt(v, x1, x2, ...) (defined in init.c), called directly
  with this loop's own (x1[i], x2[j]) rather than through the g_i_rad/g_j
  globals used by the radiation module's call path - those globals would
  be stale here (see init.c's UserDefOpacitiesAt() docstring).
*/
/* /////////////////////////////////////////////////////////////////////////// */

#include "pluto.h"
#include "modifications.h"

/* ********************************************************************* */
void ComputeUserVar (const Data *d, Grid *grid)
/*!
 * Compute user-defined diagnostic output variables across the active domain.
 *********************************************************************** */
{
  int i, j, k;
  int nv;
  double vi[NVAR];
  double nu1, nu2;
  double ***nu;
  double ***Te;
  double ***diskfrac;

#if PHYSICS == MHD
  double J[3], eta[3];
  double ***num;
#endif

#if RADIATION_VAR_OPACITIES
  double ***kappa_abs;
  double ***kappa_scat;
#endif

  double *x1 = grid->x[IDIR];
  double *x2 = grid->x[JDIR];
  double *x3 = grid->x[KDIR];

  /* Retrieve pointers for registered user-defined output arrays */
  nu       = GetUserVar("nu");
  Te       = GetUserVar("Te");
  diskfrac = GetUserVar("diskfrac");

#if PHYSICS == MHD
  num      = GetUserVar("num");
#endif

#if RADIATION_VAR_OPACITIES
  kappa_abs  = GetUserVar("kappa_abs");
  kappa_scat = GetUserVar("kappa_scat");
#endif

  DOM_LOOP(k, j, i) {
    /* Copy cell-centered primitive variables into local array */
    for (nv = 0; nv < NVAR; nv++) {
      vi[nv] = d->Vc[nv][k][j][i];
    }

    /* Kinematic viscosity calculation */
    Visc_nu(vi, x1[i], x2[j], x3[k], &nu1, &nu2);
    nu[k][j][i] = nu1;

#if PHYSICS == MHD
    /* Magnetic resistivity calculation */
    Resistive_eta(vi, x1[i], x2[j], x3[k], J, eta);
    num[k][j][i] = eta[0];
#endif

    /* Temperature ratio P / RHO */
    Te[k][j][i] = vi[PRS] / vi[RHO];

    /* Disk fraction computation */
    diskfrac[k][j][i] = DiskFraction(vi, x1[i], x2[j]);

#if RADIATION_VAR_OPACITIES
    /* Absorption/scattering opacities, evaluated at this cell's own
       (x1[i], x2[j]) directly - NOT via g_i_rad/g_j, which belong to
       the radiation module's own call path and would be stale here. */
    {
      double a_op, s_op;
      UserDefOpacitiesAt(vi, x1[i], x2[j], &a_op, &s_op);
      kappa_abs[k][j][i]  = a_op;
      kappa_scat[k][j][i] = s_op;
    }
#endif
  }
}

/* ********************************************************************* */
void ChangeOutputVar ()
/*!
 * Enable or disable user-defined diagnostic output variables for file output.
 *********************************************************************** */
{
  /* Enable standard user-defined variables for output (.dbl / VTK) */
  SetOutputVar("nu", DBL_OUTPUT, YES);
  
#if PHYSICS == MHD
  SetOutputVar("num", DBL_OUTPUT, YES);
#endif

  SetOutputVar("Te", DBL_OUTPUT, YES);
  SetOutputVar("diskfrac", DBL_OUTPUT, YES);

#if RADIATION_VAR_OPACITIES
  SetOutputVar("kappa_abs",  DBL_OUTPUT, YES);
  SetOutputVar("kappa_scat", DBL_OUTPUT, YES);
#endif

#if PARTICLES
  /* Optional particle output variable configuration */
  /* SetOutputVar("energy", PARTICLES_FLT_OUTPUT, NO); */
  /* SetOutputVar("x1",     PARTICLES_FLT_OUTPUT, NO); */
  /* SetOutputVar("vx1",    PARTICLES_FLT_OUTPUT, NO); */
#endif
}
