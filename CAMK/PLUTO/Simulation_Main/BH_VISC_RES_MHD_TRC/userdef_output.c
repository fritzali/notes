/* /////////////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief User-defined output routines for PLUTO simulations.

  Computes diagnostic user-defined output variables (kinematic viscosity nu,
  magnetic diffusivity num, temperature ratio Te, and disk fraction diskfrac)
  and registers them for file output (.dbl / VTK).

  \author A. Mignone (mignone@ph.unito.it)
  \date   Sep 2012

  \modified M. Cemeljic (miki@camk.edu.pl)
  \date   Apr 2021 / Modified 2024 (ccm)
*/
/* /////////////////////////////////////////////////////////////////////////// */

#include "pluto.h"

/* Optional forward declaration for disk fraction computation */
extern double DiskFraction(double *v, double x1, double x2);

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

#if PARTICLES
  /* Optional particle output variable configuration */
  /* SetOutputVar("energy", PARTICLES_FLT_OUTPUT, NO); */
  /* SetOutputVar("x1",     PARTICLES_FLT_OUTPUT, NO); */
  /* SetOutputVar("vx1",    PARTICLES_FLT_OUTPUT, NO); */
#endif
}
