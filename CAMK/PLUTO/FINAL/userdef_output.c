#include "pluto.h"
#include "modifications.h"

/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 *
 *
 ***************************************************************** */
{
  int i, j, k;
//ccm--210421
  int nv;
  double vi[NVAR];
  double nu1, nu2;
  double ***nu;
  double ***Te;
  double ***diskfrac;   //cfa

#if PHYSICS == MHD
  double J[3], eta[3];
  double ***num;
#endif

#if RADIATION_VAR_OPACITIES
  double ***kappa_abs;   //cfa
  double ***kappa_scat;
#endif

  double *x1 = grid->x[IDIR];
  double *x2 = grid->x[JDIR];
  double *x3 = grid->x[KDIR];

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
//ccm

//cfa--classifier state with every output, for restarts (init.c)
  SaveDiskProfiles();

  DOM_LOOP(k, j, i) {

//ccm--210421
    for (nv = 0; nv < NVAR; nv++) {
      vi[nv] = d->Vc[nv][k][j][i];
    }

    Visc_nu(vi, x1[i], x2[j], x3[k], &nu1, &nu2);
    nu[k][j][i] = nu1;

#if PHYSICS == MHD
    Resistive_eta(vi, x1[i], x2[j], x3[k], J, eta);
    num[k][j][i] = eta[0];
#endif

    Te[k][j][i] = vi[PRS] / vi[RHO];

//cfa
    diskfrac[k][j][i] = DiskFraction(vi, x1[i], x2[j]);

#if RADIATION_VAR_OPACITIES
//cfa--opacities at this cell, g_i_rad/g_j would be stale here
    {
      double a_op, s_op;
      UserDefOpacitiesAt(vi, x1[i], x2[j], &a_op, &s_op);
      kappa_abs[k][j][i]  = a_op;
      kappa_scat[k][j][i] = s_op;
    }
#endif
//ccm
  }
}
/* ************************************************************* */
void ChangeOutputVar ()
/* 
 *
 * 
 *************************************************************** */
{
//ccm--210421--additional variables output
  SetOutputVar("nu", DBL_OUTPUT, YES);

#if PHYSICS == MHD
  SetOutputVar("num", DBL_OUTPUT, YES);
#endif

  SetOutputVar("Te", DBL_OUTPUT, YES);
  SetOutputVar("diskfrac", DBL_OUTPUT, YES);   //cfa

#if RADIATION_VAR_OPACITIES
  SetOutputVar("kappa_abs",  DBL_OUTPUT, YES);
  SetOutputVar("kappa_scat", DBL_OUTPUT, YES);
#endif

#if PARTICLES
  //SetOutputVar ("energy",PARTICLES_FLT_OUTPUT, NO);
//  SetOutputVar ("x1",    PARTICLES_FLT_OUTPUT, NO);
  //SetOutputVar ("vx1",   PARTICLES_FLT_OUTPUT, NO);
#endif
}
