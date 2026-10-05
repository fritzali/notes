/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Build the right hand side for the resistivity operator.

  Compute the one-dimensional right hand side for the resistive
  (magnetic diffusion) operator in the direction given by ::g_dir.

  Physically, this operator advances the induction equation's
  diffusive term,
  \f[
     \pd{\vec{B}}{t} = -\nabla\times\left(\eta\,\vec{J}\right),
  \f]
  together with the associated Ohmic/Joule heating of the gas. The
  update is written in conservative, flux-divergence form: a resistive
  flux \c res_flx (computed by ResistiveFlux(), using the diffusivity
  from Resistive_eta()) is differenced across each cell to produce the
  increment to the magnetic field (and, if energy is evolved, to the
  total energy).

  =====================================================================
  With radiation enabled, the ad hoc sink that radiated away the Ohmic
  heating of disk material has been removed. The heat stays in the gas
  and is exchanged with the radiation field by the radiative transfer
  module, which is what cools the disk. The file is kept as a local copy
  so that it lines up with the non-radiative setups, which still carry
  that sink, and is otherwise identical to the PLUTO original.
  =====================================================================

  \authors A. Mignone (andrea.mignone@unito.it)\n

  \b References

  \date   July 10, 2019
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ********************************************************************* */
void ResistiveRHS (const Data *d, Data_Arr dU, double **dcoeff,
                    double **aflux, double dt, int beg, int end, Grid *grid)
/*!
 * \param [in]   d           pointer to PLUTO Data structure
 * \param [out]  dU          a 4D array containing conservative variables
 *                           increment
 * \param [out]  dcoeff      1D array of diffusion coefficients
 * \param [out]  aflux       pointer to 2D array for AMR re-fluxing
 *                           operations
 * \param [in]   dt          the current time-step
 * \param [in]  beg,end      initial and final interface indices
 * \param [in]  grid         pointer to Grid structure.
 *********************************************************************** */
{
  int i, j, k, nv;
  double *x1 = grid->x[IDIR], *x1p = grid->xr[IDIR], *dx1 = grid->dx[IDIR];
  double *x2 = grid->x[JDIR], *x2p = grid->xr[JDIR], *dx2 = grid->dx[JDIR];
  double *x3 = grid->x[KDIR], *x3p = grid->xr[KDIR], *dx3 = grid->dx[KDIR];
  double *sp = grid->sp;
  double A, dtdV, dtdl, q, rhs[NVAR];
  static double **res_flx, **fxA;
  intList var_list;

  #if HAVE_ENERGY
  var_list.nvar = 4;
  var_list.indx[i=0] = BX1;
  var_list.indx[++i] = BX2;
  var_list.indx[++i] = BX3;
  var_list.indx[++i] = ENG;
  #else
  var_list.nvar = 3;
  var_list.indx[i=0] = BX1;
  var_list.indx[++i] = BX2;
  var_list.indx[++i] = BX3;
  #endif

/* --------------------------------------------------------
   0. Allocate memory
   -------------------------------------------------------- */

  if (res_flx == NULL){
    res_flx = ARRAY_2D(NMAX_POINT, NVAR, double);
    fxA     = ARRAY_2D(NMAX_POINT, NVAR, double);
  }

/* --------------------------------------------------------
   1. Add resistivity flux and source terms to
      total flux (sweep->flux) and total source terms.
   -------------------------------------------------------- */

  i = g_i; j = g_j; k = g_k;
  ResistiveFlux (d->Vc, d->J, res_flx, dcoeff, beg-1, end, grid);

  if (g_dir == IDIR){

  /* --------------------------------------------------------
     1a. Compute fluxes & sources in the X1 direction.

         The resistive flux is multiplied by the appropriate
         interface area element A (or, for the B_theta/B_phi
         components in spherical geometry, by the cylindrical-like
         radius x1p, following PLUTO's convention for constructing
         a divergence-free update of the induction equation on a
         non-Cartesian mesh).
     -------------------------------------------------------- */

    for (i = beg-1; i <= end; i++){
      A = grid->A[IDIR][k][j][i];

      #if GEOMETRY != SPHERICAL
      FOR_EACH(nv, &var_list) fxA[i][nv] = A*res_flx[i][nv];
      #endif
      #if GEOMETRY == SPHERICAL
      fxA[i][iBR]   = A*res_flx[i][iBR];
      fxA[i][iBTH]  = x1p[i]*res_flx[i][iBTH];
      fxA[i][iBPHI] = x1p[i]*res_flx[i][iBPHI];
      #if HAVE_ENERGY
      fxA[i][ENG] = A*res_flx[i][ENG];
      #endif
      #endif
    }

    /* -- 1b. Build rhs in the X1-direction --

           Standard finite-volume flux differencing: the increment
           to each conserved variable is dt/dV times the difference
           of the (area-weighted) flux across the cell's left and
           right interfaces. In spherical geometry, B_r has no
           diffusive flux divergence contribution here (rhs[iBR]=0,
           since ResistiveFlux updates it through the transverse
           components instead), while B_theta and B_phi pick up an
           extra 1/x1 factor from the curvilinear divergence.
       -- */

    for (i = beg; i <= end; i++){
      dtdV = dt/grid->dV[k][j][i];
      dtdl = dt/dx1[i];
      #if GEOMETRY == SPHERICAL
      q = dtdl/x1[i];
      rhs[iBR]   = 0.0;
      rhs[iBTH]  = q*(fxA[i][iBTH]  - fxA[i-1][iBTH]);
      rhs[iBPHI] = q*(fxA[i][iBPHI] - fxA[i-1][iBPHI]);
      #if HAVE_ENERGY
      rhs[ENG] = dtdV*(fxA[i][ENG] - fxA[i-1][ENG]);
      #endif
      #else
      FOR_EACH(nv, &var_list) rhs[nv] = dtdV*(fxA[i][nv] - fxA[i-1][nv]);
      #if (GEOMETRY == POLAR) || (GEOMETRY == CYLINDRICAL)
      rhs[iBPHI] = dtdl*(res_flx[i][iBPHI] - res_flx[i-1][iBPHI]);
      #endif
      #endif

      FOR_EACH(nv, &var_list) dU[k][j][i][nv] += rhs[nv];
    }

  }else if (g_dir == JDIR){

  /* --------------------------------------------------------
     2a. Compute fluxes & sources in the X2 direction
     -------------------------------------------------------- */

    for (j = beg-1; j <= end; j++){
      A = grid->A[JDIR][k][j][i];
      FOR_EACH(nv, &var_list) fxA[j][nv] = A*res_flx[j][nv];
    }

  /* -- 2b. Build rhs in the X2-direction --

         Same flux-differencing logic as in the X1 sweep; in
         spherical geometry B_phi again requires an extra metric
         factor (dx_dl, the arclength-to-coordinate-increment ratio)
         to account for the curvilinear divergence in theta.
     -- */

    double **dx_dl = grid->dx_dl[JDIR];
    for (j = beg; j <= end; j++){
      dtdV = dt/grid->dV[k][j][i];
      FOR_EACH(nv, &var_list) rhs[nv] = dtdV*(fxA[j][nv] - fxA[j-1][nv]);
      #if GEOMETRY == SPHERICAL
      dtdl = dt/dx2[j]*dx_dl[j][i];
      rhs[iBPHI] = dtdl*(res_flx[j][iBPHI] - res_flx[j-1][iBPHI]);
      #endif

      FOR_EACH(nv, &var_list) dU[k][j][i][nv] += rhs[nv];
    }

  }else if (g_dir == KDIR){

  /* --------------------------------------------------------
     3a. Compute fluxes & sources in the X3 direction
     -------------------------------------------------------- */

    for (k = beg-1; k <= end; k++){
      A = grid->A[KDIR][k][j][i];
      FOR_EACH(nv, &var_list)  fxA[k][nv] = A*res_flx[k][nv];
    }

  /* -- 3b. Build rhs in the X3-direction -- */

    for (k = beg; k <= end; k++){
      dtdV = dt/grid->dV[k][j][i];
      FOR_EACH(nv, &var_list) rhs[nv] = dtdV*(fxA[k][nv] - fxA[k-1][nv]);

      FOR_EACH(nv, &var_list) dU[k][j][i][nv] += rhs[nv];
    }
  }

  #ifdef CHOMBO
  StoreAMRFlux (res_flx, aflux,-1, BX1, BX1+2, beg-1, end,grid);
  #if HAVE_ENERGY
  StoreAMRFlux (res_flx, aflux,-1, ENG, ENG, beg-1, end, grid);
  #endif
  #endif
}
