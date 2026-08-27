/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Build the right hand side for the viscosity operator.

  Compute the one-dimensional right hand side for the viscous
  operator in the direction given by ::g_dir.

  Physically, this operator advances the momentum equation's viscous
  stress-divergence term,
  \f[
     \pd{(\rho\vec v)}{t} = \nabla\cdot\tens\Pi \,,
  \f]
  together with the associated viscous heating of the gas. As in
  res_rhs.c for the resistive operator, the update is written in
  conservative flux-divergence form using a viscous flux \c ViF and
  source term \c ViS (from ViscousFlux(), which internally uses the
  coefficients from Visc_nu()); a local, explicit viscous-heating
  correction is then subtracted from the internal energy using the
  viscous momentum flux components at the cell (see step "cooling"
  below), capped so it cannot remove more than a fixed fraction of the
  local thermal energy in a single step - the direct analog of the
  Ohmic cooling safeguard in ResistiveRHS().

  \authors A. Mignone (andrea.mignone@unito.it)\n

 \b References

  \date   Aug 21, 2019
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include "modifications.h"

/* ********************************************************************* */
void ViscousRHS (const Data *d, Data_Arr dU, double *dcoeff,
                  double **aflux, double dt, int beg, int end, Grid *grid)
/*!
 * \param [in]   d           pointer to PLUTO Data structure
 * \param [out]  dU          a 4D array containing conservative variables
 *                           increment
 * \param [out]  dcoeff      1D array of diffusion coefficients
 * \param [out]  aflux       pointer to 2D array for AMR re-fluxing
 *                           operations
 * \param [in]   dt          the current time-step
 * \param [in]   beg,end     initial and final interface indices
 * \param [in]   grid        pointer to Grid structure.
 *********************************************************************** */
{
  int i, j, k, nv;
  double *x1  = grid->x[IDIR],   *x2 = grid->x[JDIR],   *x3  = grid->x[KDIR];
  double *dx1 = grid->dx[IDIR], *dx2 = grid->dx[JDIR], *dx3  = grid->dx[KDIR];
  double *x1p = grid->xr[IDIR], *x2p = grid->xr[JDIR],  *x3p = grid->xr[KDIR];
  double *x1m = grid->xl[IDIR], *x2m = grid->xl[JDIR],  *x3m = grid->xl[KDIR];
  double *s   = grid->s, *sp = grid->sp;
  double A, dtdV, wp, w;
  double rhs[NVAR];
  static double **ViF, **ViS, **fxA, **src;
  double vc[NVAR], vi[NVAR];    /* Center and interface values */

  /* -- Viscous heating correction: locally evaluated nu1 and the
        associated capped energy sink; see step-by-step comments
        below for the physical/numerical rationale. -- */
  double nu1, nu2, rr, tt, st, rr2, rr4, st2, st4;
  double nu_floor = 1e-25;      /* below this, nu1 is treated as zero  */
  double cost = 2;              /* O(1) prefactor for the viscous heat */
  double cool_cutoff = 1e-1;    /* max fraction of thermal energy      */
                                 /* removable in a single step         */
  double cool;
  intList var_list;
  #if HAVE_ENERGY
  var_list.nvar = 4;
  var_list.indx[i=0] = MX1;
  var_list.indx[++i] = MX2;
  var_list.indx[++i] = MX3;
  var_list.indx[++i] = ENG;
  #else
  var_list.nvar = 3;
  var_list.indx[i=0] = MX1;
  var_list.indx[++i] = MX2;
  var_list.indx[++i] = MX3;
  #endif

/* --------------------------------------------------------
   0. Initialize flux and src to zero
   --------------------------------------------------------- */

  if (ViF == NULL){
    ViF     = ARRAY_2D(NMAX_POINT, NVAR, double);
    ViS     = ARRAY_2D(NMAX_POINT, NVAR, double);
    fxA     = ARRAY_2D(NMAX_POINT, NVAR, double);
    src     = ARRAY_2D(NMAX_POINT, NVAR, double);
  }

  for (i = beg; i <= end; i++) NVAR_LOOP(nv) {
    fxA[i][nv] = src[i][nv] = 0.0;
    ViS[i][nv] = 0.0;  /* This resets ViS[ENG] = 0, since it's not
                        * always done in viscous flux */
  }

  i = g_i; j = g_j; k = g_k;
  ViscousFlux (d, ViF, ViS, dcoeff, beg-1, end, grid);

  #ifdef FARGO
  double **wA = FARGO_Velocity();
  #endif

  if (g_dir == IDIR){

  /* --------------------------------------------------------
     1. Compute fluxes & sources in the X1 direction
     -------------------------------------------------------- */

    for (i = beg-1; i <= end; i++){
      A = grid->A[IDIR][k][j][i];
      FOR_EACH(nv, &var_list){
        fxA[i][nv] = A*ViF[i][nv];
        src[i][nv] = ViS[i][nv];
      }

    /* -- 1a. Correct energy flux in rotating frame --

           In a frame rotating at g_OmegaZ, the work done by viscous
           stresses on the azimuthal (background) rotation must be
           added back to the energy flux: E_flux += Omega*R * M_phi_flux,
           consistent with the frame-transformation of the energy
           equation.
       -- */

      #if HAVE_ENERGY && (ROTATING_FRAME == YES)
      #if GEOMETRY == POLAR || GEOMETRY == CYLINDRICAL
      wp = g_OmegaZ*x1p[i];
      #elif GEOMETRY == SPHERICAL
      wp = g_OmegaZ*x1p[i]*s[j];
      #endif
      fxA[i][ENG] += wp*fxA[i][iMPHI];
      #endif

      #ifdef iMPHI
      fxA[i][iMPHI] *= fabs(x1p[i]);
      #endif
    }

  /* -- 1b. Build rhs in the X1-direction and update --

         Standard finite-volume flux differencing plus the local
         viscous source term (src, from ViscousFlux, e.g. curvature
         terms in the stress divergence not captured by simple flux
         differencing).
     -- */

    for (i = beg; i <= end; i++){
      dtdV = dt/grid->dV[k][j][i];
      FOR_EACH(nv, &var_list){
        rhs[nv] = dtdV*(fxA[i][nv] - fxA[i-1][nv]) + dt*src[i][nv];
      }

      /* -- 1c. Explicit viscous heating correction.

             Evaluate the local dynamical viscosity nu1 at the cell
             center, then estimate the local viscous heating rate
             from the momentum-flux components of the viscous flux
             (ViF[iMR], ViF[iMTH], ViF[iMPHI]) and remove it from
             the internal energy:

               cool = cost*tracer*dt * 0.5*(ViF_R^2+ViF_TH^2+ViF_PHI^2)/nu1

             This is the viscous-heating analog of the Ohmic cooling
             correction in ResistiveRHS(): applied only where nu1
             exceeds a numerical floor (nu_floor), and only if the
             resulting energy removal stays below a fixed fraction
             (cool_cutoff) of the local thermal energy p/(gamma-1),
             to guard against over-cooling / negative pressures from
             locally large viscous stresses.
         -- */

      NVAR_LOOP(nv) {
        vc[nv] = d->Vc[nv][k][j][i];
      }
      Visc_nu(vc, x1[i], x2[j], x3[k], &nu1, &nu2);
      if (nu1 > nu_floor) {
//      cool = cost*vc[TRC]*dt*0.5*(pow(ViF[i][iMR],2)+pow(ViF[i][iMTH],2)+pow(ViF[i][iMPHI],2))/(nu1);
  /* Replaced vc[TCR] with the DiskFraction(vc, x1, x2) reconstruction. */
        cool = cost*DiskFraction(vc, x1[i], x2[j])*dt*0.5*(pow(ViF[i][iMR],2)+pow(ViF[i][iMTH],2)
                                     +pow(ViF[i][iMPHI],2))/(nu1);
        if (fabs(cool) < cool_cutoff*vc[PRS]/(g_gamma-1)) {
          rhs[ENG] -= cool;
        }
      }

      #ifdef iMPHI
      rhs[iMPHI] /= fabs(x1[i]);
      #endif

    /* -- 1d. Correct energy rhs in rotating frame or fargo --

           Subtracts the work done against the (background/FARGO)
           azimuthal velocity w so that only the residual velocity's
           kinetic energy is evolved explicitly - standard practice
           when using a rotating frame or the FARGO orbital-advection
           algorithm.
       -- */

      #if HAVE_ENERGY
      w = 0.0;

      #ifdef FARGO
      #if GEOMETRY == SPHERICAL
      w = wA[j][i];
      #else
      w = wA[k][i];
      #endif
      #endif

      #if ROTATING_FRAME == YES
      #if GEOMETRY == POLAR || GEOMETRY == CYLINDRICAL
      #elif GEOMETRY == SPHERICAL
      w += g_OmegaZ*x1[i]*s[j];
      #endif
      #endif /* ROTATING_FRAME == YES */

      #if GEOMETRY == CARTESIAN
      rhs[ENG] -= w*rhs[MX2];
      #else
      rhs[ENG] -= w*rhs[iMPHI];
      #endif
      #endif /* HAVE_ENERGY */

    /* -- 1e. Update -- */

      FOR_EACH(nv, &var_list) dU[k][j][i][nv] += rhs[nv];
    }

  }else if (g_dir == JDIR){

  /* --------------------------------------------------------
     2. Compute fluxes & sources in the X2 direction
     -------------------------------------------------------- */

    for (j = beg-1; j <= end; j++){
      A = grid->A[JDIR][k][j][i];
      FOR_EACH(nv, &var_list){
        fxA[j][nv] = A*ViF[j][nv];
        src[j][nv] = ViS[j][nv];
      }

    /* -- 2a. Correct energy flux in rotating frame (spherical only) -- */

      #if (GEOMETRY == SPHERICAL) && HAVE_ENERGY && (ROTATING_FRAME == YES)
      wp = g_OmegaZ*x1[i]*sp[j];
      fxA[j][ENG] += wp*fxA[j][iMPHI];
      #endif

      #if (GEOMETRY == SPHERICAL)
      fxA[j][iMPHI] *= fabs(sp[j]);
      #endif
    }

  /* -- 2b. Build rhs in the X2-direction and update -- */

    for (j = beg; j <= end; j++){
      dtdV = dt/grid->dV[k][j][i];
      FOR_EACH(nv, &var_list){
        rhs[nv] = dtdV*(fxA[j][nv] - fxA[j-1][nv]) + dt*src[j][nv];
      }

      /* -- 2c. Viscous heating correction (X2 momentum-flux component) --
             Same construction as step 1c, using the local nu1 and the
             viscous momentum fluxes at the X2 interface index j.
         -- */

      NVAR_LOOP(nv) {
        vc[nv] = d->Vc[nv][k][j][i];
      }
      Visc_nu(vc, x1[i], x2[j], x3[k], &nu1, &nu2);
      if (nu1 > nu_floor) {
//      cool = cost*vc[TRC]*dt*0.5*(pow(ViF[j][iMR],2)+pow(ViF[j][iMTH],2)+pow(ViF[j][iMPHI],2))/(nu1);
  /* Replaced vc[TCR] with the DiskFraction(vc, x1, x2) reconstruction. */
        cool = cost*DiskFraction(vc, x1[i], x2[j])*dt*0.5*(pow(ViF[j][iMR],2)+pow(ViF[j][iMTH],2)
                                     +pow(ViF[j][iMPHI],2))/(nu1);
        if (fabs(cool) < cool_cutoff*vc[PRS]/(g_gamma-1)) {
          rhs[ENG] -= cool;
        }
      }

      #if (GEOMETRY == SPHERICAL)
      rhs[iMPHI] /= fabs(s[j]);
      #endif

    /* -- 2d. Correct energy rhs in rotating frame or fargo -- */

      #if HAVE_ENERGY
      w = 0.0;

      #ifdef FARGO
      #if GEOMETRY == SPHERICAL
      w = wA[j][i];
      #else
      w = wA[k][i];
      #endif
      #endif

      #if ROTATING_FRAME == YES
      #if GEOMETRY == SPHERICAL
      w += g_OmegaZ*x1[i]*s[j];
      #endif
      #endif /* ROTATING_FRAME == YES */

      #if GEOMETRY == CARTESIAN
      rhs[ENG] -= w*rhs[MX2];
      #else
      rhs[ENG] -= w*rhs[iMPHI];
      #endif
      #endif /* HAVE_ENERGY */

    /* -- 2e. Update -- */

      FOR_EACH(nv, &var_list) dU[k][j][i][nv] += rhs[nv];
    }

  }else if (g_dir == KDIR){

  /* --------------------------------------------------------
     3a. Compute fluxes & sources in the X3 direction
     -------------------------------------------------------- */

    for (k = beg-1; k <= end; k++){
      A = grid->A[KDIR][k][j][i];
      FOR_EACH(nv, &var_list){
        fxA[k][nv] = A*ViF[k][nv];
        src[k][nv] = ViS[k][nv];
      }
    }

  /* -- 3b. Build rhs in the X3-direction -- */

    for (k = beg; k <= end; k++){
      dtdV = dt/grid->dV[k][j][i];

      FOR_EACH(nv, &var_list){
        rhs[nv] = dtdV*(fxA[k][nv] - fxA[k-1][nv]) + dt*src[k][nv];
      }

      /* -- 3c. Viscous heating correction (X3 momentum-flux component) --
             Same construction as steps 1c/2c, using the local nu1 and
             the viscous momentum fluxes at the X3 interface index k.
         -- */

      NVAR_LOOP(nv) {
        vc[nv] = d->Vc[nv][k][j][i];
      }
      Visc_nu(vc, x1[i], x2[j], x3[k], &nu1, &nu2);
      if (nu1 > nu_floor) {
//      cool = cost*vc[TRC]*dt*0.5*(pow(ViF[k][iMR],2)+pow(ViF[k][iMTH],2)+pow(ViF[k][iMPHI],2))/(nu1);
  /* Replaced vc[TCR] with the DiskFraction(vc, x1, x2) reconstruction. */
        cool = cost*DiskFraction(vc, x1[i], x2[j])*dt*0.5*(pow(ViF[k][iMR],2)+pow(ViF[k][iMTH],2)
                                     +pow(ViF[k][iMPHI],2))/(nu1);
        if (fabs(cool) < cool_cutoff*vc[PRS]/(g_gamma-1)) {
          rhs[ENG] -= cool;
        }
      }

      #if HAVE_ENERGY
      w = 0.0;

      #ifdef FARGO
      #if GEOMETRY == SPHERICAL
      w = wA[j][i];
      #else
      w = wA[k][i];
      #endif
      #endif

      #if GEOMETRY == CARTESIAN
      rhs[ENG] -= w*rhs[MX2];
      #else
      rhs[ENG] -= w*rhs[iMPHI];
      #endif

      #endif /* HAVE_ENERGY */

      FOR_EACH(nv, &var_list)  dU[k][j][i][nv] += rhs[nv];
    }
  }

/* --------------------------------------------------------
   4. Store AMR fluxes
   -------------------------------------------------------- */

  #ifdef CHOMBO
  StoreAMRFlux (ViF, aflux, -1, MX1, MX3, beg-1, end, grid);
  #if HAVE_ENERGY
  StoreAMRFlux (ViF, aflux, -1, ENG, ENG, beg-1, end, grid);
  #endif
  #endif
}
