/* /////////////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Problem initialization file for accretion disk-corona simulations.

  Provides problem initialization routines for PLUTO MHD, supporting both 
  rotating conductive stellar surface boundaries and black hole event horizon
  absorbing diode boundaries.

  \author A. Mignone (mignone@ph.unito.it)
  \date Sep 2012

  \modified M. Cemeljic (miki@camk.edu.pl)
  \date Jul 2020 / Modified 2024 (ccm, csk)

  Based on appendix of "Atlas" paper, Cemeljic, 2019, A&A, 624, A31.
*/
/* /////////////////////////////////////////////////////////////////////////// */

#include "pluto.h"
#include "modifications.h"

#ifdef PARALLEL
 #include <mpi.h>
#endif

/* ---------------------------------------------------------------------
 * Inner Boundary Condition Selector (X1_BEG)
 * ---------------------------------------------------------------------
 * Select the desired inner boundary physics regime:
 *   BOUNDARY_STAR : Rotating conductive stellar surface
 *   BOUNDARY_BH   : Black hole event horizon absorbing diode boundary
 * --------------------------------------------------------------------- */
#define BOUNDARY_STAR 1
#define BOUNDARY_BH   2

#ifndef INNER_BOUNDARY
 #define INNER_BOUNDARY BOUNDARY_BH  /* Set to BOUNDARY_STAR or BOUNDARY_BH */
#endif

/* ---------------------------------------------------------------------
 * Adaptive disk versu corona classifier state
 *
 * g_rBinEdgesLog[]     : logarithmic edges of NBINS_CORONA bins spanning the
 *                        global radial domain (fixed at InitDomain() call).
 * g_rhoCoronaProfile[] : temporally smoothed, multicore reduced average density
 *                        of corona weighted cells per radial bin.
 * g_coronaProfileInit  : guards against use before InitCoronaProfile()
 *                        has run (e.g. a restart path that skips
 *                        InitDomain() call); in this case DiskFraction()
 *                        falls back to the original static analytic profile.
 * --------------------------------------------------------------------- */
static double g_rBinEdgesLog[NBINS_CORONA + 1];
static double g_rhoCoronaProfile[NBINS_CORONA];
static int    g_coronaProfileInit = 0;

/* ********************************************************************* */
static double GetCoronaRefDensity (double x1)
/*!
 * Linear interpolation (in log density, log radius space) of the
 * tabulated corona reference profile at input radius. Falls back to
 * clamped edge bins outside the tabulated range (e.g. within ghost
 * zones just past the physical domain edge).
 *********************************************************************** */
{
  double lr, lmin, lmax, dl, lc0, s, frac;
  double lrho0, lrho1, lrho;
  int    ib0, ib1;

  lr   = log10(x1);
  lmin = g_rBinEdgesLog[0];
  lmax = g_rBinEdgesLog[NBINS_CORONA];
  dl   = (lmax - lmin) / (double)NBINS_CORONA;
  lc0  = lmin + 0.5 * dl;                 /* center of first bin */

  s    = (lr - lc0) / dl;                 /* fractional bin coordinate */
  ib0  = (int)floor(s);
  frac = s - ib0;
  ib1  = ib0 + 1;

  ib0 = (ib0 < 0) ? 0 : (ib0 > NBINS_CORONA - 1 ? NBINS_CORONA - 1 : ib0);
  ib1 = (ib1 < 0) ? 0 : (ib1 > NBINS_CORONA - 1 ? NBINS_CORONA - 1 : ib1);

  lrho0 = log10(MAX(g_rhoCoronaProfile[ib0], 1.e-30));
  lrho1 = log10(MAX(g_rhoCoronaProfile[ib1], 1.e-30));
  lrho  = lrho0 + frac * (lrho1 - lrho0);

  return pow(10.0, lrho);
}

/* ********************************************************************* */
void InitCoronaProfile (Grid *grid)
/*!
 * One time setup of the log radial bin edges spanning the global
 * radial domain, seeded with the original analytic corona profile
 * (RHOC*r^-1.5) so that DiskFraction() behaves sensibly before the
 * first multicore reduced average is available.
 *
 * NOTE: g_domBeg[]/g_domEnd[] are assumed to hold the global physical
 * domain boundaries (set by PLUTO's grid setup from pluto.ini). If
 * your PLUTO version exposes these under different names, substitute
 * grid->xl_glob[IDIR][0] / grid->xr_glob[IDIR][grid->np_int_glob[IDIR]-1].
 *********************************************************************** */
{
  int    b;
  double lmin, lmax, lc, r;

  lmin = log10(g_domBeg[IDIR]);
  lmax = log10(g_domEnd[IDIR]);

  for (b = 0; b <= NBINS_CORONA; b++) {
    g_rBinEdgesLog[b] = lmin + (lmax - lmin) * (double)b / (double)NBINS_CORONA;
  }

  for (b = 0; b < NBINS_CORONA; b++) {
    lc = 0.5 * (g_rBinEdgesLog[b] + g_rBinEdgesLog[b + 1]);
    r  = pow(10.0, lc);
    g_rhoCoronaProfile[b] = g_inputParam[RHOC] * pow(r, -1.5);
  }

  g_coronaProfileInit = 1;
}

/* ********************************************************************* */
void UpdateCoronaProfile (const Data *d, Grid *grid)
/*!
 * Recompute the radially-binned corona density profile:
 *
 *   1. Soft classify each active cell using the current (previous
 *      update) profile via DiskFraction(); weight = 1 - diskFraction.
 *   2. Accumulate weight*rho*dV and weight*dV into the cell's radial bin.
 *   3. Sum across all ranks -> every rank ends up with the same global
 *      profile regardless of domain decomposition in X1/X2/X3.
 *   4. Exponential moving average in time to avoid abrupt jumps.
 *
 * Bins with zero accumulated volume this step (no corona weighted
 * cells found) retain their previous value.
 *********************************************************************** */
{
  static double sum_loc[NBINS_CORONA], vol_loc[NBINS_CORONA];
  static double sum_glob[NBINS_CORONA], vol_glob[NBINS_CORONA];

  int    i, j, k, b;
  double *x1 = grid->x[IDIR];
  double r, dV, w, vi[NVAR];

  if (!g_coronaProfileInit) InitCoronaProfile(grid);   /* restart safety net */

  for (b = 0; b < NBINS_CORONA; b++) sum_loc[b] = vol_loc[b] = 0.0;

  DOM_LOOP(k, j, i) {
    r = x1[i];

    b = (int)((log10(r) - g_rBinEdgesLog[0])
              / (g_rBinEdgesLog[NBINS_CORONA] - g_rBinEdgesLog[0])
              * NBINS_CORONA);
    if (b < 0) b = 0;
    if (b >= NBINS_CORONA) b = NBINS_CORONA - 1;

    vi[RHO] = d->Vc[RHO][k][j][i];              /* only RHO is read by DiskFraction */
    w  = 1.0 - DiskFraction(vi, r, grid->x[JDIR][j]);   /* soft corona weight, uses OLD profile */
    dV = grid->dV[k][j][i];

    sum_loc[b] += w * vi[RHO] * dV;
    vol_loc[b] += w * dV;
  }

#ifdef PARALLEL
  MPI_Allreduce(sum_loc, sum_glob, NBINS_CORONA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(vol_loc, vol_glob, NBINS_CORONA, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  for (b = 0; b < NBINS_CORONA; b++) { sum_glob[b] = sum_loc[b]; vol_glob[b] = vol_loc[b]; }
#endif

  for (b = 0; b < NBINS_CORONA; b++) {
    if (vol_glob[b] > 0.0) {
      double new_val = sum_glob[b] / vol_glob[b];
      g_rhoCoronaProfile[b] = CORONA_EMA_ALPHA * new_val
                             + (1.0 - CORONA_EMA_ALPHA) * g_rhoCoronaProfile[b];
    }
    /* else: no corona-weighted cells in this bin this step -> keep previous value */
  }
}

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*!
 * Primary initialization of primitive variables across the domain.
 *********************************************************************** */
{
  double coeff, eps2, pc, rcyl;
  double lambda;

  rcyl = x1 * sin(x2);
  eps2 = g_inputParam[EPS] * g_inputParam[EPS];
  coeff = 0.4 / eps2 * (1.0 / x1 - (1.0 - 2.5 * eps2) / rcyl);
  lambda = 2.2 / (1.0 + 2.56 * g_inputParam[ALPHAV] * g_inputParam[ALPHAV]);

  /* -------------------------------------------------------------------
     1. Initial non-rotating adiabatic corona in hydrostatic equilibrium
     ------------------------------------------------------------------- */
  v[RHO] = g_inputParam[RHOC] * pow(x1, -1.5);
  v[PRS] = 0.4 * g_inputParam[RHOC] * pow(x1, -2.5);
  pc     = v[PRS];

  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;

  /* -------------------------------------------------------------------
     2. Keplerian adiabatic disk (Kluzniak & Kita 2000)
     ------------------------------------------------------------------- */
  v[PRS] = eps2 * pow(coeff, 2.5);

  if (v[PRS] >= pc && rcyl > g_inputParam[RD]) {
    v[RHO] = pow(coeff, 1.5);
    v[VX1] = -g_inputParam[ALPHAV] / sin(x2) * eps2 * (10.0 - (32.0 / 3.0)
             * lambda * g_inputParam[ALPHAV] * g_inputParam[ALPHAV]
             - lambda * (5.0 - 1.0 / (eps2 * tan(x2) * tan(x2)))) / sqrt(rcyl);
    v[VX3] = (sqrt(1.0 - 2.5 * eps2) + (2.0 / 3.0) * eps2
             * g_inputParam[ALPHAV] * g_inputParam[ALPHAV]
             * lambda * (1.0 - 1.2 / (eps2 * tan(x2) * tan(x2)))) / sqrt(rcyl);
    v[TRC] = 1.0;     /* Disk tracer */
  } else {
    v[PRS] = 0.4 * g_inputParam[RHOC] * pow(x1, -2.5);
    v[TRC] = 0.0;     /* Corona tracer */
  }

  /* -------------------------------------------------------------------
     3. Initial Magnetic Field Setup
     ------------------------------------------------------------------- */
#if PHYSICS == MHD
#if BACKGROUND_FIELD == YES
  v[BX1] = 0.0;
  v[BX2] = 0.0;
  v[BX3] = 0.0;

  v[AX1] = 0.0;  
  v[AX2] = 0.0;
  v[AX3] = 0.0;
#else
  /* Cell-centered dipole field default */
  v[BX1] = 2.0 * g_inputParam[MU] * cos(x2) / (x1 * x1 * x1);
  v[BX2] = g_inputParam[MU] * sin(x2) / (x1 * x1 * x1);
  v[BX3] = 0.0;

  /* Vector potential for Constrained Transport */
  #if MHD_FORMULATION == CONSTRAINED_TRANSPORT
   v[AX1] = 0.0;
   v[AX2] = 0.0;
   v[AX3] = g_inputParam[MU] * sin(x2) / (x1 * x1);
  #endif
#endif
#endif /* PHYSICS == MHD */
}

/* ********************************************************************* */
double DiskFraction (double *v, double x1, double x2)
/*!
 * Continuous [0,1] replacement for the true tracer used to gate anomalous
 * viscosity/resistivity to disk material.
 *
 * Compares the local density against a radially binned, domain reduced,
 * temporally smoothed running average of the density found in cells
 * identified as corona like (see UpdateCoronaProfile() function), rather
 * than against the fixed initial analytic profile. This removes the
 * horizon side disk truncation that occurs once the inner disk density
 * evolves below the static analytic corona reference.
 *
 * The hard threshold is replaced by a logistic function in
 * log(rho / (CORONA_THRESH_FAC * rho_corona(r))), giving a smooth
 * spatial transition (in lieu of tracer diffusion) of width
 * CORONA_SIGMOID_WIDTH dex.
 *********************************************************************** */
{
  double rho_ref, ratio, arg;

  (void) x2;   /* angle dependence is folded into the radial profile;
                  kept in the signature for call-site compatibility  */

  if (!g_coronaProfileInit) {
    /* Fallback: profile not built yet (very first call, or a restart
       path that bypasses InitDomain() call). Use the original static
       analytic corona density so behavior degrades gracefully. */
    rho_ref = g_inputParam[RHOC] * pow(x1, -1.5);
  } else {
    rho_ref = GetCoronaRefDensity(x1);
  }

  ratio = v[RHO] / (CORONA_THRESH_FAC * rho_ref);
  arg   = log10(MAX(ratio, 1.e-30)) / CORONA_SIGMOID_WIDTH;
  arg   = MIN(MAX(arg, -50.0), 50.0);      /* guard exp() over/underflow */

  return 1.0 / (1.0 + exp(-arg));
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
{
  InitCoronaProfile(grid);
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
{}

#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 * Static, curl-free background magnetic field options.
 *********************************************************************** */
{
  /* Option A: Black Hole power-law field (Zhu & Stone 2018, Mishra et al. 2019) 
     ccm / MC Aug 2019
  
  double mu   = g_inputParam[MU];
  double rmin = 1.0 * g_inputParam[RD];
  double mm   = -1.25; // Mishra: -1.25, Zhu & Stone: -2.25

  if (x1 <= rmin) {
    B0[0] = mu * cos(x2) * pow(rmin, mm) * (1.0 + sin(x2));
    B0[1] = -mu * sin(x2) * pow(rmin, mm);
  } else {
    B0[0] = mu * pow(x1 * sin(x2), mm) * cos(x2) * (1.0 + sin(x2));
    B0[1] = -mu * pow(x1 * sin(x2), mm) * sin(x2);
  }
  B0[2] = 0.0; 
  */

  /* Option B: Stellar Dipole Field */
  B0[0] = 2.0 * g_inputParam[MU] * cos(x2) / (x1 * x1 * x1);
  B0[1] = g_inputParam[MU] * sin(x2) / (x1 * x1 * x1);
  B0[2] = 0.0;                             

  /* Option C: Quadrupole Field 
  B0[0] = 1.5 * g_inputParam[MU] * (3.0 * cos(x2) * cos(x2) - 1.0) / (x1 * x1 * x1 * x1);
  B0[1] = 3.0 * g_inputParam[MU] * cos(x2) * sin(x2) / (x1 * x1 * x1 * x1);
  B0[2] = 0.0;
  */       

  /* Option D: Octupole Field 
  B0[0] = 2.0 * g_inputParam[MU] * (5.0 * cos(x2) * cos(x2) * cos(x2) - 3.0 * cos(x2)) / (x1 * x1 * x1 * x1 * x1);
  B0[1] = 0.5 * g_inputParam[MU] * (15.0 * cos(x2) * cos(x2) * sin(x2) - 3.0 * sin(x2)) / (x1 * x1 * x1 * x1 * x1);
  B0[2] = 0.0;
  */         
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid)
/*!
 * User-defined boundary condition dispatcher.
 *********************************************************************** */
{
  int i, j, k;
  double *x1, *x2, *x3, *r;
  double a1, a2, a, rcyl, eps2, coeff, lambda;
  double dvar1dr, dvar2dr, dvardr;
  double cs2, dden, dfact;
  
  RBox dom_box;

  /* -----------------------------------------------------------------
     side == 0: Active Domain Internal Density Floor & Conservative Update
     ----------------------------------------------------------------- */
  if (side == 0) {    
    x1 = grid->xgc[IDIR];
    x2 = grid->xgc[JDIR];
    x3 = grid->xgc[KDIR];

    TOT_LOOP(k, j, i) {
      int convert_to_cons = 0;

      /* Domain density floor enforcement */
      if (d->Vc[RHO][k][j][i] < g_inputParam[DFLOOR]) {
        dden = d->Vc[RHO][k][j][i];
        cs2  = g_gamma * d->Vc[PRS][k][j][i] / d->Vc[RHO][k][j][i];

        d->Vc[RHO][k][j][i] = g_inputParam[DFLOOR];     
        dfact = dden / d->Vc[RHO][k][j][i];

        /* Preserve local speed of sound and momentum */
        d->Vc[PRS][k][j][i] = cs2 * d->Vc[RHO][k][j][i] / g_gamma;
        d->Vc[VX1][k][j][i] *= dfact;
        d->Vc[VX2][k][j][i] *= dfact;
        d->Vc[VX3][k][j][i] *= dfact;

        /* Reset coronal tracer */
        if (x2[j] < 0.5 * CONST_PI - atan(3.0 * g_inputParam[EPS]) ||
            x2[j] > 0.5 * CONST_PI + atan(3.0 * g_inputParam[EPS])) {
          d->Vc[TRC][k][j][i] = 0.0;
        }

        convert_to_cons = 1;
      }

      /* Recompute conservative variables if primitives changed */
      if (convert_to_cons) {
        RBoxDefine(i, i, j, j, k, k, CENTER, &dom_box);
        PrimToCons3D(d->Vc, d->Uc, &dom_box, grid);
      }
    }
  }

  /* -----------------------------------------------------------------
     side == X1_BEG: Inner Radial Boundary (Toggle Selected)
     ----------------------------------------------------------------- */
  if (side == X1_BEG) {
    if (box->vpos == CENTER) {

#if INNER_BOUNDARY == BOUNDARY_BH
      /* =============================================================
         OPTION 1: BLACK HOLE EVENT HORIZON ABSORBING DIODE BOUNDARY
         ============================================================= */
      X1_BEG_LOOP(k, j, i) {
        /* Copy primitive variables from first computational cell IBEG */
        d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][IBEG];
        d->Vc[PRS][k][j][i] = d->Vc[PRS][k][j][IBEG];
        d->Vc[VX1][k][j][i] = d->Vc[VX1][k][j][IBEG];
        d->Vc[VX2][k][j][i] = d->Vc[VX2][k][j][IBEG];
        d->Vc[VX3][k][j][i] = d->Vc[VX3][k][j][IBEG];
        d->Vc[TRC][k][j][i] = d->Vc[TRC][k][j][IBEG];

#if PHYSICS == MHD
        d->Vc[BX1][k][j][i] = d->Vc[BX1][k][j][IBEG];
        d->Vc[BX2][k][j][i] = d->Vc[BX2][k][j][IBEG];
        d->Vc[BX3][k][j][i] = d->Vc[BX3][k][j][IBEG];
#endif

        /* Diode condition: allow inflow (v_r <= 0), block outflow (v_r > 0) */
        if (d->Vc[VX1][k][j][i] > 0.0) {
          d->Vc[VX1][k][j][i] = 0.0;
        }

        /* Ghost cell floor enforcement */
        if (d->Vc[RHO][k][j][i] < g_inputParam[DFLOOR]) {
          d->Vc[RHO][k][j][i] = g_inputParam[DFLOOR];
        }
        if (d->Vc[PRS][k][j][i] < g_inputParam[DFLOOR] * 1.0e-3) {
          d->Vc[PRS][k][j][i] = g_inputParam[DFLOOR] * 1.0e-3;
        }
      }

#elif INNER_BOUNDARY == BOUNDARY_STAR
      /* =============================================================
         OPTION 2: STELLAR SURFACE BOUNDARY (Rotating / Conductive)
         ============================================================= */
      x2 = grid->x[JDIR];

      X1_BEG_LOOP(k, j, i) {
        d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][IBEG];
        d->Vc[PRS][k][j][i] = d->Vc[PRS][k][j][IBEG];
        d->Vc[TRC][k][j][i] = 0.0;

        /* Fixed velocity field at stellar surface */
        d->Vc[VX1][k][j][i] = 0.0;
        d->Vc[VX2][k][j][i] = 0.0;
        d->Vc[VX3][k][j][i] = g_inputParam[OMEGA_STAR] * grid->x[IDIR][IBEG] * sin(x2[j]);

#if PHYSICS == MHD
        /* Poloidal field matched, toroidal field set to zero (E_phi = 0) */
        d->Vc[BX1][k][j][i] = d->Vc[BX1][k][j][IBEG];
        d->Vc[BX2][k][j][i] = d->Vc[BX2][k][j][IBEG];
        d->Vc[BX3][k][j][i] = 0.0;
#endif

        if (d->Vc[RHO][k][j][i] < g_inputParam[DFLOOR]) {
          d->Vc[RHO][k][j][i] = g_inputParam[DFLOOR];
        }
      }
#endif

    }
  }

  /* -----------------------------------------------------------------
     side == X1_END: Outer Radial Boundary (Disk Injection & Corona)
     ----------------------------------------------------------------- */
  if (side == X1_END) {
    r  = grid->x[IDIR];
    x1 = grid->x[IDIR];
    x2 = grid->x[JDIR];

    if (box->vpos == CENTER) {
      BOX_LOOP(box, k, j, i) {
        d->Vc[TRC][k][j][i] = d->Vc[TRC][k][j][IEND];

        /* Logarithmic extrapolation of density */
        a1 = log10(d->Vc[RHO][k][j][IEND]   / d->Vc[RHO][k][j][IEND-1]) / log10(r[IEND]   / r[IEND-1]);
        a2 = log10(d->Vc[RHO][k][j][IEND-1] / d->Vc[RHO][k][j][IEND-2]) / log10(r[IEND-1] / r[IEND-2]);
        a  = VANLEER_LIMITER(a1, a2);
        a  = MIN(a, 0.0);

        d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][i-1] * pow(r[i] / r[i-1], a); 
        d->Vc[PRS][k][j][i] = d->Vc[PRS][k][j][IEND] * pow(d->Vc[RHO][k][j][i] / d->Vc[RHO][k][j][IEND], g_gamma);

        /* Outflow condition for poloidal velocity components */
        d->Vc[VX1][k][j][i] = d->Vc[VX1][k][j][IEND]; 
        d->Vc[VX2][k][j][i] = d->Vc[VX2][k][j][IEND]; 

#if PHYSICS == MHD
        /* Van Leer extrapolation for toroidal magnetic field */
        dvar1dr = (d->Vc[BX3][k][j][IEND]   - d->Vc[BX3][k][j][IEND-1]) / (r[IEND]   - r[IEND-1]);
        dvar2dr = (d->Vc[BX3][k][j][IEND-1] - d->Vc[BX3][k][j][IEND-2]) / (r[IEND-1] - r[IEND-2]);
        dvardr  = VANLEER_LIMITER(dvar1dr, dvar2dr);
        d->Vc[BX3][k][j][i] = d->Vc[BX3][k][j][i-1] + dvardr * (r[i] - r[i-1]);
#endif

        /* MinMod extrapolation for toroidal velocity */
        dvar1dr = (d->Vc[VX3][k][j][IEND]   - d->Vc[VX3][k][j][IEND-1]) / (r[IEND]   - r[IEND-1]);
        dvar2dr = (d->Vc[VX3][k][j][IEND-1] - d->Vc[VX3][k][j][IEND-2]) / (r[IEND-1] - r[IEND-2]);
        dvardr  = MINMOD_LIMITER(dvar1dr, dvar2dr);
        d->Vc[VX3][k][j][i] = d->Vc[VX3][k][j][i-1] + dvardr * (r[i] - r[i-1]);

        /* Re-inject initial Kluźniak & Kita disk profile in equatorial region */
        rcyl   = x1[i] * sin(x2[j]);
        eps2   = g_inputParam[EPS] * g_inputParam[EPS];
        coeff  = (g_gamma - 1.0) / g_gamma / eps2 * (1.0 / x1[i] - (1.0 - eps2 * g_gamma / (g_gamma - 1.0)) / rcyl);
        coeff  = MAX(coeff, 0.0);
        lambda = 2.2 / (1.0 + 2.56 * g_inputParam[ALPHAV] * g_inputParam[ALPHAV]);

        if (x2[j] >= 0.5 * CONST_PI - atan(1.25 * g_inputParam[EPS]) &&
            x2[j] <= 0.5 * CONST_PI + atan(1.25 * g_inputParam[EPS])) {

          d->Vc[RHO][k][j][i] = pow(coeff, 1.0 / (g_gamma - 1.0));
          
          if (d->Vc[RHO][k][j][i] == 0.0) { 
            d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][IEND];
            d->Vc[PRS][k][j][i] = eps2 * pow(coeff, g_gamma / (g_gamma - 1.0));
          }  

          d->Vc[VX1][k][j][i] = -g_inputParam[ALPHAV] / sin(x2[j]) * eps2
            * (10.0 - (32.0 / 3.0) * lambda * g_inputParam[ALPHAV] * g_inputParam[ALPHAV]
            - lambda * (5.0 - 1.0 / (eps2 * tan(x2[j]) * tan(x2[j])))) / sqrt(rcyl);

          d->Vc[VX3][k][j][i] = (sqrt(1.0 - 2.5 * eps2) + (2.0 / 3.0) * eps2
            * g_inputParam[ALPHAV] * g_inputParam[ALPHAV]
            * lambda * (1.0 - 1.2 / (eps2 * tan(x2[j]) * tan(x2[j])))) / sqrt(rcyl);
        }

        /* Prevent coronal back-inflow */
        if (x2[j] <= 0.5 * CONST_PI - atan(3.0 * g_inputParam[EPS]) ||
            x2[j] >= 0.5 * CONST_PI + atan(3.0 * g_inputParam[EPS])) {
          if (d->Vc[VX1][k][j][i] < 0.0) {
            d->Vc[VX1][k][j][i] = 0.0;
            d->Vc[VX2][k][j][i] = 0.0;
          }
        }
      }
    }
  }
}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Radial acceleration vector for point-mass central potential.
 *********************************************************************** */
{
  g[IDIR] = -1.0 / (x1 * x1);
  g[JDIR] = 0.0;
  g[KDIR] = 0.0; 
}

/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Central gravitational potential options.
 * Select regime in definitions.h: #define BODY_FORCE POTENTIAL
 *********************************************************************** */
{
  /* ccm -- Select active potential: */

  /* Option 1: Newtonian potential */
  /* return -1.0 / x1; */

  /* Option 2: Paczyński-Wiita potential (ccm) */
  return -1.0 / (x1 - 2.0);

  /* Option 3: Kluźniak-Lee potential (csk / ccm) */
  /* return -(1.0 / 6.0) * (exp(6.0 / x1) - 1.0); */

  /* Option 4: Kluźniak-Nordström pseudo-potential */
  /* double q = 1.25; */
  /* return -1.0 / x1 + 0.5 * (q * q / (x1 * x1)); */
}
#endif
