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
 * Adaptive disk versus corona classifier state
 *
 * g_rBinEdgesLog[]     : logarithmic edges of NBINS_PROFILE bins spanning the
 *                        global radial domain (fixed at InitDomain() call).
 * g_rhoCoronaProfile[] : temporally smoothed, multicore reduced average density
 *                        of corona weighted cells per radial bin.
 * g_rhoDiskProfile[]   : same, but for disk weighted cells - a running
 *                        estimate of "what disk density actually looks like
 *                        at radius r", tracked alongside the corona
 *                        reference so DiskFraction() can classify by where
 *                        the cell sits between the two evolving references,
 *                        instead of by a fixed multiplicative threshold on
 *                        the corona value alone. This lets the classifier
 *                        track a genuinely decaying inner disk instead of
 *                        drifting away from it.
 * g_profilesInit  : guards against use before InitProfiles()
 *                        has run (e.g. a restart path that skips
 *                        InitDomain() call); in this case DiskFraction()
 *                        falls back to the original static analytic profile
 *                        for the corona side, and the analytic Keplerian
 *                        disk profile for the disk side.
 * --------------------------------------------------------------------- */
static double g_rBinEdgesLog[NBINS_PROFILE + 1];
static double g_rhoCoronaProfile[NBINS_PROFILE];
static double g_rhoDiskProfile[NBINS_PROFILE];
static int    g_diskProfileValid[NBINS_PROFILE];  /* has bin b ever seen a
                                                      meaningful amount of
                                                      real disk material? */
static int    g_profilesInit = 0;
static int    g_profilesLive = 0;  /* has UpdateProfiles() run at least once
                                      on real simulation data? Until then,
                                      DiskFraction() returns the exact t=0
                                      classification (see InitDiskCell()),
                                      matching what the removed tr1 tracer
                                      would have been, instead of using the
                                      profile-based sigmoid. */
 
/* ********************************************************************* */
static double InterpLogProfile (double *profile, double x1)
/*!
 * Shared linear interpolation (in log density, log radius space) of a
 * tabulated radial profile at input radius. Used for both the corona and
 * disk reference profiles. Falls back to clamped edge bins outside the
 * tabulated range (e.g. within ghost zones just past the physical domain
 * edge).
 *********************************************************************** */
{
  double lr, lmin, lmax, dl, lc0, s, frac;
  double lrho0, lrho1, lrho;
  int    ib0, ib1;
 
  lr   = log10(x1);
  lmin = g_rBinEdgesLog[0];
  lmax = g_rBinEdgesLog[NBINS_PROFILE];
  dl   = (lmax - lmin) / (double)NBINS_PROFILE;
  lc0  = lmin + 0.5 * dl;                 /* center of first bin */
 
  s    = (lr - lc0) / dl;                 /* fractional bin coordinate */
  ib0  = (int)floor(s);
  frac = s - ib0;
  ib1  = ib0 + 1;
 
  ib0 = (ib0 < 0) ? 0 : (ib0 > NBINS_PROFILE - 1 ? NBINS_PROFILE - 1 : ib0);
  ib1 = (ib1 < 0) ? 0 : (ib1 > NBINS_PROFILE - 1 ? NBINS_PROFILE - 1 : ib1);
 
  lrho0 = log10(MAX(profile[ib0], 1.e-30));
  lrho1 = log10(MAX(profile[ib1], 1.e-30));
  lrho  = lrho0 + frac * (lrho1 - lrho0);
 
  return pow(10.0, lrho);
}
 
/* ********************************************************************* */
static double GetCoronaRefDensity (double x1)
/*!
 * Interpolated corona reference density at input radius.
 *********************************************************************** */
{
  return InterpLogProfile(g_rhoCoronaProfile, x1);
}
 
/* ********************************************************************* */
static double GetDiskRefDensity (double x1)
/*!
 * Interpolated disk reference density at input radius.
 *********************************************************************** */
{
  return InterpLogProfile(g_rhoDiskProfile, x1);
}
 
/* ********************************************************************* */
static double AnalyticCoronaDensity (double r)
/*!
 * Static analytic corona density profile (RHOC*r^-1.5), used to seed the
 * corona reference profile and as a restart-safety fallback.
 *********************************************************************** */
{
  return g_inputParam[RHOC] * pow(r, -1.5);
}
 
/* ********************************************************************* */
static double AnalyticDiskDensity (double r)
/*!
 * Midplane (rcyl = r) analytic Keplerian-disk density, same construction
 * as the torus profile in Init(). Used to seed the disk reference profile
 * and as a restart-safety fallback.
 *********************************************************************** */
{
  double eps2, coeff;
 
  eps2  = g_inputParam[EPS] * g_inputParam[EPS];
  coeff = 0.4 / eps2 * (1.0 / r - (1.0 - 2.5 * eps2) / r);
  coeff = MAX(coeff, 0.0);
 
  return pow(coeff, 1.5);
}
 
/* ********************************************************************* */
static void PatchUnvalidatedDiskBins (void)
/*!
 * Fill any bin that has never seen a meaningful amount of real disk
 * material (g_diskProfileValid[b] == 0) by copying the nearest already
 * validated bin's current g_rhoDiskProfile value. Searches outward
 * (increasing r) first, since real disk material first appears from
 * larger radii and spreads/settles inward (e.g. after ISCO truncation
 * at t=0), falling back to an inward search if nothing outward is valid
 * either. Leaves the bin's placeholder value untouched only if no bin
 * anywhere is validated yet (before any disk material exists at all).
 *********************************************************************** */
{
  int b, bb, src;
 
  for (b = 0; b < NBINS_PROFILE; b++) {
    if (g_diskProfileValid[b]) continue;
 
    src = -1;
    for (bb = b + 1; bb < NBINS_PROFILE; bb++) {
      if (g_diskProfileValid[bb]) { src = bb; break; }
    }
    if (src < 0) {
      for (bb = b - 1; bb >= 0; bb--) {
        if (g_diskProfileValid[bb]) { src = bb; break; }
      }
    }
    if (src >= 0) g_rhoDiskProfile[b] = g_rhoDiskProfile[src];
  }
}
 
/* ********************************************************************* */
void InitProfiles (Grid *grid)
/*!
 * One time setup of the log radial bin edges spanning the global
 * radial domain. Seeds both reference profiles with their respective
 * original analytic values (corona: RHOC*r^-1.5; disk: the initial
 * Keplerian-disk density coeff^1.5 from Init(), evaluated on the
 * midplane where rcyl = r) so that DiskFraction() behaves sensibly
 * before the first multicore reduced averages are available.
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
 
  for (b = 0; b <= NBINS_PROFILE; b++) {
    g_rBinEdgesLog[b] = lmin + (lmax - lmin) * (double)b / (double)NBINS_PROFILE;
  }
 
  for (b = 0; b < NBINS_PROFILE; b++) {
    lc = 0.5 * (g_rBinEdgesLog[b] + g_rBinEdgesLog[b + 1]);
    r  = pow(10.0, lc);
 
    g_rhoCoronaProfile[b] = AnalyticCoronaDensity(r);
 
    /* The analytic torus solution is only meaningful in the actual disk
       region (rcyl > RD, matching the truncation in Init()); inside
       that radius (e.g. within ISCO for a black hole run) there is no
       disk at t=0 by construction, so don't seed a fabricated "disk"
       density there - mark the bin unvalidated instead and let
       UpdateProfiles()'s nearest-valid-bin fallback fill it in
       once real disk material actually exists at some radius. */
    if (r > g_inputParam[RD]) {
      g_rhoDiskProfile[b]   = AnalyticDiskDensity(r);
      if (g_rhoDiskProfile[b] <= 0.0) g_rhoDiskProfile[b] = g_rhoCoronaProfile[b];
      g_diskProfileValid[b] = 1;
    } else {
      g_rhoDiskProfile[b]   = g_rhoCoronaProfile[b];  /* placeholder only */
      g_diskProfileValid[b] = 0;
    }
  }
 
  g_profilesInit = 1;
  PatchUnvalidatedDiskBins();   /* so the very first DiskFraction() call,
                                    before any UpdateProfiles(), does
                                    not see an unphysical inside-RD seed */
}
 
/* ********************************************************************* */
void UpdateProfiles (const Data *d, Grid *grid)
/*!
 * Recompute the radially-binned corona AND disk density profiles:
 *
 *   1. Soft classify each active cell using the current (previous
 *      update) profile via DiskFraction(); f = diskFraction.
 *   2. Accumulate (1-f)*rho*dV / (1-f)*dV into the corona bin, and
 *      f*rho*dV / f*dV into the disk bin, for the cell's radial bin.
 *   3. Sum across all ranks -> every rank ends up with the same global
 *      profiles regardless of domain decomposition in X1/X2/X3.
 *   4. Exponential moving average in time to avoid abrupt jumps, applied
 *      independently to each profile.
 *
 * Bins with zero accumulated volume this step (no weighted cells found
 * for that profile) retain their previous value. Because both profiles
 * are built from the same (previous-step) DiskFraction() classification,
 * they stay mutually consistent: the disk reference tracks the actual
 * evolving disk (including a decaying inner disk) instead of remaining
 * pinned to the static initial analytic profile.
 *********************************************************************** */
{
  static double sumC_loc[NBINS_PROFILE], volC_loc[NBINS_PROFILE];
  static double sumC_glob[NBINS_PROFILE], volC_glob[NBINS_PROFILE];
  static double sumD_loc[NBINS_PROFILE], volD_loc[NBINS_PROFILE];
  static double sumD_glob[NBINS_PROFILE], volD_glob[NBINS_PROFILE];
  static double volTot_loc[NBINS_PROFILE], volTot_glob[NBINS_PROFILE];
 
  int    i, j, k, b;
  double *x1 = grid->x[IDIR];
  double r, dV, f, wc, wd, vi[NVAR];
 
  if (!g_profilesInit) InitProfiles(grid);   /* restart safety net */
 
  for (b = 0; b < NBINS_PROFILE; b++) {
    sumC_loc[b] = volC_loc[b] = 0.0;
    sumD_loc[b] = volD_loc[b] = 0.0;
    volTot_loc[b] = 0.0;
  }
 
  DOM_LOOP(k, j, i) {
    r = x1[i];
 
    b = (int)((log10(r) - g_rBinEdgesLog[0])
              / (g_rBinEdgesLog[NBINS_PROFILE] - g_rBinEdgesLog[0])
              * NBINS_PROFILE);
    if (b < 0) b = 0;
    if (b >= NBINS_PROFILE) b = NBINS_PROFILE - 1;
 
    vi[RHO] = d->Vc[RHO][k][j][i];              
    vi[VX3] = d->Vc[VX3][k][j][i];               /* Added: VX3 (v_phi) now read by DiskFraction */
    
    f  = DiskFraction(vi, r, grid->x[JDIR][j]);  /* soft disk weight, uses OLD profiles */
    wc = 1.0 - f;
    wd = f;
    dV = grid->dV[k][j][i];
 
    sumC_loc[b] += wc * vi[RHO] * dV;
    volC_loc[b] += wc * dV;
 
    sumD_loc[b] += wd * vi[RHO] * dV;
    volD_loc[b] += wd * dV;
 
    volTot_loc[b] += dV;
  }
 
#ifdef PARALLEL
  MPI_Allreduce(sumC_loc, sumC_glob, NBINS_PROFILE, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(volC_loc, volC_glob, NBINS_PROFILE, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(sumD_loc, sumD_glob, NBINS_PROFILE, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(volD_loc, volD_glob, NBINS_PROFILE, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(volTot_loc, volTot_glob, NBINS_PROFILE, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  for (b = 0; b < NBINS_PROFILE; b++) {
    sumC_glob[b] = sumC_loc[b]; volC_glob[b] = volC_loc[b];
    sumD_glob[b] = sumD_loc[b]; volD_glob[b] = volD_loc[b];
    volTot_glob[b] = volTot_loc[b];
  }
#endif
 
  /* -- Corona profile: unchanged, still just holds its previous value
        when unpopulated this step (corona material is present nearly
        everywhere so this is not the problematic case). -- */
  for (b = 0; b < NBINS_PROFILE; b++) {
    if (volC_glob[b] > 0.0) {
      double new_val = sumC_glob[b] / volC_glob[b];
      g_rhoCoronaProfile[b] = CORONA_EMA_ALPHA * new_val
                             + (1.0 - CORONA_EMA_ALPHA) * g_rhoCoronaProfile[b];
    }
  }
 
  /* -- Disk profile: only accept this step's estimate as "real" if the
        disk-weighted volume is a non-negligible fraction of the bin's
        total volume (guards against a single transient infalling cell
        validating a bin that is otherwise all corona/vacuum, e.g. just
        inside ISCO at early times). Once a bin is validated it keeps
        being updated normally (EMA) even if a later step sees f~0
        there again (retains its last known value, exactly as before).
     -- */
  for (b = 0; b < NBINS_PROFILE; b++) {
    int meaningful = (volTot_glob[b] > 0.0)
                   && (volD_glob[b] > DISK_BIN_VOLFRAC_MIN * volTot_glob[b]);
 
    if (meaningful) {
      double new_val = sumD_glob[b] / volD_glob[b];
      g_rhoDiskProfile[b] = CORONA_EMA_ALPHA * new_val
                           + (1.0 - CORONA_EMA_ALPHA) * g_rhoDiskProfile[b];
      g_diskProfileValid[b] = 1;
    }
  }
 
  /* -- Patch never-yet-validated bins (e.g. inside ISCO before any
        material has fallen in and settled into a disk-like state) by
        copying the nearest already-valid bin's current value - see
        PatchUnvalidatedDiskBins(). -- */
  PatchUnvalidatedDiskBins();
 
  g_profilesLive = 1;   /* from here on DiskFraction() uses the adaptive
                            sigmoid; this call has classified at least one
                            full DOM_LOOP pass through real cell data. */
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
static int InitDiskCell (double x1, double x2)
/*!
 * Reproduces, bit-for-bit, the branch condition Init() uses to decide
 * disk vs corona at problem setup (the same test that used to set
 * v[TRC] = 1.0 vs 0.0). Kept as a single source of truth so DiskFraction()
 * can match the former tracer exactly at t=0 without ever reading TRC
 * (which is slated for removal).
 *
 * Returns 1 if (x1, x2) is disk by the Init() criterion, 0 otherwise.
 *********************************************************************** */
{
  double rcyl, eps2, coeff, pc, p_disk;
 
  rcyl  = x1 * sin(x2);
  eps2  = g_inputParam[EPS] * g_inputParam[EPS];
  coeff = 0.4 / eps2 * (1.0 / x1 - (1.0 - 2.5 * eps2) / rcyl);
 
  pc     = 0.4 * g_inputParam[RHOC] * pow(x1, -2.5);   /* corona pressure, Init() step 1 */
  p_disk = eps2 * pow(coeff, 2.5);                     /* disk pressure, Init() step 2 */
 
  return (p_disk >= pc && rcyl > g_inputParam[RD]) ? 1 : 0;
}
 
/* ********************************************************************* */
double DiskFraction (double *v, double x1, double x2)
/*!
 * Continuous [0,1] replacement for the true tracer used to gate anomalous
 * viscosity/resistivity to disk material.
 *
 * At t=0 (before UpdateProfiles() has ever run on real simulation data),
 * this returns exactly 0.0 or 1.0 according to InitDiskCell(), i.e. it
 * matches what the removed tr1 tracer would have been initialized to,
 * without ever reading it.
 *
 * Once UpdateProfiles() has accumulated at least one step of real
 * disk/corona-weighted data (g_profilesLive == 1), DiskFraction() uses
 * two multiplicative criteria to decide if a cell is disk material:
 * 
 * 1. DENSITY: Tracks TWO running radial references - a corona density 
 *    profile and a disk density profile (see UpdateProfiles()) - and
 *    classifies by where the cell sits between them.
 * 
 * 2. ROTATION (New): Compares the local azimuthal velocity v[VX3] against
 *    the expected Keplerian velocity at the cylindrical radius.
 *
 * Both criteria use independent logistic functions (sigmoids), and their 
 * outputs are multiplied to yield the final continuous fraction [0, 1].
 *********************************************************************** */
{
  double rho_ref_c, rho_ref_d, log_span, position, arg_den, arg_rot;
  double den_factor, rot_factor, rcyl, v_kep, rot_frac;
 
  if (!g_profilesLive) {
    /* Exact t=0 match to the former tracer: no profile, no sigmoid,
       just the same hard boolean Init() used. Covers both the
       very-first-call case (g_profilesInit == 0) and the case right
       after InitProfiles() has seeded the analytic profiles but no
       real cell has been classified through UpdateProfiles() yet. */
    return (double) InitDiskCell(x1, x2);
  }
 
  if (!g_profilesInit) {
    /* Should not normally happen (g_profilesLive implies g_profilesInit),
       but guard anyway: fall back to the original static analytic
       corona/disk densities so behavior degrades gracefully. */
    rho_ref_c = AnalyticCoronaDensity(x1);
    rho_ref_d = AnalyticDiskDensity(x1);
    if (rho_ref_d <= 0.0) rho_ref_d = rho_ref_c;
  } else {
    rho_ref_c = GetCoronaRefDensity(x1);
    rho_ref_d = GetDiskRefDensity(x1);
  }
 
  /* --- 1. DENSITY SIGMOID --- */
  log_span = log10(MAX(rho_ref_d, 1.e-30)) - log10(MAX(rho_ref_c, 1.e-30));
  log_span = (fabs(log_span) < 1.e-12) ? 1.e-12 : log_span;  /* guard degenerate span */
 
  position = (log10(MAX(v[RHO], 1.e-30)) - log10(MAX(rho_ref_c, 1.e-30))) / log_span;
 
  arg_den = (position - (1.0 - 1.0/CORONA_THRESH_FAC)) / CORONA_SIGMOID_WIDTH;
  arg_den = MIN(MAX(arg_den, -50.0), 50.0);      /* guard exp() over/underflow */
  den_factor = 1.0 / (1.0 + exp(-arg_den));

  /* --- 2. ROTATION SIGMOID --- */
  rcyl  = x1 * sin(x2);
  /* Ideal Newtonian Keplerian velocity matching the 1/sqrt(rcyl) injection setup */
  v_kep = 1.0 / sqrt(MAX(rcyl, 1.e-12)); 
  
  /* Measure absolute rotation fraction against Keplerian */
  rot_frac = fabs(v[VX3]) / v_kep;
  
  arg_rot = (rot_frac - ROTATION_THRESH_FAC) / ROTATION_SIGMOID_WIDTH;
  arg_rot = MIN(MAX(arg_rot, -50.0), 50.0);
  rot_factor = 1.0 / (1.0 + exp(-arg_rot));

  /* Multiply factors: only rapidly rotating, dense material yields f ~ 1 */
  return den_factor * rot_factor;
}
 
/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
{
  InitProfiles(grid);
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*!
 * Called by PLUTO at the cadence set by the "analysis" entry in
 * pluto.ini. Used here to periodically refresh the corona/disk
 * radial reference profiles consumed by DiskFraction().
 *********************************************************************** */
{
  UpdateProfiles(d, grid);
}

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
