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

  FURTHER MODIFIED (opacity/DiskFraction consistency pass):
  - UserDefOpacities() (fixed signature, called from the radiation module's
    implicit step in rad_step.c) now gates opacities using the same
    DiskFraction() classifier used elsewhere for viscosity/resistivity,
    instead of an independent hard-cut condition on v[]. Since
    UserDefOpacities()'s signature carries no grid position, this relies on
    two new globals: g_i_rad (radial index, written by the patched local
    copy of rad_step.c immediately before each call) and g_radGrid (grid
    pointer, set once in InitDomain()). PLUTO's own g_j global (see
    globals.h) supplies the poloidal index; it is already valid at the
    point rad_step.c calls UserDefOpacities().
  - A second entry point, UserDefOpacitiesAt(v, x1, x2, abs, scat), takes
    x1, x2 explicitly. userdef_output.c's diagnostic loop calls this
    directly with its own loop's real (x1[i], x2[j]) rather than going
    through the g_i_rad/g_j globals, which would otherwise be stale
    (left over from whatever cell the last radiation step visited, not
    the cell the diagnostic loop is currently on). UserDefOpacities()
    itself is now a thin wrapper around UserDefOpacitiesAt().
  - v[TRC] is no longer read anywhere in opacity or viscosity/resistivity
    gating; DiskFraction() is the sole disk/corona classifier used for
    physics. TRC continues to be set/evolved elsewhere exactly as before,
    for comparison/diagnostic purposes only.
*/
/* /////////////////////////////////////////////////////////////////////////// */

#include "pluto.h"
#include "modifications.h"

#ifdef PARALLEL
 #include <mpi.h>
#endif

/* ---------------------------------------------------------------------
 * Inner Boundary Condition Selector
 * ---------------------------------------------------------------------
 * select the desired inner boundary physics regime:
 *   BOUNDARY_STAR : rotating conductive stellar surface
 *   BOUNDARY_BH   : black hole event horizon absorbing diode boundary
 * --------------------------------------------------------------------- */
#define BOUNDARY_STAR 1
#define BOUNDARY_BH   2

#ifndef INNER_BOUNDARY
 #define INNER_BOUNDARY BOUNDARY_BH   // set to BOUNDARY_STAR or BOUNDARY_BH
#endif

/* ---------------------------------------------------------------------
 * Radiation Opacity Call Support
 * --------------------------------------------------------------------- */
int          g_i_rad   = 0;      // current radial grid index written by patched radiation step for opacities
static Grid *g_radGrid = NULL;   // grid pointer set once during domain initialization to get physical coordinates

/* ---------------------------------------------------------------------
 * Adaptive Disk Versus Corona Classifier State
 * ---------------------------------------------------------------------
 * tracks corona and disk density to account for evolving
 * structure in deciding on classification
 * --------------------------------------------------------------------- */
static double g_rBinEdgesLog[NBINS_PROFILE + 1];   // logarithmic edges of profile bins along the global radial domain
static double g_rhoCoronaProfile[NBINS_PROFILE];   // temporally smoothed, multicore reduced, corona weighted average density per bin
static double g_rhoDiskProfile[NBINS_PROFILE];     // temporally smoothed, multicode reduced, disk weighted average density per bin
static int    g_diskProfileValid[NBINS_PROFILE];   // check if meaningful amount of disk material has been in cell before extrapolating

static int    g_profilesInit = 0;   /* falls back to exact initial profile instead of computed sigmoid until */
static int    g_profilesLive = 0;   /* update has run at least once to match initial tracer behavior         */
 
/* ********************************************************************* */
static double InterpLogProfile (double *profile, double x1)
/*!
 * Shared linear interpolation in logarithmic density and radius space of
 * a tabulated radial profile at given radius. Used for both the corona
 * and disk reference profiles. Falls back to clamped edge bins outside
 * the tabulated range, within ghost zones just past the domain edge.
 *********************************************************************** */
{
  double lr, lmin, lmax, dl, lc0, s, frac;
  double lrho0, lrho1, lrho;
  int    ib0, ib1;
 
  lr   = log10(x1);
  lmin = g_rBinEdgesLog[0];
  lmax = g_rBinEdgesLog[NBINS_PROFILE];
  dl   = (lmax - lmin) / (double)NBINS_PROFILE;
  lc0  = lmin + 0.5 * dl;                 // center of first bin
 
  s    = (lr - lc0) / dl;                 // fractional bin coordinates
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
 * Static analytic corona density profile, used to seed the
 * corona reference profile and as a restart safety fallback.
 *********************************************************************** */
{
  return g_inputParam[RHOC] * pow(r, -1.5);
}
 
/* ********************************************************************* */
static double AnalyticDiskDensity (double r)
/*!
 * Midplane analytic disk density, same construction as the initial
 * torus profile. Used to seed the disk reference profile and as a
 * restart safety fallback.
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
 * material by copying the nearest already validated current density
 * value. Searches outward first, since real disk material appears
 * from larger radii and settles inward, falling back to an inward
 * search if nothing outward is valid either. Leaves the bin placeholder
 * value untouched only if no bin anywhere is validated yet.
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
 * One time setup of the logarithmic radial bin edges spanning the global
 * radial domain. Seeds both reference profiles with their respective
 * original analytic values so that DiskFraction behaves sensibly
 * before the first multicore reduced averages are available.
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
  lambda = 2.2 / (1.0 + 2.56 * g_inputParam[BETAV] * g_inputParam[BETAV]);

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
    v[VX1] = -g_inputParam[BETAV] / sin(x2) * eps2 * (10.0 - (32.0 / 3.0)
             * lambda * g_inputParam[BETAV] * g_inputParam[BETAV]
             - lambda * (5.0 - 1.0 / (eps2 * tan(x2) * tan(x2)))) / sqrt(rcyl);
    v[VX3] = (sqrt(1.0 - 2.5 * eps2) + (2.0 / 3.0) * eps2
             * g_inputParam[BETAV] * g_inputParam[BETAV]
             * lambda * (1.0 - 1.2 / (eps2 * tan(x2) * tan(x2)))) / sqrt(rcyl);
    v[TRC] = 1.0;     /* Disk tracer - diagnostic only, see file header */
  } else {
    v[PRS] = 0.4 * g_inputParam[RHOC] * pow(x1, -2.5);
    v[TRC] = 0.0;     /* Corona tracer - diagnostic only, see file header */
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

  #if RADIATION
  v[ENR] = Blackbody(GetTemperature(v[RHO],v[PRS]))*v[TRC] ;
  v[FR1] = 0.;
  v[FR2] = 0.;
  v[FR3] = 0.;
  #endif /* RADIATION */
}

/* ********************************************************************* */
static int InitDiskCell (double x1, double x2)
/*!
 * Reproduces, bit-for-bit, the branch condition Init() uses to decide
 * disk vs corona at problem setup (the same test that used to set
 * v[TRC] = 1.0 vs 0.0). Kept as a single source of truth so DiskFraction()
 * can match the former tracer exactly at t=0 without ever reading TRC
 * (which is diagnostic-only, see file header).
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
 * Continuous [0,1] disk/corona classifier used to gate anomalous
 * viscosity, resistivity, AND (see UserDefOpacities()/UserDefOpacitiesAt()
 * below) opacities to disk material.
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
  g_radGrid = grid;   /* [ADDED] needed by UserDefOpacities() below */
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

        /* Reset coronal tracer (diagnostic only, see file header) */
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
        lambda = 2.2 / (1.0 + 2.56 * g_inputParam[BETAV] * g_inputParam[BETAV]);

        if (x2[j] >= 0.5 * CONST_PI - atan(1.25 * g_inputParam[EPS]) &&
            x2[j] <= 0.5 * CONST_PI + atan(1.25 * g_inputParam[EPS])) {

          d->Vc[RHO][k][j][i] = pow(coeff, 1.0 / (g_gamma - 1.0));
          
          if (d->Vc[RHO][k][j][i] == 0.0) { 
            d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][IEND];
            d->Vc[PRS][k][j][i] = eps2 * pow(coeff, g_gamma / (g_gamma - 1.0));
          }  

          d->Vc[VX1][k][j][i] = -g_inputParam[BETAV] / sin(x2[j]) * eps2
            * (10.0 - (32.0 / 3.0) * lambda * g_inputParam[BETAV] * g_inputParam[BETAV]
            - lambda * (5.0 - 1.0 / (eps2 * tan(x2[j]) * tan(x2[j])))) / sqrt(rcyl);

          d->Vc[VX3][k][j][i] = (sqrt(1.0 - 2.5 * eps2) + (2.0 / 3.0) * eps2
            * g_inputParam[BETAV] * g_inputParam[BETAV]
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
 *
 * NOTE: this implements the Newtonian -1/r^2 force only. It is currently
 * unused while BODY_FORCE == POTENTIAL (BodyForcePotential() below is the
 * active potential and implements Paczynski-Wiita). If BODY_FORCE is ever
 * switched to VECTOR, this function would need to be updated to match
 * -d/dr[-1/(r-2)] rather than silently reverting to Newtonian gravity.
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
 *
 * Active option: Paczynski-Wiita, singularity at x1 = 2 (i.e. r = 2 r_g =
 * r_s, the Schwarzschild radius), consistent with UNIT_LENGTH = r_g.
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

#if RADIATION_VAR_OPACITIES
        //source:  Carroll (1996), pages 274–276
        #include <math.h>

        #define G_BF 1.0        // Bound-Free Gaunt factor
        #define G_FF 1.0        // Free-Free Gaunt factor
        #define T_FACTOR 1e1   // Bound-Free correction factor, typically 1 < T < 100

        #define C_BF 4.34e25 // Kramer's law's bound-free constant, CGS
        #define C_FF 3.68e22 // Kramer's law's free-free constant,  CGS

        #define X     0.8    // Hydrogen mass fraction
        #define Z    0.00    // Metallicity

    //scattering constants (CGS)
    const double K_BF = 1;//C_BF * G_BF * Z * (1.0  + X) / T_FACTOR;
    const double K_FF = 1;//C_FF * (1.0  - Z) * (1.0  + X);
    const double K_ES = 1;//0.2  * (1.0 + X);


/* ********************************************************************* */
void UserDefOpacitiesAt(double *v, double x1, double x2, double *abs, double *scat)
/*!
 * Core opacity evaluation, taking x1, x2 explicitly. Called from two
 * places:
 *   1. UserDefOpacities() below, which is the fixed-signature entry point
 *      the radiation module calls from rad_step.c; it supplies x1, x2 via
 *      the g_i_rad / g_j globals.
 *   2. userdef_output.c's diagnostic loop directly, with its own loop's
 *      real (x1[i], x2[j]) - NOT via the globals, which would be stale
 *      there (left over from whichever cell the last radiation implicit
 *      step visited, not the cell the diagnostic loop is currently on).
 *
 * Opacities are gated by DiskFraction(v, x1, x2), the same continuous
 * disk/corona classifier used for viscosity and resistivity, in place of
 * the previous independent hard-cut condition on v[] (and no longer
 * reads v[TRC], which is diagnostic-only - see file header).
 *********************************************************************** */
{
    double rho = v[RHO];
    double T   = GetTemperature(v[RHO], v[PRS]);
    double f;

    double kappa_es   = K_ES;                        /* cm^2/g; previously
                                    pre-multiplied by v[TRC] - f (below)
                                    now plays that role instead */
    double kappa_ffbf = (K_BF + K_FF) * rho * pow(T, -3.5) ;  // cm^2/g

    f = DiskFraction(v, x1, x2);

    /* NOTE: verify against your original K_ES calibration - the previous
       form was *scat = fact*rho*(K_ES*v[TRC]), a hard 0/1 gate; this is
       now f*rho*K_ES with f continuous in [0,1]. If K_ES's numerical
       value was tuned assuming the old hard-cut form, its effective
       normalization may need revisiting now that the gate is smooth. */
    *scat = f * rho * kappa_es;        // Thomson scattering only
    *abs  = f * rho * kappa_ffbf;      // free-free + bound-free true absorption
}

/* ********************************************************************* */
void UserDefOpacities(double *v, double *abs, double *scat)
/*!
 * Fixed-signature entry point required by radiation.h / called from the
 * radiation module's implicit step (rad_step.c). Supplies x1, x2 via
 * g_i_rad (set by the locally-patched rad_step.c immediately before this
 * call, inside RadImplicitNR()'s per-cell loop) and g_j (PLUTO's own
 * global, already valid at that point - see globals.h). See
 * UserDefOpacitiesAt() above for the actual opacity evaluation.
 *********************************************************************** */
{
    double x1 = g_radGrid->x[IDIR][g_i_rad];
    double x2 = g_radGrid->x[JDIR][g_j];

    UserDefOpacitiesAt(v, x1, x2, abs, scat);
}
#endif /* RADIATION_VAR_OPACITIES */
