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

  ========================================================================================
  ========================================================================================
  Further modified to include a replacement for the tracer and fixes to radiation module.
  Opacity code was slightly extended to expose position to DiskFrac call. Radiation runs
  in physical units: the module constants are set from the code units, the disk starts
  in local thermodynamic equilibrium and the radial boundaries handle radiation fields.
  Compton energy exchange between gas and radiation enters as an effective absorption
  opacity, the gas temperature is floored at 1e4 K, and matter and radiation decouple
  where the gas moves faster than the nonrelativistic radiation module can handle.
  ========================================================================================
  ========================================================================================
*/
/* /////////////////////////////////////////////////////////////////////////// */

#include "pluto.h"
#include "modifications.h"

#ifdef PARALLEL
 #include <mpi.h>
#endif

/* ---------------------------------------------------------------------
 * Physical Units
 * ---------------------------------------------------------------------
 * Code units are fixed in definitions.h through BH_MASS:
 *
 *   UNIT_LENGTH   = R_g = G M / c^2
 *   UNIT_VELOCITY = c   = sqrt(G M / UNIT_LENGTH)
 *   UNIT_DENSITY  = free normalization, sets the disk density in g/cm^3
 *
 * so G M = 1 in code units, which is exactly what the potential, the
 * corona and the Kluzniak & Kita disk below assume. The Keplerian speed at
 * r = 1 (the "stellar surface" of the original setup, here r = R_g) is
 * therefore identical to the speed of light, and all velocities written
 * by the initialization are consistently in units of c.
 *
 * The radiation module does not know about any of this by itself, all of
 * its constants default to unity. They are set here following the
 * convention of PLUTO's own radiative disk tests, after which:
 *
 *   GetTemperature(rho, prs)   returns the gas temperature in Kelvin
 *   Blackbody(T [K])           returns a_R T^4 in code energy density
 *   opacities                  are per unit mass, in code units
 * --------------------------------------------------------------------- */
#define X_MASSFRAC 0.70    // hydrogen mass fraction
#define Z_MASSFRAC 0.02    // metallicity
#define Y_MASSFRAC (1.0 - X_MASSFRAC - Z_MASSFRAC)   // helium mass fraction, derived

#define T_OPAC_MIN 1.0e4   // Kelvin, lower validity limit of Kramers opacities used below
#define T_GAS_MIN  1.0e4   // Kelvin, gas temperature floor, no recombination or molecular physics below

#define BETA_DECOUPLE 0.5  // v/c above which matter and radiation are decoupled, see UserDefOpacitiesAt
#define BETA_WIDTH    0.05 // width of that transition in v/c

/* ********************************************************************* */
static double MeanMolWeight (void)
/*!
 * Mean molecular weight for a fully ionized gas of the composition above.
 *********************************************************************** */
{
  return 1.0 / (2.0*X_MASSFRAC + 0.75*Y_MASSFRAC + 0.5*Z_MASSFRAC);
}

/* ********************************************************************* */
static void SetPhysicalUnits (void)
/*!
 * Set the radiation module constants in code units. Cheap and idempotent,
 * called from Init(), which PLUTO runs on restarts as well.
 *********************************************************************** */
{
#if RADIATION
  g_radiationConst = 4.0 * CONST_sigma / CONST_c
                     / (UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY);   // a_R / (rho0 v0^2), per K^4
  g_idealGasConst  = MeanMolWeight() * CONST_amu / CONST_kB
                     * UNIT_VELOCITY * UNIT_VELOCITY;                    // T [K] = g_idealGasConst * p / rho
  g_radC           = CONST_c / UNIT_VELOCITY;                            // exactly 1 for UNIT_VELOCITY = c
  g_reducedC       = g_inputParam[REDUCED_C] * g_radC;                   // must stay above all gas speeds
#endif
}

#if RADIATION
/* ********************************************************************* */
static double LTETemperature (double rho, double ptot)
/*!
 * Temperature at which gas and radiation in local thermodynamic equilibrium
 * together provide the total pressure, i.e. the root of
 *
 *   rho T / g_idealGasConst + g_radiationConst T^4 / 3 = ptot
 *
 * The left side is increasing and convex in T. Both the pure gas and the
 * pure radiation estimates lie above the root, so Newton started from the
 * smaller of the two converges monotonically from above.
 *********************************************************************** */
{
  double T, Tg, Tr, f, df;
  int    n;

  Tg = g_idealGasConst * ptot / rho;
  Tr = pow(3.0 * ptot / g_radiationConst, 0.25);
  T  = MIN(Tg, Tr);

  for (n = 0; n < 100; n++) {
    f  = rho * T / g_idealGasConst + g_radiationConst * T*T*T*T / 3.0 - ptot;
    df = rho / g_idealGasConst + 4.0 / 3.0 * g_radiationConst * T*T*T;
    T -= f / df;
    if (fabs(f) < 1.e-12 * ptot) break;
  }

  return T;
}
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
 * search if nothing outward is valid either. Leaves the bin
 * placeholder untouched only if no other bin is validated.
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
 
    if (r > g_inputParam[RD]) {
      g_rhoDiskProfile[b]   = AnalyticDiskDensity(r);
      if (g_rhoDiskProfile[b] <= 0.0) g_rhoDiskProfile[b] = g_rhoCoronaProfile[b];
      g_diskProfileValid[b] = 1;
    } else {
      g_rhoDiskProfile[b]   = g_rhoCoronaProfile[b];  // placeholder only in truncated disk part
      g_diskProfileValid[b] = 0;
    }
  }
 
  g_profilesInit = 1;
  PatchUnvalidatedDiskBins();   // ensures the first DiskFraction call does not see unphysical inner disk seed
}
 
/* ********************************************************************* */
void UpdateProfiles (const Data *d, Grid *grid)
/*!
 * Recompute the radially binned corona and disk density profiles:
 *
 *   1. Soft classify each active cell using the current profile from the
 *      previous update via DiskFraction: f = diskfrac
 *   2. Accumulate into the corona bins: (1-f)*rho*dV / (1-f)*dV
 *      Accumulate into the disk bins: f*rho*dV / f*dV
 *   3. Sum across all ranks, so all end up with the same global
 *      profiles regardless of domain decomposition in X1/X2/X3.
 *   4. Exponential moving average in time to avoid abrupt jumps
 *      applied independently to each profile.
 *
 * Bins which accumulated zero volume this step due to no weighted cells
 * being found for that profile retain their previous value. Because
 * profiles are built from the same DiskFraction classification,
 * they stay mutually consistent and track their evolution instead
 * of relying on a fixed initial profile.
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
 
  if (!g_profilesInit) InitProfiles(grid);   // restart safety net
 
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
 
    vi[RHO] = d->Vc[RHO][k][j][i];               // density variable
    vi[VX3] = d->Vc[VX3][k][j][i];               // rotation variable
    
    f  = DiskFraction(vi, r, grid->x[JDIR][j]);  // soft disk weight using old profiles
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
 
  /* Corona profile unchanged, still just holds its previous
   * value when unpopulated this step since it is present
   * nearly everywhere.
   */
  for (b = 0; b < NBINS_PROFILE; b++) {
    if (volC_glob[b] > 0.0) {
      double new_val = sumC_glob[b] / volC_glob[b];
      g_rhoCoronaProfile[b] = CORONA_EMA_ALPHA * new_val
                             + (1.0 - CORONA_EMA_ALPHA) * g_rhoCoronaProfile[b];
    }
  }
 
  /* Disk profile only accept this estimate step as real if the
   * disk weighted volume is a non negligible fraction of the total
   * total bin volume. Once a bin is validated it keeps being updated
   * normally even if a later step sees close to no disk material there,
   * then retaining its last known value, exactly as before.
   */
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
 
  /* Patch never validated bins by copying the nearest already
   * valid current value.
   */
  PatchUnvalidatedDiskBins();
 
  g_profilesLive = 1;   // from here on DiskFraction uses the adaptive sigmoid
}

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*!
 * Primary initialization of primitive variables across the domain.
 *********************************************************************** */
{
  double coeff, eps2, pc, rcyl;
  double lambda;

  SetPhysicalUnits();

  rcyl = x1 * sin(x2);
  eps2 = g_inputParam[EPS] * g_inputParam[EPS];
  coeff = 0.4 / eps2 * (1.0 / x1 - (1.0 - 2.5 * eps2) / rcyl);
  lambda = 2.2 / (1.0 + 2.56 * g_inputParam[BETAV] * g_inputParam[BETAV]);

  /* -------------------------------------------------------------------
     1. Static adiabatic corona in hydrostatic equilibrium
     ------------------------------------------------------------------- */
  v[RHO] = g_inputParam[RHOC] * pow(x1, -1.5);
  v[PRS] = 0.4 * g_inputParam[RHOC] * pow(x1, -2.5);
  pc     = v[PRS];

  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;

  /* -------------------------------------------------------------------
     2. Keplerian adiabatic disk after Kluzniak & Kita
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
    v[TRC] = 1.0;     // purely diagnostic disk tracer, does not feed back into physics anywhere
  } else {
    v[PRS] = 0.4 * g_inputParam[RHOC] * pow(x1, -2.5);
    v[TRC] = 0.0;     // purely diagnostic corona tracer, does not feed back into pysics anywhere
  }

  /* -------------------------------------------------------------------
     3. Magnetic dipole field setup
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
  /* cell centered default dipole */
  v[BX1] = 2.0 * g_inputParam[MU] * cos(x2) / (x1 * x1 * x1);
  v[BX2] = g_inputParam[MU] * sin(x2) / (x1 * x1 * x1);
  v[BX3] = 0.0;

  /* vector potential for constrained transport */
  #if MHD_FORMULATION == CONSTRAINED_TRANSPORT
   v[AX1] = 0.0;
   v[AX2] = 0.0;
   v[AX3] = g_inputParam[MU] * sin(x2) / (x1 * x1);
  #endif
#endif
#endif /* PHYSICS == MHD */

  #if RADIATION
  /* -------------------------------------------------------------------
     4. Radiation field in local thermodynamic equilibrium

     The analytic disk pressure above is the total pressure that holds
     the disk in equilibrium. With radiation enabled it is shared between
     gas and radiation at a common temperature, which for AGN parameters
     makes the disk radiation pressure dominated, as expected. Putting it
     all into the gas instead would give T ~ 1e9-1e10 K and, through a_R T^4,
     a radiation energy many orders of magnitude above the gas energy.

     The corona is transparent by construction, since opacities are gated
     by DiskFraction, so it remains purely adiabatic gas and starts with
     the minimum radiation energy the module allows.
     ------------------------------------------------------------------- */
  if (DiskFraction(v, x1, x2) > 0.5) {
    double T = LTETemperature(v[RHO], v[PRS]);
    v[PRS] = v[RHO] * T / g_idealGasConst;
    v[ENR] = Blackbody(T);
  } else {
    v[ENR] = RADIATION_MIN_ERAD;
  }
  v[FR1] = 0.;
  v[FR2] = 0.;
  v[FR3] = 0.;
  #endif /* RADIATION */
}

/* ********************************************************************* */
static int InitDiskCell (double x1, double x2)
/*!
 * Reproduces exactly the branch condition the initialization uses to decide
 * disk versus corona at problem setup. Kept as a single source of truth so
 * DiskFraction can initially match the former tracer exactly without ever
 * reading the tracer.
 *********************************************************************** */
{
  double rcyl, eps2, coeff, pc, p_disk;
 
  rcyl  = x1 * sin(x2);
  eps2  = g_inputParam[EPS] * g_inputParam[EPS];
  coeff = 0.4 / eps2 * (1.0 / x1 - (1.0 - 2.5 * eps2) / rcyl);
 
  pc     = 0.4 * g_inputParam[RHOC] * pow(x1, -2.5);   // corona pressure
  p_disk = eps2 * pow(coeff, 2.5);                     // disk pressure
 
  return (p_disk >= pc && rcyl > g_inputParam[RD]) ? 1 : 0;
}
 
/* ********************************************************************* */
double DiskFraction (double *v, double x1, double x2)
/*!
 * Continuous [0,1] disk/corona classifier used to gate anomalous
 * viscosity, resistivity, and opacities to disk material.
 *
 * At initialization, this returns exactly 0.0 or 1.0 like the tracer.
 *
 * Once updating the profiles has accumulated at least one step of real
 * disk/corona weighted data, DiskFraction uses two multiplicative criteria
 * to decide if a cell is disk material:
 * 
 * 1. DENSITY: Tracks two running radial references, a corona density 
 *    profile and a disk density profile,  and classifies by where the
 *    cell sits between them.
 * 
 * 2. ROTATION: Compares the local azimuthal velocity against
 *    the expected Keplerian velocity at the cylindrical radius.
 *
 * Both criteria use independent logistic functions, and their outputs
 * are multiplied to yield the final continuous fraction [0, 1].
 *********************************************************************** */
{
  double rho_ref_c, rho_ref_d, log_span, position, arg_den, arg_rot;
  double den_factor, rot_factor, rcyl, v_kep, rot_frac;
 
  if (!g_profilesLive) {
    /* Exact initial match to the former tracer, same boolean used before.
       Covers cases of both the very first call and the one right after
       analytic profiles have been seeded initially but no real cell has
       been classified through updating profiles yet. */
    return (double) InitDiskCell(x1, x2);
  }
 
  if (!g_profilesInit) {
    /* Should not normally happen, but guard anyway by falling back to the
       original static analytic corona/disk densities so behavior degrades
       gracefully. */
    rho_ref_c = AnalyticCoronaDensity(x1);
    rho_ref_d = AnalyticDiskDensity(x1);
    if (rho_ref_d <= 0.0) rho_ref_d = rho_ref_c;
  } else {
    rho_ref_c = GetCoronaRefDensity(x1);
    rho_ref_d = GetDiskRefDensity(x1);
  }
 
  /* --- 1. DENSITY SIGMOID --- */
  log_span = log10(MAX(rho_ref_d, 1.e-30)) - log10(MAX(rho_ref_c, 1.e-30));
  log_span = (fabs(log_span) < 1.e-12) ? 1.e-12 : log_span;  // guard degenerate span
 
  position = (log10(MAX(v[RHO], 1.e-30)) - log10(MAX(rho_ref_c, 1.e-30))) / log_span;
 
  arg_den = (position - (1.0 - 1.0/CORONA_THRESH_FAC)) / CORONA_SIGMOID_WIDTH;
  arg_den = MIN(MAX(arg_den, -50.0), 50.0);      // guard exponential overflow/underflow
  den_factor = 1.0 / (1.0 + exp(-arg_den));

  /* --- 2. ROTATION SIGMOID --- */
  rcyl  = x1 * sin(x2);
  v_kep = 1.0 / sqrt(MAX(rcyl, 1.e-12));    // ideal Newtonian Keplerian velocity matching injection step
  
  rot_frac = fabs(v[VX3]) / v_kep;    // measure absolute rotation fraction against Keplerian
  
  arg_rot = (rot_frac - ROTATION_THRESH_FAC) / ROTATION_SIGMOID_WIDTH;
  arg_rot = MIN(MAX(arg_rot, -50.0), 50.0);
  rot_factor = 1.0 / (1.0 + exp(-arg_rot));

  return den_factor * rot_factor;    // multiply factors to only classify rotating and dense material as disk
}
 
/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
{
  g_radGrid = grid;      // grid needed for user defined opacities
  InitProfiles(grid);    // initial seeding of profiles
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*!
 * Called by PLUTO at the cadence set before compilation. Used here to
 * periodically refresh the corona/disk radial reference profiles
 * consumed by DiskFraction.
 *********************************************************************** */
{
  UpdateProfiles(d, grid);
}

#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 * Static, curl free background magnetic field options.
 *********************************************************************** */
{
  /* Option A: black hole powerlaw field
  
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

  /* Option B: stellar dipole field */
  B0[0] = 2.0 * g_inputParam[MU] * cos(x2) / (x1 * x1 * x1);
  B0[1] = g_inputParam[MU] * sin(x2) / (x1 * x1 * x1);
  B0[2] = 0.0;                             

  /* Option C: quadrupole field 
  B0[0] = 1.5 * g_inputParam[MU] * (3.0 * cos(x2) * cos(x2) - 1.0) / (x1 * x1 * x1 * x1);
  B0[1] = 3.0 * g_inputParam[MU] * cos(x2) * sin(x2) / (x1 * x1 * x1 * x1);
  B0[2] = 0.0;
  */       

  /* Option D: octupole field 
  B0[0] = 2.0 * g_inputParam[MU] * (5.0 * cos(x2) * cos(x2) * cos(x2) - 3.0 * cos(x2)) / (x1 * x1 * x1 * x1 * x1);
  B0[1] = 0.5 * g_inputParam[MU] * (15.0 * cos(x2) * cos(x2) * sin(x2) - 3.0 * sin(x2)) / (x1 * x1 * x1 * x1 * x1);
  B0[2] = 0.0;
  */         
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid)
/*!
 * Custom user defined boundary condition dispatcher.
 *********************************************************************** */
{
  int i, j, k;
  double *x1, *x2, *x3, *r;
  double a1, a2, a, rcyl, eps2, coeff, lambda;
  double dvar1dr, dvar2dr, dvardr;
  double cs2, dden, dfact;
  
  RBox dom_box;

  /* -----------------------------------------------------------------
     Active Domain Internal Density Floor & Conservative Update
     ----------------------------------------------------------------- */
  if (side == 0) {    
    x1 = grid->xgc[IDIR];
    x2 = grid->xgc[JDIR];
    x3 = grid->xgc[KDIR];

    TOT_LOOP(k, j, i) {
      int convert_to_cons = 0;

      /* domain density floor enforcement */
      if (d->Vc[RHO][k][j][i] < g_inputParam[DFLOOR]) {
        dden = d->Vc[RHO][k][j][i];
        cs2  = g_gamma * d->Vc[PRS][k][j][i] / d->Vc[RHO][k][j][i];

        d->Vc[RHO][k][j][i] = g_inputParam[DFLOOR];     
        dfact = dden / d->Vc[RHO][k][j][i];

        /* preserve local speed of sound and momentum */
        d->Vc[PRS][k][j][i] = cs2 * d->Vc[RHO][k][j][i] / g_gamma;
        d->Vc[VX1][k][j][i] *= dfact;
        d->Vc[VX2][k][j][i] *= dfact;
        d->Vc[VX3][k][j][i] *= dfact;

        /* threshold to reset purely diagnostic tracer for corona */
        if (x2[j] < 0.5 * CONST_PI - atan(3.0 * g_inputParam[EPS]) ||
            x2[j] > 0.5 * CONST_PI + atan(3.0 * g_inputParam[EPS])) {
          d->Vc[TRC][k][j][i] = 0.0;
        }

        convert_to_cons = 1;
      }

#if RADIATION
      /* gas temperature floor, below which the model has no valid physics;
         stops surface layers that lose their radiation into the transparent
         corona from cooling catastrophically and collapsing */
      if (GetTemperature(d->Vc[RHO][k][j][i], d->Vc[PRS][k][j][i]) < T_GAS_MIN) {
        d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i] * T_GAS_MIN / g_idealGasConst;
        convert_to_cons = 1;
      }
#endif

      /* recompute conservative variables if primitives changed */
      if (convert_to_cons) {
        RBoxDefine(i, i, j, j, k, k, CENTER, &dom_box);
        PrimToCons3D(d->Vc, d->Uc, &dom_box, grid);
      }
    }
  }

  /* -----------------------------------------------------------------
     Inner Radial Boundary (Toggle Selected)
     ----------------------------------------------------------------- */
  if (side == X1_BEG) {
    if (box->vpos == CENTER) {

#if INNER_BOUNDARY == BOUNDARY_BH
      /* =============================================================
         OPTION 1: BLACK HOLE EVENT HORIZON ABSORBING DIODE BOUNDARY
         ============================================================= */
      X1_BEG_LOOP(k, j, i) {
        /* copy primitive variables from first computational cell */
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

        /* diode condition allows inflow but block outflow */
        if (d->Vc[VX1][k][j][i] > 0.0) {
          d->Vc[VX1][k][j][i] = 0.0;
        }

        /* ghost cell floor enforcement */
        if (d->Vc[RHO][k][j][i] < g_inputParam[DFLOOR]) {
          d->Vc[RHO][k][j][i] = g_inputParam[DFLOOR];
        }
        if (d->Vc[PRS][k][j][i] < g_inputParam[DFLOOR] * 1.0e-3) {
          d->Vc[PRS][k][j][i] = g_inputParam[DFLOOR] * 1.0e-3;
        }

#if RADIATION
        /* radiation may only stream into the horizon, never out of it */
        d->Vc[ENR][k][j][i] = d->Vc[ENR][k][j][IBEG];
        d->Vc[FR1][k][j][i] = MIN(d->Vc[FR1][k][j][IBEG], 0.0);
        d->Vc[FR2][k][j][i] = d->Vc[FR2][k][j][IBEG];
        d->Vc[FR3][k][j][i] = d->Vc[FR3][k][j][IBEG];
#endif
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

        /* fixed velocity field at stellar surface */
        d->Vc[VX1][k][j][i] = 0.0;
        d->Vc[VX2][k][j][i] = 0.0;
        d->Vc[VX3][k][j][i] = g_inputParam[OMEGA_STAR] * grid->x[IDIR][IBEG] * sin(x2[j]);

#if PHYSICS == MHD
        /* poloidal field matched, toroidal field set to zero */
        d->Vc[BX1][k][j][i] = d->Vc[BX1][k][j][IBEG];
        d->Vc[BX2][k][j][i] = d->Vc[BX2][k][j][IBEG];
        d->Vc[BX3][k][j][i] = 0.0;
#endif

        if (d->Vc[RHO][k][j][i] < g_inputParam[DFLOOR]) {
          d->Vc[RHO][k][j][i] = g_inputParam[DFLOOR];
        }

#if RADIATION
        /* zero gradient radiation field at the stellar surface */
        d->Vc[ENR][k][j][i] = d->Vc[ENR][k][j][IBEG];
        d->Vc[FR1][k][j][i] = d->Vc[FR1][k][j][IBEG];
        d->Vc[FR2][k][j][i] = d->Vc[FR2][k][j][IBEG];
        d->Vc[FR3][k][j][i] = d->Vc[FR3][k][j][IBEG];
#endif
      }
#endif

    }
  }

  /* -----------------------------------------------------------------
     Outer Radial Boundary (Disk Injection & Corona)
     ----------------------------------------------------------------- */
  if (side == X1_END) {
    r  = grid->x[IDIR];
    x1 = grid->x[IDIR];
    x2 = grid->x[JDIR];

    if (box->vpos == CENTER) {
      BOX_LOOP(box, k, j, i) {
        d->Vc[TRC][k][j][i] = d->Vc[TRC][k][j][IEND];

#if RADIATION
        /* radiation leaves the domain freely, nothing streams back in */
        d->Vc[ENR][k][j][i] = d->Vc[ENR][k][j][IEND];
        d->Vc[FR1][k][j][i] = MAX(d->Vc[FR1][k][j][IEND], 0.0);
        d->Vc[FR2][k][j][i] = d->Vc[FR2][k][j][IEND];
        d->Vc[FR3][k][j][i] = d->Vc[FR3][k][j][IEND];
#endif

        /* logarithmic extrapolation of density */
        a1 = log10(d->Vc[RHO][k][j][IEND]   / d->Vc[RHO][k][j][IEND-1]) / log10(r[IEND]   / r[IEND-1]);
        a2 = log10(d->Vc[RHO][k][j][IEND-1] / d->Vc[RHO][k][j][IEND-2]) / log10(r[IEND-1] / r[IEND-2]);
        a  = VANLEER_LIMITER(a1, a2);
        a  = MIN(a, 0.0);

        d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][i-1] * pow(r[i] / r[i-1], a); 
        d->Vc[PRS][k][j][i] = d->Vc[PRS][k][j][IEND] * pow(d->Vc[RHO][k][j][i] / d->Vc[RHO][k][j][IEND], g_gamma);

        /* outflow condition for poloidal velocity components */
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

        /* Reinject initial Kluźniak & Kita disk profile in equatorial region */
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

        /* prevent coronal backwards inflow */
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
 * Radial acceleration vector for point mass central potential.
 * This currently implements the Newtonian force only and is
 * currently unused.  
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
 * Select regime in definitions.
 *********************************************************************** */
{
  /* Option 1: Newtonian potential */
  // return -1.0 / x1;

  /* Option 2: Paczyński-Wiita potential */
  return -1.0 / (x1 - 2.0);

  /* Option 3: Kluźniak-Lee potential */
  // return -(1.0 / 6.0) * (exp(6.0 / x1) - 1.0);

  /* Option 4: Kluźniak-Nordström pseudo potential */
  // double q = 1.25;
  // return -1.0 / x1 + 0.5 * (q * q / (x1 * x1));
}
#endif

#if RADIATION_VAR_OPACITIES
/* ---------------------------------------------------------------------
 * Kramers and Electron Scattering Opacities (CGS)
 * ---------------------------------------------------------------------
 *   kappa_ff = C_FF g_ff (1 - Z) (1 + X)   rho T^-3.5   cm^2/g
 *   kappa_bf = C_BF g_bf Z (1 + X) / t     rho T^-3.5   cm^2/g
 *   kappa_es = 0.2 (1 + X)                              cm^2/g
 * --------------------------------------------------------------------- */
#define G_BF     1.0       // bound-free Gaunt factor
#define G_FF     1.0       // free-free Gaunt factor
#define T_FACTOR 1.0e1     // bound-free guillotine factor, typically 1 < t < 100

#define C_BF 4.34e25       // Kramers bound-free constant, CGS
#define C_FF 3.68e22       // Kramers free-free constant, CGS

static const double K_BF = C_BF * G_BF * Z_MASSFRAC * (1.0 + X_MASSFRAC) / T_FACTOR;
static const double K_FF = C_FF * G_FF * (1.0 - Z_MASSFRAC) * (1.0 + X_MASSFRAC);
static const double K_ES = 0.2 * (1.0 + X_MASSFRAC);

/* ********************************************************************* */
static double ComptonOpacity (double *v, double kappa_es)
/*!
 * Effective absorption opacity reproducing Compton energy exchange.
 *
 * PLUTO treats scattering as coherent, so electron scattering moves no
 * energy between gas and radiation. In an AGN disk this Compton exchange
 * is what keeps the gas at the radiation temperature, Kramers absorption
 * falls as T^-3.5 and alone cannot. The nonrelativistic Compton rate
 *
 *   L_C = rho kappa_es c E_r 4 k (T_g - T_r) / (m_e c^2) ,  T_r = (E_r/a_R)^(1/4)
 *
 * is matched exactly by an absorption term rho kappa_c c (a_R T_g^4 - E_r)
 * with
 *
 *   kappa_c = kappa_es 4 k T_r^4 / (m_e c^2 (T_g + T_r) (T_g^2 + T_r^2))
 *
 * after cancelling (T_g - T_r) analytically, so it is finite and positive
 * everywhere: kappa_es k T / (m_e c^2) in equilibrium, falling as T_g^-3
 * for hot gas, which keeps the exchange linear in T_g as it should be.
 * Needs v[ENR] to be current, see the patch in rad_step.c.
 *********************************************************************** */
{
  double Tg, Tr;

  Tg = GetTemperature(v[RHO], v[PRS]);
  Tr = pow(MAX(v[ENR], 0.0) / g_radiationConst, 0.25);

  if (Tg + Tr <= 0.0) return 0.0;

  return MIN(kappa_es * 4.0 * CONST_kB * Tr*Tr*Tr*Tr
             / (CONST_me * CONST_c * CONST_c * (Tg + Tr) * (Tg*Tg + Tr*Tr)),
             kappa_es);   // only reached for T_r above ~1e9 K, keeps scattering non negative
}

/* ********************************************************************* */
void UserDefOpacitiesAt(double *v, double x1, double x2, double *abs, double *scat)
/*!
 * Core opacity evaluation, taking position parameters explicitly. Called
 * from two places:
 *   1. For user defined opacities, the fixed signature entry point that
 *      the radiation module calls from radiation steps, it supplies x1,x2
 *      via the g_i_rad/g_j globals.
 *   2. For user defined ouputs, its diagnostic loop directly uses own
 *      real x1[i],x2[j] instead of the globals, which would be stale
 *      there, left over from whichever cell the last radiation implicit
 *      step visited, not the cell the diagnostic loop is currently on.
 *
 * PLUTO multiplies the returned coefficients by the density itself, see
 * e.g. dt*rho0*g_reducedC*abs_op in rad_step.c, so they are opacities per
 * unit mass. Converting cm^2/g to code units uses UNIT_DENSITY*UNIT_LENGTH,
 * so that kappa*rho*dx is the same optical depth in both systems.
 *
 * Opacities are gated by DiskFraction the same way viscosity and
 * resistivity are, in place of the previous tracer.
 *********************************************************************** */
{
  double rho_cgs, T, kappa_ffbf, kappa_es, kappa_c, f, beta;

  rho_cgs = v[RHO] * UNIT_DENSITY;                                // g/cm^3
  T       = MAX(GetTemperature(v[RHO], v[PRS]), T_OPAC_MIN);      // Kelvin

  kappa_ffbf = (K_BF + K_FF) * rho_cgs * pow(T, -3.5);            // cm^2/g
  kappa_es   = K_ES;                                              // cm^2/g
  kappa_c    = ComptonOpacity(v, kappa_es);                       // cm^2/g

  f = DiskFraction(v, x1, x2);

  /* The nonrelativistic radiation module assumes v << c, its v^2/c^2
     source corrections turn into anti damping once v approaches c, and
     with Newtonian dynamics in the Paczynski-Wiita potential the plunging
     gas next to the horizon reaches several c. Matter and radiation are
     therefore decoupled smoothly above BETA_DECOUPLE, which only affects
     the plunging region inside a few R_g. */
  beta = sqrt(v[VX1]*v[VX1] + v[VX2]*v[VX2] + v[VX3]*v[VX3]) / g_radC;
  f   *= 1.0 / (1.0 + exp(MIN((beta - BETA_DECOUPLE) / BETA_WIDTH, 50.0)));

  /* The Compton part only exchanges energy, so it is moved from scattering
     to absorption, leaving the total extinction at exactly Thomson. */
  *abs  = f * (kappa_ffbf + kappa_c)  * UNIT_DENSITY * UNIT_LENGTH;   // free-free + bound-free + Compton
  *scat = f * (kappa_es   - kappa_c)  * UNIT_DENSITY * UNIT_LENGTH;   // remaining coherent Thomson
}

/* ********************************************************************* */
void UserDefOpacities(double *v, double *abs, double *scat)
/*!
 * Fixed signature entry point required by radiation module, called from
 * the implicit step. Supplies x1 and x2 via g_i (set by the locally patched
 * radiation step immediately before this call) and g_j (PLUTO global,
 * already valid at that point).
 *********************************************************************** */
{
    double x1 = g_radGrid->x[IDIR][g_i_rad];
    double x2 = g_radGrid->x[JDIR][g_j];

    UserDefOpacitiesAt(v, x1, x2, abs, scat);
}
#endif /* RADIATION_VAR_OPACITIES */
