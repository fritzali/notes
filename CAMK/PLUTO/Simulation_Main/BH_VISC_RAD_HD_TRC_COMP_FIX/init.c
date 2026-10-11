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
  where the gas moves faster than the nonrelativistic radiation module can handle. The
  Compton rate is held fixed above 1e9 K, and absorption and Compton vanish below
  DiskFraction 0.05, while Thomson scattering applies to all gas.
  The DiskFraction classifier state is stored with every output and restored on restart.
  Its reference profiles are smoothed with a time constant, their logarithmic span is
  kept above a minimum so that bins cannot invert, and they can optionally be built
  from the density structure instead of the classification itself.
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

#define BETA_DECOUPLE 0.9  // v/c above which matter and radiation are decoupled, see UserDefOpacitiesAt
#define BETA_WIDTH    0.03 // width of that transition in v/c

#define T_COMPTON_MAX 1.0e9   // Kelvin, gas temperature above which the Compton rate is held fixed, see ComptonOpacity
#define F_OPAC_MIN    0.05    // DiskFraction below which all opacities vanish, see UserDefOpacitiesAt

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
static double ProfileWeight (void)
/*!
 * Weight of the newest estimate in the exponential moving average of the
 * reference profiles, for the time that passes between two updates. With
 * a time constant instead of a fixed weight per update, the memory of the
 * profiles is the same at every resolution and update cadence. g_dt is the
 * time step about to be taken, which restart.out stores, so a restarted
 * run uses exactly the same weight as the original one.
 *********************************************************************** */
{
  Runtime *rt = RuntimeGet();
  double   dt_upd;

  dt_upd = (rt->anl_dn > 0) ? rt->anl_dn * g_dt : rt->anl_dt;
  return 1.0 - exp(-MAX(dt_upd, 0.0) / DF_EMA_TAU);
}

/* ********************************************************************* */
static void RepairInvertedBins (void)
/*!
 * Replace the reference densities of every bin whose disk reference does
 * not exceed the corona reference by at least DF_SPAN_MIN dex. Such a bin
 * cannot separate disk from corona, and with the self referencing profiles
 * an inverted bin keeps itself inverted: dense gas classified as corona
 * raises the corona reference, thin gas classified as disk lowers the disk
 * reference. Both values are copied from the nearest bin with a sufficient
 * span, searching outward first, like PatchUnvalidatedDiskBins(). Without
 * any such bin, e.g. while the whole disk is disrupted, the disk reference
 * is raised to the minimum span instead. Changes of the number of repaired
 * bins are reported in the log.
 *********************************************************************** */
{
  static int n_last = 0;
  int    b, bb, src, n_bad = 0;
  int    bad[NBINS_PROFILE];
  double span;

  for (b = 0; b < NBINS_PROFILE; b++) {
    span   = log10(MAX(g_rhoDiskProfile[b], 1.e-30)) - log10(MAX(g_rhoCoronaProfile[b], 1.e-30));
    bad[b] = (span < DF_SPAN_MIN);
    n_bad += bad[b];
  }

  for (b = 0; b < NBINS_PROFILE; b++) {
    if (!bad[b]) continue;

    src = -1;
    for (bb = b + 1; bb < NBINS_PROFILE; bb++) {
      if (!bad[bb]) { src = bb; break; }
    }
    if (src < 0) {
      for (bb = b - 1; bb >= 0; bb--) {
        if (!bad[bb]) { src = bb; break; }
      }
    }

    if (src >= 0) {
      g_rhoCoronaProfile[b] = g_rhoCoronaProfile[src];
      g_rhoDiskProfile[b]   = g_rhoDiskProfile[src];
    } else {
      g_rhoDiskProfile[b]   = g_rhoCoronaProfile[b] * pow(10.0, DF_SPAN_MIN);
    }
  }

  if (n_bad != n_last) {
    print ("> RepairInvertedBins(): %d of %d bins below the minimum span at t = %.6e\n",
           n_bad, NBINS_PROFILE, g_time);
    n_last = n_bad;
  }
}

/* ********************************************************************* */
void UpdateProfiles (const Data *d, Grid *grid)
/*!
 * Recompute the radially binned corona and disk density profiles. The
 * reference densities of each bin are estimated in one of two ways,
 * chosen by DISKFRAC_REFERENCE in modifications.h:
 *
 *   DF_REF_SELF  Soft classify each active cell with DiskFraction using
 *                the current profiles, f = diskfrac, and accumulate
 *                (1-f)*rho*dV / (1-f)*dV for the corona and
 *                f*rho*dV / f*dV for the disk. The profiles follow the
 *                disk material wherever it is, including lifted gas.
 *
 *   DF_REF_GEOM  Take the maximum density of the bin for the disk and the
 *                volume weighted geometric mean density of the cells
 *                above z/R = DF_GEOM_ZCOR for the corona. Independent of
 *                the classification, so it cannot feed back on itself.
 *
 * Both are summed across all ranks, so all end up with the same global
 * profiles regardless of domain decomposition in X1/X2/X3, then smoothed
 * by an exponential moving average with time constant DF_EMA_TAU, and
 * finally checked for bins without a usable span, see RepairInvertedBins().
 *
 * Bins without any contribution this step retain their previous value.
 *********************************************************************** */
{
  static double sumC_loc[NBINS_PROFILE], volC_loc[NBINS_PROFILE];
  static double sumC_glob[NBINS_PROFILE], volC_glob[NBINS_PROFILE];
  static double sumD_loc[NBINS_PROFILE], volD_loc[NBINS_PROFILE];
  static double sumD_glob[NBINS_PROFILE], volD_glob[NBINS_PROFILE];
  static double volTot_loc[NBINS_PROFILE], volTot_glob[NBINS_PROFILE];

  int    i, j, k, b;
  double *x1 = grid->x[IDIR];
  double r, dV, alpha, new_val;
#if DISKFRAC_REFERENCE == DF_REF_SELF
  double f, wc, wd, vi[NVAR];
#else
  double rho, theta;
#endif

  if (!g_profilesInit) InitProfiles(grid);   // restart safety net
  alpha = ProfileWeight();

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

    dV = grid->dV[k][j][i];

#if DISKFRAC_REFERENCE == DF_REF_SELF
    vi[RHO] = d->Vc[RHO][k][j][i];               // density variable
    vi[VX3] = d->Vc[VX3][k][j][i];               // rotation variable

    f  = DiskFraction(vi, r, grid->x[JDIR][j]);  // soft disk weight using old profiles
    wc = 1.0 - f;
    wd = f;

    sumC_loc[b] += wc * vi[RHO] * dV;
    volC_loc[b] += wc * dV;

    sumD_loc[b] += wd * vi[RHO] * dV;
    volD_loc[b] += wd * dV;
#else
    rho   = d->Vc[RHO][k][j][i];
    theta = grid->x[JDIR][j];

    if (cos(theta) > DF_GEOM_ZCOR * sin(theta)) {   // z/R = cot(theta) above DF_GEOM_ZCOR
      sumC_loc[b] += log(MAX(rho, 1.e-30)) * dV;
      volC_loc[b] += dV;
    }
    sumD_loc[b] = MAX(sumD_loc[b], rho);           // maximum density of the bin
#endif

    volTot_loc[b] += dV;
  }

#ifdef PARALLEL
  MPI_Allreduce(sumC_loc, sumC_glob, NBINS_PROFILE, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(volC_loc, volC_glob, NBINS_PROFILE, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
 #if DISKFRAC_REFERENCE == DF_REF_SELF
  MPI_Allreduce(sumD_loc, sumD_glob, NBINS_PROFILE, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
 #else
  MPI_Allreduce(sumD_loc, sumD_glob, NBINS_PROFILE, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
 #endif
  MPI_Allreduce(volD_loc, volD_glob, NBINS_PROFILE, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(volTot_loc, volTot_glob, NBINS_PROFILE, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  for (b = 0; b < NBINS_PROFILE; b++) {
    sumC_glob[b] = sumC_loc[b]; volC_glob[b] = volC_loc[b];
    sumD_glob[b] = sumD_loc[b]; volD_glob[b] = volD_loc[b];
    volTot_glob[b] = volTot_loc[b];
  }
#endif

  /* Corona profile, holds its previous value when unpopulated this step. */
  for (b = 0; b < NBINS_PROFILE; b++) {
    if (volC_glob[b] > 0.0) {
#if DISKFRAC_REFERENCE == DF_REF_SELF
      new_val = sumC_glob[b] / volC_glob[b];
#else
      new_val = exp(sumC_glob[b] / volC_glob[b]);
#endif
      g_rhoCoronaProfile[b] = alpha * new_val + (1.0 - alpha) * g_rhoCoronaProfile[b];
    }
  }

  /* Disk profile. In the self referencing mode an estimate only counts if
   * the disk weighted volume is a non negligible fraction of the total bin
   * volume. Once a bin is validated it keeps being updated normally even
   * if a later step sees close to no disk material there, then retaining
   * its last known value. In the geometric mode every populated bin counts.
   */
  for (b = 0; b < NBINS_PROFILE; b++) {
#if DISKFRAC_REFERENCE == DF_REF_SELF
    int meaningful = (volTot_glob[b] > 0.0)
                   && (volD_glob[b] > DISK_BIN_VOLFRAC_MIN * volTot_glob[b]);
#else
    int meaningful = (volTot_glob[b] > 0.0);
#endif

    if (meaningful) {
#if DISKFRAC_REFERENCE == DF_REF_SELF
      new_val = sumD_glob[b] / volD_glob[b];
#else
      new_val = sumD_glob[b];
#endif
      g_rhoDiskProfile[b] = alpha * new_val + (1.0 - alpha) * g_rhoDiskProfile[b];
      g_diskProfileValid[b] = 1;
    }
  }

  /* Patch never validated bins by copying the nearest already
   * valid current value, then repair bins without a usable span.
   */
  PatchUnvalidatedDiskBins();
  RepairInvertedBins();

  g_profilesLive = 1;   // from here on DiskFraction uses the adaptive sigmoid
}

/* ---------------------------------------------------------------------
 * Classifier State Across Restarts
 * ---------------------------------------------------------------------
 * The reference profiles carry memory through their temporal smoothing,
 * so they are not a function of the current state alone. To restart
 * with exactly the classification the run had, one record per output
 * is appended to DISKFRAC_FILE in the output directory, holding the
 * time, the two flags and all bin values. On restart, the record with
 * the restart time is read back.
 * --------------------------------------------------------------------- */
#define DISKFRAC_FILE "diskfrac.out"

/* ********************************************************************* */
static int AnalysisDue (void)
/*!
 * Same condition as CheckForAnalysis() in main.c, evaluated for the
 * current time, time step and step number.
 *********************************************************************** */
{
  Runtime *rt = RuntimeGet();
  double   t  = g_time, tnext = g_time + g_dt;
  int      due_dt, due_dn;

  due_dt = ((int)(tnext / rt->anl_dt) - (int)(t / rt->anl_dt))
           || g_stepNumber == 0 || fabs(t - rt->tstop) < 1.e-9;
  due_dt = due_dt && (rt->anl_dt > 0.0);
  due_dn = (rt->anl_dn > 0) && (g_stepNumber % rt->anl_dn == 0);

  return due_dt || due_dn;
}

/* ********************************************************************* */
void SaveDiskProfiles (void)
/*!
 * Append the current classifier state to DISKFRAC_FILE. Called from
 * ComputeUserVar(), which PLUTO runs once per step right before writing
 * outputs, so there is one record per output time. The profiles are
 * identical on all ranks after the reduction in UpdateProfiles(), so
 * only rank 0 writes. A fresh run, writing its initial output at step
 * zero, starts a new file, while a restarted run appends.
 *********************************************************************** */
{
  static double t_last = -1.0;
  char  fname[512];
  FILE *fp;
  int   b;

  if (g_time == t_last) return;
  t_last = g_time;
  if (prank != 0) return;

  snprintf(fname, sizeof(fname), "%s/%s", RuntimeGet()->output_dir, DISKFRAC_FILE);
  fp = fopen(fname, g_stepNumber == 0 ? "w" : "a");
  if (fp == NULL) {
    print ("! SaveDiskProfiles(): cannot open %s\n", fname);
    return;
  }

  fprintf(fp, "%.17e %d %d %d", g_time, g_profilesInit, g_profilesLive, NBINS_PROFILE);
  for (b = 0; b < NBINS_PROFILE; b++) fprintf(fp, " %.17e", g_rhoCoronaProfile[b]);
  for (b = 0; b < NBINS_PROFILE; b++) fprintf(fp, " %.17e", g_rhoDiskProfile[b]);
  for (b = 0; b < NBINS_PROFILE; b++) fprintf(fp, " %d", g_diskProfileValid[b]);
  fprintf(fp, "\n");
  fclose(fp);
}

/* ********************************************************************* */
static int ReadDiskProfiles (double t, double *rc, double *rd, int *valid, int *flags)
/*!
 * Search DISKFRAC_FILE for a record at time t and copy it into the
 * given arrays. Restart times are restored exactly from restart.out
 * and records are written with 17 significant digits, so the match is
 * exact up to a tiny tolerance. After repeated restarts several records
 * can share a time, the last one belongs to the latest run and wins.
 * Returns 1 when a record was found.
 *********************************************************************** */
{
  char   fname[512];
  FILE  *fp;
  double tt, crc[NBINS_PROFILE], crd[NBINS_PROFILE];
  int    cvalid[NBINS_PROFILE], cinit, clive, nb, b, ok, found = 0;

  snprintf(fname, sizeof(fname), "%s/%s", RuntimeGet()->output_dir, DISKFRAC_FILE);
  fp = fopen(fname, "r");
  if (fp == NULL) return 0;

  while (fscanf(fp, "%lf %d %d %d", &tt, &cinit, &clive, &nb) == 4) {
    if (nb != NBINS_PROFILE) {
      print ("! ReadDiskProfiles(): record with %d bins, expected %d\n", nb, NBINS_PROFILE);
      break;
    }
    ok = 1;
    for (b = 0; b < nb && ok; b++) ok = (fscanf(fp, "%lf", &crc[b]) == 1);
    for (b = 0; b < nb && ok; b++) ok = (fscanf(fp, "%lf", &crd[b]) == 1);
    for (b = 0; b < nb && ok; b++) ok = (fscanf(fp, "%d", &cvalid[b]) == 1);
    if (!ok) break;

    if (fabs(tt - t) <= 1.e-12 * MAX(1.0, fabs(t))) {
      for (b = 0; b < nb; b++) {
        rc[b] = crc[b];
        rd[b] = crd[b];
        valid[b] = cvalid[b];
      }
      flags[0] = cinit;
      flags[1] = clive;
      found = 1;
    }
  }
  fclose(fp);
  return found;
}

/* ********************************************************************* */
static void RestoreDiskProfiles (const Data *d, Grid *grid)
/*!
 * Restore the classifier state when the run was restarted. Called at the
 * top of UserDefBoundary(), which runs before anything else of a step
 * uses DiskFraction.
 *
 * Startup() already calls the boundaries once at t = 0, before the
 * restart file is read, so the check waits for the first call with
 * g_time > 0. A fresh run gets there after its first step, with step
 * number one, while a restarted run starts from the restored step
 * number of its output, so the two are told apart by that number.
 *
 * In the main loop the output is written before Analysis() updates the
 * profiles in the same step, so a stored record is one update behind
 * what the original run used next. A restarted run also skips output and
 * analysis on its first iteration. The update is therefore repeated here
 * on the restored state, whenever Analysis() was due at that step, which
 * reproduces the original profiles exactly.
 *********************************************************************** */
{
  static int checked = 0;
  int    flags[2] = {0, 0}, found = 0;

  if (checked || g_time <= 0.0) return;
  checked = 1;
  if (g_stepNumber <= 1) return;           // fresh run, nothing to restore
  if (!g_profilesInit) InitProfiles(grid);   // bin edges, before the stored values are read

  if (prank == 0) {
    found = ReadDiskProfiles(g_time, g_rhoCoronaProfile, g_rhoDiskProfile, g_diskProfileValid, flags);
  }
#ifdef PARALLEL
  MPI_Bcast(&found, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (found) {
    MPI_Bcast(g_rhoCoronaProfile, NBINS_PROFILE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(g_rhoDiskProfile,   NBINS_PROFILE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(g_diskProfileValid, NBINS_PROFILE, MPI_INT,    0, MPI_COMM_WORLD);
    MPI_Bcast(flags, 2, MPI_INT, 0, MPI_COMM_WORLD);
  }
#endif

  if (!found) {
    print ("! RestoreDiskProfiles(): no record for t = %.10e in %s, "
           "re-seeding from the analytic profiles\n", g_time, DISKFRAC_FILE);
    return;
  }

  g_profilesInit = flags[0];
  g_profilesLive = flags[1];
  print ("> RestoreDiskProfiles(): restored classifier state at t = %.10e\n", g_time);

  if (AnalysisDue()) UpdateProfiles(d, grid);
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
  log_span = MAX(log_span, DF_SPAN_MIN);   // never let a degenerate or inverted span flip the density criterion
 
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

  RestoreDiskProfiles(d, grid);   // only acts once, on the first step after a restart

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
 *
 * The rate above is nonrelativistic and only valid for k T_g << m_e c^2,
 * i.e. T_g well below ~6e9 K. Hot corona gas starts near 1e11 K, where it
 * would overestimate the exchange enormously. Above T_COMPTON_MAX the rate
 * is therefore held at its value for T_g = T_COMPTON_MAX, which requires
 *
 *   kappa_c = kappa_es 4 k T_r^4 (T_max - T_r) / (m_e c^2 (T_g^4 - T_r^4))
 *
 * so that rho kappa_c c (a_R T_g^4 - E_r) = L_C(T_max). Merely capping T_g
 * inside the formula above would instead let the absorption form grow as
 * T_g^4. Both branches agree at T_g = T_COMPTON_MAX.
 *********************************************************************** */
{
  double Tg, Tr, Tr4, mc2, kappa_c;

  Tg  = GetTemperature(v[RHO], v[PRS]);
  Tr  = pow(MAX(v[ENR], 0.0) / g_radiationConst, 0.25);
  Tr4 = Tr*Tr*Tr*Tr;
  mc2 = CONST_me * CONST_c * CONST_c;

  if (Tg + Tr <= 0.0) return 0.0;

  if (Tg <= T_COMPTON_MAX) {
    kappa_c = kappa_es * 4.0 * CONST_kB * Tr4 / (mc2 * (Tg + Tr) * (Tg*Tg + Tr*Tr));
  } else if (Tr < T_COMPTON_MAX) {
    kappa_c = kappa_es * 4.0 * CONST_kB * Tr4 * (T_COMPTON_MAX - Tr)
              / (mc2 * (Tg*Tg*Tg*Tg - Tr4));    // Tg > T_COMPTON_MAX > Tr, so positive
  } else {
    kappa_c = 0.0;                              // radiation itself above the cap, no valid rate
  }

  return MIN(kappa_c, kappa_es);   // keeps scattering non negative
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
 * Only the opacities that exchange energy between gas and radiation, Kramers
 * absorption and Compton, are gated by DiskFraction, in place of the previous
 * tracer, and vanish entirely below F_OPAC_MIN. Thomson scattering only
 * transfers momentum and applies to all ionized gas, so it is not gated: when
 * it was, radiation could only push on gas classified as disk, and gas lifted
 * across the classification edge lost its support and fell back, piling up
 * into a dense, Rayleigh-Taylor unstable sheet at the disk surface.
 *********************************************************************** */
{
  double rho_cgs, T, kappa_ffbf, kappa_es, kappa_c, fd, gv, beta;

  rho_cgs = v[RHO] * UNIT_DENSITY;                                // g/cm^3
  T       = MAX(GetTemperature(v[RHO], v[PRS]), T_OPAC_MIN);      // Kelvin

  kappa_ffbf = (K_BF + K_FF) * rho_cgs * pow(T, -3.5);            // cm^2/g
  kappa_es   = K_ES;                                              // cm^2/g
  kappa_c    = ComptonOpacity(v, kappa_es);                       // cm^2/g

  /* The nonrelativistic radiation module assumes v << c, its v^2/c^2
     source corrections turn into anti damping once v approaches c, and
     with Newtonian dynamics in the Paczynski-Wiita potential the plunging
     gas next to the horizon reaches several c. Matter and radiation are
     therefore decoupled smoothly above BETA_DECOUPLE. The threshold has to
     stay well above the orbital speed, which is already 0.6 c at R = 6 in
     this potential, so that only the plunging region is affected. These
     terms scale with the total opacity, so scattering is included. */
  beta = sqrt(v[VX1]*v[VX1] + v[VX2]*v[VX2] + v[VX3]*v[VX3]) / g_radC;
  gv   = 1.0 / (1.0 + exp(MIN((beta - BETA_DECOUPLE) / BETA_WIDTH, 50.0)));

  /* DiskFraction never reaches exactly zero in the corona, the rotation
     sigmoid alone leaves about 1e-3 there. Combined with the hot corona
     gas this was enough to Compton cool the whole corona, which then lost
     its pressure support and drained into the hole. The energy exchange
     therefore vanishes below F_OPAC_MIN, and the remaining range is stretched
     back onto [0,1] so that disk material keeps its full coupling. */
  fd = DiskFraction(v, x1, x2);
  fd = (fd > F_OPAC_MIN) ? (fd - F_OPAC_MIN) / (1.0 - F_OPAC_MIN) : 0.0;

  /* The Compton part only exchanges energy, so it is moved from scattering
     to absorption where it applies, leaving the total extinction at exactly
     Thomson plus Kramers everywhere. */
  *abs  = gv * fd * (kappa_ffbf + kappa_c)      * UNIT_DENSITY * UNIT_LENGTH;   // free-free + bound-free + Compton, disk only
  *scat = gv * (kappa_es - fd * kappa_c)        * UNIT_DENSITY * UNIT_LENGTH;   // coherent Thomson, all gas
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
