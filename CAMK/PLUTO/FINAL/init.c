/* ///////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Contains basic functions for problem initialization.

  The init.c file collects most of the user-supplied functions useful
  for problem configuration.
  It is automatically searched for by the makefile.

  \author A. Mignone (mignone@ph.unito.it)
  \date   Sep 10, 2012

 \modified by M. Cemeljic (miki@camk.edu.pl)
 \date August 2025
 Setup for pluto 4.4.3 by Miki, plus Sukalpa Kundu's cooling. Based
 on the setup described in Appendix of "Atlas" paper, Cemeljic, 2019,
 A&A, 624, A31.

 \modified by F. A. Agildere (fritz.agildere@udo.edu)
 \date October 2026
 Black hole version: absorbing (diode) inner boundary, Paczynski-Wiita
 potential, units set by BH_MASS. Radiation (M1, nonrel.) with Kramers,
 Thomson and Compton as effective absorption. The passive tracer does
 not follow the disc once radiation is on, so viscosity and opacities
 use DiskFraction() instead, which classifies disc vs corona from
 density and rotation, against radial reference profiles updated in
 Analysis(). Tunables in modifications.h.
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include "modifications.h"

#ifdef PARALLEL
 #include <mpi.h>
#endif

/* cfa--code units from definitions.h: UNIT_LENGTH=R_g, UNIT_VELOCITY=c,
   so GM=1 as in the corona and KK disc below. Constants of the radiation
   module are set in SetPhysicalUnits(), afterwards GetTemperature()
   gives T in K and Blackbody() a_R T^4 in code units. */
#define X_MASSFRAC 0.70    // hydrogen
#define Z_MASSFRAC 0.02    // metals
#define Y_MASSFRAC (1.0 - X_MASSFRAC - Z_MASSFRAC)

#define T_OPAC_MIN 1.0e4   // K, lower limit for Kramers
#define T_GAS_MIN  1.0e4   // K, gas temperature floor

#define BETA_DECOUPLE 0.9  // v/c where matter and radiation decouple
#define BETA_WIDTH    0.03

#define T_COMPTON_MAX 1.0e9   // K, Compton rate fixed above this T_g
#define F_OPAC_MIN    0.05    // no absorption/Compton below this DiskFraction

/* ************************************************************** */
static double MeanMolWeight (void)
/* fully ionized gas */
{
  return 1.0 / (2.0*X_MASSFRAC + 0.75*Y_MASSFRAC + 0.5*Z_MASSFRAC);
}

/* ************************************************************** */
static void SetPhysicalUnits (void)
/* cfa--radiation constants in code units, called from Init() so it
   is also set on restart */
{
#if RADIATION
  g_radiationConst = 4.0 * CONST_sigma / CONST_c
                     / (UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY);
  g_idealGasConst  = MeanMolWeight() * CONST_amu / CONST_kB
                     * UNIT_VELOCITY * UNIT_VELOCITY;
  g_radC           = CONST_c / UNIT_VELOCITY;    // =1 here
  g_reducedC       = g_inputParam[REDUCED_C] * g_radC;
#endif
}

#if RADIATION
/* ************************************************************** */
static double LTETemperature (double rho, double ptot)
/* cfa--T where gas + radiation pressure in LTE give ptot,
   rho T/g_idealGasConst + a_R T^4/3 = ptot. Newton from the smaller
   of the pure gas and pure radiation guess, converges from above. */
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

/* cfa--choice of inner boundary: star surface or BH diode */
#define BOUNDARY_STAR 1
#define BOUNDARY_BH   2

#ifndef INNER_BOUNDARY
 #define INNER_BOUNDARY BOUNDARY_BH
#endif

/* cfa--radial index of the cell in the implicit radiation step, set in
   rad_step.c, so the opacities can call DiskFraction() at (x1,x2) */
int          g_i_rad   = 0;
static Grid *g_radGrid = NULL;

/* cfa--DiskFraction reference profiles, NBINS_PROFILE log bins in r,
   same on all ranks after the reduction in UpdateProfiles() */
static double g_rBinEdgesLog[NBINS_PROFILE + 1];
static double g_rhoCoronaProfile[NBINS_PROFILE];
static double g_rhoDiskProfile[NBINS_PROFILE];
static int    g_diskProfileValid[NBINS_PROFILE];   // bin has seen disc material

static int    g_profilesInit = 0;
static int    g_profilesLive = 0;   // 0: DiskFraction = initial disc cell (as tracer)

/* ************************************************************** */
static double InterpLogProfile (double *profile, double x1)
/* linear in log r and log rho, edge bins outside */
{
  double lr, lmin, lmax, dl, lc0, s, frac;
  double lrho0, lrho1, lrho;
  int    ib0, ib1;

  lr   = log10(x1);
  lmin = g_rBinEdgesLog[0];
  lmax = g_rBinEdgesLog[NBINS_PROFILE];
  dl   = (lmax - lmin) / (double)NBINS_PROFILE;
  lc0  = lmin + 0.5 * dl;

  s    = (lr - lc0) / dl;
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

static double GetCoronaRefDensity (double x1)
{
  return InterpLogProfile(g_rhoCoronaProfile, x1);
}

static double GetDiskRefDensity (double x1)
{
  return InterpLogProfile(g_rhoDiskProfile, x1);
}

/* analytic corona and midplane disc densities, same as in Init(),
   used to seed the profiles */
static double AnalyticCoronaDensity (double r)
{
  return g_inputParam[RHOC] * pow(r, -1.5);
}

static double AnalyticDiskDensity (double r)
{
  double eps2, coeff;

  eps2  = g_inputParam[EPS] * g_inputParam[EPS];
  coeff = 0.4 / eps2 * (1.0 / r - (1.0 - 2.5 * eps2) / r);
  coeff = MAX(coeff, 0.0);

  return pow(coeff, 1.5);
}

/* ************************************************************** */
static void PatchUnvalidatedDiskBins (void)
/* cfa--bins without disc material so far copy the nearest valid bin,
   searching outward first (disc comes in from larger r) */
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

/* ************************************************************** */
void InitProfiles (Grid *grid)
/* cfa--bin edges over the global radial domain, profiles seeded with
   the analytic densities */
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
      g_rhoDiskProfile[b]   = g_rhoCoronaProfile[b];   // inside RD, no disc yet
      g_diskProfileValid[b] = 0;
    }
  }

  g_profilesInit = 1;
  PatchUnvalidatedDiskBins();
}

/* ************************************************************** */
static double ProfileWeight (void)
/* cfa--EMA weight 1-exp(-dt/tau) for the time between two updates,
   same memory at any resolution. Was a fixed 0.05 per update before. */
{
  Runtime *rt = RuntimeGet();
  double   dt_upd;

  dt_upd = (rt->anl_dn > 0) ? rt->anl_dn * g_dt : rt->anl_dt;
  return 1.0 - exp(-MAX(dt_upd, 0.0) / DF_EMA_TAU);
}

/* ************************************************************** */
static void RepairInvertedBins (void)
/* cfa--bins with log(rho_d/rho_c) < DF_SPAN_MIN copy both references
   from the nearest good bin (outward first), or get rho_d raised to
   the minimum span. Self referencing profiles otherwise keep an
   inverted bin inverted (seen at r~2.4-3.4 in the 120x90 run). */
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

/* ************************************************************** */
void UpdateProfiles (const Data *d, Grid *grid)
/* cfa--new reference densities per bin.
   DF_REF_SELF: (1-f) and f weighted mean density, f from the current
                profiles, so they follow the disc material.
   DF_REF_GEOM: max density of the bin for the disc, geometric mean
                above z/R=DF_GEOM_ZCOR for the corona (less tested).
   Summed over all ranks, then EMA, patch and repair. Bins without
   contributions keep their old value. */
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

  if (!g_profilesInit) InitProfiles(grid);
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
    vi[RHO] = d->Vc[RHO][k][j][i];   // only rho and vphi are used
    vi[VX3] = d->Vc[VX3][k][j][i];

    f  = DiskFraction(vi, r, grid->x[JDIR][j]);
    wc = 1.0 - f;
    wd = f;

    sumC_loc[b] += wc * vi[RHO] * dV;
    volC_loc[b] += wc * dV;

    sumD_loc[b] += wd * vi[RHO] * dV;
    volD_loc[b] += wd * dV;
#else
    rho   = d->Vc[RHO][k][j][i];
    theta = grid->x[JDIR][j];

    if (cos(theta) > DF_GEOM_ZCOR * sin(theta)) {   // z/R > DF_GEOM_ZCOR
      sumC_loc[b] += log(MAX(rho, 1.e-30)) * dV;
      volC_loc[b] += dV;
    }
    sumD_loc[b] = MAX(sumD_loc[b], rho);
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

/* disc estimate counts only if its f weighted volume is a non
   negligible part of the bin (self mode) */
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

  PatchUnvalidatedDiskBins();
  RepairInvertedBins();

  g_profilesLive = 1;
}

/* cfa--the profiles have memory (EMA), so they are written with every
   output and read back on restart, one record per output time */
#define DISKFRAC_FILE "diskfrac.out"

/* ************************************************************** */
static int AnalysisDue (void)
/* same test as CheckForAnalysis() in main.c */
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

/* ************************************************************** */
void SaveDiskProfiles (void)
/* cfa--called from ComputeUserVar(), i.e. once per output. Only rank 0
   writes, new file at step 0, appended after a restart. */
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

/* ************************************************************** */
static int ReadDiskProfiles (double t, double *rc, double *rd, int *valid, int *flags)
/* cfa--record at time t, the last one wins after several restarts */
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

/* ************************************************************** */
static void RestoreDiskProfiles (const Data *d, Grid *grid)
/* cfa--called at the top of UserDefBoundary(). Startup() calls the
   boundaries once at t=0 before reading the restart, so wait for the
   first call with t>0. Fresh run is then at step 1, a restart at the
   step of its output. The stored record is one update behind (output
   is written before Analysis()), so the update is repeated here when
   it was due, which gives back the same profiles. */
{
  static int checked = 0;
  int    flags[2] = {0, 0}, found = 0;

  if (checked || g_time <= 0.0) return;
  checked = 1;
  if (g_stepNumber <= 1) return;
  if (!g_profilesInit) InitProfiles(grid);

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

/* ************************************************************** */
void Init (double *v, double x1, double x2, double x3)
/*
 *
 **************************************************************** */
{
  double coeff, eps2, pc, rcyl;
  double lambda;

  SetPhysicalUnits();

  rcyl = x1 * sin(x2);
  eps2 = g_inputParam[EPS] * g_inputParam[EPS];
  coeff = 0.4 / eps2 * (1.0 / x1 - (1.0 - 2.5 * eps2) / rcyl);
  lambda = 2.2 / (1.0 + 2.56 * g_inputParam[BETAV] * g_inputParam[BETAV]);

/* initial non-rotating adiabatic corona in hydrostatic equilibrium  */
  v[RHO] = g_inputParam[RHOC] * pow(x1, -1.5);
  v[PRS] = 0.4 * g_inputParam[RHOC] * pow(x1, -2.5);
  pc     = v[PRS];

  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;

/* Keplerian adiabatic disk in vertical pressure equilibrium with the
   adiabatic corona, as given by Kluzniak & Kita (2000) */
  v[PRS] = eps2 * pow(coeff, 2.5);

  if (v[PRS] >= pc && rcyl > g_inputParam[RD]) {
    v[RHO] = pow(coeff, 1.5);
    v[VX1] = -g_inputParam[BETAV] / sin(x2) * eps2 * (10.0 - (32.0 / 3.0)
             * lambda * g_inputParam[BETAV] * g_inputParam[BETAV]
             - lambda * (5.0 - 1.0 / (eps2 * tan(x2) * tan(x2)))) / sqrt(rcyl);
    v[VX3] = (sqrt(1.0 - 2.5 * eps2) + (2.0 / 3.0) * eps2
             * g_inputParam[BETAV] * g_inputParam[BETAV]
             * lambda * (1.0 - 1.2 / (eps2 * tan(x2) * tan(x2)))) / sqrt(rcyl);
    v[TRC] = 1.0;     /* Track the disc material */
  } else {
    v[PRS] = 0.4 * g_inputParam[RHOC] * pow(x1, -2.5);
    v[TRC] = 0.0;     /* Track the corona */
  }
//cfa--TRC is only a diagnostic now, the physics uses DiskFraction()

#if PHYSICS == MHD
#if BACKGROUND_FIELD == YES
  v[BX1] = 0.0;
  v[BX2] = 0.0;
  v[BX3] = 0.0;

  v[AX1] = 0.0;
  v[AX2] = 0.0;
  v[AX3] = 0.0;
#else
//ccm--dipole stellar field without background field
  v[BX1] = 2.0 * g_inputParam[MU] * cos(x2) / (x1 * x1 * x1);
  v[BX2] = g_inputParam[MU] * sin(x2) / (x1 * x1 * x1);
  v[BX3] = 0.0;

  #if MHD_FORMULATION == CONSTRAINED_TRANSPORT
   v[AX1] = 0.0;
   v[AX2] = 0.0;
   v[AX3] = g_inputParam[MU] * sin(x2) / (x1 * x1);
  #endif
#endif
#endif//if MHD loop

  #if RADIATION
//cfa--the KK disc pressure is the total one, shared by gas and radiation
// at a common T (LTE), so the disc starts radiation pressure dominated.
// Corona starts without radiation (transparent, opacities gated off).
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
  #endif
}

/* ************************************************************** */
static int InitDiskCell (double x1, double x2)
/* cfa--same condition as for the disc in Init(), so DiskFraction()
   starts exactly as the tracer did */
{
  double rcyl, eps2, coeff, pc, p_disk;

  rcyl  = x1 * sin(x2);
  eps2  = g_inputParam[EPS] * g_inputParam[EPS];
  coeff = 0.4 / eps2 * (1.0 / x1 - (1.0 - 2.5 * eps2) / rcyl);

  pc     = 0.4 * g_inputParam[RHOC] * pow(x1, -2.5);
  p_disk = eps2 * pow(coeff, 2.5);

  return (p_disk >= pc && rcyl > g_inputParam[RD]) ? 1 : 0;
}

/* ************************************************************** */
double DiskFraction (double *v, double x1, double x2)
/* cfa--replaces the tracer. f in [0,1] as product of two sigmoids:
   density position s between log rho_c(r) and log rho_d(r), midpoint
   1-1/CORONA_THRESH_FAC (=1/6) of the way, and |v_phi| against the
   Newtonian Keplerian speed R^-1/2, midpoint 0.5. Only dense and
   rotating gas counts as disc. */
{
  double rho_ref_c, rho_ref_d, log_span, position, arg_den, arg_rot;
  double den_factor, rot_factor, rcyl, v_kep, rot_frac;

  if (!g_profilesLive) {
    return (double) InitDiskCell(x1, x2);   // before the first update
  }

  if (!g_profilesInit) {
    rho_ref_c = AnalyticCoronaDensity(x1);   // should not happen
    rho_ref_d = AnalyticDiskDensity(x1);
    if (rho_ref_d <= 0.0) rho_ref_d = rho_ref_c;
  } else {
    rho_ref_c = GetCoronaRefDensity(x1);
    rho_ref_d = GetDiskRefDensity(x1);
  }

  log_span = log10(MAX(rho_ref_d, 1.e-30)) - log10(MAX(rho_ref_c, 1.e-30));
  log_span = MAX(log_span, DF_SPAN_MIN);   // inverted span cannot flip the criterion

  position = (log10(MAX(v[RHO], 1.e-30)) - log10(MAX(rho_ref_c, 1.e-30))) / log_span;

  arg_den = (position - (1.0 - 1.0/CORONA_THRESH_FAC)) / CORONA_SIGMOID_WIDTH;
  arg_den = MIN(MAX(arg_den, -50.0), 50.0);
  den_factor = 1.0 / (1.0 + exp(-arg_den));

  rcyl  = x1 * sin(x2);
  v_kep = 1.0 / sqrt(MAX(rcyl, 1.e-12));

  rot_frac = fabs(v[VX3]) / v_kep;

  arg_rot = (rot_frac - ROTATION_THRESH_FAC) / ROTATION_SIGMOID_WIDTH;
  arg_rot = MIN(MAX(arg_rot, -50.0), 50.0);
  rot_factor = 1.0 / (1.0 + exp(-arg_rot));

  return den_factor * rot_factor;
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*!
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
  g_radGrid = grid;   //cfa--grid for the opacities
  InitProfiles(grid);
}

/* **************************************************************** */
void Analysis (const Data *d, Grid *grid)
/*
 *
 * PURPOSE
 *
 *   Perform some pre-processing data
 *
 * ARGUMENTS
 *
 *   d:      the PLUTO Data structure.
 *   grid:   pointer to array of GRID structures
 *
 **************************************************************** */
{
  UpdateProfiles(d, grid);   //cfa--every step, see pluto.ini
}

#if PHYSICS == MHD
/* ************************************************************** */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*
 *
 * PURPOSE
 *
 *   Define the component of a static, curl-free background
 *   magnetic field.
 *
 **************************************************************** */
{
/* power law field for BH (Mishra: mm=-1.25, Zhu & Stone: -2.25)
  double mu   = g_inputParam[MU];
  double rmin = 1.0 * g_inputParam[RD];
  double mm   = -1.25;

  if (x1 <= rmin) {
    B0[0] = mu * cos(x2) * pow(rmin, mm) * (1.0 + sin(x2));
    B0[1] = -mu * sin(x2) * pow(rmin, mm);
  } else {
    B0[0] = mu * pow(x1 * sin(x2), mm) * cos(x2) * (1.0 + sin(x2));
    B0[1] = -mu * pow(x1 * sin(x2), mm) * sin(x2);
  }
  B0[2] = 0.0;
*/

/* dipole */
  B0[0] = 2.0 * g_inputParam[MU] * cos(x2) / (x1 * x1 * x1);
  B0[1] = g_inputParam[MU] * sin(x2) / (x1 * x1 * x1);
  B0[2] = 0.0;

/* quadrupole
  B0[0] = 1.5 * g_inputParam[MU] * (3.0 * cos(x2) * cos(x2) - 1.0) / (x1 * x1 * x1 * x1);
  B0[1] = 3.0 * g_inputParam[MU] * cos(x2) * sin(x2) / (x1 * x1 * x1 * x1);
  B0[2] = 0.0;
*/

/* octupole
  B0[0] = 2.0 * g_inputParam[MU] * (5.0 * cos(x2) * cos(x2) * cos(x2) - 3.0 * cos(x2)) / (x1 * x1 * x1 * x1 * x1);
  B0[1] = 0.5 * g_inputParam[MU] * (15.0 * cos(x2) * cos(x2) * sin(x2) - 3.0 * sin(x2)) / (x1 * x1 * x1 * x1 * x1);
  B0[2] = 0.0;
*/
}
#endif

/* ************************************************************** */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid)
/*
 *
 * PURPOSE
 *
 *   Define user-defined boundary conditions.
 *
 **************************************************************** */
{
  int i, j, k;
  double *x1, *x2, *x3, *r;
  double a1, a2, a, rcyl, eps2, coeff, lambda;
  double dvar1dr, dvar2dr, dvardr;
  double cs2, dden, dfact;

  RBox dom_box;

  RestoreDiskProfiles(d, grid);   //cfa--does something only once, after a restart

  if (side == 0) {
/* -- Impose conditions for the solution inside domain -- */
/*MC130421: changing a primitive variable in an active zone requires
recomputing the conservative variables in that zone as well, by the call
to PrimToCons3D() --see bottom of page 51 Userguide */
    x1 = grid->xgc[IDIR];
    x2 = grid->xgc[JDIR];
    x3 = grid->xgc[KDIR];

    TOT_LOOP(k, j, i) {
      int convert_to_cons = 0;

//ccm--density floor, keep cs and momentum. cfa--whole domain, not only near the star
      if (d->Vc[RHO][k][j][i] < g_inputParam[DFLOOR]) {
        dden = d->Vc[RHO][k][j][i];
        cs2  = g_gamma * d->Vc[PRS][k][j][i] / d->Vc[RHO][k][j][i];

        d->Vc[RHO][k][j][i] = g_inputParam[DFLOOR];
        dfact = dden / d->Vc[RHO][k][j][i];

        d->Vc[PRS][k][j][i] = cs2 * d->Vc[RHO][k][j][i] / g_gamma;
        d->Vc[VX1][k][j][i] *= dfact;
        d->Vc[VX2][k][j][i] *= dfact;
        d->Vc[VX3][k][j][i] *= dfact;

//ccm--ensure that TRC=0.0 in corona
        if (x2[j] < 0.5 * CONST_PI - atan(3.0 * g_inputParam[EPS]) ||
            x2[j] > 0.5 * CONST_PI + atan(3.0 * g_inputParam[EPS])) {
          d->Vc[TRC][k][j][i] = 0.0;
        }

        convert_to_cons = 1;
      }

#if RADIATION
//cfa--gas temperature floor, no recombination/molecules in the model.
// Without it the surface layers cooled to ~1e3 K and collapsed.
      if (GetTemperature(d->Vc[RHO][k][j][i], d->Vc[PRS][k][j][i]) < T_GAS_MIN) {
        d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i] * T_GAS_MIN / g_idealGasConst;
        convert_to_cons = 1;
      }
#endif

      if (convert_to_cons) {
        RBoxDefine(i, i, j, j, k, k, CENTER, &dom_box);
        PrimToCons3D(d->Vc, d->Uc, &dom_box, grid);
      }
    }
  }

  if (side == X1_BEG) {  /* -- X1_BEG boundary -- */
    if (box->vpos == CENTER) {

#if INNER_BOUNDARY == BOUNDARY_BH
/*############# BH: absorbing diode, cfa #####################*/
      X1_BEG_LOOP(k, j, i) {
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

//cfa--only inflow into the hole
        if (d->Vc[VX1][k][j][i] > 0.0) {
          d->Vc[VX1][k][j][i] = 0.0;
        }

        if (d->Vc[RHO][k][j][i] < g_inputParam[DFLOOR]) {
          d->Vc[RHO][k][j][i] = g_inputParam[DFLOOR];
        }
        if (d->Vc[PRS][k][j][i] < g_inputParam[DFLOOR] * 1.0e-3) {
          d->Vc[PRS][k][j][i] = g_inputParam[DFLOOR] * 1.0e-3;
        }

#if RADIATION
//cfa--radiation also only streams in
        d->Vc[ENR][k][j][i] = d->Vc[ENR][k][j][IBEG];
        d->Vc[FR1][k][j][i] = MIN(d->Vc[FR1][k][j][IBEG], 0.0);
        d->Vc[FR2][k][j][i] = d->Vc[FR2][k][j][IBEG];
        d->Vc[FR3][k][j][i] = d->Vc[FR3][k][j][IBEG];
#endif
      }

#elif INNER_BOUNDARY == BOUNDARY_STAR
/*############# rotating stellar surface #####################*/
      x2 = grid->x[JDIR];

      X1_BEG_LOOP(k, j, i) {
        d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][IBEG];
        d->Vc[PRS][k][j][i] = d->Vc[PRS][k][j][IBEG];
        d->Vc[TRC][k][j][i] = 0.0;

        d->Vc[VX1][k][j][i] = 0.0;
        d->Vc[VX2][k][j][i] = 0.0;
        d->Vc[VX3][k][j][i] = g_inputParam[OMEGA_STAR] * grid->x[IDIR][IBEG] * sin(x2[j]);

#if PHYSICS == MHD
        d->Vc[BX1][k][j][i] = d->Vc[BX1][k][j][IBEG];
        d->Vc[BX2][k][j][i] = d->Vc[BX2][k][j][IBEG];
        d->Vc[BX3][k][j][i] = 0.0;
#endif

        if (d->Vc[RHO][k][j][i] < g_inputParam[DFLOOR]) {
          d->Vc[RHO][k][j][i] = g_inputParam[DFLOOR];
        }

#if RADIATION
        d->Vc[ENR][k][j][i] = d->Vc[ENR][k][j][IBEG];   //cfa--zero gradient
        d->Vc[FR1][k][j][i] = d->Vc[FR1][k][j][IBEG];
        d->Vc[FR2][k][j][i] = d->Vc[FR2][k][j][IBEG];
        d->Vc[FR3][k][j][i] = d->Vc[FR3][k][j][IBEG];
#endif
      }
#endif

    }
  }

  if (side == X1_END) {  /* -- X1_END boundary -- */
    r  = grid->x[IDIR];
    x1 = grid->x[IDIR];
    x2 = grid->x[JDIR];

    if (box->vpos == CENTER) {
      BOX_LOOP(box, k, j, i) {
        d->Vc[TRC][k][j][i] = d->Vc[TRC][k][j][IEND];

#if RADIATION
//cfa--radiation leaves freely, nothing streams back in
        d->Vc[ENR][k][j][i] = d->Vc[ENR][k][j][IEND];
        d->Vc[FR1][k][j][i] = MAX(d->Vc[FR1][k][j][IEND], 0.0);
        d->Vc[FR2][k][j][i] = d->Vc[FR2][k][j][IEND];
        d->Vc[FR3][k][j][i] = d->Vc[FR3][k][j][IEND];
#endif

/*############# Logarithmic extrapolation of rho #####################*/
        a1 = log10(d->Vc[RHO][k][j][IEND]   / d->Vc[RHO][k][j][IEND-1]) / log10(r[IEND]   / r[IEND-1]);
        a2 = log10(d->Vc[RHO][k][j][IEND-1] / d->Vc[RHO][k][j][IEND-2]) / log10(r[IEND-1] / r[IEND-2]);
        a  = VANLEER_LIMITER(a1, a2);
        a  = MIN(a, 0.0);

        d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][i-1] * pow(r[i] / r[i-1], a);
        d->Vc[PRS][k][j][i] = d->Vc[PRS][k][j][IEND] * pow(d->Vc[RHO][k][j][i] / d->Vc[RHO][k][j][IEND], g_gamma);

        d->Vc[VX1][k][j][i] = d->Vc[VX1][k][j][IEND];
        d->Vc[VX2][k][j][i] = d->Vc[VX2][k][j][IEND];

#if PHYSICS == MHD
/*############# Toroidal magnetic field, VanLeer #################*/
        dvar1dr = (d->Vc[BX3][k][j][IEND]   - d->Vc[BX3][k][j][IEND-1]) / (r[IEND]   - r[IEND-1]);
        dvar2dr = (d->Vc[BX3][k][j][IEND-1] - d->Vc[BX3][k][j][IEND-2]) / (r[IEND-1] - r[IEND-2]);
        dvardr  = VANLEER_LIMITER(dvar1dr, dvar2dr);
        d->Vc[BX3][k][j][i] = d->Vc[BX3][k][j][i-1] + dvardr * (r[i] - r[i-1]);
#endif

/*############# vphi, MinMod #################*/
        dvar1dr = (d->Vc[VX3][k][j][IEND]   - d->Vc[VX3][k][j][IEND-1]) / (r[IEND]   - r[IEND-1]);
        dvar2dr = (d->Vc[VX3][k][j][IEND-1] - d->Vc[VX3][k][j][IEND-2]) / (r[IEND-1] - r[IEND-2]);
        dvardr  = MINMOD_LIMITER(dvar1dr, dvar2dr);
        d->Vc[VX3][k][j][i] = d->Vc[VX3][k][j][i-1] + dvardr * (r[i] - r[i-1]);

//ccm--inject the KK disc at the outer boundary, around the equator
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

//ccm--no backward inflow of the corona
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
/* ************************************************************************ */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{
  g[IDIR] = -1.0 / (x1 * x1);
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}

/* ************************************************************************ */
double BodyForcePotential(double x1, double x2, double x3)
/*
 *
 *
 *
 *************************************************************************** */
{
//  return -1.0 / x1;   /* Newtonian */

  return -1.0 / (x1 - 2.0);   //ccm--Paczynski-Wiita

//  return -(1.0 / 6.0) * (exp(6.0 / x1) - 1.0);   //csk/ccm--Kluzniak-Lee

//  double q = 1.25;   /* Kluzniak-Nordstrom */
//  return -1.0 / x1 + 0.5 * (q * q / (x1 * x1));
}
#endif

#if RADIATION_VAR_OPACITIES
/* cfa--Kramers ff+bf and Thomson in cm^2/g,
     kappa_ff = C_FF g_ff (1-Z)(1+X) rho T^-3.5
     kappa_bf = C_BF g_bf Z(1+X)/t   rho T^-3.5
     kappa_es = 0.2 (1+X)                        */
#define G_BF     1.0
#define G_FF     1.0
#define T_FACTOR 1.0e1     // guillotine factor

#define C_BF 4.34e25
#define C_FF 3.68e22

static const double K_BF = C_BF * G_BF * Z_MASSFRAC * (1.0 + X_MASSFRAC) / T_FACTOR;
static const double K_FF = C_FF * G_FF * (1.0 - Z_MASSFRAC) * (1.0 + X_MASSFRAC);
static const double K_ES = 0.2 * (1.0 + X_MASSFRAC);

/* ************************************************************** */
static double ComptonOpacity (double *v, double kappa_es)
/* cfa--Compton heating/cooling as an effective absorption opacity.
   Scattering is coherent in PLUTO, so it never exchanges energy.
   With
     kappa_c = kappa_es 4kT_r^4 / (m_e c^2 (T_g+T_r)(T_g^2+T_r^2))
   the absorption term rho kappa_c c (a_R T_g^4 - E_r) is exactly the
   nonrel. rate rho kappa_es c E_r 4k(T_g-T_r)/(m_e c^2). Needs ENR to
   be current in the implicit step, see rad_step.c.
   Above T_COMPTON_MAX the rate is held at its value there (corona is
   at ~1e11 K, far outside the validity), capping T_g inside the first
   form would let it grow as T_g^4 instead. */
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
              / (mc2 * (Tg*Tg*Tg*Tg - Tr4));
  } else {
    kappa_c = 0.0;
  }

  return MIN(kappa_c, kappa_es);   // scattering stays >= 0
}

/* ************************************************************** */
void UserDefOpacitiesAt(double *v, double x1, double x2, double *abs, double *scat)
/* cfa--opacities per unit mass in code units, with the position given
   explicitly. Also used for the kappa output in userdef_output.c. */
{
  double rho_cgs, T, kappa_ffbf, kappa_es, kappa_c, fd, gv, beta;

  rho_cgs = v[RHO] * UNIT_DENSITY;
  T       = MAX(GetTemperature(v[RHO], v[PRS]), T_OPAC_MIN);

  kappa_ffbf = (K_BF + K_FF) * rho_cgs * pow(T, -3.5);
  kappa_es   = K_ES;
  kappa_c    = ComptonOpacity(v, kappa_es);

/* decouple matter and radiation above ~0.9c. Plunging gas moves at
   several c in the PW potential, and the v^2/c^2 terms of the nonrel.
   module then blow up. Orbital speed is already 0.61c at R=6, so the
   threshold has to stay well above that. */
  beta = sqrt(v[VX1]*v[VX1] + v[VX2]*v[VX2] + v[VX3]*v[VX3]) / g_radC;
  gv   = 1.0 / (1.0 + exp(MIN((beta - BETA_DECOUPLE) / BETA_WIDTH, 50.0)));

/* energy exchange only in the disc, off below F_OPAC_MIN (the corona
   keeps f~1e-3 and was Compton cooled away otherwise) */
  fd = DiskFraction(v, x1, x2);
  fd = (fd > F_OPAC_MIN) ? (fd - F_OPAC_MIN) / (1.0 - F_OPAC_MIN) : 0.0;

/* Compton moved from scattering to absorption, total stays Thomson.
   Thomson itself acts on all gas (gating it gave a dense RT-like
   sheet at the disc surface). */
  *abs  = gv * fd * (kappa_ffbf + kappa_c)      * UNIT_DENSITY * UNIT_LENGTH;
  *scat = gv * (kappa_es - fd * kappa_c)        * UNIT_DENSITY * UNIT_LENGTH;
}

/* ************************************************************** */
void UserDefOpacities(double *v, double *abs, double *scat)
/* called by the radiation module, (x1,x2) from g_i_rad and g_j */
{
    double x1 = g_radGrid->x[IDIR][g_i_rad];
    double x2 = g_radGrid->x[JDIR][g_j];

    UserDefOpacitiesAt(v, x1, x2, abs, scat);
}
#endif
