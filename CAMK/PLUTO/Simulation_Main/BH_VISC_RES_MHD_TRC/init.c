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
 * START: Lagrangian Tracer Replacement from Physical Properties
 * --------------------------------------------------------------------- */

/* Ensure iVPHI points to azimuthal velocity component in spherical geometry */
#ifndef iVPHI
  #define iVPHI VX3
#endif

/* Global persistent arrays for cached lookup */
static double ***disk_frac_cache = NULL;  /* EMA-relaxed value actually used by physics */
static double **raw_mask         = NULL;  /* instantaneous, hard-thresholded classification */
static double **smooth_mask      = NULL;  /* spatially blurred version of raw_mask (this step) */
static double *local_rho_buf     = NULL;
static double *global_rho_buf    = NULL;

/* Set to 1 once InitDomain has primed disk_frac_cache with a first value,
 * so the very first Analysis() call blends against a real field instead
 * of zero-initialized memory (which would otherwise be seen as "all corona"
 * and instantaneously kill dissipation everywhere on step 1). */
static int disk_frac_primed = 0;

/* -------------------------------------------------------------------
 * DESIGN NOTE (read before touching this file again):
 *
 * The original version of this cache recomputed a HARD boolean
 * (rot_flag && rho_flag) every call to Analysis() and fed it straight
 * into disk_frac_cache, which in turn gates viscosity/opacity. That
 * is a discontinuous function of the solution feeding back into the
 * coefficients of the PDE being solved every step: a cell straddling
 * the threshold flips 0<->1 from one step to the next, dissipation
 * switches on/off like a delta-function forcing term, and the
 * resulting shock drives local density toward DFLOOR, which can push
 * neighbouring cells across their own threshold on the next step.
 * That is a real, first-order stability defect, not a corner case.
 *
 * Fix implemented below: disk_frac_cache is an EXPONENTIAL MOVING
 * AVERAGE (EMA) of the hard/blurred classification, updated once per
 * call with a fixed relaxation rate EMA_ALPHA. The physics-facing
 * value can now only change by a bounded amount per call, so it can
 * no longer act as an impulsive forcing term. Do not read raw_mask or
 * smooth_mask directly from viscosity/opacity code -- always go
 * through DiskFraction(), which returns the relaxed value.
 * ------------------------------------------------------------------- */
#define EMA_ALPHA 0.1   /* relaxation rate; ~1/EMA_ALPHA calls to reach steady state */

/* MPI Gathering Buffers */
#ifdef PARALLEL
static MPI_Comm theta_comm = MPI_COMM_NULL;
static int theta_comm_init = 0;
static int *recv_counts    = NULL;
static int *disps          = NULL;
#endif

/* -------------------------------------------------------------------
 * Helper: Double comparison for quicksort median calculation
 * ------------------------------------------------------------------- */
static int CompareDoubles (const void *a, const void *b)
{
  double da = *(const double *)a;
  double db = *(const double *)b;
  return (da > db) - (da < db);
}

/* -------------------------------------------------------------------
 * Helper: Index clamping for boundary domain padding
 * ------------------------------------------------------------------- */
static inline int ClampIndex (int val, int min_val, int max_val)
{
  if (val < min_val) return min_val;
  if (val > max_val) return max_val;
  return val;
}

/* -------------------------------------------------------------------
 * STEP 2: Smooth Median Generator & Rotation Veto (MPI-Safe Timestep Hook)
 * ------------------------------------------------------------------- */
void UpdateDiskFractionCache (const Data *d, Grid *grid)
{
  int i, j, k;

  /* 1. Persistent memory allocations */
  if (disk_frac_cache == NULL) {
    disk_frac_cache = ARRAY_3D(NX3_TOT, NX2_TOT, NX1_TOT, double);
  }
  if (raw_mask == NULL) {
    raw_mask = ARRAY_2D(NX2_TOT, NX1_TOT, double);
  }
  if (smooth_mask == NULL) {
    smooth_mask = ARRAY_2D(NX2_TOT, NX1_TOT, double);
  }
  if (local_rho_buf == NULL) {
    local_rho_buf = ARRAY_1D(NX2_TOT, double);
  }

  /* 2. Setup MPI Sub-Communicator along Theta (JDIR) for each Radial (IDIR) Slice */
  int global_n2 = NX2;
#ifdef PARALLEL
  if (!theta_comm_init) {
    /* Group ranks sharing the same radial rank index (g_domBeg[IDIR]) */
    MPI_Comm_split(MPI_COMM_WORLD, g_domBeg[IDIR], prank, &theta_comm);
    theta_comm_init = 1;
  }

  int theta_proc_size;
  MPI_Comm_size(theta_comm, &theta_proc_size);

  if (recv_counts == NULL) {
    recv_counts = (int *)malloc(theta_proc_size * sizeof(int));
    disps       = (int *)malloc(theta_proc_size * sizeof(int));
  }

  int local_n2 = NX2; /* Local active theta cells */
  MPI_Allgather(&local_n2, 1, MPI_INT, recv_counts, 1, MPI_INT, theta_comm);

  disps[0] = 0;
  for (int r = 1; r < theta_proc_size; r++) {
    disps[r] = disps[r - 1] + recv_counts[r - 1];
  }
  global_n2 = disps[theta_proc_size - 1] + recv_counts[theta_proc_size - 1];

  if (global_rho_buf == NULL) {
    global_rho_buf = ARRAY_1D(global_n2, double);
  }
#else
  global_rho_buf = local_rho_buf;
#endif

  /* Classification parameters */
  const double median_factor = 5.0;
  const double vphi_frac     = 0.5;
  const int blur_radius      = 1; /* 3x3 stencil window (radius=2 gives 5x5, not 3x3 -- fixed) */

  /* -----------------------------------------------------------------
   * Part A: Compute Global Radial Medians and Apply Rotation Veto
   * ----------------------------------------------------------------- */
  for (i = 0; i < NX1_TOT; i++) {
    double R = grid->x[IDIR][i];

    /* Extract local polar density slice for radial cell i, azimuthally
     * averaged over all local KBEG..KEND planes. Sampling KBEG alone
     * assumes axisymmetry; in a 3D run with turbulence or warping,
     * a single phi-plane is not representative of the shell and will
     * misclassify material at other azimuths. */
    int local_idx = 0;
    for (j = JBEG; j <= JEND; j++) {
      double rho_sum = 0.0;
      for (k = KBEG; k <= KEND; k++) {
        rho_sum += d->Vc[RHO][k][j][i];
      }
      local_rho_buf[local_idx++] = rho_sum / (double)(KEND - KBEG + 1);
    }

    /* Gather full polar density vector across theta sub-communicator */
#ifdef PARALLEL
    MPI_Allgatherv(local_rho_buf, NX2, MPI_DOUBLE,
                   global_rho_buf, recv_counts, disps, MPI_DOUBLE, theta_comm);
#endif

    /* Compute median density across global theta span */
    qsort(global_rho_buf, global_n2, sizeof(double), CompareDoubles);
    double med_rho;
    if (global_n2 % 2 == 0) {
      med_rho = 0.5 * (global_rho_buf[global_n2 / 2 - 1] + global_rho_buf[global_n2 / 2]);
    } else {
      med_rho = global_rho_buf[global_n2 / 2];
    }

    /* Apply thresholding and rotation veto across active and ghost cells.
     * Both rho and vphi are azimuthally averaged over local phi planes
     * for the same reason as above -- consistency with the median that
     * was computed from azimuthally-averaged density. */
    for (j = 0; j < NX2_TOT; j++) {
      double theta = grid->x[JDIR][j];
      double Rcyl  = R * sin(theta);

      double vK = (Rcyl > 1.0e-12) ? (1.0 / sqrt(Rcyl)) : 0.0;

      double rho_avg = 0.0, vphi_avg = 0.0;
      for (k = KBEG; k <= KEND; k++) {
        rho_avg  += d->Vc[RHO][k][j][i];
        vphi_avg += d->Vc[iVPHI][k][j][i];
      }
      rho_avg  /= (double)(KEND - KBEG + 1);
      vphi_avg /= (double)(KEND - KBEG + 1);

      int rot_flag = (vphi_avg > vphi_frac * vK);
      int rho_flag = (rho_avg > median_factor * med_rho);

      if (rot_flag && rho_flag) {
        raw_mask[j][i] = 1.0;
      } else {
        raw_mask[j][i] = 0.0;
      }
    }
  }

  /* -----------------------------------------------------------------
   * Part B: Spatial 2D Moving-Box Blur (3x3 stencil, blur_radius=1)
   * ----------------------------------------------------------------- */
  for (j = JBEG; j <= JEND; j++) {
    for (i = IBEG; i <= IEND; i++) {
      double sum = 0.0;
      int count = 0;

      for (int dj = -blur_radius; dj <= blur_radius; dj++) {
        for (int di = -blur_radius; di <= blur_radius; di++) {
          int sample_i = ClampIndex(i + di, 0, NX1_TOT - 1);
          int sample_j = ClampIndex(j + dj, 0, NX2_TOT - 1);

          sum += raw_mask[sample_j][sample_i];
          count++;
        }
      }

      smooth_mask[j][i] = sum / (double)count;
    }
  }

  /* -----------------------------------------------------------------
   * Part C: Temporal Relaxation (EMA) -- THE critical stability fix.
   *
   * disk_frac_cache is never overwritten with the instantaneous
   * smooth_mask. Instead it is relaxed toward it:
   *
   *   F^{n+1} = alpha * smooth_mask + (1 - alpha) * F^n
   *
   * so viscosity/opacity gated by DiskFraction() can only change by a
   * bounded amount (~EMA_ALPHA) per call, never jump discontinuously.
   * On the very first call (t=0, from InitDomain), there is no
   * meaningful F^n yet, so we seed disk_frac_cache directly from
   * smooth_mask instead of blending against zero-initialized memory
   * (blending against zero would read as "instant corona everywhere"
   * and kill dissipation on step 1).
   * ----------------------------------------------------------------- */
  for (j = JBEG; j <= JEND; j++) {
    for (i = IBEG; i <= IEND; i++) {
      double relaxed;
      if (!disk_frac_primed) {
        relaxed = smooth_mask[j][i];
      } else {
        relaxed = EMA_ALPHA * smooth_mask[j][i]
                + (1.0 - EMA_ALPHA) * disk_frac_cache[KBEG][j][i];
      }
      for (k = KBEG; k <= KEND; k++) {
        disk_frac_cache[k][j][i] = relaxed;
      }
    }
  }
  disk_frac_primed = 1;

  /* -----------------------------------------------------------------
   * Part D: Pad Ghost Zones for Boundary Queries
   * ----------------------------------------------------------------- */
  for (k = 0; k < NX3_TOT; k++) {
    int clamped_k = ClampIndex(k, KBEG, KEND);
    for (j = 0; j < NX2_TOT; j++) {
      int clamped_j = ClampIndex(j, JBEG, JEND);
      for (i = 0; i < NX1_TOT; i++) {
        int clamped_i = ClampIndex(i, IBEG, IEND);
        disk_frac_cache[k][j][i] = disk_frac_cache[clamped_k][clamped_j][clamped_i];
      }
    }
  }
}

/* -------------------------------------------------------------------
 * STEP 3: Fast O(1) Lookup Function
 * ------------------------------------------------------------------- */
double DiskFraction (double *v, double x1, double x2)
{
  if (disk_frac_cache != NULL) {
    int i = ClampIndex(g_i, 0, NX1_TOT - 1);
    int j = ClampIndex(g_j, 0, NX2_TOT - 1);
    int k = ClampIndex(g_k, 0, NX3_TOT - 1);
    return disk_frac_cache[k][j][i];
  }
  return 0.0;
}

/* -------------------------------------------------------------------
 * STEP 1: Integration Hooks inside init.c
 * ------------------------------------------------------------------- */

/* Called BEFORE time-stepping starts (Cold-Start Guard at t=0) */
void InitDomain (Data *d, Grid *grid)
{
  /* 1. Force immediate boundary fill and MPI exchange for initial condition Vc */
  Boundary(d, X1_BEG, grid);
  Boundary(d, X1_END, grid);
  Boundary(d, X2_BEG, grid);
  Boundary(d, X2_END, grid);
#if DIMENSIONS == 3
  Boundary(d, X3_BEG, grid);
  Boundary(d, X3_END, grid);
#endif

  /* 2. Safely populate disk fraction lookup cache using valid ghost data */
  UpdateDiskFractionCache(d, grid);
}

/* Called at the end of every timestep.
 *
 * KNOWN LIMITATION -- read before relying on this for fast disk
 * evolution: PLUTO does not expose a user hook that runs before the
 * first RK/CT sub-stage of a step, only Analysis() (end of step) and
 * InitDomain() (t=0 only). That means disk_frac_cache, and therefore
 * every viscosity/opacity call that reads DiskFraction(), is always
 * one full step stale relative to the sub-stages currently being
 * integrated -- Vc has already advanced by the time this runs.
 *
 * With EMA_ALPHA ~ 0.08 the cache changes slowly enough per call that
 * one step of staleness is not itself destabilizing (it was the
 * *instantaneous hard* reclassification, not the staleness, that
 * caused the original explosion). But if the disk boundary is
 * expected to move on a timescale comparable to a few timesteps, this
 * lag will visibly trail the true disk/corona interface. If that
 * becomes a problem, the correct fix is not to shrink EMA_ALPHA (that
 * reintroduces the impulsive-forcing failure mode) but to update the
 * cache from inside a lower-level hook that runs pre-stage, e.g.
 * wrapping RightHandSide() or the viscosity source term itself --
 * confirm the exact call graph against the PLUTO version in use
 * before doing that; do not guess at an API that isn't verified. */
void Analysis (const Data *d, Grid *grid)
{
  UpdateDiskFractionCache(d, grid);
}

/* ---------------------------------------------------------------------
 * END: Lagrangian Tracer Replacement from Physical Properties
 * --------------------------------------------------------------------- */

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
    v[TRC] = 1.0;     /* Disk tracer (kept for diagnostics/comparison; no
                          longer used to gate viscosity/resistivity) */
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
