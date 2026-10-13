#ifndef MODIFICATIONS_H
#define MODIFICATIONS_H

/* cfa--DiskFraction tunables, see init.c */

#define NBINS_PROFILE         90      // log radial bins of the reference profiles
#define DF_EMA_TAU             6.0    // smoothing time constant of the profiles, in t_g

#define DF_SPAN_MIN            0.5    // min. log10(rho_d/rho_c) span, in dex

/* reference profiles: SELF from the classification itself (default),
   GEOM from bin max. density and corona mean above DF_GEOM_ZCOR */
#define DF_REF_SELF            0
#define DF_REF_GEOM            1
#define DISKFRAC_REFERENCE     DF_REF_SELF
#define DF_GEOM_ZCOR           0.5    // z/R

#define CORONA_THRESH_FAC      1.2    // density midpoint at 1-1/1.2 of the log span

#define CORONA_SIGMOID_WIDTH   0.3

#define DISK_BIN_VOLFRAC_MIN   1.e-3  // min. disc volume fraction for a bin to count

#define ROTATION_THRESH_FAC    0.5    // midpoint at v_phi = 0.5 v_kep
#define ROTATION_SIGMOID_WIDTH 0.1

double DiskFraction (double *v, double x1, double x2);
void   SaveDiskProfiles (void);

#if RADIATION_VAR_OPACITIES
void UserDefOpacitiesAt(double *v, double x1, double x2, double *abs, double *scat);
#endif

#endif
