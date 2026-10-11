#ifndef MODIFICATIONS_H
#define MODIFICATIONS_H

/* ---------------------------------------------------------------------------------
   Tunables for the adaptive disk versus corona classifier.
   Kept as compile time constants to avoid touching definitions enumeration.
   --------------------------------------------------------------------------------- */

#define NBINS_PROFILE         90        /* number of logarithmic radial bins */
#define DF_EMA_TAU             6.0      /* time constant of the temporal smoothing
                                           of both reference profiles, in t_g */
#define DF_SPAN_MIN            0.5      /* minimum logarithmic density span between
                                           disk and corona reference, in dex */

/* Reference profiles of the density criterion, chosen at compile time:
     DF_REF_SELF  disk and corona weighted bin means, built from the current
                  classification itself
     DF_REF_GEOM  bin maximum density for the disk and geometric mean density
                  above z/R = DF_GEOM_ZCOR for the corona, independent of the
                  classification
   The classification of each cell works the same way in both modes. */
#define DF_REF_SELF            0
#define DF_REF_GEOM            1
#define DISKFRAC_REFERENCE     DF_REF_SELF
#define DF_GEOM_ZCOR           0.5      /* z/R above which cells define the corona
                                           reference in the DF_REF_GEOM mode */

#define CORONA_THRESH_FAC      1.2      /* transition midpoint, in terms of
                                           position along the corona to disk
                                           logarithmic density span */
#define CORONA_SIGMOID_WIDTH   0.3      /* transition width, , in terms of
                                           position along the corona to disk
                                           logarithmic density span */
#define DISK_BIN_VOLFRAC_MIN   1.e-3    /* minimum disk weighted volume
                                           fraction of total bin volume
                                           this step, below which the bin
                                           is not considered to have real
                                           disk material */

#define ROTATION_THRESH_FAC    0.5      /* midpoint fraction of local Keplerian velocity (0.5 means v_phi = 0.5 * v_kep) */
#define ROTATION_SIGMOID_WIDTH 0.1      /* transition width for the rotation sigmoid */

double DiskFraction (double *v, double x1, double x2);
void   SaveDiskProfiles (void);

#if RADIATION_VAR_OPACITIES
void UserDefOpacitiesAt(double *v, double x1, double x2, double *abs, double *scat);
#endif

#endif /* MODIFICATIONS_H */
