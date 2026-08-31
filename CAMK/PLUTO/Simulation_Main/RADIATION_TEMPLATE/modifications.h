#ifndef MODIFICATIONS_H
#define MODIFICATIONS_H

/* ---------------------------------------------------------------------------------
   Tunables for the adaptive disk versus corona classifier.
   Kept as compile time constants to avoid touching definitions enumeration.
   --------------------------------------------------------------------------------- */

#define NBINS_PROFILE         90        /* number of logarithmic radial bins */
#define CORONA_AVG_FREQ       10        /* update every number of steps */
#define CORONA_EMA_ALPHA       0.05     /* temporal smoothing weight */
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

#if RADIATION_VAR_OPACITIES
void UserDefOpacitiesAt(double *v, double x1, double x2, double *abs, double *scat);
#endif

#endif /* MODIFICATIONS_H */
