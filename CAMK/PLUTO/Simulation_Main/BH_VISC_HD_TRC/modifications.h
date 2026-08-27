#ifndef MODIFICATIONS_H
#define MODIFICATIONS_H

/* ---------------------------------------------------------------------------------
   Tunables for the adaptive disk versus corona classifier.
   Kept as compile-time constants to avoid touching "definitions.h" enumeration.
   --------------------------------------------------------------------------------- */

#define NBINS_PROFILE          90       /* number of logarithmic radial bins */
#define CORONA_AVG_FREQ       10        /* update every number of steps */
#define CORONA_EMA_ALPHA       0.05     /* temporal smoothing weight */
#define CORONA_THRESH_FAC      1.2      /* transition midpoint, in terms of
                                           position along the corona->disk
                                           log-density span (see DiskFraction) */
#define CORONA_SIGMOID_WIDTH   0.3      /* transition width, same units */
#define DISK_BIN_VOLFRAC_MIN   1.e-4    /* minimum disk-weighted volume
                                           fraction of a bin's total volume
                                           this step, below which the bin
                                           is NOT considered to have "real"
                                           disk material (guards against a
                                           single transient infalling cell
                                           validating an ISCO-side bin) */

#define ROTATION_THRESH_FAC    0.5      /* Midpoint fraction of local Keplerian velocity (e.g. 0.5 means v_phi = 0.5 * v_kep) */
#define ROTATION_SIGMOID_WIDTH 0.1      /* Transition width for the rotation sigmoid */

double DiskFraction (double *v, double x1, double x2);

#endif /* MODIFICATIONS_H */
