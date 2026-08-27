#ifndef MODIFICATIONS_H
#define MODIFICATIONS_H

/* ---------------------------------------------------------------------------------
   Tunables for the adaptive disk versus corona classifier.
   Kept as compile-time constants to avoid touching "definitions.h" enumeration.
   --------------------------------------------------------------------------------- */

#define NBINS_CORONA          90       /* number of logarithmic radial bins */
#define CORONA_AVG_FREQ       10       /* update every number of steps */
#define CORONA_EMA_ALPHA       0.05    /* temporal smoothing weight */
#define CORONA_THRESH_FAC      1.05    /* transition midpoint, in terms of
                                           position along the corona->disk
                                           log-density span (see DiskFraction) */
#define CORONA_SIGMOID_WIDTH   0.25    /* transition width, same units */
#define DISK_BIN_VOLFRAC_MIN   1.e-3   /* minimum disk-weighted volume
                                           fraction of a bin's total volume
                                           this step, below which the bin
                                           is NOT considered to have "real"
                                           disk material (guards against a
                                           single transient infalling cell
                                           validating an ISCO-side bin) */

double DiskFraction (double *v, double x1, double x2);

#endif /* MODIFICATIONS_H */
