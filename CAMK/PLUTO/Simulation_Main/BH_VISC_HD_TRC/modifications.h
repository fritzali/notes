#ifndef MODIFICATIONS_H
#define MODIFICATIONS_H

/* ---------------------------------------------------------------------------------
   Tunables for the adaptive disk versus corona classifier.
   Kept as compile-time constants to avoid touching "definitions.h" enumeration.
   --------------------------------------------------------------------------------- */

#define NBINS_CORONA          90       /* number of logarithmic radial bins */
#define CORONA_AVG_FREQ       10       /* update every number of steps */
#define CORONA_EMA_ALPHA       0.05    /* temporal smoothing weight */
#define CORONA_THRESH_FAC      1.5     /* disk to corona density ratio threshold */
#define CORONA_SIGMOID_WIDTH   0.5     /* transition width in logarithmic units */

double DiskFraction (double *v, double x1, double x2);

#endif /* MODIFICATIONS_H */
