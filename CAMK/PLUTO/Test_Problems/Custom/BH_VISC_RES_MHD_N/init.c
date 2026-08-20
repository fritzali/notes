/* /////////////////////////////////////////////////////////////////////////// */
/*!
  \file
  \brief Contains basic functions for problem initialization.

  The init.c file collects most of the user-supplied functions useful
  for problem configuration.
  It is automatically searched for by the makefile.

  \author A. Mignone (mignone@ph.unito.it)
  \date Sep 2012

  \modified M. Cemeljic (miki@camk.edu.pl)
  \date Jul 2020

  Stripped minimal version of ideal MHD setup for teaching purposes.
  Based on appendix of "Atlas" paper, Cemeljic, 2019, A&A, 624, A31
*/
/* /////////////////////////////////////////////////////////////////////////// */
#include "pluto.h"

/* ************************************************************** */
void Init (double *v, double x1, double x2, double x3)
/*
 *
 **************************************************************** */
{
  double coeff, eps2, pc, rcyl;
  double br,bth;
  double lambda;
  double xhi2, Rco;

  rcyl=x1*sin(x2);
  eps2=g_inputParam[EPS]*g_inputParam[EPS];
  coeff=2./5./eps2*(1./x1-(1.-5./2.*eps2)/rcyl);
  lambda=11./5./(1.+64./25.*g_inputParam[ALPHAV]*g_inputParam[ALPHAV]);

/* initial non-rotating adiabatic corona in hydrostatic equilibrium  */
  v[RHO] = g_inputParam[RHOC]*pow(x1,-3./2.);
  v[PRS] = 2./5.*g_inputParam[RHOC]*pow(x1,-5./2.);

  pc=v[PRS];

  v[VX1] = 0.0;
  v[VX2] = 0.0;
  v[VX3] = 0.0;

/* Keplerian adiabatic disk in vertical pressure equilibrium with the
   adiabatic corona, as given by Kluzniak & Kita (2000) */

  v[PRS]=eps2*pow(coeff,5./2.);

    if (v[PRS] >= pc && rcyl > g_inputParam[RD])
      {v[RHO] = pow(coeff,3./2.);
       v[VX1] = -g_inputParam[ALPHAV]/sin(x2)*eps2*(10.-32./3.
       *lambda*g_inputParam[ALPHAV]*g_inputParam[ALPHAV]
       -lambda*(5.-1./(eps2*tan(x2)*tan(x2))))/sqrt(rcyl);
       v[VX3] = (sqrt(1.-5./2.*eps2)+2./3.*eps2
       *g_inputParam[ALPHAV]*g_inputParam[ALPHAV]
       *lambda*(1.-6./(5.*eps2*tan(x2)*tan(x2))))/sqrt(rcyl);
       v[TRC] = 1.0;     /* Track the disc material */
       }
    else
      {
       v[PRS]=2./5.*g_inputParam[RHOC]*pow(x1,-5./2.);
       v[TRC] = 0.0;     /* Track the corona */
      }
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
{}


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
{}


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
 *
 * ARGUMENTS
 *
 *   x1, x2, x3  (IN)    coordinates
 *
 *   B0         (OUT)    vector component of the background field.
 *
 *
 **************************************************************** */
{
/* BH case, Zhu&Stone (2018), Mishra et al.(2019) 
   MC Aug 2019    */

// double mu,rmin,mm;
// 
// mu=g_inputParam[MU];
// rmin=1.*g_inputParam[RD];
// 
// mm=-5./4.;//5./4.;//BMishra -5/4, Zhu&Stone -9/4
// 
// if(x1<=rmin){
// B0[0]=mu*cos(x2)*pow(rmin,mm)*(1.+sin(x2));
// B0[1]=-mu*sin(x2)*pow(rmin,mm);
// }else{
// B0[0]=mu*pow(x1*sin(x2),mm)*cos(x2)*(1.+sin(x2));
// B0[1]=-mu*pow(x1*sin(x2),mm)*sin(x2);
// }
// B0[2] = 0.0; 


/* dipole */
   B0[0] = 2.*g_inputParam[MU]*cos(x2)/(x1*x1*x1);
   B0[1] = g_inputParam[MU]*sin(x2)/(x1*x1*x1);
   B0[2] = 0.0;                             
/* */

/* quadrupole 
  B0[0] = 3.0/2.0*g_inputParam[MU]*(3.0*cos(x2)*cos(x2)-1.0)/(x1*x1*x1*x1);
  B0[1] = 3.0*g_inputParam[MU]*cos(x2)*sin(x2)/(x1*x1*x1*x1);
  B0[2] = 0.0;
*/       

/* octupole 
  B0[0] = 2.0*g_inputParam[MU]*(5.0*cos(x2)*cos(x2)*cos(x2)-3.0*cos(x2))/(x1*x1*x1*x1*x1);
  B0[1] = 0.5*g_inputParam[MU]*(15.0*cos(x2)*cos(x2)*sin(x2)-3.0*sin(x2))/(x1*x1*x1*x1*x1);
  B0[2] = 0.0;
 */         
}
#endif


/* ************************************************************** */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid)
/* 
 * PURPOSE:
 *   User-defined boundary conditions for accretion disk around 
 *   a non-rotating black hole in Paczyński-Wiita potential.
 *
 **************************************************************** */
{
  int i, j, k;
  double *x1, *x2, *x3, *r;
  double a1, a2, a, rcyl, eps2, coeff, lambda;
  double dvar1dr, dvar2dr, dvardr;
  double cs2, dden, dfact;
  
  RBox dom_box;

  /* -----------------------------------------------------------------
     side == 0: Interior active domain checks (density floor & entropy)
     ----------------------------------------------------------------- */
  if (side == 0) {    
    x1 = grid->xgc[IDIR];
    x2 = grid->xgc[JDIR];
    x3 = grid->xgc[KDIR];

    TOT_LOOP(k,j,i) {
      int convert_to_cons = 0;

      /* Density floor check inside domain */
      if (d->Vc[RHO][k][j][i] < g_inputParam[DFLOOR]) {
        dden = d->Vc[RHO][k][j][i];
        cs2  = g_gamma * d->Vc[PRS][k][j][i] / d->Vc[RHO][k][j][i];

        /* Reset density to floor value */
        d->Vc[RHO][k][j][i] = g_inputParam[DFLOOR];     
        dfact = dden / d->Vc[RHO][k][j][i];

        /* Modify pressure to preserve local sound speed */
        d->Vc[PRS][k][j][i] = cs2 * d->Vc[RHO][k][j][i] / g_gamma;

        /* Rescale velocities to conserve momentum */
        d->Vc[VX1][k][j][i] *= dfact;
        d->Vc[VX2][k][j][i] *= dfact;
        d->Vc[VX3][k][j][i] *= dfact;

        /* Clear tracer in coronal region to prevent numerical artifacts */
        if (x2[j] < 0.5 * CONST_PI - atan(3.0 * g_inputParam[EPS]) ||
            x2[j] > 0.5 * CONST_PI + atan(3.0 * g_inputParam[EPS])) {
          d->Vc[TRC][k][j][i] = 0.0;
        }

        convert_to_cons = 1;
      }

      /* Recompute conservative variables if primitives were modified */
      if (convert_to_cons) {
        RBoxDefine(i, i, j, j, k, k, CENTER, &dom_box);
        PrimToCons3D(d->Vc, d->Uc, &dom_box, grid);
      }
    }
  }

  /* -----------------------------------------------------------------
     side == X1_BEG: Event Horizon Absorbing Boundary (r = 2.1)
     ----------------------------------------------------------------- */
  if (side == X1_BEG) {
    if (box->vpos == CENTER) {
      X1_BEG_LOOP(k,j,i) {
        /* Standard zero-gradient copy from the first active cell (IBEG) */
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

        /* DIODE CONDITION: Allow inflow into the black hole (v_r <= 0),
           strictly prevent material or waves from re-entering grid (v_r > 0) */
        if (d->Vc[VX1][k][j][i] > 0.0) {
          d->Vc[VX1][k][j][i] = 0.0;
        }

        /* Floor enforcement in ghost cells */
        if (d->Vc[RHO][k][j][i] < g_inputParam[DFLOOR]) {
          d->Vc[RHO][k][j][i] = g_inputParam[DFLOOR];
        }
        if (d->Vc[PRS][k][j][i] < g_inputParam[DFLOOR] * 1.0e-3) {
          d->Vc[PRS][k][j][i] = g_inputParam[DFLOOR] * 1.0e-3;
        }
      }
    }
  }

  /* -----------------------------------------------------------------
     side == X1_END: Outer Radial Boundary (Disk Feeding & Corona)
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

        /* Outflow condition for poloidal velocities */
        d->Vc[VX1][k][j][i] = d->Vc[VX1][k][j][IEND]; 
        d->Vc[VX2][k][j][i] = d->Vc[VX2][k][j][IEND]; 

#if PHYSICS == MHD
        /* Van Leer extrapolation for toroidal magnetic field */
        dvar1dr = (d->Vc[BX3][k][j][IEND]   - d->Vc[BX3][k][j][IEND-1]) / (r[IEND]   - r[IEND-1]);
        dvar2dr = (d->Vc[BX3][k][j][IEND-1] - d->Vc[BX3][k][j][IEND-2]) / (r[IEND-1] - r[IEND-2]);
        dvardr  = VANLEER_LIMITER(dvar1dr, dvar2dr);
        d->Vc[BX3][k][j][i] = d->Vc[BX3][k][j][i-1] + dvardr * (r[i] - r[i-1]);
#endif

        /* Linear extrapolation of azimuthal velocity */
        dvar1dr = (d->Vc[VX3][k][j][IEND]   - d->Vc[VX3][k][j][IEND-1]) / (r[IEND]   - r[IEND-1]);
        dvar2dr = (d->Vc[VX3][k][j][IEND-1] - d->Vc[VX3][k][j][IEND-2]) / (r[IEND-1] - r[IEND-2]);
        dvardr  = MINMOD_LIMITER(dvar1dr, dvar2dr);
        d->Vc[VX3][k][j][i] = d->Vc[VX3][k][j][i-1] + dvardr * (r[i] - r[i-1]);

        /* Re-inject initial Kluźniak & Kita disk analytical profile in equatorial belt */
        rcyl   = x1[i] * sin(x2[j]);
        eps2   = g_inputParam[EPS] * g_inputParam[EPS];
        coeff  = (g_gamma - 1.0) / g_gamma / eps2 * (1.0 / x1[i] - (1.0 - eps2 * g_gamma / (g_gamma - 1.0)) / rcyl);
        coeff  = MAX(coeff, 0.0);
        lambda = 11.0 / 5.0 / (1.0 + 64.0 / 25.0 * g_inputParam[ALPHAV] * g_inputParam[ALPHAV]);

        if (x2[j] >= 0.5 * CONST_PI - atan(1.25 * g_inputParam[EPS]) &&
            x2[j] <= 0.5 * CONST_PI + atan(1.25 * g_inputParam[EPS])) {

          d->Vc[RHO][k][j][i] = pow(coeff, 1.0 / (g_gamma - 1.0));
          
          if (d->Vc[RHO][k][j][i] == 0.0) { 
            d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][IEND];
            d->Vc[PRS][k][j][i] = eps2 * pow(coeff, g_gamma / (g_gamma - 1.0));
          }  

          d->Vc[VX1][k][j][i] = -g_inputParam[ALPHAV] / sin(x2[j]) * eps2
            * (10.0 - 32.0 / 3.0 * lambda * g_inputParam[ALPHAV] * g_inputParam[ALPHAV]
            - lambda * (5.0 - 1.0 / (eps2 * tan(x2[j]) * tan(x2[j])))) / sqrt(rcyl);

          d->Vc[VX3][k][j][i] = (sqrt(1.0 - 2.5 * eps2) + 2.0 / 3.0 * eps2
            * g_inputParam[ALPHAV] * g_inputParam[ALPHAV]
            * lambda * (1.0 - 6.0 / (5.0 * eps2 * tan(x2[j]) * tan(x2[j])))) / sqrt(rcyl);
        }

        /* Prevent inflow from the high-latitude coronal region */
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

  /* -----------------------------------------------------------------
     Polar boundaries (X2_BEG, X2_END) & Azimuthal boundaries (X3)
     (Handled automatically by PLUTO via pluto.ini settings)
     ----------------------------------------------------------------- */
}


#if BODY_FORCE != NO

/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 * \param [in] v  pointer to a cell-centered vector of primitive 
 *                variables
 * \param [out] g acceleration vector
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 *********************************************************************** */
{
    g[IDIR] = -1.0/x1/x1;
    g[JDIR] = 0.0;
    g[KDIR] = 0.0; 

}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the gravitational potential as function of the coordinates.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * 
 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
 *
 *********************************************************************** */
{
 //ccm--060824--uncomment the wanted pseudo-potential
 //set in definitions.h #define BODY_FORCE POTENTIAL 
// double q;
// q=1.25;
return -1./x1;//Newtonian
//return -1./(x1-2.);//0.0;//ccm--Paczynski-Wiita  
//return -(1.0/6.0)*(exp(6.0/x1)-1.0);//KluzniakLee
//return -1./x1+0.5*(q*q/x1/x1);//KluzniakNordstrom

//  return 0.0;
}
#endif

