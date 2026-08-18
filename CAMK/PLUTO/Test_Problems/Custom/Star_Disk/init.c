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
{
}

