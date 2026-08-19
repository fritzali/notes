#include "pluto.h"

/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 *
 *
 ***************************************************************** */
{
  int i, j, k;  
//ccm--210421
  int nv;
  double vi[NVAR];
  double nu1, nu2;
  double ***nu;
#if PHYSICS == MHD
  double J[3],eta[3];
  double ***num;
#endif
  double ***Te;  
  double *x1 = grid->x[IDIR];
  double *x2 = grid->x[JDIR];
  double *x3 = grid->x[KDIR];   
  nu=GetUserVar("nu");
#if PHYSICS == MHD  
  num=GetUserVar("num");
#endif
  Te=GetUserVar("Te");
//ccm

  DOM_LOOP(k,j,i){

//ccm--210421
for (nv = 0; nv < NVAR; nv++) vi[nv] = d->Vc[nv][k][j][i];
Visc_nu(vi, x1[i], x2[j], x3[k], &nu1, &nu2);
nu[k][j][i]=nu1;
#if PHYSICS == MHD
Resistive_eta (vi, x1[i], x2[j], x3[k], J, eta);
num[k][j][i]=eta[0];
#endif
Te[k][j][i]=vi[PRS]/vi[RHO];
//ccm
  }
}
/* ************************************************************* */
void ChangeOutputVar ()
/* 
 *
 * 
 *************************************************************** */
{ 
  Image *image;
//ccm--210421--additional variables output
  SetOutputVar("nu",DBL_OUTPUT,YES);
#if PHYSICS == MHD
  SetOutputVar("num",DBL_OUTPUT,YES);
#endif
  SetOutputVar("Te",DBL_OUTPUT,YES);

#if PARTICLES
  //SetOutputVar ("energy",PARTICLES_FLT_OUTPUT, NO);
//  SetOutputVar ("x1",    PARTICLES_FLT_OUTPUT, NO);
  //SetOutputVar ("vx1",   PARTICLES_FLT_OUTPUT, NO);
#endif
}
