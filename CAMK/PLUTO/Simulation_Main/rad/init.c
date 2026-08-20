/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Contains basic functions for problem initialization.

  The init.c file collects most of the user-supplied functions useful 
  for problem configuration.
  It is automatically searched for by the makefile.

  \author A. Mignone (mignone@ph.unito.it)
  \date   Sep 10, 2012

 \modified by M. Cemeljic (miki@camk.edu.pl)
 \date August 2025   
 Tilted field setup for pluto 4.4.3 by Miki, plus Sukalpa Kundu's
 cooling. Based on the setup described in Appendix of "Atlas" paper,
 Cemeljic, 2019, A&A, 624, A31. Can be used, with ksi=0 in axisymm.case
 dipole, in 2D&3D, HD and MHD, with minor changes/choices.
*/
/* ///////////////////////////////////////////////////////////////////// */
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
 double ksi;
 
  #if RADIATION
  #if RADIATION_NR
  g_reducedC = g_inputParam[REDUCED_C];
  #endif
  #endif

 ksi=0.*CONST_PI/10.;//mag. dipole tilt angle in rad

 //g_gamma=5./3.;//ccm261124--here we can change from PLUTO default g_gamma=5./3.

 rcyl=x1*sin(x2);
 eps2=g_inputParam[EPS]*g_inputParam[EPS];
// coeff=2./5./eps2*(1./x1-(1.-5./2.*eps2)/rcyl);
//ccm230824:
 coeff=(g_gamma-1.)/g_gamma/eps2*(1./x1-(1.-eps2*g_gamma/(g_gamma-1.))/rcyl); //coeff=rho_d/rho_d0

 lambda=11./5./(1.+64./25.*g_inputParam[BETAV]*g_inputParam[BETAV]); 

/* initial non-rotating adiabatic corona in hydrostatic equilibrium  */
// v[RHO] = g_inputParam[RHOC]*pow(x1,-3./2.);
// v[PRS] = 2./5.*g_inputParam[RHOC]*pow(x1,-5./2.);
 v[RHO] = g_inputParam[RHOC]*pow(x1,-1./(g_gamma-1.)); 
 v[PRS] = (g_gamma-1.)/g_gamma*g_inputParam[RHOC]*pow(x1,-g_gamma/(g_gamma-1.));

 pc=v[PRS];
 
  v[VX1] = 0.0;
  v[VX2] = 0.0;  
  v[VX3] = 0.0;

/* Keplerian adiabatic disk in vertical pressure equilibrium with the
   adiabatic corona, as given by Kluzniak & Kita (2000) */
   
  //v[PRS]=eps2*pow(coeff,5./2.);
  v[PRS]=eps2*pow(coeff,g_gamma/(g_gamma-1.));
    if (v[PRS] >= pc && rcyl > g_inputParam[RD])
      {//v[RHO] = pow(coeff,3./2.);
      v[RHO] = pow(coeff,1./(g_gamma-1.));
      v[PRS]=eps2*pow(coeff,g_gamma/(g_gamma-1.));
       v[VX1] = -g_inputParam[BETAV]/sin(x2)*eps2*(10.-32./3.
       *lambda*g_inputParam[BETAV]*g_inputParam[BETAV]
       -lambda*(5.-1./(eps2*tan(x2)*tan(x2))))/sqrt(rcyl);
       v[VX3] = (sqrt(1.-5./2.*eps2)+2./3.*eps2
       *g_inputParam[BETAV]*g_inputParam[BETAV]
       *lambda*(1.-6./(5.*eps2*tan(x2)*tan(x2))))/sqrt(rcyl); 
       v[TRC] = 1.0;     /* Track the disc material */
      }
    else
      {//v[PRS]=2./5.*g_inputParam[RHOC]*pow(x1,-5./2.);
      v[PRS]=(g_gamma-1.)/g_gamma*g_inputParam[RHOC]*pow(x1,-g_gamma/(g_gamma-1.));
       v[TRC] = 0.0;     /* Track the corona */
      }  

#if PHYSICS == MHD
#if BACKGROUND_FIELD == 1
   v[BX1] = 0.0;
   v[BX2] = 0.0;
   v[BX3] = 0.0;

   v[AX1] = 0.0;  
   v[AX2] = 0.0;
   v[AX3] = 0.0;
#else
//ccm--for tilted dipole stellar field without background field    
   v[BX1] = 2.*g_inputParam[MU]*(cos(x2)*cos(ksi)+sin(x2)*cos(x3)*sin(ksi))/(x1*x1*x1);
   v[BX2] = -g_inputParam[MU]*(sin(ksi)*cos(x2)*cos(x3)-cos(ksi)*sin(x2))/(x1*x1*x1);
   v[BX3] = g_inputParam[MU]*sin(ksi)*sin(x3)/(x1*x1*x1);
#endif

#endif//if MHD loop

  #if RADIATION
  v[ENR] = Blackbody(GetTemperature(v[RHO],v[PRS]))*v[TRC] ;
  v[FR1] = 0.;
  v[FR2] = 0.;
  v[FR3] = 0.;
  #endif
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
{/* dipole*/
   B0[0] = 2.*g_inputParam[MU]*cos(x2)/(x1*x1*x1);
   B0[1] = g_inputParam[MU]*sin(x2)/(x1*x1*x1);
   B0[2] = 0.0;                             
 /**/

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
 *
 * PURPOSE
 *
 *   Define user-defined boundary conditions.
 *
 **************************************************************** */
{
  int   i, j, k, nv, nth, thmin, thmax, jnew, nn, mm;
  double  *r,*rl,*rr,*thr,*thl,*Ar,*s,*dmu,*dx;;
  double a1, a2, a, smin, smax, Constreact;

  double dvar1dr, dvar2dr,dvardr;
  double dvar1dth, dvar2dth,dvardth;
  double db1dth, db2dth,dbdth;
  double vA, vB, Bpol, kappa, vphi_in, vpol2;
  double Bg[3], Bgi[3], Bgjp1[3], Bgjm1[3], Bgji[3];
  double vBi, Bpoli;
  double Bgo[3];
  double vBo, Bpolo, t;
  double *x1,*x2,*x3;
  double th1, th2, th1l, th2l, th1ln, th2ln, dthe, rhoth1, rhoth2;
  double rcyl, eps2, coeff, lambda, temp, Bpol2, vA_init, kk, vtot2, valfp2;
  double kapji, kapjip1, kapjp1, kapjm1, kapip2, kapjim2;
  double cs2, dden, dprs, dfact, eint, etot;
  double ksi;
  
//  EMF *emf = d->emf;
//  double ***Ex3e = emf->Ex3e;

// g_gamma=5./3.;//ccm190525--set 4./3. if we want non-adiabatic gamma in the disc
// gamc=5./3.;//ccm190525--if we want the adiabatic gamma in corona

 ksi=0.*CONST_PI/10.;
 
  if (side == 0) {    
/* -- Impose conditions for the solution inside domain near the star  -- */
/*MC130421: changing a primitive variable in an active zone requires
recomputing the conservative variables in that zone as well, by the call
to PrimToCons3D() --see bottom of page 51 Userguide */
    x1 = grid->xgc[IDIR];//gc is for conservative grid
    x2 = grid->xgc[JDIR];
    x3 = grid->xgc[KDIR];
/*### Avoid too small density near the star at the beginning and correct 
      the pressure to conserve the same sound speed for material in corona ###*/
//MC130421removed--   DOM_LOOP(k,j,i){ 
  RBox dom_box;
  TOT_LOOP(k,j,i){
    int convert_to_cons = 0;
//ccm--here set the rho(x1[i]) <= rho one cell above the star
  if (x1[i] <= x1[IBEG+1] && d->Vc[RHO][k][j][i]< g_inputParam[DFLOOR]){
//ccm--save values of rho and cs2 before reset to dfloor
     dden=d->Vc[RHO][k][j][i];
     cs2=g_gamma*d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i];
     d->Vc[RHO][k][j][i]=g_inputParam[DFLOOR];     
     dfact=dden/d->Vc[RHO][k][j][i];
//ccm--modify prs so that we stay with the same cs     
     d->Vc[PRS][k][j][i]=cs2*d->Vc[RHO][k][j][i]/g_gamma;
//ccm--to conserve momentum & energy, modify velocities     
     d->Vc[VX1][k][j][i]=dfact*d->Vc[VX1][k][j][i];
     d->Vc[VX2][k][j][i]=dfact*d->Vc[VX2][k][j][i];
     d->Vc[VX3][k][j][i]=dfact*d->Vc[VX3][k][j][i];
//ccm--ensure that TRC=0.0 in corona. Around reconnection region and outflows
//     sometimes TRC obtains spurious values, this is to prevent it. 
//  if(x2[j]<0.5*CONST_PI-atan(1.25*g_inputParam[EPS]) && d->Vc[RHO][k][j][i]<
  if(x2[j]<0.5*CONST_PI-atan(3.*g_inputParam[EPS]) && d->Vc[RHO][k][j][i]<
          1.e3*g_inputParam[DFLOOR])d->Vc[TRC][k][j][i]=0.0;
//ccm--for full [0,PI] half-plane: 
  if(x2[j]>0.5*CONST_PI+atan(3.*g_inputParam[EPS]) && d->Vc[RHO][k][j][i]<
          1.e3*g_inputParam[DFLOOR])d->Vc[TRC][k][j][i]=0.0;
 //Ex3e[0][j][IBEG-1] = 0.0;
// Ex3e[0][j][IBEG-1] = 0.0;

/*###### Avoid numerical heating in the corona ####### */
//ccm--smax is entropy at the stellar surface, smin is just a small number
 smax=8.61774;
 smin=0.01;
 d->Vc[PRS][k][j][i]=MAX(MIN(d->Vc[PRS][k][j][i],smax*pow(d->Vc[RHO][k][j][i],g_gamma))
 ,smin*pow(d->Vc[RHO][k][j][i],g_gamma));
     convert_to_cons = 1;
     }
     if (convert_to_cons) {
          RBoxDefine (i, i, j, j, k, k, CENTER, &dom_box);
          PrimToCons3D(d->Vc, d->Uc, &dom_box, grid);
       }
   } /* DOM_LOOP() */
} /* if (side == 0) */

  if (side == X1_BEG){  /* -- X1_BEG boundary -- */
#if RADIATION
      NRAD_LOOP(nv)  d->Vc[nv][k][j][i] = d->Vc[nv][k][j][IBEG];   ////// csk add radiation
#endif

r = grid->x[IDIR];
rr= grid->xr[IDIR];
rl= grid->xl[IDIR];
thr=grid->xr[JDIR];
thl=grid->xl[JDIR];
x1 = grid->x[IDIR];
x2 = grid->x[JDIR];
x3 = grid->x[KDIR]; 

/* ccm1214--in use of limiters using VanLeer
for rho and \vec{B}, and MinMod for pressure and \vec{v} */
    if (box->vpos == CENTER){
//ccmt--reverting to old loop structure, as new one misses a cell at jmin
//ccm   BOX_LOOP(box,k,j,i){
          X1_BEG_LOOP(k,j,i){

/*############# Linear extrapolation of rho #####################*/
 if (d->Vc[TRC][k][j][IBEG]> 0.01 && d->Vc[VX1][k][j][IBEG]<0.0){
 dvar1dr=(d->Vc[RHO][k][j][IBEG+2]-d->Vc[RHO][k][j][IBEG+1])/(r[IBEG+2]-r[IBEG+1]);
 dvar2dr=(d->Vc[RHO][k][j][IBEG+1]-d->Vc[RHO][k][j][IBEG])/(r[IBEG+1]-r[IBEG]);
 dvardr=VANLEER_LIMITER(dvar1dr,dvar2dr);
// dvardr=MINMOD(dvar1dr,dvar2dr);
 dvardr=MIN(dvardr,0.0);
 d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][i+1]-dvardr*(r[i+1]-r[i]); 
 d->Vc[PRS][k][j][i] = d->Vc[PRS][k][j][i+1]
      *pow(d->Vc[RHO][k][j][i]/d->Vc[RHO][k][j][i+1],g_gamma);
 d->Vc[TRC][k][j][i] =d->Vc[TRC][k][j][IBEG];
 }
 else
 {
  d->Vc[RHO][k][j][i] = g_inputParam[RHOC]*pow(r[i],-3./2.);
//ccm--put tempf in pluto.ini about few hundreds. It is so
// called effective temperature from polytropic law theory.  
 temp=2./5.;
 if(d->Vc[VX1][k][j][IBEG]>0.0) temp -= g_inputParam[TEMPF]
 *d->Vc[VX1][k][j][IBEG]*d->Vc[VX1][k][j][IBEG];
 d->Vc[PRS][k][j][i] = temp * g_inputParam[RHOC]*pow(r[i],-5./2.);
 d->Vc[TRC][k][j][i] = 0.0;
 }

//d->Vc[VX1][k][j][i]=d->Vc[VX1][k][j][IBEG];
d->Vc[VX1][k][j][i]=(d->Vc[VX1][k][j][IBEG] < 0) ? d->Vc[VX1][k][j][IBEG] : 0;  //csk;
d->Vc[VX2][k][j][i]=d->Vc[VX2][k][j][IBEG];
//ccm010720--in the MHD case, the below line is over-written at the
//end of the Special BC block.
d->Vc[VX3][k][j][i] = g_inputParam[OMG]*x1[i]*sin(x2[j]);

#if PHYSICS == MHD
//ccm240225--for NS, poloidal mag field is rotating with the stellar surface
//   d->Vc[BX1][k][j][i] = 2.*g_inputParam[MU]*(cos(x2[j])*cos(ksi)+sin(x2[j])*cos(x3[k]-g_inputParam[OMG]*g_time)*sin(ksi))/(x1[i]*x1[i]*x1[i]);
//   d->Vc[BX2][k][j][i] = -g_inputParam[MU]*(sin(ksi)*cos(x2[j])*cos(x3[k]-g_inputParam[OMG]*g_time)-cos(ksi)*sin(x2[j]))/(x1[i]*x1[i]*x1[i]);
//   d->Vc[BX2][k][j][i] = d->Vc[BX2][k][j][IBEG];//
//   d->Vc[BX3][k][j][i] = g_inputParam[MU]*sin(ksi)*sin(x3[k]-g_inputParam[OMG]*g_time)/(x1[i]*x1[i]*x1[i]);

/*############# Toroidal magnetic field BC     #################*/
// d->Vc[BX3][k][j][i] =d->Vc[BX3][k][j][i+1]*x1[i+1]/x1[i]; //ala Romanova //
//     d->Vc[BX3][k][j][i] = 0.;    // steady-state BC 

//########## Special BC on Bphi  ############
/**/
#if BACKGROUND_FIELD == 1
 if (i==IBEG-1 && j>=JBEG-1 && j<=JEND+1) {

 BackgroundField (x1[i+1],x2[j],x3[k],Bg);
 Bpol=sqrt(pow(Bg[0]+d->Vc[BX1][k][j][i+1],2)+pow(Bg[1]+d->Vc[BX2][k][j][i+1],2));
 vB= d->Vc[VX1][k][j][i+1]*(Bg[0]+d->Vc[BX1][k][j][i+1])+d->Vc[VX2][k][j][i+1]
 *(Bg[1]+d->Vc[BX2][k][j][i+1]);

 BackgroundField (x1[i],x2[j],x3[k],Bgi);
 Bpoli=sqrt(pow(Bgi[0]+d->Vc[BX1][k][j][i],2)+pow(Bgi[1]+d->Vc[BX2][k][j][i],2));
 vBi= d->Vc[VX1][k][j][i]*(Bgi[0]+d->Vc[BX1][k][j][i])+d->Vc[VX2][k][j][i]
 *(Bgi[1]+d->Vc[BX2][k][j][i]);

  db1dth=(d->Vc[BX3][k][j+1][i+1]-d->Vc[BX3][k][j][i+1])/(x2[j+1]-x2[j]);
  db2dth=(d->Vc[BX3][k][j][i+1]-d->Vc[BX3][k][j-1][i+1])/(x2[j]-x2[j-1]);
  if (x2[j]<0.) db2dth=2.*d->Vc[BX3][k][j][i+1]/(x2[j]-x2[j-1]);
  if (x2[j]>3.1415926) db2dth=-2.*d->Vc[BX3][k][j][i+1]/(x2[j+1]-x2[j]);
//ccm--in [0,pi] case, uncomment the line below--
  if (x2[j]>1.5708) db1dth=-2.*d->Vc[BX3][k][j][i+1]/(x2[j+1]-x2[j]);
//ccm--  
  dbdth=VANLEER_LIMITER(db1dth,db2dth);
//  dbdth=MINMOD_LIMITER(db1dth,db2dth);
  dvar1dth=(d->Vc[VX3][k][j+1][i+1]-d->Vc[VX3][k][j][i+1])/(x2[j+1]-x2[j]);
  dvar2dth=(d->Vc[VX3][k][j][i+1]-d->Vc[VX3][k][j-1][i+1])/(x2[j]-x2[j-1]);
  if (x2[j]<0.) dvar2dth=2.*d->Vc[VX3][k][j][i+1]/(x2[j]-x2[j-1]);
  if (x2[j]>3.1415926) dvar2dth=-2.*d->Vc[VX3][k][j][i+1]/(x2[j+1]-x2[j]);
//ccm--in [0,pi] case, uncomment the line below--
  if (x2[j]>1.5708) dvar1dth=-2.*d->Vc[VX3][k][j][i+1]/(x2[j+1]-x2[j]);
  dvardth=MINMOD_LIMITER(dvar1dth,dvar2dth);
//ccm--
 dvar1dr=(d->Vc[VX3][k][j][i+2]-d->Vc[VX3][k][j][i+1])/(r[i+2]-r[i+1]);
 dvar2dr=(d->Vc[VX3][k][j][i+1]-g_inputParam[OMG]*x1[i]*sin(x2[j])
 -vBi/pow(Bpoli,2)*d->Vc[BX3][k][j][i+1])/(r[i+1]-r[i]);
 dvardr=MINMOD_LIMITER(dvar1dr,dvar2dr);

 if(abs(d->Vc[VX3][k][j][IBEG]-g_inputParam[OMG]*x1[IBEG]*sin(x2[j])
 -d->Vc[BX3][k][j][IBEG]*vBi/pow(Bpoli,2))>0.01*g_inputParam[OMG]
 *x1[IBEG]*sin(x2[j])) {
 Constreact=0.1;
  } else {
  Constreact=1.;
   }
  d->Vc[BX3][k][j][i] = d->Vc[BX3][k][j][i+1]+ (r[i+1]-r[i])
  *(-sqrt(d->Vc[RHO][k][j][i+1])*Constreact*Bpol/(r[i+1]-r[i])
  *(g_inputParam[OMG]*x1[i+1]*sin(x2[j])+d->Vc[BX3][k][j][i+1]
  *vB/pow(Bpol,2)-d->Vc[VX3][k][j][i+1])-d->Vc[RHO][k][j][i+1]
  *d->Vc[VX1][k][j][i+1]*(dvar2dr+d->Vc[VX3][k][j][i+1]/r[i+1])
  -d->Vc[RHO][k][j][i+1]*d->Vc[VX2][k][j][i+1]/r[i+1]
  *(dvardth+d->Vc[VX3][k][j][i+1]*cos(x2[j])/sin(x2[j]))
  +(d->Vc[BX1][k][j][i+1]+Bg[0])/r[i+1]*d->Vc[BX3][k][j][i+1]
  +(d->Vc[BX2][k][j][i+1]+Bg[1])/r[i+1]*(dbdth+cos(x2[j])/sin(x2[j])
  *d->Vc[BX3][k][j][i+1]))/((d->Vc[BX1][k][j][i+1]+Bg[0])
  -d->Vc[RHO][k][j][i+1]*d->Vc[VX1][k][j][i+1]*vBi/pow(Bpoli,2));
}
#else
//ccm--if no Bg--
 if (i==IBEG-1 && j>=JBEG-1 && j<=JEND+1) {

 Bpol=sqrt(pow(d->Vc[BX1][k][j][i+1],2)+pow(d->Vc[BX2][k][j][i+1],2));
 vB= d->Vc[VX1][k][j][i+1]*(d->Vc[BX1][k][j][i+1])+d->Vc[VX2][k][j][i+1]
 *(d->Vc[BX2][k][j][i+1]);

 Bpoli=sqrt(pow(d->Vc[BX1][k][j][i],2)+pow(d->Vc[BX2][k][j][i],2));
 vBi= d->Vc[VX1][k][j][i]*(d->Vc[BX1][k][j][i])+d->Vc[VX2][k][j][i]
 *(d->Vc[BX2][k][j][i]);

  db1dth=(d->Vc[BX3][k][j+1][i+1]-d->Vc[BX3][k][j][i+1])/(x2[j+1]-x2[j]);
  db2dth=(d->Vc[BX3][k][j][i+1]-d->Vc[BX3][k][j-1][i+1])/(x2[j]-x2[j-1]);
  if (x2[j]<0.) db2dth=2.*d->Vc[BX3][k][j][i+1]/(x2[j]-x2[j-1]);
  if (x2[j]>3.1415926) db2dth=-2.*d->Vc[BX3][k][j][i+1]/(x2[j+1]-x2[j]);
//ccm--in [0,pi] case, uncomment the line below--
  if (x2[j]>1.5708) db1dth=-2.*d->Vc[BX3][k][j][i+1]/(x2[j+1]-x2[j]);
//ccm--  
  dbdth=VANLEER_LIMITER(db1dth,db2dth);
//  dbdth=MINMOD_LIMITER(db1dth,db2dth);
  dvar1dth=(d->Vc[VX3][k][j+1][i+1]-d->Vc[VX3][k][j][i+1])/(x2[j+1]-x2[j]);
  dvar2dth=(d->Vc[VX3][k][j][i+1]-d->Vc[VX3][k][j-1][i+1])/(x2[j]-x2[j-1]);
  if (x2[j]<0.) dvar2dth=2.*d->Vc[VX3][k][j][i+1]/(x2[j]-x2[j-1]);
  if (x2[j]>3.1415926) dvar2dth=-2.*d->Vc[VX3][k][j][i+1]/(x2[j+1]-x2[j]);
//ccm--in [0,pi] case, uncomment the line below--
  if (x2[j]>1.5708) dvar1dth=-2.*d->Vc[VX3][k][j][i+1]/(x2[j+1]-x2[j]);
//  dvardth=MINMOD_LIMITER(dvar1dth,dvar2dth);
  dvardth=VANLEER_LIMITER(dvar1dth,dvar2dth);
//ccm--
 dvar1dr=(d->Vc[VX3][k][j][i+2]-d->Vc[VX3][k][j][i+1])/(r[i+2]-r[i+1]);
 dvar2dr=(d->Vc[VX3][k][j][i+1]-g_inputParam[OMG]*x1[i]*sin(x2[j])
 -vBi/pow(Bpoli,2)*d->Vc[BX3][k][j][i+1])/(r[i+1]-r[i]);
// dvardr=MINMOD_LIMITER(dvar1dr,dvar2dr);
  dvardr=VANLEER_LIMITER(dvar1dr,dvar2dr);

 if(abs(d->Vc[VX3][k][j][IBEG]-g_inputParam[OMG]*x1[IBEG]*sin(x2[j])
 -d->Vc[BX3][k][j][IBEG]*vBi/pow(Bpoli,2))>0.01*g_inputParam[OMG]
 *x1[IBEG]*sin(x2[j])) {
 Constreact=0.1;
  } else {
  Constreact=1.;
   }
  d->Vc[BX3][k][j][i] = d->Vc[BX3][k][j][i+1]+ (r[i+1]-r[i])
  *(-sqrt(d->Vc[RHO][k][j][i+1])*Constreact*Bpol/(r[i+1]-r[i])
  *(g_inputParam[OMG]*x1[i+1]*sin(x2[j])+d->Vc[BX3][k][j][i+1]
  *vB/pow(Bpol,2)-d->Vc[VX3][k][j][i+1])-d->Vc[RHO][k][j][i+1]
  *d->Vc[VX1][k][j][i+1]*(dvar2dr+d->Vc[VX3][k][j][i+1]/r[i+1])
  -d->Vc[RHO][k][j][i+1]*d->Vc[VX2][k][j][i+1]/r[i+1]
  *(dvardth+d->Vc[VX3][k][j][i+1]*cos(x2[j])/sin(x2[j]))
  +(d->Vc[BX1][k][j][i+1])/r[i+1]*d->Vc[BX3][k][j][i+1]
  +(d->Vc[BX2][k][j][i+1])/r[i+1]*(dbdth+cos(x2[j])/sin(x2[j])
  *d->Vc[BX3][k][j][i+1]))/((d->Vc[BX1][k][j][i+1])
  -d->Vc[RHO][k][j][i+1]*d->Vc[VX1][k][j][i+1]*vBi/pow(Bpoli,2));
}
#endif
/**/
//ccm--after the changes above, recompute the BX3:
  dvardr=(d->Vc[BX3][k][j][IBEG]-d->Vc[BX3][k][j][IBEG-1])/(r[IBEG]-r[IBEG-1]);
  d->Vc[BX3][k][j][i] = d->Vc[BX3][k][j][i+1] -dvardr*(r[i+1]-r[i]);
  
//ccmt--at the end, set the stellar rotation, including E_phi=0 condition.
//Do not forget to set E_phi=0 in ct.c routine, too.
//ccm--180825--in pluto4.4.3 do not set it, as it breaks the MPI run.
//E=0 seems still to be well defined by the setting of dipole atop the
//star, not evolving it, using the correction for Bphi.
/*############# BC for Vphi #####################*/
//ccm281224--adding the magnetic torque at the stellar surface
//ccmt--at the end, set the stellar rotation, including E_phi=0 condition.
//If using CT, do not forget to set E_phi=0 in ct.c routine, too.
//ccm210525-->is E_phi=0 still a b.c. in the case with rotating tilted dipole?<--
//I did not do anything like this,

#if BACKGROUND_FIELD == 1
   BackgroundField (x1[i],x2[j],x3[k],Bgji);
   d->Vc[VX3][k][j][i] = g_inputParam[OMG]*x1[i]*sin(x2[j])
   +(d->Vc[VX1][k][j][i]*(Bgji[0]+d->Vc[BX1][k][j][i]) 
   +d->Vc[VX2][k][j][i]*(Bgji[1]+d->Vc[BX2][k][j][i]))   
   *(d->Vc[BX3][k][j][i])/((Bgji[0]+d->Vc[BX1][k][j][i])
   *(Bgji[0]+d->Vc[BX1][k][j][i])
   +(Bgji[1]+d->Vc[BX2][k][j][i])*(Bgji[1]+d->Vc[BX2][k][j][i]));
#else
   d->Vc[VX3][k][j][i] = g_inputParam[OMG]*x1[i]*sin(x2[j])
   +(d->Vc[VX1][k][j][i]*(d->Vc[BX1][k][j][i]) 
   +d->Vc[VX2][k][j][i]*(d->Vc[BX2][k][j][i]))   
   *(d->Vc[BX3][k][j][i])/((d->Vc[BX1][k][j][i])
   *(d->Vc[BX1][k][j][i])
   +(d->Vc[BX2][k][j][i])*(d->Vc[BX2][k][j][i]));
#endif       

#endif //end of if MHD loop 
}
  }else if (box->vpos == X2FACE){
#ifdef STAGGERED_MHD
//ccm190925--this is for the 2nd component (to BX3 from above) and for
//BX1 we do not define b.c., it will be defined by divB=0 condition.
 BOX_LOOP(box,k,j,i) d->Vs[BX2s][k][j][i] = d->Vs[BX2s][k][j][IBEG];
//ccm
//JTOT_LOOP(j) Ex3e[0][j][IBEG-1] = 0.0;

#endif
 }else if (box->vpos == X3FACE){
#ifdef STAGGERED_MHD
//ccm-not called in 2D, but yes in 3D, then def. smthing as e.g. Romanova
//ccm or the more complicated bphi condition. You have to include it here,
//not in cell centered part.
//  BOX_LOOP(box,k,j,i) d->Vs[BX3s][k][j][i] = d->Vs[BX3s][k][j][IBEG];
#endif
   } 
}

  if (side == X1_END){  /* -- X1_END boundary -- */

#if RADIATION
      NRAD_LOOP(nv)  d->Vc[nv][k][j][i] = d->Vc[nv][k][j][IEND];   ////// csk add radiation
#endif
r  = grid->x[IDIR];
x1 = grid->x[IDIR];
x2 = grid->x[JDIR];
x3 = grid->x[KDIR];

  if (box->vpos == CENTER){
        BOX_LOOP(box, k, j, i){
      d->Vc[TRC][k][j][i] =d->Vc[TRC][k][j][IEND];

/*############# Logarithmic extrapolation of rho, prs ###########*/
 a1=log10(d->Vc[RHO][k][j][IEND]/d->Vc[RHO][k][j][IEND-1])/log10(r[IEND]/r[IEND-1]);
 a2=log10(d->Vc[RHO][k][j][IEND-1]/d->Vc[RHO][k][j][IEND-2])/log10(r[IEND-1]/r[IEND-2]);
 a=VANLEER_LIMITER(a1,a2);
// a=MINMOD(a1,a2);
 a=MIN(a,0.0);
 d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][i-1]*pow(r[i]/r[i-1],a); 
//ccm--230824-- g_gamma
// d->Vc[PRS][k][j][i] = d->Vc[PRS][k][j][IEND]
//     *pow(d->Vc[RHO][k][j][i]/d->Vc[RHO][k][j][IEND],g_gamma);
 d->Vc[PRS][k][j][i] = d->Vc[PRS][k][j][IEND]
     *pow(d->Vc[RHO][k][j][i]/d->Vc[RHO][k][j][IEND],g_gamma); 
/*############# Ouflow BC for poloidal velocities #####################*/
 d->Vc[VX1][k][j][i] = d->Vc[VX1][k][j][IEND]; 
 d->Vc[VX2][k][j][i] = d->Vc[VX2][k][j][IEND]; 

#if PHYSICS == MHD
 dvar1dr=(d->Vc[BX3][k][j][IEND]-d->Vc[BX3][k][j][IEND-1])/(r[IEND]-r[IEND-1]);;
 dvar2dr=(d->Vc[BX3][k][j][IEND-1]-d->Vc[BX3][k][j][IEND-2])/(r[IEND-1]-r[IEND-2]);;
 dvardr=VANLEER_LIMITER(dvar1dr,dvar2dr);
// dvardr=MINMOD(dvar1dr,dvar2dr);
 d->Vc[BX3][k][j][i] = d->Vc[BX3][k][j][i-1]+dvardr*(r[i]-r[i-1]);

#endif
/*############# Linear extrapolation of Vphi #################*/    
 dvar1dr=(d->Vc[VX3][k][j][IEND]-d->Vc[VX3][k][j][IEND-1])/(r[IEND]-r[IEND-1]);;
 dvar2dr=(d->Vc[VX3][k][j][IEND-1]-d->Vc[VX3][k][j][IEND-2])/(r[IEND-1]-r[IEND-2]);;
 dvardr=MINMOD_LIMITER(dvar1dr,dvar2dr);
 d->Vc[VX3][k][j][i] = d->Vc[VX3][k][j][i-1]+dvardr*(r[i]-r[i-1]);

//ccm--at the disk outer boundary in R, set the initial HD values
 rcyl=x1[i]*sin(x2[j]);
 eps2=g_inputParam[EPS]*g_inputParam[EPS];
// coeff=2./5./eps2*(1./x1[i]-(1.-5./2.*eps2)/rcyl);
 coeff=(g_gamma-1.)/g_gamma/eps2*(1./x1[i]-(1.-eps2*g_gamma/(g_gamma-1.))/rcyl);
 coeff = MAX(coeff,0.0);
 lambda=11./5./(1.+64./25.*g_inputParam[BETAV]*g_inputParam[BETAV]);

//ccm--here we assume disk will puff-up for 25%, so we include this at b.c.
//ccm--for [0,pi/2] case use only the line below:
    if (x2[j] >= 0.5*CONST_PI-atan(1.25*g_inputParam[EPS])&&
//ccm--line below added for full [0,PI] half-plane
    x2[j] <= 0.5*CONST_PI+atan(1.25*g_inputParam[EPS])){
    d->Vc[RHO][k][j][i] = pow(coeff,1./(g_gamma-1.));//pow(coeff,3./2.);
    
    if(d->Vc[RHO][k][j][i] == 0.0){ d->Vc[RHO][k][j][i] = d->Vc[RHO][k][j][IEND];
    d->Vc[PRS][k][j][i] = eps2*pow(coeff,g_gamma/(g_gamma-1.));//eps2*pow(coeff,5./2.)
        }  

       d->Vc[VX1][k][j][i] = -g_inputParam[BETAV]/sin(x2[j])
       *eps2*(10.-32./3.*lambda*g_inputParam[BETAV]*g_inputParam[BETAV]
       -lambda*(5.-1./(eps2*tan(x2[j])*tan(x2[j]))))/sqrt(rcyl);
       d->Vc[VX3][k][j][i] = (sqrt(1.-5./2.*eps2)+2./3.*eps2
       *g_inputParam[BETAV]*g_inputParam[BETAV]
       *lambda*(1.-6./(5.*eps2*tan(x2[j])*tan(x2[j]))))/sqrt(rcyl); 
      }
//ccm--velocity could roll back even higher than the disk will pile-up,
//   so we prevent inflow in the corona from even higher region
   if (x2[j] <= 0.5*CONST_PI-atan(3.*g_inputParam[EPS]) ||
//ccm--line below added for full [0,PI] half-plane
   x2[j] >= 0.5*CONST_PI+atan(3.*g_inputParam[EPS]))
      if(d->Vc[VX1][k][j][i]<0.0){
      d->Vc[VX1][k][j][i] = 0.0;
      d->Vc[VX2][k][j][i] = 0.0;}

}
   }else if (box->vpos == X2FACE){
#ifdef STAGGERED_MHD
  BOX_LOOP(box,k,j,i)  d->Vs[BX2s][k][j][i] = d->Vs[BX2s][k][j][IEND];
#endif
  }else if (box->vpos == X3FACE){
#ifdef STAGGERED_MHD
//ccm not used in 2D, in 3D yes but then do better, kao MR ili kaj vec hoces
//  BOX_LOOP(box,k,j,i)  d->Vs[BX3s][k][j][i] = d->Vs[BX3s][k][j][IEND];
#endif
  }
  }  

  if (side == X2_BEG){  /* -- X2_BEG boundary -- */
    if (box->vpos == CENTER){
}else if (box->vpos == X1FACE){
      #ifdef STAGGERED_MHD
       BOX_LOOP(box, k, j, i) d->Vs[BX1s][k][j][i] = 0.0;
      #endif
   }else if (box->vpos == X3FACE){
      #ifdef STAGGERED_MHD
      #if INCLUDE_KDIR
       BOX_LOOP(box, k, j, i) d->Vs[BX3s][k][j][i] = 0.0;
//ccm
//JTOT_LOOP(j) Ex3e[0][j][IBEG-1] = 0.0;

      #endif 
      #endif
}

#if RADIATION
      NRAD_LOOP(nv)  d->Vc[nv][k][j][i] = d->Vc[nv][k][JBEG][i];  //csk add radiation
#endif
}

  if (side == X2_END){  /* -- X2_END boundary -- */
   if (box->vpos == CENTER){
     }else if (box->vpos == X1FACE){
      #ifdef STAGGERED_MHD
       BOX_LOOP(box, k, j, i) d->Vs[BX1s][k][j][i] = 0.0;
      #endif
    }else if (box->vpos == X3FACE){
      #ifdef STAGGERED_MHD
      #if INCLUDE_KDIR
       BOX_LOOP(box, k, j, i) d->Vs[BX3s][k][j][i] = 0.0;
       #endif
      #endif
 }
#if RADIATION
      NRAD_LOOP(nv)  d->Vc[nv][k][j][i] = d->Vc[nv][k][JEND][i];  //csk add radiation
#endif
  }

  if (side == X3_BEG){  /* -- X3_BEG boundary -- */
    X3_BEG_LOOP(k,j,i){
    #if RADIATION
      NRAD_LOOP(nv)  d->Vc[nv][k][j][i] = d->Vc[nv][KBEG][j][i];  //csk add radiation
    #endif
    }
  }
 
  if (side == X3_END) {  /* -- X3_END boundary -- */
    X3_END_LOOP(k,j,i){
    #if RADIATION
      NRAD_LOOP(nv)  d->Vc[nv][k][j][i] = d->Vc[nv][KEND][j][i];  //csk add radiation
    #endif

    }
  }
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



#if RADIATION_VAR_OPACITIES
        //source:  Carroll (1996), pages 274–276
        #include <math.h>

        #define G_BF 1.0        // Bound-Free Gaunt factor
        #define G_FF 1.0        // Free-Free Gaunt factor
        #define T_FACTOR 1e1   // Bound-Free correction factor, typically 1 < T < 100

        #define C_BF 4.34e25 // Kramer's law's bound-free constant, CGS
        #define C_FF 3.68e22 // Kramer's law's free-free constant,  CGS

        #define X     0.8    // Hydrogen mass fraction
        #define Z    0.00    // Metallicity

    //scattering constants (CGS)
    const double K_BF = 1;//C_BF * G_BF * Z * (1.0  + X) / T_FACTOR;
    const double K_FF = 1;//C_FF * (1.0  - Z) * (1.0  + X);
    const double K_ES = 1;//0.2  * (1.0 + X);


void UserDefOpacities(double *v, double *abs, double *scat){
    double rho = v[RHO];
    double T   = GetTemperature(v[RHO], v[PRS]);

    double kappa_es   = K_ES * v[TRC];                                        // cm^2/g
    double kappa_ffbf = (K_BF + K_FF) * rho * pow(T, -3.5) ;          // cm^2/g

    double fact=(v[VX3]>0.05 && v[RHO]>1e-4 && GetTemperature(v[RHO],v[PRS])<0.001) ? 1 : 0;

    *scat = fact * rho * kappa_es;        // Thomson scattering only
    *abs  = fact * rho * kappa_ffbf;      // free-free + bound-free true absorption
}
#endif

