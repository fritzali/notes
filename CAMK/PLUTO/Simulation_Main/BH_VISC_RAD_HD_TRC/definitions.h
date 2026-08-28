#define  PHYSICS                        HD
#define  DIMENSIONS                     2
#define  GEOMETRY                       SPHERICAL
#define  BODY_FORCE                     POTENTIAL
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        1
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            10

/* -- physics dependent declarations -- */

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  INCLUDE_LES                    NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      EXPLICIT
#define  RADIATION                      YES
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  BETAM                          0
#define  MU                             1
#define  TEMPF                          2
#define  RHOC                           3
#define  RD                             4
#define  EPS                            5
#define  OMG                            6
#define  BETAV                          7
#define  DFLOOR                         8
#define  REDUCED_C                      9

/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_DENSITY                   1.16e-8
#define  UNIT_VELOCITY                  2.998e10
#define  UNIT_LENGTH                    1.476923e13
#define  WARNING_MESSAGES               NO
#define  INTERNAL_BOUNDARY              YES
#define  SHOCK_FLATTENING               MULTID
#define  CT_EN_CORRECTION               YES
#define  VTK_VECTOR_DUMP                YES
#define  CT_EMF_AVERAGE                 UCT_HLL
#define  LIMITER                        VANLEER_LIM
#define  CHAR_LIMITING                  NO
#define  RADIATION_VAR_OPACITIES        YES
#define  RADIATION_IMPL                 RADIATION_NEWTON_NR_RAD
#define  USE_CMA                        YES
#define  RADIATION_DIFF_LIMITING        YES

/* [End] user-defined constants (do not change this line) */
