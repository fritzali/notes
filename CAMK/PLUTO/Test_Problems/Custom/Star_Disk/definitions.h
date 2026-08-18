#define  PHYSICS                             HD
#define  DIMENSIONS                          2
#define  COMPONENTS                          3
#define  GEOMETRY                            SPHERICAL
#define  BODY_FORCE                          VECTOR
#define  FORCED_TURB                         NO
#define  COOLING                             NO
#define  RECONSTRUCTION                      LINEAR
#define  TIME_STEPPING                       RK2
#define  DIMENSIONAL_SPLITTING               NO
#define  NTRACER                             1
#define  USER_DEF_PARAMETERS                 9

/* -- physics dependent declarations --- */

#define  DUST_FLUID                          NO
#define  EOS                                 IDEAL
#define  ENTROPY_SWITCH                      NO
#define  THERMAL_CONDUCTION                  NO
#define  VISCOSITY                           EXPLICIT
#define  ROTATING_FRAME                      NO

/* -- user-defined parameters (labels) -- */

#define  ALPHAM                              0
#define  MU                                  1
#define  TEMPF                               2
#define  RHOC                                3
#define  RD                                  4
#define  EPS                                 5
#define  OMG                                 6
#define  ALPHAV                              7
#define  DFLOOR                              8

/* [Beg] user-defined constants (do not change this line) */

#define  WARNING_MESSAGES                    NO
#define  INTERNAL_BOUNDARY                   YES
#define  SHOCK_FLATTENING                    MULTID
#define  CT_EN_CORRECTION                    YES
#define  UNIT_DENSITY                        8.5e-11
#define  UNIT_LENGTH                         1.392e11
#define  UNIT_VELOCITY                       2.1839e7
#define  VTK_VECTOR_DUMP                     YES
#define  CT_EMF_AVERAGE                      ARITHMETIC
#define  LIMITER                             VENLEER_LIM

/* [End] user-defined constants (do not change this line) */
