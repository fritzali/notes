#define  PHYSICS                        MHD
#define  DIMENSIONS                     2
#define  GEOMETRY                       SPHERICAL
#define  BODY_FORCE                     VECTOR
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        1
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            9

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   CONSTRAINED_TRANSPORT
#define  BACKGROUND_FIELD               YES
#define  AMBIPOLAR_DIFFUSION            NO
#define  RESISTIVITY                    EXPLICIT
#define  HALL_MHD                       NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      EXPLICIT
#define  RADIATION                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  ALPHAM                         0
#define  MU                             1
#define  TEMPF                          2
#define  RHOC                           3
#define  RD                             4
#define  EPS                            5
#define  OMG                            6
#define  ALPHAV                         7
#define  DFLOOR                         8

/* [Beg] user-defined constants (do not change this line) */

#define  WARNING_MESSAGES               NO
#define  INTERNAL_BOUNDARY              YES
#define  SHOCK_FLATTENING               MULTID
#define  CT_EN_CORRECTION               YES
#define  UNIT_DENSITY                   8.5e-11
#define  UNIT_LENGTH                    1.392e11
#define  UNIT_VELOCITY                  2.1839e7
#define  VTK_VECTOR_DUMP                YES
#define  CT_EMF_AVERAGE                 ARITHMETIC
#define  LIMITER                        VANLEER_LIM

/* [End] user-defined constants (do not change this line) */
