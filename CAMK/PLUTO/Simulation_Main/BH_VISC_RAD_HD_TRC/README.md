## PLUTO – Radiation Physics Template

In the PLUTO code, there is currently a problem where the tracer does not evolve properly if radiation physics are enabled. To fix this, I implemented
a drop in replacement called *DiskFraction* that identifies the disk with `1` and the corona with `0` by recovering the component from physical quantities.
This is then used to multiplicatively toggle the viscosity, resistivity, and opacity, the latter of which is of special interest for radiative transfer.

### Relevant Files

Changes were made to the following files:

1. `definitions.h` sets the standard simulation flags
2. `init.c` contains the main method implementation
3. `modifications.h` helps with additional definitions
4. `pluto.ini` initializes the grid and usual variables
5. `rad_step.c` now exposes the position for user defined opacities
6. `res_eta.c` replaces the tracer with *`diskfrac`*
7. `res_rhs.c` replaces the tracer with *`diskfrac`*
8. `userdef_output.c` writes additional outputs at each step
9. `visc_nu.c` replaces the tracer with *`diskfrac`*
10. `viscous_rhs.c` replaces the tracer with *`diskfrac`*

### Conceptual Implementation



> For more rigorous details, check the comments in the respective files themselves.
