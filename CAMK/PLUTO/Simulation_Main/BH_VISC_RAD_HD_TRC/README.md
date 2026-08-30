## PLUTO – Radiation Physics Template

In the PLUTO code, there is currently a problem where the tracer does not evolve properly if radiation physics are enabled. To fix this, I implemented
a drop in replacement called *DiskFraction* that identifies the disk with `1` and the corona with `0` by recovering the component from physical quantities.
This is then used to multiplicatively toggle the viscosity, resistivity, and opacity, the latter of which is of special interest for radiative transfer.

### Relevant Files

Changes were made to the following files:

- `definitions.h` sets the standard simulation flags
- `init.c` contains the main method implementation
- `modifications.h` helps with additional definitions
- `pluto.ini` initializes the grid and usual variables
- `rad_step.c` now exposes the position for user defined opacities
- `res_eta.c` replaces the tracer with *`diskfrac`*
- `res_rhs.c` replaces the tracer with *`diskfrac`*
- `userdef_output.c` writes additional outputs at each step
- `visc_nu.c` replaces the tracer with *`diskfrac`*
- `viscous_rhs.c` replaces the tracer with *`diskfrac`*

### Conceptual Implementation



> For more rigorous details, check the comments in the respective files themselves.
