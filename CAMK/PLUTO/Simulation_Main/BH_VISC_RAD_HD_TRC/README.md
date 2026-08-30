## PLUTO – Radiation Physics Template

In the PLUTO code, there is currently a problem where the tracer does not evolve properly if radiation physics are enabled. To fix this, I implemented
a drop in replacement called *DiskFraction* that identifies the disk with `1` and the corona with `0` by recovering the component from physical quantities.
This is then used to multiplicatively toggle the viscosity, resistivity, and opacity, the latter of which is of special interest for radiative transfer.

### Relevant Files

Changes were made to the following files:

- `definitions.h`
- `init.c`
- `modifications.h`
- `pluto.ini`
- `rad_step.c`
- `res_eta.c`
- `res_rhs.c`
- `userdef_output.c`
- `visc_nu.c`
- `viscous_rhs.c`

### Conceptual Implementation



> For more rigorous details, check the comments in the respective files themselves.
