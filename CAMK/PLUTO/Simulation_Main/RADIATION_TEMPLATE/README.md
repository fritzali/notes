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

`DiskFraction` replaces the tracer with a continuous `[0,1]` classifier built from two independent, multiplicatively combined logistic criteria:

The first is based on density. Rather than comparing against a rigid threshold or constructing histograms to separate distributions, a running
logarithmic radial reference profile is maintained for the disk and corona, respectively, at analysis cadence, by feeding the previous
classification back into the decision. This is achieved by `UpdateProfiles` using the current `diskfrac` to compute the weighted average
density for each component, with weights `f` and `1-f` for disk and corona, which for each radius are then temporally smoothed using an
exponential moving average. In the case of disk bins, their values are only accepted if the weighted disk volume makes up some predefined
fraction of the total bin volume. This way, bins that have never held significant amounts of disk material can copy from the nearest valid bin,
looking outward first and then inward, instead of constructing noisy profiles from empty cells. Once updated profiles have been built, the new
classification is calculated in logarithmic density space based on a sigmoid function.

The second criterion is rotational. It compares azimuthal velocity of each cell against that expected from Keplerian rotation in a Newtonian
potential at its cylindrical radius, without correcting for the Paczyński–Wiita modification. Analogous to the previous case, parses this
comparison through a sigmoid to arrive at an independent factor, which is useful especially for the exclusion of dense corona regions at
high latitude. Multiplying both factors gives a reasonable approximation for the passive tracer that is advected within the flux calculation.

At initialization, before any update has been computed, the classification artificially reproduces that of the native tracer exactly, with sharp
binary assignments between disk and corona, preventing artifacts and false labels inside the inner disk regions.



> For more details, check the comments in the respective files themselves.

### Message Passing Interface

When compiled for multiple cores, PLUTO uses the standard MPI library to split the simulation grid into pieces, a process called domain decomposition.
Each subset of the domain as well as a border of ghost cells shared with other pieces for boundary conditions gets assigned to a parallelized process or
rank, which has its own memory and no access to the whole grid. Since the radial profiles are constructed from averages across the entire domain, this
has to be treated specially by the `DiskFrac` implementation.

At each `DiskUpdate` call, all ranks compute from their local cells the local weighted radial profiles for density and volume. Next, these get reduced in
a collective operation across all ranks to construct the global profile, which is then passed back to each rank. With this identical information,
all individual ranks compute averages, smoothing, and classifications on their own, leading to a globally consistent update.
