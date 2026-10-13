## Numerical Simulations of AGN Accretion Disks in PLUTO

Code for the final presentation of my [CAMK](https://www.camk.edu.pl/en/) summer internship, summarizing the [report](../../report) on two dimensional
radiation hydrodynamic simulations of a thin viscous accretion disk around a supermassive black hole in [PLUTO](https://plutocode.ph.unito.it/),
the `DiskFraction` classifier replacing the passive tracer, and the steady states and density waves of the radiation runs.

The figures are made by `plot.py` from the simulation outputs in `../../PLUTO/Simulation_Main` with the routines of the report plotting script,
written to `content/` at slide size. Run `make` to plot and compile both the presentation and the flat handout, or `make slides` once the figures exist.
Single figure steps are selected with `SLIDE_ONLY`, for example `SLIDE_ONLY=domain,curves make plot`.
