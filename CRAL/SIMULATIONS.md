## Astrophysical Simulations – Methods, Codes and Projects

### 1. Methods

#### 1.1. Hydrodynamics

The equations of gas dynamics can be discretized on a mesh that is fixed in space, on elements that move with the gas, or on something in between. The review of [Vogelsberger 2020] calls these the Eulerian, the Lagrangian and the arbitrary Lagrangian Eulerian formulation, and the list below is ordered from the first to the second.

* **Grid:** gas flows through fixed cells and fluxes from Riemann problems act on the cell faces, which gives high order accuracy, almost no oscillations behind shocks and low numerical diffusion according to [Vogelsberger 2020] while [Hopkins 2015] lists the weaknesses of meshes that do not move, namely advection errors and mixing that grow with the bulk velocity of the gas, effects of the grid alignment, weaker conservation of angular momentum and a less natural coupling to gravity
* **Refinement:** the same scheme on a mesh that is refined wherever a criterion asks for it, which is how grid codes cover the large dynamic range of collapse and of cosmological structure
* **Moving Mesh:** a finite volume scheme on a Voronoi mesh whose generating points move with the gas, introduced by [Springel 2010] to combine the shock capturing of grids with Galilean invariance, where the mesh deforms and changes its topology without tangling
* **Meshless:** kernels define volumes and effective faces between neighbouring particles and a Riemann solver acts across them, in a finite mass and a finite volume variant, for which [Hopkins 2015] reports proper convergence, sharp shocks, good mixing and no need for artificial diffusion, good conservation of angular momentum, and more noise than on a mesh
* **Riemann Particles:** smoothed particle hydrodynamics in which the force between two particles follows from a Riemann problem as formulated by [Inutsuka 2002] so that shocks need no artificial viscosity
* **Particles:** smoothed particle hydrodynamics follows particles of fixed mass, takes densities from kernel sums and handles shocks with artificial viscosity, so that the resolution follows the mass and energy, momentum, angular momentum and entropy are conserved together as long as no viscosity acts, while the comparison of [Hopkins 2015] finds slower convergence, more noise and numerical viscosity, weaker mixing and less accurate subsonic flow than in the newer methods
* **Finite Differences:** derivatives of high order on a fixed grid, used for weakly compressible turbulence in direct and in large eddy simulations

As a rule, grid methods are strongest for shocks and idealised turbulence, and methods that follow the mass are strongest where the range of densities is very large, as in the formation of galaxies.

#### 1.2. Gravity

Collisionless matter and the self gravity of the gas need a solver for the gravitational field, and [Vogelsberger 2020] groups the methods as follows.

* **Direct:** the sum over all pairs of particles, exact but only affordable for small systems
* **Tree:** distant groups of particles are replaced by their multipole moments
* **Multipole:** the fast multipole method evaluates forces between whole tree nodes and scales linearly with the number of particles
* **Mesh:** the density is deposited on a grid and the field follows from a Fourier transform or from a multigrid relaxation, which also works on refined meshes
* **Hybrid:** a mesh for large distances combined with a tree or with direct sums for short distances

### 2. Codes

#### 2.1. Particles and Moving Meshes

* **Lineage:** the Gadget code is the root of a large family, and [Springel 2021] notes that its second version was the last one with a full description until the public fourth version, that the third version was written for the Aquarius project and never fully documented, and that Arepo and Gizmo both descend from it
* **Arepo:** the moving mesh code with magnetohydrodynamics and a hybrid tree and mesh gravity solver, public since [Weinberger 2020] and used for the Illustris runs and their successors
* **Gizmo:** the meshless finite mass and finite volume methods of [Hopkins 2015] on top of the gravity solver and parallelization of the third Gadget version
* **Swift:** an open source code with task based parallelism described by [Schaller 2024] where gravity uses the fast multipole method with an optional mesh, several flavours of smoothed particle hydrodynamics are implemented with Sphenix as the default, and moving mesh and mesh free schemes are announced for later releases
* **Phantom:** a widely used smoothed particle hydrodynamics code that runs on processors without accelerators and whose scheme [David--Cléris 2025] reimplemented
* **Trees:** Gasoline and ChaNGa combine smoothed particle hydrodynamics with tree gravity, the latter with fast multipoles

#### 2.2. Grids

* **Cosmology:** Ramses and the earlier code of Kravtsov combine an adaptive mesh with a multilevel mesh solver for gravity, and Enzo combines it with a multigrid solver
* **Flash:** an adaptive mesh code with magnetohydrodynamics that is behind several of the interstellar medium and black hole projects below
* **Athena:** the framework of [Stone 2020] refines in blocks of an octree, supports curvilinear coordinates and solves magnetohydrodynamics with and without relativity
* **Pluto:** the modular code of [Mignone 2007] offers Godunov schemes for Newtonian and relativistic flows with and without magnetic fields in several coordinate systems, and was first applied to stellar and extragalactic jets
* **Accelerators:** Cholla by [Schneider 2015] runs hydrodynamics on a static mesh entirely on graphics processors, and [David--Cléris 2025] names Idefix as a Godunov code on a fixed grid for hybrid machines next to the frameworks Parthenon and Quokka
* **Framework:** the code described by [David--Cléris 2025] is written in C++ with the SYCL standard and aims at exascale machines, was tested on a large cluster of graphics processors, started from the particle scheme of Phantom, and its [repository] now also holds a Godunov particle solver and grid solvers named after Ramses and Zeus

#### 2.3. Gravity Only

* **Scale:** according to [Springel 2021] the codes that reached the largest particle numbers are specialised gravity solvers, among them PKDGRAV with tree and fast multipole forces and HACC with a hybrid of mesh and short range forces
* **Volumes:** the table of [Vogelsberger 2020] lists their largest dark matter runs, among them one for the Euclid survey, next to the earlier runs with tree and mesh codes

#### 2.4. Relativity

* **Comparison:** for accretion onto black holes with turbulence from the magnetorotational instability [Porth 2019] compared a set of independent codes for relativistic magnetohydrodynamics, among them Athena++ and members of the Harm family, and found that their agreement improves with resolution

#### 2.5. Turbulence

* **Pencil:** the code of the [Pencil Code Collaboration 2021] uses finite differences of high order for compressible flows with magnetic fields and particles, and serves direct simulations of turbulence and dynamos

### 3. Projects

#### 3.1. Dark Matter

* **Web:** the Millennium run of [Springel 2005] followed only dark matter in a large volume, the largest run of its time, and added the formation of galaxies and quasars in post processing
* **Successors:** [Vogelsberger 2020] lists the larger volumes that followed and the zooms on single haloes such as the Aquarius project

#### 3.2. Cosmological Volumes

* **Moving Mesh:** the Illustris run and its successor both use Arepo, and the model of [Pillepich 2018] for the successor adds magnetic fields, a kinetic mode of black hole feedback at low accretion rates and revised galactic winds
* **Calibration:** [Schaye 2015] calibrated the unresolved feedback from stars and black holes to the observed stellar mass function and sizes of nearby galaxies, with an updated formulation of smoothed particle hydrodynamics, and discussed how this limits the predictive power of such runs
* **Alternatives:** the table of [Vogelsberger 2020] lists further volumes, among them the Horizon run on an adaptive mesh, Magneticum and Romulus with particles and Simba with the meshless finite mass method
* **Clusters:** the suite of [Schaye 2023] calibrates stellar and black hole feedback by machine learning to the stellar mass function and to the gas fractions of clusters, varies feedback and cosmology including the neutrino mass, and contained the largest hydrodynamical run to the present day at its time
* **Lensing:** with this suite [McCarthy 2023] found that the effect of baryons is not large enough to remove the tension between the clustering amplitude from lensing and the one from the microwave background
* **Cold Gas:** the model of [Schaye 2026] runs on Swift with the Sphenix scheme, simulates the multiphase interstellar medium without a floor on temperature or pressure, follows dust grains coupled to chemistry and cooling, includes feedback from stars before they explode, and is calibrated to the stellar mass function, the sizes and the black hole masses of nearby galaxies
* **Reionization:** [Kannan 2022] coupled radiative transfer and magnetohydrodynamics on the moving mesh to an established galaxy model and obtained reionization histories that match the observed neutral fraction and optical depth

#### 3.3. Zooms

* **Galaxies:** single galaxies are resimulated at high resolution with every class of method, and [Vogelsberger 2020] lists Auriga on the moving mesh, Latte from the Fire project with the meshless finite mass method, and several suites with smoothed particle hydrodynamics
* **Clusters:** the same table lists zooms on galaxy clusters with particles and with the moving mesh

#### 3.4. Interstellar Medium

* **Stratified Boxes:** [Walch 2015] and [Kim 2017] simulate a vertical column through a galactic disc with self gravity, magnetic fields, cooling and supernovae, the first with chemistry of molecular gas in Flash and the second with resolved supernova remnants in Athena, to study how feedback sets the phases of the gas, the star formation rate and fountain flows
* **Sonic Scale:** [Federrath 2021] ran the largest simulation of supersonic turbulence on a uniform grid with Flash, resolved the supersonic and the subsonic cascade at once, and found that the transition between them is smooth instead of sharp
* **Forcing:** [Federrath 2010] compared solenoidal and compressive driving of supersonic turbulence and found a much wider density distribution for compressive driving at the same Mach number
* **Driving:** [Brucy 2020] simulated sections of galactic discs with stellar feedback and found that feedback regulates star formation in conditions like the Milky Way but not in gas rich galaxies, where additional driving on large scales is needed, and [Brucy 2023] added magnetic fields and the dependence on the compressibility of the driving
* **Comparison:** [Price 2010] ran the same driven supersonic turbulence with a grid code and with smoothed particle hydrodynamics and found agreement in the basic statistics once the number of particles roughly equals the number of cells

#### 3.5. Black Holes

* **Cold Accretion:** [Gaspari 2013] followed gas from the scale of a massive galaxy down to the black hole on an adaptive mesh in Flash and found that with cooling, heating and turbulence cold clouds and filaments condense, collide and rain onto the black hole, which raises the accretion rate far above the classical estimate for hot gas and allows feedback to regulate itself
* **Horizon Scale:** the comparison of [Porth 2019] established that independent codes agree on the turbulent accretion flows that are used to interpret images of black holes

[Vogelsberger 2020]: https://arxiv.org/abs/1909.07976
[Hopkins 2015]: https://arxiv.org/abs/1409.7395
[Springel 2010]: https://arxiv.org/abs/0901.4107
[Inutsuka 2002]: https://arxiv.org/abs/astro-ph/0206401
[Springel 2021]: https://arxiv.org/abs/2010.03567
[Weinberger 2020]: https://arxiv.org/abs/1909.04667
[Schaller 2024]: https://arxiv.org/abs/2305.13380
[David--Cléris 2025]: https://arxiv.org/abs/2503.09713
[repository]: https://github.com/Shamrock-code/Shamrock
[Stone 2020]: https://arxiv.org/abs/2005.06651
[Mignone 2007]: https://arxiv.org/abs/astro-ph/0701854
[Schneider 2015]: https://arxiv.org/abs/1410.4194
[Porth 2019]: https://arxiv.org/abs/1904.04923
[Pencil Code Collaboration 2021]: https://arxiv.org/abs/2009.08231
[Springel 2005]: https://arxiv.org/abs/astro-ph/0504097
[Pillepich 2018]: https://arxiv.org/abs/1703.02970
[Schaye 2015]: https://arxiv.org/abs/1407.7040
[Schaye 2023]: https://arxiv.org/abs/2306.04024
[McCarthy 2023]: https://arxiv.org/abs/2309.07959
[Schaye 2026]: https://doi.org/10.1093/mnras/stag375
[Kannan 2022]: https://arxiv.org/abs/2110.00584
[Walch 2015]: https://arxiv.org/abs/1412.2749
[Kim 2017]: https://arxiv.org/abs/1612.03918
[Federrath 2021]: https://arxiv.org/abs/2011.06238
[Federrath 2010]: https://doi.org/10.1051/0004-6361/200912437
[Brucy 2020]: https://doi.org/10.3847/2041-8213/ab9830
[Brucy 2023]: https://doi.org/10.1051/0004-6361/202244915
[Price 2010]: https://doi.org/10.1111/j.1365-2966.2010.16810.x
[Gaspari 2013]: https://arxiv.org/abs/1301.3130
