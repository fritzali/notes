# Shamrock install & test log — Fritz's laptop

Machine: Arch Linux (rolling, kernel 7.2.x), AMD Ryzen 8C/16T, 32 GB RAM (no swap), gcc 16.2.1, cmake 4.4.5, Python 3.14.
Raw command outputs: `shamrock_logs/NN_*.log` (this folder). Everything below "Laptop" was run by Fritz on the laptop;
anything run by Claude in its own sandbox is labelled as such and is NOT a laptop result.

## Laptop — install (Fri 6 Nov 2026)
| Step | Command (abridged) | Result | Log |
|---|---|---|---|
| 1 | `sudo pacman -Syu --needed ... llvm20 clang20 lld boost openmp openmpi ninja` | OK; installed llvm20/clang20 20.1.8, boost 1.92 + full system upgrade (107 pkgs) | 01_pacman.log |
| 2 | `git clone --depth 1 --recurse-submodules --shallow-submodules https://github.com/Shamrock-code/Shamrock.git`, renamed to `CRAL/SHAMROCK`, all `.git` removed | OK; Shamrock commit **f09b23795efc** (5 Nov 2026, main); source tree 62 MB | 02_clone.log |
| 3 | `SHAMROCK_BUILD_NPROC=10 ./env/new-env --machine archlinux.acpp --builddir ~/work/codes/shamrock-build -- --backend omp` | OK; ninja, 10 compile jobs, acpp target omp | 03_newenv.log |
| 4 | `source ./activate` (builds AdaptiveCpp) | OK; all 12 required pacman packages present; AdaptiveCpp develop **6c923f68** (5 Nov 2026) built in ~1–2 min; added `-DSHAMROCK_CHECK_SUBMODULES_COMMIT_HASH=Off` to CMAKE_OPT (needed because .git was removed inside a parent repo) | 04_activate.log |
| 5a | `shamconfigure` | OK; SYCL impl ACPPDirect, flags `--acpp-targets='omp' -O3 -march=native` | 05_configure.log |
| 5b | `shammake` | **OK, rc=0**, 539 targets in ~13 min; `shamrock` + `shamrock_test` + Python module `pyshamrock.cpython-314` built. 992 warnings, all `-Wpass-failed=transform-warning` ("loop not vectorized", optimiser hint from AdaptiveCpp's OpenMP launcher, harmless) | 06_build.log |

Notes:
- Shamrock prints version `2025.10.0+git.424c520.main`: the `424c520` hash is Fritz's *parent* repo commit (side-effect of removing Shamrock's .git); the real Shamrock commit is f09b23795efc.
- Pitfall met: `sudo pacman ... | tee` shows no download progress bar (output not a TTY) — looked hung, wasn't.

## Claude sandbox (not the laptop)
- Attempted a reference build in Claude's Ubuntu 24.04 sandbox (2 cores, debian-generic.acpp). It was stopped mid-AdaptiveCpp-build by the sandbox; **no sandbox result exists**. Source reading (task 4) was done on a checkout of the same commit f09b23795efc.

## Laptop — verification (Fri 6 Nov 2026, ~20:35)
| Step | Command | Result | Log |
|---|---|---|---|
| 6 | `./shamrock --smi` ; `./shamrock --smi --sycl-cfg 0:0` | OK. One device: id 0 = AdaptiveCpp OpenMP, CPU, AMD Ryzen AI 7 350 w/ Radeon 860M, 16 compute units, 31.3 GB. MPI world size 1. `- Code init: DONE.` One benign warning (mem_base_addr_align reported as 8 bits, Shamrock defaults to 64) | 07_smoke.log |
| 7a | `pull_reffiles` | cloned Shamrock-code/reference-files into build dir | 08_reffiles.log |
| 7b | `./shamrock_test --test-list` (without --sycl-cfg) | exited early: needs `--sycl-cfg`; not needed, skipped | 09_testlist.log |
| 7c | `./shamrock_test --sycl-cfg 0:0 --loglevel 1 --unittest` | **415/415 unit tests passed, 0 failed, rc=0, ~81 s**. Breakdown incl. 86 `shammodels` tests (17 ramses/Godunov, 9 sph, 55 gsph, 1 amr). The `NameError: name 'b' is not defined` warning is a deliberate negative test (src/tests/test_pythonscript_handle.cpp expects that script to fail). Only the *unit* test category was run (not `--validation`/benchmarks) | 10_unittest.log |

## Laptop — Python for example scripts (Fri 6 Nov, ~20:43)
- `shamrock` embeds the **system** Python it was configured with (/usr/bin/python3, 3.14.9). Example scripts need numpy/matplotlib/Pillow.
- Attempt 1: borrowed the mamba env `notes` (conda-forge Python 3.14.3, numpy 2.3.5, matplotlib 3.10.8, Pillow 12.1.1) via `--pypath-from-bin <notes python>` (documented flag). A trivial import test worked (embedded 3.14.9 imported numpy from the notes env).
- 9A with that setup **failed at `import shamrock`** (log 12_sod_godunov_amr.log): `shamrock.utils` optionally imports numba; numba from the notes env raised `RecursionError: Stack overflow (used 25 kB) while decoding a JSON object`, process aborted (signal 6). No simulation ran.
  Diagnosis: `--pypath-from-bin` replaces the *whole* sys.path, so the embedded 3.14.9 interpreter was using conda's 3.14.3 stdlib and compiled modules — mixing two builds of CPython. The exact mechanism of the spurious stack-overflow is my (Claude's) guess, not verified. Fix: a venv built from system Python (the project's documented approach), without numba (numba is optional in shamrock's pylib: `try: from numba import njit / except ImportError`).
- Also: `rc=` was empty because that terminal was zsh (`PIPESTATUS` is bash-only).
- Attempt 2: venv from system Python: `/usr/bin/python3 -m venv ~/work/codes/shamrock-build/.pyvenv; pip install numpy matplotlib pillow` -> Python 3.14.9, numpy 2.5.4, matplotlib 3.11.2, Pillow 13.0.0 (11_pydeps.log). Runs use `--pypath-from-bin <venv python>`.

## Laptop — examples
### 9A Sod tube, Godunov (RAMSES-like) solver + AMR — `examples/tests_ci/sod_tube_godunov_amr.py` (Fri 6 Nov ~20:47)
- Command: `shamrock --sycl-cfg 0:0 --pypath-from-bin <venv python> --rscript .../sod_tube_godunov_amr.py`
- Setup (from script): gamma 1.4, (rho,P) = (1,1)|(0.125,0.1), HLL Riemann solver, minmod limiter, face time interpolation, density/mass-based AMR (crit_mass 1.2e-6), evolved to t = 0.245.
- **Result: rc=0, 212 timesteps (dt ~1.1e-3), wall time 2 min 0 s** (user 12m51s, sys 13m08s on 16 threads).
- L2 errors vs analytic Sod: rho 0.07809, vx 0.17424, P 0.12382 — all below the script's CI pass thresholds (0.07913, 0.17763, 0.12516), so the script's built-in check **passed**.
- Plot `shamrock_results/A_sod_godunov_amr.png`: rarefaction (0.71–0.98), plateaus (rho 0.426 / 0.266, P 0.303, v 0.927), contact (~1.23) and shock (~1.43) all on the analytic curves; contact smeared over ~0.05–0.08, shock over ~2–3 cells. AMR level 2 for x < ~1.25 and level 1 to the right: the criterion is mass-based, so it refines the dense gas, not the discontinuities.
- Caveats: (1) ~6,600 AdaptiveCpp "Could not find embedded pointer in kernel blob" warnings because `ACPP_DEBUG_LEVEL=0` (normally set by `activate`) was not set in this shell — diagnostic only. (2) Shamrock printed "Current python is not the expected version … Expected prefix: …/miniforge3/envs/notes": the mamba env was still on PATH (mamba deactivate is a no-op in a bash without mamba's shell hook). Imports came from the venv and the solver/analysis are C++, so the result is not affected as far as Claude can tell, but this is a hygiene issue to fix before further runs.

### 9B Sod tube, SPH — `examples/sph/run_sod.py`, RESOL=128 (Fri 6 Nov ~20:51–20:57)
- Command: `RESOL=128 shamrock --sycl-cfg 0:0 --pypath-from-bin <venv python> --rscript .../sph/run_sod.py` (ACPP_DEBUG_LEVEL=0 set; PATH still contained the mamba env -> same "not the expected version" notice as 9A).
- Setup (from script): M6 kernel, adiabatic gamma 1.4, Cullen–Dehnen (CD10) artificial-viscosity switch (alpha 0–1), periodic box, HCP lattices with 8:1 density ratio (~83k particles, my estimate from resolution 128×24×24 + 1/8), CFL 0.1, to t = 0.245 with 10 analysis snapshots.
- **Result: rc=0, wall time 5 min 42 s** (user 46m39s, sys 40m02s).
- Final L2 vs analytic: rho 1.6155e-4, vx 1.1627e-3, vy 2.99e-5, vz 1.74e-7, P 1.2483e-4 (`shamrock_results/B_sod_sph_result.json`). For reference only: the separate CI script `tests_ci/sod_tube_sph.py` (same physical setup, not identical script) stores expected values 1.61549e-4 / 1.16270e-3 / 1.24836e-4 — agreement to ~1e-5 relative, not bitwise. I did not run that CI script.
- Plot `B_sod_sph_final.png`: rho/v/P on the analytic curves; usual SPH features: small pressure blip and smeared density at the contact (~0.72), post-shock velocity noise, shock smoothed over a few smoothing lengths (~0.93). The AV parameter alpha is ~0 in smooth flow and rises to 1 at the shock (the CD10 switch working as intended).
- `B_sod_sph_l2.png` as produced by the example is **garbled** (the example draws onto the leftover animation figure; with the Agg backend it is never cleared). Data are intact; re-plot with `shamrock_results/replot_sph_l2.py` (written by Claude).
- Repeated log warning "High interface/patch volume ratio … increase the patch split crit": Shamrock performance hint about its domain decomposition; not an error.

## Code reading — where external forces / source terms enter (done by Claude on a checkout of the same commit f09b23795efc; Fritz has not yet verified by reading)
Legend: **[sure]** = read in the code; **[inferred]** = my interpretation/design guess, to be confirmed with supervisors.

### Godunov / RAMSES-like solver — `src/shammodels/ramses/` (C++ namespace `shammodels::basegodunov`)
- [sure] An external-force framework already exists. Config: `include/shammodels/ramses/config/ExtForceConfig.hpp` — `std::variant` with currently one force type, `ConstantGravity`, plus a time-quadrature switch `coupling` = `Midpoint` (default) | `Trapezoidal` (`config/enum_ExtForceCoupling.hpp`). Python: `cfg.add_ext_force_constant_gravity(g=...)`, `cfg.set_ext_force_coupling(ExtForceCoupling.Trapezoidal)` (bindings in `src/pyRamsesModel.cpp`).
- [sure] Wiring: `Solver::init_solver_graph()` in `src/Solver.cpp` (block starting `if (solver_config.has_ext_forces() || solver_config.is_gravity_on())`, ~line 1574) builds an "External forces" OperationSequence appended to the solver sequence:
  1. `modules::ResetExtForceAcc` — zero the per-cell acceleration field `a_ext` (one per time level: t^{n+1/2} for midpoint; t^n and t^{n+1} for trapezoidal);
  2. one `AddExtForce*` node per configured force accumulates into `a_ext`: `AddExtForceConstantGravity` (`a_ext += g`), `AddExtForceField` (adds a per-cell acceleration field — used for self-gravity);
  3. `modules::ApplyExtForceGas` (`src/modules/ApplyExtForceGas.cpp`) kicks the conservative state at t^n with time-centred quantities: rho_half = rho + dt/2·∂t rho, (rho v)_half = rho v + dt/2·∂t(rho v); then rho v += dt·rho_half·a and rho E += dt·a·(rho v)_half. `ApplyExtForceGasTrapezoidal` / `ApplyExtForceDust*` are the variants.
- [sure] In `Solver::evolve_once()` the solver sequence (fluxes -> ∂tU, then the external-force kick of U^n in place) is evaluated first, then `modules::TimeIntegrator::forward_euler(dt)` adds dt·∂tU. So the source term is an additive, time-centred kick on the conservative variables, combined with the flux update in the same step.
- [sure] Known open issue in the code: TODO in `ApplyExtForceGas.cpp` — the midpoint energy update uses (rho v)_half without the kick, giving an internal-energy drift of ~rho a² dt²/2 per step; the CI test `examples/tests_ci/ramses_ext_force_coupling.py` checks exactly this drift (and that trapezoidal conserves e_int in free fall). Relevant to "biases introduced by the forcing".
- [sure] Design comment in Solver.cpp: midpoint needs forces at t^{n+1/2} (velocity-dependent forces would need half-step velocities); trapezoidal needs t^n and t^{n+1}, and a force at t^{n+1} may only depend on rho^{n+1}.
- [inferred] A turbulence forcing would most naturally be a new variant in `ExtForceConfig` plus a node that writes a(x,t) into `a_ext` (pattern: `AddExtForceConstantGravity` / `AddExtForceField`). How an OU (stochastic) acceleration maps onto midpoint vs trapezoidal levels is an open design question.

### SPH solver — `src/shammodels/sph/`
- [sure] Config shared with other models: `src/shammodels/common/include/shammodels/common/ExtForceConfig.hpp` — variants `PointMass`, `PN_PW` (Paczyński–Wiita), `LenseThirring`, `ShearingBoxForce`, `VerticalDiscPotential`, `VelocityDissipation` (f = −η v). Python: `cfg.add_ext_force_point_mass`, `_paczynski_wiita`, `_lense_thirring`, `_shearing_box`, `_velocity_dissipation`, `_vertical_disc_potential` (`src/pySPHModel.cpp`).
- [sure] Implementation: `modules::ExternalForces` (`include/.../modules/ExternalForces.hpp`, `src/modules/ExternalForces.cpp`), two entry points:
  - `compute_ext_forces_indep_v()` — called once per step in `Solver::evolve_once()` right after the leapfrog predictor (the solver sequence); zeroes and fills the per-particle field `axyz_ext` with velocity-independent forces at the predicted positions (point mass, PW, LT, vertical potential, inertial part of shearing box).
  - `add_ext_forces()` — called at the end of `Solver::update_derivs()`, after the hydro accelerations; does `axyz += axyz_ext` and adds velocity-dependent forces (shearing-box non-inertial part, `VelocityDissipation`). `update_derivs` is inside the corrector loop, so this can run more than once per step; then the leapfrog corrector kicks v with the new axyz.
  - Force kernels are solver-graph nodes in `src/shammodels/common/include/shammodels/common/modules/AddForce*.hpp` (e.g. `AddForceVelocityDissipation`: `axyz_ext -= eta * vxyz`).
- [inferred] An SPH turbulence forcing would most naturally be a velocity-independent force evaluated at particle positions in `compute_ext_forces_indep_v()` (new variant + a node like the `AddForce*` ones).

### Other facts
- [sure] No existing turbulence/OU forcing: no hit for ornstein/uhlenbeck/turbulen*/stirring in src/examples/doc (only unrelated "forcing" in GPU comm code). GitHub keyword search (turbulence, forcing) in Shamrock-code/Shamrock issues/PRs: no related item. Cannot rule out unpublished branches.
- [sure] `AGENTS.md`: bitwise run-to-run determinism is a "hard invariant". [inferred] A stochastic forcing must therefore use seeded, reproducible random numbers, identical on all MPI ranks.
- [not checked] Whether the timestep (CFL) accounts for external accelerations in either solver; whether Godunov ext forces are exercised with AMR refinement (the CI test uses a fixed grid).

### 9C Kelvin–Helmholtz, Godunov solver — `examples/ramses/run_kh.py` (Fri 6 Nov ~21:00–21:06)
- Command: `shamrock --sycl-cfg 0:0 --pypath-from-bin <venv python> --rscript .../ramses/run_kh.py` (unmodified gallery script).
- Setup (from script): 3D, 32³ blocks × 2³ cells = 64³ cells, no AMR, HLLC (PVRS) Riemann solver, minmod limiter, face time interpolation, gamma 5/3, periodic in x,y,z; rho 1 | 3.375, P = 3.5 both sides, slip velocity 1 (Mach 0.41 / 0.76 as printed), velocity perturbation sin(2π·2x) with Gaussian envelope around y = 0.5; evolved to t = 3.175 with 25 outputs.
- **Result: rc=0, 2032 timesteps (dt ≈ 1.5e-3), wall time 5 min 28 s** (user 42m46s, sys 41m18s). Output 979 MB (mostly 25 VTK dumps) in ~/work/codes/shamrock-runs/kh_ramses/_to_trash.
- Plots `shamrock_results/rho_periodic_0005/0012/0025.png` + `rho_periodic.gif`: t = 0.66 flat interface with slight modulation; t = 1.59 two bumps growing at x ≈ 0.25, 0.75; t = 3.175 two rolled-up KH billows, as expected for the seeded n = 2 mode. The second (periodic-wrap) shear layer at y = 0/1 is unperturbed in the setup and shows only a thin mixed strip.
- Caveats: this shell again lacked `ACPP_DEBUG_LEVEL=0` (~16,000 AdaptiveCpp accessor warnings in the log) and still had the mamba env on PATH ("not the expected version" notice). Neither changes the C++ solver; both are hygiene.
- Claude's runtime estimates for comparison: predicted ~2,000 steps (correct) and 10–25 min (actual 5.5 min — per-step cost overestimated).

## Optional: toy OU forcing check (standalone NumPy, NOT Shamrock) — `CRAL/ou_toy/ou_forcing_toy.py`
- Written by Claude. Fourier-space OU forcing (Eswaran & Pope / Federrath et al. 2010 style): 43 modes (half space, 1<|k|<3, parabolic amplitude peaked at |k|=2), Helmholtz projection with solenoidal weight zeta, seeded RNG. Two integrators: exact OU update (valid for any dt) and Euler–Maruyama.
- **Run by Claude in its sandbox** (Python 3.13, numpy 2.5.3; seed 42, T=1, 400 T sampled) — not yet on the laptop:
  - exact, dt/T=0.05: fitted T = 1.010 (coeffs) / 0.995 (field); variance ratio 1.006
  - exact, dt/T=0.5 : fitted T = 0.981 / 0.990; variance ratio 0.998
  - Euler, dt/T=0.05: fitted T = 0.975 / 0.962 (predicted 0.975); variance ratio 1.021 (pred. 1.026)
  - Euler, dt/T=0.5 : fitted T = 0.721 / 0.757 (predicted 0.721); variance ratio 1.345 (pred. 1.333)
  - zeta=1: max |k·f̂|/(|k||f̂|) = 7.7e-16 (divergence-free to round-off); zeta=0: curl = 0 exactly.
  - Scatter of a few % in fitted T is statistical (finite sample); the Euler–Maruyama bias (shorter T, larger variance when dt is not ≪ T) is systematic and matches the analytic formulas.
