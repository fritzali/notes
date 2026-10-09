## PLUTO – Radiation with Compton Coupling, Classifier Fixes

Same setup as `BH_VISC_RAD_HD_TRC_COMP_LONG/`, the 120×90 run, with fixes to the DiskFraction classifier, see
[Classifier Fixes](#classifier-fixes) at the end. The sections before it describe the setup it builds on.

### Setup

Same setup as `BH_VISC_RAD_HD_TRC/`, with changes to the thermal coupling between gas and radiation. Without them, viscous heating drives the disk gas
far above the radiation temperature, from about `5e5 K` to `1e10–1e11 K` in the inner disk, because the only channel between gas and
radiation is Kramers absorption, which falls as `T^-3.5`. The overpressured gas then inflates bubbles that repeatedly break up the disk.

### Changes Relative to `BH_VISC_RAD_HD_TRC/`

1. `init.c` adds Compton energy exchange as an effective absorption opacity
2. `init.c` floors the gas temperature at `1e4 K`
3. `init.c` decouples matter and radiation where the gas moves faster than about `0.9 c`
4. `rad_step.c` keeps the radiation energy density of the primitive state current inside the implicit iteration
5. `definitions.h` lowers `UNIT_DENSITY` to `5.8e-10`, so that the disk can radiate its viscous heating, see below
6. `pluto.ini` lowers `BETAV` to `0.5` for the same reason, see below
7. `init.c` holds the Compton rate fixed above a gas temperature of `1e9 K`
8. `init.c` switches absorption and Compton off below `DiskFraction = 0.05`, so that the corona is thermally decoupled
9. `init.c` applies Thomson scattering to all gas instead of gating it by DiskFraction

### Compton Coupling

PLUTO treats scattering as coherent, so electron scattering never moves energy between gas and radiation. In an AGN disk, however,
Compton scattering is what keeps the gas at the radiation temperature. The nonrelativistic Compton exchange rate

    L_C = rho kappa_es c E_r 4 k (T_g - T_r) / (m_e c^2)

is reproduced exactly by an absorption term `rho kappa_c c (a_R T_g^4 - E_r)` with

    kappa_c = kappa_es 4 k T_r^4 / (m_e c^2 (T_g + T_r) (T_g^2 + T_r^2))

where the factor `T_g - T_r` has been cancelled analytically, so `kappa_c` is finite and positive everywhere. This part is moved from
scattering to absorption, so the total extinction stays exactly Thomson and only the energy exchange changes. Since `kappa_c` depends on
the radiation energy, `rad_step.c` now refreshes `v[ENR]` together with `v[PRS]` at every Newton iteration of the implicit step.

### Temperature Floor

The transparent corona lets radiation leave the disk surface freely, and with nothing heating them, surface layers cooled to about
`1e3 K` and collapsed into dense sheets. Below `1e4 K` the model has no valid physics anyway, as neither recombination nor molecular
opacities are included, so the gas temperature is floored there in the same pass that enforces the density floor.

### Decoupling Near the Horizon

With Newtonian dynamics in the Paczyński–Wiita potential, free fall gives `v^2 = 2 / (r - 2)` in units of `c`, so the plunging gas next
to the inner boundary at `r = 2.1` moves at several `c`. The nonrelativistic radiation module assumes `v << c`, and its `v^2/c^2` source
corrections turn into anti-damping once the gas approaches `c`. With Compton coupling the plunging gas carries enough radiation for this
to produce a NaN within about 12 time units. Opacities are therefore multiplied by a logistic factor that switches matter–radiation
coupling off above `v = 0.9 c`, with width `0.03 c`. The threshold must stay well above the orbital speed, which in this potential is
already `0.61 c` at `R = 6` and `0.47 c` at `R = 8`. A first try at `0.5 c` decoupled the whole inner disk, which then heated up just as
without Compton. At `0.9 c` only the plunging region inside `r ~ 4` is affected, and the gas dynamics there stay untouched. `BH_VISC_RAD_HD_TRC/` has the same superluminal plunge, and it survives there only because that region holds
much less radiation.

### Test Runs and Disk Parameters

Serial test runs to `t = 200` on the default grid, compared against `BH_VISC_RAD_HD_TRC/`, show that the Compton coupling works as
intended: the midplane gas stays at `1e6–1e7 K` and the disk remains radiation pressure dominated by `1e3–1e4`, where previously gas and
radiation pressure equalized within `t ~ 20`. The disk nevertheless inflates radiation-driven cavities from `t ~ 50` and is disrupted by
`t ~ 150`. This follows from the parameters rather than the numerics: with `BETAV = 1` and `UNIT_DENSITY = 1.16e-8`, viscous heating in
the initial disk exceeds what radiative diffusion can carry out by a factor `12–31` at `R = 8–20`, which corresponds to strongly
super-Eddington accretion. The radiation energy of the disk then doubles within about `40–60` time units at `R = 8–10`, so it cannot
stay thin.

Heating over cooling scales with `BETAV * UNIT_DENSITY`. This setup therefore uses `UNIT_DENSITY = 5.8e-10`, twenty times lower than
`BH_VISC_RAD_HD_TRC/`, which brings both into rough balance while keeping the disk optically thick, at `tau ~ 100` from midplane to
surface. The disk midplane then starts at `2–3e5 K`, and in the test run the disk stays thin and intact through `t = 200`.

### Full Run and Follow-Up Changes

The full run to `t = 1000` with the changes above kept the disk coherent throughout, with the midplane gas at `1e6–1e7 K`. Two issues
remained, which the following changes address. They have not been tested yet.

The inner disk still developed radiation-supported bubbles over a few hundred time units. At `UNIT_DENSITY = 5.8e-10` the heating
excess was still about `1.6` at `R = 8` and `1.2` at `R = 10`. Since this ratio scales with `BETAV * UNIT_DENSITY`, `BETAV` is now `0.5`,
which brings it to about `0.8` and `0.6` there. `BETAV` also sets the radial velocity of the analytic disk, so the initial and reinjected
disk change consistently with it. The outer disk is now somewhat overcooled, at a ratio of about `0.3` at `R = 20`.

The corona was Compton cooled and drained into the hole by `t ~ 300`. DiskFraction bottoms out near `1.5e-3` in the corona instead of
zero, because the rotation sigmoid never vanishes completely. At the corona temperatures near `1e11 K`, where `k T` is about `46 m_e c^2`
and the nonrelativistic rate is far outside its validity, even this leak drained the thermal energy of the corona by `100–1000` times
within `20–50` time units. Two changes address this:

- Opacities now vanish below `DiskFraction = 0.05`, and the range above is stretched back onto `[0,1]`, so disk material keeps its full
  opacity. Viscosity and resistivity are unaffected.
- Above `T_g = 1e9 K` the Compton rate is held at its value for `T_g = 1e9 K`. This needs
  `kappa_c = kappa_es 4 k T_r^4 (T_max - T_r) / (m_e c^2 (T_g^4 - T_r^4))`, which matches the previous expression at `T_g = T_max`.
  Simply capping `T_g` inside the previous expression would let the absorption form grow as `T_g^4` instead.

The corona drained in `BH_VISC_RAD_HD_TRC/` as well, though more slowly. Its mass fell from `2.8` to `0.66` code units by `t = 200`.
This happens because the corona is built in Newtonian hydrostatic equilibrium but evolves in the stronger Paczyński–Wiita potential, and
the outer boundary blocks coronal inflow, so it is never replenished. None of the changes here address that.

### Ungated Thomson Scattering

The run with the changes above showed a thin, dense layer along the disk surface from `t ~ 40`. It was up to `15` times denser than
the midplane, sat exactly at the DiskFraction boundary, and wrinkled and lifted off in a Rayleigh–Taylor-like way. The disk interior is
supported almost entirely by radiation, with radiation pressure exceeding gas pressure by `1e3–1e4`, but radiation could only push on gas
classified as disk, because all opacities were multiplied by DiskFraction. Gas lifted across the classification edge lost its opacity and
fell back, while gas just below was still pushed up, so material converged and piled up at the edge. Meanwhile the interior heated and
inflated, which inverted the vertical density profile. Once the corona had drained, after `t ~ 400`, the disk became smooth and
resembled `BH_VISC_HD/`.

Thomson scattering only transfers momentum, and its opacity per unit mass is the same for all ionized gas, so there is no physical reason
to restrict it to the disk. Only the energy exchange, Kramers absorption and Compton, needs gating to keep the hot corona from being
Compton cooled. Scattering is therefore no longer gated, while absorption and Compton keep the DiskFraction gating and the `0.05` cutoff.
The velocity decoupling above `0.9 c` still applies to all of them, because the unstable `v^2/c^2` terms scale with the total opacity.
These changes have not been tested yet. The corona, now with a Thomson optical depth of about `1–3`, may be pushed outward by radiation,
and the `kappa_scat` output is no longer zero there.

### Classifier Fixes

DiskFraction classifies each cell by two logistic criteria: its density relative to a corona and a disk reference density at
its radius, and its rotation relative to the Keplerian velocity. The reference densities are kept per logarithmic radial bin and
updated during the run.

#### What Went Wrong in the 120×90 Run

The stored classifier state in `BH_VISC_RAD_HD_TRC_COMP_LONG/diskfrac.out` shows that the density criterion inverted. When the
disk was disrupted around `t ≈ 4650`, the disk and corona references converged at all radii, as expected when the density
structure is mixed. Everywhere else they separated again as the disk was rebuilt, but in the bins at `r ≈ 2.4–3.4` the disk
reference stayed `0.9–1.15` dex below the corona reference until the end of the run at `t = 4e4`.

- With a negative span between the references, gas thinner than the corona reference counts as disk and denser gas as corona.
  The guard only caught a span of exactly zero.
- The references were built from the current classification, so dense gas classified as corona raised the corona reference and
  thin gas classified as disk lowered the disk reference. The inversion kept itself in place.
- The rotation criterion could not compensate. Gas plunging inside the innermost stable orbit keeps about its specific angular
  momentum there, `ℓ ≈ 3.7`, so its rotation relative to the Newtonian Keplerian velocity, `ℓ/sqrt(R) ≈ 2–2.5` at `r ≈ 3`, passes
  the threshold of `0.5` for disk and corona gas alike.
- The result was a box of cells at `r ≈ 2.3–3.6` up to about `45°` from the midplane classified as disk, while the dense
  plunging stream inside it had no viscosity, absorption or Compton coupling.

In addition, `pluto.ini` runs `Analysis()` every step and each update replaced `5 %` of the profiles, so their memory was only
about 20 steps, `0.12 t_g` at 120×90 and `0.17 t_g` at 90×45. `CORONA_AVG_FREQ` was defined but never used.

#### Changes

1. **Minimum span.** `DiskFraction()` uses at least `DF_SPAN_MIN = 0.5` dex for the logarithmic span between the references,
   so a degenerate or inverted span can no longer flip the density criterion.
2. **Repair of inverted bins.** After every update, `RepairInvertedBins()` replaces both reference densities of each bin whose
   span is below `DF_SPAN_MIN` with those of the nearest bin with a sufficient span, searching outward first. If no bin has one,
   e.g. while the whole disk is disrupted, the disk reference is raised to the minimum span above the corona reference. This
   removes the contaminated values, so an inversion cannot persist. `pluto.0.log` reports whenever the number of repaired bins
   changes.
3. **Smoothing by a time constant.** Each update weights the new estimate by `1 - exp(-Δt/DF_EMA_TAU)`, with `Δt` the time
   between two updates and `DF_EMA_TAU = 6 t_g`. The memory of the profiles is now the same at every resolution and update
   cadence. `CORONA_EMA_ALPHA` and the unused `CORONA_AVG_FREQ` are removed.
4. **Choice of the reference profiles**, set at compile time by `DISKFRAC_REFERENCE` in `modifications.h`:
   - `DF_REF_SELF`, the default: as before, the references are the disk and corona weighted mean densities of each bin. They
     follow the disk material wherever it is, also when it is lifted to higher altitudes.
   - `DF_REF_GEOM`: the disk reference is the maximum density of each bin and the corona reference the volume weighted geometric
     mean density of the cells above `z/R = DF_GEOM_ZCOR = 0.5`. These do not depend on the classification, so they cannot feed
     back on it.

   In both modes each cell is classified the same way, by its own density between the two references and its rotation. Lifted
   disk material with disk-like density and rotation is classified as disk either way; only the references differ.

The record format of `diskfrac.out` is unchanged, so restarts keep restoring the classifier state exactly, as long as the mode
and the tunables stay the same between segments. The weight of the repeated update after a restart uses the time step stored in
`restart.out`, which is the one the original run used.

These changes are compiled but not yet tested in a run. To check a run:

- `pluto.0.log` lists the repairs.
- The span `log10(ρ_d/ρ_c)` stays at or above `0.5` in every record of `diskfrac.out`.
- The snapshots no longer show the classified box at `r ≈ 2.3–3.6`.
