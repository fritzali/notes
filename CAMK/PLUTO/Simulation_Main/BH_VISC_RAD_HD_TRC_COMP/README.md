## PLUTO – Radiation with Compton Coupling

Same setup as `BH_VISC_RAD_HD_TRC/`, with changes to the thermal coupling between gas and radiation. Without them, viscous heating drives the disk gas
far above the radiation temperature, from about `5e5 K` to `1e10–1e11 K` in the inner disk, because the only channel between gas and
radiation is Kramers absorption, which falls as `T^-3.5`. The overpressured gas then inflates bubbles that repeatedly break up the disk.

### Changes Relative to `BH_VISC_RAD_HD_TRC/`

1. `init.c` adds Compton energy exchange as an effective absorption opacity
2. `init.c` floors the gas temperature at `1e4 K`
3. `init.c` decouples matter and radiation where the gas moves faster than about `0.9 c`
4. `rad_step.c` keeps the radiation energy density of the primitive state current inside the implicit iteration

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

Heating over cooling scales with `BETAV * UNIT_DENSITY`. Lowering the density to `UNIT_DENSITY = 5.8e-10` brings both into rough
balance while keeping the disk optically thick, at `tau ~ 100` from midplane to surface, and in that run the disk stays thin and intact
through `t = 200`.
