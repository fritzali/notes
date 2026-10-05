## PLUTO – Radiation with Compton Coupling

Same setup as `BH_VISC_RAD_HD_TRC/`, with two additions to the thermal physics. Without them, viscous heating drives the disk gas
far above the radiation temperature, from about `5e5 K` to `1e10–1e11 K` in the inner disk, because the only channel between gas and
radiation is Kramers absorption, which falls as `T^-3.5`. The overpressured gas then inflates bubbles that repeatedly break up the disk.

### Changes Relative to `BH_VISC_RAD_HD_TRC/`

1. `init.c` adds Compton energy exchange as an effective absorption opacity
2. `init.c` floors the gas temperature at `1e4 K`
3. `rad_step.c` keeps the radiation energy density of the primitive state current inside the implicit iteration

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
