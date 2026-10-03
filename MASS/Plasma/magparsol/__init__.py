"""
magparsol
==========
A Python library for single-particle orbital simulations of astrophysical
and space plasmas.

Quick start
-----------
::

    from magparsol import EarthDipole, BorisC, dipole_initial_conditions, R_EARTH

    field    = EarthDipole()
    particle = dipole_initial_conditions()
    sim      = BorisC(particle, field, dt=1e-4, t_max=6.0, store_dt=1e-3)
    history  = sim.run()

    sim.plot_trajectory_3d(history, earth_sphere=True, ax_lim=4)
    sim.plot_energy(history)

Public API
----------
Constants   : C, Q_E, M_E, M_P, EPS0, MU0, K_B, R_EARTH, DIPOLE_MOMENT,
              DIPOLE_TILT_DEG, B_FLOOR
Fields      : FieldModel, UniformB, UniformEB, CyclotronWaveField,
              EarthDipole, CustomField
Particles   : ParticleState, single_particle, random_ensemble,
              maxwellian_ensemble, relativistic_thermal_ensemble,
              dipole_initial_conditions
Integrators : BorisA, BorisB, BorisC, RKnonrel, RKrel
Diagnostics : TrajectoryHistory, relative_energy_error, gyrofrequency,
              gyroperiod, gyroradius, check_dt_resolution, suggest_dt
Plotting    : plot_trajectory_3d, plot_trajectory_2d, plot_energy,
              plot_speed, LivePlotter
Radiation   : radiated_power, total_radiated_energy, spectrum_fft,
              spectrum_retarded, ensemble_spectrum
Field lines : plot_field_lines
Animation   : plot_overview, make_panel_gif, make_overview_gif,
              mean_gyrofrequency
Style       : style.use() applies the package look (palette, legends, dpi)
"""

# ── Constants ─────────────────────────────────────────────────────────────────
from magparsol.constants import (
    C, Q_E, M_E, M_P, EPS0, MU0, K_B,
    R_EARTH, DIPOLE_MOMENT, DIPOLE_TILT_DEG,
    B_FLOOR,
)

# ── Fields ────────────────────────────────────────────────────────────────────
from magparsol.fields import (
    FieldModel,
    UniformB,
    UniformEB,
    CyclotronWaveField,
    EarthDipole,
    CustomField,
)

# ── Particles ─────────────────────────────────────────────────────────────────
from magparsol.particles import (
    ParticleState,
    single_particle,
    random_ensemble,
    maxwellian_ensemble,
    relativistic_thermal_ensemble,
    dipole_initial_conditions,
)

# ── Integrators ───────────────────────────────────────────────────────────────
from magparsol.integrators import (
    RKnonrel,
    RKrel,
    BorisA,
    BorisB,
    BorisC,
)

# ── Diagnostics ───────────────────────────────────────────────────────────────
from magparsol.diagnostics import (
    TrajectoryHistory,
    relative_energy_error,
    gyrofrequency,
    gyroperiod,
    gyroradius,
    check_dt_resolution,
    suggest_dt,
)

# ── Plotting ──────────────────────────────────────────────────────────────────
from magparsol.plotting import (
    plot_trajectory_3d,
    plot_trajectory_2d,
    plot_energy,
    plot_speed,
    LivePlotter,
)

# ── Radiation ─────────────────────────────────────────────────────────────────
from magparsol.radiation import (
    radiated_power,
    total_radiated_energy,
    spectrum_fft,
    spectrum_retarded,
    ensemble_spectrum,
)

# ── Field lines ───────────────────────────────────────────────────────────────
from magparsol.fieldlines import plot_field_lines

# ── Animation / overview ──────────────────────────────────────────────────────
from magparsol.animation import (
    plot_overview,
    make_panel_gif,
    make_overview_gif,
    mean_gyrofrequency,
)

from magparsol import style

__version__ = "0.4.0"

__all__ = [
    # constants
    "C", "Q_E", "M_E", "M_P", "EPS0", "MU0", "K_B",
    "R_EARTH", "DIPOLE_MOMENT", "DIPOLE_TILT_DEG", "B_FLOOR",
    # fields
    "FieldModel", "UniformB", "UniformEB", "CyclotronWaveField",
    "EarthDipole", "CustomField",
    # particles
    "ParticleState", "single_particle", "random_ensemble",
    "maxwellian_ensemble", "relativistic_thermal_ensemble",
    "dipole_initial_conditions",
    # integrators
    "RKnonrel", "RKrel", "BorisA", "BorisB", "BorisC",
    # diagnostics
    "TrajectoryHistory", "relative_energy_error",
    "gyrofrequency", "gyroperiod", "gyroradius",
    "check_dt_resolution", "suggest_dt",
    # plotting
    "plot_trajectory_3d", "plot_trajectory_2d",
    "plot_energy", "plot_speed", "LivePlotter",
    # radiation
    "radiated_power", "total_radiated_energy",
    "spectrum_fft", "spectrum_retarded", "ensemble_spectrum",
    # field lines
    "plot_field_lines",
    # animation
    "plot_overview", "make_panel_gif", "make_overview_gif",
    "mean_gyrofrequency",
    # style
    "style",
]
