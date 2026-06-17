"""
plasma_box
==========
A Python library for single-particle orbital simulations of astrophysical
and space plasmas.

Quick start
-----------
::

    from plasma_box import (
        EarthDipole, BorisC,
        dipole_initial_conditions,
    )
    from plasma_box.constants import R_EARTH

    field    = EarthDipole()
    particle = dipole_initial_conditions()
    sim      = BorisC(particle, field, dt=1e-4, t_max=6.0, store_dt=1e-3)
    history  = sim.run()

    sim.plot_trajectory_3d(history, earth_sphere=True)
    sim.plot_energy(history)

Public API
----------
Constants : see ``plasma_box.constants``
Fields    : UniformB, UniformEB, CyclotronWaveField, EarthDipole, CustomField
Particles : ParticleState, single_particle, random_ensemble, dipole_initial_conditions
Integrators: RKNonrel, RKRelativistic, BorisA, BorisB, BorisC
Diagnostics: TrajectoryHistory, relative_energy_error, gyroperiod, gyrofrequency,
             gyroradius, check_dt_resolution, suggest_dt
Plotting  : plot_trajectory_3d, plot_trajectory_2d, plot_energy, plot_speed, LivePlotter
"""

# ── Constants ─────────────────────────────────────────────────────────────────
from plasma_box.constants import (
    C, Q_E, M_E, M_P,
    R_EARTH, DIPOLE_MOMENT, DIPOLE_TILT_DEG,
    B_FLOOR,
)

# ── Fields ────────────────────────────────────────────────────────────────────
from plasma_box.fields import (
    FieldModel,
    UniformB,
    UniformEB,
    CyclotronWaveField,
    EarthDipole,
    CustomField,
)

# ── Particles ─────────────────────────────────────────────────────────────────
from plasma_box.particles import (
    ParticleState,
    single_particle,
    random_ensemble,
    dipole_initial_conditions,
)

# ── Integrators ───────────────────────────────────────────────────────────────
from plasma_box.integrators.runge_kutta import (
    RKNonrel,
    RKRelativistic,
)
from plasma_box.integrators.boris import (
    BorisA,
    BorisB,
    BorisC,
)

# ── Diagnostics ───────────────────────────────────────────────────────────────
from plasma_box.diagnostics import (
    TrajectoryHistory,
    relative_energy_error,
    gyrofrequency,
    gyroperiod,
    gyroradius,
    check_dt_resolution,
    suggest_dt,
)

# ── Plotting ──────────────────────────────────────────────────────────────────
from plasma_box.plotting import (
    plot_trajectory_3d,
    plot_trajectory_2d,
    plot_energy,
    plot_speed,
    LivePlotter,
)

__version__ = "0.1.0"

__all__ = [
    # constants
    "C", "Q_E", "M_E", "M_P",
    "R_EARTH", "DIPOLE_MOMENT", "DIPOLE_TILT_DEG", "B_FLOOR",
    # fields
    "FieldModel", "UniformB", "UniformEB", "CyclotronWaveField",
    "EarthDipole", "CustomField",
    # particles
    "ParticleState", "single_particle", "random_ensemble",
    "dipole_initial_conditions",
    # integrators
    "RKNonrel", "RKRelativistic", "BorisA", "BorisB", "BorisC",
    # diagnostics
    "TrajectoryHistory", "relative_energy_error",
    "gyrofrequency", "gyroperiod", "gyroradius",
    "check_dt_resolution", "suggest_dt",
    # plotting
    "plot_trajectory_3d", "plot_trajectory_2d",
    "plot_energy", "plot_speed", "LivePlotter",
]
