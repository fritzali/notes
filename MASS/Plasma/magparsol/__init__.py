"""
magparsol: charged particle orbits in prescribed electromagnetic fields.

Example::

    from magparsol import EarthDipole, BorisC, dipole_initial_conditions

    field = EarthDipole()
    sim = BorisC(dipole_initial_conditions(), field, dt=1e-4, t_max=6.0,
                 store_dt=1e-3)
    history = sim.run()
    sim.plot_trajectory_3d(history, earth_sphere=True, ax_lim=4)
    sim.plot_energy(history)

Everything listed in ``__all__`` is importable from the package root:

* constants (``constants``)
* field models (``fields``)
* particle states and initial conditions (``particles``)
* integrators: Boris A, B, C and Runge–Kutta (``integrators``)
* trajectory history and gyration diagnostics (``diagnostics``)
* trajectory, energy and speed plots, live plotting (``plotting``)
* radiated power and spectra (``radiation``)
* field line plots (``fieldlines``)
* overview figure and GIF animations (``animation``)
* the shared plot style, ``style.use()`` (``style``)
"""

from magparsol.constants import (
    C, Q_E, M_E, M_P, EPS0, MU0, K_B,
    R_EARTH, DIPOLE_MOMENT, DIPOLE_TILT_DEG,
    B_FLOOR,
)

from magparsol.fields import (
    FieldModel,
    UniformB,
    UniformEB,
    CyclotronWaveField,
    EarthDipole,
    CustomField,
)

from magparsol.particles import (
    ParticleState,
    single_particle,
    random_ensemble,
    maxwellian_ensemble,
    relativistic_thermal_ensemble,
    dipole_initial_conditions,
)

from magparsol.integrators import (
    RKnonrel,
    RKrel,
    BorisA,
    BorisB,
    BorisC,
)

from magparsol.diagnostics import (
    TrajectoryHistory,
    relative_energy_error,
    gyrofrequency,
    gyroperiod,
    gyroradius,
    check_dt_resolution,
    suggest_dt,
)

from magparsol.plotting import (
    plot_trajectory_3d,
    plot_trajectory_2d,
    plot_energy,
    plot_speed,
    LivePlotter,
)

from magparsol.radiation import (
    radiated_power,
    total_radiated_energy,
    spectrum_fft,
    spectrum_retarded,
    ensemble_spectrum,
)

from magparsol.fieldlines import plot_field_lines

from magparsol.animation import (
    plot_overview,
    make_panel_gif,
    make_overview_gif,
    mean_gyrofrequency,
)

from magparsol import style

__version__ = "0.1.5"

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
