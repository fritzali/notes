"""
examples.py
-----------
Usage examples for plasma_box covering all five original use cases.

Run individual examples or all of them:

    python examples.py              # runs all
    python examples.py cyclotron    # runs one example by keyword
"""

import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")   # non-interactive backend for headless/test runs
                        # replace with "TkAgg" or remove for interactive use
import matplotlib.pyplot as plt

# ── import everything from the library ────────────────────────────────────────
import sys, os
sys.path.insert(0, os.path.dirname(__file__))

from plasma_box import (
    # constants
    Q_E, M_P, M_E, C, R_EARTH,
    # fields
    UniformB, UniformEB, CyclotronWaveField, EarthDipole, CustomField,
    # particles
    single_particle, random_ensemble, dipole_initial_conditions,
    # integrators
    RKNonrel, RKRelativistic, BorisA, BorisB, BorisC,
    # diagnostics
    gyroperiod, suggest_dt, check_dt_resolution,
    # plotting
    plot_trajectory_3d, plot_trajectory_2d, plot_energy, plot_speed,
)


# ─────────────────────────────────────────────────────────────────────────────
# 1.  CYCLOTRON RESONANCE  (Boris-B, uniform B + oscillating E)
# ─────────────────────────────────────────────────────────────────────────────

def example_cyclotron_resonance():
    """Proton in cyclotron resonance: uniform B + sinusoidal E at omega_c.

    Matches CiklotronskaRezonanca/NEREL_PUTANJE_VERLE.py.
    Boris-B is used (= non-relativistic Boris-Verlet at v << c).
    """
    print("\n── Example 1: Cyclotron Resonance ──────────────────────────────")

    B_mag = 3e-10   # 3 µG [T]
    field = CyclotronWaveField(
        B=(0.0, 0.0, B_mag),
        q=Q_E, m=M_P,
        E_amp=5e-5,
        E_axis=1,    # oscillating Ey
    )

    # Print the resonant frequency for reference
    omega_c = field.omega_c
    T_c = 2 * np.pi / omega_c
    print(f"  Cyclotron frequency   ω_c = {omega_c:.4e} rad/s")
    print(f"  Cyclotron period      T_c = {T_c:.4e} s")

    particle = single_particle(
        q=Q_E, m=M_P,
        r0=[0.0, 0.0, 0.0],
        v0=[0.0, 1e4, 5e3],    # small initial velocity [m/s]
    )

    # dt << T_c for accuracy
    dt = T_c / 200.0
    t_max = 200.0 * T_c      # 200 gyro-periods

    sim = BorisB(particle.copy(), field, dt=dt, t_max=t_max,
                 store_dt=T_c, relativistic=False)
    print(f"  dt = {dt:.4e} s,  t_max = {t_max:.4e} s")
    history = sim.run(progress_every=0)
    print(f"  Stored {len(history)} snapshots.")

    # The particle should gain energy over time (resonant acceleration)
    fig, axes = plot_energy(history, sim.state.m, relativistic=False,
                            title="Cyclotron Resonance – Energy")
    fig.savefig("out_cyclotron_energy.png", dpi=120)
    print("  Saved: out_cyclotron_energy.png")

    fig, axes = plot_trajectory_2d(history, length_unit=1.0, unit_label="m",
                                   title="Cyclotron Resonance – Trajectory")
    fig.savefig("out_cyclotron_traj.png", dpi=120)
    print("  Saved: out_cyclotron_traj.png")
    plt.close("all")


# ─────────────────────────────────────────────────────────────────────────────
# 2.  NON-RELATIVISTIC RK4 – EARTH DIPOLE FIELD
# ─────────────────────────────────────────────────────────────────────────────

def example_rk_nonrel_dipole():
    """Non-relativistic proton in Earth's dipole field using RK4.

    Matches RungeKuttaNerelativisticki/DKP_mdipol.py.
    """
    print("\n── Example 2: Non-relativistic RK4, Earth Dipole ───────────────")

    field = EarthDipole()
    particle = single_particle(
        q=Q_E, m=M_P,
        r0=[0.0, -7.85 * R_EARTH, -1.53 * R_EARTH],
        v0=[0.0, 1.8e6, 1.8e6],
    )

    dt = 0.001
    t_max = 5000.0
    print(f"  dt = {dt} s,  t_max = {t_max} s")

    sim = RKNonrel(particle.copy(), field, dt=dt, t_max=t_max,
                   store_dt=1.0, relativistic=False)
    history = sim.run(progress_every=100_000)
    print(f"  Stored {len(history)} snapshots.")

    fig, ax = plot_trajectory_3d(
        history,
        length_unit=R_EARTH, unit_label=r"$R_E$",
        ax_lim=10.0, earth_sphere=True, color="green",
        title="Non-rel RK4 – Earth Dipole",
    )
    fig.savefig("out_rk_nonrel_dipole.png", dpi=120)
    print("  Saved: out_rk_nonrel_dipole.png")

    fig2, axes2 = plot_energy(history, sim.state.m, relativistic=False,
                              title="Non-rel RK4 – Energy conservation")
    fig2.savefig("out_rk_nonrel_energy.png", dpi=120)
    print("  Saved: out_rk_nonrel_energy.png")
    plt.close("all")


# ─────────────────────────────────────────────────────────────────────────────
# 3.  NON-RELATIVISTIC RK4 – ADAPTIVE STEP, UNIFORM B FIELD
# ─────────────────────────────────────────────────────────────────────────────

def example_rk_nonrel_adaptive():
    """Non-relativistic proton in uniform B – adaptive RK45 step size.

    Demonstrates the adaptive Dormand-Prince controller.
    """
    print("\n── Example 3: Non-relativistic Adaptive RK45, Uniform B ────────")

    B_mag = 3e-10
    field = UniformB(B=(0.0, 0.0, B_mag))
    particle = random_ensemble(N=1, q=Q_E, m=M_P, r0=[0.0, 0.0, 0.0],
                               v_max=0.01 * C, seed=42)

    Tc = gyroperiod(Q_E, M_P, B_mag)
    print(f"  Gyro-period T_c = {Tc:.4e} s")

    sim = RKNonrel(
        particle.copy(), field,
        dt=Tc / 50.0,           # initial step
        t_max=10.0 * Tc,
        adaptive=True,
        rtol=1e-8, atol=1e-11,
        store_dt=Tc / 20.0,
        relativistic=False,
    )
    history = sim.run()
    print(f"  Stored {len(history)} snapshots.")

    fig, axes = plot_trajectory_2d(history, length_unit=1.0, unit_label="m",
                                   title="Adaptive RK45 – Uniform B")
    fig.savefig("out_rk_adaptive_traj.png", dpi=120)
    print("  Saved: out_rk_adaptive_traj.png")

    fig2, _ = plot_energy(history, sim.state.m, relativistic=False,
                          title="Adaptive RK45 – Energy")
    fig2.savefig("out_rk_adaptive_energy.png", dpi=120)
    print("  Saved: out_rk_adaptive_energy.png")
    plt.close("all")


# ─────────────────────────────────────────────────────────────────────────────
# 4.  RELATIVISTIC RK4 – EARTH DIPOLE
# ─────────────────────────────────────────────────────────────────────────────

def example_rk_relativistic():
    """Relativistic proton in Earth's dipole field using RK4.

    Matches RungeKuttaRelativisticka/REL_PUTANJE_RK4.py.
    Initial velocity ~0.616c at 60° pitch angle.
    """
    print("\n── Example 4: Relativistic RK4, Earth Dipole ───────────────────")

    field = EarthDipole()
    particle = dipole_initial_conditions()   # 2.5 R_E, v~0.616c

    speed0 = np.linalg.norm(particle.v[0])
    gamma0 = 1.0 / np.sqrt(1.0 - (speed0/C)**2)
    print(f"  Initial speed  v/c = {speed0/C:.4f},  γ = {gamma0:.4f}")

    dt = 1e-4
    t_max = 6.0
    sim = RKRelativistic(particle.copy(), field, dt=dt, t_max=t_max,
                         store_dt=1e-3, relativistic=True)
    print(f"  dt = {dt} s,  t_max = {t_max} s")
    history = sim.run(progress_every=10_000)
    print(f"  Stored {len(history)} snapshots.")

    fig, ax = plot_trajectory_3d(
        history,
        length_unit=R_EARTH, unit_label=r"$R_E$",
        ax_lim=10.0, earth_sphere=True, color="darkorange",
        title="Relativistic RK4 – Earth Dipole",
    )
    fig.savefig("out_rk_rel_traj.png", dpi=120)
    print("  Saved: out_rk_rel_traj.png")

    fig2, _ = plot_energy(history, sim.state.m, relativistic=True,
                          title="Relativistic RK4 – Energy conservation")
    fig2.savefig("out_rk_rel_energy.png", dpi=120)
    print("  Saved: out_rk_rel_energy.png")
    plt.close("all")


# ─────────────────────────────────────────────────────────────────────────────
# 5.  RELATIVISTIC BORIS-C – EARTH DIPOLE (main target use case)
# ─────────────────────────────────────────────────────────────────────────────

def example_boris_c_earth():
    """Relativistic Boris-C proton in Earth's dipole field.

    This is the primary new contribution of the library.
    Uses the Zenitani & Umeda (2018) analytic rotation (Eqs. 11-12).
    """
    print("\n── Example 5: Relativistic Boris-C, Earth Dipole ───────────────")

    field = EarthDipole()
    particle = dipole_initial_conditions()

    speed0 = np.linalg.norm(particle.v[0])
    gamma0 = 1.0 / np.sqrt(1.0 - (speed0/C)**2)
    print(f"  Initial speed  v/c = {speed0/C:.4f},  γ = {gamma0:.4f}")

    dt = 1e-4
    t_max = 6.0
    sim = BorisC(particle.copy(), field, dt=dt, t_max=t_max,
                 store_dt=1e-3, relativistic=True)
    print(f"  dt = {dt} s,  t_max = {t_max} s")
    history = sim.run(progress_every=10_000)
    print(f"  Stored {len(history)} snapshots.")

    fig, ax = plot_trajectory_3d(
        history,
        length_unit=R_EARTH, unit_label=r"$R_E$",
        ax_lim=10.0, earth_sphere=True, color="royalblue",
        title="Relativistic Boris-C – Earth Dipole",
    )
    fig.savefig("out_boris_c_traj.png", dpi=120)
    print("  Saved: out_boris_c_traj.png")

    fig2, _ = plot_energy(history, sim.state.m, relativistic=True,
                          title="Boris-C – Energy conservation")
    fig2.savefig("out_boris_c_energy.png", dpi=120)
    print("  Saved: out_boris_c_energy.png")
    plt.close("all")


# ─────────────────────────────────────────────────────────────────────────────
# 6.  BORIS-A/B/C COMPARISON ON SAME INITIAL CONDITIONS
# ─────────────────────────────────────────────────────────────────────────────

def example_boris_comparison():
    """Compare Boris-A, Boris-B, Boris-C energy conservation.

    All three integrators run identical initial conditions in uniform B.
    Energy deviation highlights the phase error in Boris-B at large dt.
    """
    print("\n── Example 6: Boris-A/B/C Comparison, Uniform B ────────────────")

    B_mag = 3e-10
    field = UniformB(B=(0.0, 0.0, B_mag))
    v0 = 0.5 * C   # mildly relativistic

    base_particle = single_particle(
        q=Q_E, m=M_P,
        r0=[0.0, 0.0, 0.0],
        v0=[0.0, v0, 0.0],
    )

    Tc = gyroperiod(Q_E, M_P, B_mag)
    print(f"  Gyro-period T_c = {Tc:.4e} s")

    # Use a relatively large dt to show differences between A/B/C
    dt = Tc / 20.0
    t_max = 50.0 * Tc
    store_dt = Tc / 10.0

    results = {}
    for cls, label in [(BorisA, "Boris-A"), (BorisB, "Boris-B"), (BorisC, "Boris-C")]:
        sim = cls(base_particle.copy(), field, dt=dt, t_max=t_max,
                  store_dt=store_dt, relativistic=True)
        history = sim.run()
        results[label] = (history, sim.state.m)
        print(f"  {label}: {len(history)} snapshots")

    # Plot relative energy error for all three
    from plasma_box.diagnostics import relative_energy_error
    fig, ax = plt.subplots(figsize=(11, 5))
    colors = {"Boris-A": "steelblue", "Boris-B": "crimson", "Boris-C": "seagreen"}
    for label, (history, m) in results.items():
        err = relative_energy_error(history, m, relativistic=True)[:, 0] * 100.0
        ax.plot(history.t / Tc, err, label=label, color=colors[label], linewidth=1.2)
    ax.set_xlabel(r"Time $[T_c]$")
    ax.set_ylabel(r"$\Delta\varepsilon\,/\,\varepsilon_0$ [%]")
    ax.set_title(f"Boris A/B/C Energy Error Comparison  (dt = Tc/20, v=0.5c)")
    ax.axhline(0, color="k", linewidth=0.6, linestyle="--")
    ax.legend()
    ax.grid(True, alpha=0.35)
    plt.tight_layout()
    fig.savefig("out_boris_comparison.png", dpi=120)
    print("  Saved: out_boris_comparison.png")
    plt.close("all")


# ─────────────────────────────────────────────────────────────────────────────
# 7.  CUSTOM FIELD (legacy adapter demo)
# ─────────────────────────────────────────────────────────────────────────────

def example_custom_field():
    """Demonstrate wrapping a user-defined legacy scalar field function."""
    print("\n── Example 7: Custom Field Adapter ─────────────────────────────")

    # User's legacy function with the old (x, y, z, t) -> 6-tuple signature
    def my_field(x, y, z, t):
        Bz = 1e-9
        return 0.0, 0.0, Bz, 0.0, 0.0, 0.0

    field = CustomField(my_field, vector_api=False)
    B_mag = 1e-9
    particle = single_particle(q=Q_E, m=M_P, r0=[0,0,0], v0=[0, 1e5, 0])
    Tc = gyroperiod(Q_E, M_P, B_mag)
    dt = Tc / 100.0
    t_max = 5.0 * Tc

    sim = BorisC(particle.copy(), field, dt=dt, t_max=t_max,
                 store_dt=dt*5, relativistic=False)
    history = sim.run()
    print(f"  Custom field: {len(history)} snapshots stored.")

    fig, axes = plot_trajectory_2d(history, length_unit=1.0, unit_label="m",
                                   title="Custom Field (legacy adapter) – Boris-C")
    fig.savefig("out_custom_field.png", dpi=120)
    print("  Saved: out_custom_field.png")
    plt.close("all")


# ─────────────────────────────────────────────────────────────────────────────
# Dispatch
# ─────────────────────────────────────────────────────────────────────────────

EXAMPLES = {
    "cyclotron":  example_cyclotron_resonance,
    "rk_nonrel":  example_rk_nonrel_dipole,
    "adaptive":   example_rk_nonrel_adaptive,
    "rk_rel":     example_rk_relativistic,
    "boris_c":    example_boris_c_earth,
    "comparison": example_boris_comparison,
    "custom":     example_custom_field,
}

if __name__ == "__main__":
    if len(sys.argv) > 1:
        key = sys.argv[1].lower()
        if key not in EXAMPLES:
            print(f"Unknown example '{key}'. Choose from: {list(EXAMPLES)}")
            sys.exit(1)
        EXAMPLES[key]()
    else:
        for fn in EXAMPLES.values():
            fn()
    print("\nDone.")
