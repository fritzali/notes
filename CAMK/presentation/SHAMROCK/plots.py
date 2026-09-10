from program import plot_density, plot_velocity, plot_field, plot_tracer

plot_density('../../PLUTO/Simulation_Main/BH_VISC_HD/', time=0.0)
plot_density('../../PLUTO/Simulation_Main/BH_VISC_HD/', time=250.0)
plot_density('../../PLUTO/Simulation_Main/BH_VISC_HD/', time=500.0)
plot_density('../../PLUTO/Simulation_Main/BH_VISC_HD/', time=750.0)
plot_density('../../PLUTO/Simulation_Main/BH_VISC_HD_TRC/', time=0.0)
plot_density('../../PLUTO/Simulation_Main/BH_VISC_HD_TRC/', time=250.0)
plot_density('../../PLUTO/Simulation_Main/BH_VISC_HD_TRC/', time=500.0)
plot_density('../../PLUTO/Simulation_Main/BH_VISC_HD_TRC/', time=750.0)
plot_tracer('../../PLUTO/Simulation_Main/BH_VISC_HD_TRC/', time=0.0, mode='truth')
