
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate

boltzmann_constant = 1.380649e-16
radiation_constant = 7.56573325e-15
proton_mass = 1.67262192369e-24
speed_of_light = 2.99792458e10

# --------------------------------------------------------------------
#   Problem parameters : use same values from the input file
# --------------------------------------------------------------------

T_fluid = 1000
T_rad   = 3000
rho_fluid = 1.0e-11

gamma = 1.4
kappa = 1.0
mean_molecular_weight = 1.0

base_dir = "/vast/home/yoonsoo/flastro/build/"
# base_dir = "/vast/home/yoonsoo/cdss_data/HeatingCooling/"


# --------------------------------------------------------------------
#   Load simulation data
# --------------------------------------------------------------------
t = []
e_fluid = []
e_rad = []

for i in 1 * np.arange(41):
    file_name = base_dir + "output-%05d-1D-0.raw" % i
    data = np.loadtxt(file_name).T
    time, cell_idx, x, density, pressure, velocity, totalE, radE = data

    idx = len(time) // 2

    t += [time[idx]]
    e_fluid += [totalE[idx]]
    e_rad += [radE[idx]]

e_fluid = np.array(e_fluid)
e_rad   = np.array(e_rad)

temperature_rad = (e_rad / radiation_constant) ** (0.25)
temperature_fluid = (e_fluid * proton_mass * (gamma - 1.0)) / (boltzmann_constant * density[idx])


# --------------------------------------------------------------------
#   Solve ODEs to get a reference solution
# --------------------------------------------------------------------
particle_mass = mean_molecular_weight * proton_mass

def ode_rhs(t, y):
    e = y[0]
    E = y[1]
    qdot = speed_of_light * kappa * rho_fluid * (
        E - radiation_constant * (
            (e * (gamma - 1) * particle_mass) / (rho_fluid * boltzmann_constant)
        )**4)
    return [qdot, -qdot]

e0 = rho_fluid * boltzmann_constant * T_fluid / ((gamma-1) * particle_mass)
E0 = radiation_constant * T_rad**4

result = scipy.integrate.solve_ivp(ode_rhs, t_span=[t[0], t[-1]], y0=[e0, E0],
                                   method='Radau', rtol=1e-8, atol=1e-8)
Tf_ref = result.y[0] * (gamma - 1) * particle_mass / (rho_fluid * boltzmann_constant)
Tr_ref = (result.y[1] / radiation_constant)**0.25

# --------------------------------------------------------------------
#   Plot results
# --------------------------------------------------------------------

plt.figure(figsize=(5, 5), facecolor='w', dpi=100)

# plt.plot(t, e_fluid, 'o', label=r'$E_\mathrm{fluid}$')
# plt.plot(t, e_rad, 'o', label=r'$E_\mathrm{rad}$')

plt.plot(result.t, Tf_ref, color='k', linewidth=0.4, alpha=0.5)
plt.plot(result.t, Tr_ref, color='k', linewidth=0.4, alpha=0.5)

plt.plot(t, temperature_fluid, 'o', label=r'$T_\mathrm{fluid}$',
         markersize=2.5)
plt.plot(t, temperature_rad, 'o', label=r'$T_\mathrm{rad}$',
         markersize=2.5)

plt.xscale('log')
plt.xlim(0.2, 1.1 * t[-1])
# plt.yscale('log')
# plt.ylim(0, 1.05 * max(np.max(temperature_fluid), np.max(temperature_rad)))

plt.ylabel('Temperature [K]', labelpad=10)

plt.xlabel(r'$t$', fontsize=11)
plt.legend()

plt.savefig("heating_and_cooling.png", dpi=120, bbox_inches='tight')
# plt.show()
