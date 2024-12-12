
import numpy as np
import matplotlib.pyplot as plt

base_dir = "/vast/home/yoonsoo/flastro/build_catalyst/"
# base_dir = "/vast/home/yoonsoo/cdss_data/HeatingCooling/"

t = []
e_fluid = []
e_rad = []

plt.figure(figsize=(8, 4), facecolor='w', dpi=100)

# D = 1.0
def analytic_solution(x, t):
    # return 2.0 + np.exp(-4 * np.pi**2 * t) * np.sin(2 * np.pi * x)  # sine wave
    return np.exp(- (x-10.0)**2 / (4 * (t + 1.0))) / np.sqrt(t + 1.0)

c = 0
x_dense = np.linspace(0, 20, 200)
# x_dense = np.linspace(0, 1, 200)

for i in 1 * np.arange(2):
    file_name = base_dir + "output-%05d-1D-0.raw" % i
    data = np.loadtxt(file_name).T
    time, cell_idx, x, density, pressure, velocity, totalE, radE = data

    plt.plot(x_dense, analytic_solution(x_dense, time[0]), '--', linewidth=0.6, color='C%d' % c)
    plt.plot(x, radE, 'o', label=r'$t = %.1g$' % time[0],
             markersize=3, color='C%d' % c)

    c += 1
    # plt.ylim(0, 1.05 * max(np.max(temperature_fluid), np.max(temperature_rad)))

plt.xlabel(r'$x$')
plt.ylabel('$E$')

plt.legend()

# plt.savefig("Diffusion_dt10_lev7_lowest_4.png", dpi=100, bbox_inches='tight')

plt.show()

