
import numpy as np
import matplotlib.pyplot as plt

base_dir = "/vast/home/yoonsoo/cdss_data/RadAdvection/"

t = []
e_fluid = []
e_rad = []

for i in 200 * np.arange(1, 7):
    file_name = base_dir + "output-%05d-1D-0.raw" % i
    data = np.loadtxt(file_name).T
    time, cell_idx, density, pressure, velocity, momentum, totalE, radE = data

    idx = len(time) // 2

    t += [time[idx]]
    e_fluid += [totalE[idx]]
    e_rad += [radE[idx]]

    # Make plot
    plt.figure(figsize=(6, 4), facecolor='w', dpi=100)

    plt.title("t = %.2f" % time[0])

    plt.plot(cell_idx, density, label=r'$\rho$')
    plt.plot(cell_idx, radE, label=r'$E_\mathrm{rad}$')

    plt.legend()

    plt.savefig("%d.png" % i, dpi=100, bbox_inches='tight')

