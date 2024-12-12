
import numpy as np
import matplotlib.pyplot as plt

base_dir = "/vast/home/yoonsoo/flastro/build_catalyst/"

t = []
e_fluid = []
e_rad = []

# for i in np.arange(0, 100 10):
for i in [100, 200]:
    file_name = base_dir + "output-%05d-1D-0.raw" % i
    data = np.loadtxt(file_name).T
    time, cell_idx, coord_x, density, pressure, velocity, totalE, radE = data

    # Make plot
    plt.figure(figsize=(6, 4), facecolor='w', dpi=100)

    plt.title("t = %.2f" % time[0])

    plt.plot(coord_x, density, label=r'$\rho$')
    # plt.plot(cell_idx, radE, label=r'$E_\mathrm{rad}$')

    plt.ylim(0, 1.1)

    plt.legend()

    plt.savefig("sod-HLL-%d.png" % i, dpi=100, bbox_inches='tight')
    plt.close()


# plt.figure(figsize=(7, 4), facecolor='w', dpi=100)

# for file_name in ["sod-hllc.raw", "sod-hll.raw"]:
#     data = np.loadtxt(file_name).T
#     time, cell_idx, coord_x, density, pressure, velocity, totalE, radE = data

#     plt.plot(coord_x, density)
#     # plt.plot(cell_idx, radE, label=r'$E_\mathrm{rad}$')

# plt.ylim(0, 0.5)

# plt.xlim(0.58, 0.7)
# # plt.ylim(0.2, 0.5)

#     plt.savefig("rk-fixed-%d.png" % i, dpi=100, bbox_inches='tight')
#     plt.close()
