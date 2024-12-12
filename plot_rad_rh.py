
import numpy as np
import matplotlib.pyplot as plt

base_dir = "/vast/home/yoonsoo/flastro/build/"

# Make plot
plt.figure(figsize=(6, 4), facecolor='w', dpi=100)

time_index = [0, 100000, 200000]

for i in time_index:
    file_name = base_dir + "output-%05d-1D-0.raw" % i
    data = np.loadtxt(file_name).T
    time, cell_idx, coord_x, density, pressure, velocity, totalE, radE = data

    # plt.plot(coord_x, density, label=r'$t = %.2g$' % time[0])
    # plt.plot(coord_x, radE, label=r'$t = %.2g$' % time[0])
    plt.plot(coord_x, pressure, label=r'$t = %.2g$' % time[0])

# plt.ylim(0, 1.1)

plt.legend()

# plt.savefig("sod-HLL-%d.png" % i, dpi=100, bbox_inches='tight')
# plt.close()

plt.show()
