
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import cm, ticker

plt.rcParams.update({"text.usetex": True})

# base_dir = "/vast/home/yoonsoo/flastro/build/"
base_dir = "/vast/home/yoonsoo/flastro/build_gpu/"
# base_dir = "/vast/home/yoonsoo/cdss_data/RadiativeKHI_GPU/"

######################################################

# time_index = 500 * np.arange(16, 41)
time_index=[1000]
time_index = 100 * np.arange(0, 101)


levels = [9, 9]

rho0 = 1e-6

Nx = 2**levels[0]
Ny = 2**levels[0]
######################################################

for i in time_index:

    print("i = ", i)

    filename="output-%05d-2D-0.raw" % i

    sim_data = np.loadtxt(base_dir + filename).T
    time, ix, iy, x, y, density, pressure, vx, vy, totalE, radE = sim_data

    reshaped_sim_data = [dat.reshape(Nx, Ny) for dat in sim_data]
    time, ix, iy, x, y, density, pressure, vx, vy, totalE, radE = reshaped_sim_data

    # 
    plt.figure(figsize=(6, 6))

    cs = plt.contourf(x / 1e7, y / 1e7,
                    density / rho0,
                    levels=np.linspace(0.9, 2.1, 200),
                    extend='both')
    cbar = plt.colorbar(cs, ticks=[1, 2], shrink=0.4, pad=0.12,
                        label=r'$\rho/\rho_0$', orientation='horizontal')
    
    # cs = plt.contourf(x / 1e7, y / 1e7,
                    # pressure + radE/3)
    # cs = plt.contourf(x / 1e7, y / 1e7, totalE)
    # cs = plt.contourf(x / 1e7, y / 1e7, radE)
    # cbar = plt.colorbar(cs)

    plt.xlabel(r'$x / L$', fontsize=11)
    plt.ylabel(r'$y / L$', fontsize=11)

    plt.title(r'$t = %.2f$' % time[0,0], fontsize=12, pad=15)

    # plt.ylim(0.5, 1.0)
    plt.gca().set_aspect('equal')

    plt.tight_layout()

    plt.savefig("RadKHI_GPU_%05d.png" % i, dpi=150, bbox_inches='tight')

    # plt.show()
    plt.close()
