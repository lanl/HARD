
import numpy as np
import matplotlib.pyplot as plt

base_dir = "/vast/home/yoonsoo/flastro/build/"


# Make plot
plt.figure(figsize=(5, 4), facecolor='w', dpi=100)

x = [0.01000, 0.10000, 0.17783,0.31623,0.45000,0.50000,0.56234,0.75000,1.00000,1.33352, 1.77828, 3.16228, 5.62341]

dsol_tau_01 = [0.09403, 0.09326, 0.09128, 0.08230, 0.06086, 0.04766, 0.03171, 0.00755, 0.00064]
tsol_tau_01 = [0.09531, 0.09531, 0.09532, 0.09529, 0.0882, 0.0476, 0.0037]

dsol_tau_1 = [0.5035, 0.4971, 0.4830, 0.4374, 0.3665, 0.3327, 0.2902, 0.1887, 0.1015, 0.0406, 0.01011, 0.00003]
tsol_tau_1 = [0.6430, 0.63585, 0.6195, 0.5618, 0.4471, 0.3580, 0.2537, 0.1143, 0.03648, 0.00291]
              
dsol_tau_10 = [1.86, 1.854, 1.828, 1.748, 1.628, 1.5723, 1.500, 1.29758, 1.0601, 0.79696, 0.5298, 0.1218, 0.00445]
tsol_tau_10 = [2.2357, 2.21944, 2.1834, 2.0644, 1.8607, 1.7317, 1.5749, 1.2739, 0.9878, 0.70822, 0.45016, 0.09673, 0.00375]


for i in [10]:
    file_name = base_dir + "output-%05d-1D-0.raw" % i
    data = np.loadtxt(file_name).T
    time, cell_idx, coord_x, density, pressure, velocity, totalE, radE = data

    plt.plot(coord_x, radE, color='k')

    plt.plot(x[:-4], dsol_tau_01, 'o', label="Exact diffusion", color='C0', markersize=5)
    plt.plot(x[:-6], tsol_tau_01, 'x', label="Exact transport", color='C0', markersize=5)
    
    # plt.plot(x[:-1], dsol_tau_1, 'o', label="Exact diffusion", color='C0', markersize=5)
    # plt.plot(x[:-3], tsol_tau_1, 'x', label="Exact transport", color='C0', markersize=5)

    # plt.plot(x[:], dsol_tau_10, 'o', label="Exact diffusion", color='C0', markersize=5)
    # plt.plot(x[:], tsol_tau_10, 'x', label="Exact transport", color='C0', markersize=5)


plt.legend()

plt.xlabel(r'$x$')
plt.ylabel(r'$E$')

plt.xlim(0, 6.0)

plt.savefig("Su-Olson_tau01.png", dpi=120, bbox_inches='tight')

plt.show()

# plt.close()
