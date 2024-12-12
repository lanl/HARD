
import numpy as np
import matplotlib.pyplot as plt

base_dir = "/vast/home/yoonsoo/flastro/build_catalyst/"

t_f = 4.0
v = 0.2

def analytic_solution(x, t):
    return 1.0 + 0.2 * np.sin(2.0 * np.pi * (x - 0.2 * t))


error_norm = []
levels = np.array([4, 5, 6, 7, 8, 9])

for level in levels:
    file_name = base_dir + "lev%d.raw" % level
    data = np.loadtxt(file_name).T
    time, cell_idx, coord_x, density, pressure, velocity, totalE, radE = data

    error = density - analytic_solution(coord_x, time[3])
    error_norm += [np.sqrt(np.sum(error**2) / len(time))]

error_norm = np.array(error_norm)

# Make plot
plt.figure(figsize=(5, 5), facecolor='w', dpi=100)

plt.plot(levels, error_norm, '-o')

lev_dense = np.linspace(3.5, 9.5)
err_dense = 0.8 / (2**lev_dense)**2
plt.plot(lev_dense, err_dense, '--', color='k', linewidth=0.6, label='2nd order')

plt.xlabel('Grid level')
plt.ylabel('Error')

plt.yscale('log')

plt.legend()

plt.savefig('Conv-SineWave-New.png', dpi=120, bbox_inches='tight')
# plt.show()

print(np.log2(error_norm[:-1] / error_norm[1:]))
