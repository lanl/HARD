
import numpy as np
import matplotlib.pyplot as plt

from matplotlib import cm, ticker

base_dir = "/vast/home/yoonsoo/flastro/build/"
# base_dir = "/vast/home/yoonsoo/cdss_data/Shadow2D/"
# base_dir = "/vast/home/yoonsoo/cdss_data/Diffusion2D/"

def diffusion_2d_solution(x,y,t):
    return 2.0 + np.exp(-8.0 * np.pi**2 * t) * np.sin(2*np.pi*x) * np.sin(2*np.pi*y)

lev = 7
x0, y0 = 0.4, 0.5
cylinder_radius = 0.15

mpi_rank = 20

for i in 100 * np.arange(0, 101):
    print("i = ", i)

    file_name = base_dir + "output-%05d-2D-0.raw" % i
    data = np.loadtxt(file_name).T

    time, ix, iy, x, y, density, pressure, vx, vy, totalE, radE = data

    # Make plot
    plt.figure(figsize=(6, 5), facecolor='w', dpi=100)

    plt.title("t = %.2g" % time[0])
        
    Nx, Ny = 2**lev + 1, 2**lev + 1

    # Z = density.reshape(Nx, Ny)
    Z = radE.reshape(Nx, Ny)
    Zmin = np.min(Z)
    Zmax = 0.9 * np.max(Z)
    # Zmin, Zmax = 0.03, 1.0

    cs = plt.contourf(x.reshape(Nx, Ny), y.reshape(Nx, Ny), Z,
                    #   levels=np.logspace(-1.7, 0, 51), locator=ticker.LogLocator(), 
                    levels = np.linspace(Zmin, Zmax, 101),
                extend='both')
    plt.contour(x.reshape(Nx, Ny), y.reshape(Nx, Ny), Z,
                levels=np.linspace(Zmin, Zmax, 15), colors='w', linewidths=0.4)
    
    # residual = radE - diffusion_2d_solution(x, y, time)
    # cs = plt.contourf(x.reshape(Nx, Ny), y.reshape(Nx, Ny), np.abs(residual.reshape(Nx, Ny)),
                    #   levels=np.linspace(0, 1e-2), extend='both')

    # cbar = plt.gcf().colorbar(cs, shrink=0.6)
    # cbar.set_ticks([1e-2, 1e-1, 1.0])
    plt.gca().set_aspect('equal')

    circle1 = plt.Circle((x0, y0), cylinder_radius, facecolor='gray', zorder=2)
    plt.gca().add_patch(circle1)

    plt.savefig("Shadow2D_%04d.png" % i, dpi=150, bbox_inches='tight')
    # plt.show()
    plt.close()
