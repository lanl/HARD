import glob
import os

import matplotlib.pyplot as plt
import numpy as np
from acoustic_solution import Acoustic


def main() -> None:
    """
    Print the acoustic wave output L2 errors
    """

    gamma = 1.4
    r0 = 1.0
    p0 = 1.0

    # Get output files
    output_dir = "../build/output/"
    files = glob.glob(os.path.join(output_dir, "*.raw"))
    files.sort()

    first_file = True
    for file in files:
        t_arr, x, r, p, u = np.loadtxt(file, usecols=(0, 2, 3, 4, 5)).T
        t = t_arr[0]

        # Set up initial conditions
        if first_file:
            acoustic = Acoustic(x, gamma, r0, p0, r, p, u)
            first_file = False

        solution = acoustic(x, t)

        plt.plot(x, r, "b-", label="density")
        plt.plot(x, solution.density, "b--")

        plt.plot(x, p, "g-", label="pressure")
        plt.plot(x, solution.pressure, "g--")

        # plt.plot(x, u, "k-", label="velocity")
        # plt.plot(x, solution.velocity, "k--")

        plt.legend()

        figname = os.path.split(file)[-1].split('.')[0] + ".pdf"
        plt.savefig(os.path.join("figs", figname))
        plt.clf()


if __name__ == "__main__":
    main()
