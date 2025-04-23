import glob
import os

import numpy as np
from acoustic_solution import Acoustic
from numpy.typing import NDArray


def l2(analytical: NDArray | float, numerical: NDArray | float, dx: float,
       divide: bool = True) -> float:
    if divide:
        return np.sqrt(dx *
                       np.sum(((numerical - analytical) / analytical) ** 2))
    else:
        return np.sqrt(dx * np.sum((numerical - analytical) ** 2))


def main() -> None:
    """
    Print the acoustic wave output errors
    """

    gamma = 1.4
    r0 = 1.0
    p0 = 1.0

    # Get output files
    output_dir = "../build/output/"
    files = glob.glob(os.path.join(output_dir, "*.raw"))
    files.sort()

    # Extract data from first file for the solution
    t_arr, x, r, p, u = np.loadtxt(files[0], usecols=(0, 2, 3, 4, 5)).T
    acoustic = Acoustic(x, gamma, r0, p0, r, p, u)
    dx = x[1] - x[0]

    # Get last file data for the errors
    t_arr, x, r, p, u = np.loadtxt(files[-1], usecols=(0, 2, 3, 4, 5)).T
    t = t_arr[0]
    solution = acoustic(x, t)

    print("Relative L2 errors for:\n")
    print(f"Density: {l2(solution.density, r, dx)}")
    print(f"Pressure: {l2(solution.pressure, p, dx)}")
    print(f"Velocity: {l2(solution.velocity, u, dx, divide=False)}")
    print()
    print(f"dx = {dx}")
    print(f"tf = {t}")


if __name__ == "__main__":
    main()
