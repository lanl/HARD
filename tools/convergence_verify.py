import glob
import os
from collections.abc import Callable

import matplotlib.pyplot as plt
import numpy as np
from acoustic_solution import Acoustic
from verify_lib import (compute_l2_error, find_last_output, parse_cli,
                        parse_config, wrapFunction)


def main() -> None:
    """
    Test convergence of the code
    """

    # Slope tolerance
    tol = 1e-1

    # Convergence order
    ord = 2

    # Get all files
    yaml_file, out_dir, _, make_plot = parse_cli(get_file=False)
    problem, gamma, x0, x1, problem_dict = parse_config(yaml_file)

    # Find every example and build the dx and L2 error arrays
    if out_dir is None:
        pattern = "output_*"
    else:
        pattern = os.path.join(out_dir, "output_*")

    dirs = glob.glob(pattern)
    dirs.sort()
    dx = []
    l2err = []

    acoustic_instance: wrapFunction
    rho_exact: Callable

    first_loop = True
    for dir in dirs:
        last_output = find_last_output(dir=dir)
        assert last_output is not None

        out_tuple = np.loadtxt(last_output, usecols=(0, 2, 3, 4, 5)).T

        # Extract physical quantities from tuple
        t_arr, x_num, rho_num, p_num, u_num = out_tuple
        time = t_arr[0]
        dx.append(x_num[1] - x_num[0])

        # Instantiate our solution class in the first loop
        if first_loop:
            acoustic_instance = wrapFunction(
                Acoustic(gamma, x0, x1, problem_dict), time)
            rho_exact = acoustic_instance.density

            first_loop = False

        l2err.append(compute_l2_error(x_num, rho_num, rho_exact))

    # Find slope of logarithms for convergence test
    logE = np.log(l2err)
    logDX = np.log(dx)
    slope = (logE[1:] - logE[:-1]) / (logDX[1:] - logDX[:-1])

    if make_plot:
        def plot_order(dx, l2err, ord):
            plt.plot(dx, np.array(dx) ** ord *
                     l2err[0] / dx[0] ** ord, "--", label=f"order {ord}")

        plt.plot(dx, l2err, "o-")
        plot_order(dx, l2err, 3)
        plot_order(dx, l2err, 2)
        plot_order(dx, l2err, 1)

        plt.xlabel("dx")
        plt.ylabel("l2-error")
        plt.xscale("log")
        plt.yscale("log")
        plt.legend()
        plt.savefig("convergence.pdf")

    # Assert that slope is expected order
    try:
        assert abs(abs(np.mean(slope)) - ord) < tol
    except AssertionError:
        s = f"Mean of slope is {np.mean(slope)}"
        print(s)
        raise


if __name__ == "__main__":
    main()
