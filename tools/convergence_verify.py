import glob
import os
from collections.abc import Callable

import matplotlib.pyplot as plt
import numpy as np
from acoustic_solution import Acoustic
from verify_lib import (compute_l1_error_fvm, find_last_output, parse_cli,
                        parse_config, wrapFunction)


def main() -> None:
    """
    Test convergence of the code
    """

    # Slope tolerance
    tol = 1e-1

    # Convergence order
    ord = 3

    # Get all files
    yaml_file, out_dir, _, make_plot = parse_cli(get_file=False)
    problem, gamma, x0, x1, problem_dict = parse_config(yaml_file)

    # Find every example and build the dx and L1 error arrays
    if out_dir is None:
        pattern = "output_*"
    else:
        pattern = os.path.join(out_dir, "output_*")

    dirs = glob.glob(pattern)
    dirs.sort(key=lambda x: int(x.split("_")[-1]))
    dx = []
    l1err = []

    acoustic_instance: wrapFunction
    u_exact: Callable

    first_loop = True
    for dir in dirs:
        last_output = find_last_output(dir=dir)
        assert last_output is not None

        out_tuple = np.loadtxt(last_output, usecols=(0, 2, 3, 4, 5)).T

        # Extract physical quantities from tuple
        t_arr, x_num, rho_num, p_num, u_num = out_tuple
        time = t_arr[0]

        # Instantiate our solution class in the first loop
        if first_loop:
            acoustic_instance = wrapFunction(
                Acoustic(gamma, x0, x1, problem_dict), time)
            u_exact = acoustic_instance.velocity

            first_loop = False

        dx.append(x_num[1] - x_num[0])
        l1err.append(compute_l1_error_fvm(x_num, u_num, u_exact))

    if make_plot:
        def plot_order(dx, l1err, ord):
            plt.plot(dx, np.array(dx) ** ord *
                     l1err[0] / dx[0] ** ord, "--", label=f"order {ord}")

        plt.plot(dx, l1err, "o-")
        plot_order(dx, l1err, 3)
        plot_order(dx, l1err, 2)
        plot_order(dx, l1err, 1)

        plt.xscale("log")
        plt.yscale("log")

        plt.xlabel("dx")
        plt.ylabel("L1 error")
        plt.title("Error convergence")

        plt.legend()
        plt.savefig("convergence.pdf")

    # Find slope of logarithms for convergence test
    logE = np.log(l1err)
    logDX = np.log(dx)
    slope = (logE[1:] - logE[:-1]) / (logDX[1:] - logDX[:-1])

    # Assert that slope is expected order
    try:
        assert abs(abs(np.mean(slope)) - ord) < tol
    except AssertionError:
        s = f"Mean of slope is {np.mean(slope)}"
        print(s)
        raise


if __name__ == "__main__":
    main()
