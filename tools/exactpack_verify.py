import glob
import itertools
import os
import sys
from collections.abc import Callable

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
import yaml
from acoustic_solution import Acoustic
from numpy.typing import NDArray

# Threshold for passing L2 error tests
tolerance = 1.0e-0


class wrapFunction(object):
    def __init__(self, solver: Callable, extract: list[str]) -> None:
        """
        Create the wrapper for each attribute
        """

        for name in extract:
            self.__doWrapping(solver, name)

    def __doWrapping(self, f: Callable, name: str) -> None:
        """
        Transform the ExactPack array tuple into a Callable tuple
        """

        if isinstance(f, Acoustic):
            def wrap(x, t):
                return f(x, t)._asdict()[name]
        else:
            def wrap(x, t):
                if type(x) in [float, int]:
                    return f([x], t)[name]
                else:
                    return f(x, t)[name]

        self.__dict__[name] = wrap


def color_text(text, status):
    if status == "PASS":
        return f"\033[92m{text}\033[0m"  # green
    elif status == "FAIL":
        return f"\033[91m{text}\033[0m"  # red
    else:
        return text


def first_and_last_output(pattern="output-*-1D-0.raw"):
    files = glob.glob(pattern)
    if not files:
        print("No matching output files found.")
        return None

    # Sort the files so they are in order
    files.sort()

    # Return first and last
    return files[0], files[-1]


def parse_config(yaml_file):
    with open(yaml_file, 'r') as f:
        config = yaml.safe_load(f)
    gamma = float(config.get('gamma', 1.4))
    x0 = float(config['coords'][0][0])
    x1 = float(config['coords'][1][0])
    problem = config['problem']
    return problem, gamma, x0, x1


def sedov_analytic_solution(x, t, gamma):
    from exactpack.solvers.sedov import Sedov

    names = ["density", "pressure", "velocity"]

    solver = Sedov(gamma=gamma, geometry=1, eblast=0.0673185)
    result = wrapFunction(solver, extract=names)

    return result.density, result.pressure, result.velocity


def riemann_analytic_solution(x, t, gamma, left_state, right_state):
    try:
        from exactpack.solvers.riemann.ep_riemann import IGEOS_Solver
    except ImportError:
        print("Error: ExactPack Riemann solver not found.")
        sys.exit(1)

    names = ["density", "pressure", "velocity"]

    rho_l, u_l, p_l = left_state
    rho_r, u_r, p_r = right_state

    solver = IGEOS_Solver(
        rl=rho_l, ul=u_l, pl=p_l, gl=gamma,
        rr=rho_r, ur=u_r, pr=p_r, gr=gamma,
        xmin=min(x), xd0=0.5 * (min(x) + max(x)), xmax=max(x), t=t
    )

    sol = wrapFunction(solver, extract=names)
    return sol.density, sol.pressure, sol.velocity


def acoustic_analytic_solution(x, x_analytic, t, gamma, r0, p0, init_r, init_p,
                               init_u) -> tuple[Callable, Callable, Callable]:
    """
    Return the acoustic analytical solution
    """

    names = ["density", "pressure", "velocity"]

    result = wrapFunction(
        Acoustic(x, gamma, r0, p0, init_r, init_p, init_u), extract=names)

    return result.density, result.pressure, result.velocity


def compute_l2_error(x_num: NDArray, numerical: NDArray,
                     analytical: Callable) -> float:

    # For every cell, calculate the "volume" integral
    dx = x_num[1] - x_num[0]
    error = 0
    for i, x in enumerate(x_num):
        x0 = x[i] - dx * 0.5
        x1 = x[i] + dx * 0.5
        error += (numerical[i] - quad(analytical, x0, x1) / dx) ** 2

    return np.sqrt(dx * error)


def plot_comparison(x_num, num_vals, x_exact, exact_vals,
                    quantity, time, tag):
    plt.figure()
    plt.plot(x_exact, exact_vals, label="Analytic", linestyle="--")
    plt.plot(x_num, num_vals, label="Simulation", marker='o',
             linestyle='none', markersize=4)
    plt.xlabel("x")
    plt.ylabel(quantity)
    plt.title(f"{quantity} at t = {time:.4f}")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    filename = f"{quantity.lower()}_comparison_{tag}.png"
    plt.savefig(filename)
    plt.close()
    print(f"Saved plot: {filename}")


def main():
    if len(sys.argv) < 2:
        print("Usage: python compare_solution.py sod.yaml \
        [output.raw] [--plot]")
        sys.exit(1)

    yaml_file = sys.argv[1]
    first_raw_file = None
    last_raw_file = None
    make_plot = False

    for arg in sys.argv[2:]:
        if arg.endswith(".raw"):
            last_raw_file = arg
        elif arg == "--plot":
            make_plot = True

    # If no file is passed, select the latest available output
    if first_raw_file is None:
        first_raw_file, last_raw_file = first_and_last_output()
        if None in [first_raw_file, last_raw_file]:
            print("No raw file found and none provided.")
            sys.exit(1)
        else:
            print(f"Auto-selected input file: {last_raw_file}")

    # Read in the yaml file and the last raw file
    problem, gamma, x0, x1 = parse_config(yaml_file)
    out_tuple = np.loadtxt(last_raw_file, usecols=(0, 2, 3, 4, 5)).T

    # Extract physical quantities from tuple
    t_arr, x_num, rho_num, p_num, u_num = out_tuple
    time = t_arr[0]
    dx = x_num[1] - x_num[0]

    x_analytic = np.linspace(x0, x1, 1000)

    if problem == "sod":
        left_state = (1.0, 0.0, 1.0)
        right_state = (0.125, 0.0, 0.1)
        rho_exact, p_exact, u_exact = riemann_analytic_solution(
            x_analytic, time, gamma, left_state, right_state)

    elif problem == "leblanc":
        left_state = (1.0, 0.0, 0.1)
        right_state = (1e-3, 0.0, 1e-10)
        rho_exact, p_exact, u_exact = riemann_analytic_solution(
            x_analytic, time, gamma, left_state, right_state)

    elif problem == "rankine-hugoniot":
        left_state = (1.0, 0.0, 1.0)
        right_state = (0.25, 0.0, 0.1795)
        rho_exact, p_exact, u_exact = riemann_analytic_solution(
            x_analytic, time, gamma, left_state, right_state)

    elif problem == "sedov":
        rho_exact, p_exact, u_exact = sedov_analytic_solution(
            x_analytic, time, gamma)

    elif problem == "acoustic-wave":
        # Find initial configuration
        out_tuple = np.loadtxt(first_raw_file, usecols=(0, 2, 3, 4, 5)).T
        t_arr, x0, rho0, p0, u0 = out_tuple

        rho_exact, p_exact, u_exact = acoustic_analytic_solution(
            x0, x_analytic, time, gamma, 1.0, 1.0, rho0, p0, u0)

    else:
        print(f"Unsupported problem type '{problem}'")
        sys.exit(1)

    rho_ref = rho_exact(x_num, time)
    p_ref = p_exact(x_num, time)
    u_ref = u_exact(x_num, time)

    if make_plot:
        tag = os.path.splitext(os.path.basename(last_raw_file))[0].replace(
            "output-", "").replace(".raw", "")
        plot_comparison(x_num, rho_num, x_analytic, rho_exact,
                        "Density", time, tag)
        plot_comparison(x_num, p_num, x_analytic, p_exact,
                        "Pressure", time, tag)
        plot_comparison(x_num, u_num, x_analytic, u_exact,
                        "Velocity", time, tag)

    # For Sedov: only compare for x > 0.1
    if problem == "sedov":
        sedov_cutoff = 0.2
        mask = x_num > sedov_cutoff
        x_num = x_num[mask]
        rho_num = rho_num[mask]
        p_num = p_num[mask]
        u_num = u_num[mask]
        rho_ref = rho_ref[mask]
        p_ref = p_ref[mask]
        u_ref = u_ref[mask]

    err_rho = compute_l2_error(rho_num, rho_ref, dx)
    if problem in ["leblanc", "sedov"]:
        # NOTE: The leblanc and sedov problems can have 0 pressure,
        # so we do not use the relative L2 error in this case
        err_p = compute_l2_error(p_num, p_ref, dx, divide=False)
    else:
        err_p = compute_l2_error(p_num, p_ref, dx)
    err_u = compute_l2_error(u_num, u_ref, dx, divide=False)

    status_rho = "PASS" if err_rho <= tolerance else "FAIL"
    status_p = "PASS" if err_p <= tolerance else "FAIL"
    status_u = "PASS" if err_u <= tolerance else "FAIL"

    print(f"Problem type: {problem}")
    print(f"Compared solution at t = {time:.4f}")
    print("L2 Errors and Status:")
    print(f"  Density  : {err_rho:.6e} \
    [{color_text(status_rho, status_rho)}]")
    print(f"  Pressure : {err_p:.6e} \
    [{color_text(status_p, status_p)}]")
    print(f"  Velocity : {err_u:.6e} \
    [{color_text(status_u, status_u)}]")

    if any(e > tolerance for e in (err_rho, err_p, err_u)):
        s = "Test FAILED: relative L2 error exceeds tolerance of"
        s += f" {tolerance:.2e}"
        print(s)
        sys.exit(1)
    else:
        print("Test PASSED")
        sys.exit(0)


if __name__ == "__main__":
    main()
