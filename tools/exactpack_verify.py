import os
import sys
from collections import namedtuple
from collections.abc import Callable

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray
from verify_lib import (acoustic_analytic_solution, compute_l2_error,
                        parse_cli, parse_config, wrapFunction)

# Threshold for passing L2 error tests
Tolerance = namedtuple("Tolerance", ["rho", "p", "u"])
tolerance = Tolerance(1e-2, 1e-2, 5e-1)


def color_text(text: str, status: str) -> str:
    if status == "PASS":
        return f"\033[92m{text}\033[0m"  # green
    elif status == "FAIL":
        return f"\033[91m{text}\033[0m"  # red
    else:
        return text


def sedov_analytic_solution(x: NDArray, t: float, gamma: float) -> tuple[
    Callable,
    Callable,
    Callable
]:
    from exactpack.solvers.sedov import Sedov

    solver = Sedov(gamma=gamma, geometry=1, eblast=0.0673185)
    result = wrapFunction(solver, t)

    return result.density, result.pressure, result.velocity


def riemann_analytic_solution(x: NDArray, t: float, gamma: float,
                              left_state: tuple[float, float, float],
                              right_state: tuple[float, float, float]
                              ) -> tuple[Callable, Callable, Callable]:

    try:
        from exactpack.solvers.riemann.ep_riemann import IGEOS_Solver
    except ImportError:
        sys.exit("Error: ExactPack Riemann solver not found.")

    rho_l, u_l, p_l = left_state
    rho_r, u_r, p_r = right_state

    solver = IGEOS_Solver(
        rl=rho_l, ul=u_l, pl=p_l, gl=gamma,
        rr=rho_r, ur=u_r, pr=p_r, gr=gamma,
        xmin=min(x), xd0=0.5 * (min(x) + max(x)), xmax=max(x), t=t
    )

    sol = wrapFunction(solver, t)
    return sol.density, sol.pressure, sol.velocity


def plot_comparison(x: NDArray, num_vals: NDArray, exact_vals: NDArray,
                    quantity: str, time: float, tag: str, problem: str) -> None:
    plt.figure()
    plt.plot(x, exact_vals, label="Analytic", linestyle="--")
    plt.plot(x, num_vals, label="Simulation", marker='o',
             linestyle='none', markersize=4)
    plt.xlabel("x")
    plt.ylabel(quantity)
    plt.title(f"{quantity} at t = {time:.4f}")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    filename = f"{problem}_{quantity.lower()}_comparison_{tag}.png"
    plt.savefig(filename)
    plt.close()
    print(f"Saved plot: {filename}")


def main() -> None:

    # Get the values
    yaml_file, out_dir, raw_file, make_plot = parse_cli()
    assert raw_file is not None

    problem, gamma, x0, x1, problem_dict = parse_config(yaml_file)

    # Read in the last raw file
    out_tuple = np.loadtxt(raw_file, usecols=(0, 2, 3, 4, 5)).T

    # Extract physical quantities from tuple
    t_arr, x_arr, rho_num, p_num, u_num = out_tuple
    time = t_arr[0]

    if problem == "sod":
        left_state = (1.0, 0.0, 1.0)
        right_state = (0.125, 0.0, 0.1)
        rho_exact, p_exact, u_exact = riemann_analytic_solution(
            x_arr, time, gamma, left_state, right_state)

    elif problem == "leblanc":
        left_state = (1.0, 0.0, 0.1)
        right_state = (1e-3, 0.0, 1e-10)
        rho_exact, p_exact, u_exact = riemann_analytic_solution(
            x_arr, time, gamma, left_state, right_state)

    elif problem == "rankine-hugoniot":
        left_state = (1.0, 0.0, 1.0)
        right_state = (0.25, 0.0, 0.1795)
        rho_exact, p_exact, u_exact = riemann_analytic_solution(
            x_arr, time, gamma, left_state, right_state)

    elif problem == "sedov":
        rho_exact, p_exact, u_exact = sedov_analytic_solution(
            x_arr, time, gamma)

    elif problem == "acoustic-wave":
        rho_exact, p_exact, u_exact = acoustic_analytic_solution(
            gamma, time, x0, x1, problem_dict)

    else:
        sys.exit(f"Unsupported problem type '{problem}'")

    rho_ref = rho_exact(x_arr)
    p_ref = p_exact(x_arr)
    u_ref = u_exact(x_arr)

    if make_plot:
        assert raw_file is not None
        tag = os.path.splitext(os.path.basename(raw_file))[0].replace(
            "output-", "").replace(".raw", "")
        plot_comparison(x_arr, rho_num, rho_ref, "Density", time, tag, problem)
        plot_comparison(x_arr, p_num, p_ref, "Pressure", time, tag, problem)
        plot_comparison(x_arr, u_num, u_ref, "Velocity", time, tag, problem)

    # For Sedov: only compare for x > 0.1
    if problem == "sedov":
        sedov_cutoff = 0.2
        mask = x_arr > sedov_cutoff
        x_arr = x_arr[mask]
        rho_num = rho_num[mask]
        p_num = p_num[mask]
        u_num = u_num[mask]
        rho_ref = rho_ref[mask]
        p_ref = p_ref[mask]
        u_ref = u_ref[mask]

    err_rho = compute_l2_error(x_arr, rho_num, rho_exact)
    err_p = compute_l2_error(x_arr, p_num, p_exact)
    err_u = compute_l2_error(x_arr, u_num, u_exact)

    status_rho = "PASS" if err_rho <= tolerance.rho else "FAIL"
    status_p = "PASS" if err_p <= tolerance.p else "FAIL"
    status_u = "PASS" if err_u <= tolerance.u else "FAIL"

    print(f"Problem type: {problem}")
    print(f"Compared solution at t = {time:.4f}")
    print("L2 Errors and Status:")
    print(f"  Density  : {err_rho:.6e} \
    [{color_text(status_rho, status_rho)}]")
    print(f"  Pressure : {err_p:.6e} \
    [{color_text(status_p, status_p)}]")
    print(f"  Velocity : {err_u:.6e} \
    [{color_text(status_u, status_u)}]")

    if any(e > t for e, t in zip((err_rho, err_p, err_u),
                                 (tolerance.rho, tolerance.p, tolerance.u))):

        s = "Test FAILED: relative L2 error exceeds tolerance of"
        s += f" {tolerance:.2e}"
        sys.exit(s)
    else:
        print("Test PASSED")


if __name__ == "__main__":
    main()
