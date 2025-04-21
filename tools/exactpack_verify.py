import yaml
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import sys
import os
import glob
import re


# Threshold for passing L2 error tests
tolerance = 5.0e-1


def color_text(text, status):
    if status == "PASS":
        return f"\033[92m{text}\033[0m"  # green
    elif status == "FAIL":
        return f"\033[91m{text}\033[0m"  # red
    else:
        return text


def find_latest_output(pattern="output-*-1D-0.raw"):
    files = glob.glob(pattern)
    if not files:
        print("No matching output files found.")
        return None
    # Extract numeric ID from each filename
    files_with_step = []
    for f in files:
        match = re.search(r"output-(\d+)-1D-0\.raw", f)
        if match:
            step = int(match.group(1))
            files_with_step.append((step, f))
    if not files_with_step:
        print("No valid files matching pattern.")
        return None
    # Return file with max step number
    return max(files_with_step)[1]


def parse_config(yaml_file):
    with open(yaml_file, 'r') as f:
        config = yaml.safe_load(f)
    gamma = float(config.get('gamma', 1.4))
    x0 = float(config['coords'][0][0])
    x1 = float(config['coords'][1][0])
    problem = config['problem']
    return problem, gamma, x0, x1


def read_raw_output(raw_file):
    with open(raw_file, 'r') as f:
        lines = f.readlines()

    data = []
    for line in lines:
        if line.startswith('#') or not line.strip():
            continue
        values = list(map(float, line.strip().split()))
        data.append(values)

    arr = np.array(data)
    x = arr[:, 2]
    dx = abs(x[0]-x[1])
    density = arr[:, 3]
    pressure = arr[:, 4]
    velocity = arr[:, 5]
    time = arr[0, 0]
    return time, x, density, pressure, velocity, dx


def sedov_analytic_solution(x, t, gamma):
    from exactpack.solvers.sedov import Sedov
    solver = Sedov(gamma=gamma, geometry=1, eblast=0.0673185)
    result = solver(x, t)
    return result.density, result.pressure, result.velocity


def riemann_analytic_solution(x, t, gamma, left_state, right_state):
    try:
        from exactpack.solvers.riemann.ep_riemann import IGEOS_Solver
    except ImportError:
        print("Error: ExactPack Riemann solver not found.")
        sys.exit(1)

    rho_l, u_l, p_l = left_state
    rho_r, u_r, p_r = right_state

    solver = IGEOS_Solver(
        rl=rho_l, ul=u_l, pl=p_l, gl=gamma,
        rr=rho_r, ur=u_r, pr=p_r, gr=gamma,
        xmin=min(x), xd0=0.5 * (min(x) + max(x)), xmax=max(x), t=t
    )

    sol = solver(x, t)
    return sol.density, sol.pressure, sol.velocity


def compute_l2_error(numerical, analytical, dx):
    return np.sqrt(dx * np.sum((numerical - analytical)** 2 ))


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
    raw_file = None
    make_plot = False

    for arg in sys.argv[2:]:
        if arg.endswith(".raw"):
            raw_file = arg
        elif arg == "--plot":
            make_plot = True

    # If no file is passed, select the latest available output
    if raw_file is None:
        raw_file = find_latest_output()
        if raw_file is None:
            print("No raw file found and none provided.")
            sys.exit(1)
        else:
            print(f"Auto-selected input file: {raw_file}")

    problem, gamma, x0, x1 = parse_config(yaml_file)
    time, x_num, rho_num, p_num, u_num, dx = read_raw_output(raw_file)

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

    else:
        print(f"Unsupported problem type '{problem}'")
        sys.exit(1)

    rho_interp = interp1d(x_analytic, rho_exact, bounds_error=False,
                          fill_value="extrapolate")
    p_interp = interp1d(x_analytic, p_exact, bounds_error=False,
                        fill_value="extrapolate")
    u_interp = interp1d(x_analytic, u_exact, bounds_error=False,
                        fill_value="extrapolate")

    rho_ref = rho_interp(x_num)
    p_ref = p_interp(x_num)
    u_ref = u_interp(x_num)

    if make_plot:
        tag = os.path.splitext(os.path.basename(raw_file))[0].replace(
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
    err_p = compute_l2_error(p_num, p_ref, dx)
    err_u = compute_l2_error(u_num, u_ref, dx)

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
        print("Test FAILED: L2 error exceeds tolerance of "+str(tolerance))
        sys.exit(1)
    else:
        print("Test PASSED")
        sys.exit(0)


if __name__ == "__main__":
    main()
