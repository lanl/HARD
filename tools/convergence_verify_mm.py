import glob
import os
from collections.abc import Callable
from importlib import util as imputil
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FixedFormatter, FixedLocator, NullLocator
from contact_move_solution import ContactMove
from verify_lib import (
    compute_l1_error_fvm,
    find_last_output,
    parse_cli,
    wrapFunction,
)


def _parse_mm_config(
    file_path: str,
) -> tuple[str, float, float, np.ndarray, np.ndarray, dict[str, Any]]:
    """
    Parse a multi-material config file for CONTACT_move verification.
    """

    module_name = "contact_move_config"
    spec = imputil.spec_from_file_location(module_name, file_path)
    assert spec is not None
    assert spec.loader is not None

    module = imputil.module_from_spec(spec)
    spec.loader.exec_module(module)

    config = module.config
    eos_list = config.get("eos_list", [])
    if len(eos_list) < 2:
        raise ValueError("CONTACT_move config must define at least two EOS entries.")

    gamma0 = float(eos_list[0].get("gamma", 1.4))
    gamma1 = float(eos_list[1].get("gamma", 1.4))
    x0 = np.array(config["coords"][0])
    x1 = np.array(config["coords"][1])
    problem = str(config["problem"])
    problem_dict = config.get("problem_parameters", {})

    return problem, gamma0, gamma1, x0, x1, problem_dict


def _compute_masked_l1_error(
    x_num: np.ndarray,
    rho_num: np.ndarray,
    rho_exact: Callable,
    mask: np.ndarray,
) -> float:
    """
    Compute a normalized L1 error on a masked subset of the domain.
    """

    x_sel = x_num[mask]
    rho_sel = rho_num[mask]

    if x_sel.size == 0:
        raise ValueError("Smooth-region mask selected no cells.")

    err = np.abs(rho_sel - rho_exact(x_sel))
    return float(np.mean(err))


def _save_convergence_plot(
    grid_sizes: np.ndarray,
    full_err: np.ndarray,
    smooth_err: np.ndarray,
    selected_err: np.ndarray,
    ord: int,
    dim: int,
    field_name: str,
) -> None:
    fig, ax = plt.subplots(figsize=(7, 5))

    def plot_order(n_vals, err_vals, order, color: str | None = None):
        n_arr = np.asarray(n_vals, dtype=float)
        ref = err_vals[0] * (n_arr[0] / n_arr) ** order
        if color is None:
            ax.plot(
                n_arr,
                ref,
                "--",
                label=f"order {order}",
                linewidth=1.5,
            )
        else:
            ax.plot(
                n_arr,
                ref,
                "--",
                label=f"order {order}",
                linewidth=1.5,
                color=color,
            )

    colors_6 = [
        "#000000",
        "#c1272d",
        "#0000a7",
        "#eecc16",
        "#008176",
        "#b3b3b3",
    ]

    ax.plot(
        grid_sizes,
        full_err,
        "o-",
        linewidth=2,
        color=colors_6[0],
        label="full",
    )
    ax.plot(
        grid_sizes,
        smooth_err,
        "s-",
        linewidth=2,
        color=colors_6[4],
        label="smooth",
    )
    plot_order(grid_sizes, selected_err, 0.5, color=colors_6[1])
    plot_order(grid_sizes, selected_err, 1, color=colors_6[2])
    plot_order(grid_sizes, selected_err, 2, color=colors_6[3])
    plot_order(grid_sizes, selected_err, 3, color=colors_6[4])
    plot_order(grid_sizes, selected_err, 4)

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlabel("N")
    ax.set_ylabel(f"{field_name.capitalize()} L1 error (-)")

    ax.xaxis.set_major_locator(FixedLocator(grid_sizes.astype(float)))
    ax.xaxis.set_major_formatter(FixedFormatter([str(val) for val in grid_sizes]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.offsetText.set_visible(False)
    ax.set_xlim(grid_sizes[0] * 0.95, grid_sizes[-1] * 1.05)

    ax.tick_params(
        axis="x", which="major", top=True, bottom=True, labeltop=False, pad=6, labelsize=8
    )
    ax.tick_params(axis="x", which="minor", bottom=False, top=False, labelbottom=False)
    ax.tick_params(axis="y", which="major", right=True, labelright=False)
    ax.tick_params(axis="y", which="minor", right=True)

    ax.legend()
    fig.tight_layout()

    plot_path = Path.cwd() / f"contact_move_order_{ord}_convergence-{dim}D_{field_name}.pdf"
    fig.savefig(plot_path, bbox_inches="tight")
    plt.close(fig)
    print(f"\n   --> saved {field_name} plot: {plot_path}")


def main() -> None:
    """
    Test convergence of the CONTACT_move problem using density and pressure L1 errors.
    """

    tol = 0.1
    target_time = 0.1

    config_file, dim, out_dir, _, _, make_plot = parse_cli(get_file=False)
    problem, gamma0, gamma1, x0, x1, problem_dict = _parse_mm_config(config_file)

    ord = int(problem_dict["convergence_order"])
    exclusion_width = float(problem_dict.get("exclusion_width", 0.1))
    convergence_metric = str(problem_dict.get("convergence_metric", "full")).lower()

    if out_dir is None:
        pattern = "output_*"
    else:
        pattern = os.path.join(out_dir, "output_*")

    dirs = [dirpath for dirpath in glob.glob(pattern) if os.path.isdir(dirpath)]
    dirs.sort(key=lambda x: int(os.path.basename(x).split("_")[-1]))

    print("\n   --> using output directories:")
    for dirpath in dirs:
        print(f"       {dirpath}")
    dx = []
    rho_l1err_full = []
    rho_l1err_smooth = []
    p_l1err_full = []
    p_l1err_smooth = []
    u_l1err_full = []
    u_l1err_smooth = []

    cont_move_instance: wrapFunction
    rho_exact: Callable
    p_exact: Callable
    u_exact: Callable
    interface_position = None

    first_loop = True
    for dir in dirs:
        last_output, combined = find_last_output(dir=dir)
        assert last_output is not None
        print(f"\n   --> reading file: {last_output}")

        out_tuple = np.loadtxt(
            last_output,
            comments=["time", "#"],
            delimiter=",",
            usecols=(0, 2, 3, 4, 10),
        ).T

        if combined:
            os.remove(last_output)

        t_arr, x_num, rho_num, p_num, u_num = out_tuple
        time = t_arr[0]

        if not np.isclose(time, target_time):
            raise AssertionError(
                f"Expected CONTACT_move output at t={target_time}, found t={time}."
            )

        if first_loop:
            cont_move = ContactMove(gamma0, gamma1, x0, x1, problem_dict, dim=dim)
            cont_move_instance = wrapFunction(
                cont_move,
                time,
                ["density", "pressure", "velocity", "material"],
            )
            rho_exact = cont_move_instance.density
            p_exact = cont_move_instance.pressure
            u_exact = cont_move_instance.velocity
            interface_position = cont_move._interface_position(time)
            first_loop = False

        assert interface_position is not None

        smooth_mask = np.abs(x_num - interface_position) > exclusion_width

        current_dx = x_num[1] - x_num[0]
        current_rho_l1_full = compute_l1_error_fvm(x_num, rho_num, rho_exact, dim)
        current_rho_l1_smooth = _compute_masked_l1_error(
            x_num, rho_num, rho_exact, smooth_mask
        )
        current_p_l1_full = compute_l1_error_fvm(x_num, p_num, p_exact, dim)
        current_p_l1_smooth = _compute_masked_l1_error(
            x_num, p_num, p_exact, smooth_mask
        )
        current_u_l1_full = compute_l1_error_fvm(x_num, u_num, u_exact, dim)
        current_u_l1_smooth = _compute_masked_l1_error(
            x_num, u_num, u_exact, smooth_mask
        )

        dx.append(current_dx)
        rho_l1err_full.append(current_rho_l1_full)
        rho_l1err_smooth.append(current_rho_l1_smooth)
        p_l1err_full.append(current_p_l1_full)
        p_l1err_smooth.append(current_p_l1_smooth)
        u_l1err_full.append(current_u_l1_full)
        u_l1err_smooth.append(current_u_l1_smooth)

        # print(f"   --> dx = {current_dx}")
        # print(f"   --> density full L1 = {current_rho_l1_full}")
        # print(f"   --> density smooth L1 = {current_rho_l1_smooth}")
        # print(f"   --> pressure full L1 = {current_p_l1_full}")
        # print(f"   --> pressure smooth L1 = {current_p_l1_smooth}")
        # print(f"   --> velocity full L1 = {current_u_l1_full}")
        # print(f"   --> velocity smooth L1 = {current_u_l1_smooth}")

    l1err = rho_l1err_full if convergence_metric == "full" else rho_l1err_smooth

    if make_plot:

        dx_arr = np.asarray(dx, dtype=float)
        grid_sizes = np.rint(1.0 / dx_arr).astype(int)
        sort_idx = np.argsort(grid_sizes)
        grid_sizes = grid_sizes[sort_idx]

        rho_l1err_full_arr = np.asarray(rho_l1err_full, dtype=float)[sort_idx]
        rho_l1err_smooth_arr = np.asarray(rho_l1err_smooth, dtype=float)[sort_idx]
        rho_l1err_arr = np.asarray(l1err, dtype=float)[sort_idx]

        p_l1err_full_arr = np.asarray(p_l1err_full, dtype=float)[sort_idx]
        p_l1err_smooth_arr = np.asarray(p_l1err_smooth, dtype=float)[sort_idx]
        p_l1err_arr = (
            np.asarray(p_l1err_full, dtype=float)[sort_idx]
            if convergence_metric == "full"
            else np.asarray(p_l1err_smooth, dtype=float)[sort_idx]
        )

        u_l1err_full_arr = np.asarray(u_l1err_full, dtype=float)[sort_idx]
        u_l1err_smooth_arr = np.asarray(u_l1err_smooth, dtype=float)[sort_idx]
        u_l1err_arr = (
            np.asarray(u_l1err_full, dtype=float)[sort_idx]
            if convergence_metric == "full"
            else np.asarray(u_l1err_smooth, dtype=float)[sort_idx]
        )

        _save_convergence_plot(
            grid_sizes,
            rho_l1err_full_arr,
            rho_l1err_smooth_arr,
            rho_l1err_arr,
            ord,
            dim,
            "density",
        )
        _save_convergence_plot(
            grid_sizes,
            p_l1err_full_arr,
            p_l1err_smooth_arr,
            p_l1err_arr,
            ord,
            dim,
            "pressure",
        )
        _save_convergence_plot(
            grid_sizes,
            u_l1err_full_arr,
            u_l1err_smooth_arr,
            u_l1err_arr,
            ord,
            dim,
            "velocity",
        )

    logE = np.log(l1err)
    logDX = np.log(dx)
    slope = (logE[1:] - logE[:-1]) / (logDX[1:] - logDX[:-1])

    min_expected_slope = 0.5
    mean_slope = abs(float(np.mean(slope)))

    try:
        print("\n   --> minimum required density slope = ", min_expected_slope)
        print("\n   --> measured density slope = ", mean_slope)
        assert mean_slope >= min_expected_slope
    except AssertionError:
        s = (
            f"\n   --> mean slope for CONTACT_move density ({convergence_metric}) is "
            f"{mean_slope}"
        )
        print(s)
        print(f"\n   --> density full-domain L1 errors: {rho_l1err_full}")
        print(f"\n   --> density smooth-region L1 errors: {rho_l1err_smooth}")
        print(f"\n   --> pressure full-domain L1 errors: {p_l1err_full}")
        print(f"\n   --> pressure smooth-region L1 errors: {p_l1err_smooth}")
        print(f"\n   --> velocity full-domain L1 errors: {u_l1err_full}")
        print(f"\n   --> velocity smooth-region L1 errors: {u_l1err_smooth}\n")
        raise


if __name__ == "__main__":
    main()
