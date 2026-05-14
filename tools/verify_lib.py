import glob
import os
import subprocess
import sys
from collections.abc import Callable
from importlib import util as imputil
from typing import Any

import numpy as np
from numpy.typing import NDArray


class wrapFunction(object):
    # Declaring the attributes to provide mypy with enough information
    density: Callable
    pressure: Callable
    velocity: Callable
    velocity_x: Callable
    velocity_y: Callable
    velocity_z: Callable

    def __init__(self, solver: Callable, t: float, extract: list[str]) -> None:
        """
        Create the wrapper for each attribute
        """

        for name in extract:
            self.__doWrapping(solver, t, name)

    def __doWrapping(self, f: Callable, t: float, name: str) -> None:
        """
        Transform the ExactPack array tuple into a Callable tuple
        """

        def wrap(coordinates):
            return f(coordinates, t)[name]

        self.__dict__[name] = wrap


def parse_cli(
    get_file: bool = True,
) -> tuple[str, int, str | None, str | None, bool, bool]:
    """
    Parse command line input
    """

    if len(sys.argv) < 2:
        s = f"Usage: python {sys.argv[0]} <config_file> <dim>"
        s += " [output_dir | output.csv]"
        s += " [--plot]"
        sys.exit(s)

    # Read the config file
    config_file = sys.argv[1]

    # Get the number of dimensions
    dim = int(sys.argv[2])
    assert dim in [1, 2, 3]

    csv_file = None
    combined = False
    out_dir = None
    make_plot = False

    for arg in sys.argv[3:]:
        if arg.endswith(".csv"):
            csv_file = arg
        elif arg == "--plot":
            make_plot = True
        else:
            out_dir = arg

    if get_file:
        # If no file is passed, select the latest available output
        if csv_file is None:
            pattern = f"output-{dim}D-?-*.csv"
            csv_file, combined = find_last_output(pattern=pattern, dir=out_dir)
            assert csv_file is not None

            print(f"Auto-selected input file: {csv_file}")

    return config_file, dim, out_dir, csv_file, combined, make_plot


def simple_quad(f: Callable, x0: NDArray, x1: NDArray, deg: int = 10) -> float:
    """
    Use a Gauss-Legendre quadrature of degree deg
    """

    points, weights = np.polynomial.legendre.leggauss(deg)

    def transform(x):
        return ((x1 - x0) * x + x1 + x0) * 0.5

    return np.sum(f(transform(points)) * weights) * 0.5 * (x1 - x0)


def get_dx(array: NDArray) -> float | None:
    """
    Get the minimum dx from array. Assume that the coordinates are ordered.
    """

    x0 = array[0]
    for x in array:
        dx = x - x0

        if dx > 0:
            return dx

    return None


def compute_l1_error_fvm(
    x_num: NDArray, numerical: NDArray, analytical: Callable, dim: int
) -> float:

    assert dim in [1, 2, 3]

    # For every cell, calculate the "volume" integral
    error = 0

    # NOTE: This calculation assumes a homogeneous grid
    if dim == 1:
        dx = x_num[1] - x_num[0]
        for i, x in enumerate(x_num):
            x1 = x + dx * 0.5
            x0 = x - dx * 0.5

            error += abs(dx * numerical[i] - simple_quad(analytical, x0, x1))
    elif dim == 2:
        dx = get_dx(x_num[0])
        dy = get_dx(x_num[1])

        assert dx is not None
        assert dy is not None
        for i, (x, y) in enumerate(zip(x_num[0], x_num[1])):
            x1 = x + dx * 0.5
            y1 = y + dy * 0.5

            x0 = x - dx * 0.5
            y0 = y - dy * 0.5

            error += abs(
                dx * dy * numerical[i]
                - simple_quad(
                    lambda y: simple_quad(
                        lambda x: analytical([x, y]), x0, x1
                    ),
                    y0,
                    y1,
                )
            )
    else:
        dx = get_dx(x_num[0])
        dy = get_dx(x_num[1])
        dz = get_dx(x_num[2])

        assert dx is not None
        assert dy is not None
        assert dz is not None
        for i, (x, y, z) in enumerate(zip(x_num[0], x_num[1], x_num[2])):
            x1 = x + dx * 0.5
            y1 = y + dy * 0.5
            z1 = z + dz * 0.5

            x0 = x - dx * 0.5
            y0 = y - dy * 0.5
            z0 = z - dz * 0.5

            error += abs(
                dx * dy * dz * numerical[i]
                - simple_quad(
                    lambda z: simple_quad(
                        lambda y: simple_quad(
                            lambda x: analytical([x, y, z]), x0, x1
                        ),
                        y0,
                        y1,
                    ),
                    z0,
                    z1,
                )
            )

    return error


def find_last_output(
    pattern: str = "output-?D-?-*.csv", dir: str | None = None
) -> tuple[str | None, bool]:
    """
    Return the last csv file and a boolean that indicates where the
    files were combined or not.
    """

    if dir is not None:
        pattern = os.path.join(dir, pattern)

    files = glob.glob(pattern)
    if not files:
        print("No matching output files found.")
        return None, False

    # Sort the files so they are in order
    files.sort()

    # If the last file is not a "*D-0-*.csv" file, this means
    # we need to combine them
    # We need to take every Nth file, where N = len(files) / (M + 1)
    # where M is the number in *D-M-*.csv
    mm = int(files[-1].split("-")[-2])
    if mm == 0:
        return files[-1], False
    else:
        nn = int(len(files) / (mm + 1))
        combine = files[nn - 1 :: nn]

        # Combine the files
        new_file = combine[-1].replace(f"D-{mm}-", f"D-{mm + 1}-")
        with open(new_file, "w") as fwrite:
            subprocess.run(["cat"] + combine, stdout=fwrite)
        return new_file, True


def parse_config(
    file_path: str,
) -> tuple[str, float, NDArray, NDArray, dict[str, Any]]:

    module_name = "custom_module"

    spec = imputil.spec_from_file_location(module_name, file_path)
    assert spec is not None
    assert spec.loader is not None

    module = imputil.module_from_spec(spec)
    spec.loader.exec_module(module)

    gamma = float(module.config.get("gamma", 1.4))
    x0 = np.array(module.config["coords"][0])
    x1 = np.array(module.config["coords"][1])
    problem = module.config["problem"]
    problem_dict = module.config.get("problem_parameters")

    return problem, gamma, x0, x1, problem_dict
