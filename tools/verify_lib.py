import argparse
import glob
import os
import subprocess
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

    parser = argparse.ArgumentParser()
    parser.add_argument("config_file", help="The config file")
    parser.add_argument("dim", type=int, help="Dimension of the problem")
    parser.add_argument(
        "-o", "--output", help="Optional output file or directory"
    )
    parser.add_argument(
        "--plot", action="store_true", help="If set, produces a plot"
    )
    args = parser.parse_args()

    # Assert the number of dimensions
    assert args.dim in [1, 2, 3]

    csv_file = None
    combined = False
    out_dir = None

    if args.output is not None:
        if args.output.endswith(".csv"):
            csv_file = args.output
        else:
            out_dir = args.output

    if get_file:
        # If no file is passed, select the latest available output
        if csv_file is None:
            pattern = f"output-{args.dim}D-?-*.csv"
            csv_file, combined = find_last_output(pattern=pattern, dir=out_dir)
            assert csv_file is not None

            print(f"Auto-selected input file: {csv_file}")

    return args.config_file, args.dim, out_dir, csv_file, combined, args.plot


def simple_quad(f: Callable, x0: NDArray, x1: NDArray, deg: int = 10) -> float:
    """
    Use a Gauss-Legendre quadrature of degree deg
    """

    points, weights = np.polynomial.legendre.leggauss(deg)

    def transform(x):
        return ((x1 - x0) * x + x1 + x0) * 0.5

    return np.sum(f(transform(points)) * weights) * 0.5 * (x1 - x0)


def simple_quad_2d(
    f: Callable, x0: float, x1: float, y0: float, y1: float, deg: int = 10
) -> float:
    """
    Use a tensor-product Gauss-Legendre quadrature of degree deg
    over the rectangle [x0, x1] x [y0, y1]
    """

    points, weights = np.polynomial.legendre.leggauss(deg)

    xs = ((x1 - x0) * points + x1 + x0) * 0.5
    ys = ((y1 - y0) * points + y1 + y0) * 0.5
    xg, yg = np.meshgrid(xs, ys, indexing="ij")
    wg = np.outer(weights, weights).ravel()

    return (
        np.sum(f([xg.ravel(), yg.ravel()]) * wg)
        * 0.25
        * (x1 - x0)
        * (y1 - y0)
    )


def simple_quad_3d(
    f: Callable,
    x0: float,
    x1: float,
    y0: float,
    y1: float,
    z0: float,
    z1: float,
    deg: int = 10,
) -> float:
    """
    Use a tensor-product Gauss-Legendre quadrature of degree deg
    over the box [x0, x1] x [y0, y1] x [z0, z1]
    """

    points, weights = np.polynomial.legendre.leggauss(deg)

    xs = ((x1 - x0) * points + x1 + x0) * 0.5
    ys = ((y1 - y0) * points + y1 + y0) * 0.5
    zs = ((z1 - z0) * points + z1 + z0) * 0.5
    xg, yg, zg = np.meshgrid(xs, ys, zs, indexing="ij")
    wg = np.outer(np.outer(weights, weights), weights).ravel()

    return (
        np.sum(f([xg.ravel(), yg.ravel(), zg.ravel()]) * wg)
        * 0.125
        * (x1 - x0)
        * (y1 - y0)
        * (z1 - z0)
    )


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
                - simple_quad_2d(analytical, x0, x1, y0, y1)
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
                - simple_quad_3d(analytical, x0, x1, y0, y1, z0, z1)
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
