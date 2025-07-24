import glob
import os
import sys
from collections.abc import Callable
from typing import Any

import numpy as np
import yaml
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


def parse_cli(get_file: bool = True
              ) -> tuple[str, int, str | None, str | None, bool]:
    """
    Parse command line input
    """

    if len(sys.argv) < 2:
        s = f"Usage: python {sys.argv[0]} <config_file> <dim>"
        s += " [output_dir | output.csv]"
        s += " [--plot]"
        sys.exit(s)

    # Read in the yaml file
    yaml_file = sys.argv[1]
    problem, gamma, x0, x1, problem_dict = parse_config(yaml_file)

    # Get the number of dimensions
    dim = int(sys.argv[2])
    assert dim in [1, 2, 3]

    csv_file = None
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
            pattern = f"output-{dim}D-0-*.csv"
            csv_file = find_last_output(pattern=pattern, dir=out_dir)
            assert csv_file is not None

            print(f"Auto-selected input file: {csv_file}")

    return yaml_file, dim, out_dir, csv_file, make_plot


def simple_quad(f: Callable, x0: NDArray, x1: NDArray, deg: int = 10) -> float:
    """
    Use a Gauss-Legendre quadrature of degree deg
    """

    points, weights = np.polynomial.legendre.leggauss(deg)
    def transform(x): return ((x1 - x0) * x + x1 + x0) * 0.5

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


def compute_l1_error_fvm(x_num: NDArray, numerical: NDArray,
                         analytical: Callable, dim: int) -> float:

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
        dy = x_num[1][1] - x_num[1][0]

        assert dx is not None
        for i, (x, y) in enumerate(zip(x_num[0], x_num[1])):
            x1 = x + dx * 0.5
            x0 = x - dx * 0.5
            y1 = y + dy * 0.5
            y0 = y - dy * 0.5

            error += abs(dx * dy * numerical[i] - simple_quad(
                lambda y: simple_quad(
                    lambda x: analytical([x, y]),
                    x0, x1),
                y0, y1))
    else:
        dx = get_dx(x_num[0])
        dy = get_dx(x_num[1])
        dz = x_num[2][1] - x_num[2][0]

        assert dx is not None
        assert dy is not None
        for i, (x, y, z) in enumerate(zip(x_num[0], x_num[1], x_num[2])):
            x1 = x + dx * 0.5
            x0 = x - dx * 0.5
            y1 = y + dy * 0.5
            y0 = y - dy * 0.5
            z1 = z + dz * 0.5
            z0 = z - dz * 0.5

            error += abs(dx * dy * dz * numerical[i] - simple_quad(
                lambda z: simple_quad(
                    lambda y: simple_quad(
                        lambda x: analytical([x, y, z]), x0, x1),
                    y0, y1),
                z0, z1))

    return error


def find_last_output(pattern: str = "output-?D-0-*.csv",
                     dir: str | None = None) -> str | None:

    if dir is not None:
        pattern = os.path.join(dir, pattern)

    files = glob.glob(pattern)
    if not files:
        print("No matching output files found.")
        return None

    # Sort the files so they are in order
    files.sort()

    # Return first and last
    return files[-1]


def parse_config(yaml_file: str
                 ) -> tuple[str, float, NDArray, NDArray, dict[str, Any]]:
    with open(yaml_file, 'r') as f:
        config = yaml.safe_load(f)
    gamma = float(config.get('gamma', 1.4))
    x0 = np.array(config['coords'][0])
    x1 = np.array(config['coords'][1])
    problem = config['problem']
    problem_dict = config.get("problem_parameters")

    return problem, gamma, x0, x1, problem_dict
