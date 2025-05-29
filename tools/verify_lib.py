import glob
import os
import sys
from collections.abc import Callable

import numpy as np
import yaml
from acoustic_solution import Acoustic
from numpy.typing import NDArray


class wrapFunction(object):

    # Declaring the attributes to provide mypy
    # with enough information
    density: Callable
    pressure: Callable
    velocity: Callable

    def __init__(self, solver: Callable, t: float,
                 extract: list[str] = ["density", "pressure", "velocity"]
                 ) -> None:
        """
        Create the wrapper for each attribute
        """

        for name in extract:
            self.__doWrapping(solver, t, name)

    def __doWrapping(self, f: Callable, t: float, name: str) -> None:
        """
        Transform the ExactPack array tuple into a Callable tuple
        """

        if isinstance(f, Acoustic):
            def wrap(x):
                return f(x, t)._asdict()[name]
        else:
            def wrap(x):
                return f(x, t)[name]

        self.__dict__[name] = wrap


def parse_cli(get_file: bool = True) -> tuple[
        str,
        str | None,
        str | None,
        bool
]:
    """
    Parse command line input
    """

    if len(sys.argv) < 2:
        s = f"Usage: python {sys.argv[0]} <config_file>"
        s += " [output_dir | output.raw]"
        s += " [--plot]"
        sys.exit(s)

    # Read in the yaml file
    yaml_file = sys.argv[1]
    problem, gamma, x0, x1, problem_dict = parse_config(yaml_file)

    raw_file = None
    out_dir = None
    make_plot = False

    for arg in sys.argv[2:]:
        if arg.endswith(".raw"):
            raw_file = arg
        elif arg == "--plot":
            make_plot = True
        else:
            out_dir = arg

    if get_file:

        # If no file is passed, select the latest available output
        if raw_file is None:
            raw_file = find_last_output(dir=out_dir)
            assert raw_file is not None

            print(f"Auto-selected input file: {raw_file}")

    return yaml_file, out_dir, raw_file, make_plot


def simple_quad(f: Callable, x0: float, x1: float, deg=10) -> float:
    """
    Use a Gauss-Legendre quadrature of degree deg
    """

    points, weights = np.polynomial.legendre.leggauss(deg)
    def transform(x): return ((x1 - x0) * x + x1 + x0) * 0.5

    return np.sum(f(transform(points)) * weights) * 0.5 * (x1 - x0)


def compute_l1_error_fvm(x_num: NDArray, numerical: NDArray,
                         analytical: Callable) -> float:

    # For every cell, calculate the "volume" integral
    dx = x_num[1] - x_num[0]
    error = 0

    for i, x in enumerate(x_num):
        x1 = x + dx * 0.5
        x0 = x - dx * 0.5

        error += abs(dx * numerical[i] - simple_quad(analytical, x0, x1))

    return error


def find_last_output(pattern: str = "output-*-1D-0.raw",
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
                 ) -> tuple[str, float, float, float, dict[str, str]]:
    with open(yaml_file, 'r') as f:
        config = yaml.safe_load(f)
    gamma = float(config.get('gamma', 1.4))
    x0 = float(config['coords'][0][0])
    x1 = float(config['coords'][1][0])
    problem = config['problem']
    problem_dict = config.get("problem_parameters")

    return problem, gamma, x0, x1, problem_dict


def acoustic_analytic_solution(gamma: float, t: float, x0: float, x1: float,
                               problem_dict: dict[str, str]) -> tuple[Callable,
                                                                      Callable,
                                                                      Callable
                                                                      ]:
    """
    Return the acoustic analytical solution
    """

    names = ["density", "pressure", "velocity"]

    result = wrapFunction(
        Acoustic(gamma, x0, x1, problem_dict), t, extract=names)

    return result.density, result.pressure, result.velocity
