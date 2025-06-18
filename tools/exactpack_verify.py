import os
import sys
from collections.abc import Callable

import matplotlib.pyplot as plt
import numpy as np
from acoustic_solution import Acoustic
from numpy.typing import NDArray
from verify_lib import (compute_l1_error_fvm, parse_cli, parse_config,
                        wrapFunction)


class ProblemData(object):
    """
    Define the quantities for each problem
    """

    def __init__(self, name: str):
        """
        name - problem name
        """

        # Define names for Riemann problems
        self.riemann_problems = ["sod", "leblanc", "rankine-hugoniot"]

        # Default values
        self.name = name
        self.usecols = [0, 2, 3, 4, 9]
        self.extract = ["density", "pressure", "velocity"]
        self.labels = ["Density", "Pressure", "Velocity"]
        self.tolerances = [1e-2, 1e-2, 5e-1]
        self.function: Callable

        # Define specifics per problem and analytical solutions
        if self.name in self.riemann_problems:
            self.function = self.__riemann_analytic_solution
        elif self.name == "sedov":
            self.function = self.__sedov_analytic_solution
        elif self.name == "acoustic-wave":
            self.function = self.__acoustic_analytic_solution
        else:
            sys.exit(f"Unsupported problem type '{name}'")

        # For riemann problems, define left and right states
        if self.name in self.riemann_problems:
            if self.name == "sod":
                self.left_state = (1.0, 0.0, 1.0)
                self.right_state = (0.125, 0.0, 0.1)
            elif self.name == "leblanc":
                self.left_state = (1.0, 0.0, 0.1)
                self.right_state = (1e-3, 0.0, 1e-10)
            elif self.name == "rankine-hugoniot":
                self.left_state = (1.0, 0.0, 1.0)
                self.right_state = (0.25, 0.0, 0.1795)

    def __acoustic_analytic_solution(self, gamma: float, t: float, x0: float,
                                     x1: float, problem_dict: dict[str, str]
                                     ) -> tuple[Callable, Callable, Callable]:
        """
        Return the acoustic analytical solution
        """

        result = wrapFunction(
            Acoustic(gamma, x0, x1, problem_dict), t, self.extract)

        return result.density, result.pressure, result.velocity

    def __sedov_analytic_solution(self, x: NDArray, t: float, gamma: float
                                  ) -> tuple[Callable, Callable, Callable]:

        from exactpack.solvers.sedov import Sedov

        solver = Sedov(gamma=gamma, geometry=1, eblast=0.0673185)
        result = wrapFunction(solver, t, self.extract)

        return result.density, result.pressure, result.velocity

    def __riemann_analytic_solution(self, x: NDArray, t: float, gamma: float,
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

        sol = wrapFunction(solver, t, self.extract)
        return sol.density, sol.pressure, sol.velocity


class Problem(object):
    """
    Class to hold each type of problem parameters
    """

    def __init__(self, yaml_file: str) -> None:
        """
        yaml_file - the yaml_file with the problem configuration
        """

        name, gamma, x0, x1, problem_dict = parse_config(yaml_file)

        self.name = name
        self.gamma = gamma
        self.x0 = x0
        self.x1 = x1
        self.problem_dict = problem_dict

        # Initialize quantities
        self.data = ProblemData(self.name)

        # Variables to be defined later
        self.time: float
        self.x_arr: NDArray
        self.numerical: list[NDArray]
        self.references: list[NDArray]
        self.analytical: list[Callable]
        self.errors: list[float]
        self.status: list[str]

    def load_csv_file(self, csv_file: str) -> None:
        """
        Load the last csv file
        """
        self.csv_file = csv_file
        csv_out = np.loadtxt(csv_file, delimiter=",",
                             skiprows=1, usecols=self.data.usecols).T

        # Extract time and x_arr
        t_index = self.data.usecols.index(0)
        x_index = self.data.usecols.index(2)
        self.time = csv_out[t_index][0]
        self.x_arr = csv_out[x_index]

        # Extract the other quantities
        self.numerical = [x for i, x in enumerate(
            csv_out) if i not in [t_index, x_index]]

    def get_exact_solutions(self) -> None:
        """
        Get the exact solutions depending on the problem
        """

        assert self.x_arr is not None
        assert self.time is not None

        # Create input interator
        if self.name in self.data.riemann_problems:
            input = [self.x_arr, self.time, self.gamma,
                     self.data.left_state, self.data.right_state]
        elif self.name == "sedov":
            input = [self.x_arr, self.time, self.gamma]
        elif self.name == "acoustic-wave":
            input = [self.gamma, self.time,
                     self.x0, self.x1, self.problem_dict]
        elif self.name == "su-olson":
            input = [self.x_arr, self.time, self.data.extract]

        # Save the exact output functions
        self.analytical = self.data.function(*input)
        self.references = [f(self.x_arr) for f in self.analytical]

    def __plot_comparison(self, num_vals: NDArray, exact_vals: NDArray,
                          label: str, tag: str) -> None:
        """
        Plot the numerical output and the reference
        """

        plt.figure()
        plt.plot(self.x_arr, exact_vals, label="Analytic", linestyle="--")
        plt.plot(self.x_arr, num_vals, label="Simulation", marker='o',
                 linestyle='none', markersize=4)

        plt.xlabel("x")
        plt.ylabel(label)
        plt.title(f"{label} at t = {self.time:.4f}")

        plt.legend()
        plt.grid(True)
        plt.tight_layout()

        filename = f"{label.lower()}_comparison_{tag}.png"
        plt.savefig(filename)
        plt.close()

        print(f"Saved plot: {filename}")

    def make_plot(self) -> None:
        """
        Create plot with references
        """

        # Return the references
        tag = os.path.splitext(os.path.basename(self.csv_file))[0].replace(
            "output-", "").replace(".csv", "")

        for lab, num, ref in zip(self.data.labels, self.numerical,
                                 self.references):
            self.__plot_comparison(num, ref, lab, tag)

    def sedov_cutoff(self, cutoff: float = 0.2) -> None:
        """
        Remove the values up to cutoff for the sedov solution
        """

        mask = self.x_arr > cutoff
        self.x_arr = self.x_arr[mask]
        self.numerical = [x[mask] for x in self.numerical]
        self.references = [x[mask] for x in self.references]

    def compute_errors(self) -> None:
        """
        Compute all the errors
        """

        self.errors = [compute_l1_error_fvm(self.x_arr, num, exact)
                       for num, exact in zip(self.numerical, self.analytical)]

    def __color_text(self, text: str, status: str) -> str:
        if status == "PASS":
            return f"\033[92m{text}\033[0m"  # green
        elif status == "FAIL":
            return f"\033[91m{text}\033[0m"  # red
        else:
            return text

    def check_error_status(self) -> None:
        """
        Get status for every error
        """

        tols = self.data.tolerances
        failed = [False if x < y else True for x, y in zip(self.errors, tols)]
        statuses = ["FAIL" if x else "PASS" for x in failed]

        print(f"Problem type: {self.name}")
        print(f"Compared solution at t = {self.time:.4f}")
        print("L1 Errors and Status:")
        for label, err, stat in zip(self.data.labels, self.errors, statuses):
            print(f"  {label}  : {err:.6e} [{self.__color_text(stat, stat)}]")

        if any(failed):
            s = "Test FAILED: relative L1 error exceeds tolerance of "
            s += " ".join([f"{tol:.2e}" for tol in tols])
            sys.exit(s)
        else:
            print("Test PASSED")


def main() -> None:

    # Get the values
    yaml_file, _, csv_file, make_plot = parse_cli()
    assert csv_file is not None

    # Instantiate problem object
    problem = Problem(yaml_file)

    # Load csv file
    problem.load_csv_file(csv_file)

    # Get exact solutions
    problem.get_exact_solutions()

    if make_plot:
        problem.make_plot()

    # For Sedov: only compare for x > 0.2
    if problem == "sedov":
        problem.sedov_cutoff(cutoff=0.2)

    # Compute errors
    problem.compute_errors()

    # Error status
    problem.check_error_status()


if __name__ == "__main__":
    main()
