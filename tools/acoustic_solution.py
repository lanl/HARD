from collections import namedtuple

import numpy as np
from numpy.typing import NDArray

Solution = namedtuple("Solution", ["density", "pressure", "velocity"])


class Acoustic(object):
    """
    Class for acoustic wave solution
    """

    def __init__(self, gamma: float, x0: float, x1: float,
                 problem_dict: dict[str, str]) -> None:

        self.gamma = gamma
        self.x0 = x0
        self.x1 = x1
        self.problem_dict = problem_dict

        self.r0 = float(self.problem_dict["r0"])
        self.p0 = float(self.problem_dict["p0"])
        self.cs = np.sqrt(gamma * self.p0 / self.r0)

        self.amplitude = float(self.problem_dict["amplitude"])
        self.scale = 2 * np.pi * float(self.problem_dict["scale"])

    def __perturbation(self):
        """
        Perturbation function shape
        """

        return lambda x: np.sin(self.scale * x) * self.amplitude

    def __call__(self, x: NDArray | float, t: float) -> Solution:
        """
        Return named tuple with what the density, pressure and velocity are at
        time t

        Assume periodic boundary conditions
        """

        # Take the initial solution and transport it by cs * t, assuming
        # periodic boundary conditions

        perturbation = self.__perturbation()
        span = abs(self.x1 - self.x0)

        density = self.r0
        pressure = self.p0

        # Positive movement
        x_init = (x - self.cs * t - self.x0) % span + self.x0

        density += perturbation(x_init) * 0.5
        pressure += perturbation(x_init) * 0.5 * self.cs ** 2
        velocity = perturbation(x_init) * 0.5 * self.cs

        # Entropic mode (static density)
        density += perturbation(x)

        # Negative movement
        x_init = (x + self.cs * t - self.x0) % span + self.x0

        density -= perturbation(x_init) * 0.5
        pressure -= perturbation(x_init) * 0.5 * self.cs ** 2
        velocity += perturbation(x_init) * 0.5 * self.cs

        return Solution(density, pressure, velocity)
