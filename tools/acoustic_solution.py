from collections.abc import Callable
from typing import Any

import numpy as np
from numpy.typing import NDArray


class Acoustic(object):

    def __init__(self, gamma: float, x0: NDArray, x1: NDArray,
                 problem_dict: dict[str, Any], dim: int) -> None:

        self.gamma = gamma
        self.x0 = x0[:dim]
        self.x1 = x1[:dim]
        self.span = np.abs(self.x1 - self.x0)
        self.problem_dict = problem_dict
        self.dim = dim

        self.r0 = float(self.problem_dict["r0"])
        self.p0 = float(self.problem_dict["p0"])
        self.cs = np.sqrt(gamma * self.p0 / self.r0)

        self.amplitude = float(self.problem_dict["amplitude"])
        self.scale = 2 * np.pi * \
            np.array(self.problem_dict["scale"])[:dim] / self.span

        # Scale norm
        self.scale_norm = np.sqrt(self.scale.dot(self.scale))

        # Reshape arrays so that they can operate correctly
        # with self.x0 and self.scale
        if dim != 1:
            self.x0 = self.x0.reshape((dim, 1))
            self.scale = self.scale.reshape((dim, 1))
            self.span = self.span.reshape((dim, 1))

        # Directions
        self.directions = self.scale / self.scale_norm

    def __perturbation(self) -> Callable:
        """
        Perturbation function shape
        """

        if self.dim == 1:
            return lambda coords: np.sin(self.scale * coords) * self.amplitude
        else:
            return lambda coords: np.sin(np.sum(self.scale * coords,
                                                axis=0)) * self.amplitude

    def __call__(self, coordinates: list[NDArray] | list[float], t: float
                 ) -> dict[str, NDArray]:
        '''
        Evaluate the perturbation shape
        '''

        # Take the initial solution and transport it by cs * t, assuming
        # periodic boundary conditions

        perturbation = self.__perturbation()
        density = self.r0
        pressure = self.p0

        # Positive movement
        shifted_coords = (coordinates - t * self.cs * self.directions -
                          self.x0) % self.span + self.x0

        density += perturbation(shifted_coords) * 0.5
        pressure += perturbation(shifted_coords) * 0.5 * self.cs ** 2
        velocity = perturbation(shifted_coords) * 0.5 * self.cs

        # Entropic mode (static density)
        density += perturbation(coordinates)

        # Negative movement
        shifted_coords = (coordinates + t * self.cs * self.directions -
                          self.x0) % self.span + self.x0

        density -= perturbation(shifted_coords) * 0.5
        pressure -= perturbation(shifted_coords) * 0.5 * self.cs ** 2
        velocity += perturbation(shifted_coords) * 0.5 * self.cs

        solution = {
            "density": density,
            "pressure": pressure,
            "velocity": velocity
        }

        if self.dim > 1:
            solution["velocity_x"] = velocity * self.directions[0]
            solution["velocity_y"] = velocity * self.directions[1]
        if self.dim == 3:
            solution["velocity_z"] = velocity * self.directions[2]

        return solution
