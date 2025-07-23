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

        self.r0 = float(self.problem_dict["r0"])
        self.p0 = float(self.problem_dict["p0"])
        self.cs = np.sqrt(gamma * self.p0 / self.r0)

        self.amplitude = float(self.problem_dict["amplitude"])
        self.scale = 2 * np.pi * \
            np.array(self.problem_dict["scale"])[:dim] / self.span


class Acoustic_1D(Acoustic):
    """
    Class for acoustic wave solution for 1D
    """

    def __init__(self, *args):
        super().__init__(*args, dim=1)

    def __perturbation(self):
        """
        Perturbation function shape
        """

        return lambda x: np.sin(self.scale * x) * self.amplitude

    def __call__(self, x: NDArray | float, t: float) -> dict[str, NDArray]:
        """
        Return dictionary with what the density, pressure and velocity are at
        time t

        Assume periodic boundary conditions
        """

        # Take the initial solution and transport it by cs * t, assuming
        # periodic boundary conditions

        perturbation = self.__perturbation()
        density = self.r0
        pressure = self.p0

        # Positive movement
        x_init = (x - self.cs * t - self.x0) % self.span + self.x0

        density += perturbation(x_init) * 0.5
        pressure += perturbation(x_init) * 0.5 * self.cs ** 2
        velocity = perturbation(x_init) * 0.5 * self.cs

        # Entropic mode (static density)
        density += perturbation(x)

        # Negative movement
        x_init = (x + self.cs * t - self.x0) % self.span + self.x0

        density -= perturbation(x_init) * 0.5
        pressure -= perturbation(x_init) * 0.5 * self.cs ** 2
        velocity += perturbation(x_init) * 0.5 * self.cs

        solution = {
            "density": density,
            "pressure": pressure,
            "velocity": velocity
        }

        return solution


class Acoustic_2D(Acoustic):
    """
    Class for acoustic wave solution for 2D
    """

    def __init__(self, *args):
        super().__init__(*args, dim=2)

    def __perturbation(self):
        """
        Perturbation function shape
        """

        return lambda x, y: np.sin(self.scale[0] * x +
                                   self.scale[1] * y) * self.amplitude

    def __call__(self, x: NDArray | float, y: NDArray | float, t: float
                 ) -> dict[str, NDArray]:
        """
        Return dictionary with what the density, pressure and velocity are at
        time t

        Assume periodic boundary conditions
        """

        # Take the initial solution and transport it by cs * t, assuming
        # periodic boundary conditions

        perturbation = self.__perturbation()
        density = self.r0
        pressure = self.p0

        x_init = x * 0

        # Positive movement
        x_init = (x - self.cs * t - self.x0[0]) % self.span[0] + self.x0[0]
        y_init = (y - self.cs * t - self.x0[1]) % self.span[1] + self.x0[1]

        density += perturbation(x_init, y_init) * 0.5
        pressure += perturbation(x_init, y_init) * 0.5 * self.cs ** 2

        velocity = perturbation(x_init, y_init) * 0.5 * self.cs

        # Entropic mode (static density)
        density += perturbation(x, y)

        # Negative movement
        x_init = (x + self.cs * t - self.x0[0]) % self.span[0] + self.x0[0]
        y_init = (y + self.cs * t - self.x0[1]) % self.span[1] + self.x0[1]

        density -= perturbation(x_init, y_init) * 0.5
        pressure -= perturbation(x_init, y_init) * 0.5 * self.cs ** 2
        velocity += perturbation(x_init, y_init) * 0.5 * self.cs

        # Scale norm
        norm_scale = self.scale.dot(self.scale)

        solution = {
            "density": density,
            "pressure": pressure,
            "velocity_x": velocity * self.scale[0] / norm_scale,
            "velocity_y": velocity * self.scale[1] / norm_scale
        }

        return solution


class Acoustic_3D(Acoustic):
    """
    Class for acoustic wave solution for 3D
    """

    def __init__(self, *args):
        super().__init__(*args, dim=3)

    def __perturbation(self):
        """
        Perturbation function shape
        """

        return lambda x, y, z: np.sin(self.scale[0] * x +
                                      self.scale[1] * y +
                                      self.scale[2] * z) * self.amplitude

    def __call__(self, x: NDArray | float, y: NDArray | float,
                 z: NDArray | float, t: float) -> dict[str, NDArray]:
        """
        Return dictionary with what the density, pressure and velocity are at
        time t

        Assume periodic boundary conditions
        """

        # Take the initial solution and transport it by cs * t, assuming
        # periodic boundary conditions

        perturbation = self.__perturbation()
        density = self.r0
        pressure = self.p0

        x_init = x * 0

        # Positive movement
        x_init = (x - self.cs * t - self.x0[0]) % self.span[0] + self.x0[0]
        y_init = (y - self.cs * t - self.x0[1]) % self.span[1] + self.x0[1]
        z_init = (z - self.cs * t - self.x0[2]) % self.span[2] + self.x0[2]

        density += perturbation(x_init, y_init, z_init) * 0.5
        pressure += perturbation(x_init, y_init, z_init) * 0.5 * self.cs ** 2

        velocity = perturbation(x_init, y_init, z_init) * 0.5 * self.cs

        # Entropic mode (static density)
        density += perturbation(x, y, z)

        # Negative movement
        x_init = (x + self.cs * t - self.x0[0]) % self.span[0] + self.x0[0]
        y_init = (y + self.cs * t - self.x0[1]) % self.span[1] + self.x0[1]
        z_init = (z + self.cs * t - self.x0[2]) % self.span[2] + self.x0[2]

        density -= perturbation(x_init, y_init, z_init) * 0.5
        pressure -= perturbation(x_init, y_init, z_init) * 0.5 * self.cs ** 2
        velocity += perturbation(x_init, y_init, z_init) * 0.5 * self.cs

        # Scale norm
        norm_scale = self.scale.dot(self.scale)

        solution = {
            "density": density,
            "pressure": pressure,
            "velocity_x": velocity * self.scale[0] / norm_scale,
            "velocity_y": velocity * self.scale[1] / norm_scale,
            "velocity_z": velocity * self.scale[2] / norm_scale
        }

        return solution
