from collections import namedtuple

import numpy as np
import scipy as scp
from numpy.typing import NDArray

Solution = namedtuple("Solution", ["density", "pressure", "velocity"])


class Acoustic(object):
    """
    Class for acoustic wave solution
    """

    def __init__(self, grid: NDArray, gamma: float, r0: float, p0: float,
                 init_r: NDArray, init_p: NDArray, init_u: NDArray) -> None:
        self.gamma = gamma
        self.r0 = r0
        self.p0 = p0
        self.cs = np.sqrt(gamma * p0 / r0)
        self.init_r = scp.interpolate.interp1d(
            grid, init_r, fill_value="extrapolate")
        self.init_p = scp.interpolate.interp1d(
            grid, init_p, fill_value="extrapolate")
        self.init_u = scp.interpolate.interp1d(
            grid, init_u, fill_value="extrapolate")

    def __call__(self, x: NDArray | float, t: float) -> Solution:
        """
        Return named tuple with what the density, pressure and velocity are at
        time t

        Assume periodic boundary conditions
        """

        # Take the initial solution and transport it by cs * t, assuming
        # periodic boundary conditions

        density = self.r0
        pressure = self.p0

        # Positive movement
        x_init = x - self.cs * t
        while np.min(x_init) < 0:
            x_init = np.where(x_init < 0, x_init + 1, x_init)

        density += (self.init_r(x_init) - self.r0) * 0.5
        pressure += (self.init_r(x_init) - self.r0) * 0.5 * self.gamma
        velocity = self.init_u(x_init) * 0.5

        # Entropic mode (static density)
        density += self.init_r(x) - self.r0

        # Negative movement
        x_init = x + self.cs * t
        while np.max(x_init) > 1:
            x_init = np.where(x_init > 1, x_init - 1, x_init)

        density += (self.r0 - self.init_r(x_init)) * 0.5
        pressure += (self.r0 - self.init_r(x_init)) * 0.5 * self.gamma
        velocity += self.init_u(x_init) * 0.5

        return Solution(density, pressure, velocity)
