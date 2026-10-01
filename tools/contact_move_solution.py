from typing import Any

import numpy as np
from numpy.typing import NDArray


class ContactMove(object):
    def __init__(
        self,
        gamma0: float,
        gamma1: float,
        x0: NDArray,
        x1: NDArray,
        problem_dict: dict[str, Any],
        dim: int,
    ) -> None:
        self.gamma0 = gamma0
        self.gamma1 = gamma1
        self.x0 = x0[:dim]
        self.x1 = x1[:dim]
        self.problem_dict = problem_dict
        self.dim = dim

        # Exact CONTACT_move test parameters from shock_tube_mm.hh
        self.rL = 1.0
        self.uL = 1.0
        self.pL = 1.0

        self.rR = 0.125
        self.uR = 1.0
        self.pR = 1.0

        self.x_interface0 = 0.5
        self.material_left = 0
        self.material_right = 1

    def _x_coordinates(
        self, coordinates: list[NDArray] | list[float] | NDArray
    ) -> NDArray:
        """
        Extract the x-coordinate array from 1D/2D/3D coordinate input.
        """

        coords = np.asarray(coordinates)
        if self.dim == 1:
            return coords
        return np.asarray(coordinates[0])

    def _interface_position(self, t: float) -> float:
        """
        Exact contact position advected at constant speed.
        """

        return self.x_interface0 + self.uL * t

    def __call__(
        self, coordinates: list[NDArray] | list[float] | NDArray, t: float
    ) -> dict[str, NDArray]:
        """
        Evaluate the exact moving-contact solution.

        For CONTACT_move, the contact is simply advected in x with speed u = 1.0.
        At t = 0.2 the discontinuity is therefore located at x = 0.7.
        """

        x_coords = self._x_coordinates(coordinates)
        x_interface = self._interface_position(t)
        left_mask = x_coords < x_interface

        density = np.where(left_mask, self.rL, self.rR)
        pressure = np.where(left_mask, self.pL, self.pR)
        velocity = np.where(left_mask, self.uL, self.uR)
        material = np.where(left_mask, self.material_left, self.material_right)

        solution = {
            "density": density,
            "pressure": pressure,
            "velocity": velocity,
            "material": material,
        }

        if self.dim > 1:
            zeros = np.zeros_like(velocity)
            solution["velocity_x"] = velocity
            solution["velocity_y"] = zeros
        if self.dim == 3:
            solution["velocity_z"] = np.zeros_like(velocity)

        return solution
