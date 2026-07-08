import unittest

import numpy as np

from verify_lib import compute_l1_error_fvm


def cell_integral_power(bounds: list[tuple[float, float]]) -> float:
    result = 1.0
    for lower, upper in bounds:
        result *= (upper**2 - lower**2) * 0.5
    return result


class ComputeL1ErrorFVMTest(unittest.TestCase):
    def test_2d_uses_full_tensor_product_for_cell_integrals(self) -> None:
        centers = np.array([0.25, 0.75])
        xg, yg = np.meshgrid(centers, centers, indexing="ij")
        x_num = np.array([xg.ravel(), yg.ravel()])

        numerical = []
        for x, y in zip(x_num[0], x_num[1]):
            bounds = [
                (x - 0.25, x + 0.25),
                (y - 0.25, y + 0.25),
            ]
            numerical.append(cell_integral_power(bounds) / 0.25)

        error = compute_l1_error_fvm(
            x_num,
            np.array(numerical),
            lambda xy: xy[0] * xy[1],
            2,
        )

        self.assertAlmostEqual(error, 0.0)

    def test_3d_uses_full_tensor_product_for_cell_integrals(self) -> None:
        centers = np.array([0.25, 0.75])
        xg, yg, zg = np.meshgrid(centers, centers, centers, indexing="ij")
        x_num = np.array([xg.ravel(), yg.ravel(), zg.ravel()])

        numerical = []
        for x, y, z in zip(x_num[0], x_num[1], x_num[2]):
            bounds = [
                (x - 0.25, x + 0.25),
                (y - 0.25, y + 0.25),
                (z - 0.25, z + 0.25),
            ]
            numerical.append(cell_integral_power(bounds) / 0.125)

        error = compute_l1_error_fvm(
            x_num,
            np.array(numerical),
            lambda xyz: xyz[0] * xyz[1] * xyz[2],
            3,
        )

        self.assertAlmostEqual(error, 0.0)


if __name__ == "__main__":
    unittest.main()
