import unittest

import numpy as np

from verify_lib import compute_l1_error_fvm, simple_quad_2d, simple_quad_3d


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


class SimpleQuadTensorProductTest(unittest.TestCase):
    def test_2d_asymmetric_box_and_integrand(self) -> None:
        got = simple_quad_2d(
            lambda c: np.sin(c[0] + 2.0 * c[1]), 0.0, 1.0, 0.0, 2.0
        )
        expected = (np.sin(1.0) - np.sin(5.0) + np.sin(4.0)) * 0.5
        self.assertAlmostEqual(got, expected)

    def test_3d_asymmetric_box_and_integrand(self) -> None:
        got = simple_quad_3d(
            lambda c: np.sin(c[0] + 2.0 * c[1]) * np.exp(0.5 * c[2]),
            0.0,
            1.0,
            0.0,
            2.0,
            0.0,
            3.0,
        )
        expected = (np.sin(1.0) - np.sin(5.0) + np.sin(4.0)) * (
            np.exp(1.5) - 1.0
        )
        self.assertAlmostEqual(got, expected)


if __name__ == "__main__":
    unittest.main()
