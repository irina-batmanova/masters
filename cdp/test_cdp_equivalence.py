#!/usr/bin/env sage
from sage.all import *
from piecewise_affine_function import AffineFunction, PiecewiseAffineFunction
from cdp import CDP


class TestCDPEquality:

    def test_1d_equal(self):
        base1 = Polyhedron(vertices=[[-1], [1]])
        # y = 1 + x, x \in [-1, 0]
        f_11 = AffineFunction([1, 1], Polyhedron(vertices=[[-1], [0]]))
        # y = 1 - x, x \in [0, 1]
        f_12 = AffineFunction([1, -1], Polyhedron(vertices=[[0], [1]]))
        f_1 = PiecewiseAffineFunction([f_11, f_12])
        # y = 1/2 + 1/2x, x \in [-1, 1]
        f_2 = PiecewiseAffineFunction([AffineFunction([1 / 2, 1 / 2], Polyhedron(vertices=[[-1], [1]]))])
        cdp1 = CDP([f_1, f_2], base1)
        base2 = Polyhedron(vertices=[[-1], [1]])
        # y = 2, x \in [-1, 0]
        y_11 = AffineFunction([2, 0], Polyhedron(vertices=[[-1], [0]]))
        # y = 2 - 2x, x \in [0, 1]
        y_12 = AffineFunction([2, -2], Polyhedron(vertices=[[0], [1]]))
        y_1 = PiecewiseAffineFunction([y_11, y_12])
        # y = -1/2 + 1/2x, x \in [-1, 1]
        y_2 = PiecewiseAffineFunction([AffineFunction([-1 / 2, 1 / 2], Polyhedron(vertices=[[-1], [1]]))])
        cdp2 = CDP([y_1, y_2], base2)
        assert cdp1.equal(cdp2)
        assert cdp2.equal(cdp1)
        assert cdp1.equal(cdp1)

    def test_1d_not_equal(self):
        base1 = Polyhedron(vertices=[[-1], [1]])
        # y = 1 + x, x \in [-1, 0]
        f_11 = AffineFunction([1, 1], Polyhedron(vertices=[[-1], [0]]))
        # y = 1 - x, x \in [0, 1]
        f_12 = AffineFunction([1, -1], Polyhedron(vertices=[[0], [1]]))
        f_1 = PiecewiseAffineFunction([f_11, f_12])
        # y = 1/2 + 1/2x, x \in [-1, 1]
        f_2 = PiecewiseAffineFunction([AffineFunction([1 / 2, 1 / 2], Polyhedron(vertices=[[-1], [1]]))])
        cdp1 = CDP([f_1, f_2], base1)
        base2 = Polyhedron(vertices=[[-1], [1]])
        # y = 2, x \in [-1, 0]
        y_11 = AffineFunction([2, 0], Polyhedron(vertices=[[-1], [0]]))
        # y = 2 - x, x \in [0, 1]
        y_12 = AffineFunction([2, -1], Polyhedron(vertices=[[0], [1]]))
        y_1 = PiecewiseAffineFunction([y_11, y_12])
        # y = -1/2 + 1/2x, x \in [-1, 1]
        y_2 = PiecewiseAffineFunction([AffineFunction([-1 / 2, 1 / 2], Polyhedron(vertices=[[-1], [1]]))])
        cdp2 = CDP([y_1, y_2], base2)
        assert not cdp1.equal(cdp2)
        assert not cdp2.equal(cdp1)

    def test_2d_equal(self):
        base1 = Polyhedron(vertices=[[1, 0], [0, -1], [-1, 0], [0, 1]])
        f_11 = AffineFunction([1, -1, -1], Polyhedron(vertices=[[0, 0], [1, 0], [0, 1]]))
        f_12 = AffineFunction([1, -1, 1], Polyhedron(vertices=[[0, 0], [1, 0], [0, -1]]))
        f_13 = AffineFunction([1, 1, -1], Polyhedron(vertices=[[0, 0], [-1, 0], [0, 1]]))
        f_14 = AffineFunction([1, 1, 1], Polyhedron(vertices=[[0, 0], [-1, 0], [0, -1]]))
        f_1 = PiecewiseAffineFunction([f_11, f_12, f_13, f_14])
        f_21 = AffineFunction([1, 0, 0], Polyhedron(vertices=[[0, 1], [1, 0], [0, -1]]))
        f_22 = AffineFunction([1, 1, 1], Polyhedron(vertices=[[0, 0], [0, -1], [-1, 0]]))
        f_23 = AffineFunction([1, 1, -1], Polyhedron(vertices=[[0, 0], [-1, 0], [0, 1]]))
        f_2 = PiecewiseAffineFunction([f_21, f_22, f_23])
        cdp1 = CDP([f_1, f_2], base1)
        cdp2 = deepcopy(cdp1)
        cdp2.shear([1, -1], [1, 1])
        phi = linear_transformation(matrix(ZZ, [[-1, 0], [0, -1]]))
        cdp2.transform_base(phi)
        assert cdp1.equal(cdp2)


test = TestCDPEquality()
test.test_1d_equal()
test.test_1d_not_equal()
test.test_2d_equal()

# TODO: why this transfomation breaks everything?
# phi = linear_transformation(matrix(ZZ, [[-1, 1], [0, -1]]))
