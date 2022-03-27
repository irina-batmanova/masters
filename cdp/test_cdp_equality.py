#!/usr/bin/env sage
from sage.all import *
from piecewise_affine_function import AffineFunction, PiecewiseAffineFunction
from cdp import CDP


class TestCDPEquality:

    def test_equal(self):
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


test = TestCDPEquality()
test.test_equal()
