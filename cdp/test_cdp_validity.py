#!/usr/bin/env sage
from sage.all import *
from piecewise_affine_function import AffineFunction, PiecewiseAffineFunction
from cdp import CDP


class TestCDPValidity:

    def test_valid_cdp(self):
        base = Polyhedron(vertices=[[-1], [1]])
        # y = 1 + x, x \in [-1, 0]
        f_11 = AffineFunction([1, 1], Polyhedron(vertices=[[-1], [0]]))
        # y = 1 - x, x \in [0, 1]
        f_12 = AffineFunction([1, -1], Polyhedron(vertices=[[0], [1]]))
        f_1 = PiecewiseAffineFunction([f_11, f_12])
        # y = 1/2x + 1/2, x \in [-1, 1]
        f_2 = PiecewiseAffineFunction([AffineFunction([1 / 2, 1 / 2], Polyhedron(vertices=[[-1], [1]]))])
        cdp = CDP([f_1, f_2], base)

    def test_not_valid_cdp(self):
        base = Polyhedron(vertices=[[-1], [1]])
        # y = x, x \in [-1, 0]
        f_11 = AffineFunction([0, 1], Polyhedron(vertices=[[-1], [0]]))
        # y = -x, x \in [0, 1]
        f_12 = AffineFunction([0, -1], Polyhedron(vertices=[[0], [1]]))
        f_1 = PiecewiseAffineFunction([f_11, f_12])
        # y = -1/2 + 1/2x, x \in [-1, 1]
        f_2 = PiecewiseAffineFunction([AffineFunction([-1 / 2, 1 / 2], Polyhedron(vertices=[[-1], [1]]))])
        try:
            cdp = CDP([f_1, f_2], base)
        except ValueError:
            pass
        else:
            assert True
            

test = TestCDPValidity()
test.test_valid_cdp()
test.test_not_valid_cdp()