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

    def test_equality(self):
        base = Polyhedron(vertices=[[-1], [1]])
        f_11 = AffineFunction([1, 1], Polyhedron(vertices=[[-1], [0]]))
        f_12 = AffineFunction([1, -1], Polyhedron(vertices=[[0], [1]]))
        f_1 = PiecewiseAffineFunction([f_11, f_12])
        f_2 = PiecewiseAffineFunction([AffineFunction([1 / 2, 1 / 2], Polyhedron(vertices=[[-1], [1]]))])
        cdp1 = CDP([f_1, f_2], base)
        base2 = Polyhedron(vertices=[[-1], [1]])
        g_11 = AffineFunction([1, 1], Polyhedron(vertices=[[-1], [0]]))
        g_12 = AffineFunction([1, -1], Polyhedron(vertices=[[0], [1]]))
        g_1 = PiecewiseAffineFunction([g_11, g_12])
        g_2 = PiecewiseAffineFunction([AffineFunction([1 / 2, 1 / 2], Polyhedron(vertices=[[-1], [1]]))])
        cdp2 = CDP([g_2, g_1], base2)
        assert cdp1 == cdp2

    def test_inequality(self):
        base = Polyhedron(vertices=[[-1], [1]])
        f_11 = AffineFunction([1, 1], Polyhedron(vertices=[[-1], [0]]))
        f_12 = AffineFunction([1, -1], Polyhedron(vertices=[[0], [1]]))
        f_1 = PiecewiseAffineFunction([f_11, f_12])
        f_2 = PiecewiseAffineFunction([AffineFunction([1 / 2, 1 / 2], Polyhedron(vertices=[[-1], [1]]))])
        cdp1 = CDP([f_1, f_2], base)
        base2 = Polyhedron(vertices=[[-1], [2]])
        g_11 = AffineFunction([1, 1], Polyhedron(vertices=[[-1], [0]]))
        g_12 = AffineFunction([1, -1], Polyhedron(vertices=[[0], [2]]))
        g_1 = PiecewiseAffineFunction([g_11, g_12])
        g_2 = PiecewiseAffineFunction([AffineFunction([1 / 2, 1 / 2], Polyhedron(vertices=[[-1], [2]]))])
        cdp2 = CDP([g_2, g_1], base2)
        assert not cdp1 == cdp2
            

test = TestCDPValidity()
test.test_valid_cdp()
test.test_not_valid_cdp()
test.test_equality()
test.test_inequality()
