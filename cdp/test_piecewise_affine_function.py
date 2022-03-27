#!/usr/bin/env sage
from sage.all import *

from piecewise_affine_function import AffineFunction, PiecewiseAffineFunction


class TestPiecewiseAffineFunctionStr:

    def test_1_piece(self):
        # y = 1 + x, x \in [-1, 1]
        f_1 = AffineFunction([1, 1], Polyhedron(vertices=[[-1], [1]]))
        f = PiecewiseAffineFunction([f_1])
        assert f.__str__() == 'Piecewise affine function:\nAffine function 1 + x_1 with domain [(-1), (1)]'

    def test_2_pieces(self):
        # y = 1 + x, x \in [-1, 0]
        f_1 = AffineFunction([1, 1], Polyhedron(vertices=[[-1], [0]]))
        # y = 1 - x, x \in [0, 1]
        f_2 = AffineFunction([1, -1], Polyhedron(vertices=[[0], [1]]))
        f = PiecewiseAffineFunction([f_1, f_2])
        assert f.__str__() == 'Piecewise affine function:\nAffine function 1 + x_1 with domain [(-1), (0)]\n' \
                              'Affine function 1 - x_1 with domain [(0), (1)]'


test = TestPiecewiseAffineFunctionStr()
test.test_1_piece()
test.test_2_pieces()