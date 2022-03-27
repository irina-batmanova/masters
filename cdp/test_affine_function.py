#!/usr/bin/env sage
from sage.all import *

from piecewise_affine_function import AffineFunction

# TODO: How to use pytest with sage? Calling all test funcs manually is obviously bad
# (and it's impossible to use fixtures)


one_dim_polyhedron = Polyhedron([[-1], [1]])
one_dim_domain = '[(-1), (1)]'
two_dim_polyhedron = Polyhedron([[0, 0], [2, 2], [0, 2], [2, 0]])
two_dim_domain = '[(0, 0), (0, 2), (2, 0), (2, 2)]'


class TestAffineFunctionStr:

    def test_const_str(self):
        f = AffineFunction([-1, 0], one_dim_polyhedron)
        assert f.__str__() == f'Affine function -1 with domain {one_dim_domain}'

    def test_single_variable_str(self):
        f = AffineFunction([-1, 1], one_dim_polyhedron)
        assert f.__str__() == f'Affine function -1 + x_1 with domain {one_dim_domain}'

    def test_2variables_str(self):
        f = AffineFunction([-2, -2, 1], two_dim_polyhedron)
        assert f.__str__() == f'Affine function -2 - 2x_1 + x_2 with domain {two_dim_domain}'

    def test_skip_variable(self):
        f = AffineFunction([-2, 0, 1], two_dim_polyhedron)
        assert f.__str__() == f'Affine function -2 + x_2 with domain {two_dim_domain}'

    def test_zero_func(self):
        f = AffineFunction([0, 0], one_dim_polyhedron)
        assert f.__str__() == f'Affine function 0 with domain {one_dim_domain}'


class TestAffineFunctionValue:

    def test_const_value(self):
        f = AffineFunction([-1, 0], one_dim_polyhedron)
        # f = -1, f(-1) = -1
        assert f.value([-1]) == -1

    def test_single_variable(self):
        f = AffineFunction([1, -1], one_dim_polyhedron)
        # f = 1 - x
        assert f.value([-1]) == 2

    def test_skip_variable(self):
        f = AffineFunction([-2, 0, 1], two_dim_polyhedron)
        # f = -2 + x_2
        assert f.value([1, 1]) == -1


test = TestAffineFunctionStr()
test.test_const_str()
test.test_single_variable_str()
test.test_2variables_str()
test.test_skip_variable()
test.test_zero_func()

test_val = TestAffineFunctionValue()
test_val.test_const_value()
test_val.test_single_variable()
test_val.test_skip_variable()