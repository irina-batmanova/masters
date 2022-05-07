#!/usr/bin/env sage
from sage.all import *
from cdp import generate_cdp_from_polytope, CDP
from piecewise_affine_function import AffineFunction, PiecewiseAffineFunction


def test_generate_cdp_pyramid():
    poly = Polyhedron(vertices=[[1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0],
                                [0, 0, -1]])
    cdp = generate_cdp_from_polytope(poly)
    print(cdp)
    f_1 = AffineFunction(coefficients=[1, -1, -1], domain=Polyhedron(vertices=[[0, 0], [0, 1], [1, 0]]))
    f_2 = AffineFunction(coefficients=[1, 1, -1], domain=Polyhedron(vertices=[[0, 0], [0, 1], [-1, 0]]))
    f = PiecewiseAffineFunction([f_1, f_2])
    g = deepcopy(f)
    want_cdp = CDP(base=Polyhedron(vertices=[[-1, 0], [1, 0], [0, 1]]), psi_list=[f, g])
    assert cdp == want_cdp


def test_generate_cdp_square():
    poly = Polyhedron(vertices=[[-2, 0], [0, 2], [1, 2], [2, 1], [2, -2], [-2, -2]])
    cdp = generate_cdp_from_polytope(poly)
    print(cdp)
    # assert False


test_generate_cdp_pyramid()
test_generate_cdp_square()


