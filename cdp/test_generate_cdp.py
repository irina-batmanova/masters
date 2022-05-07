#!/usr/bin/env sage
from sage.all import *
from cdp import generate_cdp_from_polytope


def test_generate_cdp_pyramid():
    poly = Polyhedron(vertices=[[1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0],
                                [0, 0, -1]])
    cdp = generate_cdp_from_polytope(poly)
    print(cdp)
    base_verts = [v.vector() for v in cdp.base.vertices()]

    assert len(base_verts) == 3


test_generate_cdp_pyramid()

