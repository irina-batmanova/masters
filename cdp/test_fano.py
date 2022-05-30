#!/usr/bin/env sage
from sage.all import *
from cdp import facet_is_at_height_one


def test_facet_height_one():
    p = Polyhedron(vertices=[[1, 0], [0, 1], [-1, 0], [0, -1]])
    assert facet_is_at_height_one(p.facets()[1].vertices()) is True
    assert facet_is_at_height_one(p.facets()[0].vertices()) is True


def test_facet_height_one_3d():
    p = Polyhedron(vertices=[[1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 0]])
    assert facet_is_at_height_one(p.facets()[0].vertices()) is True


def test_facet_not_height_one():
    p = Polyhedron(vertices=[[2, 0], [0, 2], [-2, 0], [0, -2]])
    assert facet_is_at_height_one(p.facets()[1].vertices()) is False


test_facet_height_one()
test_facet_height_one_3d()
test_facet_not_height_one()