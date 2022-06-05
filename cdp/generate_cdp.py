#!/usr/bin/env sage
from sage.all import *
from sympy.matrices import Matrix
from sympy import Rational, lcm
from piecewise_affine_function import PiecewiseAffineFunction, AffineFunction
from cdp import CDP


def generate_cdp_from_polytope(poly: Polyhedron):
    # Walk over all facets, if facet's normal projection to x_n is positive -
    # this facet if a part of psi_1 graph, otherwise - part of -psi_2 grapg

    def plane_from_points(points):
        # (p_11, ..., p_1n), ..., (p_n1, ..., p_nn)
        # p_11, p_21, ..., p_n1
        #         ...             = P
        # p_1n, p_2n, ..., p_nn
        # AP = b, b = (1, ..., 1) => A = bP^-1
        # 1 - a_1*x_1 - ... - a_n*x_n = 0, return [1, -a1, ..., -a_n]
        points = points[:len(points[0])]
        b = Matrix([[1 for i in range(len(points[0]))]])
        P = Matrix(points).T
        A = b * P.inv()
        return [Rational(1)] + [-a for a in list(A.row(0))]

    def integer_coefs_from_rational(coefs):
        qs = [c.q for c in coefs]
        l = lcm(qs)
        coefs = [c * l for c in coefs]
        return [c / coefs[-1] for c in coefs]

    def piecewise_from_facets(facets, invert_coefs=False):
        pieces = []
        for facet in facets:
            # Find a plane equation with rational coefficients from facet's vertices
            # and multiply by least common multiple of denominators to obtain integer
            # coefficients.
            r = integer_coefs_from_rational(
                plane_from_points([p.vector() for p in facet]))
            if invert_coefs:
                r = [-c for c in r]
            pieces.append(AffineFunction(coefficients=r[:-1], domain=Polyhedron(
                vertices=[p.vector()[:-1] for p in facet])))
        return PiecewiseAffineFunction(affine_pieces=pieces)

    psi1_facets = []
    psi2_facets = []
    for facet in poly.facets():
        # Find out whether a facet belongs to psi1 or psi2
        coef = facet.normal_cone(direction='outer').rays_list()[0][-1]
        if coef > 0:
            psi1_facets.append(facet.vertices())
        elif coef < 0:
            psi2_facets.append(facet.vertices())
    psi1 = piecewise_from_facets(psi1_facets, invert_coefs=True)
    psi2 = piecewise_from_facets(psi2_facets)
    return CDP(psi_list=[psi1, psi2], base=Polyhedron(
        vertices=[v.vector()[:-1] for v in poly.vertices()]))
