#!/usr/bin/env sage
from sage.all import *
from piecewise_affine_function import AffineFunction, PiecewiseAffineFunction
from cdp import CDP


# TODO: remove copy-paste polyhedron


class TestCDPTransform:

    def test_shear(self):
        base = Polyhedron(vertices=[[-1], [1]])
        # y = 1 + x, x \in [-1, 0]
        f_11 = AffineFunction([1, 1], Polyhedron(vertices=[[-1], [0]]))
        # y = 1 - x, x \in [0, 1]
        f_12 = AffineFunction([1, -1], Polyhedron(vertices=[[0], [1]]))
        f_1 = PiecewiseAffineFunction([f_11, f_12])
        # y = 1/2x + 1/2, x \in [-1, 1]
        f_2 = PiecewiseAffineFunction([AffineFunction([1 / 2, 1 / 2], Polyhedron(vertices=[[-1], [1]]))])
        cdp = CDP([f_1, f_2], base)
        cdp.shear([-1, 1], [-1])
        assert cdp.psi_list[0].affine_pieces[0].coefs == [1, 2]
        assert cdp.psi_list[1].affine_pieces[0].coefs == [0.5, -0.5]

    def test_translate(self):
        base = Polyhedron(vertices=[[-1], [1]])
        # y = 1 + x, x \in [-1, 0]
        f_11 = AffineFunction([1, 1], Polyhedron(vertices=[[-1], [0]]))
        # y = 1 - x, x \in [0, 1]
        f_12 = AffineFunction([1, -1], Polyhedron(vertices=[[0], [1]]))
        f_1 = PiecewiseAffineFunction([f_11, f_12])
        # y = 1/2x + 1/2, x \in [-1, 1]
        f_2 = PiecewiseAffineFunction([AffineFunction([1 / 2, 1 / 2], Polyhedron(vertices=[[-1], [1]]))])
        cdp = CDP([f_1, f_2], base)
        cdp.shear([-1, 1], [-1])
        cdp.translate([1, -1])
        assert cdp.psi_list[0].affine_pieces[0].coefs == [2, 2]
        assert cdp.psi_list[1].affine_pieces[0].coefs == [-0.5, -0.5]

    def test_transform_base(self):
        base = Polyhedron(vertices=[[-1], [1]])
        # y = 1 + x, x \in [-1, 0]
        f_11 = AffineFunction([1, 1], Polyhedron(vertices=[[-1], [0]]))
        # y = 1 - x, x \in [0, 1]
        f_12 = AffineFunction([1, -1], Polyhedron(vertices=[[0], [1]]))
        f_1 = PiecewiseAffineFunction([f_11, f_12])
        # y = 1/2x + 1/2, x \in [-1, 1]
        f_2 = PiecewiseAffineFunction([AffineFunction([1 / 2, 1 / 2], Polyhedron(vertices=[[-1], [1]]))])
        cdp = CDP([f_1, f_2], base)
        cdp.shear([-1, 1], [-1])
        cdp.translate([1, -1])
        A = matrix(ZZ, [[-1]])
        phi = linear_transformation(A)
        cdp.transform_base(phi)
        assert cdp.psi_list[0].affine_pieces[0].coefs == [2, -2]
        domain_verts = [list(vert) for vert in cdp.psi_list[0].affine_pieces[0].domain.vertices()]
        assert domain_verts == [[0], [1]]
        assert cdp.psi_list[1].affine_pieces[0].coefs == [-0.5, 0.5]


test = TestCDPTransform()
test.test_shear()
test.test_translate()
test.test_transform_base()

