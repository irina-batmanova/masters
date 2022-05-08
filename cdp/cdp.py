#!/usr/bin/env sage
from sage.all import *
from typing import List
from piecewise_affine_function import PiecewiseAffineFunction, AffineFunction
from collections import defaultdict
from itertools import permutations
import numpy as np
from sympy.matrices import Matrix
from sympy import Rational, lcm


class CDP:
    def __init__(self, psi_list: List[PiecewiseAffineFunction], base: Polyhedron):
        # Check that sum of psi is non negative on borders of domains
        # (check polytop vertices of domains laying inside of function scope)
        for psi in psi_list:
            for piece in psi.affine_pieces:
                for vert in piece.domain.vertices():
                    if vert in base:
                        try:
                            s = sum([ps.value(vert.vector()) for ps in psi_list])
                        except ValueError:
                            continue
                        else:
                            # TODO: На внутренних точках должно быть строго > 0,
                            # здесь надо проверить, что точка не находится на грани base
                            if s < 0:
                                raise ValueError(
                                    f'Not a valid CDP - sum of psi is {s} on {vert}')
        # Check that sum of psi is non negative on base vertices
        for vert in base.vertices():
            try:
                s = sum([ps.value(vert) for ps in psi_list])
            except ValueError:
                raise ValueError('Not a valid CDP: psi is not defined on base')
            else:
                if s < 0:
                    raise ValueError(
                        f'Not a valid CDP - sum of psi is {s} on {vert}')
        self.psi_list = psi_list
        self.base = base
        self.n = len(self.base.vertices()[0].vector())
        self.k = len(self.base.vertices())
        # TODO: is there a better solution for vertices traversal?
        self.base_adjacency_map = defaultdict(set)

    def __str__(self):
        psi_str = "\n".join([str(psi) for psi in self.psi_list])
        return f'CDP object, psi list: {psi_str},\nbase: {self.base.vertices()}'

    def __eq__(self, other):
        """
        Note: the order of the functions in psi_list is not important for this implementation.
        """
        if not self.base == other.base:
            return False
        if not len(self.psi_list) == len(other.psi_list):
            return False
        for psi in self.psi_list:
            if psi not in other.psi_list:
                return False
        return True

    def _build_adjacency_map(self):
        for i in range(self.k):
            for j in range(i + 1, self.k):
                if self.base.vertices()[i].is_incident(self.base.vertices()[j]):
                    self.base_adjacency_map[i].add(j)
                    self.base_adjacency_map[j].add(i)

    def transform_base(self, phi: linear_transformation):
        vertices = []
        for vert in self.base.vertices():
            vertices.append(phi(vert.vector()))
        self.base = Polyhedron(vertices=vertices)
        try:
            inv = phi.inverse().matrix()
        except ZeroDivisionError:
            raise ValueError(f'phi is not invertible')
        inv = np.array([np.array(row) for row in inv])
        # print("invert ", inv)
        for i in range(len(self.psi_list)):
            self.psi_list[i].transform(phi, inv)

    def translate(self, alpha_list: List[int]):
        if len(alpha_list) != len(self.psi_list):
            raise ValueError(f'Length of alpha_list is {len(alpha_list)}, '
                             f'should be {len(self.psi_list)}')
        if sum(alpha_list) != 0:
            raise ValueError('Sum of coefficients should be 0')
        for idx, psi in enumerate(self.psi_list):
            for piece in psi.affine_pieces:
                piece.coefs[0] += alpha_list[idx]

    def shear(self, beta_list: List[int], v: List[int]):
        if len(beta_list) != len(self.psi_list):
            raise ValueError(f'Length of beta_list is {len(beta_list)}, '
                             f'should be {len(self.psi_list)}')
        if sum(beta_list) != 0:
            raise ValueError('Sum of coefficients should be 0')
        m = len(self.psi_list[0].affine_pieces[0].coefs) - 1
        if len(v) != m:
            raise ValueError(f'Wrong dimension of v: {len(v)}, should be {m}')
        for idx, psi in enumerate(self.psi_list):
            for piece in psi.affine_pieces:
                for j, coef in enumerate(v):
                    piece.coefs[j + 1] += coef * beta_list[idx]

    def _vert_permutation_is_valid(self, perm):
        """Check that vertices permutation saves incidence
        """
        for vert, inc_list in self.base_adjacency_map.items():
            for other_vert in inc_list:
                if not perm[other_vert] in self.base_adjacency_map[perm[vert]]:
                    return False
        return True

    def _get_transform_matrix(self, points, point_images):
        left = np.matmul(points.T, points)
        try:
            inverse = np.linalg.inv(left)
        except np.linalg.LinAlgError:
            raise ValueError(f'Base of the CDP should be non-degenerate '
                             f'(dimension {self.n} is provided, but '
                             'the real dimension is lower)')
        A = np.matmul(np.matmul(point_images, points), inverse)
        return A

    def _domains_match(self, psi1: PiecewiseAffineFunction, psi2: PiecewiseAffineFunction) -> bool:
        if len(psi1.affine_pieces) != len(psi2.affine_pieces):
            return False
        return set([p.domain for p in psi1.affine_pieces]) == set([p.domain for p in psi2.affine_pieces])

    def _get_equivalence_classes(self, other_psi_list):
        used = [False for _ in range(len(other_psi_list))]
        classes = [[] for _ in range(len(self.psi_list))]
        for i, psi1 in enumerate(self.psi_list):
            for j, psi2 in enumerate(other_psi_list):
                if self._domains_match(psi1, psi2):
                    used[j] = True
                    classes[i].append(j)
        if len([c for c in classes if len(c) == 0]) > 0 or len([u for u in used if not u]) > 0:
            return False, []
        return True, classes

    def _list_mappings(self, classes):
        dicts_res = self._list_mappings_recursive(classes, dict(), idx=0)
        res = []
        for small_res in dicts_res:
            res.append([0] * len(classes))
            for k, v in small_res.items():
                res[-1][v] = k
        return res

    def _list_mappings_recursive(self, classes, current_dict, idx=0):
        res = []
        for elem in classes[idx]:
            if elem in current_dict:
                continue
            if idx + 1 == len(classes):
                d = deepcopy(current_dict)
                d[elem] = idx
                return [d]
            d = deepcopy(current_dict)
            d[elem] = idx
            res += self._list_mappings_recursive(classes, d, idx + 1)
        return res

    def _can_be_translated(self, mapping, other_psi_list):
        alpha_list = []
        for i, j in enumerate(mapping):
            can, alpha = self.psi_list[i].can_be_translated(other_psi_list[j])
            if not can:
                return False
            alpha_list.append(alpha)
        if sum(alpha_list) != 0:
            return False
        return True

    def _can_be_sheared(self, mapping, other_psi_list):
        for i, j in enumerate(mapping):
            if not self.psi_list[i].cat_be_sheared(other_psi_list[j]):
                return False
        return True

    def equal(self, other_cdp):
        """
        Suppose that zero functions are already removed
        """
        # Check that self.base is convertible to other_cdp.base with some phi
        # Check that psi lists split into matching equivalence classes
        # Check that a constant for translation exists
        # Check that a vector for shearing exists
        if len(self.base.vertices()) != len(other_cdp.base.vertices()):
            return False
        if len(self.psi_list) != len(other_cdp.psi_list):
            return False
        G = np.array([np.array(vert.vector()) for vert in other_cdp.base.vertices()])
        G = G.T
        for perm in permutations([i for i in range(self.k)]):
            if not self._vert_permutation_is_valid(perm):
                continue
            V = np.array([np.array(self.base.vertices()[i].vector()) for i in perm])
            # A transforms base of one CDP to the base of another
            A = self._get_transform_matrix(V, G)
            A = linear_transformation(matrix(QQ, A))
            cdp_after_base_transform = deepcopy(self)
            try:
                cdp_after_base_transform.transform_base(A)
            except ValueError:
                continue
            splits, classes = cdp_after_base_transform._get_equivalence_classes(other_cdp.psi_list)
            if not splits:
                continue
            mappings = cdp_after_base_transform._list_mappings(classes)
            for m in mappings:
                can = cdp_after_base_transform._can_be_translated(m, other_cdp.psi_list)
                if not can:
                    continue
                can = cdp_after_base_transform._can_be_sheared(m, other_cdp.psi_list)
                if not can:
                    continue
                return True
        return False


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
            r = integer_coefs_from_rational(plane_from_points([p.vector() for p in facet]))
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
    return CDP(psi_list=[psi1, psi2], base=Polyhedron(vertices=[v.vector()[:-1] for v in poly.vertices()]))

