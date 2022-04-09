#!/usr/bin/env sage
from sage.all import *
from typing import List
from piecewise_affine_function import PiecewiseAffineFunction
from collections import defaultdict
from itertools import permutations
import numpy as np


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

    def equal(self, other_cdp):
        # Check that self.base is convertible to other_cdp.base with some phi
        # Check that sum(psi_i) stays the same on whole base (check on vertices)
        if len(self.base.vertices()) != len(other_cdp.base.vertices()):
            print("lengths are different")
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
            all_sums_are_equal = True
            for vert in cdp_after_base_transform.base.vertices():
                sum1 = sum([p.value(vert) for p in cdp_after_base_transform.psi_list])
                sum2 = sum([p.value(vert) for p in other_cdp.psi_list])
                if sum1 != sum2:
                    all_sums_are_equal = False
            if all_sums_are_equal:
                return True
        return False
