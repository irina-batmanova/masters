#!/usr/bin/env sage
from typing import List
from sage.all import *
import numpy as np


class AffineFunction:
    def __init__(self, coefficients: List[float], domain: Polyhedron):
        dim = len(domain.vertices()[0].vector())
        if len(coefficients) != dim + 1:
            raise ValueError(f'Domain dimension {len(coefficients) - 1} '
                             f'and coefficients list dimension {dim} do not match')
        self.coefs = list(coefficients)
        self.domain = domain
        self.dim = dim

    def __eq__(self, other):
        return self.coefs == other.coefs and self.domain == other.domain

    def value(self, x: List[int]):
        if x not in self.domain:
            raise ValueError(f'{x} is not in function domain {self.domain.vertices()}')
        if len(x) != self.dim:
            raise ValueError(f'Wrong dimension of x: {len(x)}')
        return sum([a * b for a, b in zip(self.coefs[1:], x)]) + self.coefs[0]

    def _coef_str(self, i):
        if self.coefs[i] == 0:
            return ''
        if self.coefs[i] == 1:
            return f' + x_{i}'
        elif self.coefs[i] == -1:
            return f' - x_{i}'
        elif self.coefs[i] > 0:
            return f' + {self.coefs[i]}x_{i}'
        return f' - {-self.coefs[i]}x_{i}'

    def __str__(self):
        coefs_repr = []
        if self.coefs[0] != 0 or len(self.coefs) == 2 and self.coefs[1] == 0:
            coefs_repr.append(str(self.coefs[0]))
        for i in range(1, len(self.coefs)):
            coefs_repr.append(self._coef_str(i))
        coefs_repr = "".join(coefs_repr)
        return f'Affine function {coefs_repr} with domain ' \
               f'{[vert.vector() for vert in self.domain.vertices()]}'


class PiecewiseAffineFunction:
    def __init__(self, affine_pieces: List[AffineFunction]):
        # TODO: check regularity
        # TODO: check convexity
        # TODO: проверить, что значения в вершинах многогранников одинаковые на смежных
        # кусках и целые
        self.affine_pieces = affine_pieces

    def __eq__(self, other):
        if not len(self.affine_pieces) == len(other.affine_pieces):
            return False
        for piece in self.affine_pieces:
            if not piece in other.affine_pieces:
                return False
        return True

    def value(self, x: List[int]):
        for piece in self.affine_pieces:
            if x in piece.domain:
                return piece.value(x)
        raise ValueError(f'{x} is not in function domain')

    # TODO: аргументы похожи по смыслу, но имеют разный тип - поправить
    def transform(self, phi, phi_inverse):
        for j in range(len(self.affine_pieces)):
            res = np.matmul(self.affine_pieces[j].coefs[1:], phi_inverse)
            self.affine_pieces[j].coefs = [self.affine_pieces[j].coefs[0]]
            self.affine_pieces[j].coefs.extend(res)
            vertices = []
            for vert in self.affine_pieces[j].domain.vertices():
                vertices.append(phi(vert.vector()))
            self.affine_pieces[j].domain = Polyhedron(vertices=vertices)

    def _domains_mapping(self, other_psi):
        if len(self.affine_pieces) != len(other_psi.affine_pieces):
            return False, 0
        domains_mapping = [0 for _ in range(len(other_psi.affine_pieces))]
        for i, piece in enumerate(self.affine_pieces):
            for j, other_piece in enumerate(other_psi.affine_pieces):
                if piece.domain == other_piece.domain:
                    domains_mapping[i] = j
        return domains_mapping

    def can_be_translated(self, other_psi):
        domains_mapping = self._domains_mapping(other_psi)
        alpha = other_psi.affine_pieces[domains_mapping[0]].coefs[0] - self.affine_pieces[0].coefs[0]
        for i in range(1, len(domains_mapping)):
            a = other_psi.affine_pieces[domains_mapping[i]].coefs[0] - self.affine_pieces[i].coefs[0]
            if a != alpha:
                return False, 0
        return True, alpha

    def cat_be_sheared(self, other_psi):
        domains_mapping = self._domains_mapping(other_psi)
        n = len(other_psi.affine_pieces[domains_mapping[0]].coefs)
        v = [other_psi.affine_pieces[domains_mapping[0]].coefs[i] - self.affine_pieces[0].coefs[i]
             for i in range(1, n)]
        for j in range(1, len(domains_mapping)):
            other_v = [other_psi.affine_pieces[domains_mapping[j]].coefs[i] - self.affine_pieces[j].coefs[i]
                for i in range(1, n)]
            if other_v != v:
                return False
        return True

    def __str__(self):
        resp = '\n'.join([str(piece) for piece in self.affine_pieces])
        return 'Piecewise affine function:\n' + resp

