#!/usr/bin/env sage
from typing import List
from sage.all import *


class AffineFunction:
    def __init__(self, coefficients: List[int], domain: Polyhedron):
        dim = len(domain.vertices()[0].vector())
        if len(coefficients) != dim + 1:
            raise ValueError(f'Domain dimension {len(coefficients) - 1} '
                             f'and coefficients list dimension {dim} do not match')
        self.coefs = list(coefficients)
        self.domain = domain
        self.dim = dim

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

    def value(self, x: List[int]):
        for piece in self.affine_pieces:
            if x in piece.domain:
                return piece.value(x)
        raise ValueError(f'{x} is not in function domain')

    # TODO: аргументы похожи по смыслу, но имеют разный тип - поправить
    def transform(self, phi, phi_inverse):
        for j in range(len(self.affine_pieces)):
            res = self.affine_pieces[j].coefs[1:] * phi_inverse
            res = res.tolist()[0]
            res = [self.affine_pieces[j].coefs[0]] + res
            self.affine_pieces[j].coefs = res
            vertices = []
            for vert in self.affine_pieces[j].domain.vertices():
                vertices.append(phi(vert.vector()))
            self.affine_pieces[j].domain = Polyhedron(vertices=vertices)

    def __str__(self):
        resp = '\n'.join([str(piece) for piece in self.affine_pieces])
        return 'Piecewise affine function:\n' + resp

