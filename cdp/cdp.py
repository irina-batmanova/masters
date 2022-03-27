#!/usr/bin/env sage
from sage.all import *

from piecewise_affine_function import AffineFunction


f = AffineFunction([-1, 1], Polyhedron(vertices=[[-1], [1]]))
print(f)
