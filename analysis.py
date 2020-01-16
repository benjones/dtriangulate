#!/usr/local/bin/python3

'''Do some algebraic analysis of what happens for the bad cases'''


import sympy as sp

ax, ay, bx, by, alpha, t = sp.symbols('ax ay bx by alpha t')

cx = t*ax + (1 - t)*bx + alpha*(by - ay)
cy = t*ay + (1 - t)*by + alpha*(ax - bx)

det = (ax - cx)*(by - cy) - (ay - cy)*(bx - cx)

print(det)

print(sp.simplify(det))
