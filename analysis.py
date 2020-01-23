#!/usr/local/bin/python3

'''Do some algebraic analysis of what happens for the bad cases'''


import sympy as sp
'''
ax, ay, bx, by, alpha, t = sp.symbols('ax ay bx by alpha t')

cx = t*ax + (1 - t)*bx + alpha*(by - ay)
cy = t*ay + (1 - t)*by + alpha*(ax - bx)

det = (ax - cx)*(by - cy) - (ay - cy)*(bx - cx)

print(det)

print(sp.simplify(det))
'''

'''
ax, ay, bx, by, cx, cy, dx, dy = sp.symbols("ax ay bx by cx cy dx dy")

adx = ax - dx
ady = ay - dy
bdx = bx - dx
bdy = by - dy
cdx = cx - dx
cdy = cy - dy

al = adx*adx + ady*ady
bl = bdx*bdx + bdy*bdy
cl = cdx*cdx + cdy*cdy


acf = al*(bdx*cdy - bdy*cdx)
bcf = bl*(ady*cdx - adx*cdy)
ccf = cl*(adx*bdy - ady*bdx)

det = acf + bcf + ccf

print(det)
print()
print(sp.simplify(sp.expand(det)))
print()
print(sp.latex(sp.collect(sp.collect(sp.collect(sp.collect(sp.collect(sp.collect(sp.simplify(sp.expand(det)), ax), ay), bx), by), cx), cy)))
print()
print(sp.collect(sp.simplify(sp.expand(det)), [ax - dx, ay - dy]))


print(sp.collect.__doc__)
'''

adx, ady, bdx, bdy, cdx, cdy = sp.symbols('adx ady bdx bdy cdx cdy')

det = (adx*adx + ady*ady)*(bdx*cdy - bdy*cdx) +(bdx*bdx + bdy*bdy)*(ady*cdx - adx*cdy) + (cdx*cdx + cdy*cdy)*(adx*bdy - ady*bdx)

print(sp.latex(sp.collect(sp.expand(det), [adx, ady, bdx, bdy, cdx, cdy])))
