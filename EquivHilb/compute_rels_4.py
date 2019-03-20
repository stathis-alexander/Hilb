#!/usr/bin/python
# This is a template file for a running code using the CE_alg library.

from CE_alg import *
from partitions import *

nak1111 = [1,1,1,1,1]
nak211 = [-24*V,-4*(U+3*V),-8*(U+V),-4*(3*U+V),-24*U]
nak31 = [24*V**2,6*V*(U+V),12*U*V,6*U*(U+V),24*U**2]
nak22 = [24*V**2,8*U*V,8*(U**2-U*V+V**2),8*U*V,24*U**2]
nak4 = [-24*V**3,-8*U*V**2,-4*U*V*(U+V),-8*U**2*V,-24*U**3]

def mult(one,two):
   return [simplify(a*b) for a,b in zip(one,two)]

chob = Matrix([nak1111,nak211,nak31,nak22,nak4])
for entry in chob:
   print(entry)

change_of_basis = simplify(chob.inv())

product1 = mult(nak211,nak211)

rel1 = Matrix([product1]) * change_of_basis
print(simplify(rel1))

