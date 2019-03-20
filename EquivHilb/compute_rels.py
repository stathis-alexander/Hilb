#!/usr/bin/python
# This is a template file for a running code using the CE_alg library.

from CE_alg import *
from partitions import *

nak111 = [1,1,1]
nak21 = [-6*V,-2*(U+V),-6*U]
nak3 = [6*V**2,3*U*V,6*U**2]

def mult(one,two):
   return [simplify(a*b) for a,b in zip(one,two)]

chob = Matrix([nak111,nak21,nak3])
change_of_basis = simplify(chob.inv())

product1 = mult(nak21,nak21)
print(product1)

rel1 = Matrix([product1]) * change_of_basis
print(simplify(rel1))

