#!/usr/bin/python
# This is a template file for a running code using the CE_alg library.

from CE_alg import *
from partitions import *

matrix_file = "rels4"

nak1111 = [1,1,1,1,1]
nak211 =  [-24*V,-4*(U+3*V),-8*(U+V),-4*(3*U+V),-24*U]
nak31 =   [24*V**2,6*V*(U+V),12*U*V,6*U*(U+V),24*U**2]
nak22 =   [24*V**2,8*U*V,8*(U**2-U*V+V**2),8*U*V,24*U**2]
nak4 =    [-24*V**3,-8*U*V**2,-4*U*V*(U+V),-8*U**2*V,-24*U**3]

fps = [[4],[3,1],[2,2],[2,1,1],[1,1,1,1]]

chob = Matrix([nak1111,nak211,nak31,nak22,nak4])
change_of_basis = load_or_compute(chob,matrix_file)

chs = compute_ch(fps)

for i in range(len(chs)):
  print("ch"+str(i+1)+":")
  print(simplify( Matrix([chs[i]]) * change_of_basis ))
