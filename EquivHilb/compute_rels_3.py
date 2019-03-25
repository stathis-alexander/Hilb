#!/usr/bin/python
# This is a template file for a running code using the CE_alg library.

from CE_alg import *
from partitions import *

matrix_file = "rels3"

fps = [[3],[2,1],[1,1,1]]

nak111 = [1,1,1]
nak21 = [-6*V,-2*(U+V),-6*U]
nak3 = [6*V**2,3*U*V,6*U**2]

chob = Matrix([nak111,nak21,nak3])
change_of_basis = load_or_compute(chob,matrix_file)

chs = compute_ch(fps)
print_chs(chs,fps)

for i in range(len(chs)):
  print("ch"+str(i+1)+":")
  print(simplify( Matrix([chs[i]]) * change_of_basis ))
