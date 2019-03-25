#!/usr/bin/python
# This is a template file for a running code using the CE_alg library.

from CE_alg import *
from partitions import *

matrix_file = "rels2"

nak11 = [1,1]
nak2 = [-2*V,-2*U]

fps = [[2],[1,1]]

chob = Matrix([nak11,nak2])
change_of_basis = load_or_compute(chob,matrix_file)

chs = compute_ch(fps)
print_chs(chs,fps)

for i in range(len(chs)):
  print("ch"+str(i+1)+":")
  print(simplify( Matrix([chs[i]]) * change_of_basis ))


