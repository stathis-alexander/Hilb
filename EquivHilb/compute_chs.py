#!/usr/bin/python
# compute the chern characters in terms of nak basis for n points (up to 5) using CE's algorithm at fixed points

from CE_alg import *
from partitions import *

for num_points in range(2,7):
  matrix_file = "rels" + str(num_points)

  nak = NAKS()

  chob = Matrix([nak.nak_nofp(x) for x in [passtr(P) for P in PARTITIONS[num_points]]])
  change_of_basis = load_or_compute(chob,matrix_file)

  chs = compute_ch(PARTITIONS[num_points])
  for i in range(len(chs)):
     print("ch"+str(i+1)+":")
     print(simplify( Matrix([chs[i]]) * change_of_basis ))

