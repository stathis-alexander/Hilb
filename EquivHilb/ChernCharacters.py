#!/usr/bin/python
# compute the chern characters in terms of nak basis for n points (up to 5) using CE's algorithm at fixed points

from lib.CE_alg import *
from lib.partitions import *

nstart = 4
nstop = 4
prec = 6

# print a list the chern classes at each fixed point given the chs and fixed points in a list
def print_chs(chs,fps):

  lengths=[len(fp) for fp in fps]
  for ch in chs:
    for x in ch:
      lengths += [len(str(x))]
  tabs = int(math.ceil(max(lengths)/8))

  print(end="\t")
  for x in fps:
    space = int(math.floor(len(x)/8))
    print(str(x), end=("\t"*(tabs-space)))
  print()

  for i,ch in enumerate(chs):
    print(str(i+1),end="\t")
    for x in ch:
      space = int(math.floor(len(str(x))/8))
      print(str(x), end=("\t"*(tabs-space)))
    print()

# compute the restrictions of the chern classes of the tautological bundles at each fixed point of a list of fixed points given as partitions
# up to precision degree
def compute_ch(list_of_fixedpoints,precision):
  chs = []
  for fp in list_of_fixedpoints:
    roots = []
  
    for j in range(len(fp)):
      i = fp[j]
      while i > 0:
        roots.append((i-1)*V + j * U)
        i -= 1

    ch_classes = [0] * (precision+1)
    for x in roots:
      if x == 0:
        continue
      
      expansion = poly(series(exp(x),x,0,precision+1), domain="QQ")
      while expansion:
        term = LT(expansion)
        expansion = expansion - term
      
        if not term.is_Order:
          deg = total_degree(term)
          ch_classes[deg] += term

    chs.append(ch_classes)
  chs = [list(x) for x in list(zip(*chs))[1:]]
  return chs

#print and generate

for num_points in range(nstart,nstop+1):
  print("Printing Chern Characters for "+str(num_points)+" points.")
  
  matrix_file = "rels" + str(num_points)

  nak = NAKS()

  chob = Matrix([nak.nak_nofp(x) for x in [passtr(P) for P in PARTITION(num_points)]])
  change_of_basis = load_or_compute(chob,matrix_file)

  chs = compute_ch(PARTITION(num_points),prec)
  for i in range(len(chs)):
     print("ch"+str(i+1)+":")
     coefs = list(simplify( Matrix([chs[i]]) * change_of_basis ))
     print([x for x in zip(coefs,PARTITION(num_points))])

