#!/usr/bin/python
# A general library for computing in equivariant cohomology of the hilbert scheme.` 

from sympy import *
import math
import pickle
from fractions import Fraction

# Some Global Variables
PARTITIONS = {
  1: [[1]],
  2: [[2],[1,1]],
  3: [[3],[2,1],[1,1,1]],
  4: [[4],[3,1],[2,2],[2,1,1],[1,1,1,1]],
  5: [[5],[4,1],[3,2],[3,1,1],[2,2,1],[2,1,1,1],[1,1,1,1,1]],
  6: [[6],[5,1],[4,2],[4,1,1],[3,3],[3,2,1],[3,1,1,1],[2,2,2],[2,2,1,1],[2,1,1,1,1],[1,1,1,1,1,1]]
}

def passtr(part):
  return "".join(str(x) for x in part)

X, Y, U, V = symbols("X Y U V")
# mult two fp elements together
def mults(s,one):
  return [s * x for x in one]

def mult(one,two):
  return [simplify(a*b) for a,b in zip(one,two)]

def add(one,two):
  return [simplify(a+b) for a,b in zip(one,two)]

# compute inverse matrix or load it from file using pickle 
def load_or_compute(chob, filename):
  fn = "./pickle_stores/" + filename
  
  try:
    with open(fn,"rb") as f:
      change_of_basis = pickle.loads(f.read())
      print("Matrix loaded from "+filename+".")
  except IOError as err:
    print("No saved matrix data in "+filename+".")
    change_of_basis = cancel(chob.inv())
    with open(fn,"wb") as f:
      f.write(pickle.dumps(change_of_basis))
      print("Matrix data saved in "+filename+".")
  return change_of_basis

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
def compute_ch(list_of_fixedpoints):
  chs = []
  for fp in list_of_fixedpoints:
    roots = []
  
    for j in range(len(fp)):
      i = fp[j]
      while i > 0:
        roots.append((i-1)*V + j * U)
        i -= 1

    ch_classes = [0] * 6
    for x in roots:
      if x == 0:
        continue
      
      expansion = poly(series(exp(x),x), domain="QQ")
      while expansion:
        term = LT(expansion)
        expansion = expansion - term
      
        if not term.is_Order:
          deg = total_degree(term)
          ch_classes[deg] += term

    chs.append(ch_classes)

  chs = [list(x) for x in list(zip(*chs))[1:-1]]

  return chs
