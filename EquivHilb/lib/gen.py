#!/usr/bin/python
# A general library for computing in equivariant cohomology of the hilbert scheme.` 

from sympy import *
import math
import pickle
from fractions import Fraction

from lib.partitions import PARTITION

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

