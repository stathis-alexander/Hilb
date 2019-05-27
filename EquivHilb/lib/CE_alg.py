#
# The following functions compute the action of the operators of Nakajima on the fixed point basis
# as described in "On the equivariant cohomology of the Hilbert scheme of points in the plane" by
# Chaput and Evain. Published in Annales de L'Institut Fourier, Tome 65, p 1201-1250, 2015.
#
# A partition of a nonnegative integer N is a nonincreasing list of positive integers whose sum is N
# A monomial is represented by an (a,b) pair such that the monomial is X^aY^b
#
#

from scipy.misc import comb
from sympy import *

from lib.gen import *
from lib.partitions import PARTITION, subpartition

import math


class NAKS:
  _NEW_NAKS = False
  _NAKS = {"" : [], "1" : [(1,[1])]}

  def __init__(self):
    try:
      with open("./pickle_stores/naks","rb") as f:
        self._NAKS = pickle.loads(f.read())
    except IOError:
      pass

  def nak(self, part):
    if part in self._NAKS:
      return self._NAKS[part]

    partition = part[:]
    while partition and partition not in self._NAKS:
      if partition[0] == "-":
        partition = partition[2:]
      else:
        partition = partition[1:]

    if not partition:
      start = []
      degrees = reversed([int(x) for x in part])
      n = 0
    else:
      start = self._NAKS[partition]
      degrees = reversed([int(x) for x in part[:part.rfind(partition)]])
      n = sum([int(x) for x in partition])

    for degree in degrees:
      n = n + degree

      start = q(degree, start, PARTITION(n))
      partition = str(degree) + partition
      self._NAKS[partition] = start
    
    self._NEW_NAKS = True
    return start

  def nak_nofp(self, part):
    return [x[0] for x in self.nak(part)]

  # this should be a destructor and not a method, but destructors in python MAKE ABSOLUTELY NO SENSE and DO NOT WORK AS EXPECTED, so it's not
  def store(self):
    if self._NEW_NAKS:
      try:
        with open("./naks","wb") as f:
          f.write(pickle.dumps(self._NAKS))
      except IOError:
        pass
      
# _w
# input: an (x,y) corresponding to a monomial
# output: the weight in U,V of the monomial
def _w(x):
  return x[0] * U + x[1] * V

def _generate_M(part,part2,M):
  if part == part2:
    yield M

  if len(part) != len(part2):
    yield from _generate_M(part+[1],part2,M+[(len(part),0)])
  
  if part:
    k = part[0] + 1
    for i in range(len(part)):
      if part[i] < part2[i] and (part[i] < k):
        yield from _generate_M(part[:i]+[part[i]+1]+part[i+1:],part2,M+[(i,part[i])])
      k = part[i]

# _P_m
# input: the weight of the operator, n, and the standard skew young diagram, M
# output: the polynomial P_M of M
def _P_m(n,M):
  Pm = 0

  for i in range(1,n+1):
    summand = (-1)**(i-1) * (comb(n-1,i-1,exact=True))
    for j in range(len(M)):
      if (i-1) != j:
        summand *= _w(M[j])
    Pm += summand
  
  Pm *= Fraction((-1)**(n-1), int((math.factorial(n-1))))

  return simplify(Pm)

# _q_struct_const
# input: the weight of the operator, i, a partition, part upon which it acts, and a partition, part_2,
#        which is obtained by adding i boxes to part
# output: the structure constant of part_2 for q_i applied to part
def _q_struct_const(n,part,part_2):
  output = 0

  for M in _generate_M(part,part_2,[]):
    summand = _P_m(n,M)
    
    nu = part[:]
    for x in M:
      mu = []
      if len(nu) == x[0]:
        mu = nu + [1]
      else:
        mu = nu[:x[0]] + [nu[x[0]] + 1] + nu[x[0] + 1:]
      summand *= _q1_struct_const(nu,mu)
      nu = mu

    output = simplify(output + summand)
  
  return simplify(output)

# q
# input: the weight of the operator, i, and a partition, part, upon which it acts
# output: the expression in the fixed point bais of the action of q_i on part
def q(n,part,fps):
  return _qlist(n,part,fps)

# _clefts
# input: a partition, part
# output: an ordered list containing the cleft monomials associated to the partition
# the list is ordered such that c_{i+1}/c_{i} is divisible by Y
def _clefts(part):
  c = []
  if not part:
    return [(0,0)]

  k = part[0]+1
  for i in range(len(part)):
    if part[i] < k:
      c.append((i,part[i]))
      k = part[i]
  c.append((len(part),0))
  return c

# _new_box
# input: a partition, part, and a partition, part_2, obtained from adding one box to part
# output: a list of clefts and the index into the clefts list of the added box
def _new_box(part,part_2):
  c = _clefts(part)

  if len(part) < len(part_2):
    return c,len(c)-1
  for i in range(len(part_2)):
    if (i,part_2[i]-1) in c:
      return c,(c.index((i,part_2[i]-1)))
  
  raise ValueError("part_2",part_2,"not obtained from adding a box to part",part)

# _qlist
# input: a positive integer, a list of tuples, the first entry a coefficient and the second a paritition, a list of partitions of length pos int + len partitioons in second input
# output: q(i) applied to every element of the input

def _qlist(i,nak_list,partitions):
  if not nak_list:
    return [(_q_struct_const(i,[],lam),lam) for lam in partitions]
  
  outs = []
  for x in nak_list:
    if i < -1:
      raise ValueError("No known algorithm for q_n with n < -1.")
    elif i == -1:
      outs.append([(simplify(x[0] * _qminus1_struct_const(part,x[1])),part) for part in partitions if subpartition(part,x[1])])
    elif i == 1:
      outs.append([(simplify(x[0] * _q1_struct_const(x[1],part)),part) for part in partitions if subpartition(x[1], part)])
    else:
      outs.append([(simplify(x[0] * _q_struct_const(i,x[1],part)),part) for part in partitions if subpartition(x[1], part)])

  outclass = []
  for out in outs:
    for y in out:
      inlist = False
      for i in range(len(outclass)):
        if y[1] == outclass[i][1]:
          inlist = True
          coef = simplify(y[0] + outclass[i][0])
          outclass.pop(i)
          outclass.insert(i,(coef,y[1]))
          break

      if not inlist:
        outclass.append(y)

  return outclass

#### The following is deprecated and has been replaced by _qlist()
# _q1list
# input: a list of tuples, the first entry a coefficient and the second a partition
# output: _q1 applied to every element of the input
#def _q1list(nak_list,fps):
#  outs = []
#  for x in nak_list:
#    outs.append([(x[0] * y[0],y[1])for y in _q1(x[1],fps)])
#
#  outclass = []
#  for out in outs:
#    for y in out:
#      inlist = False
#      for i in range(len(outclass)):
#        if y[1] == outclass[i][1]:
#          inlist = True
#          coef = simplify(y[0] + outclass[i][0])
#          outclass.pop(i)
#          outclass.insert(i,(coef,y[1]))
#          break
#
#      if not inlist:
#        outclass.append(y)
#
#  return outclass

# _q1
# input: a partition, part
# output: the expression in the fixed point basis of the action of q_1 on part
def _q1(part,fps):
  return [(_q1_struct_const(part,x),x) for x in fps]

# _q1_struct_const
# input: a partition, part, and a partition, part_2, obtained by adding one box to part
# output: the structure constant for part_2 of q_1 applied to part 
def _q1_struct_const(part, part_2):
  return simplify(_coker(part,part_2)/_ker(part,part_2))

# _qminus1
# input: a partition, mu, a list of fixed points of the hilbert scheme of one fewer points
# output: the expression in the fixed point basis of the action of qminus1 on mu
def _qminus1(mu,fps):
  return [(_qminus1_struct_const(x,mu),x) for x in fps]

# _qminus1_struct_const
# input: a partition, mu, and a partition, lam, obtained by removing a box from mu 
# output: the structure constant for part_2 of q_-1 applied to part. 
def _qminus1_struct_const(lam,mu):
  return simplify(_tangent(lam)/_tangent(mu)*_coker(lam,mu)/_ker(lam,mu))

# _tangent
# input: a partition, part
# output: the product of the weights of the tangent space to the fixed point corresponding to part
def _tangent(part):
  output = 1
  for i in range(len(part)):
    for j in range(part[i]):
      l = [x for x in part[(i+1):] if j < x]
      
      a = len(l)
      b = part[i] - 1 - j
      output *= (a*U - (b+1)*V)*(-(a+1)*U+b*V)

  return output

# _coker
# input: a partition, part, and a partition, part_2, obtained by adding one box to part
# output: the numerator of the structure constant for part_2 of q_1 applied to part 
def _coker(part,part_2):
  c,k = _new_box(part,part_2)

  output = 1
  for i in range(len(c)):
    if i != k:
      output = output * ((c[i][0] - c[k][0] -1)*U + (c[i][1]-c[k][1]-1)*V)
  return simplify(output)

# _ker
# input: a partition, part, and a partition, part_2, obtained by adding one box to part
# output: the denominator of the structure constant for part_2 of q_1 applied to part 
def _ker(part,part_2):
  c,k = _new_box(part,part_2)

  output = 1
  for i in range(len(c)):
    if i > k:
      output = output * (-(-c[i-1][1] + c[k][1] + 1)*V + (c[i][0]-c[k][0]-1)*U)
    if i < k:
      output = output * (-(c[k][0] - c[i+1][0] + 1)*U + (c[i][1]-c[k][1]-1)*V)

  return simplify(output)



