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
from partitions import *
import math

X, Y, U, V = symbols("X Y U V")

# w
# input: an (x,y) corresponding to a monomial
# output: the weight in U,V of the monomial
def w(x):
  return x[0] * U + x[1] * V

def generate_M(part,part2,M):
  if part == part2:
    yield M

  if len(part) != len(part2):
    yield from generate_M(part+[1],part2,M+[(len(part),0)])
  
  if part:
    k = part[0] + 1
    for i in range(len(part)):
      if part[i] < part2[i] and (part[i] < k):
        yield from generate_M(part[:i]+[part[i]+1]+part[i+1:],part2,M+[(i,part[i])])
      k = part[i]

# P_m
# input: the weight of the operator, n, and the standard skew young diagram, M
# output: the polynomial P_M of M
def P_m(n,M):
  Pm = 0

  for i in range(1,n+1):
    summand = (-1)**(i-1) * (comb(n-1,i-1,exact=True))
    for j in range(len(M)):
      if (i-1) != j:
        summand *= w(M[j])
    Pm += summand
  
  Pm *= (-1)**(n-1) / int((math.factorial(n-1)))

  return simplify(Pm)

# q_struct_const
# input: the weight of the operator, i, a partition, part upon which it acts, and a partition, part_2,
#        which is obtained by adding i boxes to part
# output: the structure constant of part_2 for q_i applied to part
def q_struct_const(n,part,part_2):
  output = 0

  for M in generate_M(part,part_2,[]):
    summand = P_m(n,M)
    nu = part[:]
    for x in M:
      mu = []
      if len(nu) == x[0]:
        mu = nu + [1]
      else:
        mu = nu[:x[0]] + [nu[x[0]] + 1] + nu[x[0] + 1:]
      summand *= q1_struct_const(nu,mu)
      nu = mu

    output+= summand
  
  return simplify(output)

# q
# input: the weight of the operator, i, and a partition, part, upon which it acts
# output: the expression in the fixed point bais of the action of q_i on part
def q(n,part):
# generate new partitions
# compute structure constant for each partition pair
  return []

# clefts
# input: a partition, part
# output: an ordered list containing the cleft monomials associated to the partition
# the list is ordered such that c_{i+1}/c_{i} is divisible by Y
def clefts(part):
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

# new_box
# input: a partition, part, and a partition, part_2, obtained from adding one box to part
# output: a list of clefts and the index into the clefts list of the added box
def new_box(part,part_2):
  c = clefts(part)

  if len(part) < len(part_2):
    return c,len(c)-1
  for i in range(len(part_2)):
    if (i,part_2[i]-1) in c:
      return c,(c.index((i,part_2[i]-1)))
  
  raise ValueError("part_2",part_2,"not obtained from adding a box to part",part)

# q1
# input: a partition, part
# output: the expression in the fixed point basis of the action of q_1 on part
def q1(part):
  return [(q1_struct_const(part,x),x) for x in gen_partitions_1(part)]

# q1_struct_const
# input: a partition, part, and a partition, part_2, obtained by adding one box to part
# output: the structure constant for part_2 of q_1 applied to part 
def q1_struct_const(part, part_2):
  return simplify(coker(part,part_2)/ker(part,part_2))


# coker
# input: a partition, part, and a partition, part_2, obtained by adding one box to part
# output: the numerator of the structure constant for part_2 of q_1 applied to part 
def coker(part,part_2):
  c,k = new_box(part,part_2)

  output = 1
  for i in range(len(c)):
    if i != k:
      output = output * ((c[i][0] - c[k][0] -1)*U + (c[i][1]-c[k][1]-1)*V)
  return simplify(output)

# coker
# input: a partition, part, and a partition, part_2, obtained by adding one box to part
# output: the denominator of the structure constant for part_2 of q_1 applied to part 
def ker(part,part_2):
  c,k = new_box(part,part_2)

  output = 1
  for i in range(len(c)):
    if i > k:
      output = output * (-(c[i-1][1] - c[k][1] + 1)*V + (c[i][0]-c[k][0]-1)*U)
    if i < k:
      output = output * (-(c[k][0] - c[i+1][0] + 1)*U + (c[i][1]-c[k][1]-1)*V)

  return simplify(output)



