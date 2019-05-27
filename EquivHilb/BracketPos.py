#!/usr/bin/python
#
# This file computes the bracket [q_n',q_m](q_1^k(1)) in terms of thei
# Nakajima classes. 1 is the vacuum vector. See 3.10 in Lehn's Chern 
# Classes paper.
#
# This is optimized for positive n,m, and is slightly faster than
# BracketNeg. 
#

from lib.CE_alg import *
from lib.partitions import PARTITION

### Assume n,m > 0
n = 2
m = 1
k = 2

num_points = n + m + k
matrix_file = "rels" + str(num_points)

nak = NAKS()

chob = Matrix([nak.nak_nofp(x) for x in [passtr(P) for P in PARTITION(num_points)]])
change_of_basis = load_or_compute(chob,matrix_file)

A = nak.nak_nofp("2" + "1" * (n+m+k-2))
B = nak.nak_nofp(str(max(n,m)) + str(min(n,m)) + "1" * k)
C = nak.nak_nofp("2" + "1" * (m+k-2))
D = nak.nak_nofp(str(m) + "1" * k)
E = nak.nak_nofp("2" + "1" * (n+k-2))
F = nak.nak_nofp(str(n) + "1" * k)

ints = sorted([m,n,2])
G = nak.nak_nofp(str(ints[0]) + str(ints[1]) + str(ints[2]) + "1" * (k-2))

sums = [mult(A,B),mult(C,D),mult(E,F),G]

sums[1] = [x[0] for x in q(n,zip(sums[1],PARTITION(m+k)),PARTITION(num_points))]
sums[2] = [x[0] for x in q(m,zip(sums[2],PARTITION(n+k)),PARTITION(num_points))]

if n+m+k-2 < 0:
  sums[0] = [0] * len(PARTITION(num_points))
else:
  sums[0] = simplify( Matrix([sums[0]]) * change_of_basis )
if m+k-2 < 0 :
  sums[1] = [0] * len(PARTITION(num_points))
else:  
  sums[1] = simplify( Matrix([sums[1]]) * change_of_basis )
if n+k-2 < 0 :
  sums[2] = [0] * len(PARTITION(num_points))
else:  
  sums[2] = simplify( Matrix([sums[2]]) * change_of_basis )
if k-1 <= 0:
  sums[3] = [0] * len(PARTITION(num_points))
else:
  sums[3] = simplify( Matrix([sums[3]]) * change_of_basis )


sums[0] = mults(Fraction(1,math.factorial(n+m+k-2)), sums[0])
sums[1] = mults(Fraction(-1,math.factorial(m+k-2)), sums[1])
sums[2] = mults(Fraction(-1,math.factorial(n+k-2)), sums[2])
sums[3] = mults(k * (k-1), sums[3])

for s in sums:
  print(s)

output = add(sums[0],add(sums[1],add(sums[2],sums[3])))
print("n=",n,"m=",m,"k=",k)
print([x for x in zip(mults(Fraction(-1,2), output),PARTITION(num_points))])

