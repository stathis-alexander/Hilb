#!/usr/bin/python
# compute the chern characters in terms of nak basis for n points (up to 5) using CE's algorithm at fixed points

from CE_alg import *
from partitions import *

## assume n > 0

n = 1
m = -1
k = 5
num_points = n + m + k
matrix_file = "rels" + str(num_points)

nak = NAKS()

chob = Matrix([nak.nak_nofp(x) for x in [passtr(P) for P in PARTITIONS[num_points]]])
change_of_basis = load_or_compute(chob,matrix_file)

A = [x[0] for x in q(n,q(m,nak.nak("1" * k),PARTITIONS[k+m]),PARTITIONS[n+m+k])]
B = nak.nak_nofp("2" + "1" * (n+m+k-2))

C = [x[0] for x in q(n,q(m,nak.nak("2" + "1" * (k-2)),PARTITIONS[k+m]),PARTITIONS[n+m+k])]

D = [x[0] for x in q(m,nak.nak("1" * k),PARTITIONS[k+m])]
D = mult(nak.nak_nofp("2" + "1" * (k+m-2)),D)
D = [x[0] for x in q(n,zip(D,PARTITIONS[k+m]),PARTITIONS[n+m+k])]
E = mult(nak.nak_nofp("2" + "1" * (k+n-2)),nak.nak_nofp(str(n) + "1" * k))
E = [x[0] for x in q(m,zip(E,PARTITIONS[k+n]),PARTITIONS[n+m+k])]

sums = [mult(A,B),D,E,C]

if n+m+k-2 < 0:
  sums[0] = [0] * len(PARTITIONS[num_points])
else:
  sums[0] = simplify( Matrix([sums[0]]) * change_of_basis )
  sums[0] = mults(Fraction(1,math.factorial(n+m+k-2)), sums[0])
if m+k-2 < 0 :
  sums[1] = [0] * len(PARTITIONS[num_points])
else:  
  sums[1] = simplify( Matrix([sums[1]]) * change_of_basis )
  sums[1] = mults(Fraction(-1,math.factorial(m+k-2)), sums[1])
if n+k-2 < 0 :
  sums[2] = [0] * len(PARTITIONS[num_points])
else:  
  sums[2] = simplify( Matrix([sums[2]]) * change_of_basis )
  sums[2] = mults(Fraction(-1,math.factorial(n+k-2)), sums[2])
if k-1 <= 0:
  sums[3] = [0] * len(PARTITIONS[num_points])
else:
  sums[3] = simplify( Matrix([sums[3]]) * change_of_basis )
  sums[3] = mults(k * (k-1), sums[3])

for s in sums:
  print(s)

output = add(sums[0],add(sums[1],add(sums[2],sums[3])))
print("n=",n,"m=",m,"k=",k)
print([x for x in zip(mults(Fraction(-1,2), output),PARTITIONS[num_points])])

