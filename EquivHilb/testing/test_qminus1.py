#!/usr/bin/python
# This is a template file for a running code using the CE_alg library.

from CE_alg import *
from partitions import *

lam = [1,1]
mu = [3,2]

#print(coker(lam,mu))
#print(ker(lam,mu))
#print(tangent(lam))
#print(tangent(mu))
#print(qminus1_struct_const([1,1],[2,1]))

print(q(-1,[(1,mu)],PARTITIONS[4]))
