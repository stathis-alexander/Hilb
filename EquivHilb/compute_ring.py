#!/usr/bin/python
# This is a template file for a running code using the CE_alg library.

from CE_alg import *
from partitions import *

nak11 = [(2,[2]),(2,[1,1])]
nak2 = [(-2*V,[2]),(-2*U,[1,1])]

def product(class1,class2):
 return [(simplify(x[0]*y[0]),x[1]) for x,y in zip(class1,class2)]

print(product(nak2,nak2))
