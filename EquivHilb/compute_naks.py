#!/usr/bin/python
# This is a template file for a running code using the CE_alg library.

from CE_alg import *
from partitions import *

print("Compute Nakajima classes in terms of fix classes.")
print("Hilb(2):")

partitions = [[2],[1,1]]

nak11 = q1([1])
nak2 = [(q_struct_const(2,[],lam),lam) for lam in partitions]
print("nak(1,1) =",nak11)
print("nak(2) =",nak2)

print("Hilb(3):")

partitions = [[3],[2,1],[1,1,1]]

nak111 = q1list(nak11)
nak21 = q1list(nak2)
nak3 = [(q_struct_const(3,[],lam),lam) for lam in partitions]

print("nak(1,1,1) =",nak111)
print("nak(2,1) =",nak21)
print("nak(3) =", nak3)

print("Hilb(4):")
partitions = [[4],[3,1],[2,2],[2,1,1],[1,1,1,1]]

nak1111 = q1list(nak111)
nak31 = q1list(nak3)
nak22 = qlist(2,nak2,partitions)
nak211 = q1list(nak21)
nak4 = [(q_struct_const(4,[],lam),lam) for lam in partitions]

print("nak(1,1,1,1) =",nak1111)
print("nak(3,1) =",nak31)
print("nak(2,2) =",nak22)
print("nak(2,1,1) =",nak211)
print("nak(4) =",nak4)
