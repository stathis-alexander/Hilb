#!/usr/bin/python
from CE_alg import *
from partitions import *

import math

#Check commutation of q3,q1 ---

q3_111 = q_struct_const(3,[],[1,1,1])
q3_21 = q_struct_const(3,[],[2,1])
q3_3 = q_struct_const(3,[],[3])

print(q3_111)
print(q3_21)
print(q3_3)

print("--")

q1_21 = q1([2,1])
q1_3 = q1([3])

print(q1([1,1,1]))
print(q1([2,1]))
print(q1([3]))
print("--")
print(simplify(q3_3 * q1_3[1][0] + q1_21[0][0] * q3_21))
print(q_struct_const(3,[1],[3,1]))
print("--")
print(simplify(q3_3 * q1_3[1][0] + q1_21[0][0] * q3_21 - q_struct_const(3,[1],[3,1])) == 0)

