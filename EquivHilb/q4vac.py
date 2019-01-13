#!/usr/bin/python
from CE_alg import *
from partitions import *

import math

#Compute q4(vacuum) ---

parts4 = [[4],[3,1],[2,2],[2,1,1],[1,1,1,1]]

for par in parts4:
  print(par,":")
  print(q_struct_const(4,[],par))

