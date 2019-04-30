#!/usr/bin/python
# this thing is going to precompute a bunch of nak classes in the fixed point basis

from CE_alg import *
from partitions import *

nak = NAKS()


for i in range(1,7):
  for part in PARTITIONS[i]:
    partstr = "".join(str(x) for x in part) 
    print(partstr + ": " +str(nak.nak(partstr)))

nak.store()
