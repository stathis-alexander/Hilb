#!/usr/bin/python
#
# This script computes the Nakajima classes in the Hilbert scheme of 
# n points on the affine plane in the fixed point basis. 
#
#

from CE_alg import *
from partitions import *

nstart = 1
nstop = 6

# Load NAKS class from CE_Alg.
nak = NAKS()

# This will print all the Nakajima classes in the fixed point basis
# for n=1 to n=6. 
for i in range(nstart,nstop+1):
  for part in PARTITIONS[i]:
    partstr = passtr(part)
    print(partstr + ": " +str(nak.nak(partstr)))

# NAKS "descructor"
nak.store()

