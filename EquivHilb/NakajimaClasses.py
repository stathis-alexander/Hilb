#!/usr/bin/python
#
# This script computes the Nakajima classes in the Hilbert scheme of 
# n points on the affine plane in the fixed point basis. 
#
#

from CE_alg import *
from partitions import PARTITION

nstart = 7
nstop = 7

# Load NAKS class from CE_Alg.
nak = NAKS()

# This will print all the Nakajima classes in the fixed point basis
# for n=nstart to n=nstop. 
for i in range(nstart,nstop+1):
  for part in PARTITION(i):
    partstr = passtr(part)
    print(partstr + ": " +str(nak.nak(partstr)))

# NAKS "descructor"
nak.store()

