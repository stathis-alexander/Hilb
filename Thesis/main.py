#!/usr/bin/python

from compcodim import *
from pieri import *

# change me, the first list is the partition to use, the second is a number in the partition and describes which line to intersect with the fixed line from H
data = [[2,2,1],2]

# don't touch me
output = intersect_with_H(data)

# don't touch me, I'm for printing
for i in range(len(output[0])):
  print(output[1][i],output[0][i])
  
#output is the coefficient, followed by a list with a partition describing the locus, these partitions describe moving lines with the entry + 1 many points on them, and an integer which is the length of schemes in the locus supported at a general point
