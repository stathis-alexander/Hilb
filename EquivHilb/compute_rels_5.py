#!/usr/bin/python
# This is a template file for a running code using the CE_alg library.

from CE_alg import *
from partitions import *

matrix_file = "rels5"

nak11111 = [1,1,1,1,1,1,1]
nak2111 =  [-120*V,-12*U-72*V,-24*U-48*V,-36*U-36*V,-48*U-24*V,-72*U-12*V,-120*U]
nak221 =   [120*V**2,24*V*(U+V),8*U**2+8*U*V+24*V**2,40*U*V,24*U**2+8*U*V+8*V**2,24*U*(U+V),120*U**2]
nak311 =  [120*V**3,V*(18*U+48*V),V*(36*U+12*V),12*U**2+24*U*V+12*V**2,U*(12*U+36*V),U*(48*U+18*V),120*U**2]
nak32 =   [-120*V**3,-30*U*V**2,-12*V*(U**2+V**2),-12*U*V*(U+V),-12*U*(U**2+V**2),-30*U**2*V,-120*U**3]
nak41 =   [-120*V**3,-24*V**2*(U+V),-U*V*(4*U+28*V),-16*U*V*(U+V),-U*V*(28*U+4*V),-24*U**2*(U+V),-120*U**3]
nak5 =    [120*V**4,30*U*V**3,10*U*V**2*(U+V),20*U**2*V**2,10*U**2*V*(U+V),30*U**3*V,120*U**4] 

fps = [[5],[4,1],[3,2],[3,1,1],[2,2,1],[2,1,1,1],[1,1,1,1,1]]

chob = Matrix([nak11111,nak2111,nak221,nak311,nak32,nak41,nak5])
change_of_basis = load_or_compute(chob,matrix_file)

chs = compute_ch(fps)

for i in range(len(chs)):
  print("ch"+str(i+1)+":")
  print(simplify( Matrix([chs[i]]) * change_of_basis ))
