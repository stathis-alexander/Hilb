#!/usr/bin/python
# compute the chern characters in terms of nak basis for n points (up to 5) using CE's algorithm at fixed points

from CE_alg import *
from partitions import *

num_points = 4

# 2 point classes
nak11 = [1,1]
nak2 = [-2*V,-2*U]
fp2 = [[2],[1,1]]

chob = Matrix([nak11,nak2])
matrix_file = "rels2"

# 3 point classes
nak111 = [1,1,1]
nak21 = [-6*V,-2*(U+V),-6*U]
nak3 = [6*V**2,3*U*V,6*U**2]
fp3 = [[3],[2,1],[1,1,1]]

# 4 point classes
nak1111 = [1,1,1,1,1]
nak211 =  [-24*V,-4*(U+3*V),-8*(U+V),-4*(3*U+V),-24*U]
nak31 =   [24*V**2,6*V*(U+V),12*U*V,6*U*(U+V),24*U**2]
nak22 =   [24*V**2,8*U*V,8*(U**2-U*V+V**2),8*U*V,24*U**2]
nak4 =    [-24*V**3,-8*U*V**2,-4*U*V*(U+V),-8*U**2*V,-24*U**3]
fp4 = [[4],[3,1],[2,2],[2,1,1],[1,1,1,1]]

# 5 point classes
nak11111 = [1,1,1,1,1,1,1]
nak2111 =  [-120*V,-12*U-72*V,-24*U-48*V,-36*U-36*V,-48*U-24*V,-72*U-12*V,-120*U]
nak221 =   [120*V**2,24*V*(U+V),8*U**2+8*U*V+24*V**2,40*U*V,24*U**2+8*U*V+8*V**2,24*U*(U+V),120*U**2]
nak311 =  [120*V**2,V*(18*U+48*V),V*(36*U+12*V),12*U**2+24*U*V+12*V**2,U*(12*U+36*V),U*(48*U+18*V),120*U**2]
nak32 =   [-120*V**3,-30*U*V**2,-12*V*(U**2+V**2),-12*U*V*(U+V),-12*U*(U**2+V**2),-30*U**2*V,-120*U**3]
nak41 =   [-120*V**3,-24*V**2*(U+V),-U*V*(4*U+28*V),-16*U*V*(U+V),-U*V*(28*U+4*V),-24*U**2*(U+V),-120*U**3]
nak5 =    [120*V**4,30*U*V**3,10*U*V**2*(U+V),20*U**2*V**2,10*U**2*V*(U+V),30*U**3*V,120*U**4] 
fp5 = [[5],[4,1],[3,2],[3,1,1],[2,2,1],[2,1,1,1],[1,1,1,1,1]]

if num_points == 3:
  matrix_file = "rels3"
  chob = Matrix([nak111,nak21,nak3])
elif num_points == 4:
  matrix_file = "rels4"
  chob = Matrix([nak1111,nak211,nak31,nak22,nak4])
elif num_points == 5:
  matrix_file = "rels5"
  chob = Matrix([nak11111,nak2111,nak221,nak311,nak32,nak41,nak5])

change_of_basis = load_or_compute(chob,matrix_file)
nak = NAKS()

print([(y,z) for y,z in zip(nak311, [x[0] for x in nak.nak("311")])] )


#p1 = mult(nak211,nak211)
#p2 = mult(nak2,nak11)
#p3 = mult(nak21,nak21)

#p2 = [simplify(x[0]) for x in qlist(2,zip(p2,fp2),fp4)]
#p3 = [simplify(x[0]) for x in q1list(zip(p3,fp3))]

#print(simplify( Matrix([p1]) * change_of_basis ))
#print(simplify( Matrix([p2]) * change_of_basis ))
#print(simplify( Matrix([p3]) * change_of_basis ))


