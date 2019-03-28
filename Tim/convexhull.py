#!/usr/bin/python
import numpy as np
from scipy.linalg import null_space

o2_2 = [1,1,0,1,0]
o2_11 = [2,0,0,2,1]
o3_2 = [1,2,1,2,0]
o3_11 = [2,1,0,4,7]
t1_2 = [5,5,1,7,5]
t1_11 = [7,3,0,9,12]
e25_2 = [35,25,3,43,40]
e25_11 = [40,20,1,47,54]
e35_2 = [35,30,6,52,50]
e35_11 = [40,25,3,58,69]
e513_2 = [247,156,15,294,294]
e513_11 = [260,143,10,304,329]
e813_2 = [247,195,36,369,384]
e813_11 = [260,182,28,385,434]
od_2 = [0,0,1,0,0]
od_11 = [0,0,1,0,2]
td_2 = [0,0,2,0,1]
td_11 = [0,0,2,0,3]

points = [o2_2, o2_11, o3_2, o3_11, t1_2, t1_11, e25_2, e25_11, e35_2, e35_11, e513_2, e513_11, e813_2, e813_11, od_2, od_11, td_2, td_11]

for i,point in enumerate(points):
  print(str(i)+": "+str(point))

extremal = []

for i in range(len(points)):
  for j in range(i+1,len(points)):
    for k in range(j+1,len(points)):
      for l in range(k+1,len(points)):
        A = np.array([points[i],points[j],points[k],points[l]])
        normal = null_space(A)
        if normal.shape[1] != 1:
          continue
        prod = [np.asscalar(np.dot(x,normal)) for x in points]
        if all([x >= 0 for x in prod]) or all([x <= 0 for x in prod]):
          extremal.append([i,j,k,l])

print(extremal)

output = set()
for x in extremal:
  for y in x:
    output.add(y)

print(output)
