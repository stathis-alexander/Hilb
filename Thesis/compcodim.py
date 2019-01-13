from itertools import permutations, combinations

####
## A Simple Tree Data Structure
## 

class Tree():
   def __init__(self, parent, horiz, vert):
      self.parent = parent
      self.children = []
      self.horiz = horiz
      self.vert = vert
      self.numLeaf = 0

# generates mixed partitions for n
def mixed_partition(n):
   partition = [0 for i in range(n+1)]
   k = 1
   y = n - 1
   while k != 0:
      x = partition[k - 1] + 1
      k -= 1
      while 2*x <= y:
         partition[k] = x
         y -= x
         k += 1
      l = k + 1
      while x <= y:
         partition[k] = x
         partition[l] = y
         if k + 2 <= 3:
            output = list(reversed(partition[:k+2]))
            output = output + [0] * (1-k)
            yield output
         x += 1
         y -= 1
      partition[k] = x + y
      y = x + y - 1
      if k+1 <= 3:
         output = list(reversed(partition[:k+1]))
         output = output + [0] * (2-k)
         yield output

# Code from http://jeromekelleher.net/partitions.php cite later
# modified above to only emit partitions of length m if specified
def partition(n, m = None):
   partition = [0 for i in range(n+1)]
   k = 1
   y = n - 1
   while k != 0:
      x = partition[k - 1] + 1
      k -= 1
      while 2*x <= y:
         partition[k] = x
         y -= x
         k += 1
      l = k + 1
      while x <= y:
         partition[k] = x
         partition[l] = y
         if not m or (k + 2) == m:
            yield list(reversed(partition[:k + 2]))
         x += 1
         y -= 1
      partition[k] = x + y
      y = x + y - 1
      if not m or k+1 == m:
         yield list(reversed(partition[:k + 1]))

# A diagram is a list containing exactly three members:
# 0: A nonincreasing sequence of positive integers representing incidence conditions on horizontal lines
# 1: A nonincreasing sequence of positive integers representing incidence conditions on vertical lines
# 2: a set of (i,j) pairs corresponding to when the ith horizontal line must meet the jth vertical line
# it's assumed that the sum of the entries in the vertical list is equal the sum of the entries in the horizontal list
def generate_tree( diagram ):
   deg_tuple1 = diagram[1]
   deg_tuple2 = diagram[0]
   point_set = diagram[2]
   
   t = Tree( None, deg_tuple1, deg_tuple2 )
   generate_tree_worker( t, point_set, 0 )
   return t
   
def generate_tree_worker( t, point_set, height ):
   
   if all( i == 0 for i in t.horiz):
      if all( degree == 0 for degree in t.vert ):
         return 1
      else:
         return 0

   index_list = [i for i, e in enumerate(t.vert) if e != 0]
   line = next((i for i,e in enumerate(t.horiz) if e != 0), None)
   required_points = {(x,y) for (x,y) in point_set if x == line}
   choices = t.horiz[line]
   if choices > len( index_list ):
      return 0
   
   number_viable = 0
   child_deg_tuple_zero = list(t.horiz)
   child_deg_tuple_zero[line] = 0
   subsequences = combinations(index_list, choices)
   
   for subsequence in subsequences:
      if any(y not in subsequence for (x,y) in required_points):
         continue

      new_list = list( t.vert )
      for i in subsequence:
         new_list[i] -= 1
      
      child = Tree(t, child_deg_tuple_zero, new_list)
      t.children.append(child)
      number_viable += generate_tree_worker(child, point_set, height+1)
   
   t.numLeaf = number_viable;
   return number_viable
# this function takes as input two mixed partitions and returns all the diagrams to compute the interseciton number for those corresponding basis elements
def generate_diagrams(sigma, tau):
   sigma_fixed = sigma[0] + sigma[1]
   tau_fixed = tau[0] + tau[1]

   list_of_diagrams = []

   point_set_sigma = set( (len(sigma_fixed)+i,i) for i in range(0,len(tau[0])) )
   point_set_tau = set( (i,len(tau_fixed)+i) for i in range(0, len(sigma[0])) )

   point_set = point_set_sigma | point_set_tau

   all_permutations_sigma = {x for x in permutations(sigma[2])}
   all_permutations_tau = {x for x in permutations(tau[2])}

   for c_tau in all_permutations_tau:
      for c_sigma in all_permutations_sigma:
         final_list_sigma = sigma_fixed + list(c_sigma)
         final_list_tau = tau_fixed + list(c_tau)
         list_of_diagrams.append([final_list_tau,final_list_sigma,point_set])

   return list_of_diagrams
# this takes as input two mixed partitions and computes the intersection number. it is assumed that the length of the first element of the first mixed partition is the length of the third element of the second, and vice versa
def compute_intersection_number( mixed_partition_1, mixed_partition_2 ):
  list_of_diagrams = generate_diagrams( mixed_partition_1, mixed_partition_2 )
  number = 0
  for diagram in list_of_diagrams:
    t = generate_tree( diagram )
    number += t.numLeaf
  return number

def print_diagram(diagram):
   horiz_labels = [str(x) for x in diagram[0]]
   vert_labels = [str(x) for x in diagram[1]]
   point_set = diagram[2]
   col_size = max([len(x) for x in horiz_labels] + [len(x) for x in vert_labels])
   print(" ".ljust(col_size), " ".join([x.ljust(col_size) for x in vert_labels]))
   for i in range(0,len(horiz_labels)):
      string = ""
      hold = [x for (x,y) in point_set if y == i]
      if hold:
         if hold[0] == 0:
            string = "*"
         else:
            string = "*".rjust(hold[0]*(col_size+1)  + 1)
      print (horiz_labels[i].ljust(col_size), string) 

def print_diagrams( list_of_diagrams ):
   for diagram in list_of_diagrams:
      print_diagram( diagram )

def print_dfs(tree):
   discovered = []
   stack = []
   stack.append((tree,0))
   while stack:
      nodepair = stack.pop()
      node = nodepair[0]
      depth = nodepair[1]
      if node not in discovered:
         print (str(depth) + ":: " + str(node.horiz) + "::" + str(node.vert) + "::" + str(node.numLeaf))
         discovered.append(node)
         for child in node.children:
            stack.append((child,depth + 1))

def generate_basis(N, codim):
   list_of_basis = []
   fixed_range = codim // 2 + 1
   for mixed in mixed_partition(N):
      perms = {x for x in permutations(mixed)}
      for perm in perms:
         for fixed in range(0, fixed_range):
            moving = fixed + N - codim
            on_lines = N - moving - fixed
            if perm[0] < fixed or perm[2] < moving or (not fixed and perm[0]) or (not moving and perm[2]):
               continue
            for first in partition(perm[0], fixed):
               if len(first) == 1 and first[0] == 0:
                  first = []
               for second in partition(perm[1]):
                  if len(second) == 1 and second[0] == 0:
                     second = []
                  for third in partition(perm[2], moving):
                     if len(third) == 1 and third[0] == 0:
                        third = []
                     list_of_basis.append([first, second, third])
   return list_of_basis

def order_basis(list_of_basis):
   print(list_of_basis.sort())
   
def key_fcn(basis):
   new_base = [] 
   for entry in basis:
      a = sum(entry)
      new_base.append([a] + entry)
   return new_base

def compute_intersection_matrix(N, d):
   basis = generate_basis(N,d)
   opposite_basis = [list(reversed(x)) for x in basis]
   
   basis = sorted(basis, key=key_fcn)
   opposite_basis = sorted(opposite_basis, key=key_fcn)
   
   print(" ----- Basis Codim", d ,"(rows) ----- ")
   for (i,b) in enumerate(basis):
      print(i, ":: ", b, sep = "")
   print(" ----- Basis Codim", 2*N-d ,"(cols) ----- ")
   for (i,b) in enumerate(opposite_basis):
      print(i, ":: ", b, sep="")

   matrix = []

   for b in basis:
      row = []
      for o in opposite_basis:
         if len(b[0]) == len(o[2]) and len(b[2]) == len(o[0]):
            row.append(compute_intersection_number(b,o))
         else:
            row.append(0)
      matrix.append(row)

   return matrix
