from numpy import *
from compcodim import *
from itertools import permutations

###
# Input is a list of two entries: [partition,entry_to_intersect_by_line]
# A partition is a nonincreasing list of integers, each integer decsribes the number of points on a moving line
# entry_to_intersect_by_line is an integer in the partition which describes which moving line to intersect with the fixed line from H
###
def intersect_with_H(left):
  
  components = []
  coefficients = []
  
  [comps,coefs] = first_degeneration(left)
  
  while comps:
    
    comp = comps.pop(0)
    coef = coefs.pop(0)
    
    components.append(comp)
    coefficients.append(coef)
      
    if len(comp) == 3:
      continue
    
    [temp_comps,temp_coefs] = second_degeneration(comp)
    
    del temp_comps[0]
    del temp_coefs[0]
    
    for a in temp_comps:
      i = temp_comps.index(a)
      if a in comps:
        j = comps.index(a)
        if coefs[j] == temp_coefs[i] * coef:
          del coefs[j]
          del comps[j]
        else:
          coefs[j] -= temp_coefs[i] * coef
      else:
        comps.append(a)
        coefs.append(-1 * temp_coefs[i] * coef)

  return [components,coefficients]

###
# Degeneration One
###

def first_degeneration(left):
  entry = left[1]
  part = list(left[0])
  
  entry_indices = [i for i,e in enumerate(part) if e == entry]
  
  no_ones = [x for x in part if x != 1]
  
  sigma = list(part)
  sigma.remove(entry)
  
  components = [[0,[entry],sigma]]
  coefficients = [entry]
  
  part[:] = [x - 1 for x in part]
  entry -= 1
  
  component_groups = []
  
  for k in range(1,len(no_ones)+1):
    group = []
    choices = combinations(range(len(no_ones)),k)
    
    for choice in choices:
      if any([i in choice for i in entry_indices]):
        temp = list(part)

        for i in choice:
          temp[i] = temp[i] - 1
        
        if all([temp[i] >= temp[i+1] for i in range(len(temp)-1)]):
          group.append([temp,k])
       
    component_groups.append(group)

  for lams in component_groups:
    taus = generate_taus(lams)
    
    sols = degeneration1_intersection_solutions(left, taus)
    matrx = generate_intersection_matrix(lams, taus)
    
    # Issue is in intersect_original_locus, of course. *sigh*
    
    coefs = determine_coefficients(matrx,sols)
    
    components += lams
    coefficients += coefs
    
  return [components,coefficients]


# left is the same format with different meaning [partition,an_entry_of_the_partition_to_meet]
def intersect_original_locus(left, tau):
  entry = left[1]
  part = list(left[0])
  
  tau_part = list(tau[0])
  tau_len = tau[1]
  
  perms = {x for x in permutations(part)}
  choices = len(perms)
  
  entries = [x for x,y in enumerate(part) if y == entry]
  
  part = [x-1 for x in part]
  entry -= 1
    
  diagram_0 = part
  diagram_1 = tau_part + [tau_len]
  diagram_2 = {(len(tau_part),entries[0])}
  
  diagram = [diagram_0,diagram_1,diagram_2]

  t = generate_tree(diagram)
  
  number = t.numLeaf * choices * len(entries)

  return number

def degeneration1_intersection_solutions(left, taus):
  return [intersect_original_locus(left, tau) for tau in taus]


###
# Degeneration Two
###

def second_degeneration(left):
  special_pt = left[1]
  part = list(left[0])
  
  no_zeros = [entry for entry in part if entry != 0]
  
  components = [[part,special_pt]]
  coefficients = [1]
  
  component_groups = []
  
  for k in range(1,len(no_zeros)+1):
    group = []
    choices = combinations(range(len(no_zeros)),k)
    
    for choice in choices:
      temp = list(part)
      for i in choice:
        temp[i] = temp[i] - 1
        
      if all([temp[i] >= temp[i+1] for i in range(len(temp)-1)]):
        group.append([temp,special_pt+k])
        
    component_groups.append(group)
  
  for lams in component_groups:
    taus = generate_taus(lams)
    
    sols = degeneration2_intersection_solutions(left, taus)
    matrx = generate_intersection_matrix(lams, taus)
    
    coefs = determine_coefficients(matrx,sols)
    
    components += lams
    coefficients += coefs
    
  return [components,coefficients]
  
# this computes a single intersection number from the left hand side with a tau class, both given as [partition,length_at_special_pt]
def intersect_left(left, tau):
  special_pt = tau[1]
  left_spc_pt = left[1]

  lam_part = list(left[0])
  tau_part = list(tau[0])
  
  length = special_pt - left_spc_pt
  
  tau_part.append(length)
  
  left_component = [[],[],[x+1 for x in lam_part]]
  corr_class = [[1]* len(lam_part),tau_part,[]]
  
  number = compute_intersection_number(left_component, corr_class)
  
  return number
  
# this takes a list of taus and computes the intersection of this component with the left hand side
def degeneration2_intersection_solutions(left, taus):
  return [intersect_left(left,tau) for tau in taus]

###
# General Purpose
###

def generate_taus(lams):
  taus = []
  
  for lam in lams:
    lam_part = list(lam[0])
    lam_spc_pt = lam[1]
    
    no_zeros = [entry for entry in lam_part if entry != 0]
    
    tau = [conjugate(no_zeros),lam_spc_pt]
    taus.append(tau)
  
  return taus

# returns conjugate of a partition
# taken from: https://www.ics.uci.edu/~eppstein/PADS/IntegerPartitions.py
def conjugate(p):
    result = []
    j = len(p)
    if j <= 0:
        return result
    while True:
        result.append(j)
        while len(result) >= p[j-1]:
            j -= 1
            if j == 0:
                return result

# lam and tau are two lists [partition,length_at_special_pt]
def intersect_components(lam, tau):
  special_pt = tau[1]
  lam_spc_pt = lam[1]
  
  if special_pt != lam_spc_pt:
    return 0
    
  lam_part = list(lam[0])
  tau_part = list(tau[0])
  
  component = [[],[],[x+1 for x in lam_part]]
  corr_class = [[1]* len(lam_part),tau_part,[]]
  
  number = compute_intersection_number(component, corr_class)
    
  return number

# this takes a list of lambdas and a list of taus and computes the matrix to solve
def generate_intersection_matrix(lams, taus):
  matrx = []
  for tau in taus:
    temp = [intersect_components(lam,tau) for lam in lams]
    matrx.append(temp)
  return matrx
  
# this solves the system Ax = b for A from generate intersection matrix and b from generate_intersection_solutions
def determine_coefficients(matrx, sols):
  output = linalg.solve(array(matrx),array(sols))
  
  coefs = [int(round(out)) for out in output]
  
  return coefs



