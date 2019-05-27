#
# a partition is a nonincreasing list of positive integers
# lists of partitions are always assumed to be in REVERSE
# LEXICOGRAPHIC ORDER
#

# We use these a lot, so might as well hard code them. 
_PARTITIONS = {
  1: [[1]],
  2: [[2],[1,1]],
  3: [[3],[2,1],[1,1,1]],
  4: [[4],[3,1],[2,2],[2,1,1],[1,1,1,1]],
  5: [[5],[4,1],[3,2],[3,1,1],[2,2,1],[2,1,1,1],[1,1,1,1,1]],
  6: [[6],[5,1],[4,2],[4,1,1],[3,3],[3,2,1],[3,1,1,1],[2,2,2],[2,2,1,1],[2,1,1,1,1],[1,1,1,1,1,1]]
}

# A very efficient pythonic way of generating partitions
# not sure where I got this, but I def didn't come up with it
def _gen_partitions(n):
    """
    Generate all integer partitions of n
    in lexicographic order
    """
    if n == 0:
        yield []
        return
    for p in _gen_partitions(n - 1):
        yield p+[1]
        if len(p) == 1 or (len(p) > 1 and p[-1] < p[-2]):
            p[-1] += 1
            yield p

# PARTITION(n)
# Input: an integer n
# Output: all partitions of n in REVERSE LEXICOGRAPHIC ORDER
# We want to be able to handle a call for partitions of n > 6 just in case, but we want to use the hard coded ones for n <= 6. 
def PARTITION(n):
  if(n <= 0):
    return 0
  if(n <= 6):
    return _PARTITIONS[n]
  
  return [x for x in reversed([y for y in _gen_partitions(n)])]


# DEPRECATED
# def gen_partitions_1(inner):
#  partitions = []
#
#  for i in range(len(inner)):
#    if i == 0 or inner[i] != inner[i-1]:
#      partitions.append(inner[:i] + [inner[i]+1] + inner[i+1:])
# 
#  partitions.append(list(inner)+[1])
#
#  return partitions

# subpartition
# input: two partitions as lists of positive integers
# output: True if first partition is a subpartition of the second
def subpartition(part1,part2):
  if len(part1) > len(part2):
    return False
  
  return all([part1[i] <= part2[i] for i in range(len(part1))])
