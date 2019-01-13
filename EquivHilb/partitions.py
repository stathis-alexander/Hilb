#
# a partition is a nonincreasing list of positive integers
#

# a generating function for all partitions of n
def gen_partitions(n):
    """
    Generate all integer partitions of n
    in lexicographic order
    """
    if n == 0:
        yield []
        return
    for p in gen_partitions(n - 1):
        yield p+[1]
        if len(p) == 1 or (len(p) > 1 and p[-1] < p[-2]):
            p[-1] += 1
            yield p


def gen_partitions_1(inner):
  partitions = []

  for i in range(len(inner)):
    if i == 0 or inner[i] != inner[i-1]:
      partitions.append(inner[:i] + [inner[i]+1] + inner[i+1:])
 
  partitions.append(list(inner)+[1])

  return partitions





