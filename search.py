from itertools import chain, combinations
import numpy

def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))

def phi_index(p, i):
  """ phi(p^i), with boundary conditions for 0 """
  if i == 0:
    return 1
  return (p-1)*pow(p, i-1)

class SoMat(object):
  """ a matrix with possible duplications of rows """
  
  def __init__(self, mat, row_mults):
    self.mat = mat
    self.row_mults = row_mults

  def __str__(self):
    return('\n'.join(['\t'.join([str(cell) for cell in row]) + '\t (x{})'.format(self.row_mults[ind])
                      for (ind, row) in enumerate(self.mat)]))
    
  def __mul__(self, other):
    n1 = len(self.mat) # rows
    n2 = len(other.mat) # rows
    m1 = len(self.mat[0]) # cols (equals rows, but separate in case we want to generalize)
    m2 = len(other.mat[0])
    new_mat = []
    new_row_mults = []
    for i in range(n1):
      for k in range(n2):
        new_mat.append([])
        new_row_mults.append(self.row_mults[i] * other.row_mults[k])
        for j in range(m1):
          for l in range(m2):
            new_mat[-1].append(self.mat[i][j] * other.mat[k][l])
    return SoMat(new_mat, new_row_mults)

  @classmethod
  def _single_mat_func(cls, p, i, j):
    if j >= i + 2:
      return 0
    if j == i + 1:
      return -pow(p, i)
    else:
      return phi_index(p, j)
    
  @classmethod
  def SingleDimension(cls, p, N):
    new_mat = []
    new_row_mults = []
    for i in range(N+1):
      new_mat.append([])
      new_row_mults.append(phi_index(p, N-i))
      for j in range(N+1):
        new_mat[-1].append(cls._single_mat_func(p, i, j))
    return cls(new_mat, new_row_mults)

  @classmethod
  def FromTuples(cls, ft):
    """ example ft: [(2,2), (3,5)] means we want 2^2 3^5"""
    cur_mat = cls.SingleDimension(2, 0)
    for (a, b) in ft:
      new_mat = cls.SingleDimension(a, b)
      cur_mat *= new_mat
    return cur_mat

def deflate(vect, mults):
  """ takes a vector with multiplicities and makes a bigger one with no multiplicities"""
  ourdict = {}
  for (i,v) in enumerate(vect):
    if not (v in ourdict):
      ourdict[v] = mults[i]
    else:
      ourdict[v] += mults[i]
  return ourdict

def mat_deflate(mat, rows, mults):
  """
  given a matrix with *column* multiplicities, pick particular rows, add them, and then deflate.

  for speed of computation we transpose everything, so in the original problem we sum columns
  but for our problem we sum rows

  [mults] would then be column multiplicities (but in our problem they are row multiplicities)
  
  In our problem, this gets the column vector sum to a format that can be compared to another
  to see if there's a clash.
  """
  col_count = len(mat)
  vect = numpy.array([0 for i in range(col_count)])
  for a in rows:
    vect += numpy.array(mat[a])
  dvect = deflate(vect, mults)
  return tuple(sorted(list(dvect.items()), key = lambda x: x[0]))

def create_dvects(somat):
  rows = len(somat.mat)
  cols = len(somat.mat[0])
  zmat = tuple(zip(*somat.mat))
  result_dict = {}
  for A in powerset(range(1, cols)):
    dvect = mat_deflate(zmat, A, somat.row_mults)
    result_dict[tuple(A)] = dvect
  return result_dict

def search_clash(somat, verbose=False):
  """
  Given a matrix, see if any pair of subsets of columns have sum equal (as sets)
  """
  cols = len(somat.mat[0])
  dvects = create_dvects(somat)
  keys = [tuple(k) for k in powerset(range(1, cols))]
  for i in range(len(keys) - 1):
    for j in range(i + 1, len(keys)):
      A = keys[i]
      B = keys[j]
      dvect = dvects[A]
      dvect2 = dvects[B]
      if dvect == dvect2:
        if A == B and not verbose:
          continue
        print("clash found: {} and {}".format(A, B))
        print("  " + str(dvect))
        print("  " + str(dvect2))

def search_clash_jit(somat):
  # just-in-time version
  rows = len(somat.mat)
  cols = len(somat.mat[0])
  zmat = tuple(zip(*somat.mat))
  result_dict = {}
  reverse_dict = {}
  for A in powerset(range(1, cols)):
    dvect = mat_deflate(zmat, A, somat.row_mults)
    result_dict[tuple(A)] = dvect
    if dvect in reverse_dict:
      B = reverse_dict[dvect]
      print(f"Clash found for {A} and {B}: {dvect}")
    reverse_dict[dvect] = tuple(A)
  
def prime_factors(n):
  i = 2
  factors = {}
  while i * i <= n:
    if n % i:
      i += 1
    else:
      n //= i
      if i in factors:
        factors[i] = factors[i]+1
      else:
        factors[i] = 1
  if n > 1:
    if n in factors:
      factors[n] = factors[n] + 1
    factors[n] = 1
  return factors

def factorization_tuple(n):
  factors = prime_factors(n)
  return sorted(list(factors.items()), key = lambda x: x[0])

def mega_search(initial_N=1, upper_bound=2000, row_bound=24):
  """
  Searching N from [initial_N] to [upper_bound], skipping any N where the number of rows
  is too large.

  In particular, for N=p^2q^2, p < 11, we only need to search up to 7^2*(7^3)^2
  """
  for N in range(initial_N, upper_bound):
    ft = factorization_tuple(N)
    somat = SoMat.FromTuples(ft)
    rows = len(somat.mat)
    print(f"N: {N} number of rows: {rows}")
    if rows >= row_bound:
      print("  skipped!")
      continue
    search_clash_jit(somat)
    N += 1

# For the Paper #

def p2q2_search(p_bound, q_bound):
  f = open("primes.txt", "r")
  primes = [int(line) for line in f]
  for i in range(len(primes)-1):
    if primes[i] >= p_bound:
      break
    for j in range(i+1, len(primes)):
      if primes[j] >= q_bound:
        break
      p = primes[i]
      q = primes[j]
      N = p*p*q*q
      somat = SoMat.FromTuples([(p,2), (q,2)])
      rows = len(somat.mat)
      assert rows == 9
      search_clash(somat)
