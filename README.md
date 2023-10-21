# README

This repository is for python code that explores ICGs (integral circulant graphs) as presented in the paper "On Isospectral Integral Circulant Graphs". The discussion here assumes the context given in the paper.

`primes.txt` was copied from [The first 1000 and 10000 primes](https://www.di-mgt.com.au/primes1000.txt). It contains the first 1000 primes up to 7919.

## Examples of Usage

A **clash** (not defined in the paper, but is implied by the definitions involving "cospectral pair") is when a matrix contains $2$ column sets that sum to vectors which are equal under permutation.

First, we create a matrix with 2 columns that are unequal but are equal under permutation, to make sure that our `search_clash` function correctly finds the clash.

```
In [2]: import search                                                           

In [7]: test_clash_mat = search.SoMat([[1, 0, 3], [4, 1, 1], [5, 3, 0]], [2,1,2]
   ...: )                                                                       

In [8]: search.search_clash(test_clash_mat)                                     
clash found: (1,) and (2,)
  ((0, 2), (1, 1), (3, 2))
  ((0, 2), (1, 1), (3, 2))

In [9]: print(test_clash_mat)                                                   
1	0	3	 (x2)
4	1	1	 (x1)
5	3	0	 (x2)
```

## Mega Search

The `mega_search` is a general search for $N$ in a specified range with an upper bound on rows (we ran out of memory when the number of rows got too high). Recall that the number of rows is $\tau(N)$, the number of divisors of $N$.

We ran `mega_search` up through $30000$ with an upper bound of $24$ rows and found no counterexamples. So this would skip $N$ of form such as $p^3q^2r$. This can be replicated as

```
import search
search.mega_search(initial_N=1, upper_bound=30000, row_bound=24)
```

## Proof of Result Stated in Paper

The paper claims that when $p < 11$, $q < p^3$, there are no clashes when $N = p^2q^2$.

We search through all $N=p^2q^2$ with $p < 11$, $p < q$, $q < 7000$ in `primes.txt`. As 7000 is much bigger than $7^3$, this search is more than sufficient to prove our claim. It can be replicated as

```
import search                                                           
search.p2q2_search(11, 7000) 
```
