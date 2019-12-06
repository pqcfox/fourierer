import numpy as np
import cmath

def get_factors(n):
	pairwise_factors = []
	for i in range(int(cmath.sqrt(n))):
		if n % i == 0:
			factors.append((i, n/i))
	return factors

# find the prime factorization using sqrt(len(X)) property of primes
def find_primes(n):
	primes = set()
	while n % 2 == 0:
		primes.add(2)
		n /= 2
	for i in range(3, int(cmath.sqrt(n))+1, 2):
		while n % i == 0:
			primes.add(i)
			n /= i
	if n > 2:
		primes.add(n)
	return primes

def find_relative_primes(pairwise_factors)
	# start search from end to get optimal factorization
	for i in range(len(pairwise_factors) - 1, 0, -1):
		pair = pairwise_factors[i]
		set_primes_1 = find_primes(pair[0])
		set_primes_2 = find_primes(pair[1])
		if len(set_primes_1.intersection(set_primes_2)) == 0:
			return pair
	raise ValueError('length of X cannot be used for prime factor FFT')

"""
Solves the 1D FFT using the Good-Thomas algorithm.
Input:
	X = 1D numpy array of length N, where N can be factored into at least two relatively prime factors
Output:
	result = transformed 1D numpy array
"""
def prime_fft(X):
	pairwise_factors = get_factors(len(X))
	N_1, N_2 = find_relative_primes(pairwise_factors)

	# reindexing for our new arrays of N_1, N_2 
	# need to solve for p,q,r,s using Chinese Remainder Theorem
	

	# now we can compute a 2D DFT, which is a series of different 1D DFTs


# compare to numpy 1d fft