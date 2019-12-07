import numpy as np
import cmath
import math

def get_factors(n):
	pairwise_factors = []
	for i in range(1, int(math.sqrt(n)) + 1):
		if n % i == 0:
			pairwise_factors.append((i, int(n/i)))
	return pairwise_factors

# find the prime factorization using sqrt(len(X)) property of primes
def find_primes(n):
	primes = set()
	while n % 2 == 0:
		primes.add(2)
		n /= 2
	for i in range(3, int(math.sqrt(n))+1, 2):
		while n % i == 0:
			primes.add(int(i))
			n /= i
	if n > 2:
		primes.add(int(n))
	return primes

def find_relative_primes(pairwise_factors):
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
	N = len(X)
	pairwise_factors = get_factors(N)
	N_1, N_2 = find_relative_primes(pairwise_factors)

	# reindexing for our new arrays of N_1, N_2 
	input_index = np.zeros((N_1, N_2))
	i = 0
	j = 0
	for n in range(N):	
		input_index[i % N_1][j % N_2] = n
		i += 1
		j += 1
	output_index = np.zeros((N_1, N_2))
	for k_1 in range(N_1):
		for k_2 in range(N_2):
			output_index[k_1][k_2] = (k_1 * N_2 + k_2 * N_1) % N

	# 2D DFT
	output = np.zeros((N_1, N_2), dtype=complex)
	for k_1 in range(N_1):
		for k_2 in range(N_2):
			# very inefficient DFT
			output[k_1][k_2] = sum([sum([X[int(input_index[i][j])] * cmath.exp(-2.0j * np.pi * i * k_1 / N_1) for i in range(N_1)]) * cmath.exp(-2.0j * np.pi * j * k_2 / N_2) for j in range(N_2)])

	result = np.zeros(N, dtype=complex)
	for k_1 in range(N_1):
		for k_2 in range(N_2):
			result[int(output_index[k_1][k_2])] = output[k_1][k_2]
	return result

X = np.tile(np.array([0,1,2,3,4], dtype=complex), 6)
transformed = prime_fft(X)
# compare to numpy 1d fft
np_result = np.fft.fft(X)

for i in range(len(transformed)):
	print(transformed[i], np_result[i])

