import random
import sys
from statsmodels.sandbox.stats.runs import runstest_1samp
import numpy as np
import math

####################################### Pseudo Random Number generatiors ##########################################

'''
	Linear Congruential PRNG
'''
class LinearCongruential(object):
	def __init__(self, m, a, c, seed):
		self.a = a
		self.m = m
		self.c = c
		self.seed = seed

	def generate(self):
		x = self.seed
		while True:
			yield x
			x = (self.a * x + self.c) % self.m


'''
	Blum Blum Shub PRNG
'''
class BlumBlumShub(object):

	'''Initialize the seed value'''
	def __init__(self, seed, n = 70891 * 85257, random = False):
		self.seed = seed		
		if random:
			self.n = self.generateN()

	'''
		Returns a list of prime numbers from [2, n] using
		Sieve of Eratosthenes algorithm
	'''
	def get_primes(self, n):
		primes = []
		non_primes = []
		for i in range(2, n + 1):
			if i not in non_primes:
				primes.append(i)
				for j in range(i*i, n + 1, i):
					non_primes.append(j)

		return primes

	'''Checks if the two numbers are coprime or not'''
	def check_coprime(self, n1, n2):
		while n2 != 0:
			n1, n2 = n2, (n1 % n2)

		return n1 == 1

	'''
		Method used to generate n, which in turn is the product of the 
		prime numbers p and q.
		1. checks congruent 3 modulo 4
		2. checks if the value is greater than a threshold
		3. Also checks if p * q is coprime to the seed value
	'''
	def generateN(self):
		threshold = 5000
		primes = self.get_primes(10000)
		p, q = -1, -1
		while True:
			p = random.choice(primes)
			if p % 4 == 3 and p > threshold:
				break

		while True:
			q = random.choice(primes)
			if q % 4 == 3 and q > threshold:
				if p != q and self.check_coprime(self.seed, p * q):
					break

		return p * q

	'''
		Method used to generate a pseudo random number
	'''
	def generate(self):
		x = self.seed
		while True:
			yield x

			x = (x ** 2) % self.n

############################### Randomness tests ###################################
def RunsTest(seq):
	C = False
	if len(seq) >= 50:
		C = True
	stats = runstest_1samp(seq, correction = C)
	return {'z-score': stats[0], 'p-score': stats[1]}

def KSTest(seq, div):
	seq.sort()
	seq = [x / div for x in seq]
	N = len(seq)
	D_plus, D_minus = -10000000000000000, -10000000000000000

	for i in range(1, N + 1):
	    x = i / N - seq[i-1]
	    D_plus = max(D_plus, x)

	for i in range(1, N + 1):
		y = (i - 1) / N
		y = seq[i - 1] - y
		D_minus = max(D_minus, y)

	ans = max(math.sqrt(N) * D_plus, math.sqrt(N) * D_minus)
	return ans

######################################### main ###################################################
seed = int(input("Enter the seed value for the two PRNGs: "))
n = int(input("Enter the length of random number sequence: "))

LC = LinearCongruential(2**16, 101427, 321, seed)
bbs = BlumBlumShub(seed, random = True)

generator_LC = LC.generate()
generator_bbs = bbs.generate()

LC_sequence, bbs_sequence = [], []
for i in range(n):
	LC_sequence.append(next(generator_LC))
	bbs_sequence.append(next(generator_bbs))

print('Linear Congruential PRNG: ', LC_sequence)
print("##############################################################################")
print("Runs Test:")
print(RunsTest(LC_sequence))
print()
print("Kolmogorov-Smirnov Test: ")
print('D: ', KSTest(LC_sequence, 2**16), '}')
print()

print('BlumBlumShub PRNG: ', bbs_sequence)
print("##############################################################################")
print("Runs Test:")
print(RunsTest(bbs_sequence))
print()
print("Kolmogorov-Smirnov Test: ")
print('{D: ', KSTest(bbs_sequence, bbs.n), '}')
print()

