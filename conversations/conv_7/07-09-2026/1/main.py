from sympy import isprime
from sympy import primerange
import math

def macmahon_M1(n):
	"""M1(n): sum of multiplicities for partitions into 1 part size"""
	M1 = 0
	for s1 in range(1, n+1):
		if n % s1 == 0:
			m1 = n // s1
			M1 += m1
	return M1

def macmahon_M2(n):
	"""M2(n): sum of products of multiplicities for partitions into 2 distinct part sizes"""
	M2 = 0
	# generate all pairs of distinct part sizes
	for s1 in range(1, n):
		for s2 in range(s1+1, n+1):
			# solve m1*s1 + m2*s2 = n
			max_m1 = n // s1
			for m1 in range(1, max_m1+1):
				remainder = n - m1*s1
				if remainder > 0 and remainder % s2 == 0:
					m2 = remainder // s2
					M2 += m1 * m2
	return M2

def macmahon_M3(n, M1=None, M2=None):
	"""Compute M3 using the Ono polynomial"""
	if M1 is None:
		M1 = macmahon_M1(n)
	if M2 is None:
		M2 = macmahon_M2(n, M1)
	a = 3*n**3 - 13*n**2 + 18*n - 8
	b = 12*n**2 - 120*n + 212
	c = 960
	M3 = (a*M1 + b*M2)/c
	return M3

def get_R(n, M1, M2):
	return M1 * (n**2 - 3*n + 2) - (8 * M2)

# Generate all primes less than 100
primes = list(primerange(2, 100))

print("n	p	q	S	M1	M2	R	R_mod_n")

for i, p in enumerate(primes):
	for q in primes[i+1:]:
		n = p * q
		S = p + q

		if(n < 1000):
			M1_real = macmahon_M1(n)
			M2_real = macmahon_M2(n)

			x = (p**3 + 1) * (q**3 + 1)
			
			M1 = (p + 1) * (q + 1)
			M2 = (x - (2*n*M1) + M1) // 8

			R = get_R(n, M1, M2)
			
			print(f"{n}	{p}	{q}	{S}	{M1_real}	{M1}	{M2_real}	{M2}	{R}")
			#print(f"{n}	{p}	{q}	{S}	{M1}	{M2}	{R}")