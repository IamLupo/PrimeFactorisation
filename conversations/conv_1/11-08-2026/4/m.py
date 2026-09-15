from math import prod

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

def macmahon_M3_direct(n):
    total = 0

    for s1 in range(1, n):
        for s2 in range(s1 + 1, n):
            for s3 in range(s2 + 1, n + 1):

                for m1 in range(1, n // s1 + 1):
                    for m2 in range(1, n // s2 + 1):

                        rem = n - m1*s1 - m2*s2

                        if rem > 0 and rem % s3 == 0:
                            m3 = rem // s3
                            total += m1 * m2 * m3

    return total

def sigma_k_from_factors(factors, k):
    """
    factors = dict {prime: exponent}
    computes sigma_k(n)
    """

    terms = []

    for p, a in factors.items():
        numerator = p**(k * (a + 1)) - 1
        denominator = p**k - 1
        terms.append(numerator // denominator)

    return prod(terms)

def calc_M1_old(primes):
	M1 = 1

	for p in primes:
		M1 *= p + 1

	return M1

def calc_M1(primes):
	M1 = 1

	for base, power in primes.items():
		x = 1
		for p in range(power):
			x += pow(base, p + 1)
		M1 *= x

	return M1

def calc_M2(n, primes, M1):
	M2 = 1
	x = 1

	for base, power in primes.items():
		y = 1

		for p in range(power):
			y += pow(base, (p + 1) * 3)

		x *= y

	return (x - (2 * n * M1) + M1) // 8

def calc_x(primes):
    x = 1

    for base, power in primes.items():
        y = 1
        for p in range(power):
            y += pow(base, (p + 1) * 3)
        x *= y

    return x

def get_R(n, M1, M2):
    return M1 * (n**2 - 3*n + 2) - (8 * M2)

def calc_xy(n, p, q):
    if(n % 4 == 3):
        x = (q - p + 6) // 4
        y = (p + q) // 4
    else:
        x = (3 * p - q + 6) // 4
        y = (3 * p + q) // 4

    return x, y