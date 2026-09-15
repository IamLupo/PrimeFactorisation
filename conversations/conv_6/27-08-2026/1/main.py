import random
import sympy


# ------------------------------------------------------------
# 1. Generate random primes
# ------------------------------------------------------------

p = sympy.randprime(700, 1000)
q = sympy.randprime(700, 1000)

n = p * q

print("=" * 70)
print("RANDOM SEMIPRIME")
print("=" * 70)

print(f"p = {p}")
print(f"q = {q}")
print(f"n = {n}")
print()


# ------------------------------------------------------------
# 2. Fingerprint function
# ------------------------------------------------------------

def fingerprint(x): return [ x % 2, x % 3, x % 5, x % 7, x % 11, x % 13 ]


target = fingerprint(n)

print("Target fingerprint:")
print(f"F({n}) = {target}")
print()


# ------------------------------------------------------------
# 3. Search from 10 through n
# ------------------------------------------------------------

matches = []
prime_matches = []

for x in range(10, n + 1):

    if fingerprint(x) == target:

        matches.append(x)

        if sympy.isprime(x):
            prime_matches.append(x)


# ------------------------------------------------------------
# 4. Results
# ------------------------------------------------------------

print("=" * 70)
print("RESULTS")
print("=" * 70)

print(f"Numbers checked       : {n - 9}")
print(f"Fingerprint matches   : {len(matches)}")
print(f"Prime matches         : {len(prime_matches)}")
print()

print("First 100 matching values:")
print(matches[:100])

print()

print("Prime matching values:")
print(prime_matches)

print()

print("Original factors:")
print(f"p = {p}")
print(f"q = {q}")

print()

print("Are p and q in the matching prime set?")
print(f"p in matches: {p in prime_matches}")
print(f"q in matches: {q in prime_matches}")