import sympy
import math


# ============================================================
# TWO MODULI SETS
# ============================================================

MODULI1 = [3, 5, 7, 11, 13]
MODULI2 = [7, 11, 13, 17]


def fingerprint(x, moduli):
    return [x % m for m in moduli]


# ============================================================
# RANDOM SEMIPRIME
# ============================================================

p = sympy.randprime(100, 1000)
q = sympy.randprime(100, 1000)

n = p * q


fp1 = fingerprint(n, MODULI1)
fp2 = fingerprint(n, MODULI2)


print("=" * 70)
print("RANDOM SEMIPRIME")
print("=" * 70)

print(f"p = {p}")
print(f"q = {q}")
print(f"n = {n}")
print()


# ============================================================
# FINGERPRINT 1
# ============================================================

print("=" * 70)
print("FINGERPRINT 1")
print("=" * 70)

print(f"MODULI1 = {MODULI1}")
print(f"F1(n)   = {fp1}")
print()


# ============================================================
# FINGERPRINT 2
# ============================================================

print("=" * 70)
print("FINGERPRINT 2")
print("=" * 70)

print(f"MODULI2 = {MODULI2}")
print(f"F2(n)   = {fp2}")
print()


# ============================================================
# SEARCH
# ============================================================

matches1 = []
matches2 = []
matches_both = []


for x in range(10, n + 1):

    f1 = fingerprint(x, MODULI1)
    f2 = fingerprint(x, MODULI2)

    if f1 == fp1:
        matches1.append(x)

    if f2 == fp2:
        matches2.append(x)

    if f1 == fp1 and f2 == fp2:
        matches_both.append(x)


# ============================================================
# RESULTS
# ============================================================

print("=" * 70)
print("RESULTS")
print("=" * 70)

print(f"Numbers checked      : {n - 9}")
print()

print(f"Fingerprint 1 matches: {len(matches1)}")
print(f"Fingerprint 2 matches: {len(matches2)}")
print(f"Both fingerprints    : {len(matches_both)}")
print()


# ============================================================
# PRINT MATCHES
# ============================================================

print("First 100 matches for F1:")
print(matches1[:100])

print()

print("First 100 matches for F2:")
print(matches2[:100])

print()

print("First 100 matches for BOTH:")
print(matches_both[:100])

print()


# ============================================================
# PRIME MATCHES
# ============================================================

prime_matches1 = [x for x in matches1 if sympy.isprime(x)]
prime_matches2 = [x for x in matches2 if sympy.isprime(x)]
prime_matches_both = [x for x in matches_both if sympy.isprime(x)]


print("=" * 70)
print("PRIME MATCHES")
print("=" * 70)

print(f"Prime matches F1    : {len(prime_matches1)}")
print(f"Prime matches F2    : {len(prime_matches2)}")
print(f"Prime matches BOTH  : {len(prime_matches_both)}")

print()

print("Prime matches BOTH:")
print(prime_matches_both)

print()


# ============================================================
# CHECK ACTUAL FACTORS
# ============================================================

print("=" * 70)
print("ACTUAL FACTORS")
print("=" * 70)

print(f"p = {p}")
print(f"q = {q}")

print()

print(f"p in F1 matches   : {p in matches1}")
print(f"q in F1 matches   : {q in matches1}")

print(f"p in F2 matches   : {p in matches2}")
print(f"q in F2 matches   : {q in matches2}")

print(f"p in BOTH        : {p in matches_both}")
print(f"q in BOTH        : {q in matches_both}")
