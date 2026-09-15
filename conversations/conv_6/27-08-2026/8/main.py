import sympy
import math


# ============================================================
# SETTINGS
# ============================================================

MODULI1 = [3, 5, 7, 11, 13]
MODULI2 = [17, 19, 23, 29, 31]

PRIME_LIMIT = 5000


# ============================================================
# FINGERPRINT
# ============================================================

def fingerprint(x, moduli):
    return tuple(x % m for m in moduli)


# ============================================================
# RANDOM SEMIPRIME
# ============================================================

p = sympy.randprime(1000, PRIME_LIMIT)
q = sympy.randprime(1000, PRIME_LIMIT)

n = p * q

target1 = fingerprint(n, MODULI1)
target2 = fingerprint(n, MODULI2)


# ============================================================
# HEADER
# ============================================================

print("=" * 80)
print("TWO-FINGERPRINT FACTOR-PAIR SEARCH")
print("=" * 80)

print(f"p = {p}")
print(f"q = {q}")
print(f"n = {n}")
print()

print(f"MODULI1 = {MODULI1}")
print(f"F1(n)   = {target1}")
print()

print(f"MODULI2 = {MODULI2}")
print(f"F2(n)   = {target2}")
print()


# ============================================================
# PRIME LIST
# ============================================================

primes = list(sympy.primerange(2, PRIME_LIMIT + 1))

print(f"Prime candidates <= {PRIME_LIMIT}: {len(primes)}")
print()


# ============================================================
# SEARCH
#
# We deliberately DO NOT require:
#
#     candidate_p * candidate_q == n
#
# We only require both fingerprints to match.
# ============================================================

matches1 = []
matches2 = []
matches_both = []


for i, a in enumerate(primes):

    for b in primes[i:]:

        product = a * b

        fp1 = fingerprint(product, MODULI1)
        fp2 = fingerprint(product, MODULI2)

        if fp1 == target1:
            matches1.append((a, b, product))

        if fp2 == target2:
            matches2.append((a, b, product))

        if fp1 == target1 and fp2 == target2:
            matches_both.append((a, b, product))


# ============================================================
# RESULTS
# ============================================================

print("=" * 80)
print("RESULTS")
print("=" * 80)

print(f"F1 matching prime pairs   : {len(matches1)}")
print(f"F2 matching prime pairs   : {len(matches2)}")
print(f"BOTH matching prime pairs : {len(matches_both)}")
print()


# ============================================================
# CHECK WHETHER ORIGINAL PAIR APPEARS
# ============================================================

actual_pair = tuple(sorted((p, q)))

found = False

for a, b, product in matches_both:

    if (a, b) == actual_pair:
        found = True
        break


print(f"Actual factor pair found: {found}")
print()


# ============================================================
# PRINT ALL BOTH-MATCHING PAIRS
# ============================================================

print("=" * 80)
print("PAIRS MATCHING BOTH FINGERPRINTS")
print("=" * 80)

for a, b, product in matches_both:

    marker = ""

    if (a, b) == actual_pair:
        marker = " <-- ACTUAL FACTORS"

    print(
        f"{a:6} * {b:6} = {product:12}"
        f"{marker}"
    )


# ============================================================
# DISTANCE FROM n
# ============================================================

print()
print("=" * 80)
print("DISTANCE FROM n")
print("=" * 80)

for a, b, product in matches_both:

    print(
        f"{a:6} * {b:6} = {product:12}   "
        f"difference = {product - n:+12}"
    )


# ============================================================
# FINGERPRINT PERIOD
# ============================================================

M1 = math.prod(MODULI1)
M2 = math.prod(MODULI2)

L = math.lcm(M1, M2)

print()
print("=" * 80)
print("PERIOD")
print("=" * 80)

print(f"M1 = {M1}")
print(f"M2 = {M2}")
print(f"LCM = {L}")
print()

print(
    "Every integer x satisfying both fingerprints obeys:"
)

print(
    f"x = n + k*{L}"
)
