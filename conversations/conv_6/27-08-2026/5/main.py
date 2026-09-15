import sympy
import math


# ============================================================
# TWO FINGERPRINT SYSTEMS
# ============================================================

MODULI1 = [3, 5, 7, 11, 13]
MODULI2 = [7, 11, 13, 17]

MODULI_ALL = sorted(set(MODULI1 + MODULI2))

M1 = math.prod(MODULI1)
M2 = math.prod(MODULI2)
M = math.prod(MODULI_ALL)


def fingerprint(x, moduli):
    return [x % m for m in moduli]


# ============================================================
# RANDOM SEMIPRIME
# ============================================================

p = sympy.randprime(100, 1000)
q = sympy.randprime(100, 1000)

n = p * q


print("=" * 70)
print("RANDOM SEMIPRIME")
print("=" * 70)

print(f"p = {p}")
print(f"q = {q}")
print(f"n = {n}")
print()

print(f"MODULI1 = {MODULI1}")
print(f"F1(n)   = {fingerprint(n, MODULI1)}")
print()

print(f"MODULI2 = {MODULI2}")
print(f"F2(n)   = {fingerprint(n, MODULI2)}")
print()

print(f"Combined moduli M = {MODULI_ALL}")
print(f"Combined period M = {M}")
print(f"n mod M = {n % M}")
print()


# ============================================================
# ACTUAL FACTORS
# ============================================================

print("=" * 70)
print("ACTUAL FACTORS")
print("=" * 70)

print(f"p mod M = {p % M}")
print(f"q mod M = {q % M}")

print(
    f"(p * q) mod M = {(p * q) % M}"
)

print()


# ============================================================
# TEST EVERY PRIME p <= sqrt(n)
# ============================================================

limit = math.isqrt(n)

prime_candidates = list(sympy.primerange(2, limit + 1))

compatible = []


for candidate_p in prime_candidates:

    # p must be invertible modulo M.
    if math.gcd(candidate_p, M) != 1:
        continue

    rp = candidate_p % M

    # From:
    #
    #     p*q = n (mod M)
    #
    # therefore:
    #
    #     q = n * p^-1 (mod M)

    rq = (n * pow(rp, -1, M)) % M

    # ========================================================
    # Now q must be:
    #
    #     q = rq + k*M
    #
    # and also:
    #
    #     q > p
    #
    #     p*q = n
    # ========================================================

    # Search q values congruent to rq.
    #
    # Since q <= n/p, only a few values can exist.

    q_min = candidate_p
    q_max = n // candidate_p

    if rq == 0:
        k_start = 0
    else:
        k_start = (q_min - rq + M - 1) // M

    k_end = (q_max - rq) // M

    for k in range(k_start, k_end + 1):

        candidate_q = rq + k * M

        if candidate_q < candidate_p:
            continue

        if candidate_p * candidate_q != n:
            continue

        if not sympy.isprime(candidate_q):
            continue

        compatible.append(
            (
                candidate_p,
                candidate_q,
                rp,
                rq,
                k
            )
        )


# ============================================================
# RESULTS
# ============================================================

print("=" * 70)
print("FACTOR SEARCH")
print("=" * 70)

print(f"Prime p candidates <= sqrt(n): {len(prime_candidates)}")
print(f"Compatible factor pairs       : {len(compatible)}")
print()

for cp, cq, rp, rq, k in compatible:

    print(
        f"p = {cp:6}  "
        f"q = {cq:6}  "
        f"p mod M = {rp:6}  "
        f"q mod M = {rq:6}  "
        f"k = {k}"
    )


# ============================================================
# SHOW THE RESIDUE RELATION FOR EVERY PRIME
# ============================================================

print()
print("=" * 70)
print("PRIME -> REQUIRED q RESIDUE")
print("=" * 70)

for candidate_p in prime_candidates:

    if math.gcd(candidate_p, M) != 1:
        continue

    rp = candidate_p % M
    rq = (n * pow(rp, -1, M)) % M

    marker = ""

    if candidate_p == p:
        marker = " <-- ACTUAL p"

    print(
        f"p={candidate_p:5}  "
        f"required q mod {M} = {rq:6}"
        f"{marker}"
    )
