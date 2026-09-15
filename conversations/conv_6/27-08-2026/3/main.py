import sympy
import math


MODULI = [3, 5, 7, 11, 13]
M = math.prod(MODULI)


def fingerprint(x):
    return [x % m for m in MODULI]


# ============================================================
# RANDOM SEMIPRIME
# ============================================================

p = sympy.randprime(100, 1000)
q = sympy.randprime(100, 1000)

n = p * q

target = n % M

print("=" * 70)
print("RANDOM SEMIPRIME")
print("=" * 70)

print(f"p = {p}")
print(f"q = {q}")
print(f"n = {n}")
print()

print(f"M = {M}")
print(f"n mod M = {target}")
print(f"F(n) = {fingerprint(n)}")
print()


# ============================================================
# PRIME CANDIDATES
# ============================================================

limit = math.isqrt(n)

prime_candidates = list(sympy.primerange(2, limit + 1))


print("=" * 70)
print("PRIME CANDIDATES")
print("=" * 70)

print(f"sqrt(n)               = {math.sqrt(n):.3f}")
print(f"Prime candidates      = {len(prime_candidates)}")
print()


# ============================================================
# TEST EACH POSSIBLE PRIME p'
# ============================================================

compatible = []

for candidate_p in prime_candidates:

    # Candidate must be invertible modulo M
    if math.gcd(candidate_p, M) != 1:
        continue

    inverse = pow(candidate_p, -1, M)

    candidate_q_residue = (target * inverse) % M

    # q must satisfy:
    #
    # q = candidate_q_residue + M*k
    #
    # But because q > candidate_p and candidate_p*q = n,
    # we can calculate the exact q directly.
    if n % candidate_p != 0:
        continue

    candidate_q = n // candidate_p

    # Verify the residue relationship
    if candidate_q % M != candidate_q_residue:
        continue

    if sympy.isprime(candidate_q):
        compatible.append(
            (
                candidate_p,
                candidate_q,
                candidate_p % M,
                candidate_q % M
            )
        )


# ============================================================
# RESULTS
# ============================================================

print("=" * 70)
print("COMPATIBLE PRIME FACTOR PAIRS")
print("=" * 70)

print(f"Compatible prime pairs: {len(compatible)}")
print()

for pair in compatible:
    cp, cq, rp, rq = pair

    print(
        f"p = {cp:6}  "
        f"q = {cq:6}  "
        f"rp = {rp:6}  "
        f"rq = {rq:6}"
    )


# ============================================================
# NEW EXPERIMENT:
# PRIME RESIDUE CLASSES BELOW sqrt(n)
# ============================================================

print()
print("=" * 70)
print("PRIME RESIDUE CLASSES")
print("=" * 70)

prime_residue_candidates = []

for prime in prime_candidates:

    if math.gcd(prime, M) != 1:
        continue

    rp = prime % M
    rq = (target * pow(rp, -1, M)) % M

    prime_residue_candidates.append(
        (prime, rp, rq)
    )

print(
    f"Prime candidates with valid residue: "
    f"{len(prime_residue_candidates)}"
)

print()

for item in prime_residue_candidates[:100]:
    prime, rp, rq = item

    marker = ""

    if prime == p:
        marker = " <-- ACTUAL p"

    if prime == q:
        marker = " <-- ACTUAL q"

    print(
        f"prime={prime:5}  "
        f"rp={rp:5}  "
        f"required rq={rq:5}"
        f"{marker}"
    )
