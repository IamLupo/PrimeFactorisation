import sympy
import math


# ============================================================
# SETTINGS
# ============================================================

MODULI = [3, 5, 7, 11, 13]

M = math.prod(MODULI)


# ============================================================
# FINGERPRINT
# ============================================================

def fingerprint(x):
    return [x % m for m in MODULI]


# ============================================================
# RANDOM SEMIPRIME
# ============================================================

p = sympy.randprime(100, 1000)
q = sympy.randprime(100, 1000)

n = p * q

target_fp = fingerprint(n)
target_mod = n % M


print("=" * 70)
print("RANDOM SEMIPRIME")
print("=" * 70)

print(f"p = {p}")
print(f"q = {q}")
print(f"n = {n}")
print()

print("Moduli:")
print(MODULI)

print(f"Combined modulus M = {M}")
print()

print("Target fingerprint:")
print(f"F({n}) = {target_fp}")
print()

print(f"n mod M = {target_mod}")
print()


# ============================================================
# ACTUAL FACTOR RESIDUES
# ============================================================

rp_actual = p % M
rq_actual = q % M

print("=" * 70)
print("ACTUAL FACTOR RESIDUES")
print("=" * 70)

print(f"p mod M = {rp_actual}")
print(f"q mod M = {rq_actual}")

print(
    f"({rp_actual} * {rq_actual}) mod {M} = "
    f"{(rp_actual * rq_actual) % M}"
)

print()


# ============================================================
# SOLVE:
#
#     rp * rq = target_mod (mod M)
#
# For each rp, solve the linear congruence for rq.
# ============================================================

def solve_linear_congruence(a, b, m):
    """
    Solve:

        a*x = b (mod m)

    Returns all solutions x in [0, m-1].
    """

    g = math.gcd(a, m)

    # No solution
    if b % g != 0:
        return []

    # Reduce equation
    a2 = a // g
    b2 = b // g
    m2 = m // g

    # Now gcd(a2, m2) == 1
    inv = pow(a2, -1, m2)

    x0 = (b2 * inv) % m2

    # There are exactly g solutions modulo m
    solutions = []

    for k in range(g):
        x = x0 + k * m2
        solutions.append(x)

    return solutions


# ============================================================
# FIND COMPATIBLE FACTOR RESIDUES
# ============================================================

matches = []

for rp in range(M):

    rq_solutions = solve_linear_congruence(
        rp,
        target_mod,
        M
    )

    for rq in rq_solutions:

        # Safety check
        assert (rp * rq) % M == target_mod

        matches.append((rp, rq))


# ============================================================
# RESULTS
# ============================================================

print("=" * 70)
print("RESULTS")
print("=" * 70)

print(f"Total possible (rp,rq) pairs : {M * M:,}")
print(f"Compatible pairs             : {len(matches):,}")
print()

actual_found = (rp_actual, rq_actual) in matches
swapped_found = (rq_actual, rp_actual) in matches

print(f"Actual pair found            : {actual_found}")
print(f"Swapped pair found           : {swapped_found}")
print()


# ============================================================
# FIRST MATCHES
# ============================================================

print("=" * 70)
print("FIRST 100 COMPATIBLE PAIRS")
print("=" * 70)

for i, (rp, rq) in enumerate(matches[:100], 1):

    marker = ""

    if (rp, rq) == (rp_actual, rq_actual):
        marker = "  <-- ACTUAL"

    elif (rp, rq) == (rq_actual, rp_actual):
        marker = "  <-- ACTUAL SWAPPED"

    print(f"{i:4}: rp = {rp:5}, rq = {rq:5}{marker}")


# ============================================================
# UNIQUE RESIDUES
# ============================================================

rp_values = sorted(set(rp for rp, rq in matches))
rq_values = sorted(set(rq for rp, rq in matches))

print()
print("=" * 70)
print("UNIQUE RESIDUES")
print("=" * 70)

print(f"Unique rp values: {len(rp_values):,}")
print(f"Unique rq values: {len(rq_values):,}")

print()

print("First 100 possible rp:")
print(rp_values[:100])

print()

print("First 100 possible rq:")
print(rq_values[:100])