import sympy
import math


# ============================================================
# SETTINGS
# ============================================================

MODULI_SETS = [
    [3],
    [3, 5],
    [3, 5, 7],
    [3, 5, 7, 11],
    [3, 5, 7, 11, 13],

    [17],
    [17, 19],
    [17, 19, 23],
    [17, 19, 23, 29],
    [17, 19, 23, 29, 31],

    [3, 5, 7, 11, 13, 17, 19, 23, 29, 31],
]

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

print("=" * 80)
print("INCREMENTAL FINGERPRINT FACTOR EXPERIMENT")
print("=" * 80)

print(f"p = {p}")
print(f"q = {q}")
print(f"n = {n}")
print()

primes = list(sympy.primerange(2, PRIME_LIMIT + 1))

print(f"Prime candidates <= {PRIME_LIMIT}: {len(primes)}")
print()


# ============================================================
# TEST EACH MODULUS SET
# ============================================================

for moduli in MODULI_SETS:

    M = math.prod(moduli)

    print()
    print("=" * 80)
    print(f"MODULI = {moduli}")
    print("=" * 80)

    print(f"M = {M:,}")

    target = fingerprint(n, moduli)

    print(f"F(n) = {target}")
    print()

    compatible_p = []
    compatible_pairs = []

    # --------------------------------------------------------
    # For each possible prime p:
    #
    #     p*q = n (mod M)
    #
    # therefore:
    #
    #     q = n * p^-1 (mod M)
    # --------------------------------------------------------

    for candidate_p in primes:

        if math.gcd(candidate_p, M) != 1:
            continue

        rp = candidate_p % M

        rq = (n * pow(rp, -1, M)) % M

        # q is represented by:
        #
        #     q = rq + k*M
        #
        # For this small experiment, M may be larger than q,
        # in which case rq itself is the only possible q
        # below PRIME_LIMIT.

        possible_q = []

        k_min = 0

        if rq == 0:
            k_min = 0

        while True:

            candidate_q = rq + k_min * M

            if candidate_q > PRIME_LIMIT:
                break

            if candidate_q >= 2 and candidate_q in primes:
                possible_q.append(candidate_q)

            k_min += 1

        if possible_q:

            compatible_p.append(candidate_p)

            for candidate_q in possible_q:

                if fingerprint(
                    candidate_p * candidate_q,
                    moduli
                ) == target:

                    compatible_pairs.append(
                        (
                            candidate_p,
                            candidate_q,
                            candidate_p * candidate_q
                        )
                    )

    # --------------------------------------------------------
    # RESULTS
    # --------------------------------------------------------

    print(f"Possible p values     : {len(compatible_p)}")
    print(f"Possible prime pairs  : {len(compatible_pairs)}")

    print()

    if compatible_pairs:

        for a, b, product in compatible_pairs:

            marker = ""

            if {a, b} == {p, q}:
                marker = " <-- ACTUAL"

            print(
                f"{a:6} * {b:6} = "
                f"{product:12}"
                f"{marker}"
            )
    else:

        print("No compatible prime pairs.")


# ============================================================
# DIRECT TWO-STAGE TEST
# ============================================================

print()
print("=" * 80)
print("TWO-STAGE TEST")
print("=" * 80)

M1 = math.prod([3, 5, 7, 11, 13])
M2 = math.prod([17, 19, 23, 29, 31])

print(f"M1 = {M1:,}")
print(f"M2 = {M2:,}")
print()

stage1 = []
stage2 = []

for candidate_p in primes:

    # --------------------------------------------------------
    # FIRST FINGERPRINT
    # --------------------------------------------------------

    if math.gcd(candidate_p, M1) != 1:
        continue

    rq1 = (n * pow(candidate_p, -1, M1)) % M1

    # --------------------------------------------------------
    # SECOND FINGERPRINT
    # --------------------------------------------------------

    if math.gcd(candidate_p, M2) != 1:
        continue

    rq2 = (n * pow(candidate_p, -1, M2)) % M2

    stage1.append(
        (candidate_p, rq1)
    )

    stage2.append(
        (candidate_p, rq2)
    )


print(f"Candidate p values tested: {len(primes)}")
print(f"Stage-1 results           : {len(stage1)}")
print(f"Stage-2 results           : {len(stage2)}")

print()

print("Actual p:")
print(p)

print()

for candidate_p, rq1 in stage1:

    if candidate_p == p:

        print(
            f"Actual p = {candidate_p}"
        )

        print(
            f"Required q mod M1 = {rq1}"
        )

        break

for candidate_p, rq2 in stage2:

    if candidate_p == p:

        print(
            f"Actual p = {candidate_p}"
        )

        print(
            f"Required q mod M2 = {rq2}"
        )

        break


# ============================================================
# CRT COMBINATION OF THE TWO q RESIDUES
# ============================================================

print()
print("=" * 80)
print("CRT COMBINATION FOR ACTUAL p")
print("=" * 80)

rq1 = q % M1
rq2 = q % M2

print(f"q mod M1 = {rq1}")
print(f"q mod M2 = {rq2}")

print()

# Since gcd(M1, M2) = 1:
#
# q = rq1 (mod M1)
# q = rq2 (mod M2)

# q = rq1 + M1*t
#
# M1*t = rq2-rq1 (mod M2)

t = (
    (rq2 - rq1)
    * pow(M1, -1, M2)
) % M2

q_reconstructed = rq1 + M1 * t

M_combined = M1 * M2

print(f"Combined modulus = {M_combined:,}")
print(f"Reconstructed q   = {q_reconstructed}")
print(f"Actual q          = {q}")
print()

print(
    f"Correct = "
    f"{q_reconstructed == q}"
)
