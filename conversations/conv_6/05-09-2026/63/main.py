# ================================================================
# START EXPERIMENT 146
# Symmetric CRT residue coupling
#
# Goal:
#
#   Use identical radix systems for p and q:
#
#       R = [17, 43, 59]
#
#   so both factors have residues modulo
#
#       M = 17*43*59 = 43129.
#
#   For every possible p residue P (mod M),
#
#       q = n * P^{-1} (mod M)
#
#   is forced.
#
#   We then ask:
#
#       1. How many p residues are compatible?
#       2. How many p/q residue pairs remain?
#       3. Can the exact factor be recovered by lifting the
#          residue classes?
#
# No carry-bound interval is used.
#
# This is a control for the hypothesis that the important
# information is the coupling between the two factor residue
# systems rather than the individual carry matrices.
# ================================================================

import math
import random
import time


# ------------------------------------------------
# Configuration
# ------------------------------------------------

RADICES = [17, 43, 59]

M = 1

for r in RADICES:
    M *= r

BIT_SIZES = [30, 36, 42, 48, 54]


# ================================================================
# Miller-Rabin
# ================================================================

def is_probable_prime(n):

    if n < 2:
        return False

    small_primes = [
        2, 3, 5, 7, 11, 13, 17,
        19, 23, 29, 31, 37
    ]

    for p in small_primes:

        if n == p:
            return True

        if n % p == 0:
            return False

    d = n - 1
    s = 0

    while d % 2 == 0:
        d //= 2
        s += 1

    for a in [2, 3, 5, 7, 11, 13, 17]:

        if a >= n:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        witness_passed = False

        for _ in range(s - 1):

            x = (x * x) % n

            if x == n - 1:
                witness_passed = True
                break

        if not witness_passed:
            return False

    return True


def random_prime(bits):

    while True:

        x = random.getrandbits(bits)

        x |= 1
        x |= (1 << (bits - 1))

        if is_probable_prime(x):
            return x


def random_semiprime(bits):

    p_bits = bits // 2
    q_bits = bits - p_bits

    while True:

        p = random_prime(p_bits)
        q = random_prime(q_bits)

        if p != q:
            return p, q


# ================================================================
# CRT helpers
# ================================================================

def crt_pair(a1, m1, a2, m2):

    """
    Solve

        x = a1 mod m1
        x = a2 mod m2

    for coprime m1,m2.
    """

    inv = pow(m1, -1, m2)

    t = (
        (a2 - a1)
        * inv
    ) % m2

    x = a1 + m1 * t

    return x % (m1 * m2)


def residue_vector_to_modulus(residues):

    x = residues[0]
    modulus = RADICES[0]

    for i in range(1, len(RADICES)):

        x = crt_pair(
            x,
            modulus,
            residues[i],
            RADICES[i]
        )

        modulus *= RADICES[i]

    return x % modulus


def residue_vector(x):

    return [
        x % r
        for r in RADICES
    ]


# ================================================================
# Random factor test
# ================================================================

def invertible_mod_M(x):

    return math.gcd(x, M) == 1


# ================================================================
# Build all invertible residue classes
# ================================================================

def build_residue_table():

    table = {}

    for p_res in range(M):

        if not invertible_mod_M(p_res):
            continue

        table[p_res] = p_res

    return table


# ================================================================
# Symmetric residue coupling
# ================================================================

def coupled_residues(n):

    """

    For every invertible residue

        p = P (mod M),

    compute the forced residue

        q = n * P^{-1} (mod M).

    """

    candidates = []

    for p_res in range(M):

        if not invertible_mod_M(p_res):
            continue

        q_res = (
            n
            * pow(
                p_res,
                -1,
                M
            )
        ) % M

        candidates.append(
            (
                p_res,
                q_res
            )
        )

    return candidates


# ================================================================
# Factor lifting
# ================================================================

def lift_residue(
    residue,
    modulus,
    n
):

    """
    Enumerate p = residue + t*modulus below sqrt(n).

    This is the only lifting step.
    """

    root = math.isqrt(n)

    if residue == 0:

        first = modulus

    else:

        first = residue

        if first < 2:
            first += modulus

    if first > root:
        return []

    count = (
        (root - first) // modulus
    ) + 1

    return [
        first + i * modulus
        for i in range(count)
    ]


# ================================================================
# Main
# ================================================================

print("=" * 72)
print("START EXPERIMENT 146")
print("Symmetric CRT residue coupling")
print("=" * 72)
print()

print("RADICES =", RADICES)
print("M       =", M)
print()

# ------------------------------------------------
# Static information.
# ------------------------------------------------

invertible_count = sum(
    1
    for x in range(M)
    if math.gcd(x, M) == 1
)

print(
    "INVERTIBLE RESIDUES MOD M =",
    invertible_count
)

print()


for bits in BIT_SIZES:

    start = time.perf_counter()

    print("=" * 72)
    print(
        f"GENERATING {bits}-BIT SEMIPRIME"
    )
    print("=" * 72)
    print()

    p, q = random_semiprime(bits)
    n = p * q

    print("-" * 72)
    print("n bits =", n.bit_length())
    print("n      =", n)
    print("true p =", p)
    print("true q =", q)
    print()

    print(
        "true p mod M =",
        p % M
    )

    print(
        "true q mod M =",
        q % M
    )

    print()

    # ------------------------------------------------------------
    # Direct modular coupling.
    # ------------------------------------------------------------

    if math.gcd(p, M) != 1:

        print(
            "TRUE p IS NOT INVERTIBLE MOD M"
        )

        print(
            "This case is unsuitable for the inverse coupling."
        )

        print()

        continue

    forced_q_res = (
        n
        * pow(
            p % M,
            -1,
            M
        )
    ) % M

    print(
        "FORCED q RESIDUE FROM TRUE p =",
        forced_q_res
    )

    print(
        "ACTUAL q RESIDUE             =",
        q % M
    )

    print(
        "TRUE MODULAR COUPLING:",
        "PASS"
        if forced_q_res == q % M
        else "FAIL"
    )

    print()

    # ------------------------------------------------------------
    # Enumerate every admissible p residue.
    # ------------------------------------------------------------

    coupling_start = time.perf_counter()

    pairs = coupled_residues(n)

    coupling_elapsed = (
        time.perf_counter()
        - coupling_start
    )

    print(
        "COUPLED RESIDUE PAIRS =",
        len(pairs)
    )

    print(
        "COUPLING ENUMERATION RUNTIME =",
        f"{coupling_elapsed:.6f} s"
    )

    print()

    # ------------------------------------------------------------
    # Check uniqueness of q residue mapping.
    # ------------------------------------------------------------

    q_residues = {}

    for p_res, q_res in pairs:

        q_residues.setdefault(
            q_res,
            0
        )

        q_residues[q_res] += 1

    max_q_multiplicity = max(
        q_residues.values()
    )

    distinct_q_residues = len(
        q_residues
    )

    print(
        "DISTINCT q RESIDUES =",
        distinct_q_residues
    )

    print(
        "MAX q-RESIDUE MULTIPLICITY =",
        max_q_multiplicity
    )

    print()

    # ------------------------------------------------------------
    # Find the true residue pair.
    # ------------------------------------------------------------

    true_pair_present = (
        (
            p % M,
            q % M
        )
        in set(pairs)
    )

    swapped_pair_present = (
        (
            q % M,
            p % M
        )
        in set(pairs)
    )

    print(
        "TRUE RESIDUE PAIR PRESENT =",
        true_pair_present
    )

    print(
        "SWAPPED RESIDUE PAIR PRESENT =",
        swapped_pair_present
    )

    print()

    # ------------------------------------------------------------
    # Lift every p residue below sqrt(n).
    #
    # This measures exactly how much work remains after the
    # symmetric modular coupling.
    # ------------------------------------------------------------

    root = math.isqrt(n)

    total_lift_points = 0
    factor_hits = []

    lift_start = time.perf_counter()

    for p_res, q_res in pairs:

        p_points = lift_residue(
            p_res,
            M,
            n
        )

        total_lift_points += len(p_points)

        for p_candidate in p_points:

            if p_candidate <= 1:
                continue

            if n % p_candidate != 0:
                continue

            q_candidate = n // p_candidate

            if q_candidate % M != q_res:
                continue

            factor_hits.append(
                (
                    p_candidate,
                    q_candidate
                )
            )

    lift_elapsed = (
        time.perf_counter()
        - lift_start
    )

    factor_hits = sorted(
        set(
            tuple(sorted(pair))
            for pair in factor_hits
        )
    )

    print("LIFTING")
    print()

    print(
        "sqrt(n) =",
        root
    )

    print(
        "residue modulus M =",
        M
    )

    print(
        "sqrt(n) / M =",
        f"{root / M:.6f}"
    )

    print(
        "total p lift points =",
        total_lift_points
    )

    print(
        "lift runtime =",
        f"{lift_elapsed:.6f} s"
    )

    print()

    print(
        "EXACT FACTOR HITS =",
        len(factor_hits)
    )

    for pair in factor_hits:

        print(
            "  ",
            pair
        )

    print()

    # ------------------------------------------------------------
    # Compare the residue-class work with ordinary sqrt search.
    # ------------------------------------------------------------

    print("COMPLEXITY COMPARISON")
    print()

    print(
        "ordinary p candidates below sqrt =",
        root
    )

    print(
        "symmetric CRT lift candidates =",
        total_lift_points
    )

    if root:

        print(
            "reduction factor =",
            f"{root / max(total_lift_points, 1):.6f}"
        )

    print()

    # ------------------------------------------------------------
    # Show the actual residue-vector structure.
    # ------------------------------------------------------------

    print("TRUE RADIX RESIDUES")
    print()

    print(
        "p residues =",
        residue_vector(p)
    )

    print(
        "q residues =",
        residue_vector(q)
    )

    reconstructed_p_res = residue_vector_to_modulus(
        residue_vector(p)
    )

    reconstructed_q_res = residue_vector_to_modulus(
        residue_vector(q)
    )

    print(
        "CRT p residue =",
        reconstructed_p_res
    )

    print(
        "CRT q residue =",
        reconstructed_q_res
    )

    print()

    # ------------------------------------------------------------
    # Information threshold.
    # ------------------------------------------------------------

    print("MODULUS THRESHOLD")
    print()

    if M > root:

        print(
            "M > sqrt(n): "
            "each invertible residue class contains at most "
            "one p below sqrt(n)."
        )

    else:

        print(
            "M <= sqrt(n): "
            "multiple p values can occupy one residue class."
        )

        print(
            "maximum approximate lift count =",
            root // M + 1
        )

    print()

    elapsed = time.perf_counter() - start

    print(
        "TOTAL RUNTIME =",
        f"{elapsed:.6f} s"
    )

    print()
    print("-" * 72)
    print()


print("=" * 72)
print("FINISHED EXPERIMENT 146")
print("=" * 72)
