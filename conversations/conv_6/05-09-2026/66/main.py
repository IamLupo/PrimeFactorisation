# ================================================================
# START EXPERIMENT 149
# Correct discriminant residue structure:
#   S^2 - d^2 = 4n
#
# S = p + q
# d = p - q
#
# Goal:
#   Determine whether the (S,d) representation provides any
#   information reduction compared with direct (p,q) residues.
# ================================================================

import math
import random
import time


# ------------------------------------------------
# Configuration
# ------------------------------------------------

RADICES = [
    17,
    43,
    59,
    71,
    83,
    97,
]

BIT_SIZES = [
    30,
    36,
    42,
    48,
    54,
]


# ================================================================
# Miller-Rabin
# ================================================================

def is_probable_prime(n):

    if n < 2:
        return False

    small_primes = [
        2, 3, 5, 7, 11, 13, 17,
        19, 23, 29, 31, 37,
    ]

    for prime in small_primes:

        if n == prime:
            return True

        if n % prime == 0:
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
# CRT
# ================================================================

def crt_pair(a1, m1, a2, m2):
    """
    Solve:

        x = a1 (mod m1)
        x = a2 (mod m2)

    with gcd(m1,m2)=1.
    """

    inv = pow(
        m1,
        -1,
        m2
    )

    t = (
        (a2 - a1)
        * inv
    ) % m2

    x = a1 + m1 * t

    return (
        x % (m1 * m2),
        m1 * m2
    )


# ================================================================
# Local (S,d) states
# ================================================================

def local_sd_states(n, r):
    """
    Enumerate:

        S^2 - d^2 = 4n (mod r)
    """

    target = (
        4 * n
    ) % r

    states = []

    for s in range(r):

        s2 = (
            s * s
        ) % r

        for d in range(r):

            if (
                (s2 - d * d) % r
                == target
            ):

                states.append(
                    (s, d)
                )

    return states


# ================================================================
# Direct local factor states
# ================================================================

def local_factor_states(n, r):
    """
    Enumerate all invertible p residues and force q:

        q = n * p^{-1} mod r

    Then map them to:

        S = p + q
        d = p - q
    """

    states = []

    for p in range(1, r):

        if math.gcd(p, r) != 1:
            continue

        q = (
            n
            * pow(
                p,
                -1,
                r
            )
        ) % r

        s = (
            p + q
        ) % r

        d = (
            p - q
        ) % r

        states.append(
            (p, q, s, d)
        )

    return states


# ================================================================
# Local equivalence check
# ================================================================

def check_local_equivalence(n, r):

    sd_states = set(
        local_sd_states(
            n,
            r
        )
    )

    factor_states = set()

    for p, q, s, d in local_factor_states(
        n,
        r
    ):

        factor_states.add(
            (s, d)
        )

        # Because d = p-q and swapping p,q gives -d.
        factor_states.add(
            (
                s,
                (-d) % r
            )
        )

    return (
        sd_states,
        factor_states
    )


# ================================================================
# Combine local (S,d) states using CRT
# ================================================================

def combine_sd_states(n):

    current = [
        (
            0,
            0,
            1
        )
    ]

    for r in RADICES:

        local = local_sd_states(
            n,
            r
        )

        next_states = []

        for S0, d0, modulus0 in current:

            for s, d in local:

                S_new, modulus_s = crt_pair(
                    S0,
                    modulus0,
                    s,
                    r
                )

                d_new, modulus_d = crt_pair(
                    d0,
                    modulus0,
                    d,
                    r
                )

                if modulus_s != modulus_d:
                    raise RuntimeError(
                        "CRT modulus mismatch"
                    )

                next_states.append(
                    (
                        S_new,
                        d_new,
                        modulus_s
                    )
                )

        current = sorted(
            set(next_states)
        )

    return current


# ================================================================
# Direct p-residue CRT states
# ================================================================

def combined_factor_residue_states(n):

    current = [
        (
            0,
            1
        )
    ]

    for r in RADICES:

        local = []

        for p in range(1, r):

            if math.gcd(p, r) != 1:
                continue

            local.append(p)

        next_states = []

        for x, modulus in current:

            for residue in local:

                new_x, new_modulus = crt_pair(
                    x,
                    modulus,
                    residue,
                    r
                )

                next_states.append(
                    (
                        new_x,
                        new_modulus
                    )
                )

        current = sorted(
            set(next_states)
        )

    return current


# ================================================================
# Main experiment
# ================================================================

print("=" * 72)
print("START EXPERIMENT 149")
print("Correct discriminant residue structure")
print("=" * 72)
print()

print(
    "Testing S^2 - d^2 = 4n (mod r)"
)

print(
    "RADICES =",
    RADICES
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

    S_true = p + q
    d_true = p - q

    print("-" * 72)

    print(
        "n bits =",
        n.bit_length()
    )

    print(
        "n      =",
        n
    )

    print(
        "true p =",
        p
    )

    print(
        "true q =",
        q
    )

    print(
        "true S =",
        S_true
    )

    print(
        "true d =",
        d_true
    )

    print()

    # ------------------------------------------------------------
    # Global identity.
    # ------------------------------------------------------------

    lhs = (
        S_true * S_true
        - d_true * d_true
    )

    rhs = 4 * n

    print(
        "GLOBAL DISCRIMINANT IDENTITY"
    )

    print(
        "S^2 - d^2 =",
        lhs
    )

    print(
        "4n         =",
        rhs
    )

    print(
        "IDENTITY:",
        "PASS"
        if lhs == rhs
        else "FAIL"
    )

    print()

    # ------------------------------------------------------------
    # Local state comparison.
    # ------------------------------------------------------------

    product_sd = 1
    product_factor = 1

    local_equivalence = True
    local_true_present = True

    print(
        "LOCAL STATE COUNTS"
    )

    print()

    for r in RADICES:

        sd_states, factor_states = \
            check_local_equivalence(
                n,
                r
            )

        sd_count = len(sd_states)
        factor_count = len(factor_states)

        product_sd *= sd_count
        product_factor *= factor_count

        true_sd = (
            S_true % r,
            d_true % r
        )

        contains_true = (
            true_sd in sd_states
        )

        if not contains_true:
            local_true_present = False

        equivalent = (
            sd_states == factor_states
        )

        if not equivalent:
            local_equivalence = False

        print(
            f"r={r:3d}: "
            f"(S,d)={sd_count:4d}, "
            f"factor={factor_count:4d}, "
            f"equivalent={equivalent}, "
            f"true_present={contains_true}"
        )

    print()

    print(
        "PRODUCT LOCAL (S,d) STATES =",
        product_sd
    )

    print(
        "PRODUCT LOCAL FACTOR STATES =",
        product_factor
    )

    print(
        "LOCAL REPRESENTATIONS IDENTICAL =",
        local_equivalence
    )

    print(
        "TRUE LOCAL STATE PRESENT =",
        local_true_present
    )

    print()

    # ------------------------------------------------------------
    # CRT combine.
    # ------------------------------------------------------------

    crt_start = time.perf_counter()

    sd_crt = combine_sd_states(
        n
    )

    crt_elapsed = (
        time.perf_counter()
        - crt_start
    )

    print(
        "(S,d) CRT COMBINATION"
    )

    print()

    print(
        "CRT states =",
        len(sd_crt)
    )

    if sd_crt:

        print(
            "CRT modulus =",
            sd_crt[0][2]
        )

    print(
        "CRT runtime =",
        f"{crt_elapsed:.6f} s"
    )

    print()

    # ------------------------------------------------------------
    # Direct p CRT.
    # ------------------------------------------------------------

    pcrt_start = time.perf_counter()

    p_crt = combined_factor_residue_states(
        n
    )

    pcrt_elapsed = (
        time.perf_counter()
        - pcrt_start
    )

    print(
        "DIRECT p-RESIDUE CRT"
    )

    print()

    print(
        "p CRT states =",
        len(p_crt)
    )

    if p_crt:

        print(
            "p CRT modulus =",
            p_crt[0][1]
        )

    print(
        "p CRT runtime =",
        f"{pcrt_elapsed:.6f} s"
    )

    print()

    # ------------------------------------------------------------
    # CRT modulus.
    # ------------------------------------------------------------

    M = 1

    for r in RADICES:
        M *= r

    # ------------------------------------------------------------
    # True CRT state.
    # ------------------------------------------------------------

    true_S_residue = (
        S_true % M
    )

    true_d_residue = (
        d_true % M
    )

    true_sd_present = (
        (
            true_S_residue,
            true_d_residue,
            M
        )
        in sd_crt
    )

    true_sd_swapped_present = (
        (
            true_S_residue,
            (-true_d_residue) % M,
            M
        )
        in sd_crt
    )

    print(
        "TRUE CRT STATE"
    )

    print()

    print(
        "S mod M =",
        true_S_residue
    )

    print(
        "d mod M =",
        true_d_residue
    )

    print(
        "true (S,d) present =",
        true_sd_present
    )

    print(
        "swapped d sign present =",
        true_sd_swapped_present
    )

    print()

    # ------------------------------------------------------------
    # Convert every (S,d) state back to (p,q).
    #
    # Since p,q are odd and M is odd:
    #
    #     p = (S+d)/2
    #     q = (S-d)/2
    #
    # is valid modulo M.
    # ------------------------------------------------------------

    inv2 = pow(
        2,
        -1,
        M
    )

    factor_pairs_from_sd = set()

    for S_residue, d_residue, modulus in sd_crt:

        p_residue = (
            (
                S_residue
                + d_residue
            )
            * inv2
        ) % modulus

        q_residue = (
            (
                S_residue
                - d_residue
            )
            * inv2
        ) % modulus

        factor_pairs_from_sd.add(
            (
                p_residue,
                q_residue
            )
        )

    # ------------------------------------------------------------
    # Build direct (p,q) residue relation.
    # ------------------------------------------------------------

    direct_factor_pairs = set()

    for p_residue, modulus in p_crt:

        q_residue = (
            n
            * pow(
                p_residue,
                -1,
                modulus
            )
        ) % modulus

        direct_factor_pairs.add(
            (
                p_residue,
                q_residue
            )
        )

    # ------------------------------------------------------------
    # Exact equivalence test.
    # ------------------------------------------------------------

    exact_equivalence = (
        factor_pairs_from_sd
        == direct_factor_pairs
    )

    print(
        "SD -> FACTOR RESIDUE CONVERSION"
    )

    print()

    print(
        "(S,d) CRT states =",
        len(sd_crt)
    )

    print(
        "(p,q) from SD states =",
        len(factor_pairs_from_sd)
    )

    print(
        "direct (p,q) states =",
        len(direct_factor_pairs)
    )

    print()

    print(
        "EXACT CRT REPRESENTATION EQUIVALENCE =",
        "PASS"
        if exact_equivalence
        else "FAIL"
    )

    print()

    # ------------------------------------------------------------
    # Information reduction.
    # ------------------------------------------------------------

    if len(direct_factor_pairs) > 0:

        reduction = (
            len(direct_factor_pairs)
            /
            max(
                len(factor_pairs_from_sd),
                1
            )
        )

    else:

        reduction = 0.0

    print(
        "INFORMATION REDUCTION"
    )

    print()

    print(
        "direct factor-residue states =",
        len(direct_factor_pairs)
    )

    print(
        "(S,d) residue states =",
        len(factor_pairs_from_sd)
    )

    print(
        "state reduction factor =",
        f"{reduction:.6f}"
    )

    print()

    # ------------------------------------------------------------
    # Difference between true factor coordinates and SD.
    # ------------------------------------------------------------

    p_from_true_sd = (
        (
            true_S_residue
            + true_d_residue
        )
        * inv2
    ) % M

    q_from_true_sd = (
        (
            true_S_residue
            - true_d_residue
        )
        * inv2
    ) % M

    print(
        "TRUE STATE RECONSTRUCTION"
    )

    print()

    print(
        "true p mod M =",
        p % M
    )

    print(
        "from S,d      =",
        p_from_true_sd
    )

    print(
        "true q mod M =",
        q % M
    )

    print(
        "from S,d      =",
        q_from_true_sd
    )

    print()

    true_reconstruction = (
        p_from_true_sd == p % M
        and
        q_from_true_sd == q % M
    )

    print(
        "TRUE RECONSTRUCTION =",
        "PASS"
        if true_reconstruction
        else "FAIL"
    )

    print()

    # ------------------------------------------------------------
    # No integer factor search.
    # ------------------------------------------------------------

    print(
        "NO INTEGER FACTOR LIFT PERFORMED."
    )

    print(
        "This experiment measures only residue information."
    )

    print()

    # ------------------------------------------------------------
    # Runtime.
    # ------------------------------------------------------------

    elapsed = (
        time.perf_counter()
        - start
    )

    print(
        "TOTAL RUNTIME =",
        f"{elapsed:.6f} s"
    )

    print()

    print("-" * 72)
    print()


print("=" * 72)
print("FINISHED EXPERIMENT 149")
print("=" * 72)