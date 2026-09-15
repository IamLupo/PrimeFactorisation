#!/usr/bin/env python3

"""
START EXPERIMENT 183

FULL RESIDUE-CLASS REDUNDANCY TEST

Goal
----

Experiments 180-182 suggest that once we choose a residue class for p,
the corresponding residue class for q is completely determined by

    pq = n.

For any modulus M coprime to n:

    p ≡ P (mod M)

implies

    q ≡ n * P^{-1} (mod M).

Therefore the pair

    (P,Q)

does not contain two independent pieces of information.

This experiment makes that statement explicit.

For several CRT modulus products M, we enumerate every unit residue

    P in (Z/MZ)^*

and define

    Q = n * P^{-1} mod M.

Then we verify:

    P*Q ≡ n (mod M)

for every P.

We also verify that adding more radices does NOT reduce the number of
possible P classes beyond excluding P values that are non-units.

Finally we compare:

    number of admissible P residues modulo M
versus
    number of integers p <= sqrt(n).

The expected result is:

    #classes ≈ phi(M)

and every one has a compatible Q.

This means small-radix information alone is not selecting the
factor residue. It merely partitions the ordinary factor search into
CRT residue classes.

A second test verifies the same phenomenon with the full separate
R1 and R2 systems:

    P_R1 = p mod product(R1)
    Q_R1 = n*P_R1^-1 mod product(R1)

and independently

    Q_R2 = q mod product(R2)
    P_R2 = n*Q_R2^-1 mod product(R2).

The cross-system reconstruction is checked for arbitrary residue
choices as well as the true factor residues.

No C-values are used.

The hidden p,q are used only for benchmark validation.

"""


import math
import random
import time


# ============================================================
# CONFIG
# ============================================================

R1 = [17, 43, 59, 71, 83]
R2 = [19, 37, 61, 73, 89]

BITS = [30, 36, 42, 48, 54]

RANDOM_TESTS = 10

RANDOM_SEED = 183


# ============================================================
# NUMBER THEORY
# ============================================================

def is_prime(n):

    if n < 2:
        return False

    small = [
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37
    ]

    for p in small:

        if n == p:
            return True

        if n % p == 0:
            return False

    d = 41
    step = 2

    while d * d <= n:

        if n % d == 0:
            return False

        d += step
        step = 6 - step

    return True


def egcd(a, b):

    if b == 0:
        return a, 1, 0

    g, x1, y1 = egcd(
        b,
        a % b,
    )

    return (
        g,
        y1,
        x1 - (a // b) * y1,
    )


def inv_mod(a, m):

    g, x, _ = egcd(
        a,
        m,
    )

    if g != 1:
        raise ValueError(
            f"{a} is not invertible mod {m}; gcd={g}"
        )

    return x % m


# ============================================================
# CRT
# ============================================================

def crt_pair(
    a1,
    m1,
    a2,
    m2,
):

    g = math.gcd(
        m1,
        m2,
    )

    if (a2 - a1) % g != 0:
        return None

    m1r = m1 // g
    m2r = m2 // g

    rhs = (
        (a2 - a1)
        // g
    ) % m2r

    if m2r == 1:

        t = 0

    else:

        t = (
            rhs
            * inv_mod(
                m1r % m2r,
                m2r,
            )
        ) % m2r

    x = (
        a1
        + m1 * t
    )

    return (
        x % (m1 * m2r),
        m1 * m2r,
    )


def crt_vector(
    residues,
    moduli,
):

    x = 0
    M = 1

    for a, m in zip(
        residues,
        moduli,
    ):

        merged = crt_pair(
            x,
            M,
            a,
            m,
        )

        if merged is None:
            return None

        x, M = merged

    return x, M


# ============================================================
# SEMIPRIME GENERATION
# ============================================================

def make_semiprime(bits):

    low = 1 << (
        bits // 2 - 1
    )

    high = 1 << (
        bits // 2 + 1
    )

    for p in range(
        low | 1,
        high,
        2,
    ):

        if not is_prime(p):
            continue

        target = (
            (1 << bits)
            // p
        )

        for delta in range(
            -1000,
            1001,
            2,
        ):

            q = target + delta

            if q <= p:
                continue

            if not is_prime(q):
                continue

            n = p * q

            if n.bit_length() == bits:

                return p, q, n

    raise RuntimeError(
        "semiprime generation failed"
    )


# ============================================================
# EULER PHI
# ============================================================

def phi_product(
    primes,
):

    result = 1

    for p in primes:

        result *= (
            p - 1
        )

    return result


# ============================================================
# RESIDUE-PAIR TEST
# ============================================================

def test_modulus(
    n,
    primes,
):

    M = math.prod(
        primes
    )

    phi_M = phi_product(
        primes
    )

    unit_count = 0

    bad = 0

    seen_q = set()

    start = time.perf_counter()

    for P in range(
        1,
        M,
    ):

        if math.gcd(
            P,
            M,
        ) != 1:
            continue

        unit_count += 1

        Q = (
            n
            * inv_mod(
                P,
                M,
            )
        ) % M

        seen_q.add(
            Q
        )

        if (
            P * Q
        ) % M != n % M:

            bad += 1

    elapsed = (
        time.perf_counter()
        - start
    )

    return {
        "M": M,
        "phi": phi_M,
        "unit_count": unit_count,
        "q_count": len(seen_q),
        "bad": bad,
        "time": elapsed,
    }


# ============================================================
# DIRECT RANDOM RESIDUE TEST
# ============================================================

def random_residue_test(
    n,
    primes,
    rng,
):

    M = math.prod(
        primes
    )

    for _ in range(
        RANDOM_TESTS
    ):

        # Pick arbitrary unit P.
        while True:

            P = rng.randrange(
                1,
                M,
            )

            if math.gcd(
                P,
                M,
            ) == 1:

                break

        Q = (
            n
            * inv_mod(
                P,
                M,
            )
        ) % M

        lhs = (
            P * Q
        ) % M

        rhs = n % M

        print(
            f"    P={P} "
            f"Q={Q} "
            f"P*Q mod M={lhs} "
            f"n mod M={rhs} "
            f"ok={lhs == rhs}"
        )


# ============================================================
# FULL R1/R2 TEST
# ============================================================

def cross_system_test(
    n,
):

    M1 = math.prod(
        R1
    )

    M2 = math.prod(
        R2
    )

    total = M1 * M2

    print()
    print(
        "FULL R1/R2 MODULUS TEST"
    )

    print(
        f"    M1 = product(R1) = {M1}"
    )

    print(
        f"    M2 = product(R2) = {M2}"
    )

    print(
        f"    M1*M2 = {total}"
    )

    # --------------------------------------------------------
    # Choose arbitrary residue classes for p modulo M1.
    # --------------------------------------------------------

    for trial in range(
        RANDOM_TESTS
    ):

        P1 = random.randrange(
            1,
            M1,
        )

        while math.gcd(
            P1,
            M1,
        ) != 1:

            P1 = random.randrange(
                1,
                M1,
            )

        Q1 = (
            n
            * inv_mod(
                P1,
                M1,
            )
        ) % M1

        print()
        print(
            f"    RANDOM R1 TRIAL "
            f"{trial + 1}"
        )

        print(
            f"        P mod M1 = "
            f"{P1}"
        )

        print(
            f"        Q mod M1 = "
            f"{Q1}"
        )

        print(
            f"        PQ=n mod M1 = "
            f"{(P1 * Q1) % M1 == n % M1}"
        )

    # --------------------------------------------------------
    # Do the same on R2.
    # --------------------------------------------------------

    for trial in range(
        RANDOM_TESTS
    ):

        Q2 = random.randrange(
            1,
            M2,
        )

        while math.gcd(
            Q2,
            M2,
        ) != 1:

            Q2 = random.randrange(
                1,
                M2,
            )

        P2 = (
            n
            * inv_mod(
                Q2,
                M2,
            )
        ) % M2

        print()
        print(
            f"    RANDOM R2 TRIAL "
            f"{trial + 1}"
        )

        print(
            f"        Q mod M2 = "
            f"{Q2}"
        )

        print(
            f"        P mod M2 = "
            f"{P2}"
        )

        print(
            f"        PQ=n mod M2 = "
            f"{(P2 * Q2) % M2 == n % M2}"
        )


# ============================================================
# RUN INSTANCE
# ============================================================

def run_instance(bits):

    print()
    print("=" * 72)
    print(
        f"START INSTANCE {bits}-BIT"
    )
    print("=" * 72)

    p, q, n = make_semiprime(
        bits
    )

    sqrt_n = math.isqrt(n)

    print(
        f"bits(n) = {n.bit_length()}"
    )

    print(
        f"sqrt(n) = {sqrt_n}"
    )

    print(
        f"hidden p = {p}"
    )

    print(
        f"hidden q = {q}"
    )

    # --------------------------------------------------------
    # Moduli.
    # --------------------------------------------------------

    modulus_sets = [
        [17],
        [17, 43],
        [17, 43, 59],
        [17, 43, 59, 71],
    ]

    print()
    print(
        "P-RESIDUE FREEDOM"
    )

    for primes in modulus_sets:

        result = test_modulus(
            n,
            primes,
        )

        M = result["M"]

        sqrt_density = (
            sqrt_n / M
        )

        print()
        print(
            f"  primes={primes}"
        )

        print(
            f"    M = "
            f"{M}"
        )

        print(
            f"    phi(M) = "
            f"{result['phi']}"
        )

        print(
            f"    enumerated units = "
            f"{result['unit_count']}"
        )

        print(
            f"    distinct Q values = "
            f"{result['q_count']}"
        )

        print(
            f"    identity failures = "
            f"{result['bad']}"
        )

        print(
            f"    sqrt(n)/M = "
            f"{sqrt_density:.3f}"
        )

        print(
            f"    time = "
            f"{result['time']:.6f}s"
        )

        print(
            f"    all units have compatible Q = "
            f"{result['bad'] == 0}"
        )

    # --------------------------------------------------------
    # Random direct verification.
    # --------------------------------------------------------

    rng = random.Random(
        RANDOM_SEED + bits
    )

    print()
    print(
        "RANDOM UNIT TESTS"
    )

    random_residue_test(
        n,
        [17, 43, 59],
        rng,
    )

    # --------------------------------------------------------
    # True residue class.
    # --------------------------------------------------------

    print()
    print(
        "TRUE FACTOR RESIDUE"
    )

    M = math.prod(
        [17, 43, 59]
    )

    P_true = p % M

    Q_true = q % M

    Q_from_n = (
        n
        * inv_mod(
            P_true,
            M,
        )
    ) % M

    print(
        f"    M = {M}"
    )

    print(
        f"    true P mod M = "
        f"{P_true}"
    )

    print(
        f"    true Q mod M = "
        f"{Q_true}"
    )

    print(
        f"    derived Q from n,P = "
        f"{Q_from_n}"
    )

    print(
        f"    exact match = "
        f"{Q_true == Q_from_n}"
    )

    # --------------------------------------------------------
    # Full R1/R2 system.
    # --------------------------------------------------------

    cross_system_test(
        n
    )

    # --------------------------------------------------------
    # CRT reconstruction of arbitrary vectors.
    # --------------------------------------------------------

    print()
    print(
        "ARBITRARY VECTOR CRT TEST"
    )

    arbitrary_A = [
        1,
        2,
        3,
        4,
        5,
    ]

    arbitrary_B = [
        1,
        2,
        3,
        4,
        5,
    ]

    A_result = crt_vector(
        arbitrary_A,
        R1,
    )

    B_result = crt_vector(
        arbitrary_B,
        R2,
    )

    if (
        A_result is not None
        and B_result is not None
    ):

        P0, M1 = A_result

        Q0, M2 = B_result

        print(
            f"    arbitrary A = "
            f"{arbitrary_A}"
        )

        print(
            f"    arbitrary B = "
            f"{arbitrary_B}"
        )

        print(
            f"    P0 = {P0}"
        )

        print(
            f"    M1 = {M1}"
        )

        print(
            f"    Q0 = {Q0}"
        )

        print(
            f"    M2 = {M2}"
        )

        # Check the factor equation modulo M1 and M2
        # using the induced opposite-side residues.

        Q_on_M1 = (
            n
            * inv_mod(
                P0,
                M1,
            )
        ) % M1

        P_on_M2 = (
            n
            * inv_mod(
                Q0,
                M2,
            )
        ) % M2

        print(
            f"    induced Q mod M1 = "
            f"{Q_on_M1}"
        )

        print(
            f"    induced P mod M2 = "
            f"{P_on_M2}"
        )

        print(
            f"    arbitrary P valid = "
            f"{math.gcd(P0, M1) == 1}"
        )

        print(
            f"    arbitrary Q valid = "
            f"{math.gcd(Q0, M2) == 1}"
        )

    print()
    print(
        "INTERPRETATION"
    )

    print(
        "    Every unit residue P mod M has exactly "
        "one compatible Q mod M."
    )

    print(
        "    Therefore pq=n does not select P mod M."
    )

    print(
        "    Adding more small radices only increases "
        "the modulus; it does not reduce the number "
        "of admissible residue classes."
    )

    print(
        "    The residue method becomes useful only "
        "when an independent source selects the "
        "correct residue class."
    )

    print()
    print(
        f"FINISHED INSTANCE {bits}-BIT"
    )


# ============================================================
# MAIN
# ============================================================

def main():

    print()
    print(
        "START EXPERIMENT 183"
    )

    print()

    print(
        f"R1 = {R1}"
    )

    print(
        f"R2 = {R2}"
    )

    print(
        f"BITS = {BITS}"
    )

    print()

    for bits in BITS:

        run_instance(
            bits
        )

    print()
    print(
        "FINISHED EXPERIMENT 183"
    )


if __name__ == "__main__":
    main()
