#!/usr/bin/env python3

import math
import random
import time


# ========================================================================
# START EXPERIMENT 126
# Direct Diophantine recovery of row quotients
# ========================================================================

print("=" * 72)
print("START EXPERIMENT 126")
print("Direct Diophantine recovery of row quotients")
print("=" * 72)


# ------------------------------------------------------------------------
# GRID
# ------------------------------------------------------------------------

R1 = [17, 43, 47]
R2 = [19, 37, 53]

print()
print("R1 =", R1)
print("R2 =", R2)


# ------------------------------------------------------------------------
# MILLER-RABIN
# ------------------------------------------------------------------------

def is_probable_prime(n):

    if n < 2:
        return False

    small = (
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37
    )

    for p in small:

        if n == p:
            return True

        if n % p == 0:
            return False

    d = n - 1
    s = 0

    while d % 2 == 0:
        d //= 2
        s += 1

    bases = (
        2,
        325,
        9375,
        28178,
        450775,
        9780504,
        1795265022,
    )

    for a in bases:

        if a % n == 0:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        for _ in range(s - 1):

            x = x * x % n

            if x == n - 1:
                break

        else:
            return False

    return True


# ------------------------------------------------------------------------
# PRIME
# ------------------------------------------------------------------------

def random_prime(bits):

    while True:

        x = random.getrandbits(bits)

        x |= 1
        x |= 1 << (bits - 1)

        if is_probable_prime(x):
            return x


# ------------------------------------------------------------------------
# SEMIPRIME
# ------------------------------------------------------------------------

def random_semiprime(bits):

    pb = bits // 2
    qb = bits - pb

    while True:

        p = random_prime(pb)
        q = random_prime(qb)

        if p == q:
            continue

        n = p * q

        if n.bit_length() == bits:
            return n, min(p, q), max(p, q)


# ------------------------------------------------------------------------
# CARRY STRUCTURE
# ------------------------------------------------------------------------

def carry_state(p, q, r1, r2):

    k, a = divmod(p, r1)
    ell, b = divmod(q, r2)

    c1 = (k * b) // r2
    c2 = (ell * a) // r1

    beta = (k * b) % r2
    alpha = (ell * a) % r1

    c3 = (
        r1 * beta
        + r2 * alpha
        + a * b
    ) // (r1 * r2)

    E = c1 + c2 + c3

    Q = (p * q) // (r1 * r2)

    return {
        "k": k,
        "ell": ell,
        "a": a,
        "b": b,
        "E": E,
        "Q": Q,
        "K": Q - E,
    }


# ------------------------------------------------------------------------
# DIRECT CONTROL
# ------------------------------------------------------------------------

def direct_factor(n):

    limit = math.isqrt(n)

    tested = 0

    for p in range(2, limit + 1):

        tested += 1

        if n % p == 0:
            return p, n // p, tested

    return None, None, tested


# ------------------------------------------------------------------------
# NEW ALGEBRA
#
# For a fixed column j:
#
#   K0 = k0 * ell
#   K1 = k1 * ell
#
# Therefore:
#
#   K0 * k1 = K1 * k0
#
# and:
#
#   r0*k0 - r1*k1 = a1-a0
#
# Substitute:
#
#   k0 = K0/d
#   k1 = K1/d
#
# where d = ell.
#
# Instead of guessing E, use the fact that
#
#   K0 = Q0-E0
#   K1 = Q1-E1.
#
# The experiment tests whether the small residue difference
#
#   delta = a1-a0
#
# plus divisibility constraints can directly recover ell.
# ------------------------------------------------------------------------

def recover_from_column(
    n,
    r0,
    r1,
    r2,
    Q0,
    Q1,
    a0,
    a1,
):
    """
    Attempt to recover ell from the exact row relation.

    We know:

        r0*k0 + a0 = r1*k1 + a1

    hence:

        r0*k0 - r1*k1 = a1-a0.

    Also:

        K0 = k0*ell
        K1 = k1*ell.

    Thus:

        r0*K0 - r1*K1
            = ell*(a1-a0).

    We search only over the possible value of
    ell inferred from the difference, using the fact
    that E values are themselves sums of carry components.
    """

    delta = a1 - a0

    if delta == 0:
        return []

    results = []

    # ------------------------------------------------------------
    # The E values are not arbitrary. For this experiment we
    # derive candidate E pairs from the requirement that
    #
    #   r0(Q0-E0)-r1(Q1-E1)
    #
    # is a multiple of delta.
    #
    # Instead of scanning a full E x E grid, use a single E0
    # and solve the corresponding E1 modulo |delta|.
    #
    # We cap E using a bound based on sqrt(n).
    # ------------------------------------------------------------

    s = math.isqrt(n)

    max_k = s // r0 + 1
    max_l = s // r2 + 1

    EMAX = max_k + max_l + 2

    step = abs(delta)

    for E0 in range(EMAX + 1):

        K0 = Q0 - E0

        if K0 <= 0:
            continue

        base = r0 * K0 - r1 * Q1

        # Need:
        #
        #   base + r1*E1
        #
        # divisible by delta.
        #
        # Solve:
        #
        #   r1*E1 == -base (mod delta)
        #
        if math.gcd(r1, step) != 1:

            # General modular case.
            g = math.gcd(r1, step)

            if (-base) % g != 0:
                continue

            rr = r1 // g
            mm = step // g
            bb = (-base // g) % mm

            if mm == 1:
                E1_first = 0
            else:
                inv = pow(rr, -1, mm)
                E1_first = (bb * inv) % mm

            E1_step = mm

        else:

            inv = pow(
                r1,
                -1,
                step
            )

            E1_first = (
                (-base)
                * inv
            ) % step

            E1_step = step

        E1 = E1_first

        while E1 <= EMAX:

            K1 = Q1 - E1

            if K1 > 0:

                numerator = (
                    r0 * K0
                    - r1 * K1
                )

                if numerator % delta == 0:

                    ell = numerator // delta

                    if ell > 0:

                        if K0 % ell == 0 and K1 % ell == 0:

                            k0 = K0 // ell
                            k1 = K1 // ell

                            p0 = r0 * k0 + a0
                            p1 = r1 * k1 + a1

                            if p0 == p1:

                                p = p0

                                if 1 < p < n:

                                    if n % p == 0:

                                        q = n // p

                                        results.append(
                                            (
                                                p,
                                                q,
                                                ell,
                                                k0,
                                                k1,
                                                E0,
                                                E1,
                                            )
                                        )

            E1 += E1_step

    return results


# ------------------------------------------------------------------------
# TEST CASE
# ------------------------------------------------------------------------

random.seed(126)

TEST_BITS = [
    30,
    36,
    42,
]

for bits in TEST_BITS:

    print()
    print("=" * 72)
    print("GENERATING", bits, "-BIT SEMIPRIME")
    print("=" * 72)

    n, true_p, true_q = random_semiprime(bits)

    print()
    print("-" * 72)
    print("n bits =", n.bit_length())
    print("n      =", n)
    print("true p =", true_p)
    print("true q =", true_q)

    # ------------------------------------------------------------
    # Direct control
    # ------------------------------------------------------------

    t0 = time.perf_counter()

    pd, qd, tested = direct_factor(n)

    t1 = time.perf_counter()

    direct_time = t1 - t0

    print()
    print("DIRECT SQRT(n) SCAN")
    print("p tested =", tested)
    print("runtime  =", f"{direct_time:.6f}", "s")

    # ------------------------------------------------------------
    # Use all columns.
    # ------------------------------------------------------------

    print()
    print("DIophantine CARRY RECOVERY")

    t0 = time.perf_counter()

    all_candidates = []

    for r2 in R2:

        Q0 = n // (R1[0] * r2)
        Q1 = n // (R1[1] * r2)

        # --------------------------------------------------------
        # Enumerate only the small residues a0,a1.
        # --------------------------------------------------------

        for a0 in range(R1[0]):

            for a1 in range(R1[1]):

                candidates = recover_from_column(
                    n,
                    R1[0],
                    R1[1],
                    r2,
                    Q0,
                    Q1,
                    a0,
                    a1,
                )

                all_candidates.extend(candidates)

    t1 = time.perf_counter()

    carry_time = t1 - t0

    # ------------------------------------------------------------
    # Deduplicate
    # ------------------------------------------------------------

    unique = {}

    for candidate in all_candidates:

        key = (
            candidate[0],
            candidate[1],
        )

        unique[key] = candidate

    all_candidates = list(unique.values())

    recovered = any(
        {x[0], x[1]}
        == {true_p, true_q}
        for x in all_candidates
    )

    print()
    print("candidate recoveries =", len(all_candidates))
    print("runtime =", f"{carry_time:.6f}", "s")
    print("TRUE FACTORIZATION FOUND =", recovered)

    # ------------------------------------------------------------
    # Print candidates
    # ------------------------------------------------------------

    if all_candidates:

        print()
        print("RECOVERED CANDIDATES")

        for x in all_candidates[:20]:

            p, q, ell, k0, k1, E0, E1 = x

            exact = (
                {p, q}
                == {true_p, true_q}
            )

            print(
                "p =", p,
                "q =", q,
                "ell =", ell,
                "k0 =", k0,
                "k1 =", k1,
                "E0 =", E0,
                "E1 =", E1,
                "EXACT =", exact,
            )

    # ------------------------------------------------------------
    # Compare
    # ------------------------------------------------------------

    print()
    print("TIMING")

    if carry_time > 0:

        print(
            "direct / carry =",
            f"{direct_time / carry_time:.3f}x"
        )

    # ------------------------------------------------------------
    # Full verification
    # ------------------------------------------------------------

    if recovered:

        print()
        print("FULL GRID VERIFICATION")

        x = next(
            x
            for x in all_candidates
            if {x[0], x[1]}
            == {true_p, true_q}
        )

        p, q = x[:2]

        verified = True

        for r1 in R1:

            for r2 in R2:

                c = carry_state(
                    p,
                    q,
                    r1,
                    r2,
                )

                if (
                    c["Q"] - c["E"]
                    != c["k"] * c["ell"]
                ):
                    verified = False

        print(
            "carry identity verified =",
            verified
        )


# ========================================================================
# FINISHED EXPERIMENT 126
# ========================================================================

print()
print("=" * 72)
print("FINISHED EXPERIMENT 126")
print("=" * 72)
