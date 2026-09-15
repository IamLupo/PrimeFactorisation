#!/usr/bin/env python3

import math
import random
import time


# ============================================================
# START EXPERIMENT 107
# Direct p-space search
# ============================================================

R1 = [17, 43]
R2 = [19, 47]

R10, R11 = R1
R20, R21 = R2

R00 = R10 * R20

RANDOM_SEED = 107


# ============================================================
# PRIME TEST
# ============================================================

def is_probable_prime(n):

    if n < 2:
        return False

    small = (
        2, 3, 5, 7, 11, 13, 17,
        19, 23, 29, 31, 37
    )

    for p in small:
        if n % p == 0:
            return n == p

    d = n - 1
    s = 0

    while d % 2 == 0:
        d //= 2
        s += 1

    for a in (
        2, 3, 5, 7,
        11, 13, 17
    ):

        if a >= n:
            continue

        x = pow(a, d, n)

        if x == 1 or x == n - 1:
            continue

        for _ in range(s - 1):

            x = (x * x) % n

            if x == n - 1:
                break

        else:
            return False

    return True


def random_prime(bits):

    while True:

        p = random.getrandbits(bits)

        p |= (1 << (bits - 1))
        p |= 1

        if is_probable_prime(p):
            return p


def make_semiprime(bits):

    factor_bits = bits // 2

    while True:

        p = random_prime(factor_bits)
        q = random_prime(factor_bits)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        n = p * q

        if n.bit_length() >= bits - 1:
            return n, p, q


# ============================================================
# CARRY
# ============================================================

def carries(r1, r2, k, l, a, b):

    c1 = (k * b) // r2
    c2 = (l * a) // r1

    beta = (k * b) % r2
    alpha = (l * a) % r1

    c3 = (
        r1 * beta
        + r2 * alpha
        + a * b
    ) // (r1 * r2)

    return c1, c2, c3, c1 + c2 + c3


# ============================================================
# FIRST-CELL SOLVER
# ============================================================

def solve_l0_fixed_t(
    Q,
    r1,
    r2,
    k,
    a,
    b,
    t,
):
    """
    Solve the first-cell carry equation with

        l = t + r1*m

    for one fixed residue class t = l mod r1.

    Returns l or None.
    """

    c1 = (k * b) // r2
    beta = (k * b) % r2

    alpha = (t * a) % r1

    c3 = (
        r1 * beta
        + r2 * alpha
        + a * b
    ) // (r1 * r2)

    floor_t = (t * a) // r1

    numerator = (
        Q
        - k * t
        - c1
        - floor_t
        - c3
    )

    denominator = k * r1 + a

    if numerator < 0:
        return None

    if numerator % denominator != 0:
        return None

    m = numerator // denominator

    if m < 0:
        return None

    l = t + r1 * m

    if l <= 0:
        return None

    # Exact verification.

    _, _, _, Ecarry = carries(
        r1,
        r2,
        k,
        l,
        a,
        b,
    )

    Equot = Q - k * l

    if Ecarry != Equot:
        return None

    return l


# ============================================================
# MODULAR q FORCING
# ============================================================

def forced_q_residue(n, p, modulus):

    g = math.gcd(p, modulus)

    if g != 1:
        return None

    inv = pow(
        p,
        -1,
        modulus,
    )

    return (
        n * inv
    ) % modulus


def decode_q_residue(q0):

    # q = 19*l + b
    #
    # q0 = 19*t + b
    #
    # where t = l mod 17.

    b = q0 % R20
    t = q0 // R20

    return t, b


# ============================================================
# SEARCH
# ============================================================

def search(n):

    Q00 = n // R00

    p_limit = math.isqrt(n)

    counters = {
        "p_tested": 0,
        "p_modular_forced": 0,
        "p_modular_failed": 0,
        "first_cell_candidates": 0,
        "bucket_survivors": 0,
        "mod_R00_survivors": 0,
        "exact": 0,
    }

    bucket_survivors = []
    congruence_survivors = []

    start = time.perf_counter()

    # --------------------------------------------------------
    # Since n is odd and is a product of two odd primes,
    # only odd p need to be examined.
    #
    # The factor 2 case is handled separately.
    # --------------------------------------------------------

    p_start = 3

    if n % 2 == 0:

        p_start = 2

    else:

        # ----------------------------------------------------
        # Special small-prime cases.
        # ----------------------------------------------------

        for small in (17, 19):

            if small <= p_limit and n % small == 0:

                q = n // small

                if small <= q:

                    bucket_survivors.append(
                        (small, q, 0)
                    )

                    congruence_survivors.append(
                        (small, q, 0)
                    )

        if p_start == 3:
            pass

    # --------------------------------------------------------
    # Direct p scan.
    # --------------------------------------------------------

    if p_start == 2 and n % 2 == 0:

        p_range = range(
            2,
            p_limit + 1
        )

    else:

        p_range = range(
            3,
            p_limit + 1,
            2
        )

    for p in p_range:

        counters["p_tested"] += 1

        # Avoid the special primes whose inverses modulo R00
        # do not exist. They were handled above.
        if math.gcd(p, R00) != 1:

            counters["p_modular_failed"] += 1
            continue

        # ----------------------------------------------------
        # p directly determines k0 and a0.
        # ----------------------------------------------------

        k0 = p // R10
        a0 = p % R10

        # ----------------------------------------------------
        # Modularly force q mod R00.
        # ----------------------------------------------------

        q0 = forced_q_residue(
            n,
            p,
            R00,
        )

        if q0 is None:

            counters["p_modular_failed"] += 1
            continue

        counters["p_modular_forced"] += 1

        # q0 gives:
        #
        #   b0 = q mod 19
        #   t  = l0 mod 17
        #

        t, b0 = decode_q_residue(q0)

        # ----------------------------------------------------
        # First-cell quotient.
        # ----------------------------------------------------

        l0 = solve_l0_fixed_t(
            Q00,
            R10,
            R20,
            k0,
            a0,
            b0,
            t,
        )

        if l0 is None:
            continue

        counters["first_cell_candidates"] += 1

        # ----------------------------------------------------
        # Construct q.
        # ----------------------------------------------------

        q = R20 * l0 + b0

        if p > q:
            continue

        product = p * q

        # ----------------------------------------------------
        # First-cell bucket.
        # ----------------------------------------------------

        if product // R00 != n // R00:
            continue

        counters["bucket_survivors"] += 1

        error = product - n

        bucket_survivors.append(
            (p, q, error)
        )

        # ----------------------------------------------------
        # Congruence closure.
        # ----------------------------------------------------

        if product % R00 != n % R00:
            continue

        counters["mod_R00_survivors"] += 1

        congruence_survivors.append(
            (p, q, error)
        )

        if product == n:
            counters["exact"] += 1

    counters["runtime"] = (
        time.perf_counter() - start
    )

    return (
        counters,
        bucket_survivors,
        congruence_survivors,
    )


# ============================================================
# REPORT
# ============================================================

def report(
    n,
    true_p,
    true_q,
    counters,
    bucket,
    congruence,
):

    print()
    print("-" * 72)

    print(
        "n bits =",
        n.bit_length()
    )

    print()
    print(
        "n      =",
        n
    )

    print(
        "true p =",
        true_p
    )

    print(
        "true q =",
        true_q
    )

    print()
    print("GRID")

    print(
        "r1 =",
        R1
    )

    print(
        "r2 =",
        R2
    )

    print(
        "R00 =",
        R00
    )

    print()
    print("SEARCH RESULTS")

    for key, value in counters.items():

        if key == "runtime":

            print(
                f"{key:<24} = "
                f"{value:.6f} s"
            )

        else:

            print(
                f"{key:<24} = "
                f"{value}"
            )

    # --------------------------------------------------------
    # Error bound.
    # --------------------------------------------------------

    if bucket:

        max_error = max(
            abs(error)
            for _, _, error in bucket
        )

        print()
        print("BUCKET ERROR")

        print(
            "max |pq-n| =",
            max_error
        )

        print(
            "R00 - 1    =",
            R00 - 1
        )

    # --------------------------------------------------------
    # Survivors.
    # --------------------------------------------------------

    print()
    print("CONGRUENCE SURVIVORS")

    for p, q, error in congruence[:20]:

        print(
            "  "
            f"p={p} "
            f"q={q} "
            f"pq-n={error}"
        )

    if len(congruence) > 20:

        print(
            "  ..."
            f"{len(congruence) - 20}"
            " more"
        )

    print()
    print("FINAL VERIFICATION")

    exact = [
        (p, q)
        for p, q, error in congruence
        if p * q == n
    ]

    print(
        "exact factorizations =",
        len(exact)
    )

    print(
        "true factorization found =",
        (true_p, true_q) in exact
    )


# ============================================================
# MAIN
# ============================================================

def main():

    random.seed(RANDOM_SEED)

    print("=" * 72)
    print("START EXPERIMENT 107")
    print("Direct p-space search")
    print("=" * 72)

    print()
    print("r1 =", R1)
    print("r2 =", R2)
    print("R00 =", R00)

    for bits in (30, 36, 42):

        print()
        print("=" * 72)
        print(
            f"GENERATING {bits}-BIT SEMIPRIME"
        )
        print("=" * 72)

        n, p, q = make_semiprime(bits)

        (
            counters,
            bucket,
            congruence,
        ) = search(n)

        report(
            n,
            p,
            q,
            counters,
            bucket,
            congruence,
        )

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 107")
    print("=" * 72)


if __name__ == "__main__":
    main()
