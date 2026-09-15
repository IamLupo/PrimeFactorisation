#!/usr/bin/env python3

import math
import random
import time


# ============================================================
# START EXPERIMENT 106
# Early modular forcing
# ============================================================

R1 = [17, 43]
R2 = [19, 47]

R00 = R1[0] * R2[0]

RANDOM_SEED = 106


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
# FIXED-RESIDUE l SOLVER
# ============================================================

def solve_l0_fixed_t(Q, r1, r2, k, a, b, t):

    """
    Solve

        Q - k*l = carry_E(k,l,a,b)

    but with

        l = t + r1*m

    for one specific residue class t.

    This avoids scanning all t.
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
        b
    )

    Equot = Q - k * l

    if Ecarry != Equot:
        return None

    return l


# ============================================================
# MODULAR FORCING
# ============================================================

def forced_q_residue(n, p, modulus):

    """
    Solve

        p*q = n (mod modulus)

    for q.

    Returns q0 in [0, modulus-1], or None if impossible.
    """

    g = math.gcd(p, modulus)

    if n % g != 0:
        return None

    if g != 1:
        return None

    p_inv = pow(
        p,
        -1,
        modulus
    )

    return (n * p_inv) % modulus


def decode_q_residue(q0, r2, r1):

    """
    q = r2*l + b

    with

        0 <= b < r2

    and q known modulo r1*r2.

    Write

        q0 = r2*t + b

    where

        0 <= t < r1.

    Then t = l mod r1.
    """

    b = q0 % r2

    t = q0 // r2

    return t, b


# ============================================================
# SEARCH
# ============================================================

def search(n):

    r10, r11 = R1
    r20, r21 = R2

    Q00 = n // R00

    p_limit = math.isqrt(n)

    k0_max = p_limit // r10

    counters = {
        "k0_tested": 0,

        "a_tests": 0,

        "modular_forced": 0,
        "modular_failed": 0,

        "l_candidates": 0,

        "bucket_survivors": 0,
        "mod_R00_survivors": 0,
    }

    bucket = []
    congruence = []

    start = time.perf_counter()

    for k0 in range(1, k0_max + 1):

        counters["k0_tested"] += 1

        for a0 in range(r10):

            counters["a_tests"] += 1

            # ------------------------------------------------
            # p is now known.
            # ------------------------------------------------

            p = r10 * k0 + a0

            # ------------------------------------------------
            # Since we ultimately want p*q=n mod R00,
            # force q modulo R00.
            #
            # For a large prime p, gcd(p,R00)=1 almost
            # always. Handle the exceptional case explicitly.
            # ------------------------------------------------

            q0 = forced_q_residue(
                n,
                p,
                R00
            )

            if q0 is None:

                counters["modular_failed"] += 1

                continue

            counters["modular_forced"] += 1

            # ------------------------------------------------
            # q0 = r20*t + b0
            #
            # Thus b0 and l0 mod r10 are now FORCED.
            # ------------------------------------------------

            t, b0 = decode_q_residue(
                q0,
                r20,
                r10
            )

            # ------------------------------------------------
            # Solve the first-cell carry equation only in
            # this forced residue class.
            # ------------------------------------------------

            l0 = solve_l0_fixed_t(
                Q00,
                r10,
                r20,
                k0,
                a0,
                b0,
                t
            )

            if l0 is None:
                continue

            counters["l_candidates"] += 1

            # ------------------------------------------------
            # Construct q.
            # ------------------------------------------------

            q = r20 * l0 + b0

            if p > q:
                continue

            # ------------------------------------------------
            # First-cell bucket.
            # ------------------------------------------------

            product = p * q

            if product // R00 != n // R00:
                continue

            counters["bucket_survivors"] += 1

            error = product - n

            bucket.append(
                (p, q, error)
            )

            # ------------------------------------------------
            # Congruence closure.
            #
            # This should now be automatic because q was
            # forced modulo R00, but explicitly verify it.
            # ------------------------------------------------

            if product % R00 != n % R00:
                continue

            counters["mod_R00_survivors"] += 1

            congruence.append(
                (p, q, error)
            )

    counters["runtime"] = (
        time.perf_counter() - start
    )

    return counters, bucket, congruence


# ============================================================
# REPORT
# ============================================================

def report(
    n,
    true_p,
    true_q,
    counters,
    bucket,
    congruence
):

    print()
    print("-" * 72)

    print("n bits =", n.bit_length())

    print()
    print("n      =", n)
    print("true p =", true_p)
    print("true q =", true_q)

    print()
    print("GRID")

    print("r1  =", R1)
    print("r2  =", R2)

    print()
    print("R00 =", R00)

    print()
    print("SEARCH RESULTS")

    print(
        "k0 tested            =",
        counters["k0_tested"]
    )

    print(
        "a tests              =",
        counters["a_tests"]
    )

    print(
        "modular forced       =",
        counters["modular_forced"]
    )

    print(
        "modular failed       =",
        counters["modular_failed"]
    )

    print(
        "l candidates         =",
        counters["l_candidates"]
    )

    print(
        "bucket survivors     =",
        counters["bucket_survivors"]
    )

    print(
        "mod-R00 survivors    =",
        counters["mod_R00_survivors"]
    )

    print(
        "runtime              = %.6f s"
        % counters["runtime"]
    )

    # --------------------------------------------------------
    # Bucket error.
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
            "R00-1      =",
            R00 - 1
        )

    # --------------------------------------------------------
    # Congruence candidates.
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
            f"{len(congruence)-20}"
            " more"
        )

    # --------------------------------------------------------
    # Final verification.
    # --------------------------------------------------------

    exact = [
        (p, q)
        for p, q, error in congruence
        if p * q == n
    ]

    print()
    print("FINAL VERIFICATION")

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
    print("START EXPERIMENT 106")
    print("Early modular forcing")
    print("=" * 72)

    print()
    print("r1 =", R1)
    print("r2 =", R2)
    print("R00 =", R00)

    for bits in (30, 36, 42):

        print()
        print("=" * 72)
        print(f"GENERATING {bits}-BIT SEMIPRIME")
        print("=" * 72)

        n, p, q = make_semiprime(bits)

        counters, bucket, congruence = search(n)

        report(
            n,
            p,
            q,
            counters,
            bucket,
            congruence
        )

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 106")
    print("=" * 72)


if __name__ == "__main__":
    main()

