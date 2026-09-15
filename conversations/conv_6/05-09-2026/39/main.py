#!/usr/bin/env python3

import math
import random
import time


# ============================================================
# START EXPERIMENT 122
# Reciprocal two-sided interval search
# ============================================================

R00 = 15
R10 = 145
R01 = 93
R11 = 899

M = math.lcm(R00, R11)

RANDOM_SEED = 122


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
        2, 3, 5, 7, 11, 13, 17
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

        p |= 1 << (bits - 1)
        p |= 1

        if is_probable_prime(p):
            return p


def make_semiprime(bits):

    fb = bits // 2

    while True:

        p = random_prime(fb)
        q = random_prime(fb)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        n = p * q

        if n.bit_length() >= bits - 1:
            return n, p, q


# ============================================================
# QUOTIENT BUCKET
# ============================================================

def partner_interval(n, x, R):

    """
    If x is one factor and y is the other, and

        floor(x*y/R) = floor(n/R),

    then

        R*Q <= x*y < R*(Q+1).

    Hence

        ceil(R*Q/x) <= y
        <= floor((R*(Q+1)-1)/x).
    """

    Q = n // R

    low = (
        R * Q + x - 1
    ) // x

    high = (
        R * (Q + 1) - 1
    ) // x

    return low, high


# ============================================================
# MODULAR PARTNER
# ============================================================

def partner_residue(n, x, R):

    g = math.gcd(x, R)

    if n % g != 0:
        return None

    solutions = solve_linear_congruence(
        x,
        n,
        R
    )

    return solutions


def solve_linear_congruence(a, b, m):

    g = math.gcd(a, m)

    if b % g != 0:
        return []

    aa = a // g
    bb = b // g
    mm = m // g

    if mm == 1:
        return [0]

    inv = pow(
        aa % mm,
        -1,
        mm
    )

    x0 = (
        bb * inv
    ) % mm

    return [
        x0 + t * mm
        for t in range(g)
    ]


# ============================================================
# INTERVAL + RESIDUE INTERSECTION
# ============================================================

def first_in_interval_with_residue(
    low,
    high,
    residue,
    modulus
):

    if low > high:
        return None

    delta = (
        residue - low
    ) % modulus

    x = low + delta

    if x > high:
        return None

    return x


# ============================================================
# P-DIRECTION
# ============================================================

def search_from_p(n):

    limit = math.isqrt(n)

    p_tested = 0
    q_candidates = 0
    exact = []

    start = time.perf_counter()

    for p in range(
        3,
        limit + 1,
        2
    ):

        p_tested += 1

        # ----------------------------------------------------
        # q modulo R00
        # ----------------------------------------------------

        q00_states = partner_residue(
            n,
            p,
            R00
        )

        if not q00_states:
            continue

        # ----------------------------------------------------
        # q modulo R11
        # ----------------------------------------------------

        q11_states = partner_residue(
            n,
            p,
            R11
        )

        if not q11_states:
            continue

        # ----------------------------------------------------
        # Second quotient interval.
        # ----------------------------------------------------

        q_low, q_high = partner_interval(
            n,
            p,
            R11
        )

        for q00 in q00_states:

            for q11 in q11_states:

                # ------------------------------------------------
                # CRT:
                #
                # q == q00 mod 15
                # q == q11 mod 899
                # ------------------------------------------------

                combined = crt(
                    q00,
                    R00,
                    q11,
                    R11
                )

                if combined is None:
                    continue

                q_res, modulus = combined

                q = first_in_interval_with_residue(
                    q_low,
                    q_high,
                    q_res,
                    modulus
                )

                while q is not None and q <= q_high:

                    q_candidates += 1

                    if p <= q:

                        if p * q == n:

                            exact.append(
                                (p, q)
                            )

                            return {
                                "p_tested": p_tested,
                                "q_candidates":
                                    q_candidates,
                                "exact": exact,
                                "time":
                                    time.perf_counter()
                                    - start,
                            }

                    q += modulus

                    if q > q_high:
                        break

    return {
        "p_tested": p_tested,
        "q_candidates":
            q_candidates,
        "exact": exact,
        "time":
            time.perf_counter()
            - start,
    }


# ============================================================
# CRT
# ============================================================

def crt(a1, m1, a2, m2):

    g = math.gcd(m1, m2)

    if (a2 - a1) % g != 0:
        return None

    m1r = m1 // g
    m2r = m2 // g

    rhs = (
        a2 - a1
    ) // g

    if m2r == 1:
        t = 0
    else:

        inv = pow(
            m1r % m2r,
            -1,
            m2r
        )

        t = (
            rhs * inv
        ) % m2r

    x = (
        a1
        + m1 * t
    )

    modulus = (
        m1 * m2r
    )

    return (
        x % modulus,
        modulus
    )


# ============================================================
# Q-DIRECTION
# ============================================================

def search_from_q(n):

    limit = math.isqrt(n)

    q_tested = 0
    p_candidates = 0
    exact = []

    start = time.perf_counter()

    for q in range(
        3,
        limit + 1,
        2
    ):

        q_tested += 1

        p00_states = partner_residue(
            n,
            q,
            R00
        )

        if not p00_states:
            continue

        p11_states = partner_residue(
            n,
            q,
            R11
        )

        if not p11_states:
            continue

        p_low, p_high = partner_interval(
            n,
            q,
            R11
        )

        for p00 in p00_states:

            for p11 in p11_states:

                combined = crt(
                    p00,
                    R00,
                    p11,
                    R11
                )

                if combined is None:
                    continue

                p_res, modulus = combined

                p = first_in_interval_with_residue(
                    p_low,
                    p_high,
                    p_res,
                    modulus
                )

                while p is not None and p <= p_high:

                    p_candidates += 1

                    if p <= q:

                        if p * q == n:

                            exact.append(
                                (p, q)
                            )

                            return {
                                "q_tested": q_tested,
                                "p_candidates":
                                    p_candidates,
                                "exact": exact,
                                "time":
                                    time.perf_counter()
                                    - start,
                            }

                    p += modulus

                    if p > p_high:
                        break

    return {
        "q_tested": q_tested,
        "p_candidates":
            p_candidates,
        "exact": exact,
        "time":
            time.perf_counter()
            - start,
    }


# ============================================================
# REPORT
# ============================================================

def report(
    n,
    true_p,
    true_q,
    pside,
    qside
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
    print("MODULI")

    print(
        "R00 =",
        R00
    )

    print(
        "R11 =",
        R11
    )

    print(
        "M   =",
        M
    )

    print()
    print("P-DIRECTION")

    print(
        "p tested      =",
        pside["p_tested"]
    )

    print(
        "q candidates  =",
        pside["q_candidates"]
    )

    print(
        "exact         =",
        len(pside["exact"])
    )

    print(
        "runtime       = %.6f s"
        % pside["time"]
    )

    print()
    print("Q-DIRECTION")

    print(
        "q tested      =",
        qside["q_tested"]
    )

    print(
        "p candidates  =",
        qside["p_candidates"]
    )

    print(
        "exact         =",
        len(qside["exact"])
    )

    print(
        "runtime       = %.6f s"
        % qside["time"]
    )

    print()
    print("P-DIRECTION EXACT")

    for p, q in pside["exact"]:

        print(
            "  "
            f"{p} * {q} = {p*q}"
        )

    print()
    print("Q-DIRECTION EXACT")

    for p, q in qside["exact"]:

        print(
            "  "
            f"{p} * {q} = {p*q}"
        )

    print()
    print(
        "TRUE FOUND FROM P =",
        (
            true_p,
            true_q
        ) in pside["exact"]
    )

    print(
        "TRUE FOUND FROM Q =",
        (
            true_p,
            true_q
        ) in qside["exact"]
    )


# ============================================================
# MAIN
# ============================================================

def main():

    random.seed(
        RANDOM_SEED
    )

    print("=" * 72)
    print("START EXPERIMENT 122")
    print("Reciprocal two-sided interval search")
    print("=" * 72)

    print()
    print(
        "R00 =",
        R00
    )

    print(
        "R11 =",
        R11
    )

    print(
        "M   =",
        M
    )

    for bits in (
        30,
        36,
        42,
        48,
    ):

        print()
        print("=" * 72)
        print(
            f"GENERATING {bits}-BIT SEMIPRIME"
        )
        print("=" * 72)

        n, p, q = make_semiprime(
            bits
        )

        pside = search_from_p(
            n
        )

        qside = search_from_q(
            n
        )

        report(
            n,
            p,
            q,
            pside,
            qside
        )

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 122")
    print("=" * 72)


if __name__ == "__main__":
    main()
