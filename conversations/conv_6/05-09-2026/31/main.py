#!/usr/bin/env python3

import math
import random
import time


# ============================================================
# START EXPERIMENT 114
# Difference-space factor recovery
# ============================================================

R1 = [3, 29]
R2 = [5, 31]

r10, r11 = R1
r20, r21 = R2

R00 = r10 * r20
R11 = r11 * r21
M = math.lcm(R00, R11)

RANDOM_SEED = 114


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
# VERIFY THE COMPLETE CARRY GRID
# ============================================================

def verify_grid(n, p, q):

    if p > q:
        p, q = q, p

    values = []

    for r1 in R1:

        k = p // r1
        a = p % r1

        for r2 in R2:

            l = q // r2
            b = q % r2

            Q = n // (r1 * r2)

            E_quot = Q - k * l

            E_carry = carries(
                r1,
                r2,
                k,
                l,
                a,
                b
            )[3]

            values.append(
                E_quot == E_carry
            )

    return all(values)


# ============================================================
# FERMAT SEARCH
# ============================================================

def fermat_factor(n):

    # p <= q
    #
    # q-p = d
    #
    # s = p+q
    #
    # s^2 - d^2 = 4n

    a = math.isqrt(n)

    if a * a < n:
        a += 1

    iterations = 0

    while True:

        iterations += 1

        b2 = a * a - n

        b = math.isqrt(b2)

        if b * b == b2:

            p = a - b
            q = a + b

            if p > 1 and q > 1 and p * q == n:

                return p, q, iterations

        a += 1


# ============================================================
# FERMAT + CARRY VERIFICATION
# ============================================================

def fermat_with_carry(n):

    start = time.perf_counter()

    p, q, iterations = fermat_factor(n)

    elapsed = (
        time.perf_counter()
        - start
    )

    carry_ok = verify_grid(
        n,
        p,
        q
    )

    return {
        "p": p,
        "q": q,
        "iterations": iterations,
        "time": elapsed,
        "carry_ok": carry_ok,
    }


# ============================================================
# DIRECT p SCAN BASELINE
# ============================================================

def direct_scan(n):

    limit = math.isqrt(n)

    tested = 0

    start = time.perf_counter()

    for p in range(
        3,
        limit + 1,
        2
    ):

        tested += 1

        if n % p == 0:

            q = n // p

            elapsed = (
                time.perf_counter()
                - start
            )

            return {
                "p": p,
                "q": q,
                "tested": tested,
                "time": elapsed,
            }

    return {
        "p": None,
        "q": None,
        "tested": tested,
        "time": (
            time.perf_counter()
            - start
        ),
    }


# ============================================================
# CARRY-MODULAR DIRECT SCAN
# ============================================================

def carry_modular_scan(n):

    """
    This is the Experiment-113 style baseline.

    For each p:

        pq = n mod R00
        pq = n mod R11

    are combined to force q mod M.

    Then the second quotient bucket limits q.

    Finally the exact product is checked.
    """

    limit = math.isqrt(n)

    p_tested = 0
    q_candidates = 0

    start = time.perf_counter()

    for p in range(
        3,
        limit + 1,
        2
    ):

        p_tested += 1

        if math.gcd(p, M) != 1:
            continue

        # q modulo the two moduli

        q1 = (
            n
            * pow(
                p,
                -1,
                R00
            )
        ) % R00

        q2 = (
            n
            * pow(
                p,
                -1,
                R11
            )
        ) % R11

        # ----------------------------------------------------
        # CRT
        # ----------------------------------------------------

        # q = q1 (mod R00)
        #
        # q = q2 (mod R11)

        inv = pow(
            R00,
            -1,
            R11
        )

        t = (
            (q2 - q1)
            * inv
        ) % R11

        q_res = (
            q1
            + R00 * t
        ) % M

        # ----------------------------------------------------
        # Second quotient bucket.
        # ----------------------------------------------------

        Q11 = n // R11

        q_low = (
            R11 * Q11 + p - 1
        ) // p

        q_high = (
            R11 * (Q11 + 1) - 1
        ) // p

        # First q in interval with q=q_res mod M.

        q = (
            q_low
            + (q_res - q_low) % M
        )

        while q <= q_high:

            q_candidates += 1

            if p <= q and p * q == n:

                elapsed = (
                    time.perf_counter()
                    - start
                )

                return {
                    "p": p,
                    "q": q,
                    "p_tested": p_tested,
                    "q_candidates": q_candidates,
                    "time": elapsed,
                }

            q += M

    return {
        "p": None,
        "q": None,
        "p_tested": p_tested,
        "q_candidates": q_candidates,
        "time": (
            time.perf_counter()
            - start
        ),
    }


# ============================================================
# REPORT
# ============================================================

def report(
    n,
    true_p,
    true_q,
    direct,
    modular,
    fermat,
):

    gap = true_q - true_p

    sqrt_n = math.isqrt(n)

    print()
    print("-" * 72)

    print(
        "n bits =",
        n.bit_length()
    )

    print()
    print("n      =", n)
    print("true p =", true_p)
    print("true q =", true_q)

    print()
    print("FACTOR GAP")

    print(
        "q-p     =",
        gap
    )

    print(
        "sqrt(n) =",
        sqrt_n
    )

    print(
        "gap/sqrt(n) =",
        gap / sqrt_n
    )

    print()
    print("DIRECT P-SCAN")

    print(
        "p tested =",
        direct["tested"]
    )

    print(
        "time     = %.6f s"
        % direct["time"]
    )

    print(
        "found    =",
        direct["p"] == true_p
    )

    print()
    print("TWO-CELL MODULAR SCAN")

    print(
        "p tested       =",
        modular["p_tested"]
    )

    print(
        "q candidates   =",
        modular["q_candidates"]
    )

    print(
        "time           = %.6f s"
        % modular["time"]
    )

    print(
        "found          =",
        modular["p"] == true_p
    )

    print()
    print("FERMAT")

    print(
        "iterations     =",
        fermat["iterations"]
    )

    print(
        "time           = %.6f s"
        % fermat["time"]
    )

    print(
        "p              =",
        fermat["p"]
    )

    print(
        "q              =",
        fermat["q"]
    )

    print(
        "carry grid OK   =",
        fermat["carry_ok"]
    )

    print(
        "found           =",
        (
            fermat["p"] == true_p
            and
            fermat["q"] == true_q
        )
    )

    print()
    print("FERMAT ITERATION CHECK")

    expected_bound = (
        gap * gap
    )

    print(
        "(q-p)^2 =",
        expected_bound
    )

    # --------------------------------------------------------
    # Speed ratios
    # --------------------------------------------------------

    print()
    print("TIMING")

    if fermat["time"] > 0:

        print(
            "direct / Fermat = %.3fx"
            % (
                direct["time"]
                / fermat["time"]
            )
        )

        print(
            "modular / Fermat = %.3fx"
            % (
                modular["time"]
                / fermat["time"]
            )
        )


# ============================================================
# MAIN
# ============================================================

def main():

    random.seed(
        RANDOM_SEED
    )

    print("=" * 72)
    print("START EXPERIMENT 114")
    print("Difference-space factor recovery")
    print("=" * 72)

    print()
    print("r1 =", R1)
    print("r2 =", R2)
    print("R00 =", R00)
    print("R11 =", R11)
    print("M   =", M)

    for bits in (
        30,
        36,
        42,
        48,
        54,
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

        # ----------------------------------------------------
        # Direct baseline
        # ----------------------------------------------------

        direct = direct_scan(
            n
        )

        # ----------------------------------------------------
        # Experiment-113 style search
        # ----------------------------------------------------

        modular = carry_modular_scan(
            n
        )

        # ----------------------------------------------------
        # Fermat
        # ----------------------------------------------------

        fermat = fermat_with_carry(
            n
        )

        report(
            n,
            p,
            q,
            direct,
            modular,
            fermat,
        )

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 114")
    print("=" * 72)


if __name__ == "__main__":
    main()

