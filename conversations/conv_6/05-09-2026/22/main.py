#!/usr/bin/env python3

import math
import random
import time


# ============================================================
# START EXPERIMENT 105
# Larger first-cell scale
# ============================================================

R1 = [17, 43]
R2 = [19, 47]

R00 = R1[0] * R2[0]

RANDOM_SEED = 105


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
# ALGEBRAIC l SOLVER
# ============================================================

def solve_l0(Q, r1, r2, k, a, b):

    solutions = []

    c1 = (k * b) // r2
    beta = (k * b) % r2

    denominator = k * r1 + a

    for t in range(r1):

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

        if numerator < 0:
            continue

        if numerator % denominator != 0:
            continue

        m = numerator // denominator

        if m < 0:
            continue

        l = t + r1 * m

        if l <= 0:
            continue

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
            continue

        solutions.append(l)

    return solutions


# ============================================================
# SEARCH
# ============================================================

def search(n):

    r10, r11 = R1
    r20, r21 = R2

    Q00 = n // (r10 * r20)

    counters = {
        "k0_tested": 0,
        "residue_tests": 0,
        "first_cell_solutions": 0,

        "bucket_survivors": 0,
        "mod_R00_survivors": 0,

        "time_l_solver": 0.0,
        "time_total": 0.0,
    }

    bucket_survivors = []
    congruence_survivors = []

    p_limit = math.isqrt(n)

    k0_max = p_limit // r10

    search_start = time.perf_counter()

    for k0 in range(1, k0_max + 1):

        counters["k0_tested"] += 1

        for a0 in range(r10):

            for b0 in range(r20):

                counters["residue_tests"] += 1

                t0 = time.perf_counter()

                l_values = solve_l0(
                    Q00,
                    r10,
                    r20,
                    k0,
                    a0,
                    b0,
                )

                counters["time_l_solver"] += (
                    time.perf_counter() - t0
                )

                counters["first_cell_solutions"] += len(
                    l_values
                )

                for l0 in l_values:

                    p = r10 * k0 + a0
                    q = r20 * l0 + b0

                    if p > q:
                        continue

                    product = p * q

                    # ------------------------------------------------
                    # Quotient-bucket condition.
                    # ------------------------------------------------

                    if product // R00 != n // R00:
                        continue

                    counters["bucket_survivors"] += 1

                    error = product - n

                    bucket_survivors.append(
                        (p, q, error)
                    )

                    # ------------------------------------------------
                    # Congruence closure.
                    # ------------------------------------------------

                    if product % R00 != n % R00:
                        continue

                    counters["mod_R00_survivors"] += 1

                    congruence_survivors.append(
                        (p, q, error)
                    )

    counters["time_total"] = (
        time.perf_counter()
        - search_start
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

    print("n bits =", n.bit_length())

    print()
    print("n      =", n)
    print("true p =", true_p)
    print("true q =", true_q)

    print()
    print("GRID")
    print("r1 =", R1)
    print("r2 =", R2)

    print()
    print("R00 =", R00)

    print()
    print("SEARCH RESULTS")

    print(
        "k0 tested             =",
        counters["k0_tested"]
    )

    print(
        "residue tests         =",
        counters["residue_tests"]
    )

    print(
        "first-cell solutions  =",
        counters["first_cell_solutions"]
    )

    print(
        "bucket survivors      =",
        counters["bucket_survivors"]
    )

    print(
        "mod-R00 survivors     =",
        counters["mod_R00_survivors"]
    )

    print(
        "l-solver time         = %.6f s"
        % counters["time_l_solver"]
    )

    print(
        "total search time     = %.6f s"
        % counters["time_total"]
    )

    # --------------------------------------------------------
    # Bucket theorem check.
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
    # Congruence survivors.
    # --------------------------------------------------------

    print()
    print("CONGRUENCE SURVIVORS")

    for p, q, error in congruence[:25]:

        print(
            "  "
            f"p={p} "
            f"q={q} "
            f"pq-n={error}"
        )

    if len(congruence) > 25:

        print(
            "  ... "
            f"{len(congruence) - 25} more"
        )

    # --------------------------------------------------------
    # Exact verification.
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
    print("START EXPERIMENT 105")
    print("Larger first-cell scale")
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
    print("FINISHED EXPERIMENT 105")
    print("=" * 72)


if __name__ == "__main__":
    main()
