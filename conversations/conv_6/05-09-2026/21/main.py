#!/usr/bin/env python3

import math
import random
import time


# ============================================================
# START EXPERIMENT 104
# Quotient bucket + congruence closure
# ============================================================

R1 = [3, 29]
R2 = [5, 31]

R00 = R1[0] * R2[0]


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

    for a in (2, 3, 5, 7, 11, 13, 17):

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
# CARRY EQUATION
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
# EXACT l SOLVER
# ============================================================

def solve_l0(Q, r1, r2, k, a, b):

    solutions = []

    c1 = (k * b) // r2
    beta = (k * b) % r2

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

        denominator = k * r1 + a

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

        # Verify first-cell carry equation.

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
    Q10 = n // (r11 * r20)
    Q01 = n // (r10 * r21)
    Q11 = n // (r11 * r21)

    p_limit = math.isqrt(n)

    k0_max = p_limit // r10

    counters = {
        "k0_tested": 0,
        "residue_tests": 0,
        "first_cell_solutions": 0,

        "bucket_survivors": 0,

        "mod_R00_survivors": 0,
        "prime_survivors": 0,

        "exact": 0,
    }

    bucket_survivors = []
    congruence_survivors = []
    prime_survivors = []

    for k0 in range(1, k0_max + 1):

        counters["k0_tested"] += 1

        for a0 in range(r10):

            for b0 in range(r20):

                counters["residue_tests"] += 1

                l0_values = solve_l0(
                    Q00,
                    r10,
                    r20,
                    k0,
                    a0,
                    b0
                )

                counters["first_cell_solutions"] += len(
                    l0_values
                )

                for l0 in l0_values:

                    p = r10 * k0 + a0
                    q = r20 * l0 + b0

                    if p > q:
                        continue

                    # ------------------------------------------------
                    # FIRST-CELL QUOTIENT BUCKET
                    # ------------------------------------------------

                    product = p * q

                    if product // R00 != n // R00:
                        continue

                    counters["bucket_survivors"] += 1

                    bucket_survivors.append(
                        (p, q, product - n)
                    )

                    # ------------------------------------------------
                    # CONGRUENCE CLOSURE
                    #
                    # Because:
                    #
                    # |pq-n| < R00
                    #
                    # and
                    #
                    # pq == n (mod R00),
                    #
                    # we must have pq=n.
                    # ------------------------------------------------

                    if product % R00 != n % R00:
                        continue

                    counters["mod_R00_survivors"] += 1

                    congruence_survivors.append(
                        (p, q, product - n)
                    )

                    # ------------------------------------------------
                    # PRIME CHECK
                    #
                    # This is additional information, not required
                    # for the congruence closure.
                    # ------------------------------------------------

                    if not is_probable_prime(p):
                        continue

                    if not is_probable_prime(q):
                        continue

                    counters["prime_survivors"] += 1

                    prime_survivors.append(
                        (p, q)
                    )

                    # Exact check is ONLY a final measurement.
                    if product == n:
                        counters["exact"] += 1

    return {
        "counters": counters,
        "bucket_survivors": bucket_survivors,
        "congruence_survivors": congruence_survivors,
        "prime_survivors": prime_survivors,
    }


# ============================================================
# REPORT
# ============================================================

def report(n, true_p, true_q, result, elapsed):

    counters = result["counters"]

    print()
    print("-" * 72)

    print("n bits =", n.bit_length())

    print()
    print("n      =", n)
    print("true p =", true_p)
    print("true q =", true_q)

    print()
    print("R00 =", R00)

    print()
    print("TRUE PRODUCT")
    print("p*q =", true_p * true_q)

    print()
    print("TRUE MOD R00")
    print("n mod R00  =", n % R00)
    print(
        "(p*q) mod R00 =",
        (true_p * true_q) % R00
    )

    print()
    print("SEARCH RESULTS")

    for key, value in counters.items():

        print(
            f"{key:<24} = {value}"
        )

    print(
        "runtime                  = %.6f s"
        % elapsed
    )

    # --------------------------------------------------------
    # Demonstrate the key theorem numerically.
    # --------------------------------------------------------

    if result["bucket_survivors"]:

        maximum_error = max(
            abs(error)
            for _, _, error
            in result["bucket_survivors"]
        )

        print()
        print("BUCKET ERROR")

        print(
            "max |pq-n| among bucket survivors =",
            maximum_error
        )

        print(
            "R00 - 1                         =",
            R00 - 1
        )

    print()
    print("MOD-R00 SURVIVORS")

    for p, q, error in result["congruence_survivors"][:20]:

        print(
            "  "
            f"p={p} "
            f"q={q} "
            f"pq-n={error}"
        )

    if len(result["congruence_survivors"]) > 20:

        print(
            "  ... "
            f"{len(result['congruence_survivors']) - 20}"
            " more"
        )

    print()
    print("PRIME SURVIVORS")

    for p, q in result["prime_survivors"][:20]:

        print(
            "  "
            f"p={p} "
            f"q={q} "
            f"product={p*q}"
        )


# ============================================================
# MAIN
# ============================================================

def main():

    random.seed(104)

    print("=" * 72)
    print("START EXPERIMENT 104")
    print("Quotient bucket + congruence closure")
    print("=" * 72)

    for bits in (30, 36, 42):

        print()
        print("=" * 72)
        print(f"GENERATING {bits}-BIT SEMIPRIME")
        print("=" * 72)

        n, p, q = make_semiprime(bits)

        start = time.perf_counter()

        result = search(n)

        elapsed = (
            time.perf_counter()
            - start
        )

        report(
            n,
            p,
            q,
            result,
            elapsed
        )

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 104")
    print("=" * 72)


if __name__ == "__main__":
    main()
