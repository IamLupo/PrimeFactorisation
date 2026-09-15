#!/usr/bin/env python3

import math
import random
import time


# ============================================================
# START EXPERIMENT 121
# Tiny-q interval intersection
# ============================================================

R00 = 15
R11 = 899

RANDOM_SEED = 121


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
# SECOND-CELL INTERVAL
# ============================================================

def second_q_interval(n, p):

    Q11 = n // R11

    low = (
        R11 * Q11 + p - 1
    ) // p

    high = (
        R11 * (Q11 + 1) - 1
    ) // p

    return low, high


# ============================================================
# FIRST MODULAR CONDITION
# ============================================================

def q_mod_15(n, p):

    if math.gcd(p, R00) != 1:
        return None

    return (
        n
        * pow(
            p,
            -1,
            R00
        )
    ) % R00


# ============================================================
# SECOND MODULAR CONDITION
# ============================================================

def product_mod_899_ok(n, p, q):

    return (
        (p * q) % R11
        ==
        n % R11
    )


# ============================================================
# SEARCH
# ============================================================

def search(n):

    limit = math.isqrt(n)

    counters = {
        "p_tested": 0,
        "p_mod15": 0,

        "interval_empty": 0,
        "interval_width_total": 0,
        "interval_width_1": 0,
        "interval_width_2": 0,
        "interval_width_gt2": 0,

        "mod15_candidates": 0,
        "mod899_pass": 0,

        "exact": 0,
    }

    exact = []

    start = time.perf_counter()

    # --------------------------------------------------------
    # Search p directly.
    # --------------------------------------------------------

    for p in range(
        3,
        limit + 1,
        2
    ):

        counters["p_tested"] += 1

        q15 = q_mod_15(
            n,
            p
        )

        if q15 is None:
            continue

        counters["p_mod15"] += 1

        # ----------------------------------------------------
        # Very narrow second-cell interval.
        # ----------------------------------------------------

        q_low, q_high = second_q_interval(
            n,
            p
        )

        if q_low > q_high:

            counters["interval_empty"] += 1

            continue

        width = (
            q_high
            - q_low
            + 1
        )

        counters["interval_width_total"] += width

        if width == 1:
            counters["interval_width_1"] += 1
        elif width == 2:
            counters["interval_width_2"] += 1
        else:
            counters["interval_width_gt2"] += 1

        # ----------------------------------------------------
        # Intersect the interval with
        #
        # q == q15 (mod 15).
        # ----------------------------------------------------

        delta = (
            q15 - q_low
        ) % R00

        q = q_low + delta

        while q <= q_high:

            counters["mod15_candidates"] += 1

            # ------------------------------------------------
            # Second-cell modulus.
            # ------------------------------------------------

            if product_mod_899_ok(
                n,
                p,
                q
            ):

                counters["mod899_pass"] += 1

                # ------------------------------------------------
                # Exact final verification.
                # ------------------------------------------------

                if p * q == n:

                    counters["exact"] += 1

                    exact.append(
                        (p, q)
                    )

                    elapsed = (
                        time.perf_counter()
                        - start
                    )

                    return {
                        "counters": counters,
                        "exact": exact,
                        "time": elapsed,
                    }

            q += R00

    elapsed = (
        time.perf_counter()
        - start
    )

    return {
        "counters": counters,
        "exact": exact,
        "time": elapsed,
    }


# ============================================================
# DIRECT BASELINE
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

            return {
                "tested": tested,
                "p": p,
                "q": q,
                "time":
                    time.perf_counter()
                    - start,
            }

    return {
        "tested": tested,
        "p": None,
        "q": None,
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
    result,
    direct
):

    c = result["counters"]

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
    print("MODULI")

    print("R00 =", R00)
    print("R11 =", R11)

    print()
    print("DIRECT CONTROL")

    print(
        "p tested =",
        direct["tested"]
    )

    print(
        "runtime  = %.6f s"
        % direct["time"]
    )

    print()
    print("EXPERIMENT 121")

    print(
        "p tested             =",
        c["p_tested"]
    )

    print(
        "p passing mod 15     =",
        c["p_mod15"]
    )

    print(
        "interval empty       =",
        c["interval_empty"]
    )

    print(
        "interval width = 1   =",
        c["interval_width_1"]
    )

    print(
        "interval width = 2   =",
        c["interval_width_2"]
    )

    print(
        "interval width > 2   =",
        c["interval_width_gt2"]
    )

    print(
        "total q interval size =",
        c["interval_width_total"]
    )

    print(
        "mod-15 q candidates  =",
        c["mod15_candidates"]
    )

    print(
        "mod-899 pass         =",
        c["mod899_pass"]
    )

    print(
        "exact                =",
        c["exact"]
    )

    print(
        "runtime              = %.6f s"
        % result["time"]
    )

    print()
    print("EXACT")

    for p, q in result["exact"]:

        print(
            "  "
            f"{p} * {q} = {p*q}"
        )

    print()
    print(
        "TRUE FOUND =",
        (
            true_p,
            true_q
        ) in result["exact"]
    )


# ============================================================
# MAIN
# ============================================================

def main():

    random.seed(
        RANDOM_SEED
    )

    print("=" * 72)
    print("START EXPERIMENT 121")
    print("Tiny-q interval intersection")
    print("=" * 72)

    print()

    print("R00 =", R00)
    print("R11 =", R11)

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

        direct = direct_scan(
            n
        )

        result = search(
            n
        )

        report(
            n,
            p,
            q,
            result,
            direct
        )

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 121")
    print("=" * 72)


if __name__ == "__main__":
    main()
