#!/usr/bin/env python3

import math
import random
import time


# ============================================================
# START EXPERIMENT 117
# Precomputed modular residue sieve
# ============================================================

R1 = [3, 29]
R2 = [5, 31]

R00 = R1[0] * R2[0]      # 15
R11 = R1[1] * R2[1]      # 899

# gcd(15,899)=1
M = R00 * R11             # 13485

RANDOM_SEED = 117


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
# PRECOMPUTE p -> q MOD M
# ============================================================

def build_residue_table(n):
    """
    For every x modulo M:

        gcd(x,M)=1

    precompute

        q = n*x^(-1) mod M.

    The expensive modular inverses are therefore performed
    only M times, instead of once per p candidate.
    """

    table = [-1] * M

    start = time.perf_counter()

    n_mod = n % M

    for x in range(M):

        if math.gcd(x, M) != 1:
            continue

        table[x] = (
            n_mod
            * pow(
                x,
                -1,
                M
            )
        ) % M

    elapsed = (
        time.perf_counter()
        - start
    )

    return table, elapsed


# ============================================================
# SECOND-CELL q INTERVAL
# ============================================================

def q_interval(n, p):

    Q11 = n // R11

    low = (
        R11 * Q11 + p - 1
    ) // p

    high = (
        R11 * (Q11 + 1) - 1
    ) // p

    return low, high


# ============================================================
# SEARCH
# ============================================================

def search(n, table):

    limit = math.isqrt(n)

    p_tested = 0
    table_hits = 0
    bucket_candidates = 0
    exact = []

    start = time.perf_counter()

    # Only odd p for odd semiprimes.
    for p in range(
        3,
        limit + 1,
        2
    ):

        p_tested += 1

        # ----------------------------------------------------
        # O(1) lookup instead of pow(p,-1,M).
        # ----------------------------------------------------

        q_res = table[p % M]

        if q_res < 0:
            continue

        table_hits += 1

        # ----------------------------------------------------
        # Second-cell quotient bucket.
        # ----------------------------------------------------

        q_low, q_high = q_interval(
            n,
            p
        )

        # Find first q >= q_low such that
        #
        # q == q_res (mod M).

        q = (
            q_low
            + (q_res - q_low) % M
        )

        while q <= q_high:

            bucket_candidates += 1

            # p <= q is our orientation.
            if p <= q:

                if p * q == n:

                    exact.append(
                        (p, q)
                    )

                    elapsed = (
                        time.perf_counter()
                        - start
                    )

                    return {
                        "p_tested": p_tested,
                        "table_hits": table_hits,
                        "bucket_candidates":
                            bucket_candidates,
                        "exact": exact,
                        "time": elapsed,
                    }

            q += M

    elapsed = (
        time.perf_counter()
        - start
    )

    return {
        "p_tested": p_tested,
        "table_hits": table_hits,
        "bucket_candidates":
            bucket_candidates,
        "exact": exact,
        "time": elapsed,
    }


# ============================================================
# OLD EXPERIMENT-113 STYLE SEARCH
# ============================================================

def old_search(n):

    limit = math.isqrt(n)

    p_tested = 0
    bucket_candidates = 0
    exact = []

    start = time.perf_counter()

    for p in range(
        3,
        limit + 1,
        2
    ):

        p_tested += 1

        if math.gcd(p, M) != 1:
            continue

        q_res = (
            n
            * pow(
                p,
                -1,
                M
            )
        ) % M

        q_low, q_high = q_interval(
            n,
            p
        )

        q = (
            q_low
            + (q_res - q_low) % M
        )

        while q <= q_high:

            bucket_candidates += 1

            if p <= q and p * q == n:

                exact.append(
                    (p, q)
                )

                elapsed = (
                    time.perf_counter()
                    - start
                )

                return {
                    "p_tested": p_tested,
                    "bucket_candidates":
                        bucket_candidates,
                    "exact": exact,
                    "time": elapsed,
                }

            q += M

    return {
        "p_tested": p_tested,
        "bucket_candidates":
            bucket_candidates,
        "exact": exact,
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
    table_time,
    old,
    new
):

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
    print("MODULUS")

    print("R00 =", R00)
    print("R11 =", R11)
    print("M   =", M)

    print()
    print("TABLE BUILD")

    print(
        "time = %.9f s"
        % table_time
    )

    print()
    print("OLD EXPERIMENT 113")

    print(
        "p tested       =",
        old["p_tested"]
    )

    print(
        "q candidates   =",
        old["bucket_candidates"]
    )

    print(
        "exact          =",
        len(old["exact"])
    )

    print(
        "runtime        = %.6f s"
        % old["time"]
    )

    print()
    print("EXPERIMENT 117")

    print(
        "p tested       =",
        new["p_tested"]
    )

    print(
        "table hits     =",
        new["table_hits"]
    )

    print(
        "q candidates   =",
        new["bucket_candidates"]
    )

    print(
        "exact          =",
        len(new["exact"])
    )

    print(
        "runtime        =",
        "%.6f s"
        % new["time"]
    )

    print()
    print("TOTAL WITH TABLE")

    total_new = (
        table_time
        + new["time"]
    )

    print(
        "runtime        =",
        "%.6f s"
        % total_new
    )

    if new["time"] > 0:

        print()
        print("SEARCH SPEEDUP")

        print(
            "old/new = %.3fx"
            % (
                old["time"]
                / new["time"]
            )
        )

    if total_new > 0:

        print(
            "old/(table+new) = %.3fx"
            % (
                old["time"]
                / total_new
            )
        )

    print()
    print("RESULT")

    for p, q in new["exact"]:

        print(
            "  "
            f"{p} * {q} = {p*q}"
        )

    print()
    print(
        "TRUE FOUND =",
        (true_p, true_q)
        in new["exact"]
    )


# ============================================================
# MAIN
# ============================================================

def main():

    random.seed(
        RANDOM_SEED
    )

    print("=" * 72)
    print("START EXPERIMENT 117")
    print("Precomputed modular residue sieve")
    print("=" * 72)

    print()
    print("r1 =", R1)
    print("r2 =", R2)
    print("R00 =", R00)
    print("R11 =", R11)
    print("M   =", M)

    # --------------------------------------------------------
    # Build table once.
    # --------------------------------------------------------

    # Use a representative n for the residue table.
    #
    # IMPORTANT:
    # The table depends on n mod M, so it must be rebuilt
    # for every new n.
    # --------------------------------------------------------

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

        n, true_p, true_q = make_semiprime(
            bits
        )

        table, table_time = build_residue_table(
            n
        )

        print()
        print(
            "Residue table built."
        )

        old = old_search(
            n
        )

        new = search(
            n,
            table
        )

        report(
            n,
            true_p,
            true_q,
            table_time,
            old,
            new
        )

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 117")
    print("=" * 72)


if __name__ == "__main__":
    main()
