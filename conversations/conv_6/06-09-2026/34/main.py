#!/usr/bin/env python3

import math
import random
import time


# ============================================================================
# START EXPERIMENT 190
#
# SEGMENTED FERMAT QUADRATIC-RESIDUE SIEVE
#
# Fermat:
#
#     x = (p + q) / 2
#     y = (q - p) / 2
#
#     x^2 - y^2 = n
#
# Therefore:
#
#     x^2 = n + y^2
#
# A necessary condition is:
#
#     n + y^2
#
# must be a quadratic residue modulo every chosen prime.
#
# Experiment 189 constructed a global CRT wheel.
#
# Experiment 190 instead works on SEGMENTS of y:
#
#     [start, start + BLOCK_SIZE)
#
# and marks impossible residues directly.
#
# This avoids constructing a potentially huge CRT residue table.
# ============================================================================


R1 = [17, 43, 59, 71, 83]
R2 = [19, 37, 61, 73, 89]

FILTER_PRIMES = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]

DEPTHS = [1, 2, 3, 4, 5, 6, 7]

BLOCK_SIZE = 1 << 20

# Maximum y candidates examined by each method.
MAX_Y = 25_000_000


# ============================================================================
# DETERMINISTIC 64-BIT MILLER-RABIN
# ============================================================================

MR_BASES = [
    2,
    325,
    9375,
    28178,
    450775,
    9780504,
    1795265022,
]


def is_prime(n: int) -> bool:
    if n < 2:
        return False

    small_primes = (
        2, 3, 5, 7, 11, 13, 17, 19,
        23, 29, 31, 37
    )

    for p in small_primes:
        if n == p:
            return True
        if n % p == 0:
            return False

    d = n - 1
    s = 0

    while (d & 1) == 0:
        s += 1
        d >>= 1

    for a in MR_BASES:
        if a % n == 0:
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


def next_prime(n: int) -> int:
    if n <= 2:
        return 2

    if n % 2 == 0:
        n += 1

    while not is_prime(n):
        n += 2

    return n


# ============================================================================
# CONTROLLED-BALANCED TEST INSTANCES
# ============================================================================

def make_instance(bits: int, target_y: int):
    """
    Construct p,q with:

        y = (q-p)/2

    close to target_y.

    This lets us test larger integer sizes without requiring
    hundreds of millions of Fermat iterations.
    """
    half = bits // 2

    # Start slightly below 2^(bits/2).
    base = (1 << half) - 1_000_003

    p = base

    if p % 2 == 0:
        p -= 1

    while not is_prime(p):
        p -= 2

    q_start = p + 2 * target_y

    if q_start <= p:
        raise ValueError("target_y must be positive")

    q = next_prime(q_start)

    n = p * q

    # Retry with a slightly different p if product has wrong size.
    while n.bit_length() != bits:
        base -= 100_003
        p = base

        if p % 2 == 0:
            p -= 1

        while not is_prime(p):
            p -= 2

        q = next_prime(p + 2 * target_y)
        n = p * q

    return p, q


# ============================================================================
# FERMAT TEST
# ============================================================================

def fermat_test_y(n: int, y: int):
    value = n + y * y

    x = math.isqrt(value)

    if x * x != value:
        return None

    p = x - y
    q = x + y

    if p > 1 and q > 1 and p * q == n:
        return p, q

    return None


# ============================================================================
# ORDINARY FERMAT
# ============================================================================

def ordinary_fermat(n: int, max_y: int):
    start = time.perf_counter()

    for y in range(max_y + 1):
        result = fermat_test_y(n, y)

        if result is not None:
            elapsed = time.perf_counter() - start
            return y, y + 1, result, elapsed

    elapsed = time.perf_counter() - start

    return None, max_y + 1, None, elapsed


# ============================================================================
# QUADRATIC RESIDUE TABLE
# ============================================================================

def build_qr_table(prime: int):
    return {
        (x * x) % prime
        for x in range(prime)
    }


# ============================================================================
# PRECOMPUTE BAD RESIDUES
# ============================================================================

def build_bad_residues(n: int, prime: int):
    """
    Return residues r modulo prime for which

        n + r^2

    is NOT a quadratic residue modulo prime.
    """
    qr = build_qr_table(prime)

    bad = []

    for r in range(prime):
        if (n + r * r) % prime not in qr:
            bad.append(r)

    return bad


# ============================================================================
# SEGMENTED QR SIEVE
# ============================================================================

def segmented_fermat(
    n: int,
    primes,
    max_y: int,
    block_size: int
):
    """
    Search y = 0..max_y using segmented quadratic-residue marking.

    The segment bytearray contains:

        1 = currently possible
        0 = impossible

    For each prime we mark all y congruence classes that violate

        n + y^2 = quadratic residue (mod prime)

    The bytearray.find(1, ...) call is implemented in C and avoids
    walking every byte in Python.
    """
    start_time = time.perf_counter()

    bad_by_prime = {
        p: build_bad_residues(n, p)
        for p in primes
    }

    build_time = time.perf_counter() - start_time

    tested_y = 0
    blocks = 0
    surviving = 0

    search_start = time.perf_counter()

    block_start = 0

    while block_start <= max_y:

        block_end = min(
            block_start + block_size,
            max_y + 1
        )

        block_len = block_end - block_start

        # 1 = alive
        # 0 = eliminated
        alive = bytearray(b"\x01") * block_len

        # ------------------------------------------------------------
        # MARK IMPOSSIBLE RESIDUES
        # ------------------------------------------------------------

        for prime in primes:

            bad = bad_by_prime[prime]

            for residue in bad:
                if residue < block_start:
                    first = (
                        block_start
                        + ((residue - block_start) % prime)
                    )
                else:
                    first = residue

                if first >= block_end:
                    continue

                idx = first - block_start

                count = (
                    (block_len - 1 - idx) // prime
                ) + 1

                alive[idx:block_len:prime] = (
                    b"\x00" * count
                )

        # ------------------------------------------------------------
        # TEST SURVIVING Y VALUES
        # ------------------------------------------------------------

        pos = 0

        while True:
            pos = alive.find(1, pos)

            if pos < 0:
                break

            y = block_start + pos

            surviving += 1

            result = fermat_test_y(n, y)

            tested_y += 1

            if result is not None:
                elapsed = time.perf_counter() - search_start

                return {
                    "found_y": y,
                    "result": result,
                    "tested_survivors": tested_y,
                    "surviving_total": surviving,
                    "blocks": blocks + 1,
                    "build_time": build_time,
                    "search_time": elapsed,
                    "total_time": build_time + elapsed,
                }

            pos += 1

        blocks += 1

        block_start = block_end

    elapsed = time.perf_counter() - search_start

    return {
        "found_y": None,
        "result": None,
        "tested_survivors": tested_y,
        "surviving_total": surviving,
        "blocks": blocks,
        "build_time": build_time,
        "search_time": elapsed,
        "total_time": build_time + elapsed,
    }


# ============================================================================
# INSTANCE
# ============================================================================

def run_instance(bits: int, p: int, q: int):

    n = p * q

    true_y = (q - p) // 2
    true_x = (q + p) // 2

    print()
    print("=" * 72)
    print(f"START INSTANCE {bits}-BIT")
    print("=" * 72)

    print(f"bits(n) = {n.bit_length()}")
    print(f"p = {p}")
    print(f"q = {q}")
    print(f"true x = {true_x}")
    print(f"true y = {true_y}")
    print(f"gap = {q - p}")

    # ------------------------------------------------------------------------
    # ORDINARY
    # ------------------------------------------------------------------------

    print()
    print("ORDINARY FERMAT")

    ordinary_y, ordinary_tests, ordinary_result, ordinary_time = (
        ordinary_fermat(n, MAX_Y)
    )

    print(f"    tests = {ordinary_tests}")
    print(f"    time = {ordinary_time:.6f}s")
    print(f"    found y = {ordinary_y}")
    print(f"    result = {ordinary_result}")

    # ------------------------------------------------------------------------
    # SEGMENTED DEPTHS
    # ------------------------------------------------------------------------

    for depth in DEPTHS:

        primes = FILTER_PRIMES[:depth]

        print()
        print(f"SEGMENTED DEPTH {depth}")
        print(f"    primes = {primes}")
        print(f"    block size = {BLOCK_SIZE}")

        result = segmented_fermat(
            n,
            primes,
            MAX_Y,
            BLOCK_SIZE
        )

        found_y = result["found_y"]
        found = result["result"]

        print(f"    blocks = {result['blocks']}")
        print(f"    survivors tested = {result['tested_survivors']}")
        print(f"    build time = {result['build_time']:.6f}s")
        print(f"    search time = {result['search_time']:.6f}s")
        print(f"    total time = {result['total_time']:.6f}s")
        print(f"    found y = {found_y}")
        print(f"    FOUND = {found}")

        if found is not None:
            correct = (
                found == (p, q)
                or found == (q, p)
            )

            print(f"    correct = {correct}")

            if ordinary_time > 0:
                print(
                    "    speedup vs ordinary = "
                    f"{ordinary_time / result['total_time']:.3f}x"
                )

    # ------------------------------------------------------------------------
    # EXACT IDENTITY
    # ------------------------------------------------------------------------

    print()
    print("EXACT FERMAT IDENTITY")

    print(f"    x = {true_x}")
    print(f"    y = {true_y}")
    print(f"    x²-y² = {true_x * true_x - true_y * true_y}")
    print(f"    n = {n}")
    print(
        f"    identity = "
        f"{(true_x * true_x - true_y * true_y) == n}"
    )

    print()
    print(f"FINISHED INSTANCE {bits}-BIT")


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("START EXPERIMENT 190")
    print()
    print("Segmented Fermat quadratic-residue sieve")
    print(f"BLOCK_SIZE = {BLOCK_SIZE}")
    print(f"MAX_Y = {MAX_Y}")
    print()

    # ------------------------------------------------------------------------
    # Existing-style balanced instances.
    #
    # These use the same factor pairs from the previous experiments so the
    # results remain directly comparable.
    # ------------------------------------------------------------------------

    existing = {
        30: (16451, 64271),
        36: (131267, 522517),
        42: (1049089, 4191259),
        48: (8390069, 33547589),
    }

    for bits in [30, 36, 42, 48]:
        p, q = existing[bits]

        run_instance(bits, p, q)

    # ------------------------------------------------------------------------
    # New larger controlled-gap instances.
    #
    # Keep y below MAX_Y so all methods can actually reach the factor.
    # ------------------------------------------------------------------------

    new_specs = [
        (54, 10_000_000),
        (60, 15_000_000),
        (66, 20_000_000),
    ]

    for bits, target_y in new_specs:

        p, q = make_instance(bits, target_y)

        run_instance(bits, p, q)

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 190")
    print("=" * 72)


if __name__ == "__main__":
    main()
