#!/usr/bin/env python3

import math
import time


# ============================================================================
# START EXPERIMENT 191
#
# QUADRATIC-RESIDUE SIEVE SCALING STUDY
#
# Experiment 190 showed that the segmented QR sieve gives large practical
# speedups for Fermat factorization.
#
# Experiment 191 asks two separate questions:
#
#   1. How does the surviving y-density behave as n gets larger?
#
#   2. Does that density remain useful when the factors become progressively
#      farther apart?
#
# We therefore separate:
#
#   A) DENSITY MEASUREMENT
#      Measure the fraction of y-values surviving the modular sieve over a
#      fixed interval. This is cheap even for large integers.
#
#   B) ACTUAL FERMAT RECOVERY
#      For smaller instances, actually search until the hidden y is reached.
#
# Important:
# The sieve itself never uses p or q. The hidden factors are used only to
# generate controlled test numbers and to verify the result.
#
# ============================================================================


FILTER_PRIMES = [
    3, 5, 7, 11, 13, 17, 19, 23, 29, 31
]

DEPTHS = list(range(1, len(FILTER_PRIMES) + 1))

BLOCK_SIZE = 1 << 20

# Number of y values used for the density experiment.
DENSITY_RANGE = 1 << 22       # 4,194,304 values

# Actual Fermat search is deliberately capped.
MAX_SEARCH_Y = 25_000_000

# Actual recovery only for these bit sizes.
RECOVERY_BITS = [48, 54, 60, 66]

# Factor-distance ratios:
#
#     y = (q-p)/2
#
# expressed relative to p.
#
# Small ratio => close factors => Fermat is favorable.
# Large ratio => distant factors => Fermat becomes expensive.
GAP_RATIOS = [
    0.001,
    0.005,
    0.01,
    0.025,
    0.05,
    0.10,
]


# ============================================================================
# DETERMINISTIC MILLER-RABIN
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
        d >>= 1
        s += 1

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


def previous_prime(n: int) -> int:

    if n <= 2:
        return 2

    if n % 2 == 0:
        n -= 1

    while not is_prime(n):
        n -= 2

    return n


def next_prime(n: int) -> int:

    if n <= 2:
        return 2

    if n % 2 == 0:
        n += 1

    while not is_prime(n):
        n += 2

    return n


# ============================================================================
# CONTROLLED INSTANCE GENERATION
# ============================================================================

def make_instance(bits: int, gap_ratio: float):

    """
    Construct:

        p ~ 2^(bits/2)
        y ~ gap_ratio * p
        q = p + 2y

    and move q to the next prime.

    This gives controlled factor-distance ratios.
    """

    half_bits = bits // 2

    center = 1 << half_bits

    # Start slightly below 2^(bits/2).
    p0 = center - 1_000_003

    p = previous_prime(p0)

    target_y = max(
        1,
        int(round(gap_ratio * p))
    )

    q_target = p + 2 * target_y

    q = next_prime(q_target)

    n = p * q

    # Ensure the resulting n has the requested bit length.
    while n.bit_length() != bits:

        p -= 1009

        p = previous_prime(p)

        target_y = max(
            1,
            int(round(gap_ratio * p))
        )

        q = next_prime(p + 2 * target_y)

        n = p * q

    return p, q


# ============================================================================
# FERMAT VALUES
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
# BUILD QR CONDITIONS
# ============================================================================

def build_bad_residues(n: int, prime: int):

    """
    Compute residues r modulo prime such that

        n + r^2

    is NOT a quadratic residue modulo prime.
    """

    qr = {
        (x * x) % prime
        for x in range(prime)
    }

    bad = []

    for r in range(prime):

        if (n + r * r) % prime not in qr:
            bad.append(r)

    return bad


# ============================================================================
# SEGMENTED DENSITY MEASUREMENT
# ============================================================================

def measure_density(
    n: int,
    primes,
    max_y: int,
    block_size: int
):
    """
    Measure QR survivor density over y=0..max_y.

    No square-root / isqrt operations are performed here.

    This isolates the sieve itself.
    """

    bad_by_prime = {
        p: build_bad_residues(n, p)
        for p in primes
    }

    survivors = 0
    total = max_y + 1

    start_time = time.perf_counter()

    block_start = 0

    while block_start <= max_y:

        block_end = min(
            block_start + block_size,
            max_y + 1
        )

        block_len = block_end - block_start

        alive = bytearray(
            b"\x01"
        ) * block_len

        for prime in primes:

            for residue in bad_by_prime[prime]:

                first = (
                    block_start
                    + ((residue - block_start) % prime)
                )

                if first >= block_end:
                    continue

                idx = first - block_start

                count = (
                    (block_len - 1 - idx) // prime
                ) + 1

                alive[idx:block_len:prime] = (
                    b"\x00" * count
                )

        survivors += alive.count(1)

        block_start = block_end

    elapsed = time.perf_counter() - start_time

    density = survivors / total

    return survivors, total, density, elapsed


# ============================================================================
# TRUE-Y SURVIVAL TEST
# ============================================================================

def true_y_survives(n: int, y: int, primes):

    """
    Verify directly that the actual hidden y satisfies all QR conditions.
    """

    for prime in primes:

        qr = {
            (x * x) % prime
            for x in range(prime)
        }

        if (n + y * y) % prime not in qr:
            return False

    return True


# ============================================================================
# ACTUAL SEGMENTED FERMAT SEARCH
# ============================================================================

def segmented_fermat(
    n: int,
    primes,
    target_y: int,
    max_y: int,
    block_size: int
):

    """
    Search y=0..max_y using the segmented QR sieve.

    Return:

        found_y
        result
        survivors_tested
        elapsed
    """

    bad_by_prime = {
        p: build_bad_residues(n, p)
        for p in primes
    }

    survivors_tested = 0

    start_time = time.perf_counter()

    block_start = 0

    while block_start <= max_y:

        block_end = min(
            block_start + block_size,
            max_y + 1
        )

        block_len = block_end - block_start

        alive = bytearray(
            b"\x01"
        ) * block_len

        for prime in primes:

            for residue in bad_by_prime[prime]:

                first = (
                    block_start
                    + ((residue - block_start) % prime)
                )

                if first >= block_end:
                    continue

                idx = first - block_start

                count = (
                    (block_len - 1 - idx) // prime
                ) + 1

                alive[idx:block_len:prime] = (
                    b"\x00" * count
                )

        pos = 0

        while True:

            pos = alive.find(1, pos)

            if pos < 0:
                break

            y = block_start + pos

            survivors_tested += 1

            result = fermat_test_y(n, y)

            if result is not None:

                elapsed = time.perf_counter() - start_time

                return (
                    y,
                    result,
                    survivors_tested,
                    elapsed
                )

            pos += 1

        block_start = block_end

    elapsed = time.perf_counter() - start_time

    return (
        None,
        None,
        survivors_tested,
        elapsed
    )


# ============================================================================
# ORDINARY FERMAT
# ============================================================================

def ordinary_fermat(
    n: int,
    max_y: int
):

    start_time = time.perf_counter()

    for y in range(max_y + 1):

        result = fermat_test_y(n, y)

        if result is not None:

            elapsed = time.perf_counter() - start_time

            return (
                y,
                result,
                y + 1,
                elapsed
            )

    elapsed = time.perf_counter() - start_time

    return (
        None,
        None,
        max_y + 1,
        elapsed
    )


# ============================================================================
# DENSITY INSTANCE
# ============================================================================

def run_density_instance(
    bits: int,
    gap_ratio: float
):

    p, q = make_instance(bits, gap_ratio)

    n = p * q

    x = (p + q) // 2
    y = (q - p) // 2

    print()
    print("=" * 72)
    print(f"START DENSITY INSTANCE")
    print("=" * 72)

    print(f"bits(n) = {n.bit_length()}")
    print(f"gap ratio y/p = {gap_ratio}")
    print(f"p = {p}")
    print(f"q = {q}")
    print(f"x = {x}")
    print(f"true y = {y}")
    print(f"gap = {q - p}")
    print(f"density range = {DENSITY_RANGE}")

    for depth in DEPTHS:

        primes = FILTER_PRIMES[:depth]

        survivors, total, density, elapsed = (
            measure_density(
                n,
                primes,
                DENSITY_RANGE - 1,
                BLOCK_SIZE
            )
        )

        survives = true_y_survives(
            n,
            y,
            primes
        )

        expected_survivors = (
            int(round((y + 1) * density))
        )

        print()
        print(f"DEPTH {depth}")
        print(f"    primes = {primes}")
        print(f"    survivors in sample = {survivors}")
        print(f"    sample total = {total}")
        print(f"    density = {density:.12f}")
        print(
            f"    estimated survivors to true y = "
            f"{expected_survivors}"
        )
        print(f"    true y survives = {survives}")
        print(f"    sieve time = {elapsed:.6f}s")

    print()
    print("EXACT IDENTITY")
    print(f"    x²-y² = {x * x - y * y}")
    print(f"    n = {n}")
    print(
        f"    identity = "
        f"{x * x - y * y == n}"
    )

    print()
    print("FINISHED DENSITY INSTANCE")


# ============================================================================
# RECOVERY INSTANCE
# ============================================================================

def run_recovery_instance(
    bits: int,
    gap_ratio: float
):

    p, q = make_instance(bits, gap_ratio)

    n = p * q

    y = (q - p) // 2

    print()
    print("=" * 72)
    print("START RECOVERY INSTANCE")
    print("=" * 72)

    print(f"bits(n) = {n.bit_length()}")
    print(f"gap ratio y/p = {gap_ratio}")
    print(f"p = {p}")
    print(f"q = {q}")
    print(f"true y = {y}")
    print(f"MAX_SEARCH_Y = {MAX_SEARCH_Y}")

    print()
    print("ORDINARY FERMAT")

    ordinary_y, ordinary_result, ordinary_tests, ordinary_time = (
        ordinary_fermat(
            n,
            MAX_SEARCH_Y
        )
    )

    print(f"    found y = {ordinary_y}")
    print(f"    result = {ordinary_result}")
    print(f"    tests = {ordinary_tests}")
    print(f"    time = {ordinary_time:.6f}s")

    for depth in DEPTHS:

        primes = FILTER_PRIMES[:depth]

        print()
        print(f"SEGMENTED DEPTH {depth}")
        print(f"    primes = {primes}")

        start = time.perf_counter()

        found_y, result, survivors, search_time = (
            segmented_fermat(
                n,
                primes,
                y,
                MAX_SEARCH_Y,
                BLOCK_SIZE
            )
        )

        total_time = time.perf_counter() - start

        correct = (
            result == (p, q)
            or result == (q, p)
        )

        print(f"    survivors tested = {survivors}")
        print(f"    found y = {found_y}")
        print(f"    FOUND = {result}")
        print(f"    correct = {correct}")
        print(f"    search time = {search_time:.6f}s")
        print(f"    total time = {total_time:.6f}s")

        if result is not None and ordinary_time > 0:

            print(
                f"    speedup = "
                f"{ordinary_time / total_time:.3f}x"
            )

    print()
    print("EXACT IDENTITY")

    x = (p + q) // 2

    print(f"    x = {x}")
    print(f"    y = {y}")
    print(f"    x²-y² = {x * x - y * y}")
    print(f"    n = {n}")
    print(
        f"    identity = "
        f"{x * x - y * y == n}"
    )

    print()
    print("FINISHED RECOVERY INSTANCE")


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("START EXPERIMENT 191")
    print()
    print("Quadratic-residue Fermat sieve scaling study")
    print()
    print(f"FILTER_PRIMES = {FILTER_PRIMES}")
    print(f"DENSITY_RANGE = {DENSITY_RANGE}")
    print(f"BLOCK_SIZE = {BLOCK_SIZE}")
    print(f"MAX_SEARCH_Y = {MAX_SEARCH_Y}")
    print()

    # ------------------------------------------------------------------------
    # PART A:
    # DENSITY ACROSS BIT SIZES
    #
    # This is cheap even when the true Fermat y is enormous.
    # ------------------------------------------------------------------------

    print()
    print("#" * 72)
    print("PART A — DENSITY SCALING")
    print("#" * 72)

    density_bits = [
        48,
        54,
        60,
        66,
        72,
        78,
    ]

    # Use several factor distances.
    #
    # To keep runtime manageable, each bit size is tested at three distances.
    density_ratios = [
        0.001,
        0.01,
        0.05,
    ]

    for bits in density_bits:

        for gap_ratio in density_ratios:

            run_density_instance(
                bits,
                gap_ratio
            )

    # ------------------------------------------------------------------------
    # PART B:
    # ACTUAL FERMAT RECOVERY
    #
    # Only run instances whose y is deliberately below MAX_SEARCH_Y.
    # ------------------------------------------------------------------------

    print()
    print("#" * 72)
    print("PART B — ACTUAL FERMAT RECOVERY")
    print("#" * 72)

    # These ratios produce controlled y values.
    recovery_ratios = [
        0.00001,
        0.00005,
    ]

    for bits in RECOVERY_BITS:

        for gap_ratio in recovery_ratios:

            run_recovery_instance(
                bits,
                gap_ratio
            )

    # ------------------------------------------------------------------------
    # FINAL
    # ------------------------------------------------------------------------

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 191")
    print("=" * 72)


if __name__ == "__main__":
    main()
