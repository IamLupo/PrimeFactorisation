#!/usr/bin/env python3

import math
import time


# ============================================================================
# EXPERIMENT 189
#
# Correct Fermat quadratic-residue wheel.
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
# A necessary modular condition is:
#
#     n + y^2 must be a quadratic residue modulo every wheel prime.
#
# IMPORTANT:
# We build the wheel directly for y.
# There is NO factor-dependent information in the wheel.
#
# This fixes Experiment 188, where the 4*y^2 scaling was applied incorrectly.
# ============================================================================

R1 = [17, 43, 59, 71, 83]
R2 = [19, 37, 61, 73, 89]

BITS = [30, 36, 42, 48, 54]

FILTER_PRIMES = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
WHEEL_DEPTHS = [1, 2, 3, 4, 5, 6, 7]

MAX_GAP_TESTS = 50_000_000
MAX_WHEEL_CLASSES = 5_000_000


# ============================================================================
# TEST INSTANCES
# ============================================================================

FACTORS = {
    30: (16451, 64271),
    36: (131267, 522517),
    42: (1049089, 4191259),
    48: (8390069, 33547589),
    54: (67112971, 268418039),
}


# ============================================================================
# BASIC FERMAT SEARCH
# ============================================================================

def fermat_test_y(n: int, y: int):
    """
    Test one Fermat y-value.

        x^2 = n + y^2

    If x^2 is a square, return (x-y, x+y).
    """
    value = n + y * y
    x = math.isqrt(value)

    if x * x == value:
        p = x - y
        q = x + y

        if p > 1 and q > 1 and p * q == n:
            return p, q

    return None


def ordinary_fermat(n: int, max_tests: int):
    """
    Standard Fermat search:

        y = 0, 1, 2, ...

    Return:
        (tests, result, elapsed)
    """
    start = time.perf_counter()

    for y in range(max_tests):
        result = fermat_test_y(n, y)

        if result is not None:
            elapsed = time.perf_counter() - start
            return y + 1, result, elapsed

    elapsed = time.perf_counter() - start
    return max_tests, None, elapsed


# ============================================================================
# QUADRATIC RESIDUE TABLE
# ============================================================================

def quadratic_residues(prime: int):
    """
    Return the set of quadratic residues modulo prime.
    """
    return {x * x % prime for x in range(prime)}


# ============================================================================
# BUILD CORRECT Y-WHEEL
# ============================================================================

def build_y_wheel(n: int, primes):
    """
    Construct allowed residues of y modulo M = product(primes).

    A residue y is retained iff

        n + y^2

    is a quadratic residue modulo every wheel prime.

    By the Chinese remainder theorem, testing modulo each prime
    independently is sufficient.
    """
    modulus = 1
    allowed = [0]

    for p in primes:
        residues = quadratic_residues(p)

        new_allowed = []

        for base in allowed:
            for r in range(p):
                candidate = base + modulus * r

                if (n + candidate * candidate) % p in residues:
                    new_allowed.append(candidate)

        modulus *= p
        allowed = new_allowed

        if len(allowed) > MAX_WHEEL_CLASSES:
            return None, None, None

    allowed.sort()

    return modulus, allowed, len(allowed)


# ============================================================================
# WHEEL SEARCH
# ============================================================================

def wheel_fermat(n: int, modulus: int, allowed, max_tests: int):
    """
    Enumerate y values in increasing order using the periodic wheel.

    This avoids the expensive min() / active-progression approach from
    the earlier implementation.
    """
    start = time.perf_counter()

    tests = 0
    y = None

    # We need y values starting at 0.
    period = modulus

    while tests < max_tests:
        for residue in allowed:
            y_candidate = y_base + residue if False else None

        break

    # Explicit period enumeration.
    period_base = 0

    while tests < max_tests:
        for residue in allowed:
            y = period_base + residue

            # The wheel itself is already guaranteed to satisfy the
            # quadratic-residue conditions.
            result = fermat_test_y(n, y)

            tests += 1

            if result is not None:
                elapsed = time.perf_counter() - start
                return tests, result, elapsed, y

            if tests >= max_tests:
                break

        period_base += period

    elapsed = time.perf_counter() - start
    return tests, None, elapsed, None


# ============================================================================
# MAIN
# ============================================================================

def run_instance(bits: int, p: int, q: int):

    n = p * q

    print()
    print("=" * 72)
    print(f"START INSTANCE {bits}-BIT")
    print("=" * 72)

    sqrt_n = math.isqrt(n)

    if sqrt_n * sqrt_n < n:
        sqrt_n += 1

    true_x = (p + q) // 2
    true_y = (q - p) // 2

    print(f"bits(n) = {n.bit_length()}")
    print(f"sqrt(n) = {math.isqrt(n)}")
    print(f"hidden p = {p}")
    print(f"hidden q = {q}")
    print(f"true y = {true_y}")
    print(f"true gap = {q - p}")

    # ------------------------------------------------------------------------
    # ORDINARY FERMAT
    # ------------------------------------------------------------------------

    print()
    print("ORDINARY FERMAT SEARCH")

    ordinary_tests, ordinary_result, ordinary_time = ordinary_fermat(
        n,
        MAX_GAP_TESTS
    )

    print(f"    tests = {ordinary_tests}")
    print(f"    time = {ordinary_time:.6f}s")
    print(f"    found = {ordinary_result}")

    # ------------------------------------------------------------------------
    # WHEELS
    # ------------------------------------------------------------------------

    for depth in WHEEL_DEPTHS:

        primes = FILTER_PRIMES[:depth]

        print()
        print(f"WHEEL DEPTH {depth}")
        print(f"    primes = {primes}")

        build_start = time.perf_counter()

        modulus, allowed, count = build_y_wheel(n, primes)

        build_time = time.perf_counter() - build_start

        if modulus is None:
            print("    wheel build aborted: class limit exceeded")
            continue

        density = count / modulus

        true_residue = true_y % modulus
        true_survives = true_residue in set(allowed)

        print(f"    modulus = {modulus}")
        print(f"    allowed y residues = {count}")
        print(f"    y density = {density:.12f}")
        print(f"    true y residue = {true_residue}")
        print(f"    true y survives = {true_survives}")
        print(f"    wheel build time = {build_time:.6f}s")

        if not true_survives:
            print("    ERROR: true factor y was filtered out")
            continue

        wheel_tests, wheel_result, wheel_time, found_y = wheel_fermat(
            n,
            modulus,
            allowed,
            MAX_GAP_TESTS
        )

        print(f"    wheel tests = {wheel_tests}")
        print(f"    wheel search time = {wheel_time:.6f}s")
        print(f"    found y = {found_y}")
        print(f"    FOUND = {wheel_result}")
        print(f"    correct = {wheel_result == (p, q) or wheel_result == (q, p)}")

        if wheel_result is not None and ordinary_time > 0:
            print(
                f"    speedup vs ordinary = "
                f"{ordinary_time / wheel_time:.3f}x"
            )

    # ------------------------------------------------------------------------
    # EXACT IDENTITY
    # ------------------------------------------------------------------------

    print()
    print("EXACT FERMAT IDENTITY")

    x = (p + q) // 2
    y = (q - p) // 2

    print(f"    x = {x}")
    print(f"    y = {y}")
    print(f"    x²-y² = {x*x - y*y}")
    print(f"    n = {n}")
    print(f"    identity = {(x*x - y*y) == n}")

    print()
    print(f"FINISHED INSTANCE {bits}-BIT")


def main():

    print("START EXPERIMENT 189")
    print()
    print("Correct Fermat quadratic-residue wheel")
    print()

    for bits in BITS:
        p, q = FACTORS[bits]
        run_instance(bits, p, q)

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 189")
    print("=" * 72)


if __name__ == "__main__":
    main()