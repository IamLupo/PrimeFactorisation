#!/usr/bin/env python3

import math
import time


# ============================================================================
# START EXPERIMENT 192
#
# ANALYTICAL QUADRATIC-RESIDUE DENSITY + PRIME-SELECTION OPTIMIZATION
#
# Experiment 191 showed that the QR-sieve density is highly structured.
#
# For an odd prime l not dividing n, the Fermat condition is
#
#     n + y^2 = x^2 (mod l)
#
# so n + y^2 must be a quadratic residue modulo l.
#
# Let chi_l(-n) be the Legendre symbol.
#
# The exact number of allowed y residues modulo l is:
#
#     A_l = (l + chi_l(-n)) / 2
#
# hence exact density:
#
#     rho_l = A_l / l
#
# and for a set S of distinct primes:
#
#     rho(S) = product rho_l
#
# because the congruence conditions combine by CRT.
#
# This experiment does FOUR things:
#
#   PART A
#       Compute the analytical density prime-by-prime.
#
#   PART B
#       Compare analytical density against measured segmented density.
#
#   PART C
#       Automatically select primes using a cost/value optimization.
#
#   PART D
#       Benchmark the automatically selected wheels against fixed wheels.
#
# The hidden factors are used ONLY to generate/verify instances.
# The QR conditions themselves never use the factors.
# ============================================================================


FILTER_PRIMES = [
    3, 5, 7, 11, 13, 17, 19, 23,
    29, 31, 37, 41, 43, 47, 53, 59
]

# Prime pool used for automatic optimization.
OPT_PRIMES = FILTER_PRIMES

# Analytical / empirical density sample.
DENSITY_RANGE = 1 << 22

# Segmented sieve block.
BLOCK_SIZE = 1 << 20

# Actual search limit.
MAX_SEARCH_Y = 25_000_000

# Budget is expressed in approximate "marking cost".
#
# A prime p requires approximately O(p) work to construct its bad-residue
# table, and then O(number of bad residue classes) marking work per block.
#
# We use p as a simple conservative proxy.
MAX_WHEEL_COST = 120

# Benchmark the automatically selected wheel plus common fixed wheels.
FIXED_WHEELS = [
    [3],
    [3, 5],
    [3, 5, 7],
    [3, 5, 7, 11],
    [3, 5, 7, 11, 13],
    [3, 5, 7, 11, 13, 17],
    [3, 5, 7, 11, 13, 17, 19],
]

# Existing controlled instances from previous experiments.
INSTANCES = {
    48: (8390069, 33547589),
    54: (124517461, 144517463),
    60: (1058841403, 1088841421),
    66: (8569934017, 8609934041),
}


# ============================================================================
# MILLER-RABIN
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
        23, 29, 31, 37, 41, 43, 47
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


# ============================================================================
# LEGENDRE SYMBOL
# ============================================================================

def legendre_symbol(a: int, p: int) -> int:
    """
    Return:

        +1 if a is a nonzero quadratic residue mod p
        -1 if a is a non-residue mod p
         0 if a == 0 mod p
    """

    a %= p

    if a == 0:
        return 0

    value = pow(a, (p - 1) // 2, p)

    if value == 1:
        return 1

    if value == p - 1:
        return -1

    raise RuntimeError(
        f"Unexpected Legendre-symbol value {value}"
    )


# ============================================================================
# ANALYTICAL PRIME DENSITY
# ============================================================================

def prime_density(n: int, p: int):
    """
    Exact density of allowed y modulo p.

    For p not dividing n:

        A = (p + chi(-n)) / 2

    where chi is the Legendre symbol.

    If p divides n, the formula has the chi=0 special case and the exact
    count is handled explicitly.
    """

    chi = legendre_symbol(-n, p)

    if n % p == 0:
        # Then x^2 == y^2 mod p.
        # Exactly the solutions y == +/- x mod p are represented, giving
        # (p+1)/2 distinct y values because the zero class is shared.
        allowed = (p + 1) // 2

    else:
        allowed = (p + chi) // 2

    density = allowed / p

    return chi, allowed, density


# ============================================================================
# EXACT PRODUCT DENSITY
# ============================================================================

def analytical_density(n: int, primes):

    density = 1.0
    rows = []

    for p in primes:

        chi, allowed, rho = prime_density(n, p)

        density *= rho

        rows.append(
            {
                "prime": p,
                "chi": chi,
                "allowed": allowed,
                "density": rho,
            }
        )

    return density, rows


# ============================================================================
# BAD RESIDUES FOR EMPIRICAL CHECK
# ============================================================================

def build_bad_residues(n: int, p: int):

    quadratic_residues = {
        (x * x) % p
        for x in range(p)
    }

    bad = []

    for r in range(p):

        if (n + r * r) % p not in quadratic_residues:
            bad.append(r)

    return bad


# ============================================================================
# EMPIRICAL DENSITY
# ============================================================================

def measured_density(
    n: int,
    primes,
    max_y: int,
    block_size: int
):

    bad_by_prime = {
        p: build_bad_residues(n, p)
        for p in primes
    }

    survivors = 0
    total = max_y + 1

    start = time.perf_counter()

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

        for p in primes:

            for residue in bad_by_prime[p]:

                first = (
                    block_start
                    + ((residue - block_start) % p)
                )

                if first >= block_end:
                    continue

                idx = first - block_start

                count = (
                    (block_len - 1 - idx) // p
                ) + 1

                alive[idx:block_len:p] = (
                    b"\x00" * count
                )

        survivors += alive.count(1)

        block_start = block_end

    elapsed = time.perf_counter() - start

    return survivors, total, survivors / total, elapsed


# ============================================================================
# AUTOMATIC PRIME OPTIMIZATION
# ============================================================================

def optimize_primes(n: int, candidate_primes, budget: int):
    """
    0/1 knapsack.

    Each prime has:

        cost  = p

        value = -log2(rho_p)

    where rho_p is the exact surviving density.

    Thus value measures the number of bits of candidate reduction gained.

    The optimizer maximizes total modular reduction subject to the cost
    budget.
    """

    items = []

    for p in candidate_primes:

        chi, allowed, rho = prime_density(n, p)

        if rho <= 0:
            continue

        value = -math.log2(rho)

        cost = p

        items.append(
            (
                p,
                cost,
                value,
                rho,
                chi,
                allowed,
            )
        )

    # dp[c] = (value, selected primes)
    dp = [
        (0.0, [])
        for _ in range(budget + 1)
    ]

    for item in items:

        p, cost, value, rho, chi, allowed = item

        if cost > budget:
            continue

        for c in range(budget, cost - 1, -1):

            old_value, old_set = dp[c - cost]

            new_value = old_value + value

            if new_value > dp[c][0]:

                dp[c] = (
                    new_value,
                    old_set + [p]
                )

    best_value = -1.0
    best_set = []

    for value, selected in dp:

        if value > best_value:

            best_value = value
            best_set = selected

    best_set.sort()

    density, _ = analytical_density(
        n,
        best_set
    )

    return best_set, density, best_value


# ============================================================================
# GREEDY PRIME OPTIMIZER
# ============================================================================

def greedy_primes(n: int, candidate_primes, budget: int):

    remaining = budget
    selected = []

    candidates = []

    for p in candidate_primes:

        chi, allowed, rho = prime_density(n, p)

        # Score = modular reduction per estimated construction cost.
        score = (-math.log2(rho)) / p

        candidates.append(
            (
                score,
                p,
                rho,
                chi,
                allowed,
            )
        )

    candidates.sort(reverse=True)

    for score, p, rho, chi, allowed in candidates:

        if p > remaining:
            continue

        selected.append(p)
        remaining -= p

    selected.sort()

    density, _ = analytical_density(
        n,
        selected
    )

    return selected, density


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
# SEGMENTED FERMAT SEARCH
# ============================================================================

def segmented_fermat(
    n: int,
    primes,
    max_y: int,
    block_size: int
):

    bad_by_prime = {
        p: build_bad_residues(n, p)
        for p in primes
    }

    survivors = 0

    start = time.perf_counter()

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

        for p in primes:

            for residue in bad_by_prime[p]:

                first = (
                    block_start
                    + ((residue - block_start) % p)
                )

                if first >= block_end:
                    continue

                idx = first - block_start

                count = (
                    (block_len - 1 - idx) // p
                ) + 1

                alive[idx:block_len:p] = (
                    b"\x00" * count
                )

        pos = 0

        while True:

            pos = alive.find(1, pos)

            if pos < 0:
                break

            y = block_start + pos

            survivors += 1

            result = fermat_test_y(n, y)

            if result is not None:

                elapsed = time.perf_counter() - start

                return (
                    y,
                    result,
                    survivors,
                    elapsed
                )

            pos += 1

        block_start = block_end

    elapsed = time.perf_counter() - start

    return (
        None,
        None,
        survivors,
        elapsed
    )


# ============================================================================
# PART A — ANALYTICAL DENSITY
# ============================================================================

def run_analytical_density(n: int):

    print()
    print("=" * 72)
    print("PART A — ANALYTICAL PRIME DENSITIES")
    print("=" * 72)

    for p in OPT_PRIMES:

        chi, allowed, rho = prime_density(n, p)

        bits_reduction = -math.log2(rho)

        print(
            f"prime={p:2d} "
            f"chi(-n)={chi:+d} "
            f"allowed={allowed:2d}/{p:<2d} "
            f"rho={rho:.12f} "
            f"reduction={bits_reduction:.6f} bits"
        )


# ============================================================================
# PART B — ANALYTICAL VS EMPIRICAL
# ============================================================================

def run_density_comparison(n: int):

    print()
    print("=" * 72)
    print("PART B — ANALYTICAL VS EMPIRICAL")
    print("=" * 72)

    for depth in range(1, len(FILTER_PRIMES) + 1):

        primes = FILTER_PRIMES[:depth]

        analytical, _ = analytical_density(
            n,
            primes
        )

        survivors, total, measured, elapsed = (
            measured_density(
                n,
                primes,
                DENSITY_RANGE - 1,
                BLOCK_SIZE
            )
        )

        error = measured - analytical

        print()
        print(f"DEPTH {depth}")
        print(f"    primes = {primes}")
        print(f"    analytical density = {analytical:.12f}")
        print(f"    measured density   = {measured:.12f}")
        print(f"    absolute error     = {error:+.12e}")
        print(f"    survivors          = {survivors}")
        print(f"    sample total       = {total}")
        print(f"    measurement time   = {elapsed:.6f}s")


# ============================================================================
# PART C — AUTOMATIC OPTIMIZATION
# ============================================================================

def run_optimization(n: int):

    print()
    print("=" * 72)
    print("PART C — AUTOMATIC PRIME SELECTION")
    print("=" * 72)

    print(f"candidate primes = {OPT_PRIMES}")
    print(f"cost budget = {MAX_WHEEL_COST}")

    optimal, optimal_density, optimal_bits = optimize_primes(
        n,
        OPT_PRIMES,
        MAX_WHEEL_COST
    )

    greedy, greedy_density = greedy_primes(
        n,
        OPT_PRIMES,
        MAX_WHEEL_COST
    )

    print()
    print("KNAPSACK OPTIMUM")
    print(f"    primes = {optimal}")
    print(f"    total cost = {sum(optimal)}")
    print(f"    analytical density = {optimal_density:.12f}")
    print(f"    reduction = {optimal_bits:.6f} bits")
    print(f"    reduction factor = {1.0 / optimal_density:.3f}x")

    print()
    print("GREEDY")
    print(f"    primes = {greedy}")
    print(f"    total cost = {sum(greedy)}")
    print(f"    analytical density = {greedy_density:.12f}")
    print(
        f"    reduction = "
        f"{-math.log2(greedy_density):.6f} bits"
    )
    print(
        f"    reduction factor = "
        f"{1.0 / greedy_density:.3f}x"
    )

    return optimal, greedy


# ============================================================================
# PART D — BENCHMARK
# ============================================================================

def run_benchmark(
    bits: int,
    p: int,
    q: int
):

    n = p * q
    true_y = (q - p) // 2

    print()
    print("=" * 72)
    print(f"START BENCHMARK {bits}-BIT")
    print("=" * 72)

    print(f"p = {p}")
    print(f"q = {q}")
    print(f"true y = {true_y}")

    if true_y > MAX_SEARCH_Y:

        print(
            "    SKIPPED: true y exceeds MAX_SEARCH_Y"
        )

        print(
            f"FINISHED BENCHMARK {bits}-BIT"
        )

        return

    # ------------------------------------------------------------------------
    # Ordinary reference
    # ------------------------------------------------------------------------

    start = time.perf_counter()

    for y in range(true_y + 1):

        result = fermat_test_y(n, y)

        if result is not None:
            break

    ordinary_time = time.perf_counter() - start

    print()
    print("ORDINARY")

    print(f"    tests = {true_y + 1}")
    print(f"    time = {ordinary_time:.6f}s")
    print(f"    result = {result}")

    # ------------------------------------------------------------------------
    # Build automatic sets
    # ------------------------------------------------------------------------

    optimal, optimal_density, _ = optimize_primes(
        n,
        OPT_PRIMES,
        MAX_WHEEL_COST
    )

    greedy, greedy_density = greedy_primes(
        n,
        OPT_PRIMES,
        MAX_WHEEL_COST
    )

    benchmark_wheels = []

    for wheel in FIXED_WHEELS:

        if wheel not in benchmark_wheels:
            benchmark_wheels.append(wheel)

    if optimal not in benchmark_wheels:
        benchmark_wheels.append(optimal)

    if greedy not in benchmark_wheels:
        benchmark_wheels.append(greedy)

    # ------------------------------------------------------------------------
    # Bench wheels
    # ------------------------------------------------------------------------

    for wheel in benchmark_wheels:

        analytical, _ = analytical_density(
            n,
            wheel
        )

        start = time.perf_counter()

        found_y, found, survivors, search_time = (
            segmented_fermat(
                n,
                wheel,
                MAX_SEARCH_Y,
                BLOCK_SIZE
            )
        )

        total_time = time.perf_counter() - start

        correct = (
            found == (p, q)
            or found == (q, p)
        )

        print()
        print(f"WHEEL {wheel}")
        print(f"    cost = {sum(wheel)}")
        print(f"    analytical density = {analytical:.12f}")
        print(
            f"    predicted survivor count = "
            f"{int((true_y + 1) * analytical)}"
        )
        print(f"    survivors tested = {survivors}")
        print(f"    found y = {found_y}")
        print(f"    result = {found}")
        print(f"    correct = {correct}")
        print(f"    time = {total_time:.6f}s")

        if ordinary_time > 0:

            print(
                f"    speedup = "
                f"{ordinary_time / total_time:.3f}x"
            )

    print()
    print(f"OPTIMUM SET = {optimal}")
    print(f"OPTIMUM DENSITY = {optimal_density:.12f}")
    print(f"GREEDY SET = {greedy}")
    print(f"GREEDY DENSITY = {greedy_density:.12f}")

    print()
    print(f"FINISHED BENCHMARK {bits}-BIT")


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("START EXPERIMENT 192")
    print()
    print("Analytical QR density + automatic prime optimization")
    print()
    print(f"OPT_PRIMES = {OPT_PRIMES}")
    print(f"DENSITY_RANGE = {DENSITY_RANGE}")
    print(f"BLOCK_SIZE = {BLOCK_SIZE}")
    print(f"MAX_SEARCH_Y = {MAX_SEARCH_Y}")
    print(f"MAX_WHEEL_COST = {MAX_WHEEL_COST}")

    # ------------------------------------------------------------------------
    # Use the 48-bit instance as the detailed analytical reference.
    # ------------------------------------------------------------------------

    reference_p, reference_q = INSTANCES[48]
    reference_n = reference_p * reference_q

    print()
    print("=" * 72)
    print("REFERENCE INSTANCE")
    print("=" * 72)

    print(f"p = {reference_p}")
    print(f"q = {reference_q}")
    print(f"n = {reference_n}")
    print(
        f"true y = "
        f"{(reference_q - reference_p) // 2}"
    )

    run_analytical_density(reference_n)
    run_density_comparison(reference_n)

    optimal, greedy = run_optimization(
        reference_n
    )

    # ------------------------------------------------------------------------
    # Benchmark multiple bit sizes.
    #
    # Automatic selection is recomputed for EVERY n because the Legendre
    # symbols depend on n mod p.
    # ------------------------------------------------------------------------

    for bits in [48, 54, 60, 66]:

        p, q = INSTANCES[bits]

        run_benchmark(
            bits,
            p,
            q
        )

    print()
    print("=" * 72)
    print("FINISHED EXPERIMENT 192")
    print("=" * 72)


if __name__ == "__main__":
    main()
