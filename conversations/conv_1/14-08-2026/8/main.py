#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 51
SUM-DISCRIMINANT SIEVE / CYCLOTOMIC VS RANDOM MODULUS CONTROL
NO CSV OUTPUT
==============================================================================

NEW ALGEBRAIC VARIABLE:

    s = p + q
    d = |p - q|

Then

    s^2 - 4n = d^2

because

    (p+q)^2 - 4pq = (p-q)^2.

For a prime modulus ell, any true factorization must satisfy

    s^2 - 4n  is a quadratic residue mod ell.

Therefore we can sieve possible integer sums s BEFORE doing the
expensive exact perfect-square test.

This experiment compares:

    A. ordinary sum scan
       enumerate every feasible even s

    B. cyclotomic discriminant sieve
       progressively require
           s^2 - 4n
       to be a quadratic residue modulo the non-exceptional
       cyclotomic prime factors of F(r)=r^2+r+1

    C. random-prime control sieve
       same procedure using an equal number of unrelated primes

The purpose is NOT to claim factoring from s alone.

The purpose is to determine whether the cyclotomic family produces
special restrictions on s beyond generic quadratic-residue filtering.

No pair enumeration is performed.
No candidate-pair list is constructed.
"""

from __future__ import annotations

import math
import random
import time
from bisect import bisect_left

# ============================================================================
# CONFIGURATION
# ============================================================================

SEED = 20260814

P_MIN = 2_000_000
P_MAX = 4_200_000

TARGETS = [
    (3318013, 4042603),
    (2129167, 3402323),
    (2224517, 3978749),
    (3685051, 4020281),
    (2399627, 2452649),
    (2593039, 2996527),
    (2149859, 2772097),
    (2060543, 2514401),
    (2675423, 2883973),
    (2828887, 3960137),
    (3497381, 3793241),
    (2193509, 4011353),
]

# Non-exceptional prime factors of F(r)=r^2+r+1 from Experiments 46-50.
CYCLOTOMIC_PRIMES = [
    7, 13, 19, 31, 37, 61, 67, 79, 127, 307, 331, 631, 1723
]

# Equal-size random/control family used in Experiment 48.
CONTROL_PRIMES = [
    673, 1109, 2309, 2351, 4421, 4561, 4759,
    6211, 7879, 7951, 8273, 8689, 9781
]

# Run all moduli by default.
PREFIXES = list(range(1, len(CYCLOTOMIC_PRIMES) + 1))

# ============================================================================
# BASIC UTILITIES
# ============================================================================

def sieve_primes(limit: int) -> list[int]:
    """Return all primes < limit."""
    if limit < 2:
        return []

    sieve = bytearray(b"\x01") * limit
    sieve[0:2] = b"\x00\x00"

    root = math.isqrt(limit - 1)
    for p in range(2, root + 1):
        if sieve[p]:
            start = p * p
            sieve[start:limit:p] = b"\x00" * (
                ((limit - 1 - start) // p) + 1
            )

    return [i for i, flag in enumerate(sieve) if flag]


def closest_prime_order(primes: list[int], n: int) -> list[int]:
    """Prime order by increasing distance to sqrt(n), interleaving sides."""
    root = math.isqrt(n)
    idx = bisect_left(primes, root)

    out = []
    lo = idx - 1
    hi = idx

    while lo >= 0 or hi < len(primes):
        if lo < 0:
            out.append(primes[hi])
            hi += 1
            continue

        if hi >= len(primes):
            out.append(primes[lo])
            lo -= 1
            continue

        dl = abs(primes[lo] - root)
        dh = abs(primes[hi] - root)

        if dl <= dh:
            out.append(primes[lo])
            lo -= 1
        else:
            out.append(primes[hi])
            hi += 1

    return out


def exact_factor_search(
    primes_ordered: list[int],
    n: int,
    p_expected: int,
    q_expected: int,
) -> tuple[int | None, int, float]:
    """Ordinary closest-prime baseline."""
    t0 = time.perf_counter()
    tests = 0

    limit = math.isqrt(n)

    for p in primes_ordered:
        if p > limit:
            break

        tests += 1

        if n % p == 0:
            q = n // p
            ok = (
                (p == p_expected and q == q_expected)
                or
                (p == q_expected and q == p_expected)
            )
            return p, tests, time.perf_counter() - t0

    return None, tests, time.perf_counter() - t0


def quadratic_residue_mask(modulus: int) -> bytearray:
    """
    mask[x] = 1 iff x is a quadratic residue mod modulus.

    Includes zero.
    moduli in this experiment are prime.
    """
    mask = bytearray(modulus)
    for a in range(modulus):
        mask[(a * a) % modulus] = 1
    return mask


def allowed_sum_residues(n: int, ell: int, qr_mask: bytearray) -> set[int]:
    """
    Return residues a mod ell such that

        a^2 - 4n

    is a quadratic residue mod ell.
    """
    n4 = (4 * n) % ell
    allowed = set()

    for a in range(ell):
        delta = (a * a - n4) % ell
        if qr_mask[delta]:
            allowed.add(a)

    return allowed


def sum_range(p_min: int, p_max: int) -> tuple[int, int]:
    """
    Since all primes in the interval are odd:

        s=p+q

    is even, with

        2*p_min <= s < 2*p_max.
    """
    s_min = 2 * p_min
    s_max = 2 * (p_max - 1)

    if s_min % 2:
        s_min += 1
    if s_max % 2:
        s_max -= 1

    return s_min, s_max


def recover_from_sum(
    s: int,
    n: int,
    prime_set: set[int],
) -> tuple[int, int] | None:
    """
    Exact test:

        d^2 = s^2 - 4n

    and

        p=(s-d)/2
        q=(s+d)/2.

    Then require both factors to be primes in the population.
    """
    d2 = s * s - 4 * n
    if d2 < 0:
        return None

    d = math.isqrt(d2)
    if d * d != d2:
        return None

    if (s - d) & 1:
        return None

    p = (s - d) // 2
    q = (s + d) // 2

    if p not in prime_set or q not in prime_set:
        return None

    return (p, q)


def plain_sum_scan(
    n: int,
    s_min: int,
    s_max: int,
    prime_set: set[int],
    expected: tuple[int, int],
) -> tuple[int | None, int, float]:
    """
    Enumerate every feasible even s and perform an exact discriminant test.
    """
    t0 = time.perf_counter()
    exact_tests = 0

    for s in range(s_min, s_max + 1, 2):
        exact_tests += 1
        recovered = recover_from_sum(s, n, prime_set)

        if recovered is not None:
            p, q = recovered
            ep, eq = expected

            if {p, q} == {ep, eq}:
                return s, exact_tests, time.perf_counter() - t0

    return None, exact_tests, time.perf_counter() - t0


def sieve_sum_candidates(
    n: int,
    s_min: int,
    s_max: int,
    moduli: list[int],
    residue_tables: dict[int, tuple[bytearray, set[int]]],
    expected: tuple[int, int],
) -> tuple[
    int | None,
    list[tuple[int, int, int, float]],
    int,
    int,
    float,
]:
    """
    Sequentially filter feasible even sums.

    Returns:
        recovered_s
        stage statistics
        final number of candidates
        exact square tests
        total filtering+exact time
    """
    t0 = time.perf_counter()

    candidates = list(range(s_min, s_max + 1, 2))

    stage_rows = []

    for stage, ell in enumerate(moduli, 1):
        before = len(candidates)
        allowed = residue_tables[ell][1]

        # Keep only s whose residue is compatible with:
        #
        #     s^2 - 4n  quadratic residue mod ell.
        candidates = [s for s in candidates if (s % ell) in allowed]

        after = len(candidates)
        retained = after / before if before else 0.0
        reduction = before / after if after else float("inf")

        stage_rows.append(
            (
                stage,
                ell,
                before,
                after,
                retained,
                reduction,
                time.perf_counter() - t0,
            )
        )

    exact_tests = 0
    recovered_s = None

    for s in candidates:
        exact_tests += 1

        recovered = recover_from_sum(
            s,
            n,
            GLOBAL_PRIME_SET,
        )

        if recovered is None:
            continue

        if set(recovered) == set(expected):
            recovered_s = s
            break

    total = time.perf_counter() - t0

    return (
        recovered_s,
        stage_rows,
        len(candidates),
        exact_tests,
        total,
    )


# ============================================================================
# GLOBAL PRIME POPULATION
# ============================================================================

GLOBAL_PRIME_SET: set[int] = set()


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    global GLOBAL_PRIME_SET

    random.seed(SEED)

    print("=" * 78)
    print("KAPPA EXPERIMENT 51")
    print("SUM-DISCRIMINANT SIEVE / CYCLOTOMIC VS RANDOM CONTROL")
    print("NO CSV OUTPUT")
    print("=" * 78)

    print()
    print("random seed      =", SEED)
    print("targets          =", len(TARGETS))
    print(f"prime interval   = [{P_MIN:,}, {P_MAX:,})")
    print()

    # ------------------------------------------------------------------------
    # PRIME POPULATION
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("1. PRIME POPULATION")
    print("=" * 78)

    t0 = time.perf_counter()
    primes = sieve_primes(P_MAX)
    primes = [p for p in primes if p >= P_MIN]

    GLOBAL_PRIME_SET = set(primes)

    generation_time = time.perf_counter() - t0

    print(f"prime population = {len(primes):,}")
    print(f"generation time  = {generation_time:.4f}s")

    # ------------------------------------------------------------------------
    # TARGETS
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. TARGETS")
    print("=" * 78)

    for i, (p, q) in enumerate(TARGETS, 1):
        print(f"target {i:2d}: p={p} q={q} n={p*q}")

    # ------------------------------------------------------------------------
    # SUM RANGE
    # ------------------------------------------------------------------------

    s_min, s_max = sum_range(P_MIN, P_MAX)
    total_sum_population = ((s_max - s_min) // 2) + 1

    print()
    print("=" * 78)
    print("3. FEASIBLE SUM RANGE")
    print("=" * 78)

    print(f"s minimum       = {s_min:,}")
    print(f"s maximum       = {s_max:,}")
    print(f"even sums       = {total_sum_population:,}")
    print()
    print("For every target:")
    print()
    print("    d^2 = s^2 - 4n")
    print("    p   = (s-d)/2")
    print("    q   = (s+d)/2")

    # ------------------------------------------------------------------------
    # MODULI
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. MODULUS FAMILIES")
    print("=" * 78)

    print("cyclotomic =", CYCLOTOMIC_PRIMES)
    print("control    =", CONTROL_PRIMES)

    # ------------------------------------------------------------------------
    # PRECOMPUTE QUADRATIC RESIDUES
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. QUADRATIC-RESIDUE TABLES")
    print("=" * 78)

    all_moduli = sorted(set(CYCLOTOMIC_PRIMES + CONTROL_PRIMES))

    residue_tables: dict[int, tuple[bytearray, set[int]]] = {}

    t0 = time.perf_counter()

    for ell in all_moduli:
        qr_mask = quadratic_residue_mask(ell)

        # The allowed s residues depend on n, so only the QR mask
        # is precomputed globally.
        residue_tables[ell] = (qr_mask, set())

    print(
        "precomputed QR tables for",
        len(all_moduli),
        "prime moduli in",
        f"{time.perf_counter() - t0:.4f}s",
    )

    # ------------------------------------------------------------------------
    # BASELINE
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. BASELINE COMPARISON")
    print("=" * 78)

    baseline_rows = []

    total_plain_exact = 0
    total_plain_time = 0.0

    for i, (p, q) in enumerate(TARGETS, 1):
        n = p * q

        recovered_s, exact_tests, elapsed = plain_sum_scan(
            n,
            s_min,
            s_max,
            GLOBAL_PRIME_SET,
            (p, q),
        )

        total_plain_exact += exact_tests
        total_plain_time += elapsed

        correct = recovered_s == p + q

        baseline_rows.append(
            (
                recovered_s,
                exact_tests,
                elapsed,
                correct,
            )
        )

        print(
            f"target {i:2d}: "
            f"exact sum tests={exact_tests:8,d} "
            f"time={elapsed:.6f}s "
            f"correct={correct}"
        )

    print()
    print(f"baseline total exact tests = {total_plain_exact:,}")
    print(f"baseline total time         = {total_plain_time:.6f}s")

    # ------------------------------------------------------------------------
    # CYCLOTOMIC SIEVE
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. CYCLOTOMIC SUM-DISCRIMINANT SIEVE")
    print("=" * 78)

    cyclo_summary = []

    total_cyclo_time = 0.0
    total_cyclo_exact = 0
    recovered_count = 0

    for i, (p, q) in enumerate(TARGETS, 1):
        n = p * q

        # Build target-specific allowed residue sets.
        target_tables = {}

        for ell in CYCLOTOMIC_PRIMES:
            qr_mask = residue_tables[ell][0]
            allowed = allowed_sum_residues(n, ell, qr_mask)
            target_tables[ell] = (qr_mask, allowed)

        (
            recovered_s,
            stage_rows,
            final_candidates,
            exact_tests,
            elapsed,
        ) = sieve_sum_candidates(
            n,
            s_min,
            s_max,
            CYCLOTOMIC_PRIMES,
            target_tables,
            (p, q),
        )

        total_cyclo_time += elapsed
        total_cyclo_exact += exact_tests

        correct = recovered_s == p + q
        recovered_count += int(correct)

        cyclo_summary.append(
            (
                final_candidates,
                exact_tests,
                elapsed,
                correct,
            )
        )

        print()
        print("-" * 78)
        print(f"TARGET {i:2d}  n={n}")
        print("-" * 78)

        print(
            "stage modulus before after retained reduction"
        )

        for (
            stage,
            ell,
            before,
            after,
            retained,
            reduction,
            _,
        ) in stage_rows:
            print(
                f"{stage:5d} "
                f"{ell:8d} "
                f"{before:8,d} "
                f"{after:8,d} "
                f"{retained:9.6f} "
                f"{reduction:10.3f}x"
            )

        print(
            f"final candidates = {final_candidates:,}"
        )
        print(
            f"exact square tests = {exact_tests:,}"
        )
        print(
            f"recovered s       = {recovered_s}"
        )
        print(
            f"true s            = {p+q}"
        )
        print(
            f"correct           = {correct}"
        )
        print(
            f"total time        = {elapsed:.6f}s"
        )

    # ------------------------------------------------------------------------
    # CONTROL SIEVE
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. RANDOM-PRIME CONTROL SUM-DISCRIMINANT SIEVE")
    print("=" * 78)

    control_summary = []

    total_control_time = 0.0
    total_control_exact = 0
    control_recovered = 0

    for i, (p, q) in enumerate(TARGETS, 1):
        n = p * q

        target_tables = {}

        for ell in CONTROL_PRIMES:
            qr_mask = residue_tables[ell][0]
            allowed = allowed_sum_residues(n, ell, qr_mask)
            target_tables[ell] = (qr_mask, allowed)

        (
            recovered_s,
            stage_rows,
            final_candidates,
            exact_tests,
            elapsed,
        ) = sieve_sum_candidates(
            n,
            s_min,
            s_max,
            CONTROL_PRIMES,
            target_tables,
            (p, q),
        )

        total_control_time += elapsed
        total_control_exact += exact_tests

        correct = recovered_s == p + q
        control_recovered += int(correct)

        control_summary.append(
            (
                final_candidates,
                exact_tests,
                elapsed,
                correct,
            )
        )

        print()
        print("-" * 78)
        print(f"TARGET {i:2d}  n={n}")
        print("-" * 78)

        print(
            "stage modulus before after retained reduction"
        )

        for (
            stage,
            ell,
            before,
            after,
            retained,
            reduction,
            _,
        ) in stage_rows:
            print(
                f"{stage:5d} "
                f"{ell:8d} "
                f"{before:8,d} "
                f"{after:8,d} "
                f"{retained:9.6f} "
                f"{reduction:10.3f}x"
            )

        print(
            f"final candidates = {final_candidates:,}"
        )
        print(
            f"exact square tests = {exact_tests:,}"
        )
        print(
            f"recovered s       = {recovered_s}"
        )
        print(
            f"true s            = {p+q}"
        )
        print(
            f"correct           = {correct}"
        )
        print(
            f"total time        = {elapsed:.6f}s"
        )

    # ------------------------------------------------------------------------
    # SUMMARY STATISTICS
    # ------------------------------------------------------------------------

    def median(values: list[float]) -> float:
        values = sorted(values)
        n = len(values)
        if n % 2:
            return values[n // 2]
        return (values[n // 2 - 1] + values[n // 2]) / 2.0

    baseline_tests = [float(x[1]) for x in baseline_rows]
    baseline_times = [float(x[2]) for x in baseline_rows]

    cyclo_candidates = [float(x[0]) for x in cyclo_summary]
    cyclo_exact = [float(x[1]) for x in cyclo_summary]
    cyclo_times = [float(x[2]) for x in cyclo_summary]

    control_candidates = [float(x[0]) for x in control_summary]
    control_exact = [float(x[1]) for x in control_summary]
    control_times = [float(x[2]) for x in control_summary]

    # ------------------------------------------------------------------------
    # FINAL SUMMARY
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. GLOBAL SUMMARY")
    print("=" * 78)

    print()
    print(
        "method                         "
        "median final   median exact   median time"
    )
    print("-" * 78)

    print(
        f"{'plain sum scan':30s} "
        f"{total_sum_population:12,d} "
        f"{median(baseline_tests):14,.1f} "
        f"{median(baseline_times):.6f}s"
    )

    print(
        f"{'cyclotomic sieve':30s} "
        f"{median(cyclo_candidates):12,.1f} "
        f"{median(cyclo_exact):14,.1f} "
        f"{median(cyclo_times):.6f}s"
    )

    print(
        f"{'random control sieve':30s} "
        f"{median(control_candidates):12,.1f} "
        f"{median(control_exact):14,.1f} "
        f"{median(control_times):.6f}s"
    )

    print()
    print(
        "cyclotomic final-candidate reduction = "
        f"{total_sum_population / median(cyclo_candidates):.3f}x"
    )
    print(
        "control final-candidate reduction     = "
        f"{total_sum_population / median(control_candidates):.3f}x"
    )

    print()
    print(
        "cyclotomic exact-test reduction = "
        f"{median(baseline_tests) / median(cyclo_exact):.3f}x"
    )
    print(
        "control exact-test reduction    = "
        f"{median(baseline_tests) / median(control_exact):.3f}x"
    )

    print()
    print(
        "cyclotomic recovery = "
        f"{recovered_count}/{len(TARGETS)}"
    )
    print(
        "control recovery    = "
        f"{control_recovered}/{len(TARGETS)}"
    )

    # ------------------------------------------------------------------------
    # ALGEBRAIC DIAGNOSTIC
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. ALGEBRAIC DIAGNOSTIC")
    print("=" * 78)

    print()
    print("The tested necessary condition is:")
    print()
    print("    d^2 = s^2 - 4n")
    print()
    print("For every modulus ell:")
    print()
    print("    s^2 - 4n must be a quadratic residue mod ell.")
    print()
    print("The important comparison is:")
    print()
    print("    CYCLOTOMIC:")
    print("        ell | r^2+r+1")
    print("        ell != 3")
    print("        r^3 = 1 mod ell")
    print()
    print("    CONTROL:")
    print("        unrelated prime moduli")
    print()
    print("Interpretation:")
    print()
    print("    If cyclotomic and control produce similar survivor")
    print("    densities, the useful effect is probably generic")
    print("    quadratic-residue filtering.")
    print()
    print("    If cyclotomic consistently produces materially fewer")
    print("    surviving sums than the control family, then the")
    print("    cubic-root structure is contributing additional")
    print("    information about the feasible sum s=p+q.")
    print()
    print("    If exact square tests collapse but runtime does not,")
    print("    the mathematical filter is strong but the implementation")
    print("    still needs optimization.")
    print()
    print("No prime-pair enumeration is performed.")
    print("No CSV files are produced.")

    print()
    print("=" * 78)
    print("EXPERIMENT 51 COMPLETE")
    print("=" * 78)

    total_runtime = generation_time + total_plain_time + \
        total_cyclo_time + total_control_time

    print(f"approx measured runtime = {total_runtime:.4f}s")


if __name__ == "__main__":
    main()

