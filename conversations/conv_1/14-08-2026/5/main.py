#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 48
CYCLOTOMIC SIGNATURE VS RANDOM-PRIME MODULUS CONTROL
NO CSV OUTPUT
==============================================================================

QUESTION
--------
Is the useful filtering caused by the special structure

    F(r) = r^2 + r + 1 = Phi_3(r)

or is it simply the accumulation of ordinary modular constraints?

This experiment compares:

    A) cyclotomic prime factors of F(r)
    B) unrelated control primes

For each cumulative stage it measures:

    - combined modulus M
    - phi(M)
    - observed compatible-prime fraction
    - theoretical 1/phi(M)
    - ratio observed / (1/phi(M))
    - exact divisibility work
    - correctness

NO prime-pair enumeration.
NO candidate-pair lists.
NO CSV OUTPUT.

Requires:
    sympy
"""

from __future__ import annotations

import math
import random
import statistics
import time
from dataclasses import dataclass
from typing import Iterable

import sympy


# ============================================================================
# CONFIGURATION
# ============================================================================

SEED = 20260814

PRIME_LO = 2_000_000
PRIME_HI = 4_200_000

TARGETS = 12

# Use enough stages to see accumulation without making the experiment huge.
CONTROL_MODULI_COUNT = 13

# Cyclotomic r values whose F(r) prime factors were identified previously.
R_VALUES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]

# We exclude ell=3 from the main cyclotomic comparison because it is the
# exceptional factor with r == 1 (mod 3), order 1 rather than order 3.
EXCLUDE_EXCEPTIONAL_3 = True

# Exact-search benchmark only for selected late stages.
EXACT_BENCHMARK_STAGES = {3, 5, 8, 10, 13}


# ============================================================================
# FORMATTING
# ============================================================================

def line(ch: str = "-", n: int = 78) -> None:
    print(ch * n)


def fmt_int(x: int) -> str:
    return f"{x:,}"


def fmt_float(x: float, digits: int = 9) -> str:
    return f"{x:.{digits}g}"


# ============================================================================
# PRIME POPULATION
# ============================================================================

def generate_prime_population() -> list[int]:
    """
    Generate every prime in the chosen interval.

    Uses sympy.isprime as requested.
    """
    return list(sympy.primerange(PRIME_LO, PRIME_HI))


# ============================================================================
# TARGET GENERATION
# ============================================================================

def generate_targets(
    primes: list[int],
    count: int,
    rng: random.Random,
) -> list[tuple[int, int, int]]:
    """
    Generate distinct semiprime targets from the prime population.

    The two factors are chosen from the same population.
    """
    targets = []
    seen = set()

    while len(targets) < count:
        p = rng.choice(primes)
        q = rng.choice(primes)

        if p == q:
            continue

        a, b = sorted((p, q))
        if (a, b) in seen:
            continue

        seen.add((a, b))
        targets.append((a, b, a * b))

    return targets


# ============================================================================
# CYCLOTOMIC FAMILY
# ============================================================================

def F(r: int) -> int:
    return r * r + r + 1


def cyclotomic_prime_inventory(
    r_values: Iterable[int],
) -> list[int]:
    """
    Collect distinct prime factors ell of F(r).
    """
    factors = set()

    for r in r_values:
        value = F(r)
        factorization = sympy.factorint(value)

        for ell in factorization:
            ell = int(ell)

            if EXCLUDE_EXCEPTIONAL_3 and ell == 3:
                continue

            factors.add(ell)

    return sorted(factors)


# ============================================================================
# CONTROL MODULI
# ============================================================================

def choose_control_primes(
    exclude: set[int],
    count: int,
    rng: random.Random,
) -> list[int]:
    """
    Choose unrelated prime moduli.

    We deliberately choose primes not appearing in the cyclotomic family.

    To make the comparison reasonably fair, the control primes are chosen
    from a broad logarithmic scale rather than all being tiny primes.
    """
    candidates = [
        p
        for p in list(sympy.primerange(5, 10000))
        if p not in exclude
    ]

    rng.shuffle(candidates)

    selected = []

    for p in candidates:
        # Avoid accidental repeated small structure by spacing selections.
        if all(math.gcd(p, q) == 1 for q in selected):
            selected.append(p)

        if len(selected) == count:
            break

    if len(selected) < count:
        raise RuntimeError("Could not construct control prime family.")

    return sorted(selected)


# ============================================================================
# MODULUS / PHI
# ============================================================================

def lcm_list(values: list[int]) -> int:
    result = 1

    for value in values:
        result = math.lcm(result, value)

    return result


def phi_from_factorization(moduli: list[int]) -> int:
    """
    phi(LCM(moduli)), computed exactly from the distinct prime factors.
    """
    M = lcm_list(moduli)
    fac = sympy.factorint(M)

    phi = M

    for p in fac:
        phi = phi // p * (p - 1)

    return int(phi)


# ============================================================================
# RESIDUE SIGNATURE INDEX
# ============================================================================

@dataclass
class SignatureIndex:
    moduli: tuple[int, ...]
    signatures: set[tuple[int, ...]]
    M: int
    phi: int


def build_signature_index(
    primes: list[int],
    moduli: list[int],
) -> SignatureIndex:
    """
    Build the set of occupied prime residue signatures.

    For moduli [m1, ..., mk], a prime q receives:

        (q mod m1, ..., q mod mk)

    Only the signature itself is stored, not the primes belonging to it.

    This keeps memory modest.
    """
    signatures = set()

    for p in primes:
        signatures.add(tuple(p % m for m in moduli))

    M = lcm_list(moduli)
    phi = phi_from_factorization(moduli)

    return SignatureIndex(
        moduli=tuple(moduli),
        signatures=signatures,
        M=M,
        phi=phi,
    )


# ============================================================================
# N-ONLY COMPATIBILITY
# ============================================================================

def compatible_prime_count(
    primes: list[int],
    n: int,
    index: SignatureIndex,
) -> tuple[int, float]:
    """
    Count primes p for which there exists a prime q in the population
    satisfying

        p*q == n (mod every selected modulus).

    For each modulus m:

        q == n * p^(-1) mod m

    because every searchable prime is coprime to these small moduli.

    The implementation checks the target q-signature against the occupied
    prime-signature set.
    """
    moduli = index.moduli
    signatures = index.signatures

    survivors = 0

    for p in primes:
        wanted = []

        valid = True

        for m in moduli:
            g = math.gcd(p, m)

            if g != 1:
                # Generalized condition: if p is not invertible, fall back
                # to direct residue compatibility.
                #
                # This branch is mainly defensive for the control family.
                n_mod = n % m
                p_mod = p % m

                # Search existence of q is handled by building a direct
                # residue set only in the special non-unit case.
                valid = False
                break

            inv = pow(p, -1, m)
            wanted.append((n * inv) % m)

        if not valid:
            continue

        if tuple(wanted) in signatures:
            survivors += 1

    fraction = survivors / len(primes)

    return survivors, fraction


# ============================================================================
# EXACT SEARCH
# ============================================================================

def ordinary_prime_scan(
    ordered_primes: list[int],
    n: int,
) -> tuple[int | None, int]:
    """
    Baseline scan of primes in supplied order.
    """
    tests = 0

    for p in ordered_primes:
        if n % p == 0:
            q = n // p

            if sympy.isprime(q):
                return p, tests + 1

        tests += 1

    return None, tests


def crt_filtered_exact_scan(
    ordered_primes: list[int],
    primes: list[int],
    n: int,
    index: SignatureIndex,
) -> tuple[int | None, int]:
    """
    Scan only primes whose residue signature is compatible with n.
    """
    tests = 0
    moduli = index.moduli
    signatures = index.signatures

    for p in ordered_primes:
        wanted = []

        valid = True

        for m in moduli:
            if math.gcd(p, m) != 1:
                valid = False
                break

            inv = pow(p, -1, m)
            wanted.append((n * inv) % m)

        if not valid:
            continue

        if tuple(wanted) not in signatures:
            continue

        tests += 1

        if n % p == 0:
            q = n // p

            if sympy.isprime(q):
                return p, tests

    return None, tests


# ============================================================================
# OBSERVATION
# ============================================================================

def print_family_table(
    name: str,
    moduli: list[int],
) -> None:
    print()
    print(f"{name}")
    line()

    print(
        f"{'stage':>5} "
        f"{'added':>8} "
        f"{'M':>16} "
        f"{'bits':>5} "
        f"{'phi(M)':>16} "
        f"{'phi/M':>12}"
    )

    current = []

    for i, m in enumerate(moduli, 1):
        current.append(m)

        M = lcm_list(current)
        phi = phi_from_factorization(current)

        print(
            f"{i:5d} "
            f"{m:8d} "
            f"{M:16,d} "
            f"{M.bit_length():5d} "
            f"{phi:16,d} "
            f"{phi / M:12.6g}"
        )


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    rng = random.Random(SEED)
    total_start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 48")
    print("CYCLOTOMIC SIGNATURE VS RANDOM-PRIME MODULUS CONTROL")
    print("NO CSV OUTPUT")
    print("=" * 78)

    print(f"random seed      = {SEED}")
    print(f"targets          = {TARGETS}")
    print(f"prime interval   = [{PRIME_LO:,}, {PRIME_HI:,})")
    print(f"R values         = {R_VALUES}")

    # ------------------------------------------------------------------
    # PRIME POPULATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. PRIME POPULATION")
    print("=" * 78)

    t0 = time.perf_counter()

    primes = generate_prime_population()

    population_time = time.perf_counter() - t0

    print(f"prime population = {len(primes):,}")
    print(f"generation time  = {population_time:.4f}s")

    # ------------------------------------------------------------------
    # TARGETS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. TARGETS")
    print("=" * 78)

    targets = generate_targets(primes, TARGETS, rng)

    for i, (p, q, n) in enumerate(targets, 1):
        print(
            f"target {i:2d}: "
            f"p={p} q={q} n={n}"
        )

    # ------------------------------------------------------------------
    # CYCLOTOMIC FAMILY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. CYCLOTOMIC PRIME-FACTOR FAMILY")
    print("=" * 78)

    cyclo_moduli = cyclotomic_prime_inventory(R_VALUES)

    print("distinct non-exceptional cyclotomic prime factors:")
    print(cyclo_moduli)

    print()
    print("F(r) factorization:")

    for r in R_VALUES:
        print(
            f"r={r:2d} "
            f"F(r)={F(r):5d} "
            f"factors={sympy.factorint(F(r))}"
        )

    # ------------------------------------------------------------------
    # CONTROL FAMILY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. CONTROL PRIME FAMILY")
    print("=" * 78)

    control_exclude = set(cyclo_moduli)

    control_moduli = choose_control_primes(
        exclude=control_exclude,
        count=min(CONTROL_MODULI_COUNT, len(cyclo_moduli)),
        rng=rng,
    )

    # To compare stage-by-stage, use the same number of moduli.
    cyclo_moduli = cyclo_moduli[:len(control_moduli)]

    print(f"cyclotomic count = {len(cyclo_moduli)}")
    print(f"control count    = {len(control_moduli)}")
    print()
    print(f"cyclotomic = {cyclo_moduli}")
    print(f"control    = {control_moduli}")

    # ------------------------------------------------------------------
    # MODULUS TABLES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. CUMULATIVE MODULUS STRUCTURE")
    print("=" * 78)

    print_family_table(
        "CYCLOTOMIC FAMILY",
        cyclo_moduli,
    )

    print_family_table(
        "CONTROL FAMILY",
        control_moduli,
    )

    # ------------------------------------------------------------------
    # BUILD PREFIX INDEXES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. SIGNATURE INDEX CONSTRUCTION")
    print("=" * 78)

    cyclo_indexes: dict[int, SignatureIndex] = {}
    control_indexes: dict[int, SignatureIndex] = {}

    for stage in range(1, len(cyclo_moduli) + 1):
        t0 = time.perf_counter()

        idx = build_signature_index(
            primes,
            cyclo_moduli[:stage],
        )

        elapsed = time.perf_counter() - t0

        cyclo_indexes[stage] = idx

        print(
            f"cyclotomic stage={stage:2d} "
            f"M={idx.M:>20,d} "
            f"phi(M)={idx.phi:>20,d} "
            f"signatures={len(idx.signatures):>8,d} "
            f"time={elapsed:.4f}s"
        )

    for stage in range(1, len(control_moduli) + 1):
        t0 = time.perf_counter()

        idx = build_signature_index(
            primes,
            control_moduli[:stage],
        )

        elapsed = time.perf_counter() - t0

        control_indexes[stage] = idx

        print(
            f"control     stage={stage:2d} "
            f"M={idx.M:>20,d} "
            f"phi(M)={idx.phi:>20,d} "
            f"signatures={len(idx.signatures):>8,d} "
            f"time={elapsed:.4f}s"
        )

    # ------------------------------------------------------------------
    # SHARED PRIME ORDER
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. BASELINE PRIME ORDER")
    print("=" * 78)

    sqrt_ordered_primes = []

    for n_data in targets:
        pass

    # Since each target has a different sqrt, create the order per target.
    # This is cheap compared with exact search.
    def closest_order(n: int) -> list[int]:
        root = math.isqrt(n)

        return sorted(
            primes,
            key=lambda p: abs(p - root),
        )

    # ------------------------------------------------------------------
    # MAIN COMPARISON
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. OBSERVED COMPATIBLE-PRIME DENSITY")
    print("=" * 78)

    cyclo_results: dict[tuple[int, int], int] = {}
    control_results: dict[tuple[int, int], int] = {}

    for target_id, (_, _, n) in enumerate(targets, 1):
        print()
        print(f"TARGET {target_id}: n={n}")

        print()
        print(
            f"{'stage':>5} "
            f"{'cyclo':>10} "
            f"{'cyclo frac':>13} "
            f"{'1/phi':>13} "
            f"{'cyclo ratio':>13} "
            f"{'control':>10} "
            f"{'control frac':>13} "
            f"{'control ratio':>13}"
        )
        line()

        for stage in range(1, len(cyclo_moduli) + 1):
            c_idx = cyclo_indexes[stage]
            r_idx = control_indexes[stage]

            cyclo_count, cyclo_fraction = compatible_prime_count(
                primes,
                n,
                c_idx,
            )

            control_count, control_fraction = compatible_prime_count(
                primes,
                n,
                r_idx,
            )

            cyclo_expected = 1.0 / c_idx.phi
            control_expected = 1.0 / r_idx.phi

            cyclo_ratio = (
                cyclo_fraction / cyclo_expected
                if cyclo_expected > 0
                else float("inf")
            )

            control_ratio = (
                control_fraction / control_expected
                if control_expected > 0
                else float("inf")
            )

            cyclo_results[(target_id, stage)] = cyclo_count
            control_results[(target_id, stage)] = control_count

            print(
                f"{stage:5d} "
                f"{cyclo_count:10,d} "
                f"{cyclo_fraction:13.6g} "
                f"{cyclo_expected:13.6g} "
                f"{cyclo_ratio:13.6g} "
                f"{control_count:10,d} "
                f"{control_fraction:13.6g} "
                f"{control_ratio:13.6g}"
            )

    # ------------------------------------------------------------------
    # SUMMARY OVER TARGETS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. MEDIAN DENSITY SUMMARY")
    print("=" * 78)

    print(
        f"{'stage':>5} "
        f"{'cyclo median':>15} "
        f"{'cyclo frac':>13} "
        f"{'cyclo / 1phi':>14} "
        f"{'control median':>15} "
        f"{'control frac':>13} "
        f"{'control / 1phi':>15}"
    )
    line()

    for stage in range(1, len(cyclo_moduli) + 1):
        c_values = [
            cyclo_results[(target_id, stage)]
            for target_id in range(1, TARGETS + 1)
        ]

        r_values = [
            control_results[(target_id, stage)]
            for target_id in range(1, TARGETS + 1)
        ]

        c_median = statistics.median(c_values)
        r_median = statistics.median(r_values)

        c_fraction = c_median / len(primes)
        r_fraction = r_median / len(primes)

        c_expected = 1.0 / cyclo_indexes[stage].phi
        r_expected = 1.0 / control_indexes[stage].phi

        c_ratio = c_fraction / c_expected
        r_ratio = r_fraction / r_expected

        print(
            f"{stage:5d} "
            f"{c_median:15,.1f} "
            f"{c_fraction:13.6g} "
            f"{c_ratio:14.6g} "
            f"{r_median:15,.1f} "
            f"{r_fraction:13.6g} "
            f"{r_ratio:15.6g}"
        )

    # ------------------------------------------------------------------
    # EXACT RECOVERY BENCHMARK
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. EXACT RECOVERY BENCHMARK")
    print("=" * 78)

    print(
        "Selected stages are tested against ordinary closest-prime scanning."
    )

    baseline_times = []
    cyclo_times = []
    control_times = []

    for target_id, (_, _, n) in enumerate(targets, 1):
        ordered = closest_order(n)

        # Baseline
        t0 = time.perf_counter()

        baseline_factor, baseline_tests = ordinary_prime_scan(
            ordered,
            n,
        )

        baseline_elapsed = time.perf_counter() - t0

        baseline_times.append(baseline_elapsed)

        if target_id <= 3:
            print()
            print(
                f"target {target_id:2d}: "
                f"baseline tests={baseline_tests:,} "
                f"time={baseline_elapsed:.6f}s "
                f"correct={baseline_factor is not None}"
            )

        for stage in sorted(EXACT_BENCHMARK_STAGES):
            if stage > len(cyclo_moduli):
                continue

            # Cyclotomic
            t0 = time.perf_counter()

            c_factor, c_tests = crt_filtered_exact_scan(
                ordered,
                primes,
                n,
                cyclo_indexes[stage],
            )

            c_elapsed = time.perf_counter() - t0
            cyclo_times.append(c_elapsed)

            # Control
            t0 = time.perf_counter()

            r_factor, r_tests = crt_filtered_exact_scan(
                ordered,
                primes,
                n,
                control_indexes[stage],
            )

            r_elapsed = time.perf_counter() - t0
            control_times.append(r_elapsed)

            if target_id <= 3:
                print(
                    f"  stage {stage:2d} "
                    f"cyclo tests={c_tests:8,d} "
                    f"time={c_elapsed:.6f}s "
                    f"correct={c_factor is not None}"
                )
                print(
                    f"           "
                    f"control tests={r_tests:8,d} "
                    f"time={r_elapsed:.6f}s "
                    f"correct={r_factor is not None}"
                )

    # ------------------------------------------------------------------
    # FINAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("11. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        """
This experiment tests the central hypothesis:

    Is the observed KAPPA filtering special to
    F(r)=r^2+r+1,

or is it simply the normal behaviour of accumulating
more modular constraints?

Interpretation:

    cyclotomic ratio near control ratio
        ->
        no evidence that the Phi_3 structure itself
        provides extra filtering.

    cyclotomic ratio materially different from control
        ->
        investigate the algebraic cyclotomic structure further.

Important distinction:

    1/phi(M)

is the natural unit-group reference, while

    1/M

is generally the wrong baseline when searchable primes
are all units modulo M.

The experiment therefore reports both the measured density
and its ratio to 1/phi(M).

The exact-search benchmark then asks whether any density
advantage translates into actual factor-search speed.

No prime-pair enumeration is performed.
No CSV files are produced.
"""
    )

    total_elapsed = time.perf_counter() - total_start

    print("=" * 78)
    print("EXPERIMENT 48 COMPLETE")
    print("=" * 78)
    print(f"total runtime = {total_elapsed:.3f}s")


if __name__ == "__main__":
    main()

