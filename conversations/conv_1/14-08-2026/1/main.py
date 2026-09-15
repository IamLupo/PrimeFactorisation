#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 44
DIRECT RESIDUE-CLASS RETRIEVAL / CRT SEARCH
NO CSV OUTPUT
==============================================================================

Goal
----
Experiment 43 showed:

    ~177x fewer divisibility tests at prefix 7

but was still slower because it scanned the entire prime population and
performed a Python-level residue filter on every prime.

Experiment 44 removes that scan.

Instead:

    prime -> residue class
             |
             v
    residue class -> primes

For each target n:

    1. compute the generalized CRT modulus M
    2. determine the compatible residue classes for p
    3. directly retrieve the primes belonging to those classes
    4. test n % p only on the retrieved primes

This measures whether the large mathematical filtering advantage can
be converted into an actual wall-clock advantage.

No candidate-pair list is built.
No CSV files are produced.
Only console output is generated.
"""

from __future__ import annotations

import math
import random
import time
from collections import defaultdict
from statistics import median
from typing import Dict, Iterable, List, Sequence, Tuple

from sympy import isprime


# =============================================================================
# CONFIGURATION
# =============================================================================

SEED = 20260814

TARGETS = 100

PRIME_LO = 2_000_000
PRIME_HI = 4_200_000

PREFIXES = [6, 7]

R_VALUES = [
    2, 3, 5, 7, 11, 13, 17, 19,
    23, 29, 31, 37, 41, 43, 47,
]

PRINT_TARGETS = True

# Maximum number of retrieved candidate primes to print for one target.
PRINT_MAX_CANDIDATES = 12


# =============================================================================
# MODULUS CONSTRUCTION
# =============================================================================

def build_moduli(r_values: Sequence[int]) -> List[int]:
    """
    Reconstruct the modulus family used throughout the experiments.

    m(r) = r^2 + r + 1
    """
    return [r * r + r + 1 for r in r_values]


def lcm_many(values: Sequence[int]) -> int:
    result = 1
    for value in values:
        result = math.lcm(result, value)
    return result


# =============================================================================
# PRIME GENERATION
# =============================================================================

def generate_primes(lo: int, hi: int) -> List[int]:
    primes: List[int] = []

    # Handle the even prime explicitly.
    if lo <= 2 < hi:
        primes.append(2)

    start = max(lo, 3)
    if start % 2 == 0:
        start += 1

    for x in range(start, hi, 2):
        if isprime(x):
            primes.append(x)

    return primes


# =============================================================================
# TARGET GENERATION
# =============================================================================

def generate_targets(
    primes: Sequence[int],
    count: int,
    seed: int,
) -> List[Tuple[int, int, int]]:
    """
    Choose distinct semiprime targets from the supplied prime population.

    The factors are selected from the central portion of the interval to
    avoid trivial boundary effects.
    """
    rng = random.Random(seed)

    if len(primes) < 2:
        raise RuntimeError("Prime population is too small.")

    targets: List[Tuple[int, int, int]] = []
    seen_pairs = set()

    # Prefer factors away from the interval boundaries.
    low_index = max(0, len(primes) // 20)
    high_index = min(len(primes) - 1, len(primes) - len(primes) // 20 - 1)

    available = primes[low_index:high_index + 1]

    while len(targets) < count:
        p = rng.choice(available)
        q = rng.choice(available)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        key = (p, q)
        if key in seen_pairs:
            continue

        seen_pairs.add(key)
        targets.append((p, q, p * q))

    return targets


# =============================================================================
# RESIDUE / CRT HELPERS
# =============================================================================

def residue_vector(x: int, moduli: Sequence[int]) -> Tuple[int, ...]:
    return tuple(x % m for m in moduli)


def merge_congruence(
    a1: int,
    m1: int,
    a2: int,
    m2: int,
) -> Tuple[int, int]:
    """
    Solve:

        x = a1 (mod m1)
        x = a2 (mod m2)

    for non-coprime moduli.

    Returns:
        (x0, lcm(m1, m2))

    Raises ValueError if inconsistent.
    """
    g = math.gcd(m1, m2)

    if (a2 - a1) % g != 0:
        raise ValueError("Inconsistent generalized CRT system.")

    m1g = m1 // g
    m2g = m2 // g

    # m1g and m2g are coprime.
    t = ((a2 - a1) // g) * pow(m1g, -1, m2g)
    t %= m2g

    x = a1 + m1 * t
    mod = m1 * m2g

    return x % mod, mod


def generalized_crt(
    residues: Sequence[int],
    moduli: Sequence[int],
) -> Tuple[int, int]:
    if len(residues) != len(moduli):
        raise ValueError("Residue/modulus length mismatch.")

    x = residues[0] % moduli[0]
    mod = moduli[0]

    for a, m in zip(residues[1:], moduli[1:]):
        x, mod = merge_congruence(x, mod, a % m, m)

    return x, mod


def compatible_p_classes(
    n: int,
    moduli: Sequence[int],
) -> Tuple[int, Dict[int, int]]:
    """
    Find all residue classes x modulo M for which there exists y satisfying:

        x*y = n (mod each selected modulus)

    For unit x, y is determined.

    Because the moduli are small for this experiment, we construct the
    compatible classes by enumerating the unit residues modulo each modulus
    and merging the resulting congruence systems.

    Returns:
        M
        dict {canonical residue x mod M: canonical q residue mod M}
    """
    M = lcm_many(moduli)

    # We work modulo M.
    #
    # Any prime in our population is coprime to every modulus in this
    # particular family, so x is a unit modulo each m.
    #
    # For each possible residue x modulo M represented by its residues
    # modulo the individual moduli, q is n*x^{-1} modulo every modulus.

    #
    # Generate compatible residue vectors progressively.
    #
    classes: Dict[int, int] = {0: 0}

    current_mod = 1

    for m in moduli:
        new_classes: Dict[int, int] = {}

        # Unit residues modulo m.
        unit_residues = [
            a for a in range(m)
            if math.gcd(a, m) == 1
        ]

        # For each prior partial x, combine with a new residue a.
        #
        # We don't actually need to combine against every q independently:
        # q residue follows from the inverse relation.
        #
        # Store pairs (x,q) in the dictionary using a packed tuple first.
        partial_items = list(classes.items())

        for x_old, q_old in partial_items:
            for a in unit_residues:
                # q residue modulo m
                q_a = (n % m) * pow(a, -1, m) % m

                # Merge x residues.
                try:
                    x_new, merged_mod = merge_congruence(
                        x_old,
                        current_mod,
                        a,
                        m,
                    )
                except ValueError:
                    continue

                # Merge q residues.
                try:
                    q_new, _ = merge_congruence(
                        q_old,
                        current_mod,
                        q_a,
                        m,
                    )
                except ValueError:
                    continue

                new_classes[x_new] = q_new

        classes = new_classes
        current_mod = math.lcm(current_mod, m)

    return M, classes


# =============================================================================
# DIRECT INDEX
# =============================================================================

class DirectResidueIndex:
    """
    Direct mapping:

        residue modulo M -> tuple/list of primes

    No scanning is needed after the index exists.
    """

    def __init__(self, primes: Sequence[int], M: int):
        self.M = M
        self.groups: Dict[int, Tuple[int, ...]] = {}

        buckets: Dict[int, List[int]] = defaultdict(list)

        for p in primes:
            if math.gcd(p, M) != 1:
                continue

            buckets[p % M].append(p)

        # Convert to tuples to make accidental mutation impossible.
        self.groups = {
            residue: tuple(values)
            for residue, values in buckets.items()
        }

    def get(self, residue: int) -> Tuple[int, ...]:
        return self.groups.get(residue, ())

    def occupied(self) -> int:
        return len(self.groups)

    def max_class_size(self) -> int:
        if not self.groups:
            return 0
        return max(len(v) for v in self.groups.values())


# =============================================================================
# BASELINE
# =============================================================================

def closest_order(
    primes: Sequence[int],
    n: int,
) -> List[int]:
    root = math.isqrt(n)

    # Sorting 147k primes repeatedly would add unnecessary overhead.
    #
    # For the benchmark, build the order once per target.
    return sorted(
        primes,
        key=lambda p: (abs(p - root), p),
    )


def ordinary_scan(
    primes: Sequence[int],
    n: int,
) -> Tuple[int | None, int, float]:
    ordered = closest_order(primes, n)

    tests = 0
    start = time.perf_counter()

    for p in ordered:
        tests += 1

        if n % p == 0:
            q = n // p
            if q >= p and isprime(q):
                elapsed = time.perf_counter() - start
                return p, tests, elapsed

    elapsed = time.perf_counter() - start
    return None, tests, elapsed


# =============================================================================
# DIRECT CRT SEARCH
# =============================================================================

def build_target_q_class_map(
    n: int,
    moduli: Sequence[int],
) -> Tuple[int, Dict[int, int]]:
    """
    Wrapper around generalized CRT class construction.

    Returns:
        M
        {p_residue: q_residue}
    """
    return compatible_p_classes(n, moduli)


def direct_crt_search(
    n: int,
    index: DirectResidueIndex,
    p_to_q: Dict[int, int],
    ordered_primes: Sequence[int],
) -> Tuple[int | None, int, int, int, float]:
    """
    Retrieve only compatible residue classes directly.

    We do not scan the entire prime population.

    We also preserve closest-to-sqrt ordering by assigning a rank to the
    retrieved primes instead of sorting the entire population.
    """
    start = time.perf_counter()

    # Direct retrieval phase.
    candidate_set = set()

    for p_residue, q_residue in p_to_q.items():
        # We only need to retrieve the p-side classes.
        #
        # q_residue is retained for diagnostics and verification.
        for p in index.get(p_residue):
            candidate_set.add(p)

    retrieval_count = len(candidate_set)

    # Rank only the retrieved candidates according to the shared
    # closest-to-sqrt ordering.
    rank = {p: i for i, p in enumerate(ordered_primes)}

    candidates = sorted(
        candidate_set,
        key=rank.__getitem__,
    )

    tests = 0

    for p in candidates:
        tests += 1

        if n % p != 0:
            continue

        q = n // p

        if q >= 2 and isprime(q):
            elapsed = time.perf_counter() - start
            return p, retrieval_count, tests, len(candidates), elapsed

    elapsed = time.perf_counter() - start
    return None, retrieval_count, tests, len(candidates), elapsed


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    total_start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 44")
    print("DIRECT RESIDUE-CLASS RETRIEVAL / CRT SEARCH")
    print("NO CSV OUTPUT")
    print("=" * 78)
    print()

    print(f"random seed      = {SEED}")
    print(f"targets          = {TARGETS}")
    print(f"prime interval   = [{PRIME_LO:,}, {PRIME_HI:,})")
    print(f"prefixes         = {PREFIXES}")
    print(f"R values         = {R_VALUES}")
    print()

    # -------------------------------------------------------------------------
    # MODULI
    # -------------------------------------------------------------------------

    all_moduli = build_moduli(R_VALUES)

    print("-" * 78)
    print("1. MODULUS INVENTORY")
    print("-" * 78)

    for r, m in zip(R_VALUES, all_moduli):
        print(f"r={r:>3} m={m:>10,}")

    print()

    prefix_moduli: Dict[int, List[int]] = {}
    prefix_M: Dict[int, int] = {}

    print("-" * 78)
    print("2. GENERALIZED CRT PREFIX MODULI")
    print("-" * 78)

    for prefix in PREFIXES:
        mods = all_moduli[:prefix]
        M = lcm_many(mods)

        prefix_moduli[prefix] = mods
        prefix_M[prefix] = M

        print(
            f"prefix={prefix:>2} "
            f"M={M:>14,} "
            f"bits={M.bit_length():>2} "
            f"moduli={mods}"
        )

    print()

    # -------------------------------------------------------------------------
    # PRIME POPULATION
    # -------------------------------------------------------------------------

    print("-" * 78)
    print("3. PRIME POPULATION")
    print("-" * 78)

    start = time.perf_counter()
    primes = generate_primes(PRIME_LO, PRIME_HI)
    prime_time = time.perf_counter() - start

    print(f"prime population = {len(primes):,}")
    print(f"generation time  = {prime_time:.6f}s")
    print()

    # -------------------------------------------------------------------------
    # TARGETS
    # -------------------------------------------------------------------------

    print("-" * 78)
    print("4. TARGETS")
    print("-" * 78)

    targets = generate_targets(primes, TARGETS, SEED)

    for i, (p, q, n) in enumerate(targets, 1):
        print(
            f"target {i:>3}: "
            f"p={p} q={q} n={n}"
        )

    print()

    # -------------------------------------------------------------------------
    # SHARED CLOSEST ORDER
    # -------------------------------------------------------------------------

    print("-" * 78)
    print("5. SHARED CLOSEST-PRIME ORDER")
    print("-" * 78)
    print(
        "The same closest-to-sqrt order is used by the baseline and "
        "direct-retrieval searches."
    )
    print()

    ordered_targets: Dict[int, List[int]] = {}

    for _, _, n in targets:
        ordered_targets[n] = closest_order(primes, n)

    # -------------------------------------------------------------------------
    # DIRECT INDEX CONSTRUCTION
    # -------------------------------------------------------------------------

    print("-" * 78)
    print("6. DIRECT RESIDUE-CLASS INDEX")
    print("-" * 78)

    indexes: Dict[int, DirectResidueIndex] = {}

    index_total = 0.0

    for prefix in PREFIXES:
        M = prefix_M[prefix]

        start = time.perf_counter()
        index = DirectResidueIndex(primes, M)
        elapsed = time.perf_counter() - start

        indexes[prefix] = index
        index_total += elapsed

        print(
            f"prefix={prefix:>2} "
            f"M={M:>14,} "
            f"occupied={index.occupied():>8,} "
            f"max_class={index.max_class_size():>4} "
            f"time={elapsed:.6f}s"
        )

    print(f"total index construction = {index_total:.6f}s")
    print()

    # -------------------------------------------------------------------------
    # BASELINE
    # -------------------------------------------------------------------------

    print("-" * 78)
    print("7. ORDINARY CLOSEST-PRIME BASELINE")
    print("-" * 78)

    baseline_times: List[float] = []
    baseline_tests: List[int] = []
    baseline_correct = 0

    for i, (true_p, true_q, n) in enumerate(targets, 1):
        ordered = ordered_targets[n]

        start = time.perf_counter()
        tests = 0
        found = None

        for p in ordered:
            tests += 1

            if n % p == 0:
                q = n // p

                if q >= p and isprime(q):
                    found = (p, q)
                    break

        elapsed = time.perf_counter() - start

        ok = (
            found is not None
            and set(found) == {true_p, true_q}
        )

        if ok:
            baseline_correct += 1

        baseline_times.append(elapsed)
        baseline_tests.append(tests)

        print(
            f"target {i:>3}: "
            f"tests={tests:>7,} "
            f"time={elapsed:.6f}s "
            f"correct={ok}"
        )

    print()

    # -------------------------------------------------------------------------
    # DIRECT RETRIEVAL
    # -------------------------------------------------------------------------

    results = {}

    for prefix in PREFIXES:
        print("-" * 78)
        print(f"8. DIRECT RETRIEVAL PREFIX {prefix}")
        print("-" * 78)

        M = prefix_M[prefix]
        index = indexes[prefix]

        times: List[float] = []
        tests_list: List[int] = []
        retrieved_list: List[int] = []
        correct_count = 0

        for i, (true_p, true_q, n) in enumerate(targets, 1):
            ordered = ordered_targets[n]

            start = time.perf_counter()

            M2, p_to_q = build_target_q_class_map(
                n,
                prefix_moduli[prefix],
            )

            if M2 != M:
                raise RuntimeError(
                    f"CRT modulus mismatch: {M2} != {M}"
                )

            found, retrieved, tests, candidate_count, elapsed_search = (
                direct_crt_search(
                    n,
                    index,
                    p_to_q,
                    ordered,
                )
            )

            total_elapsed = time.perf_counter() - start

            ok = (
                found is not None
                and set((found, n // found)) == {true_p, true_q}
            )

            if ok:
                correct_count += 1

            times.append(total_elapsed)
            tests_list.append(tests)
            retrieved_list.append(retrieved)

            print(
                f"target {i:>3}: "
                f"classes={len(p_to_q):>7,} "
                f"retrieved={retrieved:>6,} "
                f"tests={tests:>7,} "
                f"time={total_elapsed:.6f}s "
                f"correct={ok}"
            )

            if PRINT_TARGETS and found is not None:
                print(
                    f"             recovered=({found}, {n // found})"
                )

        results[prefix] = {
            "times": times,
            "tests": tests_list,
            "retrieved": retrieved_list,
            "correct": correct_count,
        }

        print()
        print(
            f"prefix {prefix}: "
            f"median_tests={median(tests_list):,.1f} "
            f"median_retrieved={median(retrieved_list):,.1f} "
            f"median_time={median(times):.6f}s "
            f"batch={sum(times):.6f}s "
            f"recovered={correct_count}/{TARGETS}"
        )
        print()

        # Print a few examples for the first target only.
        if PRINT_TARGETS and targets:
            true_p, true_q, n = targets[0]
            M2, p_to_q = build_target_q_class_map(
                n,
                prefix_moduli[prefix],
            )

            candidate_preview: List[int] = []

            for residue in p_to_q:
                for p in index.get(residue):
                    candidate_preview.append(p)
                    if len(candidate_preview) >= PRINT_MAX_CANDIDATES:
                        break

                if len(candidate_preview) >= PRINT_MAX_CANDIDATES:
                    break

            candidate_preview.sort()

            print(
                f"target 1 direct-retrieval preview "
                f"(first {len(candidate_preview)} candidates):"
            )
            print(f"  {candidate_preview}")
            print()

    # -------------------------------------------------------------------------
    # FINAL COMPARISON
    # -------------------------------------------------------------------------

    print("=" * 78)
    print("9. FINAL COMPARISON")
    print("=" * 78)

    print()
    print(
        f"{'method':32s}"
        f"{'median tests':>16s}"
        f"{'median time':>16s}"
        f"{'batch time':>16s}"
    )
    print("-" * 78)

    baseline_batch = sum(baseline_times)

    print(
        f"{'ordinary closest-prime scan':32s}"
        f"{median(baseline_tests):>16,.1f}"
        f"{median(baseline_times):>16.6f}"
        f"{baseline_batch:>16.6f}"
    )

    for prefix in PREFIXES:
        data = results[prefix]

        print(
            f"{('direct CRT prefix ' + str(prefix)):32s}"
            f"{median(data['tests']):>16,.1f}"
            f"{median(data['times']):>16.6f}"
            f"{sum(data['times']):>16.6f}"
        )

    print()

    # -------------------------------------------------------------------------
    # SPEEDUP / TEST REDUCTION
    # -------------------------------------------------------------------------

    print("=" * 78)
    print("10. TEST REDUCTION AND RETRIEVAL EFFECT")
    print("=" * 78)

    base_med_tests = median(baseline_tests)

    for prefix in PREFIXES:
        data = results[prefix]

        med_tests = median(data["tests"])
        med_retrieved = median(data["retrieved"])
        med_time = median(data["times"])

        test_reduction = (
            base_med_tests / med_tests
            if med_tests > 0 else float("inf")
        )

        time_speed = (
            median(baseline_times) / med_time
            if med_time > 0 else float("inf")
        )

        print(
            f"prefix {prefix}: "
            f"tests={base_med_tests:,.1f} -> {med_tests:,.1f} "
            f"({test_reduction:.3f}x reduction), "
            f"retrieved={med_retrieved:,.1f}, "
            f"time speed={time_speed:.3f}x"
        )

    print()

    # -------------------------------------------------------------------------
    # INDEX AMORTIZATION
    # -------------------------------------------------------------------------

    print("=" * 78)
    print("11. INDEX AMORTIZATION")
    print("=" * 78)

    base_median_time = median(baseline_times)

    for prefix in PREFIXES:
        direct_median_time = median(results[prefix]["times"])

        if direct_median_time < base_median_time:
            per_target_saving = base_median_time - direct_median_time

            break_even = (
                index_total / per_target_saving
                if per_target_saving > 0 else float("inf")
            )
        else:
            per_target_saving = 0.0
            break_even = math.inf

        print(
            f"prefix {prefix}: "
            f"baseline median={base_median_time:.6f}s "
            f"direct median={direct_median_time:.6f}s "
            f"break-even targets="
            f"{'NEVER' if math.isinf(break_even) else f'{break_even:.1f}'}"
        )

    print()

    # -------------------------------------------------------------------------
    # CORRECTNESS
    # -------------------------------------------------------------------------

    print("=" * 78)
    print("12. CORRECTNESS")
    print("=" * 78)

    print(
        f"ordinary baseline: "
        f"{baseline_correct}/{TARGETS}"
    )

    for prefix in PREFIXES:
        print(
            f"direct CRT prefix {prefix}: "
            f"{results[prefix]['correct']}/{TARGETS}"
        )

    print()

    # -------------------------------------------------------------------------
    # INTERPRETATION
    # -------------------------------------------------------------------------

    print("=" * 78)
    print("13. INTERPRETATION")
    print("=" * 78)

    print()
    print(
        "Experiment 43 still scanned every prime before the residue filter.\n"
        "Experiment 44 removes that scan from the CRT candidate-retrieval\n"
        "phase by directly indexing residue class -> prime."
    )

    print()
    print(
        "The key comparison is now:"
    )

    print()
    print(
        "    ordinary prime scan\n"
        "        versus\n"
        "    direct compatible-class retrieval"
    )

    print()
    print(
        "If prefix 7 now approaches or beats the ordinary scan in wall-clock\n"
        "time, that would show that the previous slowdown was primarily due\n"
        "to scanning/filtering overhead."
    )

    print()
    print(
        "If prefix 7 is still slower, then the remaining cost is concentrated\n"
        "in constructing the target-specific compatible residue classes and\n"
        "retrieving/ranking the corresponding primes."
    )

    print()
    print(
        "The experiment deliberately reports all three quantities:\n"
        "    retrieved primes\n"
        "    actual n % p tests\n"
        "    elapsed time"
    )

    print()
    print("No candidate-pair list is constructed.")
    print("No CSV files are produced.")

    total_elapsed = time.perf_counter() - total_start

    print()
    print("=" * 78)
    print("EXPERIMENT 44 COMPLETE")
    print("=" * 78)
    print(f"total runtime = {total_elapsed:.6f}s")


if __name__ == "__main__":
    main()

