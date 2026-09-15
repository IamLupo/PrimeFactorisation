#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 47
CYCLOTOMIC PRIME-FACTOR FILTER CONTRIBUTION
NO CSV OUTPUT
==============================================================================

Goal
----
Experiment 46 established:

    F(r) = r^2 + r + 1 = Phi_3(r)
    F(r) | r^3 - 1

and showed that, apart from ell=3, every prime ell | F(r)
has ord_ell(r) = 3.

Experiment 47 breaks the modulus family into its DISTINCT PRIME
FACTORS and asks:

    1. Which prime factors are genuinely independent?
    2. What residue restriction does each prime ell impose?
    3. How much candidate-prime reduction does each ell provide?
    4. How much additional reduction appears when independent ell's
       are accumulated?
    5. Is the observed ~4x multiplier over 1/M related to the
       cyclotomic prime structure?

This experiment does NOT enumerate prime pairs.

It works with:
    * the fixed prime population
    * n-only residue constraints
    * distinct prime factors of F(r)
    * generalized CRT / residue intersections

All exact factor recovery is a control only.

NO CSV OUTPUT.
Only console output.
"""

from __future__ import annotations

import math
import random
import time
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Set, Tuple

import sympy


# ============================================================================
# CONFIGURATION
# ============================================================================

SEED = 20260814

PRIME_LOW = 2_000_000
PRIME_HIGH = 4_200_000

PRIME_POOL_SIZE = 5000

TARGETS = 12

R_VALUES = [
    2, 3, 5, 7, 11, 13, 17,
    19, 23, 29, 31, 37, 41, 43, 47
]

PREFIXES = [3, 4, 5, 6, 7]

# We only need one exact prime population for this experiment.
# This keeps runtime and memory modest.
assert PRIME_LOW < PRIME_HIGH


# ============================================================================
# DATA TYPES
# ============================================================================

@dataclass(frozen=True)
class Target:
    p: int
    q: int
    n: int


@dataclass(frozen=True)
class CyclotomicFactor:
    ell: int
    rs: Tuple[int, ...]
    order3_count: int
    exceptional: bool


# ============================================================================
# BASIC MATH
# ============================================================================

def F(r: int) -> int:
    return r * r + r + 1


def prime_population(low: int, high: int) -> List[int]:
    """
    Generate every prime in [low, high).

    sympy.isprime is intentionally used here, matching the experiment
    convention established earlier.
    """
    primes: List[int] = []
    for x in range(low | 1, high, 2):
        if sympy.isprime(x):
            primes.append(x)
    return primes


def choose_targets(
    rng: random.Random,
    primes: Sequence[int],
    count: int,
) -> List[Target]:
    """
    Generate distinct semiprime targets from the population.
    """
    targets: List[Target] = []
    seen_pairs: Set[Tuple[int, int]] = set()

    while len(targets) < count:
        p = primes[rng.randrange(len(primes))]
        q = primes[rng.randrange(len(primes))]

        if p == q:
            continue

        a, b = sorted((p, q))
        pair = (a, b)

        if pair in seen_pairs:
            continue

        seen_pairs.add(pair)
        targets.append(Target(a, b, a * b))

    return targets


def multiplicative_order(a: int, modulus: int) -> int:
    """
    Small-modulus multiplicative order using sympy.
    """
    return int(sympy.n_order(a % modulus, modulus))


def lcm_many(values: Iterable[int]) -> int:
    result = 1
    for value in values:
        result = math.lcm(result, value)
    return result


# ============================================================================
# GENERALIZED CRT
# ============================================================================

def generalized_pair_residues(
    n: int,
    modulus: int,
) -> Set[Tuple[int, int]]:
    """
    Enumerate unit residue pairs (a,b) modulo modulus such that:

        a*b == n (mod modulus)

    This routine is only intended for the small prime factors ell of
    F(r), not for the giant LCMs used in earlier experiments.

    For each a coprime to modulus:

        b = n * a^{-1} mod modulus.
    """
    result: Set[Tuple[int, int]] = set()

    n_mod = n % modulus

    for a in range(modulus):
        if math.gcd(a, modulus) != 1:
            continue

        b = (n_mod * pow(a, -1, modulus)) % modulus

        if (a * b) % modulus == n_mod:
            result.add((a, b))

    return result


def allowed_prime_residues(
    n: int,
    ell: int,
) -> Set[int]:
    """
    Return residues a mod ell for which there exists a unit b mod ell
    satisfying a*b = n mod ell.

    Since all searchable primes in our interval are > ell and therefore
    nonzero mod ell unless ell itself lies in the interval (it does not),
    we effectively work with unit residues.
    """
    n_mod = n % ell
    allowed: Set[int] = set()

    for a in range(1, ell):
        if math.gcd(a, ell) != 1:
            continue

        b = (n_mod * pow(a, -1, ell)) % ell

        if (a * b) % ell == n_mod:
            allowed.add(a)

    return allowed


# ============================================================================
# CYCLOTOMIC STRUCTURE
# ============================================================================

def factor_family() -> Dict[int, List[int]]:
    """
    Return prime-factor -> list of r values for which ell | F(r).
    """
    family: Dict[int, List[int]] = defaultdict(list)

    for r in R_VALUES:
        factors = sympy.factorint(F(r))
        for ell in factors:
            family[int(ell)].append(r)

    return dict(family)


def build_factor_inventory() -> Dict[int, CyclotomicFactor]:
    family = factor_family()
    inventory: Dict[int, CyclotomicFactor] = {}

    for ell, rs in sorted(family.items()):
        order3 = 0
        for r in rs:
            if ell == 3:
                continue

            if math.gcd(r, ell) == 1 and multiplicative_order(r, ell) == 3:
                order3 += 1

        inventory[ell] = CyclotomicFactor(
            ell=ell,
            rs=tuple(rs),
            order3_count=order3,
            exceptional=(ell == 3),
        )

    return inventory


# ============================================================================
# PRIME RESIDUE INDEX
# ============================================================================

def build_residue_index(
    primes: Sequence[int],
    ell: int,
) -> Dict[int, int]:
    counts: Dict[int, int] = defaultdict(int)

    for p in primes:
        counts[p % ell] += 1

    return dict(counts)


# ============================================================================
# TARGET-LEVEL FILTER MEASUREMENTS
# ============================================================================

def count_factor_candidates(
    primes: Sequence[int],
    n: int,
    ell: int,
    residue_index: Dict[int, int],
) -> int:
    """
    Count searchable prime p for which there exists a prime residue q
    compatible with:

        p*q = n mod ell.

    This is a prime-level filtering count, not pair enumeration.
    """
    n_mod = n % ell
    count = 0

    for p in primes:
        a = p % ell

        if a == 0:
            continue

        if math.gcd(a, ell) != 1:
            continue

        b = (n_mod * pow(a, -1, ell)) % ell

        if residue_index.get(b, 0) > 0:
            count += 1

    return count


def exact_recovery_scan(
    primes: Sequence[int],
    n: int,
    target_p: int,
    target_q: int,
) -> Tuple[bool, int]:
    """
    Exact control: ordinary closest-to-sqrt prime scan.
    """
    root = math.isqrt(n)

    ordered = sorted(
        primes,
        key=lambda p: (abs(p - root), p),
    )

    tests = 0

    for p in ordered:
        if p * p > n:
            continue

        tests += 1

        if n % p == 0:
            q = n // p
            ok = (
                sympy.isprime(q)
                and {p, q} == {target_p, target_q}
            )
            return ok, tests

    return False, tests


def factor_residue_hit(
    p: int,
    n: int,
    ell: int,
) -> bool:
    """
    Check whether p has an ell-residue compatible with some unit q.
    """
    a = p % ell

    if a == 0:
        return False

    if math.gcd(a, ell) != 1:
        return False

    b = (n % ell) * pow(a, -1, ell) % ell

    return True if b is not None else False


# ============================================================================
# INDEPENDENCE ANALYSIS
# ============================================================================

def incremental_lcm_factors(
    inventory: Dict[int, CyclotomicFactor],
) -> List[Tuple[int, int, int, int]]:
    """
    For each distinct cyclotomic prime ell, report:

        ell
        previous LCM
        new LCM
        incremental multiplier
    """
    rows = []
    current = 1

    for ell in sorted(inventory):
        previous = current
        current = math.lcm(current, ell)
        multiplier = current // previous

        rows.append((ell, previous, current, multiplier))

    return rows


def pairwise_factor_overlap(
    inventory: Dict[int, CyclotomicFactor],
) -> List[Tuple[int, int]]:
    """
    Counts pairs of ell values that occur in the same F(r).
    """
    rows: List[Tuple[int, int]] = []

    ell_values = sorted(inventory)

    for i in range(len(ell_values)):
        a = ell_values[i]
        a_rs = set(inventory[a].rs)

        for j in range(i + 1, len(ell_values)):
            b = ell_values[j]
            b_rs = set(inventory[b].rs)

            if a_rs & b_rs:
                rows.append((a, b))

    return rows


# ============================================================================
# EXPERIMENT
# ============================================================================

def main() -> None:
    random.seed(SEED)
    rng = random.Random(SEED)

    total_start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 47")
    print("CYCLOTOMIC PRIME-FACTOR FILTER CONTRIBUTION")
    print("NO CSV OUTPUT")
    print("=" * 78)

    print()
    print("=" * 78)
    print("1. CYCLOTOMIC MODULUS FAMILY")
    print("=" * 78)

    for r in R_VALUES:
        value = F(r)
        factors = sympy.factorint(value)

        print(
            f"r={r:>3} "
            f"F(r)={value:>6} "
            f"factors={factors}"
        )

    # ----------------------------------------------------------------------
    # Prime population
    # ----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. PRIME POPULATION")
    print("=" * 78)

    start = time.perf_counter()

    primes = prime_population(PRIME_LOW, PRIME_HIGH)

    elapsed = time.perf_counter() - start

    print(f"prime population = {len(primes):,}")
    print(f"generation time  = {elapsed:.4f}s")

    # ----------------------------------------------------------------------
    # Targets
    # ----------------------------------------------------------------------

    targets = choose_targets(rng, primes, TARGETS)

    print()
    print("=" * 78)
    print("3. TARGETS")
    print("=" * 78)

    for i, target in enumerate(targets, 1):
        print(
            f"target {i:>2}: "
            f"p={target.p} "
            f"q={target.q} "
            f"n={target.n}"
        )

    # ----------------------------------------------------------------------
    # Cyclotomic factor inventory
    # ----------------------------------------------------------------------

    inventory = build_factor_inventory()

    print()
    print("=" * 78)
    print("4. DISTINCT CYCLOTOMIC PRIME FACTOR INVENTORY")
    print("=" * 78)

    print(
        f"{'ell':>7} "
        f"{'ell%3':>7} "
        f"{'r values':<22} "
        f"{'order-3':>8} "
        f"{'exception':>10}"
    )
    print("-" * 78)

    for ell, info in sorted(inventory.items()):
        print(
            f"{ell:>7} "
            f"{ell % 3:>7} "
            f"{str(list(info.rs)):<22} "
            f"{info.order3_count:>8} "
            f"{str(info.exceptional):>10}"
        )

    print()
    print(f"distinct ell = {len(inventory)}")

    # ----------------------------------------------------------------------
    # Prime residue indices
    # ----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. PRIME RESIDUE INDEXES")
    print("=" * 78)

    residue_indexes: Dict[int, Dict[int, int]] = {}

    index_start = time.perf_counter()

    for ell in sorted(inventory):
        residue_indexes[ell] = build_residue_index(primes, ell)

        occupied = len(residue_indexes[ell])
        max_class = max(residue_indexes[ell].values())

        print(
            f"ell={ell:>5} "
            f"occupied={occupied:>5} "
            f"max_class={max_class:>6}"
        )

    print(
        f"index build time = "
        f"{time.perf_counter() - index_start:.4f}s"
    )

    # ----------------------------------------------------------------------
    # Algebraic classification of the factors
    # ----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. CYCLOTOMIC PRIME CLASSIFICATION")
    print("=" * 78)

    count_ell3 = 0
    count_ell1 = 0
    order3_checks = 0
    order3_passes = 0

    for ell, info in sorted(inventory.items()):
        if ell == 3:
            count_ell3 += 1
        elif ell % 3 == 1:
            count_ell1 += 1

        for r in info.rs:
            if ell == 3:
                continue

            order3_checks += 1

            if multiplicative_order(r, ell) == 3:
                order3_passes += 1

    print(f"ell = 3 factors              = {count_ell3}")
    print(f"ell ≡ 1 (mod 3) factors      = {count_ell1}")
    print(f"order-3 checks                = {order3_checks}")
    print(f"order-3 passes                = {order3_passes}")
    print(
        "order-3 non-exceptional pass = "
        f"{order3_passes == order3_checks}"
    )

    # ----------------------------------------------------------------------
    # LCM independence
    # ----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. INDEPENDENT PRIME-FACTOR LCM GROWTH")
    print("=" * 78)

    print(
        f"{'step':>4} "
        f"{'ell':>7} "
        f"{'new multiplier':>15} "
        f"{'LCM':>24} "
        f"{'bits':>5}"
    )
    print("-" * 78)

    current_lcm = 1

    for step, ell in enumerate(sorted(inventory), 1):
        previous = current_lcm
        current_lcm = math.lcm(current_lcm, ell)
        multiplier = current_lcm // previous

        print(
            f"{step:>4} "
            f"{ell:>7} "
            f"{multiplier:>15} "
            f"{current_lcm:>24} "
            f"{current_lcm.bit_length():>5}"
        )

    # ----------------------------------------------------------------------
    # Pairwise overlap
    # ----------------------------------------------------------------------

    overlap = pairwise_factor_overlap(inventory)

    print()
    print("=" * 78)
    print("8. PRIME-FACTOR REUSE")
    print("=" * 78)

    if not overlap:
        print("no shared F(r) factor pairs")
    else:
        print(
            f"factor pairs sharing at least one r = {len(overlap)}"
        )

        for a, b in overlap:
            shared_rs = sorted(
                set(inventory[a].rs) & set(inventory[b].rs)
            )
            print(
                f"ell={a:>5} "
                f"ell={b:>5} "
                f"shared r={shared_rs}"
            )

    # ----------------------------------------------------------------------
    # Per-target single-factor contribution
    # ----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. SINGLE PRIME-FACTOR FILTER CONTRIBUTION")
    print("=" * 78)

    print(
        "For each target and each cyclotomic prime ell, count searchable"
    )
    print(
        "prime p values for which a compatible prime residue q exists."
    )
    print()

    global_single_counts: Dict[int, List[int]] = defaultdict(list)

    # Do not print every target x ell entry.
    # Instead print compact per-target summaries.
    for target_id, target in enumerate(targets, 1):
        rows = []

        for ell in sorted(inventory):
            count = count_factor_candidates(
                primes,
                target.n,
                ell,
                residue_indexes[ell],
            )
            rows.append((ell, count))
            global_single_counts[ell].append(count)

        best = sorted(rows, key=lambda x: x[1])

        print(
            f"target {target_id:>2}: "
            f"best ell={best[0][0]:>5} "
            f"candidates={best[0][1]:>8,}; "
            f"worst ell={best[-1][0]:>5} "
            f"candidates={best[-1][1]:>8,}"
        )

    # ----------------------------------------------------------------------
    # Incremental accumulation of DISTINCT ell factors
    # ----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. INCREMENTAL INDEPENDENT-ELL FILTER")
    print("=" * 78)

    ordered_ells = sorted(inventory)

    # For each target keep a set of p candidates.
    # This is only ~147k elements, so memory is safe.
    #
    # Important:
    # We start from all primes and progressively impose:
    #
    #      p*q = n (mod ell)
    #
    # for each DISTINCT ell.
    #
    # This directly measures the contribution of independent
    # cyclotomic prime moduli.

    global_stage_survivors: Dict[int, List[int]] = defaultdict(list)
    global_stage_reduction: Dict[int, List[float]] = defaultdict(list)

    for target_id, target in enumerate(targets, 1):
        candidates: Set[int] = set(primes)

        print()
        print(
            f"TARGET {target_id:>2} "
            f"n={target.n}"
        )
        print(
            f"{'step':>4} "
            f"{'ell':>7} "
            f"{'before':>10} "
            f"{'after':>10} "
            f"{'retained':>10} "
            f"{'reduction':>12}"
        )
        print("-" * 78)

        for step, ell in enumerate(ordered_ells, 1):
            before = len(candidates)
            residue_index = residue_indexes[ell]
            n_mod = target.n % ell

            next_candidates: Set[int] = set()

            for p in candidates:
                a = p % ell

                if a == 0:
                    continue

                # All primes in our search interval exceed ell.
                # Nevertheless, retain the gcd check for correctness.
                if math.gcd(a, ell) != 1:
                    continue

                b = (n_mod * pow(a, -1, ell)) % ell

                if residue_index.get(b, 0) > 0:
                    next_candidates.add(p)

            candidates = next_candidates
            after = len(candidates)

            retained = after / before if before else 0.0
            reduction = (
                before / after
                if after
                else float("inf")
            )

            print(
                f"{step:>4} "
                f"{ell:>7} "
                f"{before:>10,} "
                f"{after:>10,} "
                f"{retained:>10.6f} "
                f"{reduction:>11.3f}x"
            )

            global_stage_survivors[ell].append(after)
            global_stage_reduction[ell].append(reduction)

        true_present = (
            target.p in candidates
            and target.q in primes
        )

        print(
            f"final survivors = {len(candidates):,} "
            f"true p retained = {target.p in candidates} "
            f"recovered-control={true_present}"
        )

    # ----------------------------------------------------------------------
    # Compare cyclotomic factor accumulation with raw F(r) accumulation
    # ----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("11. CYCLOTOMIC VS RAW-F(r) INFORMATION")
    print("=" * 78)

    raw_lcm = 1

    print(
        f"{'r':>4} "
        f"{'F(r)':>8} "
        f"{'raw LCM':>22} "
        f"{'bits':>5} "
        f"{'new prime factors':<20}"
    )
    print("-" * 78)

    seen_factor_primes: Set[int] = set()

    for r in R_VALUES:
        value = F(r)
        factors = set(sympy.factorint(value).keys())

        previous = raw_lcm
        raw_lcm = math.lcm(raw_lcm, value)

        new_factors = sorted(factors - seen_factor_primes)
        seen_factor_primes.update(factors)

        print(
            f"{r:>4} "
            f"{value:>8} "
            f"{raw_lcm:>22} "
            f"{raw_lcm.bit_length():>5} "
            f"{str(new_factors):<20}"
        )

    # ----------------------------------------------------------------------
    # Compact global summary
    # ----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("12. GLOBAL SUMMARY")
    print("=" * 78)

    print(
        f"{'ell':>7} "
        f"{'median single-filter candidates':>34} "
        f"{'median final-stage survivors':>30}"
    )
    print("-" * 78)

    for ell in ordered_ells:
        single_values = global_single_counts[ell]
        stage_values = global_stage_survivors[ell]

        median_single = (
            sorted(single_values)[len(single_values) // 2]
        )
        median_stage = (
            sorted(stage_values)[len(stage_values) // 2]
        )

        print(
            f"{ell:>7} "
            f"{median_single:>34,} "
            f"{median_stage:>30,}"
        )

    # ----------------------------------------------------------------------
    # Exact recovery control
    # ----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("13. EXACT RECOVERY CONTROL")
    print("=" * 78)

    successes = 0

    for target in targets:
        ok, tests = exact_recovery_scan(
            primes,
            target.n,
            target.p,
            target.q,
        )

        if ok:
            successes += 1

        print(
            f"n={target.n} "
            f"baseline_tests={tests:>7,} "
            f"correct={ok}"
        )

    print()
    print(f"correctly recovered = {successes}/{len(targets)}")

    # ----------------------------------------------------------------------
    # Explicit mathematical diagnostic
    # ----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("14. MATHEMATICAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "For ell != 3 with ell | F(r):"
    )
    print(
        "    r^2 + r + 1 = 0 (mod ell)"
    )
    print(
        "    r^3 = 1 (mod ell)"
    )
    print(
        "    ord_ell(r) = 3"
    )
    print(
        "    therefore ell = 1 (mod 3)"
    )

    print()
    print(
        "This experiment asks whether the independent prime factors ell"
    )
    print(
        "provide a stronger N-only filter than treating each F(r) as a"
    )
    print(
        "single composite modulus."
    )

    print()
    print(
        "The critical comparison is:"
    )
    print(
        "    raw F(r) structure"
    )
    print(
        "versus"
    )
    print(
        "    distinct cyclotomic prime factors ell"
    )

    print()
    print(
        "No prime-pair enumeration is performed."
    )
    print(
        "No CSV files are produced."
    )

    total_elapsed = time.perf_counter() - total_start

    print()
    print("=" * 78)
    print("EXPERIMENT 47 COMPLETE")
    print("=" * 78)
    print(
        f"total runtime = "
        f"{total_elapsed:.4f}s "
        f"({total_elapsed / 60:.2f} min)"
    )


if __name__ == "__main__":
    main()

