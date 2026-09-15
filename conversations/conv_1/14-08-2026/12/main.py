#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 55
DIRECT SUM-RESIDUE RECONSTRUCTION / MODULAR LIFTING TREE
NO CSV OUTPUT
==============================================================================

PURPOSE
-------
The previous experiments established that the discriminant sieve

    d^2 = s^2 - 4n

is extremely effective, but we were still finding s by scanning the
2,200,000 possible even sums.

This experiment asks a different question:

    Can we reconstruct the surviving s values directly from their
    modular residue classes, rather than scanning every s?

We build the admissible residue set incrementally:

    S_0
      -> S_1 mod m_1
      -> S_2 mod (m_1*m_2)
      -> ...
      -> S_k

and compare the amount of residue information against the size of the
original interval.

IMPORTANT:
    This is a diagnostic experiment.
    It deliberately stops the modular tree before memory growth becomes
    excessive and compares several modulus orderings.

The experiment measures:

    1. residue-tree growth
    2. interval-compatible residue counts
    3. point at which the cumulative modulus exceeds the search interval
    4. whether direct modular reconstruction becomes competitive
    5. whether a better modulus ordering materially changes the result

No candidate-pair enumeration.
No CSV output.
==============================================================================
"""

from __future__ import annotations

import math
import random
import time
from dataclasses import dataclass
from typing import Iterable, Sequence


SEED = 20260814

LO = 2_000_000
HI = 4_200_000

R_VALUES = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
    31, 37, 41, 43, 47
]

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

SEARCH_S_MIN = LO * 2
SEARCH_S_MAX = HI * 2 - 2
EVEN_SUM_COUNT = (SEARCH_S_MAX - SEARCH_S_MIN) // 2 + 1

# Prevent accidental runaway memory use.
MAX_TREE_STATES = 2_500_000


# ---------------------------------------------------------------------------
# Basic arithmetic
# ---------------------------------------------------------------------------

def primes_upto(n: int) -> list[int]:
    """Simple sieve."""
    if n < 2:
        return []

    sieve = bytearray(b"\x01") * (n + 1)
    sieve[:2] = b"\x00\x00"

    limit = int(math.isqrt(n))
    for p in range(2, limit + 1):
        if sieve[p]:
            start = p * p
            sieve[start:n + 1:p] = b"\x00" * (
                ((n - start) // p) + 1
            )

    return [i for i, flag in enumerate(sieve) if flag]


def isqrt_exact(x: int) -> int:
    """Integer square root with an exact-square test."""
    if x < 0:
        return -1
    r = math.isqrt(x)
    return r if r * r == x else -1


def factor_is_valid(p: int, q: int, n: int) -> bool:
    return p * q == n and p >= LO and p < HI and q >= LO and q < HI


def exact_sum_recovery(n: int, candidate_sums: Iterable[int]) -> tuple[int | None, int]:
    """
    Test candidate sums exactly.

    For each s:
        d^2 = s^2 - 4n
        p = (s-d)/2
        q = (s+d)/2

    Returns:
        (recovered_sum, number_of_exact_tests)
    """
    tests = 0

    for s in candidate_sums:
        d2 = s * s - 4 * n
        if d2 < 0:
            continue

        d = isqrt_exact(d2)
        if d < 0:
            tests += 1
            continue

        if ((s - d) & 1) or ((s + d) & 1):
            tests += 1
            continue

        p = (s - d) // 2
        q = (s + d) // 2

        tests += 1

        if factor_is_valid(p, q, n):
            return s, tests

    return None, tests


# ---------------------------------------------------------------------------
# Quadratic-residue admissible sum residues
# ---------------------------------------------------------------------------

def qr_table(modulus: int) -> list[bool]:
    """
    QR[a] is True when a is a quadratic residue mod modulus.

    We include zero as a residue.
    """
    table = [False] * modulus

    for x in range(modulus):
        table[(x * x) % modulus] = True

    return table


def allowed_sum_residues(n: int, modulus: int, qr: Sequence[bool]) -> list[int]:
    """
    Return residues s mod modulus satisfying

        s^2 - 4n is a quadratic residue mod modulus.
    """
    c = (4 * n) % modulus

    out = []

    for s in range(modulus):
        if qr[(s * s - c) % modulus]:
            out.append(s)

    return out


# ---------------------------------------------------------------------------
# CRT
# ---------------------------------------------------------------------------

def crt_pair(a: int, m: int, b: int, p: int) -> tuple[int, int]:
    """
    CRT for coprime moduli.

    x = a mod m
    x = b mod p

    returns x mod m*p.
    """
    # x = a + m*k
    # m*k = b-a mod p
    k = ((b - a) * pow(m, -1, p)) % p

    new_mod = m * p
    x = (a + m * k) % new_mod

    return x, new_mod


# ---------------------------------------------------------------------------
# Direct modular lifting
# ---------------------------------------------------------------------------

@dataclass
class LiftStage:
    stage: int
    modulus: int
    residue_count: int
    interval_hits: int
    state_time: float
    skipped: bool = False


def interval_hits_for_residues(
    residues: Sequence[int],
    modulus: int,
    s_min: int,
    s_max: int,
) -> int:
    """
    Count actual even sums in the interval satisfying x == r mod modulus.

    The interval is small enough that this can be computed analytically
    without enumerating every point.
    """
    count = 0

    for r in residues:
        # Find smallest x >= s_min with x == r mod modulus.
        delta = (r - s_min) % modulus
        first = s_min + delta

        if first > s_max:
            continue

        count += (s_max - first) // modulus + 1

    return count


def lift_residue_tree(
    n: int,
    moduli: Sequence[int],
    state_limit: int = MAX_TREE_STATES,
) -> tuple[list[LiftStage], list[int] | None]:
    """
    Build the admissible residue classes by CRT lifting.

    At each stage:

        x satisfies all previous QR constraints
        and additionally satisfies the new QR constraint.

    The number of states may explode. Once state_limit is exceeded,
    the experiment records the stage and stops.

    Returned final residues are modulo the final cumulative modulus.
    """
    residues = [0]
    modulus = 1

    stages: list[LiftStage] = []

    for stage_no, p in enumerate(moduli, start=1):
        t0 = time.perf_counter()

        qr = qr_table(p)
        allowed = allowed_sum_residues(n, p, qr)

        # For the first modulus, simply use the allowed residues.
        if modulus == 1:
            new_residues = allowed
            new_modulus = p
        else:
            new_modulus = modulus * p

            # Predict state count before constructing it.
            predicted = len(residues) * len(allowed)

            if predicted > state_limit:
                elapsed = time.perf_counter() - t0

                stages.append(
                    LiftStage(
                        stage=stage_no,
                        modulus=new_modulus,
                        residue_count=predicted,
                        interval_hits=0,
                        state_time=elapsed,
                        skipped=True,
                    )
                )

                break

            new_residues = []

            inv_modulus = pow(modulus, -1, p)

            # Optimized CRT loop.
            for a in residues:
                amod = a % p

                for b in allowed:
                    k = ((b - amod) * inv_modulus) % p
                    x = a + modulus * k
                    new_residues.append(x)

        residues = new_residues
        modulus = new_modulus

        interval_hits = interval_hits_for_residues(
            residues,
            modulus,
            SEARCH_S_MIN,
            SEARCH_S_MAX,
        )

        elapsed = time.perf_counter() - t0

        stages.append(
            LiftStage(
                stage=stage_no,
                modulus=modulus,
                residue_count=len(residues),
                interval_hits=interval_hits,
                state_time=elapsed,
                skipped=False,
            )
        )

    return stages, residues


# ---------------------------------------------------------------------------
# Ordering heuristics
# ---------------------------------------------------------------------------

def expected_single_modulus_survival(modulus: int) -> float:
    """
    For quadratic-residue filtering on a prime modulus, approximately half
    the residues survive.

    This function is intentionally simple: we are comparing ordering
    strategies, not modeling exact local statistics.
    """
    if modulus <= 2:
        return 1.0
    return 0.5


def cumulative_modulus(moduli: Sequence[int]) -> int:
    out = 1
    for m in moduli:
        out *= m
    return out


def print_ordering_diagnostic(name: str, moduli: Sequence[int]) -> None:
    M = cumulative_modulus(moduli)

    print()
    print(f"ORDERING: {name}")
    print(f"moduli = {list(moduli)}")
    print(f"final formal modulus bits = {M.bit_length()}")
    print(f"first modulus product exceeding search interval = ", end="")

    prod = 1
    found = False
    for i, m in enumerate(moduli, start=1):
        prod *= m
        if prod > EVEN_SUM_COUNT and not found:
            print(f"stage {i} (M={prod:,})")
            found = True

    if not found:
        print("never")


# ---------------------------------------------------------------------------
# Target diagnostics
# ---------------------------------------------------------------------------

def print_target_header(index: int, p: int, q: int) -> int:
    n = p * q
    true_s = p + q

    print()
    print("=" * 78)
    print(f"TARGET {index:2d}")
    print("=" * 78)
    print(f"p={p} q={q} n={n} true_s={true_s}")

    return n


def print_stages(stages: Sequence[LiftStage]) -> None:
    print()
    print(
        "stage modulus                           "
        "residues       interval_hits       time"
    )
    print("-" * 78)

    for st in stages:
        if st.skipped:
            print(
                f"{st.stage:5d} {st.modulus:>32,} "
                f"{st.residue_count:>14,} "
                f"{'STOPPED':>16} "
                f"{st.state_time:>10.6f}s"
            )
        else:
            print(
                f"{st.stage:5d} {st.modulus:>32,} "
                f"{st.residue_count:>14,} "
                f"{st.interval_hits:>16,} "
                f"{st.state_time:>10.6f}s"
            )


# ---------------------------------------------------------------------------
# Main experiment
# ---------------------------------------------------------------------------

def main() -> None:
    random.seed(SEED)

    print("=" * 78)
    print("KAPPA EXPERIMENT 55")
    print("DIRECT SUM-RESIDUE RECONSTRUCTION / MODULAR LIFTING TREE")
    print("NO CSV OUTPUT")
    print("=" * 78)

    print()
    print("random seed            =", SEED)
    print("prime interval         =", f"[{LO:,}, {HI:,})")
    print("candidate sums         =", f"{EVEN_SUM_COUNT:,}")
    print("sum interval           =", f"[{SEARCH_S_MIN:,}, {SEARCH_S_MAX:,}]")

    # ----------------------------------------------------------------------
    # Prime population
    # ----------------------------------------------------------------------
    t0 = time.perf_counter()
    primes = primes_upto(HI - 1)
    primes = [p for p in primes if p >= LO]
    prime_generation = time.perf_counter() - t0

    print()
    print("prime population       =", f"{len(primes):,}")
    print("generation time        =", f"{prime_generation:.6f}s")

    # ----------------------------------------------------------------------
    # Modulus families
    # ----------------------------------------------------------------------
    cyclo = [r * r + r + 1 for r in R_VALUES]

    # The previous experiments found these prime factors.
    # We deliberately include only the non-exceptional prime factors.
    base_cyclo = [
        7, 13, 19, 31, 37, 61, 67,
        79, 127, 307, 331, 631, 1723
    ]

    # Add 3 separately if desired; it is a degenerate/order-1 factor.
    base_cyclo_no3 = base_cyclo

    orderings = {
        "ascending": tuple(sorted(base_cyclo_no3)),
        "descending": tuple(sorted(base_cyclo_no3, reverse=True)),
        "small-first": tuple(base_cyclo_no3),
        "large-first": tuple(reversed(base_cyclo_no3)),
    }

    print()
    print("base cyclotomic primes =", base_cyclo_no3)

    for name, mods in orderings.items():
        print_ordering_diagnostic(name, mods)

    # ----------------------------------------------------------------------
    # Main target loop
    # ----------------------------------------------------------------------
    all_summary = []

    for index, (p_true, q_true) in enumerate(TARGETS, start=1):
        n = print_target_header(index, p_true, q_true)

        # --------------------------------------------------------------
        # Baseline full sum scan
        # --------------------------------------------------------------
        t0 = time.perf_counter()

        baseline_sums = range(
            SEARCH_S_MIN,
            SEARCH_S_MAX + 1,
            2,
        )

        baseline_s, baseline_tests = exact_sum_recovery(
            n,
            baseline_sums,
        )

        baseline_time = time.perf_counter() - t0

        print()
        print("BASELINE")
        print(
            f"recovered_s={baseline_s} "
            f"tests={baseline_tests:,} "
            f"time={baseline_time:.6f}s "
            f"correct={baseline_s == p_true + q_true}"
        )

        target_summary = {
            "baseline_tests": baseline_tests,
            "baseline_time": baseline_time,
            "orders": {},
        }

        # --------------------------------------------------------------
        # Modular lifting order comparisons
        # --------------------------------------------------------------
        for name, mods in orderings.items():
            print()
            print("-" * 78)
            print(f"MODULAR RECONSTRUCTION: {name}")

            t0 = time.perf_counter()

            stages, final_residues = lift_residue_tree(
                n,
                mods,
                state_limit=MAX_TREE_STATES,
            )

            total_tree_time = time.perf_counter() - t0

            print_stages(stages)

            completed = bool(stages) and not stages[-1].skipped

            if final_residues is None or not completed:
                print()
                print(
                    "TREE STATUS: STOPPED BEFORE COMPLETE RECONSTRUCTION"
                )
                print(
                    "This is itself a useful result: direct CRT lifting "
                    "becomes combinatorially expensive."
                )

                target_summary["orders"][name] = {
                    "completed": False,
                    "tree_time": total_tree_time,
                }

                continue

            # ----------------------------------------------------------
            # Convert final residue classes into interval candidates.
            # ----------------------------------------------------------
            final_modulus = stages[-1].modulus

            interval_candidates = []

            for r in final_residues:
                # Smallest value in interval with this residue.
                delta = (r - SEARCH_S_MIN) % final_modulus
                first = SEARCH_S_MIN + delta

                if first <= SEARCH_S_MAX and (first & 1) == 0:
                    interval_candidates.append(first)

                # If modulus is below interval, additional lifts exist.
                # They are enumerated here only for diagnostic purposes.
                if final_modulus <= SEARCH_S_MAX - SEARCH_S_MIN:
                    x = first + final_modulus

                    while x <= SEARCH_S_MAX:
                        if (x & 1) == 0:
                            interval_candidates.append(x)
                        x += final_modulus

            interval_candidates = sorted(set(interval_candidates))

            print()
            print(
                f"INTERVAL CANDIDATES = {len(interval_candidates):,}"
            )

            t1 = time.perf_counter()

            recovered_s, exact_tests = exact_sum_recovery(
                n,
                interval_candidates,
            )

            exact_time = time.perf_counter() - t1

            total_time = total_tree_time + exact_time

            print(
                f"RECOVERY "
                f"s={recovered_s} "
                f"exact_tests={exact_tests:,} "
                f"tree_time={total_tree_time:.6f}s "
                f"exact_time={exact_time:.6f}s "
                f"total={total_time:.6f}s "
                f"correct={recovered_s == p_true + q_true}"
            )

            target_summary["orders"][name] = {
                "completed": True,
                "tree_time": total_tree_time,
                "candidate_count": len(interval_candidates),
                "exact_tests": exact_tests,
                "total_time": total_time,
            }

        all_summary.append(target_summary)

    # ----------------------------------------------------------------------
    # Global summary
    # ----------------------------------------------------------------------
    print()
    print("=" * 78)
    print("GLOBAL SUMMARY")
    print("=" * 78)

    def median(values: Sequence[float]) -> float:
        if not values:
            return float("nan")
        values = sorted(values)
        n = len(values)
        if n & 1:
            return values[n // 2]
        return (values[n // 2 - 1] + values[n // 2]) / 2

    baseline_tests = [
        x["baseline_tests"]
        for x in all_summary
    ]

    baseline_times = [
        x["baseline_time"]
        for x in all_summary
    ]

    print()
    print(
        f"baseline median exact tests = "
        f"{median(baseline_tests):,.1f}"
    )
    print(
        f"baseline median time        = "
        f"{median(baseline_times):.6f}s"
    )

    print()
    print(
        "ordering              completed   median tree time   "
        "median candidates   median exact"
    )
    print("-" * 78)

    for name in orderings:
        completed_rows = [
            x["orders"][name]
            for x in all_summary
            if x["orders"].get(name, {}).get("completed")
        ]

        if not completed_rows:
            print(
                f"{name:20s} {'0/12':>10s} "
                f"{'--':>18s} {'--':>19s} {'--':>15s}"
            )
            continue

        tree_times = [
            x["tree_time"] for x in completed_rows
        ]

        candidates = [
            x["candidate_count"] for x in completed_rows
        ]

        exacts = [
            x["exact_tests"] for x in completed_rows
        ]

        print(
            f"{name:20s} "
            f"{len(completed_rows)}/12"
            f"{median(tree_times):18.6f}s "
            f"{median(candidates):19,.1f} "
            f"{median(exacts):15,.1f}"
        )

    # ----------------------------------------------------------------------
    # Mathematical interpretation
    # ----------------------------------------------------------------------
    print()
    print("=" * 78)
    print("MATHEMATICAL DIAGNOSTIC")
    print("=" * 78)

    print(
        """
The reconstruction tree is testing whether the QR information can be
converted directly into a sparse set of CRT residue classes for s.

At each prime ell:

    s^2 - 4n must be a quadratic residue mod ell.

If constraints were fully independent, approximately half of all residue
classes would survive at every new prime.

Therefore after k independent primes we expect roughly

    2.2e6 / 2^k

candidate sums.

The experiment asks whether CRT lifting actually provides a practical
reconstruction mechanism.

Three outcomes are especially important:

A. SUCCESS
   The modular tree remains small enough to reach the full modulus and
   reconstruct s directly with very few interval candidates.

B. COMBINATORIAL WALL
   The number of CRT residue states explodes before the modulus becomes
   useful. This means the information exists but is not naturally
   reconstructible through naive CRT lifting.

C. ORDERING EFFECT
   One modulus ordering is substantially better than another. That would
   suggest an adaptive reconstruction strategy.

A negative result is still valuable: it would show that the QR sieve is
best understood as a fast filtering mechanism rather than a direct CRT
decoder.

No prime-pair enumeration.
No CSV output.
"""
    )

    print("=" * 78)
    print("EXPERIMENT 55 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

