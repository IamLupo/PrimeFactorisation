#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 52
SEGMENTED QUADRATIC-RESIDUE SIEVE VS EXPLICIT FILTERING
NO CSV OUTPUT
==============================================================================

Goal
----
Experiment 51 established that the powerful effect is the discriminant:

    d^2 = s^2 - 4n

and that both the cyclotomic and random-prime modulus families behave
approximately like generic quadratic-residue filters.

Experiment 52 changes the implementation model.

Instead of:

    for every candidate s:
        for every modulus ell:
            test whether s^2 - 4n is a QR mod ell

we build a true periodic sieve:

    for each modulus ell:
        compute all admissible residues s mod ell
        strike out all forbidden positions in the interval

This measures whether the quadratic-residue information can be turned
into a genuinely sieve-like algorithm rather than repeated Python-level
membership tests.

We compare:

    A. ordinary exact square scan
    B. explicit QR filtering
    C. segmented QR sieve
    D. segmented QR sieve + exact square verification

The experiment also compares:

    cyclotomic modulus family
    random-prime control family

No prime-pair enumeration.
No CSV files.
==============================================================================
"""

from __future__ import annotations

import math
import random
import statistics
import time
from dataclasses import dataclass
from typing import Iterable


# =============================================================================
# CONFIGURATION
# =============================================================================

SEED = 20260814

TARGETS = 12

P_MIN = 2_000_000
P_MAX = 4_200_000

# Same factor pairs used throughout the preceding experiments.
KNOWN_PAIRS = [
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

CYCLOTOMIC_MODULI = [
    7, 13, 19, 31, 37, 61, 67, 79,
    127, 307, 331, 631, 1723
]

CONTROL_MODULI = [
    673, 1109, 2309, 2351, 4421, 4561, 4759,
    6211, 7879, 7951, 8273, 8689, 9781
]

# Test subsets as well as full depth.
DEPTHS = [1, 2, 3, 4, 5, 7, 9, 11, 13]

# Sieve should work on a compact bytearray segment.
# Larger values reduce Python loop overhead but consume more memory.
SEGMENT_SIZE = 1_000_000


# =============================================================================
# DATA STRUCTURES
# =============================================================================

@dataclass
class Target:
    index: int
    p: int
    q: int
    n: int
    true_s: int


@dataclass
class Result:
    name: str
    target: int
    candidates: int
    exact_tests: int
    found_s: int | None
    elapsed: float


# =============================================================================
# PRIME GENERATION
# =============================================================================

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
            sieve[start:limit:p] = b"\x00" * (((limit - 1 - start) // p) + 1)

    return [i for i, flag in enumerate(sieve) if flag]


# =============================================================================
# TARGETS
# =============================================================================

def build_targets() -> list[Target]:
    targets: list[Target] = []

    for i, (p, q) in enumerate(KNOWN_PAIRS[:TARGETS], start=1):
        if not (P_MIN <= p < P_MAX and P_MIN <= q < P_MAX):
            raise ValueError(f"Target outside prime interval: {(p, q)}")

        n = p * q
        s = p + q

        targets.append(
            Target(
                index=i,
                p=p,
                q=q,
                n=n,
                true_s=s,
            )
        )

    return targets


# =============================================================================
# DISCRIMINANT DOMAIN
# =============================================================================

def feasible_even_sum_bounds() -> tuple[int, int]:
    """
    For p,q in [P_MIN,P_MAX):

        s_min = 2*P_MIN
        s_max = 2*(P_MAX-1)

    Since p,q are odd primes here, s is even.
    """
    return 2 * P_MIN, 2 * (P_MAX - 1)


# =============================================================================
# QR TABLES
# =============================================================================

def build_qr_residue_table(modulus: int) -> bytearray:
    """
    qr[r] = 1 iff r is a quadratic residue modulo modulus.

    All moduli used here are odd primes.
    """
    qr = bytearray(modulus)

    for x in range(modulus):
        qr[(x * x) % modulus] = 1

    return qr


def build_qr_tables(moduli: Iterable[int]) -> dict[int, bytearray]:
    return {m: build_qr_residue_table(m) for m in moduli}


# =============================================================================
# BASIC EXACT CHECK
# =============================================================================

def exact_square_root(x: int) -> int | None:
    if x < 0:
        return None

    r = math.isqrt(x)
    if r * r == x:
        return r

    return None


def factor_from_sum(n: int, s: int) -> tuple[int, int] | None:
    """
    Given s=p+q, solve:

        d^2 = s^2 - 4n
        p = (s-d)/2
        q = (s+d)/2
    """
    disc = s * s - 4 * n

    if disc < 0:
        return None

    d = math.isqrt(disc)

    if d * d != disc:
        return None

    if (s - d) & 1:
        return None

    p = (s - d) // 2
    q = (s + d) // 2

    if p * q != n:
        return None

    return p, q


# =============================================================================
# METHOD A: PLAIN EXACT SUM SCAN
# =============================================================================

def plain_sum_scan(target: Target) -> tuple[int, int, float]:
    """
    Scan every even s and test whether s^2-4n is a square.
    """
    s_min, s_max = feasible_even_sum_bounds()

    tests = 0
    start = time.perf_counter()

    for s in range(s_min, s_max + 1, 2):
        tests += 1
        result = factor_from_sum(target.n, s)

        if result is not None:
            return s, tests, time.perf_counter() - start

    raise RuntimeError("Plain sum scan failed")


# =============================================================================
# METHOD B: EXPLICIT QR FILTER
# =============================================================================

def explicit_qr_filter(
    target: Target,
    moduli: list[int],
    qr_tables: dict[int, bytearray],
) -> tuple[int, int, int, float]:
    """
    Reference implementation corresponding to Experiment 51.

    It checks every candidate s against every modulus.
    """
    s_min, s_max = feasible_even_sum_bounds()

    candidates = 0
    exact_tests = 0

    start = time.perf_counter()

    for s in range(s_min, s_max + 1, 2):
        good = True

        s2 = s * s
        four_n = 4 * target.n

        for m in moduli:
            residue = (s2 - four_n) % m

            if not qr_tables[m][residue]:
                good = False
                break

        if not good:
            continue

        candidates += 1
        exact_tests += 1

        if factor_from_sum(target.n, s) is not None:
            return (
                candidates,
                exact_tests,
                s,
                time.perf_counter() - start,
            )

    raise RuntimeError("Explicit QR filter failed")


# =============================================================================
# PERIODIC SIEVE PRECOMPUTATION
# =============================================================================

def admissible_sum_residues(
    n: int,
    modulus: int,
    qr: bytearray,
) -> list[int]:
    """
    Return residues a mod modulus for which:

        a^2 - 4n

    is a quadratic residue modulo modulus.

    These are the only residues that may survive the sieve.
    """
    four_n_mod = (4 * n) % modulus

    allowed = []

    for a in range(modulus):
        disc_residue = ((a * a) - four_n_mod) % modulus

        if qr[disc_residue]:
            allowed.append(a)

    return allowed


def admissible_parity_adjusted_residues(
    n: int,
    modulus: int,
    qr: bytearray,
) -> set[int]:
    """
    Same condition, but retain residues represented by even s.

    The sieve itself operates on k where:

        s = s_min + 2k.

    So we transform the condition into:

        (s_min + 2k)^2 - 4n is QR mod modulus.

    Since modulus is odd, 2 has an inverse.
    """
    s_min, _ = feasible_even_sum_bounds()

    inv2 = pow(2, -1, modulus)

    out: set[int] = set()

    for s_residue in admissible_sum_residues(n, modulus, qr):
        k = ((s_residue - (s_min % modulus)) * inv2) % modulus
        out.add(k)

    return out


# =============================================================================
# METHOD C: SEGMENTED PERIODIC QR SIEVE
# =============================================================================

def sieve_segment(
    low_k: int,
    high_k: int,
    n: int,
    moduli: list[int],
    qr_tables: dict[int, bytearray],
) -> bytearray:
    """
    Create a live mask for k in [low_k, high_k).

    k represents even sums:

        s = s_min + 2k.

    mask[k-low_k] = 1 survives all QR tests.

    Instead of calculating a modular square for every candidate,
    this routine walks only the forbidden residue classes for each
    modulus and strikes them out periodically.
    """
    length = high_k - low_k
    live = bytearray(b"\x01") * length

    s_min, _ = feasible_even_sum_bounds()

    for modulus in moduli:
        allowed = admissible_parity_adjusted_residues(
            n, modulus, qr_tables[modulus]
        )

        allowed_set = set(allowed)

        # Small-prime moduli have relatively few residue classes.
        # We strike forbidden classes rather than compute a fresh square.
        #
        # Number of forbidden classes = modulus - len(allowed).
        forbidden = [
            r for r in range(modulus)
            if r not in allowed_set
        ]

        for r in forbidden:
            first = low_k + ((r - low_k) % modulus)

            for k in range(first, high_k, modulus):
                live[k - low_k] = 0

    return live


def segmented_qr_sieve(
    target: Target,
    moduli: list[int],
    qr_tables: dict[int, bytearray],
) -> tuple[int, int, float, int]:
    """
    Full segmented sieve over the feasible even-sum range.

    Returns:

        first recovered s
        number of survivors
        elapsed
        number of exact square tests
    """
    s_min, s_max = feasible_even_sum_bounds()

    # We enumerate even sums through:
    #
    #   s = s_min + 2k
    #
    # where k = 0,...,K-1.
    k_count = ((s_max - s_min) // 2) + 1

    survivors_total = 0
    exact_tests = 0

    start = time.perf_counter()

    for low_k in range(0, k_count, SEGMENT_SIZE):
        high_k = min(low_k + SEGMENT_SIZE, k_count)

        live = sieve_segment(
            low_k,
            high_k,
            target.n,
            moduli,
            qr_tables,
        )

        for offset, flag in enumerate(live):
            if not flag:
                continue

            survivors_total += 1

            k = low_k + offset
            s = s_min + 2 * k

            exact_tests += 1

            result = factor_from_sum(target.n, s)

            if result is not None:
                return (
                    s,
                    survivors_total,
                    time.perf_counter() - start,
                    exact_tests,
                )

    raise RuntimeError("Segmented sieve failed")


# =============================================================================
# FASTER PERIODIC VERSION
# =============================================================================

def precompute_sieve_patterns(
    n: int,
    moduli: list[int],
    qr_tables: dict[int, bytearray],
) -> dict[int, tuple[int, list[int], list[int]]]:
    """
    Precompute forbidden residue classes for k.

    For each modulus:

        modulus
        allowed residues
        forbidden residues

    This isolates all algebraic work from the main segmented loop.
    """
    patterns = {}

    s_min, _ = feasible_even_sum_bounds()

    for modulus in moduli:
        allowed = admissible_parity_adjusted_residues(
            n, modulus, qr_tables[modulus]
        )
        allowed_set = set(allowed)

        forbidden = [
            r for r in range(modulus)
            if r not in allowed_set
        ]

        patterns[modulus] = (
            modulus,
            allowed,
            forbidden,
        )

    return patterns


def segmented_qr_sieve_precomputed(
    target: Target,
    moduli: list[int],
    qr_tables: dict[int, bytearray],
) -> tuple[int, int, float, int]:
    """
    Same sieve but moves all residue-class generation outside the segment loop.
    """
    s_min, s_max = feasible_even_sum_bounds()

    k_count = ((s_max - s_min) // 2) + 1

    patterns = precompute_sieve_patterns(
        target.n,
        moduli,
        qr_tables,
    )

    survivors_total = 0
    exact_tests = 0

    start = time.perf_counter()

    for low_k in range(0, k_count, SEGMENT_SIZE):
        high_k = min(low_k + SEGMENT_SIZE, k_count)

        length = high_k - low_k
        live = bytearray(b"\x01") * length

        for modulus, _, forbidden in patterns.values():
            for r in forbidden:
                first = low_k + ((r - low_k) % modulus)

                for k in range(first, high_k, modulus):
                    live[k - low_k] = 0

        for offset, flag in enumerate(live):
            if not flag:
                continue

            survivors_total += 1

            k = low_k + offset
            s = s_min + 2 * k

            exact_tests += 1

            result = factor_from_sum(target.n, s)

            if result is not None:
                return (
                    s,
                    survivors_total,
                    time.perf_counter() - start,
                    exact_tests,
                )

    raise RuntimeError("Precomputed segmented sieve failed")


# =============================================================================
# CONTROL VALIDATION
# =============================================================================

def verify_true_sum_survives(
    target: Target,
    moduli: list[int],
    qr_tables: dict[int, bytearray],
) -> bool:
    s = target.true_s

    for m in moduli:
        disc_residue = (s * s - 4 * target.n) % m

        if not qr_tables[m][disc_residue]:
            return False

    return True


# =============================================================================
# REPORTING
# =============================================================================

def median(values: list[float]) -> float:
    return statistics.median(values)


def main() -> None:
    random.seed(SEED)

    print("=" * 78)
    print("KAPPA EXPERIMENT 52")
    print("SEGMENTED QUADRATIC-RESIDUE SIEVE VS EXPLICIT FILTERING")
    print("NO CSV OUTPUT")
    print("=" * 78)
    print()

    print("random seed      =", SEED)
    print("targets          =", TARGETS)
    print(
        "prime interval   = "
        f"[{P_MIN:,}, {P_MAX:,})"
    )
    print()

    # -------------------------------------------------------------------------
    # PRIME POPULATION
    # -------------------------------------------------------------------------

    t0 = time.perf_counter()
    primes = sieve_primes(P_MAX)
    generation_time = time.perf_counter() - t0

    prime_population = [
        p for p in primes
        if P_MIN <= p < P_MAX
    ]

    print("-" * 78)
    print("1. PRIME POPULATION")
    print("-" * 78)
    print("prime population =", f"{len(prime_population):,}")
    print("generation time  =", f"{generation_time:.6f}s")
    print()

    # -------------------------------------------------------------------------
    # TARGETS
    # -------------------------------------------------------------------------

    targets = build_targets()

    print("-" * 78)
    print("2. TARGETS")
    print("-" * 78)

    for t in targets:
        print(
            f"target {t.index:2d}: "
            f"p={t.p} q={t.q} n={t.n} s={t.true_s}"
        )

    print()

    # -------------------------------------------------------------------------
    # DOMAIN
    # -------------------------------------------------------------------------

    s_min, s_max = feasible_even_sum_bounds()
    candidate_count = ((s_max - s_min) // 2) + 1

    print("-" * 78)
    print("3. EVEN-SUM DOMAIN")
    print("-" * 78)
    print("s minimum       =", f"{s_min:,}")
    print("s maximum       =", f"{s_max:,}")
    print("candidate sums  =", f"{candidate_count:,}")
    print("representation  = s = s_min + 2k")
    print()

    # -------------------------------------------------------------------------
    # MODULI
    # -------------------------------------------------------------------------

    print("-" * 78)
    print("4. MODULUS FAMILIES")
    print("-" * 78)
    print("cyclotomic =", CYCLOTOMIC_MODULI)
    print("control    =", CONTROL_MODULI)
    print()

    # -------------------------------------------------------------------------
    # QR TABLES
    # -------------------------------------------------------------------------

    all_moduli = sorted(set(CYCLOTOMIC_MODULI + CONTROL_MODULI))

    t0 = time.perf_counter()
    qr_tables = build_qr_tables(all_moduli)
    qr_time = time.perf_counter() - t0

    print("-" * 78)
    print("5. QR TABLE PREPARATION")
    print("-" * 78)
    print("moduli =", len(all_moduli))
    print("time   =", f"{qr_time:.6f}s")
    print()

    # -------------------------------------------------------------------------
    # TRUE-SUM CROSS CHECK
    # -------------------------------------------------------------------------

    print("-" * 78)
    print("6. TRUE-SUM SURVIVAL CHECK")
    print("-" * 78)

    for family_name, moduli in [
        ("cyclotomic", CYCLOTOMIC_MODULI),
        ("control", CONTROL_MODULI),
    ]:
        ok = all(
            verify_true_sum_survives(t, moduli, qr_tables)
            for t in targets
        )

        print(
            f"{family_name:12s}: "
            f"{12 if ok else 0}/{len(targets)} true sums survive"
        )

    print()

    # -------------------------------------------------------------------------
    # BASELINE
    # -------------------------------------------------------------------------

    baseline_tests = []
    baseline_times = []

    print("-" * 78)
    print("7. PLAIN EXACT SUM SCAN")
    print("-" * 78)

    for target in targets:
        s, tests, elapsed = plain_sum_scan(target)

        baseline_tests.append(tests)
        baseline_times.append(elapsed)

        print(
            f"target {target.index:2d}: "
            f"tests={tests:9,d} "
            f"time={elapsed:.6f}s "
            f"correct={s == target.true_s}"
        )

    print(
        "\nbaseline total tests =",
        f"{sum(baseline_tests):,}"
    )
    print(
        "baseline total time  =",
        f"{sum(baseline_times):.6f}s"
    )
    print()

    # -------------------------------------------------------------------------
    # EXPLICIT FILTER
    # -------------------------------------------------------------------------

    explicit_summary: dict[str, dict[str, list[float]]] = {}

    for family_name, moduli in [
        ("cyclotomic", CYCLOTOMIC_MODULI),
        ("control", CONTROL_MODULI),
    ]:
        print("-" * 78)
        print(f"8. EXPLICIT QR FILTER: {family_name.upper()}")
        print("-" * 78)

        exact_tests = []
        times = []
        survivors = []

        for target in targets:
            cand, tests, found_s, elapsed = explicit_qr_filter(
                target,
                moduli,
                qr_tables,
            )

            exact_tests.append(tests)
            times.append(elapsed)
            survivors.append(cand)

            print(
                f"target {target.index:2d}: "
                f"survivors={cand:7,d} "
                f"exact_tests={tests:7,d} "
                f"time={elapsed:.6f}s "
                f"correct={found_s == target.true_s}"
            )

        explicit_summary[family_name] = {
            "exact_tests": exact_tests,
            "times": times,
            "survivors": survivors,
        }

        print(
            f"{family_name} median survivors = "
            f"{median(survivors):,.1f}"
        )
        print(
            f"{family_name} median exact tests = "
            f"{median(exact_tests):,.1f}"
        )
        print(
            f"{family_name} median time = "
            f"{median(times):.6f}s"
        )
        print()

    # -------------------------------------------------------------------------
    # SEGMENTED SIEVE
    # -------------------------------------------------------------------------

    sieve_summary: dict[str, dict[str, list[float]]] = {}

    for family_name, moduli in [
        ("cyclotomic", CYCLOTOMIC_MODULI),
        ("control", CONTROL_MODULI),
    ]:
        print("-" * 78)
        print(f"9. SEGMENTED QR SIEVE: {family_name.upper()}")
        print("-" * 78)

        exact_tests = []
        times = []
        survivors = []

        for target in targets:
            s, found_count, elapsed, tests = segmented_qr_sieve_precomputed(
                target,
                moduli,
                qr_tables,
            )

            exact_tests.append(tests)
            times.append(elapsed)
            survivors.append(found_count)

            print(
                f"target {target.index:2d}: "
                f"survivors_seen={found_count:7,d} "
                f"exact_tests={tests:7,d} "
                f"time={elapsed:.6f}s "
                f"correct={s == target.true_s}"
            )

        sieve_summary[family_name] = {
            "exact_tests": exact_tests,
            "times": times,
            "survivors": survivors,
        }

        print(
            f"{family_name} median survivors seen = "
            f"{median(survivors):,.1f}"
        )
        print(
            f"{family_name} median exact tests = "
            f"{median(exact_tests):,.1f}"
        )
        print(
            f"{family_name} median time = "
            f"{median(times):.6f}s"
        )
        print()

    # -------------------------------------------------------------------------
    # DEPTH ANALYSIS
    # -------------------------------------------------------------------------

    print("-" * 78)
    print("10. DEPTH ANALYSIS")
    print("-" * 78)

    for family_name, moduli in [
        ("cyclotomic", CYCLOTOMIC_MODULI),
        ("control", CONTROL_MODULI),
    ]:
        print()
        print(f"{family_name.upper()}")

        for depth in DEPTHS:
            selected = moduli[:depth]

            counts = []

            for target in targets:
                patterns = precompute_sieve_patterns(
                    target.n,
                    selected,
                    qr_tables,
                )

                # Instead of running the full sieve, estimate survivor count
                # using one complete residue period when feasible.
                #
                # We use LCM here, but cap construction because the product of
                # 13 distinct primes can become enormous.
                period = 1
                safe = True

                for m in selected:
                    period = math.lcm(period, m)

                    if period > 2_000_000:
                        safe = False
                        break

                if safe:
                    valid_k = 0

                    for k in range(period):
                        ok = True

                        for m, _, forbidden in patterns.values():
                            if (k % m) in set(forbidden):
                                ok = False
                                break

                        if ok:
                            valid_k += 1

                    fraction = valid_k / period
                    estimated = round(candidate_count * fraction)
                    counts.append(estimated)

                else:
                    # Monte Carlo estimate when the complete period is too
                    # large. Fixed seed keeps the experiment reproducible.
                    rng = random.Random(SEED + depth)
                    samples = 100_000
                    hits = 0

                    forbidden_sets = {
                        m: set(forbidden)
                        for m, _, forbidden in patterns.values()
                    }

                    for _ in range(samples):
                        k = rng.randrange(candidate_count)

                        ok = True
                        for m, forbidden in forbidden_sets.items():
                            if (k % m) in forbidden:
                                ok = False
                                break

                        if ok:
                            hits += 1

                    counts.append(
                        round(candidate_count * hits / samples)
                    )

            print(
                f"depth={depth:2d} "
                f"median estimated survivors="
                f"{median(counts):10,.1f}"
            )

    # -------------------------------------------------------------------------
    # FINAL COMPARISON
    # -------------------------------------------------------------------------

    print()
    print("-" * 78)
    print("11. FINAL COMPARISON")
    print("-" * 78)

    print(
        f"{'method':30s}"
        f"{'median survivors':>18s}"
        f"{'median exact':>18s}"
        f"{'median time':>16s}"
    )
    print("-" * 82)

    print(
        f"{'plain exact scan':30s}"
        f"{candidate_count:18,d}"
        f"{median(baseline_tests):18,.1f}"
        f"{median(baseline_times):16.6f}"
    )

    for family_name in ["cyclotomic", "control"]:
        data = explicit_summary[family_name]

        print(
            f"{'explicit QR ' + family_name:30s}"
            f"{median(data['survivors']):18,.1f}"
            f"{median(data['exact_tests']):18,.1f}"
            f"{median(data['times']):16.6f}"
        )

        data2 = sieve_summary[family_name]

        print(
            f"{'segmented sieve ' + family_name:30s}"
            f"{median(data2['survivors']):18,.1f}"
            f"{median(data2['exact_tests']):18,.1f}"
            f"{median(data2['times']):16.6f}"
        )

    # -------------------------------------------------------------------------
    # DIAGNOSTIC
    # -------------------------------------------------------------------------

    print()
    print("-" * 78)
    print("12. DIAGNOSTIC")
    print("-" * 78)

    cyclo_tests = median(sieve_summary["cyclotomic"]["exact_tests"])
    control_tests = median(sieve_summary["control"]["exact_tests"])
    baseline_med = median(baseline_tests)

    cyclo_time = median(sieve_summary["cyclotomic"]["times"])
    control_time = median(sieve_summary["control"]["times"])
    baseline_time = median(baseline_times)

    print()
    print("baseline median exact tests =",
          f"{baseline_med:,.1f}")

    print("cyclotomic median exact tests =",
          f"{cyclo_tests:,.1f}")

    print("control median exact tests =",
          f"{control_tests:,.1f}")

    print()
    print(
        "cyclotomic exact-test reduction =",
        f"{baseline_med / cyclo_tests:.2f}x"
    )

    print(
        "control exact-test reduction =",
        f"{baseline_med / control_tests:.2f}x"
    )

    print()
    print(
        "baseline median time =",
        f"{baseline_time:.6f}s"
    )

    print(
        "cyclotomic sieve median time =",
        f"{cyclo_time:.6f}s"
    )

    print(
        "control sieve median time =",
        f"{control_time:.6f}s"
    )

    print()
    print("The decisive question is now:")
    print()
    print("    Can periodic QR sieving reduce runtime,")
    print("    not merely the number of exact square tests?")
    print()
    print("A second question is:")
    print()
    print("    Does the cyclotomic family retain any measurable")
    print("    advantage over an equal-sized random-prime family")
    print("    once both are implemented as genuine periodic sieves?")
    print()
    print("If both families have approximately the same survivor")
    print("density and the random family is equally fast or faster,")
    print("the x^2+x+1 structure is probably not providing the")
    print("computational advantage.")
    print()
    print("If the segmented sieve is substantially faster than the")
    print("explicit filter, we have isolated an implementation path")
    print("that deserves scaling experiments.")
    print()
    print("If the sieve remains slower than the plain scan despite")
    print("thousands-fold reduction in exact square tests, then the")
    print("main bottleneck is candidate marking / memory traffic,")
    print("not the arithmetic test itself.")
    print()

    print("=" * 78)
    print("EXPERIMENT 52 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

