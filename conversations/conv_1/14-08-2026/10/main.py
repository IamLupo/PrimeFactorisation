#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 53
BITSET / VECTORIZED DISCRIMINANT SIEVE
PYTHON FILTER COST VS BULK MASK OPERATIONS
NO CSV OUTPUT
==============================================================================

Goal
----
Experiment 52 showed:

    ~10^6 -> ~10^2 candidate sums

but the QR filter was slower than the plain scan because the filtering
itself is expensive in Python.

Experiment 53 therefore keeps exactly the same mathematical condition:

    d^2 = s^2 - 4n

and asks whether the filtering can be made cheap with bulk operations.

Three implementations are benchmarked:

    A. plain Python exact scan
    B. NumPy vectorized QR filtering
    C. NumPy packed-bit / boolean-mask intersection

Both cyclotomic and random-control modulus families are tested.

Important:
    This experiment does NOT claim the cyclotomic family is special.
    The primary question is implementation speed.

No prime-pair enumeration.
No CSV output.
==============================================================================

"""

from __future__ import annotations

import math
import random
import statistics
import time
from dataclasses import dataclass
from typing import Iterable

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit(
        "Experiment 53 requires NumPy.\n"
        "Install with: python -m pip install numpy"
    ) from exc


# =============================================================================
# CONFIG
# =============================================================================

SEED = 20260814

TARGETS = 12

PRIME_MIN = 2_000_000
PRIME_MAX = 4_200_000

R_VALUES = [
    2, 3, 5, 7, 11, 13, 17, 19,
    23, 29, 31, 37, 41, 43, 47
]

CYCLOTOMIC = [
    7, 13, 19, 31, 37, 61, 67,
    79, 127, 307, 331, 631, 1723
]

CONTROL = [
    673, 1109, 2309, 2351, 4421, 4561, 4759,
    6211, 7879, 7951, 8273, 8689, 9781
]

RUN_CONTROL = True
RUN_SEGMENTED = True
RUN_PACKED = True

# Number of times to repeat the vectorized benchmark.
VECTOR_REPEATS = 3

# Segment size for the segmented implementation.
SEGMENT_SIZE = 250_000

# Limit output size.
PRINT_PER_TARGET = True


# =============================================================================
# DATA
# =============================================================================

@dataclass(frozen=True)
class Target:
    p: int
    q: int

    @property
    def n(self) -> int:
        return self.p * self.q

    @property
    def s(self) -> int:
        return self.p + self.q


# We reuse the established 12-target corpus from Experiments 45-52.
TARGET_DATA = [
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


# =============================================================================
# PRIME GENERATION
# =============================================================================

def generate_primes(lo: int, hi: int) -> list[int]:
    """
    Simple segmented sieve.
    """
    if hi <= lo:
        return []

    limit = math.isqrt(hi - 1) + 1

    base = bytearray(b"\x01") * limit
    base[0:2] = b"\x00\x00"

    for p in range(2, math.isqrt(limit - 1) + 1):
        if base[p]:
            start = p * p
            base[start:limit:p] = b"\x00" * (
                ((limit - 1 - start) // p) + 1
            )

    small_primes = [i for i, flag in enumerate(base) if flag]

    segment_size = 1_000_000
    primes: list[int] = []

    for left in range(lo, hi, segment_size):
        right = min(left + segment_size, hi)
        block = bytearray(b"\x01") * (right - left)

        for p in small_primes:
            if p * p >= right:
                break

            start = max(p * p, ((left + p - 1) // p) * p)
            if start >= right:
                continue

            block[start - left : right - left : p] = b"\x00" * (
                ((right - 1 - start) // p) + 1
            )

        primes.extend(
            left + i
            for i, flag in enumerate(block)
            if flag and left + i >= 2
        )

    return primes


# =============================================================================
# QR TABLES
# =============================================================================

def quadratic_residue_table(m: int) -> np.ndarray:
    """
    Boolean QR table:
        qr[x] = True iff x is a quadratic residue mod m.

    m is prime.
    """
    table = np.zeros(m, dtype=np.bool_)

    # 0 is a square.
    table[0] = True

    x = np.arange(1, m, dtype=np.int64)
    table[(x * x) % m] = True

    return table


def prepare_qr_tables(moduli: Iterable[int]) -> dict[int, np.ndarray]:
    tables: dict[int, np.ndarray] = {}

    t0 = time.perf_counter()

    for m in moduli:
        tables[m] = quadratic_residue_table(m)

    elapsed = time.perf_counter() - t0
    print(f"QR table preparation = {elapsed:.6f}s")

    return tables


# =============================================================================
# SUM DOMAIN
# =============================================================================

def build_sum_domain() -> tuple[int, int, np.ndarray]:
    """
    All possible even sums p+q where

        p,q >= 2,000,000
        p,q <  4,200,000

    Hence:

        4,000,000 <= s <= 8,399,998
        s even
    """
    s_min = 2 * PRIME_MIN
    s_max = 2 * (PRIME_MAX - 1)

    sums = np.arange(
        s_min,
        s_max + 1,
        2,
        dtype=np.int64,
    )

    return s_min, s_max, sums


# =============================================================================
# BASELINE
# =============================================================================

def baseline_sum_scan(
    n: int,
    sums: np.ndarray,
) -> tuple[int | None, int]:
    """
    Exact Python scan over all possible sums.

    For each s:
        d^2 = s^2 - 4n
    and d must be a nonnegative perfect square.

    Then:
        p = (s-d)/2
        q = (s+d)/2

    Since sqrt computation is performed only after checking d2 >= 0,
    this is the direct baseline used throughout this line of experiments.
    """
    four_n = 4 * n
    tests = 0

    for s in sums.tolist():
        d2 = s * s - four_n
        if d2 < 0:
            continue

        tests += 1

        d = math.isqrt(d2)
        if d * d != d2:
            continue

        p = (s - d) // 2
        q = (s + d) // 2

        if p * q == n:
            return s, tests

    return None, tests


# =============================================================================
# NUMPY VECTOR FILTER
# =============================================================================

def vectorized_filter(
    n: int,
    sums: np.ndarray,
    moduli: list[int],
    qr_tables: dict[int, np.ndarray],
) -> np.ndarray:
    """
    Build a boolean mask over all sums.

    A sum survives modulus m iff:

        s^2 - 4n

    is a quadratic residue modulo m.

    Everything is done using NumPy arrays.
    """
    mask = np.ones(sums.size, dtype=np.bool_)

    for m in moduli:
        residues = ((sums % m) ** 2 - (4 * n) % m) % m
        mask &= qr_tables[m][residues]

        if not mask.any():
            break

    return sums[mask]


def exact_check_candidates(
    n: int,
    candidates: np.ndarray,
) -> tuple[int | None, int]:
    """
    Exact square test over already-filtered candidate sums.
    """
    four_n = 4 * n
    tests = 0

    for s in candidates.tolist():
        d2 = int(s) * int(s) - four_n
        if d2 < 0:
            continue

        tests += 1

        d = math.isqrt(d2)
        if d * d != d2:
            continue

        p = (int(s) - d) // 2
        q = (int(s) + d) // 2

        if p * q == n:
            return int(s), tests

    return None, tests


# =============================================================================
# PACKED-BIT IMPLEMENTATION
# =============================================================================

def packed_filter(
    n: int,
    sums: np.ndarray,
    moduli: list[int],
    qr_tables: dict[int, np.ndarray],
) -> np.ndarray:
    """
    Boolean-mask intersection followed by packbits.

    We keep the logical mask in NumPy and compress it only between stages.

    This benchmark is primarily intended to measure whether reducing the
    memory footprint helps the repeated-mask workload.
    """
    mask = np.ones(sums.size, dtype=np.bool_)

    for m in moduli:
        residues = ((sums % m) ** 2 - (4 * n) % m) % m
        mask &= qr_tables[m][residues]

        if not mask.any():
            break

    packed = np.packbits(mask)
    unpacked = np.unpackbits(packed)[:sums.size]

    return sums[unpacked]


# =============================================================================
# SEGMENTED NUMPY IMPLEMENTATION
# =============================================================================

def segmented_filter(
    n: int,
    sums: np.ndarray,
    moduli: list[int],
    qr_tables: dict[int, np.ndarray],
    segment_size: int = SEGMENT_SIZE,
) -> np.ndarray:
    """
    Process the sum domain in chunks.

    This avoids materializing every intermediate candidate array at once.
    """
    survivors: list[np.ndarray] = []

    for start in range(0, sums.size, segment_size):
        stop = min(start + segment_size, sums.size)
        seg = sums[start:stop]

        mask = np.ones(seg.size, dtype=np.bool_)

        for m in moduli:
            residues = ((seg % m) ** 2 - (4 * n) % m) % m
            mask &= qr_tables[m][residues]

            if not mask.any():
                break

        if mask.any():
            survivors.append(seg[mask])

    if not survivors:
        return np.empty(0, dtype=np.int64)

    return np.concatenate(survivors)


# =============================================================================
# BIT-PACK MICROBENCHMARK
# =============================================================================

def mask_only_benchmark(
    n: int,
    sums: np.ndarray,
    moduli: list[int],
    qr_tables: dict[int, np.ndarray],
    repeats: int = VECTOR_REPEATS,
) -> tuple[float, float, int]:
    """
    Benchmark just the QR filtering operation, excluding exact square tests.
    """
    normal_times = []
    packed_times = []

    normal_count = 0

    for _ in range(repeats):
        t0 = time.perf_counter()
        candidates = vectorized_filter(
            n,
            sums,
            moduli,
            qr_tables,
        )
        normal_times.append(time.perf_counter() - t0)
        normal_count = int(candidates.size)

        t0 = time.perf_counter()
        _ = packed_filter(
            n,
            sums,
            moduli,
            qr_tables,
        )
        packed_times.append(time.perf_counter() - t0)

    return (
        statistics.median(normal_times),
        statistics.median(packed_times),
        normal_count,
    )


# =============================================================================
# TARGET RUNNER
# =============================================================================

def run_target(
    index: int,
    target: Target,
    sums: np.ndarray,
    cyclo_tables: dict[int, np.ndarray],
    control_tables: dict[int, np.ndarray],
):
    n = target.n

    print()
    print("-" * 78)
    print(
        f"TARGET {index:2d} "
        f"p={target.p} q={target.q} "
        f"n={n} s={target.s}"
    )
    print("-" * 78)

    # -------------------------------------------------------------------------
    # A. Baseline
    # -------------------------------------------------------------------------
    t0 = time.perf_counter()
    baseline_s, baseline_tests = baseline_sum_scan(n, sums)
    baseline_time = time.perf_counter() - t0

    baseline_ok = baseline_s == target.s

    print(
        f"baseline:"
        f" tests={baseline_tests:9,d}"
        f" time={baseline_time:.6f}s"
        f" recovered={baseline_s}"
        f" correct={baseline_ok}"
    )

    # -------------------------------------------------------------------------
    # B. Vectorized cyclotomic
    # -------------------------------------------------------------------------
    t0 = time.perf_counter()
    cyclo_candidates = vectorized_filter(
        n,
        sums,
        CYCLOTOMIC,
        cyclo_tables,
    )
    filter_time = time.perf_counter() - t0

    t1 = time.perf_counter()
    recovered_s, exact_tests = exact_check_candidates(
        n,
        cyclo_candidates,
    )
    exact_time = time.perf_counter() - t1

    total_time = filter_time + exact_time
    correct = recovered_s == target.s

    print()
    print("vectorized cyclotomic:")
    print(
        f"  candidates={cyclo_candidates.size:8,d}"
        f" filter={filter_time:.6f}s"
        f" exact_tests={exact_tests:6,d}"
        f" exact={exact_time:.6f}s"
        f" total={total_time:.6f}s"
        f" correct={correct}"
    )

    # -------------------------------------------------------------------------
    # C. Packed-bit cyclotomic
    # -------------------------------------------------------------------------
    packed_candidates = np.empty(0, dtype=np.int64)
    packed_filter_time = None
    packed_exact_time = None
    packed_total_time = None
    packed_correct = False

    if RUN_PACKED:
        t0 = time.perf_counter()
        packed_candidates = packed_filter(
            n,
            sums,
            CYCLOTOMIC,
            cyclo_tables,
        )
        packed_filter_time = time.perf_counter() - t0

        t1 = time.perf_counter()
        packed_s, packed_tests = exact_check_candidates(
            n,
            packed_candidates,
        )
        packed_exact_time = time.perf_counter() - t1
        packed_total_time = (
            packed_filter_time + packed_exact_time
        )
        packed_correct = packed_s == target.s

        print()
        print("packed-bit cyclotomic:")
        print(
            f"  candidates={packed_candidates.size:8,d}"
            f" filter={packed_filter_time:.6f}s"
            f" exact_tests={packed_tests:6,d}"
            f" exact={packed_exact_time:.6f}s"
            f" total={packed_total_time:.6f}s"
            f" correct={packed_correct}"
        )

    # -------------------------------------------------------------------------
    # D. Segmented cyclotomic
    # -------------------------------------------------------------------------
    if RUN_SEGMENTED:
        t0 = time.perf_counter()
        segmented_candidates = segmented_filter(
            n,
            sums,
            CYCLOTOMIC,
            cyclo_tables,
        )
        segmented_filter_time = time.perf_counter() - t0

        t1 = time.perf_counter()
        segmented_s, segmented_tests = exact_check_candidates(
            n,
            segmented_candidates,
        )
        segmented_exact_time = time.perf_counter() - t1

        segmented_total_time = (
            segmented_filter_time + segmented_exact_time
        )

        print()
        print("segmented cyclotomic:")
        print(
            f"  candidates={segmented_candidates.size:8,d}"
            f" filter={segmented_filter_time:.6f}s"
            f" exact_tests={segmented_tests:6,d}"
            f" exact={segmented_exact_time:.6f}s"
            f" total={segmented_total_time:.6f}s"
            f" correct={segmented_s == target.s}"
        )

    # -------------------------------------------------------------------------
    # E. Control
    # -------------------------------------------------------------------------
    if RUN_CONTROL:
        t0 = time.perf_counter()
        control_candidates = vectorized_filter(
            n,
            sums,
            CONTROL,
            control_tables,
        )
        control_filter_time = time.perf_counter() - t0

        t1 = time.perf_counter()
        control_s, control_tests = exact_check_candidates(
            n,
            control_candidates,
        )
        control_exact_time = time.perf_counter() - t1

        control_total = (
            control_filter_time + control_exact_time
        )

        print()
        print("vectorized control:")
        print(
            f"  candidates={control_candidates.size:8,d}"
            f" filter={control_filter_time:.6f}s"
            f" exact_tests={control_tests:6,d}"
            f" exact={control_exact_time:.6f}s"
            f" total={control_total:.6f}s"
            f" correct={control_s == target.s}"
        )

    return {
        "baseline_tests": baseline_tests,
        "baseline_time": baseline_time,
        "cyclo_candidates": int(cyclo_candidates.size),
        "cyclo_filter_time": filter_time,
        "cyclo_exact_tests": exact_tests,
        "cyclo_exact_time": exact_time,
        "cyclo_total_time": total_time,
        "packed_candidates": int(packed_candidates.size),
        "packed_total_time": packed_total_time,
    }


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    random.seed(SEED)

    print("=" * 78)
    print("KAPPA EXPERIMENT 53")
    print("BITSET / VECTORIZED DISCRIMINANT SIEVE")
    print("PYTHON FILTER COST VS BULK MASK OPERATIONS")
    print("NO CSV OUTPUT")
    print("=" * 78)

    # -------------------------------------------------------------------------
    # 1. Targets
    # -------------------------------------------------------------------------
    targets = [Target(p, q) for p, q in TARGET_DATA]

    print()
    print("-" * 78)
    print("1. TARGETS")
    print("-" * 78)

    for i, t in enumerate(targets, 1):
        print(
            f"target {i:2d}: "
            f"p={t.p} q={t.q} "
            f"n={t.n} s={t.s}"
        )

    # -------------------------------------------------------------------------
    # 2. Prime population
    # -------------------------------------------------------------------------
    print()
    print("-" * 78)
    print("2. PRIME POPULATION")
    print("-" * 78)

    t0 = time.perf_counter()
    primes = generate_primes(PRIME_MIN, PRIME_MAX)
    prime_generation_time = time.perf_counter() - t0

    print(f"prime population = {len(primes):,}")
    print(f"generation time  = {prime_generation_time:.6f}s")

    # Keep variable alive explicitly for experimental control.
    _ = primes

    # -------------------------------------------------------------------------
    # 3. Sum domain
    # -------------------------------------------------------------------------
    print()
    print("-" * 78)
    print("3. EVEN-SUM DOMAIN")
    print("-" * 78)

    s_min, s_max, sums = build_sum_domain()

    print(f"s minimum       = {s_min:,}")
    print(f"s maximum       = {s_max:,}")
    print(f"candidate sums  = {sums.size:,}")
    print("representation  = s = s_min + 2k")

    # -------------------------------------------------------------------------
    # 4. Modulus families
    # -------------------------------------------------------------------------
    print()
    print("-" * 78)
    print("4. MODULUS FAMILIES")
    print("-" * 78)

    print("cyclotomic =", CYCLOTOMIC)
    print("control    =", CONTROL)

    # -------------------------------------------------------------------------
    # 5. QR preparation
    # -------------------------------------------------------------------------
    print()
    print("-" * 78)
    print("5. QR TABLE PREPARATION")
    print("-" * 78)

    all_moduli = sorted(set(CYCLOTOMIC + CONTROL))
    qr_tables = prepare_qr_tables(all_moduli)

    cyclo_tables = {
        m: qr_tables[m]
        for m in CYCLOTOMIC
    }

    control_tables = {
        m: qr_tables[m]
        for m in CONTROL
    }

    # -------------------------------------------------------------------------
    # 6. Structural sanity checks
    # -------------------------------------------------------------------------
    print()
    print("-" * 78)
    print("6. TRUE-SUM SURVIVAL CHECK")
    print("-" * 78)

    for i, target in enumerate(targets, 1):
        n = target.n

        cyclo_ok = target.s in vectorized_filter(
            n,
            sums,
            CYCLOTOMIC,
            cyclo_tables,
        )

        control_ok = target.s in vectorized_filter(
            n,
            sums,
            CONTROL,
            control_tables,
        )

        print(
            f"target {i:2d}: "
            f"cyclotomic={cyclo_ok} "
            f"control={control_ok}"
        )

    # -------------------------------------------------------------------------
    # 7. Run targets
    # -------------------------------------------------------------------------
    results = []

    for i, target in enumerate(targets, 1):
        results.append(
            run_target(
                i,
                target,
                sums,
                cyclo_tables,
                control_tables,
            )
        )

    # -------------------------------------------------------------------------
    # 8. Aggregate summary
    # -------------------------------------------------------------------------
    print()
    print("=" * 78)
    print("8. GLOBAL SUMMARY")
    print("=" * 78)

    baseline_median_tests = statistics.median(
        r["baseline_tests"] for r in results
    )

    baseline_median_time = statistics.median(
        r["baseline_time"] for r in results
    )

    cyclo_median_candidates = statistics.median(
        r["cyclo_candidates"] for r in results
    )

    cyclo_median_exact_tests = statistics.median(
        r["cyclo_exact_tests"] for r in results
    )

    cyclo_median_total_time = statistics.median(
        r["cyclo_total_time"] for r in results
    )

    print(
        f"{'method':32s}"
        f"{'median candidates':>20s}"
        f"{'median exact':>16s}"
        f"{'median time':>16s}"
    )
    print("-" * 84)

    print(
        f"{'plain Python baseline':32s}"
        f"{'-':>20s}"
        f"{baseline_median_tests:16,.1f}"
        f"{baseline_median_time:16.6f}"
    )

    print(
        f"{'NumPy vectorized cyclotomic':32s}"
        f"{cyclo_median_candidates:20,.1f}"
        f"{cyclo_median_exact_tests:16,.1f}"
        f"{cyclo_median_total_time:16.6f}"
    )

    if RUN_PACKED:
        valid = [
            r["packed_total_time"]
            for r in results
            if r["packed_total_time"] is not None
        ]

        packed_median_time = statistics.median(valid)

        print(
            f"{'NumPy packed-bit cyclotomic':32s}"
            f"{cyclo_median_candidates:20,.1f}"
            f"{cyclo_median_exact_tests:16,.1f}"
            f"{packed_median_time:16.6f}"
        )

    # -------------------------------------------------------------------------
    # 9. Speed ratios
    # -------------------------------------------------------------------------
    print()
    print("-" * 78)
    print("9. SPEED DIAGNOSTIC")
    print("-" * 78)

    print(
        f"plain baseline median time       "
        f"= {baseline_median_time:.6f}s"
    )
    print(
        f"vectorized cyclotomic median     "
        f"= {cyclo_median_total_time:.6f}s"
    )

    if cyclo_median_total_time > 0:
        print(
            f"baseline / vectorized ratio      "
            f"= {baseline_median_time / cyclo_median_total_time:.3f}x"
        )

    if RUN_PACKED:
        valid = [
            r["packed_total_time"]
            for r in results
            if r["packed_total_time"] is not None
        ]
        packed_median_time = statistics.median(valid)

        print(
            f"packed-bit cyclotomic median     "
            f"= {packed_median_time:.6f}s"
        )

        if packed_median_time > 0:
            print(
                f"baseline / packed ratio          "
                f"= {baseline_median_time / packed_median_time:.3f}x"
            )

    # -------------------------------------------------------------------------
    # 10. Mathematical reduction
    # -------------------------------------------------------------------------
    print()
    print("-" * 78)
    print("10. MATHEMATICAL REDUCTION")
    print("-" * 78)

    print(
        f"median baseline exact tests      = "
        f"{baseline_median_tests:,.1f}"
    )
    print(
        f"median cyclotomic exact tests    = "
        f"{cyclo_median_exact_tests:,.1f}"
    )

    if cyclo_median_exact_tests > 0:
        print(
            f"exact-test reduction             = "
            f"{baseline_median_tests / cyclo_median_exact_tests:,.3f}x"
        )

    print()
    print("The mathematical filter is unchanged:")
    print("    d^2 = s^2 - 4n")
    print()
    print("Only the computational representation changes.")

    # -------------------------------------------------------------------------
    # 11. Decision criterion
    # -------------------------------------------------------------------------
    print()
    print("-" * 78)
    print("11. DECISION CRITERION")
    print("-" * 78)

    print(
        "SUCCESS CONDITION:"
    )
    print(
        "    bulk/vectorized filtering becomes faster than"
    )
    print(
        "    the plain Python exact-sum scan while preserving"
    )
    print(
        "    12/12 exact recovery."
    )

    print()
    print(
        "If vectorization is still slower:"
    )
    print(
        "    the next bottleneck is NumPy arithmetic/modulo cost."
    )

    print()
    print(
        "If vectorization becomes faster:"
    )
    print(
        "    the QR discriminant sieve becomes computationally"
    )
    print(
        "    interesting, and larger target batches should be tested."
    )

    print()
    print(
        "If packed-bit is faster than ordinary vectorization:"
    )
    print(
        "    memory bandwidth / mask representation is important."
    )

    print()
    print(
        "If cyclotomic and control have similar runtime:"
    )
    print(
        "    the advantage is generic quadratic-residue filtering,"
    )
    print(
        "    not yet evidence for special F(x)=x^2+x+1 structure."
    )

    print()
    print("=" * 78)
    print("EXPERIMENT 53 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()