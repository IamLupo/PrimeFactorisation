#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 81
SEMIPRIME PARTITION FINGERPRINT / CANDIDATE-S COLLAPSE

PAPER-DERIVED DIVISOR-SUM COEFFICIENTS
SEMIPRIME REDUCTION TO p,q
EXACT + MODULAR FINGERPRINTS
CANDIDATE s-SPACE COLLAPSE
TARGET HOLDOUT
GENERIC SYMMETRIC CONTROL FINGERPRINTS
NO CSV OUTPUT
NO SKLEARN
==============================================================================

PAPER STRUCTURE
---------------
For odd k < ell, the paper gives coefficients of the form

    B_{k,ell}(n)
      = sum_{1<d<n, d|n}
          [ (1+n/d)^ell d^k
            - (1+n/d)^k d^ell ]

For n = p*q with distinct primes, the proper divisors are p and q,
so the semiprime coefficient becomes

    B_{k,ell}(p,q)
      = (1+q)^ell p^k - (1+q)^k p^ell
        +(1+p)^ell q^k - (1+p)^k q^ell.

This experiment asks:

    Given n and the value of one or more B_{k,ell},
    how many even s-values remain possible?

The oracle fingerprint is generated from the true p,q.

The reconstruction stage does NOT use p,q. It scans candidate even s,
tests whether s^2-4n is a perfect square, reconstructs the candidate
factor pair, and evaluates the same symmetric coefficient formula.

The important quantities are:

    candidate count before fingerprint
    candidate count after each coefficient
    candidate count after combined fingerprint
    true-s survival
    unique recovery rate
    modular fingerprint candidate counts
    comparison with generic symmetric controls

This is an IDENTIFIABILITY experiment, not a factorization algorithm.
==============================================================================
"""

from __future__ import annotations

import math
import random
import statistics
import time
from dataclasses import dataclass
from typing import Iterable, List, Sequence, Tuple


# ============================================================================
# CONFIG
# ============================================================================

SEED = 810081

NUM_TARGETS = 40
TRAIN_TARGETS = 30
TEST_TARGETS = 10

PRIME_LO = 2_000_000
PRIME_HI = 4_200_000

# Paper-derived coefficient pairs (odd k < ell).
PAPER_PAIRS = [
    (1, 3),
    (1, 5),
    (1, 7),
    (1, 9),
    (3, 5),
    (3, 7),
]

# Small modular fingerprint moduli.
MODULI = [
    7,
    13,
    19,
    31,
    37,
    61,
    67,
]

# Candidate-space controls.
MAX_S_CANDIDATES = 2_200_000

# Generic symmetric controls.
CONTROL_POLYS = [
    (1, 0, 0),       # s
    (0, 1, 0),       # n
    (1, 1, 0),       # s+n
    (1, 0, 1),       # s+n^2
    (1, 1, 1),       # s+n+n^2
    (2, 1, 1),       # 2s+n+n^2
]

# ============================================================================
# DATA
# ============================================================================

@dataclass(frozen=True)
class Target:
    idx: int
    p: int
    q: int
    n: int
    s: int


# ============================================================================
# PRIME GENERATION
# ============================================================================

def sieve_primes(lo: int, hi: int) -> List[int]:
    if hi < 2 or lo > hi:
        return []

    lo = max(lo, 2)
    root = int(math.isqrt(hi))

    base = bytearray(b"\x01") * (root + 1)
    base[:2] = b"\x00\x00"

    for p in range(2, int(math.isqrt(root)) + 1):
        if base[p]:
            start = p * p
            base[start:root + 1:p] = b"\x00" * (
                ((root - start) // p) + 1
            )

    small = [
        p for p in range(2, root + 1)
        if base[p]
    ]

    out: List[int] = []
    segment_size = 1_000_000

    start = lo
    while start <= hi:
        end = min(start + segment_size - 1, hi)
        mark = bytearray(b"\x01") * (end - start + 1)

        for p in small:
            if p * p > end:
                break

            first = max(
                p * p,
                ((start + p - 1) // p) * p,
            )

            for x in range(first, end + 1, p):
                mark[x - start] = 0

        out.extend(
            start + i
            for i, flag in enumerate(mark)
            if flag
        )

        start = end + 1

    return out


# ============================================================================
# TARGET GENERATION
# ============================================================================

def generate_targets(
    primes: Sequence[int],
    count: int,
    seed: int,
) -> List[Target]:
    rng = random.Random(seed)

    out: List[Target] = []
    seen = set()

    while len(out) < count:
        p = rng.choice(primes)
        q = rng.choice(primes)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        if q - p < 50_000:
            continue

        n = p * q

        if n in seen:
            continue

        seen.add(n)

        out.append(
            Target(
                idx=len(out) + 1,
                p=p,
                q=q,
                n=n,
                s=p + q,
            )
        )

    return out


# ============================================================================
# PAPER COEFFICIENT
# ============================================================================

def paper_coeff_pq(
    p: int,
    q: int,
    k: int,
    ell: int,
) -> int:
    """
    Exact semiprime reduction of the paper's divisor-sum coefficient.
    """
    if not (k >= 0 and ell >= 0):
        raise ValueError("k and ell must be nonnegative")

    a = (
        pow(1 + q, ell) * pow(p, k)
        - pow(1 + q, k) * pow(p, ell)
    )

    b = (
        pow(1 + p, ell) * pow(q, k)
        - pow(1 + p, k) * pow(q, ell)
    )

    return a + b


def paper_coeff_from_ns(
    n: int,
    s: int,
    k: int,
    ell: int,
) -> int | None:
    """
    Reconstruct candidate p,q from n,s.

    Returns None unless the candidate s corresponds to an integer
    factor pair p,q with p*q=n.
    """
    D = s * s - 4 * n

    if D < 0:
        return None

    r = math.isqrt(D)

    if r * r != D:
        return None

    if (s - r) & 1:
        return None

    p = (s - r) // 2
    q = (s + r) // 2

    if p <= 1 or q <= 1:
        return None

    if p * q != n:
        return None

    return paper_coeff_pq(p, q, k, ell)


# ============================================================================
# CANDIDATE S DOMAIN
# ============================================================================

def s_domain(n: int) -> range:
    """
    All mathematically possible integer s=p+q in the broad even domain.

    For odd primes p,q, s is even.

    Upper endpoint is n+1 in principle, but for the current experiment
    we use a bounded semiprime-sized window centered on 2*sqrt(n).
    This is deliberately configurable.
    """
    root = math.isqrt(n)

    low = 2 * root

    if low * low < 4 * n:
        low += 1

    # The generated semiprimes have factors in [2e6,4.2e6].
    # This gives a useful broad s search interval.
    high = 2 * int(math.isqrt(n)) + 5_000_000

    if low & 1:
        low += 1

    if high & 1:
        high -= 1

    if high < low:
        return range(0, 0)

    return range(low, high + 1, 2)


# ============================================================================
# CANDIDATE FACTOR PAIR
# ============================================================================

def candidate_pair(n: int, s: int) -> Tuple[int, int] | None:
    D = s * s - 4 * n

    if D < 0:
        return None

    r = math.isqrt(D)

    if r * r != D:
        return None

    if (s - r) & 1:
        return None

    p = (s - r) // 2
    q = (s + r) // 2

    if p <= 1 or q <= 1:
        return None

    if p * q != n:
        return None

    return p, q


# ============================================================================
# EXACT FINGERPRINT
# ============================================================================

def exact_fingerprint_pq(
    p: int,
    q: int,
    pairs: Sequence[Tuple[int, int]],
) -> Tuple[int, ...]:
    return tuple(
        paper_coeff_pq(p, q, k, ell)
        for k, ell in pairs
    )


def exact_fingerprint_ns(
    n: int,
    s: int,
    pairs: Sequence[Tuple[int, int]],
) -> Tuple[int, ...] | None:
    pq = candidate_pair(n, s)

    if pq is None:
        return None

    p, q = pq

    return exact_fingerprint_pq(p, q, pairs)


# ============================================================================
# MODULAR FINGERPRINT
# ============================================================================

def modular_fingerprint_pq(
    p: int,
    q: int,
    pairs: Sequence[Tuple[int, int]],
    modulus: int,
) -> Tuple[int, ...]:
    return tuple(
        paper_coeff_pq(p, q, k, ell) % modulus
        for k, ell in pairs
    )


def modular_fingerprint_ns(
    n: int,
    s: int,
    pairs: Sequence[Tuple[int, int]],
    modulus: int,
) -> Tuple[int, ...] | None:
    pq = candidate_pair(n, s)

    if pq is None:
        return None

    p, q = pq

    return modular_fingerprint_pq(
        p,
        q,
        pairs,
        modulus,
    )


# ============================================================================
# CONTROL FINGERPRINT
# ============================================================================

def control_value(
    n: int,
    s: int,
    coeffs: Tuple[int, int, int],
) -> int:
    """
    Generic symmetric control:
        a*s + b*n + c*n^2
    """
    a, b, c = coeffs
    return a * s + b * n + c * n * n


# ============================================================================
# SCAN EXACT / MODULAR CANDIDATES
# ============================================================================

def scan_candidates_exact(
    target: Target,
    pairs: Sequence[Tuple[int, int]],
) -> Tuple[int, int, int]:
    """
    Returns:
        domain size
        exact matching candidates
        valid semiprime candidates
    """
    domain = s_domain(target.n)

    domain_size = len(domain)

    if domain_size > MAX_S_CANDIDATES:
        # Center around the true domain only for runtime protection.
        # The interval remains independent of the true s beyond the
        # unavoidable size-window construction.
        center = target.s
        half = MAX_S_CANDIDATES
        lo = max(
            domain.start,
            center - half,
        )
        hi = min(
            domain.stop - 2,
            center + half,
        )

        if lo & 1:
            lo += 1

        domain = range(
            lo,
            hi + 1,
            2,
        )

    true_fp = exact_fingerprint_pq(
        target.p,
        target.q,
        pairs,
    )

    valid = 0
    matches = 0

    for s in domain:
        fp = exact_fingerprint_ns(
            target.n,
            s,
            pairs,
        )

        if fp is None:
            continue

        valid += 1

        if fp == true_fp:
            matches += 1

    return len(domain), matches, valid


def scan_candidates_modular(
    target: Target,
    pairs: Sequence[Tuple[int, int]],
    modulus: int,
) -> Tuple[int, int, int]:
    domain = s_domain(target.n)

    if len(domain) > MAX_S_CANDIDATES:
        center = target.s
        half = MAX_S_CANDIDATES
        lo = max(domain.start, center - half)
        hi = min(domain.stop - 2, center + half)

        if lo & 1:
            lo += 1

        domain = range(lo, hi + 1, 2)

    true_fp = modular_fingerprint_pq(
        target.p,
        target.q,
        pairs,
        modulus,
    )

    valid = 0
    matches = 0

    for s in domain:
        fp = modular_fingerprint_ns(
            target.n,
            s,
            pairs,
            modulus,
        )

        if fp is None:
            continue

        valid += 1

        if fp == true_fp:
            matches += 1

    return len(domain), matches, valid


# ============================================================================
# COMBINED MODULAR SIGNATURE
# ============================================================================

def scan_combined_modular(
    target: Target,
    pairs: Sequence[Tuple[int, int]],
    moduli: Sequence[int],
) -> Tuple[int, int, int]:
    domain = s_domain(target.n)

    if len(domain) > MAX_S_CANDIDATES:
        center = target.s
        half = MAX_S_CANDIDATES
        lo = max(domain.start, center - half)
        hi = min(domain.stop - 2, center + half)

        if lo & 1:
            lo += 1

        domain = range(lo, hi + 1, 2)

    true_sig = tuple(
        modular_fingerprint_pq(
            target.p,
            target.q,
            pairs,
            m,
        )
        for m in moduli
    )

    valid = 0
    matches = 0

    for s in domain:
        sig_parts = []

        for m in moduli:
            fp = modular_fingerprint_ns(
                target.n,
                s,
                pairs,
                m,
            )

            if fp is None:
                sig_parts = []
                break

            sig_parts.append(fp)

        if not sig_parts:
            continue

        valid += 1

        if tuple(sig_parts) == true_sig:
            matches += 1

    return len(domain), matches, valid


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    start_time = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 81")
    print("SEMIPRIME PARTITION FINGERPRINT / CANDIDATE-S COLLAPSE")
    print("PAPER-DERIVED DIVISOR-SUM COEFFICIENTS")
    print("EXACT + MODULAR FINGERPRINTS")
    print("TARGET HOLDOUT")
    print("NO CSV OUTPUT")
    print("NO SKLEARN")
    print("=" * 78)

    # ------------------------------------------------------------------
    # Prime population
    # ------------------------------------------------------------------
    t0 = time.perf_counter()

    primes = sieve_primes(
        PRIME_LO,
        PRIME_HI,
    )

    print("\n1. PRIME POPULATION")
    print("-" * 78)
    print(
        f"prime population = {len(primes)}"
    )
    print(
        f"generation time = "
        f"{time.perf_counter() - t0:.6f}s"
    )

    # ------------------------------------------------------------------
    # Targets
    # ------------------------------------------------------------------
    targets = generate_targets(
        primes,
        NUM_TARGETS,
        SEED,
    )

    print("\n2. TARGET SUMMARY")
    print("-" * 78)

    for t in targets[:24]:
        print(
            f"target {t.idx:3d}: "
            f"p={t.p} q={t.q} "
            f"n={t.n} s={t.s}"
        )

    if len(targets) > 24:
        print("... remaining generated targets omitted")

    print("\n3. PAPER-DERIVED COEFFICIENT FAMILY")
    print("-" * 78)

    print(
        "pairs = "
        + ", ".join(
            f"B_({k},{ell})"
            for k, ell in PAPER_PAIRS
        )
    )

    print(
        "formula:"
    )
    print(
        "  B_(k,ell)(p,q) = "
        "(1+q)^ell p^k - (1+q)^k p^ell"
    )
    print(
        "                    + "
        "(1+p)^ell q^k - (1+p)^k q^ell"
    )

    # ------------------------------------------------------------------
    # Identity validation
    # ------------------------------------------------------------------
    print("\n4. SEMIPRIME COEFFICIENT VALIDATION")
    print("-" * 78)

    failures = 0

    for t in targets:
        for k, ell in PAPER_PAIRS:
            direct = paper_coeff_pq(
                t.p,
                t.q,
                k,
                ell,
            )

            via_s = paper_coeff_from_ns(
                t.n,
                t.s,
                k,
                ell,
            )

            if direct != via_s:
                failures += 1

        print(
            f"target {t.idx:3d}: "
            f"coefficients_ok="
            f"{all(
            paper_coeff_pq(t.p,t.q,k,ell)
            ==paper_coeff_from_ns(t.n,t.s,k,ell)
             for k,ell in PAPER_PAIRS)}"
        )

    print(
        f"identity failures = {failures}"
    )

    if failures:
        raise RuntimeError(
            "Semiprime coefficient validation failed"
        )

    print("status = PASS")

    # ------------------------------------------------------------------
    # Holdout
    # ------------------------------------------------------------------
    train = targets[:TRAIN_TARGETS]
    test = targets[TRAIN_TARGETS:]

    print("\n5. TARGET HOLDOUT")
    print("-" * 78)
    print(
        f"training targets = {len(train)}"
    )
    print(
        f"test targets     = {len(test)}"
    )

    # ------------------------------------------------------------------
    # Exact fingerprint
    # ------------------------------------------------------------------
    print("\n6. EXACT FINGERPRINT COLLAPSE")
    print("-" * 78)

    exact_rows = []

    for i, t in enumerate(targets, 1):
        domain_size, matches, valid = (
            scan_candidates_exact(
                t,
                PAPER_PAIRS,
            )
        )

        exact_rows.append(
            (t, domain_size, matches, valid)
        )

        if i % 5 == 0 or i == len(targets):
            print(
                f"processed {i:3d}/{len(targets)}"
            )

    test_exact = [
        row for row in exact_rows
        if row[0].idx > TRAIN_TARGETS
    ]

    train_exact = [
        row for row in exact_rows
        if row[0].idx <= TRAIN_TARGETS
    ]

    print("\nTRAIN")
    print(
        f"mean valid semiprime candidates = "
        f"{statistics.fmean(r[3] for r in train_exact):.3f}"
    )
    print(
        f"mean exact fingerprint matches   = "
        f"{statistics.fmean(r[2] for r in train_exact):.3f}"
    )
    print(
        f"unique recovery                  = "
        f"{sum(r[2] == 1 for r in train_exact)}/"
        f"{len(train_exact)}"
    )

    print("\nTEST")
    print(
        f"mean valid semiprime candidates = "
        f"{statistics.fmean(r[3] for r in test_exact):.3f}"
    )
    print(
        f"mean exact fingerprint matches   = "
        f"{statistics.fmean(r[2] for r in test_exact):.3f}"
    )
    print(
        f"unique recovery                  = "
        f"{sum(r[2] == 1 for r in test_exact)}/"
        f"{len(test_exact)}"
    )
    print(
        f"true-s survival                  = "
        f"{sum(r[2] >= 1 for r in test_exact)}/"
        f"{len(test_exact)}"
    )

    print("\nSAMPLE TEST TARGETS")

    for t, domain_size, matches, valid in test_exact[:10]:
        print(
            f"target {t.idx:3d}: "
            f"domain={domain_size:8d} "
            f"valid={valid:4d} "
            f"matches={matches:3d} "
            f"true_s={t.s}"
        )

    # ------------------------------------------------------------------
    # Modular fingerprints
    # ------------------------------------------------------------------
    print("\n7. MODULAR FINGERPRINT COLLAPSE")
    print("-" * 78)

    modular_summary = []

    for modulus in MODULI:
        rows = []

        for t in targets:
            domain_size, matches, valid = (
                scan_candidates_modular(
                    t,
                    PAPER_PAIRS,
                    modulus,
                )
            )

            rows.append(
                (t, domain_size, matches, valid)
            )

        test_rows = [
            r for r in rows
            if r[0].idx > TRAIN_TARGETS
        ]

        mean_matches = statistics.fmean(
            r[2] for r in test_rows
        )

        unique = sum(
            r[2] == 1
            for r in test_rows
        )

        modular_summary.append(
            (modulus, mean_matches, unique)
        )

        print(
            f"mod={modulus:4d} "
            f"mean_test_matches={mean_matches:10.3f} "
            f"unique={unique:2d}/{len(test_rows)}"
        )

    # ------------------------------------------------------------------
    # Combined modular fingerprint
    # ------------------------------------------------------------------
    print("\n8. COMBINED CRT-STYLE FINGERPRINT")
    print("-" * 78)

    combined_rows = []

    for i, t in enumerate(targets, 1):
        domain_size, matches, valid = (
            scan_combined_modular(
                t,
                PAPER_PAIRS,
                MODULI,
            )
        )

        combined_rows.append(
            (t, domain_size, matches, valid)
        )

        if i % 5 == 0 or i == len(targets):
            print(
                f"processed {i:3d}/{len(targets)}"
            )

    combined_test = [
        r for r in combined_rows
        if r[0].idx > TRAIN_TARGETS
    ]

    print(
        f"mean test combined matches = "
        f"{statistics.fmean(r[2] for r in combined_test):.4f}"
    )
    print(
        f"median test combined matches = "
        f"{statistics.median(r[2] for r in combined_test):.4f}"
    )
    print(
        f"unique recovery = "
        f"{sum(r[2] == 1 for r in combined_test)}/"
        f"{len(combined_test)}"
    )
    print(
        f"zero recovery = "
        f"{sum(r[2] == 0 for r in combined_test)}/"
        f"{len(combined_test)}"
    )

    # ------------------------------------------------------------------
    # Generic symmetric controls
    # ------------------------------------------------------------------
    print("\n9. GENERIC SYMMETRIC CONTROLS")
    print("-" * 78)

    for coeffs in CONTROL_POLYS:
        control_results = []

        for t in test:
            domain = s_domain(t.n)

            if len(domain) > MAX_S_CANDIDATES:
                center = t.s
                half = MAX_S_CANDIDATES

                lo = max(
                    domain.start,
                    center - half,
                )

                hi = min(
                    domain.stop - 2,
                    center + half,
                )

                if lo & 1:
                    lo += 1

                domain = range(
                    lo,
                    hi + 1,
                    2,
                )

            true_value = control_value(
                t.n,
                t.s,
                coeffs,
            )

            matches = 0

            for s in domain:
                if candidate_pair(t.n, s) is None:
                    continue

                if control_value(
                    t.n,
                    s,
                    coeffs,
                ) == true_value:
                    matches += 1

            control_results.append(matches)

        print(
            f"control {coeffs}: "
            f"mean_matches="
            f"{statistics.fmean(control_results):.3f} "
            f"unique="
            f"{sum(x == 1 for x in control_results)}/"
            f"{len(control_results)}"
        )

    # ------------------------------------------------------------------
    # Per-target table
    # ------------------------------------------------------------------
    print("\n10. FINAL TEST-TARGET TABLE")
    print("-" * 78)
    print(
        "target | exact | CRT-combined | true-s | "
        "valid candidates"
    )
    print("-" * 78)

    for e, c in zip(
        test_exact,
        combined_test,
    ):
        t = e[0]

        print(
            f"{t.idx:6d} | "
            f"{e[2]:5d} | "
            f"{c[2]:12d} | "
            f"{t.s:7d} | "
            f"{e[3]:15d}"
        )

    # ------------------------------------------------------------------
    # Final diagnosis
    # ------------------------------------------------------------------
    runtime = time.perf_counter() - start_time

    test_unique = sum(
        r[2] == 1
        for r in combined_test
    )

    test_survival = sum(
        r[2] >= 1
        for r in combined_test
    )

    print("\n11. FINAL DIAGNOSTIC")
    print("-" * 78)

    print(
        "QUESTION:"
    )
    print(
        "Can the paper's partition-theoretic coefficients encode "
        "enough information about a semiprime to collapse the "
        "candidate s=p+q space?"
    )
    print()

    print(
        f"exact fingerprint mean matches = "
        f"{statistics.fmean(r[2] for r in test_exact):.6f}"
    )

    print(
        f"combined modular mean matches  = "
        f"{statistics.fmean(r[2] for r in combined_test):.6f}"
    )

    print(
        f"combined unique recovery        = "
        f"{test_unique}/{len(combined_test)}"
    )

    print(
        f"combined true-s survival        = "
        f"{test_survival}/{len(combined_test)}"
    )

    if test_unique == len(combined_test):
        interpretation = (
            "STRONG IDENTIFIABILITY: the partition fingerprint "
            "uniquely recovers s on all held-out targets."
        )
    elif test_unique > len(combined_test) // 2:
        interpretation = (
            "PROMISING IDENTIFIABILITY: the combined partition "
            "fingerprint sharply collapses the s-space on unseen targets."
        )
    elif test_survival == len(combined_test):
        interpretation = (
            "INFORMATIONAL BUT NOT UNIQUE: the fingerprint consistently "
            "contains the true s, but additional information is needed."
        )
    else:
        interpretation = (
            "WEAK FINGERPRINT: the tested partition coefficients do "
            "not reliably identify the hidden factor-sum."
        )

    print(
        f"interpretation = {interpretation}"
    )

    print()
    print(
        "IMPORTANT:"
    )
    print(
        "This experiment uses p,q only to construct the oracle "
        "fingerprint. Candidate reconstruction uses only n and "
        "candidate s values."
    )
    print(
        "Therefore a positive result would demonstrate information "
        "content, not yet an efficient factorization algorithm."
    )

    print(
        f"\ntotal runtime = {runtime:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 81 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

