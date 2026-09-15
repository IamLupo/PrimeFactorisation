#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 95
PAPER PRIME-DETECTOR RESOLVENT SPECTRUM
f_(k,l) FAMILY: (1,3), (1,5), (1,7), (3,5), (3,7), (5,7)

SEMIPRIME SPECIALIZATION:
    n = p*q
    s = p+q

PAPER-DERIVED COEFFICIENT:
    B_(k,l)(p,q)
      = (1+q)^l p^k - (1+q)^k p^l
        + (1+p)^l q^k - (1+p)^k q^l

The candidate-side computation NEVER uses p or q.

Instead B_(k,l) is reduced to a symmetric polynomial in (n,s).

QUESTIONS:
  1. Which paper detector gives the strongest s-collapse?
  2. Does f_(1,3) outperform E1/H6?
  3. Do multiple paper detectors add information?
  4. Does the effect survive unseen targets?
  5. Do cyclotomic moduli outperform matched controls?

NO FACTOR-PAIR CHECK DURING SIEVE
STRICT TARGET HOLDOUT
NO CSV
NO SKLEARN
==============================================================================
"""

from __future__ import annotations

import math
import statistics
import time
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np


# ============================================================================
# PARAMETERS
# ============================================================================

NUM_TARGETS = 40
TRAIN_TARGETS = 30
TEST_TARGETS = 10

# Same scale as earlier experiments.
PRIME_MIN = 2_000_000
PRIME_MAX = 4_200_000

# Every even s in this range is tested.
S_MIN = 4_000_000
S_MAX = 8_400_000

# Paper-derived detector family.
DETECTORS: Tuple[Tuple[int, int], ...] = (
    (1, 3),
    (1, 5),
    (1, 7),
    (3, 5),
    (3, 7),
    (5, 7),
)

# Cyclotomic-style moduli used in previous experiments.
C3 = (7, 13, 19)
C5 = (7, 13, 19, 31, 37)
C7 = (7, 13, 19, 31, 37, 61, 67)

# Matched controls.
R3 = (673, 4561, 4759)
R5 = (673, 4561, 4759, 6211, 7879)
R7 = (673, 4561, 4759, 6211, 7879, 7951, 8689)


# ============================================================================
# TARGET
# ============================================================================

@dataclass(frozen=True)
class Target:
    p: int
    q: int
    n: int
    s: int


# ============================================================================
# PRIME GENERATION
# ============================================================================

def sieve_primes(lo: int, hi: int) -> List[int]:
    """
    Return primes in [lo, hi].

    A standard bytearray sieve is sufficient for the prime population.
    """
    if hi < 2:
        return []

    flags = bytearray(b"\x01") * (hi + 1)
    flags[0:2] = b"\x00\x00"

    limit = int(math.isqrt(hi))

    for p in range(2, limit + 1):
        if flags[p]:
            start = p * p
            flags[start : hi + 1 : p] = b"\x00" * (
                ((hi - start) // p) + 1
            )

    return [p for p in range(max(2, lo), hi + 1) if flags[p]]


def generate_targets(primes: Sequence[int], count: int) -> List[Target]:
    """
    Deterministic target generation.

    We avoid repeated p,q pairs and keep p < q.
    """
    targets: List[Target] = []
    seen = set()

    # Deterministic pseudo-random-looking stride pattern.
    plen = len(primes)

    i = 17
    j = 113

    while len(targets) < count:
        p = primes[i % plen]
        q = primes[j % plen]

        if p > q:
            p, q = q, p

        if p != q and (p, q) not in seen:
            seen.add((p, q))
            targets.append(Target(p, q, p * q, p + q))

        i += 137
        j += 593

    return targets


# ============================================================================
# DIRECT PAPER COEFFICIENT
# ============================================================================

def B_direct(p: int, q: int, k: int, ell: int) -> int:
    """
    Oracle construction using p,q.

    This is used ONLY for validation/oracle target generation.
    """
    return (
        (1 + q) ** ell * p**k
        - (1 + q) ** k * p**ell
        + (1 + p) ** ell * q**k
        - (1 + p) ** k * q**ell
    )


# ============================================================================
# SYMMETRIC POLYNOMIAL CONSTRUCTION
# ============================================================================

def build_B_symmetric_polynomial(
    k: int,
    ell: int,
) -> Dict[Tuple[int, int], int]:
    """
    Build the symmetric monomial representation

        B(p,q) = sum c[a,b] * (p^a q^b + p^b q^a)

    for a>b, and

        c[a,a] * p^a q^a

    for a=b.

    We never use p,q in candidate evaluation.
    """

    terms: Dict[Tuple[int, int], int] = {}

    def add(a: int, b: int, coeff: int) -> None:
        terms[(a, b)] = terms.get((a, b), 0) + coeff

    # +(1+q)^ell p^k
    for b in range(ell + 1):
        add(k, b, math.comb(ell, b))

    # -(1+q)^k p^ell
    for b in range(k + 1):
        add(ell, b, -math.comb(k, b))

    # +(1+p)^ell q^k
    for a in range(ell + 1):
        add(a, k, math.comb(ell, a))

    # -(1+p)^k q^ell
    for a in range(k + 1):
        add(a, ell, -math.comb(k, a))

    # Group swapped monomials.
    grouped: Dict[Tuple[int, int], int] = {}

    processed = set()

    for (a, b), coeff in list(terms.items()):
        if (a, b) in processed:
            continue

        if a == b:
            c = coeff
            grouped[(a, b)] = grouped.get((a, b), 0) + c
            processed.add((a, b))
            continue

        ca = terms.get((a, b), 0)
        cb = terms.get((b, a), 0)

        if ca != cb:
            raise ArithmeticError(
                f"B({k},{ell}) failed symmetry at ({a},{b}): "
                f"{ca} != {cb}"
            )

        u, v = max(a, b), min(a, b)
        grouped[(u, v)] = grouped.get((u, v), 0) + ca

        processed.add((a, b))
        processed.add((b, a))

    return {
        key: value
        for key, value in grouped.items()
        if value != 0
    }


def power_sum_polynomial_value(
    r: int,
    n: int,
    s: int,
) -> int:
    """
    R_r = p^r + q^r for p+q=s, pq=n.

    R_0 = 2
    R_1 = s
    R_r = s R_(r-1) - n R_(r-2)
    """
    if r == 0:
        return 2
    if r == 1:
        return s

    r0 = 2
    r1 = s

    for _ in range(2, r + 1):
        r0, r1 = r1, s * r1 - n * r0

    return r1


def B_from_ns(
    n: int,
    s: int,
    polynomial: Dict[Tuple[int, int], int],
) -> int:
    """
    Evaluate the symmetric B_(k,l) polynomial using only n and s.
    No factor-pair test is performed.
    """
    total = 0

    for (a, b), coeff in polynomial.items():
        if a == b:
            total += coeff * (n**a)
        else:
            lo = b
            diff = a - b
            total += coeff * (n**lo) * power_sum_polynomial_value(
                diff, n, s
            )

    return total


# ============================================================================
# MODULAR VERSION
# ============================================================================

def power_sum_mod(
    r: int,
    nmod: int,
    smod: int,
    mod: int,
) -> int:
    if r == 0:
        return 2 % mod

    if r == 1:
        return smod % mod

    r0 = 2 % mod
    r1 = smod % mod

    for _ in range(2, r + 1):
        r0, r1 = r1, (smod * r1 - nmod * r0) % mod

    return r1


def B_from_ns_mod(
    n: int,
    s: int,
    mod: int,
    polynomial: Dict[Tuple[int, int], int],
) -> int:
    nmod = n % mod
    smod = s % mod

    total = 0

    # Cache R_j for this (n,s,mod).
    max_diff = 0
    for (a, b) in polynomial:
        if a != b:
            max_diff = max(max_diff, a - b)

    R = [0] * (max_diff + 1)

    if max_diff >= 0:
        R[0] = 2 % mod
    if max_diff >= 1:
        R[1] = smod

    for r in range(2, max_diff + 1):
        R[r] = (smod * R[r - 1] - nmod * R[r - 2]) % mod

    for (a, b), coeff in polynomial.items():
        c = coeff % mod

        if a == b:
            total = (total + c * pow(nmod, a, mod)) % mod
        else:
            total = (
                total
                + c
                * pow(nmod, b, mod)
                * R[a - b]
            ) % mod

    return total


# ============================================================================
# MODULAR LOOKUP TABLES
# ============================================================================

def allowed_residue_table(
    n: int,
    mod: int,
    polynomial: Dict[Tuple[int, int], int],
    target_residue: int,
) -> np.ndarray:
    """
    allowed[r] = True iff B(n,r) == target B residue mod mod.
    """
    allowed = np.zeros(mod, dtype=np.bool_)

    for r in range(mod):
        value = B_from_ns_mod(
            n=n,
            s=r,
            mod=mod,
            polynomial=polynomial,
        )

        if value == target_residue:
            allowed[r] = True

    return allowed


# ============================================================================
# SIEVE
# ============================================================================

S_VALUES = np.arange(
    S_MIN,
    S_MAX + 1,
    2,
    dtype=np.int64,
)


def modular_survivors(
    target: Target,
    polynomial: Dict[Tuple[int, int], int],
    moduli: Sequence[int],
) -> np.ndarray:
    """
    Return even s candidates surviving all modular detector congruences.

    Candidate side knows only n and s.
    """
    mask = np.ones(len(S_VALUES), dtype=np.bool_)

    for mod in moduli:
        oracle = B_direct(
            target.p,
            target.q,
            infer_pair_from_poly(polynomial)[0],
            infer_pair_from_poly(polynomial)[1],
        ) % mod

        allowed = allowed_residue_table(
            target.n,
            mod,
            polynomial,
            oracle,
        )

        residues = np.mod(S_VALUES, mod)
        mask &= allowed[residues]

        if not mask.any():
            break

    return S_VALUES[mask]


def infer_pair_from_poly(
    polynomial: Dict[Tuple[int, int], int],
) -> Tuple[int, int]:
    """
    Identify the paper pair from the constructed polynomial.

    Only used for oracle validation, never candidate reconstruction.
    """
    # The six target families have distinct maximum degrees.
    # We can instead compare against known generated families.
    for pair in DETECTORS:
        candidate = build_B_symmetric_polynomial(*pair)
        if candidate == polynomial:
            return pair

    raise ValueError("Unknown detector polynomial")


# ============================================================================
# EXACT CHECK
# ============================================================================

def exact_matches(
    target: Target,
    detector: Tuple[int, int],
    polynomial: Dict[Tuple[int, int], int],
    candidates: np.ndarray,
) -> List[int]:
    oracle = B_direct(
        target.p,
        target.q,
        detector[0],
        detector[1],
    )

    matches: List[int] = []

    for s_val in candidates.tolist():
        value = B_from_ns(
            target.n,
            int(s_val),
            polynomial,
        )

        if value == oracle:
            matches.append(int(s_val))

    return matches


# ============================================================================
# VALIDATION
# ============================================================================

def validate_detector(
    targets: Sequence[Target],
    detector: Tuple[int, int],
    polynomial: Dict[Tuple[int, int], int],
) -> None:
    failures = 0

    for t in targets:
        direct = B_direct(
            t.p,
            t.q,
            detector[0],
            detector[1],
        )

        reconstructed = B_from_ns(
            t.n,
            t.s,
            polynomial,
        )

        if direct != reconstructed:
            failures += 1

    print(
        f"  f_{detector} identity failures = {failures}"
    )

    if failures:
        raise ArithmeticError(
            f"Detector {detector} failed symmetric reconstruction."
        )


# ============================================================================
# ONE DETECTOR FAMILY
# ============================================================================

@dataclass
class DetectorResult:
    detector: Tuple[int, int]
    mean_survivors: float
    median_survivors: float
    unique: int
    true_survival: int
    exact_target_matches: List[int]


def evaluate_detector_family(
    targets: Sequence[Target],
    detector: Tuple[int, int],
    polynomial: Dict[Tuple[int, int], int],
    moduli: Sequence[int],
    label: str,
) -> DetectorResult:

    counts: List[int] = []
    exact_counts: List[int] = []
    exact_solutions: List[int] = []

    unique = 0
    true_survival = 0

    total = len(targets)

    print(
        f"    detector f_{detector}, moduli={list(moduli)}"
    )

    for idx, target in enumerate(targets, 1):
        candidates = modular_survivors(
            target=target,
            polynomial=polynomial,
            moduli=moduli,
        )

        matches = exact_matches(
            target=target,
            detector=detector,
            polynomial=polynomial,
            candidates=candidates,
        )

        counts.append(len(candidates))
        exact_counts.append(len(matches))

        if len(matches) == 1:
            unique += 1

        if target.s in candidates:
            true_survival += 1

        exact_solutions.extend(matches)

        if idx in {1, 5, 10, total}:
            print(
                f"      processed {idx:2d}/{total}: "
                f"modular={len(candidates):5d} "
                f"exact={len(matches):3d}"
            )

    mean_survivors = statistics.mean(counts)
    median_survivors = statistics.median(counts)

    print(
        f"      mean modular survivors = {mean_survivors:.3f}"
    )
    print(
        f"      median modular survivors = {median_survivors:.3f}"
    )
    print(
        f"      unique exact recovery = {unique}/{total}"
    )
    print(
        f"      true modular survival = {true_survival}/{total}"
    )

    return DetectorResult(
        detector=detector,
        mean_survivors=mean_survivors,
        median_survivors=median_survivors,
        unique=unique,
        true_survival=true_survival,
        exact_target_matches=exact_counts,
    )


# ============================================================================
# INTERSECTION OF DETECTORS
# ============================================================================

def evaluate_detector_intersection(
    targets: Sequence[Target],
    detector_polynomials: Dict[
        Tuple[int, int],
        Dict[Tuple[int, int], int],
    ],
    detectors: Sequence[Tuple[int, int]],
    moduli: Sequence[int],
    label: str,
) -> None:

    print(
        f"\n    INTERSECTION: {[f'f{d}' for d in detectors]}"
    )

    survivor_counts: List[int] = []
    unique = 0
    true_survival = 0

    for idx, target in enumerate(targets, 1):
        mask = np.ones(len(S_VALUES), dtype=np.bool_)

        for detector in detectors:
            polynomial = detector_polynomials[detector]

            oracle = B_direct(
                target.p,
                target.q,
                detector[0],
                detector[1],
            )

            for mod in moduli:
                oracle_mod = oracle % mod

                allowed = allowed_residue_table(
                    target.n,
                    mod,
                    polynomial,
                    oracle_mod,
                )

                mask &= allowed[np.mod(S_VALUES, mod)]

                if not mask.any():
                    break

        candidates = S_VALUES[mask]

        survivor_counts.append(len(candidates))

        # Exact intersection check.
        exact: List[int] = []

        for s_val in candidates.tolist():
            ok = True

            for detector in detectors:
                polynomial = detector_polynomials[detector]

                oracle = B_direct(
                    target.p,
                    target.q,
                    detector[0],
                    detector[1],
                )

                candidate_value = B_from_ns(
                    target.n,
                    int(s_val),
                    polynomial,
                )

                if candidate_value != oracle:
                    ok = False
                    break

            if ok:
                exact.append(int(s_val))

        if len(exact) == 1:
            unique += 1

        if target.s in candidates:
            true_survival += 1

        if idx in {1, 5, 10, len(targets)}:
            print(
                f"      processed {idx:2d}/{len(targets)}: "
                f"modular={len(candidates):5d} "
                f"exact={len(exact):3d}"
            )

    print(
        f"      mean survivors = "
        f"{statistics.mean(survivor_counts):.3f}"
    )
    print(
        f"      median survivors = "
        f"{statistics.median(survivor_counts):.3f}"
    )
    print(
        f"      unique exact recovery = "
        f"{unique}/{len(targets)}"
    )
    print(
        f"      true modular survival = "
        f"{true_survival}/{len(targets)}"
    )


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:

    t0 = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 95")
    print("PAPER PRIME-DETECTOR RESOLVENT SPECTRUM")
    print("f_(k,l) SEMIPRIME SPECIALIZATION")
    print("LOW-WEIGHT MACMAHON / EISENSTEIN FAMILY")
    print("STRICT TARGET HOLDOUT")
    print("CYCLOTOMIC VS CONTROL")
    print("NO FACTOR-PAIR CHECK DURING SIEVE")
    print("NO CSV")
    print("NO SKLEARN")
    print("=" * 78)

    print("\n1. PARAMETERS")
    print("-" * 78)
    print(f"targets          = {NUM_TARGETS}")
    print(f"training targets = {TRAIN_TARGETS}")
    print(f"test targets     = {TEST_TARGETS}")
    print(f"s domain         = [{S_MIN}, {S_MAX}]")
    print(f"even candidates   = {len(S_VALUES)}")
    print(f"detectors        = {DETECTORS}")
    print(f"C3               = {C3}")
    print(f"C5               = {C5}")
    print(f"C7               = {C7}")
    print(f"R3               = {R3}")
    print(f"R5               = {R5}")
    print(f"R7               = {R7}")

    # ------------------------------------------------------------------------
    # PRIME POPULATION
    # ------------------------------------------------------------------------

    print("\n2. PRIME POPULATION")
    print("-" * 78)

    start = time.perf_counter()

    primes = sieve_primes(
        PRIME_MIN,
        PRIME_MAX,
    )

    print(f"prime population = {len(primes)}")
    print(
        f"generation time = "
        f"{time.perf_counter() - start:.6f}s"
    )

    targets = generate_targets(
        primes,
        NUM_TARGETS,
    )

    for i, target in enumerate(targets[:24], 1):
        print(
            f"target {i:3d}: "
            f"p={target.p} "
            f"q={target.q} "
            f"n={target.n} "
            f"s={target.s}"
        )

    if NUM_TARGETS > 24:
        print("... remaining generated targets omitted")

    train = targets[:TRAIN_TARGETS]
    test = targets[TRAIN_TARGETS:]

    print("\n3. TARGET HOLDOUT")
    print("-" * 78)
    print(f"training targets = {len(train)}")
    print(f"test targets     = {len(test)}")

    # ------------------------------------------------------------------------
    # BUILD POLYNOMIALS
    # ------------------------------------------------------------------------

    detector_polynomials = {}

    print("\n4. SEMIPRIME SYMMETRIC RESOLVENT VALIDATION")
    print("-" * 78)

    for detector in DETECTORS:
        print(
            f"\nf_{detector}"
        )

        poly = build_B_symmetric_polynomial(
            detector[0],
            detector[1],
        )

        detector_polynomials[detector] = poly

        validate_detector(
            targets,
            detector,
            poly,
        )

        print(
            f"  symmetric monomials = {len(poly)}"
        )

    # ------------------------------------------------------------------------
    # FACTORIZATION / DEGREE DIAGNOSTIC
    # ------------------------------------------------------------------------

    print("\n5. LOW-WEIGHT STRUCTURE")
    print("-" * 78)

    for detector in DETECTORS:

        poly = detector_polynomials[detector]

        max_degree = max(
            a + b
            for (a, b) in poly
        )

        max_s_degree = max(
            a
            for (a, b) in poly
        )

        print(
            f"f_{detector}: "
            f"max total degree={max_degree:2d}, "
            f"max p-power={max_s_degree:2d}, "
            f"terms={len(poly):2d}"
        )

    # ------------------------------------------------------------------------
    # SINGLE DETECTOR SPECTRUM
    # ------------------------------------------------------------------------

    configs = [
        ("C3", C3),
        ("C5", C5),
        ("C7", C7),
        ("R3", R3),
        ("R5", R5),
        ("R7", R7),
    ]

    all_results = {}

    print("\n6. SINGLE-DETECTOR SPECTRUM")
    print("=" * 78)

    for label, moduli in configs:

        print(
            f"\n{'-' * 78}"
        )
        print(
            f"{label}: moduli={list(moduli)}"
        )
        print(
            f"{'-' * 78}"
        )

        train_results = {}
        test_results = {}

        # Only print every detector's detailed timings on test.
        for detector in DETECTORS:

            print(
                f"\n  TRAIN {label}"
            )

            result_train = evaluate_detector_family(
                train,
                detector,
                detector_polynomials[detector],
                moduli,
                label,
            )

            print(
                f"\n  TEST {label}"
            )

            result_test = evaluate_detector_family(
                test,
                detector,
                detector_polynomials[detector],
                moduli,
                label,
            )

            train_results[detector] = result_train
            test_results[detector] = result_test

        all_results[label] = (
            train_results,
            test_results,
        )

    # ------------------------------------------------------------------------
    # SUMMARY
    # ------------------------------------------------------------------------

    print("\n7. CROSS-DETECTOR SUMMARY")
    print("-" * 78)

    print(
        "label detector        "
        "train_mean  test_mean  "
        "train_unique test_unique"
    )

    for label, (train_results, test_results) in all_results.items():

        for detector in DETECTORS:

            tr = train_results[detector]
            te = test_results[detector]

            print(
                f"{label:5s} "
                f"f{detector!s:10s} "
                f"{tr.mean_survivors:10.3f} "
                f"{te.mean_survivors:10.3f} "
                f"{tr.unique:4d}/{len(train):<4d} "
                f"{te.unique:4d}/{len(test):<4d}"
            )

    # ------------------------------------------------------------------------
    # INTERSECTIONS
    # ------------------------------------------------------------------------

    print("\n8. DETECTOR INTERSECTION TEST")
    print("-" * 78)

    # The most important question:
    # do different paper detectors add information beyond f_(1,3)?
    #
    # We test on the strongest cyclotomic configuration first.

    evaluate_detector_intersection(
        test,
        detector_polynomials,
        detectors=((1, 3), (1, 5)),
        moduli=C5,
        label="C5",
    )

    evaluate_detector_intersection(
        test,
        detector_polynomials,
        detectors=((1, 3), (1, 7)),
        moduli=C5,
        label="C5",
    )

    evaluate_detector_intersection(
        test,
        detector_polynomials,
        detectors=((1, 3), (3, 5)),
        moduli=C5,
        label="C5",
    )

    evaluate_detector_intersection(
        test,
        detector_polynomials,
        detectors=((1, 3), (1, 5), (1, 7)),
        moduli=C7,
        label="C7",
    )

    evaluate_detector_intersection(
        test,
        detector_polynomials,
        detectors=DETECTORS,
        moduli=C7,
        label="C7",
    )

    # ------------------------------------------------------------------------
    # E1 BASELINE
    # ------------------------------------------------------------------------

    print("\n9. E1 / H6 BASELINE")
    print("-" * 78)

    def E1_value(n: int, s: int) -> int:
        return s * ((n + 1) ** 2 - s ** 2)

    def E1_mod(n: int, s: int, mod: int) -> int:
        nm = n % mod
        sm = s % mod
        return sm * ((nm + 1) ** 2 - sm ** 2) % mod

    def evaluate_E1(
        subset: Sequence[Target],
        moduli: Sequence[int],
    ) -> Tuple[float, float, int, int]:

        counts = []
        unique = 0
        survival = 0

        for target in subset:

            oracle = E1_direct = E1_value(
                target.n,
                target.s,
            )

            mask = np.ones(
                len(S_VALUES),
                dtype=np.bool_,
            )

            for mod in moduli:
                allowed = np.zeros(mod, dtype=np.bool_)
                target_residue = E1_direct % mod

                for r in range(mod):
                    if E1_mod(
                        target.n,
                        r,
                        mod,
                    ) == target_residue:
                        allowed[r] = True

                mask &= allowed[np.mod(S_VALUES, mod)]

            candidates = S_VALUES[mask]

            exact = [
                int(x)
                for x in candidates.tolist()
                if E1_value(
                    target.n,
                    int(x),
                ) == oracle
            ]

            counts.append(len(candidates))

            if len(exact) == 1:
                unique += 1

            if target.s in candidates:
                survival += 1

        return (
            statistics.mean(counts),
            statistics.median(counts),
            unique,
            survival,
        )

    for label, moduli in (
        ("C3", C3),
        ("C5", C5),
        ("C7", C7),
    ):

        tr = evaluate_E1(train, moduli)
        te = evaluate_E1(test, moduli)

        print(
            f"{label}: "
            f"train mean={tr[0]:.3f} "
            f"test mean={te[0]:.3f} "
            f"test unique={te[2]}/{len(test)} "
            f"test survival={te[3]}/{len(test)}"
        )

    # ------------------------------------------------------------------------
    # FINAL DIAGNOSTIC
    # ------------------------------------------------------------------------

    print("\n" + "=" * 78)
    print("10. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        """
The experiment isolates the paper's low-weight prime-detecting family.

The main question is NOT whether the coefficients detect primes;
the paper already establishes that property.

The question is whether the semiprime specialization of those same
coefficients produces a better factor-recovery resolvent than E1/H6.

A strong result would look like:

  * one low-weight f_(k,l) has substantially fewer s survivors;
  * the effect survives unseen targets;
  * it beats the E1/H6 baseline;
  * multiple paper detectors compress further when intersected;
  * controls do not show an equally strong advantage.

A particularly interesting outcome would be:

    f_(1,3) or another low-weight detector
        -> very low-degree polynomial in (n,s)
        -> strong modular collapse
        -> unique unseen s

without requiring a factor-pair test.

That would identify a concrete paper-derived resolvent worth trying
to connect to an N-only computation.

A failure is also useful:
if all f_(k,l) behave like E1/H6 or their intersections merely
reproduce the same information, then the entire low-weight
prime-detector family is probably not giving a fundamentally new
semiprime resolvent.
"""
    )

    print(
        f"\ntotal runtime = "
        f"{time.perf_counter() - t0:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 95 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

