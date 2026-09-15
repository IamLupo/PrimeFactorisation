#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
EXPERIMENT 505
==============================================================================
LEAVE-ONE-OUT FERMAT MIDPOINT LAW AUDIT
==============================================================================

GOAL
------------------------------------------------------------------------------
Test whether the Fermat midpoint

    A = (p+q)/2

can be reconstructed from N-only arithmetic / continued-fraction
observables by a genuinely determined low-complexity exact law.

IMPORTANT
------------------------------------------------------------------------------
Experiment 504 produced a misleading "linear hit" because the linear
system was underdetermined and SymPy returned free parameters tau0, tau1, ...

This experiment rejects such solutions.

A candidate is accepted only when:

  1. coefficients are uniquely determined;
  2. no free symbolic parameters remain;
  3. the law is exact on the training set;
  4. it predicts every held-out example exactly.

The factor pair is used ONLY to construct ground-truth A for validation.

No trial division or factor search is used to generate candidates.
==============================================================================
"""

from __future__ import annotations

import math
from fractions import Fraction
from itertools import combinations
from typing import Dict, List, Tuple

import sympy as sp


# ============================================================================
# TEST DATA
# ============================================================================

# Deliberately larger and more varied than Experiment 504.
# These are only validation pairs; candidate construction uses N-derived data.
PAIRS = [
    (50387, 282589),
    (1009, 10007),
    (10007, 1000003),
    (100003, 100019),
    (2000003, 3000017),
    (50021, 50047),
    (300007, 900001),

    (101, 1009),
    (211, 1009),
    (1009, 2003),
    (10007, 10009),
    (100003, 100019),
    (100003, 100129),
    (200003, 300007),
    (500009, 700001),
]


# ============================================================================
# INTEGER HELPERS
# ============================================================================

def ceil_sqrt(n: int) -> int:
    r = math.isqrt(n)
    return r if r * r == n else r + 1


def sqrt_cf(n: int, max_terms: int = 32) -> List[int]:
    a0 = math.isqrt(n)

    if a0 * a0 == n:
        return [a0]

    m = 0
    d = 1
    a = a0

    cf = [a0]

    for _ in range(max_terms - 1):
        m = d * a - m
        d = (n - m * m) // d
        a = (a0 + m) // d
        cf.append(a)

        if a == 2 * a0:
            break

    return cf


def convergents(cf: List[int]) -> List[Tuple[int, int]]:
    out = []

    p0, p1 = 0, 1
    q0, q1 = 1, 0

    for a in cf:
        p = a * p1 + p0
        q = a * q1 + q0

        out.append((p, q))

        p0, p1 = p1, p
        q0, q1 = q1, q

    return out


# ============================================================================
# N-ONLY FEATURES
# ============================================================================

def features_from_N(n: int) -> Dict[str, int]:
    floor_s = math.isqrt(n)
    ceil_s = ceil_sqrt(n)

    cf = sqrt_cf(n)
    conv = convergents(cf)

    f: Dict[str, int] = {
        "N": n,
        "sqrt_floor": floor_s,
        "sqrt_ceil": ceil_s,
        "r0": ceil_s * ceil_s - n,
        "floor_gap": n - floor_s * floor_s,
        "cf_len": len(cf),
    }

    # First few partial quotients.
    for i in range(min(8, len(cf))):
        f[f"cf{i}"] = cf[i]

    # First few convergents and signed residues.
    for i, (h, k) in enumerate(conv[:8]):
        residue = h * h - n * k * k

        f[f"h{i}"] = h
        f[f"k{i}"] = k
        f[f"res{i}"] = residue
        f[f"absres{i}"] = abs(residue)

        s = math.isqrt(abs(residue))
        f[f"res_sqrt{i}"] = s

        if s * s == abs(residue):
            f[f"res_square{i}"] = 1
        else:
            f[f"res_square{i}"] = 0

    return f


# ============================================================================
# DATASET
# ============================================================================

DATA: List[Dict[str, int]] = []

for p, q in PAIRS:
    n = p * q
    A = (p + q) // 2

    row = features_from_N(n)

    # Ground truth only.
    row["A_true"] = A
    row["B_true"] = abs(p - q) // 2
    row["d_true"] = A - row["sqrt_ceil"]

    DATA.append(row)


# ============================================================================
# FEATURE SETS
# ============================================================================

BASIC = [
    "N",
    "sqrt_floor",
    "sqrt_ceil",
    "r0",
    "floor_gap",
    "cf0",
    "cf1",
    "cf2",
    "cf3",
    "cf4",
    "cf5",
    "cf6",
    "cf7",
]

CONV = [
    "h0",
    "k0",
    "res0",
    "absres0",
    "h1",
    "k1",
    "res1",
    "absres1",
    "h2",
    "k2",
    "res2",
    "absres2",
    "h3",
    "k3",
    "res3",
    "absres3",
]

FEATURES = [
    x for x in BASIC + CONV
    if all(x in row for row in DATA)
]


# ============================================================================
# UNIQUE EXACT LINEAR SOLVER
# ============================================================================

def solve_unique_linear(
    train: List[Dict[str, int]],
    feature_names: List[str],
    target: str,
):
    """
    Solve

        target = c0 + c1*f1 + ... + ck*fk

    over Q.

    Return:
        coefficients if UNIQUE,
        None otherwise.

    A free parameter means the model is rejected.
    """

    M = []
    y = []

    for row in train:
        M.append(
            [sp.Integer(1)]
            + [sp.Integer(row[f]) for f in feature_names]
        )
        y.append(sp.Integer(row[target]))

    M = sp.Matrix(M)
    y = sp.Matrix(y)

    # Augmented system.
    aug = M.row_join(y)

    rank_M = M.rank()
    rank_aug = aug.rank()

    if rank_M != rank_aug:
        return None

    unknowns = len(feature_names) + 1

    # Unique solution requires full column rank.
    if rank_M != unknowns:
        return None

    try:
        sol = list(sp.linsolve((M, y)))
    except Exception:
        return None

    if len(sol) != 1:
        return None

    coeffs = list(sol[0])

    # Explicitly reject any residual symbolic/free parameters.
    free_symbols = set()

    for c in coeffs:
        free_symbols |= c.free_symbols

    if free_symbols:
        return None

    return coeffs


# ============================================================================
# EVALUATION
# ============================================================================

def evaluate_linear(
    row: Dict[str, int],
    feature_names: List[str],
    coeffs,
) -> sp.Rational:

    value = coeffs[0]

    for c, f in zip(coeffs[1:], feature_names):
        value += c * row[f]

    return sp.simplify(value)


# ============================================================================
# LEAVE-ONE-OUT TEST
# ============================================================================

def leave_one_out(
    feature_names: List[str],
    target: str = "A_true",
):
    accepted = []

    print()
    print(f"FEATURE FAMILY: {feature_names}")

    for holdout in range(len(DATA)):
        train = [
            row for i, row in enumerate(DATA)
            if i != holdout
        ]

        test = DATA[holdout]

        coeffs = solve_unique_linear(
            train,
            feature_names,
            target,
        )

        if coeffs is None:
            print(
                f"  holdout={holdout:02d}: "
                f"NO UNIQUE EXACT MODEL"
            )
            continue

        pred = evaluate_linear(
            test,
            feature_names,
            coeffs,
        )

        truth = sp.Integer(test[target])
        ok = (pred == truth)

        print(
            f"  holdout={holdout:02d}: "
            f"prediction={pred} "
            f"truth={truth} "
            f"PASS={ok}"
        )

        if ok:
            accepted.append(holdout)

    return accepted


# ============================================================================
# TARGETED ONE-FEATURE SEARCH
# ============================================================================

def exact_one_feature_search():
    print()
    print("[5] SINGLE-FEATURE EXACT SEARCH")
    print("-" * 78)

    hits = []

    for f in FEATURES:
        good = True

        for row in DATA:
            x = row[f]
            A = row["A_true"]

            # Only test direct identity f == A.
            if x != A:
                good = False
                break

        if good:
            hits.append(f)

    if not hits:
        print("  no feature equals A on the full dataset")
    else:
        for f in hits:
            print(f"  HIT: A = {f}")

    return hits


# ============================================================================
# SMALL AFFINE SEARCH
# ============================================================================

def affine_feature_search():
    print()
    print("[6] TWO-FEATURE AFFINE SEARCH")
    print("-" * 78)

    hits = []

    # A = c0 + c1*f
    for f in FEATURES:
        coeffs = solve_unique_linear(
            DATA,
            [f],
            "A_true",
        )

        if coeffs is None:
            continue

        predictions = [
            evaluate_linear(row, [f], coeffs)
            for row in DATA
        ]

        if all(
            pred == row["A_true"]
            for pred, row in zip(predictions, DATA)
        ):
            hits.append((f, coeffs))

    if not hits:
        print("  no exact one-feature affine law")
    else:
        for f, coeffs in hits:
            print(
                f"  HIT: A = {coeffs[0]} + "
                f"({coeffs[1]})*{f}"
            )

    return hits


# ============================================================================
# MAIN
# ============================================================================

print("EXPERIMENT 505 START")
print("=" * 78)
print("LEAVE-ONE-OUT FERMAT MIDPOINT LAW AUDIT")
print("=" * 78)


# ============================================================================
# 1. DATASET
# ============================================================================

print()
print("[1] DATASET")
print("-" * 78)

for i, row in enumerate(DATA):
    print(
        f"  {i:02d}: "
        f"N={row['N']} "
        f"ceil={row['sqrt_ceil']} "
        f"A={row['A_true']} "
        f"d={row['d_true']}"
    )


# ============================================================================
# 2. FERMAT SANITY
# ============================================================================

print()
print("[2] FERMAT SANITY")
print("-" * 78)

failures = 0

for row in DATA:
    A = row["A_true"]
    B = row["B_true"]
    n = row["N"]

    ok = (A * A - n == B * B)

    if not ok:
        failures += 1

print(f"  A^2-N=B^2 failures = {failures}")


# ============================================================================
# 3. FEATURE RANK
# ============================================================================

print()
print("[3] FEATURE MATRIX RANK")
print("-" * 78)

for count in [1, 2, 3, 4, 5, 6, 8, 10]:

    names = FEATURES[:count]

    M = sp.Matrix([
        [sp.Integer(1)] + [
            sp.Integer(row[f])
            for f in names
        ]
        for row in DATA
    ])

    print(
        f"  features={count:02d} "
        f"rank={M.rank():02d} "
        f"columns={M.cols:02d}"
    )


# ============================================================================
# 4. EXPLICITLY TEST 504'S FAILURE MODE
# ============================================================================

print()
print("[4] FREE-PARAMETER REJECTION TEST")
print("-" * 78)

names = FEATURES[:10]

M = sp.Matrix([
    [sp.Integer(1)] + [
        sp.Integer(row[f])
        for f in names
    ]
    for row in DATA[:5]
])

y = sp.Matrix([
    row["A_true"]
    for row in DATA[:5]
])

solution = sp.linsolve((M, y))

print(f"  equations = {M.rows}")
print(f"  unknowns  = {M.cols}")
print(f"  rank      = {M.rank()}")

if solution is sp.EmptySet:
    print("  no solution")
else:
    sol = next(iter(solution))
    free = set()

    for x in sol:
        free |= x.free_symbols

    print(f"  free symbols = {sorted(map(str, free))}")

    if free:
        print("  UNDERDETERMINED -> REJECT")
    else:
        print("  UNIQUE -> candidate may be tested")


# ============================================================================
# 5. SINGLE FEATURE
# ============================================================================

single_hits = exact_one_feature_search()


# ============================================================================
# 6. AFFINE
# ============================================================================

affine_hits = affine_feature_search()


# ============================================================================
# 7. LEAVE-ONE-OUT WITH SMALL FEATURE FAMILIES
# ============================================================================

print()
print("[7] LEAVE-ONE-OUT EXACT VALIDATION")
print("-" * 78)

families = [
    ["sqrt_ceil"],
    ["sqrt_floor"],
    ["r0"],
    ["cf0"],
    ["cf1"],
    ["cf2"],
    ["sqrt_ceil", "r0"],
    ["sqrt_ceil", "cf1"],
    ["sqrt_ceil", "cf1", "cf2"],
    ["sqrt_ceil", "r0", "cf1"],
    ["sqrt_ceil", "r0", "cf1", "cf2"],
]

loo_results = []

for family in families:
    accepted = leave_one_out(
        family,
        "A_true",
    )

    loo_results.append(
        (family, accepted)
    )


# ============================================================================
# 8. OFFSET TARGET
# ============================================================================

print()
print("[8] LEAVE-ONE-OUT SEARCH FOR FERMAT OFFSET d")
print("-" * 78)

offset_results = []

for family in families:
    accepted = leave_one_out(
        family,
        "d_true",
    )

    offset_results.append(
        (family, accepted)
    )


# ============================================================================
# 9. RESULTS
# ============================================================================

print()
print("[9] FINAL RESULTS")
print("-" * 78)

full = len(DATA)

print(f"  dataset size = {full}")
print(f"  single-feature A hits = {len(single_hits)}")
print(f"  affine A hits         = {len(affine_hits)}")

surviving_A = [
    family
    for family, accepted in loo_results
    if len(accepted) == full
]

surviving_d = [
    family
    for family, accepted in offset_results
    if len(accepted) == full
]

print()
print("  LOO-surviving A laws:")

if not surviving_A:
    print("    NONE")
else:
    for family in surviving_A:
        print(f"    {family}")

print()
print("  LOO-surviving d laws:")

if not surviving_d:
    print("    NONE")
else:
    for family in surviving_d:
        print(f"    {family}")


# ============================================================================
# 10. INTERPRETATION
# ============================================================================

print()
print("[10] INTERPRETATION")
print("-" * 78)

if surviving_A:
    print(
        """
A candidate exact law survived leave-one-out validation.
This is the first result in this line of experiments that
would justify deeper symbolic investigation.
"""
    )
else:
    print(
        """
No tested low-complexity N-only law survived leave-one-out
validation.

This is importantly stronger than Experiment 504:

  * underdetermined solutions are rejected;
  * free tau parameters are rejected;
  * interpolation on the complete dataset is not sufficient;
  * held-out exact prediction is required.

The Fermat identity itself remains valid, but no tested
continued-fraction feature family has yet supplied A=(p+q)/2.
"""
    )

print()
print("===============================================================================")
print("EXPERIMENT 505 FINISHED")
print("===============================================================================")
