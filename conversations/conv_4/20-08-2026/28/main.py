#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
EXPERIMENT 504
==============================================================================
FERMAT MIDPOINT FROM N-ONLY CONTINUED-FRACTION SIGNATURES
==============================================================================

TARGET
------------------------------------------------------------------------------
For N = p*q define

    A = (p+q)/2
    B = (q-p)/2

so that

    A^2 - N = B^2.

The missing upstream quantity is therefore

    A

or equivalently

    d = A - ceil(sqrt(N)).

This experiment asks:

    Can d, or A, be predicted from N-only continued-fraction data?

The factor pair p,q is used ONLY as ground truth for validation.

No factor search is performed.

==============================================================================
"""

from __future__ import annotations

import math
from fractions import Fraction
from typing import Iterable

import sympy as sp


# ============================================================================
# SYMBOLS
# ============================================================================

N, A, B, C, R, D, Q = sp.symbols(
    "N A B C R D Q",
    integer=True
)


# ============================================================================
# TEST CASES
# ============================================================================

TEST_PAIRS = [
    (50387, 282589),
    (1009, 10007),
    (10007, 1000003),
    (100003, 100019),
    (2000003, 3000017),
    (50021, 50047),
    (300007, 900001),
]


# ============================================================================
# BASIC INTEGER SQRT
# ============================================================================

def isqrt(n: int) -> int:
    return math.isqrt(n)


def ceil_sqrt(n: int) -> int:
    r = math.isqrt(n)
    return r if r * r == n else r + 1


# ============================================================================
# CONTINUED FRACTION OF sqrt(N)
# ============================================================================

def sqrt_cf(n: int, max_terms: int = 128) -> list[int]:
    """
    Continued fraction coefficients of sqrt(n).

    For nonsquare n:
        sqrt(n) = [a0; a1, a2, ...]
    """
    a0 = math.isqrt(n)

    if a0 * a0 == n:
        return [a0]

    m = 0
    d = 1
    a = a0

    coeffs = [a0]

    for _ in range(max_terms - 1):
        m = d * a - m
        d = (n - m * m) // d
        a = (a0 + m) // d
        coeffs.append(a)

        # sqrt(N) has periodic CF.
        if a == 2 * a0:
            break

    return coeffs


# ============================================================================
# CONVERGENTS
# ============================================================================

def convergents(cf: list[int]) -> list[tuple[int, int]]:
    """
    Return convergents h/k.
    """
    out: list[tuple[int, int]] = []

    p_minus_2, p_minus_1 = 0, 1
    q_minus_2, q_minus_1 = 1, 0

    for a in cf:
        p = a * p_minus_1 + p_minus_2
        q = a * q_minus_1 + q_minus_2

        out.append((p, q))

        p_minus_2, p_minus_1 = p_minus_1, p
        q_minus_2, q_minus_1 = q_minus_1, q

    return out


# ============================================================================
# GROUND-TRUTH FERMAT DATA
# ============================================================================

def fermat_data(p: int, q: int) -> dict[str, int]:
    n = p * q
    a = (p + q) // 2
    b = abs(p - q) // 2

    c = ceil_sqrt(n)

    assert a * a - n == b * b

    return {
        "p": p,
        "q": q,
        "N": n,
        "A": a,
        "B": b,
        "C": c,
        "d": a - c,
        "r0": c * c - n,
        "delta": (p - q) ** 2,
    }


# ============================================================================
# FEATURE EXTRACTION
# ============================================================================

def extract_features(n: int, cf: list[int], convs: list[tuple[int, int]]) -> dict[str, int]:
    c = ceil_sqrt(n)
    floor_r = math.isqrt(n)

    features: dict[str, int] = {
        "N": n,
        "sqrt_floor": floor_r,
        "sqrt_ceil": c,
        "r0": c * c - n,
        "ceil_minus_floor": c - floor_r,
        "cf_len": len(cf),
    }

    # First partial quotients.
    for i, a in enumerate(cf[:12]):
        features[f"cf{i}"] = a

    # Small convergent features.
    for i, (h, k) in enumerate(convs[:12]):
        residue = h * h - n * k * k

        features[f"h{i}"] = h
        features[f"k{i}"] = k
        features[f"res{i}"] = residue

        # Absolute residue and signed distance from a square.
        features[f"absres{i}"] = abs(residue)

        sq = math.isqrt(abs(residue))
        features[f"res_sqrt{i}"] = sq
        features[f"res_square{i}"] = int(sq * sq == abs(residue))

    return features


# ============================================================================
# SMALL INTEGER FEATURE SEARCH
# ============================================================================

def normalize(expr: sp.Expr) -> sp.Expr:
    return sp.factor(sp.cancel(expr))


def polynomial_feature_search(
    rows: list[tuple[dict[str, int], int]],
    feature_names: list[str],
    target_name: str,
    max_features: int = 6,
) -> list[tuple[str, sp.Expr]]:
    """
    Search exact formulas of the form

        target = c0 + sum ci * feature_i

    using rational coefficients.

    This is deliberately linear: the point is to detect unexpectedly
    simple transfer laws before considering nonlinear families.
    """
    selected = feature_names[:max_features]

    if not selected:
        return []

    symbols = [sp.Symbol(name) for name in selected]

    # Matrix over Q:
    # [1, f1, f2, ...] * coeffs = target
    M = []
    y = []

    for feats, target in rows:
        M.append(
            [sp.Integer(1)]
            + [sp.Integer(feats[name]) for name in selected]
        )
        y.append(sp.Integer(target))

    M = sp.Matrix(M)
    y = sp.Matrix(y)

    try:
        solution = sp.linsolve((M, y))
    except Exception:
        return []

    if solution is sp.EmptySet:
        return []

    results: list[tuple[str, sp.Expr]] = []

    for sol in solution:
        if len(sol) != len(selected) + 1:
            continue

        expr = sol[0]

        for coeff, sym in zip(sol[1:], symbols):
            expr += coeff * sym

        # Verify exactly against all data.
        ok = True

        for feats, target in rows:
            value = expr.subs(
                {
                    sym: feats[name]
                    for sym, name in zip(symbols, selected)
                }
            )

            if sp.simplify(value - target) != 0:
                ok = False
                break

        if ok:
            results.append(("linear", normalize(expr)))

    return results


# ============================================================================
# SPECIALIZED TRANSFER SEARCH
# ============================================================================

def search_simple_relations(rows: list[dict[str, int]]) -> None:
    """
    Search targeted expressions for the Fermat offset d.
    """
    print()
    print("[8] TARGETED FERMAT-OFFSET SEARCH")
    print("-" * 78)

    hits: list[tuple[str, list[int]]] = []

    for i, row in enumerate(rows):
        n = row["N"]
        c = row["sqrt_ceil"]
        d = row["target_d"]

        candidates = {
            "r0": row["r0"],
            "r0/2": row["r0"] // 2,
            "sqrt(r0)": row["r0_sqrt"],
            "cf0": row["cf0"],
            "cf1": row.get("cf1", 0),
            "cf2": row.get("cf2", 0),
            "cf3": row.get("cf3", 0),
        }

        print(f"  case {i + 1}: d={d}")

        for name, value in candidates.items():
            if value == d:
                print(f"    HIT: d = {name}")
                hits.append((name, [i + 1]))

    if not hits:
        print("  no direct single-feature hits")
    else:
        print(f"  direct hits = {len(hits)}")


# ============================================================================
# RESIDUE-SQUARE SEARCH
# ============================================================================

def square_residue_observables(
    n: int,
    convs: list[tuple[int, int]],
    target_a: int,
) -> list[dict[str, int]]:
    out = []

    for idx, (h, k) in enumerate(convs):
        residue = h * h - n * k * k

        if residue < 0:
            continue

        r = math.isqrt(residue)

        if r * r != residue:
            continue

        # Factor-exposure channel.
        g1 = math.gcd(h - r, n)
        g2 = math.gcd(h + r, n)

        out.append(
            {
                "idx": idx,
                "h": h,
                "k": k,
                "r": r,
                "g1": g1,
                "g2": g2,
                "is_fermat_A": int(h == target_a and k == 1),
                "h_minus_A": h - target_a,
                "k_minus_1": k - 1,
            }
        )

    return out


# ============================================================================
# MAIN
# ============================================================================

print("EXPERIMENT 504 START")
print("=" * 78)
print("FERMAT MIDPOINT FROM N-ONLY CONTINUED-FRACTION SIGNATURES")
print("=" * 78)


# ============================================================================
# 1. FERMAT IDENTITIES
# ============================================================================

print()
print("[1] FERMAT MIDPOINT IDENTITIES")
print("-" * 78)

failures = 0

for p, q in TEST_PAIRS:
    data = fermat_data(p, q)

    n = data["N"]
    a = data["A"]
    b = data["B"]

    pass_a = (a * a - n == b * b)
    pass_n = (a * a - b * b == n)

    if not pass_a or not pass_n:
        failures += 1

    print(
        f"  ({p},{q}) "
        f"A={a} B={b} "
        f"A^2-N=B^2={b*b} "
        f"PASS={pass_a and pass_n}"
    )

print(f"  IDENTITY FAILURES = {failures}")


# ============================================================================
# 2. N-ONLY FEATURE EXTRACTION
# ============================================================================

print()
print("[2] N-ONLY CONTINUED-FRACTION FEATURES")
print("-" * 78)

rows: list[dict[str, int]] = []

for p, q in TEST_PAIRS:
    data = fermat_data(p, q)
    n = data["N"]

    cf = sqrt_cf(n, max_terms=32)
    convs = convergents(cf)

    feats = extract_features(n, cf, convs)

    feats["target_A"] = data["A"]
    feats["target_d"] = data["d"]
    feats["target_B"] = data["B"]

    # r0 square root if possible.
    r0 = feats["r0"]
    r0_sqrt = math.isqrt(r0)
    feats["r0_sqrt"] = r0_sqrt if r0_sqrt * r0_sqrt == r0 else -1

    rows.append(feats)

    print(f"  ({p},{q})")
    print(f"    N           = {n}")
    print(f"    ceil(sqrtN) = {data['C']}")
    print(f"    A           = {data['A']}")
    print(f"    d=A-ceil    = {data['d']}")
    print(f"    r0          = {data['r0']}")
    print(f"    CF prefix   = {cf[:8]}")


# ============================================================================
# 3. DISTANCE FROM sqrt(N)
# ============================================================================

print()
print("[3] FERMAT OFFSET PROFILE")
print("-" * 78)

for row in rows:
    print(
        f"  N={row['N']}"
        f"  ceil={row['sqrt_ceil']}"
        f"  A={row['target_A']}"
        f"  d={row['target_d']}"
        f"  r0={row['r0']}"
    )


# ============================================================================
# 4. CONVERGENT SIGNATURE SEARCH
# ============================================================================

print()
print("[4] CONVERGENT SIGNATURE SEARCH")
print("-" * 78)

for p, q in TEST_PAIRS:
    data = fermat_data(p, q)
    n = data["N"]

    cf = sqrt_cf(n, max_terms=32)
    convs = convergents(cf)

    print(f"  ({p},{q})")

    found = False

    for idx, (h, k) in enumerate(convs[:16]):
        residue = h * h - n * k * k

        if residue >= 0:
            r = math.isqrt(residue)

            if r * r == residue:
                found = True

                g1 = math.gcd(h - r, n)
                g2 = math.gcd(h + r, n)

                print(
                    f"    idx={idx}"
                    f" h={h}"
                    f" k={k}"
                    f" r={r}"
                    f" gcds=({g1},{g2})"
                    f" A-match={h == data['A'] and k == 1}"
                )

    if not found:
        print("    no exact square residues in first 16 convergents")


# ============================================================================
# 5. DIRECT FEATURE -> A SEARCH
# ============================================================================

print()
print("[5] LINEAR N-ONLY A SEARCH")
print("-" * 78)

candidate_features = [
    "N",
    "sqrt_floor",
    "sqrt_ceil",
    "r0",
    "cf0",
    "cf1",
    "cf2",
    "cf3",
    "cf4",
    "cf5",
    "k0",
    "k1",
    "k2",
    "k3",
    "res0",
    "res1",
    "res2",
    "res3",
]

candidate_features = [
    name
    for name in candidate_features
    if all(name in row for row in rows)
]

hits = polynomial_feature_search(
    [(row, row["target_A"]) for row in rows],
    candidate_features,
    "target_A",
    max_features=min(10, len(candidate_features)),
)

if not hits:
    print("  no exact linear A representation found")
else:
    for kind, expr in hits:
        print(f"  HIT [{kind}]: A = {expr}")


# ============================================================================
# 6. DIRECT FEATURE -> d SEARCH
# ============================================================================

print()
print("[6] LINEAR N-ONLY FERMAT-OFFSET SEARCH")
print("-" * 78)

offset_features = [
    "N",
    "sqrt_floor",
    "sqrt_ceil",
    "r0",
    "cf0",
    "cf1",
    "cf2",
    "cf3",
    "cf4",
    "cf5",
    "k0",
    "k1",
    "k2",
    "k3",
    "res0",
    "res1",
    "res2",
    "res3",
]

offset_features = [
    name
    for name in offset_features
    if all(name in row for row in rows)
]

hits_d = polynomial_feature_search(
    [(row, row["target_d"]) for row in rows],
    offset_features,
    "target_d",
    max_features=min(10, len(offset_features)),
)

if not hits_d:
    print("  no exact linear d representation found")
else:
    for kind, expr in hits_d:
        print(f"  HIT [{kind}]: d = {expr}")


# ============================================================================
# 7. SMALL NORMALIZED SIGNATURE SEARCH
# ============================================================================

print()
print("[7] NORMALIZED SIGNATURE SEARCH")
print("-" * 78)

signature_hits = 0

for row in rows:
    n = row["N"]
    c = row["sqrt_ceil"]
    d = row["target_d"]
    r0 = row["r0"]

    signatures = {
        "(r0+1)/2": Fraction(r0 + 1, 2),
        "r0/c": Fraction(r0, c),
        "r0/(2*c)": Fraction(r0, 2 * c),
        "d/c": Fraction(d, c),
        "d/r0": Fraction(d, r0) if r0 != 0 else None,
        "r0/(2*d+1)": (
            Fraction(r0, 2 * d + 1)
            if 2 * d + 1 != 0
            else None
        ),
    }

    hits_here = []

    for name, value in signatures.items():
        if value == d:
            hits_here.append(name)

    if hits_here:
        signature_hits += len(hits_here)

        print(
            f"  N={n}"
            f" d={d}"
            f" hits={hits_here}"
        )

if signature_hits == 0:
    print("  no normalized direct hits")


# ============================================================================
# 8. SQUARE-DISTANCE TEST
# ============================================================================

print()
print("[8] SQUARE-DISTANCE TEST")
print("-" * 78)

square_hits = 0

for row in rows:
    n = row["N"]
    a = row["target_A"]
    b = row["target_B"]

    x = a * a - n
    root = math.isqrt(x)

    passed = root * root == x and root == b

    if passed:
        square_hits += 1

    print(
        f"  N={n}"
        f"  A^2-N={x}"
        f"  sqrt={root}"
        f"  B={b}"
        f"  PASS={passed}"
    )


# ============================================================================
# 9. SUMMARY
# ============================================================================

print()
print("[9] EXPERIMENT STATUS")
print("-" * 78)

print(f"  Fermat identity failures          = {failures}")
print(f"  linear A hits                     = {len(hits)}")
print(f"  linear d hits                     = {len(hits_d)}")
print(f"  normalized signature hits         = {signature_hits}")
print(f"  square-distance confirmations     = {square_hits}")

print()
print("MAIN QUESTION")
print("-" * 78)

print(
    """
Can an N-only continued-fraction observable predict

    A = (p+q)/2

without first recovering p or q?

The experiment separates:

    N
     -> CF / sqrt(N) observables
     -> candidate A
     -> A^2-N
     -> factor gap
     -> p,q

from the already-known downstream identity

    A^2-N = ((p-q)/2)^2.

Only exact relations surviving all test cases are reported as hits.
"""
)

print()
print("=" * 78)
print("EXPERIMENT 504 FINISHED")
print("=" * 78)
