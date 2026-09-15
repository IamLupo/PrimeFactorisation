#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
==============================================================================
KAPPA EXPERIMENT 122
C/D TOP-EDGE LAW
TOP N-DEGREE COEFFICIENTS ONLY
EXACT SYMBOLIC ARITHMETIC
NO FULL TENSOR SEARCH
NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
SAFE SYMPY DOMAINS
==============================================================================

Purpose
-------
Experiment 121 established the exact decomposition

    Q(N,S) = C(N,X) + (2S-1)D(N,X),
    X = S(S-1).

But constructing and validating the complete symbolic C/D tensor is expensive.

Experiment 122 asks a narrower and more promising question:

    What is the exact law for the TOP N-DEGREE coefficient
    of C and D?

and, when available,

    what is the exact law for the SECOND N-DEGREE coefficient?

The search is performed in

    r = (k-1)/2
    m = (ell-1)/2

using exact rational interpolation.

The analysis deliberately does NOT:
    * reconstruct the whole C/D tensor for every coefficient;
    * perform factorization;
    * perform full LOO over every monomial slot;
    * repeatedly invoke SymPy polynomial division.

The primary tests are:

    1. top-N coefficient law in (r,m)
    2. second-N coefficient law in (r,m)
    3. leave-one-ell-out validation
    4. forward ell holdout

A positive result would identify the highest-degree kernel directly.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import sympy as sp


# ============================================================================
# PARAMETERS
# ============================================================================

MAX_ELL = 17
FORWARD_ELL = 19

MAX_BIV_DEG = 4

# Only the odd detector weights.
K_VALUES = None


# ============================================================================
# SYMBOLS
# ============================================================================

N, S, X = sp.symbols("N S X")
r_sym, m_sym = sp.symbols("r m")


# ============================================================================
# DATA
# ============================================================================

@dataclass(frozen=True)
class EdgeRecord:
    k: int
    ell: int
    r: int
    m: int

    # C top edge
    C_top_degree: int
    C_top_X_degree: int
    C_top_coeffs: Tuple[sp.Rational, ...]

    # C second N layer
    C_second_degree: Optional[int]
    C_second_X_degree: Optional[int]
    C_second_coeffs: Tuple[sp.Rational, ...]

    # D top edge
    D_top_degree: Optional[int]
    D_top_X_degree: Optional[int]
    D_top_coeffs: Tuple[sp.Rational, ...]

    # D second N layer
    D_second_degree: Optional[int]
    D_second_X_degree: Optional[int]
    D_second_coeffs: Tuple[sp.Rational, ...]


# ============================================================================
# POWER SUMS
# ============================================================================

def power_sum_table(max_power: int) -> List[sp.Expr]:
    P = [sp.Integer(0)] * (max_power + 1)

    P[0] = sp.Integer(2)

    if max_power >= 1:
        P[1] = S

    for j in range(2, max_power + 1):
        P[j] = sp.expand(S * P[j - 1] - N * P[j - 2])

    return P


# ============================================================================
# SYMMETRIC MONOMIAL
# ============================================================================

def symmetric_monomial(
    a: int,
    b: int,
    P: Sequence[sp.Expr],
) -> sp.Expr:

    if a < b:
        a, b = b, a

    return sp.expand(N**b * P[a - b])


# ============================================================================
# RAW DETECTOR NUMERATOR
# ============================================================================

def raw_detector_numerator(k: int, ell: int) -> sp.Expr:

    P = power_sum_table(max(k, ell))

    F = sp.Integer(0)

    # p^k (1+q)^ell + q^k (1+p)^ell
    for t in range(ell + 1):
        F += sp.binomial(ell, t) * symmetric_monomial(k, t, P)

    # - p^ell (1+q)^k - q^ell (1+p)^k
    for t in range(k + 1):
        F -= sp.binomial(k, t) * symmetric_monomial(ell, t, P)

    return sp.expand(F)


# ============================================================================
# EXACT Q
# ============================================================================

def quotient_polynomial(k: int, ell: int) -> sp.Expr:

    F = raw_detector_numerator(k, ell)

    poly = sp.Poly(
        F,
        S,
        domain=sp.ZZ.poly_ring(N),
    )

    divisor = sp.Poly(
        S + 1,
        S,
        domain=sp.ZZ.poly_ring(N),
    )

    q, rem = sp.div(poly, divisor)

    remainder = sp.expand(rem.as_expr())

    if remainder != 0:
        raise ArithmeticError(
            f"Nonzero quotient remainder ({k},{ell}): "
            f"{sp.factor(remainder)}"
        )

    return sp.expand(q.as_expr())


# ============================================================================
# QUADRATIC REDUCTION
# ============================================================================

def reduce_mod_quadratic(expr: sp.Expr) -> sp.Expr:

    poly = sp.Poly(
        sp.expand(expr),
        S,
        domain=sp.QQ.poly_ring(N, X),
    )

    modulus = sp.Poly(
        S**2 - S - X,
        S,
        domain=sp.QQ.poly_ring(N, X),
    )

    return sp.expand(poly.rem(modulus).as_expr())


def branch_decompose(Q: sp.Expr) -> Tuple[sp.Expr, sp.Expr]:

    R = reduce_mod_quadratic(Q)

    poly = sp.Poly(
        R,
        S,
        domain=sp.QQ.poly_ring(N, X),
    )

    A = sp.expand(poly.coeff_monomial(S**0))
    B = sp.expand(poly.coeff_monomial(S**1))

    D = sp.expand(B / 2)
    C = sp.expand(A + B / 2)

    return C, D


# ============================================================================
# EDGE EXTRACTION
# ============================================================================

def polynomial_x_coefficients(expr: sp.Expr) -> Dict[int, sp.Expr]:
    """
    Treat expr as polynomial in X with coefficients in QQ[N].
    Return x_degree -> coefficient.
    """

    if expr == 0:
        return {}

    poly = sp.Poly(
        sp.expand(expr),
        X,
        domain=sp.QQ.poly_ring(N),
    )

    result: Dict[int, sp.Expr] = {}

    for (x_degree,), coeff in poly.terms():
        result[x_degree] = sp.expand(coeff)

    return result


def n_coefficients(expr: sp.Expr) -> Dict[int, sp.Expr]:
    """
    Return N-degree -> coefficient in QQ[X].
    """

    if expr == 0:
        return {}

    poly = sp.Poly(
        sp.expand(expr),
        N,
        domain=sp.QQ.poly_ring(X),
    )

    result: Dict[int, sp.Expr] = {}

    for (n_degree,), coeff in poly.terms():
        result[n_degree] = sp.expand(coeff)

    return result


def edge_layer(expr: sp.Expr) -> Tuple[
    Optional[int],
    Optional[int],
    Tuple[sp.Rational, ...]
]:
    """
    Extract the highest N-degree layer.

    Return:
        N degree,
        X degree,
        coefficients ordered X^0, X^1, ...

    We do not factor the polynomial.
    """

    if expr == 0:
        return None, None, tuple()

    nmap = n_coefficients(expr)

    if not nmap:
        return None, None, tuple()

    top_n = max(nmap)
    layer = sp.Poly(
        sp.expand(nmap[top_n]),
        X,
        domain=sp.QQ,
    )

    if layer.is_zero:
        return None, None, tuple()

    top_x = layer.degree()

    coeffs = tuple(
        sp.Rational(layer.coeff_monomial(X**j))
        for j in range(top_x + 1)
    )

    return top_n, top_x, coeffs


def second_layer(expr: sp.Expr) -> Tuple[
    Optional[int],
    Optional[int],
    Tuple[sp.Rational, ...]
]:
    """
    Extract the second-highest N-degree layer.
    """

    if expr == 0:
        return None, None, tuple()

    nmap = n_coefficients(expr)

    if len(nmap) < 2:
        return None, None, tuple()

    degrees = sorted(nmap)

    second_n = degrees[-2]

    layer = sp.Poly(
        sp.expand(nmap[second_n]),
        X,
        domain=sp.QQ,
    )

    if layer.is_zero:
        return None, None, tuple()

    top_x = layer.degree()

    coeffs = tuple(
        sp.Rational(layer.coeff_monomial(X**j))
        for j in range(top_x + 1)
    )

    return second_n, top_x, coeffs


# ============================================================================
# DATASET
# ============================================================================

def build_dataset(max_ell: int) -> List[EdgeRecord]:

    records: List[EdgeRecord] = []

    for ell in range(3, max_ell + 1, 2):

        ks = (
            list(range(1, ell, 2))
            if K_VALUES is None
            else [k for k in K_VALUES if k < ell]
        )

        for k in ks:

            Q = quotient_polynomial(k, ell)
            C, D = branch_decompose(Q)

            ct_n, ct_x, ct_coeffs = edge_layer(C)
            cs_n, cs_x, cs_coeffs = second_layer(C)

            dt_n, dt_x, dt_coeffs = edge_layer(D)
            ds_n, ds_x, ds_coeffs = second_layer(D)

            records.append(
                EdgeRecord(
                    k=k,
                    ell=ell,
                    r=(k - 1) // 2,
                    m=(ell - 1) // 2,

                    C_top_degree=(
                        int(ct_n) if ct_n is not None else -1
                    ),
                    C_top_X_degree=(
                        int(ct_x) if ct_x is not None else -1
                    ),
                    C_top_coeffs=ct_coeffs,

                    C_second_degree=(
                        int(cs_n) if cs_n is not None else None
                    ),
                    C_second_X_degree=(
                        int(cs_x) if cs_x is not None else None
                    ),
                    C_second_coeffs=cs_coeffs,

                    D_top_degree=(
                        int(dt_n) if dt_n is not None else None
                    ),
                    D_top_X_degree=(
                        int(dt_x) if dt_x is not None else None
                    ),
                    D_top_coeffs=dt_coeffs,

                    D_second_degree=(
                        int(ds_n) if ds_n is not None else None
                    ),
                    D_second_X_degree=(
                        int(ds_x) if ds_x is not None else None
                    ),
                    D_second_coeffs=ds_coeffs,
                )
            )

    return records


# ============================================================================
# EDGE VALIDATION
# ============================================================================

def print_degree_profile(records: Sequence[EdgeRecord]) -> None:

    print()
    print("=" * 78)
    print("EDGE DEGREE PROFILE")
    print("=" * 78)

    by_kind = {
        "C-top": [r.C_top_degree for r in records],
        "C-second": [
            r.C_second_degree
            for r in records
            if r.C_second_degree is not None
        ],
        "D-top": [
            r.D_top_degree
            for r in records
            if r.D_top_degree is not None
        ],
        "D-second": [
            r.D_second_degree
            for r in records
            if r.D_second_degree is not None
        ],
    }

    for name, vals in by_kind.items():
        print(
            f"{name:10s}: "
            f"min={min(vals) if vals else '-':>3} "
            f"max={max(vals) if vals else '-':>3} "
            f"distinct={sorted(set(vals))}"
        )


def print_edge_samples(records: Sequence[EdgeRecord]) -> None:

    print()
    print("=" * 78)
    print("EDGE SAMPLES")
    print("=" * 78)

    for rec in records[:20]:

        ctop = rec.C_top_coeffs
        dtop = rec.D_top_coeffs

        print(
            f"({rec.k},{rec.ell}) "
            f"r={rec.r} m={rec.m}"
        )

        print(
            f"  C-top: degree={rec.C_top_degree} "
            f"X-degree={rec.C_top_X_degree} "
            f"coeffs={ctop}"
        )

        print(
            f"  D-top: degree={rec.D_top_degree} "
            f"X-degree={rec.D_top_X_degree} "
            f"coeffs={dtop}"
        )


# ============================================================================
# SIMPLE BIVARIATE BASIS
# ============================================================================

def monomial_basis(degree: int) -> List[Tuple[int, int]]:
    basis = []

    for i in range(degree + 1):
        for j in range(degree + 1 - i):
            basis.append((i, j))

    return basis


def exact_fit(
    points: Sequence[Tuple[int, int]],
    values: Sequence[sp.Rational],
    degree: int,
) -> Optional[sp.Expr]:

    basis = monomial_basis(degree)

    if len(points) < len(basis):
        return None

    A = sp.Matrix([
        [
            sp.Integer(r) ** i * sp.Integer(m) ** j
            for i, j in basis
        ]
        for r, m in points
    ])

    b = sp.Matrix(values)

    try:
        result = sp.linsolve((A, b))
    except Exception:
        return None

    tuples = list(result)

    if len(tuples) != 1:
        return None

    sol = tuples[0]

    if any(v.free_symbols for v in sol):
        return None

    out = sp.Integer(0)

    for coeff, (i, j) in zip(sol, basis):
        out += coeff * r_sym**i * m_sym**j

    return sp.expand(out)


def evaluate(poly: sp.Expr, r: int, m: int) -> sp.Rational:
    return sp.Rational(
        sp.expand(
            poly.subs(
                {
                    r_sym: r,
                    m_sym: m,
                }
            )
        )
    )


# ============================================================================
# COEFFICIENT LAW SEARCH
# ============================================================================

def extract_component(
    record: EdgeRecord,
    kind: str,
    x_index: int,
) -> Optional[sp.Rational]:

    if kind == "C-top":
        vals = record.C_top_coeffs
    elif kind == "C-second":
        vals = record.C_second_coeffs
    elif kind == "D-top":
        vals = record.D_top_coeffs
    elif kind == "D-second":
        vals = record.D_second_coeffs
    else:
        raise ValueError(kind)

    if x_index >= len(vals):
        return sp.Rational(0)

    return sp.Rational(vals[x_index])


def search_component_law(
    records: Sequence[EdgeRecord],
    kind: str,
    x_index: int,
    max_degree: int,
) -> Tuple[Optional[sp.Expr], Optional[int], int]:

    usable = []

    for rec in records:
        value = extract_component(rec, kind, x_index)

        if value is not None:
            usable.append((rec.r, rec.m, value))

    if len(usable) < 6:
        return None, None, 0

    points = [(r, m) for r, m, _ in usable]
    values = [v for _, _, v in usable]

    for degree in range(max_degree + 1):
        poly = exact_fit(points, values, degree)

        if poly is None:
            continue

        failures = 0

        for r, m, actual in usable:
            if evaluate(poly, r, m) != actual:
                failures += 1

        if failures == 0:
            return poly, degree, len(usable)

    return None, None, len(usable)


# ============================================================================
# LOO ELL TEST
# ============================================================================

def loo_test(
    records: Sequence[EdgeRecord],
    kind: str,
    x_index: int,
    degree: int,
) -> Tuple[bool, int]:

    unique_m = sorted(set(rec.m for rec in records))

    total = 0
    failures = 0

    for holdout_m in unique_m:

        train = [
            rec for rec in records
            if rec.m != holdout_m
            and extract_component(rec, kind, x_index) is not None
        ]

        test = [
            rec for rec in records
            if rec.m == holdout_m
            and extract_component(rec, kind, x_index) is not None
        ]

        if not train or not test:
            continue

        points = [(r.r, r.m) for r in train]
        values = [
            extract_component(r, kind, x_index)
            for r in train
        ]

        poly = exact_fit(points, values, degree)

        if poly is None:
            failures += len(test)
            total += len(test)
            continue

        for rec in test:
            actual = extract_component(
                rec,
                kind,
                x_index,
            )

            predicted = evaluate(
                poly,
                rec.r,
                rec.m,
            )

            total += 1

            if predicted != actual:
                failures += 1

    return failures == 0, failures


# ============================================================================
# FORWARD ELL TEST
# ============================================================================

def forward_holdout(
    train_records: Sequence[EdgeRecord],
    future_records: Sequence[EdgeRecord],
    kind: str,
    x_index: int,
    degree: int,
) -> Tuple[bool, int]:

    train = [
        rec for rec in train_records
        if extract_component(rec, kind, x_index) is not None
    ]

    future = [
        rec for rec in future_records
        if extract_component(rec, kind, x_index) is not None
    ]

    if not train or not future:
        return False, 0

    points = [(r.r, r.m) for r in train]

    values = [
        extract_component(r, kind, x_index)
        for r in train
    ]

    poly = exact_fit(points, values, degree)

    if poly is None:
        return False, len(future)

    failures = 0

    for rec in future:
        actual = extract_component(
            rec,
            kind,
            x_index,
        )

        predicted = evaluate(
            poly,
            rec.r,
            rec.m,
        )

        if predicted != actual:
            failures += 1

    return failures == 0, failures


# ============================================================================
# RUN EDGE SEARCH
# ============================================================================

def run_edge_search(
    records: Sequence[EdgeRecord],
) -> None:

    print()
    print("=" * 78)
    print("TOP-EDGE LAW SEARCH")
    print("=" * 78)

    # We only inspect X powers 0..3.
    # Higher X powers are a secondary question.
    for kind in (
        "C-top",
        "C-second",
        "D-top",
        "D-second",
    ):

        print()
        print(kind)

        for x_index in range(4):

            poly, degree, count = search_component_law(
                records,
                kind,
                x_index,
                MAX_BIV_DEG,
            )

            if poly is None:
                continue

            loo_ok, loo_failures = loo_test(
                records,
                kind,
                x_index,
                degree,
            )

            print(
                f"  X^{x_index}: "
                f"degree={degree} "
                f"points={count} "
                f"LOO={'PASS' if loo_ok else 'FAIL'} "
                f"failures={loo_failures}"
            )

            print(
                f"      law = {sp.factor(poly)}"
            )


# ============================================================================
# FORWARD HOLDOUT
# ============================================================================

def run_forward_test() -> None:

    print()
    print("=" * 78)
    print("FORWARD ELL HOLDOUT")
    print("=" * 78)

    train = build_dataset(MAX_ELL)

    future_all = build_dataset(FORWARD_ELL)

    future = [
        rec for rec in future_all
        if rec.ell == FORWARD_ELL
    ]

    print(
        f"training records = {len(train)}"
    )
    print(
        f"future ell={FORWARD_ELL} records = {len(future)}"
    )

    for kind in (
        "C-top",
        "C-second",
        "D-top",
        "D-second",
    ):

        for x_index in range(3):

            # First find a law on the training data.
            poly, degree, count = search_component_law(
                train,
                kind,
                x_index,
                MAX_BIV_DEG,
            )

            if poly is None:
                continue

            ok, failures = forward_holdout(
                train,
                future,
                kind,
                x_index,
                degree,
            )

            print(
                f"{kind} X^{x_index}: "
                f"degree={degree} "
                f"forward={'PASS' if ok else 'FAIL'} "
                f"failures={failures}"
            )

            if ok:
                print(
                    f"    law = {sp.factor(poly)}"
                )


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:

    start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 122")
    print("C/D TOP-EDGE LAW")
    print("TOP N-DEGREE COEFFICIENTS ONLY")
    print("EXACT SYMBOLIC ARITHMETIC")
    print("NO FULL TENSOR SEARCH")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("SAFE SYMPY DOMAINS")
    print("=" * 78)

    print()
    print("1. BUILDING EDGE DATASET")
    print("-------------------------")

    t0 = time.perf_counter()

    records = build_dataset(MAX_ELL)

    print(
        f"records = {len(records)}"
    )
    print(
        f"build time = {time.perf_counter() - t0:.6f}s"
    )

    print_degree_profile(records)

    print_edge_samples(records)

    t0 = time.perf_counter()

    run_edge_search(records)

    print()
    print(
        f"edge-search time = "
        f"{time.perf_counter() - t0:.6f}s"
    )

    print()
    print("2. FORWARD HOLDOUT")
    print("------------------")

    # To avoid rebuilding MAX_ELL data twice unnecessarily in the main run,
    # use a small explicit forward dataset here.
    future_records = build_dataset(FORWARD_ELL)
    future_records = [
        rec for rec in future_records
        if rec.ell == FORWARD_ELL
    ]

    print(
        f"future ell={FORWARD_ELL} records = "
        f"{len(future_records)}"
    )

    for kind in (
        "C-top",
        "C-second",
        "D-top",
        "D-second",
    ):

        for x_index in range(3):

            poly, degree, _ = search_component_law(
                records,
                kind,
                x_index,
                MAX_BIV_DEG,
            )

            if poly is None:
                continue

            ok, failures = forward_holdout(
                records,
                future_records,
                kind,
                x_index,
                degree,
            )

            print(
                f"{kind} X^{x_index}: "
                f"degree={degree} "
                f"forward={'PASS' if ok else 'FAIL'} "
                f"failures={failures}"
            )

            if ok:
                print(
                    f"    law = {sp.factor(poly)}"
                )

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    print()
    print("Experiment 122 deliberately ignores the bulk of the C/D tensor.")

    print()
    print("The primary object is the upper N-degree edge:")
    print()
    print("    C_top(X), C_second(X)")
    print("    D_top(X), D_second(X)")
    print()
    print("and whether their coefficients admit exact laws in")
    print()
    print("    r=(k-1)/2, m=(ell-1)/2.")

    print()
    print("A strong result is:")
    print()
    print("    exact low-degree law")
    print("    + leave-one-ell-out PASS")
    print("    + ell=19 forward PASS")
    print()
    print("That would justify deriving the top-edge kernel analytically")
    print("before touching the full coefficient tensor again.")

    print()
    print(
        f"total runtime = "
        f"{time.perf_counter() - start:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 122 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

