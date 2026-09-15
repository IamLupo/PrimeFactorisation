#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
==============================================================================
KAPPA EXPERIMENT 123
C/D SCALAR EDGE LAWS
TOP N-COEFFICIENTS AT X=0
SECOND N-LAYER AT X=0
EXACT LOW-DEGREE (r,m) LAW SEARCH
NO FULL TENSOR SEARCH
NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
SAFE SYMPY DOMAINS
==============================================================================

Definitions
-----------

For

    Q(N,S) = C(N,X) + (2S-1) D(N,X),
    X = S(S-1),

define

    L_C(r,m) = leading coefficient in N of C(N,0)
    L_D(r,m) = leading coefficient in N of D(N,0)

and, when present,

    M_C(r,m) = second coefficient in N of C(N,0)
    M_D(r,m) = second coefficient in N of D(N,0).

The experiment searches exact polynomial laws in

    r = (k-1)/2
    m = (ell-1)/2

and validates them by:

    * exact fit on all training weights
    * leave-one-ell-out
    * forward ell holdout

No search is performed on zero X-layers.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import List, Optional, Sequence, Tuple

import sympy as sp


# ============================================================================
# PARAMETERS
# ============================================================================

TRAIN_MAX_ELL = 23
FORWARD_ELL = 25
MAX_BIV_DEG = 6


N, S, X = sp.symbols("N S X")
r_sym, m_sym = sp.symbols("r m")


# ============================================================================
# RECORD
# ============================================================================

@dataclass(frozen=True)
class Record:
    k: int
    ell: int
    r: int
    m: int

    c_degree: int
    c_top: sp.Rational
    c_second: sp.Rational

    d_degree: Optional[int]
    d_top: sp.Rational
    d_second: sp.Rational


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
# RAW DETECTOR
# ============================================================================

def raw_detector_numerator(k: int, ell: int) -> sp.Expr:

    P = power_sum_table(max(k, ell))

    F = sp.Integer(0)

    for t in range(ell + 1):
        F += sp.binomial(ell, t) * symmetric_monomial(k, t, P)

    for t in range(k + 1):
        F -= sp.binomial(k, t) * symmetric_monomial(ell, t, P)

    return sp.expand(F)


# ============================================================================
# EXACT QUOTIENT
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

def reduce_quadratic(expr: sp.Expr) -> sp.Expr:

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

    return sp.expand(
        poly.rem(modulus).as_expr()
    )


def branch_decompose(Q: sp.Expr) -> Tuple[sp.Expr, sp.Expr]:

    R = reduce_quadratic(Q)

    poly = sp.Poly(
        R,
        S,
        domain=sp.QQ.poly_ring(N, X),
    )

    a0 = sp.expand(
        poly.coeff_monomial(S**0)
    )

    a1 = sp.expand(
        poly.coeff_monomial(S**1)
    )

    C = sp.expand(a0 + a1 / 2)
    D = sp.expand(a1 / 2)

    return C, D


# ============================================================================
# SCALAR N-EDGE
# ============================================================================

def scalar_n_layers(
    expr: sp.Expr,
) -> Tuple[int, sp.Rational, sp.Rational]:

    """
    Evaluate at X=0 first, then extract the top two N coefficients.

    Returns:

        degree_N,
        leading coefficient,
        second coefficient

    This directly avoids all irrelevant X^j=0 layers.
    """

    e0 = sp.expand(expr.subs(X, 0))

    if e0 == 0:
        return -1, sp.Rational(0), sp.Rational(0)

    poly = sp.Poly(
        e0,
        N,
        domain=sp.QQ,
    )

    degree = int(poly.degree())

    top = sp.Rational(
        poly.coeff_monomial(N**degree)
    )

    if degree >= 1:
        second = sp.Rational(
            poly.coeff_monomial(N**(degree - 1))
        )
    else:
        second = sp.Rational(0)

    return degree, top, second


# ============================================================================
# DATASET
# ============================================================================

def build_dataset(max_ell: int) -> List[Record]:

    rows: List[Record] = []

    for ell in range(3, max_ell + 1, 2):

        for k in range(1, ell, 2):

            Q = quotient_polynomial(k, ell)

            C, D = branch_decompose(Q)

            cdeg, ctop, csecond = scalar_n_layers(C)

            ddeg_raw, dtop, dsecond = scalar_n_layers(D)

            ddeg = (
                None
                if ddeg_raw < 0
                else ddeg_raw
            )

            rows.append(
                Record(
                    k=k,
                    ell=ell,
                    r=(k - 1) // 2,
                    m=(ell - 1) // 2,

                    c_degree=cdeg,
                    c_top=ctop,
                    c_second=csecond,

                    d_degree=ddeg,
                    d_top=dtop,
                    d_second=dsecond,
                )
            )

    return rows


# ============================================================================
# POLYNOMIAL BASIS
# ============================================================================

def basis_monomials(
    degree: int,
) -> List[Tuple[int, int]]:

    out = []

    for i in range(degree + 1):
        for j in range(degree + 1 - i):
            out.append((i, j))

    return out


def exact_fit(
    points: Sequence[Tuple[int, int]],
    values: Sequence[sp.Rational],
    degree: int,
) -> Optional[sp.Expr]:

    basis = basis_monomials(degree)

    if len(points) < len(basis):
        return None

    A = sp.Matrix([
        [
            sp.Integer(r) ** i *
            sp.Integer(m) ** j
            for i, j in basis
        ]
        for r, m in points
    ])

    b = sp.Matrix(values)

    try:
        solution_set = sp.linsolve((A, b))
    except Exception:
        return None

    sols = list(solution_set)

    if len(sols) != 1:
        return None

    sol = sols[0]

    if any(v.free_symbols for v in sol):
        return None

    result = sp.Integer(0)

    for coeff, (i, j) in zip(sol, basis):
        result += coeff * r_sym**i * m_sym**j

    return sp.expand(result)


def evaluate(
    poly: sp.Expr,
    r: int,
    m: int,
) -> sp.Rational:

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
# LAW SEARCH
# ============================================================================

def field_value(
    rec: Record,
    field: str,
) -> sp.Rational:

    value = getattr(rec, field)

    if value is None:
        raise ValueError(
            f"Field {field} missing."
        )

    return sp.Rational(value)


def search_law(
    rows: Sequence[Record],
    field: str,
    max_degree: int,
) -> Tuple[
    Optional[sp.Expr],
    Optional[int],
    int,
]:

    points = [
        (rec.r, rec.m)
        for rec in rows
    ]

    values = [
        field_value(rec, field)
        for rec in rows
    ]

    if len(rows) < 6:
        return None, None, len(rows)

    for degree in range(max_degree + 1):

        poly = exact_fit(
            points,
            values,
            degree,
        )

        if poly is None:
            continue

        failures = 0

        for rec in rows:

            actual = field_value(
                rec,
                field,
            )

            predicted = evaluate(
                poly,
                rec.r,
                rec.m,
            )

            if predicted != actual:
                failures += 1

        if failures == 0:
            return poly, degree, len(rows)

    return None, None, len(rows)


# ============================================================================
# LOO
# ============================================================================

def loo(
    rows: Sequence[Record],
    field: str,
    degree: int,
) -> Tuple[int, int]:

    total = 0
    failures = 0

    unique_m = sorted(
        set(rec.m for rec in rows)
    )

    for held_m in unique_m:

        train = [
            rec
            for rec in rows
            if rec.m != held_m
        ]

        test = [
            rec
            for rec in rows
            if rec.m == held_m
        ]

        if not train or not test:
            continue

        poly = exact_fit(
            [(rec.r, rec.m) for rec in train],
            [
                field_value(rec, field)
                for rec in train
            ],
            degree,
        )

        if poly is None:
            failures += len(test)
            total += len(test)
            continue

        for rec in test:

            actual = field_value(
                rec,
                field,
            )

            predicted = evaluate(
                poly,
                rec.r,
                rec.m,
            )

            total += 1

            if predicted != actual:
                failures += 1

    return total, failures


# ============================================================================
# FORWARD HOLDOUT
# ============================================================================

def forward_test(
    train: Sequence[Record],
    future: Sequence[Record],
    field: str,
    degree: int,
) -> Tuple[int, int]:

    poly = exact_fit(
        [(rec.r, rec.m) for rec in train],
        [
            field_value(rec, field)
            for rec in train
        ],
        degree,
    )

    if poly is None:
        return len(future), len(future)

    total = 0
    failures = 0

    for rec in future:

        actual = field_value(
            rec,
            field,
        )

        predicted = evaluate(
            poly,
            rec.r,
            rec.m,
        )

        total += 1

        if predicted != actual:
            failures += 1

    return total, failures


# ============================================================================
# PRINT EDGE TABLE
# ============================================================================

def print_edge_table(
    rows: Sequence[Record],
) -> None:

    print()
    print("=" * 78)
    print("SCALAR EDGE TABLE")
    print("=" * 78)

    for rec in rows:

        print(
            f"(k,ell)=({rec.k},{rec.ell}) "
            f"r={rec.r} m={rec.m} | "
            f"C deg={rec.c_degree} "
            f"L={rec.c_top} "
            f"M={rec.c_second} | "
            f"D deg={rec.d_degree} "
            f"L={rec.d_top} "
            f"M={rec.d_second}"
        )


# ============================================================================
# RUN SEARCH
# ============================================================================

def run_search(
    rows: Sequence[Record],
) -> None:

    targets = [
        ("C-top", "c_top"),
        ("C-second", "c_second"),
        ("D-top", "d_top"),
        ("D-second", "d_second"),
    ]

    print()
    print("=" * 78)
    print("EXACT SCALAR EDGE LAW SEARCH")
    print("=" * 78)

    for label, field in targets:

        print()
        print(label)

        poly, degree, count = search_law(
            rows,
            field,
            MAX_BIV_DEG,
        )

        if poly is None:
            print(
                f"  no exact law of degree <= {MAX_BIV_DEG}"
            )
            continue

        total, failures = loo(
            rows,
            field,
            degree,
        )

        print(
            f"  degree = {degree}"
        )
        print(
            f"  points = {count}"
        )
        print(
            f"  LOO = "
            f"{'PASS' if failures == 0 else 'FAIL'} "
            f"({failures}/{total})"
        )
        print(
            f"  law = {sp.factor(poly)}"
        )


# ============================================================================
# FORWARD
# ============================================================================

def run_forward(
    train: Sequence[Record],
    future: Sequence[Record],
) -> None:

    print()
    print("=" * 78)
    print("FORWARD ELL HOLDOUT")
    print("=" * 78)

    targets = [
        ("C-top", "c_top"),
        ("C-second", "c_second"),
        ("D-top", "d_top"),
        ("D-second", "d_second"),
    ]

    for label, field in targets:

        poly, degree, _ = search_law(
            train,
            field,
            MAX_BIV_DEG,
        )

        if poly is None:
            print(
                f"{label}: "
                f"no training law"
            )
            continue

        total, failures = forward_test(
            train,
            future,
            field,
            degree,
        )

        print(
            f"{label}: "
            f"degree={degree} "
            f"forward="
            f"{'PASS' if failures == 0 else 'FAIL'} "
            f"({failures}/{total})"
        )

        if failures == 0:
            print(
                f"  law = {sp.factor(poly)}"
            )


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:

    t_start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 123")
    print("C/D SCALAR EDGE LAWS")
    print("TOP N-COEFFICIENTS AT X=0")
    print("SECOND N-LAYER AT X=0")
    print("EXACT LOW-DEGREE (r,m) LAW SEARCH")
    print("NO FULL TENSOR SEARCH")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("SAFE SYMPY DOMAINS")
    print("=" * 78)

    print()
    print("1. BUILDING DATASET")
    print("--------------------")

    t0 = time.perf_counter()

    train = build_dataset(
        TRAIN_MAX_ELL
    )

    print(
        f"records = {len(train)}"
    )

    print(
        f"build time = "
        f"{time.perf_counter() - t0:.6f}s"
    )

    print_edge_table(train)

    run_search(train)

    print()
    print("2. FORWARD HOLDOUT")
    print("------------------")

    t0 = time.perf_counter()

    future_all = build_dataset(
        FORWARD_ELL
    )

    future = [
        rec
        for rec in future_all
        if rec.ell == FORWARD_ELL
    ]

    print(
        f"future ell = {FORWARD_ELL}"
    )
    print(
        f"future records = {len(future)}"
    )

    run_forward(
        train,
        future,
    )

    print(
        f"forward build time = "
        f"{time.perf_counter() - t0:.6f}s"
    )

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    print()
    print("Experiment 123 intentionally reduces the problem to four scalar")
    print("functions:")
    print()
    print("    L_C(r,m)")
    print("    M_C(r,m)")
    print("    L_D(r,m)")
    print("    M_D(r,m)")
    print()
    print("where L is the highest-N coefficient and M the next one at X=0.")
    print()
    print("This directly tests whether the complicated C/D tensor has a")
    print("simple combinatorial upper edge.")
    print()
    print("A useful positive result is:")
    print()
    print("    exact polynomial law")
    print("    + LOO PASS")
    print("    + forward ell PASS.")
    print()
    print("A failure is also useful: it means the upper edge is not a simple")
    print("low-degree polynomial in (r,m), and we should test binomial/factorial")
    print("normalizations rather than larger polynomial degrees.")

    print()
    print(
        f"total runtime = "
        f"{time.perf_counter() - t_start:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 123 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

