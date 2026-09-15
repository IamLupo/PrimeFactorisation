#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 126
DIRECT POST-BOUNDARY LEADING-LAW CERTIFICATE
EXACT FINITE BINOMIAL KERNEL -> Q(N,S) -> ONE NEWTON COEFFICIENT
NO FULL C/D TENSOR
NO FULL NEWTON TENSOR
NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
SAFE SYMPY DOMAINS

TARGET:

    r = (k-1)/2
    m = (ell-1)/2
    r = k + d,  d > 0
    j = ell - 1 - r

Experiment 112R suggested the exact leading coefficient law

    L_d(a,b) = 2*(-1)^((d-1)/2) * (a+b+1)                    d odd

    L_d(a,b) = (-1)^(d/2) * (a+b+1)*(2*a - 2*b + d + 1)       d even

with

    a = r
    b = m

This experiment does NOT construct the C/D tensor.

It constructs Q directly from the finite binomial kernel, extracts only
the required Newton coefficient c_j(N), determines its ACTUAL leading
coefficient, and tests the law exactly.

It also records cancellation exceptions separately.

==============================================================================
"""

from __future__ import annotations

import sys
import time
from dataclasses import dataclass
from typing import Dict, List, Tuple

import sympy as sp


# -----------------------------------------------------------------------------
# SYMBOLS / DOMAINS
# -----------------------------------------------------------------------------

S = sp.Symbol("S")
N = sp.Symbol("N")

QQN = sp.QQ.frac_field(N)

X = sp.Symbol("X")


# -----------------------------------------------------------------------------
# CONFIGURATION
# -----------------------------------------------------------------------------

TRAIN_MAX_ELL = 23
FORWARD_ELL_VALUES = (25, 27)

D_VALUES = tuple(range(1, 11))

# Keep output compact.
PRINT_WITNESSES = True


# -----------------------------------------------------------------------------
# BASIC POWER-SUM KERNEL
#
# P_n = p^n + q^n
#
# P_0 = 2
# P_1 = S
# P_n = S P_{n-1} - N P_{n-2}
# -----------------------------------------------------------------------------

_POWER_SUM_CACHE: Dict[int, sp.Expr] = {
    0: sp.Integer(2),
    1: S,
}


def power_sum(n: int) -> sp.Expr:
    if n < 0:
        raise ValueError("negative power-sum index")

    if n in _POWER_SUM_CACHE:
        return _POWER_SUM_CACHE[n]

    for i in range(2, n + 1):
        if i not in _POWER_SUM_CACHE:
            _POWER_SUM_CACHE[i] = sp.expand(
                S * _POWER_SUM_CACHE[i - 1]
                - N * _POWER_SUM_CACHE[i - 2]
            )

    return _POWER_SUM_CACHE[n]


# -----------------------------------------------------------------------------
# SYMMETRIC MONOMIAL PAIR
#
# p^a q^b + q^a p^b
#
# = N^min(a,b) * (p^|a-b| + q^|a-b|)
# -----------------------------------------------------------------------------

def symmetric_monomial_pair(a: int, b: int) -> sp.Expr:
    if a < 0 or b < 0:
        raise ValueError("negative monomial exponent")

    t = min(a, b)
    delta = abs(a - b)

    return sp.expand(
        N**t * power_sum(delta)
    )


# -----------------------------------------------------------------------------
# DIRECT FINITE BINOMIAL KERNEL
#
# F_(k,ell)
# =
# p^k(1+q)^ell + q^k(1+p)^ell
# - p^ell(1+q)^k - q^ell(1+p)^k
#
# expressed directly in S,N.
# -----------------------------------------------------------------------------

def raw_kernel(k: int, ell: int) -> sp.Expr:
    first = sp.Integer(0)
    second = sp.Integer(0)

    for t in range(ell + 1):
        first += sp.binomial(ell, t) * symmetric_monomial_pair(
            k, t
        )

    for t in range(k + 1):
        second += sp.binomial(k, t) * symmetric_monomial_pair(
            ell, t
        )

    return sp.expand(first - second)


# -----------------------------------------------------------------------------
# EXACT QUOTIENT BY S+1
# -----------------------------------------------------------------------------

def quotient_polynomial(k: int, ell: int) -> sp.Expr:
    F = raw_kernel(k, ell)

    num = sp.Poly(F, S, domain=QQN)
    den = sp.Poly(S + 1, S, domain=QQN)

    Qpoly, rem = sp.div(num, den)

    rem_expr = sp.expand(rem.as_expr())

    if rem_expr != 0:
        raise ArithmeticError(
            f"Nonzero remainder for ({k},{ell}): {rem_expr}"
        )

    return sp.expand(Qpoly.as_expr())


# -----------------------------------------------------------------------------
# NEWTON / POWER-SUM BASIS
#
# P_j = p^j+q^j
#
# Q(S,N) = sum_j c_j(N) P_j(S,N)
#
# We only extract one target coefficient c_target.
# -----------------------------------------------------------------------------

_BASIS_CACHE: Dict[int, sp.Expr] = {}


def newton_basis(j: int) -> sp.Expr:
    if j in _BASIS_CACHE:
        return _BASIS_CACHE[j]

    P = power_sum(j)
    _BASIS_CACHE[j] = sp.expand(P)
    return _BASIS_CACHE[j]


def extract_newton_coefficient(
    Q: sp.Expr,
    target_j: int,
) -> sp.Expr:
    """
    Triangular extraction from the top S-degree down to target_j.

    P_j is monic in S^j, so the current top S coefficient is exactly
    the Newton coefficient at that degree.
    """

    poly = sp.Poly(Q, S, domain=QQN)
    degree = poly.degree()

    if degree < target_j:
        raise ArithmeticError(
            f"Q degree {degree} < target j {target_j}"
        )

    residual = sp.expand(Q)
    coeffs: Dict[int, sp.Expr] = {}

    for j in range(degree, target_j - 1, -1):
        poly_res = sp.Poly(residual, S, domain=QQN)
        c_j = sp.cancel(poly_res.coeff_monomial(S**j))

        coeffs[j] = sp.expand(c_j)

        residual = sp.expand(
            residual - c_j * newton_basis(j)
        )

    # We deliberately do NOT require the remainder below target_j to vanish.
    # Only the requested coefficient is being extracted.

    return sp.expand(coeffs[target_j])


# -----------------------------------------------------------------------------
# DIRECT COEFFICIENT OBJECT
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class Record:
    k: int
    ell: int
    d: int
    r: int
    m: int
    j: int
    coefficient: sp.Expr
    actual_degree: int | None
    actual_lead: sp.Expr
    nominal_degree: int
    expected_lead: sp.Expr
    cancellation: bool


# -----------------------------------------------------------------------------
# POST-BOUNDARY PARAMETERS
# -----------------------------------------------------------------------------

def post_boundary_parameters(
    k: int,
    ell: int,
    d: int,
) -> Tuple[int, int, int, int, int]:
    """
    r = k+d
    m = (ell-1)/2
    j = ell-1-r

    Returns:
        r, m, j, nominal_degree, a
    """

    r = k + d
    m = (ell - 1) // 2
    j = ell - 1 - r

    if k % 2 != 1 or ell % 2 != 1:
        raise ValueError("k and ell must be odd")

    if r > ell - 1:
        raise ValueError("post-boundary index exceeds Newton range")

    # Experiment 112R degree profile:
    #
    # d=1,2 -> degree k
    # d=3,4 -> degree k+1
    # d=5,6 -> degree k+2
    #
    # Hence:
    #
    # nominal degree = k + floor((d-1)/2)
    nominal_degree = k + (d - 1) // 2

    return r, m, j, nominal_degree, r


# -----------------------------------------------------------------------------
# CLOSED LEADING-LAW CONJECTURE
# -----------------------------------------------------------------------------

def expected_lead(k: int, ell: int, d: int) -> sp.Expr:
    """
    Experiment-112R conjectured law.

    a = (k-1)/2
    b = (ell-1)/2
    """

    a = sp.Integer((k - 1) // 2)
    b = sp.Integer((ell - 1) // 2)

    common = a + b + 1

    if d % 2 == 1:
        sign = (-1) ** ((d - 1) // 2)
        return sp.expand(
            2 * sign * common
        )

    sign = (-1) ** (d // 2)

    return sp.expand(
        sign * common * (
            2 * a - 2 * b + d + 1
        )
    )


# -----------------------------------------------------------------------------
# ACTUAL DEGREE / LEAD
# -----------------------------------------------------------------------------

def polynomial_degree_and_lead(
    expr: sp.Expr,
) -> Tuple[int | None, sp.Expr]:

    poly = sp.Poly(sp.expand(expr), N, domain=sp.QQ)

    if poly.is_zero:
        return None, sp.Integer(0)

    deg = int(poly.degree())
    lead = sp.expand(poly.LC())

    return deg, lead


# -----------------------------------------------------------------------------
# RECORD
# -----------------------------------------------------------------------------

def coefficient_record(
    k: int,
    ell: int,
    d: int,
) -> Record:

    r, m, j, nominal_degree, _ = post_boundary_parameters(
        k, ell, d
    )

    Q = quotient_polynomial(k, ell)

    coefficient = extract_newton_coefficient(Q, j)

    actual_degree, actual_lead = polynomial_degree_and_lead(
        coefficient
    )

    exp_lead = expected_lead(k, ell, d)

    cancellation = (
        actual_degree != nominal_degree
        or sp.expand(actual_lead - exp_lead) != 0
    )

    return Record(
        k=k,
        ell=ell,
        d=d,
        r=r,
        m=m,
        j=j,
        coefficient=coefficient,
        actual_degree=actual_degree,
        actual_lead=actual_lead,
        nominal_degree=nominal_degree,
        expected_lead=exp_lead,
        cancellation=cancellation,
    )


# -----------------------------------------------------------------------------
# TRAINING / HOLDOUT GRID
# -----------------------------------------------------------------------------

def available_weights(max_ell: int) -> List[Tuple[int, int, int]]:
    rows: List[Tuple[int, int, int]] = []

    for ell in range(3, max_ell + 1, 2):
        for k in range(1, ell, 2):
            for d in D_VALUES:
                r = k + d
                j = ell - 1 - r

                if r <= ell - 1 and j >= 0:
                    rows.append((k, ell, d))

    return rows


# -----------------------------------------------------------------------------
# BUILD DATASET
# -----------------------------------------------------------------------------

def build_dataset(
    max_ell: int,
) -> List[Record]:

    rows: List[Record] = []

    for k, ell, d in available_weights(max_ell):
        rows.append(
            coefficient_record(k, ell, d)
        )

    return rows


# -----------------------------------------------------------------------------
# EXACT VALIDATION OF THE LEADING LAW
# -----------------------------------------------------------------------------

def validate_leading_law(
    rows: List[Record],
    label: str,
) -> int:

    failures = 0
    nondegenerate = 0
    cancellation_count = 0

    print()
    print("=" * 78)
    print(f"{label}")
    print("=" * 78)

    for rec in rows:
        expected = sp.expand(rec.expected_lead)

        if rec.actual_degree != rec.nominal_degree:
            cancellation_count += 1
            continue

        nondegenerate += 1

        if sp.expand(rec.actual_lead - expected) != 0:
            failures += 1

            if PRINT_WITNESSES:
                print(
                    f"FAIL ({rec.k},{rec.ell}), d={rec.d}: "
                    f"actual_deg={rec.actual_degree} "
                    f"lead={rec.actual_lead} "
                    f"expected={expected}"
                )

    print()
    print(f"rows                 = {len(rows)}")
    print(f"nondegenerate rows   = {nondegenerate}")
    print(f"cancellation rows    = {cancellation_count}")
    print(f"law failures         = {failures}")

    return failures


# -----------------------------------------------------------------------------
# CANCELLATION CLASSIFICATION
# -----------------------------------------------------------------------------

def print_cancellations(
    rows: List[Record],
) -> None:

    cancels = [
        rec for rec in rows
        if rec.actual_degree != rec.nominal_degree
    ]

    print()
    print("=" * 78)
    print("CANCELLATION LOCUS")
    print("=" * 78)

    if not cancels:
        print("none")
        return

    for rec in cancels:
        print(
            f"({rec.k},{rec.ell}) d={rec.d}: "
            f"actual_degree={rec.actual_degree}, "
            f"nominal_degree={rec.nominal_degree}, "
            f"lead={rec.actual_lead}, "
            f"expected={rec.expected_lead}"
        )


# -----------------------------------------------------------------------------
# SAMPLE CERTIFICATES
# -----------------------------------------------------------------------------

def print_sample_certificates(
    rows: List[Record],
) -> None:

    print()
    print("=" * 78)
    print("EXACT CERTIFICATE SAMPLES")
    print("=" * 78)

    wanted = [
        (1, 5, 1),
        (3, 7, 1),
        (5, 9, 1),
        (1, 7, 2),
        (3, 9, 2),
        (5, 11, 2),
        (1, 9, 3),
        (3, 11, 4),
        (5, 13, 5),
        (1, 13, 6),
        (3, 17, 8),
        (1, 23, 10),
    ]

    index = {
        (rec.k, rec.ell, rec.d): rec
        for rec in rows
    }

    for key in wanted:
        rec = index.get(key)

        if rec is None:
            continue

        print(
            f"({rec.k},{rec.ell}) "
            f"d={rec.d} r={rec.r} j={rec.j}"
        )
        print(
            f"  c_j(N)       = {sp.factor(rec.coefficient)}"
        )
        print(
            f"  degree        = {rec.actual_degree}"
        )
        print(
            f"  actual lead   = {rec.actual_lead}"
        )
        print(
            f"  expected lead = {rec.expected_lead}"
        )
        print(
            f"  match         = "
            f"{sp.expand(rec.actual_lead - rec.expected_lead) == 0}"
        )
        print()


# -----------------------------------------------------------------------------
# FORWARD HOLDOUT
# -----------------------------------------------------------------------------

def forward_holdout(
    train_max_ell: int,
    future_ells: Tuple[int, ...],
) -> None:

    print()
    print("=" * 78)
    print("FORWARD ELL HOLDOUT")
    print("=" * 78)

    total = 0
    failures = 0
    cancellations = 0

    for ell in future_ells:
        rows = build_dataset(ell)

        rows = [
            rec for rec in rows
            if rec.ell == ell
        ]

        ell_fail = 0
        ell_cancel = 0

        for rec in rows:
            total += 1

            if rec.actual_degree != rec.nominal_degree:
                cancellations += 1
                ell_cancel += 1
                continue

            if sp.expand(rec.actual_lead - rec.expected_lead) != 0:
                failures += 1
                ell_fail += 1

        print(
            f"ell={ell}: rows={len(rows)} "
            f"law_failures={ell_fail} "
            f"cancellations={ell_cancel}"
        )

    print()
    print(f"forward rows       = {total}")
    print(f"forward failures   = {failures}")
    print(f"forward cancellations = {cancellations}")

    if failures:
        print("STATUS = FAIL")
    else:
        print("STATUS = PASS")


# -----------------------------------------------------------------------------
# RATIO CERTIFICATE
#
# Verify the formula's induced ratios in b for fixed k,d where both
# numerator and denominator are nonzero.
# This is not a fit; it is a consistency certificate.
# -----------------------------------------------------------------------------

def ratio_certificate(rows: List[Record]) -> None:

    print()
    print("=" * 78)
    print("RATIO CERTIFICATE")
    print("=" * 78)

    by_kd: Dict[Tuple[int, int], List[Record]] = {}

    for rec in rows:
        by_kd.setdefault(
            (rec.k, rec.d),
            []
        ).append(rec)

    checks = 0
    failures = 0

    for (k, d), group in sorted(by_kd.items()):

        group = sorted(
            group,
            key=lambda rec: rec.m
        )

        for left, right in zip(group, group[1:]):

            # Only compare nondegenerate actual-leading terms.
            if (
                left.actual_degree != left.nominal_degree
                or right.actual_degree != right.nominal_degree
            ):
                continue

            L1 = sp.expand(left.actual_lead)
            L2 = sp.expand(right.actual_lead)

            if L1 == 0:
                continue

            observed = sp.cancel(L2 / L1)

            expected_ratio = sp.cancel(
                right.expected_lead /
                left.expected_lead
            )

            checks += 1

            if sp.simplify(observed - expected_ratio) != 0:
                failures += 1

                if PRINT_WITNESSES:
                    print(
                        f"FAIL k={k} d={d}: "
                        f"m={left.m}->{right.m}, "
                        f"observed={observed}, "
                        f"expected={expected_ratio}"
                    )

    print(f"ratio checks    = {checks}")
    print(f"ratio failures  = {failures}")

    if failures == 0:
        print("ratio status    = PASS")
    else:
        print("ratio status    = FAIL")


# -----------------------------------------------------------------------------
# SUMMARY BY d
# -----------------------------------------------------------------------------

def summary_by_d(
    rows: List[Record],
) -> None:

    print()
    print("=" * 78)
    print("SUMMARY BY d")
    print("=" * 78)

    for d in D_VALUES:

        group = [
            rec for rec in rows
            if rec.d == d
        ]

        if not group:
            continue

        nondeg = [
            rec for rec in group
            if rec.actual_degree == rec.nominal_degree
        ]

        fails = [
            rec for rec in nondeg
            if sp.expand(
                rec.actual_lead - rec.expected_lead
            ) != 0
        ]

        print(
            f"d={d:2d}: rows={len(group):3d} "
            f"nondeg={len(nondeg):3d} "
            f"cancel={len(group)-len(nondeg):3d} "
            f"fail={len(fails):3d}"
        )


# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------

def main() -> None:

    t0 = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 126")
    print("DIRECT POST-BOUNDARY LEADING-LAW CERTIFICATE")
    print("EXACT FINITE BINOMIAL KERNEL -> ONE NEWTON COEFFICIENT")
    print("NO FULL C/D TENSOR")
    print("NO FULL NEWTON TENSOR")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("SAFE SYMPY DOMAINS")
    print("=" * 78)

    print()
    print("1. TRAINING DATASET")
    print("-" * 78)

    t1 = time.perf_counter()

    train_rows = build_dataset(TRAIN_MAX_ELL)

    build_time = time.perf_counter() - t1

    print(f"training max ell = {TRAIN_MAX_ELL}")
    print(f"training rows    = {len(train_rows)}")
    print(f"build time       = {build_time:.6f}s")

    # -------------------------------------------------------------------------
    # Exact law test
    # -------------------------------------------------------------------------

    failures = validate_leading_law(
        train_rows,
        "2. EXACT POST-BOUNDARY LEADING LAW"
    )

    print_cancellations(train_rows)
    summary_by_d(train_rows)
    print_sample_certificates(train_rows)
    ratio_certificate(train_rows)

    # -------------------------------------------------------------------------
    # Forward holdout
    # -------------------------------------------------------------------------

    t2 = time.perf_counter()

    forward_holdout(
        TRAIN_MAX_ELL,
        FORWARD_ELL_VALUES,
    )

    forward_time = time.perf_counter() - t2

    print()
    print(f"forward build/test time = {forward_time:.6f}s")

    # -------------------------------------------------------------------------
    # Final diagnostic
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    if failures == 0:
        print(
            """
The Experiment-112R post-boundary leading law survives the
direct finite-kernel/Newton-coefficient calculation on the
training grid, apart from explicitly classified cancellation
cases.

The important point is that this script does not construct the
C/D tensor and does not fit a polynomial law.

It verifies the leading coefficient directly from the detector
kernel.

Next mathematical target:

    prove the closed form symbolically from the finite binomial sum.

The ratio certificate provides an additional hypergeometric-style
consistency check.
"""
        )
    else:
        print(
            """
The proposed leading law does NOT survive the direct kernel
certificate.

Inspect the printed counterexamples before doing any further
symbolic expansion.
"""
        )

    total_time = time.perf_counter() - t0

    print()
    print(f"total runtime = {total_time:.6f}s")
    print("=" * 78)
    print("EXPERIMENT 126 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

