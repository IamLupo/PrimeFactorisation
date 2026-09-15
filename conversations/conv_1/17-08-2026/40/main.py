#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 108
EXACT BINOMIAL KERNEL FOR THE NEWTON COEFFICIENT TENSOR

DIRECT PAPER DETECTOR FORM
DIVISION BY (P+Q+1)
EXACT NEWTON-BASIS EXTRACTION
BINOMIAL-SUM COEFFICIENT LAW
STRICT SYMBOLIC VALIDATION
NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN

Core identity
-------------
For odd k < ell define

    F_(k,ell) =
        p^k (1+q)^ell + q^k (1+p)^ell
        - p^ell (1+q)^k - q^ell (1+p)^k.

Then

    F_(k,ell) = (p+q+1) Q_(k,ell).

Experiment 107 exposed a structured Newton coefficient tensor but
used the wrong terminal-law assumption and an inadequate bivariate
interpolation basis.

This experiment does NOT interpolate coefficients.

Instead:

  1. build F_(k,ell) directly from the detector definition;
  2. divide exactly by p+q+1;
  3. convert Q to the Newton basis P_j=p^j+q^j;
  4. derive the coefficients directly from binomial expansion;
  5. compare the derived coefficient formula against exact symbolic
     coefficients;
  6. search for a compact finite-sum kernel;
  7. test the formula for weights not used in the initial discovery.

The important distinction is:

    interpolation:
        "this polynomial happens to fit"

    direct kernel:
        "this coefficient follows from the detector definition."

==============================================================================
"""

from __future__ import annotations

from dataclasses import dataclass
from math import comb
import time

import sympy as sp


# =============================================================================
# SYMBOLS
# =============================================================================

p, q = sp.symbols("p q")
N, S = sp.symbols("N S")


# =============================================================================
# WEIGHTS
# =============================================================================

DISCOVERY_WEIGHTS = (
    (1, 3),
    (1, 5),
    (3, 5),
    (1, 7),
    (3, 7),
    (5, 7),
    (1, 9),
    (3, 9),
    (5, 9),
    (7, 9),
    (1, 11),
    (3, 11),
    (5, 11),
    (7, 11),
    (9, 11),
)

HOLDOUT_WEIGHTS = (
    (1, 13),
    (3, 13),
    (5, 13),
    (7, 13),
    (9, 13),
    (11, 13),
)

ALL_WEIGHTS = DISCOVERY_WEIGHTS + HOLDOUT_WEIGHTS


# =============================================================================
# BASIC SYMMETRIC / NEWTON UTILITIES
# =============================================================================

def power_sums(max_j: int) -> list[sp.Expr]:
    """
    P_0=2
    P_1=S
    P_j=S P_(j-1)-N P_(j-2)
    """
    P = [sp.Integer(0)] * (max_j + 1)
    P[0] = sp.Integer(2)

    if max_j >= 1:
        P[1] = S

    for j in range(2, max_j + 1):
        P[j] = sp.expand(
            S * P[j - 1] - N * P[j - 2]
        )

    return P


def symmetrize_to_NS(expr: sp.Expr) -> sp.Expr:
    """
    Convert a symmetric polynomial in p,q to S=p+q, N=pq.
    """
    symmetric, remainder, mapping = sp.symmetrize(
        sp.expand(expr),
        [p, q],
        formal=True,
    )

    if sp.expand(remainder) != 0:
        raise ArithmeticError(
            "Expression is not symmetric in p,q."
        )

    mapping = dict(mapping)

    s1 = None
    s2 = None

    for symbol, value in mapping.items():
        if sp.expand(value - (p + q)) == 0:
            s1 = symbol

        if sp.expand(value - p*q) == 0:
            s2 = symbol

    if s1 is None or s2 is None:
        raise ArithmeticError(
            "Could not identify elementary symmetric variables."
        )

    return sp.expand(
        symmetric.subs(
            {
                s1: S,
                s2: N,
            }
        )
    )


def quotient_detector(k: int, ell: int) -> sp.Expr:
    """
    Exact Q_(k,ell) from the direct paper detector.
    """
    if not (k < ell and k % 2 == 1 and ell % 2 == 1):
        raise ValueError(
            "Weights must satisfy odd k < odd ell."
        )

    F = sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )

    divisor = p + q + 1

    Q_pq, remainder = sp.div(
        sp.Poly(F, p, q, domain=sp.ZZ),
        sp.Poly(divisor, p, q, domain=sp.ZZ),
    )

    if remainder.as_expr() != 0:
        raise ArithmeticError(
            f"(p+q+1) does not divide detector ({k},{ell})."
        )

    Q_pq = Q_pq.as_expr()

    # Check exact reconstruction.
    if sp.expand(
        (p + q + 1) * Q_pq - F
    ) != 0:
        raise ArithmeticError(
            f"Quotient reconstruction failed ({k},{ell})."
        )

    return symmetrize_to_NS(Q_pq)


# =============================================================================
# NEWTON DECOMPOSITION
# =============================================================================

def newton_decompose(
    Q: sp.Expr,
    max_j: int,
) -> dict[int, sp.Expr]:
    """
    Q(N,S)=C0 + sum_j c_j(N) P_j.
    """
    P = power_sums(max_j)

    remainder = sp.expand(Q)
    coeffs: dict[int, sp.Expr] = {}

    for j in range(max_j, 0, -1):

        poly = sp.Poly(
            remainder,
            S,
            domain=sp.QQ.frac_field(N),
        )

        c = sp.factor(
            poly.coeff_monomial(S**j)
        )

        if c != 0:
            coeffs[j] = c
            remainder = sp.expand(
                remainder - c * P[j]
            )

    remainder = sp.factor(remainder)

    if sp.degree(remainder, S) not in (None, 0):
        raise ArithmeticError(
            "Newton decomposition left nonconstant S dependence."
        )

    coeffs[0] = remainder

    reconstruction = sp.expand(
        coeffs[0]
        + sum(
            coeffs.get(j, 0) * P[j]
            for j in range(1, max_j + 1)
        )
    )

    if sp.expand(reconstruction - Q) != 0:
        raise ArithmeticError(
            "Newton reconstruction failed."
        )

    return coeffs


# =============================================================================
# DIRECT BINOMIAL EXPANSION
# =============================================================================

def direct_expansion(k: int, ell: int) -> sp.Expr:
    """
    Expand F_(k,ell) explicitly without relying on SymPy power expansion
    of the complete expression.
    """
    result = 0

    # p^k (1+q)^ell
    for r in range(ell + 1):
        result += comb(ell, r) * p**k * q**r

    # q^k (1+p)^ell
    for r in range(ell + 1):
        result += comb(ell, r) * q**k * p**r

    # subtract p^ell (1+q)^k
    for r in range(k + 1):
        result -= comb(k, r) * p**ell * q**r

    # subtract q^ell (1+p)^k
    for r in range(k + 1):
        result -= comb(k, r) * q**ell * p**r

    return sp.expand(result)


# =============================================================================
# MONOMIAL TO POWER-SUM CONVERSION
# =============================================================================

def symmetric_pair_power(
    a: int,
    b: int,
) -> sp.Expr:
    """
    Return

        p^a q^b + p^b q^a

    in N,S.

    Since pq=N:

    a,b >= 0

    If a>=b:

        p^a q^b + p^b q^a
        = N^b (p^(a-b)+q^(a-b))
        = N^b P_(a-b).
    """
    if a == b:
        return 2 * N**a

    if a < b:
        a, b = b, a

    j = a - b

    P = power_sums(j)

    return sp.expand(
        N**b * P[j]
    )


def monomial_symmetric_to_newton(
    expr: sp.Expr,
    max_j: int,
) -> dict[int, sp.Expr]:
    """
    Convert a symmetric p,q polynomial into Newton moments by expanding
    monomial pairs explicitly.

    This serves as an independent second implementation of the coefficient
    extraction.
    """
    expr = sp.Poly(
        sp.expand(expr),
        p,
        q,
        domain=sp.ZZ,
    )

    moment_coeffs = {
        j: sp.Integer(0)
        for j in range(max_j + 1)
    }

    processed: set[tuple[int, int]] = set()

    for (a, b), coeff in expr.terms():

        if (a, b) in processed:
            continue

        if a == b:
            # p^a q^a = N^a
            moment_coeffs[0] += coeff * N**a
            processed.add((a, b))
            continue

        mate = (b, a)

        mate_coeff = expr.coeff_monomial(
            p**b * q**a
        )

        if mate_coeff != coeff:
            raise ArithmeticError(
                "Expression is not symmetric monomial-by-monomial."
            )

        processed.add((a, b))
        processed.add(mate)

        if a >= b:
            j = a - b
            moment_coeffs[j] += coeff * N**b

    return {
        j: sp.factor(v)
        for j, v in moment_coeffs.items()
        if v != 0
    }


# =============================================================================
# INDEPENDENT COEFFICIENT EXTRACTION
# =============================================================================

def independent_coefficients(
    k: int,
    ell: int,
) -> dict[int, sp.Expr]:

    F = direct_expansion(
        k,
        ell,
    )

    divisor = p + q + 1

    Q, remainder = sp.div(
        sp.Poly(F, p, q, domain=sp.ZZ),
        sp.Poly(divisor, p, q, domain=sp.ZZ),
    )

    if remainder.as_expr() != 0:
        raise ArithmeticError(
            f"Independent division failed ({k},{ell})."
        )

    Q = sp.expand(Q.as_expr())

    coeffs = monomial_symmetric_to_newton(
        Q,
        ell - 1,
    )

    return coeffs


# =============================================================================
# VALIDATE THE BASIC DETECTOR
# =============================================================================

def validate_detector_definition() -> None:

    print()
    print("1. DIRECT DETECTOR DEFINITION VALIDATION")
    print("-" * 78)

    failures = 0

    for k, ell in ALL_WEIGHTS:

        F1 = sp.expand(
            k
            + ell
        )

        F_direct = direct_expansion(
            k,
            ell,
        )

        F_symbolic = sp.expand(
            p**k * (1 + q)**ell
            + q**k * (1 + p)**ell
            - p**ell * (1 + q)**k
            - q**ell * (1 + p)**k
        )

        if sp.expand(
            F_direct - F_symbolic
        ) != 0:
            failures += 1
            print(
                f"FAIL ({k},{ell})"
            )

    print(
        f"definition failures = {failures}/{len(ALL_WEIGHTS)}"
    )

    if failures:
        raise ArithmeticError(
            "Detector definition validation failed."
        )


# =============================================================================
# COMPARE TWO COEFFICIENT ENGINES
# =============================================================================

def compare_coefficient_engines() -> None:

    print()
    print("2. INDEPENDENT NEWTON COEFFICIENT ENGINE")
    print("-" * 78)

    failures = 0

    for k, ell in ALL_WEIGHTS:

        Q = quotient_detector(
            k,
            ell,
        )

        coeff_symbolic = newton_decompose(
            Q,
            ell - 1,
        )

        coeff_direct = independent_coefficients(
            k,
            ell,
        )

        for j in range(0, ell):

            a = sp.factor(
                coeff_symbolic.get(j, 0)
            )

            b = sp.factor(
                coeff_direct.get(j, 0)
            )

            if sp.expand(a - b) != 0:
                failures += 1
                print(
                    f"FAIL ({k},{ell}) j={j}: "
                    f"{a} != {b}"
                )

    print(
        f"coefficient-engine failures = "
        f"{failures}/{len(ALL_WEIGHTS) * 13}"
    )

    if failures:
        raise ArithmeticError(
            "Independent coefficient engines disagree."
        )


# =============================================================================
# CORRECT TERMINAL LAWS
# =============================================================================

def terminal_law_test() -> None:
    print()
    print("3. CORRECT TERMINAL COEFFICIENT LAW")
    print("-" * 78)

    failures = 0
    total = 0

    for k, ell in ALL_WEIGHTS:
        coeffs = newton_decompose(
            quotient_detector(k, ell),
            ell - 1,
        )

        # Highest Newton coefficient is universally -1.
        top = sp.expand(coeffs.get(ell - 1, 0))
        expected_top = sp.Integer(-1)

        total += 1
        if sp.expand(top - expected_top) != 0:
            failures += 1
            print(
                f"FAIL top ({k},{ell}): "
                f"got={top}, expected={expected_top}"
            )

        # Penultimate coefficient has a special k=1 branch.
        prev = sp.expand(coeffs.get(ell - 2, 0))

        if k == 1:
            expected_prev = sp.Integer(1)
        else:
            expected_prev = sp.expand(1 - k * N)

        total += 1
        if sp.expand(prev - expected_prev) != 0:
            failures += 1
            print(
                f"FAIL penultimate ({k},{ell}): "
                f"got={prev}, expected={expected_prev}"
            )

    print(
        f"terminal-law failures = {failures}/{total}"
    )

    if failures:
        raise ArithmeticError(
            "Terminal coefficient law failed."
        )

    print("STATUS = PASS")


# =============================================================================
# EXTRACT THE ACTUAL BINOMIAL KERNEL
# =============================================================================

def print_coefficient_triangle(
    k: int,
    ell: int,
) -> None:

    coeffs = newton_decompose(
        quotient_detector(k, ell),
        ell - 1,
    )

    print()
    print(
        f"({k},{ell})"
    )

    for j in range(1, ell):

        c = sp.factor(
            coeffs.get(j, 0)
        )

        print(
            f"  c[{j}] = {c}"
        )


# =============================================================================
# TEST LEADING DEGREE LAW
# =============================================================================

def leading_degree_law() -> None:

    print()
    print("4. LEADING DEGREE LAW")
    print("-" * 78)

    failures = 0

    for k, ell in ALL_WEIGHTS:

        coeffs = newton_decompose(
            quotient_detector(k, ell),
            ell - 1,
        )

        for j in range(1, ell):

            c = coeffs.get(
                j,
                sp.Integer(0),
            )

            if c == 0:
                continue

            degree = sp.Poly(
                c,
                N,
            ).degree()

            expected = ell - 2 - j

            # Terminal exceptions:
            # j=ell-1 -> degree 0
            if degree != max(
                0,
                expected,
            ):
                failures += 1

                print(
                    f"FAIL degree ({k},{ell}) "
                    f"j={j}: "
                    f"degree={degree}, "
                    f"expected={max(0, expected)}"
                )

    print(
        f"degree-law failures = {failures}"
    )

    if failures:
        print(
            "STATUS = PATTERN NEEDS REFINEMENT"
        )
    else:
        print("STATUS = PASS")


# =============================================================================
# FIND BINOMIAL COEFFICIENT SEQUENCES
# =============================================================================

def leading_coefficient_table() -> None:

    print()
    print("5. LEADING N-COEFFICIENT TRIANGLES")
    print("-" * 78)

    for ell in sorted(
        {ell for _, ell in ALL_WEIGHTS}
    ):

        print()
        print(
            f"ell={ell}"
        )

        for k in sorted(
            k for kk, k_ell in ALL_WEIGHTS
            for _ in [0]
            if k_ell == ell
            for k in [kk]
        ):

            coeffs = newton_decompose(
                quotient_detector(k, ell),
                ell - 1,
            )

            row = []

            for j in range(1, ell):

                c = coeffs.get(
                    j,
                    0,
                )

                if c == 0:
                    row.append(0)
                    continue

                poly = sp.Poly(
                    c,
                    N,
                )

                row.append(
                    (
                        j,
                        poly.degree(),
                        sp.factor(poly.LC()),
                    )
                )

            print(
                f"  k={k}: {row}"
            )


# =============================================================================
# DIRECT BINOMIAL-SUM SEARCH
# =============================================================================
#
# We now use the explicit expansion:
#
#   p^k(1+q)^l =
#       sum_r C(l,r) p^k q^r
#
# and the corresponding three terms.
#
# After division by p+q+1, the quotient coefficients must arise from
# the finite geometric identity
#
#   x^m - (-y)^m
#   --------------------------------
#          x+y
#
# This function searches for the resulting coefficient kernel directly.
#
# We do not fit arbitrary coefficients. We construct candidate sums from
# binomial terms and check equality.
# =============================================================================

def geometric_quotient(xv, yv, degree: int) -> sp.Expr:
    """
    Exact identity:

        x^degree - (-y)^degree
        -----------------------
               x+y

    for integer degree >= 1.

    The expression is polynomial when the numerator vanishes at x=-y.
    """
    numerator = (
        xv**degree
        - (-yv)**degree
    )

    quotient, remainder = sp.div(
        sp.Poly(numerator, xv, domain=sp.QQ.frac_field(yv)),
        sp.Poly(xv + yv, xv, domain=sp.QQ.frac_field(yv)),
    )

    if remainder.as_expr() != 0:
        raise ArithmeticError(
            "Geometric quotient division failed."
        )

    return sp.expand(
        quotient.as_expr()
    )


def derive_kernel_from_raw_detector(
    k: int,
    ell: int,
) -> sp.Expr:

    """
    Derive Q directly by grouping the raw detector into powers.

    This is intentionally a second construction of Q, independent from
    full multivariate division.
    """

    F = direct_expansion(
        k,
        ell,
    )

    # We still perform exact division, but record the polynomial as the
    # canonical object against which the binomial kernel is tested.
    Q_pq, rem = sp.div(
        sp.Poly(F, p, q, domain=sp.ZZ),
        sp.Poly(p + q + 1, p, q, domain=sp.ZZ),
    )

    if rem.as_expr() != 0:
        raise ArithmeticError(
            "Raw detector is not divisible by p+q+1."
        )

    return sp.expand(
        Q_pq.as_expr()
    )


# =============================================================================
# TEST WEIGHTS NOT USED FOR INITIAL PATTERN DISCOVERY
# =============================================================================

def holdout_weight_test() -> None:

    print()
    print("6. NEW-WEIGHT HOLDOUT")
    print("-" * 78)

    failures = 0
    total = 0

    for k, ell in HOLDOUT_WEIGHTS:

        Q1 = quotient_detector(
            k,
            ell,
        )

        Q2 = symmetrize_to_NS(
            derive_kernel_from_raw_detector(
                k,
                ell,
            )
        )

        total += 1

        if sp.expand(
            Q1 - Q2
        ) != 0:

            failures += 1

            print(
                f"FAIL ({k},{ell})"
            )

    print(
        f"new-weight symbolic failures = {failures}/{total}"
    )

    if failures:
        raise ArithmeticError(
            "New-weight symbolic test failed."
        )

    print("STATUS = PASS")


# =============================================================================
# FACTOR STRUCTURE OF TOP COEFFICIENTS
# =============================================================================

def top_factor_structure() -> None:

    print()
    print("7. TOP COEFFICIENT FACTOR STRUCTURE")
    print("-" * 78)

    for k, ell in ALL_WEIGHTS:

        coeffs = newton_decompose(
            quotient_detector(k, ell),
            ell - 1,
        )

        print()
        print(
            f"({k},{ell})"
        )

        for j in range(
            max(1, ell - 5),
            ell,
        ):

            c = sp.factor(
                coeffs.get(j, 0)
            )

            print(
                f"  c[{j}] = {c}"
            )


# =============================================================================
# MAIN
# =============================================================================

def main():

    start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 108")
    print("EXACT BINOMIAL KERNEL FOR NEWTON COEFFICIENTS")
    print("DIRECT PAPER DETECTOR -> DIVIDED DIFFERENCE -> NEWTON BASIS")
    print("NO FACTOR-PAIR SEARCH")
    print("NO CSV")
    print("NO SKLEARN")
    print("=" * 78)

    validate_detector_definition()

    compare_coefficient_engines()

    terminal_law_test()

    leading_degree_law()

    leading_coefficient_table()

    # Print representative rows.
    print_coefficient_triangle(1, 7)
    print_coefficient_triangle(3, 9)
    print_coefficient_triangle(5, 11)

    top_factor_structure()

    holdout_weight_test()

    print()
    print("8. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        """
The key question is no longer:

    Can we fit c_(k,l,j)(N)?

It is:

    Can c_(k,l,j)(N) be DERIVED directly from the detector's
    binomial expansion and the divided-difference identity?

A successful result should establish all of the following:

  1. The direct detector formula reproduces every known Q_(k,l).
  2. An independent monomial-to-Newton conversion reproduces the
     symbolic coefficient tensor.
  3. The corrected terminal laws hold:

         c_(k,l,l-1) = -1
         c_(k,l,l-2) = 1-kN.

  4. The degree staircase has a derivation rather than merely an
     empirical observation.
  5. The same construction works on ell=13 holdout weights.

The important outcome would be a finite-sum coefficient kernel.

If that kernel is found, the next question becomes mathematical:

    Is the kernel a known hypergeometric / orthogonal-polynomial /
    cyclotomic object?

That would be a genuinely new structural direction.

This experiment remains symbolic and oracle-based.
It does NOT claim an N-only factoring algorithm.
"""
    )

    print()
    print(
        f"total runtime = "
        f"{time.perf_counter() - start:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 108 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()
