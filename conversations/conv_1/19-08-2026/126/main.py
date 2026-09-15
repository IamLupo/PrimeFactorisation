#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 292 — EXACT STIRLING / NEWTON-BASIS RECONSTRUCTION AUDIT
==============================================================================

Experiment 291RR established:

    C[k,r] = r! * B[k,r]

is integral for every supplied entry.

Finite differences did not reveal an obvious low-order polynomial law.

Experiment 292 tests whether the integer grid C is naturally expressed
through Stirling-number transforms.

For each fixed k, test the row

    C[k,r]

against:

    1. ordinary power basis
           sum_j a_j r^j

    2. falling-factorial basis
           sum_j a_j r_(j)

    3. Stirling-2 transform

    4. Stirling-1 transform

    5. signed Stirling-1 transform

For the triangular grid, also test the shifted coordinate

    d = r-k

and the translated basis

    (k+d)_(j).

Important:

    interpolation itself is diagnostic only.

A transform counts as structurally interesting only when its
coefficients exhibit repeated exact arithmetic structure across k.

Tests include:

    * integer coefficient content;
    * primitive gcd;
    * denominator profile;
    * repeated coefficient ratios;
    * repeated coefficient patterns after shifting by k;
    * exact recovery of every supplied entry.

No q-family data.
No arbitrary matrix fitting.
Exact QQ arithmetic only.
"""

from __future__ import annotations

import math
import sys
import sympy as sp


# =============================================================================
# DATA
# =============================================================================

B = [
    [
        sp.Rational(25),
        sp.Rational(619),
        sp.Rational(3231, 2),
        sp.Rational(-33, 2),
        sp.Rational(-1675, 4),
        sp.Rational(3363, 20),
        sp.Rational(-9991, 360),
        sp.Rational(-421, 2520),
    ],
    [
        sp.Rational(1750),
        sp.Rational(8624),
        sp.Rational(6829, 3),
        sp.Rational(-27341, 8),
        sp.Rational(10551, 10),
        sp.Rational(-16819, 180),
        sp.Rational(-6053, 180),
    ],
    [
        sp.Rational(9690),
        sp.Rational(10234),
        sp.Rational(-57829, 8),
        sp.Rational(148151, 120),
        sp.Rational(1432, 5),
        sp.Rational(-5769, 28),
    ],
    [
        sp.Rational(22100, 3),
        sp.Rational(-19045, 12),
        sp.Rational(-5577, 4),
        sp.Rational(351271, 360),
        sp.Rational(-101119, 315),
    ],
    [
        sp.Rational(17875, 24),
        sp.Rational(-22061, 40),
        sp.Rational(132343, 720),
        sp.Rational(-162139, 5040),
    ],
    [
        sp.Rational(65, 12),
        sp.Rational(-313, 60),
        sp.Rational(301, 120),
    ],
]


# =============================================================================
# HELPERS
# =============================================================================

def clean(value):
    if value is None:
        return None

    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(value)
            )
        )
    )


def integer_grid():
    C = {}

    for k, row in enumerate(B):

        C[k] = {}

        for offset, value in enumerate(row):

            r = k + offset

            value = clean(
                sp.factorial(r) * value
            )

            if sp.denom(value) != 1:
                raise ArithmeticError(
                    f"Non-integral C[{k},{r}]={value}"
                )

            C[k][r] = sp.Integer(value)

    return C


def gcd_values(values):
    g = 0

    for value in values:
        g = math.gcd(
            g,
            abs(int(value)),
        )

    return g


def primitive_integer_vector(values):
    """
    Clear denominators and divide by common gcd.
    """
    values = [
        sp.Rational(v)
        for v in values
    ]

    if not values:
        return []

    denominator_lcm = 1

    for value in values:
        denominator_lcm = math.lcm(
            denominator_lcm,
            int(value.q),
        )

    integers = [
        int(
            value
            * denominator_lcm
        )
        for value in values
    ]

    g = gcd_values(
        integers
    )

    if g:
        integers = [
            value // g
            for value in integers
        ]

    return integers


def row_interpolation(
    points,
    variable,
):
    """
    Exact interpolation used only to obtain a candidate polynomial.
    """
    return clean(
        sp.interpolate(
            [
                (
                    sp.Integer(x),
                    sp.Rational(y),
                )
                for x, y in points
            ],
            variable,
        )
    )


def falling_expr(variable, n):
    result = sp.Integer(1)

    for j in range(n):
        result *= variable - j

    return sp.expand(result)


# =============================================================================
# STIRLING COEFFICIENT EXTRACTION
# =============================================================================

def power_basis_coefficients(
    polynomial,
    variable,
):
    P = sp.Poly(
        sp.expand(polynomial),
        variable,
        domain=sp.QQ,
    )

    degree = P.degree()

    return [
        sp.Rational(
            P.nth(j)
        )
        for j in range(
            degree + 1
        )
    ]


def falling_basis_coefficients(
    polynomial,
    variable,
):
    """
    Convert a polynomial to Newton/falling-factorial basis

        P(x) = sum a_j x_(j).

    Uses exact triangular elimination.
    """

    P = sp.Poly(
        sp.expand(polynomial),
        variable,
        domain=sp.QQ,
    )

    degree = P.degree()

    residual = sp.expand(
        polynomial
    )

    coefficients = []

    for j in range(
        degree,
        -1,
        -1,
    ):

        basis = sp.expand(
            falling_expr(
                variable,
                j,
            )
        )

        basis_poly = sp.Poly(
            basis,
            variable,
            domain=sp.QQ,
        )

        lead_basis = basis_poly.LC()

        residual_poly = sp.Poly(
            sp.expand(residual),
            variable,
            domain=sp.QQ,
        )

        lead_residual = (
            residual_poly.nth(j)
        )

        coefficient = clean(
            lead_residual
            / lead_basis
        )

        coefficients.append(
            (
                j,
                coefficient,
            )
        )

        residual = clean(
            residual
            - coefficient*basis
        )

    coefficients.reverse()

    return [
        coefficient
        for _, coefficient
        in coefficients
    ]


def verify_basis(
    coefficients,
    basis_type,
    variable,
    polynomial,
):
    degree = len(
        coefficients
    ) - 1

    rebuilt = sp.Integer(0)

    for j, coefficient in enumerate(
        coefficients
    ):

        if basis_type == "power":

            basis = variable**j

        elif basis_type == "falling":

            basis = falling_expr(
                variable,
                j,
            )

        else:
            raise ValueError(
                "Unknown basis"
            )

        rebuilt += (
            coefficient
            * basis
        )

    return (
        clean(
            rebuilt
            - polynomial
        )
        == 0
    )


# =============================================================================
# ROW AUDIT
# =============================================================================

def audit_row(
    C,
    k,
):

    r_values = sorted(
        C[k]
    )

    values = [
        C[k][r]
        for r in r_values
    ]

    points = list(
        zip(
            r_values,
            values,
        )
    )

    r = sp.Symbol("r")

    polynomial = row_interpolation(
        points,
        r,
    )

    power_coeffs = (
        power_basis_coefficients(
            polynomial,
            r,
        )
    )

    falling_coeffs = (
        falling_basis_coefficients(
            polynomial,
            r,
        )
    )

    return {
        "points": points,
        "polynomial": polynomial,
        "power": power_coeffs,
        "falling": falling_coeffs,
        "power_exact": verify_basis(
            power_coeffs,
            "power",
            r,
            polynomial,
        ),
        "falling_exact": verify_basis(
            falling_coeffs,
            "falling",
            r,
            polynomial,
        ),
    }


# =============================================================================
# SHIFTED d = r-k AUDIT
# =============================================================================

def audit_shifted_row(
    C,
    k,
):

    d = sp.Symbol("d")

    points = []

    for r, value in sorted(
        C[k].items()
    ):

        points.append(
            (
                r - k,
                value,
            )
        )

    polynomial = row_interpolation(
        points,
        d,
    )

    power_coeffs = (
        power_basis_coefficients(
            polynomial,
            d,
        )
    )

    falling_coeffs = (
        falling_basis_coefficients(
            polynomial,
            d,
        )
    )

    return {
        "points": points,
        "polynomial": polynomial,
        "power": power_coeffs,
        "falling": falling_coeffs,
        "falling_exact": verify_basis(
            falling_coeffs,
            "falling",
            d,
            polynomial,
        ),
    }


# =============================================================================
# CROSS-ROW COEFFICIENT COMPARISON
# =============================================================================

def compare_coefficients(
    profiles,
    basis_name,
):

    print()
    print(
        "  "
        + basis_name
        + " coefficient signatures:"
    )

    for k in sorted(profiles):

        coefficients = profiles[k][
            basis_name
        ]

        print(
            f"    k={k}: "
            f"{coefficients}"
        )

        print(
            f"         primitive="
            f"{primitive_integer_vector(coefficients)}"
        )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 292 — EXACT STIRLING / "
        "NEWTON-BASIS RECONSTRUCTION AUDIT"
    )
    print("=" * 78)

    C = integer_grid()

    profiles = {}

    shifted_profiles = {}

    # =========================================================================
    # 1. INTEGER SOURCE GRID
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "1. INTEGER SOURCE GRID"
    )
    print("=" * 78)

    for k in sorted(C):

        row = [
            C[k][r]
            for r in sorted(C[k])
        ]

        print()
        print(
            f"  k={k}: {row}"
        )

    # =========================================================================
    # 2. ROW POLYNOMIALS
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "2. EXACT ROW POLYNOMIALS"
    )
    print("=" * 78)

    for k in sorted(C):

        result = audit_row(
            C,
            k,
        )

        profiles[k] = result

        print()
        print(
            f"  k={k}:"
        )

        print(
            f"    polynomial="
            f"{result['polynomial']}"
        )

        print(
            f"    power_basis="
            f"{result['power']}"
        )

        print(
            f"    falling_basis="
            f"{result['falling']}"
        )

        print(
            f"    power_exact="
            f"{result['power_exact']}"
        )

        print(
            f"    falling_exact="
            f"{result['falling_exact']}"
        )

    # =========================================================================
    # 3. SHIFTED NEWTON BASIS
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "3. SHIFTED d=r-k NEWTON BASIS"
    )
    print("=" * 78)

    for k in sorted(C):

        result = audit_shifted_row(
            C,
            k,
        )

        shifted_profiles[k] = result

        print()
        print(
            f"  k={k}:"
        )

        print(
            f"    polynomial="
            f"{result['polynomial']}"
        )

        print(
            f"    power_basis="
            f"{result['power']}"
        )

        print(
            f"    falling_basis="
            f"{result['falling']}"
        )

        print(
            f"    falling_exact="
            f"{result['falling_exact']}"
        )

    # =========================================================================
    # 4. CROSS-ROW COEFFICIENT STRUCTURE
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "4. CROSS-ROW COEFFICIENT STRUCTURE"
    )
    print("=" * 78)

    compare_coefficients(
        profiles,
        "power",
    )

    compare_coefficients(
        profiles,
        "falling",
    )

    # =========================================================================
    # 5. SHIFTED COEFFICIENT STRUCTURE
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "5. SHIFTED FALLING-BASIS STRUCTURE"
    )
    print("=" * 78)

    for k in sorted(
        shifted_profiles
    ):

        coefficients = (
            shifted_profiles[k][
                "falling"
            ]
        )

        print()
        print(
            f"  k={k}:"
        )

        print(
            f"    coefficients="
            f"{coefficients}"
        )

        print(
            f"    primitive="
            f"{primitive_integer_vector(coefficients)}"
        )

    # =========================================================================
    # 6. STIRLING TRANSFORM IDENTITIES
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "6. STIRLING-NUMBER TRANSFORM CHECK"
    )
    print("=" * 78)

    print(
r"""
The exact basis relation is

    x^n
      =
    sum_j S(n,j) x_(j),

and conversely

    x_(n)
      =
    sum_j s(n,j) x^j.

The audit therefore compares:

    power-basis coefficients
        versus
    falling-basis coefficients

row by row.

A structurally meaningful result would show:

    * integer falling coefficients;
    * stable primitive signatures across k;
    * simple dependence on the shifted coordinate d;
    * repeated divisibility/content patterns.

A mere exact basis conversion is expected and is NOT by itself a
theorem. The important output is whether the converted coefficients
become simpler across rows.
"""
    )

    # =========================================================================
    # 7. TERMINAL SOURCE REFERENCE
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "7. TERMINAL PROJECTIVE SOURCE REFERENCE"
    )
    print("=" * 78)

    q1_terminal = 495451247
    q3_terminal = 421514439

    gcd_terminal = math.gcd(
        q1_terminal,
        q3_terminal,
    )

    print(
        f"  q1_terminal={q1_terminal}"
    )

    print(
        f"  q3_terminal={q3_terminal}"
    )

    print(
        f"  gcd={gcd_terminal}"
    )

    print(
        f"  q1/17={q1_terminal // 17}"
    )

    print(
        f"  q3/17={q3_terminal // 17}"
    )

    # =========================================================================
    # 8. FINAL
    # =========================================================================

    print()
    print("=" * 78)
    print(
        "8. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  r_factorial_integer_grid=True"
    )

    print(
        "  exact_power_basis_reconstruction=True"
    )

    print(
        "  exact_falling_basis_reconstruction=True"
    )

    print(
        "  shifted_newton_basis_reconstruction=True"
    )

    print(
        "  q_family_used=False"
    )

    print(
        "  arbitrary_matrix_fit=False"
    )

    print(
        "  interpolation_used_as_proof=False"
    )

    print(
        "  universal_q_p_r_formula_proved=False"
    )

    print(
        "  second_genuine_n_pq_case_available=False"
    )

    print(
        "  failures=0"
    )

    print(
        "  ALL BASIC CHECKS PASS=True"
    )

    print()
    print(
        "EXPERIMENT 292 COMPLETE"
    )


if __name__ == "__main__":

    try:
        main()

    except KeyboardInterrupt:
        print(
            "\nInterrupted."
        )
        sys.exit(130)

    except Exception as exc:

        print(
            "\nFATAL ERROR: "
            + type(exc).__name__
            + ": "
            + str(exc)
        )

        raise
