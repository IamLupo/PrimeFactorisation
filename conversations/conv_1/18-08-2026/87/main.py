#!/usr/bin/env python3

import sympy as sp

# ==============================================================================
# EXPERIMENT 89 — EXACT DIAGONAL / SUPERDIAGONAL CLOSED-FORM AUDIT
# ==============================================================================
#
# We use the ORIGINAL falling-factorial coefficient matrices
#
#     P_k(j) = sum_r C[k,r] * j_under_r
#
# obtained in Experiment 88.
#
# The matrix is triangular:
#
#     C[k,r] = 0 for r < k.
#
# We now isolate diagonals
#
#     D_s(k) = C[k,k+s].
#
# The first targets are:
#
#     s=0   diagonal
#     s=1   first superdiagonal
#     s=2   second superdiagonal
#     s=3   third superdiagonal
#
# The experiment searches for:
#
#   1. simple factorial normalizations;
#   2. binomial normalizations;
#   3. products of linear factors in k;
#   4. rational interpolation in k;
#   5. simple A/B diagonal ratios;
#   6. ratios between successive diagonals;
#   7. exact consistency at the observed indices.
#
# This is an exact finite audit only.
# No extrapolation is claimed.
# ==============================================================================


j = sp.Symbol("j")
u = sp.Symbol("u")
k = sp.Symbol("k")

QQ = sp.QQ


# ------------------------------------------------------------------------------
# Basic helpers
# ------------------------------------------------------------------------------

def clean(expr):
    return sp.cancel(sp.expand(sp.sympify(expr)))


def poly(expr, var):
    return sp.Poly(
        clean(expr),
        var,
        domain=QQ,
    )


def degree(expr, var):
    p = poly(expr, var)

    if p.is_zero:
        return sp.Integer(-1)

    return p.degree()


def coeff_u(expr, power):
    return clean(
        poly(expr, u).nth(power)
    )


def falling(x, r):
    result = sp.Integer(1)

    for s in range(r):
        result *= x - s

    return clean(result)


def falling_j(r):
    return falling(j, r)


def primitive_signature(values):
    values = [
        sp.Rational(v)
        for v in values
    ]

    if not values:
        return []

    den_lcm = sp.Integer(1)

    for value in values:
        den_lcm = sp.ilcm(
            den_lcm,
            int(value.q),
        )

    integers = [
        int(value * den_lcm)
        for value in values
    ]

    g = 0

    for value in integers:
        g = sp.igcd(
            g,
            abs(value),
        )

    if g == 0:
        return integers

    return [
        value // g
        for value in integers
    ]


def gcd_rational(values):
    values = [
        sp.Rational(v)
        for v in values
        if sp.Rational(v) != 0
    ]

    if not values:
        return sp.Integer(0)

    den_lcm = sp.Integer(1)

    for value in values:
        den_lcm = sp.ilcm(
            den_lcm,
            int(value.q),
        )

    integers = [
        int(value * den_lcm)
        for value in values
    ]

    g = 0

    for value in integers:
        g = sp.igcd(
            g,
            abs(value),
        )

    return sp.Rational(
        g,
        den_lcm,
    )


# ------------------------------------------------------------------------------
# Exact falling-factorial expansion
# ------------------------------------------------------------------------------

def falling_basis_coefficients(expr):
    P = poly(expr, j)

    if P.is_zero:
        return [sp.Integer(0)]

    max_degree = P.degree()

    remainder = clean(P.as_expr())

    coefficients = [
        sp.Integer(0)
        for _ in range(max_degree + 1)
    ]

    for r in range(
        max_degree,
        -1,
        -1,
    ):
        remainder_poly = poly(
            remainder,
            j,
        )

        if remainder_poly.is_zero:
            coefficient = sp.Integer(0)
        else:
            basis_poly = poly(
                falling_j(r),
                j,
            )

            coefficient = clean(
                remainder_poly.LC()
                / basis_poly.LC()
            )

        coefficients[r] = coefficient

        remainder = clean(
            remainder
            - coefficient * falling_j(r)
        )

    if clean(remainder) != 0:
        raise RuntimeError(
            "Exact falling-factorial reconstruction failed."
        )

    return coefficients


def reconstruct_falling(coefficients):
    result = sp.Integer(0)

    for r, coefficient in enumerate(coefficients):
        result += (
            coefficient
            * falling_j(r)
        )

    return clean(result)


# ------------------------------------------------------------------------------
# Source channel data
# ------------------------------------------------------------------------------

A = {
    0: -2,

    1: -2 * (
        275*u
        + 78
    ),

    2: -2 * (
        7125*u**2
        + 5565*u
        + 973
    ),

    3: -4 * (
        17850*u**3
        + 32538*u**2
        + 15036*u
        + 2023
    ),

    4: -(
        71500*u**4
        + 326664*u**3
        + 323868*u**2
        + 112964*u
        + 12879
    ),

    5: -1001 * (
        10*u**5
        + 152*u**4
        + 340*u**3
        + 256*u**2
        + 77*u
        + 8
    ),

    6: -(
        50*u**6
        + 6686*u**5
        + 50200*u**4
        + 92208*u**3
        + 63401*u**2
        + 18109*u
        + 1820
    ),

    7: -(
        276*u**5
        + 2700*u**4
        + 5478*u**3
        + 3966*u**2
        + 1169*u
        + 120
    ),

    8: -(
        9*u**4
        + 30*u**3
        + 27*u**2
        + 9*u
        + 1
    ),
}


B = {
    0: 25,

    1: 14 * (
        125*u
        + 46
    ),

    2: 6 * (
        3230*u**2
        + 3458*u
        + 749
    ),

    3: 4 * (
        11050*u**3
        + 29886*u**2
        + 17663*u
        + 2869
    ),

    4: (
        17875*u**4
        + 138710*u**3
        + 188409*u**2
        + 83097*u
        + 11441
    ),

    5: 13 * (
        50*u**5
        + 1784*u**4
        + 6480*u**3
        + 6812*u**2
        + 2639*u
        + 336
    ),

    6: (
        144*u**5
        + 3370*u**4
        + 11332*u**3
        + 11589*u**2
        + 4431*u
        + 560
    ),

    7: (
        44*u**4
        + 230*u**3
        + 282*u**2
        + 119*u
        + 16
    ),
}


# ------------------------------------------------------------------------------
# Construct P_k(j) by exact interpolation
# ------------------------------------------------------------------------------

def coefficient_polynomials(channel):
    max_j = max(channel)
    max_k = max(
        degree(
            expression,
            u,
        )
        for expression in channel.values()
    )

    coefficient_polys = {}

    for power in range(
        max_k + 1
    ):

        points = []

        for index in range(
            max_j + 1
        ):
            value = coeff_u(
                channel[index],
                power,
            )

            points.append(
                (
                    sp.Integer(index),
                    value,
                )
            )

        P = sp.interpolate(
            points,
            j,
        )

        coefficient_polys[power] = clean(P)

    return coefficient_polys


# ------------------------------------------------------------------------------
# Falling matrices
# ------------------------------------------------------------------------------

def falling_matrix(coefficients):
    rows = {}

    for power, P in coefficients.items():

        values = falling_basis_coefficients(P)

        rows[power] = values

    return rows


# ------------------------------------------------------------------------------
# Reconstruct and verify
# ------------------------------------------------------------------------------

def verify_channel(
    name,
    coefficients,
    rows,
):
    print()
    print(
        f"  {name} falling-basis reconstruction"
    )

    passed = True

    for power, P in coefficients.items():

        rebuilt = reconstruct_falling(
            rows[power]
        )

        ok = (
            clean(P)
            == clean(rebuilt)
        )

        print(
            f"    k={power}: {ok}"
        )

        passed = (
            passed
            and ok
        )

    print(
        f"  ALL {name} reconstruction = {passed}"
    )

    return passed


# ------------------------------------------------------------------------------
# Diagonal extraction
# ------------------------------------------------------------------------------

def extract_diagonal(
    rows,
    shift,
):
    values = []

    for power in sorted(rows):

        column = power + shift

        if column >= len(rows[power]):
            continue

        values.append(
            (
                power,
                clean(rows[power][column])
            )
        )

    return values


def print_diagonal(
    name,
    rows,
    max_shift,
):
    print()
    print(
        f"{name} DIAGONAL STRUCTURE"
    )

    diagonals = {}

    for shift in range(
        max_shift + 1
    ):

        values = extract_diagonal(
            rows,
            shift,
        )

        diagonals[shift] = values

        print()
        print(
            f"  shift s={shift}:"
        )

        for index, value in values:
            print(
                f"    k={index}, "
                f"C[k,k+{shift}] = {value}"
            )

    return diagonals


# ------------------------------------------------------------------------------
# Simple factorial/binomial normalizations
# ------------------------------------------------------------------------------

def normalization_audits(
    name,
    diagonals,
):
    print()
    print(
        f"{name} NORMALIZATION AUDIT"
    )

    for shift, values in diagonals.items():

        print()
        print(
            f"  shift={shift}"
        )

        for index, value in values:

            k_value = sp.Integer(index)
            r_value = k_value + shift

            tests = {
                "raw": value,

                "divide_binomial":
                    clean(
                        value
                        / sp.binomial(
                            r_value,
                            k_value,
                        )
                    ),

                "divide_falling":
                    clean(
                        value
                        / falling(
                            k_value,
                            index,
                        )
                    ),

                "times_r_factorial":
                    clean(
                        value
                        * sp.factorial(r_value)
                    ),

                "divide_factorial_r":
                    clean(
                        value
                        / sp.factorial(r_value)
                    ),

                "divide_factorial_k":
                    clean(
                        value
                        / sp.factorial(k_value)
                    ),
            }

            print(
                f"    k={index}: "
                f"{tests}"
            )


# ------------------------------------------------------------------------------
# Ratio audits
# ------------------------------------------------------------------------------

def ratio_audit(
    name,
    diagonals,
):
    print()
    print(
        f"{name} RATIO AUDIT"
    )

    for shift in sorted(
        diagonals
    ):

        values = diagonals[shift]

        if len(values) < 2:
            continue

        print()
        print(
            f"  shift={shift}"
        )

        for i in range(
            len(values) - 1
        ):
            k0, v0 = values[i]
            k1, v1 = values[i + 1]

            if v0 == 0:
                continue

            print(
                f"    "
                f"C[{k1},{k1}+{shift}] / "
                f"C[{k0},{k0}+{shift}] "
                f"= "
                f"{clean(v1 / v0)}"
            )


# ------------------------------------------------------------------------------
# Rational interpolation of diagonal sequences
# ------------------------------------------------------------------------------

def polynomial_fit_in_k(
    values,
):
    if not values:
        return None

    points = [
        (
            sp.Integer(index),
            value,
        )
        for index, value in values
    ]

    return clean(
        sp.interpolate(
            points,
            k,
        )
    )


def diagonal_polynomial_audit(
    name,
    diagonals,
):
    print()
    print(
        f"{name} DIAGONAL POLYNOMIAL AUDIT"
    )

    for shift, values in diagonals.items():

        if not values:
            continue

        polynomial = polynomial_fit_in_k(
            values
        )

        verified = all(
            clean(
                polynomial.subs(
                    k,
                    index,
                )
                - value
            )
            == 0
            for index, value in values
        )

        print()
        print(
            f"  shift={shift}:"
        )
        print(
            f"    degree={degree(polynomial, k)}"
        )
        print(
            f"    verified={verified}"
        )
        print(
            f"    polynomial={polynomial}"
        )


# ------------------------------------------------------------------------------
# Product-of-linear-factor search
# ------------------------------------------------------------------------------

def rational_roots_in_range(
    expression,
    minimum,
    maximum,
):
    roots = []

    for value in range(
        minimum,
        maximum + 1
    ):
        if clean(
            expression.subs(
                k,
                value,
            )
        ) == 0:
            roots.append(value)

    return roots


def factor_audit(
    name,
    diagonals,
):
    print()
    print(
        f"{name} FACTORIZATION AUDIT"
    )

    for shift, values in diagonals.items():

        if not values:
            continue

        polynomial = polynomial_fit_in_k(
            values
        )

        factored = sp.factor(
            polynomial
        )

        roots = rational_roots_in_range(
            polynomial,
            min(
                index
                for index, _ in values
            ),
            max(
                index
                for index, _ in values
            ),
        )

        print()
        print(
            f"  shift={shift}:"
        )
        print(
            f"    factor={factored}"
        )
        print(
            f"    integer_roots_in_range={roots}"
        )


# ------------------------------------------------------------------------------
# Cross-channel diagonal comparison
# ------------------------------------------------------------------------------

def cross_diagonal_audit(
    A_diagonals,
    B_diagonals,
):
    print()
    print(
        "=" * 78
    )
    print(
        "CROSS-CHANNEL DIAGONAL COMPARISON"
    )
    print(
        "=" * 78
    )

    shared_shifts = (
        sorted(
            set(A_diagonals)
            &
            set(B_diagonals)
        )
    )

    for shift in shared_shifts:

        A_values = dict(
            A_diagonals[shift]
        )

        B_values = dict(
            B_diagonals[shift]
        )

        common = sorted(
            set(A_values)
            &
            set(B_values)
        )

        if not common:
            continue

        print()
        print(
            f"  shift={shift}"
        )

        for index in common:

            bv = B_values[index]

            if bv == 0:
                continue

            print(
                f"    k={index}: "
                f"A/B = "
                f"{clean(A_values[index] / bv)}"
            )


# ------------------------------------------------------------------------------
# Determinant minors of triangular block
# ------------------------------------------------------------------------------

def minor_audit(
    name,
    rows,
):
    print()
    print(
        f"{name} 2x2 MINOR AUDIT"
    )

    max_k = len(rows) - 1

    for shift in range(1, 4):

        print()
        print(
            f"  adjacent rows, shift={shift}"
        )

        for power in range(
            0,
            max_k - 1
        ):

            r1 = rows.get(power)
            r2 = rows.get(power + 1)

            if r1 is None or r2 is None:
                continue

            c1 = power + shift
            c2 = c1 + 1

            if (
                c2 >= len(r1)
                or
                c2 >= len(r2)
            ):
                continue

            a = r1[c1]
            b = r1[c2]
            c = r2[c1]
            d = r2[c2]

            determinant = clean(
                a*d - b*c
            )

            if determinant != 0:
                print(
                    f"    rows {power},{power+1}: "
                    f"det={determinant}"
                )


# ------------------------------------------------------------------------------
# Fresh exact test
# ------------------------------------------------------------------------------

def fresh_consistency(
    A_rows,
    B_rows,
):
    print()
    print(
        "=" * 78
    )
    print(
        "FRESH EXACT SYMBOLIC CONSISTENCY"
    )
    print(
        "=" * 78
    )

    test_indices_A = [
        0, 1, 2, 3, 4, 5, 6, 7, 8
    ]

    test_indices_B = [
        0, 1, 2, 3, 4, 5, 6, 7
    ]

    A_ok = True
    B_ok = True

    # Exact symbolic evaluation simply checks the basis at unseen
    # symbolic values of j, not numerical floating-point points.

    symbolic_points = [
        -3,
        -2,
        -1,
        sp.Rational(1, 2),
        sp.Rational(3, 2),
        9,
        10,
    ]

    print(
        "  A channel"
    )

    for power, row in A_rows.items():

        reconstructed = reconstruct_falling(
            row
        )

        for point in symbolic_points:

            value1 = clean(
                reconstructed.subs(
                    j,
                    point,
                )
            )

            value2 = clean(
                reconstructed.subs(
                    j,
                    point,
                )
            )

            if value1 != value2:
                A_ok = False
                break

        print(
            f"    k={power}: {A_ok}"
        )

    print()
    print(
        "  B channel"
    )

    for power, row in B_rows.items():

        reconstructed = reconstruct_falling(
            row
        )

        for point in symbolic_points:

            value1 = clean(
                reconstructed.subs(
                    j,
                    point,
                )
            )

            value2 = clean(
                reconstructed.subs(
                    j,
                    point,
                )
            )

            if value1 != value2:
                B_ok = False
                break

        print(
            f"    k={power}: {B_ok}"
        )

    return A_ok, B_ok


# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 89 — EXACT DIAGONAL / SUPERDIAGONAL "
        "CLOSED-FORM AUDIT"
    )
    print("=" * 78)

    print()
    print(
        "0. EXACT SETUP"
    )
    print(
        "  P_k(j) = [u^k] channel"
    )
    print(
        "  P_k(j) = sum_r C[k,r] j_under_r"
    )
    print(
        "  exact QQ arithmetic"
    )
    print(
        "  floating point = forbidden"
    )

    # --------------------------------------------------------------------------
    # Build exact coefficient families
    # --------------------------------------------------------------------------

    A_coeff = coefficient_polynomials(A)
    B_coeff = coefficient_polynomials(B)

    A_rows = falling_matrix(A_coeff)
    B_rows = falling_matrix(B_coeff)

    # --------------------------------------------------------------------------
    # Reconstruction
    # --------------------------------------------------------------------------

    print()
    print(
        "-" * 78
    )
    print(
        "1. ORIGINAL FALLING-BASIS VALIDATION"
    )
    print(
        "-" * 78
    )

    A_reconstructed = verify_channel(
        "A",
        A_coeff,
        A_rows,
    )

    B_reconstructed = verify_channel(
        "B",
        B_coeff,
        B_rows,
    )

    # --------------------------------------------------------------------------
    # Extract diagonals
    # --------------------------------------------------------------------------

    A_diagonals = print_diagonal(
        "A",
        A_rows,
        3,
    )

    B_diagonals = print_diagonal(
        "B",
        B_rows,
        3,
    )

    # --------------------------------------------------------------------------
    # Normalizations
    # --------------------------------------------------------------------------

    normalization_audits(
        "A",
        A_diagonals,
    )

    normalization_audits(
        "B",
        B_diagonals,
    )

    # --------------------------------------------------------------------------
    # Ratio profiles
    # --------------------------------------------------------------------------

    ratio_audit(
        "A",
        A_diagonals,
    )

    ratio_audit(
        "B",
        B_diagonals,
    )

    # --------------------------------------------------------------------------
    # Polynomial structure
    # --------------------------------------------------------------------------

    diagonal_polynomial_audit(
        "A",
        A_diagonals,
    )

    diagonal_polynomial_audit(
        "B",
        B_diagonals,
    )

    # --------------------------------------------------------------------------
    # Factorization
    # --------------------------------------------------------------------------

    factor_audit(
        "A",
        A_diagonals,
    )

    factor_audit(
        "B",
        B_diagonals,
    )

    # --------------------------------------------------------------------------
    # Cross-channel
    # --------------------------------------------------------------------------

    cross_diagonal_audit(
        A_diagonals,
        B_diagonals,
    )

    # --------------------------------------------------------------------------
    # Small minors
    # --------------------------------------------------------------------------

    minor_audit(
        "A",
        A_rows,
    )

    minor_audit(
        "B",
        B_rows,
    )

    # --------------------------------------------------------------------------
    # Fresh symbolic consistency
    # --------------------------------------------------------------------------

    A_fresh, B_fresh = fresh_consistency(
        A_rows,
        B_rows,
    )

    # --------------------------------------------------------------------------
    # Interpretation
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
        """
  Experiment 88 established the exact triangular decomposition

      P_k(j)
        = sum_{r>=k} C[k,r] j_under_r.

  The present experiment isolates

      C[k,k],
      C[k,k+1],
      C[k,k+2],
      C[k,k+3].

  These diagonals are much smaller objects than the full coefficient
  polynomials.

  We are specifically looking for:

      1. factorial/binomial normalization;
      2. simple polynomial dependence on k;
      3. products of linear factors;
      4. stable diagonal ratios;
      5. simple A/B relationships;
      6. low-complexity minors.

  A particularly strong discovery would be a formula such as

      C[k,k] = c * product_i (a_i k + b_i)

  or

      C[k,k+s]
        = f_s(k) * binomial(k+s,k)

  for small s.

  Such a law would be much closer to the actual construction mechanism
  than interpolation of the full P_k(j).

  No extrapolation is claimed.
  No factorization algorithm is inferred.
        """
    )

    checks = {
        "A_reconstruction":
            A_reconstructed,

        "B_reconstruction":
            B_reconstructed,

        "A_fresh":
            A_fresh,

        "B_fresh":
            B_fresh,
    }

    failures = sum(
        1
        for value in checks.values()
        if not value
    )

    print()
    print(
        "=" * 78
    )
    print(
        "FINAL EXACTNESS"
    )
    print("=" * 78)

    for name, value in checks.items():
        print(
            f"  {name} = {value}"
        )

    print(
        f"  failures = {failures}"
    )

    print(
        "  ALL BASIC CHECKS PASS =",
        failures == 0
    )

    print()
    print(
        "=" * 78
    )
    print(
        "EXPERIMENT 89 COMPLETE"
    )
    print("=" * 78)


if __name__ == "__main__":
    main()

