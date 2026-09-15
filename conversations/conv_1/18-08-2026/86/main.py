#!/usr/bin/env python3

import sympy as sp

# ==============================================================================
# EXPERIMENT 88 — EXACT ORIGINAL-COEFFICIENT STIRLING TRIANGLE AUDIT
# ==============================================================================
#
# We now analyze the ORIGINAL coefficient polynomials
#
#     P_k(j) = [u^k] A(j,u)
#     P_k(j) = [u^k] B(j,u)
#
# directly in the falling-factorial basis
#
#     j_under_r = j(j-1)...(j-r+1).
#
# Because Experiment 84 established
#
#     j_under_k | P_k(j),
#
# we expect the exact expansion to have
#
#     P_k(j) = sum_{r=k}^{m} c[k,r] j_under_r.
#
# This is the natural triangular representation.
#
# The experiment checks:
#
#   1. exact falling-basis reconstruction;
#   2. triangular support;
#   3. diagonal coefficients;
#   4. near-diagonal coefficients;
#   5. normalized diagonal ratios;
#   6. Stirling-number comparisons;
#   7. row/column gcd structure;
#   8. A/B cross-channel comparison;
#   9. exact finite consistency at a fresh symbolic value.
#
# No extrapolation.
# No floating point.
# No factorization claim.
# ==============================================================================


j = sp.Symbol("j")
u = sp.Symbol("u")

QQ = sp.QQ


# ------------------------------------------------------------------------------
# Helpers
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
    return -sp.oo if p.is_zero else p.degree()


def coeff_u(expr, k):
    return clean(
        poly(expr, u).nth(k)
    )


def falling(r):
    if r == 0:
        return sp.Integer(1)

    result = sp.Integer(1)

    for s in range(r):
        result *= j - s

    return clean(result)


def rational_gcd(values):
    values = [
        sp.Rational(v)
        for v in values
        if sp.Rational(v) != 0
    ]

    if not values:
        return sp.Integer(0)

    den_lcm = sp.Integer(1)

    for v in values:
        den_lcm = sp.ilcm(
            den_lcm,
            int(v.q),
        )

    ints = [
        int(v * den_lcm)
        for v in values
    ]

    g = 0

    for v in ints:
        g = sp.igcd(
            g,
            abs(v),
        )

    return sp.Rational(
        g,
        den_lcm,
    )


def primitive_signature(values):
    values = [
        sp.Rational(v)
        for v in values
    ]

    if not values:
        return []

    den_lcm = sp.Integer(1)

    for v in values:
        den_lcm = sp.ilcm(
            den_lcm,
            int(v.q),
        )

    ints = [
        int(v * den_lcm)
        for v in values
    ]

    g = 0

    for v in ints:
        g = sp.igcd(
            g,
            abs(v),
        )

    if g == 0:
        return ints

    return [
        v // g
        for v in ints
    ]


def falling_basis_coefficients(expr):
    """
    Exact expansion

        expr = sum_r c_r * j_under_r

    obtained by triangular extraction from highest degree downward.
    """
    P = poly(expr, j)

    if P.is_zero:
        return []

    max_degree = P.degree()
    remainder = clean(P.as_expr())

    coeffs = [
        sp.Integer(0)
        for _ in range(max_degree + 1)
    ]

    for r in range(max_degree, -1, -1):
        current = poly(
            remainder,
            j,
        )

        if current.is_zero:
            coeff = sp.Integer(0)
        else:
            basis = poly(
                falling(r),
                j,
            )

            coeff = clean(
                current.LC()
                / basis.LC()
            )

        coeffs[r] = coeff

        remainder = clean(
            remainder
            - coeff * falling(r)
        )

    if clean(remainder) != 0:
        raise RuntimeError(
            "Falling-factorial reconstruction failed."
        )

    return coeffs


def reconstruct(coeffs):
    result = sp.Integer(0)

    for r, coeff in enumerate(coeffs):
        result += (
            sp.Rational(coeff)
            * falling(r)
        )

    return clean(result)


def finite_difference(values):
    rows = [
        [
            clean(v)
            for v in values
        ]
    ]

    while len(rows[-1]) > 1:
        previous = rows[-1]

        rows.append([
            clean(
                previous[i + 1]
                - previous[i]
            )
            for i in range(
                len(previous) - 1
            )
        ])

    return rows


def newton_coefficients(values):
    table = finite_difference(values)

    return [
        clean(row[0])
        for row in table
    ]


def matrix_from_rows(rows):
    width = max(
        len(row)
        for row in rows
    )

    return sp.Matrix([
        [
            sp.Rational(v)
            for v in row
        ]
        + [
            sp.Integer(0)
        ] * (
            width - len(row)
        )
        for row in rows
    ])


# ------------------------------------------------------------------------------
# Exact source channel polynomials
# ------------------------------------------------------------------------------

A = {
    0: -2,

    1: -2 * (
        275*u + 78
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
        125*u + 46
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
# Build original coefficient polynomials P_k(j)
# ------------------------------------------------------------------------------

def build_coefficients(channel):
    max_k = max(
        degree(expr, u)
        for expr in channel.values()
    )

    max_j = max(channel)

    result = {}

    for k in range(max_k + 1):

        values = [
            coeff_u(
                channel[j0],
                k,
            )
            for j0 in range(max_j + 1)
        ]

        # Exact interpolation on all known j.
        P = sp.interpolate(
            [
                (
                    sp.Integer(j0),
                    values[j0],
                )
                for j0 in range(max_j + 1)
            ],
            j,
        )

        result[k] = clean(P)

    return result


# ------------------------------------------------------------------------------
# Full triangular falling-basis matrices
# ------------------------------------------------------------------------------

def build_falling_matrix(coefficients):
    rows = {}

    max_degree = max(
        degree(P, j)
        for P in coefficients.values()
    )

    for k, P in coefficients.items():

        full = falling_basis_coefficients(P)

        padded = (
            full
            + [
                sp.Integer(0)
            ] * (
                max_degree + 1
                - len(full)
            )
        )

        rows[k] = padded

    return rows


def verify_reconstruction(
    coefficients,
    matrix_rows,
):
    passed = True

    for k, P in coefficients.items():

        rebuilt = reconstruct(
            matrix_rows[k]
        )

        ok = (
            clean(rebuilt)
            ==
            clean(P)
        )

        print(
            f"    k={k}: {ok}"
        )

        passed = passed and ok

    return passed


def verify_triangularity(matrix_rows):
    passed = True

    for k, row in matrix_rows.items():

        bad = []

        for r in range(
            min(k, len(row))
        ):
            if row[r] != 0:
                bad.append(r)

        ok = (
            len(bad) == 0
        )

        print(
            f"    k={k}: "
            f"entries r<k zero = {ok}"
        )

        if not ok:
            print(
                f"      nonzero positions={bad}"
            )

        passed = passed and ok

    return passed


# ------------------------------------------------------------------------------
# Stirling comparison
# ------------------------------------------------------------------------------

def stirling_first_signed(n, k):
    return sp.functions.combinatorial.numbers.stirling(
        n,
        k,
        kind=1,
        signed=True,
    )


def stirling_second(n, k):
    return sp.functions.combinatorial.numbers.stirling(
        n,
        k,
        kind=2,
    )


def stirling_normalization(matrix_rows):
    """
    Compare each falling-basis coefficient c[k,r] with
    elementary Stirling/binomial quantities.

    We report ratios only when denominator is nonzero.
    """
    print()
    print(
        "  STIRLING-SECOND-KIND COMPARISON"
    )

    for k, row in matrix_rows.items():

        ratios = []

        for r, value in enumerate(row):

            if value == 0:
                continue

            S2 = stirling_second(
                r,
                k,
            )

            if S2 != 0:
                ratios.append(
                    (
                        r,
                        clean(value / S2)
                    )
                )

        print(
            f"    k={k}: {ratios}"
        )


# ------------------------------------------------------------------------------
# Diagonal / near diagonal structure
# ------------------------------------------------------------------------------

def diagonal_audit(name, matrix_rows):
    print()
    print(
        f"  {name} DIAGONAL / NEAR-DIAGONAL"
    )

    for k, row in matrix_rows.items():

        if k < len(row):
            diagonal = clean(
                row[k]
            )
        else:
            diagonal = sp.Integer(0)

        entries = []

        for delta in range(0, 4):

            r = k + delta

            if r < len(row):

                entries.append(
                    (
                        r,
                        clean(row[r])
                    )
                )

        print(
            f"    k={k}: "
            f"diagonal={diagonal} "
            f"near={entries}"
        )


# ------------------------------------------------------------------------------
# Row normalization by diagonal element
# ------------------------------------------------------------------------------

def diagonal_normalization(
    name,
    matrix_rows,
):
    print()
    print(
        f"  {name} DIAGONAL-NORMALIZED ROWS"
    )

    for k, row in matrix_rows.items():

        if k >= len(row):
            continue

        diagonal = clean(
            row[k]
        )

        if diagonal == 0:
            print(
                f"    k={k}: zero diagonal"
            )
            continue

        normalized = []

        for r in range(
            k,
            len(row)
        ):
            normalized.append(
                (
                    r,
                    clean(
                        row[r]
                        / diagonal
                    )
                )
            )

        print(
            f"    k={k}: "
            f"{normalized}"
        )


# ------------------------------------------------------------------------------
# Cross-channel triangular comparison
# ------------------------------------------------------------------------------

def cross_channel(
    A_rows,
    B_rows,
):
    print()
    print(
        "=" * 78
    )
    print(
        "3. CROSS-CHANNEL FALLING-BASIS TRIANGLE"
    )
    print(
        "=" * 78
    )

    shared = sorted(
        set(A_rows)
        &
        set(B_rows)
    )

    for k in shared:

        Arow = A_rows[k]
        Brow = B_rows[k]

        ratios = []

        for r in range(
            min(
                len(Arow),
                len(Brow),
            )
        ):

            av = Arow[r]
            bv = Brow[r]

            if av != 0 and bv != 0:
                ratios.append(
                    (
                        r,
                        clean(av / bv)
                    )
                )

        print(
            f"  k={k}: {ratios}"
        )


# ------------------------------------------------------------------------------
# Fresh symbolic consistency
# ------------------------------------------------------------------------------

def fresh_test(
    A_coeff,
    B_coeff,
    A_rows,
    B_rows,
):
    print()
    print(
        "=" * 78
    )
    print(
        "4. FRESH EXACT SYMBOLIC CONSISTENCY"
    )
    print(
        "=" * 78
    )

    fresh_j_values = [
        -2,
        -1,
        0,
        1,
        2,
        3,
        4,
        5,
        8,
    ]

    A_ok = True
    B_ok = True

    print(
        "  A channel"
    )

    for k, P in A_coeff.items():

        rebuilt = reconstruct(
            A_rows[k]
        )

        for j0 in fresh_j_values:

            direct = clean(
                P.subs(
                    j,
                    j0,
                )
            )

            test = clean(
                rebuilt.subs(
                    j,
                    j0,
                )
            )

            if direct != test:
                A_ok = False
                break

        print(
            f"    k={k}: {A_ok}"
        )

    print()
    print(
        "  B channel"
    )

    for k, P in B_coeff.items():

        rebuilt = reconstruct(
            B_rows[k]
        )

        for j0 in fresh_j_values:

            direct = clean(
                P.subs(
                    j,
                    j0,
                )
            )

            test = clean(
                rebuilt.subs(
                    j,
                    j0,
                )
            )

            if direct != test:
                B_ok = False
                break

        print(
            f"    k={k}: {B_ok}"
        )

    return A_ok, B_ok


# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 88 — EXACT ORIGINAL-COEFFICIENT "
        "STIRLING TRIANGLE AUDIT"
    )
    print("=" * 78)

    print()
    print("=" * 78)
    print(
        "0. EXACT SYMBOLIC SETUP"
    )
    print("=" * 78)

    print(
        "  P_A,k(j) = [u^k] A(j,u)"
    )
    print(
        "  P_B,k(j) = [u^k] B(j,u)"
    )
    print(
        "  basis = j_under_r"
    )
    print(
        "  arithmetic = exact QQ"
    )
    print(
        "  floating point = forbidden"
    )

    # --------------------------------------------------------------------------
    # A
    # --------------------------------------------------------------------------

    print()
    print("-" * 78)
    print(
        "A CHANNEL"
    )
    print("-" * 78)

    A_coeff = build_coefficients(A)
    A_rows = build_falling_matrix(
        A_coeff
    )

    print()
    print(
        "1. EXACT FALLING-BASIS RECONSTRUCTION"
    )

    A_reconstruction = verify_reconstruction(
        A_coeff,
        A_rows,
    )

    print(
        f"  ALL A RECONSTRUCTIONS = "
        f"{A_reconstruction}"
    )

    print()
    print(
        "2. TRIANGULARITY AUDIT"
    )

    A_triangular = verify_triangularity(
        A_rows
    )

    print(
        f"  A triangular = "
        f"{A_triangular}"
    )

    diagonal_audit(
        "A",
        A_rows,
    )

    diagonal_normalization(
        "A",
        A_rows,
    )

    print()
    print(
        "  A FALLING-BASIS MATRIX"
    )

    A_matrix = matrix_from_rows(
        [
            A_rows[k]
            for k in sorted(A_rows)
        ]
    )

    print(
        A_matrix
    )

    print()
    print(
        f"  shape={A_matrix.shape}"
    )

    print(
        f"  rank={A_matrix.rank()}"
    )

    print(
        "  column gcds:"
    )

    for c in range(
        A_matrix.cols
    ):

        print(
            f"    c={c}: "
            f"{rational_gcd(A_matrix[:, c])}"
        )

    print()
    print(
        "  primitive row signatures:"
    )

    for k, row in A_rows.items():
        print(
            f"    k={k}: "
            f"{primitive_signature(row)}"
        )

    # --------------------------------------------------------------------------
    # B
    # --------------------------------------------------------------------------

    print()
    print("-" * 78)
    print(
        "B CHANNEL"
    )
    print("-" * 78)

    B_coeff = build_coefficients(B)
    B_rows = build_falling_matrix(
        B_coeff
    )

    print()
    print(
        "1. EXACT FALLING-BASIS RECONSTRUCTION"
    )

    B_reconstruction = verify_reconstruction(
        B_coeff,
        B_rows,
    )

    print(
        f"  ALL B RECONSTRUCTIONS = "
        f"{B_reconstruction}"
    )

    print()
    print(
        "2. TRIANGULARITY AUDIT"
    )

    B_triangular = verify_triangularity(
        B_rows
    )

    print(
        f"  B triangular = "
        f"{B_triangular}"
    )

    diagonal_audit(
        "B",
        B_rows,
    )

    diagonal_normalization(
        "B",
        B_rows,
    )

    print()
    print(
        "  B FALLING-BASIS MATRIX"
    )

    B_matrix = matrix_from_rows(
        [
            B_rows[k]
            for k in sorted(B_rows)
        ]
    )

    print(
        B_matrix
    )

    print()
    print(
        f"  shape={B_matrix.shape}"
    )

    print(
        f"  rank={B_matrix.rank()}"
    )

    print(
        "  column gcds:"
    )

    for c in range(
        B_matrix.cols
    ):

        print(
            f"    c={c}: "
            f"{rational_gcd(B_matrix[:, c])}"
        )

    print()
    print(
        "  primitive row signatures:"
    )

    for k, row in B_rows.items():
        print(
            f"    k={k}: "
            f"{primitive_signature(row)}"
        )

    # --------------------------------------------------------------------------
    # Stirling comparison
    # --------------------------------------------------------------------------

    print()
    print(
        "5. STIRLING-BASIS COMPARISON"
    )

    print()
    print(
        "  A channel"
    )

    stirling_normalization(
        A_rows
    )

    print()
    print(
        "  B channel"
    )

    stirling_normalization(
        B_rows
    )

    # --------------------------------------------------------------------------
    # Cross-channel
    # --------------------------------------------------------------------------

    cross_channel(
        A_rows,
        B_rows,
    )

    # --------------------------------------------------------------------------
    # Fresh checks
    # --------------------------------------------------------------------------

    A_fresh, B_fresh = fresh_test(
        A_coeff,
        B_coeff,
        A_rows,
        B_rows,
    )

    # --------------------------------------------------------------------------
    # Interpretation
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "5. STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
        """
  The important object is now the ORIGINAL coefficient triangle

      P_k(j) = [u^k].

  Instead of dividing by the falling factorial and then studying the
  quotient, we directly expand

      P_k(j)
        = sum_{r} C[k,r] j_under_r.

  Because

      j_under_k | P_k(j),

  the matrix should be triangular in the observed range:

      C[k,r] = 0 for r < k.

  This is more informative than the quotient basis because it keeps the
  combinatorial divisibility and the remaining coefficient structure in
  the same object.

  The main tests are:

    * diagonal law:
          C[k,k]

    * near-diagonal laws:
          C[k,k+1], C[k,k+2], ...

    * column normalization;

    * Stirling-number comparisons;

    * cross-channel A/B ratios.

  A particularly interesting result would be that the diagonal or
  near-diagonal entries follow a simple product formula.

  Another strong signal would be a relation of the form

      C[k,r]
        = f(k,r) * binomial(r,k)

  or

      C[k,r]
        = f(r-k) * g(k),

  or a short expression in Stirling numbers.

  That would expose a genuinely combinatorial construction law rather
  than merely an interpolation identity.

  This experiment remains finite and exact.

  No extrapolation is claimed.

  No factorization algorithm is inferred.
        """
    )

    checks = {
        "A_falling_basis": A_reconstruction,
        "A_triangular": A_triangular,
        "B_falling_basis": B_reconstruction,
        "B_triangular": B_triangular,
        "fresh_A": A_fresh,
        "fresh_B": B_fresh,
    }

    failures = sum(
        1
        for value in checks.values()
        if not value
    )

    print()
    print("=" * 78)
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
    print("=" * 78)
    print(
        "EXPERIMENT 88 COMPLETE"
    )
    print("=" * 78)


if __name__ == "__main__":
    main()

