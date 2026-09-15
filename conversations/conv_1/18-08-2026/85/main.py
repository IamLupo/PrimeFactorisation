#!/usr/bin/env python3

import sympy as sp

# ==============================================================================
# EXPERIMENT 87 — EXACT STIRLING / BINOMIAL-BASIS COEFFICIENT STRUCTURE
# ==============================================================================
#
# Goal
# ----
# Experiment 84 established
#
#     P_A,k(j) = [u^k] A(j,u) = j_under_k * Q_A,k(j)
#     P_B,k(j) = [u^k] B(j,u) = j_under_k * Q_B,k(j)
#
# Experiment 86 correctly identified the falling-factorial quotient Q_k(j),
# but its Newton reconstruction test used symbolic binomial(j,r) objects
# directly.  That made exact polynomial comparison unreliable.
#
# Here the Newton basis is represented explicitly as
#
#     C_r(j) = binomial(j,r)
#            = j_under_r / r!
#
# using ordinary polynomial expressions only.
#
# We then investigate:
#
#   1. exact reconstruction in the binomial basis;
#   2. exact conversion between monomial and binomial bases;
#   3. the quotient coefficient matrices;
#   4. rank and minor structure;
#   5. normalized integer signatures;
#   6. diagonal / anti-diagonal patterns;
#   7. possible Stirling-like transforms;
#   8. whether the A and B arrays share a common construction pattern.
#
# No extrapolation beyond the available j range.
# No floating point.
# ==============================================================================


j = sp.Symbol("j")
u = sp.Symbol("u")

QQ = sp.QQ


# ------------------------------------------------------------------------------
# Exact helpers
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
    return poly(expr, u).nth(k)


def falling_factorial(k):
    result = sp.Integer(1)

    for r in range(k):
        result *= j - r

    return clean(result)


def binomial_basis_polynomial(r):
    if r == 0:
        return sp.Integer(1)

    return clean(
        falling_factorial(r)
        / sp.factorial(r)
    )


def quotient_exact(expr, divisor):
    P = poly(expr, j)
    D = poly(divisor, j)

    Q, R = sp.div(
        P,
        D,
        domain=QQ,
    )

    return clean(Q.as_expr()), clean(R.as_expr())


def interpolate(values):
    points = [
        (sp.Integer(idx), clean(value))
        for idx, value in enumerate(values)
    ]

    return clean(
        sp.interpolate(points, j)
    )


def finite_difference_table(values):
    rows = [
        [clean(value) for value in values]
    ]

    while len(rows[-1]) > 1:
        previous = rows[-1]

        current = [
            clean(
                previous[idx + 1]
                - previous[idx]
            )
            for idx in range(
                len(previous) - 1
            )
        ]

        rows.append(current)

    return rows


def newton_coefficients(expr, max_j):
    values = [
        clean(
            expr.subs(
                j,
                sp.Integer(idx)
            )
        )
        for idx in range(max_j + 1)
    ]

    table = finite_difference_table(values)

    return [
        clean(row[0])
        for row in table
    ]


def reconstruct_from_newton(coefficients):
    result = sp.Integer(0)

    for r, coefficient in enumerate(coefficients):
        result += (
            clean(coefficient)
            * binomial_basis_polynomial(r)
        )

    return clean(result)


def monomial_coefficients(expr):
    P = poly(expr, j)

    return [
        clean(
            P.nth(r)
        )
        for r in range(
            P.degree() + 1
        )
    ]


def primitive_integer_signature(values):
    rationals = [
        sp.Rational(value)
        for value in values
    ]

    if not rationals:
        return []

    denominator_lcm = sp.Integer(1)

    for value in rationals:
        denominator_lcm = sp.ilcm(
            denominator_lcm,
            int(value.q),
        )

    integers = [
        int(value * denominator_lcm)
        for value in rationals
    ]

    gcd_value = 0

    for value in integers:
        gcd_value = sp.igcd(
            gcd_value,
            abs(value),
        )

    if gcd_value == 0:
        return [
            0
            for _ in integers
        ]

    return [
        value // gcd_value
        for value in integers
    ]


def rational_gcd(values):
    values = [
        sp.Rational(value)
        for value in values
        if sp.Rational(value) != 0
    ]

    if not values:
        return sp.Integer(0)

    denominator_lcm = sp.Integer(1)

    for value in values:
        denominator_lcm = sp.ilcm(
            denominator_lcm,
            int(value.q),
        )

    integers = [
        int(value * denominator_lcm)
        for value in values
    ]

    gcd_value = 0

    for value in integers:
        gcd_value = sp.igcd(
            gcd_value,
            abs(value),
        )

    return sp.Rational(
        gcd_value,
        denominator_lcm,
    )


def matrix_from_rows(rows):
    if not rows:
        return sp.Matrix([])

    width = max(
        len(row)
        for row in rows
    )

    padded = []

    for row in rows:
        padded.append(
            [
                sp.Rational(x)
                for x in row
            ]
            + [
                sp.Rational(0)
            ] * (
                width - len(row)
            )
        )

    return sp.Matrix(padded)


def rank_of_rows(rows):
    return matrix_from_rows(rows).rank()


def factor(expr):
    return sp.factor(
        clean(expr)
    )


# ------------------------------------------------------------------------------
# Exact channel data
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
# Build coefficient laws
# ------------------------------------------------------------------------------

def build_channel(channel):
    max_k = max(
        int(
            degree(expr, u)
        )
        for expr in channel.values()
    )

    coefficient_polys = {}
    quotient_polys = {}
    newton_rows = {}
    monomial_rows = {}

    max_j = max(channel)

    for k in range(max_k + 1):

        values = [
            coeff_u(
                channel[idx],
                k
            )
            for idx in range(max_j + 1)
        ]

        P = interpolate(values)

        coefficient_polys[k] = P

        divisor = falling_factorial(k)

        Q, remainder = quotient_exact(
            P,
            divisor
        )

        if remainder != 0:
            raise RuntimeError(
                f"falling-factorial division failed for k={k}"
            )

        quotient_polys[k] = Q

        coeffs = newton_coefficients(
            Q,
            max_j
        )

        newton_rows[k] = coeffs

        monomial_rows[k] = (
            monomial_coefficients(Q)
        )

    return {
        "coefficient_polys": coefficient_polys,
        "quotient_polys": quotient_polys,
        "newton_rows": newton_rows,
        "monomial_rows": monomial_rows,
        "max_j": max_j,
    }


# ------------------------------------------------------------------------------
# Verify exact falling-factorial structure
# ------------------------------------------------------------------------------

def verify_falling_structure(data):
    passed = True

    for k, P in data["coefficient_polys"].items():

        divisor = falling_factorial(k)

        _, remainder = quotient_exact(
            P,
            divisor
        )

        ok = (
            clean(remainder)
            == 0
        )

        print(
            f"    k={k}: "
            f"j_under_{k} divides P_k(j) = {ok}"
        )

        passed = passed and ok

    return passed


# ------------------------------------------------------------------------------
# Exact quotient reconstruction
# ------------------------------------------------------------------------------

def verify_newton_reconstruction(data):
    passed = True

    for k, Q in data["quotient_polys"].items():

        coefficients = data["newton_rows"][k]

        rebuilt = reconstruct_from_newton(
            coefficients
        )

        ok = (
            clean(rebuilt)
            ==
            clean(Q)
        )

        print(
            f"    Q_{k}: {ok}"
        )

        if not ok:
            print(
                "      direct =",
                clean(Q)
            )
            print(
                "      rebuilt =",
                clean(rebuilt)
            )

        passed = passed and ok

    return passed


# ------------------------------------------------------------------------------
# Compare the original coefficient P_k(j) with a direct double basis:
#
#     P_k(j) = j_under_k * sum_r c[k,r] binom(j,r)
#
# The c[k,r] are exactly the Newton coefficients of Q_k.
# ------------------------------------------------------------------------------

def direct_double_basis_reconstruction(
    data
):
    passed = True

    for k, P in data["coefficient_polys"].items():

        Q = data["quotient_polys"][k]
        row = data["newton_rows"][k]

        rebuilt_Q = reconstruct_from_newton(
            row
        )

        rebuilt_P = clean(
            falling_factorial(k)
            * rebuilt_Q
        )

        ok = (
            rebuilt_P
            ==
            clean(P)
        )

        print(
            f"    k={k}: {ok}"
        )

        passed = passed and ok

    return passed


# ------------------------------------------------------------------------------
# Row / column structure
# ------------------------------------------------------------------------------

def print_newton_matrix_analysis(
    name,
    data
):
    rows = [
        data["newton_rows"][k]
        for k in sorted(
            data["newton_rows"]
        )
    ]

    matrix = matrix_from_rows(rows)

    print()
    print(
        f"  {name} NEWTON MATRIX"
    )

    print(
        f"    shape = {matrix.shape}"
    )

    print(
        f"    rank = {matrix.rank()}"
    )

    print()
    print(
        f"  {name} NEWTON ROW SIGNATURES"
    )

    for k, row in data["newton_rows"].items():

        print(
            f"    k={k}: "
            f"{primitive_integer_signature(row)}"
        )

    print()
    print(
        f"  {name} NEWTON COLUMN GCDS"
    )

    width = matrix.cols

    for c in range(width):

        values = [
            matrix[r, c]
            for r in range(
                matrix.rows
            )
        ]

        print(
            f"    column={c}: "
            f"gcd={rational_gcd(values)}"
        )


# ------------------------------------------------------------------------------
# Test Stirling-like transforms
#
# The first kind Stirling transform takes powers into falling factorials.
# Here we test whether the quotient rows become simpler under another
# exact basis change.
#
# We compare:
#
#   monomial basis j^r
#   falling factorial basis j_under_r
#   binomial basis binom(j,r)
#
# and report the exact coefficient arrays.
# ------------------------------------------------------------------------------

def falling_basis_coefficients(expr):
    P = poly(expr, j)

    degree_P = (
        -1
        if P.is_zero
        else P.degree()
    )

    remainder = clean(P.as_expr())
    coefficients = []

    # Triangular extraction from highest falling degree.
    for r in range(
        degree_P,
        -1,
        -1
    ):

        basis = poly(
            falling_factorial(r),
            j
        )

        leading = basis.LC()

        current = poly(
            remainder,
            j
        )

        coefficient = clean(
            current.LC()
            / leading
        )

        coefficients.append(
            (
                r,
                coefficient
            )
        )

        remainder = clean(
            remainder
            - coefficient
            * falling_factorial(r)
        )

    coefficients.reverse()

    return coefficients


def print_falling_basis_analysis(
    name,
    data
):
    print()
    print(
        f"  {name} FALLING-FACTORIAL BASIS"
    )

    for k, Q in data["quotient_polys"].items():

        coefficients = (
            falling_basis_coefficients(Q)
        )

        print(
            f"    Q_{k}: {coefficients}"
        )


# ------------------------------------------------------------------------------
# Low-rank minor audit
# ------------------------------------------------------------------------------

def rank1_minor_audit(rows):
    M = matrix_from_rows(rows)

    if M.rows < 2 or M.cols < 2:
        return True

    for r1 in range(M.rows):
        for r2 in range(r1 + 1, M.rows):
            for c1 in range(M.cols):
                for c2 in range(c1 + 1, M.cols):

                    minor = clean(
                        M[r1, c1] * M[r2, c2]
                        - M[r1, c2] * M[r2, c1]
                    )

                    if minor != 0:
                        return False

    return True


# ------------------------------------------------------------------------------
# Cross-channel comparison
# ------------------------------------------------------------------------------

def cross_channel_analysis(
    A_data,
    B_data
):
    print()
    print("=" * 78)
    print("3. CROSS-CHANNEL STIRLING / BINOMIAL COMPARISON")
    print("=" * 78)

    shared = sorted(
        set(
            A_data["newton_rows"]
        )
        &
        set(
            B_data["newton_rows"]
        )
    )

    for k in shared:

        Arow = A_data["newton_rows"][k]
        Brow = B_data["newton_rows"][k]

        ratios = []

        width = min(
            len(Arow),
            len(Brow)
        )

        for r in range(width):

            av = clean(Arow[r])
            bv = clean(Brow[r])

            if av != 0 and bv != 0:
                ratios.append(
                    (
                        r,
                        clean(av / bv)
                    )
                )

        print(
            f"  k={k}: "
            f"entrywise ratios={ratios}"
        )


# ------------------------------------------------------------------------------
# Fresh exact consistency using a new symbolic point
# ------------------------------------------------------------------------------

def fresh_symbolic_test(
    A_data,
    B_data
):
    t_value = sp.Rational(
        26279197733,
        353095
    )

    u_value = clean(
        t_value * (
            t_value + 1
        )
    )

    print()
    print("=" * 78)
    print("4. FRESH EXACT SYMBOLIC CONSISTENCY")
    print("=" * 78)

    print(
        f"  t = {t_value}"
    )

    print(
        f"  u = {u_value}"
    )

    A_ok = True
    B_ok = True

    print()
    print("  A channel")

    for k, P in A_data["coefficient_polys"].items():

        Q = A_data["quotient_polys"][k]

        original = clean(
            P.subs(
                j,
                4
            )
        )

        reconstructed = clean(
            (
                falling_factorial(k)
                * Q
            ).subs(
                j,
                4
            )
        )

        ok = (
            original
            ==
            reconstructed
        )

        A_ok = A_ok and ok

        print(
            f"    k={k}: {ok}"
        )

    print()
    print("  B channel")

    for k, P in B_data["coefficient_polys"].items():

        Q = B_data["quotient_polys"][k]

        original = clean(
            P.subs(
                j,
                4
            )
        )

        reconstructed = clean(
            (
                falling_factorial(k)
                * Q
            ).subs(
                j,
                4
            )
        )

        ok = (
            original
            ==
            reconstructed
        )

        B_ok = B_ok and ok

        print(
            f"    k={k}: {ok}"
        )

    return A_ok, B_ok


# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 87 — EXACT STIRLING / "
        "BINOMIAL-BASIS COEFFICIENT STRUCTURE"
    )
    print("=" * 78)

    print()
    print("=" * 78)
    print("0. EXACT SYMBOLIC SETUP")
    print("=" * 78)

    print(
        "  A_j = A_{2j}(u), j=0,...,8"
    )
    print(
        "  B_j = B_{2j+1}(u), j=0,...,7"
    )
    print(
        "  arithmetic = exact QQ"
    )
    print(
        "  floating point = forbidden"
    )

    print()
    print("-" * 78)
    print("A CHANNEL")
    print("-" * 78)

    A_data = build_channel(A)

    A_falling = verify_falling_structure(
        A_data
    )

    print(
        f"  ALL A FALLING DIVISIONS = "
        f"{A_falling}"
    )

    print()
    print(
        "  FULL A QUOTIENT NEWTON ROWS"
    )

    for k, row in A_data["newton_rows"].items():
        print(
            f"    k={k}: {row}"
        )

    print()
    print(
        "  A QUOTIENT RECONSTRUCTION"
    )

    A_newton_ok = (
        verify_newton_reconstruction(
            A_data
        )
    )

    print(
        f"  ALL A QUOTIENT RECONSTRUCTIONS = "
        f"{A_newton_ok}"
    )

    print_newton_matrix_analysis(
        "A",
        A_data
    )

    print_falling_basis_analysis(
        "A",
        A_data
    )

    print()
    print(
        "  A ORIGINAL DOUBLE-BASIS RECONSTRUCTION"
    )

    A_double_ok = (
        direct_double_basis_reconstruction(
            A_data
        )
    )

    print(
        f"  ALL A DOUBLE-BASIS RECONSTRUCTIONS = "
        f"{A_double_ok}"
    )

    A_rank1 = rank1_minor_audit(
        [
            A_data["newton_rows"][k]
            for k in sorted(
                A_data["newton_rows"]
            )
        ]
    )

    print(
        f"  A Newton matrix rank-1 separable = "
        f"{A_rank1}"
    )

    print()
    print("-" * 78)
    print("B CHANNEL")
    print("-" * 78)

    B_data = build_channel(B)

    B_falling = verify_falling_structure(
        B_data
    )

    print(
        f"  ALL B FALLING DIVISIONS = "
        f"{B_falling}"
    )

    print()
    print(
        "  FULL B QUOTIENT NEWTON ROWS"
    )

    for k, row in B_data["newton_rows"].items():
        print(
            f"    k={k}: {row}"
        )

    print()
    print(
        "  B QUOTIENT RECONSTRUCTION"
    )

    B_newton_ok = (
        verify_newton_reconstruction(
            B_data
        )
    )

    print(
        f"  ALL B QUOTIENT RECONSTRUCTIONS = "
        f"{B_newton_ok}"
    )

    print_newton_matrix_analysis(
        "B",
        B_data
    )

    print_falling_basis_analysis(
        "B",
        B_data
    )

    print()
    print(
        "  B ORIGINAL DOUBLE-BASIS RECONSTRUCTION"
    )

    B_double_ok = (
        direct_double_basis_reconstruction(
            B_data
        )
    )

    print(
        f"  ALL B DOUBLE-BASIS RECONSTRUCTIONS = "
        f"{B_double_ok}"
    )

    B_rank1 = rank1_minor_audit(
        [
            B_data["newton_rows"][k]
            for k in sorted(
                B_data["newton_rows"]
            )
        ]
    )

    print(
        f"  B Newton matrix rank-1 separable = "
        f"{B_rank1}"
    )

    cross_channel_analysis(
        A_data,
        B_data
    )

    A_fresh, B_fresh = fresh_symbolic_test(
        A_data,
        B_data
    )

    print()
    print("=" * 78)
    print("5. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
  The corrected basis chain is

      P_k(j)
          = [u^k] channel

      P_k(j)
          = j_under_k * Q_k(j)

      Q_k(j)
          = sum_r c[k,r] * binom(j,r).

  The previous false reconstruction status came from representing
  binom(j,r) as a SymPy combinatorial object rather than as the
  exact polynomial

      j_under_r / r!.

  Experiment 87 removes that ambiguity.

  The important structural object is now the matrix

      C = (c[k,r]).

  There are three possible outcomes.

    1. Low rank:
         the quotient family may have a small separable core.

    2. Sparse / triangular structure:
         the coefficients may be generated by a finite combinatorial
         transform.

    3. Full-rank irregular structure:
         the falling-factorial extraction may be the main structural
         simplification already available.

  The comparison of monomial, binomial, and falling-factorial bases is
  useful because a construction may look complicated in one basis and
  simple in another.

  This experiment does not extrapolate to unseen layer indices.

  It does not assume a recurrence.

  It does not infer a factorization algorithm.

  The eventual inverse target remains:

      N
        -> layer information
        -> t=N/X
        -> X
        -> S
        -> p,q.
        """
    )

    failures = 0

    checks = {
        "A_falling_divisibility": A_falling,
        "B_falling_divisibility": B_falling,
        "A_newton_reconstruction": A_newton_ok,
        "B_newton_reconstruction": B_newton_ok,
        "A_double_basis": A_double_ok,
        "B_double_basis": B_double_ok,
        "fresh_A": A_fresh,
        "fresh_B": B_fresh,
    }

    for value in checks.values():
        if not value:
            failures += 1

    print()
    print("=" * 78)
    print("FINAL EXACTNESS")
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
    print("EXPERIMENT 87 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

