#!/usr/bin/env python3

import sympy as sp

# ==============================================================================
# EXPERIMENT 85 — EXACT FALLING-FACTORIAL COEFFICIENT-MATRIX STRUCTURE
# ==============================================================================
#
# Goal:
#   Experiment 84 established
#
#       [u^k] A(j,u) = j_under_k * Q_A(k,j)
#       [u^k] B(j,u) = j_under_k * Q_B(k,j)
#
#   where j_under_k = j(j-1)...(j-k+1).
#
# The quotients Q_A(k,j), Q_B(k,j) still looked complicated.
#
# This experiment changes basis again:
#
#       Q(k,j) = sum_r c[k,r] * binom(j,r)
#
# and studies the exact coefficient matrices c[k,r].
#
# We test:
#
#   1. triangularity / vanishing regions;
#   2. row and column gcd structure;
#   3. diagonal and near-diagonal patterns;
#   4. whether c[k,r] factor into simple functions of k and r;
#   5. whether A/B coefficient matrices share a hidden common matrix;
#   6. exact reconstruction from the coefficient matrices.
#
# Nothing is extrapolated beyond the observed finite index range.
# No floating point is used.
#
# ==============================================================================

u = sp.Symbol("u")
j = sp.Symbol("j")

QQ = sp.QQ


# ------------------------------------------------------------------------------
# Exact symbolic helpers
# ------------------------------------------------------------------------------

def S(x):
    return sp.sympify(x)


def clean(x):
    return sp.cancel(sp.expand(S(x)))


def poly(expr, var):
    return sp.Poly(clean(expr), var, domain=QQ)


def degree(expr, var):
    p = poly(expr, var)
    return -sp.oo if p.is_zero else p.degree()


def coeff_u(expr, k):
    return poly(expr, u).nth(k)


def falling(k):
    out = sp.Integer(1)

    for r in range(k):
        out *= (j - r)

    return sp.expand(out)


def divide_j(expr, divisor):
    p = poly(expr, j)
    q = poly(divisor, j)

    quotient, remainder = sp.div(
        p,
        q,
        domain=QQ,
    )

    return (
        clean(quotient.as_expr()),
        clean(remainder.as_expr()),
    )


def interpolate(values):
    points = [
        (sp.Integer(idx), S(value))
        for idx, value in enumerate(values)
    ]

    return clean(sp.interpolate(points, j))


def finite_differences(values):
    current = [S(v) for v in values]
    rows = [current]

    while len(current) > 1:
        current = [
            clean(current[i + 1] - current[i])
            for i in range(len(current) - 1)
        ]
        rows.append(current)

    return rows


def newton_coefficients(expr, max_index):
    values = [
        clean(expr.subs(j, sp.Integer(idx)))
        for idx in range(max_index + 1)
    ]

    rows = finite_differences(values)

    coeffs = []

    for r in range(len(rows)):
        coeffs.append(
            clean(rows[r][0])
        )

    return coeffs


def binomial_basis(expr, max_index):
    coeffs = newton_coefficients(
        expr,
        max_index,
    )

    result = sp.Integer(0)

    for r, c in enumerate(coeffs):
        result += c * sp.binomial(j, r)

    return clean(result)


def primitive_integer_vector(values):
    values = [sp.Rational(v) for v in values]

    if not values:
        return []

    denominator_lcm = sp.Integer(1)

    for v in values:
        denominator_lcm = sp.ilcm(
            denominator_lcm,
            int(v.q),
        )

    ints = [
        int(v * denominator_lcm)
        for v in values
    ]

    g = 0

    for x in ints:
        g = sp.igcd(g, abs(x))

    if g == 0:
        return [0 for _ in ints]

    return [
        x // g
        for x in ints
    ]


def rational_gcd(values):
    values = [
        sp.Rational(v)
        for v in values
        if sp.Rational(v) != 0
    ]

    if not values:
        return sp.Integer(0)

    denominator_lcm = sp.Integer(1)

    for v in values:
        denominator_lcm = sp.ilcm(
            denominator_lcm,
            int(v.q),
        )

    ints = [
        int(v * denominator_lcm)
        for v in values
    ]

    g = 0

    for x in ints:
        g = sp.igcd(g, abs(x))

    if g == 0:
        return sp.Integer(0)

    return sp.Rational(g, denominator_lcm)


def factor_if_possible(expr):
    return sp.factor(clean(expr))


# ------------------------------------------------------------------------------
# Channel data
# ------------------------------------------------------------------------------

A = {
    0: sp.Integer(-2),

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
    0: sp.Integer(25),

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
# Build coefficient polynomials in j
# ------------------------------------------------------------------------------

def build_channel_polynomials(channel):
    max_k = max(
        int(degree(expr, u))
        for expr in channel.values()
    )

    coeff_polys = {}

    for k in range(max_k + 1):

        values = [
            coeff_u(
                channel[idx],
                k,
            )
            for idx in range(len(channel))
        ]

        P = interpolate(values)

        coeff_polys[k] = clean(P)

    return coeff_polys


# ------------------------------------------------------------------------------
# Remove falling factorial and construct Q_k(j)
# ------------------------------------------------------------------------------

def build_quotients(channel, coeff_polys):
    quotients = {}

    for k, P in coeff_polys.items():

        divisor = falling(k)

        Q, R = divide_j(
            P,
            divisor,
        )

        if clean(R) != 0:
            raise RuntimeError(
                f"Falling-factorial division failed for k={k}"
            )

        quotients[k] = clean(Q)

    return quotients


# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

def main():

    print("=" * 78)
    print("EXPERIMENT 85 — EXACT FALLING-FACTORIAL COEFFICIENT-MATRIX STRUCTURE")
    print("=" * 78)

    print()
    print("=" * 78)
    print("0. EXACT SYMBOLIC SETUP")
    print("=" * 78)
    print("  A_j = A_{2j}(u), j=0,...,8")
    print("  B_j = B_{2j+1}(u), j=0,...,7")
    print("  arithmetic = exact QQ")
    print("  floating point = forbidden")

    # ------------------------------------------------------------------
    # 1. Build coefficient laws
    # ------------------------------------------------------------------

    A_coeff = build_channel_polynomials(A)
    B_coeff = build_channel_polynomials(B)

    print()
    print("=" * 78)
    print("1. FALLING-FACTORIAL QUOTIENT PREPARATION")
    print("=" * 78)

    A_quot = build_quotients(
        A,
        A_coeff,
    )

    B_quot = build_quotients(
        B,
        B_coeff,
    )

    for k in sorted(A_quot):
        print(
            f"  A [u^{k}]: "
            f"degree={degree(A_coeff[k], j)} "
            f"quotient_degree={degree(A_quot[k], j)}"
        )

    print()

    for k in sorted(B_quot):
        print(
            f"  B [u^{k}]: "
            f"degree={degree(B_coeff[k], j)} "
            f"quotient_degree={degree(B_quot[k], j)}"
        )

    # ------------------------------------------------------------------
    # 2. Convert each quotient to Newton/binomial basis
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. NEWTON / BINOMIAL COEFFICIENT MATRICES")
    print("=" * 78)

    A_matrix = {}
    B_matrix = {}

    for k, Q in A_quot.items():

        coeffs = newton_coefficients(
            Q,
            8,
        )

        A_matrix[k] = coeffs

        print()
        print(
            f"  A [u^{k}] Newton coefficients:"
        )
        print(
            f"    {coeffs}"
        )

    for k, Q in B_quot.items():

        coeffs = newton_coefficients(
            Q,
            7,
        )

        B_matrix[k] = coeffs

        print()
        print(
            f"  B [u^{k}] Newton coefficients:"
        )
        print(
            f"    {coeffs}"
        )

    # ------------------------------------------------------------------
    # 3. Vanishing triangle audit
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. EXACT TRIANGULARITY AUDIT")
    print("=" * 78)

    A_triangular = True
    B_triangular = True

    print()
    print("  A-channel")

    for k, row in A_matrix.items():

        expected_zero = 0

        for r in range(k):
            expected_zero += 1

            ok = clean(row[r]) == 0

            if not ok:
                A_triangular = False

        print(
            f"    k={k}: "
            f"first {k} Newton coefficients zero = "
            f"{all(clean(row[r]) == 0 for r in range(k))}"
        )

    print()
    print("  B-channel")

    for k, row in B_matrix.items():

        ok = all(
            clean(row[r]) == 0
            for r in range(k)
        )

        B_triangular = (
            B_triangular
            and ok
        )

        print(
            f"    k={k}: "
            f"first {k} Newton coefficients zero = {ok}"
        )

    # ------------------------------------------------------------------
    # 4. Reduced matrices after deleting forced zeros
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. REDUCED NEWTON COEFFICIENT ARRAYS")
    print("=" * 78)

    print()
    print("  A reduced rows")

    A_reduced = {}

    for k, row in A_matrix.items():

        reduced = [
            clean(x)
            for x in row[k:]
        ]

        A_reduced[k] = reduced

        print(
            f"    k={k}: {reduced}"
        )

    print()
    print("  B reduced rows")

    B_reduced = {}

    for k, row in B_matrix.items():

        reduced = [
            clean(x)
            for x in row[k:]
        ]

        B_reduced[k] = reduced

        print(
            f"    k={k}: {reduced}"
        )

    # ------------------------------------------------------------------
    # 5. Column gcd structure
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. COLUMN GCD STRUCTURE")
    print("=" * 78)

    max_A_cols = max(
        len(row)
        for row in A_reduced.values()
    )

    max_B_cols = max(
        len(row)
        for row in B_reduced.values()
    )

    print()
    print("  A-channel column gcds")

    A_col_gcd = {}

    for offset in range(max_A_cols):

        vals = []

        for k, row in A_reduced.items():

            if offset < len(row):
                vals.append(
                    row[offset]
                )

        g = rational_gcd(vals)
        A_col_gcd[offset] = g

        print(
            f"    reduced-column={offset}: gcd={g}"
        )

    print()
    print("  B-channel column gcds")

    B_col_gcd = {}

    for offset in range(max_B_cols):

        vals = []

        for k, row in B_reduced.items():

            if offset < len(row):
                vals.append(
                    row[offset]
                )

        g = rational_gcd(vals)
        B_col_gcd[offset] = g

        print(
            f"    reduced-column={offset}: gcd={g}"
        )

    # ------------------------------------------------------------------
    # 6. Diagonal / near-diagonal ratios
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. DIAGONAL AND NEAR-DIAGONAL RATIO AUDIT")
    print("=" * 78)

    print()
    print("  A-channel")

    for k in sorted(A_reduced):

        row = A_reduced[k]

        if not row:
            continue

        lead = row[0]

        ratios = []

        for offset, value in enumerate(row[:4]):

            if value != 0:
                ratios.append(
                    (
                        offset,
                        clean(value / lead),
                    )
                )

        print(
            f"    k={k}: {ratios}"
        )

    print()
    print("  B-channel")

    for k in sorted(B_reduced):

        row = B_reduced[k]

        if not row:
            continue

        lead = row[0]

        ratios = []

        for offset, value in enumerate(row[:4]):

            if value != 0:
                ratios.append(
                    (
                        offset,
                        clean(value / lead),
                    )
                )

        print(
            f"    k={k}: {ratios}"
        )

    # ------------------------------------------------------------------
    # 7. Primitive integer matrix signatures
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. PRIMITIVE INTEGER SIGNATURES OF REDUCED ROWS")
    print("=" * 78)

    print()
    print("  A-channel")

    for k, row in A_reduced.items():

        signature = primitive_integer_vector(
            row
        )

        print(
            f"    k={k}: {signature}"
        )

    print()
    print("  B-channel")

    for k, row in B_reduced.items():

        signature = primitive_integer_vector(
            row
        )

        print(
            f"    k={k}: {signature}"
        )

    # ------------------------------------------------------------------
    # 8. Cross-channel comparison on shared k
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. CROSS-CHANNEL REDUCED-ROW COMPARISON")
    print("=" * 78)

    shared_k = sorted(
        set(A_reduced.keys())
        &
        set(B_reduced.keys())
    )

    cross_row_gcd = []

    for k in shared_k:

        A_row = A_reduced[k]
        B_row = B_reduced[k]

        length = min(
            len(A_row),
            len(B_row),
        )

        print()
        print(f"  k={k}")

        ratios = []

        for offset in range(length):

            a_val = A_row[offset]
            b_val = B_row[offset]

            if a_val != 0 and b_val != 0:

                ratios.append(
                    (
                        offset,
                        clean(a_val / b_val),
                    )
                )

        g = rational_gcd(
            A_row + B_row
        )

        cross_row_gcd.append(g)

        print(
            f"    entrywise ratios={ratios}"
        )
        print(
            f"    combined gcd={g}"
        )

    # ------------------------------------------------------------------
    # 9. Seek simple factorization of the reduced coefficient matrix
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. SIMPLE SEPARABLE-FORM SEARCH")
    print("=" * 78)

    print(
        """
  We test whether reduced coefficients can approximately separate as

      c[k,r] = f(k) g(r)

  exactly over QQ.

  Exact rank 1 would mean every 2x2 minor vanishes.
        """
    )

    def matrix_rank_qq(rows):
        if not rows:
            return 0

        M = sp.Matrix(
            [
                [
                    S(x)
                    for x in row
                ]
                for row in rows
            ]
        )

        return M.rank()

    A_rect = []
    B_rect = []

    for k in sorted(A_reduced):

        row = A_reduced[k]

        if row:
            A_rect.append(row)

    for k in sorted(B_reduced):

        row = B_reduced[k]

        if row:
            B_rect.append(row)

    # Pad with zeros for a common rectangle.

    def pad(rows):
        width = max(
            len(row)
            for row in rows
        )

        return [
            list(row)
            + [sp.Integer(0)] * (
                width - len(row)
            )
            for row in rows
        ]

    A_pad = pad(A_rect)
    B_pad = pad(B_rect)

    A_rank = matrix_rank_qq(A_pad)
    B_rank = matrix_rank_qq(B_pad)

    print(
        f"  A reduced coefficient matrix rank = {A_rank}"
    )
    print(
        f"  B reduced coefficient matrix rank = {B_rank}"
    )

    print(
        f"  A rank-1 separable = {A_rank == 1}"
    )
    print(
        f"  B rank-1 separable = {B_rank == 1}"
    )

    # ------------------------------------------------------------------
    # 10. Exact reconstruction from Newton matrices
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. EXACT RECONSTRUCTION FROM COEFFICIENT MATRICES")
    print("=" * 78)

    A_reconstruction_ok = True
    B_reconstruction_ok = True

    for k, Q in A_quot.items():

        rebuilt = sp.Integer(0)

        for r, coefficient in enumerate(
            A_matrix[k]
        ):
            rebuilt += (
                coefficient
                * sp.binomial(j, r)
            )

        rebuilt = clean(rebuilt)

        ok = (
            rebuilt
            ==
            clean(Q)
        )

        A_reconstruction_ok = (
            A_reconstruction_ok
            and ok
        )

        print(
            f"  A Q_{k}: {ok}"
        )

    for k, Q in B_quot.items():

        rebuilt = sp.Integer(0)

        for r, coefficient in enumerate(
            B_matrix[k]
        ):
            rebuilt += (
                coefficient
                * sp.binomial(j, r)
            )

        rebuilt = clean(rebuilt)

        ok = (
            rebuilt
            ==
            clean(Q)
        )

        B_reconstruction_ok = (
            B_reconstruction_ok
            and ok
        )

        print(
            f"  B Q_{k}: {ok}"
        )

    # ------------------------------------------------------------------
    # 11. Fresh exact semiprime consistency
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("11. FRESH EXACT SEMIPRIME CONSISTENCY")
    print("=" * 78)

    p = sp.Integer(106621)
    q = sp.Integer(246473)

    N_val = p * q
    S_val = p + q
    X_val = S_val + 1
    t_val = sp.Rational(N_val, X_val)
    u_val = clean(
        t_val * (t_val + 1)
    )

    print(f"  p = {p}")
    print(f"  q = {q}")
    print(f"  N = {N_val}")
    print(f"  S = {S_val}")
    print(f"  X = {X_val}")
    print(f"  t = {t_val}")
    print(f"  u = {u_val}")

    fresh_ok = True

    for k, P in A_coeff.items():

        predicted_Q = clean(
            sum(
                A_matrix[k][r]
                * sp.binomial(8, r)
                for r in range(
                    len(A_matrix[k])
                )
            )
        )

        # Direct coefficient polynomial at j=8
        direct_Q = clean(
            A_quot[k].subs(
                j,
                sp.Integer(8),
            )
        )

        ok = predicted_Q == direct_Q

        fresh_ok = (
            fresh_ok
            and ok
        )

        print(
            f"  A quotient k={k}: {ok}"
        )

    for k, P in B_coeff.items():

        predicted_Q = clean(
            sum(
                B_matrix[k][r]
                * sp.binomial(7, r)
                for r in range(
                    len(B_matrix[k])
                )
            )
        )

        direct_Q = clean(
            B_quot[k].subs(
                j,
                sp.Integer(7),
            )
        )

        ok = predicted_Q == direct_Q

        fresh_ok = (
            fresh_ok
            and ok
        )

        print(
            f"  B quotient k={k}: {ok}"
        )

    # ------------------------------------------------------------------
    # 12. Structural interpretation
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("12. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
  Experiment 84 established the exact factorization

      [u^k] = j_under_k * Q_k(j),

  where

      j_under_k = j(j-1)...(j-k+1).

  Experiment 85 changes basis once more and writes

      Q_k(j)
        = sum_r c[k,r] * binom(j,r).

  The important object is now the coefficient matrix c[k,r].

  A particularly strong structural signature would be a simple
  triangular array such as

      c[k,r] = f(k,r)

  with

      f(k,r) = product of a few linear factors,
      or
      f(k,r) = binomial/Stirling coefficients,
      or
      f(k,r) = a(k)b(r),
      or
      f(k,r) = a(k-r)b(k,r).

  The rank test separates a fully separable construction

      c[k,r] = f(k)g(r)

  from a genuinely two-dimensional one.

  The diagonal and near-diagonal ratios test whether the coefficients
  propagate along k with a simple rule.

  The primitive row signatures test whether rational scalings are hiding
  a small integer combinatorial array.

  This is still purely a construction audit.

  No extrapolation beyond the observed index range is claimed.

  No factorization algorithm is inferred.

  The ultimate inverse question remains

      N -> layers -> t -> X -> S -> p,q.
        """
    )

    # ------------------------------------------------------------------
    # 13. Exactness audit
    # ------------------------------------------------------------------

    checks = {
        "A_triangular": A_triangular,
        "B_triangular": B_triangular,
        "A_newton_reconstruction": A_reconstruction_ok,
        "B_newton_reconstruction": B_reconstruction_ok,
        "fresh_consistency": fresh_ok,
    }

    failures = sum(
        1
        for value in checks.values()
        if not value
    )

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
        f"  ALL BASIC CHECKS PASS = "
        f"{failures == 0}"
    )

    print()
    print("=" * 78)
    print("EXPERIMENT 85 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

