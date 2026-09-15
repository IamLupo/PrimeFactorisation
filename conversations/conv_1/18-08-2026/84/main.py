#!/usr/bin/env python3

import sympy as sp

# ==============================================================================
# EXPERIMENT 86 — EXACT QUOTIENT-MATRIX / STIRLING-BASIS STRUCTURE AUDIT
# ==============================================================================
#
# Purpose
# -------
# Experiment 84 found the exact factor
#
#     [u^k] A(j,u) = j_under_k * Q_A,k(j)
#     [u^k] B(j,u) = j_under_k * Q_B,k(j)
#
# Experiment 85 then tried a Newton/binomial analysis, but two diagnostics
# were logically mismatched:
#
#   1. the quotient Q_k(j) need not have its first k Newton coefficients zero;
#   2. reconstruction accidentally used a truncated "reduced" row instead of
#      the complete Newton row.
#
# This experiment fixes that and asks a cleaner structural question:
#
#     What is the exact coefficient structure of Q_k(j)?
#
# We examine:
#
#   A. full Newton/binomial coefficient rows;
#   B. monomial coefficient rows;
#   C. conversion between the two bases;
#   D. exact triangularity of the ORIGINAL coefficient laws;
#   E. exact reconstruction of every quotient from its FULL Newton row;
#   F. matrix ranks of the quotient coefficient arrays;
#   G. whether rows become simpler in the falling-factorial basis;
#   H. exact diagonal / near-diagonal patterns;
#   I. cross-channel comparison.
#
# No extrapolation beyond j=0,...,8 and j=0,...,7 is claimed.
# No floating point.
# ==============================================================================


u = sp.Symbol("u")
j = sp.Symbol("j")

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


def falling(k):
    out = sp.Integer(1)
    for r in range(k):
        out *= (j - r)
    return clean(out)


def quotient_by(expr, divisor):
    P = poly(expr, j)
    D = poly(divisor, j)

    Q, R = sp.div(
        P,
        D,
        domain=QQ,
    )

    return (
        clean(Q.as_expr()),
        clean(R.as_expr()),
    )


def interpolate(values):
    points = [
        (sp.Integer(i), clean(v))
        for i, v in enumerate(values)
    ]
    return clean(sp.interpolate(points, j))


def finite_difference_table(values):
    rows = [
        [clean(v) for v in values]
    ]

    while len(rows[-1]) > 1:
        prev = rows[-1]

        nxt = [
            clean(prev[i + 1] - prev[i])
            for i in range(len(prev) - 1)
        ]

        rows.append(nxt)

    return rows


def newton_coefficients_from_values(values):
    table = finite_difference_table(values)

    return [
        clean(row[0])
        for row in table
    ]


def newton_coefficients(expr, max_j):
    values = [
        clean(expr.subs(j, sp.Integer(i)))
        for i in range(max_j + 1)
    ]

    return newton_coefficients_from_values(values)


def reconstruct_from_newton(coeffs):
    out = sp.Integer(0)

    for r, c in enumerate(coeffs):
        out += c * sp.binomial(j, r)

    return clean(out)


def monomial_coefficients(expr):
    P = poly(expr, j)

    return [
        clean(P.nth(r))
        for r in range(P.degree() + 1)
    ]


def primitive_integer_signature(values):
    vals = [
        sp.Rational(v)
        for v in values
    ]

    if not vals:
        return []

    lcm = sp.Integer(1)

    for v in vals:
        lcm = sp.ilcm(
            lcm,
            int(v.q),
        )

    ints = [
        int(v * lcm)
        for v in vals
    ]

    g = 0

    for x in ints:
        g = sp.igcd(g, abs(x))

    if g == 0:
        return [0 for _ in ints]

    return [x // g for x in ints]


def rational_gcd(values):
    vals = [
        sp.Rational(v)
        for v in values
        if sp.Rational(v) != 0
    ]

    if not vals:
        return sp.Integer(0)

    lcm = sp.Integer(1)

    for v in vals:
        lcm = sp.ilcm(
            lcm,
            int(v.q),
        )

    ints = [
        int(v * lcm)
        for v in vals
    ]

    g = 0

    for x in ints:
        g = sp.igcd(g, abs(x))

    if g == 0:
        return sp.Integer(0)

    return sp.Rational(g, lcm)


def matrix_rank(rows, width=None):
    if not rows:
        return 0

    if width is None:
        width = max(len(row) for row in rows)

    padded = []

    for row in rows:
        padded.append(
            [
                sp.Rational(x)
                for x in row
            ]
            + [
                sp.Rational(0)
            ] * (width - len(row))
        )

    return sp.Matrix(padded).rank()


def factor(expr):
    return sp.factor(clean(expr))


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
# Build coefficient laws and quotient laws
# ------------------------------------------------------------------------------

def build_data(channel):
    max_k = max(
        int(degree(expr, u))
        for expr in channel.values()
    )

    coefficient_polys = {}
    quotient_polys = {}

    for k in range(max_k + 1):

        values = [
            coeff_u(
                channel[idx],
                k,
            )
            for idx in range(len(channel))
        ]

        P = interpolate(values)

        coefficient_polys[k] = clean(P)

        divisor = falling(k)

        Q, R = quotient_by(
            P,
            divisor,
        )

        if clean(R) != 0:
            raise RuntimeError(
                f"falling-factorial quotient failed for k={k}"
            )

        quotient_polys[k] = clean(Q)

    return coefficient_polys, quotient_polys


# ------------------------------------------------------------------------------
# Exact original triangularity
#
# This is where the falling factorial belongs:
#
#   coefficient P_k(j) is divisible by j_under_k.
#
# We test this directly, rather than expecting the quotient Q_k(j)
# to remain triangular.
# ------------------------------------------------------------------------------

def original_falling_triangularity(
    coefficient_polys
):
    result = True

    for k, P in coefficient_polys.items():

        divisor = falling(k)

        _, remainder = quotient_by(
            P,
            divisor,
        )

        ok = clean(remainder) == 0

        result = result and ok

        print(
            f"    k={k}: "
            f"j_under_{k} divides P_k(j) = {ok}"
        )

    return result


# ------------------------------------------------------------------------------
# Analyze a channel
# ------------------------------------------------------------------------------

def analyze_channel(
    name,
    channel,
    max_j,
):
    print()
    print("-" * 78)
    print(f"{name} CHANNEL")
    print("-" * 78)

    coefficient_polys, quotient_polys = build_data(
        channel
    )

    print()
    print("  Original falling-factorial divisibility")

    falling_ok = original_falling_triangularity(
        coefficient_polys
    )

    print(
        f"  ALL ORIGINAL FALLING DIVISIONS = "
        f"{falling_ok}"
    )

    # ------------------------------------------------------------------
    # Full quotient Newton rows
    # ------------------------------------------------------------------

    print()
    print("  FULL QUOTIENT NEWTON ROWS")

    quotient_newton = {}

    for k, Q in quotient_polys.items():

        coeffs = newton_coefficients(
            Q,
            max_j,
        )

        quotient_newton[k] = coeffs

        print(
            f"    [u^{k}] Q_k(j):"
        )
        print(
            f"      {coeffs}"
        )

    # ------------------------------------------------------------------
    # Full quotient monomial rows
    # ------------------------------------------------------------------

    print()
    print("  FULL QUOTIENT MONOMIAL ROWS")

    quotient_monomial = {}

    for k, Q in quotient_polys.items():

        coeffs = monomial_coefficients(Q)

        quotient_monomial[k] = coeffs

        print(
            f"    [u^{k}] monomial coefficients:"
        )
        print(
            f"      {coeffs}"
        )

    # ------------------------------------------------------------------
    # Reconstruction audit
    # ------------------------------------------------------------------

    print()
    print("  EXACT QUOTIENT RECONSTRUCTION")

    reconstruction_ok = True

    for k, Q in quotient_polys.items():

        rebuilt = reconstruct_from_newton(
            quotient_newton[k]
        )

        ok = (
            clean(rebuilt)
            ==
            clean(Q)
        )

        reconstruction_ok = (
            reconstruction_ok
            and ok
        )

        print(
            f"    Q_{k}: {ok}"
        )

    # ------------------------------------------------------------------
    # Degree collapse
    # ------------------------------------------------------------------

    print()
    print("  DEGREE PROFILE")

    for k, Q in quotient_polys.items():

        print(
            f"    k={k}: "
            f"degree(P_k)={degree(coefficient_polys[k], j)} "
            f"degree(Q_k)={degree(Q, j)}"
        )

    # ------------------------------------------------------------------
    # Quotient matrix rank
    # ------------------------------------------------------------------

    max_width_newton = max(
        len(row)
        for row in quotient_newton.values()
    )

    newton_rows = [
        quotient_newton[k]
        for k in sorted(quotient_newton)
    ]

    rank_newton = matrix_rank(
        newton_rows,
        max_width_newton,
    )

    print()
    print("  QUOTIENT NEWTON MATRIX")

    print(
        f"    shape = "
        f"({len(newton_rows)}, {max_width_newton})"
    )

    print(
        f"    rank = {rank_newton}"
    )

    # ------------------------------------------------------------------
    # Quotient monomial matrix rank
    # ------------------------------------------------------------------

    max_width_monomial = max(
        len(row)
        for row in quotient_monomial.values()
    )

    monomial_rows = [
        quotient_monomial[k]
        for k in sorted(quotient_monomial)
    ]

    rank_monomial = matrix_rank(
        monomial_rows,
        max_width_monomial,
    )

    print()
    print("  QUOTIENT MONOMIAL MATRIX")

    print(
        f"    shape = "
        f"({len(monomial_rows)}, {max_width_monomial})"
    )

    print(
        f"    rank = {rank_monomial}"
    )

    # ------------------------------------------------------------------
    # Row primitive signatures
    # ------------------------------------------------------------------

    print()
    print("  PRIMITIVE INTEGER SIGNATURES OF FULL NEWTON ROWS")

    for k in sorted(quotient_newton):

        signature = primitive_integer_signature(
            quotient_newton[k]
        )

        print(
            f"    k={k}: {signature}"
        )

    # ------------------------------------------------------------------
    # Exact column gcds
    # ------------------------------------------------------------------

    print()
    print("  NEWTON COLUMN GCDS")

    column_gcds = []

    for col in range(max_width_newton):

        values = []

        for row in newton_rows:

            if col < len(row):
                values.append(row[col])

        g = rational_gcd(values)

        column_gcds.append(g)

        print(
            f"    column={col}: gcd={g}"
        )

    # ------------------------------------------------------------------
    # Simple diagonal comparisons
    # ------------------------------------------------------------------

    print()
    print("  DIAGONAL / NEAR-DIAGONAL AUDIT")

    for k in sorted(quotient_newton):

        row = quotient_newton[k]

        if not row:
            continue

        base = row[0]

        ratios = []

        if base != 0:

            for r in range(
                min(4, len(row))
            ):
                if row[r] != 0:
                    ratios.append(
                        (
                            r,
                            clean(
                                row[r] / base
                            ),
                        )
                    )

        print(
            f"    k={k}: {ratios}"
        )

    # ------------------------------------------------------------------
    # Factor quotient polynomials
    # ------------------------------------------------------------------

    print()
    print("  QUOTIENT FACTORIZATIONS")

    for k, Q in quotient_polys.items():

        print(
            f"    Q_{k}(j) = {factor(Q)}"
        )

    return {
        "coefficient_polys": coefficient_polys,
        "quotient_polys": quotient_polys,
        "quotient_newton": quotient_newton,
        "quotient_monomial": quotient_monomial,
        "falling_ok": falling_ok,
        "reconstruction_ok": reconstruction_ok,
        "rank_newton": rank_newton,
        "rank_monomial": rank_monomial,
    }


# ------------------------------------------------------------------------------
# Fresh exact consistency
# ------------------------------------------------------------------------------

def fresh_consistency(
    A_data,
    B_data,
):
    # We check the reconstructed coefficient laws directly at
    # the already-known endpoint indices.
    #
    # This remains purely symbolic and does not require p,q.

    A_ok = True
    B_ok = True

    for k, P in A_data["coefficient_polys"].items():

        Q = A_data["quotient_polys"][k]
        divisor = falling(k)

        rebuilt = clean(
            divisor * Q
        )

        ok = (
            clean(rebuilt)
            ==
            clean(P)
        )

        A_ok = A_ok and ok

    for k, P in B_data["coefficient_polys"].items():

        Q = B_data["quotient_polys"][k]
        divisor = falling(k)

        rebuilt = clean(
            divisor * Q
        )

        ok = (
            clean(rebuilt)
            ==
            clean(P)
        )

        B_ok = B_ok and ok

    return A_ok, B_ok


# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

def main():

    print("=" * 78)
    print("EXPERIMENT 86 — EXACT QUOTIENT-MATRIX / STIRLING-BASIS STRUCTURE AUDIT")
    print("=" * 78)

    print()
    print("=" * 78)
    print("0. EXACT SYMBOLIC SETUP")
    print("=" * 78)
    print("  A_j = A_{2j}(u), j=0,...,8")
    print("  B_j = B_{2j+1}(u), j=0,...,7")
    print("  arithmetic = exact QQ")
    print("  floating point = forbidden")

    A_data = analyze_channel(
        "A",
        A,
        8,
    )

    B_data = analyze_channel(
        "B",
        B,
        7,
    )

    print()
    print("=" * 78)
    print("3. CROSS-CHANNEL QUOTIENT NEWTON MATRIX")
    print("=" * 78)

    A_newton = A_data["quotient_newton"]
    B_newton = B_data["quotient_newton"]

    shared_k = sorted(
        set(A_newton.keys())
        &
        set(B_newton.keys())
    )

    for k in shared_k:

        arow = A_newton[k]
        brow = B_newton[k]

        width = min(
            len(arow),
            len(brow),
        )

        ratios = []

        for r in range(width):

            av = clean(arow[r])
            bv = clean(brow[r])

            if av != 0 and bv != 0:
                ratios.append(
                    (
                        r,
                        clean(av / bv)
                    )
                )

        print(
            f"  k={k}: entrywise ratios={ratios}"
        )

    print()
    print("=" * 78)
    print("4. EXPLICIT BASIS-CORRECTION CHECK")
    print("=" * 78)

    print(
        """
  The previous run incorrectly treated the Newton coefficients
  of Q_k(j) as though they should begin with k zeros.

  That property belongs to the ORIGINAL coefficient polynomial

      P_k(j) = j_under_k Q_k(j),

  not to Q_k(j) itself.

  We therefore verify both statements separately:

      ORIGINAL:
          j_under_k | P_k(j)

      QUOTIENT:
          Q_k(j) is represented by its COMPLETE Newton row.

  No truncation is performed.
        """
    )

    print()
    print("  A original falling divisibility:")
    A_fall = original_falling_triangularity(
        A_data["coefficient_polys"]
    )

    print()
    print("  B original falling divisibility:")
    B_fall = original_falling_triangularity(
        B_data["coefficient_polys"]
    )

    print()
    print("=" * 78)
    print("5. FRESH SYMBOLIC CONSISTENCY")
    print("=" * 78)

    fresh_A, fresh_B = fresh_consistency(
        A_data,
        B_data,
    )

    print(
        f"  A quotient reconstruction = {fresh_A}"
    )

    print(
        f"  B quotient reconstruction = {fresh_B}"
    )

    print()
    print("=" * 78)
    print("6. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
  The corrected analysis separates three different levels.

  LEVEL 1
      [u^k] channel coefficient

          P_k(j)

  LEVEL 2
      exact combinatorial factor

          j_under_k

  LEVEL 3
      residual polynomial

          Q_k(j) = P_k(j) / j_under_k.

  The crucial observation from Experiment 84 was that LEVEL 2
  exists exactly for every observed coefficient.

  Experiment 86 now asks whether LEVEL 3 has a simpler basis.

  In particular, compare:

      monomial basis:
          Q_k(j) = sum a_r j^r

      Newton/binomial basis:
          Q_k(j) = sum c_r binom(j,r).

  If the Newton rows have substantially lower rank or a simple
  factor pattern, that is evidence for a combinatorial construction.

  If both bases remain high-rank and unfactored, then the falling
  factorial is still the clearest structural component currently
  exposed.

  This remains a finite exact audit only.

  No extrapolation beyond the available layer indices is claimed.
  No factorization algorithm is inferred.

  The long-term inverse question remains:

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
        "A_original_falling_factorials": A_fall,
        "B_original_falling_factorials": B_fall,
        "A_newton_reconstruction": A_data["reconstruction_ok"],
        "B_newton_reconstruction": B_data["reconstruction_ok"],
        "fresh_A": fresh_A,
        "fresh_B": fresh_B,
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
        f"  ALL EXACT CHECKS PASS = "
        f"{failures == 0}"
    )

    print()
    print("=" * 78)
    print("EXPERIMENT 86 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

