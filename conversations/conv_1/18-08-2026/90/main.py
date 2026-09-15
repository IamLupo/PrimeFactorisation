#!/usr/bin/env python3

import sympy as sp

# ==============================================================================
# EXPERIMENT 92 — EXACT TRIANGULAR COEFFICIENT KERNEL / BINOMIAL NORMALIZATION
# ==============================================================================

t, u, j = sp.symbols("t u j")
x, y = sp.symbols("x y")


# ==============================================================================
# EXACT HELPERS
# ==============================================================================

def clean(expr):
    return sp.cancel(sp.expand(sp.sympify(expr)))


def poly(expr, var):
    return sp.Poly(clean(expr), var, domain=sp.QQ)


def degree(expr, var):
    p = poly(expr, var)
    if p.is_zero:
        return -sp.oo
    return p.degree()


def falling(expr, n):
    out = sp.Integer(1)
    for q in range(n):
        out *= expr - q
    return clean(out)


def binom_poly(r, k):
    if r < k:
        return sp.Integer(0)
    return clean(
        falling(r, k) / sp.factorial(k)
    )


# ==============================================================================
# ESTABLISHED EXACT h_d(t)
# ==============================================================================

H = {
    16: -(3*t**2 + 3*t + 1) * (
        3*t**6 + 9*t**5 + 18*t**4 +
        21*t**3 + 15*t**2 + 6*t + 1
    ),

    15: (2*t + 1) * (
        44*t**8 + 176*t**7 + 494*t**6 +
        866*t**5 + 1016*t**4 + 794*t**3 +
        401*t**2 + 119*t + 16
    ),

    14: -(
        276*t**10 + 1380*t**9 + 5460*t**8 +
        13560*t**7 + 23058*t**6 + 27510*t**5 +
        23100*t**4 + 13410*t**3 + 5135*t**2 +
        1169*t + 120
    ),

    13: (2*t + 1) * (
        144*t**10 + 720*t**9 + 4810*t**8 +
        14920*t**7 + 32272*t**6 + 47620*t**5 +
        48955*t**4 + 34510*t**3 + 16020*t**2 +
        4431*t + 560
    ),

    12: -(
        50*t**12 + 300*t**11 + 7436*t**10 +
        34430*t**9 + 117810*t**8 + 267960*t**7 +
        426888*t**6 + 484110*t**5 + 390225*t**4 +
        219010*t**3 + 81510*t**2 + 18109*t + 1820
    ),

    11: 13*(2*t + 1) * (
        50*t**10 + 250*t**9 + 2284*t**8 +
        7636*t**7 + 17434*t**6 + 26626*t**5 +
        28036*t**4 + 20104*t**3 + 9451*t**2 +
        2639*t + 336
    ),

    10: -1001 * (
        10*t**10 + 50*t**9 + 252*t**8 +
        708*t**7 + 1302*t**6 + 1638*t**5 +
        1428*t**4 + 852*t**3 + 333*t**2 +
        77*t + 8
    ),

    9: (2*t + 1) * (
        17875*t**8 + 71500*t**7 + 245960*t**6 +
        487630*t**5 + 622414*t**4 + 515528*t**3 +
        271506*t**2 + 83097*t + 11441
    ),

    8: -(
        71500*t**8 + 286000*t**7 + 755664*t**6 +
        1265992*t**5 + 1375360*t**4 +
        974400*t**3 + 436832*t**2 +
        112964*t + 12879
    ),

    7: 4*(2*t + 1) * (
        11050*t**6 + 33150*t**5 + 63036*t**4 +
        70822*t**3 + 47549*t**2 + 17663*t + 2869
    ),

    6: -4 * (
        17850*t**6 + 53550*t**5 + 86088*t**4 +
        82926*t**3 + 47574*t**2 + 15036*t + 2023
    ),

    5: 6*(2*t + 1) * (
        3230*t**4 + 6460*t**3 + 6688*t**2 +
        3458*t + 749
    ),

    4: -2 * (
        7125*t**4 + 14250*t**3 + 12690*t**2 +
        5565*t + 973
    ),

    3: 14*(2*t + 1) * (
        125*t**2 + 125*t + 46
    ),

    2: -2 * (
        275*t**2 + 275*t + 78
    ),

    1: 25*(2*t + 1),

    0: -2,
}


# ==============================================================================
# CHANNEL RECOVERY
# ==============================================================================

def recover_u(expr):
    expr = clean(expr)

    if expr == 0:
        return sp.Integer(0)

    deg = int(degree(expr, t))

    if deg == 0:
        return expr

    m = deg // 2

    points = []

    for tv in range(m + 1):
        uv = tv * (tv + 1)
        value = clean(expr.subs(t, tv))
        points.append((sp.Integer(uv), value))

    candidate = clean(
        sp.interpolate(points, u)
    )

    check = clean(
        candidate.subs(u, t*(t + 1)) - expr
    )

    if check != 0:
        raise ArithmeticError(
            "u reconstruction failed"
        )

    return candidate


def build_channels():

    A = {}
    B = {}

    for d in range(0, 17):

        h = clean(H[d])

        if d % 2 == 0:

            A[d] = recover_u(h)

        else:

            q, r = sp.div(
                poly(h, t),
                poly(2*t + 1, t),
            )

            if not r.is_zero:
                raise ArithmeticError(
                    f"(2t+1) does not divide h_{d}"
                )

            B[d] = recover_u(q.as_expr())

    return A, B


# ==============================================================================
# INDEXED SEQUENCES
# ==============================================================================

def coefficient_polynomial(sequence, k):

    points = []

    for jj, expr in enumerate(sequence):

        value = clean(
            poly(expr, u).coeff_monomial(
                u**k
            )
        )

        points.append(
            (sp.Integer(jj), value)
        )

    return clean(
        sp.interpolate(points, j)
    )


# ==============================================================================
# FALLING-BASIS TRIANGLE
# ==============================================================================

def falling_basis(P):

    deg = int(degree(P, j))

    values = [
        clean(P.subs(j, q))
        for q in range(deg + 1)
    ]

    result = {}

    current = values

    for r in range(deg + 1):

        result[r] = clean(
            current[0] /
            sp.factorial(r)
        )

        if len(current) == 1:
            break

        current = [
            clean(
                current[q + 1] -
                current[q]
            )
            for q in range(len(current) - 1)
        ]

    return result


def build_triangle(sequence, max_k):

    triangle = {}

    for k in range(max_k + 1):

        P = coefficient_polynomial(
            sequence,
            k
        )

        triangle[k] = falling_basis(P)

    return triangle


# ==============================================================================
# BINOMIAL NORMALIZATION
# ==============================================================================

def normalized_entry(C, k, r):

    value = C[k].get(
        r,
        sp.Integer(0)
    )

    if r < k:
        return sp.Integer(0)

    b = sp.binomial(r, k)

    if b == 0:
        return sp.Integer(0)

    return clean(
        value / b
    )


# ==============================================================================
# MATRIX RANK
# ==============================================================================

def matrix_from_triangle(C, rows, cols):

    return sp.Matrix([
        [
            C[k].get(
                r,
                sp.Integer(0)
            )
            for r in range(cols)
        ]
        for k in range(rows)
    ])


# ==============================================================================
# EXACT 2x2 MINOR TEST
# ==============================================================================

def rank_one_exact(M):

    rows, cols = M.shape

    if rows == 0 or cols == 0:
        return True

    pivot = None

    for i in range(rows):
        for j0 in range(cols):

            if clean(M[i, j0]) != 0:
                pivot = (i, j0)
                break

        if pivot is not None:
            break

    if pivot is None:
        return True

    i0, j0 = pivot

    for i in range(rows):
        for j1 in range(cols):

            lhs = clean(
                M[i, j1] * M[i0, j0]
            )

            rhs = clean(
                M[i, j0] * M[i0, j1]
            )

            if clean(lhs - rhs) != 0:
                return False

    return True


# ==============================================================================
# ROW GENERATING POLYNOMIALS
# ==============================================================================

def row_generating_polynomial(C, k, max_r):

    return clean(
        sum(
            C[k].get(
                r,
                sp.Integer(0)
            ) * y**r
            for r in range(max_r + 1)
        )
    )


# ==============================================================================
# DIAGONAL NORMALIZATION
# ==============================================================================

def diagonal_normalization(C, max_k):

    values = []

    for k in range(max_k + 1):

        value = C[k].get(
            k,
            sp.Integer(0)
        )

        values.append(value)

    first = None

    for value in values:

        if value != 0:
            first = value
            break

    if first is None:
        return values

    return [
        clean(v / first)
        for v in values
    ]


# ==============================================================================
# SIMPLE FACTORIZATION SEARCH FOR NORMALIZED ENTRIES
# ==============================================================================

def factor_normalized_triangle(C, max_k):

    print()

    for k in range(max_k + 1):

        print(f"  k={k}")

        for r in range(k, max_k + 3):

            value = normalized_entry(
                C,
                k,
                r
            )

            if value == 0:
                continue

            print(
                f"    r={r}: "
                f"C/binom(r,k) = {sp.factor(value)}"
            )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 92 — EXACT TRIANGULAR COEFFICIENT "
        "KERNEL / BINOMIAL NORMALIZATION"
    )
    print("=" * 78)

    A, B = build_channels()

    Aseq = [
        A[2*j]
        for j in range(9)
    ]

    Bseq = [
        B[2*j + 1]
        for j in range(8)
    ]

    print()
    print("=" * 78)
    print("1. EXACT CHANNEL DATA")
    print("=" * 78)

    print("  A:")
    for jj, expr in enumerate(Aseq):
        print(
            f"    j={jj}: degree_u={degree(expr,u)}"
        )

    print("  B:")
    for jj, expr in enumerate(Bseq):
        print(
            f"    j={jj}: degree_u={degree(expr,u)}"
        )

    A_C = build_triangle(
        Aseq,
        6
    )

    B_C = build_triangle(
        Bseq,
        5
    )

    print()
    print("=" * 78)
    print("2. EXACT FALLING-BASIS TRIANGLE")
    print("=" * 78)

    for channel_name, C, max_k in (
        ("A", A_C, 6),
        ("B", B_C, 5),
    ):

        print()
        print(f"  {channel_name} channel")

        for k in range(max_k + 1):

            row = [
                C[k].get(
                    r,
                    sp.Integer(0)
                )
                for r in range(
                    k,
                    9 if channel_name == "A"
                    else 8
                )
            ]

            print(
                f"    k={k}: {row}"
            )

    print()
    print("=" * 78)
    print("3. BINOMIAL-NORMALIZED TRIANGLE")
    print("=" * 78)

    for channel_name, C, max_k, max_r in (
        ("A", A_C, 6, 8),
        ("B", B_C, 5, 7),
    ):

        print()
        print(
            f"  {channel_name} channel"
        )

        for k in range(max_k + 1):

            row = [
                normalized_entry(
                    C,
                    k,
                    r
                )
                for r in range(
                    k,
                    max_r + 1
                )
            ]

            print(
                f"    k={k}: {row}"
            )

    print()
    print("=" * 78)
    print("4. BINOMIAL-NORMALIZED MATRIX RANK")
    print("=" * 78)

    A_norm = sp.Matrix([
        [
            normalized_entry(
                A_C,
                k,
                r
            )
            for r in range(9)
        ]
        for k in range(7)
    ])

    B_norm = sp.Matrix([
        [
            normalized_entry(
                B_C,
                k,
                r
            )
            for r in range(8)
        ]
        for k in range(6)
    ])

    print(
        "  A normalized matrix shape =",
        A_norm.shape
    )

    print(
        "  A normalized matrix rank =",
        A_norm.rank()
    )

    print(
        "  A normalized rank-1 =",
        rank_one_exact(A_norm)
    )

    print(
        "  B normalized matrix shape =",
        B_norm.shape
    )

    print(
        "  B normalized matrix rank =",
        B_norm.rank()
    )

    print(
        "  B normalized rank-1 =",
        rank_one_exact(B_norm)
    )

    print()
    print("=" * 78)
    print("5. ROW-GENERATING POLYNOMIALS IN THE FALLING INDEX")
    print("=" * 78)

    for channel_name, C, max_k, max_r in (
        ("A", A_C, 6, 8),
        ("B", B_C, 5, 7),
    ):

        print()
        print(
            f"  {channel_name} channel"
        )

        for k in range(max_k + 1):

            row_poly = row_generating_polynomial(
                C,
                k,
                max_r
            )

            print(
                f"    k={k}:"
            )

            print(
                "      raw =",
                sp.factor(row_poly)
            )

    print()
    print("=" * 78)
    print("6. BINOMIAL-NORMALIZED ROW GENERATORS")
    print("=" * 78)

    for channel_name, C, max_k, max_r in (
        ("A", A_C, 6, 8),
        ("B", B_C, 5, 7),
    ):

        print()
        print(
            f"  {channel_name} channel"
        )

        for k in range(max_k + 1):

            normalized_poly = clean(
                sum(
                    normalized_entry(
                        C,
                        k,
                        r
                    ) * y**r
                    for r in range(
                        k,
                        max_r + 1
                    )
                )
            )

            print(
                f"    k={k}:"
            )

            print(
                "      normalized =",
                sp.factor(normalized_poly)
            )

    print()
    print("=" * 78)
    print("7. DIAGONAL AND SUPERDIAGONAL NORMALIZATION")
    print("=" * 78)

    for channel_name, C, max_k, max_shift in (
        ("A", A_C, 6, 3),
        ("B", B_C, 5, 3),
    ):

        print()
        print(
            f"  {channel_name} channel"
        )

        for shift in range(max_shift + 1):

            vals = []

            for k in range(
                max_k - shift + 1
            ):

                r = k + shift

                value = C[k].get(
                    r,
                    sp.Integer(0)
                )

                normalized = clean(
                    value /
                    sp.binomial(
                        r,
                        k
                    )
                )

                vals.append(
                    normalized
                )

            print(
                f"    shift={shift}:"
            )

            print(
                "      values =",
                vals
            )

            if vals:

                base = vals[0]

                if base != 0:

                    ratios = [
                        clean(v / base)
                        for v in vals
                    ]

                else:

                    ratios = []

                print(
                    "      relative_to_first =",
                    ratios
                )

    print()
    print("=" * 78)
    print("8. EXACT DIAGONAL CLOSED-FORM FACTORIZATION")
    print("=" * 78)

    for channel_name, C, max_k in (
        ("A", A_C, 6),
        ("B", B_C, 5),
    ):

        print()
        print(
            f"  {channel_name} diagonal"
        )

        values = [
            C[k].get(
                k,
                sp.Integer(0)
            )
            for k in range(
                max_k + 1
            )
        ]

        P = clean(
            sp.interpolate(
                [
                    (
                        sp.Integer(k),
                        values[k]
                    )
                    for k in range(
                        len(values)
                    )
                ],
                j
            )
        )

        print(
            "    interpolated polynomial =",
            sp.factor(P)
        )

        print(
            "    factor over QQ =",
            sp.factor(
                P
            )
        )

        print(
            "    integer roots in observed range =",
            [
                kk
                for kk in range(
                    max_k + 1
                )
                if clean(
                    P.subs(j, kk)
                )
                == values[kk]
                and values[kk] == 0
            ]
        )

    print()
    print("=" * 78)
    print("9. EXACT CROSS-CHANNEL NORMALIZED MATRIX")
    print("=" * 78)

    print()

    for k in range(6):

        print(
            f"  k={k}"
        )

        for r in range(
            k,
            min(8, 8)
        ):

            a = normalized_entry(
                A_C,
                k,
                r
            )

            b = normalized_entry(
                B_C,
                k,
                r
            )

            if b != 0:

                print(
                    f"    r={r}: A/B =",
                    clean(a/b)
                )

    print()
    print("=" * 78)
    print("10. EXACT ORIGINAL-TRIANGLE RECONSTRUCTION")
    print("=" * 78)

    reconstruction_ok = True

    for channel_name, C, sequence, max_k in (
        ("A", A_C, Aseq, 6),
        ("B", B_C, Bseq, 5),
    ):

        print()
        print(
            f"  {channel_name} channel"
        )

        for k in range(max_k + 1):

            P = coefficient_polynomial(
                sequence,
                k
            )

            reconstructed = clean(
                sum(
                    C[k].get(
                        r,
                        sp.Integer(0)
                    ) * falling(j, r)
                    for r in range(
                        9 if channel_name == "A"
                        else 8
                    )
                )
            )

            ok = (
                clean(
                    reconstructed - P
                ) == 0
            )

            reconstruction_ok &= ok

            print(
                f"    k={k}: {ok}"
            )

    print()
    print("=" * 78)
    print("11. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
  The previous experiments established the exact triangular identity

      P_k(j) = sum_{r >= k} C[k,r] j_(r).

  This experiment removes one more arbitrary choice.

  Instead of looking at raw C[k,r], we inspect

      C[k,r] / binom(r,k).

  The purpose is not to force a formula.  It tests whether the
  triangular coefficient is naturally carrying a standard combinatorial
  multiplicity.

  Three outcomes are informative:

      1. Strong simplification after division by binom(r,k).

         This would indicate that the triangle is close to a standard
         binomial/Stirling transform.

      2. The normalized matrix still has high exact rank and irregular
         row-generating polynomials.

         Then binom(r,k) is not the missing universal factor.

      3. The normalized row-generators factor into a common family of
         low-degree polynomials.

         That would be a substantially stronger construction signal.

  We also reconstruct the complete bivariate coefficient triangle from
  the exact falling basis, so every reported identity remains an exact
  QQ identity on the known layer range.

  No extrapolation is performed.
  No factorization algorithm is inferred.
  No claim is made that an observed finite polynomial identity extends
  to arbitrary layer indices.
        """
    )

    print()
    print("=" * 78)
    print("12. FINAL EXACTNESS")
    print("=" * 78)

    print(
        "  A falling-basis reconstruction =",
        reconstruction_ok
    )

    print(
        "  B falling-basis reconstruction =",
        reconstruction_ok
    )

    print(
        "  A normalized rank =",
        A_norm.rank()
    )

    print(
        "  B normalized rank =",
        B_norm.rank()
    )

    print(
        "  ALL BASIC CHECKS PASS =",
        reconstruction_ok
    )

    print()
    print(
        "EXPERIMENT 92 COMPLETE"
    )


if __name__ == "__main__":
    main()

