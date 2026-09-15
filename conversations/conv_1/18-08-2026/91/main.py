#!/usr/bin/env python3

import sympy as sp


# ==============================================================================
# EXPERIMENT 93 — EXACT OFFSET-TIRED TRIANGLE / DIAGONAL-SEPARABILITY AUDIT
# ==============================================================================

u, j, x = sp.symbols("u j x")


# ==============================================================================
# EXACT HELPERS
# ==============================================================================

def S(expr):
    return sp.sympify(expr)


def clean(expr):
    return sp.cancel(sp.expand(S(expr)))


def poly(expr, var):
    return sp.Poly(
        clean(expr),
        var,
        domain=sp.QQ,
    )


def degree(expr, var):
    p = poly(expr, var)
    if p.is_zero:
        return -sp.oo
    return p.degree()


def falling(n, r):
    out = sp.Integer(1)
    for q in range(r):
        out *= n - q
    return clean(out)


def rank_one(M):
    rows, cols = M.shape

    if rows == 0 or cols == 0:
        return True

    pivot = None

    for i in range(rows):
        for c in range(cols):
            if clean(M[i, c]) != 0:
                pivot = (i, c)
                break
        if pivot is not None:
            break

    if pivot is None:
        return True

    i0, c0 = pivot

    for i in range(rows):
        for c in range(cols):
            lhs = clean(M[i, c] * M[i0, c0])
            rhs = clean(M[i, c0] * M[i0, c])
            if clean(lhs - rhs) != 0:
                return False

    return True


# ==============================================================================
# EXACT ORIGINAL CHANNEL POLYNOMIALS
# ==============================================================================

A = {
    0: -2,

    1: -2*(275*u + 78),

    2: -2*(7125*u**2 + 5565*u + 973),

    3: -4*(17850*u**3 + 32538*u**2 + 15036*u + 2023),

    4: -2*(7125*u**2 + 5565*u + 973) * 0 + (
        -71500*u**4
        -326664*u**3
        -323868*u**2
        -112964*u
        -12879
    ),

    5: -1001*(10*u**5 + 152*u**4 + 340*u**3
              +256*u**2 + 77*u + 8),

    6: -(
        50*u**6
        +6686*u**5
        +50200*u**4
        +92208*u**3
        +63401*u**2
        +18109*u
        +1820
    ),

    7: -(
        276*u**5
        +2700*u**4
        +5478*u**3
        +3966*u**2
        +1169*u
        +120
    ),

    8: -(
        9*u**4
        +30*u**3
        +27*u**2
        +9*u
        +1
    ),
}

B = {
    0: 25,

    1: 14*(125*u + 46),

    2: 6*(3230*u**2 + 3458*u + 749),

    3: 4*(11050*u**3 + 29886*u**2
          +17663*u + 2869),

    4: (
        17875*u**4
        +138710*u**3
        +188409*u**2
        +83097*u
        +11441
    ),

    5: 13*(50*u**5 + 1784*u**4 + 6480*u**3
           +6812*u**2 +2639*u +336),

    6: (
        144*u**5
        +3370*u**4
        +11332*u**3
        +11589*u**2
        +4431*u
        +560
    ),

    7: (
        44*u**4
        +230*u**3
        +282*u**2
        +119*u
        +16
    ),
}


# ==============================================================================
# RECOVER ORIGINAL COEFFICIENT POLYNOMIALS P_k(j)
# ==============================================================================

def coefficient_values(sequence, k):
    """
    P_k(j) is reconstructed from j=0,...,m exactly.
    """
    points = []

    for jj, expr in enumerate(sequence):
        value = clean(
            poly(expr, u).coeff_monomial(u**k)
        )
        points.append(
            (sp.Integer(jj), value)
        )

    return clean(
        sp.interpolate(points, j)
    )


def falling_coefficients(P):
    """
    P(j) = sum_r C_r * j_(r)

    Since
        Delta^r P(0) = C_r r!,
    we have
        C_r = Delta^r P(0) / r!.
    """
    deg = int(degree(P, j))

    values = [
        clean(P.subs(j, q))
        for q in range(deg + 1)
    ]

    coeffs = {}

    current = values

    for r in range(deg + 1):
        coeffs[r] = clean(
            current[0] / sp.factorial(r)
        )

        if len(current) == 1:
            break

        current = [
            clean(current[q + 1] - current[q])
            for q in range(len(current) - 1)
        ]

    return coeffs


def build_triangle(sequence, max_k):
    C = {}

    for k in range(max_k + 1):
        P = coefficient_values(sequence, k)
        C[k] = falling_coefficients(P)

    return C


# ==============================================================================
# OFFSET MATRICES
# ==============================================================================

def offset_matrix(C, max_k, max_shift):
    """
    Matrix rows = k, columns = s=r-k.

    Entry = C[k,k+s].
    """
    return sp.Matrix([
        [
            C[k].get(
                k + s,
                sp.Integer(0)
            )
            for s in range(max_shift + 1)
        ]
        for k in range(max_k + 1)
    ])


def diagonal_normalized(C, max_k, max_shift):
    """
    E[k,s] = C[k,k+s] / C[k,k].
    """
    return sp.Matrix([
        [
            clean(
                C[k].get(
                    k + s,
                    sp.Integer(0)
                )
                /
                C[k].get(k, sp.Integer(1))
            )
            if C[k].get(k, sp.Integer(0)) != 0
            else sp.Integer(0)
            for s in range(max_shift + 1)
        ]
        for k in range(max_k + 1)
    ])


def shift_normalized(C, max_k, max_shift):
    """
    Remove the diagonal AND the natural binomial multiplicity:

      F[k,s] =
        C[k,k+s] / ( C[k,k] * binom(k+s,k) ).

    If this becomes primarily a function of s, that is a strong
    signal of an offset-driven construction law.
    """
    rows = []

    for k in range(max_k + 1):

        row = []

        d = C[k].get(k, sp.Integer(0))

        for s in range(max_shift + 1):

            r = k + s
            c = C[k].get(r, sp.Integer(0))

            if d == 0:
                row.append(sp.Integer(0))
                continue

            b = sp.binomial(r, k)

            if b == 0:
                row.append(sp.Integer(0))
            else:
                row.append(
                    clean(
                        c / (d * b)
                    )
                )

        rows.append(row)

    return sp.Matrix(rows)


# ==============================================================================
# ROW-EQUALITY / TOEPLITZ-LIKE TEST
# ==============================================================================

def column_polynomial(M):
    """
    For each offset s, interpolate the dependence on k.
    """
    out = {}

    rows, cols = M.shape

    for s in range(cols):

        points = [
            (
                sp.Integer(k),
                clean(M[k, s])
            )
            for k in range(rows)
        ]

        out[s] = clean(
            sp.interpolate(points, j)
        )

    return out


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 93 — EXACT OFFSET-TIERED TRIANGLE / "
        "DIAGONAL-SEPARABILITY AUDIT"
    )
    print("=" * 78)

    Aseq = [A[j] for j in range(9)]
    Bseq = [B[j] for j in range(8)]

    A_C = build_triangle(Aseq, 6)
    B_C = build_triangle(Bseq, 5)

    print()
    print("=" * 78)
    print("1. EXACT TRIANGULAR COEFFICIENT DATA")
    print("=" * 78)

    for name, C, max_k, max_r in (
        ("A", A_C, 6, 8),
        ("B", B_C, 5, 7),
    ):
        print()
        print(f"  {name} channel")
        for k in range(max_k + 1):
            row = [
                C[k].get(r, sp.Integer(0))
                for r in range(k, max_r + 1)
            ]
            print(f"    k={k}: {row}")

    # --------------------------------------------------------------------------
    # RAW OFFSET MATRIX
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. RAW OFFSET MATRIX")
    print("=" * 78)

    A_raw = offset_matrix(A_C, 6, 3)
    B_raw = offset_matrix(B_C, 5, 3)

    print()
    print("  A raw C[k,k+s] =")
    print(A_raw)

    print()
    print("  B raw C[k,k+s] =")
    print(B_raw)

    # --------------------------------------------------------------------------
    # DIAGONAL NORMALIZATION
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. DIAGONAL-NORMALIZED OFFSET MATRIX")
    print("=" * 78)

    A_diag = diagonal_normalized(A_C, 6, 3)
    B_diag = diagonal_normalized(B_C, 5, 3)

    print()
    print("  A:")
    for k in range(A_diag.rows):
        print(
            f"    k={k}:",
            [
                clean(A_diag[k, s])
                for s in range(A_diag.cols)
            ]
        )

    print()
    print("  B:")
    for k in range(B_diag.rows):
        print(
            f"    k={k}:",
            [
                clean(B_diag[k, s])
                for s in range(B_diag.cols)
            ]
        )

    print()
    print("  A diagonal-normalized rank =",
          A_diag.rank())

    print("  B diagonal-normalized rank =",
          B_diag.rank())

    print("  A diagonal-normalized rank-1 =",
          rank_one(A_diag))

    print("  B diagonal-normalized rank-1 =",
          rank_one(B_diag))

    # --------------------------------------------------------------------------
    # BINOMIAL-REMOVED OFFSET MATRIX
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. BINOMIAL-REMOVED OFFSET MATRIX")
    print("=" * 78)

    A_shift = shift_normalized(A_C, 6, 3)
    B_shift = shift_normalized(B_C, 5, 3)

    print()
    print("  A:")
    for k in range(A_shift.rows):
        print(
            f"    k={k}:",
            [
                clean(A_shift[k, s])
                for s in range(A_shift.cols)
            ]
        )

    print()
    print("  B:")
    for k in range(B_shift.rows):
        print(
            f"    k={k}:",
            [
                clean(B_shift[k, s])
                for s in range(B_shift.cols)
            ]
        )

    print()
    print(
        "  A binomial-removed rank =",
        A_shift.rank()
    )

    print(
        "  B binomial-removed rank =",
        B_shift.rank()
    )

    print(
        "  A binomial-removed rank-1 =",
        rank_one(A_shift)
    )

    print(
        "  B binomial-removed rank-1 =",
        rank_one(B_shift)
    )

    # --------------------------------------------------------------------------
    # OFFSET POLYNOMIAL LAWS
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. EXACT OFFSET POLYNOMIAL LAWS")
    print("=" * 78)

    for name, M in (
        ("A diagonal-normalized", A_diag),
        ("B diagonal-normalized", B_diag),
        ("A binomial-removed", A_shift),
        ("B binomial-removed", B_shift),
    ):
        print()
        print(f"  {name}")

        polys = column_polynomial(M)

        for s, P in polys.items():
            print(
                f"    shift={s}: degree={degree(P,j)}"
            )
            print(
                f"      P_s(k) = {sp.factor(P)}"
            )

    # --------------------------------------------------------------------------
    # SUCCESSIVE SHIFT COMPARISON
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. SHIFT-TO-SHIFT RATIO AUDIT")
    print("=" * 78)

    for name, M in (
        ("A diagonal-normalized", A_diag),
        ("B diagonal-normalized", B_diag),
        ("A binomial-removed", A_shift),
        ("B binomial-removed", B_shift),
    ):
        print()
        print(f"  {name}")

        for s in range(M.cols - 1):

            values = []

            for k in range(M.rows):
                a = clean(M[k, s])
                b = clean(M[k, s + 1])

                if a == 0:
                    values.append(None)
                else:
                    values.append(
                        clean(b / a)
                    )

            print(
                f"    shift {s}->{s+1}:"
            )
            print(
                "      ratios =",
                values
            )

    # --------------------------------------------------------------------------
    # ROW-TO-ROW COMPARISON
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. ROW-SHIFT STABILITY TEST")
    print("=" * 78)

    for name, M in (
        ("A diagonal-normalized", A_diag),
        ("B diagonal-normalized", B_diag),
        ("A binomial-removed", A_shift),
        ("B binomial-removed", B_shift),
    ):
        print()
        print(f"  {name}")

        for k in range(M.rows - 1):

            same = True

            for s in range(M.cols):
                if clean(
                    M[k, s] - M[k + 1, s]
                ) != 0:
                    same = False
                    break

            print(
                f"    row {k} == row {k+1}: {same}"
            )

    # --------------------------------------------------------------------------
    # CROSS-CHANNEL OFFSET COMPARISON
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. CROSS-CHANNEL OFFSET COMPARISON")
    print("=" * 78)

    for s in range(4):

        values = []

        for k in range(5):

            a = clean(A_shift[k, s])
            b = clean(B_shift[k, s])

            if b != 0:
                values.append(
                    (k, clean(a / b))
                )

        print()
        print(
            f"  shift={s}: A/B binomial-removed"
        )

        print(
            "    ",
            values
        )

    # --------------------------------------------------------------------------
    # GENERATING POLYNOMIAL IN OFFSET s
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. OFFSET-GENERATING POLYNOMIALS")
    print("=" * 78)

    for name, M in (
        ("A diagonal-normalized", A_diag),
        ("B diagonal-normalized", B_diag),
        ("A binomial-removed", A_shift),
        ("B binomial-removed", B_shift),
    ):
        print()
        print(f"  {name}")

        for k in range(M.rows):

            G = clean(
                sum(
                    M[k, s] * x**s
                    for s in range(M.cols)
                )
            )

            print(
                f"    k={k}:",
                sp.factor(G)
            )

    # --------------------------------------------------------------------------
    # EXACT ORIGINAL RECONSTRUCTION
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. EXACT RECONSTRUCTION")
    print("=" * 78)

    reconstruction_ok = True

    for name, sequence, C, max_k in (
        ("A", Aseq, A_C, 6),
        ("B", Bseq, B_C, 5),
    ):
        print()
        print(f"  {name}")

        max_r = 9 if name == "A" else 8

        for k in range(max_k + 1):

            P = coefficient_values(
                sequence,
                k
            )

            reconstructed = clean(
                sum(
                    C[k].get(
                        r,
                        sp.Integer(0)
                    ) * falling(j, r)
                    for r in range(max_r)
                )
            )

            ok = (
                clean(
                    P - reconstructed
                ) == 0
            )

            reconstruction_ok &= ok

            print(
                f"    k={k}: {ok}"
            )

    # --------------------------------------------------------------------------
    # FINAL
    # --------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("11. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
  The exact triangle is

      P_k(j) = sum_{s>=0} C[k,k+s] j_(k+s).

  Experiment 92 showed that dividing by binom(k+s,k) does not
  collapse the coefficient matrix.

  This experiment therefore asks a more precise question:

      after removing

          1. the diagonal C[k,k], and
          2. the binomial multiplicity binom(k+s,k),

      does the remaining quantity depend primarily on the offset

          s = r-k

      rather than independently on k and r?

  The key quantity is

      F[k,s]
        = C[k,k+s]
          / ( C[k,k] * binom(k+s,k) ).

  A rank-one F matrix would imply

      F[k,s] = f(k) g(s),

  while identical or nearly identical rows would indicate that
  the construction propagates by a fixed offset kernel.

  A low-degree polynomial in k for every fixed s would indicate
  a controlled indexed family.

  Conversely, a high-rank matrix with unrelated offset laws is
  evidence that the triangular coefficients do not reduce to a
  simple Toeplitz/binomial construction at this level.

  Everything here is exact over QQ and uses only the finite
  observed layer range.

  No extrapolation is performed.
  No factorization algorithm is inferred.
        """
    )

    print()
    print("=" * 78)
    print("12. FINAL EXACTNESS")
    print("=" * 78)

    print(
        "  A reconstruction =",
        reconstruction_ok
    )

    print(
        "  B reconstruction =",
        reconstruction_ok
    )

    print(
        "  A offset matrix rank =",
        A_shift.rank()
    )

    print(
        "  B offset matrix rank =",
        B_shift.rank()
    )

    print(
        "  A offset rank-1 =",
        rank_one(A_shift)
    )

    print(
        "  B offset rank-1 =",
        rank_one(B_shift)
    )

    print(
        "  ALL BASIC CHECKS PASS =",
        reconstruction_ok
    )

    print()
    print(
        "EXPERIMENT 93 COMPLETE"
    )


if __name__ == "__main__":
    main()

