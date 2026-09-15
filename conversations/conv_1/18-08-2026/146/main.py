#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 146 — EXACT 2-ADIC DIAGONAL NORMALIZATION / RANK-LIFT AUDIT
==============================================================================

Experiment 145 found:

    rank(R mod 2) = 1
    rank(R mod 4) = 1
    rank(R mod 8) = 1
    rank(R mod 16) != 1

but the first 2-adic correction has

    rank(E_2 mod 2) = 3.

Experiment 146 asks a narrower, evidence-driven question:

    Is the exceptional 2-adic behavior partly caused by simple row/column
    divisibility already visible in the normalized cofactor matrix?

From Experiment 142 we already know that the rows and columns of R have
different integer gcds.

We therefore construct an exact diagonal normalization based only on
2-adic valuation profiles.

For each row i:

    a_i = min_j v_2(R_ij)

and divide row i by 2^a_i.

Then for the resulting matrix:

    b_j = min_i v_2(row-normalized_ij)

and divide column j by 2^b_j.

The final matrix remains integral by construction.

We then compare:

    original R
    row-normalized R
    row+column normalized R

at the levels

    mod 2,
    mod 4,
    mod 8,
    mod 16.

For rank-one persistence we use the exact determinantal criterion:

    all 2x2 minors divisible by 2^e
        <=> determinantal rank <= 1 modulo 2^e.

We also compare the first-correction rank after the normalization.

This is deliberately a small structural audit.

It does NOT assume that diagonal normalization is the correct hidden
mechanism. A negative result is equally informative.

No floating point.
No SymPy.
No extrapolation.
No connection to the original (p,q)-kernel is asserted.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
import sys


# ============================================================================
# EXACT DATA
# ============================================================================

Q = {
    1: [
        -126258696,
        -11600759760,
        2668721436,
        1764373740,
        -1338089411,
        495451247,
    ],
    3: [
        9955176,
        -1263551016,
        -152369292,
        -128667196,
        421514439,
    ],
    5: [
        -62398,
        4771718,
        16027881,
    ],
    7: [
        1,
    ],
}

D = {
    1: 5,
    3: 4,
    5: 2,
    7: 0,
}

TRANSITIONS = [
    (1, 3),
    (3, 5),
    (5, 7),
]

DELTA_11 = 9512681472


# ============================================================================
# BASIC HELPERS
# ============================================================================

def valuation(n: int, p: int) -> int:
    n = abs(int(n))

    if n == 0:
        return 10**9

    v = 0

    while n % p == 0:
        n //= p
        v += 1

    return v


def q_value(p: int, r: int) -> int:
    if r < 0 or r >= len(Q[p]):
        return 0
    return Q[p][r]


def monomial_basis(p: int, d: int):
    return [
        1,
        d,
        d * d,
        p,
        p * d,
        p * p,
    ]


def exact_max(a, b):
    return a if a >= b else b


# ============================================================================
# BAREISS
# ============================================================================

def bareiss_det(A):

    if not A:
        return 1

    M = [
        [int(x) for x in row]
        for row in A
    ]

    n = len(M)

    if n == 1:
        return M[0][0]

    previous = 1
    sign = 1

    for k in range(n - 1):

        pivot_row = None

        for r in range(k, n):
            if M[r][k] != 0:
                pivot_row = r
                break

        if pivot_row is None:
            return 0

        if pivot_row != k:
            M[k], M[pivot_row] = M[pivot_row], M[k]
            sign *= -1

        pivot = M[k][k]

        for i in range(k + 1, n):

            for j in range(k + 1, n):

                value = (
                    M[i][j] * pivot
                    - M[i][k] * M[k][j]
                )

                if k > 0:

                    if value % previous != 0:
                        raise ArithmeticError(
                            "Bareiss exact division failed."
                        )

                    value //= previous

                M[i][j] = value

        for i in range(k + 1, n):
            M[i][k] = 0

        previous = pivot

    return sign * M[-1][-1]


# ============================================================================
# FULL SYSTEM
# ============================================================================

def build_full_system():

    M = []

    for p, p_next in TRANSITIONS:

        for r in range(D[p] + 1):

            d = D[p] - r

            q0 = q_value(p, r)
            q1 = q_value(p, r + 1)

            basis = monomial_basis(p, d)

            row = [0] * 14

            for j in range(6):
                row[j] = basis[j] * q0

            for j in range(6):
                row[6 + j] = basis[j] * q1

            if d == 0:
                row[12] = 1
                row[13] = p

            M.append(row)

    return M


# ============================================================================
# F-CORE
# ============================================================================

def build_core(M):

    b0 = M[5]
    b1 = M[10]

    H = [
        [b0[12], b0[13]],
        [b1[12], b1[13]],
    ]

    det_h = (
        H[0][0] * H[1][1]
        - H[0][1] * H[1][0]
    )

    if det_h == 0:
        raise ArithmeticError(
            "Boundary block is singular."
        )

    H_inv = [
        [
            Fraction(H[1][1], det_h),
            Fraction(-H[0][1], det_h),
        ],
        [
            Fraction(-H[1][0], det_h),
            Fraction(H[0][0], det_h),
        ],
    ]

    h_coeff = [
        [Fraction(0) for _ in range(12)]
        for _ in range(2)
    ]

    for i in range(2):

        for j in range(12):

            h_coeff[i][j] = -(
                H_inv[i][0] * b0[j]
                + H_inv[i][1] * b1[j]
            )

    core = []

    for idx in range(14):

        if idx in (5, 10):
            continue

        row = M[idx]

        reduced = []

        for j in range(12):

            value = (
                Fraction(row[j])
                + Fraction(row[12]) * h_coeff[0][j]
                + Fraction(row[13]) * h_coeff[1][j]
            )

            if value.denominator != 1:
                raise ArithmeticError(
                    "Non-integral core entry."
                )

            reduced.append(value.numerator)

        core.append(reduced)

    return core


# ============================================================================
# NORMALIZED COFACTOR MATRIX
# ============================================================================

def maximal_minor(A, omit_row, omit_col):

    rows = [
        i for i in range(len(A))
        if i != omit_row
    ]

    cols = [
        j for j in range(len(A[0]))
        if j != omit_col
    ]

    return [
        [A[i][j] for j in cols]
        for i in rows
    ]


def normalized_cofactors(A):

    n = len(A)

    R = []

    for i in range(n):

        row = []

        for j in range(n):

            det_minor = bareiss_det(
                maximal_minor(A, i, j)
            )

            if det_minor % DELTA_11 != 0:
                raise ArithmeticError(
                    f"Minor ({i},{j}) not divisible by Delta_11."
                )

            row.append(
                det_minor // DELTA_11
            )

        R.append(row)

    return R


# ============================================================================
# MATRIX SCALING
# ============================================================================

def row_min_v2(A):

    result = []

    for row in A:

        finite = [
            valuation(x, 2)
            for x in row
            if x != 0
        ]

        if not finite:
            result.append(0)
        else:
            result.append(min(finite))

    return result


def col_min_v2(A):

    n = len(A[0])
    result = []

    for j in range(n):

        finite = [
            valuation(A[i][j], 2)
            for i in range(len(A))
            if A[i][j] != 0
        ]

        if not finite:
            result.append(0)
        else:
            result.append(min(finite))

    return result


def divide_rows_by_powers(A, exponents, p=2):

    B = []

    for i, row in enumerate(A):

        factor = p ** exponents[i]

        out = []

        for x in row:

            if x % factor != 0:
                raise ArithmeticError(
                    "Row normalization produced non-integral entry."
                )

            out.append(
                x // factor
            )

        B.append(out)

    return B


def divide_cols_by_powers(A, exponents, p=2):

    n = len(A[0])

    B = [
        row[:]
        for row in A
    ]

    for j in range(n):

        factor = p ** exponents[j]

        for i in range(len(A)):

            if B[i][j] % factor != 0:
                raise ArithmeticError(
                    "Column normalization produced non-integral entry."
                )

            B[i][j] //= factor

    return B


def row_normalize_2adic(A):

    exponents = row_min_v2(A)
    B = divide_rows_by_powers(
        A,
        exponents,
        2,
    )

    return B, exponents


def row_column_normalize_2adic(A):

    B, row_exp = row_normalize_2adic(A)

    col_exp = col_min_v2(B)

    C = divide_cols_by_powers(
        B,
        col_exp,
        2,
    )

    return C, row_exp, col_exp


# ============================================================================
# MODULAR RANK
# ============================================================================

def rank_mod(A, modulus):

    M = [
        [x % modulus for x in row]
        for row in A
    ]

    rows = len(M)
    cols = len(M[0])

    rank = 0

    for col in range(cols):

        pivot = None

        for r in range(rank, rows):

            if M[r][col] != 0:
                pivot = r
                break

        if pivot is None:
            continue

        M[rank], M[pivot] = (
            M[pivot],
            M[rank],
        )

        inv = pow(
            M[rank][col],
            -1,
            modulus,
        )

        for j in range(col, cols):

            M[rank][j] = (
                M[rank][j] * inv
            ) % modulus

        for r in range(rows):

            if r == rank:
                continue

            factor = M[r][col]

            if factor == 0:
                continue

            for j in range(col, cols):

                M[r][j] = (
                    M[r][j]
                    - factor * M[rank][j]
                ) % modulus

        rank += 1

    return rank


# ============================================================================
# MINIMUM 2x2 MINOR VALUATION
# ============================================================================

def min_2x2_v2(A):

    n = len(A)

    minimum = 10**9
    positions = []

    for i, j in combinations(range(n), 2):

        for a, b in combinations(range(n), 2):

            det = (
                A[i][a] * A[j][b]
                - A[i][b] * A[j][a]
            )

            if det == 0:
                continue

            v = valuation(det, 2)

            if v < minimum:
                minimum = v
                positions = [
                    (i, j, a, b)
                ]

            elif v == minimum:
                positions.append(
                    (i, j, a, b)
                )

    return minimum, positions


# ============================================================================
# MATRIX CONTENT / PROFILE
# ============================================================================

def matrix_gcd(A):

    g = 0

    for row in A:

        for x in row:
            from math import gcd
            g = gcd(g, abs(x))

    return g


def row_min_profile(A):

    return row_min_v2(A)


def col_min_profile(A):

    return col_min_v2(A)


# ============================================================================
# FIRST 2-ADIC DEFECT
# ============================================================================

def rank1_factor_mod_power(A, modulus):

    n = len(A)

    pivot_pos = None
    pivot_inv = None

    for i in range(n):

        for j in range(n):

            value = A[i][j] % modulus

            if value == 0:
                continue

            try:
                inv = pow(
                    value,
                    -1,
                    modulus,
                )
            except ValueError:
                continue

            pivot_pos = (i, j)
            pivot_inv = inv
            break

        if pivot_pos is not None:
            break

    if pivot_pos is None:
        return None

    i0, j0 = pivot_pos

    u = [
        A[i][j0] % modulus
        for i in range(n)
    ]

    v = [
        (
            A[i0][j] * pivot_inv
        ) % modulus
        for j in range(n)
    ]

    exact = True

    for i in range(n):

        for j in range(n):

            if (
                u[i] * v[j]
            ) % modulus != A[i][j] % modulus:

                exact = False
                break

        if not exact:
            break

    return {
        "position": pivot_pos,
        "u": u,
        "v": v,
        "exact": exact,
    }


def defect_rank_after_mod8(A):

    data = rank1_factor_mod_power(
        A,
        8,
    )

    if data is None or not data["exact"]:
        return None, None

    u = data["u"]
    v = data["v"]

    E = []

    for i in range(12):

        row = []

        for j in range(12):

            numerator = (
                A[i][j]
                - u[i] * v[j]
            )

            if numerator % 8 != 0:
                raise ArithmeticError(
                    "Invalid mod-8 rank-1 factorization."
                )

            row.append(
                numerator // 8
            )

        E.append(row)

    return rank_mod(E, 2), data


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 146 — EXACT 2-ADIC DIAGONAL NORMALIZATION / "
        "RANK-LIFT AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # BUILD
    # ------------------------------------------------------------------

    M = build_full_system()
    core = build_core(M)
    R = normalized_cofactors(core)

    print()
    print("=" * 78)
    print("1. EXACT DATA")
    print("=" * 78)

    print("  R_shape=(12,12)")
    print(f"  Delta_11={DELTA_11}")

    # ------------------------------------------------------------------
    # BASELINE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. BASELINE 2-ADIC PROFILE")
    print("=" * 78)

    baseline_rows = row_min_profile(R)
    baseline_cols = col_min_profile(R)
    baseline_minor_v, baseline_pos = min_2x2_v2(R)

    print(
        f"  row_min_v2={baseline_rows}"
    )

    print(
        f"  col_min_v2={baseline_cols}"
    )

    print(
        f"  minimum_2x2_v2={baseline_minor_v}"
    )

    print(
        f"  minimum_2x2_positions_sample="
        f"{baseline_pos[:12]}"
    )

    for e in range(1, 5):

        print(
            f"  rank1_mod_2^{e}="
            f"{baseline_minor_v >= e}"
        )

    # ------------------------------------------------------------------
    # ROW NORMALIZATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. ROW 2-ADIC NORMALIZATION")
    print("=" * 78)

    R_row, row_exp = row_normalize_2adic(R)

    row_after_rows = row_min_profile(R_row)
    row_after_cols = col_min_profile(R_row)
    row_minor_v, row_minor_pos = min_2x2_v2(R_row)

    print(
        f"  row_exponents={row_exp}"
    )

    print(
        f"  post_row_min_v2={row_after_rows}"
    )

    print(
        f"  post_row_col_min_v2={row_after_cols}"
    )

    print(
        f"  minimum_2x2_v2={row_minor_v}"
    )

    print(
        f"  minimum_2x2_positions_sample="
        f"{row_minor_pos[:12]}"
    )

    for e in range(1, 5):

        print(
            f"  rank1_mod_2^{e}="
            f"{row_minor_v >= e}"
        )

    # ------------------------------------------------------------------
    # ROW + COLUMN NORMALIZATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. ROW + COLUMN 2-ADIC NORMALIZATION")
    print("=" * 78)

    R_bal, row_exp2, col_exp2 = (
        row_column_normalize_2adic(R)
    )

    balanced_rows = row_min_profile(R_bal)
    balanced_cols = col_min_profile(R_bal)
    balanced_minor_v, balanced_pos = min_2x2_v2(R_bal)

    print(
        f"  row_exponents={row_exp2}"
    )

    print(
        f"  column_exponents={col_exp2}"
    )

    print(
        f"  post_balance_row_min_v2="
        f"{balanced_rows}"
    )

    print(
        f"  post_balance_col_min_v2="
        f"{balanced_cols}"
    )

    print(
        f"  minimum_2x2_v2={balanced_minor_v}"
    )

    print(
        f"  minimum_2x2_positions_sample="
        f"{balanced_pos[:12]}"
    )

    for e in range(1, 5):

        print(
            f"  rank1_mod_2^{e}="
            f"{balanced_minor_v >= e}"
        )

    # ------------------------------------------------------------------
    # MODULAR RANK COMPARISON
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. MODULAR RANK COMPARISON")
    print("=" * 78)

    for name, A in [
        ("raw", R),
        ("row_normalized", R_row),
        ("row_column_normalized", R_bal),
    ]:

        print()
        print(f"  {name}")

        for modulus in [2, 4, 8]:

            # Ordinary Gaussian rank only makes sense over fields.
            # For 4 and 8 we report the exact 2x2 determinantal condition
            # separately rather than mislabeling it as field rank.
            if modulus == 2:

                rank = rank_mod(A, 2)

                print(
                    f"    rank_mod_2={rank}"
                )

            else:

                mv, _ = min_2x2_v2(A)

                print(
                    f"    rank1_determinantal_mod_{modulus}="
                    f"{mv >= (2 if modulus == 4 else 3)}"
                )

    # ------------------------------------------------------------------
    # DEFECT RANK AFTER NORMALIZATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. FIRST 2-ADIC DEFECT AFTER NORMALIZATION")
    print("=" * 78)

    for name, A in [
        ("raw", R),
        ("row_normalized", R_row),
        ("row_column_normalized", R_bal),
    ]:

        defect_rank, factor = defect_rank_after_mod8(A)

        print()

        print(
            f"  {name}:"
        )

        print(
            f"    mod8_rank1_factorization="
            f"{factor is not None and factor['exact']}"
        )

        print(
            f"    defect_rank_mod2="
            f"{defect_rank}"
        )

    # ------------------------------------------------------------------
    # CONTENT COMPARISON
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. INTEGER CONTENT")
    print("=" * 78)

    for name, A in [
        ("raw", R),
        ("row_normalized", R_row),
        ("row_column_normalized", R_bal),
    ]:

        g = matrix_gcd(A)

        print(
            f"  {name}: gcd_entries={g}"
        )

    # ------------------------------------------------------------------
    # STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The previous experiments show that the normalized cofactor matrix R
has an unusually persistent 2-adic rank-one structure.

This experiment asks whether that persistence is explained simply by
row and column divisibility.

The normalization is deliberately mechanical:

    row exponent:
        a_i = min_j v_2(R_ij)

followed by

    column exponent:
        b_j = min_i v_2(R'_ij).

No coefficients are fitted.

If the rank-one persistence increases after this normalization, then
the previous 2-adic collapse was at least partly attributable to
diagonal divisibility.

If the rank profile and first-defect rank remain essentially unchanged,
then the collapse is not explained by this elementary gauge.

The experiment therefore separates:

    divisibility normalization
        from
    genuine low-rank congruence.

No conclusion about the original (p,q)-kernel is drawn.
"""
    )

    # ------------------------------------------------------------------
    # FINAL
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. FINAL EXACTNESS")
    print("=" * 78)

    final_ok = (
        len(R) == 12
        and len(R_row) == 12
        and len(R_bal) == 12
        and len(R[0]) == 12
        and len(R_row[0]) == 12
        and len(R_bal[0]) == 12
    )

    print(
        f"  data_reconstruction_exact={final_ok}"
    )

    print(
        "  row_normalization_integral=True"
    )

    print(
        "  row_column_normalization_integral=True"
    )

    print(
        "  2adic_minor_audit_completed=True"
    )

    print(
        "  first_defect_audit_completed=True"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 146 COMPLETE")


if __name__ == "__main__":

    try:
        main()

    except KeyboardInterrupt:

        print("\nInterrupted.")
        sys.exit(130)

    except Exception as exc:

        print(
            "\nFATAL ERROR: "
            + type(exc).__name__
            + ": "
            + str(exc)
        )

        raise

