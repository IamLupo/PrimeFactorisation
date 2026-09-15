#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 160 — EXACT COFACTOR/ADJUGATE ORIENTATION / p-ADIC KERNEL AUDIT
==============================================================================

Experiment 159 tested C*A, where C was constructed as the cofactor matrix

    C_ij = (-1)^(i+j) det(A without row i and column j).

The correct adjugate is

    adj(A) = C^T.

Hence the exact identities are

    C^T A = det(A) I
    A C^T = det(A) I.

Experiment 160 corrects only this orientation issue.

It tests:

    1. exact C^T A identity;
    2. exact A C^T identity;
    3. the normalized identities after division by Delta_11;
    4. the distinguished normalized cofactor row/column;
    5. left-kernel and right-kernel behavior modulo powers of 2;
    6. mod-2 row/column spaces of the normalized adjugate;
    7. comparison with the mod-2 left/right nullspaces of A.

This is an orientation audit, not a new fitted model.

No recurrence fitting.
No SymPy.
No floating point.
No extrapolation.
No connection to the original (p,q)-kernel is asserted.
"""

from __future__ import annotations

from fractions import Fraction
import sys


# ============================================================================
# DATA
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
N = 12


# ============================================================================
# HELPERS
# ============================================================================

def q_value(p, r):
    if r < 0 or r >= len(Q[p]):
        return 0
    return Q[p][r]


def basis(p, d):
    return [
        1,
        d,
        d * d,
        p,
        p * d,
        p * p,
    ]


def valuation(n, p):
    n = abs(int(n))

    if n == 0:
        return float("inf")

    v = 0

    while n % p == 0:
        n //= p
        v += 1

    return v


# ============================================================================
# BAREISS
# ============================================================================

def bareiss_det(A):

    if not A:
        return 1

    M = [[int(x) for x in row] for row in A]
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
# SYSTEM
# ============================================================================

def build_system():

    M = []

    for p, pnext in TRANSITIONS:

        for r in range(D[p] + 1):

            d = D[p] - r

            q0 = q_value(p, r)
            q1 = q_value(p, r + 1)

            b = basis(p, d)

            row = [0] * 14

            for j in range(6):
                row[j] = b[j] * q0

            for j in range(6):
                row[6 + j] = b[j] * q1

            if d == 0:
                row[12] = 1
                row[13] = p

            M.append(row)

    return M


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
        raise ArithmeticError("Boundary block singular.")

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
        out = []

        for j in range(12):

            value = (
                Fraction(row[j])
                + Fraction(row[12]) * h_coeff[0][j]
                + Fraction(row[13]) * h_coeff[1][j]
            )

            if value.denominator != 1:
                raise ArithmeticError(
                    "Non-integral F-core entry."
                )

            out.append(value.numerator)

        core.append(out)

    return core


# ============================================================================
# COFACTOR MATRIX
# ============================================================================

def maximal_minor(A, omit_row, omit_col):

    rows = [i for i in range(len(A)) if i != omit_row]
    cols = [j for j in range(len(A[0])) if j != omit_col]

    return [
        [A[i][j] for j in cols]
        for i in rows
    ]


def cofactor_matrix(A):

    n = len(A)
    C = []

    for i in range(n):

        row = []

        for j in range(n):

            value = bareiss_det(
                maximal_minor(A, i, j)
            )

            if (i + j) & 1:
                value = -value

            row.append(value)

        C.append(row)

    return C


def normalize(C):

    return [
        [
            x // DELTA_11
            if x % DELTA_11 == 0
            else (_raise_division())
            for x in row
        ]
        for row in C
    ]


def _raise_division():
    raise ArithmeticError(
        "Delta_11 does not divide cofactor."
    )


# ============================================================================
# MATRIX OPERATIONS
# ============================================================================

def transpose(A):

    return [
        [
            A[i][j]
            for i in range(len(A))
        ]
        for j in range(len(A[0]))
    ]


def matmul(A, B):

    rows = len(A)
    inner = len(B)
    cols = len(B[0])

    out = [
        [0 for _ in range(cols)]
        for _ in range(rows)
    ]

    for i in range(rows):

        for k in range(inner):

            x = A[i][k]

            if x == 0:
                continue

            for j in range(cols):
                out[i][j] += x * B[k][j]

    return out


# ============================================================================
# MOD-2
# ============================================================================

def rref_mod2(A):

    M = [
        [int(x) & 1 for x in row]
        for row in A
    ]

    rows = len(M)
    cols = len(M[0])

    pivots = []
    pivot_row = 0

    for col in range(cols):

        pivot = None

        for r in range(pivot_row, rows):
            if M[r][col]:
                pivot = r
                break

        if pivot is None:
            continue

        M[pivot_row], M[pivot] = (
            M[pivot],
            M[pivot_row],
        )

        for r in range(rows):

            if r == pivot_row:
                continue

            if M[r][col]:

                for j in range(col, cols):
                    M[r][j] ^= M[pivot_row][j]

        pivots.append(col)
        pivot_row += 1

        if pivot_row == rows:
            break

    return M, pivots


def rank_mod2(A):
    _, pivots = rref_mod2(A)
    return len(pivots)


def row_space_mod2(A):

    R, _ = rref_mod2(A)

    return [
        row
        for row in R
        if any(row)
    ]


def nullspace_mod2(A):

    R, pivots = rref_mod2(A)

    n = len(A[0])

    pivot_set = set(pivots)

    free = [
        j
        for j in range(n)
        if j not in pivot_set
    ]

    basis = []

    for f in free:

        v = [0] * n
        v[f] = 1

        for i, pc in enumerate(pivots):

            if R[i][f]:
                v[pc] = 1

        basis.append(v)

    return basis


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 160 — EXACT COFACTOR/ADJUGATE ORIENTATION / "
        "p-ADIC KERNEL AUDIT"
    )
    print("=" * 78)

    M = build_system()
    A = build_core(M)
    C = cofactor_matrix(A)
    R = normalize(C)

    Adj = transpose(C)
    AdjR = transpose(R)

    det_A = bareiss_det(A)

    d12 = det_A // DELTA_11

    # ------------------------------------------------------------------
    # 1. CORRECT ADJUGATE IDENTITIES
    # ------------------------------------------------------------------

    CT_A = matmul(Adj, A)
    A_CT = matmul(A, Adj)

    normalized_left = matmul(AdjR, A)
    normalized_right = matmul(A, AdjR)

    expected = [
        [0 for _ in range(N)]
        for _ in range(N)
    ]

    for i in range(N):
        expected[i][i] = det_A

    expected_normalized = [
        [0 for _ in range(N)]
        for _ in range(N)
    ]

    for i in range(N):
        expected_normalized[i][i] = d12

    print()
    print("=" * 78)
    print("1. CORRECT ADJUGATE IDENTITIES")
    print("=" * 78)

    print(
        f"  C^T*A_exact={CT_A == expected}"
    )

    print(
        f"  A*C^T_exact={A_CT == expected}"
    )

    print(
        f"  R_adj*A_exact="
        f"{normalized_left == expected_normalized}"
    )

    print(
        f"  A*R_adj_exact="
        f"{normalized_right == expected_normalized}"
    )

    # ------------------------------------------------------------------
    # 2. DISTINGUISHED ADJUGATE ROW
    # ------------------------------------------------------------------

    r = AdjR[5]

    residual = normalized_left[5]

    print()
    print("=" * 78)
    print("2. DISTINGUISHED ADJUGATE ROW")
    print("=" * 78)

    print(
        f"  row_index=5"
    )

    print(
        f"  AdjR_row={r}"
    )

    print(
        f"  AdjR_row*A={residual}"
    )

    # ------------------------------------------------------------------
    # 3. EXACT p-ADIC PRECISION
    # ------------------------------------------------------------------

    residual_v2 = [
        valuation(x, 2)
        for x in residual
    ]

    finite = [
        x for x in residual_v2
        if x != float("inf")
    ]

    min_v = min(finite)

    non_dist = [
        residual[j]
        for j in range(N)
        if j != 5
    ]

    non_dist_min = min(
        valuation(x, 2)
        for x in non_dist
    )

    print()
    print("=" * 78)
    print("3. 2-ADIC KERNEL PRECISION")
    print("=" * 78)

    print(
        f"  d12_v2={valuation(d12,2)}"
    )

    print(
        f"  residual_min_v2={min_v}"
    )

    print(
        f"  non_distinguished_min_v2={non_dist_min}"
    )

    for e in range(1, 12):

        full_zero = all(
            x % (2 ** e) == 0
            for x in residual
        )

        non_dist_zero = all(
            x % (2 ** e) == 0
            for x in non_dist
        )

        print(
            f"  mod=2^{e}: "
            f"full_zero={full_zero} "
            f"non_dist_zero={non_dist_zero}"
        )

    # ------------------------------------------------------------------
    # 4. MOD-2 ADJUGATE VS A NULLSPACES
    # ------------------------------------------------------------------

    A2 = [
        [x & 1 for x in row]
        for row in A
    ]

    AdjR2 = [
        [x & 1 for x in row]
        for row in AdjR
    ]

    left_kernel = nullspace_mod2(A2)
    right_kernel = nullspace_mod2(
        transpose(A2)
    )

    adj_row_space = row_space_mod2(
        AdjR2
    )

    adj_col_space = row_space_mod2(
        transpose(AdjR2)
    )

    print()
    print("=" * 78)
    print("4. MOD-2 KERNEL / ADJUGATE SPACE")
    print("=" * 78)

    print(
        f"  rank_A_mod2={rank_mod2(A2)}"
    )

    print(
        f"  left_nullity={len(left_kernel)}"
    )

    print(
        f"  right_nullity={len(right_kernel)}"
    )

    print(
        f"  adjugate_rowspace_rank="
        f"{len(adj_row_space)}"
    )

    print(
        f"  adjugate_columnspace_rank="
        f"{len(adj_col_space)}"
    )

    # ------------------------------------------------------------------
    # 5. EXACT INCLUSION TESTS
    # ------------------------------------------------------------------

    def in_span(v, B):

        if not B:
            return not any(v)

        before = rank_mod2(B)
        after = rank_mod2(B + [v])

        return before == after

    print()
    print("=" * 78)
    print("5. EXACT ADJUGATE / NULLSPACE INCLUSION")
    print("=" * 78)

    all_rows_in_left = all(
        in_span(row, left_kernel)
        for row in AdjR2
        if any(row)
    )

    all_cols_in_right = all(
        in_span(col, right_kernel)
        for col in AdjR2
        if any(col)
    )

    print(
        f"  all_nonzero_adj_rows_in_left_nullspace="
        f"{all_rows_in_left}"
    )

    print(
        f"  all_nonzero_adj_columns_in_right_nullspace="
        f"{all_cols_in_right}"
    )

    # ------------------------------------------------------------------
    # 6. SUPPORT
    # ------------------------------------------------------------------

    active_rows = [
        i
        for i, row in enumerate(AdjR2)
        if any(row)
    ]

    active_cols = [
        j
        for j in range(N)
        if any(
            AdjR2[i][j]
            for i in range(N)
        )
    ]

    print()
    print("=" * 78)
    print("6. ADJUGATE SUPPORT MOD 2")
    print("=" * 78)

    print(
        f"  active_rows={active_rows}"
    )

    print(
        f"  active_columns={active_cols}"
    )

    print(
        f"  distinguished_adj_row_support="
        f"{[j for j,x in enumerate(AdjR2[5]) if x]}"
    )

    # ------------------------------------------------------------------
    # 7. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
Experiment 159 used the cofactor matrix C as though it were the
adjugate. That orientation was incorrect.

The exact convention is:

    cofactor matrix:
        C;

    adjugate:
        C^T.

Therefore the correct identities are

    C^T A = det(A) I,
    A C^T = det(A) I.

Experiment 160 tests these identities directly.

If they pass, then the previously observed mod-2 rank-one behavior
has the standard adjugate interpretation:

    every nonzero row of Adj(A) mod 2
        lies in the left kernel of A mod 2;

    every nonzero column of Adj(A) mod 2
        lies in the right kernel of A mod 2.

Likewise the normalized identity gives the exact p-adic statement

    AdjR * A = d12 I.

The terminal valuation v2(d12)=9 then gives an immediate and exact
2-adic annihilation threshold.

This is a basic algebraic identity rather than a speculative
mechanism, so it should be treated separately from the structural
interpretation of the original data.
"""
    )

    # ------------------------------------------------------------------
    # 8. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        CT_A == expected
        and A_CT == expected
        and normalized_left == expected_normalized
        and normalized_right == expected_normalized
        and len(AdjR2) == 12
        and rank_mod2(AdjR2) == 1
        and non_dist_min >= 9
    )

    print()
    print("=" * 78)
    print("8. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  adjugate_left_identity_exact="
        f"{CT_A == expected}"
    )

    print(
        f"  adjugate_right_identity_exact="
        f"{A_CT == expected}"
    )

    print(
        f"  normalized_left_identity_exact="
        f"{normalized_left == expected_normalized}"
    )

    print(
        f"  normalized_right_identity_exact="
        f"{normalized_right == expected_normalized}"
    )

    print(
        f"  mod2_adjugate_rank_one="
        f"{rank_mod2(AdjR2) == 1}"
    )

    print(
        f"  p_adic_precision_match="
        f"{non_dist_min >= 9}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 160 COMPLETE")


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

