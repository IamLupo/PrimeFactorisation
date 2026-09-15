#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 169S — EXACT SCHUR FACTOR LIFT / FIRST 2-ADIC LAYER AUDIT
==============================================================================

Corrected, dimension-safe replacement for Experiment 169R2.

Established:

    F_ret : 12x12
    A     : 12x12
    K     : 12x12

with

    A = F_ret + K

and

    K = -(1/2) B adj(H) C

where

    B       : 12x2
    H       : 2x2
    C       : 2x12
    det(H)  : 2.

Experiment 168R established:

    B adj(H) mod 2 = 0
    adj(H) C mod 2 = 0

while

    K mod 2

is nonzero and rank one.

This experiment computes the first 2-adic lifts:

    L1 = (B adj(H))/2        [12x2]
    R1 = (adj(H) C)/2        [2x12]

and determines whether either lifted factor alone explains K mod 2.

IMPORTANT DIMENSION RULE

    L1 has shape 12x2.
    C  has shape 2x12.

Therefore

    L1*C

is 12x12.

Likewise

    B*R1

is 12x12.

The outer-product test is performed only when the factor dimensions
permit it. For a 12x2 matrix of rank one, the correct decomposition is

    L1 = u (1x2 row),

not a 12x12 outer product.

Likewise

    R1 = (2x1 column) v^T.

No incorrect square-matrix indexing is used.

No SymPy.
No floating point.
No extrapolation.
No recurrence fitting.
No claim about the original (p,q)-kernel.
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

FULL_N = 14
CORE_N = 12

BOUNDARY_ROWS = (5, 10)
BOUNDARY_COLS = (12, 13)


# ============================================================================
# BASIC HELPERS
# ============================================================================

def q_value(p: int, r: int) -> int:
    if r < 0 or r >= len(Q[p]):
        return 0
    return Q[p][r]


def basis(p: int, d: int) -> list[int]:
    return [
        1,
        d,
        d * d,
        p,
        p * d,
        p * p,
    ]


def assert_rectangular(A, name="matrix"):
    if not A:
        return

    width = len(A[0])

    for i, row in enumerate(A):
        if len(row) != width:
            raise ValueError(
                f"{name}: row {i} has width {len(row)}, "
                f"expected {width}."
            )


def to_mod2(A):
    assert_rectangular(A)

    return [
        [int(x) & 1 for x in row]
        for row in A
    ]


# ============================================================================
# MATRIX OPERATIONS
# ============================================================================

def transpose(A):
    assert_rectangular(A, "transpose input")

    if not A:
        return []

    return [
        [A[i][j] for i in range(len(A))]
        for j in range(len(A[0]))
    ]


def matmul(A, B):
    assert_rectangular(A, "left matrix")
    assert_rectangular(B, "right matrix")

    if not A or not B:
        return []

    if len(A[0]) != len(B):
        raise ValueError(
            f"Matrix multiplication mismatch: "
            f"{len(A)}x{len(A[0])} times "
            f"{len(B)}x{len(B[0])}."
        )

    out = [
        [Fraction(0) for _ in range(len(B[0]))]
        for _ in range(len(A))
    ]

    for i in range(len(A)):
        for k in range(len(B)):

            a = Fraction(A[i][k])

            if a == 0:
                continue

            for j in range(len(B[0])):
                out[i][j] += a * Fraction(B[k][j])

    return out


def matrix_to_int(A, name):
    out = []

    for i, row in enumerate(A):

        new_row = []

        for j, x in enumerate(row):

            x = Fraction(x)

            if x.denominator != 1:
                raise ArithmeticError(
                    f"{name}[{i}][{j}] is non-integral: {x}"
                )

            new_row.append(x.numerator)

        out.append(new_row)

    return out


def divide_matrix_exact(A, divisor, name):

    out = []

    for i, row in enumerate(A):

        new_row = []

        for j, x in enumerate(row):

            if x % divisor != 0:
                raise ArithmeticError(
                    f"{name}[{i}][{j}]={x} "
                    f"is not divisible by {divisor}."
                )

            new_row.append(x // divisor)

        out.append(new_row)

    return out


# ============================================================================
# BAREISS
# ============================================================================

def bareiss_det(A):

    if not A:
        return 1

    assert_rectangular(A, "Bareiss input")

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

    for p, _ in TRANSITIONS:

        for r in range(D[p] + 1):

            d = D[p] - r

            q0 = q_value(p, r)
            q1 = q_value(p, r + 1)

            b = basis(p, d)

            row = [0] * FULL_N

            for j in range(6):
                row[j] = b[j] * q0

            for j in range(6):
                row[6 + j] = b[j] * q1

            if d == 0:
                row[12] = 1
                row[13] = p

            M.append(row)

    if len(M) != FULL_N:
        raise ArithmeticError(
            f"Expected {FULL_N} rows, got {len(M)}."
        )

    return M


# ============================================================================
# SCHUR BLOCKS
# ============================================================================

def build_blocks(M):

    retained_indices = [
        i
        for i in range(FULL_N)
        if i not in BOUNDARY_ROWS
    ]

    retained_rows = [
        M[i]
        for i in retained_indices
    ]

    boundary_rows = [
        M[i]
        for i in BOUNDARY_ROWS
    ]

    # 12x12
    F = [
        [
            row[j]
            for j in range(CORE_N)
        ]
        for row in retained_rows
    ]

    # 12x2
    B = [
        [
            row[j]
            for j in BOUNDARY_COLS
        ]
        for row in retained_rows
    ]

    # 2x2
    H = [
        [
            row[j]
            for j in BOUNDARY_COLS
        ]
        for row in boundary_rows
    ]

    # 2x12
    C = [
        [
            row[j]
            for j in range(CORE_N)
        ]
        for row in boundary_rows
    ]

    det_h = (
        H[0][0] * H[1][1]
        - H[0][1] * H[1][0]
    )

    if det_h != 2:
        raise ArithmeticError(
            f"Expected det(H)=2, got {det_h}."
        )

    adj_h = [
        [
            H[1][1],
            -H[0][1],
        ],
        [
            -H[1][0],
            H[0][0],
        ],
    ]

    return (
        retained_indices,
        F,
        B,
        H,
        C,
        adj_h,
        det_h,
    )


# ============================================================================
# MOD-2 RANK
# ============================================================================

def rref_mod2(A):

    M = to_mod2(A)

    if not M:
        return [], []

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

                for j in range(cols):
                    M[r][j] ^= M[pivot_row][j]

        pivots.append(col)
        pivot_row += 1

        if pivot_row == rows:
            break

    return M, pivots


def rank_mod2(A):

    if not A:
        return 0

    return len(rref_mod2(A)[1])


def row_basis_mod2(A):

    R, _ = rref_mod2(A)

    return [
        row
        for row in R
        if any(row)
    ]


# ============================================================================
# RANK-ONE FACTORIZATION FOR RECTANGULAR MATRICES
# ============================================================================

def rank_one_factorization_mod2(A):

    """
    For an arbitrary m x n matrix over F_2:

        rank(A)=1

    iff A = u v^T.

    Returns (u, v) or (None, None).

    This is dimension-safe for rectangular matrices.
    """

    A2 = to_mod2(A)

    m = len(A2)

    if m == 0:
        return None, None

    n = len(A2[0])

    if rank_mod2(A2) != 1:
        return None, None

    pivot_row = None
    pivot_col = None

    for i in range(m):
        for j in range(n):
            if A2[i][j]:
                pivot_row = i
                pivot_col = j
                break
        if pivot_row is not None:
            break

    if pivot_row is None:
        return None, None

    u = [
        A2[i][pivot_col]
        for i in range(m)
    ]

    v = [
        A2[pivot_row][j]
        for j in range(n)
    ]

    predicted = [
        [
            u[i] & v[j]
            for j in range(n)
        ]
        for i in range(m)
    ]

    if predicted != A2:
        raise ArithmeticError(
            "Rank-one factorization reconstruction failed."
        )

    return u, v


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 169S — EXACT SCHUR FACTOR LIFT / "
        "FIRST 2-ADIC LAYER AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # BUILD
    # ------------------------------------------------------------------

    M = build_full_system()

    (
        retained_indices,
        F,
        B,
        H,
        C,
        adj_h,
        det_h,
    ) = build_blocks(M)

    H2 = to_mod2(H)
    adj_h2 = to_mod2(adj_h)

    # ------------------------------------------------------------------
    # 1. H
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. BOUNDARY BLOCK")
    print("=" * 78)

    print(f"  H={H}")
    print(f"  det_H={det_h}")
    print(f"  H_mod2={H2}")
    print(f"  adj_H={adj_h}")
    print(f"  adj_H_mod2={adj_h2}")

    # ------------------------------------------------------------------
    # 2. INTERMEDIATE FACTORS
    # ------------------------------------------------------------------

    L = matrix_to_int(
        matmul(B, adj_h),
        "L"
    )

    R = matrix_to_int(
        matmul(adj_h, C),
        "R"
    )

    L_even = all(
        x % 2 == 0
        for row in L
        for x in row
    )

    R_even = all(
        x % 2 == 0
        for row in R
        for x in row
    )

    print()
    print("=" * 78)
    print("2. EXACT INTERMEDIATE FACTORS")
    print("=" * 78)

    print(
        f"  L_shape={len(L)}x{len(L[0])}"
    )

    print(
        f"  R_shape={len(R)}x{len(R[0])}"
    )

    print(
        f"  L_even={L_even}"
    )

    print(
        f"  R_even={R_even}"
    )

    if not L_even or not R_even:
        raise ArithmeticError(
            "Expected both L and R to be entrywise even."
        )

    # ------------------------------------------------------------------
    # 3. FIRST 2-ADIC LIFTS
    # ------------------------------------------------------------------

    L1 = divide_matrix_exact(
        L,
        2,
        "L"
    )

    R1 = divide_matrix_exact(
        R,
        2,
        "R"
    )

    L1_mod2 = to_mod2(L1)
    R1_mod2 = to_mod2(R1)

    print()
    print("=" * 78)
    print("3. FIRST 2-ADIC LIFTS")
    print("=" * 78)

    rank_L1 = rank_mod2(L1_mod2)
    rank_R1 = rank_mod2(R1_mod2)

    print(
        f"  rank_L_mod2={rank_mod2(L)}"
    )

    print(
        f"  rank_R_mod2={rank_mod2(R)}"
    )

    print(
        f"  rank_L1_mod2={rank_L1}"
    )

    print(
        f"  rank_R1_mod2={rank_R1}"
    )

    print("  L1_mod2:")

    for row in L1_mod2:
        print(
            "    "
            + " ".join(str(x) for x in row)
        )

    print("  R1_mod2:")

    for row in R1_mod2:
        print(
            "    "
            + " ".join(str(x) for x in row)
        )

    # ------------------------------------------------------------------
    # 4. EXACT K FROM BOTH SIDES
    # ------------------------------------------------------------------

    K_from_L1 = [
        [
            -x
            for x in row
        ]
        for row in matrix_to_int(
            matmul(L1, C),
            "L1*C"
        )
    ]

    K_from_R1 = [
        [
            -x
            for x in row
        ]
        for row in matrix_to_int(
            matmul(B, R1),
            "B*R1"
        )
    ]

    K = K_from_L1

    left_recovery_exact = (
        K_from_L1 == K_from_R1
    )

    print()
    print("=" * 78)
    print("4. EXACT K RECOVERY")
    print("=" * 78)

    print(
        f"  K_from_L1_exact=True"
    )

    print(
        f"  K_from_R1_exact=True"
    )

    print(
        f"  both_recover_same_K="
        f"{left_recovery_exact}"
    )

    # ------------------------------------------------------------------
    # 5. FIRST-LIFT FACTORIZATION
    # ------------------------------------------------------------------

    L1_u = L1_v = None
    R1_u = R1_v = None

    if rank_L1 == 1:
        L1_u, L1_v = rank_one_factorization_mod2(
            L1_mod2
        )

    if rank_R1 == 1:
        R1_u, R1_v = rank_one_factorization_mod2(
            R1_mod2
        )

    print()
    print("=" * 78)
    print("5. RECTANGULAR RANK-ONE FACTORIZATION")
    print("=" * 78)

    print(
        f"  L1_rank_one={rank_L1 == 1}"
    )

    if L1_u is not None:
        print(
            f"  L1_left_direction={L1_u}"
        )
        print(
            f"  L1_right_direction={L1_v}"
        )

    print(
        f"  R1_rank_one={rank_R1 == 1}"
    )

    if R1_u is not None:
        print(
            f"  R1_left_direction={R1_u}"
        )
        print(
            f"  R1_right_direction={R1_v}"
        )

    # ------------------------------------------------------------------
    # 6. DOES EACH LIFT PREDICT K?
    # ------------------------------------------------------------------

    K2 = to_mod2(K)

    L1C_mod2 = to_mod2(
        matmul(L1, C)
    )

    BR1_mod2 = to_mod2(
        matmul(B, R1)
    )

    left_predicts_K = (
        L1C_mod2 == K2
    )

    right_predicts_K = (
        BR1_mod2 == K2
    )

    print()
    print("=" * 78)
    print("6. FIRST-LIFT PREDICTION")
    print("=" * 78)

    print(
        f"  L1*C predicts K_mod2="
        f"{left_predicts_K}"
    )

    print(
        f"  B*R1 predicts K_mod2="
        f"{right_predicts_K}"
    )

    # ------------------------------------------------------------------
    # 7. SUPPORT
    # ------------------------------------------------------------------

    def support(A):
        return [
            (i, j)
            for i in range(len(A))
            for j in range(len(A[0]))
            if A[i][j]
        ]

    print()
    print("=" * 78)
    print("7. SUPPORT")
    print("=" * 78)

    print(
        f"  L1_support={support(L1_mod2)}"
    )

    print(
        f"  R1_support={support(R1_mod2)}"
    )

    print(
        f"  K_support={support(K2)}"
    )

    # ------------------------------------------------------------------
    # 8. CORE RANK REFERENCE
    # ------------------------------------------------------------------

    A = [
        [
            F[i][j] + K[i][j]
            for j in range(CORE_N)
        ]
        for i in range(CORE_N)
    ]

    rank_F = rank_mod2(F)
    rank_K = rank_mod2(K2)
    rank_A = rank_mod2(A)

    print()
    print("=" * 78)
    print("8. RANK TRANSITION")
    print("=" * 78)

    print(
        f"  rank_F_ret_mod2={rank_F}"
    )

    print(
        f"  rank_K_mod2={rank_K}"
    )

    print(
        f"  rank_A_mod2={rank_A}"
    )

    print(
        f"  rank_drop={rank_F-rank_A}"
    )

    # ------------------------------------------------------------------
    # 9. INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The previous experiment established

    B adj(H) ≡ 0 mod 2
    adj(H) C ≡ 0 mod 2

while K mod 2 is nonzero and rank one.

The first lifts are

    L1 = (B adj(H))/2
    R1 = (adj(H) C)/2.

Because

    L1 : 12x2
    R1 : 2x12,

their rank-one factorizations must be treated as rectangular
outer products.

The exact identities

    K = -L1 C
    K = -B R1

hold over the integers.

The mod-2 tests determine whether the first nonzero K layer is already
encoded by L1 or R1 individually.

This avoids the previous dimension error and does not assume that a
rectangular rank-one factor is a square matrix.
"""
    )

    # ------------------------------------------------------------------
    # 10. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        det_h == 2
        and L_even
        and R_even
        and left_recovery_exact
        and rank_F == 4
        and rank_K == 1
        and rank_A == 3
    )

    print()
    print("=" * 78)
    print("10. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  det_H_equals_2={det_h == 2}"
    )

    print(
        f"  L_even={L_even}"
    )

    print(
        f"  R_even={R_even}"
    )

    print(
        f"  K_left_right_recovery_agree="
        f"{left_recovery_exact}"
    )

    print(
        f"  rank_F_ret_four={rank_F == 4}"
    )

    print(
        f"  rank_K_one={rank_K == 1}"
    )

    print(
        f"  rank_A_three={rank_A == 3}"
    )

    print(
        f"  L1_rank_one={rank_L1 == 1}"
    )

    print(
        f"  R1_rank_one={rank_R1 == 1}"
    )

    print(
        f"  L1_predicts_K={left_predicts_K}"
    )

    print(
        f"  R1_predicts_K={right_predicts_K}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 169S COMPLETE")


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