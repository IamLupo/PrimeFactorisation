#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 145 — EXACT p-ADIC LIFT-DEFECT / FIRST-CORRECTION AUDIT
==============================================================================

Experiment 144 established:

    R mod 2  = rank 1
    R mod 4  = rank 1
    R mod 8  = rank 1
    R mod 16 != rank 1

    R mod 3  = rank 1
    R mod 9  != rank 1

    R mod 7  = rank 2

Experiment 145 isolates the FIRST correction responsible for those failures.

For a rank-1 factorization modulo p^e,

    R = u v^T + p^e E

with E an exact integer matrix.

We then inspect

    E mod p.

This is the first p-adic correction beyond the low-rank leading term.

For p=2 we use the exact rank-1 factorization modulo 8.
For p=3 we use the exact rank-1 factorization modulo 3.

For p=7 we use the exact rank-2 factorization modulo 7:

    R = U V + 7 E.

The questions are:

    1. Is E integral?
    2. What is rank(E mod p)?
    3. Does the correction itself have low rank?
    4. Does the first correction explain the exact loss of lifting?
    5. Are the correction ranks different for 2,3,7?

No full Smith reduction.
No floating point.
No SymPy.
No extrapolation.
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

DELTA_11 = 9512681472


# ============================================================================
# HELPERS
# ============================================================================

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


def valuation(n: int, p: int) -> int:
    n = abs(int(n))

    if n == 0:
        return 10**9

    v = 0

    while n % p == 0:
        n //= p
        v += 1

    return v


def digits(n: int) -> int:
    return len(str(abs(n))) if n else 1


# ============================================================================
# EXACT BAREISS DETERMINANT
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
                            "Bareiss division failed."
                        )
                    value //= previous

                M[i][j] = value

        for i in range(k + 1, n):
            M[i][k] = 0

        previous = pivot

    return sign * M[-1][-1]


# ============================================================================
# BUILD FULL SYSTEM
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
# BUILD F-CORE
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
                    "Non-integral F-core entry."
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

            mdet = bareiss_det(
                maximal_minor(A, i, j)
            )

            if mdet % DELTA_11 != 0:
                raise ArithmeticError(
                    "Delta_11 divisibility failed."
                )

            row.append(
                mdet // DELTA_11
            )

        R.append(row)

    return R


# ============================================================================
# MODULAR RANK
# ============================================================================

def rank_mod(A, p):

    M = [
        [x % p for x in row]
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

        M[rank], M[pivot] = M[pivot], M[rank]

        inv = pow(
            M[rank][col],
            -1,
            p,
        )

        for j in range(col, cols):
            M[rank][j] = (
                M[rank][j] * inv
            ) % p

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
                ) % p

        rank += 1

        if rank == rows:
            break

    return rank


# ============================================================================
# RANK-1 FACTOR RECOVERY
# ============================================================================

def rank1_factorization(A, modulus):

    n = len(A)

    pivot = None
    pivot_position = None

    for i in range(n):

        for j in range(n):

            a = A[i][j] % modulus

            if a == 0:
                continue

            try:
                inv = pow(a, -1, modulus)
            except ValueError:
                continue

            pivot = a
            pivot_position = (i, j)
            break

        if pivot_position is not None:
            break

    if pivot_position is None:
        return None

    i0, j0 = pivot_position
    pivot_inv = pow(pivot, -1, modulus)

    u = [
        A[i][j0] % modulus
        for i in range(n)
    ]

    v = [
        (A[i0][j] * pivot_inv) % modulus
        for j in range(n)
    ]

    exact = True

    for i in range(n):

        for j in range(n):

            lhs = A[i][j] % modulus
            rhs = (
                u[i] * v[j]
            ) % modulus

            if lhs != rhs:
                exact = False
                break

        if not exact:
            break

    return {
        "pivot": pivot,
        "position": pivot_position,
        "u": u,
        "v": v,
        "exact": exact,
    }


# ============================================================================
# RANK-2 FACTOR RECOVERY OVER F_7
# ============================================================================

def rank2_factorization_mod_prime(A, prime):

    n = len(A)

    M = [
        [x % prime for x in row]
        for row in A
    ]

    def rank_local(B):

        B = [row[:] for row in B]

        rows = len(B)
        cols = len(B[0])
        rank = 0

        for col in range(cols):

            pivot = None

            for r in range(rank, rows):

                if B[r][col] != 0:
                    pivot = r
                    break

            if pivot is None:
                continue

            B[rank], B[pivot] = (
                B[pivot],
                B[rank],
            )

            inv = pow(
                B[rank][col],
                -1,
                prime,
            )

            for j in range(col, cols):
                B[rank][j] = (
                    B[rank][j] * inv
                ) % prime

            for r in range(rows):

                if r == rank:
                    continue

                f = B[r][col]

                if f == 0:
                    continue

                for j in range(col, cols):
                    B[r][j] = (
                        B[r][j]
                        - f * B[rank][j]
                    ) % prime

            rank += 1

        return rank

    basis_rows = []

    for i in range(n):

        candidate = basis_rows + [M[i]]

        if rank_local(candidate) > len(basis_rows):

            basis_rows.append(M[i])

            if len(basis_rows) == 2:
                break

    if len(basis_rows) != 2:
        return None

    # Find two pivot columns.
    pivot_cols = None

    for a in range(n):

        for b in range(a + 1, n):

            det = (
                basis_rows[0][a] * basis_rows[1][b]
                - basis_rows[0][b] * basis_rows[1][a]
            ) % prime

            if det != 0:

                pivot_cols = (a, b)
                inv_det = pow(
                    det,
                    -1,
                    prime,
                )
                break

        if pivot_cols is not None:
            break

    if pivot_cols is None:
        return None

    a, b = pivot_cols

    coefficients = []

    for i in range(n):

        x = M[i][a]
        y = M[i][b]

        u_a = basis_rows[0][a]
        u_b = basis_rows[0][b]
        v_a = basis_rows[1][a]
        v_b = basis_rows[1][b]

        alpha = (
            (x * v_b - y * v_a)
            * inv_det
        ) % prime

        beta = (
            (u_a * y - u_b * x)
            * inv_det
        ) % prime

        coefficients.append(
            (alpha, beta)
        )

    exact = True

    for i in range(n):

        for j in range(n):

            predicted = (
                coefficients[i][0] * basis_rows[0][j]
                + coefficients[i][1] * basis_rows[1][j]
            ) % prime

            if predicted != M[i][j]:
                exact = False
                break

        if not exact:
            break

    return {
        "basis_rows": basis_rows,
        "pivot_cols": pivot_cols,
        "coefficients": coefficients,
        "exact": exact,
    }


# ============================================================================
# DEFECT BUILDERS
# ============================================================================

def rank1_defect(A, factor_data, modulus, prime):

    u = factor_data["u"]
    v = factor_data["v"]

    n = len(A)

    E = []

    for i in range(n):

        row = []

        for j in range(n):

            value = (
                A[i][j]
                - u[i] * v[j]
            )

            if value % modulus != 0:
                raise ArithmeticError(
                    "Rank-1 factor does not lift at stated modulus."
                )

            row.append(
                value // modulus
            )

        E.append(row)

    E_mod_p = [
        [
            value % prime
            for value in row
        ]
        for row in E
    ]

    return E, E_mod_p


def rank2_defect(A, factor_data, modulus, prime):

    basis = factor_data["basis_rows"]
    coeffs = factor_data["coefficients"]

    n = len(A)

    E = []

    for i in range(n):

        row = []

        for j in range(n):

            approximation = (
                coeffs[i][0] * basis[0][j]
                + coeffs[i][1] * basis[1][j]
            )

            value = (
                A[i][j]
                - approximation
            )

            if value % modulus != 0:
                raise ArithmeticError(
                    "Rank-2 factor does not lift at stated modulus."
                )

            row.append(
                value // modulus
            )

        E.append(row)

    E_mod_p = [
        [
            value % prime
            for value in row
        ]
        for row in E
    ]

    return E, E_mod_p


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 145 — EXACT p-ADIC LIFT-DEFECT / "
        "FIRST-CORRECTION AUDIT"
    )
    print("=" * 78)

    M = build_full_system()
    core = build_core(M)
    R = normalized_cofactors(core)

    print()
    print("=" * 78)
    print("1. EXACT NORMALIZED COFACTOR DATA")
    print("=" * 78)

    print("  R_shape=(12,12)")
    print(f"  Delta_11={DELTA_11}")

    # ------------------------------------------------------------------
    # 2. 2-ADIC FIRST CORRECTION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. 2-ADIC FIRST-CORRECTION AUDIT")
    print("=" * 78)

    f2 = rank1_factorization(
        R,
        8,
    )

    correction2_ok = False

    if f2 is None or not f2["exact"]:
        print(
            "  ERROR: rank-1 factorization modulo 8 "
            "could not be reconstructed."
        )
    else:

        E2, E2_mod = rank1_defect(
            R,
            f2,
            8,
            2,
        )

        rank_E2 = rank_mod(
            E2_mod,
            2,
        )

        correction2_ok = True

        print(
            f"  base_modulus=8"
        )

        print(
            f"  factorization_exact=True"
        )

        print(
            f"  defect_rank_mod_2={rank_E2}"
        )

        print(
            f"  defect_matrix_mod_2="
        )

        for row in E2_mod:
            print(
                "    "
                + " ".join(
                    str(x)
                    for x in row
                )
            )

        print()
        print(
            "  interpretation:"
        )

        print(
            "    R = u*v^T + 8E"
        )

        print(
            f"    rank(E mod 2)={rank_E2}"
        )

    # ------------------------------------------------------------------
    # 3. 3-ADIC FIRST CORRECTION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. 3-ADIC FIRST-CORRECTION AUDIT")
    print("=" * 78)

    f3 = rank1_factorization(
        R,
        3,
    )

    correction3_ok = False

    if f3 is None or not f3["exact"]:
        print(
            "  ERROR: rank-1 factorization modulo 3 "
            "could not be reconstructed."
        )
    else:

        E3, E3_mod = rank1_defect(
            R,
            f3,
            3,
            3,
        )

        rank_E3 = rank_mod(
            E3_mod,
            3,
        )

        correction3_ok = True

        print(
            "  base_modulus=3"
        )

        print(
            "  factorization_exact=True"
        )

        print(
            f"  defect_rank_mod_3={rank_E3}"
        )

        print(
            "  defect_matrix_mod_3="
        )

        for row in E3_mod:
            print(
                "    "
                + " ".join(
                    str(x)
                    for x in row
                )
            )

        print()
        print(
            "  interpretation:"
        )

        print(
            "    R = u*v^T + 3E"
        )

        print(
            f"    rank(E mod 3)={rank_E3}"
        )

    # ------------------------------------------------------------------
    # 4. 7-ADIC FIRST CORRECTION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. 7-ADIC FIRST-CORRECTION AUDIT")
    print("=" * 78)

    f7 = rank2_factorization_mod_prime(
        R,
        7,
    )

    correction7_ok = False

    if f7 is None or not f7["exact"]:
        print(
            "  ERROR: rank-2 factorization modulo 7 "
            "could not be reconstructed."
        )
    else:

        # The factorization coefficients are modulo 7.
        # We lift them to the integer representatives 0,...,6.
        E7, E7_mod = rank2_defect(
            R,
            f7,
            7,
            7,
        )

        rank_E7 = rank_mod(
            E7_mod,
            7,
        )

        correction7_ok = True

        print(
            "  base_modulus=7"
        )

        print(
            "  rank2_factorization_exact=True"
        )

        print(
            f"  defect_rank_mod_7={rank_E7}"
        )

        print(
            "  defect_matrix_mod_7="
        )

        for row in E7_mod:
            print(
                "    "
                + " ".join(
                    str(x)
                    for x in row
                )
            )

        print()
        print(
            "  interpretation:"
        )

        print(
            "    R = U*V + 7E"
        )

        print(
            f"    rank(E mod 7)={rank_E7}"
        )

    # ------------------------------------------------------------------
    # 5. CORRECTION RANK SUMMARY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. FIRST-CORRECTION RANK SUMMARY")
    print("=" * 78)

    if correction2_ok:
        print(
            f"  p=2: leading rank=1, "
            f"first defect rank={rank_E2}"
        )

    if correction3_ok:
        print(
            f"  p=3: leading rank=1, "
            f"first defect rank={rank_E3}"
        )

    if correction7_ok:
        print(
            f"  p=7: leading rank=2, "
            f"first defect rank={rank_E7}"
        )

    # ------------------------------------------------------------------
    # 6. DEFECT SUPPORT SUMMARY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. FIRST-CORRECTION SUPPORT")
    print("=" * 78)

    if correction2_ok:

        support2 = [
            (i, j)
            for i in range(12)
            for j in range(12)
            if E2_mod[i][j] != 0
        ]

        print(
            f"  p=2: nonzero_entries="
            f"{len(support2)}"
        )

        print(
            f"  p=2: support_sample="
            f"{support2[:30]}"
        )

    if correction3_ok:

        support3 = [
            (i, j)
            for i in range(12)
            for j in range(12)
            if E3_mod[i][j] != 0
        ]

        print(
            f"  p=3: nonzero_entries="
            f"{len(support3)}"
        )

        print(
            f"  p=3: support_sample="
            f"{support3[:30]}"
        )

    if correction7_ok:

        support7 = [
            (i, j)
            for i in range(12)
            for j in range(12)
            if E7_mod[i][j] != 0
        ]

        print(
            f"  p=7: nonzero_entries="
            f"{len(support7)}"
        )

        print(
            f"  p=7: support_sample="
            f"{support7[:30]}"
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
Experiment 144 showed different lift depths:

    p=2:
        rank one through 8, failure at 16;

    p=3:
        rank one modulo 3, failure modulo 9;

    p=7:
        rank two modulo 7, failure modulo 49.

Experiment 145 isolates the first correction.

For p=2:

    R = u v^T + 8E_2.

For p=3:

    R = u v^T + 3E_3.

For p=7:

    R = U V + 7E_7.

The ranks of

    E_2 mod 2,
    E_3 mod 3,
    E_7 mod 7

measure the complexity of the first obstruction beyond the leading
finite-field low-rank structure.

A low defect rank would indicate a layered p-adic structure.

A high defect rank would indicate that the low-rank phenomenon is
quickly destroyed by the first correction.

The result is purely arithmetic and does not by itself identify
the origin of the cofactor matrix.
"""
    )

    # ------------------------------------------------------------------
    # 8. FINAL
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. FINAL EXACTNESS")
    print("=" * 78)

    final_ok = (
        correction2_ok
        and correction3_ok
        and correction7_ok
    )

    print(
        f"  p2_defect_reconstruction={correction2_ok}"
    )

    print(
        f"  p3_defect_reconstruction={correction3_ok}"
    )

    print(
        f"  p7_defect_reconstruction={correction7_ok}"
    )

    print(
        "  first_correction_audit_completed=True"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 145 COMPLETE")


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
