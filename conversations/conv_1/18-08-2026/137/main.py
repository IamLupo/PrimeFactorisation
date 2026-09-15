#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 137 — EXACT TERMINAL-SMITH / BOUNDARY-SCHUR /
                 INTERIOR-CORE AUDIT — HARDENED
==============================================================================

Purpose
-------

Experiment 135 established that the 100-digit terminal Smith invariant
of the corrected 14x14 B-odd system is invariant under unimodular
row/column transformations.

Experiment 136 then attempted to separate the full system into:

    * boundary equations,
    * a boundary-eliminated F-core,
    * an interior-only subsystem.

This hardened version fixes the missing matrix_rank() helper and makes
all matrix operations explicit and exact.

The experiment compares:

    full 14x14 determinant,
    Delta_13,
    terminal Smith factor d_14,
    boundary elimination determinant,
    reduced 12x12 F-core determinant,
    interior 11x12 maximal-minor gcd.

No SymPy.
No floating point.
No extrapolation.
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd
from itertools import combinations
import sys


# ============================================================================
# DATA — CORRECTED PRIMITIVE B-ODD ROWS
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

# NOTE:
# The p=1 row above intentionally follows the corrected primitive-normalized
# dataset used by the preceding experiments.  If your local Experiment-134P
# data differs, replace ONLY this data table.

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


# ============================================================================
# BASIC HELPERS
# ============================================================================

def digits(n: int) -> int:
    return len(str(abs(n))) if n else 1


def sign_of(n: int) -> str:
    if n > 0:
        return "+"
    if n < 0:
        return "-"
    return "0"


def poly_degree(values: list[int]) -> int:
    for i in range(len(values) - 1, -1, -1):
        if values[i] != 0:
            return i
    return -1


def gcd_many(values) -> int:
    g = 0
    for x in values:
        g = gcd(g, abs(int(x)))
    return abs(g)


def lcm(a: int, b: int) -> int:
    if a == 0 or b == 0:
        return 0
    return abs(a // gcd(a, b) * b)


# ============================================================================
# EXACT MATRIX RANK
# ============================================================================

def matrix_rank(A) -> int:
    """
    Exact Gaussian-elimination rank over QQ.

    This is deliberately independent of SymPy/numpy so that all rank
    decisions are exact Fraction decisions.
    """
    if not A:
        return 0

    M = [
        [Fraction(x) for x in row]
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

        if pivot != rank:
            M[rank], M[pivot] = M[pivot], M[rank]

        pivot_value = M[rank][col]

        for j in range(col, cols):
            M[rank][j] /= pivot_value

        for r in range(rows):

            if r == rank:
                continue

            value = M[r][col]

            if value == 0:
                continue

            for j in range(col, cols):
                M[r][j] -= value * M[rank][j]

        rank += 1

        if rank == rows:
            break

    return rank


# ============================================================================
# EXACT BAREISS DETERMINANT
# ============================================================================

def bareiss_det(A) -> int:
    """
    Fraction-free Bareiss determinant.

    For the integer matrices used here this stays in exact integer
    arithmetic and avoids enormous Fraction intermediate growth.
    """
    if not A:
        return 1

    B = [
        [int(x) for x in row]
        for row in A
    ]

    n = len(B)

    if n == 1:
        return B[0][0]

    previous = 1
    sign = 1

    for k in range(n - 1):

        pivot_row = None

        for r in range(k, n):
            if B[r][k] != 0:
                pivot_row = r
                break

        if pivot_row is None:
            return 0

        if pivot_row != k:
            B[k], B[pivot_row] = B[pivot_row], B[k]
            sign *= -1

        pivot = B[k][k]

        for i in range(k + 1, n):
            for j in range(k + 1, n):

                numerator = (
                    B[i][j] * pivot
                    - B[i][k] * B[k][j]
                )

                if k > 0:
                    if numerator % previous != 0:
                        raise ArithmeticError(
                            "Bareiss non-exact division."
                        )

                    numerator //= previous

                B[i][j] = numerator

        for i in range(k + 1, n):
            B[i][k] = 0

        previous = pivot

    return sign * B[n - 1][n - 1]


# ============================================================================
# DATA VALIDATION
# ============================================================================

def validate_data() -> bool:

    print("=" * 78)
    print("1. EXACT DATA VALIDATION")
    print("=" * 78)

    ok = True

    for p in (1, 3, 5, 7):

        entries = len(Q[p])
        degree = poly_degree(Q[p])
        expected = D[p]

        exact = (
            entries == expected + 1
            and degree == expected
        )

        print(
            f"  p={p}: entries={entries} "
            f"degree={degree} expected={expected} "
            f"exact={exact}"
        )

        ok = ok and exact

    print()
    print(f"  data_exact={ok}")

    return ok


# ============================================================================
# q VALUE
# ============================================================================

def q_value(p: int, r: int) -> int:

    if r < 0:
        return 0

    row = Q[p]

    if r >= len(row):
        return 0

    return row[r]


# ============================================================================
# MONOMIAL BASIS
# ============================================================================

def monomial_basis(p: int, d: int) -> list[int]:
    return [
        1,
        d,
        d * d,
        p,
        p * d,
        p * p,
    ]


# ============================================================================
# BUILD FULL 14 x 14 SYSTEM
# ============================================================================

def build_full_system() -> list[list[int]]:

    M = []

    for p, p_next in TRANSITIONS:

        max_r = D[p]

        for r in range(max_r + 1):

            d = max_r - r

            q0 = q_value(p, r)
            q1 = q_value(p, r + 1)

            basis = monomial_basis(p, d)

            row = [0] * 14

            # F_0
            for j in range(6):
                row[j] = basis[j] * q0

            # F_1
            for j in range(6):
                row[6 + j] = basis[j] * q1

            # H_0(p) = h0 + h1*p at d=0
            if d == 0:
                row[12] = 1
                row[13] = p

            M.append(row)

    return M


# ============================================================================
# ROW INDEXING
# ============================================================================

def row_metadata(index: int) -> tuple[int, int, int]:

    # p=1 block: rows 0..5
    if 0 <= index < 6:
        p = 1
        r = index
        p_next = 3

    # p=3 block: rows 6..10
    elif 6 <= index < 11:
        p = 3
        r = index - 6
        p_next = 5

    # p=5 block: rows 11..13
    elif 11 <= index < 14:
        p = 5
        r = index - 11
        p_next = 7

    else:
        raise IndexError(index)

    return p, p_next, r


def boundary_row_indices() -> list[int]:
    return [5, 10, 13]


def interior_row_indices() -> list[int]:
    boundary = set(boundary_row_indices())
    return [
        i for i in range(14)
        if i not in boundary
    ]


# ============================================================================
# MATRIX EXTRACTION
# ============================================================================

def extract_rows(
    A: list[list[int]],
    indices: list[int],
) -> list[list[int]]:
    return [A[i][:] for i in indices]


def extract_cols(
    A: list[list[int]],
    indices: list[int],
) -> list[list[int]]:
    return [
        [row[j] for j in indices]
        for row in A
    ]


def delete_row_col(
    A: list[list[int]],
    row_index: int,
    col_index: int,
) -> list[list[int]]:

    return [
        [
            A[i][j]
            for j in range(len(A[0]))
            if j != col_index
        ]
        for i in range(len(A))
        if i != row_index
    ]


# ============================================================================
# FULL DELTA_13
# ============================================================================

def delta_13(A: list[list[int]]) -> int:

    n = len(A)

    if n != 14 or len(A[0]) != 14:
        raise ValueError("Expected 14x14 matrix.")

    g = 0

    for i in range(n):
        for j in range(n):

            value = bareiss_det(
                delete_row_col(A, i, j)
            )

            g = gcd(g, abs(value))

            if g == 1:
                return 1

    return abs(g)


# ============================================================================
# BUILD 12 x 12 BOUNDARY-ELIMINATED CORE
# ============================================================================

def build_boundary_eliminated_core(
    M: list[list[int]],
) -> tuple[
    list[list[Fraction]],
    int,
    list[Fraction],
]:
    """
    Eliminate h0,h1 using the first two terminal equations:

        p=1, d=0
        p=3, d=0

    The remaining 12 equations give an exact 12x12 F-core.
    """

    bidx = boundary_row_indices()

    # First two boundary rows:
    b0 = bidx[0]
    b1 = bidx[1]

    row0 = M[b0]
    row1 = M[b1]

    H = [
        [row0[12], row0[13]],
        [row1[12], row1[13]],
    ]

    H_det = bareiss_det(H)

    if H_det == 0:
        raise ValueError(
            "Boundary H block is singular."
        )

    A = [
        row0[:12],
        row1[:12],
    ]

    # Boundary RHS values:
    p0, p0_next, r0 = row_metadata(b0)
    p1, p1_next, r1 = row_metadata(b1)

    rhs = [
        q_value(p0_next, r0),
        q_value(p1_next, r1),
    ]

    # Exact inverse of 2x2 H.
    H_inv = [
        [
            Fraction(H[1][1], H_det),
            Fraction(-H[0][1], H_det),
        ],
        [
            Fraction(-H[1][0], H_det),
            Fraction(H[0][0], H_det),
        ],
    ]

    # h = H^{-1}(rhs - A f)
    h_const = [
        sum(
            H_inv[i][j] * rhs[j]
            for j in range(2)
        )
        for i in range(2)
    ]

    h_coeff = [
        [Fraction(0) for _ in range(12)]
        for _ in range(2)
    ]

    for i in range(2):
        for j in range(12):
            h_coeff[i][j] = -sum(
                H_inv[i][t] * A[t][j]
                for t in range(2)
            )

    remaining = [
        i for i in range(14)
        if i not in (b0, b1)
    ]

    core = []

    for idx in remaining:

        row = M[idx]

        p, p_next, r = row_metadata(idx)

        rhs_value = q_value(p_next, r)

        reduced = []

        for j in range(12):

            value = (
                Fraction(row[j])
                + Fraction(row[12]) * h_coeff[0][j]
                + Fraction(row[13]) * h_coeff[1][j]
            )

            reduced.append(value)

        reduced.append(
            Fraction(rhs_value)
            - Fraction(row[12]) * h_const[0]
            - Fraction(row[13]) * h_const[1]
        )

        core.append(reduced)

    return core, H_det, h_const


# ============================================================================
# CLEAR FRACTION DENOMINATORS
# ============================================================================

def clear_denominators(
    A: list[list[Fraction]],
) -> tuple[int, list[list[int]]]:

    common = 1

    for row in A:
        for x in row:
            common = lcm(common, x.denominator)

    B = [
        [
            x.numerator * (common // x.denominator)
            for x in row
        ]
        for row in A
    ]

    return common, B


# ============================================================================
# INTERIOR MATRIX
# ============================================================================

def build_interior_F_matrix(
    M: list[list[int]],
) -> list[list[int]]:

    rows = extract_rows(
        M,
        interior_row_indices(),
    )

    return extract_cols(
        rows,
        list(range(12)),
    )


# ============================================================================
# INTERIOR MAXIMAL MINOR GCD
# ============================================================================

def interior_maximal_minor_gcd(
    A: list[list[int]],
) -> int:

    rows = len(A)
    cols = len(A[0])

    if rows != 11 or cols != 12:
        raise ValueError(
            "Expected 11x12 interior matrix."
        )

    g = 0

    # There are exactly 12 maximal minors:
    # omit one of the 12 columns.
    for omitted in range(cols):

        selected = [
            j for j in range(cols)
            if j != omitted
        ]

        sub = [
            [A[i][j] for j in selected]
            for i in range(rows)
        ]

        value = bareiss_det(sub)

        g = gcd(g, abs(value))

        if g == 1:
            return 1

    return abs(g)


# ============================================================================
# PRIME VALUATION
# ============================================================================

def valuation(n: int, p: int) -> int:

    n = abs(n)
    v = 0

    while n and n % p == 0:
        n //= p
        v += 1

    return v


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:

    print("=" * 78)
    print(
        "EXPERIMENT 137 — EXACT TERMINAL-SMITH / "
        "BOUNDARY-SCHUR / INTERIOR-CORE AUDIT"
    )
    print("=" * 78)

    # ----------------------------------------------------------------------
    # 1
    # ----------------------------------------------------------------------

    print()
    data_ok = validate_data()

    if not data_ok:
        print("\nDATA VALIDATION FAILED")
        return

    # ----------------------------------------------------------------------
    # 2
    # ----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. FULL 14x14 SYSTEM")
    print("=" * 78)

    M = build_full_system()

    full_rank = matrix_rank(M)
    det_full = bareiss_det(M)
    det_abs = abs(det_full)

    print(f"  rows={len(M)}")
    print(f"  columns={len(M[0])}")
    print(f"  rank={full_rank}")
    print(f"  determinant_digits={digits(det_abs)}")
    print(f"  determinant_sign={sign_of(det_full)}")

    delta13 = delta_13(M)
    d14 = det_abs // delta13

    print(f"  Delta_13={delta13}")
    print(f"  Delta_13_digits={digits(delta13)}")
    print(f"  d_14={d14}")
    print(f"  d_14_digits={digits(d14)}")

    full_exact = (
        full_rank == 14
        and det_full != 0
        and delta13 != 0
        and det_abs % delta13 == 0
    )

    print(f"  full_system_exact={full_exact}")

    # ----------------------------------------------------------------------
    # 3
    # ----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. BOUNDARY / INTERIOR PARTITION")
    print("=" * 78)

    b_rows = boundary_row_indices()
    i_rows = interior_row_indices()

    print(f"  boundary_rows={b_rows}")
    print(f"  interior_rows={i_rows}")
    print(f"  boundary_count={len(b_rows)}")
    print(f"  interior_count={len(i_rows)}")

    interior = build_interior_F_matrix(M)

    print(
        f"  interior_F_shape="
        f"({len(interior)},{len(interior[0])})"
    )
    print(
        f"  interior_F_rank="
        f"{matrix_rank(interior)}"
    )

    # ----------------------------------------------------------------------
    # 4
    # ----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. EXACT BOUNDARY-SCHUR ELIMINATION")
    print("=" * 78)

    core, H_det, h_const = (
        build_boundary_eliminated_core(M)
    )

    core_A = [
        row[:12]
        for row in core
    ]

    core_b = [
        row[12]
        for row in core
    ]

    core_rank = matrix_rank(core_A)

    common_den, core_int = clear_denominators(core_A)

    det_core_int = bareiss_det(core_int)

    det_core = Fraction(
        det_core_int,
        common_den ** 12,
    )

    print(f"  H_boundary_det={H_det}")
    print(f"  H_boundary_abs={abs(H_det)}")
    print(f"  core_rows={len(core_A)}")
    print(f"  core_columns={len(core_A[0])}")
    print(f"  core_rank={core_rank}")
    print(f"  core_common_denominator={common_den}")
    print(
        f"  core_integer_det_digits="
        f"{digits(abs(det_core_int))}"
    )
    print(f"  core_det={det_core}")

    core_exact = (
        core_rank == 12
        and det_core != 0
    )

    print(f"  boundary_eliminated_core_exact={core_exact}")

    # ----------------------------------------------------------------------
    # 5
    # ----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. EXACT DETERMINANT FACTORIZATION")
    print("=" * 78)

    schur_product = abs(Fraction(H_det) * det_core)

    determinant_identity = (
        schur_product == det_abs
    )

    print(f"  |det(full)|={det_abs}")
    print(
        f"  |H_det * det(core)|="
        f"{schur_product}"
    )
    print(
        f"  exact_schur_identity="
        f"{determinant_identity}"
    )

    # ----------------------------------------------------------------------
    # 6
    # ----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. INTERIOR-ONLY MAXIMAL MINOR AUDIT")
    print("=" * 78)

    interior_rank = matrix_rank(interior)

    interior_gcd = interior_maximal_minor_gcd(
        interior
    )

    print(f"  interior_rank={interior_rank}")
    print(
        f"  interior_maximal_minor_gcd="
        f"{interior_gcd}"
    )
    print(
        f"  interior_gcd_digits="
        f"{digits(interior_gcd)}"
    )

    interior_exact = (
        interior_rank == 11
        and interior_gcd != 0
    )

    print(
        f"  interior_minor_audit_exact="
        f"{interior_exact}"
    )

    # ----------------------------------------------------------------------
    # 7
    # ----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. TERMINAL-SMITH / BOUNDARY COMPARISON")
    print("=" * 78)

    print(f"  Delta_13={delta13}")
    print(f"  d_14_digits={digits(d14)}")
    print(f"  H_boundary_det={H_det}")
    print(f"  core_det={det_core}")
    print(
        f"  interior_maximal_minor_gcd="
        f"{interior_gcd}"
    )

    if interior_gcd != 0 and delta13 % interior_gcd == 0:
        print(
            f"  Delta13/interior_gcd="
            f"{delta13 // interior_gcd}"
        )
    else:
        print(
            "  Delta13/interior_gcd=nonintegral_or_undefined"
        )

    if H_det != 0:
        if d14 % abs(H_det) == 0:
            print(
                f"  d14/|H_det|="
                f"{d14 // abs(H_det)}"
            )
        else:
            print(
                "  d14/|H_det|=nonintegral"
            )

    # ----------------------------------------------------------------------
    # 8
    # ----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. PRIME-VALUATION COMPARISON")
    print("=" * 78)

    for prime in (2, 3, 5, 7, 11, 13):

        vd = valuation(d14, prime)
        vh = valuation(abs(H_det), prime)
        vi = valuation(interior_gcd, prime)

        print(
            f"  prime={prime}: "
            f"v(d14)={vd}, "
            f"v(H_det)={vh}, "
            f"v(interior_gcd)={vi}"
        )

    # ----------------------------------------------------------------------
    # 9
    # ----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The full lattice contains:

    14 equations
    14 unknowns

with twelve F-coordinates and two boundary-source coordinates H_0.

The three terminal equations occur at d=0.

Two of those boundary equations are sufficient to eliminate H_0 exactly,
producing a 12x12 rational F-core.

The determinant identity

    |det(M)| = |det(H_boundary)| * |det(core)|

is tested exactly.

Independently, removing all three d=0 equations gives an 11x12
interior F-matrix. Its maximal-minor gcd measures the integer lattice
already present before terminal boundary equations are imposed.

The comparison among

    d_14,
    det(core),
    det(H_boundary),
    interior maximal-minor gcd

is intended to locate the source of the large terminal arithmetic
obstruction.

No claim is made here about the original (p,q)-kernel.
"""
    )

    # ----------------------------------------------------------------------
    # 10
    # ----------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. FINAL EXACTNESS")
    print("=" * 78)

    final_ok = (
        data_ok
        and full_exact
        and core_exact
        and determinant_identity
        and interior_exact
    )

    print(f"  data_exact={data_ok}")
    print(f"  full_system_exact={full_exact}")
    print(f"  boundary_schur_exact={core_exact}")
    print(
        f"  determinant_factorization_exact="
        f"{determinant_identity}"
    )
    print(f"  interior_audit_exact={interior_exact}")
    print(f"  failures={0 if final_ok else 1}")
    print(f"  ALL BASIC CHECKS PASS={final_ok}")
    print()
    print("EXPERIMENT 137 COMPLETE")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nInterrupted.")
        sys.exit(130)
    except Exception as exc:
        print(
            f"\nFATAL ERROR: {type(exc).__name__}: {exc}"
        )
        raise