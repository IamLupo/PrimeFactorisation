#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 135 — EXACT TERMINAL SMITH FACTOR /
                 UNIMODULAR INVARIANCE AUDIT
==============================================================================

Purpose
-------

Experiment 134P showed that the corrected primitive-normalized 14x14 system
has

    rank(M) = 14,

    Delta_13 = gcd(all 13x13 minors),

    d_14 = |det(M)| / Delta_13,

with d_14 carrying roughly 100 decimal digits.

Experiment 135 tests whether this large terminal invariant is intrinsic to
the integer lattice.

A unimodular row or column operation preserves the Smith normal form:

    row_i <- row_i + c row_j
    col_i <- col_i + c col_j
    row swaps
    column swaps
    multiplication of a row/column by -1

Therefore the following quantities must remain unchanged:

    rank
    |det(M)|
    Delta_13
    d_14

The experiment applies several exact unimodular transformations and checks
those quantities independently.

No SymPy.
No floating point.
No extrapolation.
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd
import sys


# ============================================================================
# CORRECTED EXPERIMENT-129 PRIMITIVE ROWS
# ============================================================================

Q = {
    1: [
        -126258696,
        -11600744760,
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


def poly_degree(row: list[int]) -> int:
    for i in range(len(row) - 1, -1, -1):
        if row[i] != 0:
            return i
    return -1


def matrix_rank(A: list[list[int]]) -> int:
    M = [
        [Fraction(x) for x in row]
        for row in A
    ]

    if not M:
        return 0

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

        pivot_value = M[rank][col]

        for j in range(col, cols):
            M[rank][j] /= pivot_value

        for r in range(rows):
            if r == rank:
                continue

            if M[r][col] == 0:
                continue

            factor = M[r][col]

            for j in range(col, cols):
                M[r][j] -= factor * M[rank][j]

        rank += 1

        if rank == rows:
            break

    return rank


# ============================================================================
# BAREISS EXACT DETERMINANT
# ============================================================================

def bareiss_det(A: list[list[int]]) -> int:
    A = [row[:] for row in A]

    n = len(A)

    if n == 0:
        return 1

    if n == 1:
        return A[0][0]

    sign = 1
    previous = 1

    for k in range(n - 1):

        pivot_row = None

        for r in range(k, n):
            if A[r][k] != 0:
                pivot_row = r
                break

        if pivot_row is None:
            return 0

        if pivot_row != k:
            A[k], A[pivot_row] = A[pivot_row], A[k]
            sign *= -1

        pivot = A[k][k]

        for i in range(k + 1, n):
            for j in range(k + 1, n):

                numerator = (
                    A[i][j] * pivot
                    - A[i][k] * A[k][j]
                )

                if k > 0:
                    if numerator % previous != 0:
                        raise ArithmeticError(
                            "Bareiss division was not exact."
                        )
                    numerator //= previous

                A[i][j] = numerator

        for i in range(k + 1, n):
            A[i][k] = 0

        previous = pivot

    return sign * A[n - 1][n - 1]


# ============================================================================
# POLYNOMIAL BASIS
# ============================================================================

def basis(p: int, d: int) -> tuple[int, ...]:
    return (
        1,
        d,
        d * d,
        p,
        p * d,
        p * p,
    )


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
# q_p(r)
# ============================================================================

def q_value(p: int, r: int) -> int:

    if r < 0:
        return 0

    row = Q[p]

    if r >= len(row):
        return 0

    return row[r]


# ============================================================================
# 14x14 SYSTEM
# ============================================================================

def build_system() -> list[list[int]]:

    M = []

    for p, p2 in TRANSITIONS:

        max_r = D[p]

        for r in range(max_r + 1):

            d = max_r - r

            q0 = q_value(p, r)
            q1 = q_value(p, r + 1)

            b = basis(p, d)

            row = [0] * 14

            # F_0 block
            for j in range(6):
                row[j] += b[j] * q0

            # F_1 block
            for j in range(6):
                row[6 + j] += b[j] * q1

            # H_0 boundary block
            if d == 0:
                row[12] += 1
                row[13] += p

            M.append(row)

    return M


# ============================================================================
# MINOR
# ============================================================================

def delete_row_col(
    A: list[list[int]],
    row_index: int,
    col_index: int,
) -> list[list[int]]:

    return [
        [
            A[i][j]
            for j in range(len(A))
            if j != col_index
        ]
        for i in range(len(A))
        if i != row_index
    ]


# ============================================================================
# DELTA_13
# ============================================================================

def delta_13(A: list[list[int]]) -> int:

    g = 0
    n = len(A)

    for i in range(n):
        for j in range(n):

            minor = delete_row_col(A, i, j)

            value = bareiss_det(minor)

            g = gcd(g, abs(value))

            if g == 1:
                return 1

    return abs(g)


# ============================================================================
# UNIMODULAR TRANSFORMATIONS
# ============================================================================

def row_add(
    A: list[list[int]],
    target: int,
    source: int,
    multiplier: int,
) -> list[list[int]]:

    B = [row[:] for row in A]

    B[target] = [
        B[target][j] + multiplier * B[source][j]
        for j in range(len(B[target]))
    ]

    return B


def col_add(
    A: list[list[int]],
    target: int,
    source: int,
    multiplier: int,
) -> list[list[int]]:

    B = [row[:] for row in A]

    for i in range(len(B)):
        B[i][target] += multiplier * B[i][source]

    return B


def row_swap(
    A: list[list[int]],
    i: int,
    j: int,
) -> list[list[int]]:

    B = [row[:] for row in A]

    B[i], B[j] = B[j], B[i]

    return B


def col_swap(
    A: list[list[int]],
    i: int,
    j: int,
) -> list[list[int]]:

    B = [row[:] for row in A]

    for row in B:
        row[i], row[j] = row[j], row[i]

    return B


def row_negate(
    A: list[list[int]],
    i: int,
) -> list[list[int]]:

    B = [row[:] for row in A]

    B[i] = [-x for x in B[i]]

    return B


# ============================================================================
# INVARIANT AUDIT
# ============================================================================

def audit_matrix(
    name: str,
    A: list[list[int]],
    base_det_abs: int,
    base_delta13: int,
    base_d14: int,
) -> bool:

    rank = matrix_rank(A)

    det = bareiss_det(A)
    det_abs = abs(det)

    delta = delta_13(A)

    if delta == 0:
        d14 = None
    else:
        d14 = det_abs // delta

    rank_ok = rank == 14
    det_ok = det_abs == base_det_abs
    delta_ok = delta == base_delta13
    d14_ok = d14 == base_d14

    exact = (
        rank_ok
        and det_ok
        and delta_ok
        and d14_ok
    )

    print(f"\n  {name}")
    print(f"    rank={rank} exact={rank_ok}")
    print(
        f"    determinant_digits={digits(det_abs)} "
        f"same_abs_det={det_ok}"
    )
    print(
        f"    Delta_13={delta} "
        f"same_Delta_13={delta_ok}"
    )
    print(
        f"    d_14={d14} "
        f"same_d_14={d14_ok}"
    )
    print(f"    invariant_exact={exact}")

    return exact


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:

    print("=" * 78)
    print(
        "EXPERIMENT 135 — EXACT TERMINAL SMITH FACTOR / "
        "UNIMODULAR INVARIANCE AUDIT"
    )
    print("=" * 78)
    print()

    data_ok = validate_data()

    if not data_ok:
        print("\nDATA VALIDATION FAILED")
        return

    M = build_system()

    print()
    print("=" * 78)
    print("2. BASE SYSTEM")
    print("=" * 78)

    base_rank = matrix_rank(M)
    base_det = bareiss_det(M)
    base_det_abs = abs(base_det)

    base_delta13 = delta_13(M)

    if base_delta13 == 0:
        print("  Delta_13=0")
        return

    base_d14 = base_det_abs // base_delta13

    base_ok = (
        base_rank == 14
        and base_det != 0
        and base_det_abs % base_delta13 == 0
    )

    print(f"  rows={len(M)}")
    print(f"  columns={len(M[0])}")
    print(f"  rank={base_rank}")
    print(f"  determinant_digits={digits(base_det_abs)}")
    print(f"  determinant_sign={sign_of(base_det)}")
    print(f"  Delta_13={base_delta13}")
    print(f"  Delta_13_digits={digits(base_delta13)}")
    print(f"  d_14={base_d14}")
    print(f"  d_14_digits={digits(base_d14)}")
    print(f"  base_exact={base_ok}")

    print()
    print("=" * 78)
    print("3. BASE PRIME VALUATIONS")
    print("=" * 78)

    def valuations(n: int, primes: list[int]) -> dict[int, int]:

        result = {}

        for prime in primes:
            x = abs(n)
            v = 0

            while x and x % prime == 0:
                x //= prime
                v += 1

            if v:
                result[prime] = v

        return result

    small_primes = [2, 3, 5, 7, 11, 13]

    print(
        f"  det valuations="
        f"{valuations(base_det_abs, small_primes)}"
    )
    print(
        f"  Delta_13 valuations="
        f"{valuations(base_delta13, small_primes)}"
    )
    print(
        f"  d_14 valuations="
        f"{valuations(base_d14, small_primes)}"
    )

    print()
    print("=" * 78)
    print("4. UNIMODULAR INVARIANCE TESTS")
    print("=" * 78)

    transformed = []

    transformed.append(
        (
            "row_add(0,1,+3)",
            row_add(M, 0, 1, 3),
        )
    )

    transformed.append(
        (
            "row_add(5,2,-7)",
            row_add(M, 5, 2, -7),
        )
    )

    transformed.append(
        (
            "col_add(3,7,+5)",
            col_add(M, 3, 7, 5),
        )
    )

    transformed.append(
        (
            "col_add(10,4,-4)",
            col_add(M, 10, 4, -4),
        )
    )

    transformed.append(
        (
            "row_swap(0,13)",
            row_swap(M, 0, 13),
        )
    )

    transformed.append(
        (
            "col_swap(1,12)",
            col_swap(M, 1, 12),
        )
    )

    transformed.append(
        (
            "row_negate(6)",
            row_negate(M, 6),
        )
    )

    # Combined operations.
    M_combo = row_add(M, 0, 1, 2)
    M_combo = col_add(M_combo, 4, 9, -3)
    M_combo = row_swap(M_combo, 2, 11)
    transformed.append(
        (
            "combined_unimodular",
            M_combo,
        )
    )

    invariance_ok = True

    for name, A in transformed:

        exact = audit_matrix(
            name,
            A,
            base_det_abs,
            base_delta13,
            base_d14,
        )

        invariance_ok = invariance_ok and exact

    print()
    print("=" * 78)
    print("5. TERMINAL-FACTOR INTERPRETATION")
    print("=" * 78)

    print(
        f"  base determinant digits = {digits(base_det_abs)}"
    )
    print(
        f"  Delta_13 digits = {digits(base_delta13)}"
    )
    print(
        f"  terminal d_14 digits = {digits(base_d14)}"
    )
    print(
        f"  determinant / d_14 = "
        f"{base_det_abs // base_d14}"
    )

    print()
    print(
        "  The determinant's magnitude is decomposed as"
    )
    print(
        "      |det(M)| = Delta_13 * d_14."
    )
    print()
    print(
        "  The unimodular tests ask whether Delta_13 and d_14 survive"
    )
    print(
        "  exact integer coordinate changes. If they do, their sizes"
    )
    print(
        "  are lattice invariants rather than artifacts of the chosen"
    )
    print(
        "  basis presentation."
    )

    print()
    print("=" * 78)
    print("6. FINAL EXACTNESS")
    print("=" * 78)

    final_ok = (
        data_ok
        and base_ok
        and invariance_ok
    )

    print(f"  data_exact={data_ok}")
    print(f"  base_system_exact={base_ok}")
    print(f"  unimodular_invariance={invariance_ok}")
    print(f"  terminal_smith_factor_invariant={invariance_ok}")
    print(f"  failures={0 if final_ok else 1}")
    print(f"  ALL BASIC CHECKS PASS={final_ok}")
    print()
    print("EXPERIMENT 135 COMPLETE")


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

