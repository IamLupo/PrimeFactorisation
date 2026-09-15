#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 139 — EXACT COMPLETE SMITH SPECTRUM OF THE F-CORE
==============================================================================

Purpose
-------

Experiment 138 established:

    full 14x14 system
        |
        | exact boundary elimination
        v
    12x12 F-core

and found

    d_14(full) = d_12(core).

Experiment 139 now determines the COMPLETE Smith invariant spectrum

    d_1 | d_2 | ... | d_12

of the 12x12 F-core.

The experiment is deliberately evidence-first.

It does NOT assume:

    * that the last invariant is a Cramer's artifact;
    * that the last invariant is fundamental;
    * that factorial/binomial normalization explains the coefficients;
    * that the system comes from the original (p,q)-kernel.

The exact Smith spectrum itself is the object under study.

The script uses integer row/column operations only.

No SymPy.
No floating point.
No extrapolation.
"""


from __future__ import annotations

from fractions import Fraction
from math import gcd
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


# ============================================================================
# BASIC HELPERS
# ============================================================================

def digits(n: int) -> int:
    return len(str(abs(n))) if n else 1


def valuation(n: int, p: int) -> int:
    n = abs(n)
    v = 0

    while n and n % p == 0:
        n //= p
        v += 1

    return v


def gcd_many(values) -> int:
    g = 0

    for x in values:
        g = gcd(g, abs(int(x)))

    return abs(g)


def row_metadata(index: int) -> tuple[int, int, int]:

    if 0 <= index < 6:
        return 1, 3, index

    if 6 <= index < 11:
        return 3, 5, index - 6

    if 11 <= index < 14:
        return 5, 7, index - 11

    raise IndexError(index)


def q_value(p: int, r: int) -> int:

    if r < 0:
        return 0

    if r >= len(Q[p]):
        return 0

    return Q[p][r]


def poly_degree(values: list[int]) -> int:

    for i in range(len(values) - 1, -1, -1):
        if values[i] != 0:
            return i

    return -1


# ============================================================================
# RANK
# ============================================================================

def matrix_rank(A) -> int:

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

        pv = M[rank][col]

        for j in range(col, cols):
            M[rank][j] /= pv

        for r in range(rows):

            if r == rank:
                continue

            factor = M[r][col]

            if factor == 0:
                continue

            for j in range(col, cols):
                M[r][j] -= factor * M[rank][j]

        rank += 1

        if rank == rows:
            break

    return rank


# ============================================================================
# EXACT BAREISS DETERMINANT
# ============================================================================

def bareiss_det(A: list[list[int]]) -> int:

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

                value = (
                    B[i][j] * pivot
                    - B[i][k] * B[k][j]
                )

                if k > 0:

                    if value % previous != 0:
                        raise ArithmeticError(
                            "Bareiss division not exact."
                        )

                    value //= previous

                B[i][j] = value

        for i in range(k + 1, n):
            B[i][k] = 0

        previous = pivot

    return sign * B[n - 1][n - 1]


# ============================================================================
# FULL SYSTEM
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


def build_full_system() -> list[list[int]]:

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
# EXACT F-CORE
# ============================================================================

def build_core(
    M: list[list[int]],
) -> tuple[list[list[int]], int]:

    # Boundary equations p=1,d=0 and p=3,d=0.
    row0 = M[5]
    row1 = M[10]

    H = [
        [row0[12], row0[13]],
        [row1[12], row1[13]],
    ]

    H_det = bareiss_det(H)

    if H_det == 0:
        raise ArithmeticError(
            "Boundary block singular."
        )

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

    _, p0_next, r0 = row_metadata(5)
    _, p1_next, r1 = row_metadata(10)

    rhs0 = q_value(p0_next, r0)
    rhs1 = q_value(p1_next, r1)

    h_const = [
        H_inv[i][0] * rhs0
        + H_inv[i][1] * rhs1
        for i in range(2)
    ]

    h_coeff = [
        [Fraction(0) for _ in range(12)]
        for _ in range(2)
    ]

    boundary_A = [
        row0[:12],
        row1[:12],
    ]

    for i in range(2):

        for j in range(12):

            h_coeff[i][j] = -(
                H_inv[i][0] * boundary_A[0][j]
                + H_inv[i][1] * boundary_A[1][j]
            )

    core_augmented = []

    for idx in range(14):

        if idx in (5, 10):
            continue

        row = M[idx]

        p, p_next, r = row_metadata(idx)

        target = Fraction(
            q_value(p_next, r)
        )

        reduced = []

        for j in range(12):

            value = (
                Fraction(row[j])
                + Fraction(row[12]) * h_coeff[0][j]
                + Fraction(row[13]) * h_coeff[1][j]
            )

            if value.denominator != 1:
                raise ArithmeticError(
                    "Core coefficient is non-integral."
                )

            reduced.append(value.numerator)

        reduced_rhs = (
            target
            - Fraction(row[12]) * h_const[0]
            - Fraction(row[13]) * h_const[1]
        )

        if reduced_rhs.denominator != 1:
            raise ArithmeticError(
                "Core RHS is non-integral."
            )

        reduced.append(reduced_rhs.numerator)

        core_augmented.append(reduced)

    return (
        [row[:12] for row in core_augmented],
        H_det,
    )


# ============================================================================
# INTEGER ROW/COLUMN OPERATIONS
# ============================================================================

def swap_rows(A, i, j):
    A[i], A[j] = A[j], A[i]


def swap_cols(A, i, j):

    for row in A:
        row[i], row[j] = row[j], row[i]


def row_add(A, target, source, q):

    if q == 0:
        return

    for j in range(len(A[0])):
        A[target][j] += q * A[source][j]


def col_add(A, target, source, q):

    if q == 0:
        return

    for i in range(len(A)):
        A[i][target] += q * A[i][source]


def row_negate(A, i):

    for j in range(len(A[0])):
        A[i][j] = -A[i][j]


def col_negate(A, j):

    for i in range(len(A)):
        A[i][j] = -A[i][j]


# ============================================================================
# PRACTICAL INTEGER SMITH REDUCTION
# ============================================================================

def smith_normal_form(A: list[list[int]]) -> list[int]:
    """
    Exact Smith normal form using integer row/column operations.

    The implementation uses repeated Euclidean reduction in the active
    submatrix. It explicitly checks the divisibility condition before
    advancing to the next diagonal position.

    Returns only the nonzero diagonal invariant factors.
    """

    M = [
        [int(x) for x in row]
        for row in A
    ]

    m = len(M)
    n = len(M[0])

    k = 0

    while k < m and k < n:

        # --------------------------------------------------------------
        # Find any nonzero entry in the active submatrix.
        # --------------------------------------------------------------

        pivot_pos = None

        for i in range(k, m):

            for j in range(k, n):

                if M[i][j] != 0:
                    pivot_pos = (i, j)
                    break

            if pivot_pos is not None:
                break

        if pivot_pos is None:
            break

        pi, pj = pivot_pos

        if pi != k:
            swap_rows(M, pi, k)

        if pj != k:
            swap_cols(M, pj, k)

        # --------------------------------------------------------------
        # Euclidean reduction.
        # --------------------------------------------------------------

        while True:

            changed = False

            # Reduce all entries in the pivot column.
            for i in range(k + 1, m):

                if M[i][k] == 0:
                    continue

                q = M[i][k] // M[k][k]

                row_add(M, i, k, -q)

                if abs(M[i][k]) < abs(M[k][k]) and M[i][k] != 0:
                    swap_rows(M, i, k)

                changed = True

            # Reduce all entries in the pivot row.
            for j in range(k + 1, n):

                if M[k][j] == 0:
                    continue

                q = M[k][j] // M[k][k]

                col_add(M, j, k, -q)

                if abs(M[k][j]) < abs(M[k][k]) and M[k][j] != 0:
                    swap_cols(M, j, k)

                changed = True

            # Continue until the pivot row/column vanish.
            active_nonzero = False

            for i in range(k + 1, m):

                if M[i][k] != 0:
                    active_nonzero = True
                    break

            if not active_nonzero:

                for j in range(k + 1, n):

                    if M[k][j] != 0:
                        active_nonzero = True
                        break

            if active_nonzero:
                continue

            # ----------------------------------------------------------
            # Pivot row and column are clean. Check divisibility of
            # every remaining active entry by the pivot.
            # ----------------------------------------------------------

            bad = None

            pivot = M[k][k]

            for i in range(k + 1, m):

                for j in range(k + 1, n):

                    if M[i][j] % pivot != 0:
                        bad = (i, j)
                        break

                if bad is not None:
                    break

            if bad is None:
                break

            # Introduce the offending entry into the pivot row.
            bi, bj = bad

            row_add(M, k, bi, 1)

            # The active pivot will now be reduced by the Euclidean loop.
            changed = True

        # --------------------------------------------------------------
        # Normalize sign.
        # --------------------------------------------------------------

        if M[k][k] < 0:
            row_negate(M, k)

        k += 1

    diag = []

    for i in range(min(m, n)):

        value = M[i][i]

        if value != 0:
            diag.append(abs(value))

    return diag


# ============================================================================
# VERIFY SMITH CHAIN
# ============================================================================

def smith_chain_ok(diag: list[int]) -> bool:

    if not diag:
        return True

    for i in range(len(diag) - 1):

        if diag[i + 1] % diag[i] != 0:
            return False

    return True


# ============================================================================
# PRODUCT CHECK
# ============================================================================

def product(values) -> int:

    out = 1

    for x in values:
        out *= x

    return out


# ============================================================================
# UNIMODULAR TEST TRANSFORMS
# ============================================================================

def transformed(A, name):

    B = [row[:] for row in A]

    if name == "row_add(0,1,+3)":
        row_add(B, 0, 1, 3)

    elif name == "row_add(5,2,-7)":
        row_add(B, 5, 2, -7)

    elif name == "col_add(3,7,+5)":
        col_add(B, 3, 7, 5)

    elif name == "col_add(10,4,-4)":
        col_add(B, 10, 4, -4)

    elif name == "row_swap(0,11)":
        swap_rows(B, 0, 11)

    elif name == "col_swap(1,10)":
        swap_cols(B, 1, 10)

    elif name == "row_negate(6)":
        row_negate(B, 6)

    else:
        raise ValueError(name)

    return B


# ============================================================================
# PRIME SPECTRUM
# ============================================================================

def mod_rank(A, p):

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

        inv = pow(M[rank][col], -1, p)

        for j in range(col, cols):
            M[rank][j] = (
                M[rank][j] * inv
            ) % p

        for r in range(rows):

            if r == rank:
                continue

            factor = M[r][col] % p

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
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 139 — EXACT COMPLETE SMITH SPECTRUM "
        "OF THE F-CORE"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. VALIDATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. EXACT DATA VALIDATION")
    print("=" * 78)

    data_exact = True

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

        data_exact = data_exact and exact

    print()
    print(f"  data_exact={data_exact}")

    if not data_exact:
        print("\nDATA VALIDATION FAILED")
        return

    # ------------------------------------------------------------------
    # 2. CORE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. EXACT 12x12 F-CORE")
    print("=" * 78)

    full = build_full_system()
    full_rank = matrix_rank(full)

    core, H_det = build_core(full)

    core_rank = matrix_rank(core)
    core_det = bareiss_det(core)

    print(f"  full_rank={full_rank}")
    print(f"  H_det={H_det}")
    print(f"  core_shape=({len(core)},{len(core[0])})")
    print(f"  core_rank={core_rank}")
    print(f"  core_det_digits={digits(abs(core_det))}")
    print(f"  core_det_sign={'+' if core_det > 0 else '-'}")

    # ------------------------------------------------------------------
    # 3. COMPLETE SMITH SPECTRUM
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. COMPLETE CORE SMITH SPECTRUM")
    print("=" * 78)

    smith = smith_normal_form(core)

    chain_ok = smith_chain_ok(smith)
    product_ok = product(smith) == abs(core_det)
    count_ok = len(smith) == 12

    print(f"  invariant_count={len(smith)}")
    print(f"  divisibility_chain={chain_ok}")
    print(f"  product_equals_abs_det={product_ok}")

    print()
    print("  invariant_factors=[")

    for i, value in enumerate(smith, start=1):

        print(
            f"    d_{i}={value}"
        )

    print("  ]")

    print()
    print("  digit_sizes=[")

    for i, value in enumerate(smith, start=1):

        print(
            f"    d_{i}: {digits(value)} digits"
        )

    print("  ]")

    # ------------------------------------------------------------------
    # 4. PRIME VALUATIONS OF EVERY SMITH FACTOR
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. COMPLETE SMITH PRIME-VALUATION TABLE")
    print("=" * 78)

    relevant_primes = [
        2, 3, 5, 7, 11, 13,
        17, 19, 23, 29, 31,
    ]

    for prime in relevant_primes:

        values = [
            valuation(d, prime)
            for d in smith
        ]

        if any(values):

            print()
            print(f"  prime={prime}")
            print(
                f"    valuations={values}"
            )
            print(
                f"    total={sum(values)}"
            )

    # ------------------------------------------------------------------
    # 5. SUPPORT OF THE SMITH SPECTRUM
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. SMITH SUPPORT PROFILE")
    print("=" * 78)

    nonunit_positions = []

    for i, value in enumerate(smith, start=1):

        if value != 1:
            nonunit_positions.append(i)

    print(
        f"  nonunit_positions={nonunit_positions}"
    )

    unit_count = sum(
        1
        for value in smith
        if value == 1
    )

    print(
        f"  unit_count={unit_count}"
    )

    # ------------------------------------------------------------------
    # 6. TERMINAL CONCENTRATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. TERMINAL CONCENTRATION AUDIT")
    print("=" * 78)

    prefix_products = []

    running = 1

    for value in smith:

        running *= value
        prefix_products.append(running)

    print(
        f"  product_first_11={prefix_products[10]}"
    )

    print(
        f"  last_factor=d_12={smith[-1]}"
    )

    ratio_last_to_prefix = Fraction(
        smith[-1],
        prefix_products[10],
    )

    print(
        f"  d12_over_product_first11="
        f"{ratio_last_to_prefix}"
    )

    # ------------------------------------------------------------------
    # 7. MODULAR RANK CROSS-CHECK
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. MODULAR RANK VS SMITH SUPPORT")
    print("=" * 78)

    test_primes = [
        2, 3, 5, 7, 11, 13,
        17, 19, 23, 29, 31,
        37, 41, 43, 47,
    ]

    for prime in test_primes:

        rank_mod = mod_rank(core, prime)

        zero_smith_count = sum(
            1
            for value in smith
            if value % prime == 0
        )

        predicted_rank = 12 - zero_smith_count

        exact = (
            rank_mod == predicted_rank
        )

        print(
            f"  prime={prime}: "
            f"rank_mod_p={rank_mod} "
            f"zero_smith_count={zero_smith_count} "
            f"predicted_rank={predicted_rank} "
            f"exact={exact}"
        )

    # ------------------------------------------------------------------
    # 8. ENTRY CONTENT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. CORE ENTRY-CONTENT CHECK")
    print("=" * 78)

    content = gcd_many(
        x
        for row in core
        for x in row
    )

    print(
        f"  gcd_all_entries={content}"
    )

    # ------------------------------------------------------------------
    # 9. UNIMODULAR SMITH INVARIANCE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. UNIMODULAR SMITH INVARIANCE")
    print("=" * 78)

    transform_names = [
        "row_add(0,1,+3)",
        "row_add(5,2,-7)",
        "col_add(3,7,+5)",
        "col_add(10,4,-4)",
        "row_swap(0,11)",
        "col_swap(1,10)",
        "row_negate(6)",
    ]

    invariance_ok = True

    for name in transform_names:

        T = transformed(core, name)

        T_smith = smith_normal_form(T)

        same_smith = T_smith == smith

        print()
        print(f"  {name}")
        print(
            f"    same_complete_smith={same_smith}"
        )

        invariance_ok = (
            invariance_ok
            and same_smith
        )

    # ------------------------------------------------------------------
    # 10. FULL VS CORE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. FULL VS CORE TERMINAL INVARIANT")
    print("=" * 78)

    full_det_abs = abs(bareiss_det(full))

    print(
        f"  full_det_digits={digits(full_det_abs)}"
    )

    print(
        f"  core_det_digits={digits(abs(core_det))}"
    )

    print(
        f"  terminal_core_factor={smith[-1]}"
    )

    print(
        f"  terminal_core_factor_digits="
        f"{digits(smith[-1])}"
    )

    # ------------------------------------------------------------------
    # 11. FINAL
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("11. FINAL EXACTNESS")
    print("=" * 78)

    final_ok = (
        data_exact
        and full_rank == 14
        and core_rank == 12
        and core_det != 0
        and count_ok
        and chain_ok
        and product_ok
        and content == 1
        and invariance_ok
    )

    print(
        f"  data_exact={data_exact}"
    )
    print(
        f"  full_rank_exact={full_rank == 14}"
    )
    print(
        f"  core_rank_exact={core_rank == 12}"
    )
    print(
        f"  complete_smith_recovered={count_ok}"
    )
    print(
        f"  smith_divisibility_chain={chain_ok}"
    )
    print(
        f"  smith_product_check={product_ok}"
    )
    print(
        f"  entry_content_is_one={content == 1}"
    )
    print(
        f"  unimodular_smith_invariance={invariance_ok}"
    )
    print(
        f"  failures={0 if final_ok else 1}"
    )
    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 139 COMPLETE")


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
