#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 150 — EXACT 2-ADIC DEFECT / TRANSITION-BLOCK PROVENANCE AUDIT
==============================================================================

Experiment 149 established that the first 2-adic obstruction

    E mod 2

has rank 3, and that deleting rows 0 and 5 reduces the rank to 1.

The 14 equations naturally originate from three transitions:

    p=1 -> 3 : rows  0..5
    p=3 -> 5 : rows  6..10
    p=5 -> 7 : rows 11..13

Experiment 150 asks whether the rank-3 obstruction is aligned with
those transition blocks.

We compute:

    * defect rank of each transition block;
    * defect support by transition;
    * exact row-type equivalence within each transition;
    * cross-transition row-type equivalence;
    * rank of cumulative transition blocks;
    * whether the complete p=5 -> 7 block vanishes;
    * whether the p=3 -> 5 contribution is contained in the same
      row direction already present in the p=1 -> 3 block.

The point is provenance inside the exact 14-equation system.

This is NOT a search for a recurrence and does not assume any hidden
mechanism.

No floating point.
No SymPy.
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


# ============================================================================
# BASIC HELPERS
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
# BUILD FULL SYSTEM
# ============================================================================

def build_full_system():

    M = []
    row_meta = []

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

            row_meta.append(
                {
                    "p": p,
                    "p_next": p_next,
                    "r": r,
                    "d": d,
                }
            )

    return M, row_meta


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
        raise ArithmeticError("Singular boundary block.")

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
# IMPORTANT ROW MAPPING
# ============================================================================

def core_row_meta():
    """
    Rows 5 and 10 of the 14x14 system are the two eliminated boundary rows.
    The corresponding rows do not appear in the 12x12 core.

    We therefore construct the exact metadata for the retained core rows.
    """

    _, meta = build_full_system()

    result = []

    for idx, info in enumerate(meta):

        if idx in (5, 10):
            continue

        result.append(info)

    return result


# ============================================================================
# NORMALIZED COFACTORS
# ============================================================================

def maximal_minor(A, omit_row, omit_col):

    rows = [i for i in range(len(A)) if i != omit_row]
    cols = [j for j in range(len(A[0])) if j != omit_col]

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
                    "Delta_11 divisibility failure."
                )

            row.append(
                mdet // DELTA_11
            )

        R.append(row)

    return R


# ============================================================================
# MOD 2 LINEAR ALGEBRA
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
            if M[r][col] == 1:
                pivot = r
                break

        if pivot is None:
            continue

        M[pivot_row], M[pivot] = M[pivot], M[pivot_row]

        for r in range(rows):

            if r == pivot_row:
                continue

            if M[r][col] == 1:

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


def row_space_basis_mod2(A):

    RREF, pivots = rref_mod2(A)

    basis = [
        row
        for row in RREF
        if any(row)
    ]

    return basis, pivots


# ============================================================================
# ROW-TYPE ANALYSIS
# ============================================================================

def row_tuple(row):
    return tuple(int(x) & 1 for x in row)


def unique_row_types(A):

    types = {}
    for i, row in enumerate(A):

        key = row_tuple(row)

        types.setdefault(
            key,
            [],
        ).append(i)

    return types


def row_type_labels(A):

    types = unique_row_types(A)

    ordered = []

    for key, indices in types.items():

        ordered.append(
            {
                "row": list(key),
                "indices": indices,
            }
        )

    ordered.sort(
        key=lambda x: x["indices"][0]
    )

    return ordered


# ============================================================================
# MATRIX UTILITIES
# ============================================================================

def submatrix(A, rows):

    return [
        A[i][:]
        for i in rows
    ]


def xor_rows(a, b):
    return [
        (x ^ y)
        for x, y in zip(a, b)
    ]


def is_zero_row(row):
    return not any(row)


def support_rows(A):

    return [
        i
        for i, row in enumerate(A)
        if any(row)
    ]


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 150 — EXACT 2-ADIC DEFECT / "
        "TRANSITION-BLOCK PROVENANCE AUDIT"
    )
    print("=" * 78)

    M, full_meta = build_full_system()
    core = build_core(M)
    R = normalized_cofactors(core)

    meta = core_row_meta()

    # Recover mod-8 rank-one structure.
    modulus = 8
    pivot_position = None
    pivot_inverse = None

    for i in range(12):

        for j in range(12):

            x = R[i][j] % modulus

            if x == 0:
                continue

            try:
                inv = pow(
                    x,
                    -1,
                    modulus,
                )
            except ValueError:
                continue

            pivot_position = (i, j)
            pivot_inverse = inv
            break

        if pivot_position is not None:
            break

    if pivot_position is None:
        raise ArithmeticError(
            "No unit pivot for mod-8 rank-one factorization."
        )

    i0, j0 = pivot_position

    u = [
        R[i][j0] % 8
        for i in range(12)
    ]

    v = [
        R[i0][j] * pivot_inverse % 8
        for j in range(12)
    ]

    E = []

    for i in range(12):

        row = []

        for j in range(12):

            diff = (
                R[i][j]
                - u[i] * v[j]
            )

            if diff % 8 != 0:
                raise ArithmeticError(
                    "Invalid mod-8 factorization."
                )

            row.append(
                (diff // 8) & 1
            )

        E.append(row)

    # ------------------------------------------------------------------
    # 1. EXACT DATA
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. EXACT TRANSITION METADATA")
    print("=" * 78)

    for i, info in enumerate(meta):

        print(
            f"  core_row={i}: "
            f"p={info['p']} -> {info['p_next']} "
            f"r={info['r']} d={info['d']}"
        )

    # ------------------------------------------------------------------
    # 2. TRANSITION ROW GROUPS
    # ------------------------------------------------------------------

    groups = {
        "p=1 -> 3": [
            i
            for i, info in enumerate(meta)
            if info["p"] == 1
        ],
        "p=3 -> 5": [
            i
            for i, info in enumerate(meta)
            if info["p"] == 3
        ],
        "p=5 -> 7": [
            i
            for i, info in enumerate(meta)
            if info["p"] == 5
        ],
    }

    print()
    print("=" * 78)
    print("2. TRANSITION-BLOCK SUPPORT")
    print("=" * 78)

    for name, rows in groups.items():

        block = submatrix(
            E,
            rows,
        )

        print(
            f"  {name}"
        )

        print(
            f"    rows={rows}"
        )

        print(
            f"    rank_mod_2={rank_mod2(block)}"
        )

        print(
            f"    nonzero_rows_relative="
            f"{support_rows(block)}"
        )

        print(
            f"    row_types="
        )

        for item in row_type_labels(block):

            print(
                f"      indices={item['indices']} "
                f"row={item['row']}"
            )

    # ------------------------------------------------------------------
    # 3. CROSS-TRANSITION ROW TYPES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. CROSS-TRANSITION ROW-TYPE EQUALITY")
    print("=" * 78)

    global_types = unique_row_types(E)

    for key, indices in sorted(
        global_types.items(),
        key=lambda pair: pair[1][0],
    ):

        annotated = []

        for i in indices:

            info = meta[i]

            annotated.append(
                (
                    i,
                    info["p"],
                    info["r"],
                    info["d"],
                )
            )

        print()
        print(
            f"  global_row_type_indices={indices}"
        )

        print(
            f"  metadata={annotated}"
        )

        print(
            f"  row={list(key)}"
        )

    # ------------------------------------------------------------------
    # 4. CUMULATIVE RANK GROWTH
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. CUMULATIVE TRANSITION-RANK AUDIT")
    print("=" * 78)

    cumulative = []

    for name in (
        "p=1 -> 3",
        "p=3 -> 5",
        "p=5 -> 7",
    ):

        cumulative.extend(
            groups[name]
        )

        print(
            f"  through {name}: "
            f"rows={cumulative} "
            f"rank={rank_mod2(submatrix(E, cumulative))}"
        )

    # ------------------------------------------------------------------
    # 5. P=5 -> 7 ZERO TEST
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. FINAL TRANSITION ZERO TEST")
    print("=" * 78)

    p57_rows = groups["p=5 -> 7"]

    p57_zero = all(
        is_zero_row(E[i])
        for i in p57_rows
    )

    print(
        f"  p=5->7 rows={p57_rows}"
    )

    print(
        f"  all_zero={p57_zero}"
    )

    # ------------------------------------------------------------------
    # 6. CROSS-TRANSITION SPAN TEST
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. SPAN-CONTINUATION AUDIT")
    print("=" * 78)

    p13_rows = groups["p=1 -> 3"]
    p35_rows = groups["p=3 -> 5"]

    basis_p13, _ = row_space_basis_mod2(
        submatrix(E, p13_rows)
    )

    basis_p35, _ = row_space_basis_mod2(
        submatrix(E, p35_rows)
    )

    combined_rank = rank_mod2(
        basis_p13 + basis_p35
    )

    print(
        f"  rank(p=1->3)={len(basis_p13)}"
    )

    print(
        f"  rank(p=3->5)={len(basis_p35)}"
    )

    print(
        f"  rank(combined)={combined_rank}"
    )

    print(
        f"  p=3->5_contained_in_p=1->3_span="
        f"{combined_rank == len(basis_p13)}"
    )

    # ------------------------------------------------------------------
    # 7. DISTINGUISHED ROW TEST
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. DISTINGUISHED ROW TEST")
    print("=" * 78)

    distinguished = [0, 5]

    remaining = [
        i
        for i in range(12)
        if i not in distinguished
    ]

    rank_remaining = rank_mod2(
        submatrix(E, remaining)
    )

    print(
        f"  distinguished_core_rows={distinguished}"
    )

    print(
        f"  rank_without_distinguished={rank_remaining}"
    )

    print(
        f"  total_rank={rank_mod2(E)}"
    )

    print(
        f"  rank_drop={rank_mod2(E) - rank_remaining}"
    )

    # ------------------------------------------------------------------
    # 8. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The first 2-adic defect is attached to the 12 retained rows of the
F-core, which inherit the three original transition blocks:

    p=1 -> 3,
    p=3 -> 5,
    p=5 -> 7.

The experiment asks whether the rank-3 obstruction is actually a
three-way phenomenon across those transitions.

The most informative possibilities are:

    * the p=5 -> 7 block is zero;
    * the p=3 -> 5 block contributes only a row direction already
      present in p=1 -> 3;
    * the p=1 -> 3 block contains the full obstruction.

If that happens, the global rank-3 defect is not evenly distributed
over the three transitions.

This does not explain why the arithmetic occurs. It only localizes
the first obstruction inside the exact finite system.

No recurrence is inferred from support alone.
No connection to the original (p,q)-kernel is asserted.
"""
    )

    # ------------------------------------------------------------------
    # 9. FINAL
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. FINAL EXACTNESS")
    print("=" * 78)

    final_ok = (
        len(meta) == 12
        and rank_mod2(E) == 3
        and len(groups["p=1 -> 3"]) == 6
        and len(groups["p=3 -> 5"]) == 5
        and len(groups["p=5 -> 7"]) == 1
    )

    print(
        f"  data_exact={len(meta) == 12}"
    )

    print(
        f"  defect_rank_exact={rank_mod2(E) == 3}"
    )

    print(
        f"  transition_partition_exact={final_ok}"
    )

    print(
        f"  row_type_audit_completed=True"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 150 COMPLETE")


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

