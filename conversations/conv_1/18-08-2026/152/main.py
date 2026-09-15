#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 152R — EXACT 2-ADIC DEFECT / LOCAL-DATA PROVENANCE AUDIT
==============================================================================

Corrected version of Experiment 152.

Experiment 151 established:

    rank(E mod 2) = 3

with exactly:

    3 nonzero defect directions
    1 zero defect direction.

Experiment 152 asked whether these directions are determined directly
by local scalar data from the row that generated them.

This corrected version keeps the original experiment intact while
fixing the metadata-key bug:

    the previous script attempted to access features[i]["d"],
    although "d" was stored only in retained_meta.

This version uses retained_meta explicitly for all p/r/d grouping.

Tests include:

    * exact defect classes;
    * single-feature separation;
    * two-feature separation;
    * q-parity grouping;
    * zero/nonzero defect classification;
    * d-only classification;
    * exact p/r/d grouping.

No fitted recurrence.
No SymPy.
No floating point.
No extrapolation.
No claim about the original (p,q)-kernel.
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


# ============================================================================
# EXACT BAREISS
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
            M[k], M[pivot_row] = (
                M[pivot_row],
                M[k],
            )
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
# BUILD 14x14 SYSTEM
# ============================================================================

def build_system():

    M = []
    meta = []

    for p, pnext in TRANSITIONS:

        for r in range(D[p] + 1):

            d = D[p] - r

            q0 = q_value(p, r)
            q1 = q_value(p, r + 1)

            b = monomial_basis(p, d)

            row = [0] * 14

            for j in range(6):
                row[j] = b[j] * q0

            for j in range(6):
                row[6 + j] = b[j] * q1

            if d == 0:
                row[12] = 1
                row[13] = p

            M.append(row)

            meta.append(
                {
                    "p": p,
                    "pnext": pnext,
                    "r": r,
                    "d": d,
                    "q0": q0,
                    "q1": q1,
                }
            )

    return M, meta


# ============================================================================
# BUILD 12x12 F-CORE
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
            "Boundary block singular."
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

        out = []

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

            out.append(value.numerator)

        core.append(out)

    return core


# ============================================================================
# NORMALIZED COFACTORS
# ============================================================================

def maximal_minor(A, omit_row, omit_col):

    rows = [
        i
        for i in range(len(A))
        if i != omit_row
    ]

    cols = [
        j
        for j in range(len(A[0]))
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
                    "Delta_11 divisibility failure."
                )

            row.append(
                det_minor // DELTA_11
            )

        R.append(row)

    return R


# ============================================================================
# MOD-8 RANK-1 DEFECT
# ============================================================================

def make_defect(R):

    n = len(R)

    pivot = None
    pivot_inverse = None

    for i in range(n):

        for j in range(n):

            x = R[i][j] % 8

            if x == 0:
                continue

            try:
                inv = pow(x, -1, 8)
            except ValueError:
                continue

            pivot = (i, j)
            pivot_inverse = inv
            break

        if pivot is not None:
            break

    if pivot is None:
        raise ArithmeticError(
            "No unit pivot modulo 8."
        )

    i0, j0 = pivot

    u = [
        R[i][j0] % 8
        for i in range(n)
    ]

    v = [
        R[i0][j] * pivot_inverse % 8
        for j in range(n)
    ]

    # Independent validation.
    exact = True

    for i in range(n):

        for j in range(n):

            if (
                u[i] * v[j] % 8
                != R[i][j] % 8
            ):
                exact = False
                break

        if not exact:
            break

    if not exact:
        raise ArithmeticError(
            "Mod-8 rank-one factorization failed."
        )

    E = []

    for i in range(n):

        row = []

        for j in range(n):

            diff = (
                R[i][j]
                - u[i] * v[j]
            )

            if diff % 8 != 0:
                raise ArithmeticError(
                    "R - uv^T is not divisible by 8."
                )

            row.append(
                (diff // 8) & 1
            )

        E.append(row)

    return E


# ============================================================================
# MOD-2 RANK
# ============================================================================

def rref_mod2(A):

    M = [
        [int(x) & 1 for x in row]
        for row in A
    ]

    rows = len(M)
    cols = len(M[0])

    pivot_cols = []
    pivot_row = 0

    for col in range(cols):

        pivot = None

        for r in range(pivot_row, rows):

            if M[r][col] == 1:
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

        pivot_cols.append(col)
        pivot_row += 1

        if pivot_row == rows:
            break

    return M, pivot_cols


def rank_mod2(A):
    _, pivots = rref_mod2(A)
    return len(pivots)


# ============================================================================
# DEFECT CLASSES
# ============================================================================

def row_tuple(row):
    return tuple(int(x) & 1 for x in row)


def row_classes(E):

    classes = {}

    for i, row in enumerate(E):

        key = row_tuple(row)

        classes.setdefault(
            key,
            [],
        ).append(i)

    return classes


# ============================================================================
# LOCAL FEATURES
# ============================================================================

def feature_dict(info):

    p = info["p"]
    r = info["r"]
    d = info["d"]
    q0 = info["q0"]
    q1 = info["q1"]

    return {
        "p_mod2": p % 2,
        "r_mod2": r % 2,
        "d_mod2": d % 2,
        "p_mod4": p % 4,
        "r_mod3": r % 3,
        "d_mod3": d % 3,
        "q0_mod2": q0 & 1,
        "q1_mod2": q1 & 1,
        "q0_xor_q1_mod2": (q0 ^ q1) & 1,
        "q0_nonzero": int(q0 != 0),
        "q1_nonzero": int(q1 != 0),
        "terminal": int(d == 0),
        "near_terminal": int(d <= 1),
    }


# ============================================================================
# FEATURE SEPARATION
# ============================================================================

def separates_exactly(
    labels,
    values,
    indices,
):
    """
    Test whether equal feature values imply equal target labels.
    """

    buckets = {}

    for i in indices:

        key = values[i]

        buckets.setdefault(
            key,
            set(),
        ).add(
            labels[i]
        )

    return all(
        len(targets) <= 1
        for targets in buckets.values()
    )


def smallest_exact_pair(
    labels,
    feature_names,
    features,
    indices,
):

    for a, b in combinations(
        feature_names,
        2,
    ):

        buckets = {}

        for i in indices:

            key = (
                features[i][a],
                features[i][b],
            )

            buckets.setdefault(
                key,
                set(),
            ).add(
                labels[i]
            )

        if all(
            len(targets) <= 1
            for targets in buckets.values()
        ):
            return a, b

    return None


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 152R — EXACT 2-ADIC DEFECT / "
        "LOCAL-DATA PROVENANCE AUDIT"
    )
    print("=" * 78)

    M, full_meta = build_system()
    core = build_core(M)
    R = normalized_cofactors(core)
    E = make_defect(R)

    retained_meta = [
        info
        for idx, info in enumerate(full_meta)
        if idx not in (5, 10)
    ]

    # ------------------------------------------------------------------
    # 1. VALIDATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. EXACT DATA VALIDATION")
    print("=" * 78)

    classes = row_classes(E)

    nonzero_classes = [
        key
        for key in classes
        if any(key)
    ]

    print(
        f"  retained_rows={len(retained_meta)}"
    )

    print(
        f"  defect_shape={len(E)}x{len(E[0])}"
    )

    print(
        f"  defect_rank={rank_mod2(E)}"
    )

    print(
        f"  distinct_total_row_types={len(classes)}"
    )

    print(
        f"  distinct_nonzero_row_types="
        f"{len(nonzero_classes)}"
    )

    # ------------------------------------------------------------------
    # 2. LOCAL DATA TABLE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. LOCAL DATA / DEFECT CLASS TABLE")
    print("=" * 78)

    sorted_classes = sorted(
        classes.items(),
        key=lambda item: item[1][0],
    )

    class_id = {}

    for cid, (_, indices) in enumerate(
        sorted_classes
    ):
        for i in indices:
            class_id[i] = cid

    feature_rows = []

    for i, info in enumerate(
        retained_meta
    ):

        f = feature_dict(info)
        feature_rows.append(f)

        print(
            f"  row={i}: "
            f"p={info['p']} "
            f"r={info['r']} "
            f"d={info['d']} "
            f"q0={info['q0']} "
            f"q1={info['q1']} "
            f"class={class_id[i]}"
        )

    # ------------------------------------------------------------------
    # 3. EXACT DEFECT CLASSES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. EXACT DEFECT CLASSES")
    print("=" * 78)

    for cid, (key, indices) in enumerate(
        sorted_classes
    ):

        locations = []

        for i in indices:

            info = retained_meta[i]

            locations.append(
                (
                    info["p"],
                    info["r"],
                    info["d"],
                )
            )

        print()
        print(
            f"  class_{cid}: rows={indices}"
        )

        print(
            f"    vector={list(key)}"
        )

        print(
            f"    locations={locations}"
        )

    # ------------------------------------------------------------------
    # 4. SINGLE FEATURE SEPARATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. SINGLE LOCAL FEATURE SEPARATION")
    print("=" * 78)

    feature_names = list(
        feature_rows[0].keys()
    )

    single_results = {}

    for name in feature_names:

        values = {
            i: feature_rows[i][name]
            for i in range(len(feature_rows))
        }

        ok = separates_exactly(
            class_id,
            values,
            range(len(feature_rows)),
        )

        single_results[name] = ok

        print(
            f"  {name}: separates_exactly={ok}"
        )

    # ------------------------------------------------------------------
    # 5. TWO-FEATURE SEPARATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. TWO-FEATURE SEPARATION SEARCH")
    print("=" * 78)

    pair = smallest_exact_pair(
        class_id,
        feature_names,
        feature_rows,
        range(len(feature_rows)),
    )

    print(
        f"  smallest_exact_feature_pair={pair}"
    )

    if pair is not None:

        a, b = pair

        buckets = {}

        for i in range(
            len(feature_rows)
        ):

            key = (
                feature_rows[i][a],
                feature_rows[i][b],
            )

            buckets.setdefault(
                key,
                [],
            ).append(i)

        for key, rows in buckets.items():

            print(
                f"    {a},{b}={key}: "
                f"rows={rows} "
                f"classes="
                f"{sorted(set(class_id[i] for i in rows))}"
            )

    # ------------------------------------------------------------------
    # 6. q-PARITY AUDIT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. q-PARITY AUDIT")
    print("=" * 78)

    q_buckets = {}

    for i, f in enumerate(
        feature_rows
    ):

        key = (
            f["q0_mod2"],
            f["q1_mod2"],
        )

        q_buckets.setdefault(
            key,
            [],
        ).append(i)

    for key, rows in q_buckets.items():

        print()
        print(
            f"  q_parity={key}: rows={rows}"
        )

        print(
            f"    classes="
            f"{sorted(set(class_id[i] for i in rows))}"
        )

    # ------------------------------------------------------------------
    # 7. ZERO/NONZERO DEFECT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. ZERO / NONZERO DEFECT AUDIT")
    print("=" * 78)

    zero_rows = [
        i
        for i, row in enumerate(E)
        if not any(row)
    ]

    nonzero_rows = [
        i
        for i, row in enumerate(E)
        if any(row)
    ]

    zero_labels = {
        i: int(i in nonzero_rows)
        for i in range(len(E))
    }

    print(
        f"  zero_defect_rows={zero_rows}"
    )

    print(
        f"  nonzero_defect_rows={nonzero_rows}"
    )

    for name in feature_names:

        values = {
            i: feature_rows[i][name]
            for i in range(len(feature_rows))
        }

        ok = separates_exactly(
            zero_labels,
            values,
            range(len(feature_rows)),
        )

        print(
            f"  {name}: "
            f"exact_zero_nonzero_separator={ok}"
        )

    # ------------------------------------------------------------------
    # 8. CORRECT d-ONLY TEST
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. BOUNDARY-COORDINATE RECONFIRMATION")
    print("=" * 78)

    d_values = {
        i: retained_meta[i]["d"]
        for i in range(len(retained_meta))
    }

    d_separates = separates_exactly(
        class_id,
        d_values,
        range(len(retained_meta)),
    )

    print(
        f"  d_alone_separates_defect_classes="
        f"{d_separates}"
    )

    # ------------------------------------------------------------------
    # 9. p/r/d COMBINATIONS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. EXACT p/r/d COMBINATION AUDIT")
    print("=" * 78)

    combinations_to_test = {
        "p": lambda info: info["p"],
        "r": lambda info: info["r"],
        "d": lambda info: info["d"],
        "(p,r)": lambda info: (
            info["p"],
            info["r"],
        ),
        "(p,d)": lambda info: (
            info["p"],
            info["d"],
        ),
        "(r,d)": lambda info: (
            info["r"],
            info["d"],
        ),
        "(p,r,d)": lambda info: (
            info["p"],
            info["r"],
            info["d"],
        ),
    }

    for name, fn in combinations_to_test.items():

        values = {
            i: fn(retained_meta[i])
            for i in range(
                len(retained_meta)
            )
        }

        ok = separates_exactly(
            class_id,
            values,
            range(len(retained_meta)),
        )

        print(
            f"  {name}: separates_exactly={ok}"
        )

    # ------------------------------------------------------------------
    # 10. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
Experiment 151 ruled out d=D(p)-r as a complete label for the first
2-adic defect direction.

Experiment 152R therefore tests the immediate local quantities from
which each row is constructed:

    p, r, d, q_p(r), q_p(r+1),

and their simple parity/nonzero reductions.

The purpose is not to find an arbitrary classifier.

The meaningful outcomes are:

    * a very small exact local separator exists;
    * q-parity explains the zero/nonzero split;
    * (p,d) or (r,d) explains the defect direction;
    * no small local description exists.

A negative result is informative because it would indicate that the
rank-3 defect is not simply a local scalar property of an individual
transition row.

No recurrence is fitted and no general theorem is inferred from the
finite classification.
"""
    )

    # ------------------------------------------------------------------
    # 11. FINAL
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("11. FINAL EXACTNESS")
    print("=" * 78)

    final_ok = (
        len(retained_meta) == 12
        and len(E) == 12
        and len(E[0]) == 12
        and len(classes) == 4
        and len(nonzero_classes) == 3
        and rank_mod2(E) == 3
    )

    print(
        f"  data_exact="
        f"{len(retained_meta) == 12 and len(E) == 12}"
    )

    print(
        f"  defect_rank_exact="
        f"{rank_mod2(E) == 3}"
    )

    print(
        f"  four_total_row_types_exact="
        f"{len(classes) == 4}"
    )

    print(
        f"  three_nonzero_directions_exact="
        f"{len(nonzero_classes) == 3}"
    )

    print(
        "  local_provenance_search_completed=True"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 152R COMPLETE")


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