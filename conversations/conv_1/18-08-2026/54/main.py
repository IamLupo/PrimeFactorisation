#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
R=5 EXPERIMENT 54
=================

SUPPORT-REMOVED COEFFICIENT-MATRIX RANK / SEPARABILITY TEST

Input:
    exact supplied R=5 discrepancy polynomials

Structure:
    E_j(k,D) = Support_j(D) * Q_j(k,D)

Questions:
    1. What is the exact coefficient-matrix rank of Q_j?
    2. Does shifting D reveal lower rank?
    3. Is any j exactly separable in k and D?
    4. Do the coefficient matrices have small-rank structure?
    5. Does j=5 exhibit special terminal behaviour?

No r=6.
No pq-kernel expansion.
No universal replacement law.
"""

from __future__ import annotations

import math
from fractions import Fraction
from itertools import combinations
from typing import Dict, List, Sequence, Tuple

import sympy as sp


# ============================================================================
# SYMBOLS
# ============================================================================

D = sp.Symbol("D")
K = sp.Symbol("K")
X = sp.Symbol("X")
Y = sp.Symbol("Y")


# ============================================================================
# GRID
# ============================================================================

K_VALUES = [3, 5, 7, 9, 11, 13]
J_VALUES = [0, 1, 2, 3, 4, 5]
D_VALUES = [6, 8, 10, 12, 14, 16]

EXPECTED_POINTS = (
    len(K_VALUES)
    * len(J_VALUES)
    * len(D_VALUES)
)


# ============================================================================
# EXACT SUPPORT STRUCTURE
# ============================================================================

SUPPORTS = {
    0: sp.Integer(1),
    1: D - 6,
    2: D - 6,
    3: (D - 8) * (D - 6),
    4: (D - 8) * (D - 6),
    5: (D - 12) * (D - 10) * (D - 8) * (D - 6),
}


D0 = {
    0: 6,
    1: 8,
    2: 8,
    3: 10,
    4: 10,
    5: 14,
}


# ============================================================================
# EXACT INPUT POLYNOMIALS
# ============================================================================

E_POLYS: Dict[Tuple[int, int], sp.Expr] = {

    # ------------------------------------------------------------------
    # j = 0
    # ------------------------------------------------------------------

    (3, 0):
        -(7363*D**5
          - 395215*D**4
          + 8629700*D**3
          - 91891700*D**2
          + 472005792*D
          - 932561280) / 960,

    (5, 0):
        (76393*D**5
         - 4305640*D**4
         + 90582740*D**3
         - 916615760*D**2
         + 4483939392*D
         - 8489617920) / 1920,

    (7, 0):
        -(15541*D**5
          - 632044*D**4
          + 11202092*D**3
          - 98824400*D**2
          + 426557568*D
          - 720147456) / 384,

    (9, 0):
        -(615341*D**5
          - 28009020*D**4
          + 533637820*D**3
          - 5012129040*D**2
          + 22981059264*D
          - 41124902400) / 3840,

    (11, 0):
        -13*(9677*D**5
             - 485334*D**4
             + 9685172*D**3
             - 94225800*D**2
             + 445818272*D
             - 820598016) / 128,

    (13, 0):
        -119*(25623*D**5
              - 1350100*D**4
              + 27793940*D**3
              - 276874160*D**2
              + 1335166912*D
              - 2495823360) / 960,

    # ------------------------------------------------------------------
    # j = 1
    # ------------------------------------------------------------------

    (3, 1):
        -(D - 6) * (
            149*D**4
            - 2683*D**3
            - 2402*D**2
            + 398140*D
            - 2180640
        ) / 96,

    (5, 1):
        -7*(D - 6) * (
            181*D**4
            - 8406*D**3
            + 150268*D**2
            - 1168344*D
            + 3311680
        ) / 48,

    (7, 1):
        (D - 6) * (
            1679*D**4
            - 96396*D**3
            + 1753904*D**2
            - 13562106*D
            + 38170160
        ) / 80,

    (9, 1):
        -11*(D - 6) * (
            11371*D**4
            - 508124*D**3
            + 8504336*D**2
            - 62038984*D
            + 166242240
        ) / 480,

    (11, 1):
        -91*(D - 6) * (
            2863*D**4
            - 136182*D**3
            + 2382328*D**2
            - 17966172*D
            + 49397920
        ) / 480,

    (13, 1):
        -17*(D - 6) * (
            531*D**4
            - 79544*D**3
            + 2118616*D**2
            - 19978264*D
            + 63103040
        ) / 120,

    # ------------------------------------------------------------------
    # j = 2
    # ------------------------------------------------------------------

    (3, 2):
        (D - 6) * (
            22707*D**4
            - 1110703*D**3
            + 19642662*D**2
            - 149850488*D
            + 416996160
        ) / 960,

    (5, 2):
        -(D - 6) * (
            2312*D**4
            - 103197*D**3
            + 1707416*D**2
            - 12332148*D
            + 32788640
        ) / 48,

    (7, 2):
        -(D - 6) * (
            15669*D**4
            - 732214*D**3
            + 12707648*D**2
            - 95638016*D
            + 263119360
        ) / 384,

    (9, 2):
        -(D - 6) * (
            21077*D**4
            - 1071324*D**3
            + 19742740*D**2
            - 155018688*D
            + 439799040
        ) / 384,

    (11, 2):
        7*(D - 6) * (
            4955*D**4
            - 182498*D**3
            + 2446376*D**2
            - 14262136*D
            + 30560960
        ) / 192,

    (13, 2):
        17*(D - 6) * (
            94421*D**4
            - 3956744*D**3
            + 61568116*D**2
            - 421903984*D
            + 1073639040
        ) / 960,

    # ------------------------------------------------------------------
    # j = 3
    # ------------------------------------------------------------------

    (3, 3):
        -7*(D - 8)*(D - 6) * (
            131*D**3
            - 5546*D**2
            + 77468*D
            - 350330
        ) / 80,

    (5, 3):
        (D - 8)*(D - 6) * (
            4051*D**3
            - 157178*D**2
            + 1973512*D
            - 8082000
        ) / 128,

    (7, 3):
        -11*(D - 8)*(D - 6) * (
            2463*D**3
            - 93638*D**2
            + 1173424*D
            - 4819640
        ) / 960,

    (9, 3):
        13*(D - 8)*(D - 6) * (
            595*D**3
            - 20198*D**2
            + 226612*D
            - 842940
        ) / 96,

    (11, 3):
        (D - 8)*(D - 6) * (
            3358*D**3
            - 117757*D**2
            + 1366672*D
            - 5251345
        ) / 8,

    (13, 3):
        17*(D - 8)*(D - 6) * (
            163617*D**3
            - 5798122*D**2
            + 67971936*D
            - 263566160
        ) / 1920,

    # ------------------------------------------------------------------
    # j = 4
    # ------------------------------------------------------------------

    (3, 4):
        -(D - 16)*(D - 8)*(D - 6) * (
            1892*D**2
            - 40889*D
            + 219490
        ) / 64,

    (5, 4):
        (D - 8)*(D - 6) * (
            8315*D**3
            - 272506*D**2
            + 2955352*D
            - 10636160
        ) / 768,

    (7, 4):
        (D - 8)*(D - 6) * (
            124993*D**3
            - 4453618*D**2
            + 52402784*D
            - 203769440
        ) / 3840,

    (9, 4):
        (D - 8)*(D - 6) * (
            5649*D**3
            - 201758*D**2
            + 2380888*D
            - 9284240
        ) / 64,

    (11, 4):
        (D - 8)*(D - 6) * (
            487679*D**3
            - 17467064*D**2
            + 206706812*D
            - 808115120
        ) / 1920,

    (13, 4):
        (D - 8)*(D - 6) * (
            2350941*D**3
            - 84327646*D**2
            + 999345048*D
            - 3911794880
        ) / 3840,

    # ------------------------------------------------------------------
    # j = 5
    # ------------------------------------------------------------------

    (3, 5):
        (D - 12)*(D - 10)*(D - 8)*(D - 6) * (
            22277*D - 371378
        ) / 960,

    (5, 5):
        -1295*(D - 14)*(D - 12)*(D - 10)*(D - 8)*(D - 6) / 96,

    (7, 5):
        7007*(D - 14)*(D - 12)*(D - 10)*(D - 8)*(D - 6) / 480,

    (9, 5):
        1377*(D - 14)*(D - 12)*(D - 10)*(D - 8)*(D - 6) / 40,

    (11, 5):
        46189*(D - 14)*(D - 12)*(D - 10)*(D - 8)*(D - 6) / 640,

    (13, 5):
        -7*(D - 12)*(D - 10)*(D - 8)*(D - 6) * (
            37989*D - 626891
        ) / 480,
}


# ============================================================================
# BASIC HELPERS
# ============================================================================

def clean(expr: sp.Expr) -> sp.Expr:
    return sp.factor(
        sp.cancel(
            sp.expand(expr)
        )
    )


def get_E(k: int, j: int) -> sp.Expr:
    return clean(
        E_POLYS[(k, j)]
    )


def get_support(j: int) -> sp.Expr:
    return SUPPORTS[j]


def get_Q(k: int, j: int) -> sp.Expr:
    """
    Exact support-removed quotient.

        E = Support * Q
    """
    E = get_E(k, j)
    S = get_support(j)

    quotient = clean(
        sp.cancel(
            E / S
        )
    )

    # Exact divisibility check.
    E_poly = sp.Poly(
        E,
        D,
        domain=sp.QQ,
    )

    S_poly = sp.Poly(
        S,
        D,
        domain=sp.QQ,
    )

    remainder = E_poly.rem(
        S_poly
    )

    if not remainder.is_zero:
        raise RuntimeError(
            f"support divisibility failure "
            f"for k={k}, j={j}"
        )

    reconstructed = clean(
        S * quotient
    )

    if clean(
        reconstructed - E
    ) != 0:
        raise RuntimeError(
            f"support quotient reconstruction failure "
            f"for k={k}, j={j}"
        )

    return quotient


def shifted_expr(
    expr: sp.Expr,
    variable: sp.Symbol,
    shift: int,
) -> sp.Expr:
    """
    Let

        variable = D - shift

    so

        D = variable + shift.
    """
    return clean(
        expr.subs(
            D,
            variable + shift,
        )
    )


def coefficient_vector(
    expr: sp.Expr,
    variable: sp.Symbol,
) -> List[sp.Rational]:
    poly = sp.Poly(
        sp.expand(expr),
        variable,
        domain=sp.QQ,
    )

    degree = poly.degree()

    if degree is None or degree == -sp.oo:
        return [sp.Rational(0)]

    return [
        sp.Rational(
            poly.nth(i)
        )
        for i in range(
            int(degree) + 1
        )
    ]


def coefficient_matrix(
    j: int,
    variable: sp.Symbol,
    shift: int,
) -> sp.Matrix:
    """
    Rows:
        k = 3,5,...,13

    Columns:
        coefficients of the shifted quotient

            Q_j(k, D)
              = q_0(k)
              + q_1(k) Z
              + ...
    """
    rows: List[List[sp.Rational]] = []

    for k in K_VALUES:
        Q = get_Q(
            k,
            j,
        )

        shifted = shifted_expr(
            Q,
            variable,
            shift,
        )

        rows.append(
            coefficient_vector(
                shifted,
                variable,
            )
        )

    width = max(
        len(row)
        for row in rows
    )

    padded_rows = []

    for row in rows:
        padded_rows.append(
            row + [
                sp.Rational(0)
            ] * (
                width - len(row)
            )
        )

    return sp.Matrix(
        padded_rows
    )


def exact_rank(
    matrix: sp.Matrix,
) -> int:
    return int(
        matrix.rank()
    )


def primitive_integer_form(
    values: Sequence[sp.Rational],
) -> Tuple[List[int], Fraction]:
    """
    Exact rational vector -> primitive integer vector + scale.

    Robust for:
      * zero vector
      * one coefficient
      * constant vectors
      * arbitrary rational vectors
    """
    vals = [
        sp.Rational(v)
        for v in values
    ]

    if not vals:
        return (
            [],
            Fraction(0, 1),
        )

    if all(
        v == 0
        for v in vals
    ):
        return (
            [0] * len(vals),
            Fraction(0, 1),
        )

    den_lcm = 1

    for v in vals:
        den_lcm = math.lcm(
            den_lcm,
            int(sp.denom(v)),
        )

    integer_values = [
        int(v * den_lcm)
        for v in vals
    ]

    nonzero_abs = [
        abs(v)
        for v in integer_values
        if v != 0
    ]

    content = nonzero_abs[0]

    for value in nonzero_abs[1:]:
        content = math.gcd(
            content,
            value,
        )

    primitive = [
        value // content
        for value in integer_values
    ]

    # Normalize sign by first nonzero entry.
    sign = 1

    for value in primitive:
        if value != 0:
            if value < 0:
                primitive = [
                    -x
                    for x in primitive
                ]
                sign = -1
            break

    scale = Fraction(
        sign * content,
        den_lcm,
    )

    return (
        primitive,
        scale,
    )


def rank_one_factorization(
    matrix: sp.Matrix,
):
    """
    If matrix has rank one, return row/column generators.
    Otherwise return None.

    This is intentionally exact.
    """
    if matrix.rank() != 1:
        return None

    # Locate one nonzero pivot entry.
    pivot = None

    for r in range(matrix.rows):
        for c in range(matrix.cols):
            if matrix[r, c] != 0:
                pivot = (
                    r,
                    c,
                )
                break

        if pivot is not None:
            break

    if pivot is None:
        return None

    r0, c0 = pivot

    column = matrix[:, c0]
    row = matrix[r0, :]

    pivot_value = matrix[
        r0,
        c0
    ]

    reconstructed = (
        column
        * row
        / pivot_value
    )

    if reconstructed != matrix:
        raise RuntimeError(
            "Internal rank-one reconstruction failure."
        )

    return (
        column,
        row,
        pivot_value,
    )


# ============================================================================
# SECTION 0
# ============================================================================

def section_0_validation() -> None:
    print("=" * 78)
    print("0. DATA VALIDATION")
    print("=" * 78)

    print(
        f"k values = {K_VALUES}"
    )
    print(
        f"j values = {J_VALUES}"
    )
    print(
        f"D values = {D_VALUES}"
    )
    print(
        f"points   = {EXPECTED_POINTS}"
    )

    expected_keys = {
        (k, j)
        for k in K_VALUES
        for j in J_VALUES
    }

    actual_keys = set(
        E_POLYS.keys()
    )

    if actual_keys != expected_keys:
        missing = sorted(
            expected_keys - actual_keys
        )
        extra = sorted(
            actual_keys - expected_keys
        )

        raise RuntimeError(
            f"grid mismatch: "
            f"missing={missing}, extra={extra}"
        )

    checked = 0

    for k in K_VALUES:
        for j in J_VALUES:
            for d in D_VALUES:
                value = clean(
                    get_E(k, j).subs(
                        D,
                        d,
                    )
                )

                if not value.is_Rational:
                    raise RuntimeError(
                        f"non-rational value "
                        f"at k={k}, j={j}, D={d}"
                    )

                checked += 1

    print(
        "grid status = OK"
    )
    print(
        f"checked points = {checked}"
    )
    print()


# ============================================================================
# SECTION 1
# ============================================================================

def section_1_support_quotients() -> None:
    print("=" * 78)
    print("1. EXACT SUPPORT-REMOVED QUOTIENT DEGREES")
    print("=" * 78)

    for j in J_VALUES:
        support = clean(
            get_support(j)
        )

        support_degree = (
            0
            if support == 1
            else int(
                sp.degree(
                    support,
                    D,
                )
            )
        )

        print(
            f"j={j}"
        )
        print(
            f"  support={support}"
        )
        print(
            f"  support_degree={support_degree}"
        )

        for k in K_VALUES:
            Q = get_Q(
                k,
                j,
            )

            degree = int(
                sp.degree(
                    Q,
                    D,
                )
            )

            print(
                f"  k={k}: "
                f"Q_degree={degree}"
            )

        print()


# ============================================================================
# SECTION 2
# ============================================================================

def section_2_rank_in_D() -> None:
    print("=" * 78)
    print("2. EXACT COEFFICIENT-MATRIX RANK IN D")
    print("=" * 78)

    for j in J_VALUES:
        matrix = coefficient_matrix(
            j,
            D,
            0,
        )

        print(
            f"j={j}: "
            f"matrix_shape={matrix.shape} "
            f"rank={exact_rank(matrix)}"
        )

        print(
            f"  matrix={matrix}"
        )

        print()


# ============================================================================
# SECTION 3
# ============================================================================

def section_3_rank_shifted_X() -> None:
    print("=" * 78)
    print("3. COEFFICIENT-MATRIX RANK AFTER X = D - 6")
    print("=" * 78)

    for j in J_VALUES:
        matrix = coefficient_matrix(
            j,
            X,
            6,
        )

        print(
            f"j={j}: "
            f"shape={matrix.shape} "
            f"rank={exact_rank(matrix)}"
        )

        print(
            f"  X-coefficient matrix={matrix}"
        )

        print()


# ============================================================================
# SECTION 4
# ============================================================================

def section_4_rank_shifted_boundary() -> None:
    print("=" * 78)
    print("4. COEFFICIENT-MATRIX RANK AFTER Y = D - D0(j)")
    print("=" * 78)

    for j in J_VALUES:
        shift = D0[j]

        matrix = coefficient_matrix(
            j,
            Y,
            shift,
        )

        print(
            f"j={j}: "
            f"D0={shift} "
            f"shape={matrix.shape} "
            f"rank={exact_rank(matrix)}"
        )

        print(
            f"  Y-coefficient matrix={matrix}"
        )

        print()


# ============================================================================
# SECTION 5
# ============================================================================

def section_5_rank_one_test() -> None:
    print("=" * 78)
    print("5. EXACT RANK-ONE / SEPARABILITY TEST")
    print("=" * 78)

    for j in J_VALUES:
        matrix = coefficient_matrix(
            j,
            Y,
            D0[j],
        )

        rank = exact_rank(
            matrix
        )

        print(
            f"j={j}: rank={rank}"
        )

        if rank == 1:
            print(
                "  RESULT: EXACT SEPARABILITY"
            )

            factorization = rank_one_factorization(
                matrix
            )

            if factorization is not None:
                column, row, pivot = factorization

                print(
                    f"  pivot={pivot}"
                )
                print(
                    f"  column_generator={column}"
                )
                print(
                    f"  row_generator={row}"
                )

        else:
            print(
                "  RESULT: NOT RANK-ONE"
            )

        print()


# ============================================================================
# SECTION 6
# ============================================================================

def section_6_minor_audit() -> None:
    print("=" * 78)
    print("6. EXACT 2x2 MINOR AUDIT")
    print("=" * 78)

    for j in J_VALUES:
        matrix = coefficient_matrix(
            j,
            Y,
            D0[j],
        )

        rows, cols = matrix.shape

        total_minors = (
            math.comb(rows, 2)
            * math.comb(cols, 2)
            if rows >= 2 and cols >= 2
            else 0
        )

        nonzero_minors = []

        for r1, r2 in combinations(
            range(rows),
            2,
        ):
            for c1, c2 in combinations(
                range(cols),
                2,
            ):
                submatrix = matrix.extract(
                    [r1, r2],
                    [c1, c2],
                )

                det = sp.factor(
                    submatrix.det()
                )

                if det != 0:
                    nonzero_minors.append(
                        (
                            r1,
                            r2,
                            c1,
                            c2,
                            det,
                        )
                    )

        print(
            f"j={j}: "
            f"matrix_shape={matrix.shape} "
            f"total_2x2_minors={total_minors} "
            f"nonzero_2x2_minors={len(nonzero_minors)}"
        )

        # Showing only a bounded sample keeps result.txt manageable.
        for item in nonzero_minors[:12]:
            r1, r2, c1, c2, det = item

            print(
                f"  rows=({r1},{r2}) "
                f"cols=({c1},{c2}) "
                f"det={det}"
            )

        if len(nonzero_minors) > 12:
            print(
                f"  ... "
                f"{len(nonzero_minors) - 12} more nonzero minors"
            )

        print()


# ============================================================================
# SECTION 7
# ============================================================================

def section_7_normalized_profiles() -> None:
    print("=" * 78)
    print("7. PRIMITIVE INTEGER NORMALIZATION OF Y-COEFFICIENT ROWS")
    print("=" * 78)

    for j in J_VALUES:
        matrix = coefficient_matrix(
            j,
            Y,
            D0[j],
        )

        print(
            f"j={j}"
        )

        for row_index, k in enumerate(
            K_VALUES
        ):
            row = [
                sp.Rational(
                    matrix[
                        row_index,
                        col,
                    ]
                )
                for col in range(
                    matrix.cols
                )
            ]

            primitive, scale = primitive_integer_form(
                row
            )

            print(
                f"  k={k}:"
            )
            print(
                f"    row={row}"
            )
            print(
                f"    primitive={primitive}"
            )
            print(
                f"    scale={scale}"
            )

        print()


# ============================================================================
# SECTION 8
# ============================================================================

def section_8_column_rank_growth() -> None:
    print("=" * 78)
    print("8. BOUNDARY-EXPANSION COLUMN-RANK GROWTH")
    print("=" * 78)

    for j in J_VALUES:
        matrix = coefficient_matrix(
            j,
            Y,
            D0[j],
        )

        print(
            f"j={j}"
        )

        for width in range(
            1,
            matrix.cols + 1,
        ):
            submatrix = matrix[
                :,
                :width,
            ]

            print(
                f"  first_{width}_columns: "
                f"rank={exact_rank(submatrix)}"
            )

        print()


# ============================================================================
# SECTION 9
# ============================================================================

def section_9_k_difference_rank() -> None:
    print("=" * 78)
    print("9. FINITE-DIFFERENCE-IN-k RANK EVOLUTION")
    print("=" * 78)

    for j in J_VALUES:
        print(
            f"j={j}"
        )

        matrix = coefficient_matrix(
            j,
            Y,
            D0[j],
        )

        current = matrix
        order = 0

        while current.rows > 0:
            print(
                f"  difference_order={order}: "
                f"shape={current.shape} "
                f"rank={exact_rank(current)}"
            )

            if current.rows == 1:
                break

            next_rows = []

            for r in range(
                current.rows - 1
            ):
                next_rows.append(
                    current.row(r + 1)
                    - current.row(r)
                )

            if not next_rows:
                break

            current = sp.Matrix.vstack(
                *next_rows
            )

            order += 1

        print()


# ============================================================================
# SECTION 10
# ============================================================================

def section_10_pairwise_row_proportionality() -> None:
    print("=" * 78)
    print("10. EXACT PAIRWISE k-ROW PROPORTIONALITY")
    print("=" * 78)

    for j in J_VALUES:
        matrix = coefficient_matrix(
            j,
            Y,
            D0[j],
        )

        print(
            f"j={j}"
        )

        for i in range(
            len(K_VALUES) - 1
        ):
            row_a = matrix.row(i)
            row_b = matrix.row(i + 1)

            proportional = False
            ratio = None

            pivot_column = None

            for c in range(
                matrix.cols
            ):
                if (
                    row_a[c] != 0
                    and row_b[c] != 0
                ):
                    pivot_column = c
                    break

            if pivot_column is not None:
                candidate_ratio = sp.cancel(
                    row_b[pivot_column]
                    / row_a[pivot_column]
                )

                if all(
                    sp.simplify(
                        row_b[c]
                        - candidate_ratio * row_a[c]
                    ) == 0
                    for c in range(
                        matrix.cols
                    )
                ):
                    proportional = True
                    ratio = candidate_ratio

            elif (
                all(
                    row_a[c] == 0
                    for c in range(
                        matrix.cols
                    )
                )
                and all(
                    row_b[c] == 0
                    for c in range(
                        matrix.cols
                    )
                )
            ):
                proportional = True
                ratio = sp.Integer(1)

            print(
                f"  "
                f"k={K_VALUES[i]} -> "
                f"k={K_VALUES[i+1]}: "
                f"proportional={proportional}"
            )

            if proportional:
                print(
                    f"    ratio={ratio}"
                )

        print()


# ============================================================================
# SECTION 11
# ============================================================================

def section_11_j5_terminal_rank() -> None:
    print("=" * 78)
    print("11. SPECIAL j=5 TERMINAL RANK AUDIT")
    print("=" * 78)

    j = 5
    shift = D0[j]

    matrix = coefficient_matrix(
        j,
        Y,
        shift,
    )

    print(
        f"D0={shift}"
    )
    print(
        f"Y=D-{shift}"
    )
    print(
        f"matrix_shape={matrix.shape}"
    )
    print(
        f"rank={exact_rank(matrix)}"
    )
    print(
        f"matrix={matrix}"
    )

    print()

    terminal_ladder = (
        (K - 5)
        * (K - 7)
        * (K - 9)
        * (K - 11)
    )

    print(
        "terminal ladder="
        f"{sp.factor(terminal_ladder)}"
    )

    for d in [14, 16]:
        values = [
            (
                sp.Integer(k),
                sp.Rational(
                    get_E(k, 5).subs(
                        D,
                        d,
                    )
                ),
            )
            for k in K_VALUES
        ]

        interpolated = clean(
            sp.interpolate(
                values,
                K,
            )
        )

        print(
            f"D={d}:"
        )
        print(
            f"  P(K)={sp.factor(interpolated)}"
        )

        remainder = sp.Poly(
            interpolated,
            K,
            domain=sp.QQ,
        ).rem(
            sp.Poly(
                terminal_ladder,
                K,
                domain=sp.QQ,
            )
        )

        survives = remainder.is_zero

        print(
            f"  ladder_survives={survives}"
        )

        if survives:
            quotient = clean(
                interpolated
                / terminal_ladder
            )

            print(
                f"  quotient="
                f"{sp.factor(quotient)}"
            )

    print()


# ============================================================================
# SECTION 12
# ============================================================================

def section_12_rank_summary() -> None:
    print("=" * 78)
    print("12. COMPACT RANK SUMMARY")
    print("=" * 78)

    print(
        "j | deg(Q) | rank(D) | rank(X=D-6) | rank(Y=D-D0)"
    )
    print("-" * 62)

    for j in J_VALUES:
        representative = get_Q(
            K_VALUES[0],
            j,
        )

        degree = int(
            sp.degree(
                representative,
                D,
            )
        )

        rank_D = exact_rank(
            coefficient_matrix(
                j,
                D,
                0,
            )
        )

        rank_X = exact_rank(
            coefficient_matrix(
                j,
                X,
                6,
            )
        )

        rank_Y = exact_rank(
            coefficient_matrix(
                j,
                Y,
                D0[j],
            )
        )

        print(
            f"{j:1d} | "
            f"{degree:6d} | "
            f"{rank_D:7d} | "
            f"{rank_X:11d} | "
            f"{rank_Y:12d}"
        )

    print()


# ============================================================================
# SECTION 13
# ============================================================================

def section_13_support_reconstruction() -> None:
    print("=" * 78)
    print("13. EXACT SUPPORT-QUOTIENT RECONSTRUCTION")
    print("=" * 78)

    tested = 0
    failures = []

    for k in K_VALUES:
        for j in J_VALUES:
            E = get_E(
                k,
                j,
            )

            S = get_support(
                j
            )

            Q = get_Q(
                k,
                j,
            )

            reconstructed = clean(
                S * Q
            )

            if clean(
                reconstructed - E
            ) != 0:
                failures.append(
                    (
                        k,
                        j,
                        reconstructed,
                        E,
                    )
                )

            tested += 1

    print(
        f"tested={tested}"
    )
    print(
        "support reconstruction failures="
        f"{len(failures)}"
    )

    if failures:
        for failure in failures[:10]:
            print(
                f"  failure={failure}"
            )

        raise RuntimeError(
            "Support-quotient reconstruction failed."
        )

    print()


# ============================================================================
# SECTION 14
# ============================================================================

def section_14_grid_reconstruction() -> None:
    print("=" * 78)
    print("14. EXACT GRID RECONSTRUCTION")
    print("=" * 78)

    tested = 0
    failures = []

    for k in K_VALUES:
        for j in J_VALUES:
            polynomial = sp.Poly(
                get_E(k, j),
                D,
                domain=sp.QQ,
            )

            for d in D_VALUES:
                expected = clean(
                    polynomial.as_expr().subs(
                        D,
                        d,
                    )
                )

                actual = clean(
                    get_E(k, j).subs(
                        D,
                        d,
                    )
                )

                if expected != actual:
                    failures.append(
                        (
                            k,
                            j,
                            d,
                            expected,
                            actual,
                        )
                    )

                tested += 1

    print(
        f"tested={tested}"
    )
    print(
        f"reconstruction failures="
        f"{len(failures)}"
    )

    if failures:
        for failure in failures[:10]:
            print(
                f"  failure={failure}"
            )

        raise RuntimeError(
            "Grid reconstruction failed."
        )

    print()


# ============================================================================
# SECTION 15
# ============================================================================

def section_15_final_diagnostic() -> None:
    print("=" * 78)
    print("15. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "The neighboring-D gcd experiment found only gcd=1 "
        "between adjacent nonzero D-slices."
    )
    print()

    print(
        "The present experiment therefore studies the full "
        "support-removed coefficient matrix."
    )
    print()

    print(
        "For each j:"
    )
    print(
        "    Q_j(k,D) = Support_j(D)^(-1) E_j(k,D)"
    )
    print()

    print(
        "The matrix rows correspond to k values."
    )
    print(
        "The columns correspond to powers of D, X=D-6, "
        "or Y=D-D0."
    )
    print()

    print(
        "Rank 1 would imply exact multiplicative "
        "k-D separability."
    )
    print()

    print(
        "Higher but small rank would indicate a finite "
        "sum of separable pieces."
    )
    print()

    print(
        "The D, X, and Y ranks are compared independently."
    )
    print()

    print(
        "The j=5 case remains isolated because its first "
        "nonzero slice and terminal ladder behave differently "
        "from j=0,...,4."
    )
    print()

    print(
        "This experiment does NOT infer a universal R=5 law."
    )
    print(
        "Any low-rank pattern still requires fresh k,D "
        "data for out-of-sample validation."
    )
    print()

    print(
        "No r=6."
    )
    print(
        "No full pq-kernel expansion."
    )
    print(
        "No replacement universal r,j formula."
    )
    print()


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    section_0_validation()

    section_1_support_quotients()

    section_2_rank_in_D()

    section_3_rank_shifted_X()

    section_4_rank_shifted_boundary()

    section_5_rank_one_test()

    section_6_minor_audit()

    section_7_normalized_profiles()

    section_8_column_rank_growth()

    section_9_k_difference_rank()

    section_10_pairwise_row_proportionality()

    section_11_j5_terminal_rank()

    section_12_rank_summary()

    section_13_support_reconstruction()

    section_14_grid_reconstruction()

    section_15_final_diagnostic()


if __name__ == "__main__":
    main()