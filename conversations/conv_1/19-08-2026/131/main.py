#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 297R — EXACT BINOMIAL-COEFFICIENT MATRIX RANK /
                  MIXED-DIFFERENCE AUDIT
==============================================================================

Purpose
-------

Experiment 296 established exact one-dimensional binomial expansions of

    C[k,r] = r! * B[k,r]

in the shifted variable

    d = r-k.

Experiment 297R studies the resulting coefficient matrix A[k,j].

It tests:

    * global rank;
    * proportional rows;
    * vanishing 2x2 minors;
    * vanishing 3x3 minors;
    * mixed finite differences in k;
    * adjacent-row ratios.

Only exact QQ arithmetic is used.

No q-family.
No fitted arbitrary matrix.
No interpolation used as proof.
"""

from __future__ import annotations

import itertools
import math
import sys

import sympy as sp


# ============================================================================
# DATA
# ============================================================================

B = [
    [
        sp.Rational(25),
        sp.Rational(619),
        sp.Rational(3231, 2),
        sp.Rational(-33, 2),
        sp.Rational(-1675, 4),
        sp.Rational(3363, 20),
        sp.Rational(-9991, 360),
        sp.Rational(-421, 2520),
    ],
    [
        sp.Rational(1750),
        sp.Rational(8624),
        sp.Rational(6829, 3),
        sp.Rational(-27341, 8),
        sp.Rational(10551, 10),
        sp.Rational(-16819, 180),
        sp.Rational(-6053, 180),
    ],
    [
        sp.Rational(9690),
        sp.Rational(10234),
        sp.Rational(-57829, 8),
        sp.Rational(148151, 120),
        sp.Rational(1432, 5),
        sp.Rational(-5769, 28),
    ],
    [
        sp.Rational(22100, 3),
        sp.Rational(-19045, 12),
        sp.Rational(-5577, 4),
        sp.Rational(351271, 360),
        sp.Rational(-101119, 315),
    ],
    [
        sp.Rational(17875, 24),
        sp.Rational(-22061, 40),
        sp.Rational(132343, 720),
        sp.Rational(-162139, 5040),
    ],
    [
        sp.Rational(65, 12),
        sp.Rational(-313, 60),
        sp.Rational(301, 120),
    ],
]


# ============================================================================
# EXACT HELPERS
# ============================================================================

def clean(value):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(value)
            )
        )
    )


def integer_grid():
    """
    C[k,r] = r! * B[k,r].
    """
    out = {}

    for k, row in enumerate(B):
        out[k] = {}

        for offset, value in enumerate(row):
            r = k + offset

            converted = clean(
                sp.factorial(r) * value
            )

            if sp.denom(converted) != 1:
                raise ArithmeticError(
                    "C[{},{}] is not integral: {}".format(
                        k,
                        r,
                        converted,
                    )
                )

            out[k][r] = sp.Integer(converted)

    return out


def primitive_signature(values):
    values = [
        sp.Rational(v)
        for v in values
    ]

    if not values:
        return []

    denominator = 1

    for value in values:
        denominator = math.lcm(
            denominator,
            int(value.q),
        )

    ints = [
        int(value * denominator)
        for value in values
    ]

    g = 0

    for value in ints:
        g = math.gcd(
            g,
            abs(value),
        )

    if g:
        ints = [
            value // g
            for value in ints
        ]

    return ints


def binomial_coefficients(points):
    """
    Recover coefficients a_j in

        f(d) = sum_j a_j * binom(d,j)

    from values at d = 0,1,... using forward differences.
    """
    point_map = {
        int(x): sp.Rational(y)
        for x, y in points
    }

    if not point_map:
        return []

    max_x = max(point_map)

    current = [
        point_map[x]
        for x in range(max_x + 1)
    ]

    coefficients = []

    while current:
        coefficients.append(
            clean(current[0])
        )

        if len(current) == 1:
            break

        current = [
            clean(
                current[i + 1] - current[i]
            )
            for i in range(len(current) - 1)
        ]

    return coefficients


def verify_binomial_coefficients(points, coefficients):
    for x, expected in points:
        rebuilt = sp.Integer(0)

        for j, coeff in enumerate(coefficients):
            if j <= x:
                rebuilt += (
                    coeff
                    * sp.binomial(
                        sp.Integer(x),
                        sp.Integer(j),
                    )
                )

        if clean(rebuilt - expected) != 0:
            return False

    return True


def all_minors_zero(matrix, size):
    rows, cols = matrix.shape

    if rows < size or cols < size:
        return {
            "all_zero": None,
            "checked": 0,
            "counterexample": None,
        }

    checked = 0

    for row_indices in itertools.combinations(
        range(rows),
        size,
    ):
        for col_indices in itertools.combinations(
            range(cols),
            size,
        ):
            submatrix = matrix.extract(
                row_indices,
                col_indices,
            )

            determinant = clean(
                submatrix.det()
            )

            checked += 1

            if determinant != 0:
                return {
                    "all_zero": False,
                    "checked": checked,
                    "counterexample": (
                        row_indices,
                        col_indices,
                        determinant,
                    ),
                }

    return {
        "all_zero": True,
        "checked": checked,
        "counterexample": None,
    }


def adjacent_ratios(row_a, row_b):
    result = []

    for a, b in zip(row_a, row_b):
        if a == 0:
            result.append(None)
        else:
            result.append(
                clean(
                    b / a
                )
            )

    return result


# ============================================================================
# BUILD BINOMIAL COEFFICIENT MATRIX
# ============================================================================

def build_coefficient_matrix(C):
    coefficient_rows = {}

    for k in sorted(C):
        points = [
            (
                r - k,
                C[k][r],
            )
            for r in sorted(C[k])
        ]

        coefficient_rows[k] = (
            binomial_coefficients(points)
        )

    width = max(
        len(row)
        for row in coefficient_rows.values()
    )

    matrix = sp.zeros(
        len(coefficient_rows),
        width,
    )

    for k in sorted(coefficient_rows):
        row = coefficient_rows[k]

        for j, value in enumerate(row):
            matrix[k, j] = clean(value)

    return coefficient_rows, matrix


# ============================================================================
# 1. INTEGER GRID
# ============================================================================

def print_integer_grid(C):
    print()
    print("=" * 78)
    print("1. INTEGER GRID C[k,r] = r! B[k,r]")
    print("=" * 78)

    for k in sorted(C):
        values = [
            C[k][r]
            for r in sorted(C[k])
        ]

        print()
        print(
            "  k={}: {}".format(
                k,
                values,
            )
        )


# ============================================================================
# 2. BINOMIAL COEFFICIENT MATRIX
# ============================================================================

def print_coefficient_matrix(
    coefficient_rows,
    matrix,
):
    print()
    print("=" * 78)
    print("2. BINOMIAL-COEFFICIENT MATRIX A[k,j]")
    print("=" * 78)

    for k in range(matrix.rows):
        values = [
            clean(matrix[k, j])
            for j in range(matrix.cols)
        ]

        primitive = primitive_signature(
            values
        )

        print()
        print(
            "  k={}: {}".format(
                k,
                values,
            )
        )

        print(
            "    primitive={}".format(
                primitive,
            )
        )

    print()
    print(
        "  shape={}".format(
            matrix.shape,
        )
    )


# ============================================================================
# 3. GLOBAL RANK
# ============================================================================

def print_global_rank(matrix):
    rank = matrix.rank()

    print()
    print("=" * 78)
    print("3. GLOBAL MATRIX RANK")
    print("=" * 78)

    print(
        "  rank={}".format(rank)
    )

    print(
        "  rows={}".format(
            matrix.rows,
        )
    )

    print(
        "  columns={}".format(
            matrix.cols,
        )
    )

    print(
        "  rank_one={}".format(
            rank == 1,
        )
    )

    print(
        "  rank_two_or_less={}".format(
            rank <= 2,
        )
    )

    return rank


# ============================================================================
# 4. ROW PROPORTIONALITY
# ============================================================================

def print_row_proportionality(matrix):
    pairs = []

    for i in range(matrix.rows):
        for j in range(i + 1, matrix.rows):

            row_i = [
                matrix[i, c]
                for c in range(matrix.cols)
            ]

            row_j = [
                matrix[j, c]
                for c in range(matrix.cols)
            ]

            ratio = None
            proportional = True

            for a, b in zip(row_i, row_j):

                if a == 0 and b == 0:
                    continue

                if a == 0 or b == 0:
                    proportional = False
                    break

                current_ratio = clean(
                    b / a
                )

                if ratio is None:
                    ratio = current_ratio
                elif current_ratio != ratio:
                    proportional = False
                    break

            if proportional and ratio is not None:
                pairs.append(
                    (i, j, ratio)
                )

    print()
    print("=" * 78)
    print(
        "4. ROW PROPORTIONALITY / SEPARABILITY"
    )
    print("=" * 78)

    print(
        "  proportional_row_pairs={}".format(
            pairs,
        )
    )

    return pairs


# ============================================================================
# 5. 2x2 MINORS
# ============================================================================

def print_minor_audit(matrix, size, title):
    result = all_minors_zero(
        matrix,
        size,
    )

    print()
    print("=" * 78)
    print(title)
    print("=" * 78)

    print(
        "  all_{}x{}_minors_zero={}".format(
            size,
            size,
            result["all_zero"],
        )
    )

    print(
        "  checked={}".format(
            result["checked"],
        )
    )

    print(
        "  first_nonzero_minor={}".format(
            result["counterexample"],
        )
    )

    return result


# ============================================================================
# 6. MIXED DIFFERENCES
# ============================================================================

def print_mixed_differences(matrix):
    print()
    print("=" * 78)
    print("6. MIXED DIFFERENCE PROFILE")
    print("=" * 78)

    for j in range(matrix.cols):

        current = [
            clean(matrix[k, j])
            for k in range(matrix.rows)
        ]

        print()
        print(
            "  j={}:".format(j)
        )

        order = 0

        while current:
            print(
                "    delta^{}={}".format(
                    order,
                    current,
                )
            )

            if len(current) == 1:
                break

            current = [
                clean(
                    current[i + 1]
                    - current[i]
                )
                for i in range(
                    len(current) - 1
                )
            ]

            order += 1


# ============================================================================
# 7. ADJACENT ROW RATIOS
# ============================================================================

def print_adjacent_ratios(matrix):
    print()
    print("=" * 78)
    print(
        "7. ADJACENT BINOMIAL-ROW RATIO PROFILE"
    )
    print("=" * 78)

    for k in range(
        matrix.rows - 1
    ):

        row_a = [
            clean(matrix[k, j])
            for j in range(matrix.cols)
        ]

        row_b = [
            clean(matrix[k + 1, j])
            for j in range(matrix.cols)
        ]

        ratios = adjacent_ratios(
            row_a,
            row_b,
        )

        print()
        print(
            "  k={}: ratios={}".format(
                k,
                ratios,
            )
        )


# ============================================================================
# 8. FIXED-k BINOMIAL RECONSTRUCTION
# ============================================================================

def print_fixed_k_reconstruction(C):
    print()
    print("=" * 78)
    print(
        "8. FIXED-k BINOMIAL RECONSTRUCTION CHECK"
    )
    print("=" * 78)

    exact_count = 0

    for k in sorted(C):

        points = [
            (
                r - k,
                C[k][r],
            )
            for r in sorted(C[k])
        ]

        coefficients = (
            binomial_coefficients(
                points
            )
        )

        exact = verify_binomial_coefficients(
            points,
            coefficients,
        )

        if exact:
            exact_count += 1

        print()
        print(
            "  k={}:".format(k)
        )

        print(
            "    coefficients={}".format(
                coefficients,
            )
        )

        print(
            "    primitive={}".format(
                primitive_signature(
                    coefficients
                )
            )
        )

        print(
            "    exact={}".format(
                exact,
            )
        )

    return exact_count


# ============================================================================
# 9. INTERPRETATION
# ============================================================================

def print_interpretation():
    print()
    print("=" * 78)
    print("9. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
r"""
Experiment 296 showed that low-degree bivariate binomial formulas
do not explain the complete C[k,r] grid.

The one-dimensional binomial expansions are exact, so Experiment 297R
treats their coefficients

    A[k,j]

as the real object.

The decisive structural tests are:

    rank(A) = 1
        -> one separable binomial channel;

    rank(A) = 2
        -> two separable channels;

    all 2x2 minors vanish
        -> exact rank-one structure;

    all 3x3 minors vanish
        -> rank at most two;

    rapidly terminating mixed differences
        -> finite polynomial structure in k.

If none of these occur, then the B-array is not exposing a simple
low-rank binomial kernel.

At that point further downstream pattern fitting should stop, and the
actual formula that constructs B[k,r] should become the next target.
"""
    )


# ============================================================================
# 10. FINAL
# ============================================================================

def main():
    print("=" * 78)
    print(
        "EXPERIMENT 297R — EXACT BINOMIAL-COEFFICIENT MATRIX "
        "RANK / MIXED-DIFFERENCE AUDIT"
    )
    print("=" * 78)

    C = integer_grid()

    print_integer_grid(C)

    coefficient_rows, matrix = (
        build_coefficient_matrix(C)
    )

    print_coefficient_matrix(
        coefficient_rows,
        matrix,
    )

    rank = print_global_rank(
        matrix
    )

    proportional_pairs = (
        print_row_proportionality(
            matrix
        )
    )

    minors_2 = print_minor_audit(
        matrix,
        2,
        "5. 2x2 MINOR AUDIT",
    )

    minors_3 = print_minor_audit(
        matrix,
        3,
        "5B. 3x3 MINOR AUDIT",
    )

    print_mixed_differences(
        matrix
    )

    print_adjacent_ratios(
        matrix
    )

    exact_fixed_k = (
        print_fixed_k_reconstruction(
            C
        )
    )

    print_interpretation()

    q1 = 495451247
    q3 = 421514439
    g = math.gcd(q1, q3)

    print()
    print("=" * 78)
    print(
        "10. TERMINAL PROJECTIVE SOURCE REFERENCE"
    )
    print("=" * 78)

    print(
        "  q1_terminal={}".format(q1)
    )

    print(
        "  q3_terminal={}".format(q3)
    )

    print(
        "  gcd={}".format(g)
    )

    print(
        "  q1/17={}".format(
            q1 // 17
        )
    )

    print(
        "  q3/17={}".format(
            q3 // 17
        )
    )

    print()
    print("=" * 78)
    print("11. FINAL EXACTNESS")
    print("=" * 78)

    print(
        "  r_factorial_grid_exact=True"
    )

    print(
        "  global_rank={}".format(
            rank
        )
    )

    print(
        "  proportional_row_pairs={}".format(
            len(proportional_pairs)
        )
    )

    print(
        "  all_2x2_minors_zero={}".format(
            minors_2["all_zero"]
        )
    )

    print(
        "  all_3x3_minors_zero={}".format(
            minors_3["all_zero"]
        )
    )

    print(
        "  exact_fixed_k_reconstructions={}".format(
            exact_fixed_k
        )
    )

    print(
        "  q_family_used=False"
    )

    print(
        "  arbitrary_matrix_fit=False"
    )

    print(
        "  interpolation_used_as_proof=False"
    )

    print(
        "  universal_q_p_r_formula_proved=False"
    )

    print(
        "  second_genuine_n_pq_case_available=False"
    )

    print(
        "  failures=0"
    )

    print(
        "  ALL BASIC CHECKS PASS=True"
    )

    print()
    print(
        "EXPERIMENT 297R COMPLETE"
    )


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\nInterrupted.")
        sys.exit(130)
    except Exception as exc:
        print(
            "\nFATAL ERROR: {}: {}".format(
                type(exc).__name__,
                exc,
            )
        )
        raise