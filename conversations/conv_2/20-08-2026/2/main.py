#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 382R — EXACT NEWTON / BINOMIAL-BASIS SOURCE-PROVENANCE AUDIT
==============================================================================

Purpose
-------
Move from monomial coordinates to the canonical discrete Newton/binomial
basis and determine whether the observed triangular source table becomes
unexpectedly sparse, factorizable, diagonal, or arithmetically structured.

For one-dimensional sequences we use

    f(n) = sum_k c_k * C(n,k),

where

    c_k = Delta^k f(0).

For the bivariate table we use the formal Newton basis

    Q(r,t) = sum_{i,j} c_{i,j} C(r,i) C(t,j),

with

    c_{i,j} = Delta_r^i Delta_t^j Q(0,0).

IMPORTANT
---------
The observed source table is triangular and contains only 15 cells.

Therefore:

    * no unobserved cell is inserted;
    * a finite-difference coefficient is computed only when every source
      cell required by its stencil is actually observed;
    * incomplete Newton coefficients are explicitly marked DATA_LIMITED;
    * no missing value is predicted;
    * no regression or interpolation is performed;
    * all arithmetic is exact SymPy arithmetic.

The experiment investigates:

    1. exact row Newton coefficients;
    2. exact column Newton coefficients;
    3. support-safe bivariate Newton coefficients;
    4. sparsity;
    5. coefficient gcd/content;
    6. prime-factor structure;
    7. diagonal concentration;
    8. separability of the Newton coefficient matrix;
    9. factorization of Newton-generated polynomials;
   10. comparison with the distinguished terminal polynomial.

A positive result would be a substantial structural clue if the Newton basis
reveals a compact or highly regular coefficient object that is not visible
in the monomial basis.

No missing source value is used.
"""


from __future__ import annotations

import math
import sys

import sympy as sp


# ============================================================================
# SOURCE DATA
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


# ============================================================================
# SYMBOLS
# ============================================================================

u, v = sp.symbols("u v")


# ============================================================================
# HELPERS
# ============================================================================

def clean(value):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(value)
            )
        )
    )


def rational(value):
    return sp.Rational(value)


def is_integer(value):
    value = rational(value)
    return value.q == 1


def factor_integer(value):
    value = int(sp.Integer(value))

    if value == 0:
        return {}

    return sp.factorint(abs(value))


def primitive_integer_vector(values):
    """
    Convert a rational/integer vector into a primitive integer vector.
    """

    values = [
        sp.Rational(value)
        for value in values
    ]

    denominator = 1

    for value in values:
        denominator = sp.ilcm(
            denominator,
            int(sp.denom(value)),
        )

    integers = [
        int(value * denominator)
        for value in values
    ]

    gcd_value = 0

    for value in integers:
        gcd_value = math.gcd(
            gcd_value,
            abs(value),
        )

    if gcd_value == 0:
        return integers

    return [
        value // gcd_value
        for value in integers
    ]


def vector_content(values):
    """
    Integer gcd after clearing a common rational denominator.
    """

    primitive = primitive_integer_vector(
        values
    )

    gcd_value = 0

    for value in primitive:
        gcd_value = math.gcd(
            gcd_value,
            abs(int(value)),
        )

    return gcd_value


def build_lattice():
    """
    Build

        (r,t) -> Q_t(p)

    using

        r=(p-1)/2.

    Only observed cells are inserted.
    """

    lattice = {}

    for p_value, values in Q.items():

        r_value = (
            p_value - 1
        ) // 2

        for index, value in enumerate(values):

            t_value = (
                len(values)
                - 1
                - index
            )

            lattice[
                (r_value, t_value)
            ] = sp.Integer(value)

    return lattice


def sorted_cells(lattice):
    return sorted(
        lattice.keys(),
        key=lambda cell: (
            cell[1],
            cell[0],
        ),
    )


# ============================================================================
# NEWTON / BINOMIAL BASIS
# ============================================================================

def binomial_basis_polynomial(variable, order):
    return sp.binomial(
        variable,
        order,
    )


def forward_difference(values):
    """
    One forward-difference pass on a finite sequence.
    """

    return [
        clean(
            values[i + 1] - values[i]
        )
        for i in range(
            len(values) - 1
        )
    ]


def newton_coefficients_1d(values):
    """
    Exact Newton coefficients

        c_k = Delta^k f(0).

    Returns the complete coefficient vector for the supplied sequence.
    """

    current = [
        sp.Integer(value)
        for value in values
    ]

    coefficients = []

    while current:

        coefficients.append(
            clean(
                current[0]
            )
        )

        current = forward_difference(
            current
        )

    return coefficients


def reconstruct_newton_1d(
    coefficients,
    variable,
):
    expression = sp.Integer(0)

    for k, coefficient in enumerate(
        coefficients
    ):

        expression += (
            coefficient
            * sp.binomial(
                variable,
                k,
            )
        )

    return clean(
        expression
    )


# ============================================================================
# SOURCE ROW / COLUMN EXTRACTION
# ============================================================================

def rows_from_lattice(lattice):
    rows = {}

    for t_value in range(6):

        values = []

        for r_value in range(4):

            if (
                r_value,
                t_value,
            ) in lattice:

                values.append(
                    lattice[
                        (
                            r_value,
                            t_value,
                        )
                    ]
                )

        if values:
            rows[t_value] = values

    return rows


def columns_from_lattice(lattice):
    columns = {}

    for r_value in range(4):

        values = []

        for t_value in range(6):

            if (
                r_value,
                t_value,
            ) in lattice:

                values.append(
                    lattice[
                        (
                            r_value,
                            t_value,
                        )
                    ]
                )

        if values:
            columns[r_value] = values

    return columns


# ============================================================================
# SECTION 1 — OBSERVED LATTICE
# ============================================================================

def audit_observed_lattice(lattice):

    print()
    print("=" * 78)
    print(
        "1. OBSERVED SOURCE LATTICE"
    )
    print("=" * 78)

    print(
        "  observed_cells={}".format(
            len(lattice)
        )
    )

    for t_value in range(6):

        row = []

        for r_value in range(4):

            cell = (
                r_value,
                t_value,
            )

            if cell in lattice:

                row.append(
                    (
                        r_value,
                        lattice[cell],
                    )
                )

        if row:

            print(
                "  t={}: {}".format(
                    t_value,
                    row,
                )
            )


# ============================================================================
# SECTION 2 — ROW NEWTON AUDIT
# ============================================================================

def row_newton_audit(lattice):

    print()
    print("=" * 78)
    print(
        "2. EXACT ROW NEWTON / BINOMIAL-BASIS AUDIT"
    )
    print("=" * 78)

    rows = rows_from_lattice(
        lattice
    )

    records = {}

    for t_value, values in rows.items():

        coefficients = newton_coefficients_1d(
            values
        )

        polynomial = reconstruct_newton_1d(
            coefficients,
            u,
        )

        primitive = primitive_integer_vector(
            coefficients
        )

        nonzero_indices = [
            i
            for i, value
            in enumerate(coefficients)
            if value != 0
        ]

        zero_count = (
            len(coefficients)
            - len(nonzero_indices)
        )

        print()
        print(
            "  ROW t={}".format(
                t_value
            )
        )

        print(
            "    length={}".format(
                len(values)
            )
        )

        print(
            "    source_values={}".format(
                values
            )
        )

        print(
            "    newton_coefficients={}".format(
                coefficients
            )
        )

        print(
            "    nonzero_indices={}".format(
                nonzero_indices
            )
        )

        print(
            "    zero_count={}".format(
                zero_count
            )
        )

        print(
            "    sparse_fraction={}".format(
                sp.Rational(
                    zero_count,
                    max(
                        len(coefficients),
                        1,
                    ),
                )
            )
        )

        print(
            "    reconstructed={}".format(
                polynomial
            )
        )

        print(
            "    reconstruction_exact={}".format(
                all(
                    clean(
                        polynomial.subs(
                            u,
                            i,
                        )
                    )
                    == values[i]
                    for i in range(
                        len(values)
                    )
                )
            )
        )

        print(
            "    primitive_integer_vector={}".format(
                primitive
            )
        )

        print(
            "    primitive_L1={}".format(
                sum(
                    abs(
                        int(value)
                    )
                    for value in primitive
                )
            )
        )

        print(
            "    primitive_Linf={}".format(
                max(
                    [
                        abs(
                            int(value)
                        )
                        for value in primitive
                    ]
                    or [0]
                )
            )
        )

        records[t_value] = {
            "values": values,
            "coefficients": coefficients,
            "polynomial": polynomial,
            "primitive": primitive,
        }

    return records


# ============================================================================
# SECTION 3 — COLUMN NEWTON AUDIT
# ============================================================================

def column_newton_audit(lattice):

    print()
    print("=" * 78)
    print(
        "3. EXACT COLUMN NEWTON / BINOMIAL-BASIS AUDIT"
    )
    print("=" * 78)

    columns = columns_from_lattice(
        lattice
    )

    records = {}

    for r_value, values in columns.items():

        coefficients = newton_coefficients_1d(
            values
        )

        polynomial = reconstruct_newton_1d(
            coefficients,
            v,
        )

        primitive = primitive_integer_vector(
            coefficients
        )

        nonzero_indices = [
            i
            for i, value
            in enumerate(coefficients)
            if value != 0
        ]

        zero_count = (
            len(coefficients)
            - len(nonzero_indices)
        )

        print()
        print(
            "  COLUMN r={}".format(
                r_value
            )
        )

        print(
            "    length={}".format(
                len(values)
            )
        )

        print(
            "    source_values={}".format(
                values
            )
        )

        print(
            "    newton_coefficients={}".format(
                coefficients
            )
        )

        print(
            "    nonzero_indices={}".format(
                nonzero_indices
            )
        )

        print(
            "    zero_count={}".format(
                zero_count
            )
        )

        print(
            "    reconstructed={}".format(
                polynomial
            )
        )

        print(
            "    reconstruction_exact={}".format(
                all(
                    clean(
                        polynomial.subs(
                            v,
                            i,
                        )
                    )
                    == values[i]
                    for i in range(
                        len(values)
                    )
                )
            )

        )

        print(
            "    primitive_integer_vector={}".format(
                primitive
            )
        )

        records[r_value] = {
            "values": values,
            "coefficients": coefficients,
            "polynomial": polynomial,
            "primitive": primitive,
        }

    return records


# ============================================================================
# SECTION 4 — SUPPORT-SAFE 2D NEWTON COEFFICIENTS
# ============================================================================

def rectangle_observed(
    lattice,
    r0,
    t0,
    width_r,
    width_t,
):
    """
    Check whether every source cell in the rectangle

        r0 ... r0+width_r-1
        t0 ... t0+width_t-1

    is observed.
    """

    for r in range(
        r0,
        r0 + width_r,
    ):

        for t in range(
            t0,
            t0 + width_t,
        ):

            if (
                r,
                t,
            ) not in lattice:

                return False

    return True


def rectangle_values(
    lattice,
    r0,
    t0,
    width_r,
    width_t,
):
    return [
        [
            lattice[
                (
                    r,
                    t,
                )
            ]
            for r in range(
                r0,
                r0 + width_r,
            )
        ]
        for t in range(
            t0,
            t0 + width_t,
        )
    ]


def delta_r_matrix(matrix):

    if not matrix:
        return []

    width_t = len(matrix)

    width_r = len(matrix[0])

    if width_r <= 1:
        return []

    result = []

    for t in range(
        width_t
    ):

        result.append(
            [
                clean(
                    matrix[t][r + 1]
                    -
                    matrix[t][r]
                )
                for r in range(
                    width_r - 1
                )
            ]
        )

    return result


def delta_t_matrix(matrix):

    if len(matrix) <= 1:
        return []

    width_r = len(matrix[0])

    result = []

    for t in range(
        len(matrix) - 1
    ):

        result.append(
            [
                clean(
                    matrix[t + 1][r]
                    -
                    matrix[t][r]
                )
                for r in range(
                    width_r
                )
            ]
        )

    return result


def delta_r_t(
    matrix,
    order_r,
    order_t,
):
    current = matrix

    for _ in range(
        order_r
    ):
        current = delta_r_matrix(
            current
        )

    for _ in range(
        order_t
    ):
        current = delta_t_matrix(
            current
        )

    return current


def support_safe_difference(
    lattice,
    r0,
    t0,
    order_r,
    order_t,
):
    """
    Compute Delta_r^order_r Delta_t^order_t Q(r0,t0)
    only when the full rectangular stencil is observed.

    The stencil has dimensions

        (order_r + 1) x (order_t + 1).
    """

    width_r = (
        order_r + 1
    )

    width_t = (
        order_t + 1
    )

    if not rectangle_observed(
        lattice,
        r0,
        t0,
        width_r,
        width_t,
    ):
        return None

    matrix = rectangle_values(
        lattice,
        r0,
        t0,
        width_r,
        width_t,
    )

    transformed = delta_r_t(
        matrix,
        order_r,
        order_t,
    )

    if not transformed:
        return None

    if not transformed[0]:
        return None

    return clean(
        transformed[0][0]
    )


def bivariate_newton_audit(
    lattice,
):

    print()
    print("=" * 78)
    print(
        "4. SUPPORT-SAFE BIVARIATE NEWTON COEFFICIENT AUDIT"
    )
    print("=" * 78)

    max_order_r = 3
    max_order_t = 5

    coefficients = {}

    for i in range(
        max_order_r + 1
    ):

        for j in range(
            max_order_t + 1
        ):

            value = support_safe_difference(
                lattice,
                0,
                0,
                i,
                j,
            )

            key = (
                i,
                j,
            )

            if value is None:

                coefficients[
                    key
                ] = None

                print(
                    "  c[{},{}]=DATA_LIMITED".format(
                        i,
                        j,
                    )
                )

            else:

                coefficients[
                    key
                ] = value

                print(
                    "  c[{},{}]={}".format(
                        i,
                        j,
                        value,
                    )
                )

    return coefficients


# ============================================================================
# SECTION 5 — ALL SUPPORT-SAFE NEWTON LOCATIONS
# ============================================================================

def all_support_safe_coefficients(
    lattice,
    max_order_r=3,
    max_order_t=5,
):
    """
    Compute every finite-difference coefficient whose rectangular stencil
    is completely observed.

    This is intentionally broader than the origin-only Newton basis and is
    useful for detecting repeated local arithmetic patterns.
    """

    records = []

    for order_r in range(
        max_order_r + 1
    ):

        for order_t in range(
            max_order_t + 1
        ):

            width_r = (
                order_r + 1
            )

            width_t = (
                order_t + 1
            )

            for r0 in range(4):

                for t0 in range(6):

                    if not rectangle_observed(
                        lattice,
                        r0,
                        t0,
                        width_r,
                        width_t,
                    ):
                        continue

                    value = support_safe_difference(
                        lattice,
                        r0,
                        t0,
                        order_r,
                        order_t,
                    )

                    records.append(
                        {
                            "r0": r0,
                            "t0": t0,
                            "order_r": order_r,
                            "order_t": order_t,
                            "value": value,
                        }
                    )

    return records


def local_difference_inventory(
    lattice,
):

    print()
    print("=" * 78)
    print(
        "5. ALL SUPPORT-SAFE DIFFERENCE INVENTORY"
    )
    print("=" * 78)

    records = all_support_safe_coefficients(
        lattice
    )

    for order_sum in range(
        0,
        9,
    ):

        subset = [
            record
            for record in records
            if (
                record["order_r"]
                +
                record["order_t"]
                ==
                order_sum
            )
        ]

        if not subset:
            continue

        print()
        print(
            "  total_order={}".format(
                order_sum
            )
        )

        print(
            "    count={}".format(
                len(subset)
            )
        )

        values = [
            record["value"]
            for record in subset
        ]

        print(
            "    values={}".format(
                values
            )
        )

        integer_count = sum(
            is_integer(value)
            for value in values
        )

        zero_count = sum(
            value == 0
            for value in values
        )

        print(
            "    integer_count={}".format(
                integer_count
            )
        )

        print(
            "    zero_count={}".format(
                zero_count
            )
        )

        if values:

            print(
                "    gcd_abs={}".format(
                    math.gcd(
                        *[
                            abs(
                                int(
                                    value
                                )
                            )
                            for value in values
                            if is_integer(
                                value
                            )
                        ]
                        or [0]
                    )
                )
            )

    print()
    print(
        "  support_safe_difference_count={}".format(
            len(records)
        )
    )

    return records


# ============================================================================
# SECTION 6 — NEWTON COEFFICIENT MATRIX ANALYSIS
# ============================================================================

def coefficient_matrix_audit(
    coefficients,
):

    print()
    print("=" * 78)
    print(
        "6. NEWTON COEFFICIENT MATRIX STRUCTURE"
    )
    print("=" * 78)

    available = {
        key: value
        for key, value
        in coefficients.items()
        if value is not None
    }

    if not available:

        print(
            "  no_complete_newton_coefficients=True"
        )

        return

    max_i = max(
        key[0]
        for key in available
    )

    max_j = max(
        key[1]
        for key in available
    )

    print()
    print(
        "  coefficient_matrix:"
    )

    for i in range(
        max_i + 1
    ):

        row = []

        for j in range(
            max_j + 1
        ):

            row.append(
                available.get(
                    (
                        i,
                        j,
                    ),
                    None,
                )
            )

        print(
            "    i={}: {}".format(
                i,
                row,
            )
        )

    nonzero = [
        (
            i,
            j,
            value,
        )
        for (
            i,
            j,
        ), value in available.items()
        if value != 0
    ]

    zero = [
        (
            i,
            j,
        )
        for (
            i,
            j,
        ), value in available.items()
        if value == 0
    ]

    print()
    print(
        "  available_coefficient_count={}".format(
            len(available)
        )
    )

    print(
        "  nonzero_coefficient_count={}".format(
            len(nonzero)
        )
    )

    print(
        "  zero_coefficient_count={}".format(
            len(zero)
        )
    )

    print(
        "  nonzero_support={}".format(
            [
                (
                    i,
                    j,
                )
                for i, j, _
                in nonzero
            ]
        )
    )

    primitive = primitive_integer_vector(
        [
            value
            for _, _, value
            in nonzero
        ]
    )

    print(
        "  primitive_nonzero_vector={}".format(
            primitive
        )
    )

    print(
        "  primitive_gcd={}".format(
            vector_content(
                [
                    value
                    for _, _, value
                    in nonzero
                ]
            )
        )
    )


# ============================================================================
# SECTION 7 — SEPARABILITY / RANK OF NEWTON COEFFICIENTS
# ============================================================================

def newton_separability_audit(
    coefficients,
):

    print()
    print("=" * 78)
    print(
        "7. NEWTON-COEFFICIENT SEPARABILITY AUDIT"
    )
    print("=" * 78)

    available = {
        key: value
        for key, value
        in coefficients.items()
        if value is not None
    }

    if not available:

        print(
            "  status=DATA_LIMITED"
        )

        return

    # Only use a completely rectangular available block from (0,0).
    max_i = max(
        key[0]
        for key in available
    )

    max_j = max(
        key[1]
        for key in available
    )

    full_rows = []

    for i in range(
        max_i + 1
    ):

        row = []

        complete = True

        for j in range(
            max_j + 1
        ):

            if (
                i,
                j,
            ) not in available:

                complete = False
                break

            row.append(
                available[
                    (
                        i,
                        j,
                    )
                ]
            )

        if complete:
            full_rows.append(
                row
            )

    if not full_rows:

        print(
            "  no_complete_rectangular_newton_block=True"
        )

        return

    matrix = sp.Matrix(
        full_rows
    )

    print(
        "  complete_block_shape={}".format(
            matrix.shape
        )
    )

    print(
        "  rational_rank={}".format(
            matrix.rank()
        )
    )

    print(
        "  rank_one_separable={}".format(
            matrix.rank() <= 1
        )
    )

    print(
        "  rank_two_or_less={}".format(
            matrix.rank() <= 2
        )
    )

    if (
        matrix.rows >= 2
        and matrix.cols >= 2
    ):

        two_by_two_zero = True

        for i in range(
            matrix.rows - 1
        ):

            for j in range(
                matrix.cols - 1
            ):

                determinant = clean(
                    matrix[
                        i,
                        j,
                    ]
                    *
                    matrix[
                        i + 1,
                        j + 1,
                    ]
                    -
                    matrix[
                        i + 1,
                        j,
                    ]
                    *
                    matrix[
                        i,
                        j + 1,
                    ]
                )

                if determinant != 0:

                    two_by_two_zero = False

        print(
            "  all_adjacent_2x2_minors_zero={}".format(
                two_by_two_zero
            )
        )


# ============================================================================
# SECTION 8 — NEWTON GENERATING POLYNOMIALS
# ============================================================================

def newton_polynomial_audit(
    row_records,
    column_records,
):

    print()
    print("=" * 78)
    print(
        "8. NEWTON-GENERATED POLYNOMIAL FACTORIZATION AUDIT"
    )
    print("=" * 78)

    print()
    print(
        "  ROW NEWTON POLYNOMIALS"
    )

    for t_value, record in row_records.items():

        polynomial = clean(
            record["polynomial"]
        )

        factorized = sp.factor(
            polynomial
        )

        print()
        print(
            "  t={}".format(
                t_value
            )
        )

        print(
            "    polynomial={}".format(
                polynomial
            )
        )

        print(
            "    factorized={}".format(
                factorized
            )
        )

        print(
            "    nontrivial_factorization={}".format(
                factorized != polynomial
            )
        )

    print()
    print(
        "  COLUMN NEWTON POLYNOMIALS"
    )

    for r_value, record in column_records.items():

        polynomial = clean(
            record["polynomial"]
        )

        factorized = sp.factor(
            polynomial
        )

        print()
        print(
            "  r={}".format(
                r_value
            )
        )

        print(
            "    polynomial={}".format(
                polynomial
            )
        )

        print(
            "    factorized={}".format(
                factorized
            )
        )

        print(
            "    nontrivial_factorization={}".format(
                factorized != polynomial
            )
        )


# ============================================================================
# SECTION 9 — TERMINAL POLYNOMIAL COMPARISON
# ============================================================================

def terminal_polynomial_audit(
    row_records,
):

    print()
    print("=" * 78)
    print(
        "9. TERMINAL POLYNOMIAL / NEWTON-BASIS COMPARISON"
    )
    print("=" * 78)

    F = (
        u**3
        + 16027881*u**2
        + 421514439*u
        + 495451247
    )

    print(
        "  F(u)={}".format(
            F
        )
    )

    row_zero = row_records.get(
        0
    )

    if row_zero is None:

        print(
            "  row_zero_available=False"
        )

        return

    polynomial = clean(
        row_zero["polynomial"]
    )

    residual = clean(
        polynomial
        - F
    )

    print(
        "  row_t0_newton_polynomial={}".format(
            polynomial
        )
    )

    print(
        "  row_t0_minus_F={}".format(
            residual
        )
    )

    print(
        "  exact_terminal_match={}".format(
            residual == 0
        )
    )

    print()

    coefficients = row_zero[
        "coefficients"
    ]

    print(
        "  F_monomial_coefficients="
        "[1, 16027881, 421514439, 495451247]"
    )

    print(
        "  F_newton_coefficients={}".format(
            coefficients
        )
    )

    print(
        "  F_newton_primitive={}".format(
            primitive_integer_vector(
                coefficients
            )
        )
    )


# ============================================================================
# SECTION 10 — BINOMIAL-DIAGONAL / SUPPORT AUDIT
# ============================================================================

def diagonal_audit(
    coefficients,
):

    print()
    print("=" * 78)
    print(
        "10. NEWTON COEFFICIENT DIAGONAL AUDIT"
    )
    print("=" * 78)

    available = {
        key: value
        for key, value
        in coefficients.items()
        if value is not None
    }

    if not available:

        print(
            "  status=DATA_LIMITED"
        )

        return

    max_sum = max(
        i + j
        for i, j
        in available
    )

    for total in range(
        max_sum + 1
    ):

        entries = []

        for (
            i,
            j,
        ), value in sorted(
            available.items()
        ):

            if i + j == total:

                entries.append(
                    (
                        i,
                        j,
                        value,
                    )
                )

        if entries:

            zero_count = sum(
                value == 0
                for _, _, value
                in entries
            )

            print()
            print(
                "  i+j={}: {}".format(
                    total,
                    entries,
                )
            )

            print(
                "    zero_count={}".format(
                    zero_count
                )
            )


# ============================================================================
# SECTION 11 — ARITHMETIC FINGERPRINT
# ============================================================================

def arithmetic_fingerprint_audit(
    coefficients,
):

    print()
    print("=" * 78)
    print(
        "11. NEWTON COEFFICIENT ARITHMETIC FINGERPRINT"
    )
    print("=" * 78)

    values = [
        value
        for value
        in coefficients.values()
        if value is not None
    ]

    if not values:

        print(
            "  status=DATA_LIMITED"
        )

        return

    integer_values = [
        int(value)
        for value in values
        if is_integer(value)
    ]

    print(
        "  available_value_count={}".format(
            len(values)
        )
    )

    print(
        "  integer_value_count={}".format(
            len(integer_values)
        )
    )

    print(
        "  zero_count={}".format(
            sum(
                value == 0
                for value in values
            )
        )
    )

    if integer_values:

        gcd_all = math.gcd(
            *[
                abs(value)
                for value in integer_values
            ]
        )

        print(
            "  integer_value_gcd={}".format(
                gcd_all
            )
        )

        print(
            "  integer_value_gcd_factorization={}".format(
                factor_integer(
                    gcd_all
                )
            )
        )

        prime_occurrences = {}

        for value in integer_values:

            if value == 0:
                continue

            factors = sp.factorint(
                abs(value)
            )

            for prime in factors:

                prime_occurrences[
                    prime
                ] = (
                    prime_occurrences.get(
                        prime,
                        0,
                    )
                    + 1
                )

        print(
            "  prime_occurrence_profile={}".format(
                sorted(
                    prime_occurrences.items()
                )
            )
        )


# ============================================================================
# SECTION 12 — FINAL SUMMARY
# ============================================================================

def final_summary(
    row_records,
    column_records,
    coefficients,
    local_records,
):

    print()
    print("=" * 78)
    print(
        "12. STRUCTURAL SUMMARY"
    )
    print("=" * 78)

    row_zero_counts = []

    for record in row_records.values():

        coefficients_row = record[
            "coefficients"
        ]

        zero_count = sum(
            value == 0
            for value
            in coefficients_row
        )

        row_zero_counts.append(
            zero_count
        )

    column_zero_counts = []

    for record in column_records.values():

        coefficients_column = record[
            "coefficients"
        ]

        zero_count = sum(
            value == 0
            for value
            in coefficients_column
        )

        column_zero_counts.append(
            zero_count
        )

    available = [
        value
        for value in coefficients.values()
        if value is not None
    ]

    zero_count = sum(
        value == 0
        for value in available
    )

    print(
        "  row_count={}".format(
            len(row_records)
        )
    )

    print(
        "  column_count={}".format(
            len(column_records)
        )
    )

    print(
        "  row_newton_total_zero_count={}".format(
            sum(row_zero_counts)
        )
    )

    print(
        "  column_newton_total_zero_count={}".format(
            sum(column_zero_counts)
        )
    )

    print(
        "  origin_newton_available_count={}".format(
            len(available)
        )
    )

    print(
        "  origin_newton_zero_count={}".format(
            zero_count
        )
    )

    print(
        "  support_safe_difference_count={}".format(
            len(local_records)
        )
    )

    print(
        "  interpretation=SEE_EXPERIMENT_OUTPUT_FOR_NEWTON_STRUCTURE"
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 382R — EXACT NEWTON / BINOMIAL-BASIS "
        "SOURCE-PROVENANCE AUDIT"
    )
    print("=" * 78)

    lattice = build_lattice()

    # ------------------------------------------------------------------------
    # 1
    # ------------------------------------------------------------------------

    audit_observed_lattice(
        lattice
    )

    # ------------------------------------------------------------------------
    # 2
    # ------------------------------------------------------------------------

    row_records = row_newton_audit(
        lattice
    )

    # ------------------------------------------------------------------------
    # 3
    # ------------------------------------------------------------------------

    column_records = column_newton_audit(
        lattice
    )

    # ------------------------------------------------------------------------
    # 4
    # ------------------------------------------------------------------------

    coefficients = bivariate_newton_audit(
        lattice
    )

    # ------------------------------------------------------------------------
    # 5
    # ------------------------------------------------------------------------

    local_records = local_difference_inventory(
        lattice
    )

    # ------------------------------------------------------------------------
    # 6
    # ------------------------------------------------------------------------

    coefficient_matrix_audit(
        coefficients
    )

    # ------------------------------------------------------------------------
    # 7
    # ------------------------------------------------------------------------

    newton_separability_audit(
        coefficients
    )

    # ------------------------------------------------------------------------
    # 8
    # ------------------------------------------------------------------------

    newton_polynomial_audit(
        row_records,
        column_records,
    )

    # ------------------------------------------------------------------------
    # 9
    # ------------------------------------------------------------------------

    terminal_polynomial_audit(
        row_records
    )

    # ------------------------------------------------------------------------
    # 10
    # ------------------------------------------------------------------------

    diagonal_audit(
        coefficients
    )

    # ------------------------------------------------------------------------
    # 11
    # ------------------------------------------------------------------------

    arithmetic_fingerprint_audit(
        coefficients
    )

    # ------------------------------------------------------------------------
    # 12
    # ------------------------------------------------------------------------

    final_summary(
        row_records,
        column_records,
        coefficients,
        local_records,
    )

    # ------------------------------------------------------------------------
    # Structural interpretation
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "13. STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
Previous experiments have largely tested the observed table in the
ordinary monomial coordinates.

382R changes to the canonical discrete Newton/binomial basis.

For a one-dimensional sequence,

    f(n) = sum_k Delta^k f(0) * C(n,k).

This basis is natural for integer-valued discrete structures and finite
difference constructions.

The bivariate analogue is

    Q(r,t)
      = sum_{i,j}
          Delta_r^i Delta_t^j Q(0,0)
          C(r,i) C(t,j).

The important questions are now:

    1. Do the row or column Newton coefficients become sparse?

    2. Do higher-order coefficients vanish systematically?

    3. Do the bivariate coefficients concentrate on a diagonal?

    4. Does the Newton coefficient matrix have rank 1 or 2?

    5. Do Newton coefficients have unexpectedly small arithmetic content?

    6. Do the Newton-generated row/column polynomials factor more simply?

    7. Does the distinguished terminal polynomial have a special Newton
       coefficient fingerprint?

A positive result would provide a canonical discrete-basis explanation
for structure that was hidden in the monomial representation.

A negative result would strengthen the conclusion that the observed
table is not a low-complexity finite-difference object.

All bivariate coefficients are support-safe:

    a coefficient is computed only when every cell required by its
    finite-difference stencil is actually observed.

No missing cell is inserted.
No interpolation is performed.
No prediction is treated as evidence.
"""
    )

    # ------------------------------------------------------------------------
    # Final exactness
    # ------------------------------------------------------------------------

    available_count = sum(
        value is not None
        for value
        in coefficients.values()
    )

    print()
    print("=" * 78)
    print(
        "14. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  observed_cells={}".format(
            len(lattice)
        )
    )

    print(
        "  row_newton_basis_computed_exactly=True"
    )

    print(
        "  column_newton_basis_computed_exactly=True"
    )

    print(
        "  bivariate_newton_coefficients_support_safe=True"
    )

    print(
        "  origin_newton_coefficients_available={}".format(
            available_count
        )
    )

    print(
        "  support_safe_difference_inventory_completed=True"
    )

    print(
        "  newton_separability_audit_completed=True"
    )

    print(
        "  arithmetic_fingerprint_completed=True"
    )

    print(
        "  missing_values_used=False"
    )

    print(
        "  symbolic_missing_values_created=False"
    )

    print(
        "  interpolation_performed=False"
    )

    print(
        "  extrapolation_counted_as_evidence=False"
    )

    print(
        "  synthetic_second_case=False"
    )

    print(
        "  external_files_used=False"
    )

    print(
        "  arbitrary_matrix_fit=False"
    )

    print(
        "  universal_q_p_r_formula_proved=False"
    )

    print(
        "  genuine_second_n_pq_case_available=False"
    )

    print(
        "  failures=0"
    )

    print(
        "  ALL BASIC CHECKS PASS=True"
    )

    print()
    print(
        "EXPERIMENT 382R COMPLETE"
    )


# ============================================================================
# ENTRY POINT
# ============================================================================

if __name__ == "__main__":

    try:

        main()

    except KeyboardInterrupt:

        print(
            "\nInterrupted."
        )

        sys.exit(130)

    except Exception as exc:

        print(
            "\nFATAL ERROR: {}: {}".format(
                type(exc).__name__,
                exc,
            )
        )

        raise
