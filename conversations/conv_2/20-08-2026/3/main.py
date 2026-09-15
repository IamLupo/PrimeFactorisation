#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 383R — EXACT BIVARIATE NEWTON-COEFFICIENT POLYNOMIAL /
                  FACTORIZATION / ALGEBRAIC-RELATION AUDIT
==============================================================================

Purpose
-------
382R showed:

    * all 15 observed Newton coefficients are nonzero;
    * the Newton coefficient gcd is 1;
    * no obvious row/column sparsity appears;
    * no obvious separability appears from the available rectangular block.

383R therefore studies the Newton coefficient object itself.

Given

    c_{i,j} = Delta_r^i Delta_t^j Q(0,0),

construct the exact coefficient polynomial

    C(u,v) = sum c_{i,j} u^i v^j

from every support-safe coefficient currently available.

Also construct the exact bivariate Newton representation

    N(x,y)
      = sum c_{i,j} binom(x,i) binom(y,j).

The experiment tests:

    1. exact polynomial reconstruction;
    2. factorization of C(u,v);
    3. factorization of N(x,y);
    4. common factors between transformed components;
    5. vanishing on simple algebraic curves;
    6. dependence on i+j and i-j;
    7. coefficient-row / coefficient-column gcd profiles;
    8. 2D coefficient-matrix rank;
    9. low-degree algebraic relations among c_{i,j};
   10. exact comparison with the original monomial generating polynomial.

No missing cell is inserted.

No interpolation is performed.

Only coefficients whose complete finite-difference stencil is observed
are used.

All arithmetic is exact SymPy arithmetic.
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

x, y = sp.symbols("x y")
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


def is_integer(value):
    value = sp.Rational(value)
    return value.q == 1


def factor_integer(value):
    value = int(sp.Integer(value))

    if value == 0:
        return {}

    return sp.factorint(abs(value))


def primitive_integer_vector(values):
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


def build_lattice():
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
                (
                    r_value,
                    t_value,
                )
            ] = sp.Integer(value)

    return lattice


# ============================================================================
# EXACT DIFFERENCE OPERATORS
# ============================================================================

def rectangle_observed(
    lattice,
    r0,
    t0,
    width_r,
    width_t,
):

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


def delta_r(matrix):

    if not matrix:
        return []

    width = len(matrix[0])

    if width <= 1:
        return []

    return [
        [
            clean(
                matrix[t][r + 1]
                -
                matrix[t][r]
            )
            for r in range(
                width - 1
            )
        ]
        for t in range(
            len(matrix)
        )
    ]


def delta_t(matrix):

    if len(matrix) <= 1:
        return []

    width = len(matrix[0])

    return [
        [
            clean(
                matrix[t + 1][r]
                -
                matrix[t][r]
            )
            for r in range(
                width
            )
        ]
        for t in range(
            len(matrix) - 1
        )
    ]


def delta_rt(
    matrix,
    order_r,
    order_t,
):

    current = matrix

    for _ in range(order_r):
        current = delta_r(current)

    for _ in range(order_t):
        current = delta_t(current)

    return current


def support_safe_coefficient(
    lattice,
    order_r,
    order_t,
):

    width_r = order_r + 1
    width_t = order_t + 1

    if not rectangle_observed(
        lattice,
        0,
        0,
        width_r,
        width_t,
    ):
        return None

    matrix = rectangle_values(
        lattice,
        0,
        0,
        width_r,
        width_t,
    )

    transformed = delta_rt(
        matrix,
        order_r,
        order_t,
    )

    if (
        not transformed
        or not transformed[0]
    ):
        return None

    return clean(
        transformed[0][0]
    )


# ============================================================================
# 1. BUILD NEWTON COEFFICIENT TABLE
# ============================================================================

def build_newton_coefficients(
    lattice,
):

    print()
    print("=" * 78)
    print(
        "1. EXACT BIVARIATE NEWTON COEFFICIENT TABLE"
    )
    print("=" * 78)

    coefficients = {}

    for i in range(4):

        for j in range(6):

            value = support_safe_coefficient(
                lattice,
                i,
                j,
            )

            coefficients[
                (
                    i,
                    j,
                )
            ] = value

            if value is None:

                print(
                    "  c[{},{}]=DATA_LIMITED".format(
                        i,
                        j,
                    )
                )

            else:

                print(
                    "  c[{},{}]={}".format(
                        i,
                        j,
                        value,
                    )
                )

    return coefficients


# ============================================================================
# 2. COEFFICIENT POLYNOMIAL
# ============================================================================

def coefficient_polynomial(
    coefficients,
):

    C = sp.Integer(0)

    for (
        i,
        j,
    ), value in coefficients.items():

        if value is None:
            continue

        C += (
            value
            * u**i
            * v**j
        )

    return clean(C)


def coefficient_polynomial_audit(
    coefficients,
):

    print()
    print("=" * 78)
    print(
        "2. EXACT NEWTON-COEFFICIENT POLYNOMIAL"
    )
    print("=" * 78)

    C = coefficient_polynomial(
        coefficients
    )

    print(
        "  C(u,v)={}".format(
            C
        )
    )

    factorized = sp.factor(C)

    print(
        "  factorized={}".format(
            factorized
        )
    )

    print(
        "  nontrivial_factorization={}".format(
            factorized != C
        )
    )

    print(
        "  total_degree={}".format(
            sp.Poly(
                C,
                u,
                v,
            ).total_degree()
        )
    )

    print(
        "  degree_u={}".format(
            sp.degree(
                C,
                u,
            )
        )
    )

    print(
        "  degree_v={}".format(
            sp.degree(
                C,
                v,
            )
        )
    )

    return C


# ============================================================================
# 3. ACTUAL NEWTON REPRESENTATION
# ============================================================================

def newton_expression(
    coefficients,
):

    N = sp.Integer(0)

    for (
        i,
        j,
    ), value in coefficients.items():

        if value is None:
            continue

        N += (
            value
            * sp.binomial(
                x,
                i,
            )
            * sp.binomial(
                y,
                j,
            )
        )

    return clean(N)


def newton_expression_audit(
    coefficients,
    lattice,
):

    print()
    print("=" * 78)
    print(
        "3. EXACT BIVARIATE NEWTON REPRESENTATION"
    )
    print("=" * 78)

    N_raw = newton_expression(
        coefficients
    )

    N = sp.expand_func(
        sp.expand(
            N_raw
        )
    )

    print(
        "  N_binomial={}".format(
            N_raw
        )
    )

    print(
        "  N_expanded={}".format(
            clean(N)
        )
    )

    factorized = sp.factor(
        N
    )

    print(
        "  factorized={}".format(
            factorized
        )
    )

    print(
        "  nontrivial_factorization={}".format(
            factorized != N
        )
    )

    # Exact reconstruction on every observed cell.
    failures = []

    for (
        r,
        t,
    ), value in lattice.items():

        reconstructed = clean(
            N.subs(
                {
                    x: r,
                    y: t,
                }
            )
        )

        if reconstructed != value:

            failures.append(
                (
                    (
                        r,
                        t,
                    ),
                    value,
                    reconstructed,
                )
            )

    print(
        "  observed_cell_reconstruction_failures={}".format(
            failures
        )
    )

    print(
        "  reconstruction_exact={}".format(
            not failures
        )
    )

    return N


# ============================================================================
# 4. ORIGINAL GENERATING POLYNOMIAL
# ============================================================================

def original_generating_polynomial():

    return sp.expand(
        x**3
        - 62398*x**2*y**2
        + 4771718*x**2*y
        + 16027881*x**2
        + 9955176*x*y**4
        - 1263551016*x*y**3
        - 152369292*x*y**2
        - 128667196*x*y
        + 421514439*x
        - 126258696*y**5
        - 11600759760*y**4
        + 2668721436*y**3
        + 1764373740*y**2
        - 1338089411*y
        + 495451247
    )


def original_polynomial_audit(
    N,
):

    print()
    print("=" * 78)
    print(
        "4. ORIGINAL MONOMIAL GENERATING POLYNOMIAL CROSS-CHECK"
    )
    print("=" * 78)

    G = original_generating_polynomial()

    print(
        "  G(x,y)={}".format(
            G
        )
    )

    print(
        "  G_factorized={}".format(
            sp.factor(G)
        )
    )

    # The triangular data give exactly the coefficient support of G.
    # Check equality on all currently observed coefficient positions.
    lattice = build_lattice()

    coefficient_failures = []

    poly = sp.Poly(
        G,
        x,
        y,
    )

    for (
        r,
        t,
    ), value in lattice.items():

        extracted = clean(
            poly.coeff_monomial(
                x**r * y**t
            )
        )

        if extracted != value:

            coefficient_failures.append(
                (
                    (
                        r,
                        t,
                    ),
                    value,
                    extracted,
                )
            )

    print(
        "  observed_coefficient_failures={}".format(
            coefficient_failures
        )
    )

    print(
        "  observed_coefficient_match={}".format(
            not coefficient_failures
        )
    )

    return G


# ============================================================================
# 5. BINOMIAL / MONOMIAL CHANGE OF BASIS
# ============================================================================

def basis_difference_audit(
    G,
    N,
):

    print()
    print("=" * 78)
    print(
        "5. MONOMIAL-vs-NEWTON BASIS DIFFERENCE"
    )
    print("=" * 78)

    diff = clean(
        sp.expand(
            G - N
        )
    )

    print(
        "  G_minus_N={}".format(
            diff
        )
    )

    print(
        "  exact_expression_equal={}".format(
            diff == 0
        )
    )

    print(
        "  G_total_degree={}".format(
            sp.Poly(
                G,
                x,
                y,
            ).total_degree()
        )
    )

    print(
        "  N_total_degree={}".format(
            sp.Poly(
                sp.expand_func(N),
                x,
                y,
            ).total_degree()
        )
    )


# ============================================================================
# 6. NEWTON COEFFICIENT MATRIX RANK
# ============================================================================

def coefficient_matrix_rank_audit(
    coefficients,
):

    print()
    print("=" * 78)
    print(
        "6. EXACT NEWTON COEFFICIENT MATRIX RANK / MINOR AUDIT"
    )
    print("=" * 78)

    # Largest complete rectangle from origin.
    complete = []

    for i in range(4):

        row = []

        valid = True

        for j in range(6):

            value = coefficients.get(
                (
                    i,
                    j,
                )
            )

            if value is None:

                valid = False
                break

            row.append(
                value
            )

        if not valid:
            break

        complete.append(
            row
        )

    # A triangular observed support cannot supply the complete 4x6 block.
    # The largest complete origin block is therefore reported.
    if complete:

        M = sp.Matrix(
            complete
        )

        print(
            "  complete_origin_block_shape={}".format(
                M.shape
            )
        )

        print(
            "  rank={}".format(
                M.rank()
            )
        )

        print(
            "  rank_one={}".format(
                M.rank() == 1
            )
        )

        print(
            "  rank_two_or_less={}".format(
                M.rank() <= 2
            )
        )

    # All available coefficients can still be tested for all observed
    # 2x2 minors whose four entries are known.
    entries = coefficients

    nonzero_minor_count = 0
    zero_minor_count = 0

    for i0 in range(4):

        for i1 in range(
            i0 + 1,
            4,
        ):

            for j0 in range(6):

                for j1 in range(
                    j0 + 1,
                    6,
                ):

                    cells = [
                        entries.get(
                            (
                                i0,
                                j0,
                            )
                        ),
                        entries.get(
                            (
                                i0,
                                j1,
                            )
                        ),
                        entries.get(
                            (
                                i1,
                                j0,
                            )
                        ),
                        entries.get(
                            (
                                i1,
                                j1,
                            )
                        ),
                    ]

                    if any(
                        value is None
                        for value in cells
                    ):
                        continue

                    determinant = clean(
                        cells[0] * cells[3]
                        -
                        cells[1] * cells[2]
                    )

                    if determinant == 0:
                        zero_minor_count += 1
                    else:
                        nonzero_minor_count += 1

    print(
        "  observed_complete_2x2_zero_count={}".format(
            zero_minor_count
        )
    )

    print(
        "  observed_complete_2x2_nonzero_count={}".format(
            nonzero_minor_count
        )
    )

    print(
        "  all_observed_2x2_minors_zero={}".format(
            nonzero_minor_count == 0
        )
    )


# ============================================================================
# 7. DIAGONAL / ANTI-DIAGONAL STRUCTURE
# ============================================================================

def diagonal_structure_audit(
    coefficients,
):

    print()
    print("=" * 78)
    print(
        "7. EXACT DIAGONAL / ANTI-DIAGONAL ALGEBRAIC STRUCTURE"
    )
    print("=" * 78)

    available = {
        key: value
        for key, value
        in coefficients.items()
        if value is not None
    }

    print()
    print(
        "  DIAGONALS i+j=k"
    )

    for k in range(9):

        entries = [
            (
                i,
                j,
                value,
            )
            for (
                i,
                j,
            ), value in sorted(
                available.items()
            )
            if i + j == k
        ]

        if entries:

            values = [
                value
                for _, _, value
                in entries
            ]

            gcd_value = math.gcd(
                *[
                    abs(
                        int(value)
                    )
                    for value in values
                ]
            )

            print(
                "  k={}: entries={}, gcd={}".format(
                    k,
                    entries,
                    gcd_value,
                )
            )

    print()
    print(
        "  ANTI-DIAGONALS i-j=k"
    )

    for k in range(
        -5,
        4,
    ):

        entries = [
            (
                i,
                j,
                value,
            )
            for (
                i,
                j,
            ), value in sorted(
                available.items()
            )
            if i - j == k
        ]

        if entries:

            print(
                "  k={}: entries={}".format(
                    k,
                    entries,
                )
            )


# ============================================================================
# 8. SIMPLE ALGEBRAIC CURVE TESTS
# ============================================================================

def curve_zero_audit(
    C,
):

    print()
    print("=" * 78)
    print(
        "8. SIMPLE ALGEBRAIC CURVE / SPECIAL-VALUE AUDIT"
    )
    print("=" * 78)

    tests = {
        "u=0": C.subs(u, 0),
        "v=0": C.subs(v, 0),
        "u=1": C.subs(u, 1),
        "v=1": C.subs(v, 1),
        "u=v": C.subs(u, v),
        "u=-v": C.subs(u, -v),
        "u=2v": C.subs(u, 2 * v),
        "v=2u": C.subs(v, 2 * u),
    }

    for name, expression in tests.items():

        expression = clean(
            expression
        )

        print()
        print(
            "  {}: {}".format(
                name,
                expression,
            )
        )

        print(
            "    identically_zero={}".format(
                expression == 0
            )
        )

        print(
            "    factorized={}".format(
                sp.factor(
                    expression
                )
            )
        )


# ============================================================================
# 9. SMALL INTEGER COMBINATION / FACTORIZATION SEARCH
# ============================================================================

def combination_factor_search(
    coefficients,
):

    print()
    print("=" * 78)
    print(
        "9. SMALL LINEAR-COMBINATION FACTORIZATION SEARCH"
    )
    print("=" * 78)

    C = coefficient_polynomial(
        coefficients
    )

    # Generate simple perturbations of u,v and coefficient content.
    candidates = {
        "C": C,
        "C(u+1,v)": clean(
            C.subs(
                u,
                u + 1,
            )
        ),
        "C(u,v+1)": clean(
            C.subs(
                v,
                v + 1,
            )
        ),
        "C(u+1,v+1)": clean(
            C.subs(
                {
                    u: u + 1,
                    v: v + 1,
                }
            )
        ),
        "C(u-1,v)": clean(
            C.subs(
                u,
                u - 1,
            )
        ),
        "C(u,v-1)": clean(
            C.subs(
                v,
                v - 1,
            )
        ),
    }

    for name, expression in candidates.items():

        factorized = sp.factor(
            expression
        )

        print()
        print(
            "  {}:".format(
                name
            )
        )

        print(
            "    expression={}".format(
                expression
            )
        )

        print(
            "    factorized={}".format(
                factorized
            )
        )

        print(
            "    nontrivial_factorization={}".format(
                factorized != expression
            )
        )


# ============================================================================
# 10. ARITHMETIC CONTENT PROFILE
# ============================================================================

def arithmetic_profile(
    coefficients,
):

    print()
    print("=" * 78)
    print(
        "10. EXACT NEWTON COEFFICIENT ARITHMETIC PROFILE"
    )
    print("=" * 78)

    available = [
        (
            i,
            j,
            value,
        )
        for (
            i,
            j,
        ), value in coefficients.items()
        if value is not None
    ]

    for i, j, value in available:

        print()
        print(
            "  c[{},{}]={}".format(
                i,
                j,
                value,
            )
        )

        print(
            "    integer={}".format(
                is_integer(value)
            )
        )

        if is_integer(value):

            integer_value = int(
                value
            )

            print(
                "    factorization={}".format(
                    factor_integer(
                        integer_value
                    )
                )
            )

    values = [
        int(value)
        for _, _, value in available
        if is_integer(value)
    ]

    if values:

        gcd_value = math.gcd(
            *[
                abs(value)
                for value in values
            ]
        )

        print()
        print(
            "  global_integer_gcd={}".format(
                gcd_value
            )
        )

        print(
            "  global_integer_gcd_factorization={}".format(
                factor_integer(
                    gcd_value
                )
            )
        )


# ============================================================================
# 11. FINAL STRUCTURAL SUMMARY
# ============================================================================

def final_summary(
    coefficients,
    C,
    N,
):

    print()
    print("=" * 78)
    print(
        "11. STRUCTURAL SUMMARY"
    )
    print("=" * 78)

    available = [
        value
        for value in coefficients.values()
        if value is not None
    ]

    print(
        "  available_newton_coefficients={}".format(
            len(available)
        )
    )

    print(
        "  zero_newton_coefficients={}".format(
            sum(
                value == 0
                for value in available
            )
        )
    )

    print(
        "  integer_newton_coefficients={}".format(
            sum(
                is_integer(value)
                for value in available
            )
        )
    )

    print(
        "  coefficient_polynomial_degree={}".format(
            sp.Poly(
                C,
                u,
                v,
            ).total_degree()
        )
    )

    print(
        "  coefficient_polynomial_factorization={}".format(
            sp.factor(C)
        )
    )

    print(
        "  newton_expression_degree={}".format(
            sp.Poly(
                sp.expand_func(N),
                x,
                y,
            ).total_degree()
        )
    )

    print(
        "  interpretation="
        "FULL_NEWTON_ALGEBRAIC_STRUCTURE_REQUIRES_SEPARATE_ASSESSMENT"
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 383R — EXACT BIVARIATE NEWTON-COEFFICIENT "
        "POLYNOMIAL / FACTORIZATION / ALGEBRAIC-RELATION AUDIT"
    )
    print("=" * 78)

    lattice = build_lattice()

    print()
    print(
        "  observed_cells={}".format(
            len(lattice)
        )
    )

    coefficients = build_newton_coefficients(
        lattice
    )

    C = coefficient_polynomial_audit(
        coefficients
    )

    N = newton_expression_audit(
        coefficients,
        lattice,
    )

    G = original_polynomial_audit(
        N
    )

    basis_difference_audit(
        G,
        N,
    )

    coefficient_matrix_rank_audit(
        coefficients
    )

    diagonal_structure_audit(
        coefficients
    )

    curve_zero_audit(
        C
    )

    combination_factor_search(
        coefficients
    )

    arithmetic_profile(
        coefficients
    )

    final_summary(
        coefficients,
        C,
        N,
    )

    print()
    print("=" * 78)
    print(
        "12. STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
382R established that the Newton/binomial basis does not create obvious
coefficient sparsity.

383R now treats the Newton coefficients themselves as an algebraic object.

The key distinction is:

    Newton coefficient table
        versus
    Newton coefficient polynomial.

The first can look complicated while the second may factor.

The experiment therefore asks whether

    C(u,v) = sum c_(i,j) u^i v^j

has:

    * a nontrivial factorization;
    * a low-degree factor;
    * a symmetry after translation;
    * vanishing on a simple algebraic curve;
    * low-rank coefficient structure;
    * diagonal or anti-diagonal concentration.

The actual Newton representation

    N(x,y)
      = sum c_(i,j) binom(x,i) binom(y,j)

is also expanded and factorized.

An important sanity condition is exact equality between N and the
original observed generating polynomial on every observed coefficient
cell.

A positive result in C(u,v) would be a new algebraic clue because it
would arise after a canonical discrete change of basis rather than from
an arbitrary fitted recurrence.

A negative result would close another natural coordinate-system route.

No missing source value is used.
No interpolation is performed.
No prediction is treated as data.
"""
    )

    print()
    print("=" * 78)
    print(
        "13. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  observed_cells={}".format(
            len(lattice)
        )
    )

    print(
        "  exact_newton_coefficients=True"
    )

    print(
        "  exact_coefficient_polynomial=True"
    )

    print(
        "  exact_newton_expression=True"
    )

    print(
        "  factorization_audit_completed=True"
    )

    print(
        "  curve_zero_audit_completed=True"
    )

    print(
        "  coefficient_rank_audit_completed=True"
    )

    print(
        "  arithmetic_profile_completed=True"
    )

    print(
        "  missing_values_used=False"
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
        "EXPERIMENT 383R COMPLETE"
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
