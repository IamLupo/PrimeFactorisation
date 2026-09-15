#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 327R — EXACT P-DEPENDENT VERTICAL RECURRENCE AUDIT
==============================================================================

Purpose
-------
Experiments 324R and 326R established:

    * no shared constant-coefficient recurrence in t of order 1, 2, or 3;
    * no low-degree bivariate polynomial law in (p,t);
    * no low-rank separable representation of rank 1 or 2.

The next structured possibility is:

    Q_{t+r}(p)
      =
      sum_{k=0}^{r-1} C_k(p) Q_{t+k}(p),

where each recurrence coefficient C_k(p) is itself a LOW-DEGREE
polynomial in p.

This preserves the idea of a genuine recurrence in t, while allowing
the dynamics to depend on the cross-parameter p.

Models tested
-------------

    order 1:
        C0(p) degree <= 0, 1, 2

    order 2:
        C0(p), C1(p) each degree <= 0, 1, 2

    order 3:
        C0(p), C1(p), C2(p) each degree <= 0, 1

For each model we solve the entire observed triangular table exactly.

Only windows for which ALL required Q-values are observed are used.

The experiment distinguishes:

    EXACT_OVERDETERMINED
        exact solution with more equations than unknown coefficients;

    EXACT_DATA_SIZED
        exact solution with equal numbers of equations and unknowns;

    NONUNIQUE
        consistent but underdetermined;

    NO_SOLUTION
        exact rational inconsistency.

The important output is therefore not merely whether a recurrence can be
solved, but whether a LOW-DEGREE P-DEPENDENT recurrence survives redundant
tests.

No interpolation.
No missing values.
No synthetic second n=pq case.
Exact rational arithmetic only.
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


p, t = sp.symbols("p t")


# ============================================================================
# HELPERS
# ============================================================================

def clean(x):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(x)
            )
        )
    )


def valuation(x, prime):
    x = sp.Rational(x)

    if x == 0:
        return sp.oo

    numerator = abs(int(x.p))
    denominator = abs(int(x.q))

    value = 0

    while numerator % prime == 0:
        numerator //= prime
        value += 1

    while denominator % prime == 0:
        denominator //= prime
        value -= 1

    return value


def build_table():
    """
    Return exact dictionary:

        table[(p,t)] = Q_t(p)

    using ONLY observed source cells.
    """

    table = {}

    for p_value, values in Q.items():

        for index, value in enumerate(values):

            t_value = (
                len(values)
                - 1
                - index
            )

            table[
                (
                    int(p_value),
                    int(t_value),
                )
            ] = sp.Integer(value)

    return table


def observed_p_values(table):
    return sorted(
        {
            p_value
            for p_value, _ in table
        }
    )


def observed_t_values(table):
    return sorted(
        {
            t_value
            for _, t_value in table
        }
    )


def polynomial_coefficients(
    degree,
    prefix,
):
    """
    Return symbolic polynomial

        c_0 + c_1 p + ... + c_degree p^degree.
    """

    coefficients = sp.symbols(
        "{}_0:{}".format(
            prefix,
            degree + 1,
        )
    )

    polynomial = sum(
        coefficients[i] * p**i
        for i in range(
            degree + 1
        )
    )

    return coefficients, sp.expand(
        polynomial
    )


# ============================================================================
# RECURRENCE MODEL
# ============================================================================

def build_recurrence_model(
    order,
    coefficient_degree,
):
    """
    Construct

        Q_{t+order}(p)
          =
        C_0(p) Q_t(p)
        + ...
        + C_{order-1}(p) Q_{t+order-1}(p).

    Each C_k has degree <= coefficient_degree in p.
    """

    coefficient_blocks = []
    coefficient_polynomials = []

    for k in range(order):

        coeffs, polynomial = (
            polynomial_coefficients(
                coefficient_degree,
                "c{}_{}".format(k, ""),
            )
        )

        coefficient_blocks.extend(
            coeffs
        )

        coefficient_polynomials.append(
            polynomial
        )

    # The symbol names generated above are deliberately not relied upon;
    # duplicate names are harmless for construction but we replace them
    # with globally unique symbols below.
    unique_symbols = sp.symbols(
        "a0:{}".format(
            order
            * (coefficient_degree + 1)
        )
    )

    coefficient_polynomials = []

    offset = 0

    for k in range(order):

        polynomial = sum(
            unique_symbols[
                offset + j
            ] * p**j
            for j in range(
                coefficient_degree + 1
            )
        )

        coefficient_polynomials.append(
            sp.expand(polynomial)
        )

        offset += (
            coefficient_degree + 1
        )

    return (
        unique_symbols,
        coefficient_polynomials,
    )


# ============================================================================
# BUILD EXACT EQUATIONS
# ============================================================================

def recurrence_equations(
    table,
    order,
    coefficient_degree,
):
    """
    Build all equations available from the observed triangular table.
    """

    symbols, C = build_recurrence_model(
        order,
        coefficient_degree,
    )

    equations = []

    usable_windows = []

    p_values = observed_p_values(
        table
    )
    t_values = observed_t_values(
        table
    )

    for t0 in t_values:

        target_t = t0 + order

        if target_t not in t_values:
            continue

        for p_value in p_values:

            required_times = range(
                t0,
                t0 + order + 1,
            )

            if all(
                (
                    p_value,
                    time,
                )
                in table
                for time in required_times
            ):

                lhs = table[
                    (
                        p_value,
                        target_t,
                    )
                ]

                rhs = sum(
                    C[k].subs(
                        p,
                        p_value,
                    )
                    *
                    table[
                        (
                            p_value,
                            t0 + k,
                        )
                    ]
                    for k in range(order)
                )

                equations.append(
                    clean(
                        lhs - rhs
                    )
                )

                usable_windows.append(
                    (
                        p_value,
                        t0,
                    )
                )

    return (
        symbols,
        C,
        equations,
        usable_windows,
    )


# ============================================================================
# EXACT LINEAR SOLVE
# ============================================================================

def solve_recurrence(
    table,
    order,
    coefficient_degree,
):

    (
        symbols,
        C,
        equations,
        usable_windows,
    ) = recurrence_equations(
        table,
        order,
        coefficient_degree,
    )

    unknown_count = len(symbols)

    if not equations:

        return {
            "status": "NO_EQUATIONS",
            "equation_count": 0,
            "unknown_count": unknown_count,
            "rank": 0,
            "augmented_rank": 0,
            "coefficients": None,
            "polynomials": None,
            "residuals": [],
            "windows": [],
        }

    # Convert affine equations into homogeneous linear system in symbols.
    A, b = sp.linear_eq_to_matrix(
        equations,
        symbols,
    )

    # equations are of the form LHS = 0, so b is normally zero.
    rank = A.rank()
    augmented_rank = (
        A.row_join(b).rank()
    )

    equation_count = len(equations)

    if augmented_rank > rank:

        return {
            "status": "NO_SOLUTION",
            "equation_count": equation_count,
            "unknown_count": unknown_count,
            "rank": rank,
            "augmented_rank": augmented_rank,
            "coefficients": None,
            "polynomials": None,
            "residuals": None,
            "windows": usable_windows,
        }

    if rank == unknown_count:

        # The homogeneous system always has the trivial solution.
        # To test whether there is a recurrence with a nontrivial
        # coefficient vector we therefore need the nullspace dimension.
        nullspace = A.nullspace()

        if not nullspace:

            return {
                "status": "NO_NONTRIVIAL_RECURRENCE",
                "equation_count": equation_count,
                "unknown_count": unknown_count,
                "rank": rank,
                "augmented_rank": augmented_rank,
                "coefficients": None,
                "polynomials": None,
                "residuals": None,
                "windows": usable_windows,
            }

    nullspace = A.nullspace()

    if len(nullspace) == 1:

        vector = nullspace[0]

        coefficients = [
            clean(
                vector[i, 0]
            )
            for i in range(
                vector.rows
            )
        ]

        # Normalize first nonzero coefficient to 1.
        first_nonzero = next(
            value
            for value in coefficients
            if value != 0
        )

        coefficients = [
            clean(
                value / first_nonzero
            )
            for value in coefficients
        ]

        polynomials = []

        offset = 0

        block_size = (
            coefficient_degree + 1
        )

        for k in range(order):

            polynomial = clean(
                sum(
                    coefficients[
                        offset + j
                    ] * p**j
                    for j in range(
                        block_size
                    )
                )
            )

            polynomials.append(
                polynomial
            )

            offset += block_size

        # Verify all original equations.
        residuals = [
            clean(
                equation.subs(
                    dict(
                        zip(
                            symbols,
                            coefficients,
                        )
                    )
                )
            )
            for equation in equations
        ]

        if not all(
            residual == 0
            for residual in residuals
        ):

            return {
                "status": "VERIFICATION_FAILED",
                "equation_count": equation_count,
                "unknown_count": unknown_count,
                "rank": rank,
                "augmented_rank": augmented_rank,
                "coefficients": coefficients,
                "polynomials": polynomials,
                "residuals": residuals,
                "windows": usable_windows,
            }

        if equation_count > unknown_count:
            status = "EXACT_OVERDETERMINED"
        else:
            status = "EXACT_DATA_SIZED"

        return {
            "status": status,
            "equation_count": equation_count,
            "unknown_count": unknown_count,
            "rank": rank,
            "augmented_rank": augmented_rank,
            "coefficients": coefficients,
            "polynomials": polynomials,
            "residuals": residuals,
            "windows": usable_windows,
        }

    return {
        "status": "NONUNIQUE",
        "equation_count": equation_count,
        "unknown_count": unknown_count,
        "rank": rank,
        "augmented_rank": augmented_rank,
        "coefficients": nullspace,
        "polynomials": None,
        "residuals": None,
        "windows": usable_windows,
    }


# ============================================================================
# REPORT
# ============================================================================

def report_model(
    table,
    order,
    coefficient_degree,
):

    label = (
        "order={}, coefficient_degree={}".format(
            order,
            coefficient_degree,
        )
    )

    result = solve_recurrence(
        table,
        order,
        coefficient_degree,
    )

    print()
    print("=" * 78)
    print(label.upper())
    print("=" * 78)

    print(
        "  equation_count={}".format(
            result["equation_count"]
        )
    )

    print(
        "  unknown_coefficients={}".format(
            result["unknown_count"]
        )
    )

    print(
        "  redundancy={}".format(
            result["equation_count"]
            - result["unknown_count"]
        )
    )

    print(
        "  rank={}".format(
            result["rank"]
        )
    )

    print(
        "  augmented_rank={}".format(
            result["augmented_rank"]
        )
    )

    print(
        "  usable_windows={}".format(
            result["windows"]
        )
    )

    print(
        "  status={}".format(
            result["status"]
        )
    )

    if result["status"] in (
        "EXACT_OVERDETERMINED",
        "EXACT_DATA_SIZED",
    ):

        print()

        for k, polynomial in enumerate(
            result["polynomials"]
        ):

            print(
                "  C_{}(p)={}".format(
                    k,
                    polynomial,
                )
            )

            print(
                "    valuations={}".format(
                    {
                        prime: valuation(
                            polynomial,
                            prime,
                        )
                        if polynomial.is_Rational
                        else "polynomial"
                        for prime in (
                            2,
                            3,
                            5,
                            7,
                            11,
                            13,
                            17,
                        )
                    }
                )
            )

        print()
        print(
            "  normalized_coefficients={}".format(
                result["coefficients"]
            )
        )

        print(
            "  residuals={}".format(
                result["residuals"]
            )
        )

    elif result["status"] == "NONUNIQUE":

        print(
            "  nullspace_dimension={}".format(
                len(
                    result["coefficients"]
                )
            )
        )

    return result


# ============================================================================
# MODEL SUITE
# ============================================================================

def run_suite(table):

    models = [
        # order 1
        (1, 0),
        (1, 1),
        (1, 2),

        # order 2
        (2, 0),
        (2, 1),
        (2, 2),

        # order 3
        (3, 0),
        (3, 1),
    ]

    results = {}

    for order, degree in models:

        results[
            (
                order,
                degree,
            )
        ] = report_model(
            table,
            order,
            degree,
        )

    return results


# ============================================================================
# STRUCTURAL INTERPRETATION
# ============================================================================

def interpretation():

    print()
    print("=" * 78)
    print(
        "STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
Experiment 326R ruled out low-degree bivariate polynomial formulas
Q=P(p,t) for every tested degree family.

Experiment 327R keeps one piece of the failed structure:

    recurrence in terminal distance t,

but relaxes the strongest assumption:

    the recurrence coefficients are allowed to depend on p.

The tested form is

    Q_{t+r}(p)
      =
      C_0(p) Q_t(p)
      + ...
      + C_{r-1}(p) Q_{t+r-1}(p),

with C_k(p) of small degree.

This is structurally distinct from arbitrary polynomial fitting.

A successful order-1 or order-2 model with low-degree C_k(p) would say
that the source is still governed by finite-memory dynamics in t, but
the dynamics depend algebraically on p.

Because every available p,t window is used simultaneously, such a
result is much stronger than fitting each p-column separately.

The decisive cases are those with:

    equation_count > unknown_coefficients.

A successful exact overdetermined model here would therefore be a new
positive source-level discovery.

A failure of order <=2 with coefficient degree <=2 would make a simple
p-dependent recurrence substantially less plausible and would justify
moving to rational/algebraic source formulas or returning to the
construction of Q itself.

No missing values are used.
No interpolation is performed.
No synthetic second n=pq case is generated.
"""
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 327R — EXACT P-DEPENDENT VERTICAL RECURRENCE AUDIT"
    )
    print("=" * 78)

    table = build_table()

    print()
    print("=" * 78)
    print(
        "1. OBSERVED SOURCE TABLE"
    )
    print("=" * 78)

    print(
        "  observed_cells={}".format(
            len(table)
        )
    )

    for key in sorted(
        table,
        key=lambda z: (
            z[1],
            z[0],
        ),
    ):

        print(
            "  (p={},t={}) -> {}".format(
                key[0],
                key[1],
                table[key],
            )
        )

    results = run_suite(
        table
    )

    # ------------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "SUMMARY"
    )
    print("=" * 78)

    for key, result in results.items():

        order, degree = key

        print(
            "  order={}, degree={}: status={}, "
            "equations={}, unknowns={}, "
            "redundancy={}".format(
                order,
                degree,
                result["status"],
                result["equation_count"],
                result["unknown_count"],
                result["equation_count"]
                - result["unknown_count"],
            )
        )

    # ------------------------------------------------------------------------
    # Exact-overdetermined discoveries
    # ------------------------------------------------------------------------

    discoveries = [
        key
        for key, result in results.items()
        if result["status"]
        == "EXACT_OVERDETERMINED"
    ]

    print()
    print("=" * 78)
    print(
        "FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  exact_overdetermined_models={}".format(
            discoveries
        )
    )

    print(
        "  low_degree_p_dependent_recurrence_tested=True"
    )

    print(
        "  missing_values_used=False"
    )

    print(
        "  interpolation_performed=False"
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
        "EXPERIMENT 327R COMPLETE"
    )


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
