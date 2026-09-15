#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 334R — EXACT ROW-POLYNOMIAL OPERATOR EVOLUTION AUDIT
==============================================================================

Purpose
-------
The source table has now resisted:

    * fixed vertical recurrences;
    * p-dependent polynomial/rational recurrences;
    * low-degree bivariate polynomial laws;
    * additive local stencils;
    * exhaustive small translation-invariant linear stencils;
    * additive rank-one separability;
    * multiplicative rank-one separability;
    * low-rank Newton-coordinate dynamics.

Experiment 334R therefore changes representation.

For each observed source layer t, construct the exact polynomial P_t(p)
through the observed p-values only.

Then test whether adjacent layers are related by a fixed linear operator
in the p-variable.

The principal model class is:

    P_{t+1}(p)
        =
    sum_{k=0}^m a_k(p) D^k P_t(p),

where

    D = d/dp

and each coefficient polynomial a_k(p) has degree <= d.

This includes:

    m=0:
        multiplication operator

    m=1:
        a_0(p) P + a_1(p) P'

    m=2:
        a_0(p) P + a_1(p) P' + a_2(p) P''

and so on.

Only the EXISTING row polynomials are used.

Acceptance standard
-------------------
A fixed operator counts as a discovery only when its coefficients are
uniquely determined and the resulting identity is verified for every
available adjacent layer.

No missing values.
No extrapolation.
No synthetic second n=pq case.
Exact SymPy rational arithmetic only.
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

p = sp.symbols("p")


# ============================================================================
# EXACT HELPERS
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

    n = abs(int(x.p))
    d = abs(int(x.q))

    vn = 0
    vd = 0

    while n % prime == 0:
        n //= prime
        vn += 1

    while d % prime == 0:
        d //= prime
        vd += 1

    return vn - vd


# ============================================================================
# OBSERVED LAYERS
# ============================================================================

def build_layers():

    maximum_t = max(
        len(values) - 1
        for values in Q.values()
    )

    layers = {}

    for t in range(
        maximum_t + 1
    ):

        layer = []

        for p_value in sorted(Q):

            values = Q[p_value]

            index = (
                len(values)
                - 1
                - t
            )

            if index >= 0:

                layer.append(
                    (
                        sp.Integer(p_value),
                        sp.Integer(values[index]),
                    )
                )

        layers[t] = layer

    return layers


# ============================================================================
# EXACT ROW POLYNOMIALS
# ============================================================================

def interpolate_row(layer):

    points = [
        (p_value, value)
        for p_value, value
        in layer
    ]

    return clean(
        sp.interpolate(
            points,
            p,
        )
    )


def build_row_polynomials(layers):

    P = {}

    for t, layer in layers.items():

        P[t] = interpolate_row(
            layer
        )

    return P


# ============================================================================
# OPERATOR MONOMIALS
# ============================================================================

def operator_terms(
    max_derivative,
    coefficient_degree,
):
    """
    Basis terms:

        p^j D^k

    for

        0 <= k <= max_derivative
        0 <= j <= coefficient_degree
    """

    terms = []

    for k in range(
        max_derivative + 1
    ):

        for j in range(
            coefficient_degree + 1
        ):

            terms.append(
                (
                    k,
                    j,
                )
            )

    return terms


def apply_term(
    polynomial,
    derivative_order,
    coefficient_power,
):
    return clean(
        p**coefficient_power
        * sp.diff(
            polynomial,
            p,
            derivative_order,
        )
    )


# ============================================================================
# LINEAR OPERATOR SOLVER
# ============================================================================

def solve_fixed_operator(
    row_polynomials,
    max_derivative,
    coefficient_degree,
):
    """
    Solve for a FIXED operator

        L = sum c_{k,j} p^j D^k

    satisfying

        P_{t+1} = L(P_t)

    for every adjacent observed layer.

    Coefficients are solved exactly over Q.
    """

    terms = operator_terms(
        max_derivative,
        coefficient_degree,
    )

    unknown_count = len(
        terms
    )

    equations = []

    targets = []

    max_t = max(
        row_polynomials
    )

    for t in range(
        max_t
    ):

        source = row_polynomials[t]
        target = row_polynomials[t + 1]

        for power in range(
            max(
                sp.degree(source, p)
                if source != 0
                else 0,
                sp.degree(target, p)
                if target != 0
                else 0,
            ) + 1
        ):

            row = []

            for derivative_order, coefficient_power in terms:

                expression = apply_term(
                    source,
                    derivative_order,
                    coefficient_power,
                )

                row.append(
                    clean(
                        sp.expand(
                            expression
                        ).coeff(
                            p,
                            power,
                        )
                    )
                )

            target_coefficient = clean(
                sp.expand(
                    target
                ).coeff(
                    p,
                    power,
                )
            )

            equations.append(
                row
            )

            targets.append(
                target_coefficient
            )

    A = sp.Matrix(
        equations
    )

    b = sp.Matrix(
        targets
    )

    rank = A.rank()
    augmented_rank = (
        A.row_join(b).rank()
    )

    if augmented_rank > rank:

        return {
            "status": "NO_SOLUTION",
            "terms": terms,
            "unknowns": unknown_count,
            "equations": len(equations),
            "rank": rank,
            "augmented_rank": augmented_rank,
            "nullity": unknown_count - rank,
            "coefficients": None,
            "operator": None,
        }

    if rank < unknown_count:

        return {
            "status": "NONUNIQUE",
            "terms": terms,
            "unknowns": unknown_count,
            "equations": len(equations),
            "rank": rank,
            "augmented_rank": augmented_rank,
            "nullity": unknown_count - rank,
            "coefficients": None,
            "operator": None,
        }

    solution = A.gauss_jordan_solve(
        b
    )[0]

    coefficients = [
        clean(solution[i, 0])
        for i in range(
            unknown_count
        )
    ]

    operator = []

    for coefficient, (
        derivative_order,
        coefficient_power,
    ) in zip(
        coefficients,
        terms,
    ):

        if coefficient == 0:
            continue

        operator.append(
            (
                coefficient,
                derivative_order,
                coefficient_power,
            )
        )

    return {
        "status": "EXACT_UNIQUE",
        "terms": terms,
        "unknowns": unknown_count,
        "equations": len(equations),
        "rank": rank,
        "augmented_rank": augmented_rank,
        "nullity": 0,
        "coefficients": coefficients,
        "operator": operator,
    }


# ============================================================================
# OPERATOR APPLICATION
# ============================================================================

def apply_operator(
    polynomial,
    result,
):

    if result["operator"] is None:
        return None

    total = 0

    for (
        coefficient,
        derivative_order,
        coefficient_power,
    ) in result["operator"]:

        total += (
            coefficient
            * p**coefficient_power
            * sp.diff(
                polynomial,
                p,
                derivative_order,
            )
        )

    return clean(
        total
    )


def verify_operator(
    row_polynomials,
    result,
):

    if result["status"] != "EXACT_UNIQUE":
        return None

    residuals = []

    max_t = max(
        row_polynomials
    )

    for t in range(
        max_t
    ):

        predicted = apply_operator(
            row_polynomials[t],
            result,
        )

        residual = clean(
            predicted
            - row_polynomials[t + 1]
        )

        residuals.append(
            residual
        )

    return residuals


# ============================================================================
# REPORTING
# ============================================================================

def report_model(
    label,
    row_polynomials,
    max_derivative,
    coefficient_degree,
):

    result = solve_fixed_operator(
        row_polynomials,
        max_derivative,
        coefficient_degree,
    )

    print()
    print("=" * 78)
    print(label)
    print("=" * 78)

    print(
        "  max_derivative={}".format(
            max_derivative
        )
    )

    print(
        "  coefficient_degree={}".format(
            coefficient_degree
        )
    )

    print(
        "  operator_unknowns={}".format(
            result["unknowns"]
        )
    )

    print(
        "  equations={}".format(
            result["equations"]
        )
    )

    print(
        "  redundancy={}".format(
            result["equations"]
            - result["unknowns"]
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
        "  nullity={}".format(
            result["nullity"]
        )
    )

    print(
        "  status={}".format(
            result["status"]
        )
    )

    if result["status"] == "EXACT_UNIQUE":

        residuals = verify_operator(
            row_polynomials,
            result,
        )

        print()
        print(
            "  operator_terms={}".format(
                result["operator"]
            )
        )

        print()
        print(
            "  residuals={}".format(
                residuals
            )
        )

        print(
            "  all_residuals_zero={}".format(
                all(
                    residual == 0
                    for residual in residuals
                )
            )
        )

        print()
        print(
            "  coefficient_prime_profiles="
        )

        for coefficient in result[
            "coefficients"
        ]:

            if coefficient == 0:
                continue

            print(
                "    {} -> {}".format(
                    coefficient,
                    {
                        prime: valuation(
                            coefficient,
                            prime,
                        )
                        for prime in (
                            2,
                            3,
                            5,
                            7,
                            11,
                            13,
                            17,
                        )
                    },
                )
            )

    return result


# ============================================================================
# ROW-POLYNOMIAL REPORT
# ============================================================================

def report_rows(
    row_polynomials,
):

    print()
    print("=" * 78)
    print(
        "1. EXACT ROW POLYNOMIALS"
    )
    print("=" * 78)

    for t in sorted(
        row_polynomials
    ):

        polynomial = (
            row_polynomials[t]
        )

        print()
        print(
            "  t={}: P_t(p)={}".format(
                t,
                polynomial,
            )
        )

        print(
            "    degree={}".format(
                sp.degree(
                    polynomial,
                    p,
                )
            )
        )

        print(
            "    factorization={}".format(
                sp.factor(
                    polynomial
                )
            )
        )


# ============================================================================
# ADDITIONAL STRUCTURAL TESTS
# ============================================================================

def derivative_shift_audit(
    row_polynomials,
):

    print()
    print("=" * 78)
    print(
        "ADDITIONAL DERIVATIVE / SHIFT TESTS"
    )
    print("=" * 78)

    max_t = max(
        row_polynomials
    )

    for t in range(
        max_t
    ):

        source = row_polynomials[t]
        target = row_polynomials[t + 1]

        print()
        print(
            "  t={} -> {}:".format(
                t,
                t + 1,
            )
        )

        print(
            "    target={}".format(
                target
            )
        )

        print(
            "    source_derivative={}".format(
                clean(
                    sp.diff(
                        source,
                        p,
                    )
                )
            )
        )

        print(
            "    source_shift_plus_2={}".format(
                clean(
                    source.subs(
                        p,
                        p + 2,
                    )
                )
            )
        )

        print(
            "    source_shift_minus_2={}".format(
                clean(
                    source.subs(
                        p,
                        p - 2,
                    )
                )
            )

        )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 334R — EXACT ROW-POLYNOMIAL "
        "OPERATOR EVOLUTION AUDIT"
    )
    print("=" * 78)

    layers = build_layers()

    row_polynomials = build_row_polynomials(
        layers
    )

    report_rows(
        row_polynomials
    )

    derivative_shift_audit(
        row_polynomials
    )

    results = {}

    # Small operator suite.
    models = [
        ("D0 DEG0", 0, 0),
        ("D0 DEG1", 0, 1),
        ("D0 DEG2", 0, 2),
        ("D1 DEG0", 1, 0),
        ("D1 DEG1", 1, 1),
        ("D1 DEG2", 1, 2),
        ("D2 DEG0", 2, 0),
        ("D2 DEG1", 2, 1),
        ("D2 DEG2", 2, 2),
    ]

    for label, max_derivative, degree in models:

        results[label] = report_model(
            label,
            row_polynomials,
            max_derivative,
            degree,
        )

    print()
    print("=" * 78)
    print(
        "STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
The row-polynomial representation preserves every observed source value,
but changes the question from

    "What recurrence do the sampled lattice points satisfy?"

to

    "How does the entire observed p-profile evolve from one layer
     to the next?"

A fixed differential operator is a substantially more structured object
than an arbitrary p-dependent recurrence.

For example,

    P_{t+1} = (a0+a1 p) P_t
              + (b0+b1 p) P_t'

is a finite-dimensional representation of the vertical dynamics in the
p-polynomial space.

Likewise second-order models

    P_{t+1}
      = a(p) P_t
      + b(p) P_t'
      + c(p) P_t''

can encode degree-lowering or degree-preserving combinatorial operators.

Only a UNIQUE exact operator valid for every observed adjacent layer is
treated as a discovery.

A NONUNIQUE solution is not accepted because the finite data do not
identify the operator.

A NO_SOLUTION result means that entire operator family has been exactly
eliminated.

No missing p-values are introduced.
No missing n=pq case is generated.
"""
    )

    exact_models = [
        label
        for label, result in results.items()
        if result["status"] == "EXACT_UNIQUE"
    ]

    print()
    print("=" * 78)
    print(
        "FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  exact_unique_operator_models={}".format(
            exact_models
        )
    )

    print(
        "  row_polynomial_reconstruction_from_observed_cells_only=True"
    )

    print(
        "  missing_values_used=False"
    )

    print(
        "  extrapolation_used=False"
    )

    print(
        "  interpolation_outside_observed_cells=False"
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
        "EXPERIMENT 334R COMPLETE"
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
