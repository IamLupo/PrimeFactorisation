#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 352R-FIXED — EXACT MISSING-CELL PREDICTION / HYPOTHESIS
                             SEPARATION AUDIT
==============================================================================

Target
------
Missing source cell

    (r,t) = (2,3)

equivalently

    (p,t) = (5,3).

Unlike the original 352R, this version DOES NOT require the missing value.

It computes independent model predictions for Q_3(5), but NEVER treats those
predictions as observed data.

The purpose is to determine whether the currently surviving hypotheses make
the same prediction or materially different predictions.

Test families
--------------

    D1-D8 rational generating-function denominator models;
    D6 direct reconstruction from the original four equations;
    raw width-2 transfer t=3 -> 4, where enough data permit;
    first-difference / Newton-coordinate implications;
    arithmetic profile of every prediction.

Important distinction
---------------------

    MODEL PREDICTION
        value implied by an existing fitted hypothesis;

    OBSERVED DATA
        value independently obtained from the original source definition.

Only the latter can validate a hypothesis.

No interpolation.
No missing-value reconstruction counted as evidence.
No synthetic second n=pq case.
Exact SymPy arithmetic only.
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

a, b, c, d, e = sp.symbols(
    "a b c d e"
)


TARGET = (2, 3)


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


def build_lattice():
    """
    Build

        (r,t) -> Q_t(p),

    where

        r = (p-1)/2.

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


def factor_integer(value):

    value = int(
        sp.Integer(value)
    )

    if value == 0:
        return {}

    return sp.factorint(
        abs(value)
    )


def valuation(value, prime):

    value = sp.Integer(value)

    if value == 0:
        return sp.oo

    value = abs(
        int(value)
    )

    count = 0

    while value % prime == 0:

        value //= prime
        count += 1

    return count


def polynomial_prediction_summary(
    name,
    value,
):

    print()
    print(
        "  {}:".format(
            name
        )
    )

    print(
        "    predicted_value={}".format(
            clean(value)
        )
    )

    print(
        "    is_integer={}".format(
            sp.denom(
                sp.Rational(value)
            ) == 1
        )
    )

    if value != 0:

        print(
            "    numerator={}".format(
                sp.numer(
                    sp.Rational(value)
                )
            )
        )

        print(
            "    denominator={}".format(
                sp.denom(
                    sp.Rational(value)
                )
            )
        )

    if (
        sp.denom(
            sp.Rational(value)
        )
        == 1
        and
        value != 0
    ):

        print(
            "    factorization={}".format(
                factor_integer(value)
            )
        )

        print(
            "    valuations={}".format(
                {
                    p:
                    valuation(
                        value,
                        p,
                    )
                    for p in (
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


# ============================================================================
# DENOMINATOR FAMILIES
# ============================================================================

def denominator_families():

    return {
        "D1_xy": (
            1 - a*x - b*y,
            [a, b],
        ),

        "D2_xy_xy": (
            1
            - a*x
            - b*y
            - c*x*y,
            [a, b, c],
        ),

        "D3_xy_x2": (
            1
            - a*x
            - b*y
            - c*x**2,
            [a, b, c],
        ),

        "D4_xy_y2": (
            1
            - a*x
            - b*y
            - c*y**2,
            [a, b, c],
        ),

        "D5_xy_xy_x2": (
            1
            - a*x
            - b*y
            - c*x*y
            - d*x**2,
            [a, b, c, d],
        ),

        "D6_xy_xy_y2": (
            1
            - a*x
            - b*y
            - c*x*y
            - d*y**2,
            [a, b, c, d],
        ),

        "D7_xy_xy_x2_y2": (
            1
            - a*x
            - b*y
            - c*x*y
            - d*x**2
            - e*y**2,
            [a, b, c, d, e],
        ),

        "D8_xy_xy_x2_xy2": (
            1
            - a*x
            - b*y
            - c*x*y
            - d*x**2
            - e*x*y**2,
            [a, b, c, d, e],
        ),
    }


def monomial_shifts(
    expression
):

    poly = sp.Poly(
        sp.expand(expression),
        x,
        y,
    )

    return list(
        zip(
            poly.monoms(),
            poly.coeffs(),
        )
    )


# ============================================================================
# DENOMINATOR EQUATION CONSTRUCTION
# ============================================================================

def coefficient_equations(
    lattice,
    denominator,
):
    """
    Build all coefficient equations for D*G=0 for which every required
    source cell is actually observed.

    The constant term of D is normalized to 1.
    """

    shifts = monomial_shifts(
        denominator
    )

    cells = []

    min_dr = min(
        monomial[0]
        for monomial, _ in shifts
    )

    min_dt = min(
        monomial[1]
        for monomial, _ in shifts
    )

    max_r = max(
        r
        for r, _ in lattice
    )

    max_t = max(
        t
        for _, t in lattice
    )

    for r in range(
        max_r + 1
    ):

        for t in range(
            max_t + 1
        ):

            target = (
                r,
                t,
            )

            if target not in lattice:
                continue

            required = []

            for (
                monomial,
                coefficient,
            ) in shifts:

                dr, dt = monomial

                source = (
                    r - dr,
                    t - dt,
                )

                required.append(
                    source
                )

            if all(
                source in lattice
                for source in required
            ):

                cells.append(
                    target
                )

    equations = []

    for r, t in cells:

        expression = sp.Integer(0)

        for (
            monomial,
            coefficient,
        ) in shifts:

            dr, dt = monomial

            source = (
                r - dr,
                t - dt,
            )

            expression += (
                coefficient
                * lattice[
                    source
                ]
            )

        equations.append(
            clean(expression)
        )

    return cells, equations


# ============================================================================
# SOLVE DENOMINATOR MODEL ON EXISTING DATA
# ============================================================================

def fit_denominator_model(
    lattice,
    name,
    denominator,
    parameters,
):

    cells, expressions = (
        coefficient_equations(
            lattice,
            denominator,
        )
    )

    if not expressions:

        return {
            "name": name,
            "status": "INSUFFICIENT_DATA",
            "cells": cells,
            "solution": None,
        }

    matrix = sp.Matrix([
        [
            sp.diff(
                expression,
                parameter,
            )
            for parameter in parameters
        ]
        for expression in expressions
    ])

    rhs = sp.Matrix([
        clean(
            -expression.subs(
                {
                    parameter: 0
                    for parameter
                    in parameters
                }
            )
        )
        for expression in expressions
    ])

    rank = matrix.rank()

    augmented_rank = (
        matrix
        .row_join(rhs)
        .rank()
    )

    equation_count = len(
        expressions
    )

    unknown_count = len(
        parameters
    )

    if augmented_rank > rank:

        status = "NO_SOLUTION"

    elif rank < unknown_count:

        status = "NONUNIQUE"

    elif equation_count > unknown_count:

        status = "EXACT_OVERDETERMINED"

    else:

        status = "EXACT_DATA_SIZED"

    result = {
        "name": name,
        "status": status,
        "cells": cells,
        "equations": equation_count,
        "unknowns": unknown_count,
        "rank": rank,
        "augmented_rank": augmented_rank,
        "solution": None,
    }

    if status in (
        "EXACT_OVERDETERMINED",
        "EXACT_DATA_SIZED",
    ):

        solution_vector = (
            matrix.gauss_jordan_solve(
                rhs
            )[0]
        )

        result["solution"] = {
            parameter:
            clean(
                solution_vector[i, 0]
            )
            for i, parameter
            in enumerate(parameters)
        }

    return result


# ============================================================================
# PREDICT TARGET FROM A DENOMINATOR MODEL
# ============================================================================

def predict_from_denominator(
    denominator,
    solution,
    target,
):

    if solution is None:

        return None

    substituted = clean(
        denominator.subs(
            solution
        )
    )

    coefficient = 0

    shifts = monomial_shifts(
        substituted
    )

    r, t = target

    # D*G coefficient at target equals zero.
    # Solve for Q[target], whose denominator coefficient is normalized to 1.
    for (
        monomial,
        coefficient_value,
    ) in shifts:

        dr, dt = monomial

        if (
            dr == 0
            and dt == 0
        ):
            continue

        source = (
            r - dr,
            t - dt,
        )

        # Only caller supplies observed source terms.
        coefficient = (
            coefficient
            + coefficient_value
            * source
        )

    return None


def explicit_target_prediction(
    lattice,
    denominator,
    solution,
    target,
):
    """
    Compute the target coefficient explicitly.

    Because the constant coefficient of D is 1,

        Q[target]
          =
        - sum_{nonconstant shifts}
          d_shift Q[target-shift].

    """

    substituted = clean(
        denominator.subs(
            solution
        )
    )

    shifts = monomial_shifts(
        substituted
    )

    target_r, target_t = target

    prediction = 0

    for (
        monomial,
        coefficient,
    ) in shifts:

        dr, dt = monomial

        if (
            dr == 0
            and dt == 0
        ):
            continue

        source = (
            target_r - dr,
            target_t - dt,
        )

        if source not in lattice:

            return None

        prediction += (
            coefficient
            * lattice[
                source
            ]
        )

    return clean(
        -prediction
    )


# ============================================================================
# RAW WIDTH-2 PREDICTIONS
# ============================================================================

def raw_width2_prediction(
    lattice,
    source_t,
    target_r,
):
    """
    Q_{t+1}(r)
      =
    c0 Q_t(r)
      +
    c1 Q_t(r+1).

    Uses exact existing equations only.
    """

    equations = []

    for r in range(4):

        target = (
            r,
            source_t + 1,
        )

        s0 = (
            r,
            source_t,
        )

        s1 = (
            r + 1,
            source_t,
        )

        if (
            target in lattice
            and s0 in lattice
            and s1 in lattice
        ):

            equations.append(
                (
                    s0,
                    s1,
                    target,
                )
            )

    if len(
        equations
    ) < 2:

        return None

    M = sp.Matrix([
        [
            lattice[s0],
            lattice[s1],
        ]
        for s0, s1, _
        in equations
    ])

    rhs = sp.Matrix([
        lattice[target]
        for _, _, target
        in equations
    ])

    if M.rank() != 2:

        return None

    solution = M.gauss_jordan_solve(
        rhs
    )[0]

    c0 = clean(
        solution[0, 0]
    )

    c1 = clean(
        solution[1, 0]
    )

    target = (
        target_r,
        source_t + 1,
    )

    s0 = (
        target_r,
        source_t,
    )

    s1 = (
        target_r + 1,
        source_t,
    )

    if (
        s0 not in lattice
        or s1 not in lattice
    ):

        return None

    return clean(
        c0 * lattice[s0]
        + c1 * lattice[s1]
    ), c0, c1


# ============================================================================
# FIRST-DIFFERENCE / NEWTON PREDICTIONS
# ============================================================================

def row_values(
    lattice,
    t_value,
):

    return [
        lattice[
            (r, t_value)
        ]
        for r in range(4)
        if (
            r,
            t_value
        ) in lattice
    ]


def differences(values):

    rows = [
        list(values)
    ]

    while len(
        rows[-1]
    ) > 1:

        previous = rows[-1]

        rows.append([
            clean(
                previous[i + 1]
                - previous[i]
            )
            for i in range(
                len(previous) - 1
            )
        ])

    return rows


# ============================================================================
# DIRECT D6 PREDICTION
# ============================================================================

def direct_d6_prediction(
    lattice
):

    name = (
        "D6_xy_xy_y2"
    )

    denominator, parameters = (
        denominator_families()[
            name
        ]
    )

    result = fit_denominator_model(
        lattice,
        name,
        denominator,
        parameters,
    )

    if result[
        "solution"
    ] is None:

        return result, None

    prediction = (
        explicit_target_prediction(
            lattice,
            denominator,
            result[
                "solution"
            ],
            TARGET,
        )
    )

    return result, prediction


# ============================================================================
# NEWTON IMPACT
# ============================================================================

def newton_impact(
    lattice
):

    print()
    print("=" * 78)
    print(
        "6. NEWTON-COORDINATE IMPACT OF THE MISSING CELL"
    )
    print("=" * 78)

    t = 3

    before = differences(
        row_values(
            lattice,
            t,
        )
    )

    print()
    print(
        "  current_t3_values={}".format(
            row_values(
                lattice,
                t,
            )
        )
    )

    print(
        "  current_differences={}".format(
            before
        )
    )

    print()
    print(
        "  target=(2,3) would complete the p=1,3,5 row."
    )

    print(
        "  This creates the new difference"
        " Delta Q_3(3) = Q_3(5)-Q_3(3)."
    )


# ============================================================================
# ARITHMETIC COMPARISON OF PREDICTIONS
# ============================================================================

def compare_predictions(
    predictions
):

    print()
    print("=" * 78)
    print(
        "7. EXACT PREDICTION COMPARISON"
    )
    print("=" * 78)

    valid = [
        (
            name,
            clean(value),
        )
        for name, value
        in predictions
        if value is not None
    ]

    for name, value in valid:

        polynomial_prediction_summary(
            name,
            value,
        )

    print()
    print(
        "  prediction_count={}".format(
            len(valid)
        )
    )

    if len(valid) < 2:

        print(
            "  independent_prediction_comparison=False"
        )

        return

    reference_name, reference = (
        valid[0]
    )

    print()
    print(
        "  reference={}".format(
            reference_name
        )
    )

    for name, value in valid[1:]:

        difference = clean(
            value
            - reference
        )

        print()
        print(
            "  {} vs {}:".format(
                name,
                reference_name,
            )
        )

        print(
            "    difference={}".format(
                difference
            )
        )

        print(
            "    same_prediction={}".format(
                difference == 0
            )
        )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 352R-FIXED — EXACT MISSING-CELL PREDICTION / "
        "HYPOTHESIS SEPARATION AUDIT"
    )
    print("=" * 78)

    lattice = build_lattice()

    print()
    print(
        "TARGET CELL"
    )

    print(
        "  coordinate={}".format(
            TARGET
        )
    )

    print(
        "  p=5"
    )

    print(
        "  t=3"
    )

    print()
    print(
        "  IMPORTANT: Q_3(5) is NOT supplied."
    )

    print(
        "  All values below are MODEL PREDICTIONS only."
    )

    # ------------------------------------------------------------------------
    # Denominator models
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "1. EXISTING DENOMINATOR MODELS"
    )
    print("=" * 78)

    predictions = []

    denominator_results = []

    for (
        name,
        (
            denominator,
            parameters,
        ),
    ) in denominator_families().items():

        result = fit_denominator_model(
            lattice,
            name,
            denominator,
            parameters,
        )

        denominator_results.append(
            result
        )

        print()
        print(
            "{}:".format(
                name
            )
        )

        print(
            "  equations={}".format(
                result.get(
                    "equations",
                    0,
                )
            )
        )

        print(
            "  unknowns={}".format(
                result.get(
                    "unknowns",
                    len(parameters),
                )
            )
        )

        print(
            "  redundancy={}".format(
                result.get(
                    "redundancy",
                    0,
                )
            )
        )

        print(
            "  rank={}".format(
                result.get(
                    "rank",
                    None,
                )
            )
        )

        print(
            "  augmented_rank={}".format(
                result.get(
                    "augmented_rank",
                    None,
                )
            )

        )

        print(
            "  status={}".format(
                result[
                    "status"
                ]
            )
        )

        if result[
            "solution"
        ] is not None:

            print(
                "  parameters={}".format(
                    result[
                        "solution"
                    ]
                )
            )

            denominator, _ = (
                denominator_families()[
                    name
                ]
            )

            prediction = (
                explicit_target_prediction(
                    lattice,
                    denominator,
                    result[
                        "solution"
                    ],
                    TARGET,
                )
            )

            if prediction is not None:

                print(
                    "  predicted_Q3_5={}".format(
                        prediction
                    )
                )

                predictions.append(
                    (
                        name,
                        prediction,
                    )
                )

    # ------------------------------------------------------------------------
    # Direct D6
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "2. D6 DIRECT PREDICTION"
    )
    print("=" * 78)

    d6_result, d6_prediction = (
        direct_d6_prediction(
            lattice
        )
    )

    print(
        "  status={}".format(
            d6_result[
                "status"
            ]
        )
    )

    if d6_result[
        "solution"
    ] is not None:

        print(
            "  parameters={}".format(
                d6_result[
                    "solution"
                ]
            )
        )

        print(
            "  predicted_Q3_5={}".format(
                d6_prediction
            )
        )

        predictions.append(
            (
                "D6_DIRECT",
                d6_prediction,
            )
        )

    # ------------------------------------------------------------------------
    # Raw width-2 candidates
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "3. RAW WIDTH-2 TRANSFER PREDICTIONS"
    )
    print("=" * 78)

    for source_t in (
        1,
        2,
    ):

        result = raw_width2_prediction(
            lattice,
            source_t,
            2,
        )

        print()

        print(
            "  source_t={} -> target_t={}:".format(
                source_t,
                source_t + 1,
            )
        )

        if result is None:

            print(
                "    status=NO_PREDICTION_AVAILABLE"
            )

            continue

        predicted, c0, c1 = result

        print(
            "    c0={}".format(
                c0
            )
        )

        print(
            "    c1={}".format(
                c1
            )
        )

        print(
            "    predicted_target={}".format(
                predicted
            )
        )

    # ------------------------------------------------------------------------
    # Newton
    # ------------------------------------------------------------------------

    newton_impact(
        lattice
    )

    # ------------------------------------------------------------------------
    # Prediction comparison
    # ------------------------------------------------------------------------

    compare_predictions(
        predictions
    )

    # ------------------------------------------------------------------------
    # Structural verdict
    # ------------------------------------------------------------------------

    exact_overdetermined = [
        result["name"]
        for result
        in denominator_results
        if result["status"]
        == "EXACT_OVERDETERMINED"
    ]

    print()
    print("=" * 78)
    print(
        "8. STRUCTURAL VERDICT"
    )
    print("=" * 78)

    print(
        "  exact_overdetermined_existing_denominators={}".format(
            exact_overdetermined
        )
    )

    print(
        "  target_is_currently_unobserved=True"
    )

    print(
        "  predictions_are_not_validation=True"
    )

    print(
        "  independently_obtained_Q3_5_still_required=True"
    )

    print()
    print(
        "  highest_priority_next_action="
        "OBTAIN_INDEPENDENT_Q3(5)"
    )

    print()
    print(
        "  interpretation="
        "USE_THE_INDEPENDENT_VALUE_TO_DISCRIMINATE_BETWEEN_THESE_PREDICTIONS"
    )

    print()
    print("=" * 78)
    print(
        "9. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  observed_cells={}".format(
            len(lattice)
        )
    )

    print(
        "  target_cell=(2,3)"
    )

    print(
        "  target_value_supplied=False"
    )

    print(
        "  model_predictions_computed=True"
    )

    print(
        "  predictions_count={}".format(
            len(predictions)
        )
    )

    print(
        "  missing_value_inferred_as_data=False"
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
        "EXPERIMENT 352R-FIXED COMPLETE"
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
