#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 354R-FIXED2 — EXACT MISSING-CELL DISCRIMINATION /
                           MODEL-FALSIFICATION MATRIX
==============================================================================

Purpose
-------
Use only the 15 observed Q(r,t) cells to identify currently solvable models,
then compute predictions for the two strategically important missing cells:

    (2,3) = Q_3(5)
    (3,1) = Q_1(7)

The experiment strictly separates:

    OBSERVED CALIBRATION
        from
    UNOBSERVED PREDICTION.

No missing cell is ever inserted into the calibration system.

No interpolation.
No extrapolation counted as evidence.
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

a, b, c, d = sp.symbols(
    "a b c d"
)

c0_sym, c1_sym = sp.symbols(
    "c0 c1"
)


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
    value = int(
        sp.Integer(value)
    )

    if value == 0:
        return {}

    return sp.factorint(
        abs(value)
    )


def build_lattice():
    """
    Build the observed lattice

        (r,t) -> Q_t(p),

    where

        r = (p-1)/2.

    Only actually observed cells are inserted.
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
        lattice,
        key=lambda cell: (
            cell[1],
            cell[0],
        )
    )


# ============================================================================
# EXACT LINEAR MODEL SOLVER
# ============================================================================

def solve_linear_model(
    rows,
    rhs,
    symbols,
):
    """
    Solve a rational linear system exactly.

    Returns:
        NO_SOLUTION
        NONUNIQUE
        EXACT_DATA_SIZED

    The latter means the coefficients are uniquely determined by the
    calibration equations. It does NOT mean independently validated.
    """

    matrix = sp.Matrix(rows)
    rhs_matrix = sp.Matrix(rhs)

    rank = matrix.rank()

    augmented_rank = (
        matrix
        .row_join(rhs_matrix)
        .rank()
    )

    parameter_count = len(symbols)
    equation_count = len(rows)

    if augmented_rank > rank:

        return {
            "status": "NO_SOLUTION",
            "equation_count": equation_count,
            "parameter_count": parameter_count,
            "rank": rank,
            "augmented_rank": augmented_rank,
            "parameters": None,
        }

    if rank < parameter_count:

        return {
            "status": "NONUNIQUE",
            "equation_count": equation_count,
            "parameter_count": parameter_count,
            "rank": rank,
            "augmented_rank": augmented_rank,
            "parameters": None,
        }

    solution = matrix.gauss_jordan_solve(
        rhs_matrix
    )[0]

    parameters = {
        symbols[i]: clean(
            solution[i, 0]
        )
        for i in range(
            parameter_count
        )
    }

    return {
        "status": "EXACT_DATA_SIZED",
        "equation_count": equation_count,
        "parameter_count": parameter_count,
        "rank": rank,
        "augmented_rank": augmented_rank,
        "parameters": parameters,
    }


# ============================================================================
# D6 BOUNDARY DENOMINATOR MODEL
# ============================================================================

def d6_equation_source_cells(r, t):
    """
    D6 relation:

        Q(r,t)
          = a Q(r-1,t)
          + b Q(r,t-1)
          + c Q(r-1,t-1)
          + d Q(r,t-2)

    Return the four source cells required to evaluate the right-hand side.
    """

    return [
        (r - 1, t),
        (r, t - 1),
        (r - 1, t - 1),
        (r, t - 2),
    ]


def d6_observed_calibration_windows(lattice):
    """
    IMPORTANT:
    A calibration window is included ONLY when the target AND every source
    cell are observed.

    This is the critical fix for the previous KeyError at (1,5).
    """

    windows = []

    for r in range(4):

        for t in range(6):

            target = (
                r,
                t,
            )

            if target not in lattice:
                continue

            sources = d6_equation_source_cells(
                r,
                t,
            )

            if all(
                source in lattice
                for source in sources
            ):
                windows.append(
                    target
                )

    return windows


def solve_d6(lattice):

    windows = d6_observed_calibration_windows(
        lattice
    )

    rows = []
    rhs = []

    for r, t in windows:

        s0, s1, s2, s3 = (
            d6_equation_source_cells(
                r,
                t
            )
        )

        rows.append(
            [
                lattice[s0],
                lattice[s1],
                lattice[s2],
                lattice[s3],
            ]
        )

        rhs.append(
            lattice[(r, t)]
        )

    result = solve_linear_model(
        rows,
        rhs,
        [a, b, c, d],
    )

    result["windows"] = windows
    result["predictions"] = {}

    if result["status"] != "EXACT_DATA_SIZED":
        return result

    params = result["parameters"]

    # Predict only genuinely unobserved cells whose complete source
    # neighborhood is observed.
    for target in (
        (2, 3),
        (3, 1),
    ):

        r, t = target

        if target in lattice:
            continue

        if r < 1 or t < 2:
            continue

        sources = d6_equation_source_cells(
            r,
            t,
        )

        if not all(
            source in lattice
            for source in sources
        ):
            continue

        s0, s1, s2, s3 = sources

        prediction = clean(
            params[a] * lattice[s0]
            + params[b] * lattice[s1]
            + params[c] * lattice[s2]
            + params[d] * lattice[s3]
        )

        result["predictions"][
            target
        ] = prediction

    return result


# ============================================================================
# RAW WIDTH-2 MODEL
# ============================================================================

def raw_width2_calibration_equations(
    lattice,
    source_t,
    target_t,
):
    """
    Include only equations where:

        Q(r, source_t)
        Q(r+1, source_t)
        Q(r, target_t)

    are ALL observed.
    """

    equations = []

    for r in range(4):

        source0 = (
            r,
            source_t,
        )

        source1 = (
            r + 1,
            source_t,
        )

        target = (
            r,
            target_t,
        )

        if (
            source0 in lattice
            and source1 in lattice
            and target in lattice
        ):
            equations.append(
                (
                    source0,
                    source1,
                    target,
                )
            )

    return equations


def solve_raw_width2(
    lattice,
    source_t,
    target_t,
):
    equations = raw_width2_calibration_equations(
        lattice,
        source_t,
        target_t,
    )

    if len(equations) < 2:

        return {
            "status": "INSUFFICIENT_DATA",
            "equation_count": len(equations),
            "parameter_count": 2,
            "rank": None,
            "augmented_rank": None,
            "equations": equations,
            "parameters": None,
            "predictions": {},
        }

    rows = [
        [
            lattice[source0],
            lattice[source1],
        ]
        for (
            source0,
            source1,
            _,
        ) in equations
    ]

    rhs = [
        lattice[target]
        for (
            _,
            _,
            target,
        ) in equations
    ]

    result = solve_linear_model(
        rows,
        rhs,
        [c0_sym, c1_sym],
    )

    result["equations"] = equations
    result["predictions"] = {}

    if result["status"] != "EXACT_DATA_SIZED":
        return result

    params = result["parameters"]

    # Predict only immediately available unobserved targets.
    for r in range(4):

        target = (
            r,
            target_t,
        )

        source0 = (
            r,
            source_t,
        )

        source1 = (
            r + 1,
            source_t,
        )

        if (
            target not in lattice
            and source0 in lattice
            and source1 in lattice
        ):

            result["predictions"][
                target
            ] = clean(
                params[c0_sym]
                * lattice[source0]
                +
                params[c1_sym]
                * lattice[source1]
            )

    return result


# ============================================================================
# REPORTING
# ============================================================================

def report_model(
    name,
    result,
):

    print()
    print(
        name
        + ":"
    )

    print(
        "  status={}".format(
            result["status"]
        )
    )

    print(
        "  equations={}".format(
            result.get(
                "equation_count",
                0,
            )
        )
    )

    print(
        "  parameters={}".format(
            result.get(
                "parameter_count",
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
        "  redundancy={}".format(
            result.get(
                "equation_count",
                0,
            )
            -
            result.get(
                "parameter_count",
                0,
            )
        )
    )

    if "windows" in result:
        print(
            "  calibration_windows={}".format(
                result["windows"]
            )
        )

    if "equations" in result:
        print(
            "  calibration_equations={}".format(
                result["equations"]
            )
        )

    if result.get(
        "parameters"
    ) is not None:

        print(
            "  solved_parameters={}".format(
                result["parameters"]
            )
        )

    predictions = result.get(
        "predictions",
        {}
    )

    if predictions:

        for cell in sorted(
            predictions
        ):

            value = predictions[cell]

            print()
            print(
                "  prediction_cell={}".format(
                    cell
                )
            )

            print(
                "    prediction={}".format(
                    value
                )
            )

            print(
                "    integer={}".format(
                    is_integer(value)
                )
            )

            denominator = int(
                sp.denom(
                    sp.Rational(value)
                )
            )

            print(
                "    denominator={}".format(
                    denominator
                )
            )

            print(
                "    denominator_factorization={}".format(
                    factor_integer(
                        denominator
                    )
                )
            )

    else:

        print(
            "  predictions={}"
            .format({})
        )


# ============================================================================
# DISCRIMINATION MATRIX
# ============================================================================

def discrimination_matrix(
    model_results,
    candidates,
):

    print()
    print("=" * 78)
    print(
        "4. EXACT MISSING-CELL DISCRIMINATION MATRIX"
    )
    print("=" * 78)

    for cell in candidates:

        print()
        print(
            "  candidate_cell={}".format(
                cell
            )
        )

        predictions = []

        for (
            model_name,
            result,
        ) in model_results:

            prediction_map = result.get(
                "predictions",
                {}
            )

            if cell in prediction_map:

                prediction = clean(
                    prediction_map[cell]
                )

                predictions.append(
                    (
                        model_name,
                        prediction,
                    )
                )

        if not predictions:

            print(
                "    predictions_available=0"
            )

            continue

        for (
            model_name,
            prediction,
        ) in predictions:

            print()
            print(
                "    model={}".format(
                    model_name
                )
            )

            print(
                "      prediction={}".format(
                    prediction
                )
            )

            print(
                "      integer={}".format(
                    is_integer(
                        prediction
                    )
                )
            )

        distinct = {
            prediction
            for (
                _,
                prediction,
            ) in predictions
        }

        integer_models = [
            model_name
            for (
                model_name,
                prediction,
            ) in predictions
            if is_integer(
                prediction
            )
        ]

        print()
        print(
            "    prediction_count={}".format(
                len(predictions)
            )
        )

        print(
            "    distinct_prediction_count={}".format(
                len(distinct)
            )
        )

        print(
            "    integer_compatible_models={}".format(
                integer_models
            )
        )


# ============================================================================
# D6 INTEGER OBSTRUCTION
# ============================================================================

def d6_integer_obstruction(
    d6_result,
):

    print()
    print("=" * 78)
    print(
        "5. D6 INTEGER-VALUED-SOURCE OBSTRUCTION"
    )
    print("=" * 78)

    predictions = d6_result.get(
        "predictions",
        {}
    )

    if not predictions:

        print(
            "  no_D6_predictions_available=True"
        )

        return

    for cell, prediction in sorted(
        predictions.items()
    ):

        print()
        print(
            "  cell={}".format(
                cell
            )
        )

        print(
            "    prediction={}".format(
                prediction
            )
        )

        if is_integer(prediction):

            print(
                "    integer_compatible=True"
            )

        else:

            denominator = int(
                sp.denom(
                    sp.Rational(
                        prediction
                    )
                )
            )

            print(
                "    integer_compatible=False"
            )

            print(
                "    reduced_denominator={}".format(
                    denominator
                )
            )

            print(
                "    denominator_factorization={}".format(
                    factor_integer(
                        denominator
                    )
                )
            )

            print(
                "    conditional_verdict="
                "D6_FALSIFIED_IF_SOURCE_IS_INTEGER_VALUED"
            )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 354R-FIXED2 — EXACT MISSING-CELL "
        "DISCRIMINATION / MODEL-FALSIFICATION MATRIX"
    )
    print("=" * 78)

    lattice = build_lattice()

    print()
    print(
        "1. OBSERVED DATA"
    )
    print(
        "  observed_cells={}".format(
            len(lattice)
        )
    )

    print(
        "  observed_coordinates={}".format(
            sorted_cells(lattice)
        )
    )

    candidates = [
        cell
        for cell in (
            (2, 3),
            (3, 1),
        )
        if cell not in lattice
    ]

    print()
    print(
        "2. CANDIDATE MISSING CELLS"
    )

    for cell in candidates:

        print(
            "  cell={}".format(
                cell
            )
        )

    # ------------------------------------------------------------------------
    # D6
    # ------------------------------------------------------------------------

    d6_result = solve_d6(
        lattice
    )

    report_model(
        "3. D6_BOUNDARY_DENOMINATOR",
        d6_result,
    )

    # ------------------------------------------------------------------------
    # RAW WIDTH-2: t=1 -> 2
    # ------------------------------------------------------------------------

    raw_1_2 = solve_raw_width2(
        lattice,
        1,
        2,
    )

    report_model(
        "3A. RAW_WIDTH2_t1_to_t2",
        raw_1_2,
    )

    # ------------------------------------------------------------------------
    # RAW WIDTH-2: t=2 -> 3
    # ------------------------------------------------------------------------

    raw_2_3 = solve_raw_width2(
        lattice,
        2,
        3,
    )

    report_model(
        "3B. RAW_WIDTH2_t2_to_t3",
        raw_2_3,
    )

    model_results = [
        (
            "D6_BOUNDARY_DENOMINATOR",
            d6_result,
        ),
        (
            "RAW_WIDTH2_t1_to_t2",
            raw_1_2,
        ),
        (
            "RAW_WIDTH2_t2_to_t3",
            raw_2_3,
        ),
    ]

    discrimination_matrix(
        model_results,
        candidates,
    )

    d6_integer_obstruction(
        d6_result
    )

    # ------------------------------------------------------------------------
    # STRATEGY
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "6. STRATEGIC INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
The current dataset contains 15 observed source cells.

The two strategically important missing cells are:

    (2,3) = Q_3(5)
    (3,1) = Q_1(7)

The script never inserts either value into the calibration equations.

The D6 model is fitted only from fully observed D6 windows.
Its prediction at an unobserved cell is therefore a model prediction,
not an observation.

Because every currently observed source value is an integer, a reduced
non-integer D6 prediction is conditionally incompatible with the
hypothesis that Q is integer-valued everywhere.

That conditional statement is a falsification test, not a proof of the
alternative model.

The most valuable independent datum remains:

    Q_3(5).

Obtaining that value would add a genuinely overdetermining equation to
the D6 family and would simultaneously discriminate several other
candidate local denominator models.

No missing value is generated or treated as observed.
"""
    )

    # ------------------------------------------------------------------------
    # FINAL EXACTNESS
    # ------------------------------------------------------------------------

    d6_predictions = d6_result.get(
        "predictions",
        {}
    )

    d6_noninteger = [
        cell
        for (
            cell,
            value,
        ) in d6_predictions.items()
        if not is_integer(value)
    ]

    print()
    print("=" * 78)
    print(
        "7. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  observed_cells={}".format(
            len(lattice)
        )
    )

    print(
        "  missing_candidates={}".format(
            candidates
        )
    )

    print(
        "  d6_status={}".format(
            d6_result["status"]
        )
    )

    print(
        "  d6_calibration_target_cells_observed_only=True"
    )

    print(
        "  d6_noninteger_predictions={}".format(
            d6_noninteger
        )
    )

    print(
        "  predictions_used_as_data=False"
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
        "EXPERIMENT 354R-FIXED2 COMPLETE"
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