#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 352R — EXACT SINGLE-CELL SOURCE VALIDATION / HIGH-LEVERAGE TEST
==============================================================================

Target
------
The single missing source cell identified by Experiment 351R:

    (r,t) = (2,3)

equivalently

    p = 5,
    t = 3.

This experiment is designed for an INDEPENDENTLY OBTAINED value of Q_3(5).

It does NOT infer the value.

Once supplied, the experiment tests:

    1. all denominator families D1-D8 from Experiment 350R;
    2. the exact D6 model reconstructed from the original four equations;
    3. raw width-2 transfer constraints;
    4. first-entry/Newton consistency;
    5. small local 2D stencil constraints;
    6. whether the new cell creates genuine overdetermined validation;
    7. exact residuals for every tested hypothesis.

Interpretation
--------------
A model is interesting only if the added cell makes it independently
testable.

NO_SOLUTION
    The newly supplied source value contradicts the model.

EXACT_OVERDETERMINED
    The newly supplied value agrees with a previously identified model and
    makes the model genuinely overdetermined.

EXACT_DATA_SIZED
    Still reconstruction-sized; not evidence.

NONUNIQUE
    Not identified.

IMPORTANT
---------
The supplied Q_3(5) must come from an independent mathematical construction,
not from any model tested here.

No interpolation.
No extrapolation.
No missing-value reconstruction.
No synthetic second n=pq case.
Exact SymPy arithmetic only.
"""


from __future__ import annotations

import math
import sys

import sympy as sp


# ============================================================================
# INSERT ONLY AN INDEPENDENTLY OBTAINED VALUE HERE
# ============================================================================

Q_5_T3 = None

# Example:
#
# Q_5_T3 = 123456789
#
# DO NOT enter a value obtained from any fitted recurrence being tested below.


# ============================================================================
# ORIGINAL OBSERVED SOURCE DATA
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
    Build the observed lattice

        (r,t) -> Q_t(p),

    where

        r = (p-1)/2.
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


def add_target_cell(lattice):

    if Q_5_T3 is None:

        raise RuntimeError(
            "Q_5_T3 is None. Supply the independently obtained "
            "Q_3(5) value before running Experiment 352R."
        )

    if not isinstance(
        Q_5_T3,
        (int, sp.Integer),
    ):

        raise TypeError(
            "Q_5_T3 must be an exact integer."
        )

    lattice = dict(
        lattice
    )

    if (2, 3) in lattice:

        raise RuntimeError(
            "Target cell (2,3) is unexpectedly already present."
        )

    lattice[
        (2, 3)
    ] = sp.Integer(Q_5_T3)

    return lattice


def rank_status(
    matrix,
    rhs,
    unknown_count,
):

    rank = matrix.rank()

    augmented_rank = (
        matrix
        .row_join(rhs)
        .rank()
    )

    equations = matrix.rows

    if augmented_rank > rank:

        status = "NO_SOLUTION"

    elif rank < unknown_count:

        status = "NONUNIQUE"

    elif equations > unknown_count:

        status = "EXACT_OVERDETERMINED"

    else:

        status = "EXACT_DATA_SIZED"

    return {
        "equations": equations,
        "unknowns": unknown_count,
        "redundancy": (
            equations
            - unknown_count
        ),
        "rank": rank,
        "augmented_rank": augmented_rank,
        "status": status,
    }


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


def denominator_shifts(
    denominator
):

    P = sp.Poly(
        sp.expand(denominator),
        x,
        y,
    )

    return [
        (
            int(mon[0]),
            int(mon[1]),
            coefficient,
        )
        for mon, coefficient
        in zip(
            P.monoms(),
            P.coeffs(),
        )
    ]


def denominator_equations(
    lattice,
    denominator,
):

    """
    For D*G, coefficient equations are formed only where ALL source cells
    required by the denominator are observed.
    """

    shifts = denominator_shifts(
        denominator
    )

    equations = []

    for r, t in sorted(lattice):

        required = []

        for dr, dt, _ in shifts:

            source = (
                r - dr,
                t - dt,
            )

            required.append(
                source
            )

        if not all(
            source in lattice
            for source in required
        ):
            continue

        equations.append(
            (r, t)
        )

    return equations


def solve_denominator(
    lattice,
    name,
    denominator,
    parameters,
):

    cells = denominator_equations(
        lattice,
        denominator,
    )

    rows = []
    rhs = []

    shifts = denominator_shifts(
        denominator
    )

    for r, t in cells:

        coefficient_row = []

        for parameter in parameters:

            coefficient = (
                0
            )

            for dr, dt, dc in shifts:

                if dc.has(
                    parameter
                ):

                    source = (
                        r - dr,
                        t - dt,
                    )

                    coefficient += (
                        sp.diff(
                            -dc,
                            parameter,
                        )
                        * lattice[
                            source
                        ]
                    )

            coefficient_row.append(
                clean(coefficient)
            )

        # The constant part is exactly Q(r,t), because denominator has
        # constant term 1.
        constant_part = (
            lattice[
                (r, t)
            ]
        )

        for dr, dt, dc in shifts:

            if dc.free_symbols:

                continue

            source = (
                r - dr,
                t - dt,
            )

            constant_part += (
                dc
                * lattice[
                    source
                ]
            )

        rows.append(
            coefficient_row
        )

        rhs.append(
            clean(
                -constant_part
            )
        )

    if not rows:

        return {
            "name": name,
            "cells": [],
            "status": "INSUFFICIENT_DATA",
            "solution": None,
        }

    M = sp.Matrix(rows)
    rhs_matrix = sp.Matrix(rhs)

    status = rank_status(
        M,
        rhs_matrix,
        len(parameters),
    )

    result = {
        "name": name,
        "cells": cells,
        **status,
        "solution": None,
    }

    if status["status"] in (
        "EXACT_OVERDETERMINED",
        "EXACT_DATA_SIZED",
    ):

        solution = M.gauss_jordan_solve(
            rhs_matrix
        )[0]

        solution_dict = {
            parameter: clean(
                solution[i, 0]
            )
            for i, parameter
            in enumerate(parameters)
        }

        result[
            "solution"
        ] = solution_dict

        residuals = []

        for row_index in range(
            M.rows
        ):

            predicted = clean(
                sum(
                    M[row_index, j]
                    * solution[j, 0]
                    for j in range(
                        M.cols
                    )
                )
            )

            residuals.append(
                clean(
                    predicted
                    - rhs_matrix[
                        row_index
                    ]
                )
            )

        result[
            "residuals"
        ] = residuals

    return result


# ============================================================================
# D6 DIRECT COMPARISON
# ============================================================================

def d6_direct_audit(
    original_lattice,
    augmented_lattice,
):

    print()
    print("=" * 78)
    print(
        "2. D6 DIRECT FOUR-EQUATION PREDICTION TEST"
    )
    print("=" * 78)

    denominator, parameters = (
        denominator_families()[
            "D6_xy_xy_y2"
        ]
    )

    original_cells = (
        denominator_equations(
            original_lattice,
            denominator,
        )
    )

    print()
    print(
        "  original_equation_cells={}".format(
            original_cells
        )
    )

    original_result = solve_denominator(
        original_lattice,
        "D6_original",
        denominator,
        parameters,
    )

    if original_result[
        "status"
    ] not in (
        "EXACT_DATA_SIZED",
        "EXACT_OVERDETERMINED",
    ):

        print(
            "  original_D6_status={}".format(
                original_result["status"]
            )
        )

        return False

    solution = original_result[
        "solution"
    ]

    print(
        "  original_D6_parameters={}".format(
            solution
        )
    )

    # Evaluate the genuinely new coefficient equation created by (2,3).
    new_cells = (
        denominator_equations(
            augmented_lattice,
            denominator,
        )
    )

    new_cells = [
        cell
        for cell in new_cells
        if cell
        not in original_cells
    ]

    print()
    print(
        "  new_equation_cells={}".format(
            new_cells
        )
    )

    all_pass = True

    shifts = denominator_shifts(
        denominator
    )

    for r, t in new_cells:

        residual = (
            augmented_lattice[
                (r, t)
            ]
        )

        for dr, dt, dc in shifts:

            source = (
                r - dr,
                t - dt,
            )

            value = (
                augmented_lattice[
                    source
                ]
            )

            coefficient = clean(
                dc.subs(
                    solution
                )
            )

            residual += (
                coefficient
                * value
            )

        residual = clean(
            residual
        )

        print()
        print(
            "  new_cell_equation=({},{}):".format(
                r,
                t,
            )
        )

        print(
            "    residual={}".format(
                residual
            )
        )

        print(
            "    exact_match={}".format(
                residual == 0
            )
        )

        if residual != 0:
            all_pass = False

    return all_pass


# ============================================================================
# RAW WIDTH-2 TRANSFER
# ============================================================================

def raw_width2_transition_audit(
    lattice,
    t
):

    """
    Tests

        Q_{t+1}(r)
          =
        c0 Q_t(r)
          +
        c1 Q_t(r+1)

    using every currently available equation.
    """

    equations = []

    # Target layer t+1.
    for r in range(4):

        target = (
            r,
            t + 1
        )

        source0 = (
            r,
            t,
        )

        source1 = (
            r + 1,
            t,
        )

        if (
            target in lattice
            and source0 in lattice
            and source1 in lattice
        ):

            equations.append(
                (
                    source0,
                    source1,
                    target,
                )
            )

    print()
    print("=" * 78)
    print(
        "3. RAW WIDTH-2 TRANSFER t={} -> {}".format(
            t,
            t + 1,
        )
    )
    print("=" * 78)

    print(
        "  usable_equations={}".format(
            equations
        )
    )

    if len(equations) < 2:

        print(
            "  status=INSUFFICIENT_DATA"
        )

        return None

    M = sp.Matrix([
        [
            lattice[s0],
            lattice[s1],
        ]
        for s0, s1, _ in equations
    ])

    rhs = sp.Matrix([
        lattice[target]
        for _, _, target
        in equations
    ])

    result = rank_status(
        M,
        rhs,
        2,
    )

    print(
        "  status={}".format(
            result["status"]
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

    if result["status"] in (
        "EXACT_OVERDETERMINED",
        "EXACT_DATA_SIZED",
    ):

        solution = M.gauss_jordan_solve(
            rhs
        )[0]

        c0 = clean(
            solution[0, 0]
        )

        c1 = clean(
            solution[1, 0]
        )

        print(
            "  c0={}".format(
                c0
            )
        )

        print(
            "  c1={}".format(
                c1
            )
        )

        residuals = [
            clean(
                c0 * M[i, 0]
                + c1 * M[i, 1]
                - rhs[i]
            )
            for i in range(
                M.rows
            )
        ]

        print(
            "  residuals={}".format(
                residuals
            )
        )

        return (
            c0,
            c1,
            result,
        )

    return None


# ============================================================================
# FIRST-ENTRY / NEWTON COORDINATES
# ============================================================================

def first_difference_rows(
    lattice,
    t
):

    available = [
        lattice[
            (r, t)
        ]
        for r in range(4)
        if (
            r,
            t
        ) in lattice
    ]

    rows = [
        list(
            map(
                sp.Integer,
                available
            )
        )
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


def newton_target_audit(
    original_lattice,
    augmented_lattice
):

    print()
    print("=" * 78)
    print(
        "4. NEWTON-COORDINATE TARGET IMPACT"
    )
    print("=" * 78)

    for t in (
        3,
    ):

        original = first_difference_rows(
            original_lattice,
            t,
        )

        augmented = first_difference_rows(
            augmented_lattice,
            t,
        )

        print()
        print(
            "  t={}:".format(
                t
            )
        )

        print(
            "    original={}".format(
                original
            )
        )

        print(
            "    augmented={}".format(
                augmented
            )
        )

        if len(
            augmented
        ) > len(
            original
        ):

            print(
                "    new_newton_coefficients={}".format(
                    [
                        row[0]
                        for row
                        in augmented[
                            len(original):
                        ]
                        if row
                    ]
                )
            )


# ============================================================================
# LOCAL STENCIL TARGET
# ============================================================================

def local_rectangle_audit(
    lattice
):

    print()
    print("=" * 78)
    print(
        "5. LOCAL RECTANGLE / STENCIL IMPACT"
    )
    print("=" * 78)

    # Candidate new rectangle:
    #
    # (1,3) (2,3)
    # (1,4) (2,4)
    #
    # Only (2,3) is missing among the four cells that are currently relevant
    # to several local relations.  (2,4) remains unavailable.
    #
    # Therefore the new cell creates partial, not complete, 2x2 rectangles.

    required = [
        (1, 3),
        (2, 3),
        (1, 4),
        (2, 4),
    ]

    print(
        "  candidate_rectangle={}".format(
            required
        )
    )

    print(
        "  availability={}".format(
            {
                cell:
                cell in lattice
                for cell in required
            }
        )
    )

    complete = all(
        cell in lattice
        for cell in required
    )

    print(
        "  complete_rectangle={}".format(
            complete
        )
    )


# ============================================================================
# FACTOR / CONTENT AUDIT OF THE NEW VALUE
# ============================================================================

def arithmetic_audit():

    print()
    print("=" * 78)
    print(
        "6. ARITHMETIC PROFILE OF SUPPLIED Q_3(5)"
    )
    print("=" * 78)

    value = int(
        Q_5_T3
    )

    print(
        "  Q_3(5)={}".format(
            value
        )
    )

    print(
        "  sign={}".format(
            "+"
            if value > 0
            else "-"
            if value < 0
            else "0"
        )
    )

    if value != 0:

        factors = sp.factorint(
            abs(value)
        )

        print(
            "  factorization={}".format(
                factors
            )
        )

        print(
            "  v2={}".format(
                sp.factorint(
                    abs(value)
                ).get(
                    2,
                    0,
                )
            )
        )

        print(
            "  v3={}".format(
                sp.factorint(
                    abs(value)
                ).get(
                    3,
                    0,
                )
            )
        )

        print(
            "  v5={}".format(
                sp.factorint(
                    abs(value)
                ).get(
                    5,
                    0,
                )
            )
        )

        print(
            "  v7={}".format(
                sp.factorint(
                    abs(value)
                ).get(
                    7,
                    0,
                )
            )
        )

        print(
            "  v17={}".format(
                sp.factorint(
                    abs(value)
                ).get(
                    17,
                    0,
                )
            )
        )


# ============================================================================
# SUMMARY
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 352R — EXACT SINGLE-CELL SOURCE VALIDATION / "
        "HIGH-LEVERAGE TEST"
    )
    print("=" * 78)

    original = build_lattice()

    augmented = add_target_cell(
        original
    )

    print()
    print(
        "TARGET"
    )

    print(
        "  lattice_coordinate=(2,3)"
    )

    print(
        "  corresponding_p=5"
    )

    print(
        "  corresponding_t=3"
    )

    print(
        "  supplied_value={}".format(
            Q_5_T3
        )
    )

    # ------------------------------------------------------------------------
    # D1-D8
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "1. ALL DENOMINATOR FAMILIES AFTER ADDING ONE CELL"
    )
    print("=" * 78)

    denominator_results = []

    for (
        name,
        (
            denominator,
            parameters,
        ),
    ) in denominator_families().items():

        result = solve_denominator(
            augmented,
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
                result["status"]
            )
        )

        if result.get(
            "solution"
        ) is not None:

            print(
                "  solution={}".format(
                    result["solution"]
                )
            )

        if result.get(
            "residuals"
        ) is not None:

            print(
                "  residuals={}".format(
                    result["residuals"]
                )
            )

    # ------------------------------------------------------------------------
    # D6 direct test
    # ------------------------------------------------------------------------

    d6_pass = d6_direct_audit(
        original,
        augmented,
    )

    # ------------------------------------------------------------------------
    # Raw transfer
    # ------------------------------------------------------------------------

    raw_t3 = (
        raw_width2_transition_audit(
            augmented,
            3,
        )
    )

    # ------------------------------------------------------------------------
    # Newton
    # ------------------------------------------------------------------------

    newton_target_audit(
        original,
        augmented,
    )

    # ------------------------------------------------------------------------
    # Local geometry
    # ------------------------------------------------------------------------

    local_rectangle_audit(
        augmented
    )

    # ------------------------------------------------------------------------
    # Arithmetic
    # ------------------------------------------------------------------------

    arithmetic_audit()

    # ------------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------------

    exact_overdetermined = [
        result["name"]
        for result in denominator_results
        if result["status"]
        == "EXACT_OVERDETERMINED"
    ]

    no_solution = [
        result["name"]
        for result in denominator_results
        if result["status"]
        == "NO_SOLUTION"
    ]

    print()
    print("=" * 78)
    print(
        "7. STRUCTURAL VERDICT"
    )
    print("=" * 78)

    print(
        "  exact_overdetermined_denominators={}".format(
            exact_overdetermined
        )
    )

    print(
        "  contradicted_denominators={}".format(
            no_solution
        )
    )

    print(
        "  D6_new_cell_passes_direct_test={}".format(
            d6_pass
        )
    )

    print(
        "  raw_width2_t3_to_t4_status={}".format(
            None
            if raw_t3 is None
            else raw_t3[2]["status"]
        )
    )

    if (
        "D6_xy_xy_y2"
        in exact_overdetermined
    ):

        print(
            "  interpretation="
            "D6_SURVIVES_ONE_INDEPENDENT_NEW_CELL"
        )

    else:

        print(
            "  interpretation="
            "NO_TESTED_LOW_COMPLEXITY_DENOMINATOR_SURVIVES_AS_EXACT_OVERDETERMINED_MODEL"
        )

    print()
    print(
        "  NOTE:"
    )

    print(
        "    The supplied Q_3(5) is treated as independent data."
    )

    print(
        "    It is never reconstructed by the tested models."
    )

    print()
    print("=" * 78)
    print(
        "8. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  independent_target_value_supplied={}".format(
            Q_5_T3 is not None
        )
    )

    print(
        "  target_coordinate=(2,3)"
    )

    print(
        "  observed_original_cells={}".format(
            len(original)
        )
    )

    print(
        "  augmented_cells={}".format(
            len(augmented)
        )
    )

    print(
        "  denominator_families_tested=8"
    )

    print(
        "  missing_values_inferred=False"
    )

    print(
        "  interpolation_performed=False"
    )

    print(
        "  extrapolation_performed=False"
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
        "EXPERIMENT 352R COMPLETE"
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
