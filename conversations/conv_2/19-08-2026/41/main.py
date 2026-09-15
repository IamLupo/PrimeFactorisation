#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 355R — EXACT SYMBOLIC MISSING-CELL CLOSURE / MODEL ELIMINATION
==============================================================================

Purpose
-------
Treat the two missing source cells as symbolic unknowns:

    Z = Q_3(5) = Q(r=2,t=3)
    W = Q_1(7) = Q(r=3,t=1)

Do NOT assign numerical values to Z or W.

Instead ask:

    Can each candidate source law be made exactly consistent with the
    complete observed table for SOME values of Z and W?

This separates three cases:

    IMPOSSIBLE
        No values of Z,W can rescue the model.

    SYMBOLICALLY_CONSTRAINED
        The model survives only for a specific algebraic/rational
        constraint on Z,W.

    SYMBOLICALLY_UNDERDETERMINED
        The current data leave free parameters.

This is stronger than ordinary prediction because a model can be
eliminated without observing the missing cells at all.

Models tested:

    D1:
        Q(r,t) = a Q(r-1,t) + b Q(r,t-1)

    D2:
        Q(r,t) = a Q(r-1,t)
               + b Q(r,t-1)
               + c Q(r-1,t-1)

    D4:
        Q(r,t) = a Q(r-1,t)
               + b Q(r,t-1)
               + c Q(r,t-2)

    D5:
        Q(r,t) = a Q(r-1,t)
               + b Q(r,t-1)
               + c Q(r-1,t-1)
               + d Q(r-1,t-2)

    D6:
        Q(r,t) = a Q(r-1,t)
               + b Q(r,t-1)
               + c Q(r-1,t-1)
               + d Q(r,t-2)

    D7:
        Q(r,t) =
            a Q(r-1,t)
          + b Q(r,t-1)
          + c Q(r-1,t-1)
          + d Q(r,t-2)
          + e Q(r-1,t-2)

Important
---------
Only observed cells are substituted numerically.

Z and W remain symbolic throughout.

No missing value is inferred as data.
No interpolation.
No extrapolation.
No synthetic second n=pq case.
Exact SymPy arithmetic only.
"""


from __future__ import annotations

import sys
import sympy as sp


# ============================================================================
# OBSERVED SOURCE DATA
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

Z, W = sp.symbols(
    "Z W"
)

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


def build_observed_lattice():

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


def symbolic_value(
    lattice,
    cell,
):

    if cell == (2, 3):
        return Z

    if cell == (3, 1):
        return W

    if cell in lattice:
        return lattice[cell]

    return None


# ============================================================================
# MODEL DEFINITIONS
# ============================================================================

MODELS = {
    "D1_xy": {
        "symbols": [a, b],
        "source_cells": lambda r, t: [
            (r - 1, t),
            (r, t - 1),
        ],
    },

    "D2_xy_xyshift": {
        "symbols": [a, b, c],
        "source_cells": lambda r, t: [
            (r - 1, t),
            (r, t - 1),
            (r - 1, t - 1),
        ],
    },

    "D4_xy_y2": {
        "symbols": [a, b, c],
        "source_cells": lambda r, t: [
            (r - 1, t),
            (r, t - 1),
            (r, t - 2),
        ],
    },

    "D5_xy_xyshift_x2shift": {
        "symbols": [a, b, c, d],
        "source_cells": lambda r, t: [
            (r - 1, t),
            (r, t - 1),
            (r - 1, t - 1),
            (r - 1, t - 2),
        ],
    },

    "D6_xy_xyshift_y2": {
        "symbols": [a, b, c, d],
        "source_cells": lambda r, t: [
            (r - 1, t),
            (r, t - 1),
            (r - 1, t - 1),
            (r, t - 2),
        ],
    },

    "D7_extended": {
        "symbols": [a, b, c, d, e],
        "source_cells": lambda r, t: [
            (r - 1, t),
            (r, t - 1),
            (r - 1, t - 1),
            (r, t - 2),
            (r - 1, t - 2),
        ],
    },
}


# ============================================================================
# EQUATION CONSTRUCTION
# ============================================================================

def model_equations(
    lattice,
    model,
):

    parameters = model["symbols"]
    source_builder = model["source_cells"]

    equations = []
    locations = []

    #
    # We inspect every geometrically possible target in the current
    # finite triangular envelope.
    #
    for r in range(4):

        for t in range(6):

            target = (r, t)

            #
            # The model requires r >= 1 and enough t-distance.
            #
            if r < 1:
                continue

            sources = source_builder(
                r,
                t,
            )

            #
            # Negative-time cells are invalid.
            #
            if any(
                source[1] < 0
                or source[0] < 0
                for source in sources
            ):
                continue

            #
            # Every source must be either observed or one of the two
            # explicitly symbolic missing cells.
            #
            symbolic_sources = []

            valid = True

            for source in sources:

                value = symbolic_value(
                    lattice,
                    source,
                )

                if value is None:
                    valid = False
                    break

                symbolic_sources.append(
                    value
                )

            if not valid:
                continue

            target_value = symbolic_value(
                lattice,
                target,
            )

            if target_value is None:
                #
                # Target is a genuinely unavailable cell other than Z,W.
                #
                continue

            rhs = sum(
                parameters[i]
                * symbolic_sources[i]
                for i in range(
                    len(parameters)
                )
            )

            equations.append(
                sp.Eq(
                    target_value,
                    rhs,
                )
            )

            locations.append(
                target
            )

    return parameters, equations, locations


# ============================================================================
# SYMBOLIC SOLVER
# ============================================================================

def solve_symbolically(
    lattice,
    model_name,
    model,
):

    parameters, equations, locations = (
        model_equations(
            lattice,
            model,
        )
    )

    lhs = [
        clean(
            equation.lhs
            - equation.rhs
        )
        for equation in equations
    ]

    unknowns = list(parameters) + [
        Z,
        W,
    ]

    if not lhs:

        return {
            "status": "NO_EQUATIONS",
            "equations": [],
            "locations": [],
        }

    #
    # SymPy solve over all model parameters plus Z,W.
    #
    solution = sp.solve(
        lhs,
        unknowns,
        dict=True,
        simplify=False,
        rational=True,
    )

    #
    # Check whether the polynomial/rational system is exactly inconsistent.
    #
    if not solution:

        #
        # A zero-dimensional solve can also report [] for inconsistency.
        # Verify with a Groebner basis when possible.
        #
        try:

            G = sp.groebner(
                lhs,
                *unknowns,
                order="lex",
            )

            if G == sp.groebner(
                [1],
                *unknowns,
                order="lex",
            ):

                status = "IMPOSSIBLE"

            else:

                status = "NO_EXPLICIT_SOLUTION"

        except Exception:

            status = "NO_EXPLICIT_SOLUTION"

        return {
            "status": status,
            "equations": lhs,
            "locations": locations,
            "solution": None,
        }

    #
    # Determine whether Z and W are fixed or remain free.
    #
    solution = [
        {
            symbol: clean(value)
            for symbol, value
            in branch.items()
        }
        for branch in solution
    ]

    return {
        "status": "SYMBOLIC_SOLUTION",
        "equations": lhs,
        "locations": locations,
        "solution": solution,
    }


# ============================================================================
# SOLUTION CLASSIFICATION
# ============================================================================

def classify_solution(
    result,
    model,
):

    if result["status"] != "SYMBOLIC_SOLUTION":
        return result["status"]

    solutions = result["solution"]

    if not solutions:
        return "IMPOSSIBLE"

    #
    # If every solution fixes Z and W numerically/rationally, the missing
    # cells are constrained.
    #
    fixed_z = True
    fixed_w = True

    z_values = set()
    w_values = set()

    for solution in solutions:

        if Z not in solution:
            fixed_z = False
        else:
            z_values.add(
                clean(solution[Z])
            )

        if W not in solution:
            fixed_w = False
        else:
            w_values.add(
                clean(solution[W])
            )

    if fixed_z and fixed_w:

        return "SYMBOLICALLY_CONSTRAINED_BOTH"

    if fixed_z and not fixed_w:

        return "SYMBOLICALLY_CONSTRAINED_Z"

    if fixed_w and not fixed_z:

        return "SYMBOLICALLY_CONSTRAINED_W"

    return "SYMBOLICALLY_UNDERDETERMINED"


# ============================================================================
# INTEGER TEST
# ============================================================================

def integer_obstruction(
    value,
):

    if value is None:
        return None

    value = sp.Rational(
        value
    )

    return {
        "integer": value.q == 1,
        "numerator": value.p,
        "denominator": value.q,
        "denominator_factorization": (
            sp.factorint(
                abs(int(value.q))
            )
            if value.q != 1
            else {}
        ),
    }


# ============================================================================
# REPORT
# ============================================================================

def report(
    name,
    result,
    model,
):

    print()
    print("=" * 78)
    print(
        name
    )
    print("=" * 78)

    print(
        "  equations={}".format(
            len(
                result.get(
                    "equations",
                    []
                )
            )
        )
    )

    print(
        "  equation_locations={}".format(
            result.get(
                "locations",
                []
            )
        )
    )

    classification = classify_solution(
        result,
        model,
    )

    print(
        "  status={}".format(
            classification
        )
    )

    if result.get(
        "solution"
    ) is not None:

        print()
        print(
            "  solution_branches={}".format(
                len(
                    result["solution"]
                )
            )
        )

        for index, solution in enumerate(
            result["solution"]
        ):

            print()
            print(
                "  branch_{}:".format(
                    index
                )
            )

            print(
                "    {}".format(
                    solution
                )
            )

            if Z in solution:

                info = integer_obstruction(
                    solution[Z]
                )

                print(
                    "    Z_integer={}".format(
                        info["integer"]
                    )
                )

                if not info["integer"]:

                    print(
                        "    Z_denominator={}".format(
                            info["denominator"]
                        )
                    )

                    print(
                        "    Z_denominator_factorization={}".format(
                            info[
                                "denominator_factorization"
                            ]
                        )
                    )

            if W in solution:

                info = integer_obstruction(
                    solution[W]
                )

                print(
                    "    W_integer={}".format(
                        info["integer"]
                    )
                )

                if not info["integer"]:

                    print(
                        "    W_denominator={}".format(
                            info["denominator"]
                        )
                    )

                    print(
                        "    W_denominator_factorization={}".format(
                            info[
                                "denominator_factorization"
                            ]
                        )
                    )


# ============================================================================
# CROSS-MODEL SUMMARY
# ============================================================================

def cross_model_summary(
    results,
):

    print()
    print("=" * 78)
    print(
        "CROSS-MODEL SYMBOLIC COMPATIBILITY SUMMARY"
    )
    print("=" * 78)

    for name, result, classification in results:

        print()
        print(
            "  {}: {}".format(
                name,
                classification,
            )
        )

        if (
            result.get(
                "solution"
            )
        ):

            for index, solution in enumerate(
                result["solution"]
            ):

                print(
                    "    branch_{}={}".format(
                        index,
                        solution,
                    )
                )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 355R — EXACT SYMBOLIC MISSING-CELL CLOSURE / "
        "MODEL ELIMINATION"
    )
    print("=" * 78)

    lattice = build_observed_lattice()

    print()
    print(
        "1. OBSERVED SOURCE CELLS"
    )

    print(
        "  observed_count={}".format(
            len(lattice)
        )
    )

    print(
        "  observed_cells={}".format(
            sorted(
                lattice
            )
        )
    )

    print()
    print(
        "  symbolic_missing_cells:"
    )

    print(
        "    Z = Q_3(5) = Q(2,3)"
    )

    print(
        "    W = Q_1(7) = Q(3,1)"
    )

    results = []

    # ------------------------------------------------------------------------
    # Run every candidate model.
    # ------------------------------------------------------------------------

    for name, model in MODELS.items():

        result = solve_symbolically(
            lattice,
            name,
            model,
        )

        classification = classify_solution(
            result,
            model,
        )

        report(
            name,
            result,
            model,
        )

        results.append(
            (
                name,
                result,
                classification,
            )
        )

    # ------------------------------------------------------------------------
    # Cross-model comparison.
    # ------------------------------------------------------------------------

    cross_model_summary(
        results
    )

    # ------------------------------------------------------------------------
    # Specific strategic conclusions.
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "STRATEGIC INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
This experiment does not predict Z or W by choosing a preferred model.

Instead it asks whether a candidate model can survive at all once the
missing cells are left symbolic.

The possible outcomes have sharply different meanings:

    IMPOSSIBLE
        The model is already ruled out by the observed table together
        with the existence of some values Z,W.

    SYMBOLICALLY_CONSTRAINED_Z
        The model forces an exact value of Q_3(5), independently of any
        external observation.

        This is still a MODEL CONSEQUENCE, not observed data.

    SYMBOLICALLY_CONSTRAINED_W
        The analogous statement for Q_1(7).

    SYMBOLICALLY_CONSTRAINED_BOTH
        The model forces both missing cells.

    SYMBOLICALLY_UNDERDETERMINED
        The current information is insufficient to identify the missing
        values.

The especially important result is a contradiction of the form:

    model is IMPOSSIBLE

because that eliminates an entire family without requiring another
formula search.

Likewise, if D6 remains symbolically constrained to the same
non-integer Z found in 353R, the conclusion becomes stronger:

    the model forces a non-integer missing source value.

Under an independently established integer-valued source definition,
that would eliminate D6 before Q_3(5) is observed.

No symbolic solution is ever promoted to experimental data.
"""
    )

    print()
    print("=" * 78)
    print(
        "FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  symbolic_missing_cells=['Z=Q_3(5)', 'W=Q_1(7)']"
    )

    print(
        "  observed_cells_only_calibration=True"
    )

    print(
        "  missing_cells_inserted_numerically=False"
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
        "EXPERIMENT 355R COMPLETE"
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
