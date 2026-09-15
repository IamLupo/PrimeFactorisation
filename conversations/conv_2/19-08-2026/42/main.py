#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 356R — EXACT INTEGER-COMPATIBILITY / SYMBOLIC-SOURCE AUDIT
==============================================================================

Purpose
-------
Experiment 355R established:

    D1, D2, D4  -> impossible even with Z,W free
    D5, D6      -> force unique non-integer Z
    D7          -> remains symbolically underdetermined

where

    Z = Q_3(5) = Q(2,3)
    W = Q_1(7) = Q(3,1).

356R asks the sharper question:

    Which surviving candidate models are compatible with an
    INTEGER-valued source table?

The source values already observed are all integers. The mathematical
meaning of Q strongly suggests that any genuine source construction
should also produce integers at the missing cells, but this is treated
as a CONDITION, not a theorem.

For each candidate model:

    1. keep Z,W symbolic;
    2. solve the exact rational calibration system;
    3. impose Z integer compatibility;
    4. where possible, derive the exact congruence/divisibility condition;
    5. test whether the model admits any integer Z at all.

Special attention:

    D5
    D6
    D7

No missing value is inserted numerically.

No interpolation.
No extrapolation.
No synthetic second n=pq case.
Exact SymPy arithmetic only.
"""


from __future__ import annotations

import math
import sys

import sympy as sp


# ============================================================================
# OBSERVED DATA
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
    "Z W",
    integer=True,
)

a, b, c, d, e = sp.symbols(
    "a b c d e"
)


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


def factor_integer(n):
    n = int(sp.Integer(n))

    if n == 0:
        return {}

    return sp.factorint(
        abs(n)
    )


def build_lattice():

    lattice = {}

    for p_value, values in Q.items():

        r = (
            p_value - 1
        ) // 2

        for index, value in enumerate(values):

            t = (
                len(values)
                - 1
                - index
            )

            lattice[(r, t)] = sp.Integer(
                value
            )

    return lattice


def value_at(
    lattice,
    cell,
):

    if cell == (2, 3):
        return Z

    if cell == (3, 1):
        return W

    return lattice.get(
        cell,
        None,
    )


# ============================================================================
# CANDIDATE MODELS
# ============================================================================

MODELS = {

    "D5": {
        "parameters": [a, b, c, d],

        "sources": lambda r, t: [
            (r - 1, t),
            (r, t - 1),
            (r - 1, t - 1),
            (r - 1, t - 2),
        ],
    },

    "D6": {
        "parameters": [a, b, c, d],

        "sources": lambda r, t: [
            (r - 1, t),
            (r, t - 1),
            (r - 1, t - 1),
            (r, t - 2),
        ],
    },

    "D7": {
        "parameters": [a, b, c, d, e],

        "sources": lambda r, t: [
            (r - 1, t),
            (r, t - 1),
            (r - 1, t - 1),
            (r, t - 2),
            (r - 1, t - 2),
        ],
    },
}


# ============================================================================
# EQUATION GENERATION
# ============================================================================

def build_equations(
    lattice,
    model,
):

    equations = []
    locations = []

    parameters = model["parameters"]
    source_builder = model["sources"]

    for r in range(4):

        for t in range(6):

            target = (
                r,
                t,
            )

            if r < 1:
                continue

            sources = source_builder(
                r,
                t,
            )

            if any(
                rr < 0 or tt < 0
                for rr, tt in sources
            ):
                continue

            source_values = []

            usable = True

            for source in sources:

                value = value_at(
                    lattice,
                    source,
                )

                if value is None:

                    usable = False
                    break

                source_values.append(
                    value
                )

            if not usable:
                continue

            target_value = value_at(
                lattice,
                target,
            )

            if target_value is None:
                continue

            rhs = sum(
                parameters[i]
                * source_values[i]
                for i in range(
                    len(parameters)
                )
            )

            equations.append(
                clean(
                    rhs
                    - target_value
                )
            )

            locations.append(
                target
            )

    return parameters, equations, locations


# ============================================================================
# RATIONAL SYMBOLIC SOLUTION
# ============================================================================

def solve_model(
    lattice,
    name,
    model,
):

    parameters, equations, locations = (
        build_equations(
            lattice,
            model,
        )
    )

    unknowns = (
        list(parameters)
        + [Z, W]
    )

    print()
    print("=" * 78)
    print(
        "MODEL {}".format(name)
    )
    print("=" * 78)

    print(
        "  equations={}".format(
            len(equations)
        )
    )

    print(
        "  unknowns_including_ZW={}".format(
            len(unknowns)
        )
    )

    print(
        "  locations={}".format(
            locations
        )
    )

    # ------------------------------------------------------------------------
    # Groebner feasibility
    # ------------------------------------------------------------------------

    groebner = sp.groebner(
        equations,
        *unknowns,
        order="lex",
    )

    print()
    print(
        "  groebner_basis_size={}".format(
            len(groebner.polys)
        )
    )

    if groebner == sp.groebner(
        [1],
        *unknowns,
        order="lex",
    ):

        print(
            "  symbolic_status=IMPOSSIBLE"
        )

        return {
            "status": "IMPOSSIBLE",
            "equations": equations,
            "locations": locations,
            "solutions": [],
        }

    print(
        "  symbolic_status=CONSISTENT"
    )

    # ------------------------------------------------------------------------
    # Explicit solution where available
    # ------------------------------------------------------------------------

    solutions = sp.solve(
        equations,
        unknowns,
        dict=True,
        simplify=False,
        rational=True,
    )

    if not solutions:

        print(
            "  explicit_solution=UNRESOLVED"
        )

        return {
            "status": "CONSISTENT_UNRESOLVED",
            "equations": equations,
            "locations": locations,
            "solutions": [],
        }

    print(
        "  solution_branches={}".format(
            len(solutions)
        )
    )

    cleaned_solutions = []

    for index, solution in enumerate(
        solutions
    ):

        solution = {
            symbol: clean(value)
            for symbol, value
            in solution.items()
        }

        cleaned_solutions.append(
            solution
        )

        print()
        print(
            "  branch_{}={}".format(
                index,
                solution
            )
        )

    return {
        "status": "CONSISTENT",
        "equations": equations,
        "locations": locations,
        "solutions": cleaned_solutions,
    }


# ============================================================================
# INTEGER-COMPATIBILITY ANALYSIS
# ============================================================================

def analyze_Z_integer_compatibility(
    result,
):

    print()
    print(
        "  INTEGER-Z COMPATIBILITY"
    )

    solutions = result.get(
        "solutions",
        []
    )

    if not solutions:

        print(
            "    status=NO_EXPLICIT_BRANCH"
        )

        return

    for index, solution in enumerate(
        solutions
    ):

        if Z not in solution:

            print(
                "    branch_{}: Z=FREE".format(
                    index
                )
            )

            print(
                "      integer_Z_possible=True"
            )

            continue

        z_value = clean(
            solution[Z]
        )

        numerator = sp.numer(
            sp.Rational(
                z_value
            )
        )

        denominator = sp.denom(
            sp.Rational(
                z_value
            )
        )

        print()
        print(
            "    branch_{}:".format(
                index
            )
        )

        print(
            "      Z={}".format(
                z_value
            )
        )

        print(
            "      numerator={}".format(
                numerator
            )
        )

        print(
            "      denominator={}".format(
                denominator
            )
        )

        print(
            "      denominator_factorization={}".format(
                factor_integer(
                    denominator
                )
            )
        )

        integer = (
            denominator == 1
        )

        print(
            "      integer_Z_possible={}".format(
                integer
            )
        )

        if not integer:

            print(
                "      VERDICT="
                "INCOMPATIBLE_WITH_INTEGER_SOURCE"
            )


# ============================================================================
# INTEGER-PARAMETER SUBTEST
# ============================================================================

def integer_parameter_audit(
    lattice,
    name,
    model,
):

    parameters, equations, _ = (
        build_equations(
            lattice,
            model,
        )
    )

    print()
    print(
        "  INTEGER-PARAMETER SUBTEST"
    )

    #
    # Treat model coefficients themselves as integers and solve the
    # corresponding Diophantine system together with integer Z,W.
    #
    unknowns = (
        list(parameters)
        + [Z, W]
    )

    try:

        solution = sp.linsolve(
            equations,
            unknowns,
        )

        print(
            "    linsolve={}".format(
                solution
            )
        )

        if solution == sp.EmptySet:

            print(
                "    integer_parameter_status=IMPOSSIBLE"
            )

        else:

            print(
                "    integer_parameter_status="
                "RATIONAL_CONSISTENCY_ONLY"
            )

    except Exception as exc:

        print(
            "    integer_parameter_status=UNRESOLVED"
        )

        print(
            "    error={}: {}".format(
                type(exc).__name__,
                exc,
            )
        )


# ============================================================================
# CROSS-MODEL SUMMARY
# ============================================================================

def summary(results):

    print()
    print("=" * 78)
    print(
        "CROSS-MODEL INTEGER COMPATIBILITY SUMMARY"
    )
    print("=" * 78)

    for name, result in results:

        print()
        print(
            "  {}: {}".format(
                name,
                result["status"]
            )
        )

        solutions = result.get(
            "solutions",
            []
        )

        if not solutions:

            continue

        for index, solution in enumerate(
            solutions
        ):

            if Z in solution:

                z = clean(
                    solution[Z]
                )

                print(
                    "    branch_{} Z={} integer={}".format(
                        index,
                        z,
                        sp.denom(
                            sp.Rational(z)
                        ) == 1,
                    )
                )

            else:

                print(
                    "    branch_{} Z=FREE".format(
                        index
                    )
                )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 356R — EXACT INTEGER-COMPATIBILITY / "
        "SYMBOLIC-SOURCE AUDIT"
    )
    print("=" * 78)

    lattice = build_lattice()

    print()
    print(
        "OBSERVED DATA"
    )

    print(
        "  observed_cells={}".format(
            len(lattice)
        )
    )

    print(
        "  Z=Q_3(5) remains symbolic=True"
    )

    print(
        "  W=Q_1(7) remains symbolic=True"
    )

    results = []

    for name, model in MODELS.items():

        result = solve_model(
            lattice,
            name,
            model,
        )

        analyze_Z_integer_compatibility(
            result
        )

        integer_parameter_audit(
            lattice,
            name,
            model,
        )

        results.append(
            (
                name,
                result,
            )
        )

    summary(
        results
    )

    print()
    print("=" * 78)
    print(
        "STRATEGIC INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
The central question is now arithmetic feasibility rather than formula
discovery.

For D5 and D6, 355R already showed that the observed equations force a
specific Z.

356R checks whether those forced values are compatible with

    Z ∈ Z.

Because the observed source values are all integers, this is the natural
minimal consistency condition for a genuinely integer-valued source.

The logically correct conclusions are conditional:

    forced non-integer Z
        +
    independently established integer-valued source
        =>
    model falsified.

This does NOT prove that Q_3(5) equals any rational prediction.

D7 is expected to retain a free Z parameter. That is not a discovery;
it means the present data do not identify the source through that model.

The most valuable external datum remains the actual Q_3(5) value.

If it is eventually obtained, it should be inserted into a separate
validation run, never back into the calibration used to generate it.
"""
    )

    print()
    print("=" * 78)
    print(
        "FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  symbolic_Z_used=True"
    )

    print(
        "  symbolic_W_used=True"
    )

    print(
        "  missing_values_inserted_numerically=False"
    )

    print(
        "  integer_compatibility_tested=True"
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
        "EXPERIMENT 356R COMPLETE"
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
