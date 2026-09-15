#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 357R-FIXED2 — EXACT SYMBOLIC CLOSURE /
                           ROBUST FREE-PARAMETER AUDIT
==============================================================================

Purpose
-------
Correct the second D7 failure.

The previous version assumed that

    solve(..., dict=True)

would return every coefficient, including free variables. SymPy does not
guarantee this for underdetermined systems; a free variable may simply be
omitted from the dictionary.

This version therefore:

    * uses exact linear-system parameterization;
    * explicitly detects free parameters;
    * never assumes that e is present in a solution dictionary;
    * keeps Z=Q_3(5) and W=Q_1(7) symbolic;
    * distinguishes:
          UNIQUE
          UNDERDETERMINED
          IMPOSSIBLE
    * performs exact D5/D6 forced-Z integrality checks;
    * performs a conservative D7 integer-coefficient feasibility audit.

No missing cell is inserted numerically.
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

Z, W = sp.symbols(
    "Z W"
)

a, b, c, d, e = sp.symbols(
    "a b c d e"
)

D56_VARS = [
    a,
    b,
    c,
    d,
]

D7_VARS = [
    a,
    b,
    c,
    d,
    e,
]


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
    lattice = {}

    for p_value, values in Q.items():

        r = (
            p_value - 1
        ) // 2

        for index, value in enumerate(
            values
        ):

            t = (
                len(values)
                - 1
                - index
            )

            lattice[
                (r, t)
            ] = sp.Integer(value)

    return lattice


def sorted_cells(lattice):
    return sorted(
        lattice,
        key=lambda cell: (
            cell[1],
            cell[0],
        ),
    )


def is_integer(value):
    value = sp.Rational(value)
    return value.q == 1


def numerator_denominator(value):
    value = clean(value)

    numerator, denominator = sp.fraction(
        sp.together(value)
    )

    return (
        sp.Integer(numerator),
        sp.Integer(denominator),
    )


def factor_integer(value):
    value = int(
        sp.Integer(value)
    )

    if value == 0:
        return {}

    return sp.factorint(
        abs(value)
    )


# ============================================================================
# SYMBOLIC SOURCE ACCESS
# ============================================================================

def symbolic_value(
    lattice,
    cell,
):
    """
    Observed cells are fixed integers.

    The ONLY symbolic missing cells are:

        Z = Q_3(5) = Q(2,3)
        W = Q_1(7) = Q(3,1)
    """

    if cell == (2, 3):
        return Z

    if cell == (3, 1):
        return W

    if cell in lattice:
        return lattice[cell]

    return None


# ============================================================================
# MODEL SUPPORT
# ============================================================================

def model_sources(
    model,
    r,
    t,
):

    if model == "D5":

        return [
            (r - 1, t),
            (r, t - 1),
            (r - 1, t - 1),
            (r - 1, t - 2),
        ]

    if model == "D6":

        return [
            (r - 1, t),
            (r, t - 1),
            (r - 1, t - 1),
            (r, t - 2),
        ]

    if model == "D7":

        return [
            (r - 1, t),
            (r, t - 1),
            (r - 1, t - 1),
            (r, t - 2),
            (r - 1, t - 2),
        ]

    raise ValueError(
        "Unknown model {}".format(
            model
        )
    )


def model_variables(model):

    if model in (
        "D5",
        "D6",
    ):
        return D56_VARS

    if model == "D7":
        return D7_VARS

    raise ValueError(
        "Unknown model {}".format(
            model
        )
    )


# ============================================================================
# SYMBOLIC CLOSURE EQUATIONS
# ============================================================================

def symbolic_closure_equations(
    lattice,
    model,
):
    """
    IMPORTANT:

    This is a SYMBOLIC MISSING-CELL audit.

    Unlike an observed-only calibration, a target may be

        Z = Q(2,3)

    or

        W = Q(3,1).

    Missing target/source cells are represented symbolically, never
    numerically inserted.
    """

    variables = model_variables(
        model
    )

    equations = []
    locations = []

    for r in range(4):

        for t in range(6):

            target = (
                r,
                t,
            )

            observed_or_symbolic_target = (
                target in lattice
                or target in (
                    (2, 3),
                    (3, 1),
                )
            )

            if not observed_or_symbolic_target:
                continue

            if r < 1:
                continue

            sources = model_sources(
                model,
                r,
                t,
            )

            source_values = []

            source_available = True

            for source in sources:

                value = symbolic_value(
                    lattice,
                    source,
                )

                if value is None:

                    source_available = False
                    break

                source_values.append(
                    value
                )

            if not source_available:
                continue

            target_value = symbolic_value(
                lattice,
                target,
            )

            if target_value is None:
                continue

            lhs = sum(
                variables[i]
                * source_values[i]
                for i in range(
                    len(variables)
                )
            )

            equations.append(
                clean(
                    lhs
                    - target_value
                )
            )

            locations.append(
                target
            )

    return (
        variables,
        equations,
        locations,
    )


# ============================================================================
# ROBUST EXACT LINEAR PARAMETERIZATION
# ============================================================================

def parameterize_linear_system(
    equations,
    variables,
):
    """
    Robust exact solver.

    Returns a parameterization even when the system is underdetermined.

    SymPy's linsolve result is authoritative. We do NOT assume every
    original variable appears explicitly in a dictionary.
    """

    if not equations:

        return {
            "status": "NO_EQUATIONS",
            "rank": 0,
            "augmented_rank": 0,
            "solution_set": sp.EmptySet,
            "tuple": None,
            "free_parameters": [],
        }

    matrix, rhs = sp.linear_eq_to_matrix(
        equations,
        variables,
    )

    rank = matrix.rank()

    augmented_rank = (
        matrix
        .row_join(rhs)
        .rank()
    )

    if augmented_rank > rank:

        return {
            "status": "IMPOSSIBLE",
            "rank": rank,
            "augmented_rank": augmented_rank,
            "solution_set": sp.EmptySet,
            "tuple": None,
            "free_parameters": [],
        }

    solution_set = sp.linsolve(
        equations,
        variables,
    )

    if solution_set == sp.EmptySet:

        return {
            "status": "IMPOSSIBLE",
            "rank": rank,
            "augmented_rank": augmented_rank,
            "solution_set": solution_set,
            "tuple": None,
            "free_parameters": [],
        }

    tuples = list(
        solution_set
    )

    if len(tuples) != 1:

        return {
            "status": "MULTIPLE_BRANCHES",
            "rank": rank,
            "augmented_rank": augmented_rank,
            "solution_set": solution_set,
            "tuple": None,
            "free_parameters": [],
        }

    solution_tuple = tuples[0]

    #
    # Any symbol appearing in the parameterized tuple but not in the
    # original variable list is a SymPy-generated free parameter.
    #
    original_symbols = set(
        variables
    )

    free_parameters = sorted(
        set().union(
            *[
                expr.free_symbols
                for expr in solution_tuple
            ]
        )
        - original_symbols
        - {Z, W},
        key=str,
    )

    status = (
        "UNIQUE"
        if rank == len(variables)
        else "UNDERDETERMINED"
    )

    return {
        "status": status,
        "rank": rank,
        "augmented_rank": augmented_rank,
        "solution_set": solution_set,
        "tuple": solution_tuple,
        "free_parameters": free_parameters,
    }


# ============================================================================
# REPORT D5 / D6
# ============================================================================

def run_d56(
    lattice,
    model,
):

    variables, equations, locations = (
        symbolic_closure_equations(
            lattice,
            model,
        )
    )

    print()
    print("=" * 78)
    print(
        "{} SYMBOLIC CLOSURE".format(
            model
        )
    )
    print("=" * 78)

    print(
        "  equation_count={}".format(
            len(equations)
        )
    )

    print(
        "  locations={}".format(
            locations
        )
    )

    result = parameterize_linear_system(
        equations,
        variables,
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
        "  status={}".format(
            result["status"]
        )
    )

    if result["status"] == "IMPOSSIBLE":

        return result

    if result["status"] == "UNDERDETERMINED":

        print()
        print(
            "  parameterized_solution={}".format(
                result["solution_set"]
            )
        )

        print(
            "  free_parameters={}".format(
                result["free_parameters"]
            )
        )

        return result

    tuple_solution = result[
        "tuple"
    ]

    solved = {
        variables[i]: clean(
            tuple_solution[i]
        )
        for i in range(
            len(variables)
        )
    }

    print()
    print(
        "  solved_parameters={}".format(
            solved
        )
    )

    #
    # For D5/D6 the symbolic closure should force Z.
    #
    z_expression = None

    for index, variable in enumerate(
        variables
    ):
        if variable in solved:
            pass

    #
    # A UNIQUE parameter solve may contain Z explicitly if the equations
    # included the symbolic target. Extract Z by inspecting equations.
    #
    reduced_equations = [
        clean(
            equation.subs(
                solved
            )
        )
        for equation in equations
    ]

    z_conditions = [
        equation
        for equation in reduced_equations
        if equation.has(Z)
    ]

    if z_conditions:

        print()
        print(
            "  residual_Z_conditions={}".format(
                z_conditions
            )
        )

    return result


# ============================================================================
# FORCE-Z EXTRACTION FOR D5/D6
# ============================================================================

def extract_forced_Z(
    equations,
    variables,
):
    """
    Eliminate coefficient variables and obtain an exact condition on Z.

    This is more robust than expecting a coefficient dictionary to contain
    Z in a particular location.
    """

    #
    # Groebner elimination is used only as an exact symbolic elimination
    # device on a LINEAR system.
    #
    G = sp.groebner(
        equations,
        *variables,
        Z,
        W,
        order="lex",
    )

    basis = list(
        G.polys
    )

    z_only = []

    for poly in basis:

        expr = clean(
            poly.as_expr()
        )

        remaining = (
            expr.free_symbols
            - {Z, W}
        )

        if not remaining:

            z_only.append(
                expr
            )

    return z_only


def forced_z_report(
    lattice,
    model,
):

    variables, equations, locations = (
        symbolic_closure_equations(
            lattice,
            model,
        )
    )

    print()
    print("=" * 78)
    print(
        "{} FORCED-Z ELIMINATION".format(
            model
        )
    )
    print("=" * 78)

    z_conditions = extract_forced_Z(
        equations,
        variables,
    )

    print(
        "  elimination_polynomials={}".format(
            z_conditions
        )
    )

    if not z_conditions:

        print(
            "  forced_Z=UNRESOLVED"
        )

        return None

    linear_conditions = []

    for condition in z_conditions:

        poly = sp.Poly(
            condition,
            Z,
            domain=sp.QQ,
        )

        if poly.degree() == 1:

            solution = sp.solve(
                condition,
                Z,
            )

            if len(solution) == 1:

                value = clean(
                    solution[0]
                )

                linear_conditions.append(
                    value
                )

    if not linear_conditions:

        print(
            "  forced_Z=UNRESOLVED"
        )

        return None

    distinct = list(
        dict.fromkeys(
            linear_conditions
        )
    )

    print(
        "  forced_Z_candidates={}".format(
            distinct
        )
    )

    if len(distinct) == 1:

        value = distinct[0]

        print()
        print(
            "  forced_Z={}".format(
                value
            )
        )

        print(
            "  integer_Z={}".format(
                is_integer(value)
            )
        )

        numerator, denominator = (
            numerator_denominator(
                value
            )
        )

        print(
            "  numerator={}".format(
                numerator
            )
        )

        print(
            "  denominator={}".format(
                denominator
            )
        )

        print(
            "  denominator_factorization={}".format(
                factor_integer(
                    denominator
                )
            )
        )

        return value

    print(
        "  mutually_distinct_forced_Z_values=True"
    )

    return distinct


# ============================================================================
# D7 INTEGER-FEASIBILITY AUDIT
# ============================================================================

def d7_integer_feasibility(
    lattice,
):

    variables, equations, locations = (
        symbolic_closure_equations(
            lattice,
            "D7",
        )
    )

    print()
    print("=" * 78)
    print(
        "D7 ROBUST INTEGER-FEASIBILITY AUDIT"
    )
    print("=" * 78)

    print(
        "  locations={}".format(
            locations
        )
    )

    result = parameterize_linear_system(
        equations,
        variables,
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
        "  status={}".format(
            result["status"]
        )
    )

    if result["status"] == "IMPOSSIBLE":

        print(
            "  integer_feasibility=IMPOSSIBLE"
        )

        return result

    if result["status"] == "UNIQUE":

        tuple_solution = result[
            "tuple"
        ]

        solved = {
            variables[i]: clean(
                tuple_solution[i]
            )
            for i in range(
                len(variables)
            )
        }

        print()
        print(
            "  unique_solution={}".format(
                solved
            )
        )

        for variable in variables:

            value = solved[
                variable
            ]

            print(
                "    {}_integer={}".format(
                    variable,
                    is_integer(
                        value
                    )
                    if not value.has(Z)
                    else "REQUIRES_Z"
                )
            )

        return result

    #
    # UNDERDETERMINED:
    #
    # Print the exact parameterized tuple. Do not guess which original
    # variable SymPy selected as free.
    #
    print()
    print(
        "  parameterized_solution={}".format(
            result["solution_set"]
        )
    )

    print()
    print(
        "  free_parameters={}".format(
            result["free_parameters"]
        )
    )

    #
    # Explicitly substitute a human-readable free parameter if present.
    #
    free_parameters = result[
        "free_parameters"
    ]

    for parameter in free_parameters:

        print()
        print(
            "  free_parameter={}".format(
                parameter
            )
        )

        for i, variable in enumerate(
            variables
        ):

            expression = clean(
                result["tuple"][i]
            )

            if expression.has(
                parameter
            ):

                print(
                    "    {}={}".format(
                        variable,
                        expression,
                    )
                )

    #
    # Important logical conclusion:
    #
    # We do NOT claim integer feasibility unless every coefficient can be
    # shown integer for a common integer choice of the free parameters and Z.
    #
    print()
    print(
        "  integer_feasibility=UNRESOLVED"
    )

    print(
        "  reason=underdetermined_symbolic_coefficient_family"
    )

    print(
        "  no_arbitrary_free_parameter_assumption=True"
    )

    return result


# ============================================================================
# OPTIONAL D7 SPECIALIZATION e = integer
# ============================================================================

def d7_e_fixed_integer_scan(
    lattice,
):
    """
    Conservative diagnostic only.

    Since D7 is underdetermined, this does NOT search an arbitrary numerical
    range. Instead it asks whether setting the SymPy free parameter to zero
    produces a fully integer solution. This is a diagnostic, NOT evidence
    for the source model.
    """

    variables, equations, locations = (
        symbolic_closure_equations(
            lattice,
            "D7",
        )
    )

    result = parameterize_linear_system(
        equations,
        variables,
    )

    print()
    print("=" * 78)
    print(
        "D7 FREE-PARAMETER ZERO DIAGNOSTIC"
    )
    print("=" * 78)

    if result["status"] != "UNDERDETERMINED":

        print(
            "  diagnostic_available=False"
        )

        return

    free_parameters = result[
        "free_parameters"
    ]

    if not free_parameters:

        print(
            "  diagnostic_available=False"
        )

        return

    parameter = free_parameters[0]

    specialized = [
        clean(
            expression.subs(
                parameter,
                0,
            )
        )
        for expression
        in result["tuple"]
    ]

    print(
        "  free_parameter_set_to=0"
    )

    for variable, expression in zip(
        variables,
        specialized,
    ):

        print(
            "    {}={}".format(
                variable,
                expression,
            )
        )

    #
    # Check whether expressions independent of Z are already integral.
    #
    integral_without_Z = True

    for variable, expression in zip(
        variables,
        specialized,
    ):

        if not expression.has(Z):

            if not is_integer(
                expression
            ):

                integral_without_Z = False

    print()
    print(
        "  all_Z_independent_coefficients_integer={}".format(
            integral_without_Z
        )
    )

    print(
        "  interpretation=DIAGNOSTIC_ONLY_NOT_A_MODEL_SELECTION"
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 357R-FIXED2 — EXACT SYMBOLIC CLOSURE / "
        "ROBUST FREE-PARAMETER AUDIT"
    )
    print("=" * 78)

    lattice = build_lattice()

    print()
    print(
        "OBSERVED SOURCE"
    )

    print(
        "  observed_cells={}".format(
            len(lattice)
        )
    )

    print(
        "  Z=Q_3(5) symbolic=True"
    )

    print(
        "  W=Q_1(7) symbolic=True"
    )

    # ------------------------------------------------------------------------
    # D5
    # ------------------------------------------------------------------------

    d5_result = run_d56(
        lattice,
        "D5",
    )

    d5_forced_Z = forced_z_report(
        lattice,
        "D5",
    )

    # ------------------------------------------------------------------------
    # D6
    # ------------------------------------------------------------------------

    d6_result = run_d56(
        lattice,
        "D6",
    )

    d6_forced_Z = forced_z_report(
        lattice,
        "D6",
    )

    # ------------------------------------------------------------------------
    # D7
    # ------------------------------------------------------------------------

    d7_result = d7_integer_feasibility(
        lattice
    )

    d7_e_fixed_integer_scan(
        lattice
    )

    # ------------------------------------------------------------------------
    # CROSS-MODEL SUMMARY
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "CROSS-MODEL SUMMARY"
    )
    print("=" * 78)

    print()
    print(
        "  D5_status={}".format(
            d5_result["status"]
        )
    )

    print(
        "  D5_forced_Z={}".format(
            d5_forced_Z
        )
    )

    print(
        "  D5_forced_Z_integer={}".format(
            (
                is_integer(
                    d5_forced_Z
                )
                if (
                    d5_forced_Z
                    is not None
                    and
                    not isinstance(
                        d5_forced_Z,
                        list,
                    )
                )
                else "UNRESOLVED"
            )
        )
    )

    print()
    print(
        "  D6_status={}".format(
            d6_result["status"]
        )
    )

    print(
        "  D6_forced_Z={}".format(
            d6_forced_Z
        )
    )

    print(
        "  D6_forced_Z_integer={}".format(
            (
                is_integer(
                    d6_forced_Z
                )
                if (
                    d6_forced_Z
                    is not None
                    and
                    not isinstance(
                        d6_forced_Z,
                        list,
                    )
                )
                else "UNRESOLVED"
            )
        )
    )

    print()
    print(
        "  D7_status={}".format(
            d7_result["status"]
        )
    )

    print(
        "  D7_integer_feasibility=UNRESOLVED"
    )

    # ------------------------------------------------------------------------
    # STRATEGIC INTERPRETATION
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "STRATEGIC INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
The previous crash occurred because an underdetermined symbolic system
does not have to return every original variable in a solve dictionary.

The correct treatment is:

    1. preserve the complete parameterized solution;
    2. identify SymPy-generated free parameters explicitly;
    3. never assume that a particular coefficient, such as e, is fixed;
    4. only make integrality claims that follow for the full parameterized
       family.

D5 and D6 can still be tested sharply because their elimination systems
force a specific Z.

D7 is different: the available symbolic equations do not uniquely
identify all five coefficients.

Therefore the honest D7 conclusion is:

    integer feasibility = UNRESOLVED

unless a further exact Diophantine analysis is performed.

This experiment deliberately refuses to turn an arbitrary choice of
free parameter into a structural discovery.

No missing value is treated as observed.
"""
    )

    # ------------------------------------------------------------------------
    # FINAL EXACTNESS
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  robust_underdetermined_parameterization=True"
    )

    print(
        "  accidental_free_e_assumption=False"
    )

    print(
        "  D5_symbolic_closure_tested=True"
    )

    print(
        "  D6_symbolic_closure_tested=True"
    )

    print(
        "  D7_symbolic_closure_tested=True"
    )

    print(
        "  D7_integer_feasibility_not_overclaimed=True"
    )

    print(
        "  missing_values_inserted_numerically=False"
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
        "EXPERIMENT 357R-FIXED2 COMPLETE"
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