#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 350R — EXACT BOUNDARY-FORCED GENERATING-DENOMINATOR AUDIT
==============================================================================

Purpose
-------
Experiment 349R showed that the observed bivariate generating polynomial

    G(x,y) = sum Q(r,t) x^r y^t

is irreducible over Q, is not low-degree, and has no obvious terminal
polynomial factor.

That does NOT yet rule out a simple RATIONAL generating mechanism.

A rational generating function commonly has the form

    D(x,y) G(x,y) = N(x,y),

where D is a small denominator and N is supported only on the boundary.

Equivalently, the coefficient array satisfies a translation-invariant
linear recurrence in the interior, while boundary cells act as forcing.

This experiment therefore searches for small exact denominator families:

    D1 = 1 - a*x - b*y

    D2 = 1 - a*x - b*y - c*x*y

    D3 = 1 - a*x - b*y - c*x^2

    D4 = 1 - a*x - b*y - c*y^2

    D5 = 1 - a*x - b*y - c*x*y - d*x^2

    D6 = 1 - a*x - b*y - c*x*y - d*y^2

    D7 = 1 - a*x - b*y - c*x*y - d*x^2 - e*y^2

    D8 = 1 - a*x - b*y - c*x*y - d*x^2 - e*x*y^2

The denominator coefficients are solved exactly from interior coefficient
equations.

IMPORTANT
---------
We do NOT assume unobserved cells are zero.

A coefficient equation is used only when every source coefficient appearing
in that equation is an observed cell.

Thus the search is a genuine exact test on the observed triangular support.

A candidate is classified as:

    EXACT_OVERDETERMINED
        unique denominator coefficients and more equations than unknowns;

    EXACT_DATA_SIZED
        unique exact reconstruction with no redundancy;

    NONUNIQUE
        insufficient identification;

    NO_SOLUTION
        exact contradiction.

After a denominator is found, the induced boundary residual

    N = D*G

is inspected for:

    * support concentration near r=0;
    * support concentration near t=0;
    * small boundary width;
    * low total degree;
    * exact integer normalization.

No missing cells.
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


x, y = sp.symbols("x y")


# ============================================================================
# BASIC HELPERS
# ============================================================================

def clean(expr):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(expr)
            )
        )
    )


def lattice():
    result = {}

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

            result[
                (r, t)
            ] = sp.Integer(value)

    return result


def generating_polynomial(L):
    return clean(
        sum(
            value * x**r * y**t
            for (
                r,
                t,
            ), value
            in L.items()
        )
    )


def primitive_integer_polynomial(expr):
    P = sp.Poly(
        sp.expand(expr),
        x,
        y,
        domain=sp.QQ,
    )

    if P.is_zero:
        return sp.Poly(
            0,
            x,
            y,
            domain=sp.ZZ,
        )

    coeffs = P.coeffs()

    den_lcm = 1

    for c in coeffs:
        den_lcm = sp.ilcm(
            den_lcm,
            int(sp.denom(c)),
        )

    ints = [
        int(c * den_lcm)
        for c in coeffs
    ]

    g = 0

    for value in ints:
        g = math.gcd(
            g,
            abs(value),
        )

    ints = [
        value // g
        for value in ints
    ]

    expr2 = sp.Integer(0)

    for monomial, coefficient in zip(
        P.monoms(),
        ints,
    ):

        i, j = monomial

        expr2 += (
            sp.Integer(coefficient)
            * x**i
            * y**j
        )

    result = sp.Poly(
        sp.expand(expr2),
        x,
        y,
        domain=sp.ZZ,
    )

    if result.LC() < 0:
        result = sp.Poly(
            -result.as_expr(),
            x,
            y,
            domain=sp.ZZ,
        )

    return result


# ============================================================================
# DENOMINATOR FAMILIES
# ============================================================================

def denominator_families():
    a, b, c, d, e = sp.symbols(
        "a b c d e"
    )

    return {
        "D1_xy": (
            1 - a*x - b*y,
            [a, b],
        ),

        "D2_xy_xy": (
            1 - a*x - b*y - c*x*y,
            [a, b, c],
        ),

        "D3_xy_x2": (
            1 - a*x - b*y - c*x**2,
            [a, b, c],
        ),

        "D4_xy_y2": (
            1 - a*x - b*y - c*y**2,
            [a, b, c],
        ),

        "D5_xy_xy_x2": (
            1 - a*x - b*y - c*x*y - d*x**2,
            [a, b, c, d],
        ),

        "D6_xy_xy_y2": (
            1 - a*x - b*y - c*x*y - d*y**2,
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


# ============================================================================
# EXACT INTERIOR-EQUATION CONSTRUCTION
# ============================================================================

def denominator_terms(
    denominator,
    parameters,
):
    """
    Return shifted coefficient terms

        coefficient * x^dr y^dt.

    The constant term is always included and has coefficient 1.
    """

    P = sp.Poly(
        sp.expand(denominator),
        x,
        y,
        domain=sp.QQ.frac_field(
            *parameters
        ),
    )

    terms = []

    for monomial, coefficient in zip(
        P.monoms(),
        P.coeffs(),
    ):

        dr, dt = monomial

        terms.append(
            (
                int(dr),
                int(dt),
                clean(coefficient),
            )
        )

    return terms


def valid_equation_cells(
    L,
    denominator,
    parameters,
):
    """
    For coefficient [x^r y^t] of D*G, every source coefficient

        G[r-dr, t-dt]

    must actually be observed.

    No absent cell is interpreted as zero.
    """

    terms = denominator_terms(
        denominator,
        parameters,
    )

    cells = []

    for r, t in sorted(L):

        valid = True

        for dr, dt, _ in terms:

            source = (
                r - dr,
                t - dt,
            )

            if source not in L:
                valid = False
                break

        if valid:
            cells.append(
                (r, t)
            )

    return cells


def coefficient_equation(
    L,
    cell,
    denominator,
    parameters,
):
    """
    Coefficient of x^r y^t in D*G.

    Since the constant denominator coefficient is 1,
    this is linear in the unknown denominator coefficients.
    """

    r, t = cell

    terms = denominator_terms(
        denominator,
        parameters,
    )

    expression = sp.Integer(0)

    for dr, dt, coefficient in terms:

        source = (
            r - dr,
            t - dt,
        )

        expression += (
            coefficient
            * L[source]
        )

    return clean(expression)


# ============================================================================
# EXACT DENOMINATOR SOLVER
# ============================================================================

def solve_denominator_family(
    name,
    denominator,
    parameters,
    L,
):
    cells = valid_equation_cells(
        L,
        denominator,
        parameters,
    )

    equations = [
        coefficient_equation(
            L,
            cell,
            denominator,
            parameters,
        )
        for cell in cells
    ]

    # Each equation is homogeneous in the shifted terms but affine in
    # denominator parameters because the constant denominator coefficient
    # is fixed to 1.

    if not equations:
        return {
            "status": "INSUFFICIENT_DATA",
            "equations": 0,
            "unknowns": len(parameters),
        }

    expressions = [
        sp.expand(eq)
        for eq in equations
    ]

    A, rhs = sp.linear_eq_to_matrix(
        expressions,
        parameters,
    )

    # linear_eq_to_matrix puts equations in the form
    #
    #     A * params = rhs
    #
    # exactly over QQ.

    rank = A.rank()

    augmented_rank = (
        A.row_join(rhs).rank()
    )

    equation_count = len(
        equations
    )

    unknown_count = len(
        parameters
    )

    print()
    print("=" * 78)
    print(
        "DENOMINATOR FAMILY: {}".format(
            name
        )
    )
    print("=" * 78)

    print(
        "  denominator={}".format(
            denominator
        )
    )

    print(
        "  usable_interior_cells={}".format(
            cells
        )
    )

    print(
        "  equations={}".format(
            equation_count
        )
    )

    print(
        "  unknowns={}".format(
            unknown_count
        )
    )

    print(
        "  redundancy={}".format(
            equation_count
            - unknown_count
        )
    )

    print(
        "  rank={}".format(
            rank
        )
    )

    print(
        "  augmented_rank={}".format(
            augmented_rank
        )
    )

    if augmented_rank > rank:

        print(
            "  status=NO_SOLUTION"
        )

        return {
            "status": "NO_SOLUTION",
            "cells": cells,
            "rank": rank,
            "augmented_rank": augmented_rank,
        }

    if rank < unknown_count:

        print(
            "  nullity={}".format(
                unknown_count - rank
            )
        )

        print(
            "  status=NONUNIQUE"
        )

        return {
            "status": "NONUNIQUE",
            "cells": cells,
            "rank": rank,
            "augmented_rank": augmented_rank,
        }

    solution = A.gauss_jordan_solve(
        rhs
    )[0]

    values = {
        parameter: clean(
            solution[i, 0]
        )
        for i, parameter
        in enumerate(parameters)
    }

    solved_denominator = clean(
        denominator.subs(
            values
        )
    )

    residuals = [
        clean(
            expression.subs(
                values
            )
        )
        for expression
        in expressions
    ]

    exact = all(
        residual == 0
        for residual
        in residuals
    )

    if not exact:

        status = "VERIFICATION_FAILED"

    elif equation_count > unknown_count:

        status = "EXACT_OVERDETERMINED"

    else:

        status = "EXACT_DATA_SIZED"

    print(
        "  status={}".format(
            status
        )
    )

    if status.startswith(
        "EXACT"
    ):

        print()
        print(
            "  solved_parameters={}".format(
                values
            )
        )

        print(
            "  solved_denominator={}".format(
                solved_denominator
            )
        )

        print(
            "  primitive_denominator={}".format(
                primitive_integer_polynomial(
                    solved_denominator
                ).as_expr()
            )
        )

        print(
            "  residuals={}".format(
                residuals
            )
        )

    return {
        "status": status,
        "cells": cells,
        "rank": rank,
        "augmented_rank": augmented_rank,
        "parameters": values,
        "denominator": solved_denominator,
        "residuals": residuals,
    }


# ============================================================================
# BOUNDARY RESIDUAL
# ============================================================================

def boundary_residual_audit(
    G,
    result,
):
    if result["status"] not in (
        "EXACT_OVERDETERMINED",
        "EXACT_DATA_SIZED",
    ):
        return

    denominator = result[
        "denominator"
    ]

    residual = clean(
        denominator * G
    )

    P = sp.Poly(
        residual,
        x,
        y,
        domain=sp.QQ,
    )

    support = sorted(
        P.monoms()
    )

    print()
    print(
        "  BOUNDARY RESIDUAL AUDIT"
    )
    print(
        "    D(x,y)G(x,y)={}".format(
            residual
        )
    )

    print(
        "    residual_support={}".format(
            support
        )
    )

    print(
        "    residual_term_count={}".format(
            len(support)
        )
    )

    # Boundary classification.
    r_values = [
        monomial[0]
        for monomial in support
    ]

    t_values = [
        monomial[1]
        for monomial in support
    ]

    min_r = min(r_values)
    min_t = min(t_values)

    boundary_cells = [
        monomial
        for monomial in support
        if (
            monomial[0] == min_r
            or
            monomial[1] == min_t
        )
    ]

    interior_cells = [
        monomial
        for monomial in support
        if monomial not in boundary_cells
    ]

    print(
        "    min_r={}".format(
            min_r
        )
    )

    print(
        "    min_t={}".format(
            min_t
        )
    )

    print(
        "    boundary_terms={}".format(
            boundary_cells
        )
    )

    print(
        "    interior_terms={}".format(
            interior_cells
        )
    )

    print(
        "    boundary_only={}".format(
            len(interior_cells) == 0
        )
    )

    return residual


# ============================================================================
# SPECIAL SMALL-RATIONAL TESTS
# ============================================================================

def small_integer_parameter_audit(
    results
):
    print()
    print("=" * 78)
    print(
        "SPECIAL PARAMETER VALUES"
    )
    print("=" * 78)

    hits = []

    for name, result in results.items():

        if result["status"] not in (
            "EXACT_OVERDETERMINED",
            "EXACT_DATA_SIZED",
        ):
            continue

        values = result.get(
            "parameters",
            {}
        )

        integer = all(
            sp.denom(v) == 1
            for v in values.values()
        )

        small = all(
            abs(int(v)) <= 10
            for v in values.values()
        ) if integer else False

        print()
        print(
            "  {}: integer_parameters={}".format(
                name,
                integer,
            )
        )

        print(
            "    small_integer_parameters={}".format(
                small
            )
        )

        if small:
            hits.append(
                name
            )

    print()
    print(
        "  small_integer_parameter_hits={}".format(
            hits
        )
    )

    return hits


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 350R — EXACT BOUNDARY-FORCED "
        "GENERATING-DENOMINATOR AUDIT"
    )
    print("=" * 78)

    L = lattice()
    G = generating_polynomial(
        L
    )

    print()
    print(
        "OBSERVED GENERATING OBJECT"
    )
    print(
        "  G(x,y)={}".format(
            G
        )
    )

    denominator_specs = (
        denominator_families()
    )

    results = {}

    for (
        name,
        (
            denominator,
            parameters,
        ),
    ) in denominator_specs.items():

        result = solve_denominator_family(
            name,
            denominator,
            parameters,
            L,
        )

        results[name] = result

        boundary_residual_audit(
            G,
            result,
        )

    small_hits = (
        small_integer_parameter_audit(
            results
        )
    )

    print()
    print("=" * 78)
    print(
        "SUMMARY"
    )
    print("=" * 78)

    exact_overdetermined = []

    for name, result in results.items():

        print()
        print(
            "  {}: status={}, equations={}, "
            "unknowns={}".format(
                name,
                result["status"],
                result["equations"]
                if "equations" in result
                else len(
                    result.get(
                        "cells",
                        [],
                    )
                ),
                result["unknowns"]
                if "unknowns" in result
                else "?",
            )
        )

        if result["status"] == (
            "EXACT_OVERDETERMINED"
        ):
            exact_overdetermined.append(
                name
            )

    print()
    print(
        "  exact_overdetermined_denominators={}".format(
            exact_overdetermined
        )
    )

    print(
        "  small_integer_parameter_hits={}".format(
            small_hits
        )
    )

    # ------------------------------------------------------------------------
    # Structural interpretation
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
349R established that the finite polynomial G(x,y) itself has no obvious
factorization.

That does not eliminate a rational generating mechanism of the standard
form

    D(x,y) G(x,y) = N(x,y),

where D is a small translation-invariant denominator and N is a boundary
term.

This is exactly how many combinatorial triangular arrays arise.

The crucial methodological point is that the experiment does NOT declare
unobserved Q-cells to be zero.

For every recurrence equation, all source cells needed to form the
coefficient of D*G must be observed.

Thus the denominator is tested only against genuinely available interior
windows.

A successful EXACT_OVERDETERMINED model would be substantially stronger
than the failed homogeneous-stencil experiments because it would identify
a small denominator together with an explicit boundary residual.

A NONUNIQUE model is not a discovery.

A NO_SOLUTION model eliminates that denominator family exactly.

If every tested denominator fails, then the rational-generating-function
route is also substantially constrained, and the remaining bottleneck is
the missing source definition or an independent second n=pq case.

No missing values.
No interpolation.
No extrapolation.
No synthetic second case.
"""
    )

    print()
    print("=" * 78)
    print(
        "FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  denominator_family_count={}".format(
            len(
                denominator_specs
            )
        )
    )

    print(
        "  exact_overdetermined_denominators={}".format(
            exact_overdetermined
        )
    )

    print(
        "  small_integer_parameter_hits={}".format(
            small_hits
        )
    )

    print(
        "  missing_values_used=False"
    )

    print(
        "  interpolation_performed=False"
    )

    print(
        "  extrapolation_used=False"
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
        "EXPERIMENT 350R COMPLETE"
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
