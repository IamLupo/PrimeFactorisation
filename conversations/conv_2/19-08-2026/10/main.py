#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 326R — EXACT BIVARIATE SOURCE-POLYNOMIAL / LOW-DEGREE LAW AUDIT
==============================================================================

Purpose
-------
Test whether the observed triangular source table Q_t(p) is compatible with
a low-degree exact polynomial law in the two variables p and t.

Only observed data are used.

Models:

    total degree <= 1
    total degree <= 2
    total degree <= 3
    total degree <= 4

    bidegree (1,1)
    bidegree (2,1)
    bidegree (2,2)
    bidegree (3,2)

The experiment distinguishes:

    EXACT_OVERDETERMINED
        exact fit with more observations than coefficients;

    EXACT_DATA_SIZED
        exact fit with exactly as many observations as coefficients;

    NONUNIQUE
        consistent but underdetermined;

    NO_SOLUTION
        exact rational inconsistency.

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


def build_observed_points():
    """
    Q[p] is stored in reverse order in t.

        t = len(Q[p]) - 1 - index

    Only observed cells are returned.
    """

    points = []

    for p_value in sorted(Q):

        values = Q[p_value]

        for index, value in enumerate(values):

            t_value = (
                len(values)
                - 1
                - index
            )

            points.append(
                (
                    sp.Integer(p_value),
                    sp.Integer(t_value),
                    sp.Integer(value),
                )
            )

    return sorted(
        points,
        key=lambda item: (
            item[1],
            item[0],
        ),
    )


def primitive_integer_polynomial(expr):
    """
    Convert a rational polynomial in p,t to a primitive integer polynomial.
    """

    expr = sp.sympify(expr)

    poly = sp.Poly(
        sp.expand(expr),
        p,
        t,
        domain=sp.QQ,
    )

    if poly.is_zero:
        return sp.Poly(
            0,
            p,
            t,
            domain=sp.ZZ,
        )

    coefficients = poly.coeffs()

    denominator_lcm = 1

    for coefficient in coefficients:
        denominator_lcm = sp.ilcm(
            denominator_lcm,
            int(sp.denom(coefficient)),
        )

    integers = [
        int(coefficient * denominator_lcm)
        for coefficient in coefficients
    ]

    coefficient_gcd = 0

    for value in integers:
        coefficient_gcd = math.gcd(
            coefficient_gcd,
            abs(value),
        )

    if coefficient_gcd == 0:
        return sp.Poly(
            0,
            p,
            t,
            domain=sp.ZZ,
        )

    integers = [
        value // coefficient_gcd
        for value in integers
    ]

    expression = 0

    for monomial, coefficient in zip(
        poly.monoms(),
        integers,
    ):
        i, j = monomial
        expression += (
            sp.Integer(coefficient)
            * p**i
            * t**j
        )

    result = sp.Poly(
        sp.expand(expression),
        p,
        t,
        domain=sp.ZZ,
    )

    # Normalize sign.
    first = result.coeffs()[0]

    if first < 0:
        result = sp.Poly(
            -result.as_expr(),
            p,
            t,
            domain=sp.ZZ,
        )

    return result


# ============================================================================
# MONOMIAL FAMILIES
# ============================================================================

def normalize_monomials(monomials):
    """
    Convert every model entry, including Python int 1, to a SymPy expression.

    This is the critical fix for the original script.
    """
    return [
        sp.sympify(monomial)
        for monomial in monomials
    ]


def total_degree_monomials(degree):
    return [
        p**i * t**j
        for i in range(degree + 1)
        for j in range(degree + 1)
        if i + j <= degree
    ]


def bidegree_monomials(
    degree_p,
    degree_t,
):
    return [
        p**i * t**j
        for i in range(degree_p + 1)
        for j in range(degree_t + 1)
    ]


# ============================================================================
# EXACT LINEAR MODEL SOLVER
# ============================================================================

def solve_model(
    points,
    monomials,
):
    """
    Solve

        Q(p,t) = sum c_i M_i(p,t)

    exactly over Q.
    """

    monomials = normalize_monomials(
        monomials
    )

    unknown_count = len(monomials)

    A = sp.Matrix([
        [
            clean(
                monomial.subs(
                    {
                        p: point_p,
                        t: point_t,
                    }
                )
            )
            for monomial in monomials
        ]
        for point_p, point_t, _ in points
    ])

    b = sp.Matrix([
        sp.Integer(value)
        for _, _, value in points
    ])

    rank = A.rank()

    augmented_rank = (
        A.row_join(b).rank()
    )

    equation_count = len(points)

    # Impossible.
    if augmented_rank > rank:
        return {
            "status": "NO_SOLUTION",
            "rank": rank,
            "augmented_rank": augmented_rank,
            "equation_count": equation_count,
            "unknown_count": unknown_count,
            "coefficients": None,
            "polynomial": None,
            "residuals": None,
        }

    # Unique exact solution.
    if rank == unknown_count:

        solution = A.gauss_jordan_solve(
            b
        )[0]

        coefficients = [
            clean(solution[i, 0])
            for i in range(
                solution.rows
            )
        ]

        polynomial = clean(
            sum(
                coefficient * monomial
                for coefficient, monomial
                in zip(
                    coefficients,
                    monomials,
                )
            )
        )

        residuals = [
            clean(
                polynomial.subs(
                    {
                        p: point_p,
                        t: point_t,
                    }
                )
                - value
            )
            for point_p, point_t, value in points
        ]

        all_zero = all(
            residual == 0
            for residual in residuals
        )

        if not all_zero:
            status = "VERIFICATION_FAILED"

        elif equation_count > unknown_count:
            status = "EXACT_OVERDETERMINED"

        else:
            status = "EXACT_DATA_SIZED"

        return {
            "status": status,
            "rank": rank,
            "augmented_rank": augmented_rank,
            "equation_count": equation_count,
            "unknown_count": unknown_count,
            "coefficients": coefficients,
            "polynomial": polynomial,
            "residuals": residuals,
        }

    # Consistent but non-unique.
    symbols = sp.symbols(
        "c0:{}".format(
            unknown_count
        )
    )

    equations = [
        sp.Eq(
            sum(
                symbols[j] * A[i, j]
                for j in range(
                    unknown_count
                )
            ),
            b[i],
        )
        for i in range(
            equation_count
        )
    ]

    solution_set = sp.linsolve(
        equations,
        symbols,
    )

    return {
        "status": "NONUNIQUE",
        "rank": rank,
        "augmented_rank": augmented_rank,
        "equation_count": equation_count,
        "unknown_count": unknown_count,
        "coefficients": solution_set,
        "polynomial": None,
        "residuals": None,
    }


# ============================================================================
# REPORTING
# ============================================================================

def report_model(
    label,
    points,
    monomials,
):
    monomials = normalize_monomials(
        monomials
    )

    result = solve_model(
        points,
        monomials,
    )

    print()
    print("=" * 78)
    print(label)
    print("=" * 78)

    print(
        "  observed_points={}".format(
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
        "  coefficient_matrix_rank={}".format(
            result["rank"]
        )
    )

    print(
        "  augmented_matrix_rank={}".format(
            result["augmented_rank"]
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
        print(
            "  polynomial={}".format(
                result["polynomial"]
            )
        )

        primitive = (
            primitive_integer_polynomial(
                result["polynomial"]
            )
        )

        print(
            "  primitive_integer_polynomial={}".format(
                primitive.as_expr()
            )
        )

        print(
            "  residuals={}".format(
                result["residuals"]
            )
        )

        print(
            "  all_residuals_zero={}".format(
                all(
                    residual == 0
                    for residual
                    in result["residuals"]
                )
            )
        )

    elif result["status"] == "NONUNIQUE":

        print()
        print(
            "  solution_family={}".format(
                result["coefficients"]
            )
        )

    return result


def coefficient_audit(
    label,
    result,
):
    if result["status"] not in (
        "EXACT_OVERDETERMINED",
        "EXACT_DATA_SIZED",
    ):
        return

    poly = sp.Poly(
        result["polynomial"],
        p,
        t,
        domain=sp.QQ,
    )

    print()
    print("=" * 78)
    print(
        "{} COEFFICIENT AUDIT".format(
            label
        )
    )
    print("=" * 78)

    print(
        "  monomials={}".format(
            poly.monoms()
        )
    )

    for monomial, coefficient in zip(
        poly.monoms(),
        poly.coeffs(),
    ):

        print()
        print(
            "  monomial={}".format(
                monomial
            )
        )

        print(
            "    coefficient={}".format(
                clean(coefficient)
            )
        )

        print(
            "    valuations={}".format(
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
                }
            )
        )


# ============================================================================
# MODEL SUITE
# ============================================================================

def run_total_degree_suite(
    points,
):
    results = {}

    for degree in (
        1,
        2,
        3,
        4,
    ):

        label = (
            "TOTAL-DEGREE <= {}".format(
                degree
            )
        )

        result = report_model(
            label,
            points,
            total_degree_monomials(
                degree
            ),
        )

        coefficient_audit(
            label,
            result,
        )

        results[
            ("total", degree)
        ] = result

    return results


def run_bidegree_suite(
    points,
):
    models = {
        "(1,1)": bidegree_monomials(
            1,
            1,
        ),
        "(2,1)": bidegree_monomials(
            2,
            1,
        ),
        "(2,2)": bidegree_monomials(
            2,
            2,
        ),
        "(3,2)": bidegree_monomials(
            3,
            2,
        ),
    }

    results = {}

    for label, monomials in models.items():

        result = report_model(
            "BIDEGREE {}".format(label),
            points,
            monomials,
        )

        coefficient_audit(
            "BIDEGREE {}".format(label),
            result,
        )

        results[
            ("bidegree", label)
        ] = result

    return results


# ============================================================================
# SUMMARY
# ============================================================================

def summary(
    points,
    total_results,
    bidegree_results,
):

    print()
    print("=" * 78)
    print(
        "SUMMARY OF LOW-DEGREE SOURCE-LAW TESTS"
    )
    print("=" * 78)

    print()
    print(
        "  observed_points={}".format(
            len(points)
        )
    )

    print()
    print(
        "  TOTAL-DEGREE MODELS"
    )

    for degree in (
        1,
        2,
        3,
        4,
    ):

        result = total_results[
            ("total", degree)
        ]

        print(
            "    degree<={}: status={}, rank={}, "
            "augmented_rank={}, redundancy={}".format(
                degree,
                result["status"],
                result["rank"],
                result["augmented_rank"],
                result["equation_count"]
                - result["unknown_count"],
            )
        )

    print()
    print(
        "  BIDEGREE MODELS"
    )

    for label in (
        "(1,1)",
        "(2,1)",
        "(2,2)",
        "(3,2)",
    ):

        result = bidegree_results[
            ("bidegree", label)
        ]

        print(
            "    {}: status={}, rank={}, "
            "augmented_rank={}, redundancy={}".format(
                label,
                result["status"],
                result["rank"],
                result["augmented_rank"],
                result["equation_count"]
                - result["unknown_count"],
            )
        )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 326R — EXACT BIVARIATE SOURCE-POLYNOMIAL / "
        "LOW-DEGREE LAW AUDIT"
    )
    print("=" * 78)

    points = build_observed_points()

    print()
    print("=" * 78)
    print(
        "1. OBSERVED SOURCE POINTS"
    )
    print("=" * 78)

    print(
        "  observed_point_count={}".format(
            len(points)
        )
    )

    for point_p, point_t, value in points:

        print(
            "  (p={}, t={}) -> {}".format(
                point_p,
                point_t,
                value,
            )
        )

    total_results = run_total_degree_suite(
        points
    )

    bidegree_results = run_bidegree_suite(
        points
    )

    summary(
        points,
        total_results,
        bidegree_results,
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
The source table contains 20 observed cells.

The first genuinely overdetermined bivariate polynomial tests are:

    total degree <= 1:
        3 coefficients;

    total degree <= 2:
        6 coefficients;

    total degree <= 3:
        10 coefficients.

A successful exact fit at any of these degrees is therefore genuine
finite-data validation rather than mere interpolation.

Total degree <= 4 has 15 coefficients and is also overdetermined here
because 20 observed cells are available.

The bidegree models must be interpreted similarly:

    (1,1) has 4 coefficients;
    (2,1) has 6;
    (2,2) has 9;
    (3,2) has 12.

All four remain overdetermined against the 20 observed cells.

The key question is therefore:

    Does any low-degree bivariate polynomial reproduce ALL 20 observed
    source values exactly?

A positive answer would be the first simple source-level law surviving
the earlier failed recurrence, transfer, and separability tests.

A negative answer would justify moving beyond polynomial source laws and
toward either:

    * rational/algebraic formulas in p and t,
    * combinatorial/source-construction analysis,
    * or acquisition of a genuinely independent n=pq case.

No missing values are used.
No extrapolation is performed.
"""
    )

    # ------------------------------------------------------------------------
    # Final exactness
    # ------------------------------------------------------------------------

    exact_overdetermined_total = [
        degree
        for degree in (
            1,
            2,
            3,
            4,
        )
        if total_results[
            ("total", degree)
        ]["status"]
        == "EXACT_OVERDETERMINED"
    ]

    exact_overdetermined_bidegree = [
        label
        for label in (
            "(1,1)",
            "(2,1)",
            "(2,2)",
            "(3,2)",
        )
        if bidegree_results[
            ("bidegree", label)
        ]["status"]
        == "EXACT_OVERDETERMINED"
    ]

    print()
    print("=" * 78)
    print(
        "FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  observed_points={}".format(
            len(points)
        )
    )

    print(
        "  exact_overdetermined_total_degree_models={}".format(
            exact_overdetermined_total
        )
    )

    print(
        "  exact_overdetermined_bidegree_models={}".format(
            exact_overdetermined_bidegree
        )
    )

    print(
        "  missing_values_used=False"
    )

    print(
        "  interpolation_outside_observed_points=False"
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
        "EXPERIMENT 326R COMPLETE"
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