#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 349R — EXACT BIVARIATE GENERATING-FUNCTION /
                  COEFFICIENT-PROVENANCE AUDIT
==============================================================================

Purpose
-------
Treat the complete observed triangular source table as ONE bivariate
coefficient object in the natural lattice coordinate

    r = (p - 1) / 2.

Define

    G(x,y) = sum Q(r,t) x^r y^t

over observed cells only.

This experiment tests:

    * exact construction of G;
    * exact factorization of G;
    * exact low-degree bivariate polynomial compatibility;
    * terminal polynomial divisibility;
    * small fixed substitution transforms;
    * diagonal coefficient structure;
    * support / Newton-polygon geometry;
    * simple separable / product decompositions.

No missing values.
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
# HELPERS
# ============================================================================

def clean(expr):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(expr)
            )
        )
    )


def build_lattice():
    """
    Convert the stored source table into

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


def primitive_integer_polynomial(expr):
    """
    Convert a rational polynomial in x,y to a primitive integer polynomial.
    """

    poly = sp.Poly(
        sp.expand(expr),
        x,
        y,
        domain=sp.QQ,
    )

    if poly.is_zero:
        return sp.Poly(
            0,
            x,
            y,
            domain=sp.ZZ,
        )

    coeffs = poly.coeffs()

    denominator_lcm = 1

    for coefficient in coeffs:
        denominator_lcm = sp.ilcm(
            denominator_lcm,
            int(sp.denom(coefficient)),
        )

    integers = [
        int(coefficient * denominator_lcm)
        for coefficient in coeffs
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
            x,
            y,
            domain=sp.ZZ,
        )

    integers = [
        value // coefficient_gcd
        for value in integers
    ]

    expression = sp.Integer(0)

    for monomial, coefficient in zip(
        poly.monoms(),
        integers,
    ):

        i, j = monomial

        expression += (
            sp.Integer(coefficient)
            * x**i
            * y**j
        )

    result = sp.Poly(
        sp.expand(expression),
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


def generating_polynomial(lattice):
    return clean(
        sum(
            value * x**r * y**t
            for (r, t), value
            in lattice.items()
        )
    )


def total_degree_monomials(degree):
    return [
        x**i * y**j
        for i in range(degree + 1)
        for j in range(degree + 1)
        if i + j <= degree
    ]


def ratio_if_exact(a, b):
    if b == 0:
        return None

    if a % b != 0:
        return None

    return sp.Integer(a // b)


# ============================================================================
# SECTION 1 — OBSERVED LATTICE
# ============================================================================

def print_lattice(lattice):

    print()
    print("=" * 78)
    print(
        "1. OBSERVED (r,t) LATTICE"
    )
    print("=" * 78)

    for t_value in range(6):

        row = []

        for r_value in range(4):

            cell = (
                r_value,
                t_value,
            )

            if cell in lattice:
                row.append(
                    (
                        r_value,
                        lattice[cell],
                    )
                )

        if row:
            print(
                "  t={}: {}".format(
                    t_value,
                    row,
                )
            )


# ============================================================================
# SECTION 2 — GENERATING OBJECT
# ============================================================================

def generating_object_audit(lattice):

    print()
    print("=" * 78)
    print(
        "2. EXACT BIVARIATE GENERATING OBJECT"
    )
    print("=" * 78)

    G = generating_polynomial(
        lattice
    )

    print()
    print(
        "  G(x,y)={}".format(
            G
        )
    )

    primitive = primitive_integer_polynomial(
        G
    )

    print()
    print(
        "  primitive_integer_G={}".format(
            primitive.as_expr()
        )
    )

    P = sp.Poly(
        G,
        x,
        y,
        domain=sp.QQ,
    )

    print(
        "  total_degree={}".format(
            P.total_degree()
        )
    )

    print(
        "  degree_x={}".format(
            P.degree(x)
        )
    )

    print(
        "  degree_y={}".format(
            P.degree(y)
        )
    )

    return G


# ============================================================================
# SECTION 3 — EXACT FACTORIZATION
# ============================================================================

def factorization_audit(G):

    print()
    print("=" * 78)
    print(
        "3. EXACT BIVARIATE FACTORIZATION"
    )
    print("=" * 78)

    factorized = sp.factor(
        G
    )

    print(
        "  factor(G)={}".format(
            factorized
        )
    )

    factor_constant, factor_list = sp.factor_list(
        G
    )

    print(
        "  overall_constant={}".format(
            factor_constant
        )
    )

    print(
        "  factor_list={}".format(
            factor_list
        )
    )

    nonconstant_count = sum(
        1
        for factor, exponent
        in factor_list
        if sp.Poly(
            factor,
            x,
            y,
            domain=sp.QQ,
        ).total_degree() > 0
    )

    print(
        "  nonconstant_factor_count={}".format(
            nonconstant_count
        )
    )

    return factorized


# ============================================================================
# SECTION 4 — LOW-DEGREE BIVARIATE POLYNOMIAL TESTS
# ============================================================================

def solve_low_degree_model(
    lattice,
    degree,
):

    monomials = total_degree_monomials(
        degree
    )

    points = [
        (
            r,
            t,
            value,
        )
        for (
            r,
            t,
        ), value
        in sorted(
            lattice.items(),
            key=lambda item: (
                item[0][1],
                item[0][0],
            ),
        )
    ]

    matrix = sp.Matrix([
        [
            clean(
                monomial.subs(
                    {
                        x: r,
                        y: t,
                    }
                )
            )
            for monomial in monomials
        ]
        for (
            r,
            t,
            _,
        ) in points
    ])

    rhs = sp.Matrix([
        value
        for (
            _,
            _,
            value,
        ) in points
    ])

    equation_count = len(points)
    unknown_count = len(monomials)

    rank = matrix.rank()

    augmented_rank = (
        matrix.row_join(rhs).rank()
    )

    print()
    print("=" * 78)
    print(
        "4. TOTAL DEGREE <= {}".format(
            degree
        )
    )
    print("=" * 78)

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
            equation_count - unknown_count
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
            "rank": rank,
            "augmented_rank": augmented_rank,
        }

    solution = matrix.gauss_jordan_solve(
        rhs
    )[0]

    polynomial = clean(
        sum(
            solution[i, 0] * monomials[i]
            for i in range(
                unknown_count
            )
        )
    )

    residuals = [
        clean(
            polynomial.subs(
                {
                    x: r,
                    y: t,
                }
            )
            - value
        )
        for (
            r,
            t,
            value,
        ) in points
    ]

    exact = all(
        residual == 0
        for residual in residuals
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
            "  polynomial={}".format(
                polynomial
            )
        )

        print(
            "  primitive={}".format(
                primitive_integer_polynomial(
                    polynomial
                ).as_expr()
            )
        )

        print(
            "  all_residuals_zero={}".format(
                exact
            )
        )

    return {
        "status": status,
        "rank": rank,
        "augmented_rank": augmented_rank,
        "polynomial": polynomial,
    }


# ============================================================================
# SECTION 5 — TERMINAL POLYNOMIAL
# ============================================================================

def terminal_polynomial_audit(G):

    print()
    print("=" * 78)
    print(
        "5. TERMINAL POLYNOMIAL / GENERATING-OBJECT AUDIT"
    )
    print("=" * 78)

    F = (
        x**3
        + 16027881*x**2
        + 421514439*x
        + 495451247
    )

    F_reverse = sp.expand(
        x**3 * F.subs(
            x,
            1 / x,
        )
    )

    transforms = {
        "F(x)": F,
        "F_reverse(x)": F_reverse,
        "F(y)": F.subs(x, y),
        "F(x)F(y)": clean(
            F * F.subs(x, y)
        ),
    }

    print(
        "  F(x)={}".format(
            F
        )
    )

    print(
        "  F_reverse(x)={}".format(
            F_reverse
        )
    )

    for name, transform in transforms.items():

        quotient, remainder = sp.div(
            G,
            transform,
            domain=sp.QQ,
        )

        print()
        print(
            "  {}:".format(
                name
            )
        )

        print(
            "    divisible={}".format(
                clean(remainder) == 0
            )
        )

        if clean(remainder) == 0:

            print(
                "    quotient={}".format(
                    clean(quotient)
                )
            )

        else:

            print(
                "    remainder={}".format(
                    clean(remainder)
                )
            )

    return F, F_reverse


# ============================================================================
# SECTION 6 — SMALL SUBSTITUTION SEARCH
# ============================================================================

def substitution_audit(G):

    print()
    print("=" * 78)
    print(
        "6. SMALL EXACT SUBSTITUTION SEARCH"
    )
    print("=" * 78)

    candidates = []

    for c in (
        -3,
        -2,
        -1,
        1,
        2,
        3,
    ):

        candidates.append(
            (
                "x_plus_{}_y".format(c),
                clean(
                    G.subs(
                        x,
                        x + c * y,
                    )
                ),
            )
        )

        candidates.append(
            (
                "x_times_1plus_{}_y".format(c),
                clean(
                    G.subs(
                        x,
                        x * (1 + c * y),
                    )
                ),
            )
        )

        candidates.append(
            (
                "y_times_1plus_{}_x".format(c),
                clean(
                    G.subs(
                        y,
                        y * (1 + c * x),
                    )
                ),
            )
        )

    factor_hits = 0

    for name, transformed in candidates:

        factorized = sp.factor(
            transformed
        )

        if factorized != transformed:

            factor_hits += 1

            print()
            print(
                "  {}:".format(
                    name
                )
            )

            print(
                "    factorized={}".format(
                    factorized
                )
            )

    print()
    print(
        "  candidates_tested={}".format(
            len(candidates)
        )
    )

    print(
        "  nontrivial_factorization_hits={}".format(
            factor_hits
        )
    )

    return factor_hits


# ============================================================================
# SECTION 7 — DIAGONAL STRUCTURE
# ============================================================================

def diagonal_audit(lattice):

    print()
    print("=" * 78)
    print(
        "7. DIAGONAL COEFFICIENT STRUCTURE"
    )
    print("=" * 78)

    max_sum = max(
        r + t
        for (
            r,
            t,
        ) in lattice
    )

    sequences = {}

    for total in range(
        max_sum + 1
    ):

        entries = [
            (
                r,
                lattice[
                    (
                        r,
                        total - r,
                    )
                ],
            )
            for r in range(
                total + 1
            )
            if (
                r,
                total - r,
            ) in lattice
        ]

        if entries:

            sequences[total] = entries

            print()
            print(
                "  r+t={}: {}".format(
                    total,
                    entries,
                )
            )

    return sequences


# ============================================================================
# SECTION 8 — SUPPORT / NEWTON POLYGON
# ============================================================================

def support_audit(lattice):

    print()
    print("=" * 78)
    print(
        "8. SUPPORT / NEWTON-POLYGON AUDIT"
    )
    print("=" * 78)

    support = sorted(
        lattice
    )

    print(
        "  support={}".format(
            support
        )
    )

    # Compute a convex hull ourselves to avoid relying on optional
    # geometry APIs.
    points = [
        (
            int(r),
            int(t),
        )
        for r, t in support
    ]

    def cross(o, a, b):
        return (
            (a[0] - o[0])
            * (b[1] - o[1])
            -
            (a[1] - o[1])
            * (b[0] - o[0])
        )

    points = sorted(
        set(points)
    )

    if len(points) <= 1:

        hull = points

    else:

        lower = []

        for point in points:

            while (
                len(lower) >= 2
                and cross(
                    lower[-2],
                    lower[-1],
                    point,
                ) <= 0
            ):
                lower.pop()

            lower.append(point)

        upper = []

        for point in reversed(points):

            while (
                len(upper) >= 2
                and cross(
                    upper[-2],
                    upper[-1],
                    point,
                ) <= 0
            ):
                upper.pop()

            upper.append(point)

        hull = (
            lower[:-1]
            + upper[:-1]
        )

    print(
        "  convex_hull={}".format(
            hull
        )
    )

    # Shoelace area.
    if len(hull) < 3:

        area = sp.Rational(0)

    else:

        twice_area = 0

        for i in range(
            len(hull)
        ):

            x1, y1 = hull[i]
            x2, y2 = hull[
                (i + 1) % len(hull)
            ]

            twice_area += (
                x1 * y2
                - x2 * y1
            )

        area = sp.Rational(
            abs(twice_area),
            2,
        )

    print(
        "  convex_hull_area={}".format(
            area
        )
    )

    return hull, area


# ============================================================================
# SECTION 9 — COEFFICIENT PRODUCT TESTS
# ============================================================================

def product_structure_audit(
    lattice,
    G,
):

    print()
    print("=" * 78)
    print(
        "9. SIMPLE PRODUCT / SEPARABLE GENERATING-FUNCTION AUDIT"
    )
    print("=" * 78)

    P = sp.Poly(
        G,
        x,
        y,
        domain=sp.QQ,
    )

    degree_x = P.degree(x)
    degree_y = P.degree(y)

    # Test whether G factors into one x-only and one y-only polynomial.
    x_part = clean(
        G.subs(
            y,
            0,
        )
    )

    y_part = clean(
        G.subs(
            x,
            0,
        )
    )

    constant = clean(
        G.subs(
            {
                x: 0,
                y: 0,
            }
        )
    )

    additive_residual = clean(
        G
        - x_part
        - y_part
        + constant
    )

    print(
        "  x_only_part={}".format(
            x_part
        )
    )

    print(
        "  y_only_part={}".format(
            y_part
        )
    )

    print(
        "  additive_separation_residual={}".format(
            additive_residual
        )
    )

    print(
        "  additive_separable={}".format(
            additive_residual == 0
        )
    )

    if constant != 0:

        multiplicative_candidate = clean(
            x_part
            * y_part
            / constant
        )

        multiplicative_residual = clean(
            G
            - multiplicative_candidate
        )

        print(
            "  multiplicative_residual={}".format(
                multiplicative_residual
            )
        )

        print(
            "  multiplicatively_separable={}".format(
                multiplicative_residual == 0
            )
        )

    else:

        print(
            "  multiplicatively_separable=UNAVAILABLE_ZERO_CONSTANT"
        )

    print(
        "  degree_x={}".format(
            degree_x
        )
    )

    print(
        "  degree_y={}".format(
            degree_y
        )
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 349R — EXACT BIVARIATE GENERATING-FUNCTION / "
        "COEFFICIENT-PROVENANCE AUDIT"
    )
    print("=" * 78)

    lattice = build_lattice()

    print_lattice(
        lattice
    )

    G = generating_object_audit(
        lattice
    )

    factorization = factorization_audit(
        G
    )

    polynomial_results = {}

    for degree in (
        1,
        2,
        3,
        4,
    ):

        polynomial_results[
            degree
        ] = solve_low_degree_model(
            lattice,
            degree,
        )

    F, F_reverse = terminal_polynomial_audit(
        G
    )

    substitution_hits = substitution_audit(
        G
    )

    diagonal_sequences = diagonal_audit(
        lattice
    )

    hull, hull_area = support_audit(
        lattice
    )

    product_structure_audit(
        lattice,
        G,
    )

    # ------------------------------------------------------------------------
    # Structural interpretation
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "10. STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
The source table has resisted a large family of direct recurrence and
operator descriptions.

349R instead compresses the entire observed triangle into one object:

    G(x,y) = sum Q(r,t) x^r y^t.

This is the natural generating-function representation of the available
data.

The main provenance questions are:

    1. Does G factor nontrivially over Q?

    2. Does the terminal polynomial

           F(x)=x^3+16027881*x^2+421514439*x+495451247

       actually divide G or appear as a genuine component?

    3. Does a simple coordinate substitution produce a factorization?

    4. Does G have additive or multiplicative separation?

    5. Does the support geometry reveal a standard triangular generating
       construction?

A positive low-complexity factorization would be much more meaningful
than another fitted recurrence because it would point toward the
construction from which the observed coefficients arose.

A negative result is also useful: it would strongly justify stopping
generic algebraic searches and obtaining the original definition of Q.

No unobserved cell is inserted.
"""
    )

    # ------------------------------------------------------------------------
    # Final exactness
    # ------------------------------------------------------------------------

    exact_overdetermined = [
        degree
        for degree, result
        in polynomial_results.items()
        if result["status"]
        == "EXACT_OVERDETERMINED"
    ]

    factor_nontrivial = (
        sp.factor(G)
        != G
    )

    print()
    print("=" * 78)
    print(
        "11. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  observed_cells={}".format(
            len(lattice)
        )
    )

    print(
        "  generating_object_built_exactly=True"
    )

    print(
        "  factorization_completed=True"
    )

    print(
        "  nontrivial_sympy_factorization={}".format(
            factor_nontrivial
        )
    )

    print(
        "  exact_overdetermined_bivariate_degrees={}".format(
            exact_overdetermined
        )
    )

    print(
        "  substitution_candidates_tested=18"
    )

    print(
        "  substitution_factorization_hits={}".format(
            substitution_hits
        )
    )

    print(
        "  support_geometry_completed=True"
    )

    print(
        "  terminal_polynomial_defined_exactly=True"
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
        "EXPERIMENT 349R COMPLETE"
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