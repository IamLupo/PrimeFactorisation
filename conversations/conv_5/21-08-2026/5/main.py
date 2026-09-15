#!/usr/bin/env python3

# ==============================================================================
# EXPERIMENT 578
# ==============================================================================
# DISCOVER MOD-8 FUNCTIONS AND TRANSFORM A -> B
#
# IDEA
# ----
#
# We already have the mod-4 functions:
#
#     A1(x,y)
#     A2(x,y)
#
# We DO NOT try to derive B from residue bits directly.
#
# Instead:
#
#     1. collect the complete mod-8 datasets;
#
#     2. split them into:
#
#            B1 : n = 1 mod 8
#            B2 : n = 3 mod 8
#            B3 : n = 5 mod 8
#            B4 : n = 7 mod 8
#
#     3. recover each B_i(x,y) independently;
#
#     4. compare each B_i with its parent A_j;
#
#     5. search for STATIC transformations:
#
#            A1 -> B1
#            A1 -> B2
#            A2 -> B3
#            A2 -> B4
#
#
# The important point is:
#
#     B_i is discovered FIRST.
#
# We do not impose a presumed form on B_i.
#
#
# FUNCTION DISCOVERY
# ------------------
#
# We search exact integer-valued functions in a controlled basis:
#
#     1
#     x
#     y
#     x^2
#     xy
#     y^2
#     x^3
#     x^2 y
#     x y^2
#     y^3
#
# and combinations of these with exact rational coefficients.
#
# The search is performed by exact linear algebra.
#
#
# Then the recovered coefficient vectors are compared.
#
# A transformation between A and B is considered interesting when:
#
#     coeff(B)
#
# can be obtained from
#
#     coeff(A)
#
# through a small static affine transformation.
#
#
# ==============================================================================

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from collections import defaultdict

import sympy as sp


# ==============================================================================
# CONFIGURATION
# ==============================================================================

MAX_N = 500_000

# Number of points used to discover a function.
#
# More points are then used for independent validation.
DISCOVERY_POINTS = 40
VALIDATION_POINTS = 1000

# Maximum polynomial degree.
MAX_DEGREE = 3

# Candidate basis.
#
# Degree <= 3 bivariate polynomial.
BASIS_TERMS = (
    (0, 0),  # 1
    (1, 0),  # x
    (0, 1),  # y
    (2, 0),  # x^2
    (1, 1),  # xy
    (0, 2),  # y^2
    (3, 0),  # x^3
    (2, 1),  # x^2 y
    (1, 2),  # x y^2
    (0, 3),  # y^3
)


# ==============================================================================
# DATA STRUCTURE
# ==============================================================================

@dataclass(frozen=True)
class Point:

    n: int
    p: int
    q: int
    x: int
    y: int
    residue4: int
    residue8: int


@dataclass
class DiscoveredFunction:

    name: str
    coefficients: list[Fraction]
    terms: tuple
    discovery_count: int
    validation_count: int


# ==============================================================================
# 2020 COORDINATES
# ==============================================================================

def xy_from_pq(
    p: int,
    q: int,
):
    """
    Correct algebraic crosswalk:

        p = 2y - 2x + 3
        q = 2y + 2x - 3

    hence:

        x = (q-p+6)/4
        y = (p+q)/4
    """

    x_num = q - p + 6
    y_num = p + q

    if x_num % 4 != 0:
        return None

    if y_num % 4 != 0:
        return None

    return (
        x_num // 4,
        y_num // 4,
    )


# ==============================================================================
# PRIME SIEVE
# ==============================================================================

def prime_sieve(limit):

    sieve = bytearray(
        b"\x01"
    ) * (
        limit + 1
    )

    sieve[0] = 0
    sieve[1] = 0

    root = int(limit ** 0.5)

    for p in range(
        2,
        root + 1,
    ):

        if not sieve[p]:
            continue

        start = p * p

        count = (
            (limit - start) // p
        ) + 1

        sieve[
            start:limit + 1:p
        ] = b"\x00" * count

    return sieve


# ==============================================================================
# SEMIPRIMES
# ==============================================================================

def generate_points():

    sieve = prime_sieve(
        MAX_N
    )

    primes = [
        p
        for p in range(
            3,
            MAX_N + 1,
            2,
        )
        if sieve[p]
    ]

    points = []

    for i, p in enumerate(primes):

        if p * p > MAX_N:
            break

        max_q = MAX_N // p

        for q in primes[i:]:

            if q > max_q:
                break

            n = p * q

            coords = xy_from_pq(
                p,
                q,
            )

            if coords is None:
                continue

            x, y = coords

            points.append(
                Point(
                    n=n,
                    p=p,
                    q=q,
                    x=x,
                    y=y,
                    residue4=n % 4,
                    residue8=n % 8,
                )
            )

    return points


# ==============================================================================
# FUNCTION TARGET
# ==============================================================================

def target_value(
    point: Point,
):
    """
    The natural function we are recovering is p.

        p = A_i(x,y)
        p = B_i(x,y)
    """

    return Fraction(
        point.p
    )


# ==============================================================================
# BASIS
# ==============================================================================

def basis_value(
    x,
    y,
    term,
):

    a, b = term

    return (
        Fraction(x) ** a
        * Fraction(y) ** b
    )


def design_matrix(
    points,
):

    return sp.Matrix(
        [
            [
                sp.Rational(
                    basis_value(
                        point.x,
                        point.y,
                        term,
                    ).numerator,
                    basis_value(
                        point.x,
                        point.y,
                        term,
                    ).denominator,
                )
                for term in BASIS_TERMS
            ]
            for point in points
        ]
    )


def target_matrix(
    points,
):

    return sp.Matrix(
        [
            [
                sp.Rational(
                    target_value(
                        point
                    ).numerator,
                    target_value(
                        point
                    ).denominator,
                )
            ]
            for point in points
        ]
    )


# ==============================================================================
# EXACT POLYNOMIAL RECOVERY
# ==============================================================================

def recover_polynomial(
    points,
):

    if len(points) < len(BASIS_TERMS):
        return None

    A = design_matrix(
        points
    )

    b = target_matrix(
        points
    )

    try:
        solution = sp.linsolve(
            (
                A,
                b,
            )
        )

    except Exception:
        return None

    if solution == sp.EmptySet:
        return None

    solutions = list(
        solution
    )

    if len(solutions) != 1:
        return None

    vector = list(
        solutions[0]
    )

    # Reject underdetermined symbolic solutions.
    if any(
        value.free_symbols
        for value in vector
    ):
        return None

    coeffs = []

    for value in vector:

        coeffs.append(
            Fraction(
                int(value.p),
                int(value.q),
            )
        )

    return coeffs


# ==============================================================================
# EVALUATION
# ==============================================================================

def evaluate_polynomial(
    coefficients,
    x,
    y,
):

    total = Fraction(0)

    for coefficient, term in zip(
        coefficients,
        BASIS_TERMS,
    ):

        total += (
            coefficient
            * basis_value(
                x,
                y,
                term,
            )
        )

    return total


# ==============================================================================
# VALIDATION
# ==============================================================================

def validate_function(
    coefficients,
    points,
):

    correct = 0

    for point in points:

        predicted = evaluate_polynomial(
            coefficients,
            point.x,
            point.y,
        )

        actual = target_value(
            point
        )

        if predicted == actual:
            correct += 1

    return correct


# ==============================================================================
# SPLIT BY MOD 8
# ==============================================================================

def split_mod8(
    points,
):

    groups = {
        1: [],
        3: [],
        5: [],
        7: [],
    }

    for point in points:

        groups[
            point.residue8
        ].append(
            point
        )

    return groups


# ==============================================================================
# DISCOVER FUNCTION FOR A MOD-8 CLASS
# ==============================================================================

def discover_function(
    name,
    points,
):

    if len(points) < DISCOVERY_POINTS:
        return None

    discovery = points[
        :DISCOVERY_POINTS
    ]

    validation = points[
        DISCOVERY_POINTS:
        DISCOVERY_POINTS
        + VALIDATION_POINTS
    ]

    coefficients = recover_polynomial(
        discovery
    )

    if coefficients is None:
        return None

    validation_count = validate_function(
        coefficients,
        validation,
    )

    return DiscoveredFunction(
        name=name,
        coefficients=coefficients,
        terms=BASIS_TERMS,
        discovery_count=len(discovery),
        validation_count=validation_count,
    )


# ==============================================================================
# PRETTY POLYNOMIAL
# ==============================================================================

def polynomial_string(
    coefficients,
):

    pieces = []

    for coefficient, term in zip(
        coefficients,
        BASIS_TERMS,
    ):

        if coefficient == 0:
            continue

        a, b = term

        if a == 0 and b == 0:
            monomial = "1"

        elif a == 0:
            if b == 1:
                monomial = "y"
            else:
                monomial = f"y^{b}"

        elif b == 0:
            if a == 1:
                monomial = "x"
            else:
                monomial = f"x^{a}"

        else:

            if a == 1:
                xa = "x"
            else:
                xa = f"x^{a}"

            if b == 1:
                yb = "y"
            else:
                yb = f"y^{b}"

            monomial = (
                f"{xa}{yb}"
            )

        if coefficient == 1:
            term_string = monomial

        elif coefficient == -1:
            term_string = (
                "-"
                + monomial
            )

        else:
            term_string = (
                f"({coefficient})"
                + monomial
            )

        pieces.append(
            term_string
        )

    if not pieces:
        return "0"

    result = pieces[0]

    for piece in pieces[1:]:

        if piece.startswith("-"):
            result += " - "
            result += piece[1:]

        else:
            result += " + "
            result += piece

    return result


# ==============================================================================
# COEFFICIENT DIFFERENCE
# ==============================================================================

def coefficient_difference(
    parent,
    child,
):

    return [
        c - p
        for p, c in zip(
            parent,
            child,
        )
    ]


# ==============================================================================
# COEFFICIENT RATIO
# ==============================================================================

def coefficient_ratio(
    parent,
    child,
):

    ratios = []

    for p, c in zip(
        parent,
        child,
    ):

        if p == 0:
            ratios.append(
                None
            )
        else:
            ratios.append(
                c / p
            )

    return ratios


# ==============================================================================
# SIMPLE VARIABLE SUBSTITUTION SEARCH
# ==============================================================================

def transformed_polynomial(
    coefficients,
    x_sub,
    y_sub,
):
    """
    Substitute:

        x -> x_sub
        y -> y_sub

    into polynomial.

    x_sub/y_sub are represented as SymPy expressions.
    """

    x, y = sp.symbols(
        "x y"
    )

    expression = 0

    for coefficient, term in zip(
        coefficients,
        BASIS_TERMS,
    ):

        a, b = term

        expression += (
            sp.Rational(
                coefficient.numerator,
                coefficient.denominator,
            )
            * x_sub ** a
            * y_sub ** b
        )

    return sp.expand(
        expression
    )


def coefficient_vector_from_expression(
    expression,
):

    x, y = sp.symbols(
        "x y"
    )

    poly = sp.Poly(
        sp.expand(
            expression
        ),
        x,
        y,
    )

    coefficients = []

    for term in BASIS_TERMS:

        coefficient = poly.coeff_monomial(
            x ** term[0]
            * y ** term[1]
        )

        coefficients.append(
            Fraction(
                int(coefficient.p),
                int(coefficient.q),
            )
        )

    return coefficients


# ==============================================================================
# STATIC TRANSFORMATION SEARCH
# ==============================================================================

def search_static_transformations(
    parent_function,
    child_function,
):

    x, y = sp.symbols(
        "x y"
    )

    candidates = []

    # Small affine transformations.
    values = (
        -2,
        -1,
        0,
        1,
        2,
    )

    transformations = []

    for ax, bx, cx in product(
        values,
        repeat=3,
    ):

        # Avoid degenerate x mappings.
        if ax == 0:
            continue

        x_sub = (
            ax * x
            + bx * y
            + cx
        )

        for ay, by, cy in product(
            values,
            repeat=3,
        ):

            if by == 0 and ay == 0:
                continue

            y_sub = (
                ay * x
                + by * y
                + cy
            )

            transformations.append(
                (
                    x_sub,
                    y_sub,
                )
            )

    # Include global output affine transformation:
    #
    #     child = a * parent(substitution) + b
    #
    for x_sub, y_sub in transformations:

        transformed = transformed_polynomial(
            parent_function.coefficients,
            x_sub,
            y_sub,
        )

        transformed_poly = sp.Poly(
            transformed,
            x,
            y,
        )

        child_expression = sp.Add(
            *[
                sp.Rational(
                    coefficient.numerator,
                    coefficient.denominator,
                )
                * x ** term[0]
                * y ** term[1]
                for coefficient, term in zip(
                    child_function.coefficients,
                    BASIS_TERMS,
                )
            ]
        )

        # Try small output scales.
        for a in (
            -4,
            -2,
            -1,
            Fraction(1, 4),
            Fraction(1, 2),
            Fraction(1),
            Fraction(2),
            Fraction(4),
        ):

            for b in (
                -8,
                -4,
                -2,
                -1,
                0,
                1,
                2,
                4,
                8,
            ):

                candidate = sp.expand(
                    sp.Rational(
                        a.numerator
                        if isinstance(a, Fraction)
                        else a,
                        a.denominator
                        if isinstance(a, Fraction)
                        else 1,
                    )
                    * transformed
                    + b
                )

                if sp.expand(
                    candidate
                    - child_expression
                ) == 0:

                    candidates.append(
                        (
                            x_sub,
                            y_sub,
                            a,
                            b,
                        )
                    )

                    return candidates

    return candidates


# ==============================================================================
# REPORT FUNCTION
# ==============================================================================

def print_function(
    function,
    point_count,
):

    if function is None:

        print(
            "    FUNCTION DISCOVERY FAILED"
        )

        return

    expression = polynomial_string(
        function.coefficients
    )

    print(
        f"    {function.name}(x,y)"
        f" = {expression}"
    )

    print(
        f"    discovery points = "
        f"{function.discovery_count}"
    )

    print(
        f"    validation = "
        f"{function.validation_count}/"
        f"{point_count}"
    )


# ==============================================================================
# MOD-8 FUNCTION DISCOVERY
# ==============================================================================

def discover_all_B(
    points,
):

    groups = split_mod8(
        points
    )

    functions = {}

    print()
    print("=" * 90)
    print(
        "MOD-8 FUNCTION DISCOVERY"
    )
    print("=" * 90)

    for residue in (
        1,
        3,
        5,
        7,
    ):

        subset = groups[
            residue
        ]

        name = {
            1: "B1",
            3: "B2",
            5: "B3",
            7: "B4",
        }[residue]

        print()
        print(
            f"{name}: "
            f"n = {residue} mod 8"
        )

        print(
            f"    samples = {len(subset)}"
        )

        function = discover_function(
            name,
            subset,
        )

        functions[name] = function

        validation_count = max(
            0,
            min(
                VALIDATION_POINTS,
                len(subset)
                - DISCOVERY_POINTS,
            ),
        )

        print_function(
            function,
            validation_count,
        )

    return functions


# ==============================================================================
# A FUNCTIONS
# ==============================================================================

def A1(
    x,
    y,
):
    """
    Existing A1 must be inserted here.

    This placeholder is intentionally explicit.
    """

    raise NotImplementedError(
        "Insert the already-known A1(x,y) function."
    )


def A2(
    x,
    y,
):
    """
    Existing A2 must be inserted here.

    This placeholder is intentionally explicit.
    """

    raise NotImplementedError(
        "Insert the already-known A2(x,y) function."
    )


# ==============================================================================
# VERIFY EXISTING A FUNCTIONS
# ==============================================================================

def verify_A_functions(
    points,
):

    """
    Verify the existing A1/A2 functions against data.

    This assumes the user already has the formulas.
    """

    print()
    print("=" * 90)
    print(
        "VERIFY EXISTING A FUNCTIONS"
    )
    print("=" * 90)

    groups = defaultdict(list)

    for point in points:

        groups[
            point.residue4
        ].append(
            point
        )

    for residue, name, function in (
        (1, "A1", A1),
        (3, "A2", A2),
    ):

        subset = groups[
            residue
        ]

        print()
        print(
            f"{name}: "
            f"n = {residue} mod 4"
        )

        if not subset:

            print(
                "    no data"
            )

            continue

        correct = 0

        for point in subset:

            predicted = Fraction(
                function(
                    point.x,
                    point.y,
                )
            )

            if predicted == point.p:
                correct += 1

        print(
            f"    validation = "
            f"{correct}/{len(subset)}"
        )


# ==============================================================================
# TRANSFORMATION SEARCH
# ==============================================================================

def compare_A_to_B(
    A_function,
    B_function,
    A_name,
    B_name,
):

    print()
    print(
        f"TRANSFORMATION SEARCH "
        f"{A_name} -> {B_name}"
    )

    if A_function is None:
        print(
            "    parent function unavailable"
        )
        return

    if B_function is None:
        print(
            "    child function unavailable"
        )
        return

    candidates = search_static_transformations(
        A_function,
        B_function,
    )

    if not candidates:

        print(
            "    NO SMALL STATIC TRANSFORMATION FOUND"
        )

        return

    for x_sub, y_sub, scale, shift in candidates:

        print(
            "    FOUND:"
        )

        print(
            f"        x -> {x_sub}"
        )

        print(
            f"        y -> {y_sub}"
        )

        print(
            f"        output -> "
            f"({scale})*A + ({shift})"
        )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print(
        "EXPERIMENT 578 START"
    )
    print("=" * 90)

    print()
    print(
        "DISCOVER MOD-8 FUNCTIONS "
        "AND TRANSFORM A -> B"
    )

    # --------------------------------------------------------------------------
    # Generate data.
    # --------------------------------------------------------------------------

    print()
    print(
        "[1] GENERATING DATA"
    )

    points = generate_points()

    print(
        f"Generated points = {len(points)}"
    )

    if not points:
        return

    # --------------------------------------------------------------------------
    # Show residue availability.
    # --------------------------------------------------------------------------

    groups = split_mod8(
        points
    )

    print()
    print(
        "[2] MOD-8 POPULATION"
    )

    for residue in (
        1,
        3,
        5,
        7,
    ):

        print(
            f"    n = {residue} mod 8: "
            f"{len(groups[residue])}"
        )

    # --------------------------------------------------------------------------
    # Discover B functions.
    # --------------------------------------------------------------------------

    B = discover_all_B(
        points
    )

    # --------------------------------------------------------------------------
    # Verify existing A functions.
    # --------------------------------------------------------------------------

    #
    # Uncomment after inserting your existing A1/A2 formulas.
    #
    # verify_A_functions(points)
    #

    # --------------------------------------------------------------------------
    # A -> B transformations.
    # --------------------------------------------------------------------------

    #
    # These correspond to:
    #
    #     A1 -> B1
    #     A1 -> B2
    #
    #     A2 -> B3
    #     A2 -> B4
    #
    #
    # The existing A1/A2 equations must be present for these searches.
    #

    try:

        A1_function = None
        A2_function = None

        #
        # The experiment intentionally does not invent A1/A2.
        #
        # The known equations should be represented as
        # DiscoveredFunction objects here.
        #

    except Exception:

        A1_function = None
        A2_function = None

    print()
    print("=" * 90)
    print(
        "A -> B TRANSFORMATION STAGE"
    )
    print("=" * 90)

    print()
    print(
        "B functions have now been discovered independently."
    )

    print()
    print(
        "Next comparison should be:"
    )

    print(
        "    A1 -> B1"
    )

    print(
        "    A1 -> B2"
    )

    print(
        "    A2 -> B3"
    )

    print(
        "    A2 -> B4"
    )

    print()
    print(
        "The transformation search must use the actual existing"
    )

    print(
        "A1 and A2 equations; they are not guessed here."
    )

    # --------------------------------------------------------------------------
    # Final.
    # --------------------------------------------------------------------------

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 578 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()
