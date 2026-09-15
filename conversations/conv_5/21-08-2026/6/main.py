#!/usr/bin/env python3

# ==============================================================================
# EXPERIMENT 579
# ==============================================================================
# MOD-8 FUNCTION RECOVERY — RANK-AWARE
#
# Goal
# ----
#
# Recover p = B_i(x,y) independently for:
#
#     B1 : n = 1 mod 8
#     B2 : n = 3 mod 8
#     B3 : n = 5 mod 8
#     B4 : n = 7 mod 8
#
# IMPORTANT:
#
# This experiment does NOT assume that the first N samples are suitable for
# interpolation.
#
# It greedily selects points which increase the rank of the polynomial design
# matrix.
#
# This prevents a data-order artifact such as:
#
#     p=3
#     x=y
#
# from making the system appear underdetermined.
#
#
# It also explicitly tests whether the recovered B_i are genuinely different.
#
# If:
#
#     B2 = B4
#
# under the current coordinate system, that means the current x,y construction
# does NOT implement the branch-dependent function family from the 2020 theory.
#
# ==============================================================================

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

import sympy as sp


# ==============================================================================
# CONFIG
# ==============================================================================

MAX_N = 500_000

# Polynomial basis.
#
# Start with degree 1 because the known coordinate relation is linear.
#
# If the actual 2020 functions are nonlinear, increase this later.
BASIS = (
    (0, 0),  # 1
    (1, 0),  # x
    (0, 1),  # y
)

VALIDATION_POINTS = 2000


# ==============================================================================
# SAMPLE
# ==============================================================================

@dataclass(frozen=True)
class Point:
    n: int
    p: int
    q: int
    x: int
    y: int


# ==============================================================================
# PRIME SIEVE
# ==============================================================================

def prime_sieve(limit):

    sieve = bytearray(
        b"\x01"
    ) * (limit + 1)

    sieve[0] = 0
    sieve[1] = 0

    root = int(limit ** 0.5)

    for p in range(2, root + 1):

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
# CURRENT 2020 CROSSWALK
# ==============================================================================

def xy_from_pq(p, q):

    x_num = q - p + 6
    y_num = p + q

    if x_num % 4 != 0:
        return None

    if y_num % 4 != 0:
        return None

    x = x_num // 4
    y = y_num // 4

    return x, y


# ==============================================================================
# GENERATE POINTS
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
                )
            )

    return points


# ==============================================================================
# BASIS EVALUATION
# ==============================================================================

def basis_row(point):

    return [
        sp.Integer(1),
        sp.Integer(point.x),
        sp.Integer(point.y),
    ]


# ==============================================================================
# RANK-AWARE SAMPLE SELECTION
# ==============================================================================

def select_independent_points(
    points,
    required_rank,
):

    selected = []
    current_matrix = None
    current_rank = 0

    # Diversity heuristic:
    #
    # First prioritize different p values, x values and y values.
    #
    # This is not mathematical filtering; it simply improves interpolation
    # conditioning.

    ordered = sorted(
        points,
        key=lambda p: (
            p.p,
            p.x,
            p.y,
        ),
    )

    # Shuffle-like deterministic reordering:
    #
    # use multiple factor branches rather than all p=3 first.
    ordered = sorted(
        points,
        key=lambda p: (
            p.p,
            p.q % 997,
            p.x,
            p.y,
        ),
    )

    for point in ordered:

        candidate = selected + [point]

        matrix = sp.Matrix(
            [
                basis_row(p)
                for p in candidate
            ]
        )

        rank = matrix.rank()

        if rank > current_rank:

            selected.append(point)
            current_rank = rank
            current_matrix = matrix

        if current_rank >= required_rank:
            break

    return selected


# ==============================================================================
# EXACT FIT
# ==============================================================================

def recover_linear_function(
    points,
):

    if len(points) < len(BASIS):
        return None

    matrix = sp.Matrix(
        [
            basis_row(p)
            for p in points
        ]
    )

    target = sp.Matrix(
        [
            p.p
            for p in points
        ]
    )

    if matrix.rank() < len(BASIS):
        return None

    solution = sp.linsolve(
        (
            matrix,
            target,
        )
    )

    solutions = list(solution)

    if len(solutions) != 1:
        return None

    vector = list(
        solutions[0]
    )

    if any(
        value.free_symbols
        for value in vector
    ):
        return None

    return [
        Fraction(
            int(value.p),
            int(value.q),
        )
        for value in vector
    ]


# ==============================================================================
# EVALUATE
# ==============================================================================

def evaluate(
    coefficients,
    point,
):

    result = Fraction(0)

    for coefficient, term in zip(
        coefficients,
        BASIS,
    ):

        a, b = term

        result += (
            coefficient
            * Fraction(point.x) ** a
            * Fraction(point.y) ** b
        )

    return result


# ==============================================================================
# VALIDATION
# ==============================================================================

def validate(
    coefficients,
    points,
):

    correct = 0

    for point in points:

        if evaluate(
            coefficients,
            point,
        ) == point.p:

            correct += 1

    return correct


# ==============================================================================
# PRETTY PRINT
# ==============================================================================

def expression(
    coefficients,
):

    parts = []

    names = (
        "1",
        "x",
        "y",
    )

    for coefficient, name in zip(
        coefficients,
        names,
    ):

        if coefficient == 0:
            continue

        if name == "1":

            parts.append(
                str(coefficient)
            )

            continue

        if coefficient == 1:
            parts.append(
                name
            )

        elif coefficient == -1:
            parts.append(
                "-" + name
            )

        else:
            parts.append(
                f"({coefficient}){name}"
            )

    if not parts:
        return "0"

    result = parts[0]

    for part in parts[1:]:

        if part.startswith("-"):

            result += " - "
            result += part[1:]

        else:

            result += " + "
            result += part

    return result


# ==============================================================================
# RESIDUE GROUPS
# ==============================================================================

def split_residue(
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
            point.n % 8
        ].append(
            point
        )

    return groups


# ==============================================================================
# DISCOVER B
# ==============================================================================

def discover_B(
    residue,
    points,
):

    print()
    print(
        f"B{ {1:1, 3:2, 5:3, 7:4}[residue] } "
        f": n = {residue} mod 8"
    )

    print(
        f"    total points = {len(points)}"
    )

    if not points:

        print(
            "    NO DATA"
        )

        return None

    # --------------------------------------------------------------
    # Rank-aware selection
    # --------------------------------------------------------------

    selected = select_independent_points(
        points,
        len(BASIS),
    )

    print(
        "    independent points =",
        len(selected),
    )

    for point in selected:

        print(
            f"        n={point.n:<8} "
            f"p={point.p:<6} "
            f"q={point.q:<6} "
            f"x={point.x:<6} "
            f"y={point.y:<6}"
        )

    coefficients = recover_linear_function(
        selected
    )

    if coefficients is None:

        print(
            "    LINEAR FUNCTION NOT RECOVERED"
        )

        return None

    print(
        "    DISCOVERED:"
    )

    print(
        f"        p = {expression(coefficients)}"
    )

    # --------------------------------------------------------------
    # Validation
    # --------------------------------------------------------------

    validation_points = points[
        len(selected):
        len(selected)
        + VALIDATION_POINTS
    ]

    if not validation_points:
        validation_points = points

    correct = validate(
        coefficients,
        validation_points,
    )

    print(
        f"    validation = "
        f"{correct}/"
        f"{len(validation_points)}"
    )

    # --------------------------------------------------------------
    # Search the whole dataset
    # --------------------------------------------------------------

    total_correct = validate(
        coefficients,
        points,
    )

    print(
        f"    FULL validation = "
        f"{total_correct}/"
        f"{len(points)}"
    )

    return coefficients


# ==============================================================================
# FUNCTION COMPARISON
# ==============================================================================

def compare_functions(
    functions,
):

    print()
    print("=" * 90)
    print(
        "B-FUNCTION COMPARISON"
    )
    print("=" * 90)

    names = sorted(
        functions
    )

    for i in range(
        len(names)
    ):

        for j in range(
            i + 1,
            len(names),
        ):

            a = names[i]
            b = names[j]

            fa = functions[a]
            fb = functions[b]

            if fa is None or fb is None:
                continue

            identical = (
                fa == fb
            )

            print(
                f"    {a} vs {b}: "
                f"{'IDENTICAL' if identical else 'DIFFERENT'}"
            )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print(
        "EXPERIMENT 579 START"
    )
    print("=" * 90)

    print()
    print(
        "RANK-AWARE MOD-8 FUNCTION RECOVERY"
    )

    print()

    points = generate_points()

    print(
        f"Generated points = {len(points)}"
    )

    groups = split_residue(
        points
    )

    print()
    print(
        "MOD-8 POPULATION"
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

    functions = {}

    print()
    print("=" * 90)
    print(
        "INDEPENDENT B-FUNCTION RECOVERY"
    )
    print("=" * 90)

    for residue in (
        1,
        3,
        5,
        7,
    ):

        name = {
            1: "B1",
            3: "B2",
            5: "B3",
            7: "B4",
        }[residue]

        functions[name] = discover_B(
            residue,
            groups[residue],
        )

    compare_functions(
        functions
    )

    print()
    print("=" * 90)
    print(
        "IMPORTANT INTERPRETATION"
    )
    print("=" * 90)

    print(
        """
This experiment deliberately separates:

    DATA AVAILABILITY

from:

    FUNCTION DISCOVERY

If B2 and B4 recover as:

    p = 2y - 2x + 3

then the CURRENT x,y coordinate system has the same p-function
on both mod-8 branches.

That does not disprove the 2020 theory.

It means this particular coordinate reconstruction has already
absorbed the residue distinction.

For the intended A -> B experiment we need the ORIGINAL 2020
branch-specific x,y definitions, i.e. the functions that produced
your existing A1 and A2.

The correct workflow is then:

    A1,A2 known
          |
          v
    collect mod-8 data
          |
          v
    independently discover B1,B2,B3,B4
          |
          v
    compare coefficient/function structure
          |
          v
    search static transformations
          |
          v
    A -> B rule

Do NOT infer the transformation before recovering B independently.
"""
    )

    print()
    print(
        "EXPERIMENT 579 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()
