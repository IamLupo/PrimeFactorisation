#!/usr/bin/env python3

# ==============================================================================
# EXPERIMENT 591
# ==============================================================================
# DISCOVER THE ACTUAL A_r -> B_r COORDINATE/FUNCTION RECURSION
#
# NO FILES
# NO WEB
#
# The previous experiment established that only certain residues admit an
# integral normalized coordinate frame at each level.
#
# We now ignore all other residues and study ONLY the active branches.
#
# For every active residue r at level z:
#
#     n mod 2^z = r
#
# we search for exact equations:
#
#     p = a*x + b*y + c
#     q = d*x + e*y + f
#
# where x,y are the normalized coordinates.
#
# We then compare:
#
#     parent residue r
#         ->
#     child residue r or r+2^z
#
# and search for the coordinate transformation:
#
#     x_child = A*x_parent + B*y_parent + C
#     y_child = D*x_parent + E*y_parent + F
#
# The main question is whether the transformation is always:
#
#     x_child = x_parent / 2
#     y_child = y_parent / 2
#
# while the FUNCTION COEFFICIENTS change predictably.
#
# ==============================================================================

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from math import isqrt


# ==============================================================================
# CONFIG
# ==============================================================================

MAX_N = 2_000_000

MIN_Z = 2
MAX_Z = 9

MAX_COEFF = 16
MAX_OFFSET = 64

MAX_EXAMPLES = 8


# ==============================================================================
# DATA
# ==============================================================================

@dataclass(frozen=True)
class Point:
    n: int
    p: int
    q: int
    x: Fraction
    y: Fraction
    residue: int


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

    root = isqrt(limit)

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

def generate_semiprimes(
    max_n,
    sieve,
):

    primes = [
        p
        for p in range(
            3,
            max_n + 1,
            2,
        )
        if sieve[p]
    ]

    result = []

    for i, p in enumerate(primes):

        if p * p > max_n:
            break

        max_q = max_n // p

        for q in primes[i:]:

            if q > max_q:
                break

            result.append(
                (
                    p * q,
                    p,
                    q,
                )
            )

    return result


# ==============================================================================
# BASE INVARIANT COORDINATES
# ==============================================================================
#
# These are the two universal linear combinations:
#
#     U = (q-p+6)/2
#     V = (p+q)/2
#
# They correspond to:
#
#     U = X
#     V = Y
#
# in the frame:
#
#     p = Y-X+3
#     q = Y+X-3.
#
# The important point is:
#
#     U,V are invariant.
#
# We only divide them by k when the resulting coordinates are integers.
#
# ==============================================================================

def invariant_coordinates(
    p,
    q,
):

    if (
        (q - p + 6) % 2
        or
        (p + q) % 2
    ):
        return None

    X = Fraction(
        q - p + 6,
        2,
    )

    Y = Fraction(
        p + q,
        2,
    )

    return X, Y


# ==============================================================================
# ACTIVE NORMALIZED POINT
# ==============================================================================
#
# x = X / 2^(z-1)
# y = Y / 2^(z-1)
#
# Only retain integral x,y.
#
# ==============================================================================

def normalized_point(
    n,
    p,
    q,
    z,
):

    coords = invariant_coordinates(
        p,
        q,
    )

    if coords is None:
        return None

    X, Y = coords

    k = 1 << (
        z - 1
    )

    x = X / k
    y = Y / k

    if (
        x.denominator != 1
        or
        y.denominator != 1
    ):
        return None

    return Point(
        n=n,
        p=p,
        q=q,
        x=x,
        y=y,
        residue=n % (1 << z),
    )


# ==============================================================================
# BUILD ACTIVE LEVELS
# ==============================================================================

def build_levels(
    semiprimes,
):

    levels = {}

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        level = defaultdict(list)

        for n, p, q in semiprimes:

            point = normalized_point(
                n,
                p,
                q,
                z,
            )

            if point is None:
                continue

            level[
                point.residue
            ].append(
                point
            )

        levels[z] = level

        total = sum(
            len(v)
            for v in level.values()
        )

        print(
            f"    z={z:<2} "
            f"active states={total:<7} "
            f"active residues="
            f"{len(level)}"
        )

    return levels


# ==============================================================================
# EXACT ONE-VARIABLE AFFINE FIT
# ==============================================================================

def fit_affine(
    rows,
):
    """
    Find:

        target = a*x + b*y + c

    with small integer coefficients.
    """

    if not rows:
        return None

    first = rows[0]

    x0 = first[0]
    y0 = first[1]
    t0 = first[2]

    for a in range(
        -MAX_COEFF,
        MAX_COEFF + 1,
    ):

        for b in range(
            -MAX_COEFF,
            MAX_COEFF + 1,
        ):

            c = (
                t0
                - a * x0
                - b * y0
            )

            if (
                c.denominator != 1
                or
                abs(c.numerator)
                > MAX_OFFSET
            ):
                continue

            valid = True

            for x, y, target in rows:

                if (
                    a * x
                    + b * y
                    + c
                    != target
                ):

                    valid = False
                    break

            if valid:
                return (
                    a,
                    b,
                    c,
                )

    return None


# ==============================================================================
# FUNCTION DISCOVERY
# ==============================================================================

def discover_functions(
    points,
):

    p_rows = [
        (
            point.x,
            point.y,
            Fraction(point.p),
        )
        for point in points
    ]

    q_rows = [
        (
            point.x,
            point.y,
            Fraction(point.q),
        )
        for point in points
    ]

    p_rule = fit_affine(
        p_rows
    )

    q_rule = fit_affine(
        q_rows
    )

    return p_rule, q_rule


# ==============================================================================
# PRINT FUNCTION
# ==============================================================================

def affine_text(
    name,
    rule,
):

    if rule is None:
        return (
            f"{name} = NO EXACT RULE"
        )

    a, b, c = rule

    return (
        f"{name} = "
        f"({a})x + "
        f"({b})y + "
        f"({c})"
    )


# ==============================================================================
# ACTIVE RESIDUE FUNCTION REPORT
# ==============================================================================

def report_functions(
    levels,
):

    print()
    print("=" * 90)
    print(
        "ACTIVE RESIDUE FUNCTION DISCOVERY"
    )
    print("=" * 90)

    discovered = {}

    for z in sorted(
        levels
    ):

        discovered[z] = {}

        print()
        print(
            f"LEVEL z={z}"
        )

        for residue in sorted(
            levels[z]
        ):

            points = levels[z][
                residue
            ]

            p_rule, q_rule = (
                discover_functions(
                    points
                )
            )

            discovered[z][
                residue
            ] = (
                p_rule,
                q_rule,
            )

            print()
            print(
                f"    residue={residue:<5} "
                f"samples={len(points)}"
            )

            print(
                "        "
                + affine_text(
                    "p",
                    p_rule,
                )
            )

            print(
                "        "
                + affine_text(
                    "q",
                    q_rule,
                )
            )

    return discovered


# ==============================================================================
# SAME-N CHILD TRANSITIONS
# ==============================================================================

def build_transitions(
    levels,
):

    transitions = defaultdict(list)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        parent_points = {}

        for residue in levels[z]:

            for point in levels[z][
                residue
            ]:

                parent_points[
                    point.n
                ] = point

        child_points = {}

        for residue in levels[z + 1]:

            for point in levels[z + 1][
                residue
            ]:

                child_points[
                    point.n
                ] = point

        for n in (
            set(parent_points)
            &
            set(child_points)
        ):

            p = parent_points[n]
            c = child_points[n]

            transitions[
                (
                    z,
                    p.residue,
                    c.residue,
                )
            ].append(
                (
                    p,
                    c,
                )
            )

    return transitions


# ==============================================================================
# COORDINATE TRANSFORMATION SEARCH
# ==============================================================================

def fit_coordinate_transform(
    pairs,
):

    if not pairs:
        return None

    rows_x = [
        (
            parent.x,
            parent.y,
            child.x,
        )
        for parent, child in pairs
    ]

    rows_y = [
        (
            parent.x,
            parent.y,
            child.y,
        )
        for parent, child in pairs
    ]

    x_rule = fit_affine(
        rows_x
    )

    y_rule = fit_affine(
        rows_y
    )

    if (
        x_rule is None
        or
        y_rule is None
    ):
        return None

    return (
        x_rule,
        y_rule,
    )


# ==============================================================================
# TRANSITION REPORT
# ==============================================================================

def report_transforms(
    transitions,
):

    print()
    print("=" * 90)
    print(
        "RESIDUE-SPECIFIC COORDINATE TRANSFORMATIONS"
    )
    print("=" * 90)

    for (
        z,
        parent_residue,
        child_residue,
    ) in sorted(
        transitions
    ):

        pairs = transitions[
            (
                z,
                parent_residue,
                child_residue,
            )
        ]

        result = fit_coordinate_transform(
            pairs
        )

        print()
        print(
            f"LEVEL {z} -> {z+1}"
        )

        print(
            f"    residue "
            f"{parent_residue}"
            f" -> "
            f"{child_residue}"
        )

        print(
            f"    states={len(pairs)}"
        )

        if result is None:

            print(
                "    NO EXACT AFFINE "
                "TRANSFORMATION"
            )

            continue

        x_rule, y_rule = result

        print(
            "    "
            + affine_text(
                "x_child",
                x_rule,
            )
        )

        print(
            "    "
            + affine_text(
                "y_child",
                y_rule,
            )
        )


# ==============================================================================
# FUNCTION COEFFICIENT TRANSITIONS
# ==============================================================================

def report_function_transitions(
    discovered,
    transitions,
):

    print()
    print("=" * 90)
    print(
        "FUNCTION TRANSITION A -> B"
    )
    print("=" * 90)

    for (
        z,
        parent_residue,
        child_residue,
    ) in sorted(
        transitions
    ):

        parent_rule = discovered[
            z
        ].get(
            parent_residue
        )

        child_rule = discovered[
            z + 1
        ].get(
            child_residue
        )

        print()
        print(
            f"{parent_residue}"
            f" mod {1 << z}"
            f" -> "
            f"{child_residue}"
            f" mod {1 << (z+1)}"
        )

        if parent_rule is None:

            print(
                "    parent rule missing"
            )

        else:

            print(
                "    PARENT:"
            )

            print(
                "        "
                + affine_text(
                    "p",
                    parent_rule[0],
                )
            )

            print(
                "        "
                + affine_text(
                    "q",
                    parent_rule[1],
                )
            )

        if child_rule is None:

            print(
                "    child rule missing"
            )

        else:

            print(
                "    CHILD:"
            )

            print(
                "        "
                + affine_text(
                    "p",
                    child_rule[0],
                )
            )

            print(
                "        "
                + affine_text(
                    "q",
                    child_rule[1],
                )
            )


# ==============================================================================
# HALVING TEST
# ==============================================================================

def halving_test(
    transitions,
):

    print()
    print("=" * 90)
    print(
        "EXACT HALVING TEST"
    )
    print("=" * 90)

    for key in sorted(
        transitions
    ):

        pairs = transitions[key]

        exact_x = all(
            child.x
            == parent.x / 2
            for parent, child
            in pairs
        )

        exact_y = all(
            child.y
            == parent.y / 2
            for parent, child
            in pairs
        )

        z, pr, cr = key

        print(
            f"z={z} "
            f"{pr}->{cr} "
            f"states={len(pairs)} "
            f"x/2={exact_x} "
            f"y/2={exact_y}"
        )


# ==============================================================================
# SIMPLE TRANSFORMATION CLASSIFICATION
# ==============================================================================

def classify_transform(
    rule,
):

    if rule is None:
        return "NONE"

    a, b, c = rule

    candidates = {
        "HALF": (
            Fraction(1, 2),
            0,
            0,
        ),
        "NEG_HALF": (
            Fraction(-1, 2),
            0,
            0,
        ),
        "Y_HALF": (
            0,
            Fraction(1, 2),
            0,
        ),
        "NEG_Y_HALF": (
            0,
            Fraction(-1, 2),
            0,
        ),
        "IDENTITY": (
            1,
            0,
            0,
        ),
        "NEG_IDENTITY": (
            -1,
            0,
            0,
        ),
    }

    current = (
        a,
        b,
        c,
    )

    for name, candidate in (
        candidates.items()
    ):

        if current == candidate:

            return name

    return "GENERAL"


# ==============================================================================
# COMPACT SUMMARY
# ==============================================================================

def compact_summary(
    transitions,
):

    print()
    print("=" * 90)
    print(
        "COMPACT TRANSITION SUMMARY"
    )
    print("=" * 90)

    for key in sorted(
        transitions
    ):

        pairs = transitions[key]

        result = fit_coordinate_transform(
            pairs
        )

        z, parent_residue, child_residue = key

        if result is None:

            print(
                f"    z={z} "
                f"{parent_residue}->{child_residue} "
                f"TRANSFORM=NONE"
            )

            continue

        x_rule, y_rule = result

        print(
            f"    z={z} "
            f"{parent_residue}->{child_residue} "
            f""
            f"x={classify_transform(x_rule)} "
            f""
            f"y={classify_transform(y_rule)}"
        )


# ==============================================================================
# EXAMPLES
# ==============================================================================

def examples(
    transitions,
):

    print()
    print("=" * 90)
    print(
        "EXAMPLE ACTIVE BRANCH TRANSITIONS"
    )
    print("=" * 90)

    shown = 0

    for key in sorted(
        transitions
    ):

        pairs = transitions[key]

        print()
        print(
            f"z={key[0]} "
            f"{key[1]}->{key[2]}"
        )

        for parent, child in pairs[
            :MAX_EXAMPLES
        ]:

            print(
                f"    n={parent.n} "
                f""
                f"({parent.x},"
                f"{parent.y})"
                f" -> "
                f"({child.x},"
                f"{child.y})"
            )

        shown += 1

        if shown >= MAX_EXAMPLES:
            break


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print(
        "EXPERIMENT 591 START"
    )
    print("=" * 90)

    print()
    print(
        "DISCOVER THE ACTUAL RESIDUE FUNCTION RECURSION"
    )

    # --------------------------------------------------------------------------
    # Generate.
    # --------------------------------------------------------------------------

    print()
    print(
        "[1] PRIME SIEVE"
    )

    sieve = prime_sieve(
        MAX_N
    )

    print(
        "[2] SEMIPRIME GENERATION"
    )

    semiprimes = generate_semiprimes(
        MAX_N,
        sieve,
    )

    print(
        f"    semiprimes={len(semiprimes)}"
    )

    # --------------------------------------------------------------------------
    # Active levels.
    # --------------------------------------------------------------------------

    print()
    print(
        "[3] ACTIVE NORMALIZED LEVELS"
    )

    levels = build_levels(
        semiprimes
    )

    # --------------------------------------------------------------------------
    # Discover functions.
    # --------------------------------------------------------------------------

    discovered = report_functions(
        levels
    )

    # --------------------------------------------------------------------------
    # Same-n transitions.
    # --------------------------------------------------------------------------

    print()
    print(
        "[4] SAME-N TRANSITIONS"
    )

    transitions = build_transitions(
        levels
    )

    for key in sorted(
        transitions
    ):

        print(
            f"    z={key[0]} "
            f"{key[1]}->{key[2]} "
            f"states={len(transitions[key])}"
        )

    # --------------------------------------------------------------------------
    # Coordinate transforms.
    # --------------------------------------------------------------------------

    report_transforms(
        transitions
    )

    # --------------------------------------------------------------------------
    # Function transitions.
    # --------------------------------------------------------------------------

    report_function_transitions(
        discovered,
        transitions,
    )

    # --------------------------------------------------------------------------
    # Halving.
    # --------------------------------------------------------------------------

    halving_test(
        transitions
    )

    # --------------------------------------------------------------------------
    # Compact.
    # --------------------------------------------------------------------------

    compact_summary(
        transitions
    )

    # --------------------------------------------------------------------------
    # Examples.
    # --------------------------------------------------------------------------

    examples(
        transitions
    )

    # --------------------------------------------------------------------------
    # Final.
    # --------------------------------------------------------------------------

    print()
    print("=" * 90)
    print(
        "FINAL AUDIT"
    )
    print("=" * 90)

    print(
        """
The purpose of Experiment 591 is to stop treating the two
factor formulas as globally valid.

Instead:

    residue r at level z
        ->
    its own normalized x,y
        ->
    its own p(x,y), q(x,y).

Then the SAME n is followed into the child residue.

The desired result is a table such as:

    parent residue
        child residue
            coordinate transform
                function coefficients

For example:

    A1 -> B1:
        x' = x/2
        y' = y/2
        p' = ...
        q' = ...

    A1 -> B2:
        x' = ...
        y' = ...
        p' = ...
        q' = ...

The important distinction is that the residue determines the
FUNCTION FRAME, while the same factor pair may have a different
normalized coordinate representation at the next level.

This directly tests the original hierarchical A -> B idea.
"""
    )

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 591 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()
