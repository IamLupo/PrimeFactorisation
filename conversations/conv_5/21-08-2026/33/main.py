#!/usr/bin/env python3

# ==============================================================================
# EXPERIMENT 590
# ==============================================================================
# DISCOVER ONE COORDINATE FRAME PER ODD RESIDUE
#
# NO FILES
# NO WEB
#
# Goal:
#
#     n mod 4
#         -> two residue classes
#
#     n mod 8
#         -> four residue classes
#
#     n mod 16
#         -> eight residue classes
#
# etc.
#
# For every residue class we construct its own coordinate frame and then
# compare the frame at z with the frame at z+1.
#
#
# IMPORTANT:
#
# We DO NOT assume that every residue uses the same A/B equation.
#
# Instead we search two elementary factor-coordinate orientations:
#
# FRAME TYPE 1
#
#     p = Y-X+3
#     q = Y+X-3
#
# FRAME TYPE 2
#
#     3p = Y+X-3
#     q  = Y-X+3
#
# For every state, both coordinate frames can be calculated.
#
# The actual residue tells us which frame belongs to the level.
#
# Then we search the transformation:
#
#     X_child = a*X_parent + b*Y_parent + c
#     Y_child = d*X_parent + e*Y_parent + f
#
# with small integer coefficients.
#
# The interesting result is not merely "x halves".
#
# We want:
#
#     residue r
#         ->
#     coordinate frame F_r
#
# and then:
#
#     F_r at level z
#         ->
#     F_s at level z+1
#
# with a STATIC transformation.
#
# ==============================================================================

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from math import isqrt


# ==============================================================================
# CONFIGURATION
# ==============================================================================

MAX_N = 2_000_000

MIN_Z = 2
MAX_Z = 9

COEFFS = (
    -3,
    -2,
    -1,
    0,
    1,
    2,
    3,
)

OFFSETS = range(
    -15,
    16,
)

MAX_EXAMPLES = 10


# ==============================================================================
# DATA
# ==============================================================================

@dataclass(frozen=True)
class State:

    n: int
    p: int
    q: int

    z: int
    residue: int

    # Coordinate frame 1
    X1: int
    Y1: int

    # Coordinate frame 2
    X2: int
    Y2: int


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
# SEMIPRIME GENERATION
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
# FRAME TYPE 1
# ==============================================================================
#
#     p = Y-X+3
#     q = Y+X-3
#
# Therefore:
#
#     2Y = p+q
#     2X = q-p+6
#
# ==============================================================================

def frame_1(
    p,
    q,
):

    if (
        (p + q) % 2 != 0
        or
        (q - p + 6) % 2 != 0
    ):
        return None

    X = (
        q - p + 6
    ) // 2

    Y = (
        p + q
    ) // 2

    if (
        Y - X + 3 != p
        or
        Y + X - 3 != q
    ):
        return None

    return X, Y


# ==============================================================================
# FRAME TYPE 2
# ==============================================================================
#
#     3p = Y+X-3
#     q  = Y-X+3
#
# Therefore:
#
#     2Y = 3p+q
#     2X = 3p-q+6
#
# ==============================================================================

def frame_2(
    p,
    q,
):

    if (
        (3 * p + q) % 2 != 0
        or
        (3 * p - q + 6) % 2 != 0
    ):
        return None

    X = (
        3 * p - q + 6
    ) // 2

    Y = (
        3 * p + q
    ) // 2

    if (
        3 * p != Y + X - 3
        or
        q != Y - X + 3
    ):
        return None

    return X, Y


# ==============================================================================
# BUILD STATE
# ==============================================================================

def build_state(
    n,
    p,
    q,
    z,
):

    modulus = 1 << z
    residue = n % modulus

    result1 = frame_1(
        p,
        q,
    )

    result2 = frame_2(
        p,
        q,
    )

    if result1 is None:
        return None

    if result2 is None:
        return None

    X1, Y1 = result1
    X2, Y2 = result2

    return State(
        n=n,
        p=p,
        q=q,
        z=z,
        residue=residue,
        X1=X1,
        Y1=Y1,
        X2=X2,
        Y2=Y2,
    )


# ==============================================================================
# BUILD LEVELS
# ==============================================================================

def build_levels(
    semiprimes,
):

    levels = {}

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        level = {}

        for n, p, q in semiprimes:

            state = build_state(
                n,
                p,
                q,
                z,
            )

            if state is not None:

                level[n] = state

        levels[z] = level

        print(
            f"    z={z:<2} "
            f"states={len(level)}"
        )

    return levels


# ==============================================================================
# RESIDUE REPORT
# ==============================================================================

def residue_report(
    levels,
):

    print()
    print("=" * 90)
    print(
        "ODD RESIDUE POPULATION"
    )
    print("=" * 90)

    for z in sorted(levels):

        counts = defaultdict(int)

        for state in levels[z].values():

            counts[
                state.residue
            ] += 1

        print()
        print(
            f"MOD {1 << z}"
        )

        for residue in sorted(
            counts
        ):

            print(
                f"    residue={residue:<5} "
                f"states={counts[residue]}"
            )


# ==============================================================================
# COORDINATE FRAME DESCRIPTION
# ==============================================================================

def frame_difference(
    state,
):

    return {
        "dX": state.X2 - state.X1,
        "dY": state.Y2 - state.Y1,
    }


# ==============================================================================
# FRAME RELATION
# ==============================================================================

def frame_relation_report(
    levels,
):

    print()
    print("=" * 90)
    print(
        "FRAME 1 -> FRAME 2 TRANSFORMATION"
    )
    print("=" * 90)

    for z in sorted(levels):

        states = levels[z]

        print()
        print(
            f"z={z}"
        )

        pairs = set()

        for state in states.values():

            pairs.add(
                (
                    state.X1,
                    state.Y1,
                    state.X2,
                    state.Y2,
                )
            )

        # Search exact affine relation:
        #
        # X2 = aX1+bY1+c
        # Y2 = dX1+eY1+f

        solutions = search_affine_2d(
            [
                (
                    s.X1,
                    s.Y1,
                    s.X2,
                    s.Y2,
                )
                for s in states.values()
            ],
            minimum_samples=5,
        )

        if solutions:

            for solution in solutions[:5]:

                print(
                    "    "
                    f"X2 = "
                    f"{solution[0]}*X1 + "
                    f"{solution[1]}*Y1 + "
                    f"{solution[2]}"
                )

                print(
                    "    "
                    f"Y2 = "
                    f"{solution[3]}*X1 + "
                    f"{solution[4]}*Y1 + "
                    f"{solution[5]}"
                )

        else:

            print(
                "    no small affine relation found"
            )


# ==============================================================================
# AFFINE 2D SEARCH
# ==============================================================================

def search_affine_2d(
    data,
    minimum_samples=1,
):

    if len(data) < minimum_samples:
        return []

    # Instead of brute-forcing 7^4*31^2 possibilities, use candidate
    # coefficient combinations and determine the constant from the first row.
    #
    # Search:
    #
    # X2 = aX1+bY1+c
    # Y2 = dX1+eY1+f
    #
    # where coefficients are in COEFFS.

    results = []

    # Use first sample for c/f.
    first = data[0]

    x1 = first[0]
    y1 = first[1]

    for a in COEFFS:

        for b in COEFFS:

            c = (
                first[2]
                - a * x1
                - b * y1
            )

            if c not in OFFSETS:
                continue

            valid = True

            for row in data:

                if (
                    a * row[0]
                    + b * row[1]
                    + c
                    != row[2]
                ):

                    valid = False
                    break

            if not valid:
                continue

            for d in COEFFS:

                for e in COEFFS:

                    f = (
                        first[3]
                        - d * x1
                        - e * y1
                    )

                    if f not in OFFSETS:
                        continue

                    valid_y = True

                    for row in data:

                        if (
                            d * row[0]
                            + e * row[1]
                            + f
                            != row[3]
                        ):

                            valid_y = False
                            break

                    if valid_y:

                        results.append(
                            (
                                a,
                                b,
                                c,
                                d,
                                e,
                                f,
                            )
                        )

    return results


# ==============================================================================
# DISCOVER RESIDUE-SPECIFIC COORDINATES
# ==============================================================================
#
# For each residue at level z:
#
#     r = n mod 2^z
#
# choose the coordinate frame that produces integer coordinates after
# division by k.
#
# Then test whether the normalized x,y values are smaller.
#
# ==============================================================================

def residue_coordinate_report(
    levels,
):

    print()
    print("=" * 90)
    print(
        "RESIDUE-SPECIFIC NORMALIZED COORDINATES"
    )
    print("=" * 90)

    for z in sorted(levels):

        k = 1 << (
            z - 1
        )

        groups = defaultdict(list)

        for state in levels[z].values():

            groups[
                state.residue
            ].append(
                state
            )

        print()
        print(
            f"LEVEL z={z}, k={k}"
        )

        for residue in sorted(
            groups
        ):

            group = groups[
                residue
            ]

            valid_1 = 0
            valid_2 = 0

            for state in group:

                if (
                    state.X1 % k == 0
                    and
                    state.Y1 % k == 0
                ):

                    valid_1 += 1

                if (
                    state.X2 % k == 0
                    and
                    state.Y2 % k == 0
                ):

                    valid_2 += 1

            print(
                f"    residue={residue:<5} "
                f"states={len(group):<7} "
                f"frame1_integral={valid_1:<7} "
                f"frame2_integral={valid_2:<7}"
            )


# ==============================================================================
# SAME-N RESIDUE TRANSITIONS
# ==============================================================================

def residue_transitions(
    levels,
):

    print()
    print("=" * 90)
    print(
        "RESIDUE -> CHILD RESIDUE"
    )
    print("=" * 90)

    transitions = defaultdict(
        lambda: defaultdict(int)
    )

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        parent = levels[z]
        child = levels[z + 1]

        for n in (
            set(parent)
            & set(child)
        ):

            parent_state = parent[n]
            child_state = child[n]

            transitions[
                (
                    z,
                    parent_state.residue,
                )
            ][
                child_state.residue
            ] += 1

    for key in sorted(
        transitions
    ):

        z, residue = key

        children = transitions[
            key
        ]

        print()
        print(
            f"z={z} "
            f"parent={residue}"
        )

        for child_residue in sorted(
            children
        ):

            print(
                f"    -> {child_residue}: "
                f"{children[child_residue]}"
            )


# ==============================================================================
# CHILD COORDINATE TRANSFORMATIONS
# ==============================================================================

def child_coordinate_transforms(
    levels,
):

    print()
    print("=" * 90)
    print(
        "RESIDUE-SPECIFIC CHILD COORDINATE TRANSFORMATIONS"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        parent = levels[z]
        child = levels[z + 1]

        grouped = defaultdict(list)

        for n in (
            set(parent)
            & set(child)
        ):

            pstate = parent[n]
            cstate = child[n]

            grouped[
                (
                    pstate.residue,
                    cstate.residue,
                )
            ].append(
                (
                    pstate,
                    cstate,
                )
            )

        print()
        print(
            f"LEVEL {z} -> {z+1}"
        )

        for transition_key in sorted(
            grouped
        ):

            parent_residue, child_residue = (
                transition_key
            )

            pairs = grouped[
                transition_key
            ]

            print()
            print(
                f"    {parent_residue} "
                f"-> {child_residue} "
                f"states={len(pairs)}"
            )

            # ------------------------------------------------------------------
            # Frame 1
            # ------------------------------------------------------------------

            data1 = [
                (
                    p.X1,
                    p.Y1,
                    c.X1,
                    c.Y1,
                )
                for p, c in pairs
            ]

            sol1 = search_affine_2d(
                data1,
                minimum_samples=min(
                    5,
                    len(data1),
                ),
            )

            if sol1:

                for a, b, c, d, e, f in sol1[:3]:

                    print(
                        "        FRAME1:"
                    )

                    print(
                        f"            "
                        f"Xc={a}*Xp+{b}*Yp+{c}"
                    )

                    print(
                        f"            "
                        f"Yc={d}*Xp+{e}*Yp+{f}"
                    )

            else:

                print(
                    "        FRAME1: "
                    "no small affine transform"
                )

            # ------------------------------------------------------------------
            # Frame 2
            # ------------------------------------------------------------------

            data2 = [
                (
                    p.X2,
                    p.Y2,
                    c.X2,
                    c.Y2,
                )
                for p, c in pairs
            ]

            sol2 = search_affine_2d(
                data2,
                minimum_samples=min(
                    5,
                    len(data2),
                ),
            )

            if sol2:

                for a, b, c, d, e, f in sol2[:3]:

                    print(
                        "        FRAME2:"
                    )

                    print(
                        f"            "
                        f"Xc={a}*Xp+{b}*Yp+{c}"
                    )

                    print(
                        f"            "
                        f"Yc={d}*Xp+{e}*Yp+{f}"
                    )

            else:

                print(
                    "        FRAME2: "
                    "no small affine transform"
                )


# ==============================================================================
# RESIDUE CHILD RULE
# ==============================================================================

def child_rule_summary(
    levels,
):

    print()
    print("=" * 90)
    print(
        "CHILD RESIDUE RULE SUMMARY"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        parent = levels[z]
        child = levels[z + 1]

        groups = defaultdict(set)

        for n in (
            set(parent)
            & set(child)
        ):

            p = parent[n]
            c = child[n]

            groups[
                p.residue
            ].add(
                c.residue
            )

        print()
        print(
            f"MOD {1 << z} -> MOD {1 << (z+1)}"
        )

        for residue in sorted(
            groups
        ):

            print(
                f"    {residue} "
                f"-> "
                f"{sorted(groups[residue])}"
            )


# ==============================================================================
# SAMPLE TRANSITIONS
# ==============================================================================

def examples(
    levels,
):

    print()
    print("=" * 90)
    print(
        "EXAMPLE RESIDUE TRANSITIONS"
    )
    print("=" * 90)

    shown = 0

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        parent = levels[z]
        child = levels[z + 1]

        for n in (
            sorted(
                set(parent)
                & set(child)
            )
        ):

            p = parent[n]
            c = child[n]

            print(
                f"    n={n} "
                f"{p.residue} -> {c.residue} "
                f""
                f"F1: "
                f"({p.X1},{p.Y1})"
                f" -> "
                f"({c.X1},{c.Y1}) "
                f""
                f"F2: "
                f"({p.X2},{p.Y2})"
                f" -> "
                f"({c.X2},{c.Y2})"
            )

            shown += 1

            if shown >= MAX_EXAMPLES:
                return


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print(
        "EXPERIMENT 590 START"
    )
    print("=" * 90)

    print()
    print(
        "DISCOVER RESIDUE-SPECIFIC FUNCTION TRANSFORMATIONS"
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
    # Build states.
    # --------------------------------------------------------------------------

    print()
    print(
        "[3] BUILD LEVELS"
    )

    levels = build_levels(
        semiprimes
    )

    # --------------------------------------------------------------------------
    # Residue populations.
    # --------------------------------------------------------------------------

    residue_report(
        levels
    )

    # --------------------------------------------------------------------------
    # Frame relationship.
    # --------------------------------------------------------------------------

    frame_relation_report(
        levels
    )

    # --------------------------------------------------------------------------
    # Residue-specific integral coordinate test.
    # --------------------------------------------------------------------------

    residue_coordinate_report(
        levels
    )

    # --------------------------------------------------------------------------
    # Parent -> child residue.
    # --------------------------------------------------------------------------

    residue_transitions(
        levels
    )

    child_rule_summary(
        levels
    )

    # --------------------------------------------------------------------------
    # Coordinate transformations.
    # --------------------------------------------------------------------------

    child_coordinate_transforms(
        levels
    )

    # --------------------------------------------------------------------------
    # Examples.
    # --------------------------------------------------------------------------

    examples(
        levels
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
Experiment 589 accidentally forced every factor pair into the
same A-coordinate frame.

Experiment 590 removes that assumption.

For every residue:

    r = n mod 2^z

we keep separate coordinate representations and explicitly measure:

    residue r
        ->
    child residue s

and:

    coordinates at z
        ->
    coordinates at z+1.

The important output is:

    RESIDUE-SPECIFIC CHILD COORDINATE TRANSFORMATIONS

A strong result would show rules such as:

    residue 3 -> residue 3:
        Xc = Xp/2
        Yc = Yp/2

    residue 3 -> residue 7:
        Xc = ...
        Yc = ...

    residue 7 -> residue 15:
        Xc = ...
        Yc = ...

etc.

That directly attacks the original hypothesis:

    A1 -> B1 or B2
    A2 -> B3 or B4

with each child having its own coordinate/function structure.

The objective is to discover those transformations rather than
assuming them.
"""
    )

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 590 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()
