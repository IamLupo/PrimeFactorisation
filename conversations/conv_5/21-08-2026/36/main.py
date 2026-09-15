#!/usr/bin/env python3

# ==============================================================================
# EXPERIMENT 594
# ==============================================================================
# FULL A/B FRAME GRAPH
#
# Keep ALL valid A/B representations for every n at every level.
#
# We explicitly test:
#
#     A -> A
#     A -> B
#     B -> A
#     B -> B
#
# for the SAME n,p,q.
#
# The main goal is to discover:
#
#     x_child = a*x_parent + b*y_parent + c
#     y_child = d*x_parent + e*y_parent + f
#
# including fractional coefficients such as 1/2, 3/2, etc.
#
# We also derive the symbolic transformation directly between the
# two coordinate frames.
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

MAX_EXAMPLES = 10


# ==============================================================================
# DATA
# ==============================================================================

@dataclass(frozen=True)
class Representation:

    n: int
    p: int
    q: int

    z: int
    residue: int

    frame: str

    x: Fraction
    y: Fraction


# ==============================================================================
# PRIME SIEVE
# ==============================================================================

def prime_sieve(limit: int):

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
# FRAME A
# ==============================================================================
#
#     p = k(y-x)+3
#     q = k(y+x)-3
#
# therefore:
#
#     x = (q-p+6)/(2k)
#     y = (p+q)/(2k)
#
# ==============================================================================

def frame_A(
    p,
    q,
    z,
):

    k = 1 << (
        z - 1
    )

    denominator = 2 * k

    x_num = (
        q
        - p
        + 6
    )

    y_num = (
        p
        + q
    )

    if (
        x_num % denominator != 0
        or
        y_num % denominator != 0
    ):
        return None

    x = Fraction(
        x_num,
        denominator,
    )

    y = Fraction(
        y_num,
        denominator,
    )

    return x, y


# ==============================================================================
# FRAME B
# ==============================================================================
#
#     3p = k(y+x)-3
#     q  = k(y-x)+3
#
# therefore:
#
#     x = (3p-q+6)/(2k)
#     y = (3p+q)/(2k)
#
# ==============================================================================

def frame_B(
    p,
    q,
    z,
):

    k = 1 << (
        z - 1
    )

    denominator = 2 * k

    x_num = (
        3 * p
        - q
        + 6
    )

    y_num = (
        3 * p
        + q
    )

    if (
        x_num % denominator != 0
        or
        y_num % denominator != 0
    ):
        return None

    x = Fraction(
        x_num,
        denominator,
    )

    y = Fraction(
        y_num,
        denominator,
    )

    return x, y


# ==============================================================================
# BUILD ALL REPRESENTATIONS
# ==============================================================================

def build_representations(
    semiprimes,
):

    levels = {}

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        table = defaultdict(list)

        modulus = 1 << z

        for n, p, q in semiprimes:

            xy_A = frame_A(
                p,
                q,
                z,
            )

            if xy_A is not None:

                x, y = xy_A

                table[n].append(
                    Representation(
                        n=n,
                        p=p,
                        q=q,
                        z=z,
                        residue=n % modulus,
                        frame="A",
                        x=x,
                        y=y,
                    )
                )

            xy_B = frame_B(
                p,
                q,
                z,
            )

            if xy_B is not None:

                x, y = xy_B

                table[n].append(
                    Representation(
                        n=n,
                        p=p,
                        q=q,
                        z=z,
                        residue=n % modulus,
                        frame="B",
                        x=x,
                        y=y,
                    )
                )

        levels[z] = table

        representation_count = sum(
            len(v)
            for v in table.values()
        )

        print(
            f"    z={z:<2} "
            f"n={len(table):<7} "
            f"representations="
            f"{representation_count}"
        )

    return levels


# ==============================================================================
# AFFINE FIT WITH FRACTIONAL COEFFICIENTS
# ==============================================================================
#
# Find:
#
#     target = a*x + b*y + c
#
# directly from exact rational data.
#
# We use three non-collinear points.
#
# ==============================================================================

def solve_3x3(
    matrix,
    vector,
):

    a = [
        [
            Fraction(v)
            for v in row
        ]
        for row in matrix
    ]

    b = [
        Fraction(v)
        for v in vector
    ]

    n = 3

    for col in range(n):

        pivot = None

        for row in range(
            col,
            n,
        ):

            if a[row][col] != 0:

                pivot = row

                break

        if pivot is None:
            return None

        if pivot != col:

            a[col], a[pivot] = (
                a[pivot],
                a[col],
            )

            b[col], b[pivot] = (
                b[pivot],
                b[col],
            )

        factor = a[col][col]

        for j in range(
            col,
            n,
        ):

            a[col][j] /= factor

        b[col] /= factor

        for row in range(n):

            if row == col:
                continue

            factor = a[row][col]

            for j in range(
                col,
                n,
            ):

                a[row][j] -= (
                    factor
                    * a[col][j]
                )

            b[row] -= (
                factor
                * b[col]
            )

    return (
        b[0],
        b[1],
        b[2],
    )


def fit_affine(
    rows,
):

    if len(rows) < 3:
        return None

    # Find 3 non-collinear points.

    base = None

    for i in range(
        len(rows)
    ):

        for j in range(
            i + 1,
            len(rows),
        ):

            for k in range(
                j + 1,
                len(rows),
            ):

                x1, y1, _ = rows[i]
                x2, y2, _ = rows[j]
                x3, y3, _ = rows[k]

                determinant = (
                    x1 * (y2 - y3)
                    +
                    x2 * (y3 - y1)
                    +
                    x3 * (y1 - y2)
                )

                if determinant != 0:

                    base = (
                        rows[i],
                        rows[j],
                        rows[k],
                    )

                    break

            if base:
                break

        if base:
            break

    if base is None:
        return None

    matrix = [
        [
            base[0][0],
            base[0][1],
            Fraction(1),
        ],
        [
            base[1][0],
            base[1][1],
            Fraction(1),
        ],
        [
            base[2][0],
            base[2][1],
            Fraction(1),
        ],
    ]

    vector = [
        base[0][2],
        base[1][2],
        base[2][2],
    ]

    rule = solve_3x3(
        matrix,
        vector,
    )

    if rule is None:
        return None

    for x, y, target in rows:

        predicted = (
            rule[0] * x
            +
            rule[1] * y
            +
            rule[2]
        )

        if predicted != target:

            return None

    return rule


# ==============================================================================
# RULE TEXT
# ==============================================================================

def rule_text(
    name,
    rule,
):

    if rule is None:

        return (
            f"{name}=NONE"
        )

    a, b, c = rule

    return (
        f"{name}="
        f"({a})x + "
        f"({b})y + "
        f"({c})"
    )


# ==============================================================================
# SAME-N TRANSITIONS
# ==============================================================================
#
# IMPORTANT:
#
# Every representation is retained.
#
# Therefore a transition may contain:
#
#     A -> A
#     A -> B
#     B -> A
#     B -> B
#
# ==============================================================================

def build_transitions(
    levels,
):

    transitions = defaultdict(list)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        common_n = (
            set(levels[z])
            &
            set(levels[z + 1])
        )

        for n in common_n:

            parents = levels[z][n]
            children = levels[z + 1][n]

            for parent in parents:

                for child in children:

                    key = (
                        z,
                        parent.residue,
                        parent.frame,
                        child.residue,
                        child.frame,
                    )

                    transitions[key].append(
                        (
                            parent,
                            child,
                        )
                    )

    return transitions


# ==============================================================================
# TRANSFORMATION DISCOVERY
# ==============================================================================

def discover_transform(
    pairs,
):

    x_rows = [
        (
            parent.x,
            parent.y,
            child.x,
        )
        for parent, child in pairs
    ]

    y_rows = [
        (
            parent.x,
            parent.y,
            child.y,
        )
        for parent, child in pairs
    ]

    x_rule = fit_affine(
        x_rows
    )

    y_rule = fit_affine(
        y_rows
    )

    return (
        x_rule,
        y_rule,
    )


# ==============================================================================
# TRANSITION REPORT
# ==============================================================================

def report_transitions(
    transitions,
):

    print()
    print("=" * 90)
    print(
        "FULL FRAME TRANSITION GRAPH"
    )
    print("=" * 90)

    for key in sorted(
        transitions
    ):

        (
            z,
            parent_residue,
            parent_frame,
            child_residue,
            child_frame,
        ) = key

        pairs = transitions[key]

        x_rule, y_rule = (
            discover_transform(
                pairs
            )
        )

        exact = (
            x_rule is not None
            and
            y_rule is not None
        )

        print()
        print(
            f"z={z}"
            f"  "
            f"{parent_residue}{parent_frame}"
            f" -> "
            f"{child_residue}{child_frame}"
        )

        print(
            f"    states={len(pairs)}"
        )

        print(
            f"    exact={exact}"
        )

        print(
            "    "
            + rule_text(
                "x_child",
                x_rule,
            )
        )

        print(
            "    "
            + rule_text(
                "y_child",
                y_rule,
            )
        )


# ==============================================================================
# TRANSITION CLASSIFICATION
# ==============================================================================

def classify(
    x_rule,
    y_rule,
):

    if (
        x_rule is None
        or
        y_rule is None
    ):
        return "NONE"

    if (
        x_rule
        == (
            Fraction(1, 2),
            Fraction(0),
            Fraction(0),
        )
        and
        y_rule
        == (
            Fraction(0),
            Fraction(1, 2),
            Fraction(0),
        )
    ):
        return "HALF"

    if (
        x_rule
        == (
            Fraction(-1, 2),
            Fraction(0),
            Fraction(0),
        )
        and
        y_rule
        == (
            Fraction(0),
            Fraction(-1, 2),
            Fraction(0),
        )
    ):
        return "NEG_HALF"

    return "GENERAL"


# ==============================================================================
# COMPACT TRANSITION SUMMARY
# ==============================================================================

def compact(
    transitions,
):

    print()
    print("=" * 90)
    print(
        "COMPACT A/B TRANSITION SUMMARY"
    )
    print("=" * 90)

    for key in sorted(
        transitions
    ):

        z, pr, pf, cr, cf = key

        x_rule, y_rule = (
            discover_transform(
                transitions[key]
            )
        )

        print(
            f"    z={z} "
            f"{pr}{pf}->{cr}{cf} "
            f"{classify(x_rule,y_rule)}"
        )


# ==============================================================================
# SYMBOLIC FRAME CHANGES
# ==============================================================================
#
# Derive B coordinates in terms of A coordinates.
#
# Starting from:
#
#     A:
#         p = k(y-x)+3
#         q = k(y+x)-3
#
# Substitute into:
#
#     B:
#         x_B=(3p-q+6)/(2k)
#         y_B=(3p+q)/(2k)
#
# ==============================================================================

def symbolic_frame_transform():

    print()
    print("=" * 90)
    print(
        "SYMBOLIC A <-> B FRAME TRANSFORMATION"
    )
    print("=" * 90)

    print()
    print(
        "FRAME A:"
    )

    print(
        "    p = k(y_A-x_A)+3"
    )

    print(
        "    q = k(y_A+x_A)-3"
    )

    print()
    print(
        "Substituting into FRAME B:"
    )

    print(
        "    x_B"
        " = (3p-q+6)/(2k)"
    )

    print(
        "    y_B"
        " = (3p+q)/(2k)"
    )

    print()
    print(
        "RESULT:"
    )

    print(
        "    x_B = "
        "(y_A - 2*x_A + 9/k) / 2"
    )

    print(
        "    y_B = "
        "(-x_A + 2*y_A + 3/k) / 2"
    )

    print()
    print(
        "So the A -> B coordinate map is:"
    )

    print(
        "    x_B = -x_A + y_A/2 + 9/(2k)"
    )

    print(
        "    y_B = -x_A/2 + y_A + 3/(2k)"
    )

    print()
    print(
        "This is a completely different transformation"
    )

    print(
        "from the level-renormalization map:"
    )

    print(
        "    x_(z+1) = x_z/2"
    )

    print(
        "    y_(z+1) = y_z/2"
    )

    print()
    print(
        "Therefore the experiment separates two operations:"
    )

    print(
        "    1. FRAME CHANGE: A <-> B"
    )

    print(
        "    2. LEVEL RENORMALIZATION: (x,y)->(x/2,y/2)"
    )


# ==============================================================================
# FRAME CHANGE TEST
# ==============================================================================

def test_frame_change(
    transitions,
):

    print()
    print("=" * 90)
    print(
        "EMPIRICAL FRAME-CHANGE TEST"
    )
    print("=" * 90)

    for key in sorted(
        transitions
    ):

        z, pr, pf, cr, cf = key

        if pf == cf:
            continue

        pairs = transitions[key]

        print()
        print(
            f"z={z} "
            f"{pf}->{cf} "
            f"residue "
            f"{pr}->{cr}"
        )

        # Expected symbolic map from A to B
        # and inverse B to A.
        #
        # We test the exact formulas empirically.

        failures = 0

        k = Fraction(
            1 << (
                z - 1
            )
        )

        for parent, child in pairs:

            if (
                pf == "A"
                and cf == "B"
            ):

                expected_x = (
                    -parent.x
                    +
                    parent.y / 2
                    +
                    Fraction(9, 2) / k
                )

                expected_y = (
                    -parent.x / 2
                    +
                    parent.y
                    +
                    Fraction(3, 2) / k
                )

            elif (
                pf == "B"
                and cf == "A"
            ):

                # Derive inverse numerically through solving
                # the same linear system.

                expected_x = (
                    -parent.x
                    +
                    parent.y / 2
                    -
                    Fraction(9, 2) / k
                )

                expected_y = (
                    -parent.x / 2
                    +
                    parent.y
                    -
                    Fraction(3, 2) / k
                )

            else:

                continue

            if (
                child.x
                != expected_x
                or
                child.y
                != expected_y
            ):

                failures += 1

        print(
            f"    states={len(pairs)}"
        )

        print(
            f"    failures={failures}"
        )


# ==============================================================================
# RESIDUE CHILD RULE
# ==============================================================================

def child_residue_rule(
    levels,
):

    print()
    print("=" * 90)
    print(
        "CHILD RESIDUE FROM PARENT FRAME"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        table = defaultdict(set)

        common_n = (
            set(levels[z])
            &
            set(levels[z + 1])
        )

        for n in common_n:

            for parent in levels[z][n]:

                for child in levels[z + 1][n]:

                    key = (
                        parent.frame,
                        parent.residue,
                        parent.x % 2,
                        parent.y % 2,
                    )

                    table[key].add(
                        (
                            child.frame,
                            child.residue,
                        )
                    )

        print()
        print(
            f"z={z}"
        )

        for key in sorted(
            table
        ):

            print(
                f"    {key} "
                f"-> "
                f"{sorted(table[key])}"
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
        "A/B TRANSITION EXAMPLES"
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
            f"{key[2]}->{key[4]} "
            f"{key[1]}->{key[3]}"
        )

        for parent, child in pairs[
            :MAX_EXAMPLES
        ]:

            print(
                f"    n={parent.n} "
                f""
                f"({parent.x},{parent.y})"
                f" -> "
                f"({child.x},{child.y})"
            )

        shown += 1

        if shown >= 12:
            break


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print(
        "EXPERIMENT 594 START"
    )
    print("=" * 90)

    print()
    print(
        "FULL A/B FRAME GRAPH "
        "WITH FRACTIONAL TRANSFORMATION SEARCH"
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

    print()
    print(
        "[2] SEMIPRIME GENERATION"
    )

    semiprimes = generate_semiprimes(
        MAX_N,
        sieve,
    )

    print(
        f"    semiprimes="
        f"{len(semiprimes)}"
    )

    # --------------------------------------------------------------------------
    # Build representations.
    # --------------------------------------------------------------------------

    print()
    print(
        "[3] BUILD ALL A/B REPRESENTATIONS"
    )

    levels = build_representations(
        semiprimes
    )

    # --------------------------------------------------------------------------
    # Symbolic frame map.
    # --------------------------------------------------------------------------

    symbolic_frame_transform()

    # --------------------------------------------------------------------------
    # Transitions.
    # --------------------------------------------------------------------------

    print()
    print(
        "[4] BUILD FULL SAME-N FRAME GRAPH"
    )

    transitions = build_transitions(
        levels
    )

    # --------------------------------------------------------------------------
    # Reports.
    # --------------------------------------------------------------------------

    report_transitions(
        transitions
    )

    compact(
        transitions
    )

    test_frame_change(
        transitions
    )

    child_residue_rule(
        levels
    )

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
Experiment 594 keeps BOTH coordinate frames.

For every level:

    FRAME A:
        p = k(y-x)+3
        q = k(y+x)-3

    FRAME B:
        3p = k(y+x)-3
        q  = k(y-x)+3

A factor pair can therefore be represented by whichever frame
is compatible with the current residue.

The experiment follows the SAME n across two consecutive levels
and explicitly measures:

    A -> A
    A -> B
    B -> A
    B -> B

The transformation search now permits exact rational coefficients.

This is important because there are two fundamentally different
transformations:

    LEVEL RENORMALIZATION:
        (x,y) -> (x/2,y/2)

    FRAME CHANGE:
        A -> B
        or
        B -> A

The original hierarchical idea requires both pieces:

    parent residue
        ->
    branch/frame selection
        ->
    coordinate transformation
        ->
    child function.

The output to look for is therefore:

    A -> B:
        exact transformation
        exact residue transition

and whether the same transformation repeats at later levels.
"""
    )

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 594 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()
