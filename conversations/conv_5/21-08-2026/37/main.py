#!/usr/bin/env python3

# ==============================================================================
# EXPERIMENT 595
# ==============================================================================
# RESIDUE REFINEMENT TREE + FRAME ASSIGNMENT
#
# Goal
# ----
#
# Separate three different structures:
#
#   1. RESIDUE REFINEMENT
#
#          r_z = n mod 2^z
#
#          r_z -> r_{z+1}
#
#   2. FRAME ASSIGNMENT
#
#          residue -> A or B
#
#   3. COORDINATE RENORMALIZATION
#
#          x_{z+1} = x_z / 2
#          y_{z+1} = y_z / 2
#
# We do NOT assume that a child residue requires a frame change.
#
# Instead we discover the complete tree:
#
#       residue + frame
#               |
#               +---- child residue + frame
#
# and derive the function attached to every residue/frame node.
#
# No files are read or written.
# No external sources are used.
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
MAX_Z = 10

MAX_EXAMPLES = 8


# ==============================================================================
# DATA STRUCTURE
# ==============================================================================

@dataclass(frozen=True)
class State:

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
    ) * (limit + 1)

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
# FRAME A
# ==============================================================================
#
#     p = k(y-x)+3
#     q = k(y+x)-3
#
# Therefore:
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

    d = 2 * k

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
        x_num % d != 0
        or
        y_num % d != 0
    ):
        return None

    return (
        Fraction(
            x_num,
            d,
        ),
        Fraction(
            y_num,
            d,
        ),
    )


# ==============================================================================
# FRAME B
# ==============================================================================
#
#     3p = k(y+x)-3
#     q  = k(y-x)+3
#
# Therefore:
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

    d = 2 * k

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
        x_num % d != 0
        or
        y_num % d != 0
    ):
        return None

    return (
        Fraction(
            x_num,
            d,
        ),
        Fraction(
            y_num,
            d,
        ),
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

        table = defaultdict(list)

        modulus = 1 << z

        for n, p, q in semiprimes:

            xy = frame_A(
                p,
                q,
                z,
            )

            if xy is not None:

                x, y = xy

                table[
                    (
                        n,
                        "A",
                    )
                ].append(
                    State(
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

            xy = frame_B(
                p,
                q,
                z,
            )

            if xy is not None:

                x, y = xy

                table[
                    (
                        n,
                        "B",
                    )
                ].append(
                    State(
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

        states = sum(
            len(values)
            for values in table.values()
        )

        print(
            f"    z={z:<2} "
            f"states={states}"
        )

    return levels


# ==============================================================================
# NODE INDEX
# ==============================================================================

def node_index(
    levels,
):

    nodes = defaultdict(list)

    for z, table in levels.items():

        for values in table.values():

            for state in values:

                key = (
                    z,
                    state.residue,
                    state.frame,
                )

                nodes[key].append(
                    state
                )

    return nodes


# ==============================================================================
# RESIDUE -> FRAME ASSIGNMENT
# ==============================================================================

def report_frame_assignment(
    nodes,
):

    print()
    print("=" * 90)
    print(
        "RESIDUE -> FRAME ASSIGNMENT"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        print()
        print(
            f"MOD {1 << z}"
        )

        modulus = 1 << z

        for residue in range(
            1,
            modulus,
            2,
        ):

            A = nodes.get(
                (
                    z,
                    residue,
                    "A",
                ),
                [],
            )

            B = nodes.get(
                (
                    z,
                    residue,
                    "B",
                ),
                [],
            )

            if not A and not B:
                continue

            total = (
                len(A)
                +
                len(B)
            )

            if A and B:

                assignment = "A+B"

            elif A:

                assignment = "A"

            else:

                assignment = "B"

            print(
                f"    residue={residue:<5} "
                f"frame={assignment:<3} "
                f"states={total}"
            )


# ==============================================================================
# SAME-N RESIDUE TRANSITIONS
# ==============================================================================

def build_node_transitions(
    levels,
):

    transitions = defaultdict(list)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        current = levels[z]
        child = levels[z + 1]

        current_n = {
            n
            for n, frame in current
        }

        child_n = {
            n
            for n, frame in child
        }

        common = (
            current_n
            &
            child_n
        )

        for n in common:

            parents = []

            children = []

            for (
                key,
                values,
            ) in current.items():

                if key[0] == n:

                    parents.extend(
                        values
                    )

            for (
                key,
                values,
            ) in child.items():

                if key[0] == n:

                    children.extend(
                        values
                    )

            for parent in parents:

                for child_state in children:

                    key = (
                        z,
                        parent.residue,
                        parent.frame,
                        child_state.residue,
                        child_state.frame,
                    )

                    transitions[
                        key
                    ].append(
                        (
                            parent,
                            child_state,
                        )
                    )

    return transitions


# ==============================================================================
# AFFINE FIT
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

    for col in range(3):

        pivot = None

        for row in range(
            col,
            3,
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
            3,
        ):

            a[col][j] /= factor

        b[col] /= factor

        for row in range(3):

            if row == col:
                continue

            factor = a[row][col]

            for j in range(
                col,
                3,
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


def affine_fit(
    rows,
):

    if len(rows) < 3:
        return None

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

                p1 = rows[i]
                p2 = rows[j]
                p3 = rows[k]

                determinant = (
                    p1[0]
                    * (
                        p2[1]
                        -
                        p3[1]
                    )
                    +
                    p2[0]
                    * (
                        p3[1]
                        -
                        p1[1]
                    )
                    +
                    p3[0]
                    * (
                        p1[1]
                        -
                        p2[1]
                    )
                )

                if determinant != 0:

                    base = (
                        p1,
                        p2,
                        p3,
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
# RESIDUE REFINEMENT GRAPH
# ==============================================================================

def report_residue_tree(
    transitions,
):

    print()
    print("=" * 90)
    print(
        "RESIDUE REFINEMENT TREE"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        print()
        print(
            f"LEVEL {z}: "
            f"mod {1 << z}"
            f" -> "
            f"mod {1 << (z + 1)}"
        )

        rows = []

        for key in sorted(
            transitions
        ):

            (
                level,
                parent_residue,
                parent_frame,
                child_residue,
                child_frame,
            ) = key

            if level != z:
                continue

            rows.append(
                key
            )

        for key in rows:

            (
                _,
                parent_residue,
                parent_frame,
                child_residue,
                child_frame,
            ) = key

            states = len(
                transitions[key]
            )

            print(
                f"    "
                f"{parent_residue}"
                f"{parent_frame}"
                f" -> "
                f"{child_residue}"
                f"{child_frame}"
                f" "
                f"states={states}"
            )


# ==============================================================================
# DISCARD CURRENT x/y BIT
# ==============================================================================
#
# Since:
#
#     x_parent = 2*x_child
#     y_parent = 2*y_child
#
# in surviving transitions, the discarded low bits are:
#
#     x_parent mod 2
#     y_parent mod 2
#
# ==============================================================================

def report_discarded_bits(
    transitions,
):

    print()
    print("=" * 90)
    print(
        "DISCARDED LOW BITS -> CHILD RESIDUE"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        table = defaultdict(set)

        for key, pairs in transitions.items():

            level = key[0]

            if level != z:
                continue

            for parent, child in pairs:

                state_key = (
                    parent.frame,
                    parent.residue,
                    int(parent.x) & 1,
                    int(parent.y) & 1,
                )

                table[
                    state_key
                ].add(
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

            values = table[key]

            print(
                f"    "
                f"{key}"
                f" -> "
                f"{sorted(values)}"
            )


# ==============================================================================
# CHILD BIT RULE
# ==============================================================================
#
# r_{z+1} = r_z + b_z * 2^z
#
# therefore the child branch is exactly the next bit.
#
# We test whether this bit is determined by the parent frame and
# discarded coordinate bits.
#
# ==============================================================================

def report_child_bits(
    transitions,
):

    print()
    print("=" * 90)
    print(
        "CHILD BIT RULE"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        table = defaultdict(set)

        for key, pairs in transitions.items():

            level = key[0]

            if level != z:
                continue

            for parent, child in pairs:

                bit = (
                    (
                        child.residue
                        -
                        parent.residue
                    )
                    >>
                    z
                ) & 1

                feature = (
                    parent.frame,
                    parent.residue,
                    int(parent.x) & 1,
                    int(parent.y) & 1,
                )

                table[
                    feature
                ].add(
                    bit
                )

        deterministic = all(
            len(values) == 1
            for values in table.values()
        )

        print()
        print(
            f"z={z} "
            f"states={len(table)} "
            f"deterministic="
            f"{deterministic}"
        )

        for feature in sorted(
            table
        ):

            print(
                f"    {feature}"
                f" -> "
                f"{sorted(table[feature])}"
            )


# ==============================================================================
# FUNCTION COEFFICIENTS
# ==============================================================================

def frame_function(
    frame,
    z,
):

    k = Fraction(
        1 << (
            z - 1
        )
    )

    if frame == "A":

        # p = k(y-x)+3
        # q = k(y+x)-3

        return (
            (-k, k, Fraction(3)),
            (k, k, Fraction(-3)),
        )

    # B:
    #
    # 3p = k(y+x)-3
    # q  = k(y-x)+3

    return (
        (
            -k / 3,
            k / 3,
            Fraction(-1),
        ),
        (
            -k,
            k,
            Fraction(3),
        ),
    )


def report_functions(
    nodes,
):

    print()
    print("=" * 90)
    print(
        "FUNCTION ATTACHED TO EACH RESIDUE"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        print()
        print(
            f"LEVEL z={z} "
            f"mod={1 << z}"
        )

        for residue in range(
            1,
            1 << z,
            2,
        ):

            found = []

            for frame in (
                "A",
                "B",
            ):

                states = nodes.get(
                    (
                        z,
                        residue,
                        frame,
                    ),
                    [],
                )

                if states:

                    found.append(
                        frame
                    )

            if not found:
                continue

            for frame in found:

                p_rule, q_rule = (
                    frame_function(
                        frame,
                        z,
                    )
                )

                print(
                    f"    residue={residue}"
                    f" frame={frame}"
                )

                print(
                    f"        p="
                    f"({p_rule[0]})x + "
                    f"({p_rule[1]})y + "
                    f"({p_rule[2]})"
                )

                print(
                    f"        q="
                    f"({q_rule[0]})x + "
                    f"({q_rule[1]})y + "
                    f"({q_rule[2]})"
                )


# ==============================================================================
# LEVEL RENORMALIZATION TEST
# ==============================================================================

def report_renormalization(
    transitions,
):

    print()
    print("=" * 90)
    print(
        "LEVEL RENORMALIZATION"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        summary = defaultdict(
            lambda: [
                0,
                0,
                0,
            ]
        )

        for key, pairs in transitions.items():

            if key[0] != z:
                continue

            (
                _,
                pr,
                pf,
                cr,
                cf,
            ) = key

            counter = summary[
                (
                    pr,
                    pf,
                    cr,
                    cf,
                )
            ]

            for parent, child in pairs:

                counter[0] += 1

                if (
                    child.x
                    ==
                    parent.x / 2
                ):

                    counter[1] += 1

                if (
                    child.y
                    ==
                    parent.y / 2
                ):

                    counter[2] += 1

        print()
        print(
            f"z={z}"
        )

        for key in sorted(
            summary
        ):

            total, x_ok, y_ok = (
                summary[key]
            )

            print(
                f"    "
                f"{key}"
                f" states={total}"
                f" x_half="
                f"{x_ok}/{total}"
                f" y_half="
                f"{y_ok}/{total}"
            )


# ==============================================================================
# FRAME CHANGE SEARCH
# ==============================================================================

def report_frame_changes(
    transitions,
):

    print()
    print("=" * 90)
    print(
        "FRAME CHANGES"
    )
    print("=" * 90)

    found = False

    for key in sorted(
        transitions
    ):

        z, pr, pf, cr, cf = key

        if pf == cf:
            continue

        found = True

        print()
        print(
            f"z={z} "
            f"{pr}{pf}"
            f" -> "
            f"{cr}{cf}"
        )

        print(
            f"    states="
            f"{len(transitions[key])}"
        )

    if not found:

        print(
            "    NO A<->B transitions "
            "observed among same-n states."
        )


# ==============================================================================
# BRANCH-DEPENDENT RECURSION
# ==============================================================================

def report_recursive_nodes(
    transitions,
):

    print()
    print("=" * 90)
    print(
        "RECURSIVE NODE MAP"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        print()
        print(
            f"z={z}"
        )

        parent_nodes = set()

        for key in transitions:

            if key[0] != z:
                continue

            parent_nodes.add(
                (
                    key[1],
                    key[2],
                )
            )

        for parent in sorted(
            parent_nodes
        ):

            children = []

            for key in sorted(
                transitions
            ):

                if key[0] != z:
                    continue

                if (
                    key[1],
                    key[2],
                ) != parent:
                    continue

                children.append(
                    (
                        key[3],
                        key[4],
                    )
                )

            print(
                f"    "
                f"{parent[0]}{parent[1]}"
                f" -> "
                f"{children}"
            )


# ==============================================================================
# EXAMPLES
# ==============================================================================

def report_examples(
    transitions,
):

    print()
    print("=" * 90)
    print(
        "EXAMPLE RECURSIVE TRANSITIONS"
    )
    print("=" * 90)

    count = 0

    for key in sorted(
        transitions
    ):

        print()
        print(
            f"z={key[0]} "
            f"{key[1]}{key[2]}"
            f" -> "
            f"{key[3]}{key[4]}"
        )

        for parent, child in transitions[key][
            :MAX_EXAMPLES
        ]:

            print(
                f"    "
                f"n={parent.n}"
                f" "
                f"({parent.x},"
                f"{parent.y})"
                f" -> "
                f"({child.x},"
                f"{child.y})"
            )

        count += 1

        if count >= 12:
            break


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print(
        "EXPERIMENT 595 START"
    )
    print("=" * 90)

    print()
    print(
        "RESIDUE REFINEMENT TREE + "
        "A/B FRAME ASSIGNMENT"
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
    # Build.
    # --------------------------------------------------------------------------

    print()
    print(
        "[3] BUILD ALL A/B LEVEL STATES"
    )

    levels = build_levels(
        semiprimes
    )

    # --------------------------------------------------------------------------
    # Nodes.
    # --------------------------------------------------------------------------

    nodes = node_index(
        levels
    )

    # --------------------------------------------------------------------------
    # Frame assignment.
    # --------------------------------------------------------------------------

    report_frame_assignment(
        nodes
    )

    # --------------------------------------------------------------------------
    # Same-n transitions.
    # --------------------------------------------------------------------------

    print()
    print(
        "[4] BUILD SAME-N RESIDUE TRANSITIONS"
    )

    transitions = (
        build_node_transitions(
            levels
        )
    )

    # --------------------------------------------------------------------------
    # Residue tree.
    # --------------------------------------------------------------------------

    report_residue_tree(
        transitions
    )

    # --------------------------------------------------------------------------
    # Function table.
    # --------------------------------------------------------------------------

    report_functions(
        nodes
    )

    # --------------------------------------------------------------------------
    # Renormalization.
    # --------------------------------------------------------------------------

    report_renormalization(
        transitions
    )

    # --------------------------------------------------------------------------
    # Low bits.
    # --------------------------------------------------------------------------

    report_discarded_bits(
        transitions
    )

    report_child_bits(
        transitions
    )

    # --------------------------------------------------------------------------
    # Frame changes.
    # --------------------------------------------------------------------------

    report_frame_changes(
        transitions
    )

    # --------------------------------------------------------------------------
    # Recursive node map.
    # --------------------------------------------------------------------------

    report_recursive_nodes(
        transitions
    )

    # --------------------------------------------------------------------------
    # Examples.
    # --------------------------------------------------------------------------

    report_examples(
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
The previous experiment showed that the dominant transitions are:

    A -> A
    B -> B

with:

    x_child = x_parent / 2
    y_child = y_parent / 2.

Experiment 595 therefore treats RESIDUE and FRAME as separate
state variables.

The object being reconstructed is:

    (z, residue, frame)
            |
            v
    (z+1, child_residue, child_frame)

while the coordinate transformation is tracked independently.

The desired hierarchy is now:

    residue/frame node
            |
            +-- child residue 0
            |
            +-- child residue 1

and each node owns its own function:

    p = A_i(x,y)
    q = A_i(x,y)

The key question is no longer:

    "Does A transform into B?"

but:

    "For a given parent frame/residue,
     which child residues occur,
     which frame do they use,
     and what static coordinate transformation connects them?"

If every node has a deterministic child-residue rule and the
function coefficients follow a fixed recurrence, that gives the
hierarchical A -> B -> C structure directly.

No external data or files are used.
"""
    )

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 595 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()
