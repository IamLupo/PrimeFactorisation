#!/usr/bin/env python3

# ==============================================================================
# EXPERIMENT 586
# ==============================================================================
# SAME-N RECURSIVE COORDINATE TRANSFORMATION
#
# NO FILES
# NO WEB
#
# We now track the SAME factor pair (n,p,q) through consecutive levels:
#
#     z -> z+1
#
# Each level has:
#
#     k = 2^(z-1)
#
# and:
#
# CASE A:
#
#     n mod 2^z = 2^z - 1
#
#     p = k(y-x) + 3
#     q = k(y+x) - 3
#
# CASE B:
#
#     otherwise
#
#     3p = k(y+x) - 3
#     q  = k(y-x) + 3
#
#
# IMPORTANT:
#
# A sample is only considered in the recursive transition if the SAME n,p,q
# has valid coordinates at BOTH z and z+1.
#
# This fixes the problem in Experiment 585 where x//2 was tested against
# samples that did not survive to the child level.
#
#
# MAIN QUESTIONS
# --------------
#
# 1. How does x change?
#
#       x_child = ?
#
# 2. How does y change?
#
#       y_child = ?
#
# 3. Does the transformation depend on the parent branch?
#
# 4. Does it depend on the child branch?
#
# 5. Does the transformation repeat at every level?
#
#
# Candidate transformations:
#
#     x/2
#     (x-1)/2
#     (x+1)/2
#     y/2
#     (y-1)/2
#     (y+1)/2
#
# More general small affine rules are also searched:
#
#     child = a*parent + b
#
# with:
#
#     a in {-2,-1,-1/2,1/2,1,2}
#
# and small integer b.
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

SHOW_EXAMPLES = True
MAX_EXAMPLES = 12


# ==============================================================================
# DATA STRUCTURE
# ==============================================================================

@dataclass(frozen=True)
class LevelState:

    n: int
    p: int
    q: int

    z: int
    k: int
    modulus: int

    branch: str

    x: int
    y: int

    u: int
    v: int


# ==============================================================================
# PRIME SIEVE
# ==============================================================================

def prime_sieve(
    limit: int,
):

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
    max_n: int,
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
# CASE A SOLVER
# ==============================================================================

def solve_case_a(
    p: int,
    q: int,
    k: int,
):

    # p = k(y-x)+3
    # q = k(y+x)-3

    if (p - 3) % k != 0:
        return None

    if (q + 3) % k != 0:
        return None

    u = (
        p - 3
    ) // k

    v = (
        q + 3
    ) // k

    if (v - u) % 2 != 0:
        return None

    if (v + u) % 2 != 0:
        return None

    x = (
        v - u
    ) // 2

    y = (
        v + u
    ) // 2

    return x, y, u, v


# ==============================================================================
# CASE B SOLVER
# ==============================================================================

def solve_case_b(
    p: int,
    q: int,
    k: int,
):

    # 3p = k(y+x)-3
    # q  = k(y-x)+3

    if (3 * p + 3) % k != 0:
        return None

    if (q - 3) % k != 0:
        return None

    v = (
        3 * p + 3
    ) // k

    u = (
        q - 3
    ) // k

    if (v - u) % 2 != 0:
        return None

    if (v + u) % 2 != 0:
        return None

    x = (
        v - u
    ) // 2

    y = (
        v + u
    ) // 2

    return x, y, u, v


# ==============================================================================
# CONSTRUCT ONE LEVEL
# ==============================================================================

def construct_level_state(
    n: int,
    p: int,
    q: int,
    z: int,
):

    k = 1 << (
        z - 1
    )

    modulus = 1 << z

    residue = n % modulus

    # --------------------------------------------------------------------------
    # Case A
    # --------------------------------------------------------------------------

    if residue == modulus - 1:

        solved = solve_case_a(
            p,
            q,
            k,
        )

        if solved is None:
            return None

        x, y, u, v = solved

        branch = "A"

        check_p = (
            k * u + 3
        )

        check_q = (
            k * v - 3
        )

    # --------------------------------------------------------------------------
    # Case B
    # --------------------------------------------------------------------------

    else:

        solved = solve_case_b(
            p,
            q,
            k,
        )

        if solved is None:
            return None

        x, y, u, v = solved

        branch = "B"

        numerator = (
            k * v - 3
        )

        if numerator % 3 != 0:
            return None

        check_p = (
            numerator // 3
        )

        check_q = (
            k * u + 3
        )

    # --------------------------------------------------------------------------
    # Exact verification.
    # --------------------------------------------------------------------------

    if check_p != p:
        raise RuntimeError(
            "p reconstruction failed"
        )

    if check_q != q:
        raise RuntimeError(
            "q reconstruction failed"
        )

    if p * q != n:
        raise RuntimeError(
            "n reconstruction failed"
        )

    if v - u != 2 * x:
        raise RuntimeError(
            "v-u != 2x"
        )

    if v + u != 2 * y:
        raise RuntimeError(
            "v+u != 2y"
        )

    return LevelState(
        n=n,
        p=p,
        q=q,
        z=z,
        k=k,
        modulus=modulus,
        branch=branch,
        x=x,
        y=y,
        u=u,
        v=v,
    )


# ==============================================================================
# BUILD LEVEL INDEX
# ==============================================================================

def build_level_index(
    semiprimes,
):

    levels = {}

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        index = {}

        for n, p, q in semiprimes:

            state = construct_level_state(
                n,
                p,
                q,
                z,
            )

            if state is None:
                continue

            # n is sufficient because p,q are uniquely the generated pair.
            index[n] = state

        levels[z] = index

        print(
            f"    z={z:<2} "
            f"valid={len(index)}"
        )

    return levels


# ==============================================================================
# TRANSITION RECORD
# ==============================================================================

@dataclass(frozen=True)
class Transition:

    parent: LevelState
    child: LevelState


# ==============================================================================
# BUILD SAME-N TRANSITIONS
# ==============================================================================

def build_transitions(
    levels,
):

    transitions = defaultdict(list)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        parent_level = levels[z]
        child_level = levels[z + 1]

        for n in (
            set(parent_level)
            & set(child_level)
        ):

            parent = parent_level[n]
            child = child_level[n]

            transitions[z].append(
                Transition(
                    parent=parent,
                    child=child,
                )
            )

    return transitions


# ==============================================================================
# TRANSITION SUMMARY
# ==============================================================================

def transition_summary(
    transitions,
):

    print()
    print("=" * 90)
    print(
        "SAME-N RECURSIVE TRANSITIONS"
    )
    print("=" * 90)

    for z in sorted(
        transitions
    ):

        items = transitions[z]

        print()
        print(
            f"z={z} -> z={z+1}"
        )

        print(
            f"    transitions={len(items)}"
        )

        pairs = CounterPairs(items)

        for key, count in sorted(
            pairs.items()
        ):

            parent_branch, child_branch = key

            print(
                f"    {parent_branch}"
                f" -> "
                f"{child_branch}: "
                f"{count}"
            )


# ==============================================================================
# COUNTER WITHOUT IMPORTING Counter
# ==============================================================================

def CounterPairs(
    transitions,
):

    counts = defaultdict(int)

    for transition in transitions:

        key = (
            transition.parent.branch,
            transition.child.branch,
        )

        counts[key] += 1

    return counts


# ==============================================================================
# TEST EXACT COORDINATE RULE
# ==============================================================================

def test_formula(
    transitions,
    name,
    function,
    attribute,
):

    mapping = defaultdict(set)

    for transition in transitions:

        parent_value = getattr(
            transition.parent,
            attribute,
        )

        child_value = getattr(
            transition.child,
            attribute,
        )

        predicted = function(
            parent_value
        )

        mapping[
            parent_value
        ].add(
            child_value - predicted
        )

    # The function is exact if the residual is always 0.
    exact = all(
        residuals == {0}
        for residuals
        in mapping.values()
    )

    return exact


# ==============================================================================
# SEARCH AFFINE TRANSFORMATION
# ==============================================================================

def search_affine_transform(
    transitions,
    attribute,
):

    """
    Search:

        child = a * parent + b

    with small rational a and integer b.
    """

    candidates_a = (
        Fraction(-2),
        Fraction(-1),
        Fraction(-1, 2),
        Fraction(1, 2),
        Fraction(1),
        Fraction(2),
    )

    candidates_b = range(
        -8,
        9,
    )

    results = []

    for a in candidates_a:

        for b in candidates_b:

            valid = True

            for transition in transitions:

                parent_value = getattr(
                    transition.parent,
                    attribute,
                )

                child_value = getattr(
                    transition.child,
                    attribute,
                )

                predicted = (
                    a * parent_value
                    + b
                )

                if predicted != child_value:

                    valid = False
                    break

            if valid:

                results.append(
                    (
                        a,
                        b,
                    )
                )

    return results


# ==============================================================================
# SEARCH BRANCH-CONDITIONED TRANSFORM
# ==============================================================================

def branch_conditioned_transforms(
    transitions,
    attribute,
):

    groups = defaultdict(list)

    for transition in transitions:

        key = (
            transition.parent.branch,
            transition.child.branch,
        )

        groups[key].append(
            transition
        )

    result = {}

    for key, group in groups.items():

        rules = search_affine_transform(
            group,
            attribute,
        )

        result[key] = rules

    return result


# ==============================================================================
# COORDINATE TRANSFORMATION REPORT
# ==============================================================================

def coordinate_transformation_report(
    transitions,
):

    print()
    print("=" * 90)
    print(
        "COORDINATE TRANSFORMATIONS"
    )
    print("=" * 90)

    for z in sorted(
        transitions
    ):

        items = transitions[z]

        if not items:
            continue

        print()
        print(
            f"LEVEL {z} -> {z+1}"
        )

        groups = branch_conditioned_transforms(
            items,
            "x",
        )

        for branch_pair in sorted(
            groups
        ):

            parent_branch, child_branch = (
                branch_pair
            )

            print()
            print(
                f"    "
                f"{parent_branch} -> "
                f"{child_branch}"
            )

            rules = groups[
                branch_pair
            ]

            if not rules:

                print(
                    "        "
                    "x_child = a*x_parent+b: NONE"
                )

            else:

                for a, b in rules:

                    print(
                        f"        "
                        f"x_child = "
                        f"({a})*x_parent "
                        f"+ ({b})"
                    )

        # ----------------------------------------------------------------------
        # y
        # ----------------------------------------------------------------------

        groups_y = branch_conditioned_transforms(
            items,
            "y",
        )

        print()

        for branch_pair in sorted(
            groups_y
        ):

            parent_branch, child_branch = (
                branch_pair
            )

            print()
            print(
                f"    "
                f"{parent_branch} -> "
                f"{child_branch}"
            )

            rules = groups_y[
                branch_pair
            ]

            if not rules:

                print(
                    "        "
                    "y_child = a*y_parent+b: NONE"
                )

            else:

                for a, b in rules:

                    print(
                        f"        "
                        f"y_child = "
                        f"({a})*y_parent "
                        f"+ ({b})"
                    )


# ==============================================================================
# TEST CANDIDATE HALF-RULES
# ==============================================================================

def half_rule_report(
    transitions,
):

    print()
    print("=" * 90)
    print(
        "HALVING / HALF-OFFSET TEST"
    )
    print("=" * 90)

    candidate_names = (
        "x/2",
        "(x-1)/2",
        "(x+1)/2",
        "y/2",
        "(y-1)/2",
        "(y+1)/2",
    )

    for z in sorted(
        transitions
    ):

        items = transitions[z]

        if not items:
            continue

        print()
        print(
            f"LEVEL {z} -> {z+1}"
        )

        groups = defaultdict(list)

        for transition in items:

            groups[
                (
                    transition.parent.branch,
                    transition.child.branch,
                )
            ].append(
                transition
            )

        for pair in sorted(groups):

            group = groups[pair]

            pb, cb = pair

            print()
            print(
                f"    {pb} -> {cb}"
            )

            tests = (
                (
                    "x/2",
                    lambda x: Fraction(x, 2),
                    "x",
                ),
                (
                    "(x-1)/2",
                    lambda x: Fraction(
                        x - 1,
                        2,
                    ),
                    "x",
                ),
                (
                    "(x+1)/2",
                    lambda x: Fraction(
                        x + 1,
                        2,
                    ),
                    "x",
                ),
                (
                    "y/2",
                    lambda y: Fraction(y, 2),
                    "y",
                ),
                (
                    "(y-1)/2",
                    lambda y: Fraction(
                        y - 1,
                        2,
                    ),
                    "y",
                ),
                (
                    "(y+1)/2",
                    lambda y: Fraction(
                        y + 1,
                        2,
                    ),
                    "y",
                ),
            )

            for name, function, attribute in tests:

                exact = True

                for transition in group:

                    parent_value = getattr(
                        transition.parent,
                        attribute,
                    )

                    child_value = getattr(
                        transition.child,
                        attribute,
                    )

                    if function(
                        parent_value
                    ) != child_value:

                        exact = False
                        break

                if exact:

                    print(
                        f"        EXACT: "
                        f"{name}"
                    )


# ==============================================================================
# SURVIVAL CONDITION
# ==============================================================================

def survival_report(
    transitions,
):

    print()
    print("=" * 90)
    print(
        "SURVIVAL / CHILD-DOMAIN CONDITIONS"
    )
    print("=" * 90)

    for z in sorted(
        transitions
    ):

        items = transitions[z]

        if not items:
            continue

        print()
        print(
            f"LEVEL {z} -> {z+1}"
        )

        for parent_branch in (
            "A",
            "B",
        ):

            parent_items = [
                t
                for t in items
                if t.parent.branch
                == parent_branch
            ]

            if not parent_items:
                continue

            print()
            print(
                f"    PARENT {parent_branch}"
            )

            # --------------------------------------------------------------
            # Parent x parity
            # --------------------------------------------------------------

            parity = defaultdict(
                int
            )

            for transition in parent_items:

                parity[
                    transition.parent.x & 1
                ] += 1

            for value in sorted(
                parity
            ):

                print(
                    f"        surviving "
                    f"x%2={value}: "
                    f"{parity[value]}"
                )

            # --------------------------------------------------------------
            # Branch transition by x parity.
            # --------------------------------------------------------------

            mapping = defaultdict(
                set
            )

            for transition in parent_items:

                mapping[
                    transition.parent.x & 1
                ].add(
                    transition.child.branch
                )

            for value in sorted(
                mapping
            ):

                print(
                    f"        x%2={value} "
                    f"-> child branches "
                    f"{sorted(mapping[value])}"
                )


# ==============================================================================
# NORMALIZED COORDINATE TRANSFORMATION
# ==============================================================================

def uv_transformation_report(
    transitions,
):

    print()
    print("=" * 90)
    print(
        "u,v TRANSFORMATION"
    )
    print("=" * 90)

    for z in sorted(
        transitions
    ):

        items = transitions[z]

        if not items:
            continue

        print()
        print(
            f"LEVEL {z} -> {z+1}"
        )

        groups = defaultdict(list)

        for transition in items:

            groups[
                (
                    transition.parent.branch,
                    transition.child.branch,
                )
            ].append(
                transition
            )

        for pair in sorted(groups):

            group = groups[pair]

            print()
            print(
                f"    {pair[0]} -> {pair[1]}"
            )

            for attribute in (
                "u",
                "v",
            ):

                rules = search_affine_transform(
                    group,
                    attribute,
                )

                if rules:

                    for a, b in rules:

                        print(
                            f"        "
                            f"{attribute}_child = "
                            f"({a})*"
                            f"{attribute}_parent "
                            f"+ ({b})"
                        )

                else:

                    print(
                        f"        "
                        f"{attribute}: "
                        f"no small affine rule"
                    )


# ==============================================================================
# EXAMPLES
# ==============================================================================

def example_transitions(
    transitions,
):

    if not SHOW_EXAMPLES:
        return

    print()
    print("=" * 90)
    print(
        "EXAMPLE SAME-N TRANSITIONS"
    )
    print("=" * 90)

    for z in sorted(
        transitions
    ):

        print()
        print(
            f"LEVEL {z} -> {z+1}"
        )

        for transition in transitions[z][
            :MAX_EXAMPLES
        ]:

            p = transition.parent
            c = transition.child

            print(
                f"    n={p.n} "
                f"{p.branch}({p.x},{p.y})"
                f" -> "
                f"{c.branch}({c.x},{c.y})"
            )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print(
        "EXPERIMENT 586 START"
    )
    print("=" * 90)

    print()
    print(
        "SAME-N RECURSIVE COORDINATE TRANSFORMATION"
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
        f"    semiprimes = "
        f"{len(semiprimes)}"
    )

    # --------------------------------------------------------------------------
    # Build levels.
    # --------------------------------------------------------------------------

    print()
    print(
        "[3] BUILD LEVEL INDEX"
    )

    levels = build_level_index(
        semiprimes
    )

    # --------------------------------------------------------------------------
    # Build same-n transitions.
    # --------------------------------------------------------------------------

    print()
    print(
        "[4] BUILD SAME-N TRANSITIONS"
    )

    transitions = build_transitions(
        levels
    )

    # --------------------------------------------------------------------------
    # Summary.
    # --------------------------------------------------------------------------

    transition_summary(
        transitions
    )

    # --------------------------------------------------------------------------
    # Coordinate transforms.
    # --------------------------------------------------------------------------

    coordinate_transformation_report(
        transitions
    )

    # --------------------------------------------------------------------------
    # Half rules.
    # --------------------------------------------------------------------------

    half_rule_report(
        transitions
    )

    # --------------------------------------------------------------------------
    # Survival conditions.
    # --------------------------------------------------------------------------

    survival_report(
        transitions
    )

    # --------------------------------------------------------------------------
    # u/v.
    # --------------------------------------------------------------------------

    uv_transformation_report(
        transitions
    )

    # --------------------------------------------------------------------------
    # Examples.
    # --------------------------------------------------------------------------

    example_transitions(
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
This experiment differs from the previous experiments in one critical way:

    ONLY THE SAME n,p,q IS COMPARED BETWEEN z AND z+1.

Therefore a transition exists only when the factor pair survives both
modular constructions.

The measured object is:

    (x_z,y_z,branch_z)
        ->
    (x_(z+1),y_(z+1),branch_(z+1))

The main target is a static transformation such as:

    x_(z+1) = x_z / 2

or:

    x_(z+1) = (x_z +/- 1) / 2

possibly conditioned on:

    parent branch
    child branch.

Likewise for y,u,v.

If the same transformation appears repeatedly at consecutive levels,
that would be direct evidence for the recursive shrinking-coordinate
mechanism you described.
"""
    )

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 586 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()
