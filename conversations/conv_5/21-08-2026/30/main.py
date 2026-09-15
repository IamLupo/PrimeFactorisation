#!/usr/bin/env python3

# ==============================================================================
# EXPERIMENT 587
# ==============================================================================
# LEVEL-INVARIANT COORDINATES
#
# NO FILES
# NO WEB
#
# Hypothesis from Experiment 586:
#
#     x_(z+1) = x_z / 2
#     y_(z+1) = y_z / 2
#
# while
#
#     k_z = 2^(z-1)
#
# doubles:
#
#     k_(z+1) = 2*k_z.
#
# Therefore:
#
#     X = k*x
#     Y = k*y
#
# should be invariant across levels.
#
#
# The experiment tests:
#
#     X_z == X_(z+1)
#     Y_z == Y_(z+1)
#
# for the SAME n,p,q.
#
# Then it tests whether the factor functions collapse to:
#
# CASE A:
#
#     p = Y-X+3
#     q = Y+X-3
#
# CASE B:
#
#     3p = Y+X-3
#     q  = Y-X+3
#
#
# If this is exact, then the apparent A/B/C/D/... hierarchy may be a
# coordinate-renormalization hierarchy rather than genuinely new algebra.
#
# ==============================================================================

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
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
# DATA
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

    X: int
    Y: int


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
# CASE A
# ==============================================================================

def solve_case_a(
    p: int,
    q: int,
    k: int,
):

    # p = k(y-x) + 3
    # q = k(y+x) - 3

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
# CASE B
# ==============================================================================

def solve_case_b(
    p: int,
    q: int,
    k: int,
):

    # 3p = k(y+x) - 3
    # q  = k(y-x) + 3

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

def construct_level(
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
    # Branch A
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
    # Branch B
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

    # --------------------------------------------------------------------------
    # Invariant coordinates.
    # --------------------------------------------------------------------------

    X = k * x
    Y = k * y

    if X != k * x:
        raise RuntimeError(
            "X invariant construction failed"
        )

    if Y != k * y:
        raise RuntimeError(
            "Y invariant construction failed"
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
        X=X,
        Y=Y,
    )


# ==============================================================================
# BUILD LEVEL INDEXES
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

            state = construct_level(
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
            f"valid={len(level)}"
        )

    return levels


# ==============================================================================
# SAME-N TRANSITIONS
# ==============================================================================

def same_n_transitions(
    levels,
):

    result = defaultdict(list)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        parent = levels[z]
        child = levels[z + 1]

        common = (
            set(parent)
            & set(child)
        )

        for n in common:

            result[z].append(
                (
                    parent[n],
                    child[n],
                )
            )

    return result


# ==============================================================================
# INVARIANT TEST
# ==============================================================================

def test_invariants(
    transitions,
):

    print()
    print("=" * 90)
    print(
        "INVARIANT TEST: X=k*x, Y=k*y"
    )
    print("=" * 90)

    global_X_failures = 0
    global_Y_failures = 0

    for z in sorted(
        transitions
    ):

        pairs = transitions[z]

        if not pairs:
            continue

        X_failures = 0
        Y_failures = 0

        for parent, child in pairs:

            if parent.X != child.X:
                X_failures += 1

            if parent.Y != child.Y:
                Y_failures += 1

        global_X_failures += X_failures
        global_Y_failures += Y_failures

        print()
        print(
            f"z={z} -> z={z+1}"
        )

        print(
            f"    transitions={len(pairs)}"
        )

        print(
            f"    X invariant failures="
            f"{X_failures}"
        )

        print(
            f"    Y invariant failures="
            f"{Y_failures}"
        )

    print()
    print(
        "GLOBAL"
    )

    print(
        f"    X failures = "
        f"{global_X_failures}"
    )

    print(
        f"    Y failures = "
        f"{global_Y_failures}"
    )


# ==============================================================================
# FUNCTION COLLAPSE TEST
# ==============================================================================

def test_function_collapse(
    levels,
):

    print()
    print("=" * 90)
    print(
        "LEVEL-INDEPENDENT FUNCTION TEST"
    )
    print("=" * 90)

    failures_A = 0
    failures_B = 0

    counts_A = defaultdict(int)
    counts_B = defaultdict(int)

    for z in levels:

        for state in levels[z].values():

            X = state.X
            Y = state.Y

            if state.branch == "A":

                predicted_p = (
                    Y - X + 3
                )

                predicted_q = (
                    Y + X - 3
                )

                counts_A[z] += 1

                if (
                    predicted_p != state.p
                    or
                    predicted_q != state.q
                ):

                    failures_A += 1

            else:

                predicted_3p = (
                    Y + X - 3
                )

                predicted_q = (
                    Y - X + 3
                )

                counts_B[z] += 1

                if (
                    predicted_3p != 3 * state.p
                    or
                    predicted_q != state.q
                ):

                    failures_B += 1

    print()
    print(
        "CASE A:"
    )

    print(
        f"    tested={sum(counts_A.values())}"
    )

    print(
        f"    failures={failures_A}"
    )

    for z in sorted(counts_A):

        print(
            f"    z={z:<2} "
            f"count={counts_A[z]}"
        )

    print()
    print(
        "CASE B:"
    )

    print(
        f"    tested={sum(counts_B.values())}"
    )

    print(
        f"    failures={failures_B}"
    )

    for z in sorted(counts_B):

        print(
            f"    z={z:<2} "
            f"count={counts_B[z]}"
        )


# ==============================================================================
# COMPARE INVARIANT COORDINATES DIRECTLY
# ==============================================================================

def invariant_collisions(
    levels,
):

    print()
    print("=" * 90)
    print(
        "INVARIANT COORDINATE COLLAPSE"
    )
    print("=" * 90)

    by_XY = defaultdict(set)

    for z in levels:

        for state in levels[z].values():

            by_XY[
                (
                    state.X,
                    state.Y,
                )
            ].add(
                (
                    state.p,
                    state.q,
                )
            )

    ambiguous = 0

    for key, factors in by_XY.items():

        if len(factors) > 1:
            ambiguous += 1

    print(
        f"unique (X,Y) states = "
        f"{len(by_XY)}"
    )

    print(
        f"ambiguous (X,Y) -> multiple (p,q) = "
        f"{ambiguous}"
    )

    if ambiguous == 0:

        print(
            "    EVERY invariant coordinate pair "
            "maps to one factor pair."
        )

    else:

        print(
            "    Some invariant coordinates "
            "map to multiple factor pairs."
        )


# ==============================================================================
# LEVEL REPRESENTATION TEST
# ==============================================================================

def representation_examples(
    transitions,
):

    print()
    print("=" * 90)
    print(
        "SAME FACTOR, DIFFERENT LEVEL REPRESENTATION"
    )
    print("=" * 90)

    shown = 0

    for z in sorted(
        transitions
    ):

        for parent, child in transitions[z]:

            print()
            print(
                f"n={parent.n}"
            )

            print(
                f"    z={parent.z}: "
                f"k={parent.k} "
                f"x={parent.x} "
                f"y={parent.y} "
                f"X={parent.X} "
                f"Y={parent.Y} "
                f"branch={parent.branch}"
            )

            print(
                f"    z={child.z}: "
                f"k={child.k} "
                f"x={child.x} "
                f"y={child.y} "
                f"X={child.X} "
                f"Y={child.Y} "
                f"branch={child.branch}"
            )

            shown += 1

            if shown >= MAX_EXAMPLES:
                return


# ==============================================================================
# SEARCH FOR LEVEL-INDEPENDENT FORM
# ==============================================================================

def factor_formula_check(
    levels,
):

    print()
    print("=" * 90)
    print(
        "DIRECT FACTOR FORMULAS IN X,Y"
    )
    print("=" * 90)

    formulas = {
        "A-p": 0,
        "A-q": 0,
        "B-3p": 0,
        "B-q": 0,
    }

    tested = {
        key: 0
        for key in formulas
    }

    for z in levels:

        for state in levels[z].values():

            if state.branch == "A":

                tested["A-p"] += 1

                if (
                    state.p
                    == state.Y
                    - state.X
                    + 3
                ):
                    formulas["A-p"] += 1

                tested["A-q"] += 1

                if (
                    state.q
                    == state.Y
                    + state.X
                    - 3
                ):
                    formulas["A-q"] += 1

            else:

                tested["B-3p"] += 1

                if (
                    3 * state.p
                    ==
                    state.Y
                    + state.X
                    - 3
                ):
                    formulas["B-3p"] += 1

                tested["B-q"] += 1

                if (
                    state.q
                    ==
                    state.Y
                    - state.X
                    + 3
                ):
                    formulas["B-q"] += 1

    for key in formulas:

        print(
            f"    {key:<5}: "
            f"{formulas[key]}/"
            f"{tested[key]}"
        )


# ==============================================================================
# SAT STATE COUNT
# ==============================================================================

def sat_state_report(
    levels,
):

    print()
    print("=" * 90)
    print(
        "INVARIANT STATE COUNT"
    )
    print("=" * 90)

    all_xy = set()

    for z in levels:

        level_xy = {
            (
                state.X,
                state.Y,
            )
            for state in levels[z].values()
        }

        all_xy |= level_xy

        print(
            f"z={z:<2} "
            f"level states={len(level_xy):<8} "
            f"cumulative invariant states="
            f"{len(all_xy)}"
        )


# ==============================================================================
# EXAMPLES
# ==============================================================================

def print_examples(
    levels,
):

    if not SHOW_EXAMPLES:
        return

    print()
    print("=" * 90)
    print(
        "EXAMPLES"
    )
    print("=" * 90)

    shown = 0

    for z in sorted(levels):

        for state in levels[z].values():

            print(
                f"    z={z} "
                f"n={state.n} "
                f"p={state.p} "
                f"q={state.q} "
                f"k={state.k} "
                f"x={state.x} "
                f"y={state.y} "
                f"X={state.X} "
                f"Y={state.Y} "
                f"branch={state.branch}"
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
        "EXPERIMENT 587 START"
    )
    print("=" * 90)

    print()
    print(
        "LEVEL-INVARIANT COORDINATES"
    )

    # --------------------------------------------------------------------------
    # Generate semiprimes.
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
    # Construct levels.
    # --------------------------------------------------------------------------

    print()
    print(
        "[3] BUILD LEVELS"
    )

    levels = build_levels(
        semiprimes
    )

    # --------------------------------------------------------------------------
    # Same-n transitions.
    # --------------------------------------------------------------------------

    print()
    print(
        "[4] SAME-N TRANSITIONS"
    )

    transitions = same_n_transitions(
        levels
    )

    for z in sorted(
        transitions
    ):

        print(
            f"    z={z} -> z={z+1}: "
            f"{len(transitions[z])}"
        )

    # --------------------------------------------------------------------------
    # Test invariant.
    # --------------------------------------------------------------------------

    test_invariants(
        transitions
    )

    # --------------------------------------------------------------------------
    # Test functions.
    # --------------------------------------------------------------------------

    test_function_collapse(
        levels
    )

    factor_formula_check(
        levels
    )

    # --------------------------------------------------------------------------
    # Collision test.
    # --------------------------------------------------------------------------

    invariant_collisions(
        levels
    )

    # --------------------------------------------------------------------------
    # State count.
    # --------------------------------------------------------------------------

    sat_state_report(
        levels
    )

    # --------------------------------------------------------------------------
    # Examples.
    # --------------------------------------------------------------------------

    representation_examples(
        transitions
    )

    print_examples(
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
Hypothesis:

    k_z = 2^(z-1)

    X = k_z*x_z
    Y = k_z*y_z

If:

    x_(z+1) = x_z/2
    y_(z+1) = y_z/2

then:

    X_(z+1) = X_z
    Y_(z+1) = Y_z.

The experiment therefore asks whether X,Y are the true
level-independent coordinates behind the apparent hierarchy.

If the function tests also pass:

    CASE A:
        p = Y-X+3
        q = Y+X-3

    CASE B:
        3p = Y+X-3
        q  = Y-X+3

at every level, then the A -> B -> C -> ... functions are
representations of the same invariant algebra under the coordinate
renormalization:

    (x,y,k)
        ->
    (x/2,y/2,2k).

That is the structure we want to establish before doing a SAT
encoding experiment.
"""
    )

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 587 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()
