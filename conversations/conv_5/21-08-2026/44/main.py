#!/usr/bin/env python3

# ==============================================================================
# EXPERIMENT 603
# ==============================================================================
# COMPLETE GLOBAL RECONSTRUCTION FROM (FRAME, X, Y)
#
# Goal
# ----
#
# Experiment 602 established:
#
#     X = 2^(z-1) * x_z
#     Y = 2^(z-1) * y_z
#
# and:
#
#     child exists at z -> z+1
#     iff
#     2^z | X
#     and
#     2^z | Y.
#
# The previous bitwise test used the wrong bit index.
#
# Correctly:
#
#     transition z -> z+1
#     depends on bit z of X and Y
#
# because:
#
#     2^z | X,Y.
#
# This experiment goes further:
#
# GIVEN ONLY:
#
#     frame
#     X
#     Y
#
# reconstruct:
#
#     x_z
#     y_z
#     p
#     q
#     n
#     residue mod 2^z
#     whether level z exists
#     whether level z+1 exists
#
# We test whether the COMPLETE observed hierarchy can be reconstructed
# without using the original level-specific formulas as input.
#
# Two frames:
#
# FRAME A
#
#     p = Y - X + 3
#     q = Y + X - 3
#
# FRAME B
#
#     3p = Y + X - 3
#     q  = Y - X + 3
#
# Therefore:
#
#     FRAME B:
#         p = (Y + X - 3) / 3
#
# The experiment tests:
#
#     n = p*q
#
# and then:
#
#     n mod 2^z
#
# against the observed residue.
#
# IMPORTANT
# =========
#
# No external files.
# No external sources.
# No output files.
# ==============================================================================


from __future__ import annotations

from dataclasses import dataclass
from collections import Counter
from math import isqrt


# ==============================================================================
# CONFIGURATION
# ==============================================================================

MAX_N = 2_000_000

MIN_Z = 2
MAX_Z = 12

SHOW_EXAMPLES = 20


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

    x: int
    y: int

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
            start:
            limit + 1:
            p
        ] = b"\x00" * count

    return sieve


# ==============================================================================
# SEMIPRIME GENERATION
# ==============================================================================

def generate_semiprimes(
    limit,
    sieve,
):

    primes = [
        p
        for p in range(
            3,
            limit + 1,
            2,
        )
        if sieve[p]
    ]

    result = []

    for i, p in enumerate(primes):

        if p * p > limit:
            break

        max_q = limit // p

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
# FRAME A REPRESENTATION
# ==============================================================================

def frame_A_coordinates(
    p,
    q,
    z,
):

    k = 1 << (
        z - 1
    )

    denominator = 2 * k

    x_num = (
        q - p + 6
    )

    y_num = (
        p + q
    )

    if x_num % denominator:
        return None

    if y_num % denominator:
        return None

    return (
        x_num // denominator,
        y_num // denominator,
    )


# ==============================================================================
# FRAME B REPRESENTATION
# ==============================================================================

def frame_B_coordinates(
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

    if x_num % denominator:
        return None

    if y_num % denominator:
        return None

    return (
        x_num // denominator,
        y_num // denominator,
    )


# ==============================================================================
# BUILD OBSERVED LEVELS
# ==============================================================================

def build_levels(
    semiprimes,
):

    levels = {
        z: {}
        for z in range(
            MIN_Z,
            MAX_Z + 1,
        )
    }

    for n, p, q in semiprimes:

        for z in range(
            MIN_Z,
            MAX_Z + 1,
        ):

            a = frame_A_coordinates(
                p,
                q,
                z,
            )

            if a is not None:

                k = 1 << (
                    z - 1
                )

                levels[z][n] = State(
                    n=n,
                    p=p,
                    q=q,
                    z=z,
                    residue=n % (
                        1 << z
                    ),
                    frame="A",
                    x=a[0],
                    y=a[1],
                    X=k * a[0],
                    Y=k * a[1],
                )

                continue

            b = frame_B_coordinates(
                p,
                q,
                z,
            )

            if b is not None:

                k = 1 << (
                    z - 1
                )

                levels[z][n] = State(
                    n=n,
                    p=p,
                    q=q,
                    z=z,
                    residue=n % (
                        1 << z
                    ),
                    frame="B",
                    x=b[0],
                    y=b[1],
                    X=k * b[0],
                    Y=k * b[1],
                )

    return levels


# ==============================================================================
# INVARIANT -> LEVEL COORDINATES
# ==============================================================================

def xy_from_invariants(
    X,
    Y,
    z,
):

    k = 1 << (
        z - 1
    )

    if X % k:
        return None

    if Y % k:
        return None

    return (
        X // k,
        Y // k,
    )


# ==============================================================================
# RECONSTRUCT FACTORS FROM (FRAME, X, Y)
# ==============================================================================

def factors_from_invariants(
    frame,
    X,
    Y,
):

    if frame == "A":

        p = (
            Y
            - X
            + 3
        )

        q = (
            Y
            + X
            - 3
        )

        return p, q

    if frame == "B":

        numerator = (
            Y
            + X
            - 3
        )

        if numerator % 3:
            return None

        p = numerator // 3

        q = (
            Y
            - X
            + 3
        )

        return p, q

    return None


# ==============================================================================
# RECONSTRUCT N
# ==============================================================================

def reconstruct_n(
    frame,
    X,
    Y,
):

    factors = factors_from_invariants(
        frame,
        X,
        Y,
    )

    if factors is None:
        return None

    p, q = factors

    return p * q


# ==============================================================================
# CORRECT SURVIVAL TEST
# ==============================================================================

def survives(
    X,
    Y,
    z,
):

    modulus = 1 << z

    return (
        X % modulus == 0
        and
        Y % modulus == 0
    )


# ==============================================================================
# CORRECT BIT TEST
#
# Transition z -> z+1 requires:
#
#     bit z of X = 0
#     bit z of Y = 0.
#
# ==============================================================================

def bit_at(
    value,
    position,
):

    return (
        abs(value) >> position
    ) & 1


def survives_by_bit(
    X,
    Y,
    z,
):

    return (
        bit_at(X, z) == 0
        and
        bit_at(Y, z) == 0
    )


# ==============================================================================
# 2-ADIC VALUATION
# ==============================================================================

def v2(value):

    if value == 0:
        return 10**9

    value = abs(value)

    result = 0

    while (
        value & 1
    ) == 0:

        value >>= 1
        result += 1

    return result


# ==============================================================================
# OBSERVED CHILD
# ==============================================================================

def child(
    levels,
    z,
    n,
):

    return levels[
        z + 1
    ].get(n)


# ==============================================================================
# TEST 1
# ==============================================================================

def test_global_invariant(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 1: GLOBAL INVARIANT"
    )
    print("=" * 90)

    known = {}

    checked = 0
    failures = 0

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        for n, state in levels[z].items():

            current = (
                state.X,
                state.Y,
            )

            if n in known:

                if known[n] != current:

                    failures += 1

                    if failures <= SHOW_EXAMPLES:

                        print(
                            f"    mismatch "
                            f"n={n} "
                            f"z={z} "
                            f"old={known[n]} "
                            f"new={current}"
                        )

            else:

                known[n] = current

            checked += 1

    print(
        f"checked={checked} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 2
# ==============================================================================

def test_coordinate_reconstruction(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 2: X,Y -> x_z,y_z"
    )
    print("=" * 90)

    total = 0
    failures = 0

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        local_failures = 0

        for state in levels[z].values():

            predicted = (
                xy_from_invariants(
                    state.X,
                    state.Y,
                    z,
                )
            )

            total += 1

            if predicted != (
                state.x,
                state.y,
            ):

                failures += 1
                local_failures += 1

        print(
            f"z={z} "
            f"failures={local_failures}"
        )

    print()
    print(
        f"GLOBAL checked={total} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 3
# ==============================================================================

def test_factor_reconstruction(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 3: (FRAME,X,Y) -> (p,q)"
    )
    print("=" * 90)

    total = 0
    failures = 0

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        local = 0

        for state in levels[z].values():

            predicted = factors_from_invariants(
                state.frame,
                state.X,
                state.Y,
            )

            total += 1

            if predicted != (
                state.p,
                state.q,
            ):

                failures += 1
                local += 1

                if local <= SHOW_EXAMPLES:

                    print(
                        f"    z={z} "
                        f"n={state.n} "
                        f"frame={state.frame} "
                        f"pred={predicted} "
                        f"actual=({state.p},{state.q})"
                    )

        print(
            f"z={z} "
            f"failures={local}"
        )

    print()
    print(
        f"GLOBAL checked={total} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 4
# ==============================================================================

def test_n_reconstruction(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 4: (FRAME,X,Y) -> n"
    )
    print("=" * 90)

    total = 0
    failures = 0

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        local = 0

        for state in levels[z].values():

            predicted = reconstruct_n(
                state.frame,
                state.X,
                state.Y,
            )

            total += 1

            if predicted != state.n:

                failures += 1
                local += 1

        print(
            f"z={z} "
            f"failures={local}"
        )

    print()
    print(
        f"GLOBAL checked={total} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 5
# ==============================================================================

def test_residue_reconstruction(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 5: (FRAME,X,Y,z) -> residue"
    )
    print("=" * 90)

    total = 0
    failures = 0

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        local = 0

        modulus = 1 << z

        for state in levels[z].values():

            reconstructed_n = reconstruct_n(
                state.frame,
                state.X,
                state.Y,
            )

            predicted = (
                reconstructed_n
                % modulus
                if reconstructed_n is not None
                else None
            )

            total += 1

            if predicted != state.residue:

                failures += 1
                local += 1

        print(
            f"z={z} "
            f"mod={modulus} "
            f"failures={local}"
        )

    print()
    print(
        f"GLOBAL checked={total} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 6
# ==============================================================================

def test_survival(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 6: DIRECT SURVIVAL"
    )
    print("=" * 90)

    total = 0
    failures = 0

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        local = 0

        for state in levels[z].values():

            actual = (
                child(
                    levels,
                    z,
                    state.n,
                )
                is not None
            )

            predicted = survives(
                state.X,
                state.Y,
                z,
            )

            total += 1

            if predicted != actual:

                failures += 1
                local += 1

        print(
            f"z={z}->{z+1} "
            f"checked={len(levels[z])} "
            f"failures={local}"
        )

    print()
    print(
        f"GLOBAL checked={total} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 7
# ==============================================================================

def test_bit_survival(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 7: CORRECTED BITWISE SURVIVAL"
    )
    print("=" * 90)

    print(
        """
For transition z -> z+1:

    survival <-> bit_z(X)=0 AND bit_z(Y)=0.

This is zero-indexed bit position z.

Examples:

    z=2 requires bit 2
    z=3 requires bit 3
    z=4 requires bit 4
    ...
"""
    )

    total = 0
    failures = 0

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        local = 0

        for state in levels[z].values():

            actual = (
                child(
                    levels,
                    z,
                    state.n,
                )
                is not None
            )

            predicted = survives_by_bit(
                state.X,
                state.Y,
                z,
            )

            total += 1

            if predicted != actual:

                failures += 1
                local += 1

                if local <= SHOW_EXAMPLES:

                    print(
                        f"    mismatch "
                        f"z={z} "
                        f"n={state.n} "
                        f"X={state.X} "
                        f"Y={state.Y} "
                        f"Xbit={bit_at(state.X,z)} "
                        f"Ybit={bit_at(state.Y,z)} "
                        f"pred={predicted} "
                        f"actual={actual}"
                    )

        print(
            f"z={z}->{z+1} "
            f"failures={local}"
        )

    print()
    print(
        f"GLOBAL checked={total} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 8
# ==============================================================================

def test_deepest_level(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 8: DEEPEST LEVEL FROM v2(X),v2(Y)"
    )
    print("=" * 90)

    failures = 0
    checked = 0

    distribution = Counter()

    examples = []

    for n, state in levels[
        MIN_Z
    ].items():

        vx = v2(
            state.X
        )

        vy = v2(
            state.Y
        )

        predicted = (
            min(vx, vy)
            + 1
        )

        actual = MIN_Z - 1

        for z in range(
            MIN_Z,
            MAX_Z + 1,
        ):

            if n in levels[z]:

                actual = z

        if predicted <= MAX_Z:

            checked += 1

            if predicted != actual:

                failures += 1

                if len(examples) < SHOW_EXAMPLES:

                    examples.append(
                        (
                            n,
                            state.p,
                            state.q,
                            state.X,
                            state.Y,
                            vx,
                            vy,
                            predicted,
                            actual,
                        )
                    )

        distribution[
            actual
        ] += 1

    print(
        f"bounded checked={checked} "
        f"failures={failures}"
    )

    if examples:

        print()
        print(
            "COUNTEREXAMPLES"
        )

        for e in examples:

            print(
                "    "
                f"n={e[0]} "
                f"p={e[1]} "
                f"q={e[2]} "
                f"X={e[3]} "
                f"Y={e[4]} "
                f"v2X={e[5]} "
                f"v2Y={e[6]} "
                f"pred={e[7]} "
                f"actual={e[8]}"
            )

    print()
    print(
        "DEPTH DISTRIBUTION"
    )

    for depth, count in sorted(
        distribution.items()
    ):

        print(
            f"    depth={depth}: "
            f"{count}"
        )


# ==============================================================================
# TEST 9
# ==============================================================================

def test_frame_invariance(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 9: FRAME INVARIANCE"
    )
    print("=" * 90)

    failures = 0
    checked = 0

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        for state in levels[z].values():

            expected = (
                "A"
                if state.n % 4 == 3
                else "B"
            )

            checked += 1

            if expected != state.frame:

                failures += 1

    print(
        f"checked={checked} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 10
# ==============================================================================

def test_minimal_state(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 10: MINIMAL STATE"
    )
    print("=" * 90)

    print(
        """
For each transition z -> z+1 we compare:

    frame + Xbit + Ybit

against

    Xbit + Ybit

because frame may be derivable entirely from n mod 4.

We test whether survival itself needs frame.
"""
    )

    groups_full = {}
    groups_xy = {}

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        for state in levels[z].values():

            actual = (
                child(
                    levels,
                    z,
                    state.n,
                )
                is not None
            )

            key_full = (
                state.frame,
                bit_at(
                    state.X,
                    z,
                ),
                bit_at(
                    state.Y,
                    z,
                ),
            )

            key_xy = (
                bit_at(
                    state.X,
                    z,
                ),
                bit_at(
                    state.Y,
                    z,
                ),
            )

            groups_full.setdefault(
                key_full,
                set()
            ).add(actual)

            groups_xy.setdefault(
                key_xy,
                set()
            ).add(actual)

    full_deterministic = all(
        len(values) == 1
        for values in groups_full.values()
    )

    xy_deterministic = all(
        len(values) == 1
        for values in groups_xy.values()
    )

    print(
        f"frame+Xbit+Ybit "
        f"deterministic={full_deterministic}"
    )

    print(
        f"Xbit+Ybit "
        f"deterministic={xy_deterministic}"
    )

    print()

    for key in sorted(
        groups_xy
    ):

        print(
            f"    "
            f"{key} -> "
            f"{sorted(groups_xy[key])}"
        )


# ==============================================================================
# TEST 11
# ==============================================================================

def test_full_reconstruction(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 11: COMPLETE GLOBAL RECONSTRUCTION"
    )
    print("=" * 90)

    total = 0
    failures = 0

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        local = 0

        for state in levels[z].values():

            total += 1

            reconstructed = reconstruct_n(
                state.frame,
                state.X,
                state.Y,
            )

            if reconstructed != state.n:

                failures += 1
                local += 1
                continue

            coordinates = (
                xy_from_invariants(
                    state.X,
                    state.Y,
                    z,
                )
            )

            if coordinates != (
                state.x,
                state.y,
            ):

                failures += 1
                local += 1
                continue

            predicted_residue = (
                reconstructed
                % (1 << z)
            )

            if predicted_residue != (
                state.residue
            ):

                failures += 1
                local += 1
                continue

        print(
            f"z={z} "
            f"failures={local}"
        )

    print()
    print(
        f"GLOBAL checked={total} "
        f"failures={failures}"
    )


# ==============================================================================
# EXAMPLES
# ==============================================================================

def print_examples(
    levels,
):

    print()
    print("=" * 90)
    print(
        "EXAMPLES OF GLOBAL RECONSTRUCTION"
    )
    print("=" * 90)

    shown = 0

    for state in levels[
        MIN_Z
    ].values():

        if shown >= SHOW_EXAMPLES:
            break

        vx = v2(
            state.X
        )

        vy = v2(
            state.Y
        )

        deepest = (
            min(vx, vy)
            + 1
        )

        p2, q2 = factors_from_invariants(
            state.frame,
            state.X,
            state.Y,
        )

        print()
        print(
            f"n={state.n}"
        )

        print(
            f"    frame={state.frame}"
        )

        print(
            f"    X={state.X} "
            f"Y={state.Y}"
        )

        print(
            f"    reconstructed "
            f"(p,q)=({p2},{q2})"
        )

        print(
            f"    v2(X)={vx} "
            f"v2(Y)={vy}"
        )

        print(
            f"    predicted deepest={deepest}"
        )

        print(
            "    levels:"
        )

        for z in range(
            MIN_Z,
            min(
                MAX_Z,
                deepest,
            ) + 1,
        ):

            xy = xy_from_invariants(
                state.X,
                state.Y,
                z,
            )

            if xy is None:
                continue

            residue = (
                state.n
                % (1 << z)
            )

            print(
                f"        z={z:<2} "
                f"x={xy[0]:<8} "
                f"y={xy[1]:<8} "
                f"residue={residue}"
            )

        shown += 1


# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

def final_summary():

    print()
    print("=" * 90)
    print(
        "FINAL STRUCTURAL SUMMARY"
    )
    print("=" * 90)

    print(
        r"""
The target global representation is:

    STATE = (FRAME, X, Y)

where X,Y do not depend on z.

FRAME is already determined at the base level:

    n mod 4 = 3  -> A
    n mod 4 = 1  -> B.

At arbitrary level z:

    k_z = 2^(z-1)

    x_z = X / k_z
    y_z = Y / k_z.

The factors are reconstructed directly:

FRAME A:

    p = Y - X + 3
    q = Y + X - 3

FRAME B:

    p = (Y + X - 3) / 3
    q = Y - X + 3.

Then:

    n = p*q.

The level residue is simply:

    r_z = n mod 2^z.

The transition:

    z -> z+1

survives iff:

    2^z | X
    2^z | Y.

Equivalently, using zero-indexed binary bits:

    bit_z(X) = 0
    bit_z(Y) = 0.

When it survives:

    x_(z+1) = x_z / 2
    y_(z+1) = y_z / 2

and FRAME is unchanged.

Therefore the hierarchy is potentially reducible to:

    one frame bit
    +
    two invariant integers X,Y
    +
    their binary divisibility structure.

The predicted deepest level is:

    deepest =
        min(v2(X),v2(Y)) + 1.

The key result of this experiment is whether
( FRAME, X, Y ) alone reproduces the entire observed
level-by-level system.
"""
    )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print(
        "EXPERIMENT 603 START"
    )
    print("=" * 90)

    print()
    print(
        "COMPLETE GLOBAL RECONSTRUCTION "
        "FROM (FRAME, X, Y)"
    )

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

    print()
    print(
        "[3] BUILD OBSERVED LEVELS"
    )

    levels = build_levels(
        semiprimes
    )

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        print(
            f"    z={z:<2} "
            f"states={len(levels[z])}"
        )

    test_global_invariant(
        levels
    )

    test_coordinate_reconstruction(
        levels
    )

    test_factor_reconstruction(
        levels
    )

    test_n_reconstruction(
        levels
    )

    test_residue_reconstruction(
        levels
    )

    test_survival(
        levels
    )

    test_bit_survival(
        levels
    )

    test_deepest_level(
        levels
    )

    test_frame_invariance(
        levels
    )

    test_minimal_state(
        levels
    )

    test_full_reconstruction(
        levels
    )

    print_examples(
        levels
    )

    final_summary()

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 603 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()
