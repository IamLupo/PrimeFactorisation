#!/usr/bin/env python3

# ==============================================================================
# EXPERIMENT 604
# ==============================================================================
# CORRECT 2-ADIC BIT INDEX + DIRECT RESIDUE LAW + SAT STATE COLLAPSE
#
# No external files.
# No external sources.
# No output files.
#
# Core invariant:
#
#     X = 2^(z-1) * x_z
#     Y = 2^(z-1) * y_z
#
# Therefore:
#
#     x_z = X >> (z-1)
#     y_z = Y >> (z-1)
#
# Transition z -> z+1 survives iff:
#
#     x_z even
#     y_z even
#
# i.e.
#
#     bit_(z-1)(X) = 0
#     bit_(z-1)(Y) = 0
#
# IMPORTANT:
#     This is z-1, not z.
#
# The experiment tests:
#
#   1. Correct bit-index survival.
#   2. Exact deepest-level formula.
#   3. Direct residue formulas from X,Y.
#   4. Whether FRAME is irrelevant to survival.
#   5. Minimal SAT transition state.
#   6. Complete reconstruction from (FRAME,X,Y).
#   7. Closed-form residue sequence for each frame.
# ==============================================================================


from __future__ import annotations

from dataclasses import dataclass
from collections import Counter
from math import isqrt


# ==============================================================================
# CONFIG
# ==============================================================================

MAX_N = 2_000_000
MIN_Z = 2
MAX_Z = 16

SHOW = 20


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
# SEMIPRIMES
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
# FRAME COORDINATES
# ==============================================================================

def frame_A_xy(
    p,
    q,
    z,
):

    k = 1 << (
        z - 1
    )

    d = 2 * k

    x_num = (
        q - p + 6
    )

    y_num = (
        p + q
    )

    if x_num % d:
        return None

    if y_num % d:
        return None

    return (
        x_num // d,
        y_num // d,
    )


def frame_B_xy(
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

    if x_num % d:
        return None

    if y_num % d:
        return None

    return (
        x_num // d,
        y_num // d,
    )


# ==============================================================================
# LEVEL BUILD
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

            a = frame_A_xy(
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

            b = frame_B_xy(
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
# BIT HELPERS
# ==============================================================================

def bit(
    value,
    position,
):

    return (
        abs(value) >> position
    ) & 1


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
# DIRECT SURVIVAL
# ==============================================================================

def direct_survival(
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
# CORRECT BIT SURVIVAL
#
# Transition z -> z+1 uses bit z-1.
# ==============================================================================

def bit_survival(
    X,
    Y,
    z,
):

    position = z - 1

    return (
        bit(X, position) == 0
        and
        bit(Y, position) == 0
    )


# ==============================================================================
# FACTOR RECONSTRUCTION
# ==============================================================================

def factors_from_XY(
    frame,
    X,
    Y,
):

    if frame == "A":

        p = (
            Y - X + 3
        )

        q = (
            Y + X - 3
        )

        return (
            p,
            q,
        )

    if frame == "B":

        numerator = (
            Y + X - 3
        )

        if numerator % 3:
            return None

        p = (
            numerator // 3
        )

        q = (
            Y - X + 3
        )

        return (
            p,
            q,
        )

    return None


# ==============================================================================
# DIRECT N
# ==============================================================================

def n_from_XY(
    frame,
    X,
    Y,
):

    factors = factors_from_XY(
        frame,
        X,
        Y,
    )

    if factors is None:
        return None

    p, q = factors

    return p * q


# ==============================================================================
# DIRECT RESIDUE FROM X,Y
#
# We deliberately calculate n from the invariant coordinates,
# then reduce modulo 2^z.
#
# A second symbolic derivation is also tested later.
# ==============================================================================

def direct_residue(
    frame,
    X,
    Y,
    z,
):

    n = n_from_XY(
        frame,
        X,
        Y,
    )

    if n is None:
        return None

    return n % (
        1 << z
    )


# ==============================================================================
# SYMBOLIC RESIDUE LAW
#
# Let:
#
#     k = 2^(z-1).
#
# FRAME A:
#
#     n = (k(y-x)+3)(k(y+x)-3)
#
#       = k^2(y^2-x^2)
#         + 6kx
#         - 9
#
# Since k^2 is divisible by 2^(z+1):
#
#     n mod 2^z = -9 mod 2^z.
#
# FRAME B:
#
#     n = [(k(y-x)+3)(k(y+x)-3)] / 3
#
#     mod 2^z:
#
#         -9 / 3
#         = -3
#
# Therefore for the ACTIVE level itself:
#
#     A: residue = -9 mod 2^z
#     B: residue = -3 mod 2^z
#
# for z >= 2.
#
# This experiment verifies this directly from X,Y.
# ==============================================================================

def symbolic_active_residue(
    frame,
    z,
):

    modulus = 1 << z

    if frame == "A":
        return (-9) % modulus

    if frame == "B":
        return (-3) % modulus

    return None


# ==============================================================================
# CHILD RESIDUE FROM THE NEXT BIT
#
# Transition z -> z+1:
#
# FRAME A:
#
#     r' = -9 + 2^z * x_z
#        mod 2^(z+1)
#
# FRAME B:
#
#     r' = -3 + 2^z * x_z
#        mod 2^(z+1)
#
# Only x_z mod 2 matters.
# ==============================================================================

def symbolic_child_residue(
    frame,
    X,
    z,
):

    xbit = bit(
        X,
        z - 1,
    )

    modulus = 1 << (
        z + 1
    )

    base = (
        -9
        if frame == "A"
        else -3
    )

    return (
        base
        + (1 << z) * xbit
    ) % modulus


# ==============================================================================
# TEST 1
# ==============================================================================

def test_correct_bit_index(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 1: CORRECTED BIT INDEX"
    )
    print("=" * 90)

    print(
        """
For z -> z+1:

    x_z = X / 2^(z-1)
    y_z = Y / 2^(z-1)

Therefore:

    x_z even <=> bit_(z-1)(X)=0
    y_z even <=> bit_(z-1)(Y)=0.
"""
    )

    failures = 0
    checked = 0

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        local = 0

        for state in levels[z].values():

            actual = (
                n_exists(
                    levels,
                    z + 1,
                    state.n,
                )
            )

            predicted = bit_survival(
                state.X,
                state.Y,
                z,
            )

            checked += 1

            if actual != predicted:

                failures += 1
                local += 1

                if local <= SHOW:

                    print(
                        f"    mismatch "
                        f"z={z} "
                        f"n={state.n} "
                        f"X={state.X} "
                        f"Y={state.Y} "
                        f"bit={z-1} "
                        f"Xbit={bit(state.X,z-1)} "
                        f"Ybit={bit(state.Y,z-1)} "
                        f"pred={predicted} "
                        f"actual={actual}"
                    )

        print(
            f"z={z}->{z+1} "
            f"failures={local}"
        )

    print()
    print(
        f"GLOBAL checked={checked} "
        f"failures={failures}"
    )


# ==============================================================================
# EXISTENCE HELPER
# ==============================================================================

def n_exists(
    levels,
    z,
    n,
):

    return (
        n in levels.get(
            z,
            {}
        )
    )


# ==============================================================================
# TEST 2
# ==============================================================================

def test_survival_equivalence(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 2: DIVISIBILITY == BITWISE SURVIVAL"
    )
    print("=" * 90)

    checked = 0
    failures = 0

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        local = 0

        for state in levels[z].values():

            a = direct_survival(
                state.X,
                state.Y,
                z,
            )

            b = bit_survival(
                state.X,
                state.Y,
                z,
            )

            checked += 1

            if a != b:

                failures += 1
                local += 1

        print(
            f"z={z}->{z+1} "
            f"failures={local}"
        )

    print()
    print(
        f"GLOBAL checked={checked} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 3
# ==============================================================================

def test_deepest(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 3: EXACT DEEPEST LEVEL"
    )
    print("=" * 90)

    checked = 0
    failures = 0

    depth_counter = Counter()

    for n, base in levels[
        MIN_Z
    ].items():

        actual = MIN_Z - 1

        for z in range(
            MIN_Z,
            MAX_Z + 1,
        ):

            if n_exists(
                levels,
                z,
                n,
            ):

                actual = z

        predicted = (
            min(
                v2(base.X),
                v2(base.Y),
            )
            + 1
        )

        # We only care about states where the predicted
        # depth falls inside the measured range.
        if predicted <= MAX_Z:

            checked += 1

            if predicted != actual:

                failures += 1

                if failures <= SHOW:

                    print(
                        f"    mismatch "
                        f"n={n} "
                        f"X={base.X} "
                        f"Y={base.Y} "
                        f"v2X={v2(base.X)} "
                        f"v2Y={v2(base.Y)} "
                        f"pred={predicted} "
                        f"actual={actual}"
                    )

        depth_counter[
            actual
        ] += 1

    print()
    print(
        f"bounded checked={checked} "
        f"failures={failures}"
    )

    print()
    print(
        "DEPTH DISTRIBUTION"
    )

    for depth, count in sorted(
        depth_counter.items()
    ):

        print(
            f"    depth={depth}: "
            f"{count}"
        )


# ==============================================================================
# TEST 4
# ==============================================================================

def test_active_residue(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 4: CLOSED-FORM ACTIVE RESIDUE"
    )
    print("=" * 90)

    checked = 0
    failures = 0

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        local = 0

        residues = {}

        for state in levels[z].values():

            predicted = symbolic_active_residue(
                state.frame,
                z,
            )

            residues.setdefault(
                state.frame,
                set()
            ).add(
                state.residue
            )

            checked += 1

            if predicted != state.residue:

                failures += 1
                local += 1

        print(
            f"z={z} "
            f"failures={local} "
            f"observed={residues}"
        )

    print()
    print(
        f"GLOBAL checked={checked} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 5
# ==============================================================================

def test_child_residue(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 5: CLOSED-FORM CHILD RESIDUE"
    )
    print("=" * 90)

    checked = 0
    failures = 0

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        local = 0

        for state in levels[z].values():

            expected = symbolic_child_residue(
                state.frame,
                state.X,
                z,
            )

            actual_child = levels[
                z + 1
            ].get(
                state.n
            )

            if actual_child is not None:

                actual = (
                    actual_child.residue
                )

                checked += 1

                if expected != actual:

                    failures += 1
                    local += 1

                    if local <= SHOW:

                        print(
                            f"    mismatch "
                            f"z={z} "
                            f"n={state.n} "
                            f"frame={state.frame} "
                            f"X={state.X} "
                            f"xbit={bit(state.X,z-1)} "
                            f"pred={expected} "
                            f"actual={actual}"
                        )

        print(
            f"z={z}->{z+1} "
            f"surviving checked={checked} "
            f"local_failures={local}"
        )

    print()
    print(
        f"GLOBAL surviving checked={checked} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 6
# ==============================================================================

def test_frame_independence_of_survival(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 6: DOES SURVIVAL NEED FRAME?"
    )
    print("=" * 90)

    outcomes = {}

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        groups = {}

        for state in levels[z].values():

            key = (
                bit(
                    state.X,
                    z - 1,
                ),
                bit(
                    state.Y,
                    z - 1,
                ),
            )

            survives_now = n_exists(
                levels,
                z + 1,
                state.n,
            )

            groups.setdefault(
                key,
                set()
            ).add(
                survives_now
            )

        deterministic = all(
            len(values) == 1
            for values in groups.values()
        )

        outcomes[z] = deterministic

        print(
            f"z={z} "
            f"Xbit+Ybit deterministic="
            f"{deterministic}"
        )

        for key in sorted(
            groups
        ):

            print(
                f"    {key} -> "
                f"{sorted(groups[key])}"
            )

    print()
    print(
        "GLOBAL:"
    )
    print(
        f"    deterministic="
        f"{all(outcomes.values())}"
    )


# ==============================================================================
# TEST 7
# ==============================================================================

def test_sat_clauses(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 7: SAT CLAUSES"
    )
    print("=" * 90)

    print(
        """
Let:

    X = bit_(z-1)(X)
    Y = bit_(z-1)(Y)
    S = child survives

Exact survival law:

    S <-> (!X & !Y)

CNF:

    (!S v !X)
    (!S v !Y)
    ( S v  X v  Y)

The frame does not occur in the survival clauses.
"""
    )

    failures = 0
    checked = 0

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        for state in levels[z].values():

            Xb = bit(
                state.X,
                z - 1,
            )

            Yb = bit(
                state.Y,
                z - 1,
            )

            S = int(
                n_exists(
                    levels,
                    z + 1,
                    state.n,
                )
            )

            clause1 = (
                (not S)
                or
                (not Xb)
            )

            clause2 = (
                (not S)
                or
                (not Yb)
            )

            clause3 = (
                S
                or
                Xb
                or
                Yb
            )

            checked += 1

            if not (
                clause1
                and
                clause2
                and
                clause3
            ):

                failures += 1

    print(
        f"checked={checked} "
        f"failures={failures}"
    )

    print()
    print(
        "CLAUSES:"
    )

    print(
        "    (!S v !X)"
    )

    print(
        "    (!S v !Y)"
    )

    print(
        "    ( S v  X v  Y)"
    )


# ==============================================================================
# TEST 8
# ==============================================================================

def test_frame_preservation(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 8: FRAME PRESERVATION"
    )
    print("=" * 90)

    checked = 0
    failures = 0

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        local = 0

        for state in levels[z].values():

            child_state = levels[
                z + 1
            ].get(
                state.n
            )

            if child_state is None:
                continue

            checked += 1

            if (
                child_state.frame
                !=
                state.frame
            ):

                failures += 1
                local += 1

                if local <= SHOW:

                    print(
                        f"    mismatch "
                        f"z={z} "
                        f"n={state.n} "
                        f"{state.frame}"
                        f"->{child_state.frame}"
                    )

        print(
            f"z={z}->{z+1} "
            f"failures={local}"
        )

    print()
    print(
        f"GLOBAL checked={checked} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 9
# ==============================================================================

def test_complete_state(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 9: COMPLETE STATE RECONSTRUCTION"
    )
    print("=" * 90)

    checked = 0
    failures = 0

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        local = 0

        for state in levels[z].values():

            checked += 1

            n2 = n_from_XY(
                state.frame,
                state.X,
                state.Y,
            )

            xy = (
                state.X // (
                    1 << (
                        z - 1
                    )
                ),
                state.Y // (
                    1 << (
                        z - 1
                    )
                ),
            )

            residue = (
                n2 % (
                    1 << z
                )
            )

            if n2 != state.n:

                failures += 1
                local += 1
                continue

            if xy != (
                state.x,
                state.y,
            ):

                failures += 1
                local += 1
                continue

            if residue != state.residue:

                failures += 1
                local += 1
                continue

        print(
            f"z={z} "
            f"failures={local}"
        )

    print()
    print(
        f"GLOBAL checked={checked} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 10
# ==============================================================================

def print_sat_transition(
):

    print()
    print("=" * 90)
    print(
        "TEST 10: MINIMAL TRANSITION"
    )
    print("=" * 90)

    print(
        """
STATE AT LEVEL z:

    F = frame
    X = bit_(z-1)(X)
    Y = bit_(z-1)(Y)

TRANSITION:

    S = !X & !Y

IF S:

    F' = F
    X' = next invariant bit
    Y' = next invariant bit

COORDINATES:

    x_(z+1) = x_z >> 1
    y_(z+1) = y_z >> 1.

Thus the actual SAT state does not need the full x_z,y_z values.
The transition only consumes one binary digit from X and Y.
"""
    )

    print(
        "    S <-> (!X & !Y)"
    )

    print(
        "    F' <-> F  when S"
    )


# ==============================================================================
# EXAMPLES
# ==============================================================================

def examples(
    levels,
):

    print()
    print("=" * 90)
    print(
        "EXAMPLES"
    )
    print("=" * 90)

    shown = 0

    for state in levels[
        MIN_Z
    ].values():

        if shown >= SHOW:
            break

        deepest = (
            min(
                v2(state.X),
                v2(state.Y),
            )
            + 1
        )

        print()
        print(
            f"n={state.n} "
            f"p={state.p} "
            f"q={state.q}"
        )

        print(
            f"    frame={state.frame} "
            f"X={state.X} "
            f"Y={state.Y}"
        )

        print(
            f"    v2X={v2(state.X)} "
            f"v2Y={v2(state.Y)} "
            f"deepest={deepest}"
        )

        print(
            "    descent:"
        )

        for z in range(
            MIN_Z,
            min(
                MAX_Z,
                deepest,
            ) + 1,
        ):

            xpos = z - 1

            Xbit = bit(
                state.X,
                xpos,
            )

            Ybit = bit(
                state.Y,
                xpos,
            )

            xy = (
                state.X >> xpos,
                state.Y >> xpos,
            )

            residue = (
                state.n
                % (
                    1 << z
                )
            )

            child_exists = (
                n_exists(
                    levels,
                    z + 1,
                    state.n,
                )
            )

            print(
                f"        z={z:<2} "
                f"x={xy[0]:<8} "
                f"y={xy[1]:<8} "
                f"bit=({Xbit},{Ybit}) "
                f"residue={residue:<6} "
                f"next={child_exists}"
            )

        shown += 1


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print(
        "EXPERIMENT 604 START"
    )
    print("=" * 90)

    print()
    print(
        "CORRECT 2-ADIC BIT INDEX + "
        "DIRECT RESIDUE LAW + SAT COLLAPSE"
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
        "[3] BUILD LEVELS"
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

    test_correct_bit_index(
        levels
    )

    test_survival_equivalence(
        levels
    )

    test_deepest(
        levels
    )

    test_active_residue(
        levels
    )

    test_child_residue(
        levels
    )

    test_frame_independence_of_survival(
        levels
    )

    test_sat_clauses(
        levels
    )

    test_frame_preservation(
        levels
    )

    test_complete_state(
        levels
    )

    print_sat_transition()

    examples(
        levels
    )

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 604 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()
