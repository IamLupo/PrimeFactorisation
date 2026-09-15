#!/usr/bin/env python3

# ==============================================================================
# EXPERIMENT 607
# ==============================================================================
# EXACT RESIDUE AUTOMATON WITH SPECIAL z=2 BASE CASE
#
# Previous experiment 606 found:
#
#   survival:
#       Alive_(z+1)
#         =
#       Alive_z & !xbit & !ybit
#
# exactly.
#
# The only residue-law failure was:
#
#   z=2 -> z=3
#
# with exactly 136377 failures.
#
# This experiment tests the exact piecewise law:
#
#   z = 2:
#
#       FRAME A:
#           r_3 = 7  if y even
#           r_3 = 3  if y odd
#
#       FRAME B:
#           r_3 = 5  if y even
#           r_3 = 1  if y odd
#
#   z >= 3:
#
#       FRAME A:
#           r_(z+1) =
#               -9 + 2^z*xbit   mod 2^(z+1)
#
#       FRAME B:
#           r_(z+1) =
#               -3 + 2^z*xbit   mod 2^(z+1)
#
# The experiment then asks whether:
#
#       (FRAME, z, xbit, ybit)
#
# is a complete residue-state description.
#
# It also tests whether the residue transition can be reduced further:
#
#       z=2:
#           FRAME + ybit
#
#       z>=3:
#           FRAME + xbit
#
# ==============================================================================


from __future__ import annotations

from dataclasses import dataclass
from collections import defaultdict, Counter
from math import isqrt


# ==============================================================================
# CONFIGURATION
# ==============================================================================

MAX_N = 2_000_000
MAX_Z = 24
EXAMPLES = 20


# ==============================================================================
# DATA STRUCTURES
# ==============================================================================

@dataclass(frozen=True)
class Pair:
    n: int
    p: int
    q: int


@dataclass(frozen=True)
class State:
    n: int
    p: int
    q: int
    frame: str
    X: int
    Y: int


# ==============================================================================
# PRIME SIEVE
# ==============================================================================

def sieve_primes(limit: int):

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

        sieve[
            start:
            limit + 1:
            p
        ] = b"\x00" * (
            ((limit - start) // p) + 1
        )

    return sieve


# ==============================================================================
# SEMIPRIME GENERATION
# ==============================================================================

def generate_semiprimes(
    limit: int,
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
                Pair(
                    n=p * q,
                    p=p,
                    q=q,
                )
            )

    return result


# ==============================================================================
# FRAME
# ==============================================================================

def frame_from_n(n: int):

    residue = n % 4

    if residue == 3:
        return "A"

    if residue == 1:
        return "B"

    raise ValueError(
        f"Unexpected n mod 4={residue}"
    )


# ==============================================================================
# GLOBAL INVARIANT COORDINATES
# ==============================================================================

def global_XY(
    p: int,
    q: int,
    frame: str,
):

    if frame == "A":

        X_num = q - p + 6
        Y_num = p + q

    elif frame == "B":

        X_num = 3 * p - q + 6
        Y_num = 3 * p + q

    else:
        raise ValueError(frame)

    assert X_num % 2 == 0
    assert Y_num % 2 == 0

    return (
        X_num // 2,
        Y_num // 2,
    )


# ==============================================================================
# V2
# ==============================================================================

def v2(value: int):

    if value == 0:
        return 10**9

    value = abs(value)

    count = 0

    while (
        value & 1
    ) == 0:

        count += 1
        value >>= 1

    return count


# ==============================================================================
# ABSOLUTE-VALUE BIT
# ==============================================================================

def bit_abs(
    value: int,
    position: int,
):

    return (
        abs(value) >> position
    ) & 1


# ==============================================================================
# ALIVE LEVEL
# ==============================================================================

def alive_at_level(
    X: int,
    Y: int,
    z: int,
):

    divisor = 1 << (
        z - 1
    )

    return (
        X % divisor == 0
        and
        Y % divisor == 0
    )


# ==============================================================================
# EXACT BASE CASE z=2 -> 3
# ==============================================================================
#
# For k=2:
#
#     A:
#       n = 4(y^2-x^2) + 12x - 9
#
#       mod 8:
#           n = 4y - 1
#             = 7 if y even
#             = 3 if y odd
#
#     B:
#       n = [same numerator] / 3
#
#       3^{-1} = 3 mod 8
#
#       n mod 8:
#           5 if y even
#           1 if y odd
#
# ==============================================================================

def residue_base_z2(
    frame: str,
    ybit: int,
):

    if frame == "A":

        return (
            7
            if ybit == 0
            else 3
        )

    if frame == "B":

        return (
            5
            if ybit == 0
            else 1
        )

    raise ValueError(frame)


# ==============================================================================
# EXACT z>=3 CHILD RESIDUE
# ==============================================================================

def residue_child_z_ge_3(
    frame: str,
    z: int,
    xbit: int,
):

    modulus = 1 << (
        z + 1
    )

    if frame == "A":

        return (
            -9
            + (1 << z) * xbit
        ) % modulus

    if frame == "B":

        return (
            -3
            + (1 << z) * xbit
        ) % modulus

    raise ValueError(frame)


# ==============================================================================
# PIECEWISE EXACT RESIDUE LAW
# ==============================================================================

def exact_child_residue(
    frame: str,
    z: int,
    xbit: int,
    ybit: int,
):

    if z == 2:

        return residue_base_z2(
            frame,
            ybit,
        )

    return residue_child_z_ge_3(
        frame,
        z,
        xbit,
    )


# ==============================================================================
# BUILD GLOBAL STATES
# ==============================================================================

def build_states(
    pairs,
):

    states = []

    for item in pairs:

        frame = frame_from_n(
            item.n
        )

        X, Y = global_XY(
            item.p,
            item.q,
            frame,
        )

        states.append(
            State(
                n=item.n,
                p=item.p,
                q=item.q,
                frame=frame,
                X=X,
                Y=Y,
            )
        )

    return states


# ==============================================================================
# TEST 1
# ==============================================================================

def test_global_states(
    states,
):

    print()
    print("=" * 90)
    print(
        "TEST 1: GLOBAL STATE CONSISTENCY"
    )
    print("=" * 90)

    failures = 0

    for state in states:

        if (
            state.frame
            != frame_from_n(state.n)
        ):
            failures += 1

        reconstructed_XY = global_XY(
            state.p,
            state.q,
            state.frame,
        )

        if reconstructed_XY != (
            state.X,
            state.Y,
        ):
            failures += 1

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 2
# ==============================================================================
# Exact piecewise residue law.
#
# IMPORTANT:
# Only evaluate the child residue for a LIVE parent level.
#
# This is a residue statement, not a survival statement.
# ==============================================================================

def test_piecewise_residue_law(
    states,
):

    print()
    print("=" * 90)
    print(
        "TEST 2: EXACT PIECEWISE CHILD RESIDUE LAW"
    )
    print("=" * 90)

    failures_by_z = Counter()
    checked_by_z = Counter()

    first_failures = []

    for state in states:

        for z in range(
            2,
            MAX_Z,
        ):

            if not alive_at_level(
                state.X,
                state.Y,
                z,
            ):
                continue

            xbit = bit_abs(
                state.X,
                z - 1,
            )

            ybit = bit_abs(
                state.Y,
                z - 1,
            )

            predicted = exact_child_residue(
                state.frame,
                z,
                xbit,
                ybit,
            )

            actual = (
                state.n
                % (
                    1 << (
                        z + 1
                    )
                )
            )

            checked_by_z[z] += 1

            if predicted != actual:

                failures_by_z[z] += 1

                if len(first_failures) < EXAMPLES:

                    first_failures.append(
                        (
                            state,
                            z,
                            xbit,
                            ybit,
                            predicted,
                            actual,
                        )
                    )

    for z in range(
        2,
        MAX_Z,
    ):

        print(
            f"z={z:2d} "
            f"checked={checked_by_z[z]} "
            f"failures={failures_by_z[z]}"
        )

    print(
        f"GLOBAL checked="
        f"{sum(checked_by_z.values())} "
        f"failures="
        f"{sum(failures_by_z.values())}"
    )

    if first_failures:

        print()
        print(
            "FIRST FAILURES"
        )

        for (
            state,
            z,
            xb,
            yb,
            pred,
            actual,
        ) in first_failures:

            print(
                f"    n={state.n} "
                f"frame={state.frame} "
                f"z={z} "
                f"bits=({xb},{yb}) "
                f"pred={pred} "
                f"actual={actual}"
            )


# ==============================================================================
# TEST 3
# ==============================================================================
# Determine the empirical transition function:
#
#     (frame, z, xbit, ybit) -> residue
#
# and prove that every observed state has one output.
# ==============================================================================

def test_full_state_residue_determinism(
    states,
):

    print()
    print("=" * 90)
    print(
        "TEST 3: FULL BOOLEAN STATE -> RESIDUE"
    )
    print("=" * 90)

    mapping = defaultdict(set)

    for state in states:

        for z in range(
            2,
            MAX_Z,
        ):

            if not alive_at_level(
                state.X,
                state.Y,
                z,
            ):
                continue

            xb = bit_abs(
                state.X,
                z - 1,
            )

            yb = bit_abs(
                state.Y,
                z - 1,
            )

            actual = (
                state.n
                % (
                    1 << (
                        z + 1
                    )
                )
            )

            mapping[
                (
                    state.frame,
                    z,
                    xb,
                    yb,
                )
            ].add(actual)

    ambiguous = [
        (
            key,
            values,
        )
        for key, values in mapping.items()
        if len(values) > 1
    ]

    print(
        f"states tested={len(mapping)}"
    )

    print(
        f"ambiguous states={len(ambiguous)}"
    )

    if ambiguous:

        print()
        print(
            "FIRST AMBIGUITIES"
        )

        for key, values in ambiguous[
            :EXAMPLES
        ]:

            print(
                f"    state={key} "
                f"residues={sorted(values)}"
            )


# ==============================================================================
# TEST 4
# ==============================================================================
# Can the z=2 output be reduced to FRAME + YBIT?
#
# Can z>=3 output be reduced to FRAME + XBIT?
# ==============================================================================

def test_minimal_residue_state(
    states,
):

    print()
    print("=" * 90)
    print(
        "TEST 4: MINIMAL RESIDUE STATE"
    )
    print("=" * 90)

    base_map = defaultdict(set)
    higher_map = defaultdict(set)

    for state in states:

        # z=2
        if alive_at_level(
            state.X,
            state.Y,
            2,
        ):

            ybit = bit_abs(
                state.Y,
                1,
            )

            residue = (
                state.n % 8
            )

            base_map[
                (
                    state.frame,
                    ybit,
                )
            ].add(residue)

        # z >= 3
        for z in range(
            3,
            MAX_Z,
        ):

            if not alive_at_level(
                state.X,
                state.Y,
                z,
            ):
                continue

            xbit = bit_abs(
                state.X,
                z - 1,
            )

            residue = (
                state.n
                % (
                    1 << (
                        z + 1
                    )
                )
            )

            higher_map[
                (
                    state.frame,
                    z,
                    xbit,
                )
            ].add(residue)

    base_ambiguous = [
        x
        for x in base_map.items()
        if len(x[1]) > 1
    ]

    higher_ambiguous = [
        x
        for x in higher_map.items()
        if len(x[1]) > 1
    ]

    print()
    print(
        "z=2:"
    )

    print(
        f"    states={len(base_map)}"
    )

    print(
        f"    ambiguous={len(base_ambiguous)}"
    )

    for key, values in sorted(
        base_map.items()
    ):

        print(
            f"    {key} -> "
            f"{sorted(values)}"
        )

    print()
    print(
        "z>=3:"
    )

    print(
        f"    states={len(higher_map)}"
    )

    print(
        f"    ambiguous={len(higher_ambiguous)}"
    )

    for key, values in sorted(
        higher_map.items()
    ):

        print(
            f"    {key} -> "
            f"{sorted(values)}"
        )


# ==============================================================================
# TEST 5
# ==============================================================================
# Survival remains independent of the residue law.
#
# Alive(z+1) iff:
#
#     current x bit = 0
#     AND
#     current y bit = 0.
# ==============================================================================

def test_survival(
    states,
):

    print()
    print("=" * 90)
    print(
        "TEST 5: SURVIVAL AUTOMATON"
    )
    print("=" * 90)

    failures = 0
    checked = 0

    for state in states:

        for z in range(
            2,
            MAX_Z,
        ):

            alive_z = alive_at_level(
                state.X,
                state.Y,
                z,
            )

            if not alive_z:
                continue

            xb = bit_abs(
                state.X,
                z - 1,
            )

            yb = bit_abs(
                state.Y,
                z - 1,
            )

            predicted = (
                xb == 0
                and
                yb == 0
            )

            actual = alive_at_level(
                state.X,
                state.Y,
                z + 1,
            )

            checked += 1

            if predicted != actual:

                failures += 1

    print(
        f"checked={checked} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 6
# ==============================================================================
# Verify the algebra directly.
#
# This makes the z=2 exception explicit instead of merely empirical.
# ==============================================================================

def test_symbolic_modular_cases():

    print()
    print("=" * 90)
    print(
        "TEST 6: SYMBOLIC MODULAR CASES"
    )
    print("=" * 90)

    print(
        """
FRAME A:

    n = (k(y-x)+3)(k(y+x)-3)

      = k^2(y^2-x^2)
        + 6kx
        - 9

FRAME B:

    n = [(k(y-x)+3)(k(y+x)-3)] / 3

      = [k^2(y^2-x^2)
        + 6kx
        - 9] / 3

For z >= 3:

    k = 2^(z-1)

    k^2 = 2^(2z-2)

    2z-2 >= z+1

    exactly when z >= 3.

Therefore:

    k^2(y^2-x^2)
        = 0 mod 2^(z+1)

and:

FRAME A:

    n
      = -9 + 3*2^z*x
      = -9 + 2^z*x     mod 2^(z+1)

FRAME B:

    n
      = -3 + 2^z*x     mod 2^(z+1).

At z=2:

    k=2

    k^2=4

    4 is NOT divisible by 8.

So the quadratic term survives and the residue depends on y parity:

FRAME A:

    n mod 8 = 7 - 4*ybit

    equivalently:
        ybit=0 -> 7
        ybit=1 -> 3

FRAME B:

    n mod 8:
        ybit=0 -> 5
        ybit=1 -> 1

This is why Experiment 606 had exactly the z=2 residue failures.
"""
    )


# ==============================================================================
# TEST 7
# ==============================================================================
# Direct SAT-style residue clauses.
#
# For z>=3, once F and Xbit are known, residue is fixed.
#
# For z=2, F and Ybit are sufficient.
# ==============================================================================

def print_sat_residue_model():

    print()
    print("=" * 90)
    print(
        "TEST 7: SAT RESIDUE MODEL"
    )
    print("=" * 90)

    print(
        """
SURVIVAL:

    S_z <-> (!X_z & !Y_z)

    CNF:

        (!S_z v !X_z)
        (!S_z v !Y_z)
        ( S_z v  X_z v  Y_z)


BASE TRANSITION z=2:

    inputs:
        F
        Y_1

    outputs:

        F=A:
            Y_1=0 -> r_3=7
            Y_1=1 -> r_3=3

        F=B:
            Y_1=0 -> r_3=5
            Y_1=1 -> r_3=1


ALL HIGHER TRANSITIONS z>=3:

    inputs:
        F
        X_(z-1)

    FRAME A:

        X=0 -> r'=-9 mod 2^(z+1)
        X=1 -> r'=2^z-9 mod 2^(z+1)

    FRAME B:

        X=0 -> r'=-3 mod 2^(z+1)
        X=1 -> r'=2^z-3 mod 2^(z+1)

The residue therefore does NOT require the full X,Y values.

The information split is:

    SURVIVAL:
        Xbit + Ybit

    BASE RESIDUE:
        FRAME + Ybit

    HIGHER RESIDUE:
        FRAME + Xbit

    FACTOR RECOVERY:
        FRAME + X + Y
"""
    )


# ==============================================================================
# TEST 8
# ==============================================================================
# Check that the active residue itself is exactly the survivor branch.
# ==============================================================================

def test_active_survivor_residue(
    states,
):

    print()
    print("=" * 90)
    print(
        "TEST 8: SURVIVING CHILD -> ACTIVE RESIDUE"
    )
    print("=" * 90)

    checked = 0
    failures = 0

    active_counts = defaultdict(
        Counter
    )

    for state in states:

        for z in range(
            2,
            MAX_Z,
        ):

            if not alive_at_level(
                state.X,
                state.Y,
                z,
            ):
                continue

            if not alive_at_level(
                state.X,
                state.Y,
                z + 1,
            ):
                continue

            child_res = (
                state.n
                % (
                    1 << (
                        z + 1
                    )
                )
            )

            expected = (
                (-9)
                if state.frame == "A"
                else (-3)
            ) % (
                1 << (
                    z + 1
                )
            )

            active_counts[
                (
                    z,
                    state.frame,
                )
            ][child_res] += 1

            checked += 1

            if child_res != expected:

                failures += 1

                if failures <= EXAMPLES:

                    print(
                        f"    mismatch "
                        f"n={state.n} "
                        f"z={z} "
                        f"frame={state.frame} "
                        f"actual={child_res} "
                        f"expected={expected}"
                    )

    print(
        f"surviving edges checked={checked} "
        f"failures={failures}"
    )

    print()
    print(
        "ACTIVE CHILD RESIDUES"
    )

    for key in sorted(
        active_counts
    ):

        print(
            f"    {key} -> "
            f"{dict(active_counts[key])}"
        )


# ==============================================================================
# TEST 9
# ==============================================================================
# COMPLETE GLOBAL REPLAY
#
# The automaton must now reproduce:
#
#   coordinates
#   residue
#   survival
#
# through the first dead edge.
# ==============================================================================

def test_complete_replay(
    states,
):

    print()
    print("=" * 90)
    print(
        "TEST 9: COMPLETE AUTOMATON REPLAY"
    )
    print("=" * 90)

    failures = 0
    checked_levels = 0

    for state in states:

        for z in range(
            2,
            MAX_Z + 1,
        ):

            alive = alive_at_level(
                state.X,
                state.Y,
                z,
            )

            if not alive:
                break

            divisor = 1 << (
                z - 1
            )

            x = state.X // divisor
            y = state.Y // divisor

            xb = x & 1
            yb = y & 1

            actual_residue = (
                state.n
                % (
                    1 << z
                )
            )

            expected_active = (
                (-9)
                if state.frame == "A"
                else (-3)
            ) % (
                1 << z
            )

            # Current level must always lie on its frame's active
            # residue while the representation exists.
            if actual_residue != expected_active:

                failures += 1

                if failures <= EXAMPLES:

                    print(
                        f"    residue mismatch "
                        f"n={state.n} "
                        f"z={z} "
                        f"actual={actual_residue} "
                        f"expected={expected_active}"
                    )

            expected_next = (
                xb == 0
                and
                yb == 0
            )

            actual_next = alive_at_level(
                state.X,
                state.Y,
                z + 1,
            )

            if z < MAX_Z:

                if expected_next != actual_next:

                    failures += 1

                    if failures <= EXAMPLES:

                        print(
                            f"    survival mismatch "
                            f"n={state.n} "
                            f"z={z}"
                        )

            checked_levels += 1

    print(
        f"checked levels={checked_levels} "
        f"failures={failures}"
    )


# ==============================================================================
# EXAMPLES
# ==============================================================================

def show_examples(
    states,
):

    print()
    print("=" * 90)
    print(
        "EXAMPLES"
    )
    print("=" * 90)

    shown = 0

    for state in states:

        if shown >= EXAMPLES:
            break

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
            f"v2Y={v2(state.Y)}"
        )

        for z in range(
            2,
            MAX_Z + 1,
        ):

            if not alive_at_level(
                state.X,
                state.Y,
                z,
            ):
                break

            divisor = 1 << (
                z - 1
            )

            x = state.X // divisor
            y = state.Y // divisor

            xb = x & 1
            yb = y & 1

            residue = (
                state.n
                % (
                    1 << z
                )
            )

            survives = (
                xb == 0
                and
                yb == 0
            )

            if z == 2:

                child_prediction = (
                    residue_base_z2(
                        state.frame,
                        yb,
                    )
                )

            else:

                child_prediction = (
                    residue_child_z_ge_3(
                        state.frame,
                        z,
                        xb,
                    )
                )

            print(
                f"    z={z:<2} "
                f"x={x:<8} "
                f"y={y:<8} "
                f"bits=({xb},{yb}) "
                f"r={residue:<6} "
                f"next_r={child_prediction:<6} "
                f"survive={survives}"
            )

        shown += 1


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print(
        "EXPERIMENT 607 START"
    )
    print("=" * 90)

    print()
    print(
        "[1] PRIME SIEVE"
    )

    sieve = sieve_primes(
        MAX_N
    )

    print()
    print(
        "[2] SEMIPRIME GENERATION"
    )

    pairs = generate_semiprimes(
        MAX_N,
        sieve,
    )

    print(
        f"    semiprimes={len(pairs)}"
    )

    print()
    print(
        "[3] BUILD GLOBAL STATES"
    )

    states = build_states(
        pairs
    )

    print(
        f"    global states={len(states)}"
    )

    test_global_states(
        states
    )

    test_piecewise_residue_law(
        states
    )

    test_full_state_residue_determinism(
        states
    )

    test_minimal_residue_state(
        states
    )

    test_survival(
        states
    )

    test_symbolic_modular_cases()

    print_sat_residue_model()

    test_active_survivor_residue(
        states
    )

    test_complete_replay(
        states
    )

    show_examples(
        states
    )

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 607 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()
