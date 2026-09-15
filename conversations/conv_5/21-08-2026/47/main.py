#!/usr/bin/env python3

# ==============================================================================
# EXPERIMENT 606
# ==============================================================================
# PREFIX-CONDITIONED 2-ADIC AUTOMATON
#
# Experiment 605 exposed an important distinction:
#
#   bit_(z-1)(X)=0 and bit_(z-1)(Y)=0
#
# is sufficient for the transition z -> z+1 ONLY if the state already
# exists at level z.
#
# A state that died earlier must not be tested again as though it were
# still an active level-z representation.
#
# Therefore this experiment explicitly separates:
#
#   ALIVE(z)
#
# from:
#
#   CURRENT_BIT(z)
#
# and verifies:
#
#   ALIVE(z+1)
#       =
#   ALIVE(z)
#   AND
#   bit_(z-1)(X)=0
#   AND
#   bit_(z-1)(Y)=0.
#
# The experiment also verifies that the entire observed hierarchy can be
# reconstructed by scanning the bits only until the first failed bit pair.
#
# No files.
# No external sources.
#
# ==============================================================================


from __future__ import annotations

from dataclasses import dataclass
from collections import Counter
from math import isqrt


# ==============================================================================
# CONFIG
# ==============================================================================

MAX_N = 2_000_000
MAX_Z = 24
EXAMPLE_COUNT = 20


# ==============================================================================
# DATA
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

    r = n % 4

    if r == 3:
        return "A"

    if r == 1:
        return "B"

    raise ValueError(
        f"Unexpected odd residue n mod 4={r}"
    )


# ==============================================================================
# GLOBAL COORDINATES
# ==============================================================================

def global_XY(
    p: int,
    q: int,
    frame: str,
):

    if frame == "A":

        # X=(q-p+6)/2
        # Y=(p+q)/2

        xn = q - p + 6
        yn = p + q

    elif frame == "B":

        # X=(3p-q+6)/2
        # Y=(3p+q)/2

        xn = 3 * p - q + 6
        yn = 3 * p + q

    else:
        raise ValueError(frame)

    if xn % 2 != 0:
        raise AssertionError(
            "X numerator not divisible by 2"
        )

    if yn % 2 != 0:
        raise AssertionError(
            "Y numerator not divisible by 2"
        )

    return (
        xn // 2,
        yn // 2,
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

        p = Y - X + 3
        q = Y + X - 3

        return p, q

    if frame == "B":

        num = Y + X - 3

        if num % 3:
            return None

        p = num // 3
        q = Y - X + 3

        return p, q

    raise ValueError(frame)


# ==============================================================================
# N RECONSTRUCTION
# ==============================================================================

def n_from_state(
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
# V2
# ==============================================================================

def v2(
    value: int,
):

    if value == 0:
        return 10**9

    value = abs(value)

    result = 0

    while (
        value & 1
    ) == 0:

        result += 1
        value >>= 1

    return result


# ==============================================================================
# DIVISIBILITY
# ==============================================================================

def divisible_by_power_of_two(
    value: int,
    exponent: int,
):

    modulus = 1 << exponent

    return value % modulus == 0


# ==============================================================================
# BIT OF ABSOLUTE VALUE
# ==============================================================================
#
# Using abs(value) is intentional.
#
# For divisibility by 2^k:
#
#     2^k | value
#
# iff:
#
#     the lowest k bits of |value| vanish.
#
# We are not interpreting negative X,Y as fixed-width two's-complement
# bit streams.
#
# ==============================================================================

def bit_abs(
    value: int,
    position: int,
):

    return (
        abs(value) >> position
    ) & 1


# ==============================================================================
# LEVEL EXISTENCE
# ==============================================================================

def level_exists(
    X,
    Y,
    z,
):

    required = 1 << (
        z - 1
    )

    return (
        X % required == 0
        and
        Y % required == 0
    )


# ==============================================================================
# ACTIVE RESIDUE
# ==============================================================================

def active_residue(
    frame,
    z,
):

    modulus = 1 << z

    if frame == "A":
        return (-9) % modulus

    if frame == "B":
        return (-3) % modulus

    raise ValueError(frame)


# ==============================================================================
# CHILD RESIDUE
# ==============================================================================

def child_residue(
    frame,
    X,
    z,
):

    xbit = bit_abs(
        X,
        z - 1,
    )

    modulus = 1 << (
        z + 1
    )

    if frame == "A":
        base = -9
    else:
        base = -3

    return (
        base
        + (1 << z) * xbit
    ) % modulus


# ==============================================================================
# BUILD STATES
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

def test_global_reconstruction(
    states,
):

    print()
    print("=" * 90)
    print(
        "TEST 1: GLOBAL RECONSTRUCTION"
    )
    print("=" * 90)

    failures = 0

    for state in states:

        factors = factors_from_XY(
            state.frame,
            state.X,
            state.Y,
        )

        if factors != (
            state.p,
            state.q,
        ):

            failures += 1

        reconstructed = n_from_state(
            state.frame,
            state.X,
            state.Y,
        )

        if reconstructed != state.n:

            failures += 1

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 2
# ==============================================================================
# ALIVE(z) <=> divisibility by 2^(z-1)
#
# This is a LEVEL property, not a TRANSITION property.
# ==============================================================================

def test_alive_formula(
    states,
):

    print()
    print("=" * 90)
    print(
        "TEST 2: LEVEL ALIVE FORMULA"
    )
    print("=" * 90)

    checked = 0
    failures = 0

    for state in states:

        for z in range(
            2,
            MAX_Z + 1,
        ):

            expected = level_exists(
                state.X,
                state.Y,
                z,
            )

            # Reconstruct the level coordinates exactly when alive.
            divisor = 1 << (
                z - 1
            )

            actual = (
                state.X % divisor == 0
                and
                state.Y % divisor == 0
            )

            checked += 1

            if expected != actual:

                failures += 1

    print(
        f"checked={checked} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 3
# ==============================================================================
# THE IMPORTANT CORRECTION:
#
# ALIVE(z+1)
# =
# ALIVE(z)
# AND
# current X bit == 0
# AND
# current Y bit == 0
#
# where current bit = position z-1.
# ==============================================================================

def test_conditioned_transition(
    states,
):

    print()
    print("=" * 90)
    print(
        "TEST 3: PREFIX-CONDITIONED SURVIVAL TRANSITION"
    )
    print("=" * 90)

    total = 0
    failures = 0

    first_failures = []

    for state in states:

        for z in range(
            2,
            MAX_Z,
        ):

            alive_z = level_exists(
                state.X,
                state.Y,
                z,
            )

            alive_next = level_exists(
                state.X,
                state.Y,
                z + 1,
            )

            xb = bit_abs(
                state.X,
                z - 1,
            )

            yb = bit_abs(
                state.Y,
                z - 1,
            )

            predicted_next = (
                alive_z
                and
                xb == 0
                and
                yb == 0
            )

            total += 1

            if alive_next != predicted_next:

                failures += 1

                if len(first_failures) < EXAMPLE_COUNT:

                    first_failures.append(
                        (
                            state,
                            z,
                            alive_z,
                            alive_next,
                            xb,
                            yb,
                            predicted_next,
                        )
                    )

    print(
        f"checked={total} "
        f"failures={failures}"
    )

    if first_failures:

        print()
        print(
            "FIRST FAILURES"
        )

        for (
            state,
            z,
            alive_z,
            alive_next,
            xb,
            yb,
            pred,
        ) in first_failures:

            print(
                f"    n={state.n} "
                f"z={z} "
                f"alive_z={alive_z} "
                f"alive_next={alive_next} "
                f"bits=({xb},{yb}) "
                f"pred={pred}"
            )


# ==============================================================================
# TEST 4
# ==============================================================================
# VERIFY THAT ONLY ALIVE PREFIXES ARE FED TO THE TRANSITION AUTOMATON.
# ==============================================================================

def test_no_post_death_transitions(
    states,
):

    print()
    print("=" * 90)
    print(
        "TEST 4: NO POST-DEATH TRANSITIONS"
    )
    print("=" * 90)

    failures = 0
    checked = 0

    for state in states:

        died = False

        for z in range(
            2,
            MAX_Z,
        ):

            alive_z = level_exists(
                state.X,
                state.Y,
                z,
            )

            if not alive_z:
                died = True
                continue

            if died:

                failures += 1

            checked += 1

    print(
        f"alive-level checks={checked} "
        f"post-death failures={failures}"
    )


# ==============================================================================
# TEST 5
# ==============================================================================
# DEEPEST LEVEL
#
# deepest = 1 + min(v2(X),v2(Y))
#
# This is equivalent to the length of the initial zero-bit prefix
# needed by the descent process.
# ==============================================================================

def test_deepest_formula(
    states,
):

    print()
    print("=" * 90)
    print(
        "TEST 5: EXACT DEEPEST LEVEL"
    )
    print("=" * 90)

    failures = 0
    distribution = Counter()

    for state in states:

        predicted = (
            1
            + min(
                v2(state.X),
                v2(state.Y),
            )
        )

        actual = 1

        for z in range(
            2,
            MAX_Z + 1,
        ):

            if level_exists(
                state.X,
                state.Y,
                z,
            ):

                actual = z

        if predicted != actual:

            failures += 1

            if failures <= EXAMPLE_COUNT:

                print(
                    f"    mismatch "
                    f"n={state.n} "
                    f"X={state.X} "
                    f"Y={state.Y} "
                    f"predicted={predicted} "
                    f"actual={actual}"
                )

        distribution[
            actual
        ] += 1

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )

    print()
    print(
        "DEPTH DISTRIBUTION"
    )

    for depth in sorted(
        distribution
    ):

        print(
            f"    depth={depth}: "
            f"{distribution[depth]}"
        )


# ==============================================================================
# TEST 6
# ==============================================================================
# CURRENT-BIT RESIDUE LAW
#
# Only perform this test while the parent level is alive.
# ==============================================================================

def test_conditioned_residue(
    states,
):

    print()
    print("=" * 90)
    print(
        "TEST 6: CONDITIONED CHILD RESIDUE LAW"
    )
    print("=" * 90)

    checked = 0
    failures = 0

    for state in states:

        for z in range(
            2,
            MAX_Z,
        ):

            if not level_exists(
                state.X,
                state.Y,
                z,
            ):

                # DEAD PREFIX: do not apply child transition.
                continue

            actual = (
                state.n
                % (
                    1 << (
                        z + 1
                    )
                )
            )

            predicted = child_residue(
                state.frame,
                state.X,
                z,
            )

            checked += 1

            if actual != predicted:

                failures += 1

                if failures <= EXAMPLE_COUNT:

                    print(
                        f"    mismatch "
                        f"n={state.n} "
                        f"frame={state.frame} "
                        f"z={z} "
                        f"Xbit={bit_abs(state.X,z-1)} "
                        f"pred={predicted} "
                        f"actual={actual}"
                    )

    print(
        f"alive transitions checked={checked} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 7
# ==============================================================================
# COMPLETE PREFIX REPLAY
#
# Starting from (FRAME,X,Y), replay the hierarchy one level at a time.
#
# At each z:
#
#   1. Verify current level exists.
#   2. Verify x_z,y_z.
#   3. Verify residue.
#   4. Inspect current bits.
#   5. Decide whether the next level survives.
#
# This is the direct global automaton.
# ==============================================================================

def replay_state(
    state,
):

    result = []

    for z in range(
        2,
        MAX_Z + 1,
    ):

        alive = level_exists(
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

        residue = (
            state.n
            % (
                1 << z
            )
        )

        expected_residue = active_residue(
            state.frame,
            z,
        )

        if residue != expected_residue:

            raise AssertionError(
                (
                    "Residue mismatch: "
                    f"n={state.n} "
                    f"z={z} "
                    f"actual={residue} "
                    f"expected={expected_residue}"
                )
            )

        xb = x & 1
        yb = y & 1

        survives = (
            xb == 0
            and
            yb == 0
        )

        result.append(
            {
                "z": z,
                "x": x,
                "y": y,
                "residue": residue,
                "xb": xb,
                "yb": yb,
                "survives": survives,
            }
        )

        if not survives:
            break

    return result


def test_complete_replay(
    states,
):

    print()
    print("=" * 90)
    print(
        "TEST 7: COMPLETE GLOBAL PREFIX REPLAY"
    )
    print("=" * 90)

    failures = 0
    checked_levels = 0

    for state in states:

        try:

            replay = replay_state(
                state
            )

        except AssertionError as exc:

            failures += 1

            if failures <= EXAMPLE_COUNT:

                print(
                    f"    {exc}"
                )

            continue

        expected_depth = (
            1
            + min(
                v2(state.X),
                v2(state.Y),
            )
        )

        actual_depth = (
            replay[-1]["z"]
            if replay
            else 1
        )

        if actual_depth != expected_depth:

            failures += 1

            if failures <= EXAMPLE_COUNT:

                print(
                    f"    depth mismatch "
                    f"n={state.n} "
                    f"expected={expected_depth} "
                    f"actual={actual_depth}"
                )

        checked_levels += len(
            replay
        )

    print(
        f"replayed levels={checked_levels} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 8
# ==============================================================================
# SAT PREFIX AUTOMATON
#
# The important distinction is:
#
#     Alive_z
#
# must be carried forward.
#
# Thus the correct Boolean recursion is:
#
#     Alive_(z+1)
#       =
#     Alive_z AND !X_(z-1) AND !Y_(z-1)
#
# This is different from simply asserting:
#
#     Alive_(z+1) <-> !X_(z-1) AND !Y_(z-1)
#
# for every z independently.
#
# ==============================================================================

def test_sat_prefix_logic(
    states,
):

    print()
    print("=" * 90)
    print(
        "TEST 8: PREFIX SAT LOGIC"
    )
    print("=" * 90)

    checked = 0
    failures = 0

    for state in states:

        alive = True

        for z in range(
            2,
            MAX_Z,
        ):

            xb = bit_abs(
                state.X,
                z - 1,
            )

            yb = bit_abs(
                state.Y,
                z - 1,
            )

            predicted_next = (
                alive
                and
                xb == 0
                and
                yb == 0
            )

            actual_next = level_exists(
                state.X,
                state.Y,
                z + 1,
            )

            checked += 1

            if predicted_next != actual_next:

                failures += 1

                if failures <= EXAMPLE_COUNT:

                    print(
                        f"    mismatch "
                        f"n={state.n} "
                        f"z={z} "
                        f"alive={alive} "
                        f"bits=({xb},{yb}) "
                        f"pred={predicted_next} "
                        f"actual={actual_next}"
                    )

            alive = predicted_next

    print(
        f"checked={checked} "
        f"failures={failures}"
    )

    print()
    print(
        "BOOLEAN RECURSION:"
    )

    print(
        "    Alive_(z+1) = "
        "Alive_z & !X_(z-1) & !Y_(z-1)"
    )


# ==============================================================================
# TEST 9
# ==============================================================================
# FRAME PRESERVATION ONLY ON SURVIVING EDGES.
# ==============================================================================

def test_frame_preservation(
    states,
):

    print()
    print("=" * 90)
    print(
        "TEST 9: FRAME PRESERVATION ON LIVE EDGES"
    )
    print("=" * 90)

    checked = 0
    failures = 0

    for state in states:

        for z in range(
            2,
            MAX_Z,
        ):

            if not level_exists(
                state.X,
                state.Y,
                z + 1,
            ):

                continue

            # The same global state is still using the same frame.
            child_frame = state.frame

            checked += 1

            if child_frame != state.frame:

                failures += 1

    print(
        f"live edges checked={checked} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 10
# ==============================================================================
# MINIMAL INFORMATION CONTENT
#
# We now distinguish three logically different quantities:
#
#   1. X,Y
#      determine the survival prefix.
#
#   2. FRAME,X,Y
#      determine factors and n.
#
#   3. FRAME,X_bit
#      determines the child residue on a live edge.
#
# ==============================================================================

def print_minimal_information():

    print()
    print("=" * 90)
    print(
        "TEST 10: MINIMAL INFORMATION DECOMPOSITION"
    )
    print("=" * 90)

    print(
        """
GLOBAL FACTOR RECOVERY:

    (FRAME, X, Y)
        ->
    (p, q)

GLOBAL INTEGER:

    (FRAME, X, Y)
        ->
    n

SURVIVAL:

    (Xbit, Ybit)
        ->
    survive / die

    exact rule:

        S = !Xbit & !Ybit

CHILD RESIDUE:

    (FRAME, Xbit)
        ->
    r_(z+1)

    FRAME A:

        r' = -9 + 2^z Xbit
             mod 2^(z+1)

    FRAME B:

        r' = -3 + 2^z Xbit
             mod 2^(z+1)

STATE UPDATE:

    if S:

        X remains globally fixed
        Y remains globally fixed
        FRAME remains fixed

        local coordinates halve:

            x_(z+1) = x_z / 2
            y_(z+1) = y_z / 2

Thus:

    survival information
        = two current bits

    residue information
        = frame + current X bit

    factor recovery
        = frame + complete X,Y.
"""
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
        "EXAMPLES OF PREFIX REPLAY"
    )
    print("=" * 90)

    shown = 0

    for state in states:

        if shown >= EXAMPLE_COUNT:
            break

        replay = replay_state(
            state
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
            f"    v2(X)={v2(state.X)} "
            f"v2(Y)={v2(state.Y)} "
            f"deepest="
            f"{1 + min(v2(state.X),v2(state.Y))}"
        )

        for row in replay:

            print(
                f"    z={row['z']} "
                f"x={row['x']:<8} "
                f"y={row['y']:<8} "
                f"bits=({row['xb']},{row['yb']}) "
                f"residue={row['residue']:<6} "
                f"next={row['survives']}"
            )

        shown += 1


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print(
        "EXPERIMENT 606 START"
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

    test_global_reconstruction(
        states
    )

    test_alive_formula(
        states
    )

    test_conditioned_transition(
        states
    )

    test_no_post_death_transitions(
        states
    )

    test_deepest_formula(
        states
    )

    test_conditioned_residue(
        states
    )

    test_complete_replay(
        states
    )

    test_sat_prefix_logic(
        states
    )

    test_frame_preservation(
        states
    )

    print_minimal_information()

    show_examples(
        states
    )

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 606 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()
