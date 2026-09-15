#!/usr/bin/env python3

# ==============================================================================
# EXPERIMENT 601
# ==============================================================================
# EXACT SURVIVAL AUTOMATON + MINIMAL STATE TEST
#
# Goal
# ----
# Experiment 600 established experimentally that:
#
#   1. x,y even -> same-frame child exists
#   2. odd x or odd y -> no child representation
#   3. the child residue is determined by frame + x parity
#
# This experiment asks whether the entire hierarchy can therefore be
# compressed into a tiny deterministic automaton.
#
# Candidate state:
#
#     (frame, x mod 2, y mod 2)
#
# Candidate transition:
#
#     if x even and y even:
#         child = same frame
#         x' = x/2
#         y' = y/2
#
#     otherwise:
#         DEAD
#
# We test:
#
#   A. exact transition law
#   B. minimality of the state variables
#   C. whether y parity is redundant once frame and x parity are known
#   D. whether frame is redundant once the parent residue is known
#   E. whether survival can be expressed directly as n modulo a larger power
#   F. whether repeated survival is equivalent to divisibility of the
#      normalized x,y coordinates
#   G. exact SAT-style clause representation
#
# No external files.
# No external sources.
# No output files.
# ==============================================================================

from __future__ import annotations

from collections import defaultdict, Counter
from dataclasses import dataclass
from math import isqrt


# ==============================================================================
# CONFIGURATION
# ==============================================================================

MAX_N = 2_000_000

MIN_Z = 2
MAX_Z = 12

SHOW_EXAMPLES = 10


# ==============================================================================
# STATE
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
# FRAME A
# ==============================================================================

def frame_A(
    p,
    q,
    z,
):

    k = 1 << (z - 1)

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


# ==============================================================================
# FRAME B
# ==============================================================================

def frame_B(
    p,
    q,
    z,
):

    k = 1 << (z - 1)

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

            a = frame_A(
                p,
                q,
                z,
            )

            if a is not None:

                levels[z][n] = State(
                    n=n,
                    p=p,
                    q=q,
                    z=z,
                    residue=n % (1 << z),
                    frame="A",
                    x=a[0],
                    y=a[1],
                )

                continue

            b = frame_B(
                p,
                q,
                z,
            )

            if b is not None:

                levels[z][n] = State(
                    n=n,
                    p=p,
                    q=q,
                    z=z,
                    residue=n % (1 << z),
                    frame="B",
                    x=b[0],
                    y=b[1],
                )

    return levels


# ==============================================================================
# CHILD
# ==============================================================================

def child_state(
    levels,
    z,
    n,
):

    return levels[
        z + 1
    ].get(n)


# ==============================================================================
# PREDICTED CHILD RESIDUE
# ==============================================================================

def predicted_child_residue(
    state,
):

    z = state.z

    modulus = 1 << (
        z + 1
    )

    bit = state.x & 1

    if z >= 3:

        if state.frame == "A":

            return (
                -9
                + (1 << z) * bit
            ) % modulus

        if state.frame == "B":

            return (
                -3
                + (1 << z) * bit
            ) % modulus

    # Exact z=2 calculation.

    k = 1 << (
        z - 1
    )

    if state.frame == "A":

        value = (
            (k * (state.y - state.x) + 3)
            *
            (k * (state.y + state.x) - 3)
        )

        return value % modulus

    numerator = (
        (k * (state.x + state.y) - 3)
        *
        (k * (state.y - state.x) + 3)
    )

    if numerator % 3:
        return None

    return (
        numerator // 3
    ) % modulus


# ==============================================================================
# AUTOMATON PREDICTION
# ==============================================================================

def automaton_prediction(
    state,
):

    """
    Candidate minimal transition.

    EVEN / EVEN:
        child survives and keeps frame.

    Anything else:
        DEAD.
    """

    xb = state.x & 1
    yb = state.y & 1

    if xb != 0 or yb != 0:

        return None

    return state.frame


# ==============================================================================
# [1] EXACT AUTOMATON TEST
# ==============================================================================

def test_automaton(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 1: EXACT SURVIVAL AUTOMATON"
    )
    print("=" * 90)

    global_failures = 0

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        tested = 0
        failures = 0

        for n, state in levels[z].items():

            child = child_state(
                levels,
                z,
                n,
            )

            predicted = (
                automaton_prediction(
                    state
                )
            )

            actual = (
                None
                if child is None
                else child.frame
            )

            tested += 1

            if predicted != actual:

                failures += 1

                if failures <= SHOW_EXAMPLES:

                    print(
                        f"    "
                        f"z={z} "
                        f"n={n} "
                        f"frame={state.frame} "
                        f"x={state.x} "
                        f"y={state.y} "
                        f"pred={predicted} "
                        f"actual={actual}"
                    )

        global_failures += failures

        print(
            f"z={z} -> z={z+1} "
            f"tested={tested} "
            f"failures={failures}"
        )

    print()
    print(
        f"GLOBAL failures={global_failures}"
    )


# ==============================================================================
# [2] STATE TABLE
# ==============================================================================

def state_table(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 2: COMPLETE PARITY STATE TABLE"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        table = defaultdict(
            Counter
        )

        for state in levels[z].values():

            child = child_state(
                levels,
                z,
                state.n,
            )

            key = (
                state.frame,
                state.x & 1,
                state.y & 1,
            )

            outcome = (
                state.frame
                if child is not None
                else "DEAD"
            )

            table[key][
                outcome
            ] += 1

        print()
        print(
            f"z={z}"
        )

        for key in sorted(table):

            print(
                f"    "
                f"state={key} "
                f"-> "
                f"{dict(table[key])}"
            )


# ==============================================================================
# [3] IS y PARITY REDUNDANT?
# ==============================================================================

def test_y_redundancy(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 3: IS y-PARITY REDUNDANT?"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        groups = defaultdict(
            set
        )

        for state in levels[z].values():

            child = child_state(
                levels,
                z,
                state.n,
            )

            outcome = (
                child is not None
            )

            key = (
                state.frame,
                state.x & 1,
            )

            groups[key].add(
                outcome
            )

        deterministic = all(
            len(values) == 1
            for values in groups.values()
        )

        print(
            f"z={z} "
            f"deterministic="
            f"{deterministic}"
        )

        if not deterministic:

            for key, values in sorted(
                groups.items()
            ):

                if len(values) > 1:

                    print(
                        f"    counterexample state={key}"
                        f" outcomes={values}"
                    )


# ==============================================================================
# [4] IS FRAME REDUNDANT GIVEN RESIDUE?
# ==============================================================================

def test_frame_from_residue(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 4: RESIDUE -> FRAME DETERMINISM"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        mapping = defaultdict(
            set
        )

        for state in levels[z].values():

            mapping[
                state.residue
            ].add(
                state.frame
            )

        deterministic = all(
            len(values) == 1
            for values in mapping.values()
        )

        print(
            f"z={z} "
            f"residue->frame "
            f"deterministic="
            f"{deterministic}"
        )

        collisions = [
            (
                residue,
                values,
            )
            for residue, values
            in mapping.items()
            if len(values) > 1
        ]

        if collisions:

            print(
                "    collisions:"
            )

            for residue, values in (
                collisions[:10]
            ):

                print(
                    f"        "
                    f"{residue}: "
                    f"{values}"
                )


# ==============================================================================
# [5] RESIDUE AUTOMATON
# ==============================================================================

def residue_automaton(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 5: RESIDUE AUTOMATON"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        mapping = defaultdict(
            set
        )

        for state in levels[z].values():

            child = child_state(
                levels,
                z,
                state.n,
            )

            if child is None:

                outcome = (
                    "DEAD"
                )

            else:

                outcome = (
                    child.frame,
                    child.residue,
                )

            mapping[
                (
                    state.frame,
                    state.residue,
                    state.x & 1,
                    state.y & 1,
                )
            ].add(
                outcome
            )

        deterministic = all(
            len(values) == 1
            for values in mapping.values()
        )

        print(
            f"z={z} "
            f"deterministic="
            f"{deterministic} "
            f"states={len(mapping)}"
        )


# ==============================================================================
# [6] DIRECT SAT CLAUSES
# ==============================================================================

def sat_clauses(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 6: SAT-STYLE STATE CLAUSES"
    )
    print("=" * 90)

    print(
        """
Boolean variables:

    F
        0 = B
        1 = A

    X
        x mod 2

    Y
        y mod 2

    S
        child survives

    C
        child frame
        0 = B
        1 = A

Candidate transition:

    S <-> (!X & !Y)

and when S is true:

    C <-> F

Equivalent CNF for survival:

    ( X \/ S )
    ( Y \/ S )
    ( !S \/ !X )
    ( !S \/ !Y )

Frame preservation:

    ( !S \/ !F \/ C )
    ( !S \/ F \/ !C )

The experiment verifies these clauses against every observed state.
"""
    )

    failures = 0
    total = 0

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        for state in levels[z].values():

            child = child_state(
                levels,
                z,
                state.n,
            )

            X = state.x & 1
            Y = state.y & 1

            F = 1 if (
                state.frame == "A"
            ) else 0

            S = 1 if (
                child is not None
            ) else 0

            C = (
                1
                if child is not None
                and child.frame == "A"
                else 0
            )

            # --------------------------------------------------------------
            # Survival clauses
            # --------------------------------------------------------------

            clauses = [

                (
                    X == 1
                    or S == 1
                ),

                (
                    Y == 1
                    or S == 1
                ),

                (
                    S == 0
                    or X == 0
                ),

                (
                    S == 0
                    or Y == 0
                ),

                (
                    S == 0
                    or F == C
                ),

            ]

            total += 1

            if not all(clauses):

                failures += 1

    print(
        f"checked={total} "
        f"failures={failures}"
    )


# ==============================================================================
# [7] SURVIVAL AS n CONGRUENCE
# ==============================================================================

def direct_n_congruence(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 7: SURVIVAL AS DIRECT n-CONGRUENCE"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        modulus = 1 << (
            z + 1
        )

        groups = defaultdict(
            Counter
        )

        for state in levels[z].values():

            child = child_state(
                levels,
                z,
                state.n,
            )

            groups[
                (
                    state.frame,
                    state.n % modulus,
                )
            ][
                child is not None
            ] += 1

        deterministic = True

        ambiguous = []

        for key, counts in groups.items():

            if (
                len(counts) > 1
            ):

                deterministic = False
                ambiguous.append(
                    (
                        key,
                        dict(counts),
                    )
                )

        print(
            f"z={z} "
            f"n mod {modulus} "
            f"deterministic="
            f"{deterministic}"
        )

        if ambiguous:

            for key, counts in (
                ambiguous[:5]
            ):

                print(
                    f"    "
                    f"{key} -> {counts}"
                )


# ==============================================================================
# [8] SURVIVAL AS FACTOR PARITY
# ==============================================================================

def factor_parity_test(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 8: SURVIVAL AS p,q MOD 2^m"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        for bits in range(
            1,
            min(
                z + 1,
                7,
            ),
        ):

            modulus = 1 << bits

            mapping = defaultdict(
                set
            )

            for state in levels[z].values():

                child = child_state(
                    levels,
                    z,
                    state.n,
                )

                key = (
                    state.frame,
                    state.p % modulus,
                    state.q % modulus,
                )

                mapping[key].add(
                    child is not None
                )

            deterministic = all(
                len(values) == 1
                for values in mapping.values()
            )

            if deterministic:

                print(
                    f"z={z} "
                    f"bits={bits} "
                    f"mod={modulus} "
                    f"DETERMINISTIC"
                )


# ==============================================================================
# [9] DEEP SURVIVAL CHAIN
# ==============================================================================

def chain_test(
    levels,
):

    print()
    print("=" * 90)
    print(
        "TEST 9: COMPLETE SURVIVAL CHAIN"
    )
    print("=" * 90)

    maximum_depth = Counter()

    for n in levels[MIN_Z]:

        deepest = MIN_Z - 1

        for z in range(
            MIN_Z,
            MAX_Z + 1,
        ):

            if n in levels[z]:

                deepest = z

        maximum_depth[
            deepest
        ] += 1

    print(
        "deepest level distribution:"
    )

    for depth, count in sorted(
        maximum_depth.items()
    ):

        print(
            f"    "
            f"depth={depth}: "
            f"{count}"
        )

    print()

    print(
        "divisibility prediction:"
    )

    failures = 0

    for n in levels[MIN_Z]:

        state = levels[
            MIN_Z
        ][n]

        deepest = (
            min(
                (
                    (
                        state.x.bit_length()
                        - (
                            abs(state.x)
                            & -abs(state.x)
                        ).bit_length()
                        + 1
                        if state.x != 0
                        else 10**9
                    ),
                    (
                        state.y.bit_length()
                        - (
                            abs(state.y)
                            & -abs(state.y)
                        ).bit_length()
                        + 1
                        if state.y != 0
                        else 10**9
                    ),
                )
            )
        )

        # The above expression is intentionally replaced below with
        # a direct valuation for clarity.

        def valuation2(v):

            if v == 0:
                return 10**9

            v = abs(v)

            result = 0

            while (
                v & 1
            ) == 0:

                result += 1
                v >>= 1

            return result

        predicted = (
            min(
                valuation2(state.x),
                valuation2(state.y),
            )
            + MIN_Z
        )

        actual = (
            max(
                z
                for z in range(
                    MIN_Z,
                    MAX_Z + 1,
                )
                if n in levels[z]
            )
        )

        # Only compare when predicted lies in our measured range.

        if predicted <= MAX_Z:

            if predicted != actual:

                failures += 1

    print(
        f"bounded-chain mismatches="
        f"{failures}"
    )


# ==============================================================================
# [10] EXAMPLES
# ==============================================================================

def examples(
    levels,
):

    print()
    print("=" * 90)
    print(
        "EXAMPLES OF AUTOMATON STATES"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        shown = 0

        print()
        print(
            f"z={z} -> z={z+1}"
        )

        for n in sorted(
            levels[z]
        ):

            if shown >= SHOW_EXAMPLES:
                break

            state = levels[z][n]

            child = child_state(
                levels,
                z,
                n,
            )

            print(
                f"    "
                f"n={n} "
                f"frame={state.frame} "
                f"residue={state.residue} "
                f"x={state.x} "
                f"y={state.y} "
                f"bits=("
                f"{state.x & 1},"
                f"{state.y & 1}) "
                f"-> "
                f"{child.frame + ':' + str(child.residue) if child else 'DEAD'}"
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
The candidate minimal hierarchy is:

    STATE_z
      =
    (frame_z, x_z mod 2, y_z mod 2)

with transition:

    x_z even AND y_z even
        ->
    SAME FRAME
        ->
    x_(z+1) = x_z/2
    y_(z+1) = y_z/2

and:

    otherwise
        ->
    DEAD.

For z >= 3 the child residue is independently determined by:

    FRAME A:
        r_(z+1)
            =
        -9 + 2^z (x_z mod 2)
        mod 2^(z+1)

    FRAME B:
        r_(z+1)
            =
        -3 + 2^z (x_z mod 2)
        mod 2^(z+1).

Thus there are potentially two distinct finite-state layers:

    1. residue generation
    2. representation survival.

The strongest possible outcome is that no information beyond

    frame
    x mod 2
    y mod 2

is required to determine descent.

That would turn the observed hierarchy into a very small deterministic
transition system suitable for a compact SAT encoding.

This experiment also checks whether some of those variables are
redundant, which is important because the smallest exact state is
more valuable than a merely correct state.
"""
    )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print(
        "EXPERIMENT 601 START"
    )
    print("=" * 90)

    print()
    print(
        "EXACT SURVIVAL AUTOMATON + MINIMAL STATE TEST"
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
            f"    z={z}"
            f" states={len(levels[z])}"
        )

    test_automaton(
        levels
    )

    state_table(
        levels
    )

    test_y_redundancy(
        levels
    )

    test_frame_from_residue(
        levels
    )

    residue_automaton(
        levels
    )

    sat_clauses(
        levels
    )

    direct_n_congruence(
        levels
    )

    factor_parity_test(
        levels
    )

    chain_test(
        levels
    )

    examples(
        levels
    )

    final_summary()

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 601 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()
