#!/usr/bin/env python3

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from math import isqrt


# ==============================================================================
# EXPERIMENT 627
# ==============================================================================
#
# CONDITIONAL SECOND-VALUATION SEARCH
#
# Experiment 626 established:
#
#   FRAME A:
#       primary = v2(n+9)
#
#   FRAME B:
#       primary = v2(n+3)
#
# but the primary valuation is not sufficient.
#
# We now condition on that primary valuation and exhaustively search for ONE
# additional valuation:
#
#       (primary_v2, v2(n+c))
#
# for c in a small offset range.
#
# Unlike the previous broad pair/triple searches, this experiment:
#
#   1. uses the theoretically preferred primary valuation;
#   2. searches only one additional observable;
#   3. verifies every candidate on the COMPLETE dataset;
#   4. prints the exact remaining collisions;
#   5. diagnoses those collisions using factor-side information.
#
# This answers:
#
#   "Can the depth be recovered from FRAME + two 2-adic valuations?"
#
# ==============================================================================


# ==============================================================================
# CONFIGURATION
# ==============================================================================

PRIME_LIMIT = 6250

OFFSET_MIN = -96
OFFSET_MAX = 96

INF = 10**9

MAX_COLLISION_EXAMPLES = 20


# ==============================================================================
# DATA STRUCTURE
# ==============================================================================

@dataclass(frozen=True)
class State:
    n: int
    p: int
    q: int


# ==============================================================================
# v2
# ==============================================================================

def v2(x: int) -> int:
    if x == 0:
        return INF

    x = abs(x)

    return (x & -x).bit_length() - 1


# ==============================================================================
# PRIME SIEVE
# ==============================================================================

def prime_sieve(limit: int) -> list[int]:

    sieve = bytearray(b"\x01") * (limit + 1)

    sieve[0] = 0
    sieve[1] = 0

    for p in range(2, isqrt(limit) + 1):

        if not sieve[p]:
            continue

        start = p * p

        count = ((limit - start) // p) + 1

        sieve[start:limit + 1:p] = b"\x00" * count

    return [
        p
        for p in range(3, limit + 1, 2)
        if sieve[p]
    ]


# ==============================================================================
# SEMIPRIME GENERATION
# ==============================================================================

def generate_states(
    primes: list[int],
) -> list[State]:

    states = []

    for i, p in enumerate(primes):

        for q in primes[i:]:

            states.append(
                State(
                    n=p * q,
                    p=p,
                    q=q,
                )
            )

    return states


# ==============================================================================
# FRAME
# ==============================================================================

def frame_of(n: int) -> str:

    r = n & 3

    if r == 3:
        return "A"

    if r == 1:
        return "B"

    raise ValueError(
        f"Unexpected odd semiprime residue n mod 4 = {r}"
    )


# ==============================================================================
# RAW FACTOR RESIDUALS
# ==============================================================================

def raw_residuals(
    frame: str,
    p: int,
    q: int,
) -> tuple[int, int]:

    if frame == "A":

        return (
            p - 3,
            q + 3,
        )

    return (
        p + 1,
        q - 3,
    )


# ==============================================================================
# TRUE DEPTH
# ==============================================================================

def depth_of(
    state: State,
) -> int:

    frame = frame_of(state.n)

    A, B = raw_residuals(
        frame,
        state.p,
        state.q,
    )

    alpha = v2(A)
    beta = v2(B)

    return (
        min(alpha, beta)
        + int(alpha == beta)
    )


# ==============================================================================
# BUILD GLOBAL ARRAYS
# ==============================================================================

def build_arrays(
    states: list[State],
):

    frames = []
    primary = []
    depths = []

    for state in states:

        frame = frame_of(state.n)

        c = (
            9
            if frame == "A"
            else 3
        )

        frames.append(frame)
        primary.append(v2(state.n + c))
        depths.append(depth_of(state))

    return (
        frames,
        primary,
        depths,
    )


# ==============================================================================
# AMBIGUITY COUNT
# ==============================================================================

def ambiguity_for_indices(
    indices: list[int],
    signatures: dict[int, tuple[int, ...]],
    depths: list[int],
):
    """
    Returns:

        ambiguous_group_count
        collision_examples

    A signature is ambiguous when the same signature maps to multiple depths.
    """

    first_depth = {}
    collisions = defaultdict(set)

    for i in indices:

        sig = signatures[i]
        d = depths[i]

        if sig not in first_depth:

            first_depth[sig] = d

        elif first_depth[sig] != d:

            collisions[sig].add(
                first_depth[sig]
            )
            collisions[sig].add(d)

    return (
        len(collisions),
        collisions,
    )


# ==============================================================================
# PRECOMPUTE VALUATIONS
# ==============================================================================

def precompute_valuations(
    states: list[State],
    offsets: list[int],
):

    columns = {}

    for c in offsets:

        columns[c] = [
            v2(state.n + c)
            for state in states
        ]

    return columns


# ==============================================================================
# TEST ONE PRIMARY / SECONDARY PAIR
# ==============================================================================

def test_pair(
    indices: list[int],
    primary: list[int],
    secondary: list[int],
    depths: list[int],
):

    signatures = {
        i: (
            primary[i],
            secondary[i],
        )
        for i in indices
    }

    ambiguous, collisions = (
        ambiguity_for_indices(
            indices,
            signatures,
            depths,
        )
    )

    return (
        ambiguous,
        collisions,
    )


# ==============================================================================
# SEARCH SECONDARY VALUATIONS
# ==============================================================================

def search_secondary_offsets(
    states: list[State],
    frames: list[str],
    primary: list[int],
    depths: list[int],
    columns: dict[int, list[int]],
):

    print()
    print("=" * 90)
    print("TEST 1: CONDITIONAL SECOND-VALUATION SEARCH")
    print("=" * 90)

    results = {}

    for frame in ("A", "B"):

        indices = [
            i
            for i, f in enumerate(frames)
            if f == frame
        ]

        primary_c = (
            9
            if frame == "A"
            else 3
        )

        rows = []

        for c in columns:

            if c == primary_c:
                continue

            ambiguous, collisions = test_pair(
                indices,
                primary,
                columns[c],
                depths,
            )

            rows.append(
                (
                    ambiguous,
                    c,
                    collisions,
                )
            )

        rows.sort(
            key=lambda x: (
                x[0],
                abs(x[1]),
            )
        )

        results[frame] = rows

        print()
        print(
            f"FRAME {frame}"
        )

        print(
            f"    primary=v2(n+{primary_c})"
        )

        print(
            "    BEST SECONDARY OFFSETS:"
        )

        for ambiguous, c, _ in rows[:20]:

            print(
                f"        c={c:4d} "
                f"ambiguous={ambiguous}"
            )

        exact = [
            row
            for row in rows
            if row[0] == 0
        ]

        print()
        print(
            f"    EXACT SECONDARY OFFSETS="
            f"{[row[1] for row in exact]}"
        )

    return results


# ==============================================================================
# GLOBAL FRAME-AWARE SEARCH
# ==============================================================================

def search_global_secondary(
    states: list[State],
    frames: list[str],
    primary: list[int],
    depths: list[int],
    columns: dict[int, list[int]],
):

    print()
    print("=" * 90)
    print("TEST 2: GLOBAL FRAME-AWARE TWO-VALUATION SIGNATURE")
    print("=" * 90)

    results = []

    for c in columns:

        # Frame itself is part of the signature.
        signatures = {
            i: (
                frames[i],
                primary[i],
                columns[c][i],
            )
            for i in range(len(states))
        }

        ambiguous, collisions = (
            ambiguity_for_indices(
                list(range(len(states))),
                signatures,
                depths,
            )
        )

        results.append(
            (
                ambiguous,
                c,
                collisions,
            )
        )

    results.sort(
        key=lambda x: (
            x[0],
            abs(x[1]),
        )
    )

    print(
        "    BEST GLOBAL SECONDARY OFFSETS:"
    )

    for ambiguous, c, _ in results[:25]:

        print(
            f"        c={c:4d} "
            f"ambiguous={ambiguous}"
        )

    exact = [
        row
        for row in results
        if row[0] == 0
    ]

    print()
    print(
        f"    EXACT GLOBAL OFFSETS="
        f"{[row[1] for row in exact]}"
    )

    return results


# ==============================================================================
# COLLISION DIAGNOSTICS
# ==============================================================================

def diagnose_collisions(
    states: list[State],
    frames: list[str],
    primary: list[int],
    depths: list[int],
    columns: dict[int, list[int]],
    results,
):

    print()
    print("=" * 90)
    print("TEST 3: COLLISION DIAGNOSTICS")
    print("=" * 90)

    # Take the best candidate for each frame.
    for frame in ("A", "B"):

        indices = [
            i
            for i, f in enumerate(frames)
            if f == frame
        ]

        primary_c = (
            9
            if frame == "A"
            else 3
        )

        candidates = [
            row
            for row in results[frame]
            if row[1] != primary_c
        ]

        if not candidates:
            continue

        ambiguous, c, collisions = (
            candidates[0]
        )

        print()
        print(
            f"FRAME {frame}"
        )

        print(
            f"    primary offset={primary_c}"
        )

        print(
            f"    secondary offset={c}"
        )

        print(
            f"    ambiguous groups={ambiguous}"
        )

        if not collisions:
            print(
                "    EXACT"
            )
            continue

        shown = 0

        for signature in sorted(collisions):

            print()
            print(
                f"    collision signature={signature}"
            )

            matching = [
                i
                for i in indices
                if (
                    primary[i],
                    columns[c][i],
                ) == signature
            ]

            by_depth = defaultdict(list)

            for i in matching:

                by_depth[
                    depths[i]
                ].append(i)

            for d in sorted(by_depth):

                print(
                    f"        depth={d}"
                )

                for i in by_depth[d][
                    :MAX_COLLISION_EXAMPLES
                ]:

                    s = states[i]

                    A, B = raw_residuals(
                        frame,
                        s.p,
                        s.q,
                    )

                    print(
                        f"            "
                        f"n={s.n:<10} "
                        f"p={s.p:<6} "
                        f"q={s.q:<6} "
                        f"A={A:<8} "
                        f"B={B:<8}"
                    )

            shown += 1

            if shown >= 5:
                break


# ==============================================================================
# TEST FACTOR-SIDE DISCRIMINATORS
# ==============================================================================

def diagnose_hidden_factor_state(
    states: list[State],
    frames: list[str],
    primary: list[int],
    depths: list[int],
):

    print()
    print("=" * 90)
    print("TEST 4: WHAT HIDDEN FACTOR VARIABLE RESOLVES THE COLLISIONS?")
    print("=" * 90)

    for frame in ("A", "B"):

        indices = [
            i
            for i, f in enumerate(frames)
            if f == frame
        ]

        print()
        print(
            f"FRAME {frame}"
        )

        # ----------------------------------------------------------------------
        # Candidate factor-side observables.
        # ----------------------------------------------------------------------

        observables = {}

        for k in range(2, 9):

            mod = 1 << k

            observables[
                f"p mod {mod}"
            ] = lambda s, mod=mod: (
                s.p % mod
            )

            observables[
                f"q mod {mod}"
            ] = lambda s, mod=mod: (
                s.q % mod
            )

            observables[
                f"(p+q) mod {mod}"
            ] = lambda s, mod=mod: (
                (s.p + s.q) % mod
            )

            observables[
                f"(p-q) mod {mod}"
            ] = lambda s, mod=mod: (
                (s.p - s.q) % mod
            )

        best = []

        for name, fn in observables.items():

            mapping = {}
            ambiguous = set()

            for i in indices:

                key = (
                    primary[i],
                    fn(states[i]),
                )

                d = depths[i]

                old = mapping.get(key)

                if old is None:

                    mapping[key] = d

                elif old != d:

                    ambiguous.add(key)

            best.append(
                (
                    len(ambiguous),
                    name,
                )
            )

        best.sort()

        print(
            "    best factor-side "
            "secondary observables:"
        )

        for ambiguous, name in best[:15]:

            print(
                f"        {name:<20} "
                f"ambiguous={ambiguous}"
            )


# ==============================================================================
# SEARCH SUM / DIFFERENCE VALUATIONS
# ==============================================================================

def search_factor_polynomial_valuations(
    states: list[State],
    frames: list[str],
    depths: list[int],
):

    print()
    print("=" * 90)
    print("TEST 5: FACTOR-SIDE VALUATION DIAGNOSTIC")
    print("=" * 90)

    for frame in ("A", "B"):

        indices = [
            i
            for i, f in enumerate(frames)
            if f == frame
        ]

        primary_c = (
            9
            if frame == "A"
            else 3
        )

        print()
        print(
            f"FRAME {frame}"
        )

        candidates = {
            "v2(A)": lambda A, B: v2(A),
            "v2(B)": lambda A, B: v2(B),
            "v2(A+B)": lambda A, B: v2(A + B),
            "v2(A-B)": lambda A, B: v2(A - B),
            "v2(gcd)": lambda A, B: v2(
                __import__("math").gcd(
                    abs(A),
                    abs(B),
                )
            ),
        }

        for name, fn in candidates.items():

            ambiguous = set()
            mapping = {}

            for i in indices:

                state = states[i]

                A, B = raw_residuals(
                    frame,
                    state.p,
                    state.q,
                )

                key = (
                    v2(state.n + primary_c),
                    fn(A, B),
                )

                d = depths[i]

                old = mapping.get(key)

                if old is None:

                    mapping[key] = d

                elif old != d:

                    ambiguous.add(key)

            print(
                f"    "
                f"{name:<12} "
                f"ambiguous={len(ambiguous)}"
            )


# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

def final_summary(
    frame_results,
    global_results,
):

    print()
    print("=" * 90)
    print("FINAL SUMMARY")
    print("=" * 90)

    exact_any = False

    for frame in ("A", "B"):

        exact = [
            row
            for row in frame_results[frame]
            if row[0] == 0
        ]

        if exact:
            exact_any = True

        print()
        print(
            f"FRAME {frame} exact secondary offsets:"
        )

        print(
            [
                row[1]
                for row in exact
            ]
        )

    global_exact = [
        row[1]
        for row in global_results
        if row[0] == 0
    ]

    if global_exact:
        exact_any = True

    print()
    print(
        "GLOBAL exact secondary offsets:"
    )

    print(
        global_exact
    )

    print()
    print(
        f"ANY EXACT TWO-VALUATION LAW = "
        f"{exact_any}"
    )

    print(
        """
Interpretation:

    If an exact offset exists for FRAME A:

        depth =
            f(
                v2(n+9),
                v2(n+c)
            )

    If an exact offset exists for FRAME B:

        depth =
            f(
                v2(n+3),
                v2(n+c)
            )

Then the factor-side information has collapsed to two
n-only 2-adic valuations.

If no exact offset exists, the collision diagnostics are
the important result: they identify which factor-side
variable survives after the primary n-only valuation.

The next search should then target that variable rather
than blindly increasing the number of offsets.
"""
    )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print("EXPERIMENT 627 START")
    print("=" * 90)
    print()
    print(
        "CONDITIONAL SECOND-VALUATION SEARCH"
    )

    # --------------------------------------------------------------------------
    # 1
    # --------------------------------------------------------------------------

    print()
    print(
        "[1] PRIME SIEVE"
    )

    primes = prime_sieve(
        PRIME_LIMIT
    )

    print(
        f"    odd primes={len(primes)}"
    )

    # --------------------------------------------------------------------------
    # 2
    # --------------------------------------------------------------------------

    print()
    print(
        "[2] SEMIPRIME GENERATION"
    )

    states = generate_states(
        primes
    )

    print(
        f"    semiprimes={len(states)}"
    )

    # --------------------------------------------------------------------------
    # 3
    # --------------------------------------------------------------------------

    print()
    print(
        "[3] BUILD GLOBAL TARGETS"
    )

    (
        frames,
        primary,
        depths,
    ) = build_arrays(
        states
    )

    # --------------------------------------------------------------------------
    # 4
    # --------------------------------------------------------------------------

    print()
    print(
        "[4] PRECOMPUTE VALUATION COLUMNS"
    )

    offsets = list(
        range(
            OFFSET_MIN,
            OFFSET_MAX + 1,
        )
    )

    print(
        f"    offsets={len(offsets)}"
    )

    columns = precompute_valuations(
        states,
        offsets,
    )

    # --------------------------------------------------------------------------
    # 5
    # --------------------------------------------------------------------------

    frame_results = search_secondary_offsets(
        states,
        frames,
        primary,
        depths,
        columns,
    )

    # --------------------------------------------------------------------------
    # 6
    # --------------------------------------------------------------------------

    global_results = search_global_secondary(
        states,
        frames,
        primary,
        depths,
        columns,
    )

    # --------------------------------------------------------------------------
    # 7
    # --------------------------------------------------------------------------

    diagnose_collisions(
        states,
        frames,
        primary,
        depths,
        columns,
        frame_results,
    )

    # --------------------------------------------------------------------------
    # 8
    # --------------------------------------------------------------------------

    diagnose_hidden_factor_state(
        states,
        frames,
        primary,
        depths,
    )

    # --------------------------------------------------------------------------
    # 9
    # --------------------------------------------------------------------------

    search_factor_polynomial_valuations(
        states,
        frames,
        depths,
    )

    # --------------------------------------------------------------------------
    # 10
    # --------------------------------------------------------------------------

    final_summary(
        frame_results,
        global_results,
    )

    print()
    print("=" * 90)
    print("EXPERIMENT 627 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
