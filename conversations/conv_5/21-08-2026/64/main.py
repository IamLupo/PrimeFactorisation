#!/usr/bin/env python3

from __future__ import annotations

from collections import defaultdict, Counter
from dataclasses import dataclass
from math import isqrt


# ==============================================================================
# EXPERIMENT 626
# ==============================================================================
#
# FRAME-CONDITIONED SINGLE-VALUATION COLLISION ANALYSIS
#
# Previous experiment:
#
#   FRAME A:
#       v2(n+9)
#
#   FRAME B:
#       v2(n+3)
#
# already gives almost deterministic depth.
#
# The previous search reported:
#
#   FRAME A:
#       depth_ambiguous = 1
#
#   FRAME B:
#       depth_ambiguous = 1
#
# Therefore we now:
#
#   1. identify the exact ambiguous valuation buckets;
#   2. print all target/depth values in those buckets;
#   3. test n mod 8,16,32,... as the smallest secondary observable;
#   4. test the odd part of (n+c);
#   5. search for the smallest exact signature;
#   6. search simple closed-form transformations of v2(n+c).
#
# This is much smaller than a triple search.
#
# ==============================================================================


# ==============================================================================
# CONFIG
# ==============================================================================

PRIME_LIMIT = 6250

INF = 10**9

MAX_RESIDUE_BITS = 16

EXAMPLE_LIMIT = 20


# ==============================================================================
# STATE
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
# ODD PART
# ==============================================================================

def odd_part(x: int) -> int:

    if x == 0:
        return 0

    x = abs(x)

    return x >> v2(x)


# ==============================================================================
# PRIME SIEVE
# ==============================================================================

def prime_sieve(
    limit: int,
) -> list[int]:

    sieve = bytearray(
        b"\x01"
    ) * (limit + 1)

    sieve[0] = 0
    sieve[1] = 0

    for p in range(
        2,
        isqrt(limit) + 1,
    ):

        if not sieve[p]:
            continue

        start = p * p

        count = (
            (limit - start)
            // p
        ) + 1

        sieve[
            start:
            limit + 1:
            p
        ] = b"\x00" * count

    return [
        p
        for p in range(
            3,
            limit + 1,
            2,
        )
        if sieve[p]
    ]


# ==============================================================================
# SEMIPRIMES
# ==============================================================================

def generate_states(
    primes: list[int],
) -> list[State]:

    states = []

    for i, p in enumerate(primes):

        for q in primes[i:]:

            states.append(
                State(
                    p * q,
                    p,
                    q,
                )
            )

    return states


# ==============================================================================
# FRAME
# ==============================================================================

def frame_from_n(
    n: int,
) -> str:

    r = n & 3

    if r == 3:
        return "A"

    if r == 1:
        return "B"

    raise ValueError(
        f"non-odd n={n}"
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

def true_depth(
    state: State,
) -> int:

    frame = frame_from_n(
        state.n
    )

    A, B = raw_residuals(
        frame,
        state.p,
        state.q,
    )

    va = v2(A)
    vb = v2(B)

    m = min(
        va,
        vb,
    )

    return (
        m
        + int(va == vb)
    )


# ==============================================================================
# BUILD TARGETS
# ==============================================================================

def build_targets(
    states: list[State],
):

    depths = []

    frames = []

    values = []

    for state in states:

        frame = frame_from_n(
            state.n
        )

        c = (
            9
            if frame == "A"
            else 3
        )

        value = v2(
            state.n + c
        )

        depth = true_depth(
            state
        )

        frames.append(frame)
        values.append(value)
        depths.append(depth)

    return (
        frames,
        values,
        depths,
    )


# ==============================================================================
# TEST ONE SIGNATURE
# ==============================================================================

def signature_ambiguity(
    indices: list[int],
    signature_values,
    depths,
):

    mapping = {}

    collisions = defaultdict(
        set
    )

    for i in indices:

        key = signature_values[i]
        d = depths[i]

        previous = mapping.get(
            key
        )

        if previous is None:

            mapping[key] = d

        elif previous != d:

            collisions[key].add(
                previous
            )
            collisions[key].add(
                d
            )

    return (
        len(collisions),
        mapping,
        collisions,
    )


# ==============================================================================
# TEST SINGLE VALUATION
# ==============================================================================

def test_single_valuation(
    states,
    frames,
    values,
    depths,
):

    print()
    print("=" * 90)
    print(
        "TEST 1: SINGLE FRAME-CONDITIONED VALUATION"
    )
    print("=" * 90)

    for frame in (
        "A",
        "B",
    ):

        indices = [
            i
            for i, f in enumerate(frames)
            if f == frame
        ]

        amb, mapping, collisions = (
            signature_ambiguity(
                indices,
                values,
                depths,
            )
        )

        print()
        print(
            f"FRAME {frame}"
        )

        print(
            f"    states={len(indices)}"
        )

        print(
            f"    ambiguous buckets={amb}"
        )

        if collisions:

            print(
                "    ambiguous valuation buckets:"
            )

            for key in sorted(
                collisions
            ):

                print(
                    f"        v={key} "
                    f"depths="
                    f"{sorted(collisions[key])}"
                )


# ==============================================================================
# SHOW AMBIGUOUS STATES
# ==============================================================================

def show_ambiguous_states(
    states,
    frames,
    values,
    depths,
):

    print()
    print("=" * 90)
    print(
        "TEST 2: EXACT AMBIGUOUS STATES"
    )
    print("=" * 90)

    for frame in (
        "A",
        "B",
    ):

        indices = [
            i
            for i, f in enumerate(frames)
            if f == frame
        ]

        buckets = defaultdict(
            list
        )

        for i in indices:

            buckets[
                values[i]
            ].append(i)

        ambiguous = []

        for value, bucket in (
            buckets.items()
        ):

            target_set = {
                depths[i]
                for i in bucket
            }

            if len(target_set) > 1:

                ambiguous.append(
                    (
                        value,
                        bucket,
                    )
                )

        print()
        print(
            f"FRAME {frame}"
        )

        if not ambiguous:

            print(
                "    no ambiguous buckets"
            )

            continue

        for value, bucket in ambiguous:

            print()
            print(
                f"    valuation={value}"
            )

            for i in bucket[
                :EXAMPLE_LIMIT
            ]:

                s = states[i]

                print(
                    f"        n={s.n:<8} "
                    f"p={s.p:<5} "
                    f"q={s.q:<5} "
                    f"depth={depths[i]}"
                )


# ==============================================================================
# TEST n MOD 2^k
# ==============================================================================

def test_n_residue_secondary(
    states,
    frames,
    values,
    depths,
):

    print()
    print("=" * 90)
    print(
        "TEST 3: v2(n+c) + n MOD 2^k"
    )
    print("=" * 90)

    global_exact = False

    for frame in (
        "A",
        "B",
    ):

        indices = [
            i
            for i, f in enumerate(frames)
            if f == frame
        ]

        c = (
            9
            if frame == "A"
            else 3
        )

        print()
        print(
            f"FRAME {frame}"
        )

        for bits in range(
            3,
            MAX_RESIDUE_BITS + 1,
        ):

            modulus = 1 << bits

            signatures = [
                (
                    values[i],
                    states[i].n
                    % modulus,
                )
                for i in indices
            ]

            amb, _, _ = (
                signature_ambiguity(
                    indices,
                    {
                        i: signatures[j]
                        for j, i
                        in enumerate(indices)
                    },
                    depths,
                )
            )

            print(
                f"    bits={bits:<2} "
                f"mod={modulus:<8} "
                f"ambiguous={amb}"
            )

            if amb == 0:

                print(
                    "        EXACT"
                )

                global_exact = True

                break

    return global_exact


# ==============================================================================
# TEST ODD PART
# ==============================================================================

def test_odd_part(
    states,
    frames,
    values,
    depths,
):

    print()
    print("=" * 90)
    print(
        "TEST 4: v2(n+c) + ODD PART"
    )
    print("=" * 90)

    exact = False

    for frame in (
        "A",
        "B",
    ):

        indices = [
            i
            for i, f in enumerate(frames)
            if f == frame
        ]

        c = (
            9
            if frame == "A"
            else 3
        )

        # signature:
        #
        #   (valuation, odd_part mod 2^k)
        #
        # The valuation itself determines the scale;
        # the odd part carries the first information beyond it.

        print()
        print(
            f"FRAME {frame}"
        )

        for bits in range(
            1,
            13,
        ):

            modulus = 1 << bits

            signatures = {
                i: (
                    values[i],
                    odd_part(
                        states[i].n + c
                    ) % modulus,
                )
                for i in indices
            }

            amb, _, _ = (
                signature_ambiguity(
                    indices,
                    signatures,
                    depths,
                )
            )

            print(
                f"    oddbits={bits:<2} "
                f"mod={modulus:<5} "
                f"ambiguous={amb}"
            )

            if amb == 0:

                print(
                    "        EXACT"
                )

                exact = True

                break

    return exact


# ==============================================================================
# SIMPLE FUNCTIONS OF v2(n+c)
# ==============================================================================

def test_simple_transformations(
    states,
    frames,
    values,
    depths,
):

    print()
    print("=" * 90)
    print(
        "TEST 5: SIMPLE FUNCTIONS OF v2(n+c)"
    )
    print("=" * 90)

    functions = {
        "v": lambda v: v,
        "v-1": lambda v: (
            v - 1
            if v < INF
            else INF
        ),
        "v+1": lambda v: (
            v + 1
            if v < INF
            else INF
        ),
        "min(v,2)": lambda v: min(v, 2),
        "min(v,3)": lambda v: min(v, 3),
        "min(v,4)": lambda v: min(v, 4),
        "min(v,5)": lambda v: min(v, 5),
        "parity(v)": lambda v: v & 1,
    }

    for frame in (
        "A",
        "B",
    ):

        indices = [
            i
            for i, f in enumerate(frames)
            if f == frame
        ]

        print()
        print(
            f"FRAME {frame}"
        )

        for name, fn in (
            functions.items()
        ):

            signature = {
                i: fn(values[i])
                for i in indices
            }

            amb, _, _ = (
                signature_ambiguity(
                    indices,
                    signature,
                    depths,
                )
            )

            print(
                f"    {name:<12} "
                f"ambiguous={amb}"
            )


# ==============================================================================
# SEARCH SECONDARY BITS
# ==============================================================================

def search_secondary_bit(
    states,
    frames,
    values,
    depths,
):

    print()
    print("=" * 90)
    print(
        "TEST 6: MINIMAL SECONDARY n-BIT"
    )
    print("=" * 90)

    for frame in (
        "A",
        "B",
    ):

        indices = [
            i
            for i, f in enumerate(frames)
            if f == frame
        ]

        c = (
            9
            if frame == "A"
            else 3
        )

        print()
        print(
            f"FRAME {frame}"
        )

        found = False

        # Test:
        #
        #   ( v2(n+c), bit_j(n) )
        #
        # and
        #
        #   ( v2(n+c), bit_j(n+c) )
        #
        for source_name, source_fn in (
            (
                "n",
                lambda n: n,
            ),
            (
                "n+c",
                lambda n: n + c,
            ),
        ):

            for bit in range(
                0,
                16,
            ):

                signature = {
                    i: (
                        values[i],
                        (
                            source_fn(
                                states[i].n
                            )
                            >> bit
                        ) & 1,
                    )
                    for i in indices
                }

                amb, _, _ = (
                    signature_ambiguity(
                        indices,
                        signature,
                        depths,
                    )
                )

                print(
                    f"    source={source_name:<3} "
                    f"bit={bit:<2} "
                    f"ambiguous={amb}"
                )

                if amb == 0:

                    print(
                        "        EXACT"
                    )

                    found = True

                    break

            if found:
                break


# ==============================================================================
# SEARCH LOCAL RESIDUE OF NORMALIZED NUMBER
# ==============================================================================

def test_normalized_residue(
    states,
    frames,
    values,
    depths,
):

    print()
    print("=" * 90)
    print(
        "TEST 7: NORMALIZED n+c ODD SIGNATURE"
    )
    print("=" * 90)

    for frame in (
        "A",
        "B",
    ):

        indices = [
            i
            for i, f in enumerate(frames)
            if f == frame
        ]

        c = (
            9
            if frame == "A"
            else 3
        )

        print()
        print(
            f"FRAME {frame}"
        )

        found = False

        for bits in range(
            1,
            13,
        ):

            modulus = 1 << bits

            signature = {}

            for i in indices:

                x = states[i].n + c
                v = values[i]

                if v >= INF:

                    normalized = 0

                else:

                    normalized = (
                        x
                        >> v
                    )

                signature[i] = (
                    v,
                    normalized % modulus,
                )

            amb, _, _ = (
                signature_ambiguity(
                    indices,
                    signature,
                    depths,
                )
            )

            print(
                f"    bits={bits:<2} "
                f"ambiguous={amb}"
            )

            if amb == 0:

                print(
                    "        EXACT"
                )

                found = True
                break


# ==============================================================================
# CROSS-FRAME UNIVERSAL SIGNATURE
# ==============================================================================

def test_universal_frame_signature(
    states,
    frames,
    values,
    depths,
):

    print()
    print("=" * 90)
    print(
        "TEST 8: UNIVERSAL FRAME-AWARE SIGNATURE"
    )
    print("=" * 90)

    # Use the frame to choose the correct shift:
    #
    #   A -> n+9
    #   B -> n+3
    #
    # Then ask whether a common formula exists.

    signature = {}

    for i, state in enumerate(
        states
    ):

        c = (
            9
            if frames[i] == "A"
            else 3
        )

        signature[i] = (
            values[i],
            state.n & 7,
        )

    amb, mapping, collisions = (
        signature_ambiguity(
            list(range(len(states))),
            signature,
            depths,
        )
    )

    print(
        f"    (v2(n+c), n mod 8)"
        f" ambiguous={amb}"
    )

    if amb:

        print(
            "    collision signatures:"
        )

        for key in sorted(
            collisions
        ):

            print(
                f"        {key} -> "
                f"{sorted(collisions[key])}"
            )


# ==============================================================================
# SUMMARY
# ==============================================================================

def final_summary(
    exact_residue,
    exact_odd,
):

    print()
    print("=" * 90)
    print(
        "FINAL SUMMARY"
    )
    print("=" * 90)

    print(
        f"""
FRAME-CONDITIONED DEPTH OBSERVATION

    FRAME A:
        primary valuation = v2(n+9)

    FRAME B:
        primary valuation = v2(n+3)

The previous experiment found only one ambiguous
valuation bucket per frame.

This experiment therefore tested:

    1. valuation alone
    2. valuation + n mod 2^k
    3. valuation + odd part
    4. valuation + a single extra bit
    5. normalized odd residue
    6. a universal frame-aware signature.

exact residue signature found = {exact_residue}
exact odd-part signature found = {exact_odd}

The important outcome is not merely finding another
formula, but determining the minimum additional
information required after:

    frame
    +
    v2(n+c).

If a one-bit refinement is exact, the depth process
has effectively collapsed to a tiny finite-state rule.
"""
    )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print(
        "EXPERIMENT 626 START"
    )
    print("=" * 90)
    print()
    print(
        "FRAME-CONDITIONED SINGLE-VALUATION COLLISION ANALYSIS"
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
        "[3] BUILD DEPTH TARGETS"
    )

    (
        frames,
        values,
        depths,
    ) = build_targets(
        states
    )

    # --------------------------------------------------------------------------
    # 4
    # --------------------------------------------------------------------------

    test_single_valuation(
        states,
        frames,
        values,
        depths,
    )

    # --------------------------------------------------------------------------
    # 5
    # --------------------------------------------------------------------------

    show_ambiguous_states(
        states,
        frames,
        values,
        depths,
    )

    # --------------------------------------------------------------------------
    # 6
    # --------------------------------------------------------------------------

    exact_residue = (
        test_n_residue_secondary(
            states,
            frames,
            values,
            depths,
        )
    )

    # --------------------------------------------------------------------------
    # 7
    # --------------------------------------------------------------------------

    exact_odd = (
        test_odd_part(
            states,
            frames,
            values,
            depths,
        )
    )

    # --------------------------------------------------------------------------
    # 8
    # --------------------------------------------------------------------------

    test_simple_transformations(
        states,
        frames,
        values,
        depths,
    )

    # --------------------------------------------------------------------------
    # 9
    # --------------------------------------------------------------------------

    search_secondary_bit(
        states,
        frames,
        values,
        depths,
    )

    # --------------------------------------------------------------------------
    # 10
    # --------------------------------------------------------------------------

    test_normalized_residue(
        states,
        frames,
        values,
        depths,
    )

    # --------------------------------------------------------------------------
    # 11
    # --------------------------------------------------------------------------

    test_universal_frame_signature(
        states,
        frames,
        values,
        depths,
    )

    # --------------------------------------------------------------------------
    # 12
    # --------------------------------------------------------------------------

    final_summary(
        exact_residue,
        exact_odd,
    )

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 626 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()
