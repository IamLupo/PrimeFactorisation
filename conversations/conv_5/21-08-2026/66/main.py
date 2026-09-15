#!/usr/bin/env python3

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from math import isqrt, gcd


# ==============================================================================
# EXPERIMENT 628
# ==============================================================================
#
# COLLISION-ONLY n-ADIC DISCRIMINATOR SEARCH
#
# EXPERIMENT 627 showed:
#
#   FRAME A:
#       primary = v2(n+9)
#
#   FRAME B:
#       primary = v2(n+3)
#
# but scanning every state for every offset is unnecessarily expensive.
#
# Here we ONLY examine primary-valuation buckets where multiple depths occur.
#
# The search is therefore conditioned on:
#
#       (frame, primary valuation)
#
# and asks:
#
#   Can a SECOND observable distinguish the colliding depths?
#
# Candidate observables:
#
#   1. v2(n+c)
#   2. (n+c) / 2^v  mod 2^k       [normalized odd part]
#   3. n mod 2^k
#   4. (n+c) mod 2^k
#
# The search terminates early when an exact discriminator is found.
#
# ==============================================================================


PRIME_LIMIT = 6250

OFFSET_MIN = -64
OFFSET_MAX = 64

ODD_BITS = 8
RESIDUE_BITS = 12

INF = 10**9


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

        sieve[start:limit + 1:p] = (
            b"\x00"
        ) * count

    return [
        p
        for p in range(3, limit + 1, 2)
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
        f"Unexpected odd n mod 4={r}"
    )


# ==============================================================================
# FACTOR RESIDUALS
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
# BUILD PRIMARY BUCKETS
# ==============================================================================

def build_primary_buckets(
    states: list[State],
):

    buckets = defaultdict(list)

    for i, state in enumerate(states):

        frame = frame_of(state.n)

        primary_c = (
            9
            if frame == "A"
            else 3
        )

        primary_v = v2(
            state.n + primary_c
        )

        depth = depth_of(state)

        buckets[
            (frame, primary_v)
        ].append(
            (
                i,
                depth,
            )
        )

    return buckets


# ==============================================================================
# COLLISION BUCKETS ONLY
# ==============================================================================

def extract_ambiguous_buckets(
    buckets,
):

    ambiguous = {}

    for key, rows in buckets.items():

        depths = {
            depth
            for _, depth in rows
        }

        if len(depths) > 1:

            ambiguous[key] = rows

    return ambiguous


# ==============================================================================
# COLLISION STATE SET
# ==============================================================================

def collect_collision_states(
    ambiguous,
):

    indices = set()

    for rows in ambiguous.values():

        for index, _ in rows:

            indices.add(index)

    return sorted(indices)


# ==============================================================================
# SECONDARY DISCRIMINATOR
# ==============================================================================

def exact_discriminator(
    indices,
    states,
    depths,
    signature_fn,
):

    mapping = {}

    for i in indices:

        sig = signature_fn(
            states[i]
        )

        d = depths[i]

        old = mapping.get(sig)

        if old is None:

            mapping[sig] = d

        elif old != d:

            return False

    return True


# ==============================================================================
# AMBIGUITY COUNT
# ==============================================================================

def ambiguity_count(
    indices,
    states,
    depths,
    signature_fn,
):

    mapping = {}
    ambiguous = set()

    for i in indices:

        sig = signature_fn(
            states[i]
        )

        d = depths[i]

        old = mapping.get(sig)

        if old is None:

            mapping[sig] = d

        elif old != d:

            ambiguous.add(sig)

    return len(ambiguous)


# ==============================================================================
# BUILD CANDIDATE INDEX GROUPS
# ==============================================================================

def candidate_collision_groups(
    rows,
    states,
    depths,
):

    groups = defaultdict(list)

    for i, _ in rows:

        state = states[i]

        groups[
            depth_of(state)
        ].append(i)

    return groups


# ==============================================================================
# PRINT PRIMARY STRUCTURE
# ==============================================================================

def print_primary_structure(
    ambiguous,
    states,
):

    print()
    print("=" * 90)
    print("TEST 1: PRIMARY COLLISION STRUCTURE")
    print("=" * 90)

    for frame in ("A", "B"):

        print()
        print(
            f"FRAME {frame}"
        )

        frame_rows = [
            (
                v,
                rows,
            )
            for (
                f,
                v
            ), rows in ambiguous.items()
            if f == frame
        ]

        frame_rows.sort()

        for v, rows in frame_rows:

            depths = sorted({
                depth
                for _, depth in rows
            })

            print(
                f"    v={v:<3} "
                f"states={len(rows):<6} "
                f"depths={depths}"
            )


# ==============================================================================
# SEARCH SECONDARY VALUATION
# ==============================================================================

def search_secondary_v2(
    ambiguous,
    states,
    depths,
):

    print()
    print("=" * 90)
    print("TEST 2: COLLISION-ONLY SECONDARY v2(n+c)")
    print("=" * 90)

    offsets = list(
        range(
            OFFSET_MIN,
            OFFSET_MAX + 1,
        )
    )

    # Only collision states.
    all_collision_indices = (
        collect_collision_states(
            ambiguous
        )
    )

    print(
        f"    collision states="
        f"{len(all_collision_indices)}"
    )

    exact = []

    best = []

    for c in offsets:

        # The primary shift is not useful as secondary.
        signature_fn = (
            lambda s, c=c: v2(
                s.n + c
            )
        )

        ambiguous_count = 0

        for key, rows in ambiguous.items():

            count = ambiguity_count(
                [i for i, _ in rows],
                states,
                depths,
                signature_fn,
            )

            ambiguous_count += count

            if count:

                # no early global calculation necessary
                pass

        best.append(
            (
                ambiguous_count,
                c,
            )
        )

        if ambiguous_count == 0:

            exact.append(c)

            print(
                f"    EXACT c={c}"
            )

            break

    best.sort()

    if not exact:

        print(
            "    no exact secondary valuation"
        )

        print(
            "    BEST CANDIDATES:"
        )

        for count, c in best[:20]:

            print(
                f"        c={c:4d} "
                f"ambiguous={count}"
            )

    return exact, best


# ==============================================================================
# SEARCH NORMALIZED ODD PART
# ==============================================================================

def normalized_odd(
    x: int,
) -> int:

    if x == 0:

        return 0

    ax = abs(x)

    v = v2(ax)

    return (
        ax
        >> v
    )


def search_normalized_odd(
    ambiguous,
    states,
    depths,
):

    print()
    print("=" * 90)
    print("TEST 3: COLLISION-ONLY NORMALIZED ODD PART")
    print("=" * 90)

    offsets = list(
        range(
            OFFSET_MIN,
            OFFSET_MAX + 1,
        )
    )

    best = []

    for c in offsets:

        signature_fn = (
            lambda s, c=c: (
                normalized_odd(
                    s.n + c
                )
                & (
                    (1 << ODD_BITS) - 1
                )
            )
        )

        total_ambiguous = 0

        for _, rows in ambiguous.items():

            total_ambiguous += (
                ambiguity_count(
                    [i for i, _ in rows],
                    states,
                    depths,
                    signature_fn,
                )
            )

        best.append(
            (
                total_ambiguous,
                c,
            )
        )

        if total_ambiguous == 0:

            print(
                f"    EXACT normalized odd "
                f"offset={c}"
            )

            return c

    best.sort()

    print(
        "    BEST NORMALIZED-ODD CANDIDATES:"
    )

    for count, c in best[:20]:

        print(
            f"        c={c:4d} "
            f"ambiguous={count}"
        )

    return None


# ==============================================================================
# SEARCH LOW RESIDUE
# ==============================================================================

def search_residue(
    ambiguous,
    states,
    depths,
):

    print()
    print("=" * 90)
    print("TEST 4: COLLISION-ONLY LOW n RESIDUE")
    print("=" * 90)

    for bits in (
        4,
        5,
        6,
        8,
        10,
        12,
        14,
    ):

        mask = (
            (1 << bits) - 1
        )

        signature_fn = (
            lambda s, mask=mask: (
                s.n & mask
            )
        )

        total = 0

        for _, rows in ambiguous.items():

            total += ambiguity_count(
                [i for i, _ in rows],
                states,
                depths,
                signature_fn,
            )

        print(
            f"    bits={bits:<2} "
            f"ambiguous={total}"
        )

        if total == 0:

            print(
                f"    EXACT n mod 2^{bits}"
            )

            return bits

    return None


# ==============================================================================
# SEARCH FACTOR-SIDE INFORMATION ON COLLISIONS
# ==============================================================================

def search_factor_side(
    ambiguous,
    states,
    depths,
):

    print()
    print("=" * 90)
    print("TEST 5: FACTOR-SIDE COLLISION SIGNATURE")
    print("=" * 90)

    candidates = {
        "p mod 8":
            lambda s: s.p & 7,

        "q mod 8":
            lambda s: s.q & 7,

        "p mod 16":
            lambda s: s.p & 15,

        "q mod 16":
            lambda s: s.q & 15,

        "(p+q) mod 16":
            lambda s: (s.p + s.q) & 15,

        "(p-q) mod 16":
            lambda s: (s.p - s.q) & 15,

        "(p+q) mod 32":
            lambda s: (s.p + s.q) & 31,

        "(p-q) mod 32":
            lambda s: (s.p - s.q) & 31,

        "p*q mod 32":
            lambda s: (s.p * s.q) & 31,
    }

    rows = []

    for name, fn in candidates.items():

        total = 0

        for _, bucket in ambiguous.items():

            total += ambiguity_count(
                [i for i, _ in bucket],
                states,
                depths,
                fn,
            )

        rows.append(
            (
                total,
                name,
            )
        )

    rows.sort()

    for count, name in rows:

        print(
            f"    {name:<20} "
            f"ambiguous={count}"
        )


# ==============================================================================
# EXPLICIT COLLISION EXAMPLES
# ==============================================================================

def print_collision_examples(
    ambiguous,
    states,
    depths,
):

    print()
    print("=" * 90)
    print("TEST 6: EXPLICIT PRIMARY COLLISIONS")
    print("=" * 90)

    shown = 0

    for (
        key,
        rows,
    ) in sorted(
        ambiguous.items()
    ):

        if shown >= 10:

            break

        frame, primary = key

        print()
        print(
            f"    FRAME={frame} "
            f"primary={primary}"
        )

        by_depth = defaultdict(list)

        for i, depth in rows:

            by_depth[
                depth
            ].append(i)

        for depth in sorted(by_depth):

            example_indices = (
                by_depth[depth][:5]
            )

            print(
                f"        depth={depth}"
            )

            for i in example_indices:

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


# ==============================================================================
# CROSS-CHECK WITH GLOBAL X,Y
# ==============================================================================

def global_xy(
    frame: str,
    p: int,
    q: int,
) -> tuple[int, int]:

    if frame == "A":

        X = (
            q - p + 6
        ) // 2

        Y = (
            p + q
        ) // 2

    else:

        X = (
            3 * p - q + 6
        ) // 2

        Y = (
            3 * p + q
        ) // 2

    return X, Y


def collision_xy_analysis(
    ambiguous,
    states,
    depths,
):

    print()
    print("=" * 90)
    print("TEST 7: PRIMARY COLLISION -> GLOBAL X,Y")
    print("=" * 90)

    mapping = {}

    ambiguous_count = 0

    for (
        key,
        rows,
    ) in ambiguous.items():

        for i, depth in rows:

            s = states[i]

            frame = key[0]

            X, Y = global_xy(
                frame,
                s.p,
                s.q,
            )

            sig = (
                v2(X),
                v2(Y),
            )

            old = mapping.get(
                (
                    key,
                    sig,
                )
            )

            if old is None:

                mapping[
                    (
                        key,
                        sig,
                    )
                ] = depth

            elif old != depth:

                ambiguous_count += 1

    print(
        f"    (primary, v2X, v2Y) "
        f"ambiguous={ambiguous_count}"
    )

    if ambiguous_count == 0:

        print(
            "    EXACT"
        )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print("EXPERIMENT 628 START")
    print("=" * 90)

    # --------------------------------------------------------------------------
    # Prime generation
    # --------------------------------------------------------------------------

    print()
    print("[1] PRIME SIEVE")

    primes = prime_sieve(
        PRIME_LIMIT
    )

    print(
        f"    odd primes={len(primes)}"
    )

    # --------------------------------------------------------------------------
    # Semiprimes
    # --------------------------------------------------------------------------

    print()
    print("[2] SEMIPRIME GENERATION")

    states = generate_states(
        primes
    )

    print(
        f"    semiprimes={len(states)}"
    )

    # --------------------------------------------------------------------------
    # Targets
    # --------------------------------------------------------------------------

    print()
    print("[3] BUILD TARGETS")

    depths = [
        depth_of(s)
        for s in states
    ]

    # --------------------------------------------------------------------------
    # Primary buckets
    # --------------------------------------------------------------------------

    print()
    print("[4] BUILD PRIMARY BUCKETS")

    buckets = build_primary_buckets(
        states
    )

    ambiguous = extract_ambiguous_buckets(
        buckets
    )

    print(
        f"    primary buckets="
        f"{len(buckets)}"
    )

    print(
        f"    ambiguous buckets="
        f"{len(ambiguous)}"
    )

    collision_indices = (
        collect_collision_states(
            ambiguous
        )
    )

    print(
        f"    collision states="
        f"{len(collision_indices)}"
    )

    # --------------------------------------------------------------------------
    # Structure
    # --------------------------------------------------------------------------

    print_primary_structure(
        ambiguous,
        states,
    )

    # --------------------------------------------------------------------------
    # Secondary valuation
    # --------------------------------------------------------------------------

    exact_v2, best_v2 = (
        search_secondary_v2(
            ambiguous,
            states,
            depths,
        )
    )

    # --------------------------------------------------------------------------
    # Odd part
    # --------------------------------------------------------------------------

    exact_odd = (
        search_normalized_odd(
            ambiguous,
            states,
            depths,
        )
    )

    # --------------------------------------------------------------------------
    # n residue
    # --------------------------------------------------------------------------

    exact_residue = (
        search_residue(
            ambiguous,
            states,
            depths,
        )
    )

    # --------------------------------------------------------------------------
    # Factor-side diagnostic
    # --------------------------------------------------------------------------

    search_factor_side(
        ambiguous,
        states,
        depths,
    )

    # --------------------------------------------------------------------------
    # Explicit examples
    # --------------------------------------------------------------------------

    print_collision_examples(
        ambiguous,
        states,
        depths,
    )

    # --------------------------------------------------------------------------
    # Global X,Y comparison
    # --------------------------------------------------------------------------

    collision_xy_analysis(
        ambiguous,
        states,
        depths,
    )

    # --------------------------------------------------------------------------
    # Summary
    # --------------------------------------------------------------------------

    print()
    print("=" * 90)
    print("FINAL SUMMARY")
    print("=" * 90)

    print()
    print(
        f"exact secondary v2 = "
        f"{exact_v2}"
    )

    print(
        f"exact normalized odd = "
        f"{exact_odd}"
    )

    print(
        f"exact low n residue = "
        f"{exact_residue}"
    )

    print()
    print(
        """
The important reduction is:

    primary = v2(n+9)    in frame A
    primary = v2(n+3)    in frame B

Only primary buckets containing multiple depths
are searched further.

Therefore this experiment avoids the O(N * offsets)
scan over the entire semiprime population.

Interpretation:

    EXACT secondary v2
        ->
    depth is recoverable from two n-valuations.

    EXACT normalized odd
        ->
    valuation + odd part is sufficient.

    EXACT low residue
        ->
    valuation plus a finite n residue determines depth.

    otherwise the remaining collision examples reveal
    which factor-side information survives the n-only
    projection.
"""
    )

    print()
    print("=" * 90)
    print("EXPERIMENT 628 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
