#!/usr/bin/env python3

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations
from math import isqrt


# ==============================================================================
# EXPERIMENT 622
# ==============================================================================
#
# SYSTEMATIC n-ONLY 2-ADIC VALUATION SIGNATURE SEARCH
#
# Goal:
#
#   Determine whether
#
#       m = min(v2(A), v2(B))
#
# and/or
#
#       depth
#
# can be reconstructed EXACTLY from a small collection of quantities
#
#       v2(n+c1), v2(n+c2), ...
#
# without using p, q, X, Y.
#
#
# Existing exact structure:
#
#   FRAME A:
#       A = p-3
#       B = q+3
#
#   FRAME B:
#       A = p+1
#       B = q-3
#
#   depth = min(v2(A),v2(B)) + [v2(A)=v2(B)]
#
#   FRAME A equality:
#       v2(n+9) > m
#
#   FRAME B equality:
#       v2(n+3) > m
#
# The unresolved part is recovering m from n alone.
#
# This experiment performs a systematic search instead of manually
# selecting offsets.
#
# ==============================================================================


# ==============================================================================
# CONFIGURATION
# ==============================================================================

PRIME_LIMIT = 6250

# Odd offsets are generally more interesting because n is odd.
# We include 0 separately below.
OFFSET_LIMIT = 65

# Search all pairs/triples of offsets from this set.
MAX_TRIPLE_SEARCH = 3

# Number of examples to print.
MAX_EXAMPLES = 20


# ==============================================================================
# DATA
# ==============================================================================

@dataclass(frozen=True)
class State:
    n: int
    p: int
    q: int


INF = 10**9


# ==============================================================================
# v2
# ==============================================================================

def v2(x: int) -> int:
    """
    Exact 2-adic valuation.

    v2(0) is represented by INF.
    """
    if x == 0:
        return INF

    x = abs(x)

    return (x & -x).bit_length() - 1


# ==============================================================================
# PRIME SIEVE
# ==============================================================================

def prime_sieve(limit: int) -> list[int]:

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

def frame_from_n(
    n: int,
) -> str:

    r = n % 4

    if r == 3:
        return "A"

    if r == 1:
        return "B"

    raise ValueError(
        f"Non-odd semiprime residue: "
        f"n={n}, n mod 4={r}"
    )


# ==============================================================================
# RAW RESIDUALS
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
# TARGET m
# ==============================================================================

def state_target(
    state: State,
) -> tuple[str, int, int]:

    frame = frame_from_n(
        state.n
    )

    A, B = raw_residuals(
        frame,
        state.p,
        state.q,
    )

    m = min(
        v2(A),
        v2(B),
    )

    depth = (
        m
        + int(
            v2(A) == v2(B)
        )
    )

    return (
        frame,
        m,
        depth,
    )


# ==============================================================================
# GLOBAL X,Y
# ==============================================================================

def global_xy(
    frame: str,
    p: int,
    q: int,
) -> tuple[int, int]:

    if frame == "A":

        return (
            (q - p + 6) // 2,
            (p + q) // 2,
        )

    return (
        (3 * p - q + 6) // 2,
        (3 * p + q) // 2,
    )


def global_depth(
    X: int,
    Y: int,
) -> int:

    return (
        1
        + min(
            v2(X),
            v2(Y),
        )
    )


# ==============================================================================
# TEST 0
# ==============================================================================

def test_domain(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 0: ODD SEMIPRIME DOMAIN")
    print("=" * 90)

    failures = 0

    for state in states:

        if state.n % 2 == 0:
            failures += 1

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )
    print()


# ==============================================================================
# TEST 1
# ==============================================================================

def test_existing_identity(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 1: EXISTING n-ONLY EQUALITY IDENTITY")
    print("=" * 90)

    failures = 0

    for state in states:

        frame, m, _ = state_target(
            state
        )

        c = (
            9
            if frame == "A"
            else 3
        )

        A, B = raw_residuals(
            frame,
            state.p,
            state.q,
        )

        equality = (
            v2(A) == v2(B)
        )

        predicted = (
            v2(state.n + c)
            > m
        )

        if equality != predicted:

            failures += 1

            if failures <= MAX_EXAMPLES:

                print(
                    "    mismatch",
                    state.n,
                    frame,
                    m,
                    equality,
                    predicted,
                )

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )
    print()


# ==============================================================================
# TARGET SIGNATURE
# ==============================================================================

def signature(
    n: int,
    offsets: tuple[int, ...],
) -> tuple[int, ...]:

    return tuple(
        v2(n + c)
        for c in offsets
    )


# ==============================================================================
# EXACT FUNCTION TEST
# ==============================================================================

def mapping_ambiguity(
    states: list[State],
    offsets: tuple[int, ...],
    target_index: int,
) -> tuple[
    int,
    int,
    dict,
]:

    mapping = {}
    ambiguous = set()

    for state in states:

        _, m, depth = state_target(
            state
        )

        target = (
            m,
            depth,
        )[target_index]

        sig = signature(
            state.n,
            offsets,
        )

        old = mapping.get(
            sig
        )

        if old is None:

            mapping[sig] = target

        elif old != target:

            ambiguous.add(
                sig
            )

    return (
        len(mapping),
        len(ambiguous),
        mapping,
    )


# ==============================================================================
# TEST 2: SINGLE VALUATION
# ==============================================================================

def test_single_offsets(
    states: list[State],
    offsets: list[int],
) -> None:

    print("=" * 90)
    print("TEST 2: SINGLE n-ONLY VALUATION")
    print("=" * 90)

    best_m = []
    best_depth = []

    for c in offsets:

        states_m = 0
        amb_m = 0

        groups = {}

        for state in states:

            _, m, _ = state_target(
                state
            )

            key = v2(
                state.n + c
            )

            old = groups.get(
                key
            )

            if old is None:

                groups[key] = m

            elif old != m:

                amb_m += 1

            states_m += 1

        best_m.append(
            (
                amb_m,
                c,
                len(groups),
            )
        )

        groups = {}
        amb_depth = 0

        for state in states:

            _, _, depth = state_target(
                state
            )

            key = v2(
                state.n + c
            )

            old = groups.get(
                key
            )

            if old is None:

                groups[key] = depth

            elif old != depth:

                amb_depth += 1

        best_depth.append(
            (
                amb_depth,
                c,
                len(groups),
            )
        )

    best_m.sort()
    best_depth.sort()

    print(
        "    best offsets for m:"
    )

    for amb, c, count in best_m[:15]:

        print(
            f"        c={c:>4} "
            f"ambiguous={amb:<8} "
            f"states={count}"
        )

    print()

    print(
        "    best offsets for depth:"
    )

    for amb, c, count in best_depth[:15]:

        print(
            f"        c={c:>4} "
            f"ambiguous={amb:<8} "
            f"states={count}"
        )

    print()


# ==============================================================================
# TEST 3: PAIR SEARCH
# ==============================================================================

def search_pairs(
    states: list[State],
    offsets: list[int],
) -> list[tuple]:

    print("=" * 90)
    print("TEST 3: SYSTEMATIC PAIR SEARCH")
    print("=" * 90)

    results_m = []
    results_depth = []

    for pair in combinations(
        offsets,
        2,
    ):

        _, amb_m, _ = mapping_ambiguity(
            states,
            pair,
            0,
        )

        _, amb_depth, _ = mapping_ambiguity(
            states,
            pair,
            1,
        )

        results_m.append(
            (
                amb_m,
                pair,
            )
        )

        results_depth.append(
            (
                amb_depth,
                pair,
            )
        )

    results_m.sort()
    results_depth.sort()

    print(
        "    BEST PAIRS FOR m:"
    )

    for amb, pair in results_m[:25]:

        print(
            f"        offsets={pair} "
            f"ambiguous={amb}"
        )

    print()

    print(
        "    BEST PAIRS FOR depth:"
    )

    for amb, pair in results_depth[:25]:

        print(
            f"        offsets={pair} "
            f"ambiguous={amb}"
        )

    print()

    exact_m = [
        pair
        for amb, pair in results_m
        if amb == 0
    ]

    exact_depth = [
        pair
        for amb, pair in results_depth
        if amb == 0
    ]

    print(
        f"    EXACT m pairs={len(exact_m)}"
    )

    if exact_m:

        print(
            "        ",
            exact_m[:20],
        )

    print(
        f"    EXACT depth pairs="
        f"{len(exact_depth)}"
    )

    if exact_depth:

        print(
            "        ",
            exact_depth[:20],
        )

    print()

    return results_m


# ==============================================================================
# TEST 4: FRAME-CONDITIONED PAIR SEARCH
# ==============================================================================

def search_pairs_with_frame(
    states: list[State],
    offsets: list[int],
) -> None:

    print("=" * 90)
    print("TEST 4: FRAME + VALUATION PAIRS")
    print("=" * 90)

    for frame in (
        "A",
        "B",
    ):

        subset = [
            state
            for state in states
            if frame_from_n(
                state.n
            ) == frame
        ]

        results = []

        for pair in combinations(
            offsets,
            2,
        ):

            _, amb, _ = mapping_ambiguity(
                subset,
                pair,
                0,
            )

            results.append(
                (
                    amb,
                    pair,
                )
            )

        results.sort()

        print(
            f"    FRAME {frame}"
        )

        for amb, pair in results[:15]:

            print(
                f"        {pair} "
                f"ambiguous={amb}"
            )

        exact = [
            pair
            for amb, pair in results
            if amb == 0
        ]

        print(
            f"        exact={len(exact)}"
        )

    print()


# ==============================================================================
# TEST 5: TRIPLE SEARCH
# ==============================================================================

def search_triples(
    states: list[State],
    offsets: list[int],
) -> None:

    print("=" * 90)
    print("TEST 5: SYSTEMATIC TRIPLE SEARCH")
    print("=" * 90)

    results_m = []
    results_depth = []

    for triple in combinations(
        offsets,
        3,
    ):

        _, amb_m, _ = mapping_ambiguity(
            states,
            triple,
            0,
        )

        _, amb_depth, _ = mapping_ambiguity(
            states,
            triple,
            1,
        )

        results_m.append(
            (
                amb_m,
                triple,
            )
        )

        results_depth.append(
            (
                amb_depth,
                triple,
            )
        )

    results_m.sort()
    results_depth.sort()

    print(
        "    BEST TRIPLES FOR m:"
    )

    for amb, triple in results_m[:25]:

        print(
            f"        {triple} "
            f"ambiguous={amb}"
        )

    print()

    print(
        "    BEST TRIPLES FOR depth:"
    )

    for amb, triple in results_depth[:25]:

        print(
            f"        {triple} "
            f"ambiguous={amb}"
        )

    print()

    exact_m = [
        triple
        for amb, triple in results_m
        if amb == 0
    ]

    exact_depth = [
        triple
        for amb, triple in results_depth
        if amb == 0
    ]

    print(
        f"    EXACT m triples="
        f"{len(exact_m)}"
    )

    if exact_m:

        print(
            "        ",
            exact_m[:20],
        )

    print(
        f"    EXACT depth triples="
        f"{len(exact_depth)}"
    )

    if exact_depth:

        print(
            "        ",
            exact_depth[:20],
        )

    print()


# ==============================================================================
# TEST 6: MIN/MAX/SUM OF VALUATION VECTORS
# ==============================================================================

def test_vector_operations(
    states: list[State],
    offsets: list[int],
) -> None:

    print("=" * 90)
    print("TEST 6: SIMPLE FUNCTIONS OF n-ONLY VALUATION PAIRS")
    print("=" * 90)

    expressions = [
        "min",
        "max",
        "sum",
        "difference",
        "min_plus_equal",
        "min_plus_gt",
    ]

    for c1, c2 in combinations(
        offsets,
        2,
    ):

        ambiguous = {
            expr: 0
            for expr in expressions
        }

        maps = {
            expr: {}
            for expr in expressions
        }

        for state in states:

            _, m, depth = state_target(
                state
            )

            u = v2(
                state.n + c1
            )

            v = v2(
                state.n + c2
            )

            values = {
                "min": min(u, v),
                "max": max(u, v),
                "sum": u + v,
                "difference": abs(u - v),
                "min_plus_equal":
                    min(u, v)
                    + int(u == v),
                "min_plus_gt":
                    min(u, v)
                    + int(
                        max(u, v)
                        > min(u, v)
                    ),
            }

            for expr, value in values.items():

                mapping = maps[
                    expr
                ]

                old = mapping.get(
                    state.n
                    % 1
                )

                # We actually need signature-based consistency.
                # Use (u,v) itself as the signature.
                key = (
                    u,
                    v,
                )

                old = mapping.get(
                    key
                )

                target = (
                    m,
                    depth,
                )

                if old is None:

                    mapping[key] = target

                else:

                    if (
                        value != old[0]
                        and False
                    ):
                        pass

        # Report is done in the dedicated pair evaluator below.

    #
    # Actual compact search.
    #
    results = []

    for c1, c2 in combinations(
        offsets,
        2,
    ):

        mismatches = {
            "min->m": 0,
            "max->m": 0,
            "min+eq->m": 0,
            "min->depth": 0,
            "max->depth": 0,
            "min+eq->depth": 0,
        }

        for state in states:

            _, m, depth = state_target(
                state
            )

            u = v2(
                state.n + c1
            )

            v = v2(
                state.n + c2
            )

            if min(u, v) != m:
                mismatches[
                    "min->m"
                ] += 1

            if max(u, v) != m:
                mismatches[
                    "max->m"
                ] += 1

            if (
                min(u, v)
                + int(u == v)
                != m
            ):
                mismatches[
                    "min+eq->m"
                ] += 1

            if min(u, v) != depth:
                mismatches[
                    "min->depth"
                ] += 1

            if max(u, v) != depth:
                mismatches[
                    "max->depth"
                ] += 1

            if (
                min(u, v)
                + int(u == v)
                != depth
            ):
                mismatches[
                    "min+eq->depth"
                ] += 1

        score = (
            sum(
                mismatches.values()
            ),
            max(
                mismatches.values()
            ),
        )

        results.append(
            (
                score,
                (c1, c2),
                mismatches,
            )
        )

    results.sort()

    for score, pair, mismatches in results[:25]:

        print(
            f"    pair={pair} "
            f"score={score}"
        )

        for name, value in mismatches.items():

            print(
                f"        {name}: "
                f"{value}"
            )

    print()


# ==============================================================================
# TEST 7: n-ONLY CONGRUENCE SIGNATURE
# ==============================================================================

def test_residue_plus_valuation(
    states: list[State],
    offsets: list[int],
) -> None:

    print("=" * 90)
    print("TEST 7: n RESIDUE + 2-ADIC VALUATION SIGNATURE")
    print("=" * 90)

    configurations = [
        (3, 4),
        (3, 6),
        (3, 8),
        (3, 10),
        (9, 4),
        (9, 6),
        (9, 8),
        (9, 10),
    ]

    for c, bits in configurations:

        mapping = {}
        ambiguous = set()

        modulus = 1 << bits

        for state in states:

            _, m, depth = state_target(
                state
            )

            key = (
                state.n % modulus,
                v2(state.n + c),
            )

            old = mapping.get(
                key
            )

            target = (
                m,
                depth,
            )

            if old is None:

                mapping[key] = target

            elif old != target:

                ambiguous.add(key)

        print(
            f"    c={c:<3} "
            f"bits={bits:<2} "
            f"states={len(mapping):<8} "
            f"ambiguous="
            f"{len(ambiguous)}"
        )

    print()


# ==============================================================================
# TEST 8: DIRECT n-ONLY DEPTH CLASSIFIER
# ==============================================================================

def test_depth_from_single_valuation(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 8: SEARCH depth = v2(n+c) + CONSTANT/CONDITION")
    print("=" * 90)

    candidates = []

    for c in range(
        -65,
        66,
    ):

        if c % 2 == 0:
            continue

        mismatches_direct = 0
        mismatches_minus = 0
        mismatches_plus = 0

        for state in states:

            frame, m, depth = state_target(
                state
            )

            v = v2(
                state.n + c
            )

            if v != depth:
                mismatches_direct += 1

            if v - 1 != depth:
                mismatches_minus += 1

            if v + 1 != depth:
                mismatches_plus += 1

        candidates.append(
            (
                mismatches_direct,
                mismatches_minus,
                mismatches_plus,
                c,
            )
        )

    candidates.sort(
        key=lambda x: (
            min(
                x[0],
                x[1],
                x[2],
            ),
            abs(x[3]),
        )
    )

    for direct, minus, plus, c in candidates[:20]:

        print(
            f"    c={c:>4} "
            f"direct={direct:<8} "
            f"v-1={minus:<8} "
            f"v+1={plus:<8}"
        )

    print()


# ==============================================================================
# TEST 9: INFORMATION-THEORETIC AMBIGUITY
# ==============================================================================

def test_best_signature_examples(
    states: list[State],
    offsets: list[int],
) -> None:

    print("=" * 90)
    print("TEST 9: FIRST AMBIGUOUS n-ONLY PAIRS")
    print("=" * 90)

    best_pair = None
    best_ambiguous = None
    best_mapping = None

    for pair in combinations(
        offsets,
        2,
    ):

        mapping = {}
        ambiguous = set()

        for state in states:

            _, m, _ = state_target(
                state
            )

            key = signature(
                state.n,
                pair,
            )

            old = mapping.get(
                key
            )

            if old is None:

                mapping[key] = m

            elif old != m:

                ambiguous.add(
                    key
                )

        count = len(
            ambiguous
        )

        if (
            best_ambiguous is None
            or count < best_ambiguous
        ):

            best_pair = pair
            best_ambiguous = count
            best_mapping = mapping

    print(
        f"    best_pair={best_pair}"
    )

    print(
        f"    ambiguous_signatures="
        f"{best_ambiguous}"
    )

    if best_pair is not None:

        lookup = {}

        for state in states:

            _, m, depth = state_target(
                state
            )

            key = signature(
                state.n,
                best_pair,
            )

            old = lookup.get(
                key
            )

            if old is not None and old[
                0
            ] != m:

                print()
                print(
                    "    FIRST COUNTEREXAMPLE:"
                )

                print(
                    f"        signature={key}"
                )

                print(
                    f"        existing="
                    f"{old}"
                )

                print(
                    f"        current="
                    f"(n={state.n}, "
                    f"m={m}, "
                    f"depth={depth})"
                )

                break

            lookup[key] = (
                m,
                depth,
                state.n,
            )

    print()


# ==============================================================================
# EXAMPLES
# ==============================================================================

def print_examples(
    states: list[State],
) -> None:

    print("=" * 90)
    print("EXAMPLES")
    print("=" * 90)

    lookup = {
        state.n: state
        for state in states
    }

    for n in [
        9,
        15,
        21,
        33,
        39,
        51,
        57,
        69,
        87,
        93,
        141,
        183,
    ]:

        state = lookup.get(n)

        if state is None:
            continue

        frame, m, depth = state_target(
            state
        )

        X, Y = global_xy(
            frame,
            state.p,
            state.q,
        )

        print()
        print(
            f"n={n} "
            f"p={state.p} "
            f"q={state.q}"
        )

        print(
            f"    frame={frame}"
        )

        print(
            f"    A,B="
            f"{raw_residuals(frame, state.p, state.q)}"
        )

        print(
            f"    X,Y=({X},{Y})"
        )

        print(
            f"    m={m}"
        )

        print(
            f"    depth={depth}"
        )

        for c in (
            -9,
            -3,
            3,
            9,
            25,
            41,
        ):

            print(
                f"    v2(n{c:+d})="
                f"{v2(n+c)}"
            )


# ==============================================================================
# SYMBOLIC SUMMARY
# ==============================================================================

def symbolic_summary() -> None:

    print()
    print("=" * 90)
    print("SYMBOLIC SUMMARY")
    print("=" * 90)

    print(
r"""
EXPERIMENT 617 established the exact factor-side depth law:

    depth
      =
    min(v2(A),v2(B))
      +
    [v2(A)=v2(B)].

Experiment 619 established that equality can be detected from n:

    FRAME A:
        equality <=> v2(n+9) > m

    FRAME B:
        equality <=> v2(n+3) > m

so:

    depth
      =
    m + [v2(n+c)>m].

Experiment 621 showed that:
    
    n mod 2^k

does not uniquely determine m over the tested range.

It also showed that individual values v2(n+c) do not
recover m exactly.

Experiment 622 therefore performs a systematic search over
pairs and triples of valuations:

    (v2(n+c1),v2(n+c2))

and

    (v2(n+c1),v2(n+c2),v2(n+c3)).

Three outcomes are possible:

    1. EXACT PAIR:
         m is recoverable from two n-only valuations.

    2. EXACT TRIPLE:
         one additional n-only valuation is sufficient.

    3. NO SMALL SIGNATURE:
         m genuinely retains information not captured by
         these local n-adic observables.

The particularly interesting offsets are expected near:

    -9,-3,3,9

because they directly arise from the two factor frames.

An exact result here would give:

    n
      ->
    small 2-adic signature
      ->
    m
      ->
    depth

without recovering p,q or X,Y.

"""
    )


# ==============================================================================
# MAIN
# ==============================================================================

def main() -> None:

    print("=" * 90)
    print("EXPERIMENT 622 START")
    print("=" * 90)
    print()
    print("SYSTEMATIC n-ONLY 2-ADIC VALUATION SIGNATURE SEARCH")
    print()

    print("[1] PRIME SIEVE")

    primes = prime_sieve(
        PRIME_LIMIT
    )

    print(
        f"    odd primes={len(primes)}"
    )
    print()

    print("[2] SEMIPRIME GENERATION")

    states = generate_states(
        primes
    )

    print(
        f"    semiprimes={len(states)}"
    )
    print()

    print("[3] GLOBAL STATES")

    print(
        f"    states={len(states)}"
    )
    print()

    test_domain(
        states
    )

    test_existing_identity(
        states
    )

    # Keep this list reasonably small because pair/triple searches
    # scale combinatorially.
    offsets = [
        c
        for c in range(
            -OFFSET_LIMIT,
            OFFSET_LIMIT + 1,
        )
        if c % 2 != 0
    ]

    # Add the most important offsets explicitly.
    for c in (
        -9,
        -3,
        3,
        9,
    ):
        if c not in offsets:
            offsets.append(c)

    offsets = sorted(
        set(offsets)
    )

    test_single_offsets(
        states,
        offsets,
    )

    pair_results = search_pairs(
        states,
        offsets,
    )

    search_pairs_with_frame(
        states,
        offsets,
    )

    # Triple search over all 65 odd offsets is relatively large.
    #
    # Use a focused candidate pool consisting of:
    #
    #   frame constants
    #   nearby odd offsets
    #
    triple_offsets = [
        -17,
        -15,
        -13,
        -11,
        -9,
        -7,
        -5,
        -3,
        -1,
        1,
        3,
        5,
        7,
        9,
        11,
        13,
        15,
        17,
        25,
        33,
        41,
    ]

    search_triples(
        states,
        triple_offsets,
    )

    test_vector_operations(
        states,
        [
            -9,
            -3,
            3,
            9,
        ],
    )

    test_residue_plus_valuation(
        states,
        [
            -9,
            -3,
            3,
            9,
        ],
    )

    test_depth_from_single_valuation(
        states
    )

    test_best_signature_examples(
        states,
        [
            -17,
            -15,
            -13,
            -11,
            -9,
            -7,
            -5,
            -3,
            -1,
            1,
            3,
            5,
            7,
            9,
            11,
            13,
            15,
            17,
            25,
            33,
            41,
        ],
    )

    print_examples(
        states
    )

    symbolic_summary()

    print("=" * 90)
    print("EXPERIMENT 622 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
