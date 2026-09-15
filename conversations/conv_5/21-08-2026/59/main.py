#!/usr/bin/env python3

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from math import isqrt


# ==============================================================================
# EXPERIMENT 621
# ==============================================================================
#
# n-ONLY 2-ADIC DEPTH / m RECOVERY
#
# Experiment 620 established exactly:
#
#     equality =
#         [v2(A) = v2(B)]
#
# and:
#
#     FRAME A:
#         equality <=> v2(n+9) > m
#
#     FRAME B:
#         equality <=> v2(n+3) > m
#
# where:
#
#     m = min(v2(A), v2(B)).
#
# Therefore:
#
#     depth =
#         m + [v2(n+c) > m]
#
# with:
#
#     c=9 for A
#     c=3 for B.
#
# However m itself was NOT recoverable from n mod 2^k alone in the
# tested lookup table.
#
# This experiment searches for a more structured n-only description.
#
# Tests:
#
#   1. equality collapse
#   2. exact depth formula
#   3. global-depth agreement
#   4. inspect v2(n+c)-m
#   5. exact n-residue ambiguity
#   6. n-only valuation candidates
#   7. combinations of v2(n+c) and n mod powers of two
#   8. direct depth signatures
#   9. search simple functions of v2(n+c)
#  10. examples
#
# ==============================================================================


# ==============================================================================
# CONFIG
# ==============================================================================

PRIME_LIMIT = 6250
MAX_BITS = 22

EXAMPLES = [
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
    111,
    123,
    129,
    141,
    159,
    177,
    183,
    201,
    213,
    219,
]


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

    for p in range(
        2,
        isqrt(limit) + 1,
    ):

        if sieve[p]:

            start = p * p

            count = (
                (limit - start) // p
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

    states: list[State] = []

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

def frame_from_n(n: int) -> str:

    r = n % 4

    if r == 3:
        return "A"

    if r == 1:
        return "B"

    raise ValueError(
        f"Unexpected odd semiprime residue: "
        f"n={n}, n mod 4={r}"
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

    if frame == "B":

        return (
            p + 1,
            q - 3,
        )

    raise ValueError(frame)


# ==============================================================================
# FRAME OFFSET
# ==============================================================================

def frame_offset(
    frame: str,
) -> int:

    return 9 if frame == "A" else 3


# ==============================================================================
# GLOBAL X,Y
# ==============================================================================

def global_xy(
    frame: str,
    p: int,
    q: int,
) -> tuple[int, int]:

    if frame == "A":

        X = (
            q
            - p
            + 6
        ) // 2

        Y = (
            p
            + q
        ) // 2

        return X, Y

    X = (
        3 * p
        - q
        + 6
    ) // 2

    Y = (
        3 * p
        + q
    ) // 2

    return X, Y


# ==============================================================================
# GLOBAL DEPTH
# ==============================================================================

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
# RAW m
# ==============================================================================

def raw_m(
    A: int,
    B: int,
) -> int:

    return min(
        v2(A),
        v2(B),
    )


# ==============================================================================
# RAW FACTOR DEPTH
# ==============================================================================

def factor_depth(
    A: int,
    B: int,
) -> int:

    alpha = v2(A)
    beta = v2(B)

    return (
        min(alpha, beta)
        + int(alpha == beta)
    )


# ==============================================================================
# TEST 1
# ==============================================================================

def test_equality_collapse(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 1: EQUALITY FLAG FROM n")
    print("=" * 90)

    failures = 0

    for state in states:

        frame = frame_from_n(
            state.n
        )

        A, B = raw_residuals(
            frame,
            state.p,
            state.q,
        )

        m = raw_m(A, B)

        actual = (
            v2(A)
            ==
            v2(B)
        )

        c = frame_offset(frame)

        predicted = (
            v2(state.n + c) > m
        )

        if actual != predicted:

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch",
                    f"n={state.n}",
                    f"frame={frame}",
                    f"A={A}",
                    f"B={B}",
                    f"m={m}",
                    f"v2(n+c)={v2(state.n+c)}",
                    f"actual={actual}",
                    f"predicted={predicted}",
                )

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )
    print()


# ==============================================================================
# TEST 2
# ==============================================================================

def test_direct_depth(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 2: DIRECT n-BASED DEPTH")
    print("=" * 90)

    failures = 0

    for state in states:

        frame = frame_from_n(
            state.n
        )

        A, B = raw_residuals(
            frame,
            state.p,
            state.q,
        )

        m = raw_m(A, B)

        c = frame_offset(frame)

        predicted = (
            m
            + int(
                v2(state.n + c) > m
            )
        )

        actual = factor_depth(
            A,
            B,
        )

        if actual != predicted:

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch",
                    f"n={state.n}",
                    f"actual={actual}",
                    f"predicted={predicted}",
                )

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )
    print()


# ==============================================================================
# TEST 3
# ==============================================================================

def test_global_depth(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 3: n-BASED DEPTH == GLOBAL DEPTH")
    print("=" * 90)

    failures = 0

    for state in states:

        frame = frame_from_n(
            state.n
        )

        A, B = raw_residuals(
            frame,
            state.p,
            state.q,
        )

        m = raw_m(A, B)

        c = frame_offset(frame)

        predicted = (
            m
            + int(
                v2(state.n + c) > m
            )
        )

        X, Y = global_xy(
            frame,
            state.p,
            state.q,
        )

        actual = global_depth(
            X,
            Y,
        )

        if actual != predicted:

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch",
                    f"n={state.n}",
                    f"frame={frame}",
                    f"X={X}",
                    f"Y={Y}",
                    f"m={m}",
                    f"actual={actual}",
                    f"predicted={predicted}",
                )

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )
    print()


# ==============================================================================
# TEST 4
# ==============================================================================

def test_offset_structure(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 4: v2(n+c) - m")
    print("=" * 90)

    data = {
        "A": {
            "equal": Counter(),
            "unequal": Counter(),
        },
        "B": {
            "equal": Counter(),
            "unequal": Counter(),
        },
    }

    failures = 0

    for state in states:

        frame = frame_from_n(
            state.n
        )

        A, B = raw_residuals(
            frame,
            state.p,
            state.q,
        )

        alpha = v2(A)
        beta = v2(B)

        m = min(
            alpha,
            beta,
        )

        relation = (
            "equal"
            if alpha == beta
            else "unequal"
        )

        delta = (
            v2(
                state.n
                + frame_offset(frame)
            )
            - m
        )

        data[frame][relation][
            delta
        ] += 1

        expected = (
            (delta > 0)
            ==
            (relation == "equal")
        )

        if not expected:
            failures += 1

    for frame in (
        "A",
        "B",
    ):

        print(
            f"    FRAME {frame}"
        )

        print(
            "        equal:",
            dict(
                sorted(
                    data[frame][
                        "equal"
                    ].items()
                )
            ),
        )

        print(
            "        unequal:",
            dict(
                sorted(
                    data[frame][
                        "unequal"
                    ].items()
                )
            ),
        )

    print(
        f"failures={failures}"
    )
    print()


# ==============================================================================
# TEST 5
# ==============================================================================

def test_exact_n_residue(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 5: m FROM n MOD 2^k")
    print("=" * 90)

    for bits in range(
        2,
        MAX_BITS + 1,
    ):

        modulus = 1 << bits

        mapping: dict[int, int] = {}
        ambiguous = set()

        for state in states:

            frame = frame_from_n(
                state.n
            )

            A, B = raw_residuals(
                frame,
                state.p,
                state.q,
            )

            m = raw_m(A, B)

            key = (
                state.n
                % modulus
            )

            old = mapping.get(key)

            if old is None:

                mapping[key] = m

            elif old != m:

                ambiguous.add(key)

        print(
            f"    bits={bits:<2} "
            f"mod={modulus:<9} "
            f"states={len(mapping):<8} "
            f"ambiguous="
            f"{len(ambiguous)}"
        )

    print()


# ==============================================================================
# TEST 6
# ==============================================================================

def test_v2_n_candidates(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 6: n-ONLY VALUATION CANDIDATES")
    print("=" * 90)

    # Candidate odd offsets.
    offsets = range(
        -63,
        64,
        2,
    )

    result = []

    for c in offsets:

        exact = 0
        delta_counts = Counter()

        for state in states:

            frame = frame_from_n(
                state.n
            )

            A, B = raw_residuals(
                frame,
                state.p,
                state.q,
            )

            m = raw_m(A, B)

            vn = v2(
                state.n + c
            )

            if vn == m:
                exact += 1

            delta_counts[
                vn - m
            ] += 1

        # Number of distinct offsets from m.
        distinct = len(
            delta_counts
        )

        result.append(
            (
                distinct,
                -exact,
                c,
                delta_counts,
            )
        )

    result.sort(
        key=lambda x: (
            x[0],
            x[1],
            abs(x[2]),
        )
    )

    for distinct, neg_exact, c, counts in result[:20]:

        print(
            f"    c={c:>4} "
            f"distinct_delta={distinct:<3} "
            f"exact_m={-neg_exact:<8} "
            f"distribution="
            f"{dict(sorted(counts.items()))}"
        )

    print()


# ==============================================================================
# TEST 7
# ==============================================================================

def test_composite_n_signature(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 7: COMPOSITE n-ONLY SIGNATURE SEARCH")
    print("=" * 90)

    # Candidate signature:
    #
    #     (
    #         frame,
    #         v2(n+c),
    #         n mod 2^b
    #     )
    #
    # We are searching for a small signature that uniquely determines m.
    #
    # Since frame is already n mod 4, the interesting information is
    # whether adding a small c gives a useful second valuation.
    #
    candidate_offsets = [
        3,
        9,
        -3,
        -9,
        1,
        5,
        7,
        11,
        13,
        15,
    ]

    for c in candidate_offsets:

        for bits in (
            4,
            5,
            6,
            7,
            8,
            9,
            10,
            12,
        ):

            modulus = 1 << bits

            mapping = {}
            ambiguous = set()

            for state in states:

                frame = frame_from_n(
                    state.n
                )

                A, B = raw_residuals(
                    frame,
                    state.p,
                    state.q,
                )

                m = raw_m(
                    A,
                    B,
                )

                signature = (
                    frame,
                    v2(
                        state.n + c
                    ),
                    state.n
                    % modulus,
                )

                old = mapping.get(
                    signature
                )

                if old is None:

                    mapping[
                        signature
                    ] = m

                elif old != m:

                    ambiguous.add(
                        signature
                    )

            if not ambiguous:

                print(
                    "    EXACT:",
                    f"c={c}",
                    f"bits={bits}",
                    f"states={len(mapping)}"
                )

                break

    print()


# ==============================================================================
# TEST 8
# ==============================================================================

def test_depth_signature(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 8: DIRECT DEPTH FROM n MOD 2^k")
    print("=" * 90)

    for bits in range(
        2,
        MAX_BITS + 1,
    ):

        modulus = 1 << bits

        mapping = {}
        ambiguous = set()

        for state in states:

            frame = frame_from_n(
                state.n
            )

            X, Y = global_xy(
                frame,
                state.p,
                state.q,
            )

            depth = global_depth(
                X,
                Y,
            )

            key = (
                state.n
                % modulus
            )

            old = mapping.get(
                key
            )

            if old is None:

                mapping[key] = depth

            elif old != depth:

                ambiguous.add(key)

        print(
            f"    bits={bits:<2} "
            f"mod={modulus:<9} "
            f"states={len(mapping):<8} "
            f"ambiguous="
            f"{len(ambiguous)}"
        )

    print()


# ==============================================================================
# TEST 9
# ==============================================================================

def test_depth_from_v2_offsets(
    states: list[State],
) -> None:

    print("=" * 90)
    print("TEST 9: DEPTH FROM SMALL n-ONLY VALUATION VECTORS")
    print("=" * 90)

    candidate_vectors = [
        (3, 9),
        (3, -3),
        (9, -3),
        (3, 5),
        (3, 7),
        (3, 11),
        (9, 1),
        (9, 5),
        (9, 13),
    ]

    for c1, c2 in candidate_vectors:

        mismatches = 0

        # Search simple expressions:
        #
        #   min(v2(n+c1),v2(n+c2))
        #   max(...)
        #   min + indicator(...)
        #
        for state in states:

            frame = frame_from_n(
                state.n
            )

            A, B = raw_residuals(
                frame,
                state.p,
                state.q,
            )

            X, Y = global_xy(
                frame,
                state.p,
                state.q,
            )

            actual = global_depth(
                X,
                Y,
            )

            u = v2(
                state.n + c1
            )

            v = v2(
                state.n + c2
            )

            candidate1 = min(
                u,
                v,
            )

            candidate2 = (
                candidate1
                + int(
                    u == v
                )
            )

            candidate3 = max(
                u,
                v,
            )

            if actual not in (
                candidate1,
                candidate2,
                candidate3,
            ):

                mismatches += 1

        print(
            f"    offsets=({c1},{c2}) "
            f"mismatches={mismatches}"
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

    for n in EXAMPLES:

        state = lookup.get(n)

        if state is None:
            continue

        frame = frame_from_n(n)

        A, B = raw_residuals(
            frame,
            state.p,
            state.q,
        )

        X, Y = global_xy(
            frame,
            state.p,
            state.q,
        )

        alpha = v2(A)
        beta = v2(B)

        m = min(
            alpha,
            beta,
        )

        c = frame_offset(
            frame
        )

        vn = v2(
            n + c
        )

        equality = (
            alpha == beta
        )

        predicted = (
            m
            + int(vn > m)
        )

        actual = global_depth(
            X,
            Y,
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
            f"    X={X} "
            f"Y={Y}"
        )

        print(
            f"    A={A} "
            f"B={B}"
        )

        print(
            f"    v2(A)="
            f"{alpha if alpha < INF else 'inf'}"
        )

        print(
            f"    v2(B)="
            f"{beta if beta < INF else 'inf'}"
        )

        print(
            f"    m={m}"
        )

        print(
            f"    c={c}"
        )

        print(
            f"    n+c={n+c}"
        )

        print(
            f"    v2(n+c)="
            f"{vn if vn < INF else 'inf'}"
        )

        print(
            f"    equality={equality}"
        )

        print(
            f"    predicted_depth="
            f"{predicted}"
        )

        print(
            f"    actual_depth="
            f"{actual}"
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
EXPERIMENT 620 established the exact identity:

    equality
      =
    [v2(A)=v2(B)]

and:

    FRAME A:
        equality
          <=>
        v2(n+9) > m

    FRAME B:
        equality
          <=>
        v2(n+3) > m

where:

    m=min(v2(A),v2(B)).

Therefore:

    depth
      =
    m + [v2(n+c)>m].

But TEST 5 of Experiment 620 showed that m is NOT uniquely
determined by n mod 2^k for the tested range.

So Experiment 621 attacks a different possibility:

    perhaps m is not a simple residue,
    but becomes recoverable from a small vector of n-only
    2-adic observables.

For example:

    v2(n+3)
    v2(n+9)
    n mod 2^k

or combinations of these.

The target is:

    (n-only signature)
        ->
    m
        ->
    depth.

An exact success would eliminate p,q and X,Y entirely
from the depth calculation.

A failure would establish that the remaining m value
contains information not represented by the tested
small n-only signatures.
"""
    )

    print()


# ==============================================================================
# MAIN
# ==============================================================================

def main() -> None:

    print("=" * 90)
    print("EXPERIMENT 621 START")
    print("=" * 90)
    print()
    print("n-ONLY 2-ADIC DEPTH / m RECOVERY")
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

    test_equality_collapse(
        states
    )

    test_direct_depth(
        states
    )

    test_global_depth(
        states
    )

    test_offset_structure(
        states
    )

    test_exact_n_residue(
        states
    )

    test_v2_n_candidates(
        states
    )

    test_composite_n_signature(
        states
    )

    test_depth_signature(
        states
    )

    test_depth_from_v2_offsets(
        states
    )

    print_examples(
        states
    )

    symbolic_summary()

    print("=" * 90)
    print("EXPERIMENT 621 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()