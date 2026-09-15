#!/usr/bin/env python3

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from random import Random
from math import isqrt


# ==============================================================================
# EXPERIMENT 624
# ==============================================================================
#
# SMALL n-ONLY 2-ADIC SIGNATURE SEARCH
#
# Experiment 623 showed:
#
#   no tested pair
#
#       (v2(n+c1), v2(n+c2))
#
# uniquely determines m.
#
# This experiment asks the next natural question:
#
#   Can a SMALL VALUATION SIGNATURE of 3 or 4 offsets determine
#   m or depth?
#
# We use:
#
#   1. a deterministic sample for fast ranking
#   2. structural/prioritized offsets
#   3. 3-tuples first
#   4. only the best tuples go to exact verification
#   5. 4-tuples only around the best 3-tuple candidates
#
# IMPORTANT:
#
#   We do NOT perform an exhaustive 4-combination search over the
#   complete offset pool.
#
# ==============================================================================


# ==============================================================================
# CONFIGURATION
# ==============================================================================

PRIME_LIMIT = 6250

SAMPLE_SIZE = 30000

TOP_TRIPLES = 40
TOP_QUADRUPLES = 25

RANDOM_SEED = 624

# Offset pool.
#
# Include:
#
#   - theoretically important constants
#   - offsets that performed well in Exp. 622/623
#   - local odd/even neighbors
#
# The pool is deliberately modest.
#
OFFSET_POOL = sorted(
    {
        -63,
        -59,
        -55,
        -51,
        -47,
        -43,
        -39,
        -35,
        -31,
        -27,
        -23,
        -19,
        -15,
        -11,
        -9,
        -7,
        -3,
        1,
        3,
        5,
        7,
        9,
        11,
        13,
        17,
        21,
        25,
        29,
        33,
        37,
        41,
        45,
        53,
        57,
        61,
    }
)

INF = 10**9


# ==============================================================================
# DATA
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

    return (
        (x & -x).bit_length()
        - 1
    )


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
        ] = (
            b"\x00"
            * count
        )

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

    r = n & 3

    if r == 3:
        return "A"

    if r == 1:
        return "B"

    raise ValueError(
        f"Expected odd semiprime, "
        f"got n={n}, n mod 4={r}"
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
# TARGETS
# ==============================================================================

def build_targets(
    states: list[State],
) -> tuple[list[int], list[int], list[int]]:

    m_values = []
    depth_values = []
    equality_values = []

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

        equality = int(
            alpha == beta
        )

        depth = (
            m + equality
        )

        m_values.append(m)
        depth_values.append(depth)
        equality_values.append(
            equality
        )

    return (
        m_values,
        depth_values,
        equality_values,
    )


# ==============================================================================
# VALUATION COLUMNS
# ==============================================================================

def build_columns(
    states: list[State],
    offsets: list[int],
) -> dict[int, list[int]]:

    print(
        "    building valuation columns..."
    )

    ns = [
        s.n
        for s in states
    ]

    columns = {}

    for c in offsets:

        columns[c] = [
            v2(n + c)
            for n in ns
        ]

    return columns


# ==============================================================================
# SAMPLE
# ==============================================================================

def make_sample(
    total: int,
    requested: int,
) -> list[int]:

    if requested >= total:

        return list(
            range(total)
        )

    rng = Random(
        RANDOM_SEED
    )

    result = rng.sample(
        range(total),
        requested,
    )

    result.sort()

    return result


# ==============================================================================
# SIGNATURE ITERATOR
# ==============================================================================

def signature_at(
    pair: tuple[int, ...],
    columns: dict[int, list[int]],
    i: int,
) -> tuple[int, ...]:

    return tuple(
        columns[c][i]
        for c in pair
    )


# ==============================================================================
# RANK A SIGNATURE
# ==============================================================================

def rank_signature(
    offsets_tuple: tuple[int, ...],
    columns: dict[int, list[int]],
    target: list[int],
    indices: list[int],
) -> tuple[int, int, int]:

    mapping = {}

    ambiguous = set()

    for i in indices:

        key = signature_at(
            offsets_tuple,
            columns,
            i,
        )

        value = target[i]

        previous = mapping.get(
            key
        )

        if previous is None:

            mapping[key] = value

        elif previous != value:

            ambiguous.add(
                key
            )

    return (
        len(ambiguous),
        len(mapping),
        len(indices),
    )


# ==============================================================================
# EXACT SIGNATURE
# ==============================================================================

def exact_signature(
    offsets_tuple: tuple[int, ...],
    columns: dict[int, list[int]],
    target: list[int],
) -> tuple[int, int, tuple | None]:

    mapping = {}

    ambiguous = set()

    first_counterexample = None

    for i in range(
        len(target)
    ):

        key = signature_at(
            offsets_tuple,
            columns,
            i,
        )

        value = target[i]

        previous = mapping.get(
            key
        )

        if previous is None:

            mapping[key] = value

        elif previous != value:

            if key not in ambiguous:

                ambiguous.add(
                    key
                )

                if (
                    first_counterexample
                    is None
                ):

                    first_counterexample = (
                        i,
                        key,
                        previous,
                        value,
                    )

    return (
        len(ambiguous),
        len(mapping),
        first_counterexample,
    )


# ==============================================================================
# TRIPLE SEARCH
# ==============================================================================

def rank_triples(
    columns,
    offsets,
    target,
    sample_indices,
) -> list[tuple[int, ...]]:

    print()
    print(
        "=" * 90
    )
    print(
        "TEST 3: FAST TRIPLE SEARCH"
    )
    print(
        "=" * 90
    )

    ranked = []

    tested = 0

    for triple in combinations(
        offsets,
        3,
    ):

        result = rank_signature(
            triple,
            columns,
            target,
            sample_indices,
        )

        ambiguous, groups, _ = result

        ranked.append(
            (
                ambiguous,
                -groups,
                triple,
            )
        )

        tested += 1

    ranked.sort()

    print(
        f"    triples tested={tested}"
    )

    print(
        "    best triples:"
    )

    for ambiguous, neg_groups, triple in ranked[:30]:

        print(
            f"        {triple} "
            f"ambiguous={ambiguous} "
            f"groups={-neg_groups}"
        )

    return [
        triple
        for _, _, triple
        in ranked[
            :TOP_TRIPLES
        ]
    ]


# ==============================================================================
# EXACT TRIPLE TEST
# ==============================================================================

def exact_triples(
    triples,
    columns,
    m_values,
    depth_values,
):

    print()
    print(
        "=" * 90
    )
    print(
        "TEST 4: EXACT TRIPLE VERIFICATION"
    )
    print(
        "=" * 90
    )

    exact_m = []
    exact_depth = []

    for triple in triples:

        m_ambiguous, m_groups, m_ce = (
            exact_signature(
                triple,
                columns,
                m_values,
            )
        )

        d_ambiguous, d_groups, d_ce = (
            exact_signature(
                triple,
                columns,
                depth_values,
            )
        )

        print()
        print(
            f"    triple={triple}"
        )

        print(
            f"        m: "
            f"ambiguous={m_ambiguous} "
            f"groups={m_groups}"
        )

        print(
            f"        depth: "
            f"ambiguous={d_ambiguous} "
            f"groups={d_groups}"
        )

        if m_ambiguous == 0:

            exact_m.append(
                triple
            )

        if d_ambiguous == 0:

            exact_depth.append(
                triple
            )

        if (
            m_ambiguous
            and m_ce
        ):

            (
                index,
                signature,
                old_value,
                new_value,
            ) = m_ce

            print(
                "        first m collision:"
            )

            print(
                f"            index={index}"
            )

            print(
                f"            signature={signature}"
            )

            print(
                f"            targets="
                f"{old_value},{new_value}"
            )

    return (
        exact_m,
        exact_depth,
    )


# ==============================================================================
# QUADRUPLE SEARCH AROUND GOOD TRIPLES
# ==============================================================================

def rank_quads_from_triples(
    triples,
    offsets,
    columns,
    target,
    sample_indices,
):

    print()
    print(
        "=" * 90
    )
    print(
        "TEST 5: LOCAL QUADRUPLE SEARCH"
    )
    print(
        "=" * 90
    )

    ranked = []

    seen = set()

    # For each strong triple, attach a fourth offset.
    #
    # This is dramatically cheaper than searching all 4-combinations.

    for triple in triples:

        remaining = [
            c
            for c in offsets
            if c not in triple
        ]

        for c in remaining:

            quad = tuple(
                sorted(
                    triple + (c,)
                )
            )

            if quad in seen:
                continue

            seen.add(
                quad
            )

            ambiguous, groups, _ = (
                rank_signature(
                    quad,
                    columns,
                    target,
                    sample_indices,
                )
            )

            ranked.append(
                (
                    ambiguous,
                    -groups,
                    quad,
                )
            )

    ranked.sort()

    print(
        f"    quadruples tested="
        f"{len(ranked)}"
    )

    print(
        "    best quadruples:"
    )

    for ambiguous, neg_groups, quad in ranked[:30]:

        print(
            f"        {quad} "
            f"ambiguous={ambiguous} "
            f"groups={-neg_groups}"
        )

    return [
        quad
        for _, _, quad in ranked[
            :TOP_QUADRUPLES
        ]
    ]


# ==============================================================================
# EXACT QUADRUPLES
# ==============================================================================

def exact_quads(
    quads,
    columns,
    m_values,
    depth_values,
):

    print()
    print(
        "=" * 90
    )
    print(
        "TEST 6: EXACT QUADRUPLE VERIFICATION"
    )
    print(
        "=" * 90
    )

    exact_m = []
    exact_depth = []

    for quad in quads:

        m_ambiguous, m_groups, m_ce = (
            exact_signature(
                quad,
                columns,
                m_values,
            )
        )

        d_ambiguous, d_groups, d_ce = (
            exact_signature(
                quad,
                columns,
                depth_values,
            )
        )

        print()
        print(
            f"    quad={quad}"
        )

        print(
            f"        m: "
            f"ambiguous={m_ambiguous} "
            f"groups={m_groups}"
        )

        print(
            f"        depth: "
            f"ambiguous={d_ambiguous} "
            f"groups={d_groups}"
        )

        if m_ambiguous == 0:

            exact_m.append(
                quad
            )

        if d_ambiguous == 0:

            exact_depth.append(
                quad
            )

        if (
            m_ambiguous
            and m_ce
        ):

            (
                index,
                signature,
                old_value,
                new_value,
            ) = m_ce

            print(
                f"        first collision: "
                f"index={index} "
                f"signature={signature} "
                f"targets="
                f"{old_value},{new_value}"
            )

    return (
        exact_m,
        exact_depth,
    )


# ==============================================================================
# EQUALITY SIGNATURE SEARCH
# ==============================================================================

def search_equality_signatures(
    columns,
    offsets,
    equality_values,
    sample_indices,
):

    print()
    print(
        "=" * 90
    )
    print(
        "TEST 7: EQUALITY-BIT TRIPLE SEARCH"
    )
    print(
        "=" * 90
    )

    ranked = []

    for triple in combinations(
        offsets,
        3,
    ):

        ambiguous, groups, _ = (
            rank_signature(
                triple,
                columns,
                equality_values,
                sample_indices,
            )
        )

        ranked.append(
            (
                ambiguous,
                -groups,
                triple,
            )
        )

    ranked.sort()

    for ambiguous, neg_groups, triple in ranked[:15]:

        print(
            f"    {triple} "
            f"ambiguous={ambiguous} "
            f"groups={-neg_groups}"
        )


# ==============================================================================
# STRUCTURAL SIGNATURES
# ==============================================================================

def structural_signatures(
    columns,
    states,
    m_values,
    depth_values,
):

    print()
    print(
        "=" * 90
    )
    print(
        "TEST 8: STRUCTURAL n-ONLY SIGNATURES"
    )
    print(
        "=" * 90
    )

    signatures = [
        (-9, -3, 3),
        (-9, -3, 9),
        (-9, 3, 9),
        (-9, -3, 61),
        (-3, 3, 61),
        (-9, 5, 61),
    ]

    targets = {
        "m": m_values,
        "depth": depth_values,
    }

    for sig in signatures:

        print()
        print(
            f"    signature={sig}"
        )

        for name, target in targets.items():

            ambiguous, groups, _ = (
                exact_signature(
                    sig,
                    columns,
                    target,
                )
            )

            print(
                f"        {name}: "
                f"ambiguous={ambiguous} "
                f"groups={groups}"
            )


# ==============================================================================
# MODULAR SIGNATURE TEST
# ==============================================================================

def modular_signature_search(
    states,
    columns,
    m_values,
    depth_values,
):

    print()
    print(
        "=" * 90
    )
    print(
        "TEST 9: n MOD 2^k + THREE VALUATIONS"
    )
    print(
        "=" * 90
    )

    candidates = [
        (-9, -3, 3),
        (-9, -3, 9),
        (-9, 3, 9),
        (-9, -3, 61),
        (-3, 3, 61),
    ]

    ns = [
        state.n
        for state in states
    ]

    for bits in (
        6,
        8,
        10,
        12,
    ):

        modulus = 1 << bits

        for sig in candidates:

            mapping = {}

            amb_m = set()
            amb_d = set()

            cols = [
                columns[c]
                for c in sig
            ]

            for i, n in enumerate(ns):

                key = (
                    n % modulus,
                    cols[0][i],
                    cols[1][i],
                    cols[2][i],
                )

                value = (
                    m_values[i],
                    depth_values[i],
                )

                previous = mapping.get(
                    key
                )

                if previous is None:

                    mapping[key] = value

                elif previous != value:

                    if previous[0] != value[0]:
                        amb_m.add(key)

                    if previous[1] != value[1]:
                        amb_d.add(key)

            print(
                f"    bits={bits:2} "
                f"sig={sig} "
                f"m_amb={len(amb_m)} "
                f"depth_amb={len(amb_d)}"
            )


# ==============================================================================
# EXAMPLES
# ==============================================================================

def examples(
    states,
    columns,
    m_values,
    depth_values,
):

    print()
    print(
        "=" * 90
    )
    print(
        "TEST 10: EXAMPLES"
    )
    print(
        "=" * 90
    )

    wanted = {
        9,
        15,
        21,
        33,
        39,
        57,
        69,
        93,
        141,
        183,
    }

    offsets = (
        -9,
        -3,
        3,
        9,
        61,
    )

    for i, state in enumerate(states):

        if state.n not in wanted:
            continue

        print()
        print(
            f"n={state.n}"
        )

        print(
            f"    p={state.p} "
            f"q={state.q}"
        )

        print(
            f"    m={m_values[i]} "
            f"depth={depth_values[i]}"
        )

        for c in offsets:

            print(
                f"    v2(n{c:+d})="
                f"{columns[c][i]}"
            )


# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

def final_summary(
    exact_m_triples,
    exact_depth_triples,
    exact_m_quads,
    exact_depth_quads,
):

    print()
    print(
        "=" * 90
    )
    print(
        "FINAL SUMMARY"
    )
    print(
        "=" * 90
    )

    print(
        f"""
The current exact law is:

    depth =
        min(v2(A),v2(B))
        + [v2(A)=v2(B)].

Experiment 623 showed that no tested PAIR of
v2(n+c) values was sufficient.

Experiment 624 extends this to:

    3-value signatures

and then:

    4-value signatures around the strongest
    3-value candidates.

Exact m triples:
    {len(exact_m_triples)}

Exact depth triples:
    {len(exact_depth_triples)}

Exact m quadruples:
    {len(exact_m_quads)}

Exact depth quadruples:
    {len(exact_depth_quads)}

Interpretation:

    exact m triple
        =>
    m is an n-only function of three valuations.

    exact depth triple
        =>
    the entire depth is an n-only function
    of three valuations.

    exact quadruple
        =>
    a slightly larger but still very small
    n-only signature is sufficient.

If all signatures fail, the evidence becomes stronger
that m is not a low-dimensional valuation signature
of the tested n+c family.
"""
    )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print(
        "EXPERIMENT 624 START"
    )
    print("=" * 90)
    print()
    print(
        "SMALL n-ONLY 2-ADIC SIGNATURE SEARCH"
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
        "[3] BUILD TARGETS"
    )

    (
        m_values,
        depth_values,
        equality_values,
    ) = build_targets(
        states
    )

    print(
        f"    states={len(states)}"
    )

    # --------------------------------------------------------------------------
    # 4
    # --------------------------------------------------------------------------

    print()
    print(
        "[4] BUILD n-ONLY VALUATIONS"
    )

    columns = build_columns(
        states,
        OFFSET_POOL,
    )

    print(
        f"    offsets="
        f"{len(OFFSET_POOL)}"
    )

    # --------------------------------------------------------------------------
    # 5
    # --------------------------------------------------------------------------

    sample_indices = make_sample(
        len(states),
        SAMPLE_SIZE,
    )

    print()
    print(
        "[5] SAMPLE"
    )

    print(
        f"    sample size="
        f"{len(sample_indices)}"
    )

    # --------------------------------------------------------------------------
    # 6
    #
    # Triple search for m.
    # --------------------------------------------------------------------------

    triples_m = rank_triples(
        columns,
        OFFSET_POOL,
        m_values,
        sample_indices,
    )

    # --------------------------------------------------------------------------
    # 7
    #
    # Exact triples.
    # --------------------------------------------------------------------------

    (
        exact_m_triples,
        exact_depth_triples,
    ) = exact_triples(
        triples_m,
        columns,
        m_values,
        depth_values,
    )

    # --------------------------------------------------------------------------
    # 8
    #
    # Quadrupes around strong m triples.
    # --------------------------------------------------------------------------

    quadruples = (
        rank_quads_from_triples(
            triples_m,
            OFFSET_POOL,
            columns,
            m_values,
            sample_indices,
        )
    )

    (
        exact_m_quads,
        exact_depth_quads,
    ) = exact_quads(
        quadruples,
        columns,
        m_values,
        depth_values,
    )

    # --------------------------------------------------------------------------
    # 9
    #
    # Search equality separately.
    # --------------------------------------------------------------------------

    search_equality_signatures(
        columns,
        OFFSET_POOL,
        equality_values,
        sample_indices,
    )

    # --------------------------------------------------------------------------
    # 10
    #
    # Structural candidates.
    # --------------------------------------------------------------------------

    structural_signatures(
        columns,
        states,
        m_values,
        depth_values,
    )

    # --------------------------------------------------------------------------
    # 11
    #
    # Modular signature.
    # --------------------------------------------------------------------------

    modular_signature_search(
        states,
        columns,
        m_values,
        depth_values,
    )

    # --------------------------------------------------------------------------
    # 12
    # --------------------------------------------------------------------------

    examples(
        states,
        columns,
        m_values,
        depth_values,
    )

    # --------------------------------------------------------------------------
    # 13
    # --------------------------------------------------------------------------

    final_summary(
        exact_m_triples,
        exact_depth_triples,
        exact_m_quads,
        exact_depth_quads,
    )

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 624 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()
