#!/usr/bin/env python3

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from random import Random
from math import isqrt


# ==============================================================================
# EXPERIMENT 625
# ==============================================================================
#
# FAST n-ONLY 2-ADIC SIGNATURE SEARCH
#
# EXPERIMENT 624 was too slow because it directly ranked all triples:
#
#     C(35,3) = 6545 triples
#
# against 30,000 sample states.
#
# That creates roughly 196 million Python-level signature evaluations.
#
# This experiment uses a staged search:
#
#     offsets
#        |
#        v
#     PAIR SEARCH
#        |
#        v
#     best pairs only
#        |
#        v
#     TRIPLE EXPANSION
#        |
#        v
#     best triples only
#        |
#        v
#     EXACT FULL-DOMAIN VERIFICATION
#
# The aim is discovery speed, not exhaustive combinatorial coverage.
#
# ==============================================================================


# ==============================================================================
# CONFIGURATION
# ==============================================================================

PRIME_LIMIT = 6250

SAMPLE_SIZE = 5000

TOP_PAIRS = 40
TOP_TRIPLES = 20

RANDOM_SEED = 625


# Smaller, structurally motivated offset pool.
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

def frame_from_n(
    n: int,
) -> str:

    r = n & 3

    if r == 3:
        return "A"

    if r == 1:
        return "B"

    raise ValueError(
        f"Invalid odd semiprime n={n}"
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
# TARGETS
# ==============================================================================

def build_targets(
    states: list[State],
) -> tuple[
    list[int],
    list[int],
    list[int],
]:

    m_values: list[int] = []
    depth_values: list[int] = []
    equality_values: list[int] = []

    for state in states:

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

        equality = int(
            va == vb
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

    return {
        c: [
            v2(n + c)
            for n in ns
        ]
        for c in offsets
    }


# ==============================================================================
# SAMPLE
# ==============================================================================

def deterministic_sample(
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
# FAST PAIR SIGNATURE
# ==============================================================================

def rank_pair(
    c1: int,
    c2: int,
    columns: dict[int, list[int]],
    target: list[int],
    indices: list[int],
) -> tuple[int, int]:

    col1 = columns[c1]
    col2 = columns[c2]

    mapping: dict[
        tuple[int, int],
        int,
    ] = {}

    ambiguous: set[
        tuple[int, int]
    ] = set()

    for i in indices:

        key = (
            col1[i],
            col2[i],
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
    )


# ==============================================================================
# EXACT PAIR
# ==============================================================================

def exact_pair(
    pair: tuple[int, int],
    columns: dict[int, list[int]],
    target: list[int],
) -> tuple[
    int,
    int,
    tuple | None,
]:

    c1, c2 = pair

    col1 = columns[c1]
    col2 = columns[c2]

    mapping: dict[
        tuple[int, int],
        int,
    ] = {}

    ambiguous = set()
    counterexample = None

    for i, value in enumerate(target):

        key = (
            col1[i],
            col2[i],
        )

        previous = mapping.get(
            key
        )

        if previous is None:

            mapping[key] = value

        elif previous != value:

            ambiguous.add(
                key
            )

            if counterexample is None:

                counterexample = (
                    i,
                    key,
                    previous,
                    value,
                )

    return (
        len(ambiguous),
        len(mapping),
        counterexample,
    )


# ==============================================================================
# PAIR SEARCH
# ==============================================================================

def search_pairs(
    columns,
    offsets,
    m_values,
    depth_values,
    equality_values,
    sample_indices,
):

    print()
    print("=" * 90)
    print(
        "TEST 1: FAST PAIR SEARCH"
    )
    print("=" * 90)

    candidates_m = []
    candidates_depth = []
    candidates_eq = []

    for pair in combinations(
        offsets,
        2,
    ):

        amb_m, groups_m = rank_pair(
            pair[0],
            pair[1],
            columns,
            m_values,
            sample_indices,
        )

        amb_d, groups_d = rank_pair(
            pair[0],
            pair[1],
            columns,
            depth_values,
            sample_indices,
        )

        amb_e, groups_e = rank_pair(
            pair[0],
            pair[1],
            columns,
            equality_values,
            sample_indices,
        )

        candidates_m.append(
            (
                amb_m,
                -groups_m,
                pair,
            )
        )

        candidates_depth.append(
            (
                amb_d,
                -groups_d,
                pair,
            )
        )

        candidates_eq.append(
            (
                amb_e,
                -groups_e,
                pair,
            )
        )

    candidates_m.sort()
    candidates_depth.sort()
    candidates_eq.sort()

    print(
        "    BEST PAIRS FOR m:"
    )

    for amb, neg_groups, pair in (
        candidates_m[:20]
    ):

        print(
            f"        {pair} "
            f"ambiguous={amb} "
            f"groups={-neg_groups}"
        )

    print()
    print(
        "    BEST PAIRS FOR depth:"
    )

    for amb, neg_groups, pair in (
        candidates_depth[:20]
    ):

        print(
            f"        {pair} "
            f"ambiguous={amb} "
            f"groups={-neg_groups}"
        )

    print()
    print(
        "    BEST PAIRS FOR equality:"
    )

    for amb, neg_groups, pair in (
        candidates_eq[:20]
    ):

        print(
            f"        {pair} "
            f"ambiguous={amb} "
            f"groups={-neg_groups}"
        )

    return (
        [
            x[2]
            for x in candidates_m[
                :TOP_PAIRS
            ]
        ],
        [
            x[2]
            for x in candidates_depth[
                :TOP_PAIRS
            ]
        ],
        [
            x[2]
            for x in candidates_eq[
                :TOP_PAIRS
            ]
        ],
    )


# ==============================================================================
# EXPAND PAIRS INTO TRIPLES
# ==============================================================================

def expand_pairs(
    pairs,
    offsets,
) -> list[tuple[int, int, int]]:

    triples = set()

    for pair in pairs:

        for c in offsets:

            if c in pair:
                continue

            triple = tuple(
                sorted(
                    (
                        pair[0],
                        pair[1],
                        c,
                    )
                )
            )

            triples.add(
                triple
            )

    return sorted(
        triples
    )


# ==============================================================================
# FAST TRIPLE RANKING
# ==============================================================================

def rank_triple(
    triple,
    columns,
    target,
    indices,
) -> tuple[int, int]:

    c1, c2, c3 = triple

    col1 = columns[c1]
    col2 = columns[c2]
    col3 = columns[c3]

    mapping = {}
    ambiguous = set()

    for i in indices:

        key = (
            col1[i],
            col2[i],
            col3[i],
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
    )


# ==============================================================================
# TRIPLE SEARCH
# ==============================================================================

def search_triples(
    candidate_pairs,
    offsets,
    columns,
    target,
    sample_indices,
):

    print()
    print("=" * 90)
    print(
        "TEST 2: TARGETED TRIPLE SEARCH"
    )
    print("=" * 90)

    triples = expand_pairs(
        candidate_pairs,
        offsets,
    )

    print(
        f"    triples generated="
        f"{len(triples)}"
    )

    ranked = []

    for triple in triples:

        amb, groups = rank_triple(
            triple,
            columns,
            target,
            sample_indices,
        )

        ranked.append(
            (
                amb,
                -groups,
                triple,
            )
        )

    ranked.sort()

    print(
        "    best triples:"
    )

    for amb, neg_groups, triple in (
        ranked[:30]
    ):

        print(
            f"        {triple} "
            f"ambiguous={amb} "
            f"groups={-neg_groups}"
        )

    return [
        triple
        for _, _, triple
        in ranked[:TOP_TRIPLES]
    ]


# ==============================================================================
# EXACT TRIPLE VERIFICATION
# ==============================================================================

def verify_triples(
    triples,
    columns,
    m_values,
    depth_values,
    equality_values,
):

    print()
    print("=" * 90)
    print(
        "TEST 3: EXACT TRIPLE VERIFICATION"
    )
    print("=" * 90)

    exact_m = []
    exact_depth = []
    exact_eq = []

    for triple in triples:

        m_amb, m_groups, m_ce = exact_signature(
            triple,
            columns,
            m_values,
        )

        d_amb, d_groups, d_ce = exact_signature(
            triple,
            columns,
            depth_values,
        )

        e_amb, e_groups, e_ce = exact_signature(
            triple,
            columns,
            equality_values,
        )

        print()
        print(
            f"    triple={triple}"
        )

        print(
            f"        m: "
            f"ambiguous={m_amb} "
            f"groups={m_groups}"
        )

        print(
            f"        depth: "
            f"ambiguous={d_amb} "
            f"groups={d_groups}"
        )

        print(
            f"        equality: "
            f"ambiguous={e_amb} "
            f"groups={e_groups}"
        )

        if m_amb == 0:
            exact_m.append(
                triple
            )

        if d_amb == 0:
            exact_depth.append(
                triple
            )

        if e_amb == 0:
            exact_eq.append(
                triple
            )

        if m_ce is not None:

            print(
                f"        first m collision="
                f"{m_ce}"
            )

    return (
        exact_m,
        exact_depth,
        exact_eq,
    )


# ==============================================================================
# EXACT SIGNATURE
# ==============================================================================

def exact_signature(
    signature,
    columns,
    target,
):

    cols = [
        columns[c]
        for c in signature
    ]

    mapping = {}
    ambiguous = set()
    counterexample = None

    for i, value in enumerate(target):

        key = tuple(
            col[i]
            for col in cols
        )

        previous = mapping.get(
            key
        )

        if previous is None:

            mapping[key] = value

        elif previous != value:

            ambiguous.add(
                key
            )

            if counterexample is None:

                counterexample = (
                    i,
                    key,
                    previous,
                    value,
                )

    return (
        len(ambiguous),
        len(mapping),
        counterexample,
    )


# ==============================================================================
# STRUCTURAL FAMILY SEARCH
# ==============================================================================

def structural_family_search(
    columns,
    m_values,
    depth_values,
):

    print()
    print("=" * 90)
    print(
        "TEST 4: THEORETICAL SIGNATURE FAMILIES"
    )
    print("=" * 90)

    families = [
        (-9, -3),
        (-9, -3, 3),
        (-9, -3, 9),
        (-9, 3, 9),
        (-9, -3, 61),
        (-9, -3, 5),
        (-9, -3, 13),
        (-3, 3, 9),
        (-3, 9, 61),
    ]

    for signature in families:

        m_amb, m_groups, _ = exact_signature(
            signature,
            columns,
            m_values,
        )

        d_amb, d_groups, _ = exact_signature(
            signature,
            columns,
            depth_values,
        )

        print(
            f"    {signature}"
        )

        print(
            f"        m: "
            f"ambiguous={m_amb} "
            f"groups={m_groups}"
        )

        print(
            f"        depth: "
            f"ambiguous={d_amb} "
            f"groups={d_groups}"
        )


# ==============================================================================
# DIRECT FRAME-SPLIT SEARCH
# ==============================================================================

def frame_split_search(
    states,
    columns,
    m_values,
    depth_values,
):

    print()
    print("=" * 90)
    print(
        "TEST 5: FRAME-SPLIT SIGNATURE SEARCH"
    )
    print("=" * 90)

    # Because frame is already known from n mod 4,
    # search signatures independently inside A/B.
    #
    # This can reveal a signature that is deterministic only
    # after the frame is supplied.

    indices_A = [
        i
        for i, s in enumerate(states)
        if (s.n & 3) == 3
    ]

    indices_B = [
        i
        for i, s in enumerate(states)
        if (s.n & 3) == 1
    ]

    candidates = [
        (-9,),
        (-3,),
        (-9, -3),
        (-9, 3),
        (-3, 3),
        (-9, -3, 3),
        (-9, -3, 9),
        (-9, 3, 9),
        (-9, -3, 61),
        (-3, 3, 61),
    ]

    for frame_name, frame_indices in (
        ("A", indices_A),
        ("B", indices_B),
    ):

        print()
        print(
            f"    FRAME {frame_name}"
        )

        for signature in candidates:

            m_amb, _, _ = exact_signature(
                signature,
                columns,
                restricted_target(
                    m_values,
                    frame_indices,
                ),
                indices=frame_indices,
            )

            d_amb, _, _ = exact_signature(
                signature,
                columns,
                restricted_target(
                    depth_values,
                    frame_indices,
                ),
                indices=frame_indices,
            )

            print(
                f"        {signature}: "
                f"m_amb={m_amb} "
                f"depth_amb={d_amb}"
            )


# ==============================================================================
# Restricted exact signature helper
# ==============================================================================

def restricted_target(
    target,
    indices,
):
    return [
        target[i]
        for i in indices
    ]


# ==============================================================================
# FRAME SPLIT EXACT SIGNATURE
# ==============================================================================

def exact_signature_indices(
    signature,
    columns,
    target,
    indices,
):

    cols = [
        columns[c]
        for c in signature
    ]

    mapping = {}
    ambiguous = set()

    for local_i, real_i in enumerate(indices):

        key = tuple(
            col[real_i]
            for col in cols
        )

        value = target[local_i]

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
    )


# Replace the earlier convenience helper with a direct implementation.
def frame_signature_test(
    signature,
    columns,
    target,
    indices,
):

    cols = [
        columns[c]
        for c in signature
    ]

    mapping = {}
    ambiguous = set()

    for real_i in indices:

        key = tuple(
            col[real_i]
            for col in cols
        )

        value = target[real_i]

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
    )


# ==============================================================================
# FIXED FRAME SEARCH
# ==============================================================================

def fixed_frame_search(
    states,
    columns,
    m_values,
    depth_values,
):

    print()
    print("=" * 90)
    print(
        "TEST 5: FRAME-SUPPLIED n-ONLY SIGNATURES"
    )
    print("=" * 90)

    indices_A = [
        i
        for i, s in enumerate(states)
        if (s.n & 3) == 3
    ]

    indices_B = [
        i
        for i, s in enumerate(states)
        if (s.n & 3) == 1
    ]

    signatures = [
        (-9,),
        (-3,),
        (-9, -3),
        (-9, 3),
        (-3, 3),
        (-9, -3, 3),
        (-9, -3, 9),
        (-9, 3, 9),
        (-9, -3, 61),
        (-3, 3, 61),
    ]

    for frame_name, indices in (
        ("A", indices_A),
        ("B", indices_B),
    ):

        print()
        print(
            f"    FRAME {frame_name}"
        )

        for signature in signatures:

            m_amb, m_groups = (
                frame_signature_test(
                    signature,
                    columns,
                    m_values,
                    indices,
                )
            )

            d_amb, d_groups = (
                frame_signature_test(
                    signature,
                    columns,
                    depth_values,
                    indices,
                )
            )

            print(
                f"        {signature}: "
                f"m_amb={m_amb} "
                f"depth_amb={d_amb}"
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
    print("=" * 90)
    print(
        "TEST 6: EXAMPLES"
    )
    print("=" * 90)

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

        frame = frame_from_n(
            state.n
        )

        print()
        print(
            f"n={state.n}"
        )

        print(
            f"    p={state.p} "
            f"q={state.q} "
            f"frame={frame}"
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
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print(
        "EXPERIMENT 625 START"
    )
    print("=" * 90)
    print()
    print(
        "FAST STAGED n-ONLY 2-ADIC SIGNATURE SEARCH"
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

    # --------------------------------------------------------------------------
    # 4
    # --------------------------------------------------------------------------

    print()
    print(
        "[4] PRECOMPUTE n-ONLY VALUATIONS"
    )

    columns = build_columns(
        states,
        OFFSET_POOL,
    )

    print(
        f"    offsets={len(OFFSET_POOL)}"
    )

    # --------------------------------------------------------------------------
    # 5
    # --------------------------------------------------------------------------

    sample_indices = deterministic_sample(
        len(states),
        SAMPLE_SIZE,
    )

    print()
    print(
        "[5] SAMPLE"
    )

    print(
        f"    sample={len(sample_indices)}"
    )

    # --------------------------------------------------------------------------
    # 6
    #
    # Pair search for m.
    # --------------------------------------------------------------------------

    (
        pairs_m,
        pairs_depth,
        pairs_eq,
    ) = search_pairs(
        columns,
        OFFSET_POOL,
        m_values,
        depth_values,
        equality_values,
        sample_indices,
    )

    # --------------------------------------------------------------------------
    # 7
    #
    # Triple expansion from best pairs for m.
    # --------------------------------------------------------------------------

    triples_m = search_triples(
        pairs_m,
        OFFSET_POOL,
        columns,
        m_values,
        sample_indices,
    )

    # --------------------------------------------------------------------------
    # 8
    #
    # Exact triples.
    # --------------------------------------------------------------------------

    exact_m, exact_depth, exact_eq = (
        verify_triples(
            triples_m,
            columns,
            m_values,
            depth_values,
            equality_values,
        )
    )

    # --------------------------------------------------------------------------
    # 9
    #
    # Structural candidates.
    # --------------------------------------------------------------------------

    structural_family_search(
        columns,
        m_values,
        depth_values,
    )

    # --------------------------------------------------------------------------
    # 10
    #
    # Frame supplied.
    # --------------------------------------------------------------------------

    fixed_frame_search(
        states,
        columns,
        m_values,
        depth_values,
    )

    # --------------------------------------------------------------------------
    # 11
    # --------------------------------------------------------------------------

    examples(
        states,
        columns,
        m_values,
        depth_values,
    )

    # --------------------------------------------------------------------------
    # 12
    # --------------------------------------------------------------------------

    print()
    print("=" * 90)
    print(
        "FINAL SUMMARY"
    )
    print("=" * 90)

    print(
        f"""
SEARCH RESULTS

Exact m triples:
    {len(exact_m)}

Exact depth triples:
    {len(exact_depth)}

Exact equality triples:
    {len(exact_eq)}

The search is deliberately staged:

    35 offsets
       ->
    pair ranking
       ->
    40 best pairs
       ->
    targeted triples
       ->
    20 exact candidates.

An exact triple for m would establish:

    m =
        f(
            v2(n+c1),
            v2(n+c2),
            v2(n+c3)
        ).

An exact triple for depth would establish:

    depth =
        f(
            v2(n+c1),
            v2(n+c2),
            v2(n+c3)
        ).

Failure of these candidates does not prove that no
n-only signature exists; it only rejects this
small structured family.
"""
    )

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 625 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()
