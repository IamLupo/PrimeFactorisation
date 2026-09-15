#!/usr/bin/env python3

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from math import isqrt
from random import Random


# ==============================================================================
# EXPERIMENT 623
# ==============================================================================
#
# FAST n-ONLY 2-ADIC SIGNATURE SEARCH
#
# Experiment 622 attempted all pairs over a large offset set and became too
# expensive.
#
# This experiment uses:
#
#     stage 1:
#         rank candidate pairs on a deterministic sample
#
#     stage 2:
#         EXACTLY verify only the best candidates against all semiprimes.
#
# Target:
#
#     m = min(v2(A), v2(B))
#
# and:
#
#     depth = m + [v2(A)=v2(B)].
#
# ==============================================================================


# ==============================================================================
# CONFIGURATION
# ==============================================================================

PRIME_LIMIT = 6250

# Sample size used for fast candidate ranking.
SAMPLE_SIZE = 25000

# Number of candidates sent to exact verification.
TOP_EXACT = 30

# Candidate offsets.
#
# The first four are structurally important:
#
#     FRAME A: n+9
#     FRAME B: n+3
#
# We also include the best offsets from Experiment 622.
#
# Additional nearby offsets are included for discovery.
OFFSET_POOL = sorted(
    {
        -59,
        -51,
        -43,
        -35,
        -27,
        -19,
        -11,
        -9,
        -3,
        3,
        5,
        9,
        13,
        21,
        29,
        37,
        41,
        45,
        53,
        61,
    }
)

MAX_EXAMPLES = 10

RANDOM_SEED = 622


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
        f"Unexpected odd-semiprime residue: "
        f"n={n} mod 4={r}"
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
# TARGET CACHE
# ==============================================================================

def build_targets(
    states: list[State],
) -> tuple[
    list[int],
    list[int],
]:

    m_values: list[int] = []
    depth_values: list[int] = []

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

        depth = (
            m
            + int(alpha == beta)
        )

        m_values.append(m)
        depth_values.append(depth)

    return (
        m_values,
        depth_values,
    )


# ==============================================================================
# PRECOMPUTE VALUATIONS
# ==============================================================================

def build_valuation_columns(
    states: list[State],
    offsets: list[int],
) -> dict[int, list[int]]:

    print(
        "    precomputing valuation columns..."
    )

    columns: dict[int, list[int]] = {}

    ns = [
        state.n
        for state in states
    ]

    for c in offsets:

        columns[c] = [
            v2(n + c)
            for n in ns
        ]

    return columns


# ==============================================================================
# RANDOM SAMPLE
# ==============================================================================

def make_sample(
    size: int,
    total: int,
) -> list[int]:

    if size >= total:

        return list(
            range(total)
        )

    rng = Random(
        RANDOM_SEED
    )

    indices = rng.sample(
        range(total),
        size,
    )

    indices.sort()

    return indices


# ==============================================================================
# PAIR SCORE
# ==============================================================================

def pair_score_sample(
    pair: tuple[int, int],
    columns: dict[int, list[int]],
    target: list[int],
    indices: list[int],
) -> tuple[int, int]:

    c1, c2 = pair

    col1 = columns[c1]
    col2 = columns[c2]

    groups: dict[
        tuple[int, int],
        int,
    ] = {}

    ambiguous = set()

    for i in indices:

        key = (
            col1[i],
            col2[i],
        )

        value = target[i]

        old = groups.get(
            key
        )

        if old is None:

            groups[key] = value

        elif old != value:

            ambiguous.add(
                key
            )

    return (
        len(ambiguous),
        len(groups),
    )


# ==============================================================================
# EXACT PAIR TEST
# ==============================================================================

def pair_exact_test(
    pair: tuple[int, int],
    columns: dict[int, list[int]],
    target: list[int],
) -> tuple[int, int, tuple | None]:

    c1, c2 = pair

    col1 = columns[c1]
    col2 = columns[c2]

    groups: dict[
        tuple[int, int],
        int,
    ] = {}

    ambiguous = set()

    first_counterexample = None

    for i in range(
        len(target)
    ):

        key = (
            col1[i],
            col2[i],
        )

        value = target[i]

        old = groups.get(
            key
        )

        if old is None:

            groups[key] = value

        elif old != value:

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
                        old,
                        value,
                    )

    return (
        len(ambiguous),
        len(groups),
        first_counterexample,
    )


# ==============================================================================
# TOP PAIRS
# ==============================================================================

def rank_pairs(
    columns: dict[int, list[int]],
    offsets: list[int],
    target: list[int],
    sample_indices: list[int],
    label: str,
) -> list[
    tuple[
        tuple[int, int],
        int,
        int,
    ]
]:

    print()
    print(
        f"    ranking pairs for {label}..."
    )

    ranked = []

    for pair in combinations(
        offsets,
        2,
    ):

        ambiguous, groups = (
            pair_score_sample(
                pair,
                columns,
                target,
                sample_indices,
            )
        )

        ranked.append(
            (
                ambiguous,
                pair,
                groups,
            )
        )

    ranked.sort(
        key=lambda x: (
            x[0],
            -x[2],
        )
    )

    print(
        f"    sample pairs tested="
        f"{len(ranked)}"
    )

    print(
        "    best sample candidates:"
    )

    for ambiguous, pair, groups in ranked[:20]:

        print(
            f"        pair={pair} "
            f"ambiguous={ambiguous} "
            f"groups={groups}"
        )

    return [
        (
            pair,
            ambiguous,
            groups,
        )
        for ambiguous, pair, groups
        in ranked[
            :TOP_EXACT
        ]
    ]


# ==============================================================================
# EXACT VERIFICATION
# ==============================================================================

def verify_top_pairs(
    candidates,
    columns,
    m_values,
    depth_values,
) -> tuple[
    list,
    list,
]:

    print()
    print(
        "    EXACT VERIFICATION"
    )

    exact_m = []
    exact_depth = []

    for pair, sample_ambiguous, _ in candidates:

        amb_m, groups_m, ce_m = pair_exact_test(
            pair,
            columns,
            m_values,
        )

        amb_d, groups_d, ce_d = pair_exact_test(
            pair,
            columns,
            depth_values,
        )

        print(
            f"        pair={pair}"
        )

        print(
            f"            m: "
            f"ambiguous={amb_m} "
            f"groups={groups_m}"
        )

        print(
            f"            depth: "
            f"ambiguous={amb_d} "
            f"groups={groups_d}"
        )

        if amb_m == 0:

            exact_m.append(
                pair
            )

        if amb_d == 0:

            exact_depth.append(
                pair
            )

        if (
            amb_m
            and ce_m
        ):

            i, key, old, value = ce_m

            print(
                f"            first m "
                f"counterexample:"
            )

            print(
                f"                index={i}"
            )

            print(
                f"                signature={key}"
            )

            print(
                f"                old_target={old}"
            )

            print(
                f"                new_target={value}"
            )

    return (
        exact_m,
        exact_depth,
    )


# ==============================================================================
# FOCUSED THEORETICAL PAIRS
# ==============================================================================

def test_theoretical_pairs(
    columns,
    m_values,
    depth_values,
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 4: THEORETICAL OFFSET PAIRS"
    )
    print("=" * 90)

    pairs = [
        (-9, -3),
        (-9, 3),
        (-9, 9),
        (-3, 3),
        (-3, 9),
        (3, 9),
    ]

    for pair in pairs:

        amb_m, groups_m, _ = pair_exact_test(
            pair,
            columns,
            m_values,
        )

        amb_d, groups_d, _ = pair_exact_test(
            pair,
            columns,
            depth_values,
        )

        print(
            f"    pair={pair}"
        )

        print(
            f"        m: "
            f"ambiguous={amb_m} "
            f"groups={groups_m}"
        )

        print(
            f"        depth: "
            f"ambiguous={amb_d} "
            f"groups={groups_d}"
        )


# ==============================================================================
# SINGLE OFFSET BASELINE
# ==============================================================================

def single_offset_baseline(
    columns,
    offsets,
    m_values,
    depth_values,
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 2: SINGLE-OFFSET BASELINE"
    )
    print("=" * 90)

    results_m = []
    results_depth = []

    for c in offsets:

        amb_m = 0
        amb_depth = 0

        groups_m = {}
        groups_depth = {}

        values = columns[c]

        for i, value in enumerate(values):

            old = groups_m.get(
                value
            )

            if old is None:
                groups_m[value] = m_values[i]

            elif old != m_values[i]:
                amb_m += 1

            old = groups_depth.get(
                value
            )

            if old is None:
                groups_depth[value] = depth_values[i]

            elif old != depth_values[i]:
                amb_depth += 1

        results_m.append(
            (
                amb_m,
                c,
            )
        )

        results_depth.append(
            (
                amb_depth,
                c,
            )
        )

    results_m.sort()
    results_depth.sort()

    print(
        "    BEST SINGLE OFFSETS FOR m:"
    )

    for amb, c in results_m[:10]:

        print(
            f"        c={c:>4} "
            f"ambiguous={amb}"
        )

    print()

    print(
        "    BEST SINGLE OFFSETS FOR depth:"
    )

    for amb, c in results_depth[:10]:

        print(
            f"        c={c:>4} "
            f"ambiguous={amb}"
        )


# ==============================================================================
# n MODULO + VALUATION
# ==============================================================================

def test_residue_signature(
    states,
    columns,
    m_values,
    depth_values,
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 5: n MOD 2^k + ONE VALUATION"
    )
    print("=" * 90)

    configs = [
        (-3, 8),
        (-3, 10),
        (-3, 12),
        (-9, 8),
        (-9, 10),
        (-9, 12),
        (3, 8),
        (3, 10),
        (3, 12),
    ]

    ns = [
        state.n
        for state in states
    ]

    for c, bits in configs:

        modulus = 1 << bits

        mapping = {}
        ambiguous_m = set()
        ambiguous_depth = set()

        values = columns[c]

        for i, n in enumerate(ns):

            key = (
                n % modulus,
                values[i],
            )

            target = (
                m_values[i],
                depth_values[i],
            )

            old = mapping.get(
                key
            )

            if old is None:

                mapping[key] = target

            elif old != target:

                if old[0] != target[0]:
                    ambiguous_m.add(
                        key
                    )

                if old[1] != target[1]:
                    ambiguous_depth.add(
                        key
                    )

        print(
            f"    c={c:>3} "
            f"bits={bits:>2} "
            f"m_ambiguous="
            f"{len(ambiguous_m):<8} "
            f"depth_ambiguous="
            f"{len(ambiguous_depth)}"
        )


# ==============================================================================
# GLOBAL DEPTH CROSSCHECK
# ==============================================================================

def test_global_crosscheck(
    states,
    m_values,
    depth_values,
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 6: n-ONLY TARGET CROSSCHECK"
    )
    print("=" * 90)

    failures = 0

    for i, state in enumerate(states):

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

        expected_m = min(
            v2(A),
            v2(B),
        )

        expected_depth = (
            1
            + min(
                v2(X),
                v2(Y),
            )
        )

        if (
            expected_m
            != m_values[i]
            or
            expected_depth
            != depth_values[i]
        ):

            failures += 1

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )


# ==============================================================================
# GLOBAL XY
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


# ==============================================================================
# EXAMPLES
# ==============================================================================

def examples(
    states,
    columns,
    m_values,
    depth_values,
) -> None:

    print()
    print("=" * 90)
    print(
        "EXAMPLES"
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

    for i, state in enumerate(states):

        if state.n not in wanted:
            continue

        frame = frame_from_n(
            state.n
        )

        A, B = raw_residuals(
            frame,
            state.p,
            state.q,
        )

        print()
        print(
            f"n={state.n} "
            f"p={state.p} "
            f"q={state.q}"
        )

        print(
            f"    frame={frame}"
        )

        print(
            f"    A={A} "
            f"B={B}"
        )

        print(
            f"    m={m_values[i]} "
            f"depth={depth_values[i]}"
        )

        for c in (
            -9,
            -3,
            3,
            9,
        ):

            print(
                f"    v2(n{c:+d})="
                f"{columns[c][i]}"
            )


# ==============================================================================
# SUMMARY
# ==============================================================================

def summary(
    exact_m,
    exact_depth,
) -> None:

    print()
    print("=" * 90)
    print(
        "FINAL SUMMARY"
    )
    print("=" * 90)

    print(
r"""
The exact factor-side law remains:

    m = min(v2(A),v2(B))

    depth =
        m + [v2(A)=v2(B)].

The open question is whether m is recoverable from a small
n-only 2-adic signature.

This experiment does NOT spend time exhaustively checking every
possible pair over the full domain.

Instead:

    1. precompute v2(n+c)
    2. rank candidate pairs on a deterministic sample
    3. exactly verify only the strongest candidates.

An exact pair in:

    EXACT m PAIRS

would establish:

    (v2(n+c1), v2(n+c2))
        ->
    m

and therefore:

    depth =
        m + equality_bit.

An exact pair in:

    EXACT DEPTH PAIRS

would be even stronger:

    (v2(n+c1), v2(n+c2))
        ->
    depth

with no explicit p,q,X,Y.

The theoretically important first candidates are:

    (-9,-3)
    (-9,3)
    (-9,9)
    (-3,3)
    (-3,9)
    (3,9).

"""
    )


# ==============================================================================
# MAIN
# ==============================================================================

def main() -> None:

    print("=" * 90)
    print(
        "EXPERIMENT 623 START"
    )
    print("=" * 90)
    print()
    print(
        "FAST TWO-STAGE n-ONLY 2-ADIC SIGNATURE SEARCH"
    )
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
    print()

    print(
        "[3] BUILD TARGETS"
    )

    m_values, depth_values = (
        build_targets(
            states
        )
    )

    print(
        f"    states={len(states)}"
    )
    print()

    offsets = OFFSET_POOL

    print(
        "[4] PRECOMPUTE n-ONLY VALUATIONS"
    )

    columns = build_valuation_columns(
        states,
        offsets,
    )

    print(
        f"    offsets={len(offsets)}"
    )
    print()

    single_offset_baseline(
        columns,
        offsets,
        m_values,
        depth_values,
    )

    sample_indices = make_sample(
        SAMPLE_SIZE,
        len(states),
    )

    print()
    print("=" * 90)
    print(
        "TEST 3: FAST PAIR RANKING"
    )
    print("=" * 90)

    print(
        f"    sample size="
        f"{len(sample_indices)}"
    )

    candidates_m = rank_pairs(
        columns,
        offsets,
        m_values,
        sample_indices,
        "m",
    )

    candidates_depth = rank_pairs(
        columns,
        offsets,
        depth_values,
        sample_indices,
        "depth",
    )

    # Merge the best candidates from both searches.
    merged = []

    seen = set()

    for pair, _, _ in (
        candidates_m
        + candidates_depth
    ):

        if pair in seen:
            continue

        seen.add(pair)
        merged.append(
            pair
        )

    merged_candidates = [
        (
            pair,
            0,
            0,
        )
        for pair in merged[:TOP_EXACT]
    ]

    print()
    print("=" * 90)
    print(
        "TEST 4: EXACT TOP-CANDIDATE VERIFICATION"
    )
    print("=" * 90)

    exact_m, exact_depth = (
        verify_top_pairs(
            merged_candidates,
            columns,
            m_values,
            depth_values,
        )
    )

    print()
    print(
        f"EXACT m PAIRS={len(exact_m)}"
    )

    for pair in exact_m:

        print(
            f"    {pair}"
        )

    print(
        f"EXACT DEPTH PAIRS="
        f"{len(exact_depth)}"
    )

    for pair in exact_depth:

        print(
            f"    {pair}"
        )

    test_theoretical_pairs(
        columns,
        m_values,
        depth_values,
    )

    test_residue_signature(
        states,
        columns,
        m_values,
        depth_values,
    )

    test_global_crosscheck(
        states,
        m_values,
        depth_values,
    )

    examples(
        states,
        columns,
        m_values,
        depth_values,
    )

    summary(
        exact_m,
        exact_depth,
    )

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 623 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()
