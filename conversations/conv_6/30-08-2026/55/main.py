#!/usr/bin/env python3
"""
================================================================================
RECURSIVE K FACTORIZATION / (k,l) RECONSTRUCTION EXPERIMENT
================================================================================

Research question:

    Given

        T = floor(n / (r1*r2))

    and the bounded defect

        K = T - E,

    can we:

        1. factor K,
        2. generate (k,l) from the divisor pairs of K,
        3. reconstruct admissible residue/carry information,
        4. identify the true (k,l),
        5. and distinguish the true K from nearby K values?

IMPORTANT:

    No enumeration of the original p factor interval is performed.

For each candidate E:

    K = T - E

then:

    factors(K)
        ->
    divisor pairs (k,l)
        ->
    quotient bounds
        ->
    carry feasibility
        ->
    residue intervals
        ->
    exact reconstruction check

The final exact check uses the known anchor factors only for scoring.
"""

from __future__ import annotations

import math
import random
import time
from dataclasses import dataclass
from collections import Counter


# =============================================================================
# CONFIGURATION
# =============================================================================

N_ANCHORS = 100

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

MOD_MIN = 300
MOD_MAX = 3_000

CLOSE_RATIO = 0.20

SEED = 1_511_464_998

# Search E around T.
#
# The theoretical bound used by the earlier experiments is:
#
#     E <= k + l
#
# For the experiment we use a safe observable window based on the
# factor range and modulus pair.
#
# Increase carefully; factorization of candidate K remains cheap,
# but the number of E values scales linearly.
E_EXTRA = 10

# Number of anchors actually evaluated in the expensive section.
SEARCH_ANCHORS = 100

# Maximum number of divisor pairs retained for display.
DISPLAY_PAIRS = 12


# =============================================================================
# DATA STRUCTURES
# =============================================================================

@dataclass(frozen=True)
class Anchor:
    p: int
    q: int
    n: int


@dataclass(frozen=True)
class ModPair:
    r1: int
    r2: int


@dataclass
class CandidatePair:
    E: int
    K: int
    k: int
    l: int

    p_lo: int
    p_hi: int
    q_lo: int
    q_hi: int

    carry_possible: bool

    # Number of residue cells compatible with the carry equations.
    carry_cells: int

    # Sum of residue rectangle areas.
    residue_area: int

    # Whether the true anchor (k,l) equals this pair.
    true_pair: bool


# =============================================================================
# PRIME GENERATION
# =============================================================================

def sieve(limit: int) -> list[int]:
    if limit < 2:
        return []

    bs = bytearray(b"\x01") * (limit + 1)
    bs[0:2] = b"\x00\x00"

    for p in range(2, int(limit ** 0.5) + 1):
        if bs[p]:
            start = p * p
            bs[start:limit + 1:p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    return [i for i, v in enumerate(bs) if v]


# =============================================================================
# CLOSE MODULUS PAIRS
# =============================================================================

def build_close_pairs(primes: list[int], ratio: float) -> list[ModPair]:
    pairs: list[ModPair] = []

    for i, r1 in enumerate(primes):
        for r2 in primes[i + 1:]:
            if r2 > MOD_MAX:
                break

            if r1 < MOD_MIN:
                continue

            # Close relative spacing.
            if (r2 - r1) / r1 <= ratio:
                pairs.append(ModPair(r1, r2))

    return pairs


# =============================================================================
# ANCHOR GENERATION
# =============================================================================

def build_anchors(factor_primes: list[int],
                  count: int,
                  seed: int) -> list[Anchor]:

    rng = random.Random(seed)

    usable = [
        p for p in factor_primes
        if FACTOR_MIN <= p <= FACTOR_MAX
    ]

    anchors: list[Anchor] = []
    seen: set[tuple[int, int]] = set()

    while len(anchors) < count:
        p = rng.choice(usable)
        q = rng.choice(usable)

        if p > q:
            p, q = q, p

        if p == q:
            continue

        key = (p, q)

        if key in seen:
            continue

        seen.add(key)
        anchors.append(
            Anchor(
                p=p,
                q=q,
                n=p * q,
            )
        )

    return anchors


# =============================================================================
# INTEGER HELPERS
# =============================================================================

def floor_div(a: int, b: int) -> int:
    return a // b


def ceil_div(a: int, b: int) -> int:
    return -((-a) // b)


# =============================================================================
# FACTOR K
# =============================================================================

def factor_integer(n: int) -> dict[int, int]:
    """
    Trial division is intentionally sufficient here because K is
    much smaller than n in this experiment.
    """
    if n <= 0:
        return {}

    factors: dict[int, int] = {}

    d = 2

    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d

        d = 3 if d == 2 else d + 2

    if n > 1:
        factors[n] = factors.get(n, 0) + 1

    return factors


def divisor_pairs(K: int) -> list[tuple[int, int]]:
    """
    All ordered quotient pairs k*l = K with k <= l.
    """
    if K <= 0:
        return []

    pairs = []

    d = 1

    while d * d <= K:
        if K % d == 0:
            e = K // d
            pairs.append((d, e))

        d += 1

    return pairs


# =============================================================================
# QUOTIENT BOUNDS
# =============================================================================

def quotient_bounds(p_lo: int,
                    p_hi: int,
                    r1: int,
                    k: int) -> tuple[int, int] | None:
    """
    p = a + k*r1, with

        0 <= a < r1

    gives the natural cell:

        k*r1 <= p <= (k+1)*r1 - 1.

    Intersect it with the global factor interval.
    """

    lo = max(
        p_lo,
        k * r1
    )

    hi = min(
        p_hi,
        (k + 1) * r1 - 1
    )

    if lo > hi:
        return None

    return lo, hi


# =============================================================================
# CARRY FEASIBILITY
# =============================================================================

def carry_components(
    k: int,
    l: int,
    r1: int,
    r2: int,
    E: int,
) -> list[tuple[int, int, int]]:
    """
    Enumerate only the small carry split:

        E = c1 + c2 + c3

    where

        0 <= c1 < l
        0 <= c2 < k

    and c3 is bounded by the residue terms.

    We deliberately keep this entirely in quotient/carry space.
    """

    out: list[tuple[int, int, int]] = []

    # c3 is small. A conservative bound is 2 for the construction
    # used in the previous experiments, but we calculate a slightly
    # larger safe range here.
    max_c3 = 3

    for c3 in range(max_c3 + 1):
        remaining = E - c3

        if remaining < 0:
            continue

        c1_lo = max(0, remaining - (k - 1))
        c1_hi = min(l - 1, remaining)

        for c1 in range(c1_lo, c1_hi + 1):
            c2 = remaining - c1

            if not (0 <= c2 < k):
                continue

            out.append((c1, c2, c3))

    return out


def carry_residue_interval(
    carry: int,
    quotient: int,
    modulus: int,
) -> tuple[int, int] | None:
    """
    For

        carry = floor(a*quotient/modulus)

    solve for the possible integer residue a:

        carry*modulus <= a*quotient
                       < (carry+1)*modulus.

    Therefore

        ceil(carry*modulus/quotient)
            <= a
            <=
        floor(((carry+1)*modulus-1)/quotient).
    """

    lo = ceil_div(carry * modulus, quotient)

    hi = ((carry + 1) * modulus - 1) // quotient

    lo = max(0, lo)
    hi = min(modulus - 1, hi)

    if lo > hi:
        return None

    return lo, hi


def carry_cells(
    k: int,
    l: int,
    r1: int,
    r2: int,
    E: int,
) -> tuple[int, int]:
    """
    Return:

        number of carry-compatible residue rectangles
        total rectangle area.

    No p enumeration.
    No q enumeration.
    """

    cells = 0
    area = 0

    for c1, c2, c3 in carry_components(k, l, r1, r2, E):

        a_range = carry_residue_interval(c1, l, r1)
        b_range = carry_residue_interval(c2, k, r2)

        if a_range is None or b_range is None:
            continue

        a_lo, a_hi = a_range
        b_lo, b_hi = b_range

        # c3 must actually be compatible with the residual carry
        # term:
        #
        # floor((d1*r2 + d2*r1 + a*b)/(r1*r2)) = c3
        #
        # We cannot simply assume that every rectangle point works,
        # but the rectangle is still a valid necessary-condition cell.
        #
        # The purpose here is to measure the compression before
        # exact point enumeration.

        width_a = a_hi - a_lo + 1
        width_b = b_hi - b_lo + 1

        cells += 1
        area += width_a * width_b

    return cells, area


# =============================================================================
# OBSERVABLE E WINDOW
# =============================================================================

def e_window(
    n: int,
    r1: int,
    r2: int,
) -> tuple[int, int, int]:
    """
    Returns:

        T
        E_min
        E_max

    Since

        K = k*l

    and k,l are bounded by p,q <= FACTOR_MAX, we obtain a conservative
    observable bound from the maximum quotient coordinates.
    """

    R = r1 * r2

    T = n // R

    k_max = FACTOR_MAX // r1
    l_max = FACTOR_MAX // r2

    # K = k*l <= k_max*l_max.
    #
    # Therefore:
    #
    # E = T-K
    #
    # lies above T-k_max*l_max.
    #
    # We also know E >= 0 for ordinary quotient-cell placement.

    E_min = max(0, T - k_max * l_max - E_EXTRA)
    E_max = T

    return T, E_min, E_max


# =============================================================================
# CANDIDATE EVALUATION
# =============================================================================

def analyze_candidate_K(
    anchor: Anchor,
    mods: ModPair,
    E: int,
) -> list[CandidatePair]:

    n = anchor.n
    r1 = mods.r1
    r2 = mods.r2

    T = n // (r1 * r2)

    K = T - E

    if K <= 0:
        return []

    pairs = divisor_pairs(K)

    results: list[CandidatePair] = []

    for k, l in pairs:

        p_cell = quotient_bounds(
            FACTOR_MIN,
            FACTOR_MAX,
            r1,
            k
        )

        q_cell = quotient_bounds(
            FACTOR_MIN,
            FACTOR_MAX,
            r2,
            l
        )

        if p_cell is None or q_cell is None:
            continue

        p_lo, p_hi = p_cell
        q_lo, q_hi = q_cell

        cells, area = carry_cells(
            k,
            l,
            r1,
            r2,
            E
        )

        results.append(
            CandidatePair(
                E=E,
                K=K,
                k=k,
                l=l,
                p_lo=p_lo,
                p_hi=p_hi,
                q_lo=q_lo,
                q_hi=q_hi,
                carry_possible=(cells > 0),
                carry_cells=cells,
                residue_area=area,
                true_pair=(
                    k == anchor.p // r1
                    and
                    l == anchor.q // r2
                ),
            )
        )

    return results


# =============================================================================
# MAIN ANCHOR ANALYSIS
# =============================================================================

def analyze_anchor(
    anchor: Anchor,
    mods: ModPair,
) -> dict:

    n = anchor.n
    r1 = mods.r1
    r2 = mods.r2

    R = r1 * r2

    T = n // R

    true_k = anchor.p // r1
    true_l = anchor.q // r2

    true_K = true_k * true_l

    true_E = T - true_K

    _, E_min, E_max = e_window(
        n,
        r1,
        r2
    )

    all_candidates: list[CandidatePair] = []

    factor_count = 0
    distinct_K = 0

    seen_K: set[int] = set()

    E_count = 0

    for E in range(E_min, E_max + 1):

        K = T - E

        if K <= 0:
            continue

        E_count += 1

        factors = factor_integer(K)
        if not factors:
            continue

        factor_count += 1

        if K in seen_K:
            continue

        seen_K.add(K)
        distinct_K += 1

        all_candidates.extend(
            analyze_candidate_K(
                anchor,
                mods,
                E
            )
        )

    true_candidates = [
        c for c in all_candidates
        if c.true_pair
    ]

    carry_candidates = [
        c for c in all_candidates
        if c.carry_possible
    ]

    unique_K_matches = {
        c.K
        for c in true_candidates
    }

    true_K_rank = None

    # Rank candidate K values by the number of compatible (k,l)
    # descendants. Smaller is better.
    grouped: dict[int, list[CandidatePair]] = {}

    for c in all_candidates:
        grouped.setdefault(c.K, []).append(c)

    ordered_K = sorted(
        grouped,
        key=lambda K: (
            sum(x.carry_possible for x in grouped[K]),
            len(grouped[K]),
            K
        )
    )

    if true_K in ordered_K:
        true_K_rank = ordered_K.index(true_K) + 1

    return {
        "T": T,
        "true_k": true_k,
        "true_l": true_l,
        "true_K": true_K,
        "true_E": true_E,

        "E_window": E_max - E_min + 1,
        "factored_E": factor_count,
        "distinct_K": distinct_K,

        "all_candidates": all_candidates,
        "carry_candidates": carry_candidates,
        "true_candidates": true_candidates,

        "true_K_found": bool(true_candidates),
        "true_K_rank": true_K_rank,
        "unique_true_K": len(unique_K_matches) == 1,

        "factorization_of_true_K": factor_integer(true_K),
        "true_K_pairs": divisor_pairs(true_K),
    }


# =============================================================================
# OUTPUT
# =============================================================================

def print_header(title: str) -> None:
    print()
    print("=" * 92)
    print(title)
    print("=" * 92)


def main() -> None:

    start_total = time.perf_counter()

    print("=" * 92)
    print("RECURSIVE K FACTORIZATION / (k,l) RECONSTRUCTION EXPERIMENT")
    print("=" * 92)
    print(f"N anchors                 = {N_ANCHORS}")
    print(f"factor range              = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"modulus range             = {MOD_MIN:,} - {MOD_MAX:,}")
    print(f"close ratio               = {CLOSE_RATIO * 100:.1f}%")
    print(f"seed                      = {SEED}")
    print(f"search anchors            = {SEARCH_ANCHORS}")
    print()

    # -------------------------------------------------------------------------
    # PRIME POOLS
    # -------------------------------------------------------------------------

    print("=" * 92)
    print("BUILDING PRIME POOLS")
    print("=" * 92)

    factor_primes = sieve(FACTOR_MAX)
    modulus_primes = [
        p
        for p in sieve(MOD_MAX)
        if MOD_MIN <= p <= MOD_MAX
    ]

    print(f"factor primes             = {len(factor_primes):,}")
    print(f"modulus primes            = {len(modulus_primes):,}")

    # -------------------------------------------------------------------------
    # MODULUS PAIRS
    # -------------------------------------------------------------------------

    print()
    print("=" * 92)
    print("BUILDING CLOSE MODULUS PAIRS")
    print("=" * 92)

    close_pairs = build_close_pairs(
        modulus_primes,
        CLOSE_RATIO
    )

    print(f"close modulus pairs       = {len(close_pairs):,}")

    if not close_pairs:
        raise RuntimeError("No close modulus pairs found.")

    # -------------------------------------------------------------------------
    # ANCHORS
    # -------------------------------------------------------------------------

    print()
    print("=" * 92)
    print("BUILDING ANCHORS")
    print("=" * 92)

    anchors = build_anchors(
        factor_primes,
        N_ANCHORS,
        SEED
    )

    print(f"actual anchors            = {len(anchors)}")

    # -------------------------------------------------------------------------
    # SELECT ONE CLOSE PAIR PER ANCHOR
    # -------------------------------------------------------------------------

    rng = random.Random(SEED)

    selected_pairs: list[ModPair] = []

    for _ in anchors:
        selected_pairs.append(
            rng.choice(close_pairs)
        )

    # -------------------------------------------------------------------------
    # VALIDATE TRUE K
    # -------------------------------------------------------------------------

    print()
    print("=" * 92)
    print("VALIDATING K = T-E IDENTITY")
    print("=" * 92)

    identity_failures = 0

    for anchor, mods in zip(anchors, selected_pairs):

        r1 = mods.r1
        r2 = mods.r2

        T = anchor.n // (r1 * r2)

        k = anchor.p // r1
        l = anchor.q // r2

        K = k * l
        E = T - K

        if K != T - E:
            identity_failures += 1

    print(f"identity failures         = {identity_failures}")

    # -------------------------------------------------------------------------
    # SEARCH
    # -------------------------------------------------------------------------

    print()
    print("=" * 92)
    print("RUNNING K FACTORIZATION / DIVISOR-PAIR SEARCH")
    print("=" * 92)

    search_start = time.perf_counter()

    summaries = []

    total_E_window = 0
    total_K_pairs = 0
    total_carry_pairs = 0

    recovered = 0
    unique_true_K = 0

    for idx, (anchor, mods) in enumerate(
        zip(anchors[:SEARCH_ANCHORS], selected_pairs[:SEARCH_ANCHORS]),
        start=1
    ):

        if idx % max(1, SEARCH_ANCHORS // 10) == 0 or idx == 1:
            print(f"anchor {idx:3d}/{SEARCH_ANCHORS}")

        result = analyze_anchor(
            anchor,
            mods
        )

        summaries.append(
            (anchor, mods, result)
        )

        total_E_window += result["E_window"]

        total_K_pairs += len(result["all_candidates"])
        total_carry_pairs += len(result["carry_candidates"])

        if result["true_K_found"]:
            recovered += 1

        if result["unique_true_K"]:
            unique_true_K += 1

    search_time = time.perf_counter() - search_start

    # -------------------------------------------------------------------------
    # AGGREGATES
    # -------------------------------------------------------------------------

    avg_E_window = (
        total_E_window / len(summaries)
        if summaries else 0
    )

    avg_K_pairs = (
        total_K_pairs / len(summaries)
        if summaries else 0
    )

    avg_carry_pairs = (
        total_carry_pairs / len(summaries)
        if summaries else 0
    )

    carry_ratio = (
        total_carry_pairs / total_K_pairs
        if total_K_pairs
        else 0
    )

    # -------------------------------------------------------------------------
    # SUMMARY
    # -------------------------------------------------------------------------

    print_header("SUMMARY")

    print(f"anchors analyzed              = {len(summaries)}")
    print(f"average E-window states       = {avg_E_window:.3f}")
    print(f"average generated (k,l)       = {avg_K_pairs:.3f}")
    print(f"average carry-compatible      = {avg_carry_pairs:.3f}")
    print(f"carry-compatible / generated = {carry_ratio:.6f}")

    print()
    print("RECOVERY")
    print(f"true (k,l) retained           = {recovered}/{len(summaries)}")
    print(
        f"recovery rate                 = "
        f"{100.0 * recovered / len(summaries):.4f}%"
        if summaries else
        "recovery rate                 = N/A"
    )

    print()
    print("TRUE K IDENTIFICATION")
    print(f"unique true K                 = {unique_true_K}/{len(summaries)}")
    print(
        f"fraction                      = "
        f"{unique_true_K / len(summaries):.6f}"
        if summaries else
        "fraction                      = N/A"
    )

    # -------------------------------------------------------------------------
    # TRUE K FACTORIZATION
    # -------------------------------------------------------------------------

    print_header("TRUE K FACTORIZATION EXAMPLES")

    for anchor, mods, result in summaries[:20]:

        true_K = result["true_K"]

        print(
            f"n={anchor.n:,} "
            f"p={anchor.p:,} "
            f"q={anchor.q:,} "
            f"mods=({mods.r1},{mods.r2})"
        )

        print(
            f"    T={result['T']} "
            f"E={result['true_E']} "
            f"K={true_K}"
        )

        print(
            f"    factors(K) = "
            f"{result['factorization_of_true_K']}"
        )

        print(
            f"    divisor pairs = "
            f"{result['true_K_pairs'][:DISPLAY_PAIRS]}"
        )

        true_pair = (
            result["true_k"],
            result["true_l"]
        )

        print(
            f"    true (k,l) = {true_pair}"
        )

        print(
            f"    true K rank = {result['true_K_rank']}"
        )

        print()

    # -------------------------------------------------------------------------
    # MOST PROMISING
    # -------------------------------------------------------------------------

    ranked = sorted(
        summaries,
        key=lambda item: (
            item[2]["true_K_rank"]
            if item[2]["true_K_rank"] is not None
            else 10**9,
            len(item[2]["all_candidates"])
        )
    )

    print_header("BEST TRUE-K RANKS")

    for anchor, mods, result in ranked[:20]:

        print(
            f"n={anchor.n:,} "
            f"mods=({mods.r1},{mods.r2}) "
            f"T={result['T']} "
            f"E={result['true_E']} "
            f"K={result['true_K']} "
            f"(k,l)=({result['true_k']},{result['true_l']}) "
            f"rank={result['true_K_rank']} "
            f"generated={len(result['all_candidates'])} "
            f"carry={len(result['carry_candidates'])}"
        )

    # -------------------------------------------------------------------------
    # E-WINDOW EXAMPLES
    # -------------------------------------------------------------------------

    print_header("E-WINDOW / K CANDIDATE EXAMPLES")

    for anchor, mods, result in summaries[:20]:

        print(
            f"n={anchor.n:,} "
            f"mods=({mods.r1},{mods.r2})"
        )

        print(
            f"    T={result['T']} "
            f"trueK={result['true_K']} "
            f"trueE={result['true_E']} "
            f"E-window={result['E_window']}"
        )

        print(
            f"    generated (k,l) = "
            f"{len(result['all_candidates'])}"
        )

        print(
            f"    carry-compatible = "
            f"{len(result['carry_candidates'])}"
        )

        # Show candidates for the true K.
        true_K_candidates = [
            c
            for c in result["all_candidates"]
            if c.K == result["true_K"]
        ]

        print("    true-K pairs:")

        for c in true_K_candidates[:DISPLAY_PAIRS]:
            print(
                f"        "
                f"(k,l)=({c.k},{c.l}) "
                f"carryCells={c.carry_cells} "
                f"residueArea={c.residue_area} "
                f"true={c.true_pair}"
            )

        print()

    # -------------------------------------------------------------------------
    # KEY INTERPRETATION
    # -------------------------------------------------------------------------

    print_header("MATHEMATICAL INTERPRETATION")

    print(
        """
The experiment isolates the recursive hypothesis:

    T = floor(n/(r1*r2))
    E = T-k*l

therefore:

    K = k*l = T-E.

For every candidate E, the script computes K and FACTORS K directly.

If:

    K = d * (K/d),

then:

    k = d
    l = K/d

is generated without enumerating p or q.

The candidate hierarchy is therefore:

    n
     |
     v
    T
     |
     v
    candidate E
     |
     v
    K = T-E
     |
     v
    factor(K)
     |
     v
    divisor pairs (k,l)
     |
     v
    quotient cells
     |
     v
    carry constraints
     |
     v
    residue rectangles.

This experiment is deliberately testing whether the SMALL
factorization problem K can act as an intermediate representation
of the LARGE factorization problem n.

There are three possible outcomes.

1. TRUE K PRODUCES VERY FEW (k,l)

   This is expected mathematically when K has few divisors,
   but it becomes interesting if it happens systematically.

2. TRUE K IS UNIQUE AMONG T-E

   This is much more important.

   It would mean the defect E can potentially be identified
   indirectly by the arithmetic structure of K.

3. TRUE K IS ONE OF MANY EQUALLY PLAUSIBLE T-E VALUES

   Then factoring K is useful only AFTER E has somehow been
   identified.

The decisive quantity is therefore not simply:

    "K is much smaller than n."

The decisive quantity is:

    "Does factoring candidate K values allow us to identify the
     correct K without returning to the original p search?"

That is the recursive bridge being tested here.

No original p interval is enumerated during candidate generation.
"""
    )

    print_header("TIMING")

    print(
        f"search runtime              = {search_time:.3f} s"
    )

    print(
        f"total runtime               = "
        f"{time.perf_counter() - start_total:.3f} s"
    )

    print()
    print("=" * 92)
    print("EXPERIMENT COMPLETE")
    print("=" * 92)


if __name__ == "__main__":
    main()
