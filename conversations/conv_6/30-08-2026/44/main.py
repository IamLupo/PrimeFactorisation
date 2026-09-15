#!/usr/bin/env python3
"""
===============================================================================
MULTI-(r1,r2) QUOTIENT-PRODUCT INTERSECTION EXPERIMENT
===============================================================================

Research target:

    T_i = floor(n / (r1_i * r2_i))

and

    k_i*l_i = T_i - E_i

with the exact bound

    0 <= E_i <= k_i + l_i.

For several independent modulus pairs, generate candidate quotient pairs
from the short E-window and intersect the implied p/q information.

This experiment asks:

    Can several independent quotient-product windows reconstruct the
    same p,q without enumerating all candidate primes?

Important:
    This is an experimental structural test, not a claim of a new
    factorization algorithm.

All final candidates are checked with:

    p*q == n

===============================================================================
"""

from __future__ import annotations

import math
import random
import statistics
import time
from dataclasses import dataclass
from typing import Iterable


# =============================================================================
# CONFIGURATION
# =============================================================================

N_ANCHORS = 100

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

MOD_MIN = 300
MOD_MAX = 3_000

NUM_PAIRS = 8

SEED = 1_511_464_998

# Number of E values inspected around T.
# The rigorous upper bound is derived from the factor range.
MAX_E_EXTRA = 0

# Only retain divisors k that can correspond to p in the factor interval.
# This avoids factoring irrelevant quotient branches.
MAX_DIVISOR_SCAN = 200_000


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


@dataclass(frozen=True)
class QuotientCandidate:
    pair_index: int
    E: int
    T: int
    product: int
    k: int
    ell: int
    p_interval_lo: int
    p_interval_hi: int
    q_interval_lo: int
    q_interval_hi: int


# =============================================================================
# ANCHOR CONSTRUCTION
# =============================================================================

def build_anchors(
    factor_primes: list[int],
    count: int,
    seed: int,
) -> list[Anchor]:

    rng = random.Random(seed)

    anchors: list[Anchor] = []

    # Prefer broadly distributed factors so that we cover both
    # low/high and balanced/unbalanced cases.
    candidates = []

    for _ in range(count * 20):
        p, q = rng.sample(factor_primes, 2)

        if p > q:
            p, q = q, p

        candidates.append((p, q))

    # Deduplicate.
    seen = set()

    for p, q in candidates:
        if (p, q) in seen:
            continue

        seen.add((p, q))
        anchors.append(Anchor(p=p, q=q, n=p * q))

        if len(anchors) >= count:
            break

    if len(anchors) < count:
        raise RuntimeError("Could not construct enough anchors.")

    return anchors


# =============================================================================
# MODULUS PAIRS
# =============================================================================

def select_modulus_pairs(
    modulus_primes: list[int],
    n: int,
    count: int,
    rng: random.Random,
) -> list[ModPair]:

    usable = []

    for _ in range(count * 100):
        r1, r2 = rng.sample(modulus_primes, 2)

        # Keep ordering deterministic.
        if r1 > r2:
            r1, r2 = r2, r1

        # Need the product below n.
        if r1 * r2 >= n:
            continue

        # Avoid pathological gcd relationships.
        if math.gcd(r1, r2) != 1:
            continue

        usable.append(ModPair(r1, r2))

        if len(usable) >= count:
            break

    # Deduplicate while preserving order.
    out: list[ModPair] = []
    seen = set()

    for x in usable:
        key = (x.r1, x.r2)

        if key not in seen:
            seen.add(key)
            out.append(x)

    if len(out) < count:
        raise RuntimeError("Could not construct enough modulus pairs.")

    return out


# =============================================================================
# QUOTIENT BOUNDS
# =============================================================================

def quotient_bounds(
    r1: int,
    r2: int,
) -> tuple[int, int, int, int]:
    """
    Since:

        p = k*r1 + a, 0 <= a < r1

    and:

        FACTOR_MIN <= p <= FACTOR_MAX,

    obtain integer bounds on k.

    Same for ell.
    """

    k_min = max(0, (FACTOR_MIN + r1 - 1) // r1)
    k_max = FACTOR_MAX // r1

    l_min = max(0, (FACTOR_MIN + r2 - 1) // r2)
    l_max = FACTOR_MAX // r2

    return k_min, k_max, l_min, l_max


def maximum_E(
    r1: int,
    r2: int,
) -> int:

    k_min, k_max, l_min, l_max = quotient_bounds(r1, r2)

    # Exact safe upper bound:
    #
    # E <= k + ell
    #
    # so globally:
    #
    # E <= k_max + l_max

    return k_max + l_max


# =============================================================================
# DIVISORS OF AN INTEGER
# =============================================================================

def divisor_pairs_in_range(
    value: int,
    k_min: int,
    k_max: int,
    l_min: int,
    l_max: int,
) -> list[tuple[int, int]]:

    if value <= 0:
        return []

    out = []

    # We only need k up to sqrt(value).
    limit = min(
        int(math.isqrt(value)),
        k_max,
        MAX_DIVISOR_SCAN,
    )

    for k in range(max(1, k_min), limit + 1):

        if value % k != 0:
            continue

        ell = value // k

        if not (l_min <= ell <= l_max):
            continue

        out.append((k, ell))

    return out


# =============================================================================
# SINGLE MODULUS-PAIR ANALYSIS
# =============================================================================

def generate_quotient_candidates(
    anchor: Anchor,
    pair_index: int,
    pair: ModPair,
) -> tuple[int, list[QuotientCandidate], int]:

    n = anchor.n
    r1 = pair.r1
    r2 = pair.r2

    T = n // (r1 * r2)

    k_min, k_max, l_min, l_max = quotient_bounds(r1, r2)

    e_max = min(
        maximum_E(r1, r2) + MAX_E_EXTRA,
        T,
    )

    candidates: list[QuotientCandidate] = []

    # The true relation is:
    #
    #     k*ell = T-E
    #
    # for some E in [0,e_max].
    for E in range(e_max + 1):

        product = T - E

        if product <= 0:
            continue

        pairs = divisor_pairs_in_range(
            product,
            k_min,
            k_max,
            l_min,
            l_max,
        )

        for k, ell in pairs:

            # Quotient cells imply intervals for p and q.
            p_lo = max(
                FACTOR_MIN,
                k * r1,
            )

            p_hi = min(
                FACTOR_MAX,
                (k + 1) * r1 - 1,
            )

            q_lo = max(
                FACTOR_MIN,
                ell * r2,
            )

            q_hi = min(
                FACTOR_MAX,
                (ell + 1) * r2 - 1,
            )

            candidates.append(
                QuotientCandidate(
                    pair_index=pair_index,
                    E=E,
                    T=T,
                    product=product,
                    k=k,
                    ell=ell,
                    p_interval_lo=p_lo,
                    p_interval_hi=p_hi,
                    q_interval_lo=q_lo,
                    q_interval_hi=q_hi,
                )
            )

    return e_max, candidates, T


# =============================================================================
# INTERSECTION HELPERS
# =============================================================================

def interval_intersection(
    a_lo: int,
    a_hi: int,
    b_lo: int,
    b_hi: int,
) -> tuple[int, int] | None:

    lo = max(a_lo, b_lo)
    hi = min(a_hi, b_hi)

    if lo > hi:
        return None

    return lo, hi


def compatible_candidates(
    left: QuotientCandidate,
    right: QuotientCandidate,
) -> bool:

    p_ok = interval_intersection(
        left.p_interval_lo,
        left.p_interval_hi,
        right.p_interval_lo,
        right.p_interval_hi,
    )

    if p_ok is None:
        return False

    q_ok = interval_intersection(
        left.q_interval_lo,
        left.q_interval_hi,
        right.q_interval_lo,
        right.q_interval_hi,
    )

    return q_ok is not None


# =============================================================================
# MULTI-PAIR INTERSECTION
# =============================================================================

def intersect_all_pairs(
    per_pair_candidates: list[list[QuotientCandidate]],
) -> tuple[list[tuple[int, int]], int]:

    """
    Maintain states:

        (p_lo,p_hi,q_lo,q_hi)

    and progressively intersect them.

    This does not assume that the true p/q are uniquely determined
    by one modulus pair.
    """

    states = [
        (
            c.p_interval_lo,
            c.p_interval_hi,
            c.q_interval_lo,
            c.q_interval_hi,
            c,
        )
        for c in per_pair_candidates[0]
    ]

    if not states:
        return [], 0

    peak_states = len(states)

    for candidates in per_pair_candidates[1:]:

        new_states = []

        for p_lo, p_hi, q_lo, q_hi, old_c in states:

            for c in candidates:

                p_int = interval_intersection(
                    p_lo,
                    p_hi,
                    c.p_interval_lo,
                    c.p_interval_hi,
                )

                if p_int is None:
                    continue

                q_int = interval_intersection(
                    q_lo,
                    q_hi,
                    c.q_interval_lo,
                    c.q_interval_hi,
                )

                if q_int is None:
                    continue

                new_states.append(
                    (
                        p_int[0],
                        p_int[1],
                        q_int[0],
                        q_int[1],
                        c,
                    )
                )

        states = new_states

        peak_states = max(peak_states, len(states))

        if not states:
            break

    # Convert interval states into integer candidate regions.
    regions = []

    for p_lo, p_hi, q_lo, q_hi, _ in states:

        regions.append(
            (
                p_lo,
                p_hi,
                q_lo,
                q_hi,
            )
        )

    return regions, peak_states


# =============================================================================
# EXACT RECONSTRUCTION
# =============================================================================

def exact_from_regions(
    n: int,
    regions: Iterable[tuple[int, int, int, int]],
) -> set[tuple[int, int]]:

    solutions: set[tuple[int, int]] = set()

    for p_lo, p_hi, q_lo, q_hi in regions:

        # Never scan huge intervals blindly.
        # We only perform exact testing when the region is small.
        width = p_hi - p_lo + 1

        if width > 1000:
            continue

        for p in range(p_lo, p_hi + 1):

            if n % p != 0:
                continue

            q = n // p

            if not (q_lo <= q <= q_hi):
                continue

            if p * q == n:
                if FACTOR_MIN <= p <= FACTOR_MAX:
                    if FACTOR_MIN <= q <= FACTOR_MAX:

                        a, b = sorted((p, q))
                        solutions.add((a, b))

    return solutions


# =============================================================================
# TRUE E CALCULATION
# =============================================================================

def true_coordinates(
    anchor: Anchor,
    pair: ModPair,
) -> tuple[int, int, int, int]:

    p = anchor.p
    q = anchor.q

    k = p // pair.r1
    ell = q // pair.r2

    T = anchor.n // (pair.r1 * pair.r2)

    E = T - k * ell

    return k, ell, T, E


# =============================================================================
# VALIDATION
# =============================================================================

def validate_pair(
    anchor: Anchor,
    pair: ModPair,
) -> bool:

    k, ell, T, E = true_coordinates(anchor, pair)

    expected = (
        T
        - (
            k * (anchor.q % pair.r2) // pair.r2
            if False
            else 0
        )
    )

    # Direct exact definition.
    return (
        0 <= E <= k + ell
        and k * ell == T - E
        and anchor.p * anchor.q == anchor.n
    )


# =============================================================================
# MAIN
# =============================================================================

def run() -> None:

    total_start = time.perf_counter()

    print("=" * 80)
    print("MULTI-(r1,r2) QUOTIENT-PRODUCT INTERSECTION EXPERIMENT")
    print("=" * 80)
    print(f"N anchors                 = {N_ANCHORS}")
    print(f"factor range              = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"modulus range             = {MOD_MIN:,} - {MOD_MAX:,}")
    print(f"independent modulus pairs = {NUM_PAIRS}")
    print(f"seed                       = {SEED}")
    print()

    all_primes = sieve(FACTOR_MAX)
    factor_primes = [
        p for p in all_primes
        if FACTOR_MIN <= p <= FACTOR_MAX
    ]

    modulus_primes = [
        p for p in sieve(MOD_MAX)
        if p >= MOD_MIN
    ]

    print("=" * 80)
    print("BUILDING PRIME POOLS")
    print("=" * 80)
    print(f"factor primes             = {len(factor_primes):,}")
    print(f"modulus primes            = {len(modulus_primes):,}")
    print()

    anchors = build_anchors(
        factor_primes,
        N_ANCHORS,
        SEED,
    )

    print("=" * 80)
    print("BUILDING ANCHORS")
    print("=" * 80)
    print(f"actual anchors             = {len(anchors)}")
    print()

    # -------------------------------------------------------------------------
    # GLOBAL VALIDATION
    # -------------------------------------------------------------------------

    print("=" * 80)
    print("VALIDATING E-BOUND")
    print("=" * 80)

    validation_failures = 0

    for anchor in anchors:

        rng = random.Random(SEED ^ anchor.n)

        pairs = select_modulus_pairs(
            modulus_primes,
            anchor.n,
            NUM_PAIRS,
            rng,
        )

        for pair in pairs:

            if not validate_pair(anchor, pair):
                validation_failures += 1

    print(f"E-bound validation failures = {validation_failures}")
    print()

    if validation_failures:
        raise RuntimeError(
            "E-bound validation failed."
        )

    # -------------------------------------------------------------------------
    # SEARCH
    # -------------------------------------------------------------------------

    print("=" * 80)
    print("RUNNING MULTI-PAIR INTERSECTION SEARCH")
    print("=" * 80)

    measurements = []

    for idx, anchor in enumerate(anchors, start=1):

        rng = random.Random(SEED ^ anchor.n ^ 0xABCDEF)

        pairs = select_modulus_pairs(
            modulus_primes,
            anchor.n,
            NUM_PAIRS,
            rng,
        )

        per_pair = []

        total_E_window = 0
        total_qp = 0
        true_found_per_pair = 0

        t0 = time.perf_counter()

        for pair_index, pair in enumerate(pairs):

            e_max, candidates, T = generate_quotient_candidates(
                anchor,
                pair_index,
                pair,
            )

            total_E_window += e_max + 1
            total_qp += len(candidates)

            # Check whether the true quotient pair exists.
            true_k, true_l, true_T, true_E = true_coordinates(
                anchor,
                pair,
            )

            if any(
                c.E == true_E
                and c.k == true_k
                and c.ell == true_l
                for c in candidates
            ):
                true_found_per_pair += 1

            per_pair.append(candidates)

        # Drop empty pair layers.
        usable_layers = [x for x in per_pair if x]

        if not usable_layers:
            regions = []
            peak_states = 0
        else:
            regions, peak_states = intersect_all_pairs(
                usable_layers
            )

        solutions = exact_from_regions(
            anchor.n,
            regions,
        )

        recovered = (
            tuple(sorted((anchor.p, anchor.q)))
            in solutions
        )

        elapsed = time.perf_counter() - t0

        measurements.append(
            {
                "E_windows": total_E_window,
                "quotient_candidates": total_qp,
                "regions": len(regions),
                "peak_intersection_states": peak_states,
                "solutions": len(solutions),
                "recovered": recovered,
                "true_pair_covered": (
                    true_found_per_pair == len(pairs)
                ),
                "time": elapsed,
            }
        )

        if idx % 10 == 0:
            print(f"anchor {idx:3d}/{len(anchors)}")

    # -------------------------------------------------------------------------
    # SUMMARY
    # -------------------------------------------------------------------------

    avg = lambda key: statistics.mean(
        m[key] for m in measurements
    )

    exact_recovered = sum(
        1 for m in measurements if m["recovered"]
    )

    print()
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)

    print(
        f"average E-window states       = "
        f"{avg('E_windows'):.3f}"
    )

    print(
        f"average quotient candidates   = "
        f"{avg('quotient_candidates'):.3f}"
    )

    print(
        f"average intersection regions  = "
        f"{avg('regions'):.3f}"
    )

    print(
        f"average peak intersection     = "
        f"{avg('peak_intersection_states'):.3f}"
    )

    print(
        f"average exact solutions       = "
        f"{avg('solutions'):.3f}"
    )

    print()
    print("=" * 80)
    print("RECOVERY")
    print("=" * 80)

    print(
        f"correctly recovered           = "
        f"{exact_recovered}/{len(anchors)}"
    )

    print(
        f"recovery rate                 = "
        f"{100.0 * exact_recovered / len(anchors):.4f}%"
    )

    # -------------------------------------------------------------------------
    # COVERAGE OF THE TRUE QUOTIENT PAIR
    # -------------------------------------------------------------------------

    coverage = sum(
        1 for m in measurements
        if m["true_pair_covered"]
    )

    print()
    print("=" * 80)
    print("TRUE QUOTIENT COVERAGE")
    print("=" * 80)

    print(
        f"all modulus pairs contained true (k,l) = "
        f"{coverage}/{len(anchors)}"
    )

    # -------------------------------------------------------------------------
    # REDUCTION
    # -------------------------------------------------------------------------

    baseline = len(factor_primes)

    print()
    print("=" * 80)
    print("SEARCH REDUCTION")
    print("=" * 80)

    print(
        f"factor-prime baseline            = {baseline:,}"
    )

    print(
        f"quotient candidates / baseline   = "
        f"{avg('quotient_candidates') / baseline:.9f}"
    )

    print(
        f"intersection regions / baseline = "
        f"{avg('regions') / baseline:.9f}"
    )

    # -------------------------------------------------------------------------
    # EXAMPLES
    # -------------------------------------------------------------------------

    print()
    print("=" * 80)
    print("EXAMPLES")
    print("=" * 80)

    for anchor, m in list(zip(anchors, measurements))[:20]:

        print(
            f"n={anchor.n:,} "
            f"p={anchor.p:,} "
            f"q={anchor.q:,} "
            f"E-window={m['E_windows']:,} "
            f"quotients={m['quotient_candidates']:,} "
            f"regions={m['regions']} "
            f"exact={m['solutions']} "
            f"recovered={m['recovered']}"
        )

    # -------------------------------------------------------------------------
    # MOST INTERESTING CASES
    # -------------------------------------------------------------------------

    best = sorted(
        zip(anchors, measurements),
        key=lambda x: (
            x[1]["regions"],
            x[1]["quotient_candidates"],
        ),
    )[:20]

    print()
    print("=" * 80)
    print("STRONGEST INTERSECTIONS")
    print("=" * 80)

    for anchor, m in best:

        print(
            f"n={anchor.n:,} "
            f"p={anchor.p:,} "
            f"q={anchor.q:,} "
            f"quotients={m['quotient_candidates']:,} "
            f"regions={m['regions']} "
            f"exact={m['solutions']}"
        )

    # -------------------------------------------------------------------------
    # TIMING
    # -------------------------------------------------------------------------

    total_elapsed = time.perf_counter() - total_start

    print()
    print("=" * 80)
    print("TIMING")
    print("=" * 80)

    print(
        f"search runtime                = "
        f"{sum(m['time'] for m in measurements):.3f} s"
    )

    print(
        f"total runtime                 = "
        f"{total_elapsed:.3f} s"
    )

    # -------------------------------------------------------------------------
    # INTERPRETATION
    # -------------------------------------------------------------------------

    print()
    print("=" * 80)
    print("MATHEMATICAL INTERPRETATION")
    print("=" * 80)

    print(
        """
For each modulus pair:

    T_i = floor(n / (r1_i*r2_i))

and:

    k_i*l_i = T_i - E_i.

The exact bound is:

    0 <= E_i <= k_i + l_i.

Therefore each pair produces a finite collection of possible
quotient products:

    T_i,
    T_i - 1,
    T_i - 2,
    ...
    T_i - E_max.

Each such product is factored only enough to obtain quotient
pairs inside the allowed factor range.

The crucial new operation is the intersection between independent
modulus pairs.

For a candidate (k_i,l_i):

    p is known to lie in

        [k_i*r1_i, (k_i+1)*r1_i - 1]

    and q is known to lie in

        [l_i*r2_i, (l_i+1)*r2_i - 1].

Multiple modulus pairs therefore produce multiple rectangles in
the (p,q) plane.

Their intersection is potentially much smaller than any single
pair's rectangle.

This is fundamentally different from the earlier experiments,
because the goal is now:

    generate quotient candidates
        ->
    intersect coordinate intervals
        ->
    obtain a small direct p,q region.

The important measurements are:

    E-window size
    quotient candidates
    intersection regions
    exact solutions.

The strongest possible result would be:

    several independent pairs
        ->
    very few quotient candidates
        ->
    one or a few tiny (p,q) regions
        ->
    exact factorization.

A negative result is equally informative.

If each additional modulus pair merely multiplies the number of
candidate quotient states before intersection, then the apparent
compression is not useful.

Most importantly, the experiment does NOT assume that the true
factors are prime until the final exact check.

The decisive test is whether the candidate regions become
substantially smaller than the original factor interval without
requiring enumeration of that interval.

Every final solution is checked with:

    p*q == n.
"""
    )


if __name__ == "__main__":
    run()
