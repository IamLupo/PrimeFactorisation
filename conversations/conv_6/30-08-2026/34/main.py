#!/usr/bin/env python3

"""
====================================================================================================
THREE-CLOSE-PRIME HYPERBOLA / RESIDUE-LIFT EXPERIMENT
====================================================================================================

This experiment is a corrected follow-up to the 1D hyperbola-strip test.

IMPORTANT CORRECTIONS
---------------------

1. Every modulus triple must satisfy:

       R = r1*r2*r3 < n

   so that:

       n = R + g

   has a genuine positive gap:

       g = n-R > 0.

2. Different anchors receive different close modulus triples.

3. We do NOT scan every integer p inside a surviving quotient cell.

4. Instead, for each quotient cell (k,l), we use the first modulus
   to produce residue candidates:

       p = k*r1 + a
       0 <= a < r1

   and:

       q = n/p.

   The modular condition

       p*q == n (mod r1)

   becomes:

       a*b == n (mod r1)

   where:

       b = q mod r1.

5. The experiment measures the number of residue candidates that
   survive the quotient-cell geometry before exact divisibility.

The goal is to determine whether the hyperbola reduction and the
close-modulus structure together produce a genuinely small candidate
set WITHOUT reverting to an integer scan through each cell.

====================================================================================================
"""

from __future__ import annotations

import math
import random
import statistics
import time
from dataclasses import dataclass
from itertools import combinations


# ================================================================================================
# CONFIGURATION
# ================================================================================================

M = 111_546_435

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

ANCHORS = 300

MODULUS_MIN = 300
MODULUS_MAX = 3_000

CLOSE_RATIO = 0.20

SEED = 1_511_464_998

PROGRESS_STEP = 25


# ================================================================================================
# DATA STRUCTURES
# ================================================================================================

@dataclass(frozen=True)
class Anchor:
    p: int
    q: int
    n: int


@dataclass
class Result:
    anchor_id: int

    p: int
    q: int
    n: int

    r1: int
    r2: int
    r3: int

    R: int
    gap: int

    k_count: int
    full_kl: int

    strip_tests: int
    product_cells: int

    residue_states: int
    r2_survivors: int
    r3_survivors: int

    exact_candidates: int
    exact_solutions: int

    recovered: bool


# ================================================================================================
# PRIME GENERATION
# ================================================================================================

def sieve(limit: int) -> list[int]:
    if limit < 2:
        return []

    a = bytearray(b"\x01") * (limit + 1)
    a[0:2] = b"\x00\x00"

    root = math.isqrt(limit)

    for p in range(2, root + 1):
        if not a[p]:
            continue

        start = p * p
        count = ((limit - start) // p) + 1

        a[start : limit + 1 : p] = b"\x00" * count

    return [
        i
        for i, flag in enumerate(a)
        if flag
    ]


# ================================================================================================
# ANCHORS
# ================================================================================================

def build_anchors(
    primes: list[int],
    count: int,
    seed: int,
) -> list[Anchor]:

    rng = random.Random(seed)

    result: list[Anchor] = []

    while len(result) < count:

        p, q = rng.sample(primes, 2)

        if p > q:
            p, q = q, p

        result.append(
            Anchor(
                p=p,
                q=q,
                n=p*q,
            )
        )

    return result


# ================================================================================================
# MODULUS SELECTION
# ================================================================================================

def choose_triple(
    n: int,
    modulus_primes: list[int],
    close_ratio: float,
) -> tuple[int, int, int]:

    """
    Choose a close triple whose product is strictly below n.

    The previous experiment accidentally selected the same largest
    triple for every anchor. Here the scale is derived from n.

    We want:

        r1*r2*r3 < n

    while keeping the primes as close as practical.
    """

    # Approximate cubic scale.
    target = max(
        MODULUS_MIN,
        int(n ** (1.0 / 3.0)),
    )

    radius = max(
        20,
        int(target * close_ratio),
    )

    candidates = [
        p
        for p in modulus_primes
        if abs(p - target) <= radius
    ]

    # Remove primes that cannot produce a useful product.
    if len(candidates) < 3:
        raise RuntimeError(
            f"Not enough modulus primes around cube-root scale {target}"
        )

    best = None
    best_score = None

    for r1, r2, r3 in combinations(candidates, 3):

        if r1 >= r2 or r2 >= r3:
            continue

        R = r1 * r2 * r3

        if R >= n:
            continue

        gap = n - R

        # Prefer large R/n while still leaving positive gap.
        ratio = R / n

        # Also prefer a tight triple.
        spread = r3 - r1

        score = (
            ratio,
            -spread,
        )

        if best_score is None or score > best_score:
            best_score = score
            best = (r1, r2, r3)

    if best is None:

        # Search globally for the largest product below n.
        global_candidates = sorted(
            modulus_primes,
            reverse=True,
        )

        best_product = -1

        for triple in combinations(
            sorted(global_candidates[:40]),
            3,
        ):

            r1, r2, r3 = triple
            R = r1*r2*r3

            if R < n and R > best_product:
                best_product = R
                best = triple

    if best is None:
        raise RuntimeError(
            f"Unable to construct R < n for n={n}"
        )

    return best


# ================================================================================================
# MODULAR INVERSE
# ================================================================================================

def inv_mod(a: int, m: int) -> int | None:
    a %= m

    if a == 0:
        return None

    try:
        return pow(a, -1, m)
    except ValueError:
        return None


# ================================================================================================
# HYPERBOLA STRIP
# ================================================================================================

def quotient_l_range(
    n: int,
    r1: int,
    r2: int,
    k: int,
) -> tuple[int, int] | None:

    p_lo = max(
        FACTOR_MIN,
        k*r1,
    )

    p_hi = min(
        FACTOR_MAX,
        (k+1)*r1 - 1,
    )

    if p_lo > p_hi:
        return None

    q_lo = n // p_hi
    q_hi = n // p_lo

    q_lo = max(
        FACTOR_MIN,
        q_lo,
    )

    q_hi = min(
        FACTOR_MAX,
        q_hi,
    )

    if q_lo > q_hi:
        return None

    l_lo = q_lo // r2
    l_hi = q_hi // r2

    return l_lo, l_hi


# ================================================================================================
# PRODUCT CELL
# ================================================================================================

def cell_bounds(
    k: int,
    l: int,
    r1: int,
    r2: int,
) -> tuple[int, int, int, int] | None:

    p_lo = max(
        FACTOR_MIN,
        k*r1,
    )

    p_hi = min(
        FACTOR_MAX,
        (k+1)*r1 - 1,
    )

    q_lo = max(
        FACTOR_MIN,
        l*r2,
    )

    q_hi = min(
        FACTOR_MAX,
        (l+1)*r2 - 1,
    )

    if p_lo > p_hi or q_lo > q_hi:
        return None

    return (
        p_lo,
        p_hi,
        q_lo,
        q_hi,
    )


def product_cell_contains_n(
    bounds: tuple[int, int, int, int],
    n: int,
) -> bool:

    p_lo, p_hi, q_lo, q_hi = bounds

    return (
        p_lo*q_lo <= n <= p_hi*q_hi
    )


# ================================================================================================
# RESIDUE STATE GENERATION
# ================================================================================================

def residue_states_for_cell(
    n: int,
    r1: int,
    r2: int,
    r3: int,
    k: int,
    l: int,
    bounds: tuple[int, int, int, int],
) -> tuple[int, int, int, set[tuple[int, int]]]:

    """
    Generate residue candidates WITHOUT scanning all p in the cell.

    p = k*r1 + a

    Because p is in the r1-sized cell, a ranges over:

        0 <= a < r1.

    For each a we derive the required q residue modulo r1:

        b == n * inverse(a) mod r1.

    Only invertible a are useful for this experiment.

    The resulting p and q are reconstructed from their quotient cells.

    The function returns:

        residue attempts
        r2 survivors
        r3 survivors
        exact candidates
    """

    p_lo, p_hi, q_lo, q_hi = bounds

    attempts = 0
    r2_hits = 0
    r3_hits = 0

    exact_candidates: set[tuple[int, int]] = set()

    n_r1 = n % r1

    # p = k*r1 + a, therefore p mod r1 = a.
    a_lo = p_lo - k*r1
    a_hi = p_hi - k*r1

    a_lo = max(0, a_lo)
    a_hi = min(r1 - 1, a_hi)

    if a_lo > a_hi:
        return (
            0,
            0,
            0,
            exact_candidates,
        )

    for a in range(a_lo, a_hi + 1):

        inv_a = inv_mod(a, r1)

        if inv_a is None:
            continue

        attempts += 1

        # q mod r1 required by p*q == n mod r1.
        b = (n_r1 * inv_a) % r1

        # q lives in the r2 quotient cell:

        # q = l*r2 + c

        # We need:
        #
        #     l*r2 + c == b (mod r1)
        #
        # Since gcd(r1,r2)=1, c is uniquely determined modulo r1.
        #
        # c is restricted to:
        #
        #     0 <= c < r2.
        #
        # We solve:
        #
        #     c == b - l*r2 (mod r1).

        c = (b - l*r2) % r1

        # If c >= r2, this particular residue state cannot fit
        # inside the q cell.
        if c >= r2:
            continue

        q_candidate = l*r2 + c

        if not (q_lo <= q_candidate <= q_hi):
            continue

        # Corresponding p candidate for this a.
        p_candidate = k*r1 + a

        if not (
            p_lo <= p_candidate <= p_hi
        ):
            continue

        r2_hits += 1

        # Exact q relationship is still not assumed.
        # Verify the second modulus independently.
        if (
            (p_candidate * q_candidate - n)
            % r2
            != 0
        ):
            continue

        # Third modulus.
        if (
            (p_candidate * q_candidate - n)
            % r3
            != 0
        ):
            continue

        r3_hits += 1

        exact_candidates.add(
            (
                p_candidate,
                q_candidate,
            )
        )

    return (
        attempts,
        r2_hits,
        r3_hits,
        exact_candidates,
    )


# ================================================================================================
# ANCHOR SEARCH
# ================================================================================================

def run_anchor(
    anchor_id: int,
    anchor: Anchor,
    mods: tuple[int, int, int],
) -> Result:

    n = anchor.n

    r1, r2, r3 = mods

    R = r1*r2*r3
    gap = n-R

    if R >= n:
        raise RuntimeError(
            f"Invalid triple: R >= n for anchor {anchor_id}"
        )

    # k values for p = k*r1+a.
    k_min = FACTOR_MIN // r1

    if k_min*r1 < FACTOR_MIN:
        k_min += 1

    k_max = FACTOR_MAX // r1

    if k_max*r1 > FACTOR_MAX:
        k_max -= 1

    k_values = list(
        range(
            max(0, k_min),
            k_max + 1,
        )
    )

    # Full quotient lattice.
    l_min = FACTOR_MIN // r2
    l_max = FACTOR_MAX // r2

    full_kl = (
        len(k_values)
        * (l_max - l_min + 1)
    )

    strip_tests = 0
    product_cells = 0
    residue_states = 0
    r2_survivors = 0
    r3_survivors = 0

    exact: set[tuple[int, int]] = set()

    for k in k_values:

        lr = quotient_l_range(
            n,
            r1,
            r2,
            k,
        )

        if lr is None:
            continue

        l_lo, l_hi = lr

        for l in range(l_lo, l_hi + 1):

            strip_tests += 1

            bounds = cell_bounds(
                k,
                l,
                r1,
                r2,
            )

            if bounds is None:
                continue

            if not product_cell_contains_n(
                bounds,
                n,
            ):
                continue

            product_cells += 1

            (
                attempts,
                r2_hits,
                r3_hits,
                candidates,
            ) = residue_states_for_cell(
                n,
                r1,
                r2,
                r3,
                k,
                l,
                bounds,
            )

            residue_states += attempts
            r2_survivors += r2_hits
            r3_survivors += r3_hits

            for p, q in candidates:

                # FINAL exact test.
                if p*q == n:
                    exact.add(
                        (min(p,q), max(p,q))
                    )

    recovered = (
        (anchor.p, anchor.q)
        in exact
    )

    return Result(
        anchor_id=anchor_id,
        p=anchor.p,
        q=anchor.q,
        n=n,
        r1=r1,
        r2=r2,
        r3=r3,
        R=R,
        gap=gap,
        k_count=len(k_values),
        full_kl=full_kl,
        strip_tests=strip_tests,
        product_cells=product_cells,
        residue_states=residue_states,
        r2_survivors=r2_survivors,
        r3_survivors=r3_survivors,
        exact_candidates=len(
            [
                x
                for x in exact
            ]
        ),
        exact_solutions=len(exact),
        recovered=recovered,
    )


# ================================================================================================
# MAIN
# ================================================================================================

def run() -> None:

    total_start = time.perf_counter()

    print("=" * 100)
    print("THREE-CLOSE-PRIME HYPERBOLA / RESIDUE-LIFT EXPERIMENT")
    print("=" * 100)
    print(f"M                         = {M:,}")
    print(f"factor range              = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"anchors                   = {ANCHORS}")
    print(f"modulus prime range       = {MODULUS_MIN:,} - {MODULUS_MAX:,}")
    print(f"close ratio               = {CLOSE_RATIO:.0%}")
    print(f"seed                      = {SEED:,}")
    print()

    # ============================================================================================
    # PRIME POOLS
    # ============================================================================================

    print("=" * 100)
    print("BUILDING PRIME POOLS")
    print("=" * 100)

    primes = sieve(FACTOR_MAX)

    factor_primes = [
        p
        for p in primes
        if FACTOR_MIN <= p <= FACTOR_MAX
    ]

    modulus_primes = [
        p
        for p in sieve(MODULUS_MAX)
        if MODULUS_MIN <= p <= MODULUS_MAX
    ]

    print(
        f"factor primes             = "
        f"{len(factor_primes):,}"
    )

    print(
        f"modulus primes            = "
        f"{len(modulus_primes):,}"
    )

    # ============================================================================================
    # ANCHORS
    # ============================================================================================

    print()
    print("=" * 100)
    print("BUILDING ANCHORS")
    print("=" * 100)

    anchors = build_anchors(
        factor_primes,
        ANCHORS,
        SEED,
    )

    print(
        f"actual anchors            = "
        f"{len(anchors):,}"
    )

    # ============================================================================================
    # MODULUS TRIPLES
    # ============================================================================================

    print()
    print("=" * 100)
    print("SELECTING ANCHOR-SPECIFIC CLOSE TRIPLES")
    print("=" * 100)

    triples = []

    for i, anchor in enumerate(
        anchors,
        1,
    ):

        triple = choose_triple(
            anchor.n,
            modulus_primes,
            CLOSE_RATIO,
        )

        triples.append(triple)

        if (
            i % PROGRESS_STEP == 0
            or i == len(anchors)
        ):
            print(
                f"anchor {i:3d}/{len(anchors)}"
            )

    # ============================================================================================
    # MODULUS VALIDATION
    # ============================================================================================

    print()
    print("=" * 100)
    print("VALIDATING MODULUS SCALE")
    print("=" * 100)

    failures = 0

    for anchor, triple in zip(
        anchors,
        triples,
    ):

        r1, r2, r3 = triple
        R = r1*r2*r3

        if not (
            r1 < r2 < r3
        ):
            failures += 1
            continue

        if R >= anchor.n:
            failures += 1
            continue

        if (
            math.gcd(r1,r2) != 1
            or math.gcd(r1,r3) != 1
            or math.gcd(r2,r3) != 1
        ):
            failures += 1
            continue

    print(
        f"modulus validation failures = "
        f"{failures}"
    )

    if failures:
        raise RuntimeError(
            "Modulus validation failed."
        )

    # ============================================================================================
    # TRIPLE DIVERSITY
    # ============================================================================================

    unique_triples = len(set(triples))

    print(
        f"unique modulus triples      = "
        f"{unique_triples}"
    )

    # ============================================================================================
    # SEARCH
    # ============================================================================================

    print()
    print("=" * 100)
    print("RUNNING HYPERBOLA / RESIDUE-LIFT SEARCH")
    print("=" * 100)

    results: list[Result] = []

    search_start = time.perf_counter()

    for i in range(
        len(anchors)
    ):

        result = run_anchor(
            i + 1,
            anchors[i],
            triples[i],
        )

        results.append(result)

        if (
            (i+1) % PROGRESS_STEP == 0
            or i+1 == len(anchors)
        ):
            print(
                f"anchor {i+1:3d}/{len(anchors)}"
            )

    search_time = (
        time.perf_counter()
        - search_start
    )

    # ============================================================================================
    # AVERAGES
    # ============================================================================================

    mean = lambda x: statistics.mean(x) if x else 0

    avg_full = mean(
        [r.full_kl for r in results]
    )

    avg_strip = mean(
        [r.strip_tests for r in results]
    )

    avg_cells = mean(
        [r.product_cells for r in results]
    )

    avg_residue = mean(
        [r.residue_states for r in results]
    )

    avg_r2 = mean(
        [r.r2_survivors for r in results]
    )

    avg_r3 = mean(
        [r.r3_survivors for r in results]
    )

    avg_exact = mean(
        [r.exact_solutions for r in results]
    )

    # ============================================================================================
    # SUMMARY
    # ============================================================================================

    print()
    print("=" * 100)
    print("SUMMARY")
    print("=" * 100)

    print(
        f"anchors analyzed            = "
        f"{len(results)}"
    )

    print(
        f"average full (k,l) pairs    = "
        f"{avg_full:,.3f}"
    )

    print(
        f"average hyperbola strips    = "
        f"{avg_strip:,.3f}"
    )

    print(
        f"average product cells       = "
        f"{avg_cells:,.3f}"
    )

    print(
        f"average residue states      = "
        f"{avg_residue:,.3f}"
    )

    print(
        f"average r2 survivors        = "
        f"{avg_r2:,.3f}"
    )

    print(
        f"average r3 survivors        = "
        f"{avg_r3:,.3f}"
    )

    print(
        f"average exact solutions     = "
        f"{avg_exact:,.3f}"
    )

    # ============================================================================================
    # REDUCTIONS
    # ============================================================================================

    print()
    print("=" * 100)
    print("SEARCH REDUCTION")
    print("=" * 100)

    if avg_full:

        print(
            f"strip / full quotient       = "
            f"{avg_strip / avg_full:.9f}"
        )

        print(
            f"strip reduction              = "
            f"{(1-avg_strip/avg_full)*100:.6f}%"
        )

    if avg_strip:

        print(
            f"cells / strip               = "
            f"{avg_cells / avg_strip:.9f}"
        )

    if avg_cells:

        print(
            f"residue / cells             = "
            f"{avg_residue / avg_cells:.9f}"
        )

        print(
            f"r3 / cells                  = "
            f"{avg_r3 / avg_cells:.9f}"
        )

    # ============================================================================================
    # ACTUAL INTEGER WORK
    # ============================================================================================

    print()
    print("=" * 100)
    print("ACTUAL INTEGER CANDIDATE WORK")
    print("=" * 100)

    print(
        """
The previous experiment scanned every integer p in every surviving
cell.

This experiment does NOT do that.

The important statistic is now:

    residue states

This counts modularly generated candidate states before exact
multiplication.

Compare:

    8,363 factor primes
    versus
    average residue states
    versus
    average r3 survivors.
"""
    )

    print(
        f"factor-prime baseline      = 8,363"
    )

    print(
        f"average residue states     = "
        f"{avg_residue:,.3f}"
    )

    print(
        f"average r3 survivors       = "
        f"{avg_r3:,.3f}"
    )

    if avg_residue:
        print(
            f"residue / prime baseline  = "
            f"{avg_residue / 8363:.9f}"
        )

    if avg_r3:
        print(
            f"r3 / prime baseline       = "
            f"{avg_r3 / 8363:.9f}"
        )

    # ============================================================================================
    # RECOVERY
    # ============================================================================================

    recovered = sum(
        r.recovered
        for r in results
    )

    print()
    print("=" * 100)
    print("RECOVERY")
    print("=" * 100)

    print(
        f"correctly recovered        = "
        f"{recovered}/{len(results)}"
    )

    print(
        f"recovery rate              = "
        f"{recovered/len(results)*100:.4f}%"
    )

    # ============================================================================================
    # GAP
    # ============================================================================================

    gaps = [
        r.gap
        for r in results
    ]

    print()
    print("=" * 100)
    print("GAP STATISTICS")
    print("=" * 100)

    print(
        f"mean gap                    = "
        f"{mean(gaps):,.3f}"
    )

    print(
        f"minimum gap                = "
        f"{min(gaps):,}"
    )

    print(
        f"maximum gap                = "
        f"{max(gaps):,}"
    )

    # Pearson correlation helper.
    def correlation(
        xs: list[float],
        ys: list[float],
    ) -> float:

        if len(xs) != len(ys) or len(xs) < 2:
            return float("nan")

        mx = statistics.mean(xs)
        my = statistics.mean(ys)

        dx = [x-mx for x in xs]
        dy = [y-my for y in ys]

        den1 = math.sqrt(
            sum(x*x for x in dx)
        )

        den2 = math.sqrt(
            sum(y*y for y in dy)
        )

        if den1 == 0 or den2 == 0:
            return float("nan")

        return (
            sum(
                x*y
                for x,y in zip(dx,dy)
            )
            / (den1*den2)
        )

    print(
        f"corr(gap, residue states)  = "
        f"{correlation(
            gaps,
            [r.residue_states for r in results]
        ):.6f}"
    )

    print(
        f"corr(gap, r3 survivors)    = "
        f"{correlation(
            gaps,
            [r.r3_survivors for r in results]
        ):.6f}"
    )

    # ============================================================================================
    # STRONGEST COLLAPSES
    # ============================================================================================

    print()
    print("=" * 100)
    print("STRONGEST RESIDUE-LIFT COLLAPSES")
    print("=" * 100)

    strongest = sorted(
        results,
        key=lambda r: (
            r.r3_survivors,
            r.residue_states,
            r.product_cells,
        ),
    )[:20]

    for r in strongest:

        print(
            f"n={r.n:,} "
            f"p={r.p:,} "
            f"q={r.q:,} "
            f"mods=({r.r1},{r.r2},{r.r3}) "
            f"R/n={r.R/r.n:.10f} "
            f"gap={r.gap:,} "
            f"strip={r.strip_tests:,} "
            f"cells={r.product_cells:,} "
            f"residue={r.residue_states:,} "
            f"r2={r.r2_survivors:,} "
            f"r3={r.r3_survivors:,} "
            f"exact={r.exact_solutions}"
        )

    # ============================================================================================
    # WEAKEST COLLAPSES
    # ============================================================================================

    print()
    print("=" * 100)
    print("WEAKEST RESIDUE-LIFT COLLAPSES")
    print("=" * 100)

    weakest = sorted(
        results,
        key=lambda r: (
            r.r3_survivors,
            r.residue_states,
            r.product_cells,
        ),
        reverse=True,
    )[:20]

    for r in weakest:

        print(
            f"n={r.n:,} "
            f"p={r.p:,} "
            f"q={r.q:,} "
            f"mods=({r.r1},{r.r2},{r.r3}) "
            f"R/n={r.R/r.n:.10f} "
            f"gap={r.gap:,} "
            f"strip={r.strip_tests:,} "
            f"cells={r.product_cells:,} "
            f"residue={r.residue_states:,} "
            f"r2={r.r2_survivors:,} "
            f"r3={r.r3_survivors:,} "
            f"exact={r.exact_solutions}"
        )

    # ============================================================================================
    # ANCHOR TABLE
    # ============================================================================================

    print()
    print("=" * 100)
    print("ANCHOR RESULTS")
    print("=" * 100)

    print(
        f"{'ID':>3} "
        f"{'p':>7} "
        f"{'q':>7} "
        f"{'r1':>5} "
        f"{'r2':>5} "
        f"{'r3':>5} "
        f"{'R/n':>10} "
        f"{'STRIP':>7} "
        f"{'CELL':>7} "
        f"{'RES':>7} "
        f"{'R3':>5} "
        f"{'EXACT':>5}"
    )

    print("-" * 100)

    for r in results:

        print(
            f"{r.anchor_id:3d} "
            f"{r.p:7,d} "
            f"{r.q:7,d} "
            f"{r.r1:5d} "
            f"{r.r2:5d} "
            f"{r.r3:5d} "
            f"{r.R/r.n:10.7f} "
            f"{r.strip_tests:7,d} "
            f"{r.product_cells:7,d} "
            f"{r.residue_states:7,d} "
            f"{r.r3_survivors:5,d} "
            f"{r.exact_solutions:5,d}"
        )

    # ============================================================================================
    # TIMING
    # ============================================================================================

    total_time = (
        time.perf_counter()
        - total_start
    )

    print()
    print("=" * 100)
    print("TIMING")
    print("=" * 100)

    print(
        f"search runtime              = "
        f"{search_time:.3f} seconds"
    )

    print(
        f"total runtime               = "
        f"{total_time:.3f} seconds"
    )

    # ============================================================================================
    # INTERPRETATION
    # ============================================================================================

    print()
    print("=" * 100)
    print("MATHEMATICAL INTERPRETATION")
    print("=" * 100)

    print(
        """
The previous hyperbola experiment established a useful geometric
fact:

    fixed k
        ->
    narrow l interval.

But that experiment eventually scanned every integer p inside each
surviving cell.

This experiment removes that final integer scan.

For:

    p = k*r1 + a

we have:

    p mod r1 = a.

Likewise:

    q = l*r2 + c.

The first modular equation is:

    p*q == n (mod r1)

so:

    a*c == n (mod r1).

For invertible a:

    c == n*a^(-1) (mod r1).

The q-cell determines:

    q = l*r2 + c.

Therefore each residue a generates at most one q residue c
compatible with the first modulus.

The search becomes:

    k
     |
     v
    hyperbola l interval
     |
     v
    product cell
     |
     v
    modular residue lift
     |
     +--> second modulus
     |
     +--> third modulus
     |
     v
    exact p*q == n

Crucially, no scan over all integers inside the cell is performed.

The decisive statistics are therefore:

    residue states
    r3 survivors

compared against:

    8,363 factor primes.

A strong result would be:

    residue states << 8,363
    r3 survivors = O(1)
    exact recovery = 100%

across independently selected modulus triples.

If this happens, the residue-lift mechanism deserves further
algebraic investigation.

If residue states remain thousands or tens of thousands per anchor,
then the hyperbola reduction is still mainly geometric.

The experiment also records:

    gap = n - r1*r2*r3

and its correlations with the candidate counts.

A meaningful correlation would justify returning to the gap/resultant
direction. A near-zero correlation would indicate that R/n is mostly
controlling the scale of the quotient geometry rather than revealing
additional structure.

Most importantly, a successful exact recovery is not sufficient by
itself. The candidate-generation cost is what determines whether
this is actually reducing factorization work.
"""
    )

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()
