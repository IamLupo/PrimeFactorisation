#!/usr/bin/env python3

"""
====================================================================================================
THREE-CLOSE-PRIME QUOTIENT-CELL DEFECT / ONE-REMAINDER ELIMINATION EXPERIMENT
====================================================================================================

FOLLOW-UP TO:
    THREE-CLOSE-PRIME HYPERBOLA / RESIDUE-LIFT EXPERIMENT

Previous result:

    ~110 hyperbola cells / anchor
    ~118,000 residue states / anchor
    ~1.7 r3 survivors / anchor
    97% recovery

The bottleneck was:

    cell
      |
      v
    enumerate a
      |
      v
    residue states

This experiment removes the a-enumeration.

For:

    p = k*r1 + a
    q = l*r2 + c

we have:

    n = (k*r1+a)(l*r2+c)

Define the exact quotient-cell defect:

    D = n - k*l*r1*r2

Then:

    D = k*r1*c + l*r2*a + a*c

and therefore:

    D - k*r1*c = a*(l*r2+c)

so:

    a = (D-k*r1*c)/(l*r2+c).

Thus, for each quotient cell (k,l), we enumerate c rather than
enumerating all a.

Because r1 and r2 are both around the same scale, we explicitly
measure:

    c states tested
    exact divisibility hits for a
    r3 survivors
    exact solutions

The experiment is NOT claimed to be faster merely because it uses
one remainder coordinate.

The important question is:

    average remainder states << 118,000 ?

and preferably:

    average remainder states << 8,363 ?

This gives a clean test of whether the quotient-cell defect itself
contains enough arithmetic structure to replace the large residue
scan.

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
# CONFIG
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
# DATA
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

    full_kl: int
    strip_cells: int
    product_cells: int

    c_states: int
    a_division_hits: int
    r2_hits: int
    r3_hits: int

    exact_solutions: int
    recovered: bool


# ================================================================================================
# SIEVE
# ================================================================================================

def sieve(limit: int) -> list[int]:

    if limit < 2:
        return []

    flags = bytearray(b"\x01") * (limit + 1)

    flags[0:2] = b"\x00\x00"

    root = math.isqrt(limit)

    for p in range(2, root + 1):

        if not flags[p]:
            continue

        start = p*p
        count = ((limit-start)//p) + 1

        flags[
            start:
            limit+1:
            p
        ] = b"\x00" * count

    return [
        i
        for i, flag in enumerate(flags)
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

    result = []

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
# CLOSE MODULUS TRIPLE
# ================================================================================================

def choose_triple(
    n: int,
    modulus_primes: list[int],
) -> tuple[int, int, int]:

    target = max(
        MODULUS_MIN,
        int(n ** (1.0/3.0)),
    )

    radius = max(
        20,
        int(target*CLOSE_RATIO),
    )

    candidates = [
        r
        for r in modulus_primes
        if abs(r-target) <= radius
    ]

    if len(candidates) < 3:

        raise RuntimeError(
            f"not enough close primes for n={n}"
        )

    best = None
    best_score = None

    for r1, r2, r3 in combinations(
        candidates,
        3,
    ):

        R = r1*r2*r3

        if R >= n:
            continue

        # Maximize R/n and secondarily minimize spread.
        score = (
            R/n,
            -(r3-r1),
        )

        if (
            best_score is None
            or score > best_score
        ):

            best_score = score
            best = (r1,r2,r3)

    if best is None:

        raise RuntimeError(
            f"could not find R<n for n={n}"
        )

    return best


# ================================================================================================
# HYPERBOLA QUOTIENT INTERVAL
# ================================================================================================

def quotient_range(
    n: int,
    r1: int,
    r2: int,
    k: int,
) -> tuple[int,int] | None:

    p_lo = max(
        FACTOR_MIN,
        k*r1,
    )

    p_hi = min(
        FACTOR_MAX,
        (k+1)*r1-1,
    )

    if p_lo > p_hi:
        return None

    q_lo = n//p_hi
    q_hi = n//p_lo

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

    l_lo = q_lo//r2
    l_hi = q_hi//r2

    return l_lo,l_hi


# ================================================================================================
# CELL BOUNDS
# ================================================================================================

def cell_bounds(
    k: int,
    l: int,
    r1: int,
    r2: int,
) -> tuple[int,int,int,int] | None:

    p_lo = max(
        FACTOR_MIN,
        k*r1,
    )

    p_hi = min(
        FACTOR_MAX,
        (k+1)*r1-1,
    )

    q_lo = max(
        FACTOR_MIN,
        l*r2,
    )

    q_hi = min(
        FACTOR_MAX,
        (l+1)*r2-1,
    )

    if p_lo > p_hi or q_lo > q_hi:
        return None

    return (
        p_lo,
        p_hi,
        q_lo,
        q_hi,
    )


def cell_contains_n(
    bounds: tuple[int,int,int,int],
    n: int,
) -> bool:

    p_lo,p_hi,q_lo,q_hi = bounds

    return (
        p_lo*q_lo <= n <= p_hi*q_hi
    )


# ================================================================================================
# DEFECT SEARCH
# ================================================================================================

def search_cell(
    n: int,
    r1: int,
    r2: int,
    r3: int,
    k: int,
    l: int,
    bounds: tuple[int,int,int,int],
) -> tuple[int,int,int,int,set[tuple[int,int]]]:

    """
    Returns:

        c states
        exact integer-division hits for a
        r2-compatible states
        r3-compatible states
        exact solutions
    """

    p_lo,p_hi,q_lo,q_hi = bounds

    # Local remainder domains.
    a_lo = p_lo-k*r1
    a_hi = p_hi-k*r1

    c_lo = q_lo-l*r2
    c_hi = q_hi-l*r2

    a_lo = max(0,a_lo)
    a_hi = min(r1-1,a_hi)

    c_lo = max(0,c_lo)
    c_hi = min(r2-1,c_hi)

    if (
        a_lo > a_hi
        or c_lo > c_hi
    ):

        return (
            0,
            0,
            0,
            0,
            set(),
        )

    D = n-k*l*r1*r2

    c_states = 0
    division_hits = 0
    r2_hits = 0
    r3_hits = 0

    exact = set()

    # ============================================================================================
    # ENUMERATE THE SMALLER REMAINDER DOMAIN
    # ============================================================================================

    if (c_hi-c_lo) <= (a_hi-a_lo):

        # --------------------------------------------------------------------
        # PRIMARY FORM:
        #
        # a = (D-k*r1*c)/(l*r2+c)
        # --------------------------------------------------------------------

        for c in range(
            c_lo,
            c_hi+1,
        ):

            c_states += 1

            denominator = l*r2+c

            if denominator <= 0:
                continue

            numerator = D-k*r1*c

            # Exact divisibility is necessary.
            if numerator % denominator != 0:
                continue

            division_hits += 1

            a = numerator//denominator

            if not (
                a_lo <= a <= a_hi
            ):
                continue

            p = k*r1+a
            q = l*r2+c

            if not (
                p_lo <= p <= p_hi
                and
                q_lo <= q <= q_hi
            ):
                continue

            r2_hits += 1

            if (
                (p*q-n) % r2
                != 0
            ):
                continue

            if (
                (p*q-n) % r3
                != 0
            ):
                continue

            r3_hits += 1

            if p*q == n:

                exact.add(
                    (
                        min(p,q),
                        max(p,q),
                    )
                )

    else:

        # --------------------------------------------------------------------
        # SYMMETRIC FORM:
        #
        # c = (D-l*r2*a)/(k*r1+a)
        #
        # This is useful when the r1 remainder domain is smaller.
        # --------------------------------------------------------------------

        for a in range(
            a_lo,
            a_hi+1,
        ):

            c_states += 1

            denominator = k*r1+a

            if denominator <= 0:
                continue

            numerator = D-l*r2*a

            if numerator % denominator != 0:
                continue

            division_hits += 1

            c = numerator//denominator

            if not (
                c_lo <= c <= c_hi
            ):
                continue

            p = k*r1+a
            q = l*r2+c

            if not (
                p_lo <= p <= p_hi
                and
                q_lo <= q <= q_hi
            ):
                continue

            r2_hits += 1

            if (
                (p*q-n) % r2
                != 0
            ):
                continue

            if (
                (p*q-n) % r3
                != 0
            ):
                continue

            r3_hits += 1

            if p*q == n:

                exact.add(
                    (
                        min(p,q),
                        max(p,q),
                    )
                )

    return (
        c_states,
        division_hits,
        r2_hits,
        r3_hits,
        exact,
    )


# ================================================================================================
# ANCHOR
# ================================================================================================

def run_anchor(
    anchor_id: int,
    anchor: Anchor,
    mods: tuple[int,int,int],
) -> Result:

    n = anchor.n

    r1,r2,r3 = mods

    R = r1*r2*r3
    gap = n-R

    if R >= n:
        raise RuntimeError(
            "R must be smaller than n"
        )

    # k domain.
    k_min = FACTOR_MIN//r1

    if k_min*r1 < FACTOR_MIN:
        k_min += 1

    k_max = FACTOR_MAX//r1

    if k_max*r1 > FACTOR_MAX:
        k_max -= 1

    k_values = range(
        max(1,k_min),
        k_max+1,
    )

    l_min = FACTOR_MIN//r2
    l_max = FACTOR_MAX//r2

    full_kl = (
        (k_max-k_min+1)
        *
        (l_max-l_min+1)
    )

    strip_cells = 0
    product_cells = 0

    c_states = 0
    division_hits = 0
    r2_hits = 0
    r3_hits = 0

    exact = set()

    for k in k_values:

        qr = quotient_range(
            n,
            r1,
            r2,
            k,
        )

        if qr is None:
            continue

        l_lo,l_hi = qr

        for l in range(
            l_lo,
            l_hi+1,
        ):

            strip_cells += 1

            bounds = cell_bounds(
                k,
                l,
                r1,
                r2,
            )

            if bounds is None:
                continue

            if not cell_contains_n(
                bounds,
                n,
            ):
                continue

            product_cells += 1

            (
                local_states,
                local_division,
                local_r2,
                local_r3,
                local_exact,
            ) = search_cell(
                n,
                r1,
                r2,
                r3,
                k,
                l,
                bounds,
            )

            c_states += local_states
            division_hits += local_division
            r2_hits += local_r2
            r3_hits += local_r3

            exact.update(
                local_exact
            )

    recovered = (
        (anchor.p,anchor.q)
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
        full_kl=full_kl,
        strip_cells=strip_cells,
        product_cells=product_cells,
        c_states=c_states,
        a_division_hits=division_hits,
        r2_hits=r2_hits,
        r3_hits=r3_hits,
        exact_solutions=len(exact),
        recovered=recovered,
    )


# ================================================================================================
# CORRELATION
# ================================================================================================

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

    den_x = math.sqrt(
        sum(
            x*x
            for x in dx
        )
    )

    den_y = math.sqrt(
        sum(
            y*y
            for y in dy
        )
    )

    if den_x == 0 or den_y == 0:
        return float("nan")

    return (
        sum(
            x*y
            for x,y in zip(dx,dy)
        )
        /
        (den_x*den_y)
    )


# ================================================================================================
# MAIN
# ================================================================================================

def run() -> None:

    total_start = time.perf_counter()

    print("="*100)
    print(
        "THREE-CLOSE-PRIME QUOTIENT-CELL DEFECT / "
        "ONE-REMAINDER ELIMINATION EXPERIMENT"
    )
    print("="*100)

    print(f"M                         = {M:,}")
    print(
        f"factor range              = "
        f"{FACTOR_MIN:,} - {FACTOR_MAX:,}"
    )
    print(f"anchors                   = {ANCHORS}")
    print(
        f"modulus prime range       = "
        f"{MODULUS_MIN:,} - {MODULUS_MAX:,}"
    )
    print(
        f"close ratio               = "
        f"{CLOSE_RATIO:.0%}"
    )
    print(f"seed                      = {SEED:,}")

    # ============================================================================================
    # PRIME POOLS
    # ============================================================================================

    print()
    print("="*100)
    print("BUILDING PRIME POOLS")
    print("="*100)

    all_primes = sieve(
        FACTOR_MAX
    )

    factor_primes = [
        p
        for p in all_primes
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
    print("="*100)
    print("BUILDING ANCHORS")
    print("="*100)

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
    print("="*100)
    print(
        "SELECTING ANCHOR-SPECIFIC CLOSE TRIPLES"
    )
    print("="*100)

    triples = []

    for i,anchor in enumerate(
        anchors,
        1,
    ):

        triple = choose_triple(
            anchor.n,
            modulus_primes,
        )

        triples.append(
            triple
        )

        if (
            i % PROGRESS_STEP == 0
            or i == len(anchors)
        ):

            print(
                f"anchor {i:3d}/{len(anchors)}"
            )

    print()
    print("="*100)
    print("VALIDATING MODULUS TRIPLES")
    print("="*100)

    failures = 0

    for anchor,triple in zip(
        anchors,
        triples,
    ):

        r1,r2,r3 = triple

        if not (
            r1 < r2 < r3
        ):
            failures += 1

        R = r1*r2*r3

        if R >= anchor.n:
            failures += 1

        if (
            math.gcd(r1,r2) != 1
            or
            math.gcd(r1,r3) != 1
            or
            math.gcd(r2,r3) != 1
        ):
            failures += 1

    print(
        f"triple validation failures = "
        f"{failures}"
    )

    if failures:
        raise RuntimeError(
            "Invalid modulus triple."
        )

    print(
        f"unique modulus triples      = "
        f"{len(set(triples))}"
    )

    # ============================================================================================
    # SEARCH
    # ============================================================================================

    print()
    print("="*100)
    print(
        "RUNNING QUOTIENT-CELL DEFECT SEARCH"
    )
    print("="*100)

    results = []

    search_start = time.perf_counter()

    for i,(anchor,triple) in enumerate(
        zip(anchors,triples),
        1,
    ):

        result = run_anchor(
            i,
            anchor,
            triple,
        )

        results.append(result)

        if (
            i % PROGRESS_STEP == 0
            or i == len(anchors)
        ):

            print(
                f"anchor {i:3d}/{len(anchors)}"
            )

    search_time = (
        time.perf_counter()
        -
        search_start
    )

    # ============================================================================================
    # AVERAGES
    # ============================================================================================

    mean = lambda xs: (
        statistics.mean(xs)
        if xs
        else 0
    )

    avg_full = mean([
        r.full_kl
        for r in results
    ])

    avg_strip = mean([
        r.strip_cells
        for r in results
    ])

    avg_cells = mean([
        r.product_cells
        for r in results
    ])

    avg_states = mean([
        r.c_states
        for r in results
    ])

    avg_division = mean([
        r.a_division_hits
        for r in results
    ])

    avg_r2 = mean([
        r.r2_hits
        for r in results
    ])

    avg_r3 = mean([
        r.r3_hits
        for r in results
    ])

    avg_exact = mean([
        r.exact_solutions
        for r in results
    ])

    # ============================================================================================
    # SUMMARY
    # ============================================================================================

    print()
    print("="*100)
    print("SUMMARY")
    print("="*100)

    print(
        f"anchors analyzed            = "
        f"{len(results)}"
    )

    print(
        f"average full (k,l) pairs    = "
        f"{avg_full:,.3f}"
    )

    print(
        f"average hyperbola cells     = "
        f"{avg_strip:,.3f}"
    )

    print(
        f"average product cells       = "
        f"{avg_cells:,.3f}"
    )

    print(
        f"average remainder states    = "
        f"{avg_states:,.3f}"
    )

    print(
        f"average exact division hits = "
        f"{avg_division:,.3f}"
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
    # REDUCTION
    # ============================================================================================

    print()
    print("="*100)
    print("SEARCH REDUCTION")
    print("="*100)

    print(
        f"cells / full quotient       = "
        f"{avg_cells/avg_full:.9f}"
    )

    print(
        f"remainder / product cells   = "
        f"{avg_states/avg_cells:.3f}"
    )

    print(
        f"r3 / remainder states       = "
        f"{avg_r3/avg_states:.9f}"
    )

    print(
        f"remainder / 8,363 baseline  = "
        f"{avg_states/8363:.9f}"
    )

    print(
        f"r3 / 8,363 baseline         = "
        f"{avg_r3/8363:.9f}"
    )

    # ============================================================================================
    # RECOVERY
    # ============================================================================================

    recovered = sum(
        r.recovered
        for r in results
    )

    print()
    print("="*100)
    print("RECOVERY")
    print("="*100)

    print(
        f"correctly recovered          = "
        f"{recovered}/{len(results)}"
    )

    print(
        f"recovery rate                = "
        f"{100*recovered/len(results):.4f}%"
    )

    # ============================================================================================
    # GAP
    # ============================================================================================

    gaps = [
        r.gap
        for r in results
    ]

    print()
    print("="*100)
    print("GAP CORRELATION")
    print("="*100)

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

    print(
        f"corr(gap,remainder states) = "
        f"{correlation(
            gaps,
            [r.c_states for r in results]
        ): .6f}"
    )

    print(
        f"corr(gap,r3 survivors)      = "
        f"{correlation(
            gaps,
            [r.r3_hits for r in results]
        ): .6f}"
    )

    # ============================================================================================
    # BEST
    # ============================================================================================

    print()
    print("="*100)
    print("STRONGEST DEFECT COLLAPSES")
    print("="*100)

    strongest = sorted(
        results,
        key=lambda r: (
            r.a_division_hits,
            r.c_states,
            r.r3_hits,
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
            f"cells={r.product_cells:,} "
            f"states={r.c_states:,} "
            f"div={r.a_division_hits:,} "
            f"r2={r.r2_hits:,} "
            f"r3={r.r3_hits:,} "
            f"exact={r.exact_solutions}"
        )

    # ============================================================================================
    # WORST
    # ============================================================================================

    print()
    print("="*100)
    print("WEAKEST DEFECT COLLAPSES")
    print("="*100)

    weakest = sorted(
        results,
        key=lambda r: (
            r.a_division_hits,
            r.c_states,
            r.r3_hits,
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
            f"cells={r.product_cells:,} "
            f"states={r.c_states:,} "
            f"div={r.a_division_hits:,} "
            f"r2={r.r2_hits:,} "
            f"r3={r.r3_hits:,} "
            f"exact={r.exact_solutions}"
        )

    # ============================================================================================
    # ANCHOR TABLE
    # ============================================================================================

    print()
    print("="*100)
    print("ANCHOR RESULTS")
    print("="*100)

    print(
        f"{'ID':>3} "
        f"{'p':>7} "
        f"{'q':>7} "
        f"{'r1':>5} "
        f"{'r2':>5} "
        f"{'r3':>5} "
        f"{'R/n':>10} "
        f"{'CELL':>6} "
        f"{'STATES':>9} "
        f"{'DIV':>6} "
        f"{'R3':>5} "
        f"{'EXACT':>5}"
    )

    print("-"*100)

    for r in results:

        print(
            f"{r.anchor_id:3d} "
            f"{r.p:7,d} "
            f"{r.q:7,d} "
            f"{r.r1:5d} "
            f"{r.r2:5d} "
            f"{r.r3:5d} "
            f"{r.R/r.n:10.7f} "
            f"{r.product_cells:6,d} "
            f"{r.c_states:9,d} "
            f"{r.a_division_hits:6,d} "
            f"{r.r3_hits:5,d} "
            f"{r.exact_solutions:5,d}"
        )

    # ============================================================================================
    # TIMING
    # ============================================================================================

    total_time = (
        time.perf_counter()
        -
        total_start
    )

    print()
    print("="*100)
    print("TIMING")
    print("="*100)

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
    print("="*100)
    print("MATHEMATICAL INTERPRETATION")
    print("="*100)

    print(
        """
The previous experiment produced:

    ~110 quotient cells
        ->
    ~118,000 residue states
        ->
    ~1.7 third-modulus survivors

The problem was the residue expansion.

This experiment starts from the exact cell identity:

    n = (k*r1+a)(l*r2+c)

and defines:

    D = n-k*l*r1*r2.

Expanding:

    D = k*r1*c + l*r2*a + a*c.

Solving for a:

    a = (D-k*r1*c)/(l*r2+c).

Therefore a is no longer enumerated.

For every remainder c there is at most one exact integer a.

The search becomes:

    quotient cell (k,l)
             |
             v
        cell defect D
             |
             v
       enumerate c
             |
             v
       exact division
             |
             v
           a
             |
             v
        r2 consistency
             |
             v
        r3 consistency
             |
             v
        exact p*q == n

The critical comparison is now:

    average remainder states
        versus
    118,013 previous residue states
        versus
    8,363 factor primes.

A strong outcome would show that the defect equation has very few
integer divisibility hits inside each quotient cell.

An especially interesting result would be:

    remainder states < 8,363
    division hits << remainder states
    r3 survivors = O(1)
    100% recovery

That would mean the quotient-cell defect contains useful arithmetic
structure beyond the geometric hyperbola strip.

If the number of remainder states is still around 100,000, then the
cell equation is not buying anything.

If the remainder-state count is only a few thousand but exact
division hits are tiny, the next experiment should eliminate the
remainder variable as well using the divisibility condition itself.

Every final candidate is verified with:

    p*q == n

so no modular false positive is counted as a successful factor.

====================================================================================================
EXPERIMENT COMPLETE
====================================================================================================
"""
    )


if __name__ == "__main__":
    run()
