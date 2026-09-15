#!/usr/bin/env python3

"""
====================================================================================================
THREE-CLOSE-PRIME 1D QUOTIENT-STRIP / HYPERBOLA ELIMINATION EXPERIMENT
====================================================================================================

Goal
----
Previous experiments still enumerated a 2D quotient lattice:

    (k,l)

This experiment attempts to collapse that lattice to one dimension.

Use:

    p = a + k*r1
    q = b + l*r2

with:

    0 <= a < r1
    0 <= b < r2

and:

    p*q = n.

For a fixed quotient k, p lies in:

    [k*r1, (k+1)*r1)

Therefore:

    q = n/p

lies inside the corresponding hyperbola interval:

    n/((k+1)*r1) < q <= n/(k*r1).

That interval determines a bounded range of possible l values.

The experiment therefore performs:

    k
    |
    v
    derive q interval
    |
    v
    derive l interval
    |
    v
    test only quotient cells touched by the hyperbola
    |
    v
    product interval
    |
    v
    exact modular residue relation
    |
    v
    exact p*q == n

The important comparison is against:

    full quotient-pair enumeration

rather than only the prime baseline.

Every reported solution is independently verified by exact multiplication.

====================================================================================================
"""

from __future__ import annotations

import math
import random
import statistics
import time
from dataclasses import dataclass
from itertools import combinations
from typing import Iterable


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

# Only these anchors are used in the expensive search section.
SEARCH_ANCHORS = 300

PROGRESS_STEP = 25


# ================================================================================================
# DATA STRUCTURES
# ================================================================================================

@dataclass(frozen=True)
class Anchor:
    p: int
    q: int
    n: int
    mods: tuple[int, int, int]


@dataclass
class AnchorResult:
    anchor_id: int
    p: int
    q: int
    n: int

    r1: int
    r2: int
    r3: int

    R: int
    gap: int

    k_range: int
    full_kl_pairs: int

    strip_tests: int
    product_cells: int
    residue_survivors: int
    r3_survivors: int

    prime_candidates: int
    exact_solutions: int

    recovered: bool


# ================================================================================================
# PRIME GENERATION
# ================================================================================================

def sieve(limit: int) -> list[int]:
    """
    Return all primes <= limit.

    Uses bytearray slicing, with the slice length computed exactly so
    Python cannot raise the extended-slice assignment error.
    """
    if limit < 2:
        return []

    a = bytearray(b"\x01") * (limit + 1)
    a[0:2] = b"\x00\x00"

    root = math.isqrt(limit)

    for p in range(2, root + 1):
        if a[p]:
            start = p * p
            count = ((limit - start) // p) + 1
            a[start : limit + 1 : p] = b"\x00" * count

    return [i for i, is_prime in enumerate(a) if is_prime]


# ================================================================================================
# CLOSE-TRIPLE SELECTION
# ================================================================================================

def choose_close_triple(
    center: int,
    modulus_primes: list[int],
    close_ratio: float,
) -> tuple[int, int, int] | None:

    if len(modulus_primes) < 3:
        return None

    max_distance = max(1, int(center * close_ratio))

    candidates = [
        p
        for p in modulus_primes
        if abs(p - center) <= max_distance
    ]

    if len(candidates) < 3:
        return None

    best = None
    best_product = -1

    for triple in combinations(candidates, 3):
        r1, r2, r3 = sorted(triple)

        if r1 * r2 * r3 > best_product:
            best_product = r1 * r2 * r3
            best = (r1, r2, r3)

    return best


def select_triples(
    anchors: list[Anchor],
    modulus_primes: list[int],
    close_ratio: float,
) -> list[tuple[int, int, int]]:

    triples: list[tuple[int, int, int]] = []

    for idx, anchor in enumerate(anchors, 1):

        # Search around the geometric scale sqrt-ish of the factors.
        center = math.isqrt(anchor.n)

        triple = choose_close_triple(
            center,
            modulus_primes,
            close_ratio,
        )

        if triple is None:
            # Fall back to the largest three modulus primes.
            triple = tuple(sorted(modulus_primes[-3:]))

        triples.append(triple)

        if idx % PROGRESS_STEP == 0 or idx == len(anchors):
            print(f"anchor {idx:3d}/{len(anchors)}", flush=True)

    return triples


# ================================================================================================
# INTEGER INTERVAL HELPERS
# ================================================================================================

def ceil_div(a: int, b: int) -> int:
    return -((-a) // b)


def floor_div(a: int, b: int) -> int:
    return a // b


def quotient_range_for_k(
    n: int,
    r1: int,
    r2: int,
    k: int,
    factor_min: int,
    factor_max: int,
) -> tuple[int, int] | None:

    p_lo = max(
        factor_min,
        k * r1,
    )

    p_hi = min(
        factor_max,
        (k + 1) * r1 - 1,
    )

    if p_lo > p_hi:
        return None

    # Since q = n/p and p>0:
    #
    # minimum q is floor(n / p_hi)
    # maximum q is floor(n / p_lo)
    #
    q_lo = n // p_hi
    q_hi = n // p_lo

    if q_lo > q_hi:
        return None

    l_lo = floor_div(q_lo, r2)
    l_hi = floor_div(q_hi, r2)

    return l_lo, l_hi


def cell_bounds(
    k: int,
    l: int,
    r1: int,
    r2: int,
    factor_min: int,
    factor_max: int,
) -> tuple[int, int, int, int] | None:

    p_lo = max(factor_min, k * r1)
    p_hi = min(factor_max, (k + 1) * r1 - 1)

    q_lo = max(factor_min, l * r2)
    q_hi = min(factor_max, (l + 1) * r2 - 1)

    if p_lo > p_hi or q_lo > q_hi:
        return None

    return p_lo, p_hi, q_lo, q_hi


def product_interval_contains_n(
    p_lo: int,
    p_hi: int,
    q_lo: int,
    q_hi: int,
    n: int,
) -> bool:

    # Everything is positive in this experiment.
    return p_lo * q_lo <= n <= p_hi * q_hi


# ================================================================================================
# MODULAR CHECKS
# ================================================================================================

def modular_cell_relation(
    n: int,
    p: int,
    q: int,
    r: int,
) -> bool:
    return (p * q - n) % r == 0


def residue_pair_search(
    n: int,
    r1: int,
    r2: int,
    k: int,
    l: int,
    p_lo: int,
    p_hi: int,
    q_lo: int,
    q_hi: int,
) -> tuple[int, list[tuple[int, int]]]:

    """
    For a surviving quotient cell, inspect possible p residues modulo r1.

    p = k*r1 + a
    q = l*r2 + b

    Since r1 is the first scale:

        p mod r1 = a.

    For each possible a in the cell, q is constrained by:

        q = n/p.

    We avoid a full factor-prime enumeration and only retain exact
    divisibility candidates.

    This function returns:

        number of integer residue attempts
        exact (p,q) candidates found
    """

    attempts = 0
    hits: list[tuple[int, int]] = []

    # The p-cell width is at most r1, so this is deliberately bounded
    # by the modulus scale rather than the complete factor interval.
    for p in range(p_lo, p_hi + 1):
        attempts += 1

        if n % p:
            continue

        q = n // p

        if not (q_lo <= q <= q_hi):
            continue

        # Exact modular verification.
        if (p * q) % r1 != n % r1:
            continue

        if (p * q) % r2 != n % r2:
            continue

        hits.append((p, q))

    return attempts, hits


# ================================================================================================
# TRUE FACTOR GENERATION
# ================================================================================================

def build_anchors(
    factor_primes: list[int],
    count: int,
    seed: int,
) -> list[Anchor]:

    rng = random.Random(seed)

    usable = factor_primes

    anchors: list[Anchor] = []

    seen: set[tuple[int, int]] = set()

    while len(anchors) < count:

        p = rng.choice(usable)
        q = rng.choice(usable)

        if p == q:
            continue

        lo, hi = sorted((p, q))

        if not (FACTOR_MIN <= lo <= FACTOR_MAX):
            continue

        if not (FACTOR_MIN <= hi <= FACTOR_MAX):
            continue

        if (lo, hi) in seen:
            continue

        seen.add((lo, hi))

        anchors.append(
            Anchor(
                p=lo,
                q=hi,
                n=lo * hi,
                mods=(0, 0, 0),
            )
        )

    return anchors


# ================================================================================================
# MAIN SEARCH
# ================================================================================================

def run_anchor(
    anchor_id: int,
    anchor: Anchor,
    mods: tuple[int, int, int],
) -> AnchorResult:

    p_true = anchor.p
    q_true = anchor.q
    n = anchor.n

    r1, r2, r3 = mods

    R = r1 * r2 * r3
    gap = n - R

    # --------------------------------------------------------------------------------------------
    # k domain
    # --------------------------------------------------------------------------------------------

    k_min = FACTOR_MIN // r1
    if k_min * r1 < FACTOR_MIN:
        k_min += 1

    k_max = FACTOR_MAX // r1

    if k_max * r1 > FACTOR_MAX:
        k_max -= 1

    k_values = list(range(max(0, k_min), k_max + 1))

    # Full 2D quotient-pair baseline.
    l_min_global = FACTOR_MIN // r2
    l_max_global = FACTOR_MAX // r2

    full_kl_pairs = (
        len(k_values)
        * (l_max_global - l_min_global + 1)
    )

    strip_tests = 0
    product_cells = 0
    residue_survivors = 0
    r3_survivors = 0
    prime_candidates = 0

    exact_pairs: set[tuple[int, int]] = set()

    # --------------------------------------------------------------------------------------------
    # Hyperbola strip search
    # --------------------------------------------------------------------------------------------

    for k in k_values:

        qr = quotient_range_for_k(
            n,
            r1,
            r2,
            k,
            FACTOR_MIN,
            FACTOR_MAX,
        )

        if qr is None:
            continue

        l_lo, l_hi = qr

        for l in range(l_lo, l_hi + 1):

            strip_tests += 1

            bounds = cell_bounds(
                k,
                l,
                r1,
                r2,
                FACTOR_MIN,
                FACTOR_MAX,
            )

            if bounds is None:
                continue

            p_lo, p_hi, q_lo, q_hi = bounds

            if not product_interval_contains_n(
                p_lo,
                p_hi,
                q_lo,
                q_hi,
                n,
            ):
                continue

            product_cells += 1

            # ------------------------------------------------------------------------------------
            # Exact search INSIDE ONLY the surviving hyperbola cell.
            #
            # This is intentionally not reported as a prime search:
            # we are measuring how much of the quotient lattice was
            # removed before touching the underlying integer cell.
            # ------------------------------------------------------------------------------------

            _, hits = residue_pair_search(
                n,
                r1,
                r2,
                k,
                l,
                p_lo,
                p_hi,
                q_lo,
                q_hi,
            )

            for p, q in hits:

                residue_survivors += 1

                # Independent third-modulus verification.
                if not modular_cell_relation(n, p, q, r3):
                    continue

                r3_survivors += 1

                # Prime/exact candidate counter.
                if p > 1 and q > 1:
                    prime_candidates += 1

                if p * q == n:
                    exact_pairs.add((min(p, q), max(p, q)))

    recovered = (p_true, q_true) in exact_pairs

    return AnchorResult(
        anchor_id=anchor_id,
        p=p_true,
        q=q_true,
        n=n,
        r1=r1,
        r2=r2,
        r3=r3,
        R=R,
        gap=gap,
        k_range=len(k_values),
        full_kl_pairs=full_kl_pairs,
        strip_tests=strip_tests,
        product_cells=product_cells,
        residue_survivors=residue_survivors,
        r3_survivors=r3_survivors,
        prime_candidates=prime_candidates,
        exact_solutions=len(exact_pairs),
        recovered=recovered,
    )


# ================================================================================================
# RUN
# ================================================================================================

def run() -> None:

    start_total = time.perf_counter()

    print("=" * 100)
    print("THREE-CLOSE-PRIME 1D QUOTIENT-STRIP / HYPERBOLA ELIMINATION EXPERIMENT")
    print("=" * 100)
    print(f"M                         = {M:,}")
    print(f"factor range              = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"anchors                   = {ANCHORS}")
    print(f"modulus prime range       = {MODULUS_MIN:,} - {MODULUS_MAX:,}")
    print(f"close ratio               = {CLOSE_RATIO:.0%}")
    print(f"seed                      = {SEED:,}")
    print()

    print("=" * 100)
    print("BUILDING PRIME POOLS")
    print("=" * 100)

    factor_primes_all = sieve(FACTOR_MAX)
    factor_primes = [
        p for p in factor_primes_all
        if FACTOR_MIN <= p <= FACTOR_MAX
    ]

    modulus_primes = [
        p
        for p in sieve(MODULUS_MAX)
        if MODULUS_MIN <= p <= MODULUS_MAX
    ]

    print(f"factor primes             = {len(factor_primes):,}")
    print(f"modulus primes            = {len(modulus_primes):,}")
    print()

    print("=" * 100)
    print("BUILDING ANCHORS")
    print("=" * 100)

    anchors = build_anchors(
        factor_primes,
        ANCHORS,
        SEED,
    )

    print(f"actual anchors            = {len(anchors):,}")
    print()

    print("=" * 100)
    print("SELECTING MAXIMUM-PRODUCT CLOSE TRIPLES")
    print("=" * 100)

    triples = select_triples(
        anchors,
        modulus_primes,
        CLOSE_RATIO,
    )

    print()
    print("=" * 100)
    print("VALIDATING SELECTED MODULI")
    print("=" * 100)

    failures = 0

    for idx, anchor in enumerate(anchors):

        r1, r2, r3 = triples[idx]

        if not (r1 < r2 < r3):
            failures += 1
            continue

        if not (
            math.gcd(r1, r2) == 1
            and math.gcd(r1, r3) == 1
            and math.gcd(r2, r3) == 1
        ):
            failures += 1
            continue

    print(f"triple validation failures = {failures}")

    if failures:
        raise RuntimeError(
            "Modulus triples failed validation."
        )

    print()

    # ============================================================================================
    # SEARCH
    # ============================================================================================

    print("=" * 100)
    print("RUNNING 1D HYPERBOLA-STRIP SEARCH")
    print("=" * 100)

    results: list[AnchorResult] = []

    search_count = min(
        SEARCH_ANCHORS,
        len(anchors),
    )

    t_search = time.perf_counter()

    for idx in range(search_count):

        result = run_anchor(
            idx + 1,
            anchors[idx],
            triples[idx],
        )

        results.append(result)

        if (idx + 1) % PROGRESS_STEP == 0 or idx + 1 == search_count:
            print(
                f"anchor {idx + 1:3d}/{search_count}",
                flush=True,
            )

    search_time = time.perf_counter() - t_search

    # ============================================================================================
    # SUMMARY
    # ============================================================================================

    print()
    print("=" * 100)
    print("SUMMARY")
    print("=" * 100)

    avg = lambda values: (
        statistics.mean(values) if values else 0.0
    )

    baseline_kl = avg(
        [r.full_kl_pairs for r in results]
    )

    avg_strip = avg(
        [r.strip_tests for r in results]
    )

    avg_cells = avg(
        [r.product_cells for r in results]
    )

    avg_residue = avg(
        [r.residue_survivors for r in results]
    )

    avg_r3 = avg(
        [r.r3_survivors for r in results]
    )

    avg_prime = avg(
        [r.prime_candidates for r in results]
    )

    avg_exact = avg(
        [r.exact_solutions for r in results]
    )

    print(f"anchors analyzed                = {len(results)}")
    print(f"average full (k,l) pairs        = {baseline_kl:,.3f}")
    print(f"average hyperbola strip tests   = {avg_strip:,.3f}")
    print(f"average product-cell survivors  = {avg_cells:,.3f}")
    print(f"average residue survivors       = {avg_residue:,.3f}")
    print(f"average r3 survivors             = {avg_r3:,.3f}")
    print(f"average prime candidates         = {avg_prime:,.3f}")
    print(f"average exact solutions          = {avg_exact:,.3f}")

    # ============================================================================================
    # REDUCTION
    # ============================================================================================

    print()
    print("=" * 100)
    print("SEARCH REDUCTION")
    print("=" * 100)

    if baseline_kl:
        print(
            "hyperbola-strip / full quotient lattice"
            f" = {avg_strip / baseline_kl:.9f}"
        )
        print(
            "hyperbola-strip reduction"
            f" = {(1.0 - avg_strip / baseline_kl) * 100:.6f}%"
        )

    if avg_strip:
        print(
            "product-cell / strip"
            f" = {avg_cells / avg_strip:.9f}"
        )
        print(
            "product-cell reduction from strip"
            f" = {(1.0 - avg_cells / avg_strip) * 100:.6f}%"
        )

    if avg_cells:
        print(
            "r3 / product-cell"
            f" = {avg_r3 / avg_cells:.9f}"
        )
        print(
            "r3 reduction from product-cell"
            f" = {(1.0 - avg_r3 / avg_cells) * 100:.6f}%"
        )

    # ============================================================================================
    # RECOVERY
    # ============================================================================================

    recovered = sum(
        1 for r in results if r.recovered
    )

    print()
    print("=" * 100)
    print("RECOVERY")
    print("=" * 100)
    print(
        f"correctly recovered              = "
        f"{recovered}/{len(results)}"
    )

    print(
        f"recovery rate                    = "
        f"{(recovered / len(results) * 100.0):.4f}%"
        if results
        else "recovery rate                    = 0.0000%"
    )

    # ============================================================================================
    # HYPERBOLA STRIP WIDTH
    # ============================================================================================

    print()
    print("=" * 100)
    print("HYPERBOLA STRIP WIDTH")
    print("=" * 100)

    strip_widths = [
        r.strip_tests / max(1, r.k_range)
        for r in results
    ]

    print(
        f"mean l states per k             = "
        f"{avg(strip_widths):,.6f}"
    )

    print(
        f"minimum l states per k          = "
        f"{min(strip_widths) if strip_widths else 0:.6f}"
    )

    print(
        f"maximum l states per k          = "
        f"{max(strip_widths) if strip_widths else 0:.6f}"
    )

    # ============================================================================================
    # MOST AGGRESSIVE COLLAPSES
    # ============================================================================================

    print()
    print("=" * 100)
    print("MOST AGGRESSIVE HYPERBOLA COLLAPSES")
    print("=" * 100)

    strongest = sorted(
        results,
        key=lambda r: (
            r.strip_tests / max(1, r.full_kl_pairs),
            r.product_cells,
        ),
    )[:20]

    for r in strongest:

        print(
            f"n={r.n:,} "
            f"p={r.p:,} "
            f"q={r.q:,} "
            f"mods=({r.r1},{r.r2},{r.r3}) "
            f"R/n={r.R / r.n:.10f} "
            f"fullKL={r.full_kl_pairs:,} "
            f"strip={r.strip_tests:,} "
            f"cells={r.product_cells:,} "
            f"residue={r.residue_survivors:,} "
            f"r3={r.r3_survivors:,} "
            f"exact={r.exact_solutions}"
        )

    # ============================================================================================
    # WEAKEST COLLAPSES
    # ============================================================================================

    print()
    print("=" * 100)
    print("WEAKEST HYPERBOLA COLLAPSES")
    print("=" * 100)

    weakest = sorted(
        results,
        key=lambda r: (
            r.strip_tests / max(1, r.full_kl_pairs),
            -r.product_cells,
        ),
        reverse=True,
    )[:20]

    for r in weakest:

        print(
            f"n={r.n:,} "
            f"p={r.p:,} "
            f"q={r.q:,} "
            f"mods=({r.r1},{r.r2},{r.r3}) "
            f"R/n={r.R / r.n:.10f} "
            f"fullKL={r.full_kl_pairs:,} "
            f"strip={r.strip_tests:,} "
            f"cells={r.product_cells:,} "
            f"residue={r.residue_survivors:,} "
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
        f"{'FULLKL':>8} "
        f"{'STRIP':>8} "
        f"{'CELL':>7} "
        f"{'RES':>6} "
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
            f"{r.full_kl_pairs:8,d} "
            f"{r.strip_tests:8,d} "
            f"{r.product_cells:7,d} "
            f"{r.residue_survivors:6,d} "
            f"{r.r3_survivors:5,d} "
            f"{r.exact_solutions:5,d}"
        )

    # ============================================================================================
    # TIMING
    # ============================================================================================

    total_time = time.perf_counter() - start_total

    print()
    print("=" * 100)
    print("TIMING")
    print("=" * 100)
    print(f"search runtime                  = {search_time:.3f} seconds")
    print(f"total runtime                   = {total_time:.3f} seconds")

    # ============================================================================================
    # INTERPRETATION
    # ============================================================================================

    print()
    print("=" * 100)
    print("MATHEMATICAL INTERPRETATION")
    print("=" * 100)

    print(
        """
Previous quotient experiments used a 2D lattice:

    (k,l)

with:

    p = a + k*r1
    q = b + l*r2.

The present experiment uses the hyperbola:

    p*q = n.

For a fixed k:

    k*r1 <= p < (k+1)*r1.

Since:

    q = n/p,

we obtain:

    n/((k+1)*r1) < q <= n/(k*r1).

Therefore the possible l values are obtained directly from
the q interval.

The search becomes:

    k
    |
    v
    hyperbola q-range
    |
    v
    l interval
    |
    v
    quotient cell
    |
    v
    product interval
    |
    v
    exact residue/divisibility test
    |
    v
    exact factorization.

The crucial statistic is:

    hyperbola-strip tests / full (k,l) pairs.

If this ratio is very small, the hyperbola itself is eliminating
most of the quotient lattice before any modular arithmetic.

The next important statistic is:

    product-cell survivors.

If the quotient strip is narrow AND product cells collapse
strongly, then the factor search is becoming effectively
one-dimensional.

However, there is a deliberate caveat:

    the final cell search still inspects integer p values
    inside surviving cells.

Therefore this experiment does NOT claim a new factoring
algorithm merely from quotient compression.

The decisive question is whether the number and width of the
surviving cells can themselves be driven close to O(1) per anchor
without reverting to essentially all integers in the factor
interval.

Every final candidate is checked with:

    p*q == n

so modular or quotient false positives cannot count as recovery.
"""
    )

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()
