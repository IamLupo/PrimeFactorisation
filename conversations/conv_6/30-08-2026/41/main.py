#!/usr/bin/env python3

import math
import random
import statistics
import time
from dataclasses import dataclass
from collections import Counter


# ==================================================================================================
# THREE-CLOSE-PRIME QUOTIENT-BLOCK / LINEAR-CONGRUENCE INVERSION EXPERIMENT
# ==================================================================================================

M = 111_546_435

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

MOD_MIN = 300
MOD_MAX = 3_000

CLOSE_RATIO = 1.20

ANCHORS = 300
SEED = 1_511_464_998


# ==================================================================================================
# DATA STRUCTURES
# ==================================================================================================

@dataclass(frozen=True)
class Anchor:
    p: int
    q: int

    @property
    def n(self) -> int:
        return self.p * self.q


@dataclass(frozen=True)
class ModTriple:
    r1: int
    r2: int
    r3: int

    @property
    def R12(self) -> int:
        return self.r1 * self.r2

    @property
    def R123(self) -> int:
        return self.r1 * self.r2 * self.r3


@dataclass
class Stats:
    quotient_blocks: int = 0

    r1_blocks: int = 0
    r12_blocks: int = 0
    r123_blocks: int = 0

    r1_candidates: int = 0
    r12_candidates: int = 0
    r123_candidates: int = 0

    exact_solutions: int = 0

    max_block_width: int = 0
    max_candidates_per_block: int = 0


# ==================================================================================================
# PRIME UTILITIES
# ==================================================================================================

def sieve(limit: int) -> list[int]:
    if limit < 2:
        return []

    bits = bytearray(b"\x01") * (limit + 1)
    bits[0] = 0
    bits[1] = 0

    root = math.isqrt(limit)

    for p in range(2, root + 1):
        if bits[p]:
            start = p * p
            bits[start:limit + 1:p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    return [i for i in range(2, limit + 1) if bits[i]]


def prime_pool(lo: int, hi: int) -> list[int]:
    return [p for p in sieve(hi) if p >= lo]


# ==================================================================================================
# GCD / MODULAR INVERSE
# ==================================================================================================

def egcd(a: int, b: int):
    if b == 0:
        return a, 1, 0

    g, x1, y1 = egcd(b, a % b)

    return (
        g,
        y1,
        x1 - (a // b) * y1,
    )


def invmod(a: int, m: int) -> int | None:
    if m == 1:
        return 0

    a %= m

    g, x, _ = egcd(a, m)

    if g != 1:
        return None

    return x % m


# ==================================================================================================
# ANCHOR GENERATION
# ==================================================================================================

def build_anchors(
    factor_primes: list[int],
    count: int,
    seed: int,
) -> list[Anchor]:
    """
    IMPORTANT:
    Each anchor is its own semiprime n = p*q.

    We do NOT use M as the factorization target.
    """

    rng = random.Random(seed)

    anchors = []

    # Try randomized distinct prime pairs.
    max_attempts = max(100_000, count * 200)

    for _ in range(max_attempts):

        p = rng.choice(factor_primes)
        q = rng.choice(factor_primes)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        if not (
            FACTOR_MIN <= p <= FACTOR_MAX
            and FACTOR_MIN <= q <= FACTOR_MAX
        ):
            continue

        anchors.append(
            Anchor(p, q)
        )

        if len(anchors) >= count:
            break

    # Remove accidental duplicates.
    anchors = list(
        dict.fromkeys(anchors)
    )

    if len(anchors) < count:
        # Deterministic fallback.
        all_pairs = []

        # Restrict the fallback so memory does not explode.
        # We only need enough unique pairs.
        sample = factor_primes[:]

        rng.shuffle(sample)

        for i in range(min(len(sample), 5000)):
            p = sample[i]

            for j in range(i + 1, min(len(sample), i + 100)):
                q = sample[j]

                if p == q:
                    continue

                pp, qq = sorted((p, q))

                if (
                    FACTOR_MIN <= pp <= FACTOR_MAX
                    and FACTOR_MIN <= qq <= FACTOR_MAX
                ):
                    all_pairs.append(
                        Anchor(pp, qq)
                    )

                    if len(all_pairs) >= count:
                        break

            if len(all_pairs) >= count:
                break

        anchors = list(
            dict.fromkeys(
                anchors + all_pairs
            )
        )

    if len(anchors) < count:
        raise RuntimeError(
            f"Could only construct {len(anchors)} anchors, "
            f"requested {count}."
        )

    return anchors[:count]


# ==================================================================================================
# CLOSE MODULUS TRIPLES
# ==================================================================================================

def build_close_triples(
    mods: list[int],
) -> list[ModTriple]:

    triples = []

    for i, r1 in enumerate(mods):

        for j in range(i + 1, len(mods)):

            r2 = mods[j]

            if r2 > r1 * CLOSE_RATIO:
                break

            for k in range(j + 1, len(mods)):

                r3 = mods[k]

                if r3 > r1 * CLOSE_RATIO:
                    break

                triples.append(
                    ModTriple(
                        r1,
                        r2,
                        r3,
                    )
                )

    return triples


def choose_triple(
    n: int,
    triples: list[ModTriple],
) -> ModTriple:

    best = None
    best_product = -1

    for t in triples:

        R = t.R123

        if R >= n:
            continue

        if R > best_product:
            best_product = R
            best = t

    if best is None:
        raise RuntimeError(
            f"No modulus triple with R<n for n={n:,}"
        )

    return best


# ==================================================================================================
# QUOTIENT BLOCK
# ==================================================================================================

def quotient_block(
    n: int,
    q: int,
    p_min: int,
    p_max: int,
):
    """
    Exact interval for:

        floor(n/p) = q

    is:

        floor(n/(q+1)) + 1 <= p <= floor(n/q)
    """

    lo = n // (q + 1) + 1
    hi = n // q

    if lo < p_min:
        lo = p_min

    if hi > p_max:
        hi = p_max

    if lo > hi:
        return None

    return lo, hi


# ==================================================================================================
# LINEAR CONGRUENCE
# ==================================================================================================

def solve_linear_congruence(
    a: int,
    b: int,
    modulus: int,
):
    """
    Solve:

        a*x = b (mod modulus)

    Return:

        (residue, reduced_modulus)

    where:

        x = residue (mod reduced_modulus)

    or None if unsatisfiable.
    """

    if modulus <= 0:
        raise ValueError("modulus must be positive")

    if modulus == 1:
        return 0, 1

    g = math.gcd(a, modulus)

    if b % g != 0:
        return None

    aa = a // g
    bb = b // g
    mm = modulus // g

    inv = invmod(
        aa,
        mm,
    )

    if inv is None:
        return None

    residue = (
        bb * inv
    ) % mm

    return residue, mm


def first_in_class(
    residue: int,
    modulus: int,
    lo: int,
):
    """
    Smallest x >= lo satisfying:

        x == residue (mod modulus)
    """

    if modulus == 1:
        return lo

    return lo + (
        (residue - lo) % modulus
    )


def class_count_in_interval(
    residue: int,
    modulus: int,
    lo: int,
    hi: int,
) -> int:

    if lo > hi:
        return 0

    first = first_in_class(
        residue,
        modulus,
        lo,
    )

    if first > hi:
        return 0

    return (
        (hi - first) // modulus
    ) + 1


# ==================================================================================================
# DIRECT QUOTIENT-BLOCK SOLVER
# ==================================================================================================

def search_anchor(
    anchor: Anchor,
    mods: ModTriple,
) -> Stats:

    n = anchor.n

    stats = Stats()

    # We only search the smaller side.
    p_min = FACTOR_MIN
    p_max = min(
        FACTOR_MAX,
        math.isqrt(n),
    )

    if p_min > p_max:
        return stats

    R1 = mods.r1
    R12 = mods.R12
    R123 = mods.R123

    # ----------------------------------------------------------------------------------------------
    # IMPORTANT:
    #
    # Distinct floor(n/p) values can be processed as quotient blocks.
    #
    # We use the standard quotient-block range decomposition:
    #
    #     q = floor(n/p)
    #
    # and jump directly from block to block.
    #
    # ----------------------------------------------------------------------------------------------

    p = p_min

    while p <= p_max:

        q = n // p

        block = quotient_block(
            n,
            q,
            p_min,
            p_max,
        )

        if block is None:
            p += 1
            continue

        lo, hi = block

        stats.quotient_blocks += 1

        width = hi - lo + 1

        if width > stats.max_block_width:
            stats.max_block_width = width

        # ==========================================================================================
        # STAGE 1: MOD r1
        #
        #     q*p == n (mod r1)
        # ==========================================================================================

        sol1 = solve_linear_congruence(
            q,
            n,
            R1,
        )

        count1 = 0

        if sol1 is not None:

            residue1, modulus1 = sol1

            count1 = class_count_in_interval(
                residue1,
                modulus1,
                lo,
                hi,
            )

            if count1:
                stats.r1_blocks += 1

        stats.r1_candidates += count1

        # ==========================================================================================
        # STAGE 2: MOD r1*r2
        #
        # Solve directly at the combined modulus.
        # ==========================================================================================

        sol2 = solve_linear_congruence(
            q,
            n,
            R12,
        )

        count2 = 0

        if sol2 is not None:

            residue2, modulus2 = sol2

            count2 = class_count_in_interval(
                residue2,
                modulus2,
                lo,
                hi,
            )

            if count2:
                stats.r12_blocks += 1

        stats.r12_candidates += count2

        # ==========================================================================================
        # STAGE 3: MOD r1*r2*r3
        #
        # Complete three-modulus condition.
        #
        #     q*p == n (mod R123)
        #
        # ==========================================================================================

        sol3 = solve_linear_congruence(
            q,
            n,
            R123,
        )

        count3 = 0

        if sol3 is not None:

            residue3, modulus3 = sol3

            count3 = class_count_in_interval(
                residue3,
                modulus3,
                lo,
                hi,
            )

            if count3:
                stats.r123_blocks += 1

        stats.r123_candidates += count3

        if count3 > stats.max_candidates_per_block:
            stats.max_candidates_per_block = count3

        # ==========================================================================================
        # EXACT VERIFICATION
        #
        # IMPORTANT:
        #
        # We DO NOT enumerate every p in the class.
        #
        # Normally count3 should be 0 or 1 because R123 is large.
        # If the progression contains multiple values, we materialize
        # only those actual congruence survivors.
        # ==========================================================================================

        if count3:

            residue3, modulus3 = sol3

            first = first_in_class(
                residue3,
                modulus3,
                lo,
            )

            last = first + (
                count3 - 1
            ) * modulus3

            candidate = first

            while candidate <= last:

                if (
                    candidate >= FACTOR_MIN
                    and candidate <= FACTOR_MAX
                    and candidate <= math.isqrt(n)
                ):

                    if n % candidate == 0:

                        other = n // candidate

                        if (
                            FACTOR_MIN <= other <= FACTOR_MAX
                            and candidate * other == n
                        ):
                            stats.exact_solutions += 1

                candidate += modulus3

        # ------------------------------------------------------------------------------------------
        # Jump directly to the beginning of the next quotient block.
        # ------------------------------------------------------------------------------------------

        p = hi + 1

    return stats


# ==================================================================================================
# MAIN
# ==================================================================================================

def run():

    total_start = time.perf_counter()

    print("=" * 100)
    print("THREE-CLOSE-PRIME QUOTIENT-BLOCK / LINEAR-CONGRUENCE INVERSION EXPERIMENT")
    print("=" * 100)

    print(f"M                         = {M:,}")
    print(f"factor range              = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"modulus prime range       = {MOD_MIN:,} - {MOD_MAX:,}")
    print(f"close ratio               = {CLOSE_RATIO:.0%}")
    print(f"anchors                   = {ANCHORS}")
    print(f"seed                      = {SEED}")

    # ==============================================================================================
    # PRIME POOLS
    # ==============================================================================================

    print("\n" + "=" * 100)
    print("BUILDING PRIME POOLS")
    print("=" * 100)

    factor_primes = prime_pool(
        FACTOR_MIN,
        FACTOR_MAX,
    )

    modulus_primes = prime_pool(
        MOD_MIN,
        MOD_MAX,
    )

    print(
        f"factor primes             = "
        f"{len(factor_primes):,}"
    )

    print(
        f"modulus primes            = "
        f"{len(modulus_primes):,}"
    )

    # ==============================================================================================
    # ANCHORS
    # ==============================================================================================

    print("\n" + "=" * 100)
    print("BUILDING ANCHORS")
    print("=" * 100)

    anchors = build_anchors(
        factor_primes,
        ANCHORS,
        SEED,
    )

    print(
        f"actual anchors            = "
        f"{len(anchors)}"
    )

    # ==============================================================================================
    # CLOSE TRIPLES
    # ==============================================================================================

    print("\n" + "=" * 100)
    print("SELECTING ANCHOR-SPECIFIC CLOSE TRIPLES")
    print("=" * 100)

    triples = build_close_triples(
        modulus_primes,
    )

    selected = []

    for i, anchor in enumerate(
        anchors,
        1,
    ):

        t = choose_triple(
            anchor.n,
            triples,
        )

        selected.append(t)

        if i % 25 == 0:
            print(
                f"anchor {i:3d}/{len(anchors)}"
            )

    print(
        f"unique modulus triples    = "
        f"{len(set(selected))}"
    )

    # ==============================================================================================
    # VALIDATION
    # ==============================================================================================

    print("\n" + "=" * 100)
    print("VALIDATING QUOTIENT-BLOCK ALGEBRA")
    print("=" * 100)

    failures = 0

    for anchor, mods in zip(
        anchors,
        selected,
    ):

        n = anchor.n
        p = anchor.p
        q = anchor.q

        # q must be floor(n/p)
        if n // p != q:
            failures += 1
            continue

        block = quotient_block(
            n,
            q,
            FACTOR_MIN,
            min(
                FACTOR_MAX,
                math.isqrt(n),
            ),
        )

        if block is None:
            failures += 1
            continue

        lo, hi = block

        if not (
            lo <= p <= hi
        ):
            failures += 1
            continue

        # Exact modular condition.
        R = mods.R123

        if (
            (p * q - n) % R
            != 0
        ):
            failures += 1

    print(
        f"identity failures         = "
        f"{failures}"
    )

    if failures:
        raise RuntimeError(
            "Quotient-block validation failed."
        )

    # ==============================================================================================
    # SEARCH
    # ==============================================================================================

    print("\n" + "=" * 100)
    print("RUNNING QUOTIENT-BLOCK SEARCH")
    print("=" * 100)

    results = []

    for i, (
        anchor,
        mods,
    ) in enumerate(
        zip(
            anchors,
            selected,
        ),
        1,
    ):

        start = time.perf_counter()

        stats = search_anchor(
            anchor,
            mods,
        )

        elapsed = (
            time.perf_counter()
            - start
        )

        p_states = min(
            FACTOR_MAX,
            math.isqrt(anchor.n),
        ) - FACTOR_MIN + 1

        prime_states = sum(
            1
            for p in factor_primes
            if FACTOR_MIN <= p <= math.isqrt(anchor.n)
        )

        recovered = (
            anchor.p * anchor.q
            == anchor.n
            and stats.exact_solutions > 0
        )

        results.append(
            (
                anchor,
                mods,
                stats,
                p_states,
                prime_states,
                elapsed,
                recovered,
            )
        )

        if i % 25 == 0:
            print(
                f"anchor {i:3d}/{len(anchors)}"
            )

    # ==============================================================================================
    # SUMMARY
    # ==============================================================================================

    print("\n" + "=" * 100)
    print("SUMMARY")
    print("=" * 100)

    avg_p_states = statistics.mean(
        r[3]
        for r in results
    )

    avg_prime_states = statistics.mean(
        r[4]
        for r in results
    )

    avg_blocks = statistics.mean(
        r[2].quotient_blocks
        for r in results
    )

    avg_r1 = statistics.mean(
        r[2].r1_candidates
        for r in results
    )

    avg_r12 = statistics.mean(
        r[2].r12_candidates
        for r in results
    )

    avg_r123 = statistics.mean(
        r[2].r123_candidates
        for r in results
    )

    avg_exact = statistics.mean(
        r[2].exact_solutions
        for r in results
    )

    avg_time = statistics.mean(
        r[5]
        for r in results
    )

    recovered = sum(
        1
        for r in results
        if r[6]
    )

    print(
        f"anchors analyzed              = "
        f"{len(results)}"
    )

    print(
        f"average integer p states      = "
        f"{avg_p_states:,.3f}"
    )

    print(
        f"average prime p states        = "
        f"{avg_prime_states:,.3f}"
    )

    print(
        f"average quotient blocks       = "
        f"{avg_blocks:,.3f}"
    )

    print(
        f"average r1 candidates         = "
        f"{avg_r1:,.3f}"
    )

    print(
        f"average r1*r2 candidates      = "
        f"{avg_r12:,.3f}"
    )

    print(
        f"average r1*r2*r3 candidates   = "
        f"{avg_r123:,.3f}"
    )

    print(
        f"average exact solutions       = "
        f"{avg_exact:,.3f}"
    )

    print(
        f"average search time            = "
        f"{avg_time * 1000:.3f} ms"
    )

    # ==============================================================================================
    # REDUCTION
    # ==============================================================================================

    print("\n" + "=" * 100)
    print("SEARCH REDUCTION")
    print("=" * 100)

    print(
        f"quotient blocks / integer p   = "
        f"{avg_blocks / avg_p_states:.9f}"
    )

    print(
        f"quotient-block reduction      = "
        f"{100.0 * (1.0 - avg_blocks / avg_p_states):.6f}%"
    )

    print(
        f"r1 / integer p                = "
        f"{avg_r1 / avg_p_states:.9f}"
    )

    print(
        f"r1*r2 / integer p             = "
        f"{avg_r12 / avg_p_states:.9f}"
    )

    print(
        f"r1*r2*r3 / integer p          = "
        f"{avg_r123 / avg_p_states:.9f}"
    )

    print(
        f"r1*r2*r3 / prime baseline     = "
        f"{avg_r123 / avg_prime_states:.9f}"
    )

    # ==============================================================================================
    # RECOVERY
    # ==============================================================================================

    print("\n" + "=" * 100)
    print("RECOVERY")
    print("=" * 100)

    print(
        f"correctly recovered           = "
        f"{recovered}/{len(results)}"
    )

    print(
        f"recovery rate                 = "
        f"{100.0 * recovered / len(results):.4f}%"
    )

    # ==============================================================================================
    # QUOTIENT BLOCK DISTRIBUTION
    # ==============================================================================================

    print("\n" + "=" * 100)
    print("QUOTIENT-BLOCK DISTRIBUTION")
    print("=" * 100)

    block_dist = Counter(
        r[2].quotient_blocks
        for r in results
    )

    for blocks, count in sorted(
        block_dist.items()
    ):

        print(
            f"quotient blocks = "
            f"{blocks:8d} "
            f"anchors = "
            f"{count:4d}"
        )

    # ==============================================================================================
    # THREE-MODULUS CANDIDATE DISTRIBUTION
    # ==============================================================================================

    print("\n" + "=" * 100)
    print("THREE-MODULUS CANDIDATE DISTRIBUTION")
    print("=" * 100)

    cdist = Counter(
        r[2].r123_candidates
        for r in results
    )

    for candidates, count in sorted(
        cdist.items()
    ):

        print(
            f"r123 candidates = "
            f"{candidates:6d} "
            f"anchors = "
            f"{count:4d}"
        )

    # ==============================================================================================
    # STRONGEST
    # ==============================================================================================

    print("\n" + "=" * 100)
    print("STRONGEST QUOTIENT-BLOCK COLLAPSES")
    print("=" * 100)

    strongest = sorted(
        results,
        key=lambda r: (
            r[2].quotient_blocks,
            r[2].r123_candidates,
        ),
    )[:20]

    for (
        anchor,
        mods,
        stats,
        p_states,
        prime_states,
        elapsed,
        recovered_flag,
    ) in strongest:

        print(
            f"n={anchor.n:,} "
            f"p={anchor.p:,} "
            f"q={anchor.q:,} "
            f"mods=({mods.r1},{mods.r2},{mods.r3}) "
            f"pStates={p_states:,} "
            f"qBlocks={stats.quotient_blocks:,} "
            f"R1={stats.r1_candidates:,} "
            f"R12={stats.r12_candidates:,} "
            f"R123={stats.r123_candidates:,} "
            f"exact={stats.exact_solutions:,}"
        )

    # ==============================================================================================
    # WEAKEST
    # ==============================================================================================

    print("\n" + "=" * 100)
    print("WEAKEST QUOTIENT-BLOCK COLLAPSES")
    print("=" * 100)

    weakest = sorted(
        results,
        key=lambda r: (
            -r[2].quotient_blocks,
            -r[2].r123_candidates,
        ),
    )[:20]

    for (
        anchor,
        mods,
        stats,
        p_states,
        prime_states,
        elapsed,
        recovered_flag,
    ) in weakest:

        print(
            f"n={anchor.n:,} "
            f"p={anchor.p:,} "
            f"q={anchor.q:,} "
            f"mods=({mods.r1},{mods.r2},{mods.r3}) "
            f"pStates={p_states:,} "
            f"qBlocks={stats.quotient_blocks:,} "
            f"R123={stats.r123_candidates:,} "
            f"exact={stats.exact_solutions:,}"
        )

    # ==============================================================================================
    # EXAMPLES
    # ==============================================================================================

    print("\n" + "=" * 100)
    print("EXAMPLES")
    print("=" * 100)

    for (
        anchor,
        mods,
        stats,
        p_states,
        prime_states,
        elapsed,
        recovered_flag,
    ) in results[:20]:

        print(
            f"n={anchor.n:,} "
            f"p={anchor.p:,} "
            f"q={anchor.q:,} "
            f"mods=({mods.r1},{mods.r2},{mods.r3}) "
            f"R={mods.R123:,} "
            f"pStates={p_states:,} "
            f"primeStates={prime_states:,} "
            f"qBlocks={stats.quotient_blocks:,} "
            f"r1={stats.r1_candidates:,} "
            f"r12={stats.r12_candidates:,} "
            f"r123={stats.r123_candidates:,} "
            f"exact={stats.exact_solutions:,} "
            f"recovered={recovered_flag} "
            f"time={elapsed * 1000:.3f}ms"
        )

    # ==============================================================================================
    # ANCHOR TABLE
    # ==============================================================================================

    print("\n" + "=" * 100)
    print("ANCHOR RESULTS")
    print("=" * 100)

    print(
        " ID "
        "         p "
        "         q "
        "    r1 "
        "    r2 "
        "    r3 "
        "  PSTATES "
        "  QBLOCK "
        "       R1 "
        "      R12 "
        "     R123 "
        " EXACT"
    )

    print("-" * 100)

    for i, (
        anchor,
        mods,
        stats,
        p_states,
        prime_states,
        elapsed,
        recovered_flag,
    ) in enumerate(
        results,
        1,
    ):

        print(
            f"{i:3d} "
            f"{anchor.p:10,d} "
            f"{anchor.q:10,d} "
            f"{mods.r1:6d} "
            f"{mods.r2:6d} "
            f"{mods.r3:6d} "
            f"{p_states:9,d} "
            f"{stats.quotient_blocks:9,d} "
            f"{stats.r1_candidates:9,d} "
            f"{stats.r12_candidates:9,d} "
            f"{stats.r123_candidates:9,d} "
            f"{stats.exact_solutions:5,d}"
        )

    # ==============================================================================================
    # MATHEMATICAL INTERPRETATION
    # ==============================================================================================

    print("\n" + "=" * 100)
    print("MATHEMATICAL INTERPRETATION")
    print("=" * 100)

    print(
        r"""
The previous failure came from treating M as though every anchor
had the same factorization target.

That is incorrect.

Each anchor is an independent:

    n = p*q.

The present experiment uses the correct anchor-specific n.

For a fixed quotient:

    q = floor(n/p),

all p values producing this q form one exact interval:

    floor(n/(q+1)) + 1 <= p <= floor(n/q).

Inside this interval q is constant.

The factor equation modulo the combined modulus

    R = r1*r2*r3

is:

    q*p == n (mod R).

This is a linear congruence in p.

Let:

    g = gcd(q,R).

A solution exists only if:

    g | n.

When that holds:

    (q/g) p == n/g (mod R/g)

and therefore:

    p == p0 (mod R/g).

So the whole quotient block is reduced to the
intersection of:

    one integer interval

with:

    one arithmetic progression.

No scan over every integer p in the block is required.

The search hierarchy is therefore:

    quotient block
          |
          v
    one linear congruence
          |
          v
    arithmetic-progression / interval intersection
          |
          v
    candidate p
          |
          v
    exact n % p
          |
          v
    factor pair.

The important distinction from the previous experiments is that
candidate p values are generated directly by solving a congruence,
rather than by enumerating p and then filtering it.

However, there is still a fundamental issue:

    the quotient blocks themselves must be visited.

Therefore this experiment does NOT establish sublinear factoring.

The decisive measurements are:

    integer p states
    quotient blocks
    three-modulus p candidates.

A promising result would be something like:

    40,000 integer p states
        ->
     4,000 quotient blocks
        ->
         1 CRT candidate.

That would establish real candidate-generation compression.

The next question would then be:

    can multiple adjacent quotient blocks be grouped into a
    larger region where the linear-congruence solution changes
    predictably?

That is the next mathematical bottleneck.

There is another special regime.

If:

    R > p_max,

then for a genuine hyperbola candidate:

    n = p*q

and the modular condition gives:

    R | (p*q-n).

If the construction uses a candidate q and p in the bounded
hyperbola interval, the resulting difference can potentially become
small enough that divisibility by R forces exact equality.

This should be tested experimentally rather than assumed.

Every generated candidate is checked with:

    p*q == n.

Only exact factor pairs count as recovery.
"""
    )

    # ==============================================================================================
    # TIMING
    # ==============================================================================================

    total = time.perf_counter() - total_start

    print("\n" + "=" * 100)
    print("TIMING")
    print("=" * 100)

    print(
        f"total runtime               = "
        f"{total:.3f} seconds"
    )

    print("\n" + "=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()