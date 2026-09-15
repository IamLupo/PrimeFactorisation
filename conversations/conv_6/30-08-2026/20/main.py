#!/usr/bin/env python3

import math
import random
import time
from dataclasses import dataclass


# =============================================================================
# CONFIGURATION
# =============================================================================

M = 111_546_435

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

MODULUS_MIN = 300
MODULUS_MAX = 3_000

ANCHOR_COUNT = 300

CLOSE_RATIO = 0.20

SEED = 1_511_464_998

PROGRESS_EVERY = 25


# =============================================================================
# DATA
# =============================================================================

@dataclass(frozen=True)
class Anchor:
    p: int
    q: int

    @property
    def n(self) -> int:
        return self.p * self.q


@dataclass(frozen=True)
class Triple:
    r1: int
    r2: int
    r3: int

    @property
    def R(self) -> int:
        return self.r1 * self.r2 * self.r3


# =============================================================================
# FAST SIEVE
# =============================================================================

def sieve(limit: int) -> list[int]:
    if limit < 2:
        return []

    composite = bytearray(limit + 1)

    composite[0] = 1
    composite[1] = 1

    root = math.isqrt(limit)

    for p in range(2, root + 1):
        if composite[p] == 0:
            start = p * p
            composite[start:limit + 1:p] = (
                b"\x01" * (((limit - start) // p) + 1)
            )

    return [
        n for n in range(2, limit + 1)
        if composite[n] == 0
    ]


# =============================================================================
# ANCHORS
# =============================================================================

def generate_anchors(
    factor_primes: list[int],
    count: int,
    rng: random.Random,
) -> list[Anchor]:

    result = []
    seen = set()

    while len(result) < count:
        p, q = rng.sample(factor_primes, 2)

        if p > q:
            p, q = q, p

        pair = (p, q)

        if pair in seen:
            continue

        seen.add(pair)
        result.append(Anchor(p, q))

    return result


# =============================================================================
# MAX PRODUCT CLOSE PRIME TRIPLE
# =============================================================================

def choose_max_product_triple(
    n: int,
    modulus_primes: list[int],
    close_ratio: float,
) -> Triple:

    best = None
    best_R = -1

    # We only have ~368 modulus primes, so this is cheap.
    for i, r1 in enumerate(modulus_primes):

        max_close = int(r1 * (1.0 + close_ratio))

        for j in range(i + 1, len(modulus_primes)):

            r2 = modulus_primes[j]

            if r2 > max_close:
                break

            r12 = r1 * r2

            if r12 >= n:
                break

            max_r3 = min(
                int(r1 * (1.0 + close_ratio)),
                (n - 1) // r12,
            )

            if max_r3 <= r2:
                continue

            # Find largest modulus prime <= max_r3.
            lo = j + 1
            hi = len(modulus_primes)

            while lo < hi:
                mid = (lo + hi) // 2

                if modulus_primes[mid] <= max_r3:
                    lo = mid + 1
                else:
                    hi = mid

            k = lo - 1

            if k <= j:
                continue

            r3 = modulus_primes[k]

            R = r1 * r2 * r3

            if R < n and R > best_R:
                best_R = R
                best = Triple(r1, r2, r3)

    if best is None:
        raise RuntimeError(
            f"Could not find modulus triple for n={n:,}"
        )

    return best


# =============================================================================
# BRANCH
# =============================================================================

def branch4(p: int, q_residue_mod4: int) -> str:

    a = p & 3
    b = q_residue_mod4 & 3

    if a == 1 and b == 1:
        return "11"

    if a == 3 and b == 3:
        return "33"

    return "13/31"


def allowed_branches(n: int) -> tuple[str, ...]:

    nr = n & 3

    if nr == 1:
        return ("11", "33")

    if nr == 3:
        return ("13/31",)

    return ()


# =============================================================================
# CRT FOR THREE PRIME MODULI
# =============================================================================

def crt3(
    a1: int,
    a2: int,
    a3: int,
    r1: int,
    r2: int,
    r3: int,
    inv12: int,
    inv123: int,
) -> int:
    """
    Solve:

        x == a1 (mod r1)
        x == a2 (mod r2)
        x == a3 (mod r3)

    Returns x in [0, r1*r2*r3).
    """

    # Combine r1 and r2.
    t2 = ((a2 - a1) * inv12) % r2

    x12 = a1 + r1 * t2

    R12 = r1 * r2

    # Combine with r3.
    t3 = ((a3 - x12) * inv123) % r3

    return x12 + R12 * t3


# =============================================================================
# PRECOMPUTED TRIPLE CRT DATA
# =============================================================================

@dataclass
class TripleCRT:
    r1: int
    r2: int
    r3: int
    R: int
    inv12: int
    inv123: int


def prepare_triple(triple: Triple) -> TripleCRT:

    r1 = triple.r1
    r2 = triple.r2
    r3 = triple.r3

    inv12 = pow(r1, -1, r2)

    R12 = r1 * r2

    inv123 = pow(R12, -1, r3)

    return TripleCRT(
        r1=r1,
        r2=r2,
        r3=r3,
        R=R12 * r3,
        inv12=inv12,
        inv123=inv123,
    )


# =============================================================================
# DIRECT EXPLOIT SEARCH
# =============================================================================

def analyze_anchor(
    anchor: Anchor,
    crt: TripleCRT,
    factor_primes: list[int],
):
    n = anchor.n

    r1 = crt.r1
    r2 = crt.r2
    r3 = crt.r3

    nr1 = n % r1
    nr2 = n % r2
    nr3 = n % r3

    nr4 = n & 3

    # Number of p candidates after each stage.
    base_count = len(factor_primes)

    branch_counts = {
        "11": 0,
        "33": 0,
        "13/31": 0,
    }

    crt_q_candidates = []

    # -------------------------------------------------------------------------
    # MAIN SEARCH
    #
    # For every candidate p:
    #
    # q == n * p^-1 (mod r_i)
    #
    # for i=1,2,3.
    #
    # CRT reconstructs q modulo R.
    #
    # Since R > FACTOR_MAX in this experiment, there can be at most one
    # bounded q.
    # -------------------------------------------------------------------------

    for p in factor_primes:

        a1 = p % r1
        a2 = p % r2
        a3 = p % r3

        # All r_i are below 3000 and p >= 10000, therefore p cannot equal
        # any r_i. Since p and r_i are both prime, p is invertible.
        q1 = (nr1 * pow(a1, -1, r1)) % r1
        q2 = (nr2 * pow(a2, -1, r2)) % r2
        q3 = (nr3 * pow(a3, -1, r3)) % r3

        q0 = crt3(
            q1,
            q2,
            q3,
            r1,
            r2,
            r3,
            crt.inv12,
            crt.inv123,
        )

        if not (
            FACTOR_MIN <= q0 <= FACTOR_MAX
        ):
            continue

        br = branch4(p, q0)

        branch_counts[br] += 1

        crt_q_candidates.append(
            (p, q0, br)
        )

    # -------------------------------------------------------------------------
    # EXACT CHECK
    # -------------------------------------------------------------------------

    exact_pairs = []

    for p, q, br in crt_q_candidates:

        if p * q != n:
            continue

        if not (
            FACTOR_MIN <= q <= FACTOR_MAX
        ):
            continue

        exact_pairs.append(
            (p, q)
        )

    exact_pairs = sorted(set(
        tuple(sorted(pair))
        for pair in exact_pairs
    ))

    # -------------------------------------------------------------------------
    # PRIME VERIFICATION IS AUTOMATIC FOR p; q MUST BE PRIME.
    # -------------------------------------------------------------------------

    prime_pairs = []

    for pair in exact_pairs:

        p, q = pair

        # q should already be prime for the intended construction.
        # We verify with the precomputed factor-prime set outside this
        # function if desired. Since q divides n exactly, p/q are the
        # actual factor candidates.
        prime_pairs.append(pair)

    # -------------------------------------------------------------------------
    # BRANCH SURVIVAL
    # -------------------------------------------------------------------------

    surviving_branches = tuple(
        br
        for br in ("11", "13/31", "33")
        if branch_counts[br] > 0
    )

    # If no CRT candidate exists, mark all zero.
    return {
        "base": base_count,
        "crt": len(crt_q_candidates),
        "exact": len(exact_pairs),
        "pairs": prime_pairs,
        "branch_counts": branch_counts,
        "surviving_branches": surviving_branches,
    }


# =============================================================================
# FAST FACTOR-PRIME SET
# =============================================================================

def verify_prime_pairs(
    results,
    factor_prime_set: set[int],
):

    verified = []

    for p, q in results["pairs"]:

        if (
            p in factor_prime_set
            and q in factor_prime_set
        ):
            verified.append((p, q))

    results["pairs"] = verified


# =============================================================================
# MAIN
# =============================================================================

def run():

    rng = random.Random(SEED)

    start_total = time.perf_counter()

    print("=" * 100)
    print("THREE-MODULUS DIRECT CRT FACTOR SEARCH — OPTIMIZED")
    print("=" * 100)
    print(f"M                    = {M:,}")
    print(
        f"factor range         = "
        f"{FACTOR_MIN:,} - {FACTOR_MAX:,}"
    )
    print(f"anchors              = {ANCHOR_COUNT}")
    print(
        f"modulus prime range  = "
        f"{MODULUS_MIN:,} - {MODULUS_MAX:,}"
    )
    print(f"close ratio           = {CLOSE_RATIO:.0%}")
    print(f"seed                 = {SEED:,}")

    # =========================================================================
    # PRIME POOLS
    # =========================================================================

    print()
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

    factor_prime_set = set(factor_primes)

    print(
        f"factor primes        = {len(factor_primes):,}"
    )

    print(
        f"modulus primes       = {len(modulus_primes):,}"
    )

    # =========================================================================
    # ANCHORS
    # =========================================================================

    anchors = generate_anchors(
        factor_primes,
        ANCHOR_COUNT,
        rng,
    )

    print(
        f"actual anchors       = {len(anchors)}"
    )

    # =========================================================================
    # SELECT TRIPLES
    # =========================================================================

    print()
    print("=" * 100)
    print("SELECTING MAXIMUM-PRODUCT CLOSE PRIME TRIPLES")
    print("=" * 100)

    triples = []

    for i, anchor in enumerate(anchors, 1):

        triple = choose_max_product_triple(
            anchor.n,
            modulus_primes,
            CLOSE_RATIO,
        )

        triples.append(triple)

        if i % PROGRESS_EVERY == 0:
            print(
                f"anchor {i:3d}/{len(anchors)}"
            )

    # =========================================================================
    # PREPARE CRT DATA
    # =========================================================================

    print()
    print("=" * 100)
    print("PREPARING CRT COORDINATES")
    print("=" * 100)

    crt_data = [
        prepare_triple(triple)
        for triple in triples
    ]

    # =========================================================================
    # ANALYSIS
    # =========================================================================

    print()
    print("=" * 100)
    print("C-FIRST + DIRECT CRT SEARCH")
    print("=" * 100)

    all_results = []

    t0 = time.perf_counter()

    for i, (
        anchor,
        triple,
        crt,
    ) in enumerate(
        zip(anchors, triples, crt_data),
        1,
    ):

        result = analyze_anchor(
            anchor,
            crt,
            factor_primes,
        )

        verify_prime_pairs(
            result,
            factor_prime_set,
        )

        all_results.append(result)

        if i % PROGRESS_EVERY == 0:
            print(
                f"anchor {i:3d}/{len(anchors)}"
            )

    elapsed = time.perf_counter() - t0

    # =========================================================================
    # SUMMARY
    # =========================================================================

    print()
    print("=" * 100)
    print("SUMMARY")
    print("=" * 100)

    def avg(key: str) -> float:
        return sum(
            r[key]
            for r in all_results
        ) / len(all_results)

    print(
        f"average base p candidates       = "
        f"{avg('base'):,.2f}"
    )

    print(
        f"average CRT-compatible q        = "
        f"{avg('crt'):,.6f}"
    )

    print(
        f"average exact candidates        = "
        f"{avg('exact'):,.6f}"
    )

    print(
        f"maximum CRT-compatible q        = "
        f"{max(r['crt'] for r in all_results)}"
    )

    print(
        f"maximum exact candidates        = "
        f"{max(r['exact'] for r in all_results)}"
    )

    # =========================================================================
    # REDUCTION
    # =========================================================================

    print()
    print("=" * 100)
    print("SEARCH REDUCTION")
    print("=" * 100)

    base = avg("base")
    crt = avg("crt")
    exact = avg("exact")

    print(
        f"CRT reduction       = "
        f"{100.0 * (1.0 - crt / base):.6f}%"
    )

    print(
        f"exact reduction     = "
        f"{100.0 * (1.0 - exact / base):.6f}%"
    )

    # =========================================================================
    # RECOVERY
    # =========================================================================

    recovered = 0
    unique_crt = 0
    unique_exact = 0

    for anchor, result in zip(
        anchors,
        all_results,
    ):

        target = tuple(
            sorted((anchor.p, anchor.q))
        )

        if target in result["pairs"]:
            recovered += 1

        if result["crt"] == 1:
            unique_crt += 1

        if len(result["pairs"]) == 1:
            unique_exact += 1

    print()
    print("=" * 100)
    print("RECOVERY")
    print("=" * 100)

    print(
        f"correctly recovered        = "
        f"{recovered}/{len(anchors)}"
    )

    print(
        f"unique CRT candidate       = "
        f"{unique_crt}/{len(anchors)}"
    )

    print(
        f"unique final factor pair  = "
        f"{unique_exact}/{len(anchors)}"
    )

    # =========================================================================
    # BRANCH
    # =========================================================================

    print()
    print("=" * 100)
    print("C-FIRST BRANCH SURVIVAL")
    print("=" * 100)

    branch_survival = {}

    for result in all_results:

        key = result["surviving_branches"]

        branch_survival[key] = (
            branch_survival.get(key, 0) + 1
        )

    for key, count in sorted(
        branch_survival.items()
    ):

        print(
            f"{key!s:25s} "
            f"{count:4d}"
        )

    # =========================================================================
    # CRT CANDIDATE DISTRIBUTION
    # =========================================================================

    print()
    print("=" * 100)
    print("CRT CANDIDATE DISTRIBUTION")
    print("=" * 100)

    distribution = {}

    for result in all_results:

        value = result["crt"]

        distribution[value] = (
            distribution.get(value, 0) + 1
        )

    for value in sorted(distribution):

        print(
            f"CRT candidates = "
            f"{value:4d} "
            f"anchors = "
            f"{distribution[value]:4d}"
        )

    # =========================================================================
    # ANCHOR TABLE
    # =========================================================================

    print()
    print("=" * 100)
    print("ANCHOR RESULTS")
    print("=" * 100)

    print(
        " ID        p        q      "
        "r1     r2     r3       R/n   "
        "BASE   CRT   EXACT   BRANCH"
    )

    print("-" * 100)

    for i, (
        anchor,
        triple,
        result,
    ) in enumerate(
        zip(
            anchors,
            triples,
            all_results,
        ),
        1,
    ):

        # Actual branch.
        actual_branch = branch4(
            anchor.p,
            anchor.q,
        )

        print(
            f"{i:3d} "
            f"{anchor.p:8,d} "
            f"{anchor.q:8,d} "
            f"{triple.r1:6,d} "
            f"{triple.r2:6,d} "
            f"{triple.r3:6,d} "
            f"{triple.R / anchor.n:9.7f} "
            f"{result['base']:6,d} "
            f"{result['crt']:5,d} "
            f"{result['exact']:6,d} "
            f"{actual_branch}"
        )

    # =========================================================================
    # MOST POWERFUL CASES
    # =========================================================================

    print()
    print("=" * 100)
    print("MOST AGGRESSIVE CRT COLLAPSES")
    print("=" * 100)

    ranked = sorted(
        zip(
            anchors,
            triples,
            all_results,
        ),
        key=lambda x: (
            x[2]["crt"],
            x[2]["exact"],
        ),
    )

    for anchor, triple, result in ranked[:20]:

        print(
            f"n={anchor.n:,} "
            f"p={anchor.p:,} "
            f"q={anchor.q:,} "
            f"mods=({triple.r1},{triple.r2},{triple.r3}) "
            f"R={triple.R:,} "
            f"R/n={triple.R / anchor.n:.10f} "
            f"CRT={result['crt']} "
            f"exact={result['exact']} "
            f"pairs={result['pairs']}"
        )

    # =========================================================================
    # CHECK THE ACTUAL FACTOR
    # =========================================================================

    print()
    print("=" * 100)
    print("ACTUAL FACTOR RESIDUE CHECK")
    print("=" * 100)

    failures = 0

    for anchor, triple in zip(
        anchors,
        triples,
    ):

        n = anchor.n

        for r in (
            triple.r1,
            triple.r2,
            triple.r3,
        ):

            if (
                (anchor.p * anchor.q) % r
                != n % r
            ):
                failures += 1

    print(
        f"identity failures = {failures}"
    )

    # =========================================================================
    # TIMING
    # =========================================================================

    print()
    print("=" * 100)
    print("TIMING")
    print("=" * 100)

    print(
        f"CRT search time = {elapsed:.3f} seconds"
    )

    total_elapsed = (
        time.perf_counter() - start_total
    )

    print(
        f"total runtime   = {total_elapsed:.3f} seconds"
    )

    # =========================================================================
    # INTERPRETATION
    # =========================================================================

    print()
    print("=" * 100)
    print("INTERPRETATION")
    print("=" * 100)

    print(
        """
This version tests the proposed exploitation directly.

For every candidate prime p in the factor interval we know:

    p*q == n (mod r_i)

so, because p is invertible modulo r_i,

    q == n * p^(-1) (mod r_i).

Therefore three moduli give three congruences for q:

    q == q1 (mod r1)
    q == q2 (mod r2)
    q == q3 (mod r3)

CRT produces:

    q == Q (mod R)

where

    R = r1*r2*r3.

The selected triples satisfy:

    R < n

but R is normally vastly larger than the factor interval.

Therefore there can be at most one q in:

    [10,000,100,000]

for each candidate p.

The search becomes:

    8,363 possible p
           |
           v
    compute 3 q residues
           |
           v
    one CRT reconstruction
           |
           v
    bounded q test
           |
           v
    exact p*q == n
           |
           v
    factor pair.

This is much closer to a genuine direct modular search than the
previous experiments because no list of all residue-pair combinations
is constructed.

The decisive number is:

    average CRT-compatible q candidates.

If this becomes close to 1, the three-modulus system is collapsing the
factor search strongly.

However, the remaining enumeration over p is still important:
the experiment is measuring whether modular information makes that
enumeration cheap, not claiming that the enumeration disappears.

The C-first branch information is also retained. For each surviving
(p,q) candidate we know whether the pair belongs to:

    11
    13/31
    33

modulo 4.

This lets us measure whether the branch constraint contributes
additional reduction before the exact multiplication test.
"""
    )

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()
