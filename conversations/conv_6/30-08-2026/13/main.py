#!/usr/bin/env python3

import math
import random
import statistics
from collections import Counter

# =============================================================================
# MAX-PRODUCT vs RANDOM CLOSE-TRIPLE EXPERIMENT
# OPTIMIZED CANDIDATE SEARCH
# =============================================================================

M = 111_546_435

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

ACTUAL_TRIALS = 300

MODULUS_MIN = 300
MODULUS_MAX = 3_000

CLOSE_RATIOS = [0.05, 0.10, 0.20]

# Much smaller because candidate search now stops at 2 matches.
RANDOM_TRIPLES_PER_ANCHOR = 20

SEED = 1_511_464_998


# =============================================================================
# SIEVE
# =============================================================================

def sieve(lo, hi):
    s = bytearray(b"\x01") * (hi + 1)
    s[0:2] = b"\x00\x00"

    for p in range(2, int(math.isqrt(hi)) + 1):
        if s[p]:
            start = p * p
            s[start:hi + 1:p] = b"\x00" * (((hi - start) // p) + 1)

    return [x for x in range(max(2, lo), hi + 1) if s[x]]


# =============================================================================
# ANCHORS
# =============================================================================

def generate_anchors(primes, count, rng):
    result = []

    while len(result) < count:
        p, q = sorted(rng.sample(primes, 2))
        result.append((p * q, p, q))

    return result


# =============================================================================
# MODULAR INVERSE CACHE
# =============================================================================

def build_inverse_tables(moduli):
    tables = {}

    for r in moduli:
        inv = {}

        for x in range(1, r):
            inv[x] = pow(x, -1, r)

        tables[r] = inv

    return tables


# =============================================================================
# PRIME RESIDUES
# =============================================================================

def build_residue_tables(primes, moduli):
    residue = {}

    for r in moduli:
        residue[r] = [p % r for p in primes]

    return residue


# =============================================================================
# BEST CLOSE TRIPLE
# =============================================================================

def best_triple(n, modulus_primes, ratio):
    limit = 1.0 + ratio

    root = n ** (1.0 / 3.0)

    lo = root / limit
    hi = root * limit

    local = [
        p for p in modulus_primes
        if lo <= p <= hi
    ]

    best = None
    best_R = -1

    L = len(local)

    for i in range(L):
        r1 = local[i]

        for j in range(i + 1, L):
            r2 = local[j]

            if r2 > r1 * limit:
                break

            base = r1 * r2

            if base >= n:
                break

            for k in range(L - 1, j, -1):
                r3 = local[k]

                if r3 > r1 * limit:
                    continue

                R = base * r3

                if R < n and R > best_R:
                    best_R = R
                    best = (r1, r2, r3, R)

                break

    return best


# =============================================================================
# RANDOM TRIPLE
# =============================================================================

def random_triple(n, modulus_primes, ratio, rng):
    limit = 1.0 + ratio

    root = n ** (1.0 / 3.0)

    lo = root / limit
    hi = root * limit

    local = [
        p for p in modulus_primes
        if lo <= p <= hi
    ]

    if len(local) < 3:
        return None

    for _ in range(100):

        r1, r2, r3 = sorted(rng.sample(local, 3))

        if r3 > r1 * limit:
            continue

        R = r1 * r2 * r3

        if R < n:
            return (r1, r2, r3, R)

    return None


# =============================================================================
# CRT FOR THREE MODULI
# =============================================================================

def crt3(a1, m1, a2, m2, a3, m3):
    # x = a1 mod m1
    # x = a2 mod m2
    # x = a3 mod m3

    m12 = m1 * m2

    x12 = ((a2 - a1) * pow(m1, -1, m2)) % m2

    x = a1 + m1 * x12

    x3 = ((a3 - x) * pow(m12, -1, m3)) % m3

    return x + m12 * x3


# =============================================================================
# FAST CANDIDATE SEARCH
# =============================================================================

def candidate_count_fast(
    n,
    p,
    q,
    triple,
    primes,
    prime_set,
    inverse_tables,
    residue_tables,
):
    """
    Return:

        0 = no compatible factor pair
        1 = exactly one compatible pair
        2 = at least two compatible pairs

    Also return whether the actual (p,q) is among them.

    We deliberately stop at two candidates because the question here
    is uniqueness, not the complete candidate count.
    """

    r1, r2, r3, R = triple

    # n residues
    n1 = n % r1
    n2 = n % r2
    n3 = n % r3

    inv1 = inverse_tables[r1]
    inv2 = inverse_tables[r2]
    inv3 = inverse_tables[r3]

    res1 = residue_tables[r1]
    res2 = residue_tables[r2]
    res3 = residue_tables[r3]

    found = 0
    actual_seen = False

    for i, a in enumerate(primes):

        a1 = res1[i]
        a2 = res2[i]
        a3 = res3[i]

        # All factor primes are > r in this experiment.
        if a1 == 0 or a2 == 0 or a3 == 0:
            continue

        b1 = (n1 * inv1[a1]) % r1
        b2 = (n2 * inv2[a2]) % r2
        b3 = (n3 * inv3[a3]) % r3

        b = crt3(
            b1, r1,
            b2, r2,
            b3, r3,
        )

        if b < FACTOR_MIN or b > FACTOR_MAX:
            continue

        if b not in prime_set:
            continue

        if a > b:
            continue

        found += 1

        if (a == p and b == q):
            actual_seen = True

        if found >= 2:
            return 2, actual_seen

    return found, actual_seen


# =============================================================================
# STATISTICS
# =============================================================================

def summarize(rows):
    counts = [x["count"] for x in rows]

    return {
        "mean_candidates": statistics.mean(counts),
        "unique_fraction": sum(c == 1 for c in counts) / len(counts),
        "ambiguous_fraction": sum(c >= 2 for c in counts) / len(counts),
        "zero_fraction": sum(c == 0 for c in counts) / len(counts),
        "n": len(rows),
    }


# =============================================================================
# MAIN
# =============================================================================

def run():

    rng = random.Random(SEED)

    factor_primes = sieve(
        FACTOR_MIN,
        FACTOR_MAX,
    )

    modulus_primes = sieve(
        MODULUS_MIN,
        MODULUS_MAX,
    )

    prime_set = set(factor_primes)

    # Only modulus primes that can occur near cube root of the tested n.
    anchors = generate_anchors(
        factor_primes,
        ACTUAL_TRIALS,
        rng,
    )

    print("=" * 100)
    print("FAST MAX-PRODUCT vs RANDOM CLOSE-TRIPLE EXPERIMENT")
    print("=" * 100)

    print(f"M                         = {M:,}")
    print(
        f"factor range              = "
        f"{FACTOR_MIN:,} - {FACTOR_MAX:,}"
    )
    print(f"actual anchors            = {ACTUAL_TRIALS:,}")
    print(
        f"modulus prime range       = "
        f"{MODULUS_MIN:,} - {MODULUS_MAX:,}"
    )
    print(f"modulus primes            = {len(modulus_primes):,}")
    print(
        f"random triples / anchor   = "
        f"{RANDOM_TRIPLES_PER_ANCHOR}"
    )
    print(f"seed                      = {SEED:,}")

    # Only use moduli relevant to factor interval / cube-root region.
    modulus_primes = [
        r for r in modulus_primes
        if r > 2
    ]

    relevant_moduli = sorted(modulus_primes)

    # These are tiny enough to precompute.
    needed = set()

    for n, _, _ in anchors:

        root = n ** (1.0 / 3.0)

        for ratio in CLOSE_RATIOS:

            lo = root / (1.0 + ratio)
            hi = root * (1.0 + ratio)

            for r in modulus_primes:
                if lo <= r <= hi:
                    needed.add(r)

    needed = sorted(needed)

    inverse_tables = build_inverse_tables(
        needed
    )

    residue_tables = build_residue_tables(
        factor_primes,
        needed,
    )

    # =====================================================================
    # EACH RATIO
    # =====================================================================

    for ratio in CLOSE_RATIOS:

        print()
        print("=" * 100)
        print(f"CLOSE RATIO = {ratio:.2%}")
        print("=" * 100)

        actual_rows = []
        random_rows = []

        # ================================================================
        # ACTUAL MAX PRODUCT
        # ================================================================

        print()
        print("ACTUAL MAXIMUM-PRODUCT TRIPLES")

        for index, (n, p, q) in enumerate(
            anchors,
            start=1,
        ):

            triple = best_triple(
                n,
                modulus_primes,
                ratio,
            )

            if triple is None:
                continue

            count, seen = candidate_count_fast(
                n,
                p,
                q,
                triple,
                factor_primes,
                prime_set,
                inverse_tables,
                residue_tables,
            )

            r1, r2, r3, R = triple

            actual_rows.append({
                "n": n,
                "p": p,
                "q": q,
                "triple": triple,
                "R": R,
                "ratio": R / n,
                "gap": n - R,
                "count": count,
                "actual_seen": seen,
            })

            if index % 50 == 0:
                print(
                    f"anchor {index:3d}/{len(anchors)}"
                )

        # ================================================================
        # RANDOM
        # ================================================================

        print()
        print("RANDOM TRIPLES")

        total_attempts = (
            ACTUAL_TRIALS
            * RANDOM_TRIPLES_PER_ANCHOR
        )

        attempt = 0

        for index, (n, p, q) in enumerate(
            anchors,
            start=1,
        ):

            for _ in range(RANDOM_TRIPLES_PER_ANCHOR):

                triple = random_triple(
                    n,
                    modulus_primes,
                    ratio,
                    rng,
                )

                if triple is None:
                    continue

                count, seen = candidate_count_fast(
                    n,
                    p,
                    q,
                    triple,
                    factor_primes,
                    prime_set,
                    inverse_tables,
                    residue_tables,
                )

                r1, r2, r3, R = triple

                random_rows.append({
                    "n": n,
                    "p": p,
                    "q": q,
                    "triple": triple,
                    "R": R,
                    "ratio": R / n,
                    "gap": n - R,
                    "count": count,
                    "actual_seen": seen,
                })

                attempt += 1

            if index % 50 == 0:
                print(
                    f"anchor {index:3d}/{len(anchors)}"
                )

        # ================================================================
        # ACTUAL SUMMARY
        # ================================================================

        a = summarize(actual_rows)

        print()
        print("=" * 100)
        print("ACTUAL MAXIMUM-PRODUCT SUMMARY")
        print("=" * 100)

        ratios = [x["ratio"] for x in actual_rows]
        gaps = [x["gap"] for x in actual_rows]

        print(
            f"anchors solved       = "
            f"{len(actual_rows):,}"
        )
        print(
            f"mean R/n             = "
            f"{statistics.mean(ratios):.10f}"
        )
        print(
            f"median R/n           = "
            f"{statistics.median(ratios):.10f}"
        )
        print(
            f"minimum R/n          = "
            f"{min(ratios):.10f}"
        )
        print(
            f"maximum R/n          = "
            f"{max(ratios):.10f}"
        )
        print(
            f"mean gap             = "
            f"{statistics.mean(gaps):,.2f}"
        )
        print(
            f"median gap           = "
            f"{statistics.median(gaps):,.2f}"
        )
        print(
            f"mean candidate count = "
            f"{a['mean_candidates']:.6f}"
        )
        print(
            f"unique pairs         = "
            f"{sum(x['count'] == 1 for x in actual_rows):,}"
        )
        print(
            f"ambiguous pairs      = "
            f"{sum(x['count'] >= 2 for x in actual_rows):,}"
        )
        print(
            f"zero pairs           = "
            f"{sum(x['count'] == 0 for x in actual_rows):,}"
        )

        # ================================================================
        # RANDOM SUMMARY
        # ================================================================

        r = summarize(random_rows)

        print()
        print("=" * 100)
        print("RANDOM TRIPLE SUMMARY")
        print("=" * 100)

        random_ratios = [
            x["ratio"]
            for x in random_rows
        ]

        random_gaps = [
            x["gap"]
            for x in random_rows
        ]

        print(
            f"random triples      = "
            f"{len(random_rows):,}"
        )

        print(
            f"mean R/n            = "
            f"{statistics.mean(random_ratios):.10f}"
        )
        print(
            f"median R/n          = "
            f"{statistics.median(random_ratios):.10f}"
        )
        print(
            f"maximum R/n         = "
            f"{max(random_ratios):.10f}"
        )
        print(
            f"mean gap            = "
            f"{statistics.mean(random_gaps):,.2f}"
        )
        print(
            f"mean candidates     = "
            f"{r['mean_candidates']:.6f}"
        )
        print(
            f"unique pairs        = "
            f"{sum(x['count'] == 1 for x in random_rows):,}"
        )
        print(
            f"ambiguous pairs     = "
            f"{sum(x['count'] >= 2 for x in random_rows):,}"
        )

        # ================================================================
        # CLOSENESS BUCKETS
        # ================================================================

        print()
        print("=" * 100)
        print("CLOSENESS-BUCKET COMPARISON")
        print("=" * 100)

        thresholds = [
            0.90,
            0.95,
            0.99,
            0.995,
            0.999,
            0.9995,
            0.9999,
        ]

        print(
            f"{'THRESH':>10}"
            f"{'RANDOM N':>12}"
            f"{'MEAN R/N':>14}"
            f"{'MEAN CAND':>14}"
            f"{'UNIQUE %':>12}"
        )

        print("-" * 68)

        for threshold in thresholds:

            subset = [
                x for x in random_rows
                if x["ratio"] >= threshold
            ]

            if not subset:
                continue

            print(
                f"{threshold:10.4f}"
                f"{len(subset):12d}"
                f"{statistics.mean(x['ratio'] for x in subset):14.10f}"
                f"{statistics.mean(x['count'] for x in subset):14.6f}"
                f"{100 * sum(x['count'] == 1 for x in subset) / len(subset):11.3f}%"
            )

        # ================================================================
        # TOP TRIPLES
        # ================================================================

        print()
        print("=" * 100)
        print("TOP ACTUAL MAX-PRODUCT TRIPLES")
        print("=" * 100)

        for x in sorted(
            actual_rows,
            key=lambda z: z["ratio"],
            reverse=True,
        )[:15]:

            r1, r2, r3, R = x["triple"]

            print(
                f"n={x['n']:,} "
                f"p={x['p']:,} "
                f"q={x['q']:,} "
                f"mods=({r1},{r2},{r3}) "
                f"R={R:,} "
                f"R/n={x['ratio']:.10f} "
                f"gap={x['gap']:,} "
                f"candidates={x['count']}"
            )

        # ================================================================
        # RANDOM TRIPLE WITH HIGHEST R/N
        # ================================================================

        print()
        print("=" * 100)
        print("BEST RANDOM TRIPLES")
        print("=" * 100)

        for x in sorted(
            random_rows,
            key=lambda z: z["ratio"],
            reverse=True,
        )[:15]:

            r1, r2, r3, R = x["triple"]

            print(
                f"n={x['n']:,} "
                f"mods=({r1},{r2},{r3}) "
                f"R={R:,} "
                f"R/n={x['ratio']:.10f} "
                f"gap={x['gap']:,} "
                f"candidates={x['count']}"
            )

    # =====================================================================
    # FINAL
    # =====================================================================

    print()
    print("=" * 100)
    print("INTERPRETATION")
    print("=" * 100)

    print(
        """
The experiment separates two effects.

1. MODULUS SIZE

    R = r1*r2*r3

    If R approaches n from below, the combined CRT residue has a very
    large modulus relative to the factor interval.

2. PARTICULAR TRIPLE

    Two triples can have similar R/n but completely different
    residue structure.

The important comparison is therefore:

    MAX-PRODUCT TRIPLE
            vs
    RANDOM TRIPLES WITH SIMILAR R/n.

The candidate search is deliberately reported as:

    0  = no compatible factor pair
    1  = unique compatible factor pair
    2  = two or more compatible pairs

because uniqueness is the important question.

The strongest result would be:

    maximum-product triples have substantially more unique
    factor identification than random triples having comparable R/n.

If the two populations behave similarly after conditioning on R/n,
then the main mechanism is simply the CRT modulus becoming larger
than the factor interval.

If maximum-product triples remain better even after matching
closeness, then the individual locations of r1,r2,r3 may contain
additional information beyond their product R.
"""
    )

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()
