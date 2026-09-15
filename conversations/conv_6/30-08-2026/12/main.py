#!/usr/bin/env python3

import math
import random
import statistics
from collections import defaultdict

# =============================================================================
# MAX-PRODUCT vs CLOSENESS-MATCHED RANDOM TRIPLE EXPERIMENT
# =============================================================================

M = 111_546_435

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

ACTUAL_TRIALS = 300

SEED = 1_511_464_998

MODULUS_MIN = 300
MODULUS_MAX = 3_000

# Prime triples are restricted by r3/r1 <= 1 + ratio.
CLOSE_RATIOS = [0.05, 0.10, 0.20]

# Random triples sampled per anchor.
RANDOM_TRIPLES_PER_ANCHOR = 100

# Closeness buckets.
#
# R/n >= threshold means the three-prime product is that close to n.
CLOSENESS_THRESHOLDS = [
    0.99,
    0.995,
    0.999,
    0.9995,
    0.9999,
    0.99999,
]


# =============================================================================
# SIEVE
# =============================================================================

def sieve(lo, hi):
    if hi < 2:
        return []

    s = bytearray(b"\x01") * (hi + 1)
    s[0:2] = b"\x00\x00"

    for p in range(2, int(math.isqrt(hi)) + 1):
        if s[p]:
            start = p * p
            s[start:hi + 1:p] = b"\x00" * (
                ((hi - start) // p) + 1
            )

    return [
        x for x in range(max(2, lo), hi + 1)
        if s[x]
    ]


# =============================================================================
# ACTUAL ANCHORS
# =============================================================================

def generate_anchors(primes, count, rng):
    out = []

    while len(out) < count:
        p, q = sorted(rng.sample(primes, 2))
        n = p * q
        out.append((n, p, q))

    return out


# =============================================================================
# BEST CLOSE TRIPLE
# =============================================================================

def best_triple(n, modulus_primes, ratio):
    """
    Find the valid triple

        r1 < r2 < r3
        r3/r1 <= 1+ratio
        r1*r2*r3 < n

    maximizing R.
    """

    limit_ratio = 1.0 + ratio
    root = n ** (1.0 / 3.0)

    local = [
        r for r in modulus_primes
        if root / limit_ratio <= r <= root * limit_ratio
    ]

    best = None
    best_R = -1

    L = len(local)

    for i in range(L):
        r1 = local[i]

        for j in range(i + 1, L):
            r2 = local[j]

            if r2 > r1 * limit_ratio:
                break

            base = r1 * r2

            if base >= n:
                break

            max_r3 = min(
                int((n - 1) // base),
                int(r1 * limit_ratio),
            )

            # local is sorted. Find largest valid r3.
            for k in range(L - 1, j, -1):
                r3 = local[k]

                if r3 > max_r3:
                    continue

                R = base * r3

                if R < n and R > best_R:
                    best_R = R
                    best = (r1, r2, r3, R)

                break

    return best


# =============================================================================
# RANDOM VALID TRIPLE
# =============================================================================

def random_triple(n, modulus_primes, ratio, rng):
    limit_ratio = 1.0 + ratio
    root = n ** (1.0 / 3.0)

    local = [
        r for r in modulus_primes
        if root / limit_ratio <= r <= root * limit_ratio
    ]

    if len(local) < 3:
        return None

    for _ in range(1000):
        r1, r2, r3 = sorted(rng.sample(local, 3))

        if r3 > r1 * limit_ratio:
            continue

        R = r1 * r2 * r3

        if R < n:
            return r1, r2, r3, R

    return None


# =============================================================================
# CRT
# =============================================================================

def crt3(a1, m1, a2, m2, a3, m3):
    """
    Solve

        x = a1 mod m1
        x = a2 mod m2
        x = a3 mod m3

    for pairwise-coprime moduli.
    """

    m12 = m1 * m2
    M = m12 * m3

    inv12 = pow(m1, -1, m2)
    x12 = (a2 - a1) % m2
    x12 = (x12 * inv12) % m2
    x = a1 + m1 * x12

    inv123 = pow(m12, -1, m3)
    x3 = (a3 - x) % m3
    x3 = (x3 * inv123) % m3

    return (x + m12 * x3) % M


# =============================================================================
# CANDIDATE COUNT
# =============================================================================

def candidate_count(
    n,
    p,
    q,
    triple,
    prime_set,
    primes,
):
    """
    Count unordered prime pairs (a,b) from the allowed factor pool
    satisfying

        a*b == n mod r

    for all three r.

    Because R=r1*r2*r3 is normally much larger than FACTOR_MAX,
    CRT turns the required b-residue into an ordinary integer
    candidate in the factor interval.
    """

    r1, r2, r3, R = triple

    count = 0
    actual_seen = False

    nr1 = n % r1
    nr2 = n % r2
    nr3 = n % r3

    for a in primes:
        if a > FACTOR_MAX:
            break

        a1 = a % r1
        a2 = a % r2
        a3 = a % r3

        # All r's are prime and a is a factor prime below 100k.
        if a1 == 0 or a2 == 0 or a3 == 0:
            continue

        b1 = (nr1 * pow(a1, -1, r1)) % r1
        b2 = (nr2 * pow(a2, -1, r2)) % r2
        b3 = (nr3 * pow(a3, -1, r3)) % r3

        b = crt3(
            b1, r1,
            b2, r2,
            b3, r3,
        )

        # CRT gives the unique residue modulo R.
        #
        # Since the factor interval is far below R for essentially
        # all triples here, b itself must be the possible factor.
        if b < FACTOR_MIN or b > FACTOR_MAX:
            continue

        if b not in prime_set:
            continue

        if a > b:
            continue

        count += 1

        if (a == p and b == q) or (a == q and b == p):
            actual_seen = True

    return count, actual_seen


# =============================================================================
# BUCKET
# =============================================================================

def closeness_bucket(ratio):
    """
    Return highest threshold satisfied.
    """

    best = None

    for t in CLOSENESS_THRESHOLDS:
        if ratio >= t:
            best = t

    return best


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

    anchors = generate_anchors(
        factor_primes,
        ACTUAL_TRIALS,
        rng,
    )

    print("=" * 100)
    print("MAX-PRODUCT vs CLOSENESS-MATCHED RANDOM TRIPLE EXPERIMENT")
    print("=" * 100)
    print(f"M                         = {M:,}")
    print(f"factor range              = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"actual anchors            = {ACTUAL_TRIALS}")
    print(f"modulus prime range       = {MODULUS_MIN:,} - {MODULUS_MAX:,}")
    print(f"modulus primes            = {len(modulus_primes)}")
    print(f"random triples / anchor   = {RANDOM_TRIPLES_PER_ANCHOR}")
    print(f"close ratios              = {CLOSE_RATIOS}")
    print(f"seed                      = {SEED:,}")

    # =====================================================================
    # TEST EACH CLOSENESS RULE
    # =====================================================================

    for close_ratio in CLOSE_RATIOS:

        print()
        print("=" * 100)
        print(f"CLOSE-RATIO = {close_ratio:.2%}")
        print("=" * 100)

        best_rows = []

        random_rows = []

        # -----------------------------------------------------------------
        # ACTUAL MAXIMUM-PRODUCT TRIPLE
        # -----------------------------------------------------------------

        print()
        print("BUILDING MAXIMUM-PRODUCT TRIPLES")

        for idx, (n, p, q) in enumerate(anchors, 1):

            triple = best_triple(
                n,
                modulus_primes,
                close_ratio,
            )

            if triple is None:
                continue

            r1, r2, r3, R = triple

            candidates, actual_seen = candidate_count(
                n,
                p,
                q,
                triple,
                prime_set,
                factor_primes,
            )

            ratio = R / n

            best_rows.append({
                "n": n,
                "p": p,
                "q": q,
                "triple": triple,
                "ratio": ratio,
                "gap": n - R,
                "candidates": candidates,
                "unique": candidates == 1,
                "actual_seen": actual_seen,
            })

            if idx % 25 == 0:
                print(f"anchor {idx:3d}/{len(anchors)}")

        # -----------------------------------------------------------------
        # RANDOM TRIPLES
        # -----------------------------------------------------------------

        print()
        print("BUILDING RANDOM TRIPLES")

        total_random = 0

        for idx, (n, p, q) in enumerate(anchors, 1):

            for _ in range(RANDOM_TRIPLES_PER_ANCHOR):

                triple = random_triple(
                    n,
                    modulus_primes,
                    close_ratio,
                    rng,
                )

                if triple is None:
                    continue

                r1, r2, r3, R = triple

                candidates, actual_seen = candidate_count(
                    n,
                    p,
                    q,
                    triple,
                    prime_set,
                    factor_primes,
                )

                ratio = R / n

                random_rows.append({
                    "n": n,
                    "p": p,
                    "q": q,
                    "triple": triple,
                    "ratio": ratio,
                    "gap": n - R,
                    "candidates": candidates,
                    "unique": candidates == 1,
                    "actual_seen": actual_seen,
                })

                total_random += 1

            if idx % 25 == 0:
                print(f"anchor {idx:3d}/{len(anchors)}")

        # =================================================================
        # MAXIMUM-PRODUCT RESULT
        # =================================================================

        print()
        print("=" * 100)
        print("MAXIMUM-PRODUCT RESULT")
        print("=" * 100)

        if best_rows:

            ratios = [x["ratio"] for x in best_rows]
            gaps = [x["gap"] for x in best_rows]
            counts = [x["candidates"] for x in best_rows]

            print(f"anchors solved       = {len(best_rows)}")
            print(f"mean R/n             = {statistics.mean(ratios):.10f}")
            print(f"median R/n           = {statistics.median(ratios):.10f}")
            print(f"minimum R/n          = {min(ratios):.10f}")
            print(f"maximum R/n          = {max(ratios):.10f}")
            print(f"mean gap             = {statistics.mean(gaps):,.2f}")
            print(f"median gap           = {statistics.median(gaps):,.2f}")
            print(f"mean candidates      = {statistics.mean(counts):.6f}")
            print(f"maximum candidates   = {max(counts)}")
            print(
                f"unique factor pair   = "
                f"{sum(x['unique'] for x in best_rows)} / "
                f"{len(best_rows)}"
            )

        # =================================================================
        # RANDOM TRIPLE RESULT
        # =================================================================

        print()
        print("=" * 100)
        print("RANDOM TRIPLE RESULT")
        print("=" * 100)

        print(f"random triples       = {total_random:,}")

        if random_rows:

            ratios = [x["ratio"] for x in random_rows]
            counts = [x["candidates"] for x in random_rows]

            print(f"mean R/n             = {statistics.mean(ratios):.10f}")
            print(f"median R/n           = {statistics.median(ratios):.10f}")
            print(f"maximum R/n          = {max(ratios):.10f}")
            print(f"mean candidates      = {statistics.mean(counts):.6f}")
            print(f"maximum candidates   = {max(counts)}")
            print(
                f"unique factor pair   = "
                f"{sum(x['unique'] for x in random_rows)} / "
                f"{len(random_rows)}"
            )

        # =================================================================
        # CLOSENESS STRATIFICATION
        # =================================================================

        print()
        print("=" * 100)
        print("RANDOM TRIPLES STRATIFIED BY R/n")
        print("=" * 100)

        print(
            f"{'THRESHOLD':>12} "
            f"{'SAMPLES':>10} "
            f"{'MEAN R/n':>14} "
            f"{'MEAN CAND':>14} "
            f"{'UNIQUE %':>12}"
        )
        print("-" * 70)

        for threshold in CLOSENESS_THRESHOLDS:

            rows = [
                x for x in random_rows
                if x["ratio"] >= threshold
            ]

            if not rows:
                print(
                    f"{threshold:12.5f} "
                    f"{0:10d}"
                )
                continue

            mean_ratio = statistics.mean(
                x["ratio"] for x in rows
            )

            mean_candidates = statistics.mean(
                x["candidates"] for x in rows
            )

            unique_fraction = (
                sum(x["unique"] for x in rows)
                / len(rows)
            )

            print(
                f"{threshold:12.5f} "
                f"{len(rows):10d} "
                f"{mean_ratio:14.10f} "
                f"{mean_candidates:14.6f} "
                f"{100*unique_fraction:11.3f}%"
            )

        # =================================================================
        # MAXIMUM-PRODUCT INSIDE RANDOM CLOSENESS DISTRIBUTION
        # =================================================================

        print()
        print("=" * 100)
        print("MAXIMUM-PRODUCT vs RANDOM CLOSENESS-MATCHED NULL")
        print("=" * 100)

        for threshold in CLOSENESS_THRESHOLDS:

            matched = [
                x for x in random_rows
                if x["ratio"] >= threshold
            ]

            actual = [
                x for x in best_rows
                if x["ratio"] >= threshold
            ]

            if not matched or not actual:
                continue

            rand_mean = statistics.mean(
                x["candidates"] for x in matched
            )

            rand_unique = (
                sum(x["unique"] for x in matched)
                / len(matched)
            )

            act_mean = statistics.mean(
                x["candidates"] for x in actual
            )

            act_unique = (
                sum(x["unique"] for x in actual)
                / len(actual)
            )

            print()
            print(f"R/n >= {threshold:.5f}")

            print(
                f"  MAX triple: "
                f"mean candidates={act_mean:.6f} "
                f"unique={100*act_unique:.3f}% "
                f"samples={len(actual)}"
            )

            print(
                f"  RANDOM:     "
                f"mean candidates={rand_mean:.6f} "
                f"unique={100*rand_unique:.3f}% "
                f"samples={len(matched)}"
            )

        # =================================================================
        # BEST TRIPLES
        # =================================================================

        print()
        print("=" * 100)
        print("TOP MAXIMUM-PRODUCT TRIPLES")
        print("=" * 100)

        for x in sorted(
            best_rows,
            key=lambda z: z["ratio"],
            reverse=True,
        )[:20]:

            r1, r2, r3, R = x["triple"]

            print(
                f"n={x['n']:,} "
                f"p={x['p']:,} q={x['q']:,} "
                f"Rmods=({r1},{r2},{r3}) "
                f"R={R:,} "
                f"R/n={x['ratio']:.10f} "
                f"gap={x['gap']:,} "
                f"candidates={x['candidates']}"
            )

    # =====================================================================
    # CONCEPTUAL CONCLUSION
    # =====================================================================

    print()
    print("=" * 100)
    print("INTERPRETATION")
    print("=" * 100)

    print(
        """
The previous experiment established that three close primes can produce
a product R extremely close to n from below.

This experiment asks a more specific question:

    Is the usefulness of (r1,r2,r3) primarily caused by R being close to n?

The comparison is:

    MAXIMUM-PRODUCT triple
            versus
    RANDOM valid triples with similar R/n.

For a triple:

    R = r1*r2*r3

we observe

    n mod r1
    n mod r2
    n mod r3.

For every possible factor a we calculate the corresponding required
factor b using:

    a*b == n (mod r1)
    a*b == n (mod r2)
    a*b == n (mod r3)

and combine those three equations with CRT.

If R is larger than the complete factor interval, the CRT residue of b
essentially determines b directly.

The crucial result is therefore not merely:

    R/n -> 1

but whether:

    larger R/n
        -> fewer surviving factor pairs

and, more importantly,

    MAXIMUM PRODUCT
        performs better than
    RANDOM triples at comparable R/n.

If the maximum-product triples do NOT outperform equally-close random
triples, then "R as close as possible to n" is probably the dominant
explanation.

If they DO outperform matched random triples, then the particular
placement of the three primes (not merely their product size) is
contributing additional structure.
"""
    )

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()
