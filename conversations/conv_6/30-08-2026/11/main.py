#!/usr/bin/env python3

import math
import random
import statistics
from collections import defaultdict

# =============================================================================
# THREE-CLOSE-PRIME MODULUS / MAXIMUM-PRODUCT EXPERIMENT
# =============================================================================

M = 111_546_435

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

ACTUAL_TRIALS = 300

# Three primes must satisfy:
#
#     r1 < r2 < r3
#
# and
#
#     r3 / r1 <= 1 + CLOSE_RATIO
#
# while
#
#     r1*r2*r3 < n
#
# The triple with maximum product is selected.
#
CLOSE_RATIOS = [
    0.05,   # within 5%
    0.10,   # within 10%
    0.20,   # within 20%
]

# Number of random close-prime triples used as a comparison.
RANDOM_TRIPLES_PER_ANCHOR = 500

SEED = 1_511_464_998

# Prime-modulus search range.
#
# Since n is around 10^8 ... 10^10 in this experiment,
# cube_root(n) is approximately 500 ... 2200.
#
MODULUS_MIN = 300
MODULUS_MAX = 3000


# =============================================================================
# PRIME GENERATION
# =============================================================================

def sieve_primes(lo, hi):
    """Return all primes in [lo, hi]."""
    if hi < 2:
        return []

    sieve = bytearray(b"\x01") * (hi + 1)
    sieve[0:2] = b"\x00\x00"

    limit = int(math.isqrt(hi))

    for p in range(2, limit + 1):
        if sieve[p]:
            start = p * p
            sieve[start:hi + 1:p] = b"\x00" * (
                ((hi - start) // p) + 1
            )

    return [
        x for x in range(max(2, lo), hi + 1)
        if sieve[x]
    ]


# =============================================================================
# FACTOR POOL
# =============================================================================

def build_factor_pool():
    return sieve_primes(FACTOR_MIN, FACTOR_MAX)


# =============================================================================
# ACTUAL SEMIPRIME ANCHORS
# =============================================================================

def generate_actual_anchors(primes, count, rng):
    anchors = []

    while len(anchors) < count:
        p = rng.choice(primes)
        q = rng.choice(primes)

        if p >= q:
            continue

        n = p * q

        anchors.append((n, p, q))

    return anchors


# =============================================================================
# PRIME RESIDUE SIGNATURE
# =============================================================================

def factor_signature(p, q, moduli):
    """
    Ordered factor fingerprint.

    We keep the smaller factor first so that (p,q) and (q,p)
    are treated identically.
    """
    if p > q:
        p, q = q, p

    return tuple(
        (p % r, q % r)
        for r in moduli
    )


# =============================================================================
# INDEX PRIME FINGERPRINTS
# =============================================================================

def build_fingerprint_index(primes, moduli):
    """
    Map:

        fingerprint -> primes

    This allows rapid determination of whether a residue vector
    identifies a unique prime.
    """
    index = defaultdict(list)

    for p in primes:
        sig = tuple(p % r for r in moduli)
        index[sig].append(p)

    return index


# =============================================================================
# FIND MAXIMUM-PRODUCT CLOSE PRIME TRIPLE
# =============================================================================

def best_close_prime_triple(n, primes, ratio):
    """
    Find

        r1 < r2 < r3

    such that

        r1*r2*r3 < n
        r3/r1 <= 1+ratio

    and maximize r1*r2*r3.

    Returns:

        (r1, r2, r3, product)

    or None.
    """

    max_ratio = 1.0 + ratio
    best = None
    best_product = -1

    # We only need primes around cube-root(n).
    #
    # The optimum under a "close primes" restriction will not
    # be arbitrarily far from cube_root(n).
    root = n ** (1.0 / 3.0)

    # Broad but still small search window.
    lower_bound = max(MODULUS_MIN, int(root / max_ratio) - 5)
    upper_bound = min(MODULUS_MAX, int(root * max_ratio) + 20)

    local = [
        p for p in primes
        if lower_bound <= p <= upper_bound
    ]

    L = len(local)

    for i in range(L):
        r1 = local[i]

        # r2 cannot be too large.
        r2_limit = r1 * max_ratio

        for j in range(i + 1, L):
            r2 = local[j]

            if r2 > r2_limit:
                break

            base = r1 * r2

            if base >= n:
                break

            # Largest possible r3 is limited by:
            #
            #   r3 < n/(r1*r2)
            #
            # and by the closeness condition.
            max_r3 = min(
                int((n - 1) // base),
                int(r1 * max_ratio)
            )

            # Find the largest prime <= max_r3.
            # Since local is sorted, walk backwards.
            for k in range(L - 1, j, -1):
                r3 = local[k]

                if r3 > max_r3:
                    continue

                product = base * r3

                if product >= n:
                    continue

                if product > best_product:
                    best_product = product
                    best = (r1, r2, r3, product)

                break

    return best


# =============================================================================
# RANDOM CLOSE-PRIME TRIPLE
# =============================================================================

def random_close_prime_triple(n, primes, ratio, rng):
    """
    Randomly sample valid close-prime triples.

    Returns:
        (r1,r2,r3,R)
    or None
    """

    max_ratio = 1.0 + ratio
    root = n ** (1.0 / 3.0)

    local = [
        p
        for p in primes
        if root / max_ratio <= p <= root * max_ratio
    ]

    if len(local) < 3:
        return None

    for _ in range(1000):
        r1, r2, r3 = sorted(rng.sample(local, 3))

        if r3 > r1 * max_ratio:
            continue

        R = r1 * r2 * r3

        if R < n:
            return r1, r2, r3, R

    return None


# =============================================================================
# CANDIDATE FACTOR MATCHING
# =============================================================================

def count_matching_prime_pairs(
    n,
    actual_p,
    actual_q,
    moduli,
    prime_signature_index,
):
    """
    Given the observations

        n mod r

    determine all prime pairs (a,b) from the factor pool whose
    factor residues are compatible with:

        a*b == n (mod r)

    for every r.

    This is the actual factor-identification test.
    """

    pairs = []

    for a in FACTOR_PRIMES:
        if a > math.isqrt(n):
            break

        ok_a = True

        b_residues = []

        for r in moduli:
            a_r = a % r
            n_r = n % r

            # a must be invertible modulo r.
            if math.gcd(a_r, r) != 1:
                ok_a = False
                break

            inv = pow(a_r, -1, r)
            b_r = (n_r * inv) % r
            b_residues.append(b_r)

        if not ok_a:
            continue

        sig = tuple(b_residues)
        candidates = prime_signature_index.get(sig, ())

        for b in candidates:
            if b < a:
                continue

            if a * b != n:
                # We are measuring possible factorization candidates
                # consistent with residues, not necessarily exact n.
                #
                # For the fingerprint experiment, any pair whose
                # residues match is a candidate.
                pass

            pairs.append((a, b))

    return pairs


# =============================================================================
# IMPORTANT:
# For a fingerprint of FACTORS, we actually want collision buckets
# among the complete prime pair signatures.
#
# Build them directly instead of scanning every prime for every anchor.
# =============================================================================

def build_pair_signature_index(primes, moduli):
    index = defaultdict(list)

    for i, p in enumerate(primes):
        for q in primes[i:]:
            sig = factor_signature(p, q, moduli)
            index[sig].append((p, q))

    return index


# =============================================================================
# MEMORY-SAFE PAIR RESOLUTION
# =============================================================================
#
# The complete pair index would be ~35 million pairs for 8363 primes,
# which is too large.
#
# Therefore we resolve each actual anchor using CRT residue constraints
# one prime at a time.
# =============================================================================

def compatible_factor_candidates(
    n,
    primes,
    moduli,
):
    """
    Find candidate unordered prime pairs (a,b) such that:

        a*b == n (mod r)

    for every r in moduli.

    This is performed incrementally so no gigantic pair table is built.
    """

    # Start with all primes as possible first factors.
    candidate_a = primes

    # Filter using the first factor coordinate.
    #
    # For each a we calculate the required residue of b.
    #
    # Instead of making a huge pair table, construct residue lookup
    # tables for b.
    residue_indices = {}

    for r in moduli:
        d = defaultdict(list)

        for p in primes:
            d[p % r].append(p)

        residue_indices[r] = d

    candidates = []

    # We enforce a <= b.
    for a in primes:
        if a > FACTOR_MAX:
            break

        valid = True
        candidate_sets = []

        for r in moduli:
            ar = a % r
            nr = n % r

            if math.gcd(ar, r) != 1:
                valid = False
                break

            br = (nr * pow(ar, -1, r)) % r

            candidate_sets.append(
                residue_indices[r].get(br, ())
            )

        if not valid:
            continue

        if not candidate_sets:
            continue

        # Intersect the smallest residue bucket first.
        candidate_sets.sort(key=len)

        base = set(candidate_sets[0])

        for s in candidate_sets[1:]:
            base.intersection_update(s)

            if not base:
                break

        for b in base:
            if a <= b:
                candidates.append((a, b))

    return candidates


# =============================================================================
# ANALYZE ONE ANCHOR
# =============================================================================

def analyze_anchor(
    n,
    p,
    q,
    moduli,
    primes,
    residue_indices,
):
    actual_sig = factor_signature(p, q, moduli)

    # Build candidates from the factor pool.
    candidate_pairs = []

    for a in primes:
        valid = True
        required = []

        for r in moduli:
            ar = a % r
            nr = n % r

            if math.gcd(ar, r) != 1:
                valid = False
                break

            br = (nr * pow(ar, -1, r)) % r
            required.append(br)

        if not valid:
            continue

        # Intersect residue classes of q.
        buckets = [
            residue_indices[r].get(br, ())
            for r, br in zip(moduli, required)
        ]

        if not buckets:
            continue

        buckets.sort(key=len)

        possible_b = set(buckets[0])

        for bucket in buckets[1:]:
            possible_b.intersection_update(bucket)

            if not possible_b:
                break

        for b in possible_b:
            if a <= b:
                candidate_pairs.append((a, b))

    unique_pairs = len(candidate_pairs)

    actual_present = tuple(sorted((p, q))) in candidate_pairs

    return {
        "actual_signature": actual_sig,
        "candidate_pairs": candidate_pairs,
        "candidate_count": unique_pairs,
        "actual_present": actual_present,
        "unique": unique_pairs == 1,
    }


# =============================================================================
# MAIN
# =============================================================================

def run():

    global FACTOR_PRIMES
    FACTOR_PRIMES = build_factor_pool()

    rng = random.Random(SEED)

    print("=" * 100)
    print("THREE-CLOSE-PRIME MAXIMUM-PRODUCT MODULUS EXPERIMENT")
    print("=" * 100)
    print(f"M                    = {M:,}")
    print(f"factor range         = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"actual anchors       = {ACTUAL_TRIALS}")
    print(f"modulus search       = {MODULUS_MIN:,} - {MODULUS_MAX:,}")
    print(f"close ratios         = {CLOSE_RATIOS}")
    print(f"random triples/anchor= {RANDOM_TRIPLES_PER_ANCHOR}")
    print(f"seed                 = {SEED:,}")

    primes = sieve_primes(MODULUS_MIN, MODULUS_MAX)

    print()
    print("=" * 100)
    print("MODULUS PRIME POOL")
    print("=" * 100)
    print(f"modulus primes = {len(primes)}")

    anchors = generate_actual_anchors(
        FACTOR_PRIMES,
        ACTUAL_TRIALS,
        rng,
    )

    print()
    print("=" * 100)
    print("ACTUAL ANCHORS")
    print("=" * 100)
    print(f"anchors = {len(anchors)}")

    all_results = {}

    for ratio in CLOSE_RATIOS:

        print()
        print("=" * 100)
        print(f"CLOSE-RATIO = {ratio:.2%}")
        print("=" * 100)

        results = []

        for idx, (n, p, q) in enumerate(anchors, 1):

            triple = best_close_prime_triple(
                n,
                primes,
                ratio,
            )

            if triple is None:
                continue

            r1, r2, r3, R = triple

            gap = n - R
            ratio_to_n = R / n

            results.append({
                "n": n,
                "p": p,
                "q": q,
                "r1": r1,
                "r2": r2,
                "r3": r3,
                "R": R,
                "gap": gap,
                "R_over_n": ratio_to_n,
            })

            if idx % 25 == 0:
                print(f"anchor {idx:3d}/{len(anchors)}")

        all_results[ratio] = results

        print()
        print("SUMMARY")
        print("-" * 100)

        values = [x["R_over_n"] for x in results]
        gaps = [x["gap"] for x in results]

        if values:
            print(f"anchors solved       = {len(results)}")
            print(f"mean R/n             = {statistics.mean(values):.10f}")
            print(f"median R/n           = {statistics.median(values):.10f}")
            print(f"minimum R/n          = {min(values):.10f}")
            print(f"maximum R/n          = {max(values):.10f}")
            print(f"mean n-R             = {statistics.mean(gaps):,.2f}")
            print(f"minimum n-R          = {min(gaps):,}")
            print(f"maximum n-R          = {max(gaps):,}")

        # ---------------------------------------------------------------------
        # PRINT SELECTED TRIPLES
        # ---------------------------------------------------------------------

        print()
        print("BEST MAX-PRODUCT TRIPLES")
        print("-" * 100)

        ranked = sorted(
            results,
            key=lambda x: x["R_over_n"],
            reverse=True,
        )

        for x in ranked[:20]:
            print(
                f"n={x['n']:,} "
                f"p={x['p']:,} q={x['q']:,} "
                f"mods=({x['r1']},{x['r2']},{x['r3']}) "
                f"R={x['R']:,} "
                f"R/n={x['R_over_n']:.10f} "
                f"gap={x['gap']:,}"
            )

    # =========================================================================
    # FACTOR FINGERPRINT RESOLUTION
    # =========================================================================

    print()
    print("=" * 100)
    print("FACTOR FINGERPRINT RESOLUTION")
    print("=" * 100)

    for ratio in CLOSE_RATIOS:

        print()
        print(f"close ratio = {ratio:.2%}")
        print("-" * 100)

        resolved = 0
        unresolved = 0
        candidate_counts = []

        for idx, x in enumerate(all_results[ratio], 1):

            moduli = (
                x["r1"],
                x["r2"],
                x["r3"],
            )

            residue_indices = {}

            for r in moduli:
                d = defaultdict(list)

                for prime in FACTOR_PRIMES:
                    d[prime % r].append(prime)

                residue_indices[r] = d

            result = analyze_anchor(
                x["n"],
                x["p"],
                x["q"],
                moduli,
                FACTOR_PRIMES,
                residue_indices,
            )

            c = result["candidate_count"]
            candidate_counts.append(c)

            if result["unique"]:
                resolved += 1
            else:
                unresolved += 1

            if idx % 25 == 0:
                print(f"anchor {idx:3d}/{len(all_results[ratio])}")

        print()
        print(f"resolved              = {resolved}")
        print(f"unresolved            = {unresolved}")

        if candidate_counts:
            print(f"mean candidates       = {statistics.mean(candidate_counts):.6f}")
            print(f"median candidates     = {statistics.median(candidate_counts):.6f}")
            print(f"maximum candidates    = {max(candidate_counts)}")

            print()
            print("CANDIDATE COUNT DISTRIBUTION")

            counts = defaultdict(int)

            for c in candidate_counts:
                counts[c] += 1

            for c in sorted(counts):
                print(
                    f"candidates={c:4d} "
                    f"anchors={counts[c]:4d}"
                )

    # =========================================================================
    # RANDOM CLOSE-TRIPLE BASELINE
    # =========================================================================

    print()
    print("=" * 100)
    print("RANDOM CLOSE-TRIPLE BASELINE")
    print("=" * 100)

    for ratio in CLOSE_RATIOS:

        random_ratios = []
        random_gaps = []
        valid = 0

        for n, _, _ in anchors:

            for _ in range(RANDOM_TRIPLES_PER_ANCHOR):

                triple = random_close_prime_triple(
                    n,
                    primes,
                    ratio,
                    rng,
                )

                if triple is None:
                    continue

                r1, r2, r3, R = triple

                valid += 1
                random_ratios.append(R / n)
                random_gaps.append(n - R)

        print()
        print(f"close ratio = {ratio:.2%}")

        if random_ratios:
            print(f"valid random triples = {valid:,}")
            print(
                f"mean random R/n       = "
                f"{statistics.mean(random_ratios):.10f}"
            )
            print(
                f"maximum random R/n    = "
                f"{max(random_ratios):.10f}"
            )
            print(
                f"mean random gap       = "
                f"{statistics.mean(random_gaps):,.2f}"
            )

    # =========================================================================
    # FINAL CONCEPTUAL CHECK
    # =========================================================================

    print()
    print("=" * 100)
    print("CONCEPTUAL CHECK")
    print("=" * 100)

    print(
        """
For each anchor n=p*q we searched for three close primes

    r1 < r2 < r3

such that

    r1*r2*r3 < n

and selected the maximum-product triple.

The main quantity is

    R/n = (r1*r2*r3)/n.

A value close to 1 means the three-prime modulus product is close
to n from below.

The factor fingerprint is

    ((p mod r1, q mod r1),
     (p mod r2, q mod r2),
     (p mod r3, q mod r3)).

The experiment asks whether those three modular observations identify
the actual factor pair inside the allowed factor-prime population.

This is deliberately different from the earlier experiments involving
the fixed modulus M.

Here r1,r2,r3 are chosen independently of M.

The important comparison is therefore:

    MAXIMUM-PRODUCT CLOSE TRIPLE
            versus
    RANDOM CLOSE TRIPLE.

A second important measurement is whether uniqueness occurs because
R is close to n, rather than merely because the individual moduli
are large.
"""
    )

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()
