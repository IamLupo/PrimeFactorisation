#!/usr/bin/env python3

import math
import random
from collections import Counter, defaultdict

# =============================================================================
# THREE-MODULUS RESIDUE-SUM / RESIDUE-PRODUCT / SYMMETRIC-POLYNOMIAL
# / QUOTIENT-COORDINATE FINGERPRINT EXPERIMENT
#
# For each actual semiprime n = p*q:
#
#   1. Select r1 < r2 < r3 as close primes with maximum product R < n.
#   2. Compute:
#
#        a1 = p mod r1
#        a2 = p mod r2
#        a3 = p mod r3
#
#   3. Compare several p-fingerprints:
#
#        RAW
#            (a1, a2, a3)
#
#        SUM/PRODUCT
#            (e1, e3)
#
#        FULL SYMMETRIC
#            (e1, e2, e3)
#
#        where
#            e1 = a1 + a2 + a3
#            e2 = a1*a2 + a1*a3 + a2*a3
#            e3 = a1*a2*a3
#
#        QUOTIENT COORDINATES
#            k1 = (p-a1)/r1
#            k2 = (p-a2)/r2
#            k3 = (p-a3)/r3
#
#   4. For each fingerprint, determine how many primes in the complete
#      factor pool share that fingerprint.
#
#   5. Repeat the same analysis for q.
#
#   6. Test whether the fingerprints can identify the actual p/q pair.
#
# This does NOT assume that the compressed quantities preserve all CRT
# information. It explicitly measures how much information is lost.
#
# IMPORTANT:
# Uniqueness is only conditional on:
#   - this factor interval
#   - this prime pool
#   - the selected three moduli
#
# =============================================================================


# =============================================================================
# CONFIGURATION
# =============================================================================

M = 111_546_435

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

ACTUAL_TRIALS = 300

MODULUS_MIN = 300
MODULUS_MAX = 3000

CLOSE_RATIO = 0.20

SEED = 1_511_464_998

# =============================================================================
# PRIME GENERATION
# =============================================================================


def sieve(limit):
    """Return all primes <= limit."""
    if limit < 2:
        return []

    is_prime = bytearray(b"\x01") * (limit + 1)
    is_prime[0:2] = b"\x00\x00"

    root = int(math.isqrt(limit))

    for p in range(2, root + 1):
        if is_prime[p]:
            start = p * p
            is_prime[start:limit + 1:p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    return [i for i in range(2, limit + 1) if is_prime[i]]


# =============================================================================
# MODULUS TRIPLE SELECTION
# =============================================================================


def build_close_triples():
    """
    Build all useful triples of modulus primes.

    A triple is accepted when:

        max(r_i) / min(r_i) <= 1 + CLOSE_RATIO

    We later select the triple whose product is the largest value < n.
    """
    primes = sieve(MODULUS_MAX)
    primes = [p for p in primes if p >= MODULUS_MIN]

    triples = []

    # The modulus range is small enough that this remains cheap.
    for i in range(len(primes)):
        r1 = primes[i]

        # Since r2/r1 and r3/r1 must remain inside the close ratio,
        # stop as soon as the ratio is too large.
        max_value = int(r1 * (1.0 + CLOSE_RATIO))

        j = i + 1

        while j < len(primes) and primes[j] <= max_value:
            r2 = primes[j]

            k = j + 1
            while k < len(primes) and primes[k] <= max_value:
                r3 = primes[k]

                R = r1 * r2 * r3

                triples.append((r1, r2, r3, R))
                k += 1

            j += 1

    return triples


def select_best_triple(n, triples):
    """
    Select the admissible triple with maximum product R < n.

    Returns:
        (r1, r2, r3, R)
    """
    best = None
    best_R = -1

    for r1, r2, r3, R in triples:
        if R < n and R > best_R:
            best_R = R
            best = (r1, r2, r3, R)

    return best


# =============================================================================
# FACTOR POOL
# =============================================================================


def build_factor_pool():
    return sieve(FACTOR_MAX)


# =============================================================================
# ACTUAL ANCHORS
# =============================================================================


def generate_actual_anchors(primes, count, rng):
    """
    Generate distinct semiprime anchors p*q with:

        FACTOR_MIN <= p,q <= FACTOR_MAX

    and p != q.

    We keep p < q for canonical representation.
    """
    anchors = []
    seen = set()

    while len(anchors) < count:
        p = rng.choice(primes)
        q = rng.choice(primes)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        key = (p, q)

        if key in seen:
            continue

        seen.add(key)
        anchors.append((p, q, p * q))

    return anchors


# =============================================================================
# FINGERPRINT FUNCTIONS
# =============================================================================


def residue_data(p, r1, r2, r3):
    a1 = p % r1
    a2 = p % r2
    a3 = p % r3

    e1 = a1 + a2 + a3
    e2 = a1 * a2 + a1 * a3 + a2 * a3
    e3 = a1 * a2 * a3

    k1 = (p - a1) // r1
    k2 = (p - a2) // r2
    k3 = (p - a3) // r3

    return {
        "raw": (a1, a2, a3),
        "sum_product": (e1, e3),
        "symmetric": (e1, e2, e3),
        "quotients": (k1, k2, k3),
        "a": (a1, a2, a3),
        "e1": e1,
        "e2": e2,
        "e3": e3,
        "k": (k1, k2, k3),
    }


# =============================================================================
# PRIME-POOL INDEX
# =============================================================================


def build_fingerprint_indices(primes, r1, r2, r3):
    """
    Build all four fingerprint indexes over the complete prime pool.

    This is O(number_of_primes), not O(number_of_prime_pairs).
    """

    indices = {
        "raw": defaultdict(list),
        "sum_product": defaultdict(list),
        "symmetric": defaultdict(list),
        "quotients": defaultdict(list),
    }

    for p in primes:
        d = residue_data(p, r1, r2, r3)

        indices["raw"][d["raw"]].append(p)
        indices["sum_product"][d["sum_product"]].append(p)
        indices["symmetric"][d["symmetric"]].append(p)
        indices["quotients"][d["quotients"]].append(p)

    return indices


# =============================================================================
# BASIC STATS
# =============================================================================


def fingerprint_bucket_stats(index):
    sizes = [len(v) for v in index.values()]

    if not sizes:
        return {
            "signatures": 0,
            "unique": 0,
            "max_bucket": 0,
            "mean_bucket": 0.0,
        }

    return {
        "signatures": len(index),
        "unique": sum(1 for x in sizes if x == 1),
        "max_bucket": max(sizes),
        "mean_bucket": sum(sizes) / len(sizes),
    }


def percentile_rank(bucket_size, index):
    """
    Percentile-like measure for an actual bucket size:
    fraction of prime fingerprints with bucket size >= bucket_size.
    """
    sizes = [len(v) for v in index.values()]

    if not sizes:
        return 0.0

    return sum(1 for x in sizes if x >= bucket_size) / len(sizes)


# =============================================================================
# PAIR CANDIDATE SEARCH
# =============================================================================


def recover_pairs_from_raw(n, primes, r1, r2, r3):
    """
    Search using the raw p fingerprint.

    This is intentionally efficient:
    only primes p matching the actual p fingerprint need to be tested.

    Since:

        p*q == n (mod r_i)

    we compute the required q residue modulo each modulus.
    """
    pd = residue_data(
        recover_pairs_from_raw.current_p,
        r1,
        r2,
        r3
    )

    _ = pd  # kept explicit for readability

    raise RuntimeError("Internal misuse")


def recover_factor_pairs(
    n,
    actual_p,
    actual_q,
    primes,
    r1,
    r2,
    r3,
    p_index
):
    """
    Recover candidate factor pairs by using the actual p fingerprint.

    For each prime p' with the same raw p fingerprint as actual_p,
    q' must satisfy:

        p' * q' == n (mod r_i)

    for all three moduli.

    We then test q' against the prime pool.
    """

    target = residue_data(actual_p, r1, r2, r3)

    p_candidates = p_index["raw"].get(target["raw"], [])

    prime_set = set(primes)

    candidates = []

    for pc in p_candidates:
        # Compute required q residues modulo each modulus.
        # Since pc is prime and r_i are different modulus primes,
        # pc is either invertible modulo r_i or exactly equal to r_i.
        #
        # In our selected range r_i can be inside factor range,
        # so explicitly handle the divisible case.

        residues = []

        valid = True

        for r in (r1, r2, r3):
            pc_mod = pc % r
            n_mod = n % r

            if pc_mod == 0:
                if n_mod != 0:
                    valid = False
                    break

                residues.append(None)
            else:
                inv = pow(pc_mod, -1, r)
                residues.append((n_mod * inv) % r)

        if not valid:
            continue

        # Test every q in the same factor pool for the required residues.
        # Because the candidate p bucket is already tiny, this remains fast.
        for qc in primes:

            if qc == pc:
                continue

            ok = True

            for r, want in zip((r1, r2, r3), residues):
                if want is not None:
                    if qc % r != want:
                        ok = False
                        break
                else:
                    if qc % r != 0:
                        ok = False
                        break

            if not ok:
                continue

            if pc * qc == n:
                pair = tuple(sorted((pc, qc)))
                if pair not in candidates:
                    candidates.append(pair)

    actual_pair = tuple(sorted((actual_p, actual_q)))

    return candidates, actual_pair


# =============================================================================
# MAIN EXPERIMENT
# =============================================================================


def run():
    rng = random.Random(SEED)

    print("=" * 100)
    print("THREE-MODULUS RESIDUE-SUM / RESIDUE-PRODUCT / SYMMETRIC")
    print("/ QUOTIENT-COORDINATE FINGERPRINT EXPERIMENT")
    print("=" * 100)

    print(f"M                         = {M:,}")
    print(f"factor range             = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"actual anchors            = {ACTUAL_TRIALS}")
    print(f"modulus prime range       = {MODULUS_MIN:,} - {MODULUS_MAX:,}")
    print(f"close ratio               = {CLOSE_RATIO:.0%}")
    print(f"seed                      = {SEED:,}")

    # -------------------------------------------------------------------------
    # PRIME POOLS
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("BUILDING PRIME POOLS")
    print("=" * 100)

    all_primes = build_factor_pool()
    factor_primes = [
        p for p in all_primes
        if FACTOR_MIN <= p <= FACTOR_MAX
    ]

    print(f"factor primes             = {len(factor_primes):,}")

    modulus_triples = build_close_triples()

    print(f"admissible modulus triples = {len(modulus_triples):,}")

    # -------------------------------------------------------------------------
    # ACTUAL ANCHORS
    # -------------------------------------------------------------------------

    anchors = generate_actual_anchors(
        factor_primes,
        ACTUAL_TRIALS,
        rng
    )

    print(f"actual anchors             = {len(anchors):,}")

    # -------------------------------------------------------------------------
    # ANALYZE EACH ANCHOR
    # -------------------------------------------------------------------------

    results = []

    for i, (p, q, n) in enumerate(anchors, 1):

        triple = select_best_triple(n, modulus_triples)

        if triple is None:
            continue

        r1, r2, r3, R = triple

        p_data = residue_data(p, r1, r2, r3)
        q_data = residue_data(q, r1, r2, r3)

        # Build a prime index specifically for this triple.
        #
        # This is O(8363) and is cheap.
        # It replaces any O(8363^2) search.
        indices = build_fingerprint_indices(
            factor_primes,
            r1,
            r2,
            r3
        )

        # -------------------------------------------------------------
        # p fingerprint bucket sizes
        # -------------------------------------------------------------

        p_bucket_raw = indices["raw"].get(
            p_data["raw"], []
        )

        p_bucket_sum_product = indices["sum_product"].get(
            p_data["sum_product"], []
        )

        p_bucket_symmetric = indices["symmetric"].get(
            p_data["symmetric"], []
        )

        p_bucket_quotients = indices["quotients"].get(
            p_data["quotients"], []
        )

        # -------------------------------------------------------------
        # q fingerprint bucket sizes
        # -------------------------------------------------------------

        q_bucket_raw = indices["raw"].get(
            q_data["raw"], []
        )

        q_bucket_sum_product = indices["sum_product"].get(
            q_data["sum_product"], []
        )

        q_bucket_symmetric = indices["symmetric"].get(
            q_data["symmetric"], []
        )

        q_bucket_quotients = indices["quotients"].get(
            q_data["quotients"], []
        )

        # -------------------------------------------------------------
        # Pair recovery from raw fingerprint
        # -------------------------------------------------------------

        pair_candidates, actual_pair = recover_factor_pairs(
            n,
            p,
            q,
            factor_primes,
            r1,
            r2,
            r3,
            indices
        )

        results.append({
            "p": p,
            "q": q,
            "n": n,
            "r1": r1,
            "r2": r2,
            "r3": r3,
            "R": R,
            "R_ratio": R / n,

            "p_data": p_data,
            "q_data": q_data,

            "p_raw_size": len(p_bucket_raw),
            "p_sum_product_size": len(p_bucket_sum_product),
            "p_symmetric_size": len(p_bucket_symmetric),
            "p_quotients_size": len(p_bucket_quotients),

            "q_raw_size": len(q_bucket_raw),
            "q_sum_product_size": len(q_bucket_sum_product),
            "q_symmetric_size": len(q_bucket_symmetric),
            "q_quotients_size": len(q_bucket_quotients),

            "pair_candidates": pair_candidates,
            "pair_count": len(pair_candidates),
        })

        if i % 25 == 0 or i == len(anchors):
            print(f"anchor {i:3d}/{len(anchors)}")

    # -------------------------------------------------------------------------
    # ASSERTIONS / SANITY
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("SANITY CHECK")
    print("=" * 100)

    identity_failures = 0
    quotient_failures = 0

    for x in results:
        p = x["p"]
        r1, r2, r3 = x["r1"], x["r2"], x["r3"]

        d = x["p_data"]

        if p % r1 != d["a"][0]:
            identity_failures += 1

        if p % r2 != d["a"][1]:
            identity_failures += 1

        if p % r3 != d["a"][2]:
            identity_failures += 1

        for r, k, a in zip(
            (r1, r2, r3),
            d["k"],
            d["a"]
        ):
            if p != k * r + a:
                quotient_failures += 1

    print(f"residue identity failures   = {identity_failures}")
    print(f"quotient identity failures  = {quotient_failures}")

    # -------------------------------------------------------------------------
    # AGGREGATES
    # -------------------------------------------------------------------------

    def mean(field):
        return sum(x[field] for x in results) / len(results)

    def resolved(field):
        return sum(
            1 for x in results
            if x[field] == 1
        )

    # -------------------------------------------------------------------------
    # P FINGERPRINT RESULTS
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("P-FINGERPRINT UNIQUENESS")
    print("=" * 100)

    print(
        f"{'FINGERPRINT':<24}"
        f"{'MEAN BUCKET':>16}"
        f"{'UNIQUE':>12}"
        f"{'AMBIGUOUS':>12}"
        f"{'MAX':>10}"
    )
    print("-" * 100)

    fingerprint_specs = [
        (
            "RAW (a1,a2,a3)",
            "p_raw_size",
        ),
        (
            "SUM + PRODUCT (e1,e3)",
            "p_sum_product_size",
        ),
        (
            "FULL SYMMETRIC (e1,e2,e3)",
            "p_symmetric_size",
        ),
        (
            "QUOTIENTS (k1,k2,k3)",
            "p_quotients_size",
        ),
    ]

    for name, field in fingerprint_specs:
        values = [x[field] for x in results]

        print(
            f"{name:<24}"
            f"{sum(values)/len(values):16.6f}"
            f"{sum(1 for v in values if v == 1):12d}"
            f"{sum(1 for v in values if v > 1):12d}"
            f"{max(values):10d}"
        )

    # -------------------------------------------------------------------------
    # Q FINGERPRINT RESULTS
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("Q-FINGERPRINT UNIQUENESS")
    print("=" * 100)

    print(
        f"{'FINGERPRINT':<24}"
        f"{'MEAN BUCKET':>16}"
        f"{'UNIQUE':>12}"
        f"{'AMBIGUOUS':>12}"
        f"{'MAX':>10}"
    )
    print("-" * 100)

    q_specs = [
        (
            "RAW (a1,a2,a3)",
            "q_raw_size",
        ),
        (
            "SUM + PRODUCT (e1,e3)",
            "q_sum_product_size",
        ),
        (
            "FULL SYMMETRIC (e1,e2,e3)",
            "q_symmetric_size",
        ),
        (
            "QUOTIENTS (k1,k2,k3)",
            "q_quotients_size",
        ),
    ]

    for name, field in q_specs:
        values = [x[field] for x in results]

        print(
            f"{name:<24}"
            f"{sum(values)/len(values):16.6f}"
            f"{sum(1 for v in values if v == 1):12d}"
            f"{sum(1 for v in values if v > 1):12d}"
            f"{max(values):10d}"
        )

    # -------------------------------------------------------------------------
    # PAIR RECOVERY
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("FACTOR-PAIR RECOVERY")
    print("=" * 100)

    pair_counts = Counter(
        x["pair_count"]
        for x in results
    )

    print(f"unique factor pair   = {pair_counts.get(1, 0):3d}")
    print(
        f"ambiguous factor pair = "
        f"{sum(v for k, v in pair_counts.items() if k > 1):3d}"
    )
    print(
        f"zero factor pair      = "
        f"{pair_counts.get(0, 0):3d}"
    )

    print()
    print("candidate distribution:")

    for count in sorted(pair_counts):
        print(
            f"    {count:2d} candidates : "
            f"{pair_counts[count]:3d} anchors"
        )

    # -------------------------------------------------------------------------
    # RATIO STATISTICS
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("SELECTED MODULUS TRIPLES")
    print("=" * 100)

    ratios = [x["R_ratio"] for x in results]

    print(f"mean R/n       = {sum(ratios)/len(ratios):.10f}")
    print(f"median R/n     = {sorted(ratios)[len(ratios)//2]:.10f}")
    print(f"minimum R/n    = {min(ratios):.10f}")
    print(f"maximum R/n    = {max(ratios):.10f}")

    # -------------------------------------------------------------------------
    # CORRELATION BETWEEN RATIO AND FINGERPRINT UNIQUENESS
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("R/n vs P-FINGERPRINT BUCKET SIZE")
    print("=" * 100)

    for threshold in (
        0.90,
        0.95,
        0.99,
        0.995,
        0.999,
    ):
        subset = [
            x for x in results
            if x["R_ratio"] >= threshold
        ]

        if not subset:
            continue

        print()
        print(f"R/n >= {threshold:.4f}")
        print(f"anchors = {len(subset)}")

        for label, field in (
            ("RAW", "p_raw_size"),
            ("SUM+PRODUCT", "p_sum_product_size"),
            ("SYMMETRIC", "p_symmetric_size"),
            ("QUOTIENTS", "p_quotients_size"),
        ):
            vals = [x[field] for x in subset]

            print(
                f"  {label:<12}"
                f" mean={sum(vals)/len(vals):.6f}"
                f" unique={sum(v == 1 for v in vals):3d}"
                f"/{len(vals)}"
            )

    # -------------------------------------------------------------------------
    # EXTREME EXAMPLES
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("EXAMPLES")
    print("=" * 100)

    print()
    print("MOST COMPRESSED SUM+PRODUCT COLLISIONS")

    examples = sorted(
        results,
        key=lambda x: (
            x["p_sum_product_size"],
            x["R_ratio"]
        ),
        reverse=True
    )

    shown = 0

    for x in examples:
        if x["p_sum_product_size"] > 1:
            print(
                f"p={x['p']:,} q={x['q']:,} n={x['n']:,} "
                f"mods=({x['r1']},{x['r2']},{x['r3']}) "
                f"S/P bucket={x['p_sum_product_size']} "
                f"raw={x['p_raw_size']} "
                f"sym={x['p_symmetric_size']} "
                f"quot={x['p_quotients_size']}"
            )

            print(
                f"    p residues={x['p_data']['a']}"
            )
            print(
                f"    (e1,e2,e3)="
                f"({x['p_data']['e1']},"
                f" {x['p_data']['e2']},"
                f" {x['p_data']['e3']})"
            )
            print(
                f"    k={x['p_data']['k']}"
            )

            shown += 1

            if shown >= 10:
                break

    if shown == 0:
        print("none")

    # -------------------------------------------------------------------------
    # TOP PAIR-AMBIGUOUS CASES
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("PAIR-AMBIGUOUS EXAMPLES")
    print("=" * 100)

    ambiguous = [
        x for x in results
        if x["pair_count"] > 1
    ]

    ambiguous.sort(
        key=lambda x: x["pair_count"],
        reverse=True
    )

    if not ambiguous:
        print("none")
    else:
        for x in ambiguous[:10]:
            print()
            print(
                f"n={x['n']:,} "
                f"actual=({x['p']:,},{x['q']:,}) "
                f"mods=({x['r1']},{x['r2']},{x['r3']})"
            )
            print(
                f"candidate pairs={x['pair_candidates']}"
            )

    # -------------------------------------------------------------------------
    # TOP MODULUS TRIPLES
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("SELECTED TRIPLE EXAMPLES")
    print("=" * 100)

    top = sorted(
        results,
        key=lambda x: x["R_ratio"],
        reverse=True
    )

    for x in top[:20]:
        print(
            f"n={x['n']:,} "
            f"mods=({x['r1']},{x['r2']},{x['r3']}) "
            f"R={x['R']:,} "
            f"R/n={x['R_ratio']:.10f} "
            f"raw={x['p_raw_size']} "
            f"sumprod={x['p_sum_product_size']} "
            f"sym={x['p_symmetric_size']} "
            f"quot={x['p_quotients_size']} "
            f"pairs={x['pair_count']}"
        )

    # -------------------------------------------------------------------------
    # INTERPRETATION
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("INTERPRETATION")
    print("=" * 100)

    print(
        """
The experiment compares four different ways of representing p.

1. RAW
       (p mod r1, p mod r2, p mod r3)

   This is the complete ordered three-residue fingerprint.

2. SUM + PRODUCT
       (e1, e3)

   where

       e1 = a1+a2+a3
       e3 = a1*a2*a3

   This is a compressed representation and can collide even when
   the original residue triple is unique.

3. FULL SYMMETRIC
       (e1,e2,e3)

   where

       e2 = a1*a2 + a1*a3 + a2*a3.

   These are the coefficients of

       (x-a1)(x-a2)(x-a3).

   Therefore they recover the multiset {a1,a2,a3}, but they do
   not preserve which residue belongs to which modulus.

4. QUOTIENT COORDINATES
       (k1,k2,k3)

   where

       p = k_i*r_i + a_i.

   These describe the position of p relative to the three moduli.

The critical comparison is:

       RAW bucket size
              versus
       SUM+PRODUCT bucket size
              versus
       SYMMETRIC bucket size
              versus
       QUOTIENT bucket size.

If SUM+PRODUCT remains almost always unique, then the two simple
quantities you proposed preserve a surprising amount of information.

If SUM+PRODUCT has many collisions but FULL SYMMETRIC is unique,
then the missing e2 term is carrying the lost information.

If the raw fingerprint is unique but all compressed forms collide,
then the information is genuinely tied to the individual residue
coordinates rather than their symmetric combinations.

The experiment also measures whether R/n is related to fingerprint
uniqueness.

A particularly interesting result would be:

    R close to n
        AND
    very small compressed fingerprint buckets.

That would indicate that the usefulness comes from the combination
of a large three-modulus CRT range and the algebraic compression.

No claim of general factorization follows from uniqueness inside
this finite prime interval. The experiment measures information
content in exactly the population being studied.
"""
    )

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()
