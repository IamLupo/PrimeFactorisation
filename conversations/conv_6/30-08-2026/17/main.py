#!/usr/bin/env python3

import random
from collections import defaultdict, Counter

# ================================================================================================
# CONFIG
# ================================================================================================

M = 111_546_435

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

MODULUS_MIN = 300
MODULUS_MAX = 3_000

CLOSE_RATIO = 0.20

ACTUAL_ANCHORS = 300
RANDOM_CONTROLS_PER_TRIPLE = 20

SEED = 1_511_464_998

PROGRESS_EVERY = 25


# ================================================================================================
# PRIME TEST
# ================================================================================================

def is_prime(n):
    if n < 2:
        return False

    if n % 2 == 0:
        return n == 2

    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2

    return True


# ================================================================================================
# PRIME POOLS
# ================================================================================================

def build_primes(lo, hi):
    return [n for n in range(lo, hi + 1) if is_prime(n)]


# ================================================================================================
# ACTUAL ANCHORS
# ================================================================================================

def build_actual_anchors(primes, count, rng):
    result = []
    seen = set()

    while len(result) < count:
        p = rng.choice(primes)
        q = rng.choice(primes)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        pair = (p, q)

        if pair in seen:
            continue

        seen.add(pair)
        result.append(pair)

    return result


# ================================================================================================
# MODULUS PRIMES
# ================================================================================================

def build_modulus_primes(lo, hi):
    return build_primes(lo, hi)


# ================================================================================================
# MAX PRODUCT CLOSE TRIPLE
# ================================================================================================

def choose_max_product_triple(n, modulus_primes, close_ratio):
    """
    Find r1 < r2 < r3 such that

        r1*r2*r3 < n

    and

        r3/r1 <= 1 + close_ratio.

    Among all such triples choose the one with maximum product.
    """

    best = None
    best_product = -1

    L = len(modulus_primes)

    for i in range(L):
        r1 = modulus_primes[i]

        max_r3 = int(r1 * (1.0 + close_ratio))

        for j in range(i + 1, L):
            r2 = modulus_primes[j]

            if r2 > max_r3:
                break

            partial = r1 * r2

            if partial >= n:
                break

            for k in range(j + 1, L):
                r3 = modulus_primes[k]

                if r3 > max_r3:
                    break

                R = partial * r3

                if R >= n:
                    break

                if R > best_product:
                    best_product = R
                    best = (r1, r2, r3)

    return best


# ================================================================================================
# FINGERPRINT
# ================================================================================================

def residues(p, mods):
    return (
        p % mods[0],
        p % mods[1],
        p % mods[2],
    )


def fingerprint_raw(p, mods):
    return residues(p, mods)


def fingerprint_sum(p, mods):
    a, b, c = residues(p, mods)
    return a + b + c


def fingerprint_product(p, mods):
    a, b, c = residues(p, mods)
    return a * b * c


def fingerprint_sum_product(p, mods):
    a, b, c = residues(p, mods)

    e1 = a + b + c
    e3 = a * b * c

    return e1, e3


def fingerprint_symmetric(p, mods):
    a, b, c = residues(p, mods)

    e1 = a + b + c
    e2 = a * b + a * c + b * c
    e3 = a * b * c

    return e1, e2, e3


# ================================================================================================
# BUILD POPULATION MAPS
# ================================================================================================

def build_maps(primes, mods):
    raw = defaultdict(list)
    sums = defaultdict(list)
    products = defaultdict(list)
    sumprod = defaultdict(list)
    symmetric = defaultdict(list)

    for p in primes:

        raw[fingerprint_raw(p, mods)].append(p)

        sums[fingerprint_sum(p, mods)].append(p)

        products[fingerprint_product(p, mods)].append(p)

        sumprod[fingerprint_sum_product(p, mods)].append(p)

        symmetric[fingerprint_symmetric(p, mods)].append(p)

    return {
        "raw": raw,
        "sum": sums,
        "product": products,
        "sumprod": sumprod,
        "symmetric": symmetric,
    }


# ================================================================================================
# BUCKET STATISTICS
# ================================================================================================

def bucket_size(mapping, p, mods):
    if mapping is None:
        raise ValueError("mapping is None")

    key = {
        "raw": fingerprint_raw,
        "sum": fingerprint_sum,
        "product": fingerprint_product,
        "sumprod": fingerprint_sum_product,
        "symmetric": fingerprint_symmetric,
    }

    raise RuntimeError("Use bucket_size_for_type()")


def bucket_size_for_type(mapping, p, mods, typ):
    if typ == "raw":
        key = fingerprint_raw(p, mods)
    elif typ == "sum":
        key = fingerprint_sum(p, mods)
    elif typ == "product":
        key = fingerprint_product(p, mods)
    elif typ == "sumprod":
        key = fingerprint_sum_product(p, mods)
    elif typ == "symmetric":
        key = fingerprint_symmetric(p, mods)
    else:
        raise ValueError(typ)

    return len(mapping[key])


# ================================================================================================
# RANDOM CONTROL SAMPLING
# ================================================================================================

def choose_random_controls(primes, excluded, count, rng):
    """
    Sample distinct primes, excluding the actual p and q.
    """

    result = []

    while len(result) < count:
        p = rng.choice(primes)

        if p in excluded:
            continue

        result.append(p)

    return result


# ================================================================================================
# MAIN
# ================================================================================================

def run():

    rng = random.Random(SEED)

    print("=" * 100)
    print("THREE-MODULUS SUM+PRODUCT FINGERPRINT VALIDATION EXPERIMENT")
    print("=" * 100)

    print(f"M                         = {M:,}")
    print(f"factor range             = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"actual anchors            = {ACTUAL_ANCHORS}")
    print(f"modulus prime range      = {MODULUS_MIN:,} - {MODULUS_MAX:,}")
    print(f"close ratio              = {CLOSE_RATIO:.0%}")
    print(f"random controls/triple   = {RANDOM_CONTROLS_PER_TRIPLE}")
    print(f"seed                      = {SEED:,}")

    # --------------------------------------------------------------------------------------------
    # PRIME POOLS
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("BUILDING PRIME POOLS")
    print("=" * 100)

    factor_primes = build_primes(
        FACTOR_MIN,
        FACTOR_MAX
    )

    modulus_primes = build_modulus_primes(
        MODULUS_MIN,
        MODULUS_MAX
    )

    factor_set = set(factor_primes)

    print(f"factor primes             = {len(factor_primes):,}")
    print(f"modulus primes            = {len(modulus_primes):,}")

    # --------------------------------------------------------------------------------------------
    # ACTUAL ANCHORS
    # --------------------------------------------------------------------------------------------

    anchors = build_actual_anchors(
        factor_primes,
        ACTUAL_ANCHORS,
        rng
    )

    print(f"actual anchors             = {len(anchors):,}")

    # --------------------------------------------------------------------------------------------
    # SELECT MODULUS TRIPLE FOR EACH ANCHOR
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("SELECTING MAXIMUM-PRODUCT CLOSE TRIPLES")
    print("=" * 100)

    selected = []

    for i, (p, q) in enumerate(anchors, 1):

        n = p * q

        mods = choose_max_product_triple(
            n,
            modulus_primes,
            CLOSE_RATIO
        )

        if mods is None:
            continue

        selected.append(
            (p, q, n, mods)
        )

        if i % PROGRESS_EVERY == 0:
            print(f"anchor {i:3d}/{len(anchors)}")

    print()
    print(f"usable anchors = {len(selected):,}")

    # --------------------------------------------------------------------------------------------
    # FINGERPRINT MAPS
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("BUILDING FULL PRIME-POPULATION FINGERPRINT MAPS")
    print("=" * 100)

    maps_by_mods = {}

    for i, (_, _, _, mods) in enumerate(selected, 1):

        if mods not in maps_by_mods:

            maps_by_mods[mods] = build_maps(
                factor_primes,
                mods
            )

        if i % PROGRESS_EVERY == 0:
            print(f"triple {i:3d}/{len(selected)}")

    print()
    print(
        f"unique modulus triples = {len(maps_by_mods):,}"
    )

    # --------------------------------------------------------------------------------------------
    # STATISTICS ACCUMULATORS
    # --------------------------------------------------------------------------------------------

    TYPES = [
        "raw",
        "sum",
        "product",
        "sumprod",
        "symmetric",
    ]

    actual_bucket_sizes = {
        t: []
        for t in TYPES
    }

    control_bucket_sizes = {
        t: []
        for t in TYPES
    }

    actual_unique = Counter()
    control_unique = Counter()

    actual_ambiguous = Counter()
    control_ambiguous = Counter()

    actual_max_bucket = Counter()
    control_max_bucket = Counter()

    # Random control bucket quantiles.
    control_distributions = {
        t: []
        for t in TYPES
    }

    # --------------------------------------------------------------------------------------------
    # ANALYSIS
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("ANALYZING ACTUAL + RANDOM PRIME POPULATION")
    print("=" * 100)

    anchor_rows = []

    for i, (p, q, n, mods) in enumerate(selected, 1):

        maps = maps_by_mods[mods]

        row = {
            "p": p,
            "q": q,
            "n": n,
            "mods": mods,
            "R": mods[0] * mods[1] * mods[2],
        }

        # ------------------------------------------------------------------------
        # Actual p
        # ------------------------------------------------------------------------

        actual_sizes = {}

        for typ in TYPES:

            size = bucket_size_for_type(
                maps[typ],
                p,
                mods,
                typ
            )

            actual_sizes[typ] = size
            actual_bucket_sizes[typ].append(size)

            actual_max_bucket[typ] = max(
                actual_max_bucket[typ],
                size
            )

            if size == 1:
                actual_unique[typ] += 1
            else:
                actual_ambiguous[typ] += 1

        row["actual_sizes"] = actual_sizes

        # ------------------------------------------------------------------------
        # Random controls
        # ------------------------------------------------------------------------

        controls = choose_random_controls(
            factor_primes,
            {p, q},
            RANDOM_CONTROLS_PER_TRIPLE,
            rng
        )

        row_control_sizes = {}

        for typ in TYPES:

            sizes = []

            for x in controls:

                size = bucket_size_for_type(
                    maps[typ],
                    x,
                    mods,
                    typ
                )

                sizes.append(size)

                control_bucket_sizes[typ].append(size)
                control_distributions[typ].append(size)

                control_max_bucket[typ] = max(
                    control_max_bucket[typ],
                    size
                )

                if size == 1:
                    control_unique[typ] += 1
                else:
                    control_ambiguous[typ] += 1

            row_control_sizes[typ] = sizes

        row["control_sizes"] = row_control_sizes

        anchor_rows.append(row)

        if i % PROGRESS_EVERY == 0:
            print(f"anchor {i:3d}/{len(selected)}")

    # --------------------------------------------------------------------------------------------
    # SUMMARY
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("FULL PRIME-POPULATION RESULTS")
    print("=" * 100)

    print(
        f"{'FINGERPRINT':28s}"
        f"{'ACT UNIQUE':>12s}"
        f"{'ACT AMBIG':>12s}"
        f"{'ACT MAX':>10s}"
        f"{'CTRL UNIQUE':>14s}"
        f"{'CTRL AMBIG':>13s}"
        f"{'CTRL MAX':>10s}"
    )

    print("-" * 105)

    for typ, label in [
        ("raw", "RAW (a1,a2,a3)"),
        ("sum", "SUM (e1)"),
        ("product", "PRODUCT (e3)"),
        ("sumprod", "SUM+PRODUCT (e1,e3)"),
        ("symmetric", "SYMMETRIC (e1,e2,e3)"),
    ]:

        print(
            f"{label:28s}"
            f"{actual_unique[typ]:12d}"
            f"{actual_ambiguous[typ]:12d}"
            f"{actual_max_bucket[typ]:10d}"
            f"{control_unique[typ]:14d}"
            f"{control_ambiguous[typ]:13d}"
            f"{control_max_bucket[typ]:10d}"
        )

    # --------------------------------------------------------------------------------------------
    # MEAN BUCKET SIZES
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("MEAN BUCKET SIZE")
    print("=" * 100)

    for typ, label in [
        ("raw", "RAW"),
        ("sum", "SUM"),
        ("product", "PRODUCT"),
        ("sumprod", "SUM+PRODUCT"),
        ("symmetric", "SYMMETRIC"),
    ]:

        actual_mean = (
            sum(actual_bucket_sizes[typ])
            / len(actual_bucket_sizes[typ])
        )

        control_mean = (
            sum(control_bucket_sizes[typ])
            / len(control_bucket_sizes[typ])
        )

        print(
            f"{label:15s}"
            f" actual={actual_mean:.6f}"
            f" control={control_mean:.6f}"
        )

    # --------------------------------------------------------------------------------------------
    # MEDIANS
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("MEDIAN BUCKET SIZE")
    print("=" * 100)

    for typ, label in [
        ("raw", "RAW"),
        ("sum", "SUM"),
        ("product", "PRODUCT"),
        ("sumprod", "SUM+PRODUCT"),
        ("symmetric", "SYMMETRIC"),
    ]:

        a = sorted(actual_bucket_sizes[typ])
        c = sorted(control_bucket_sizes[typ])

        am = a[len(a) // 2]
        cm = c[len(c) // 2]

        print(
            f"{label:15s}"
            f" actual={am}"
            f" control={cm}"
        )

    # --------------------------------------------------------------------------------------------
    # PERCENTILES OF CONTROL BUCKETS
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("RANDOM CONTROL BUCKET DISTRIBUTION")
    print("=" * 100)

    for typ, label in [
        ("raw", "RAW"),
        ("sum", "SUM"),
        ("product", "PRODUCT"),
        ("sumprod", "SUM+PRODUCT"),
        ("symmetric", "SYMMETRIC"),
    ]:

        values = sorted(control_distributions[typ])

        def percentile(frac):
            if not values:
                return 0

            index = int(frac * (len(values) - 1))
            return values[index]

        print()
        print(label)

        print(
            f"  p50  = {percentile(0.50):.3f}"
        )

        print(
            f"  p90  = {percentile(0.90):.3f}"
        )

        print(
            f"  p95  = {percentile(0.95):.3f}"
        )

        print(
            f"  p99  = {percentile(0.99):.3f}"
        )

        print(
            f"  max  = {percentile(1.00):.3f}"
        )

    # --------------------------------------------------------------------------------------------
    # ACTUAL p AGAINST RANDOM CONTROL DISTRIBUTION
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("ACTUAL p POSITION INSIDE RANDOM-PRIME BUCKET DISTRIBUTION")
    print("=" * 100)

    for typ, label in [
        ("raw", "RAW"),
        ("sum", "SUM"),
        ("product", "PRODUCT"),
        ("sumprod", "SUM+PRODUCT"),
        ("symmetric", "SYMMETRIC"),
    ]:

        controls = sorted(control_distributions[typ])

        count_at_or_below = 0
        total = len(controls)

        values = actual_bucket_sizes[typ]

        percentile_values = []

        for x in values:

            le = 0

            for c in controls:
                if c <= x:
                    le += 1

            percentile_values.append(
                le / total
            )

        mean_percentile = (
            sum(percentile_values)
            / len(percentile_values)
        )

        print(
            f"{label:15s}"
            f" mean percentile={mean_percentile:.6f}"
        )

    # --------------------------------------------------------------------------------------------
    # DIRECT SUM+PRODUCT COLLISION ANALYSIS
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("SUM+PRODUCT COLLISION ANALYSIS")
    print("=" * 100)

    collision_triples = 0
    total_collision_members = 0
    largest_collision = 0

    examples = []

    for mods, maps in maps_by_mods.items():

        mapping = maps["sumprod"]

        local_collisions = [
            (key, values)
            for key, values in mapping.items()
            if len(values) > 1
        ]

        if local_collisions:

            collision_triples += 1

            for key, values in local_collisions:

                total_collision_members += len(values)

                largest_collision = max(
                    largest_collision,
                    len(values)
                )

                if len(examples) < 20:
                    examples.append(
                        (mods, key, values)
                    )

    print(
        f"triples containing collisions = "
        f"{collision_triples}/{len(maps_by_mods)}"
    )

    print(
        f"prime members in collision buckets = "
        f"{total_collision_members:,}"
    )

    print(
        f"largest collision bucket = "
        f"{largest_collision:,}"
    )

    print()
    print("Examples:")

    if not examples:
        print("none")
    else:
        for mods, key, values in examples:

            print()
            print(f"mods={mods}")
            print(f"SUM+PRODUCT={key}")

            print(
                "primes="
                + ", ".join(
                    f"{x:,}" for x in values[:20]
                )
            )

    # --------------------------------------------------------------------------------------------
    # SUM+PRODUCT VS RAW
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("INFORMATION RETENTION: SUM+PRODUCT VS RAW")
    print("=" * 100)

    raw_loss = 0
    sumprod_loss = 0
    both_unique = 0

    for row in anchor_rows:

        raw_size = row["actual_sizes"]["raw"]
        sp_size = row["actual_sizes"]["sumprod"]

        if raw_size > 1:
            raw_loss += 1

        if sp_size > 1:
            sumprod_loss += 1

        if raw_size == 1 and sp_size == 1:
            both_unique += 1

    print(
        f"actual anchors with RAW collision       = "
        f"{raw_loss}/{len(anchor_rows)}"
    )

    print(
        f"actual anchors with SUM+PRODUCT collision = "
        f"{sumprod_loss}/{len(anchor_rows)}"
    )

    print(
        f"actual anchors unique under BOTH         = "
        f"{both_unique}/{len(anchor_rows)}"
    )

    # --------------------------------------------------------------------------------------------
    # EXAMPLES OF ACTUALS
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("ACTUAL ANCHOR EXAMPLES")
    print("=" * 100)

    for row in anchor_rows[:25]:

        mods = row["mods"]

        print(
            f"p={row['p']:,} "
            f"q={row['q']:,} "
            f"mods={mods} "
            f"R/n={row['R']/row['n']:.9f} "
            f"raw={row['actual_sizes']['raw']} "
            f"sum={row['actual_sizes']['sum']} "
            f"product={row['actual_sizes']['product']} "
            f"sumprod={row['actual_sizes']['sumprod']} "
            f"sym={row['actual_sizes']['symmetric']}"
        )

    # --------------------------------------------------------------------------------------------
    # EXTREME ACTUAL SUM+PRODUCT BUCKETS
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("LARGEST ACTUAL SUM+PRODUCT BUCKETS")
    print("=" * 100)

    sorted_rows = sorted(
        anchor_rows,
        key=lambda row: row["actual_sizes"]["sumprod"],
        reverse=True
    )

    for row in sorted_rows[:20]:

        print(
            f"p={row['p']:,} "
            f"q={row['q']:,} "
            f"mods={row['mods']} "
            f"sumprod_bucket={row['actual_sizes']['sumprod']}"
        )

        mods = row["mods"]

        key = fingerprint_sum_product(
            row["p"],
            mods
        )

        values = maps_by_mods[mods]["sumprod"][key]

        print(
            "  matching primes = "
            + ", ".join(
                f"{x:,}" for x in values[:20]
            )
        )

    # --------------------------------------------------------------------------------------------
    # RATIO / UNIQUENESS
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("R/n VS SUM+PRODUCT UNIQUENESS")
    print("=" * 100)

    buckets = [
        (0.9999, []),
        (0.9990, []),
        (0.9950, []),
        (0.9900, []),
        (0.9500, []),
        (0.9000, []),
    ]

    for row in anchor_rows:

        ratio = row["R"] / row["n"]

        for threshold, values in buckets:

            if ratio >= threshold:
                values.append(
                    row["actual_sizes"]["sumprod"]
                )

    for threshold, values in buckets:

        if not values:
            continue

        unique = sum(
            1 for x in values
            if x == 1
        )

        print()
        print(
            f"R/n >= {threshold:.4f}"
        )

        print(
            f"  anchors = {len(values)}"
        )

        print(
            f"  unique SUM+PRODUCT = "
            f"{unique}/{len(values)}"
        )

        print(
            f"  mean bucket = "
            f"{sum(values)/len(values):.6f}"
        )

        print(
            f"  maximum bucket = "
            f"{max(values)}"
        )

    # --------------------------------------------------------------------------------------------
    # FINAL
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("FINAL")
    print("=" * 100)

    print(
        """
The critical statistic is the SUM+PRODUCT bucket.

For an actual prime p:

    a1 = p mod r1
    a2 = p mod r2
    a3 = p mod r3

    S = a1 + a2 + a3
    P = a1*a2*a3

The pair

    (S,P)

is considered unique only when exactly one prime in the COMPLETE
10,000..100,000 prime population produces that same pair.

This is stronger than checking the 300 actual anchors against one
another.

The comparison against random control primes answers a second question:

    Does an actual anchor's p have an unusually small SUM+PRODUCT bucket?

Possible outcomes:

1. SUM+PRODUCT is unique for nearly every prime.

   Then the fingerprint is a powerful encoding of primes in this
   finite population.

2. SUM+PRODUCT is usually ambiguous for random primes but unique for
   actual anchors.

   That would be much more interesting and would suggest that the
   selected maximum-product triples interact with the actual factors
   in a nontrivial way.

3. SUM+PRODUCT behaves almost identically to random controls.

   Then the uniqueness is mostly explained by the size of the CRT
   fingerprint space.

4. SUM alone or PRODUCT alone performs almost as well as SUM+PRODUCT.

   Then the pair (S,P) is not where the extra information resides.

5. SUM+PRODUCT is unique but factor-pair recovery is not.

   Then the fingerprint identifies p but does not by itself identify q.

The most important result to watch is therefore NOT merely:

    300/300 unique.

It is:

    actual SUM+PRODUCT bucket distribution
        versus
    random-prime SUM+PRODUCT bucket distribution.

That tells us whether the phenomenon is special to the factors we are
studying or is simply a consequence of the large three-modulus
fingerprint space.
"""
    )

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()
