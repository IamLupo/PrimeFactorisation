#!/usr/bin/env python3

import math
import random
from collections import Counter

# ==================================================================================================
# EXACT C-CLASS ANCHOR-CENTRALITY EXPERIMENT
# ==================================================================================================

M = 111_546_435
TMAX = 25

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

TRIALS = 300
SEED = 1_511_464_998

# ==================================================================================================
# PRIME SIEVE
# ==================================================================================================

def sieve(n):
    s = bytearray(b"\x01") * (n + 1)
    s[0:2] = b"\x00\x00"

    for p in range(2, math.isqrt(n) + 1):
        if s[p]:
            start = p * p
            s[start:n + 1:p] = b"\x00" * (((n - start) // p) + 1)

    return [
        x for x in range(FACTOR_MIN, FACTOR_MAX + 1)
        if s[x]
    ]


PRIMES = sieve(FACTOR_MAX)
PRIME_SET = set(PRIMES)


# ==================================================================================================
# ANCHOR GENERATION
# ==================================================================================================

def generate_unique_anchors(count, seed):
    rng = random.Random(seed)

    anchors = []
    seen_c = set()

    while len(anchors) < count:

        p = rng.choice(PRIMES)
        q = rng.choice(PRIMES)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        C = (p * q) % M

        if C in seen_c:
            continue

        seen_c.add(C)
        anchors.append((p, q, p * q, C))

    return anchors


# ==================================================================================================
# EXACT C-CLASS
# ==================================================================================================

def build_c_class(C):
    candidates = set()

    for a in PRIMES:

        if math.gcd(a, M) != 1:
            continue

        inv = pow(a, -1, M)
        b = (C * inv) % M

        if not (FACTOR_MIN <= b <= FACTOR_MAX):
            continue

        if b not in PRIME_SET:
            continue

        if a > b:
            continue

        if (a * b) % M != C:
            continue

        candidates.add((a, b))

    return sorted(candidates)


# ==================================================================================================
# GEOMETRIC STATISTICS
# ==================================================================================================

def anchor_geometry(anchor, candidates):

    p, q = anchor

    others = [
        (a, b)
        for a, b in candidates
        if (a, b) != (p, q)
    ]

    if not others:
        return {
            "count": 0,

            "sum_abs_dp": 0.0,
            "sum_abs_dq": 0.0,
            "sum_abs_dsum": 0.0,
            "sum_abs_ddiff": 0.0,

            "sum_dist": 0.0,
            "min_dist": 0.0,
            "max_dist": 0.0,

            "sum_abs_t": 0.0,
            "mean_t": 0.0,
            "same_sign": 0,
            "opposite": 0,
        }

    abs_dp = []
    abs_dq = []
    abs_dsum = []
    abs_ddiff = []

    distances = []

    abs_t = []
    ts = []

    same_sign = 0
    opposite = 0

    for a, b in others:

        dp = a - p
        dq = b - q

        dsum = dp + dq
        ddiff = dp - dq

        delta = a * b - p * q

        if delta % (2 * M) != 0:
            raise RuntimeError("Non-integral t")

        t = delta // (2 * M)

        if abs(t) <= TMAX:
            abs_t.append(abs(t))
            ts.append(t)

        abs_dp.append(abs(dp))
        abs_dq.append(abs(dq))
        abs_dsum.append(abs(dsum))
        abs_ddiff.append(abs(ddiff))

        dist = math.hypot(dp, dq)

        distances.append(dist)

        if dp * dq > 0:
            same_sign += 1
        elif dp * dq < 0:
            opposite += 1

    return {
        "count": len(others),

        "sum_abs_dp": sum(abs_dp),
        "sum_abs_dq": sum(abs_dq),
        "sum_abs_dsum": sum(abs_dsum),
        "sum_abs_ddiff": sum(abs_ddiff),

        "sum_dist": sum(distances),
        "min_dist": min(distances),
        "max_dist": max(distances),

        "sum_abs_t": sum(abs_t),
        "mean_t": sum(ts) / len(ts) if ts else 0.0,

        "same_sign": same_sign,
        "opposite": opposite,
    }


# ==================================================================================================
# ADDITIVE STATISTICS
# ==================================================================================================

STAT_KEYS = [
    "sum_abs_dp",
    "sum_abs_dq",
    "sum_abs_dsum",
    "sum_abs_ddiff",
    "sum_dist",
    "sum_abs_t",
    "same_sign",
    "opposite",
]


# ==================================================================================================
# BUILD CLASSES
# ==================================================================================================

def build_classes(anchors):

    classes = {}

    for p, q, n, C in anchors:

        candidates = build_c_class(C)

        if (p, q) not in candidates:
            raise RuntimeError(
                f"Anchor missing from C-class: {(p, q)} C={C}"
            )

        classes[C] = candidates

    return classes


# ==================================================================================================
# ACTUAL STATISTICS
# ==================================================================================================

def actual_dataset(anchors, classes):

    rows = []

    for p, q, n, C in anchors:

        candidates = classes[C]

        g = anchor_geometry((p, q), candidates)

        rows.append({
            "p": p,
            "q": q,
            "n": n,
            "C": C,
            **g,
        })

    return rows


# ==================================================================================================
# EXACT RE-ANCHOR NULL
# ==================================================================================================

def null_dataset(classes):

    """
    Exact null:

        for each C-class,
        every candidate is equally eligible to become the anchor.

    There is NO Monte Carlo here.
    """

    per_class = {}

    for C, candidates in classes.items():

        vertex_stats = []

        for anchor in candidates:

            g = anchor_geometry(anchor, candidates)

            vertex_stats.append({
                "anchor": anchor,
                **g,
            })

        per_class[C] = vertex_stats

    return per_class


# ==================================================================================================
# AGGREGATE MEAN
# ==================================================================================================

def aggregate_mean(rows, key):

    if not rows:
        return 0.0

    return sum(float(r[key]) for r in rows) / len(rows)


# ==================================================================================================
# EXACT NULL MEAN / VARIANCE
# ==================================================================================================

def exact_null_moments(per_class, key):

    total = 0.0
    total_sq = 0.0
    choices = 0

    for vertex_stats in per_class.values():

        if not vertex_stats:
            continue

        for r in vertex_stats:

            x = float(r[key])

            total += x
            total_sq += x * x
            choices += 1

    if choices == 0:
        return 0.0, 0.0

    mean = total / choices

    var = max(
        0.0,
        total_sq / choices - mean * mean
    )

    return mean, var


# ==================================================================================================
# CLASS-CORRECT AGGREGATE NULL
# ==================================================================================================

def exact_class_null_moments(per_class, key):

    """
    Important:

    Every C-class contributes exactly ONE anchor.

    Therefore classes must be weighted equally,
    NOT by their number of candidate vertices.
    """

    class_means = []

    for vertex_stats in per_class.values():

        if not vertex_stats:
            continue

        mu = sum(
            float(r[key])
            for r in vertex_stats
        ) / len(vertex_stats)

        class_means.append(mu)

    if not class_means:
        return 0.0, 0.0

    mean = sum(class_means) / len(class_means)

    var = (
        sum((x - mean) ** 2 for x in class_means)
        / len(class_means)
    )

    return mean, var


# ==================================================================================================
# ACTUAL CLASS-LEVEL AGGREGATE
# ==================================================================================================

def actual_class_mean(rows, key):

    return aggregate_mean(rows, key)


# ==================================================================================================
# EMPIRICAL EXACT PERCENTILE
# ==================================================================================================

def exact_percentile_position(actual, values):

    """
    Mid-rank percentile.
    """

    if not values:
        return 0.5

    below = sum(v < actual for v in values)
    equal = sum(v == actual for v in values)

    return (below + 0.5 * equal) / len(values)


# ==================================================================================================
# WITHIN-CLASS RANKS
# ==================================================================================================

def anchor_ranks(anchors, classes):

    rank_rows = []

    for p, q, n, C in anchors:

        candidates = classes[C]

        tuples = []

        for candidate in candidates:

            a, b = candidate

            tuples.append({
                "anchor": candidate,
                "sum": a + b,
                "diff": abs(a - b),
                "product": a * b,
                "ratio": a / b,
            })

        original = (p, q)

        ours = next(x for x in tuples if x["anchor"] == original)

        rank_rows.append({
            "C": C,
            "anchor": original,
            "class_size": len(candidates),

            "sum_rank": 1 + sum(
                x["sum"] < ours["sum"]
                for x in tuples
            ),

            "diff_rank": 1 + sum(
                x["diff"] < ours["diff"]
                for x in tuples
            ),

            "product_rank": 1 + sum(
                x["product"] < ours["product"]
                for x in tuples
            ),

            "ratio_rank": 1 + sum(
                x["ratio"] < ours["ratio"]
                for x in tuples
            ),

            "sum_values": [x["sum"] for x in tuples],
            "diff_values": [x["diff"] for x in tuples],
            "product_values": [x["product"] for x in tuples],
            "ratio_values": [x["ratio"] for x in tuples],
        })

    return rank_rows


# ==================================================================================================
# MAIN
# ==================================================================================================

def main():

    print("=" * 100)
    print("EXACT C-CLASS ANCHOR-CENTRALITY EXPERIMENT")
    print("=" * 100)

    print(f"M                    = {M:,}")
    print(f"2M                   = {2*M:,}")
    print(f"factor range         = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"T range              = [-{TMAX}, +{TMAX}]")
    print(f"actual anchors       = {TRIALS}")
    print(f"factor primes        = {len(PRIMES):,}")
    print(f"seed                 = {SEED:,}")
    print()

    # ----------------------------------------------------------------------------------------------
    # Anchors
    # ----------------------------------------------------------------------------------------------

    print("=" * 100)
    print("GENERATING ACTUAL ANCHORS")
    print("=" * 100)

    anchors = generate_unique_anchors(TRIALS, SEED)

    print(f"anchors = {len(anchors)}")
    print()

    # ----------------------------------------------------------------------------------------------
    # Classes
    # ----------------------------------------------------------------------------------------------

    print("=" * 100)
    print("BUILDING EXACT C-CLASSES")
    print("=" * 100)

    classes = {}

    for i, (p, q, n, C) in enumerate(anchors, 1):

        candidates = build_c_class(C)

        classes[C] = candidates

        if i % 25 == 0 or i == len(anchors):

            print(
                f"C {i:3d}/{len(anchors):3d} "
                f"candidates={len(candidates):2d}"
            )

    # ----------------------------------------------------------------------------------------------
    # Actual
    # ----------------------------------------------------------------------------------------------

    actual = actual_dataset(anchors, classes)

    # ----------------------------------------------------------------------------------------------
    # Null
    # ----------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("BUILDING EXACT RE-ANCHOR NULL")
    print("=" * 100)

    null = null_dataset(classes)

    # ----------------------------------------------------------------------------------------------
    # Sanity
    # ----------------------------------------------------------------------------------------------

    identity_failures = 0
    membership_failures = 0

    for row in actual:

        p = row["p"]
        q = row["q"]
        C = row["C"]

        if (p * q) % M != C:
            identity_failures += 1

        if (p, q) not in classes[C]:
            membership_failures += 1

    print()
    print("=" * 100)
    print("SANITY CHECK")
    print("=" * 100)

    print(f"C identity failures       = {identity_failures}")
    print(f"anchor membership failures = {membership_failures}")

    # ----------------------------------------------------------------------------------------------
    # Main comparison
    # ----------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("ACTUAL vs EXACT C-CLASS NULL")
    print("=" * 100)

    print(
        f"{'STATISTIC':24s}"
        f"{'ACTUAL':>16s}"
        f"{'NULL MEAN':>16s}"
        f"{'NULL SD':>16s}"
        f"{'Z':>12s}"
    )

    print("-" * 100)

    for key in STAT_KEYS:

        actual_value = actual_class_mean(actual, key)

        null_mean, null_var = exact_class_null_moments(
            null,
            key
        )

        null_sd = math.sqrt(null_var)

        if null_sd > 0:
            z = (actual_value - null_mean) / null_sd
        else:
            z = 0.0

        print(
            f"{key:24s}"
            f"{actual_value:16.6f}"
            f"{null_mean:16.6f}"
            f"{null_sd:16.6f}"
            f"{z:12.6f}"
        )

    # ----------------------------------------------------------------------------------------------
    # Within-class geometric ranks
    # ----------------------------------------------------------------------------------------------

    ranks = anchor_ranks(anchors, classes)

    print()
    print("=" * 100)
    print("WITHIN-C-CLASS GEOMETRIC RANKS")
    print("=" * 100)

    for key, label in [
        ("sum_rank", "p+q rank"),
        ("diff_rank", "|p-q| rank"),
        ("product_rank", "pq rank"),
        ("ratio_rank", "p/q rank"),
    ]:

        ranks_values = []
        percentile_values = []

        for r in ranks:

            k = {
                "sum_rank": "sum_values",
                "diff_rank": "diff_values",
                "product_rank": "product_values",
                "ratio_rank": "ratio_values",
            }[key]

            rank = r[key]

            size = r["class_size"]

            if size > 1:
                percentile = (rank - 1) / (size - 1)
            else:
                percentile = 0.5

            ranks_values.append(rank)
            percentile_values.append(percentile)

        print()
        print(label)

        print(
            f"  mean rank       = "
            f"{sum(ranks_values)/len(ranks_values):.6f}"
        )

        print(
            f"  mean percentile = "
            f"{sum(percentile_values)/len(percentile_values):.6f}"
        )

        print(
            f"  minimum rank    = "
            f"{min(ranks_values)}"
        )

        print(
            f"  maximum rank    = "
            f"{max(ranks_values)}"
        )

    # ----------------------------------------------------------------------------------------------
    # Most central / most extreme actual anchors
    # ----------------------------------------------------------------------------------------------

    centrality_rows = []

    for row in actual:

        candidates = classes[row["C"]]

        g = anchor_geometry(
            (row["p"], row["q"]),
            candidates
        )

        centrality_rows.append(
            (
                g["sum_dist"],
                row,
                g,
                candidates
            )
        )

    centrality_rows.sort(
        key=lambda x: x[0]
    )

    print()
    print("=" * 100)
    print("MOST GEOMETRICALLY CENTRAL ACTUAL ANCHORS")
    print("=" * 100)

    for distance_sum, row, g, candidates in centrality_rows[:20]:

        print()
        print(
            f"C={row['C']:,} "
            f"anchor=({row['p']:,},{row['q']:,}) "
            f"class={len(candidates)}"
        )

        print(
            f"  sum_dist   = {g['sum_dist']:,.3f}"
        )

        print(
            f"  min_dist   = {g['min_dist']:,.3f}"
        )

        print(
            f"  max_dist   = {g['max_dist']:,.3f}"
        )

        print(
            f"  abs_dsum   = {g['sum_abs_dsum']:,.0f}"
        )

    # ----------------------------------------------------------------------------------------------
    # Most extreme actual anchors
    # ----------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("MOST GEOMETRICALLY EXTREME ACTUAL ANCHORS")
    print("=" * 100)

    for distance_sum, row, g, candidates in centrality_rows[-20:][::-1]:

        print()
        print(
            f"C={row['C']:,} "
            f"anchor=({row['p']:,},{row['q']:,}) "
            f"class={len(candidates)}"
        )

        print(
            f"  sum_dist   = {g['sum_dist']:,.3f}"
        )

        print(
            f"  min_dist   = {g['min_dist']:,.3f}"
        )

        print(
            f"  max_dist   = {g['max_dist']:,.3f}"
        )

        print(
            f"  abs_dsum   = {g['sum_abs_dsum']:,.0f}"
        )

    # ----------------------------------------------------------------------------------------------
    # Classes where centrality actually differs
    # ----------------------------------------------------------------------------------------------

    varying = []

    for C, vertex_stats in null.items():

        values = [
            s["sum_dist"]
            for s in vertex_stats
        ]

        if len(set(values)) > 1:
            varying.append(
                (
                    max(values) - min(values),
                    C,
                    vertex_stats
                )
            )

    varying.sort(reverse=True)

    print()
    print("=" * 100)
    print("C-CLASSES WITH VARYING GEOMETRIC CENTRALITY")
    print("=" * 100)

    print(
        f"classes with non-constant sum_dist = "
        f"{len(varying)}"
    )

    for spread, C, vertex_stats in varying[:30]:

        print()
        print(
            f"C={C:,} spread={spread:,.3f}"
        )

        for s in sorted(
            vertex_stats,
            key=lambda x: x["sum_dist"]
        ):

            print(
                f"  "
                f"({s['anchor'][0]:,},{s['anchor'][1]:,}) "
                f"sum_dist={s['sum_dist']:,.3f} "
                f"min={s['min_dist']:,.3f} "
                f"max={s['max_dist']:,.3f}"
            )

    # ----------------------------------------------------------------------------------------------
    # Final
    # ----------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("FINAL")
    print("=" * 100)

    print(
        "This experiment keeps every observed C = p*q mod M class "
        "fixed and changes only which member is treated as the anchor."
    )

    print()
    print(
        "The primary question is no longer collision degree."
    )

    print(
        "It is whether the ORIGINAL factor pair occupies an "
        "unusual geometric position inside its exact modular class."
    )

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    main()
