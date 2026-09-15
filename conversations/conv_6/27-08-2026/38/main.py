#!/usr/bin/env python3

import math
import random
from collections import Counter

# ==================================================================================================
# CRT COORDINATE / EXACT C-CLASS RE-ANCHORING EXPERIMENT
# ==================================================================================================

M = 111_546_435

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

TRIALS = 300
TMAX = 25

SEED = 1_511_464_998

# M is square-free:
#     3 * 5 * 7 * 11 * 13 * 17 * 19 * 23
MODULI = [3, 5, 7, 11, 13, 17, 19, 23]

assert math.prod(MODULI) == M


# ==================================================================================================
# PRIME SIEVE
# ==================================================================================================

def sieve(n):
    s = bytearray(b"\x01") * (n + 1)
    s[0] = 0
    s[1] = 0

    for p in range(2, math.isqrt(n) + 1):
        if s[p]:
            start = p * p
            s[start:n + 1:p] = b"\x00" * (
                ((n - start) // p) + 1
            )

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

        inv_a = pow(a, -1, M)
        b = (C * inv_a) % M

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
# CRT RESIDUE HELPERS
# ==================================================================================================

def cyclic_distance(x, y, r):
    """
    Minimal cyclic distance modulo r.
    """

    d = abs(x - y)

    return min(d, r - d)


def signed_cyclic_delta(x, y, r):
    """
    Signed shortest displacement x-y modulo r.

    Output is chosen in approximately [-r/2, r/2].
    """

    d = (x - y) % r

    if d > r // 2:
        d -= r

    return d


# ==================================================================================================
# CRT SIGNATURE
# ==================================================================================================

def crt_signature(pair):
    a, b = pair

    return tuple(
        (a % r, b % r)
        for r in MODULI
    )


# ==================================================================================================
# CRT ANCHOR STATISTICS
# ==================================================================================================

def crt_anchor_statistics(anchor, candidates):

    p, q = anchor

    others = [
        pair
        for pair in candidates
        if pair != anchor
    ]

    out = {
        "class_size": len(candidates),

        "sum_cyclic_a": 0.0,
        "sum_cyclic_b": 0.0,
        "sum_cyclic_total": 0.0,

        "sum_abs_signed_a": 0.0,
        "sum_abs_signed_b": 0.0,
        "sum_abs_signed_total": 0.0,

        "same_residue_count": 0,
        "opposite_residue_count": 0,

        "mod_agreement_score": 0.0,
    }

    if not others:
        return out

    for r in MODULI:

        p_r = p % r
        q_r = q % r

        for a, b in others:

            a_r = a % r
            b_r = b % r

            da = cyclic_distance(a_r, p_r, r)
            db = cyclic_distance(b_r, q_r, r)

            sda = signed_cyclic_delta(a_r, p_r, r)
            sdb = signed_cyclic_delta(b_r, q_r, r)

            out["sum_cyclic_a"] += da
            out["sum_cyclic_b"] += db
            out["sum_cyclic_total"] += da + db

            out["sum_abs_signed_a"] += abs(sda)
            out["sum_abs_signed_b"] += abs(sdb)
            out["sum_abs_signed_total"] += abs(sda) + abs(sdb)

            # CRT orientation:
            #
            # same sign means both coordinates move in the
            # same modular direction relative to the anchor.
            #
            # opposite sign means they move in different
            # modular directions.

            if sda != 0 and sdb != 0:
                if (sda > 0) == (sdb > 0):
                    out["same_residue_count"] += 1
                else:
                    out["opposite_residue_count"] += 1

            # Count coordinate agreement events.
            if a_r == p_r:
                out["mod_agreement_score"] += 1.0

            if b_r == q_r:
                out["mod_agreement_score"] += 1.0

    return out


# ==================================================================================================
# C-CLASS NULL
# ==================================================================================================

def build_null(classes):

    null = {}

    for C, candidates in classes.items():

        rows = []

        for anchor in candidates:

            stats = crt_anchor_statistics(
                anchor,
                candidates
            )

            rows.append({
                "anchor": anchor,
                **stats,
            })

        null[C] = rows

    return null


# ==================================================================================================
# ACTUAL DATASET
# ==================================================================================================

def build_actual(anchors, classes):

    rows = []

    for p, q, n, C in anchors:

        stats = crt_anchor_statistics(
            (p, q),
            classes[C]
        )

        rows.append({
            "p": p,
            "q": q,
            "n": n,
            "C": C,
            **stats,
        })

    return rows


# ==================================================================================================
# CLASS-WEIGHTED EXACT NULL
# ==================================================================================================

STAT_KEYS = [
    "sum_cyclic_a",
    "sum_cyclic_b",
    "sum_cyclic_total",
    "sum_abs_signed_a",
    "sum_abs_signed_b",
    "sum_abs_signed_total",
    "same_residue_count",
    "opposite_residue_count",
    "mod_agreement_score",
]


def actual_mean(rows, key):
    if not rows:
        return 0.0

    return sum(float(r[key]) for r in rows) / len(rows)


def exact_class_null_mean_sd(null, key):

    """
    Every C-class contributes exactly one re-anchored vertex.

    Therefore:
        1. average within each class
        2. average across classes

    This prevents classes with more vertices from dominating.
    """

    class_means = []

    for rows in null.values():

        if not rows:
            continue

        mu = sum(
            float(r[key])
            for r in rows
        ) / len(rows)

        class_means.append(mu)

    if not class_means:
        return 0.0, 0.0

    mean = sum(class_means) / len(class_means)

    var = (
        sum((x - mean) ** 2 for x in class_means)
        / len(class_means)
    )

    return mean, math.sqrt(max(0.0, var))


# ==================================================================================================
# WITHIN-CLASS RANK
# ==================================================================================================

def within_class_percentiles(anchors, classes, key):

    values = []

    for p, q, n, C in anchors:

        candidates = classes[C]

        anchor = (p, q)

        stats = []

        for candidate in candidates:

            s = crt_anchor_statistics(
                candidate,
                candidates
            )

            stats.append(
                (candidate, float(s[key]))
            )

        ours = next(
            value
            for candidate, value in stats
            if candidate == anchor
        )

        ordered = sorted(
            value
            for _, value in stats
        )

        if len(ordered) <= 1:
            percentile = 0.5
        else:
            less = sum(v < ours for v in ordered)
            equal = sum(v == ours for v in ordered)

            percentile = (
                less + 0.5 * equal
            ) / len(ordered)

        values.append(percentile)

    return values


# ==================================================================================================
# CRT SIGNATURE COUNTS
# ==================================================================================================

def signature_counts(rows):

    counter = Counter()

    for r in rows:

        pair = (r["p"], r["q"])

        counter[crt_signature(pair)] += 1

    return counter


# ==================================================================================================
# MODULUS-BY-MODULUS ANALYSIS
# ==================================================================================================

def per_modulus_statistics(anchors, classes):

    print()
    print("=" * 100)
    print("PER-MODULUS CRT ANCHOR STATISTICS")
    print("=" * 100)

    print(
        f"{'r':>4} "
        f"{'A-SAME':>10} "
        f"{'A-CROSS':>10} "
        f"{'A-DIST':>12} "
        f"{'A-AGREE':>12} "
        f"{'NULL-DIST':>12}"
    )

    for r in MODULI:

        actual_same = 0
        actual_cross = 0
        actual_dist = 0.0
        actual_agree = 0.0

        null_dist_values = []

        for p, q, n, C in anchors:

            candidates = classes[C]

            others = [
                x for x in candidates
                if x != (p, q)
            ]

            if not others:
                continue

            p_r = p % r
            q_r = q % r

            for a, b in others:

                a_r = a % r
                b_r = b % r

                da = signed_cyclic_delta(
                    a_r, p_r, r
                )

                db = signed_cyclic_delta(
                    b_r, q_r, r
                )

                actual_dist += abs(da) + abs(db)

                if da != 0 and db != 0:

                    if (da > 0) == (db > 0):
                        actual_same += 1
                    else:
                        actual_cross += 1

                if a_r == p_r:
                    actual_agree += 1

                if b_r == q_r:
                    actual_agree += 1

            # exact re-anchor null for this modulus
            for candidate in candidates:

                ca, cb = candidate

                c_r_a = ca % r
                c_r_b = cb % r

                dtotal = 0

                for x, y in candidates:

                    if (x, y) == candidate:
                        continue

                    dtotal += cyclic_distance(
                        x % r,
                        c_r_a,
                        r
                    )

                    dtotal += cyclic_distance(
                        y % r,
                        c_r_b,
                        r
                    )

                null_dist_values.append(dtotal)

        print(
            f"{r:4d} "
            f"{actual_same:10d} "
            f"{actual_cross:10d} "
            f"{actual_dist:12.3f} "
            f"{actual_agree:12.3f} "
            f"{(
                sum(null_dist_values) / len(null_dist_values)
                if null_dist_values else 0.0
            ):12.3f}"
        )


# ==================================================================================================
# MAIN
# ==================================================================================================

def main():

    print("=" * 100)
    print("CRT COORDINATE / EXACT C-CLASS RE-ANCHORING EXPERIMENT")
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

    anchors = generate_unique_anchors(
        TRIALS,
        SEED
    )

    print(f"anchors = {len(anchors)}")

    # ----------------------------------------------------------------------------------------------
    # Classes
    # ----------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("BUILDING EXACT C-CLASSES")
    print("=" * 100)

    classes = {}

    for i, (p, q, n, C) in enumerate(
        anchors,
        1
    ):

        candidates = build_c_class(C)

        classes[C] = candidates

        if i % 25 == 0 or i == len(anchors):

            print(
                f"C {i:3d}/{len(anchors):3d} "
                f"candidates={len(candidates):2d}"
            )

    # ----------------------------------------------------------------------------------------------
    # Sanity
    # ----------------------------------------------------------------------------------------------

    identity_failures = 0
    membership_failures = 0

    for p, q, n, C in anchors:

        if (p * q) % M != C:
            identity_failures += 1

        if (p, q) not in classes[C]:
            membership_failures += 1

    print()
    print("=" * 100)
    print("SANITY CHECK")
    print("=" * 100)

    print(
        f"C identity failures        = "
        f"{identity_failures}"
    )

    print(
        f"anchor membership failures = "
        f"{membership_failures}"
    )

    # ----------------------------------------------------------------------------------------------
    # Actual
    # ----------------------------------------------------------------------------------------------

    actual = build_actual(
        anchors,
        classes
    )

    # ----------------------------------------------------------------------------------------------
    # Null
    # ----------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("BUILDING EXACT CRT RE-ANCHOR NULL")
    print("=" * 100)

    null = build_null(classes)

    # ----------------------------------------------------------------------------------------------
    # Main comparison
    # ----------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("ACTUAL vs EXACT C-CLASS CRT NULL")
    print("=" * 100)

    print(
        f"{'STATISTIC':28s}"
        f"{'ACTUAL':>16s}"
        f"{'NULL MEAN':>16s}"
        f"{'NULL SD':>16s}"
        f"{'Z':>12s}"
    )

    print("-" * 100)

    for key in STAT_KEYS:

        av = actual_mean(
            actual,
            key
        )

        nm, ns = exact_class_null_mean_sd(
            null,
            key
        )

        if ns > 0:
            z = (av - nm) / ns
        else:
            z = 0.0

        print(
            f"{key:28s}"
            f"{av:16.6f}"
            f"{nm:16.6f}"
            f"{ns:16.6f}"
            f"{z:12.6f}"
        )

    # ----------------------------------------------------------------------------------------------
    # Per modulus
    # ----------------------------------------------------------------------------------------------

    per_modulus_statistics(
        anchors,
        classes
    )

    # ----------------------------------------------------------------------------------------------
    # Within-class percentile tests
    # ----------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("WITHIN-C-CLASS CRT PERCENTILES")
    print("=" * 100)

    for key, label in [
        ("sum_cyclic_a", "CRT-A cyclic distance"),
        ("sum_cyclic_b", "CRT-B cyclic distance"),
        ("sum_cyclic_total", "CRT total cyclic distance"),
        ("sum_abs_signed_total", "CRT signed displacement"),
        ("mod_agreement_score", "CRT residue agreement"),
    ]:

        percentiles = within_class_percentiles(
            anchors,
            classes,
            key
        )

        print()
        print(label)

        print(
            f"  mean percentile = "
            f"{sum(percentiles) / len(percentiles):.6f}"
        )

        print(
            f"  minimum        = "
            f"{min(percentiles):.6f}"
        )

        print(
            f"  maximum        = "
            f"{max(percentiles):.6f}"
        )

        print(
            f"  <= 0.10        = "
            f"{sum(x <= 0.10 for x in percentiles)}"
        )

        print(
            f"  >= 0.90        = "
            f"{sum(x >= 0.90 for x in percentiles)}"
        )

    # ----------------------------------------------------------------------------------------------
    # Cross-modulus signatures
    # ----------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("CRT SIGNATURE STRUCTURE")
    print("=" * 100)

    actual_signature_counter = Counter()

    for p, q, n, C in anchors:
        actual_signature_counter[
            crt_signature((p, q))
        ] += 1

    repeated = [
        (sig, count)
        for sig, count in actual_signature_counter.items()
        if count > 1
    ]

    print(
        f"unique full CRT signatures = "
        f"{len(actual_signature_counter)}"
    )

    print(
        f"repeated signatures        = "
        f"{len(repeated)}"
    )

    if repeated:

        print()
        print("Repeated signatures:")

        for sig, count in sorted(
            repeated,
            key=lambda x: (-x[1], x[0])
        )[:20]:

            print(
                f"  count={count} "
                f"{sig}"
            )

    # ----------------------------------------------------------------------------------------------
    # Largest deviations
    # ----------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("MOST EXTREME CRT ANCHOR POSITIONS")
    print("=" * 100)

    ranked = []

    for row in actual:

        candidate_rows = null[row["C"]]

        lookup = {
            x["anchor"]: x
            for x in candidate_rows
        }

        ours = lookup[(row["p"], row["q"])]

        values = [
            x["sum_cyclic_total"]
            for x in candidate_rows
        ]

        actual_value = ours["sum_cyclic_total"]

        percentile = (
            exact_percentile_position(
                actual_value,
                values
            )
        )

        ranked.append(
            (
                abs(percentile - 0.5),
                percentile,
                row
            )
        )

    ranked.sort(
        reverse=True,
        key=lambda x: x[0]
    )

    for _, percentile, row in ranked[:25]:

        print(
            f"C={row['C']:,} "
            f"anchor=({row['p']:,},{row['q']:,}) "
            f"class={row['class_size']} "
            f"percentile={percentile:.6f} "
            f"CRTdist={row['sum_cyclic_total']:.3f}"
        )

    # ----------------------------------------------------------------------------------------------
    # Final
    # ----------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("FINAL")
    print("=" * 100)

    print(
        "The C-class is held fixed."
    )

    print(
        "The only randomized quantity is which prime-pair vertex "
        "is designated as the anchor."
    )

    print()
    print(
        "This tests whether the observed factor pair has a special "
        "position in the CRT coordinate representation of its "
        "exact modular hyperbola."
    )

    print()
    print(
        "The strongest result would be a statistic that is extreme "
        "both globally and in the within-C-class percentile analysis."
    )

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    main()
