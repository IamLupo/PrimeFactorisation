#!/usr/bin/env python3

import math
import random
from collections import Counter, defaultdict

# ==================================================================================================
# EXACT C-CLASS / RE-ANCHOR GRAPH EXPERIMENT
# ==================================================================================================

M = 111_546_435
TMAX = 25

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

TRIALS = 300
SEED = 1_511_464_998

# ================================================================================================
# CONFIGURATION
# ================================================================================================

# Keep this True for the normal experiment.
# Set False if you want to use a supplied list of anchors.
GENERATE_ANCHORS = True

# ================================================================================================
# PRIME GENERATION
# ================================================================================================

def sieve(n: int):
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[0:2] = b"\x00\x00"

    limit = int(math.isqrt(n))

    for p in range(2, limit + 1):
        if sieve[p]:
            start = p * p
            sieve[start:n + 1:p] = b"\x00" * (((n - start) // p) + 1)

    return [i for i in range(FACTOR_MIN, FACTOR_MAX + 1) if sieve[i]]


PRIMES = sieve(FACTOR_MAX)
PRIME_SET = set(PRIMES)

# ================================================================================================
# MODULAR INVERSE
# ================================================================================================

def inv_mod(a: int, m: int):
    g = math.gcd(a, m)
    if g != 1:
        return None
    return pow(a, -1, m)


# ================================================================================================
# ANCHOR GENERATION
# ================================================================================================

def generate_unique_anchors(count: int, seed: int):
    """
    Generate unordered prime pairs p < q with unique C = p*q mod M.

    IMPORTANT:
    For exact comparison with your previous experiments, replace this function
    with the SAME anchor-generation routine used in those experiments if it
    differs. The random seed alone does not guarantee identical anchors if
    the sampling procedure changed.
    """
    rng = random.Random(seed)

    anchors = []
    seen_c = set()

    n_primes = len(PRIMES)

    while len(anchors) < count:
        i = rng.randrange(n_primes)
        j = rng.randrange(n_primes)

        if i == j:
            continue

        p = PRIMES[i]
        q = PRIMES[j]

        if p > q:
            p, q = q, p

        n = p * q
        C = n % M

        if C in seen_c:
            continue

        seen_c.add(C)
        anchors.append((p, q, n, C))

    return anchors


# ================================================================================================
# EXACT C-CLASS ENUMERATION
# ================================================================================================

def build_c_class(C: int):
    """
    Enumerate all unordered prime pairs (a,b), FACTOR_MIN <= a <= b <= FACTOR_MAX,
    satisfying

        a*b == C (mod M)

    Since M > FACTOR_MAX, for invertible a there is at most one possible b
    in the factor interval.
    """
    candidates = set()

    for a in PRIMES:
        inv = inv_mod(a, M)

        if inv is None:
            continue

        b = (C * inv) % M

        if b < FACTOR_MIN or b > FACTOR_MAX:
            continue

        if b not in PRIME_SET:
            continue

        if a > b:
            continue

        if (a * b) % M != C:
            continue

        candidates.add((a, b))

    return sorted(candidates)


# ================================================================================================
# VERTEX STATISTICS
# ================================================================================================

def candidate_statistics(anchor, candidates):
    """
    Treat 'anchor' as the chosen vertex.

    For every other candidate (a,b) in the same C-class:

        a*b - p*q = 2*M*t

    with |t| <= TMAX.

    Returns exact graph statistics for this vertex.
    """

    p, q = anchor
    anchor_product = p * q

    degree = 0
    same_sign = 0
    opposite_sign = 0

    sum_dp = 0
    sum_dq = 0
    sum_dsum = 0
    sum_ddiff = 0

    abs_dp = 0
    abs_dq = 0
    abs_dsum = 0
    abs_ddiff = 0

    events = []

    for a, b in candidates:
        if (a, b) == (p, q):
            continue

        delta_product = a * b - anchor_product

        if delta_product % (2 * M) != 0:
            raise RuntimeError(
                f"Non-integral t: anchor={anchor}, candidate={(a,b)}"
            )

        t = delta_product // (2 * M)

        if abs(t) > TMAX:
            continue

        dp = a - p
        dq = b - q

        dsum = dp + dq
        ddiff = dp - dq

        degree += 1

        if dp == 0 or dq == 0 or dp * dq > 0:
            same_sign += 1
        else:
            opposite_sign += 1

        sum_dp += dp
        sum_dq += dq
        sum_dsum += dsum
        sum_ddiff += ddiff

        abs_dp += abs(dp)
        abs_dq += abs(dq)
        abs_dsum += abs(dsum)
        abs_ddiff += abs(ddiff)

        events.append({
            "t": t,
            "a": a,
            "b": b,
            "dp": dp,
            "dq": dq,
            "dsum": dsum,
            "ddiff": ddiff,
        })

    return {
        "degree": degree,
        "same_sign": same_sign,
        "opposite": opposite_sign,

        "sum_dp": sum_dp,
        "sum_dq": sum_dq,
        "sum_dsum": sum_dsum,
        "sum_ddiff": sum_ddiff,

        "abs_dp": abs_dp,
        "abs_dq": abs_dq,
        "abs_dsum": abs_dsum,
        "abs_ddiff": abs_ddiff,

        "events": events,
    }


# ================================================================================================
# EXACT DISCRETE NULL
# ================================================================================================

def exact_null_distribution(classes):
    """
    Every observed C-class is retained.

    For each C-class, choose one candidate uniformly as the new anchor.

    We compute the EXACT distribution of:
        total event degree
        total same-sign events
        total opposite-sign events

    via dynamic programming.

    No Monte Carlo is used.
    """

    # DP state:
    #   (total_degree, total_same, total_opposite) -> number of equally likely selections
    #
    # We actually keep integer multiplicities. Every class contributes len(class)
    # equally likely possibilities.

    dp = {(0, 0, 0): 1}

    for C, candidates in classes.items():

        choices = []

        for anchor in candidates:
            s = candidate_statistics(anchor, candidates)

            choices.append((
                s["degree"],
                s["same_sign"],
                s["opposite"],
            ))

        new_dp = defaultdict(int)

        for (d0, ss0, oo0), multiplicity in dp.items():
            for degree, same, opposite in choices:
                new_dp[(d0 + degree,
                        ss0 + same,
                        oo0 + opposite)] += multiplicity

        dp = dict(new_dp)

    total_weight = 1
    for candidates in classes.values():
        total_weight *= len(candidates)

    return dp, total_weight


# ================================================================================================
# EXACT MOMENTS
# ================================================================================================

def exact_moments(classes):
    """
    Compute exact expectation and variance for additive statistics.
    """

    means = {
        "degree": 0.0,
        "same_sign": 0.0,
        "opposite": 0.0,

        "sum_dp": 0.0,
        "sum_dq": 0.0,
        "sum_dsum": 0.0,
        "sum_ddiff": 0.0,

        "abs_dp": 0.0,
        "abs_dq": 0.0,
        "abs_dsum": 0.0,
        "abs_ddiff": 0.0,
    }

    variances = {k: 0.0 for k in means}

    for candidates in classes.values():

        local = []

        for anchor in candidates:
            s = candidate_statistics(anchor, candidates)
            local.append(s)

        n = len(local)

        for key in means:
            vals = [float(s[key]) for s in local]

            mu = sum(vals) / n

            if n > 1:
                var = sum((x - mu) ** 2 for x in vals) / n
            else:
                var = 0.0

            means[key] += mu
            variances[key] += var

    return means, variances


# ================================================================================================
# EXACT TWO-SIDED P
# ================================================================================================

def exact_two_sided_p(dp, total_weight, observed_value, statistic_index):
    """
    Discrete exact two-sided probability defined by distance from the
    null mean.
    """

    mean = 0.0
    for state, weight in dp.items():
        mean += state[statistic_index] * weight

    mean /= total_weight

    observed_distance = abs(observed_value - mean)

    tail_weight = 0

    for state, weight in dp.items():
        distance = abs(state[statistic_index] - mean)

        if distance >= observed_distance - 1e-12:
            tail_weight += weight

    return tail_weight / total_weight


# ================================================================================================
# MAIN
# ================================================================================================

def main():

    print("=" * 100)
    print("EXACT C-CLASS / RE-ANCHOR GRAPH EXPERIMENT")
    print("=" * 100)

    print(f"M                    = {M:,}")
    print(f"2M                   = {2*M:,}")
    print(f"factor range         = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"T range              = [-{TMAX}, +{TMAX}]")
    print(f"actual trials        = {TRIALS}")
    print(f"factor primes        = {len(PRIMES):,}")
    print(f"random seed          = {SEED:,}")
    print()

    # --------------------------------------------------------------------------------------------
    # Generate anchors
    # --------------------------------------------------------------------------------------------

    if GENERATE_ANCHORS:
        anchors = generate_unique_anchors(TRIALS, SEED)
    else:
        raise RuntimeError("Supply the previous anchor list here.")

    print("=" * 100)
    print("GENERATING ACTUAL ANCHORS")
    print("=" * 100)

    print(f"anchors = {len(anchors)}")
    print()

    # --------------------------------------------------------------------------------------------
    # Build exact C classes
    # --------------------------------------------------------------------------------------------

    print("=" * 100)
    print("BUILDING EXACT C-CLASSES")
    print("=" * 100)

    classes = {}

    for i, (p, q, n, C) in enumerate(anchors, 1):

        candidates = build_c_class(C)

        if (p, q) not in candidates:
            raise RuntimeError(
                f"ANCHOR NOT IN C-CLASS: {(p,q)} C={C}"
            )

        classes[C] = candidates

        if i % 25 == 0 or i == len(anchors):
            print(
                f"C {i:3d}/{len(anchors):3d} "
                f"candidates={len(candidates):2d}"
            )

    # --------------------------------------------------------------------------------------------
    # Class statistics
    # --------------------------------------------------------------------------------------------

    size_counter = Counter(len(v) for v in classes.values())

    print()
    print("=" * 100)
    print("C-CLASS STRUCTURE")
    print("=" * 100)

    print(f"unique C classes      = {len(classes)}")
    print(
        f"minimum class size   = "
        f"{min(map(len, classes.values()))}"
    )
    print(
        f"maximum class size   = "
        f"{max(map(len, classes.values()))}"
    )
    print(
        f"mean class size      = "
        f"{sum(map(len, classes.values())) / len(classes):.6f}"
    )

    print()
    print("class-size distribution:")

    for size in sorted(size_counter):
        print(
            f"    size={size:2d} classes={size_counter[size]:4d}"
        )

    # --------------------------------------------------------------------------------------------
    # Actual vertex statistics
    # --------------------------------------------------------------------------------------------

    actual_stats = []

    for p, q, n, C in anchors:

        candidates = classes[C]

        s = candidate_statistics((p, q), candidates)

        actual_stats.append({
            "p": p,
            "q": q,
            "n": n,
            "C": C,
            **s,
        })

    # --------------------------------------------------------------------------------------------
    # Verify graph arithmetic
    # --------------------------------------------------------------------------------------------

    identity_failures = 0
    membership_failures = 0
    t_failures = 0

    for row in actual_stats:

        p = row["p"]
        q = row["q"]
        C = row["C"]

        if (p * q) % M != C:
            identity_failures += 1

        if (p, q) not in classes[C]:
            membership_failures += 1

        for e in row["events"]:
            if (e["a"] * e["b"] - p * q) != 2 * M * e["t"]:
                t_failures += 1

    print()
    print("=" * 100)
    print("EXACT GRAPH SANITY CHECK")
    print("=" * 100)

    print(f"C identity failures           = {identity_failures}")
    print(f"anchor membership failures    = {membership_failures}")
    print(f"t reconstruction failures     = {t_failures}")

    # --------------------------------------------------------------------------------------------
    # Actual totals
    # --------------------------------------------------------------------------------------------

    def actual_total(key):
        return sum(s[key] for s in actual_stats)

    actual_values = {
        key: actual_total(key)
        for key in [
            "degree",
            "same_sign",
            "opposite",
            "sum_dp",
            "sum_dq",
            "sum_dsum",
            "sum_ddiff",
            "abs_dp",
            "abs_dq",
            "abs_dsum",
            "abs_ddiff",
        ]
    }

    events = actual_values["degree"]

    # --------------------------------------------------------------------------------------------
    # Exact null moments
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("COMPUTING EXACT RE-ANCHOR NULL")
    print("=" * 100)

    null_means, null_vars = exact_moments(classes)

    # --------------------------------------------------------------------------------------------
    # Z table
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("ACTUAL vs EXACT C-CLASS NULL")
    print("=" * 100)

    print(
        f"{'STATISTIC':28s}"
        f"{'ACTUAL':>16s}"
        f"{'NULL MEAN':>16s}"
        f"{'NULL SD':>16s}"
        f"{'Z':>12s}"
    )
    print("-" * 100)

    for key in [
        "degree",
        "same_sign",
        "opposite",
        "sum_dp",
        "sum_dq",
        "sum_dsum",
        "sum_ddiff",
        "abs_dp",
        "abs_dq",
        "abs_dsum",
        "abs_ddiff",
    ]:

        actual = float(actual_values[key])
        mean = null_means[key]
        sd = math.sqrt(null_vars[key])

        if sd > 0:
            z = (actual - mean) / sd
        else:
            z = 0.0

        print(
            f"{key:28s}"
            f"{actual:16.6f}"
            f"{mean:16.6f}"
            f"{sd:16.6f}"
            f"{z:12.6f}"
        )

    # --------------------------------------------------------------------------------------------
    # Exact discrete distribution
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("BUILDING EXACT DISCRETE DISTRIBUTION")
    print("=" * 100)

    dp, total_weight = exact_null_distribution(classes)

    print(f"null states            = {len(dp):,}")
    print(f"total equally-weighted selections = {total_weight:,}")

    # --------------------------------------------------------------------------------------------
    # Exact p-values
    # --------------------------------------------------------------------------------------------

    p_degree = exact_two_sided_p(
        dp,
        total_weight,
        actual_values["degree"],
        0,
    )

    p_same = exact_two_sided_p(
        dp,
        total_weight,
        actual_values["same_sign"],
        1,
    )

    p_opp = exact_two_sided_p(
        dp,
        total_weight,
        actual_values["opposite"],
        2,
    )

    print()
    print("=" * 100)
    print("EXACT TWO-SIDED P-VALUES")
    print("=" * 100)

    print(f"event degree       = {p_degree:.10f}")
    print(f"same-sign events   = {p_same:.10f}")
    print(f"opposite events    = {p_opp:.10f}")

    # --------------------------------------------------------------------------------------------
    # Actual anchor rank inside every class
    # --------------------------------------------------------------------------------------------

    max_degree_count = 0
    strictly_max_count = 0
    min_degree_count = 0

    degree_ranks = []

    for row in actual_stats:

        candidates = classes[row["C"]]

        degrees = []

        for candidate in candidates:
            s = candidate_statistics(candidate, candidates)
            degrees.append(s["degree"])

        d = row["degree"]

        greater = sum(x > d for x in degrees)
        equal = sum(x == d for x in degrees)
        lower = sum(x < d for x in degrees)

        degree_ranks.append((greater, equal, lower, len(degrees)))

        if d == max(degrees):
            max_degree_count += 1

        if d == max(degrees) and degrees.count(d) == 1:
            strictly_max_count += 1

        if d == min(degrees):
            min_degree_count += 1

    print()
    print("=" * 100)
    print("ANCHOR POSITION INSIDE C-CLASS")
    print("=" * 100)

    print(
        f"anchors at maximum degree       = "
        f"{max_degree_count:4d}/{len(actual_stats)}"
    )

    print(
        f"anchors at unique maximum       = "
        f"{strictly_max_count:4d}/{len(actual_stats)}"
    )

    print(
        f"anchors at minimum degree       = "
        f"{min_degree_count:4d}/{len(actual_stats)}"
    )

    # --------------------------------------------------------------------------------------------
    # Degree breakdown
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("ACTUAL ANCHOR DEGREE")
    print("=" * 100)

    degree_counter = Counter(row["degree"] for row in actual_stats)

    for d in sorted(degree_counter):
        print(
            f"degree={d:2d} anchors={degree_counter[d]:4d}"
        )

    # --------------------------------------------------------------------------------------------
    # Show classes where re-anchoring actually changes degree
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("CLASSES WITH NON-CONSTANT ANCHOR DEGREE")
    print("=" * 100)

    changed = []

    for row in actual_stats:

        candidates = classes[row["C"]]

        records = []

        for candidate in candidates:
            s = candidate_statistics(candidate, candidates)

            records.append(
                (
                    candidate,
                    s["degree"],
                    s["same_sign"],
                    s["sum_dsum"],
                )
            )

        degrees = {r[1] for r in records}

        if len(degrees) > 1:
            changed.append((row, records))

    print(f"classes with varying degree = {len(changed)}")

    for row, records in changed[:50]:

        print()
        print(
            f"C={row['C']:,} "
            f"original=({row['p']:,},{row['q']:,})"
        )

        for candidate, degree, same, dsum in records:

            marker = "*"

            if candidate != (row["p"], row["q"]):
                marker = " "

            print(
                f"  {marker} ({candidate[0]:,},{candidate[1]:,}) "
                f"degree={degree:2d} "
                f"same={same:2d} "
                f"sum_dSum={dsum:+,}"
            )

    # --------------------------------------------------------------------------------------------
    # Strongest examples
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("MOST EXTREME RE-ANCHOR EFFECTS")
    print("=" * 100)

    extreme = []

    for row in actual_stats:

        candidates = classes[row["C"]]

        vertex_records = []

        for candidate in candidates:

            s = candidate_statistics(candidate, candidates)

            vertex_records.append(
                (
                    candidate,
                    s,
                )
            )

        original = (row["p"], row["q"])

        original_degree = row["degree"]

        max_degree = max(s["degree"] for _, s in vertex_records)
        min_degree = min(s["degree"] for _, s in vertex_records)

        extreme.append(
            (
                max_degree - min_degree,
                row,
                vertex_records,
            )
        )

    extreme.sort(reverse=True, key=lambda x: x[0])

    for spread, row, vertex_records in extreme[:20]:

        print()
        print(
            f"C={row['C']:,} "
            f"original=({row['p']:,},{row['q']:,}) "
            f"degree-spread={spread}"
        )

        for candidate, s in vertex_records:

            marker = "* " if candidate == (row["p"], row["q"]) else "  "

            print(
                f"{marker}"
                f"({candidate[0]:,},{candidate[1]:,}) "
                f"degree={s['degree']:2d} "
                f"same={s['same_sign']:2d} "
                f"opp={s['opposite']:2d}"
            )

    # --------------------------------------------------------------------------------------------
    # Full observed event count
    # --------------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("FINAL")
    print("=" * 100)

    print(f"actual anchors             = {len(actual_stats)}")
    print(f"actual collision events    = {events}")

    print()
    print(
        "The null keeps every observed C = p*q mod M class fixed and "
        "replaces only the chosen vertex."
    )

    print()
    print(
        "A strong anchor-specific result would require the ORIGINAL "
        "anchors to occupy unusual positions inside these finite classes."
    )

    print()
    print(
        "Most important numbers:"
    )
    print(
        f"    exact null mean degree = "
        f"{null_means['degree']:.8f}"
    )
    print(
        f"    exact null SD degree   = "
        f"{math.sqrt(null_vars['degree']):.8f}"
    )
    print(
        f"    actual degree          = "
        f"{events}"
    )
    print(
        f"    exact degree p-value   = "
        f"{p_degree:.10f}"
    )

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    main()
