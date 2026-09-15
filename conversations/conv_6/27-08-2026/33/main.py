#!/usr/bin/env python3

import math
import random
from collections import Counter, defaultdict
from statistics import mean


# =============================================================================
# CRT-HYPERBOLA / C-CONDITIONED RE-ANCHORING NULL
#
# Question:
#
#   Given C = p*q (mod M), is the original anchor (p,q) unusually positioned
#   among all prime factor pairs (a,b) satisfying
#
#       a*b = C (mod M)
#
#   when we measure collisions
#
#       a*b = p*q + 2*M*t
#
#   for |t| <= T_MAX?
#
# Null:
#   Keep every anchor's C fixed, but replace its anchor pair by a random
#   candidate pair from the SAME modular hyperbola C.
#
# This is deliberately much smaller than the previous 200k-draw experiment.
# The candidate sets are built once and the null permutations are cheap.
# =============================================================================


# =============================================================================
# CONFIGURATION
# =============================================================================

M = 111_546_435
T_MIN = -25
T_MAX = 25

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

ACTUAL_TRIALS = 300
PERMUTATIONS = 5000

SEED = 1_511_464_998

MODULI = [3, 5, 7, 11, 13, 17, 19, 23]

random.seed(SEED)


# =============================================================================
# PRIME SIEVE
# =============================================================================

def sieve(n):
    data = bytearray(b"\x01") * (n + 1)

    data[0] = 0
    data[1] = 0

    for p in range(2, math.isqrt(n) + 1):
        if data[p]:
            start = p * p
            count = ((n - start) // p) + 1
            data[start:n + 1:p] = b"\x00" * count

    return [
        x
        for x in range(FACTOR_MIN, n + 1)
        if data[x]
    ]


PRIMES = sieve(FACTOR_MAX)
PRIME_SET = set(PRIMES)


# =============================================================================
# INVERSES MOD M
# =============================================================================

INVERSE = {}

for p in PRIMES:
    if math.gcd(p, M) == 1:
        INVERSE[p] = pow(p, -1, M)


# =============================================================================
# ACTUAL ANCHORS
#
# Same deterministic construction used by the previous experiment.
# =============================================================================

def generate_actual_anchors(count):
    anchors = []
    seen = set()

    while len(anchors) < count:

        p = random.choice(PRIMES)
        q = random.choice(PRIMES)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        n = p * q

        if n in seen:
            continue

        if math.gcd(n, M) != 1:
            continue

        seen.add(n)

        anchors.append(
            {
                "p": p,
                "q": q,
                "n": n,
                "C": n % M,
            }
        )

    return anchors


# =============================================================================
# BUILD ALL PRIME POINTS ON A MOD-M HYPERBOLA
#
# For fixed C:
#
#     a*b = C (mod M)
#
# and gcd(a,M)=1, therefore
#
#     b = C * a^{-1} (mod M)
#
# is uniquely determined.
# =============================================================================

def build_candidates(C):

    candidates = set()

    for a in PRIMES:

        inv = INVERSE.get(a)

        if inv is None:
            continue

        b = (C * inv) % M

        if b < FACTOR_MIN:
            continue

        if b > FACTOR_MAX:
            continue

        if b not in PRIME_SET:
            continue

        if a == b:
            continue

        if a > b:
            a2, b2 = b, a
        else:
            a2, b2 = a, b

        candidates.add((a2, b2))

    return sorted(candidates)


# =============================================================================
# COLLISIONS FROM AN ANCHOR
# =============================================================================

def anchor_collisions(p, q, candidates):

    n = p * q

    result = []

    for a, b in candidates:

        product = a * b

        diff = product - n

        if diff == 0:
            continue

        denominator = 2 * M

        if diff % denominator != 0:
            continue

        t = diff // denominator

        if T_MIN <= t <= T_MAX:

            result.append(
                {
                    "p": p,
                    "q": q,
                    "a": a,
                    "b": b,
                    "t": t,
                    "n": n,
                    "Nt": product,
                }
            )

    return result


# =============================================================================
# DATASET STATISTICS
# =============================================================================

def normalized_factor(x):
    return (
        x - FACTOR_MIN
    ) / (
        FACTOR_MAX - FACTOR_MIN
    )


def event_statistics(events):

    if not events:
        return {
            "count": 0,
            "mean_ndp": 0.0,
            "mean_ndq": 0.0,
            "mean_dsum": 0.0,
            "mean_ddiff": 0.0,
            "same_sign": 0,
            "opposite_sign": 0,
        }

    dp = []
    dq = []
    dsum = []
    ddiff = []

    same = 0
    opposite = 0

    for e in events:

        x = e["a"] - e["p"]
        y = e["b"] - e["q"]

        dp.append(normalized_factor(e["a"]) -
                  normalized_factor(e["p"]))

        dq.append(normalized_factor(e["b"]) -
                  normalized_factor(e["q"]))

        dsum.append(x + y)
        ddiff.append(x - y)

        if x != 0 and y != 0:

            if (x > 0 and y > 0) or (x < 0 and y < 0):
                same += 1
            else:
                opposite += 1

    return {
        "count": len(events),
        "mean_ndp": mean(dp),
        "mean_ndq": mean(dq),
        "mean_dsum": mean(dsum),
        "mean_ddiff": mean(ddiff),
        "same_sign": same,
        "opposite_sign": opposite,
    }


def t_distribution(events):
    return Counter(e["t"] for e in events)


# =============================================================================
# MODULAR MUTUAL INFORMATION
# =============================================================================

def mutual_information(x, y):

    if not x:
        return 0.0

    n = len(x)

    joint = Counter(zip(x, y))
    cx = Counter(x)
    cy = Counter(y)

    result = 0.0

    for (vx, vy), count in joint.items():

        pxy = count / n
        px = cx[vx] / n
        py = cy[vy] / n

        result += pxy * math.log2(
            pxy / (px * py)
        )

    return result


def modular_mi(events, r):

    a = [e["a"] % r for e in events]
    b = [e["b"] % r for e in events]

    return mutual_information(a, b)


# =============================================================================
# EMPIRICAL P-VALUE
# =============================================================================

def two_sided_p(null_values, actual):

    if not null_values:
        return 1.0

    center = mean(null_values)

    observed = abs(actual - center)

    extreme = sum(
        abs(x - center) >= observed
        for x in null_values
    )

    return (
        (extreme + 1)
        / (len(null_values) + 1)
    )


# =============================================================================
# HEADER
# =============================================================================

print("=" * 100)
print("CRT-HYPERBOLA / C-CONDITIONED RE-ANCHORING NULL")
print("=" * 100)

print(f"M                    = {M:,}")
print(f"2M                   = {2 * M:,}")
print(
    f"factor range        = "
    f"{FACTOR_MIN:,} - {FACTOR_MAX:,}"
)
print(
    f"T range              = "
    f"[{T_MIN}, {T_MAX}]"
)
print(f"actual anchors       = {ACTUAL_TRIALS:,}")
print(f"permutations         = {PERMUTATIONS:,}")
print(f"factor primes        = {len(PRIMES):,}")
print(f"seed                 = {SEED:,}")


# =============================================================================
# GENERATE ACTUAL ANCHORS
# =============================================================================

print()
print("=" * 100)
print("GENERATING ACTUAL ANCHORS")
print("=" * 100)

anchors = generate_actual_anchors(
    ACTUAL_TRIALS
)

print(
    f"anchors = {len(anchors):,}"
)


# =============================================================================
# BUILD UNIQUE C CANDIDATE SETS
# =============================================================================

unique_C = sorted(
    set(a["C"] for a in anchors)
)

candidate_cache = {}

print()
print("=" * 100)
print("BUILDING C-CONDITIONED HYPERBOLA CANDIDATES")
print("=" * 100)

for i, C in enumerate(unique_C, 1):

    candidate_cache[C] = build_candidates(C)

    if (
        i % 25 == 0
        or i == len(unique_C)
    ):
        print(
            f"C {i:4d}/"
            f"{len(unique_C):4d} "
            f"candidates="
            f"{len(candidate_cache[C]):3d}"
        )


candidate_sizes = [
    len(candidate_cache[a["C"]])
    for a in anchors
]

print()
print(
    f"unique C values        = "
    f"{len(unique_C):,}"
)
print(
    f"minimum candidate set  = "
    f"{min(candidate_sizes)}"
)
print(
    f"maximum candidate set  = "
    f"{max(candidate_sizes)}"
)
print(
    f"mean candidate set     = "
    f"{mean(candidate_sizes):.4f}"
)


# =============================================================================
# VERIFY ACTUAL ANCHORS BELONG TO THEIR C-CLASS
# =============================================================================

anchor_membership_failures = 0

for a in anchors:

    pair = (a["p"], a["q"])

    if pair not in candidate_cache[a["C"]]:
        anchor_membership_failures += 1


print()
print("=" * 100)
print("ANCHOR MEMBERSHIP CHECK")
print("=" * 100)

print(
    f"membership failures = "
    f"{anchor_membership_failures}"
)


# =============================================================================
# ACTUAL EVENTS
# =============================================================================

actual_events = []

events_per_anchor_actual = []

for anchor_id, anchor in enumerate(anchors):

    candidates = candidate_cache[
        anchor["C"]
    ]

    events = anchor_collisions(
        anchor["p"],
        anchor["q"],
        candidates
    )

    events_per_anchor_actual.append(
        len(events)
    )

    for e in events:
        e["anchor_id"] = anchor_id

    actual_events.extend(events)


ACTUAL_STATS = event_statistics(
    actual_events
)

ACTUAL_T = t_distribution(
    actual_events
)


# =============================================================================
# DIRECT ALGEBRA CHECK
# =============================================================================

identity_failures = 0
t_failures = 0

for e in actual_events:

    if e["a"] * e["b"] != (
        e["n"] + 2 * M * e["t"]
    ):
        identity_failures += 1

    diff = (
        e["a"] * e["b"]
        - e["p"] * e["q"]
    )

    if diff // (2 * M) != e["t"]:
        t_failures += 1


print()
print("=" * 100)
print("ACTUAL LANDSCAPE")
print("=" * 100)

print(
    f"events                 = "
    f"{len(actual_events):,}"
)

print(
    f"anchors with events    = "
    f"{sum(x > 0 for x in events_per_anchor_actual):,}"
)

print(
    f"identity failures      = "
    f"{identity_failures}"
)

print(
    f"t reconstruction       = "
    f"{t_failures}"
)


# =============================================================================
# NULL PERMUTATIONS
#
# For each original anchor:
#
#   C is fixed.
#
#   A random prime pair from the SAME C candidate set becomes the new anchor.
#
# This preserves:
#   - number of anchors
#   - C distribution
#   - factor range
#   - primality
#
# It destroys:
#   - the original choice of p,q within its C-class.
# =============================================================================

print()
print("=" * 100)
print("C-CONDITIONED RE-ANCHORING NULL")
print("=" * 100)

null_event_counts = []
null_mean_ndp = []
null_mean_ndq = []
null_mean_dsum = []
null_mean_ddiff = []
null_same_fraction = []

null_mi = {
    r: []
    for r in MODULI
}

null_t_distributions = {
    t: []
    for t in range(T_MIN, T_MAX + 1)
    if t != 0
}

for permutation in range(1, PERMUTATIONS + 1):

    permutation_events = []

    for anchor in anchors:

        candidates = candidate_cache[
            anchor["C"]
        ]

        # Choose a different candidate whenever possible.
        if len(candidates) <= 1:
            new_p = anchor["p"]
            new_q = anchor["q"]
        else:
            pair = random.choice(candidates)
            new_p, new_q = pair

        events = anchor_collisions(
            new_p,
            new_q,
            candidates
        )

        permutation_events.extend(events)

    stats = event_statistics(
        permutation_events
    )

    null_event_counts.append(
        stats["count"]
    )

    null_mean_ndp.append(
        stats["mean_ndp"]
    )

    null_mean_ndq.append(
        stats["mean_ndq"]
    )

    null_mean_dsum.append(
        stats["mean_dsum"]
    )

    null_mean_ddiff.append(
        stats["mean_ddiff"]
    )

    total_sign = (
        stats["same_sign"]
        + stats["opposite_sign"]
    )

    if total_sign:
        null_same_fraction.append(
            stats["same_sign"] / total_sign
        )
    else:
        null_same_fraction.append(0.0)

    dist = t_distribution(
        permutation_events
    )

    for t in null_t_distributions:
        null_t_distributions[t].append(
            dist.get(t, 0)
        )

    for r in MODULI:
        null_mi[r].append(
            modular_mi(
                permutation_events,
                r
            )
        )

    if (
        permutation % 500 == 0
        or permutation == PERMUTATIONS
    ):
        print(
            f"permutation "
            f"{permutation:5d}/"
            f"{PERMUTATIONS}"
        )


# =============================================================================
# NULL SUMMARY
# =============================================================================

print()
print("=" * 100)
print("ACTUAL vs C-CONDITIONED NULL")
print("=" * 100)

print(
    f"{'STATISTIC':<30}"
    f"{'ACTUAL':>16}"
    f"{'NULL MEAN':>16}"
    f"{'NULL SD':>16}"
    f"{'Z':>12}"
    f"{'P(2-sided)':>16}"
)

print("-" * 100)


def print_stat(name, actual, values):

    mu = mean(values)

    if len(values) > 1:
        variance = mean(
            (x - mu) ** 2
            for x in values
        )
        sd = math.sqrt(variance)
    else:
        sd = 0.0

    if sd:
        z = (actual - mu) / sd
    else:
        z = 0.0

    p = two_sided_p(values, actual)

    print(
        f"{name:<30}"
        f"{actual:16.8f}"
        f"{mu:16.8f}"
        f"{sd:16.8f}"
        f"{z:12.4f}"
        f"{p:16.8f}"
    )


print_stat(
    "event_count",
    len(actual_events),
    null_event_counts
)

print_stat(
    "mean_ndp",
    ACTUAL_STATS["mean_ndp"],
    null_mean_ndp
)

print_stat(
    "mean_ndq",
    ACTUAL_STATS["mean_ndq"],
    null_mean_ndq
)

print_stat(
    "mean_dsum",
    ACTUAL_STATS["mean_dsum"],
    null_mean_dsum
)

print_stat(
    "mean_ddiff",
    ACTUAL_STATS["mean_ddiff"],
    null_mean_ddiff
)

actual_same_fraction = (
    ACTUAL_STATS["same_sign"]
    / len(actual_events)
    if actual_events
    else 0.0
)

print_stat(
    "same_sign_fraction",
    actual_same_fraction,
    null_same_fraction
)


# =============================================================================
# NULL QUANTILES
# =============================================================================

def quantile(values, q):

    values = sorted(values)

    if not values:
        return 0.0

    position = (
        (len(values) - 1) * q
    )

    low = math.floor(position)
    high = math.ceil(position)

    if low == high:
        return float(values[low])

    weight = position - low

    return (
        values[low] * (1 - weight)
        + values[high] * weight
    )


print()
print("=" * 100)
print("ACTUAL POSITION")
print("=" * 100)

for name, actual, values in [
    (
        "event_count",
        len(actual_events),
        null_event_counts
    ),
    (
        "mean_ndp",
        ACTUAL_STATS["mean_ndp"],
        null_mean_ndp
    ),
    (
        "mean_ndq",
        ACTUAL_STATS["mean_ndq"],
        null_mean_ndq
    ),
    (
        "mean_dsum",
        ACTUAL_STATS["mean_dsum"],
        null_mean_dsum
    ),
    (
        "mean_ddiff",
        ACTUAL_STATS["mean_ddiff"],
        null_mean_ddiff
    ),
    (
        "same_sign_fraction",
        actual_same_fraction,
        null_same_fraction
    ),
]:

    percentile = (
        sum(
            x <= actual
            for x in values
        )
        / len(values)
    )

    print(name)
    print(
        f"    actual         = "
        f"{actual:.8f}"
    )
    print(
        f"    null 2.5%      = "
        f"{quantile(values, 0.025):.8f}"
    )
    print(
        f"    null 50%       = "
        f"{quantile(values, 0.500):.8f}"
    )
    print(
        f"    null 97.5%     = "
        f"{quantile(values, 0.975):.8f}"
    )
    print(
        f"    percentile     = "
        f"{percentile:.6f}"
    )


# =============================================================================
# t DISTRIBUTION
# =============================================================================

print()
print("=" * 100)
print("t DISTRIBUTION: ACTUAL vs C-CONDITIONED NULL")
print("=" * 100)

print(
    f"{'t':>5}"
    f"{'ACT':>10}"
    f"{'NULL':>14}"
    f"{'NULL SD':>14}"
    f"{'Z':>12}"
)

for t in range(T_MIN, T_MAX + 1):

    if t == 0:
        continue

    actual_count = ACTUAL_T.get(t, 0)
    values = null_t_distributions[t]

    mu = mean(values)

    if len(values) > 1:
        sd = math.sqrt(
            mean(
                (x - mu) ** 2
                for x in values
            )
        )
    else:
        sd = 0.0

    if sd:
        z = (actual_count - mu) / sd
    else:
        z = 0.0

    print(
        f"{t:5d}"
        f"{actual_count:10d}"
        f"{mu:14.4f}"
        f"{sd:14.4f}"
        f"{z:12.4f}"
    )


# =============================================================================
# MODULAR MUTUAL INFORMATION
# =============================================================================

print()
print("=" * 100)
print("MODULAR MUTUAL INFORMATION")
print("=" * 100)

print(
    f"{'r':>5}"
    f"{'ACTUAL':>18}"
    f"{'NULL MEAN':>18}"
    f"{'NULL SD':>18}"
    f"{'Z':>12}"
    f"{'P':>14}"
)

for r in MODULI:

    actual_value = modular_mi(
        actual_events,
        r
    )

    values = null_mi[r]

    mu = mean(values)

    sd = math.sqrt(
        mean(
            (x - mu) ** 2
            for x in values
        )
    )

    if sd:
        z = (actual_value - mu) / sd
    else:
        z = 0.0

    p = two_sided_p(
        values,
        actual_value
    )

    print(
        f"{r:5d}"
        f"{actual_value:18.8f}"
        f"{mu:18.8f}"
        f"{sd:18.8f}"
        f"{z:12.4f}"
        f"{p:14.8f}"
    )


# =============================================================================
# ANCHOR DEGREE / LOCAL STRUCTURE
# =============================================================================

print()
print("=" * 100)
print("ANCHOR COLLISION DEGREE")
print("=" * 100)

degree_counts = Counter(
    events_per_anchor_actual
)

for degree in sorted(degree_counts):

    print(
        f"degree={degree:2d}"
        f" anchors={degree_counts[degree]:4d}"
    )

print()
print(
    f"mean actual degree = "
    f"{mean(events_per_anchor_actual):.6f}"
)


# =============================================================================
# C-CLASS ANALYSIS
# =============================================================================

print()
print("=" * 100)
print("C-CLASS STRUCTURE")
print("=" * 100)

c_class_sizes = []

for C in unique_C:

    candidates = candidate_cache[C]

    c_class_sizes.append(
        len(candidates)
    )

print(
    f"mean C-class size   = "
    f"{mean(c_class_sizes):.6f}"
)

print(
    f"min C-class size    = "
    f"{min(c_class_sizes)}"
)

print(
    f"max C-class size    = "
    f"{max(c_class_sizes)}"
)

print()
print(
    "candidate-set size distribution:"
)

for size, count in sorted(
    Counter(c_class_sizes).items()
):

    print(
        f"    size={size:2d} "
        f"classes={count:4d}"
    )


# =============================================================================
# EXACT C-CLASS TABLE
#
# For each observed C:
#
#   candidate pair
#   quotient k = (ab-C)/M
#   normalized location
#
# This exposes the actual geometry underlying the collision landscape.
# =============================================================================

print()
print("=" * 100)
print("OBSERVED C-CLASS CANDIDATES")
print("=" * 100)

shown = 0

for anchor in anchors:

    C = anchor["C"]

    candidates = candidate_cache[C]

    if shown >= 30:
        break

    shown += 1

    print()
    print(
        f"C={C:,} "
        f"anchor=({anchor['p']:,},{anchor['q']:,})"
    )

    anchor_n = anchor["n"]

    for a, b in candidates:

        product = a * b

        diff = product - anchor_n

        if diff % (2 * M) == 0:
            t = diff // (2 * M)
        else:
            t = None

        print(
            f"    "
            f"({a:,},{b:,}) "
            f"product={product:,} "
            f"t={t}"
        )


# =============================================================================
# FINAL ALGEBRAIC SANITY CHECK
# =============================================================================

print()
print("=" * 100)
print("FINAL SANITY CHECK")
print("=" * 100)

all_good = True

for anchor in anchors:

    candidates = candidate_cache[
        anchor["C"]
    ]

    for a, b in candidates:

        product = a * b

        if (
            product - anchor["n"]
        ) % M != 0:
            all_good = False
            break

    if not all_good:
        break


print(
    f"all candidate pairs satisfy "
    f"ab = pq (mod M): {all_good}"
)


# =============================================================================
# INTERPRETATION
# =============================================================================

print()
print("=" * 100)
print("INTERPRETATION")
print("=" * 100)

print(
    """
This experiment conditions on the exact modular class

    C = p*q (mod M).

The null does NOT replace C.
It replaces only the particular factor pair chosen as the anchor.

Therefore:

    ACTUAL
        original (p,q) anchor

    NULL
        random prime (p',q') from the same C-class

The experiment asks whether the original anchor has an unusually
high or unusual local collision structure.

This is stronger than the previous unconstrained hyperbola null,
because the residue class C is fixed for every anchor.

Important algebraic point:

    If (a,b) lies on the same modular hyperbola,

        a*b = p*q (mod M),

    then

        a*b - p*q = k*M.

Because every factor here is an odd prime and M is odd, k is even,
so

        t = k/2

    is automatically an integer.

Thus the shifted-product collision landscape is exactly the
C-conditioned prime-point geometry together with the t window.

The genuinely informative statistic is therefore whether the
ORIGINAL anchor is an unusually well-connected point of its
C-conditioned candidate set.

A significant result should survive this null.
"""
)

print()
print("=" * 100)
print("EXPERIMENT COMPLETE")
print("=" * 100)
