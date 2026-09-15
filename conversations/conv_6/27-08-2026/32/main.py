#!/usr/bin/env python3

import math
import random
from collections import Counter, defaultdict
from statistics import mean

# =============================================================================
# CRT-HYPERBOLA / MOD-M CONDITIONAL NULL EXPERIMENT
# =============================================================================

M = 111_546_435
T_MIN = -25
T_MAX = 25

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

ACTUAL_TRIALS = 300
NULL_DRAWS = 200_000

SEED = 1_511_464_998

MODULI = [3, 5, 7, 11, 13, 17, 19, 23]

random.seed(SEED)


# =============================================================================
# PRIME SIEVE
# =============================================================================

def sieve(n):
    sieve_data = bytearray(b"\x01") * (n + 1)

    sieve_data[0] = 0
    sieve_data[1] = 0

    limit = math.isqrt(n)

    for p in range(2, limit + 1):
        if sieve_data[p]:
            start = p * p
            count = ((n - start) // p) + 1
            sieve_data[start:n + 1:p] = b"\x00" * count

    return [
        x for x in range(FACTOR_MIN, n + 1)
        if sieve_data[x]
    ]


PRIMES = sieve(FACTOR_MAX)
PRIME_SET = set(PRIMES)


# =============================================================================
# MODULAR INVERSES
# =============================================================================

INVERSE = {}

for p in PRIMES:
    g = math.gcd(p, M)

    if g == 1:
        INVERSE[p] = pow(p, -1, M)


# =============================================================================
# FACTOR A NUMBER IN THE 10K-100K PRIME RANGE
# =============================================================================

def factor_in_range(n):
    """
    Return a canonical prime factor pair (p,q) with
        p*q = n
    and both factors in [FACTOR_MIN, FACTOR_MAX].

    Return None if no such pair exists.
    """

    for p in PRIMES:

        if p * p > n:
            break

        if n % p != 0:
            continue

        q = n // p

        if (
            FACTOR_MIN <= q <= FACTOR_MAX
            and q in PRIME_SET
        ):
            return (p, q)

    return None


# =============================================================================
# GENERATE ACTUAL SEMIPRIME ANCHORS
# =============================================================================

def generate_actual_anchors(count):
    """
    Generate unique semiprime anchors p*q.
    """

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

        seen.add(n)

        anchors.append(
            {
                "p": p,
                "q": q,
                "n": n,
            }
        )

    return anchors


# =============================================================================
# GENERATE ACTUAL COLLISION EVENTS
# =============================================================================

def generate_actual_events(anchors):

    events = []

    for anchor_id, anchor in enumerate(anchors, 1):

        p = anchor["p"]
        q = anchor["q"]
        n = anchor["n"]

        for t in range(T_MIN, T_MAX + 1):

            if t == 0:
                continue

            Nt = n + 2 * M * t

            if Nt <= 0:
                continue

            pair = factor_in_range(Nt)

            if pair is None:
                continue

            a, b = pair

            events.append(
                {
                    "anchor": anchor_id,
                    "p": p,
                    "q": q,
                    "n": n,
                    "t": t,
                    "Nt": Nt,
                    "a": a,
                    "b": b,
                }
            )

    return events


# =============================================================================
# EXACT MOD-M HYPERBOLA CANDIDATES
# =============================================================================

def build_hyperbola_candidates(C):
    """
    Enumerate all prime factor pairs (a,b) in the factor range satisfying

        a*b = C (mod M)

    exactly.

    Since M is much larger than the factor range, for every a there is
    one residue b = C * a^{-1} mod M.

    No rejection sampling is used.
    """

    pairs = set()

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

        x = min(a, b)
        y = max(a, b)

        pairs.add((x, y))

    return sorted(pairs)


# =============================================================================
# BASIC DATASET STATISTICS
# =============================================================================

def normalize_factor(x):
    return (
        x - FACTOR_MIN
    ) / (
        FACTOR_MAX - FACTOR_MIN
    )


def entropy(values):

    if not values:
        return 0.0

    counts = Counter(values)
    n = len(values)

    h = 0.0

    for count in counts.values():
        p = count / n
        h -= p * math.log2(p)

    return h


def mutual_information(x, y):

    if not x:
        return 0.0

    n = len(x)

    joint = Counter(zip(x, y))
    cx = Counter(x)
    cy = Counter(y)

    mi = 0.0

    for (vx, vy), count in joint.items():

        pxy = count / n
        px = cx[vx] / n
        py = cy[vy] / n

        mi += pxy * math.log2(
            pxy / (px * py)
        )

    return mi


def total_variation(x, y):

    cx = Counter(x)
    cy = Counter(y)

    nx = len(x)
    ny = len(y)

    keys = set(cx) | set(cy)

    total = 0.0

    for k in keys:

        px = cx[k] / nx
        py = cy[k] / ny

        total += abs(px - py)

    return 0.5 * total


def same_sign_fraction(pairs_a, pairs_b):

    if not pairs_a:
        return 0.0

    same = 0

    for (a, b), (p, q) in zip(pairs_a, pairs_b):

        dp = a - p
        dq = b - q

        if dp == 0 or dq == 0:
            continue

        if (dp > 0 and dq > 0) or (dp < 0 and dq < 0):
            same += 1

    return same / len(pairs_a)


# =============================================================================
# BUILD DATASET STATISTICS
# =============================================================================

def dataset_stats(pairs):

    result = {}

    if not pairs:
        return result

    A = [x for x, y in pairs]
    B = [y for x, y in pairs]

    result["mean_A"] = mean(
        normalize_factor(x)
        for x in A
    )

    result["mean_B"] = mean(
        normalize_factor(y)
        for y in B
    )

    for r in MODULI:

        rx = [x % r for x in A]
        ry = [y % r for y in B]

        result[f"mi_{r}"] = mutual_information(
            rx,
            ry
        )

        result[f"joint_entropy_{r}"] = entropy(
            list(zip(rx, ry))
        )

    signatures = []

    for a, b in pairs:

        signature = (
            tuple(a % r for r in MODULI) +
            tuple(b % r for r in MODULI)
        )

        signatures.append(signature)

    result["signature_entropy"] = entropy(signatures)

    result["signature_unique"] = len(
        set(signatures)
    )

    result["signature_concentration"] = (
        max(Counter(signatures).values())
        / len(signatures)
    )

    return result


# =============================================================================
# HEADER
# =============================================================================

print("=" * 100)
print("CRT-HYPERBOLA / MOD-M CONDITIONAL NULL EXPERIMENT")
print("=" * 100)

print(f"M                 = {M:,}")
print(f"2M                = {2 * M:,}")
print(
    f"factor range      = "
    f"{FACTOR_MIN:,} - {FACTOR_MAX:,}"
)
print(
    f"T range           = "
    f"[{T_MIN}, {T_MAX}]"
)
print(f"actual trials     = {ACTUAL_TRIALS:,}")
print(f"null draws        = {NULL_DRAWS:,}")
print(f"prime factors     = {len(PRIMES):,}")
print(f"seed              = {SEED:,}")


# =============================================================================
# ACTUAL DATA
# =============================================================================

print()
print("=" * 100)
print("GENERATING ACTUAL DATA")
print("=" * 100)

anchors = generate_actual_anchors(
    ACTUAL_TRIALS
)

actual_events = generate_actual_events(
    anchors
)

print(
    f"actual anchors = {len(anchors):,}"
)

print(
    f"actual events  = {len(actual_events):,}"
)


# =============================================================================
# MOD-M CONSTANTS
# =============================================================================

actual_C = []

for e in actual_events:

    C = (
        e["p"] * e["q"]
    ) % M

    actual_C.append(C)


unique_C = sorted(set(actual_C))

print()
print("=" * 100)
print("MOD-M LANDSCAPE")
print("=" * 100)

print(
    f"collision events       = {len(actual_events):,}"
)

print(
    f"unique C values        = {len(unique_C):,}"
)


# =============================================================================
# EXACT HYPERBOLA CANDIDATE SETS
# =============================================================================

print()
print("=" * 100)
print("BUILDING EXACT HYPERBOLA CANDIDATE SETS")
print("=" * 100)

candidate_cache = {}

for index, C in enumerate(unique_C, 1):

    candidates = build_hyperbola_candidates(C)

    candidate_cache[C] = candidates

    if (
        index % 25 == 0
        or index == len(unique_C)
    ):
        print(
            f"C {index:4d}/"
            f"{len(unique_C):4d}"
            f"   candidates={len(candidates):4d}"
        )


usable_C = [
    C for C in unique_C
    if candidate_cache[C]
]


print()
print(
    f"C with candidates = "
    f"{len(usable_C):,}"
)

print(
    f"C without candidates = "
    f"{len(unique_C) - len(usable_C):,}"
)


# =============================================================================
# VERIFY ALL ACTUAL EVENTS AGAINST MOD-M HYPERBOLA
# =============================================================================

print()
print("=" * 100)
print("ACTUAL EVENT MOD-M VERIFICATION")
print("=" * 100)

identity_failures = 0
candidate_membership_failures = 0

for e in actual_events:

    a = e["a"]
    b = e["b"]

    C = (
        e["p"] * e["q"]
    ) % M

    if (a * b - C) % M != 0:
        identity_failures += 1

    if (
        C not in candidate_cache
        or (a, b) not in candidate_cache[C]
    ):
        candidate_membership_failures += 1


print(
    f"identity failures       = "
    f"{identity_failures}"
)

print(
    f"candidate membership "
    f"failures                = "
    f"{candidate_membership_failures}"
)


# =============================================================================
# FILTER ACTUAL EVENTS TO ONES WITH A VALID NULL
# =============================================================================

usable_actual_events = []

for e in actual_events:

    C = (
        e["p"] * e["q"]
    ) % M

    if candidate_cache[C]:
        usable_actual_events.append(e)


print()
print(
    f"usable actual events = "
    f"{len(usable_actual_events):,}"
)


# =============================================================================
# CONDITIONAL HYPERBOLA NULL
# =============================================================================

print()
print("=" * 100)
print("SAMPLING CONDITIONAL HYPERBOLA NULL")
print("=" * 100)

hyperbola_pairs = []

for i in range(NULL_DRAWS):

    e = random.choice(
        usable_actual_events
    )

    C = (
        e["p"] * e["q"]
    ) % M

    candidates = candidate_cache[C]

    a, b = random.choice(candidates)

    hyperbola_pairs.append(
        {
            "a": a,
            "b": b,
            "C": C,
        }
    )

    if (
        (i + 1) % 25_000 == 0
        or i + 1 == NULL_DRAWS
    ):
        print(
            f"draw "
            f"{i + 1:7,d}/"
            f"{NULL_DRAWS:,}"
        )


# =============================================================================
# RANDOM PRIME BASELINE
# =============================================================================

print()
print("=" * 100)
print("SAMPLING RANDOM PRIME BASELINE")
print("=" * 100)

random_pairs = []

while len(random_pairs) < NULL_DRAWS:

    a = random.choice(PRIMES)
    b = random.choice(PRIMES)

    if a == b:
        continue

    if a > b:
        a, b = b, a

    random_pairs.append((a, b))


# =============================================================================
# ACTUAL PAIRS
# =============================================================================

actual_pairs = [
    (e["a"], e["b"])
    for e in usable_actual_events
]

hyper_pairs = [
    (e["a"], e["b"])
    for e in hyperbola_pairs
]


# =============================================================================
# DATASET STATISTICS
# =============================================================================

ACT = dataset_stats(actual_pairs)
HYP = dataset_stats(hyper_pairs)
RND = dataset_stats(random_pairs)


# =============================================================================
# SUMMARY STATISTICS
# =============================================================================

print()
print("=" * 100)
print("FACTOR-SIZE STATISTICS")
print("=" * 100)

print(
    f"{'STATISTIC':<32}"
    f"{'ACTUAL':>18}"
    f"{'HYPERBOLA':>18}"
    f"{'RANDOM':>18}"
)

print("-" * 100)

for key in [
    "mean_A",
    "mean_B",
    "signature_entropy",
    "signature_unique",
    "signature_concentration",
]:

    print(
        f"{key:<32}"
        f"{ACT[key]:18.8f}"
        f"{HYP[key]:18.8f}"
        f"{RND[key]:18.8f}"
    )


# =============================================================================
# MUTUAL INFORMATION
# =============================================================================

print()
print("=" * 100)
print("MODULAR MUTUAL INFORMATION")
print("=" * 100)

print(
    f"{'r':>4}"
    f"{'ACTUAL':>18}"
    f"{'HYPERBOLA':>18}"
    f"{'RANDOM':>18}"
)

for r in MODULI:

    print(
        f"{r:4d}"
        f"{ACT[f'mi_{r}']:18.8f}"
        f"{HYP[f'mi_{r}']:18.8f}"
        f"{RND[f'mi_{r}']:18.8f}"
    )


# =============================================================================
# TOTAL-VARIATION DISTANCE
# =============================================================================

print()
print("=" * 100)
print("JOINT RESIDUE DISTRIBUTION")
print("=" * 100)

print(
    f"{'r':>4}"
    f"{'TV(ACT,HYP)':>20}"
    f"{'TV(ACT,RND)':>20}"
)


for r in MODULI:

    actual_residue_pairs = [
        (a % r, b % r)
        for a, b in actual_pairs
    ]

    hyper_residue_pairs = [
        (a % r, b % r)
        for a, b in hyper_pairs
    ]

    random_residue_pairs = [
        (a % r, b % r)
        for a, b in random_pairs
    ]

    tv_h = total_variation(
        actual_residue_pairs,
        hyper_residue_pairs
    )

    tv_r = total_variation(
        actual_residue_pairs,
        random_residue_pairs
    )

    print(
        f"{r:4d}"
        f"{tv_h:20.8f}"
        f"{tv_r:20.8f}"
    )


# =============================================================================
# RESIDUE PAIR COUNTS
# =============================================================================

print()
print("=" * 100)
print("RESIDUE PAIR COUNTS: ACTUAL vs HYPERBOLA")
print("=" * 100)

for r in MODULI:

    actual_counts = Counter(
        (a % r, b % r)
        for a, b in actual_pairs
    )

    hyper_counts = Counter(
        (a % r, b % r)
        for a, b in hyper_pairs
    )

    print()
    print(f"r = {r}")

    keys = sorted(
        set(actual_counts) |
        set(hyper_counts)
    )

    for key in keys:

        print(
            f"  {key}: "
            f"actual={actual_counts[key]:6d} "
            f"hyperbola={hyper_counts[key]:6d}"
        )


# =============================================================================
# FACTOR GEOMETRY
#
# We compare each actual factor pair to the original anchor pair.
# =============================================================================

print()
print("=" * 100)
print("ACTUAL FACTOR GEOMETRY")
print("=" * 100)

by_t = defaultdict(list)

for e in usable_actual_events:
    by_t[e["t"]].append(e)


for t in sorted(by_t):

    group = by_t[t]

    dps = []
    dqs = []
    sums = []
    diffs = []

    for e in group:

        p = e["p"]
        q = e["q"]

        a = e["a"]
        b = e["b"]

        dp = a - p
        dq = b - q

        dps.append(dp)
        dqs.append(dq)

        sums.append(dp + dq)
        diffs.append(dp - dq)

    mean_dp = mean(dps)
    mean_dq = mean(dqs)
    mean_sum = mean(sums)
    mean_diff = mean(diffs)

    print(
        f"t={t:+3d} "
        f"N={len(group):4d} "
        f"mean_dp={mean_dp:+11.2f} "
        f"mean_dq={mean_dq:+11.2f} "
        f"mean_dsum={mean_sum:+12.2f} "
        f"mean_ddiff={mean_diff:+12.2f}"
    )


# =============================================================================
# DIRECT PRODUCT IDENTITY
# =============================================================================

print()
print("=" * 100)
print("PRODUCT IDENTITY")
print("=" * 100)

product_failures = 0
t_failures = 0

for e in usable_actual_events:

    p = e["p"]
    q = e["q"]

    a = e["a"]
    b = e["b"]

    t = e["t"]

    n = e["n"]
    Nt = e["Nt"]

    if a * b != Nt:
        product_failures += 1

    reconstructed_t = (
        (a * b) - (p * q)
    ) // (2 * M)

    if reconstructed_t != t:
        t_failures += 1


print(
    f"a*b = N_t failures = "
    f"{product_failures}"
)

print(
    f"t reconstruction failures = "
    f"{t_failures}"
)


# =============================================================================
# FINAL INTERPRETATION
# =============================================================================

print()
print("=" * 100)
print("INTERPRETATION")
print("=" * 100)

print(
    """
The three datasets have different meanings.

ACTUAL
    The factors of
        N_t = p*q + 2*M*t

HYPERBOLA NULL
    Prime pairs (a,b) satisfying
        a*b = p*q (mod M)

    while preserving the observed modulo-M anchor C.

RANDOM
    Unrestricted random prime pairs from the same factor range.

The important comparison is:

    ACTUAL vs HYPERBOLA

not merely:

    ACTUAL vs RANDOM.

A small ACTUAL/HYPERBOLA residue-distribution distance means
the observed residue structure is already explained by the
mod-M hyperbola and the prime/range constraints.

A persistent ACTUAL/HYPERBOLA difference would be evidence
for additional structure beyond that congruence.
"""
)

print()
print("=" * 100)
print("EXPERIMENT COMPLETE")
print("=" * 100)