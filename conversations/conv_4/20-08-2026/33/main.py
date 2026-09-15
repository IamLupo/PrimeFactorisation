#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 509
==============================================================================
DIVISOR-SPECTRUM / CANONICAL HISTORICAL-PAIR SELECTION SEARCH
==============================================================================

CORE QUESTION

Experiment 508 established:

    K = 4v - 1

and

    K = N       if N == 3 mod 4
    K = 3N      if N == 1 mod 4.

The historical coordinates satisfy

    A0 = 2y - 2x + 3
    B0 = 2y + 2x - 3

with

    A0*B0 = K.

Therefore every divisor pair

    d * e = K

produces a valid algebraic (x,y):

    y = (d + e)/4
    x = (e - d + 6)/4

whenever the required integrality conditions hold.

The research question is:

    Is there an N-only observable that selects the particular
    divisor pair corresponding to the desired p,q?

We search exact structural signatures involving:

    N
    b
    v
    g
    K
    sqrt(N)
    sqrt(K)
    divisor sum d+e
    divisor difference e-d
    parity
    residue classes
    proximity to sqrt(K)
    proximity to sqrt(N)
    modular signatures

We explicitly distinguish:

    TRIVIAL identities
    factor-selection identities
    empirical coincidences.

No floating-point arithmetic is used for exact tests.
==============================================================================

"""

from __future__ import annotations

import math
from itertools import combinations


# ============================================================================
# TEST DATA
# ============================================================================

CASES = [
    (50387, 282589),
    (1009, 10007),
    (10007, 1000003),
    (100003, 100019),
    (2000003, 3000017),
    (50021, 50047),
    (300007, 900001),
    (3, 5),
    (5, 13),
    (13, 17),
    (3, 7),
    (7, 11),
]


# ============================================================================
# BASIC N-ONLY CHANNEL
# ============================================================================

def b_value(n: int) -> int:
    return n // 2


def v_value(n: int) -> int:
    b = b_value(n)
    return (b * b) % n


def g_value(n: int) -> int:
    b = b_value(n)
    return b // 2 + 1


def K_value(n: int) -> int:
    return 4 * v_value(n) - 1


# ============================================================================
# INTEGER HELPERS
# ============================================================================

def isqrt(n: int) -> int:
    return math.isqrt(n)


def divisors(n: int) -> list[int]:
    """
    Exact positive divisors.
    """
    out = []

    r = isqrt(n)

    for d in range(1, r + 1):
        if n % d == 0:
            out.append(d)

            other = n // d

            if other != d:
                out.append(other)

    return sorted(out)


def factor_pairs(n: int):
    """
    Unique unordered positive factor pairs d <= e.
    """
    pairs = []

    for d in divisors(n):
        e = n // d

        if d > e:
            continue

        pairs.append((d, e))

    return pairs


# ============================================================================
# HISTORICAL COORDINATES
# ============================================================================

def historical_xy_from_pair(d: int, e: int):
    """
    A0=d, B0=e.

        A0 = 2y - 2x + 3
        B0 = 2y + 2x - 3

    Hence

        y = (A0+B0)/4
        x = (B0-A0+6)/4.
    """

    if (d + e) % 4 != 0:
        return None

    if (e - d + 6) % 4 != 0:
        return None

    y = (d + e) // 4
    x = (e - d + 6) // 4

    return x, y


# ============================================================================
# TRUE HISTORICAL TARGET
# ============================================================================

def true_pair_from_case(p: int, q: int):
    n = p * q

    if n % 4 == 3:
        candidates = [
            (min(p, q), max(p, q)),
        ]
    else:
        candidates = [
            (min(p, 3 * q), max(p, 3 * q)),
            (min(q, 3 * p), max(q, 3 * p)),
        ]

    K = K_value(n)

    for pair in candidates:
        if pair[0] * pair[1] == K:
            return pair

    raise RuntimeError(
        f"No historical pair found for ({p},{q})"
    )


# ============================================================================
# CANDIDATE FEATURES
# ============================================================================

def candidate_features(n: int, d: int, e: int):
    K = K_value(n)

    x_y = historical_xy_from_pair(d, e)

    x = None
    y = None

    if x_y is not None:
        x, y = x_y

    s = d + e
    diff = e - d

    floor_sqrt_n = isqrt(n)
    ceil_sqrt_n = floor_sqrt_n if floor_sqrt_n ** 2 == n else floor_sqrt_n + 1

    floor_sqrt_k = isqrt(K)
    ceil_sqrt_k = floor_sqrt_k if floor_sqrt_k ** 2 == K else floor_sqrt_k + 1

    return {
        "d": d,
        "e": e,
        "sum": s,
        "diff": diff,
        "sum2": s * s,
        "diff2": diff * diff,
        "product": d * e,
        "x": x,
        "y": y,
        "sum_over_4": s // 4 if s % 4 == 0 else None,
        "diff_over_2": diff // 2 if diff % 2 == 0 else None,

        # distance from sqrt(K)
        "sqrtK_distance": e - ceil_sqrt_k,

        # distance from sqrt(N)
        "sqrtN_distance_d": abs(d - ceil_sqrt_n),
        "sqrtN_distance_e": abs(e - ceil_sqrt_n),

        # parity
        "d_mod2": d % 2,
        "e_mod2": e % 2,

        # mod 4
        "d_mod4": d % 4,
        "e_mod4": e % 4,

        # mod 3
        "d_mod3": d % 3,
        "e_mod3": e % 3,

        # relation to N-only historical quantities
        "d_minus_v": d - v_value(n),
        "e_minus_v": e - v_value(n),
        "d_minus_g": d - g_value(n),
        "e_minus_g": e - g_value(n),

        "sum_minus_b": s - b_value(n),
        "diff_minus_b": diff - b_value(n),

        "sum_mod_n": s % n,
        "diff_mod_n": diff % n,

        "sum_mod_k": s % K,
        "diff_mod_k": diff % K,
    }


# ============================================================================
# EXACT PROPERTY TESTS
# ============================================================================

def structural_signature(n: int, d: int, e: int):
    """
    A collection of simple exact predicates.

    These are deliberately N-only once d,e are supplied as a candidate.
    """

    f = candidate_features(n, d, e)

    K = K_value(n)
    b = b_value(n)
    v = v_value(n)
    g = g_value(n)

    return {
        "d_is_prime_factor_of_N":
            n % d == 0,

        "e_is_prime_factor_of_N":
            n % e == 0,

        "one_factor_N":
            (n % d == 0 or n % e == 0),

        "contains_factor_3":
            (d % 3 == 0 or e % 3 == 0),

        "sum_div4":
            (d + e) % 4 == 0,

        "diff_plus6_div4":
            (e - d + 6) % 4 == 0,

        "sum_equals_4y":
            f["y"] is not None,

        "v_is_midpoint":
            abs(d - v) == abs(e - v),

        "g_between":
            min(d, e) <= g <= max(d, e),

        "b_between":
            min(d, e) <= b <= max(d, e),

        "sum_equals_4b":
            d + e == 4 * b,

        "diff_even":
            (e - d) % 2 == 0,

        "same_mod4":
            d % 4 == e % 4,

        "opposite_mod4":
            d % 4 != e % 4,

        "nearest_sqrtK_upper":
            e == isqrt(K) or e == isqrt(K) + 1,

        "nearest_sqrtN_upper":
            e == isqrt(n) or e == isqrt(n) + 1,

        "smallest_nontrivial_divisor":
            d != 1 and d == min(
                q for q in divisors(K)
                if q > 1
            ),

        "largest_nontrivial_divisor":
            e != K and e == max(
                q for q in divisors(K)
                if q < K
            ),

        "candidate_x_nonnegative":
            f["x"] is not None and f["x"] >= 0,

        "candidate_y_nonnegative":
            f["y"] is not None and f["y"] >= 0,

        "candidate_y_gt_x":
            f["x"] is not None and f["y"] > f["x"],

        "candidate_y_gt_2x":
            f["x"] is not None and f["y"] > 2 * f["x"],

        "candidate_y_lt_2x":
            f["x"] is not None and f["y"] < 2 * f["x"],
    }


# ============================================================================
# PRINT HEADER
# ============================================================================

print("=" * 78)
print("EXPERIMENT 509 START")
print("=" * 78)
print("DIVISOR-SPECTRUM / CANONICAL HISTORICAL-PAIR SELECTION SEARCH")
print("=" * 78)


# ============================================================================
# [1] KERNEL AUDIT
# ============================================================================

print()
print("[1] HISTORICAL KERNEL")
print("-" * 78)

failures = 0

for p, q in CASES:

    n = p * q

    v = v_value(n)
    K = K_value(n)

    expected = n if n % 4 == 3 else 3 * n

    ok = K == expected

    print(
        f"  ({p},{q}) "
        f"N={n} "
        f"v={v} "
        f"K=4v-1={K} "
        f"expected={expected} "
        f"PASS={ok}"
    )

    if not ok:
        failures += 1

print(f"  KERNEL FAILURES = {failures}")


# ============================================================================
# [2] DIVISOR SPECTRUM
# ============================================================================

print()
print("[2] DIVISOR SPECTRUM OF K=4v-1")
print("-" * 78)

for p, q in CASES:

    n = p * q
    K = K_value(n)

    pairs = factor_pairs(K)

    print(
        f"\n  ({p},{q}) "
        f"N={n} "
        f"K={K} "
        f"factor-pair count={len(pairs)}"
    )

    for d, e in pairs:

        xy = historical_xy_from_pair(d, e)

        if xy is None:
            continue

        x, y = xy

        marker = ""

        true_d, true_e = true_pair_from_case(p, q)

        if (d, e) == (true_d, true_e):
            marker = "  <-- TARGET"

        print(
            f"      d={d:<18} "
            f"e={e:<18} "
            f"d+e={d+e:<18} "
            f"e-d={e-d:<18} "
            f"x={x:<10} "
            f"y={y:<10}"
            f"{marker}"
        )


# ============================================================================
# [3] TARGET PAIR FEATURES
# ============================================================================

print()
print("[3] TARGET PAIR N-ONLY SIGNATURE")
print("-" * 78)

for p, q in CASES:

    n = p * q
    d, e = true_pair_from_case(p, q)

    f = candidate_features(n, d, e)

    print(
        f"  ({p},{q}) "
        f"d={d} "
        f"e={e} "
        f"sum={f['sum']} "
        f"diff={f['diff']} "
        f"sqrtK-distance={f['sqrtK_distance']} "
        f"x={f['x']} "
        f"y={f['y']}"
    )


# ============================================================================
# [4] PROPERTY UNIQUENESS SEARCH
# ============================================================================

print()
print("[4] SINGLE-PREDICATE TARGET SELECTION")
print("-" * 78)

property_names = [
    "d_is_prime_factor_of_N",
    "e_is_prime_factor_of_N",
    "one_factor_N",
    "contains_factor_3",
    "sum_div4",
    "diff_plus6_div4",
    "sum_equals_4y",
    "v_is_midpoint",
    "g_between",
    "b_between",
    "sum_equals_4b",
    "diff_even",
    "same_mod4",
    "opposite_mod4",
    "nearest_sqrtK_upper",
    "nearest_sqrtN_upper",
    "smallest_nontrivial_divisor",
    "largest_nontrivial_divisor",
    "candidate_x_nonnegative",
    "candidate_y_nonnegative",
    "candidate_y_gt_x",
    "candidate_y_gt_2x",
    "candidate_y_lt_2x",
]

universal_unique_hits = []

for prop in property_names:

    unique_all = True
    details = []

    for p, q in CASES:

        n = p * q
        K = K_value(n)

        target = true_pair_from_case(p, q)

        pairs = factor_pairs(K)

        matching = []

        for d, e in pairs:

            sig = structural_signature(n, d, e)

            if sig[prop]:
                matching.append((d, e))

        if len(matching) != 1 or matching[0] != target:
            unique_all = False

        details.append(
            (
                (p, q),
                target,
                matching,
            )
        )

    if unique_all:
        universal_unique_hits.append(prop)

        print(
            f"  UNIVERSAL UNIQUE SELECTOR: {prop}"
        )

print(
    f"  universal selectors found = "
    f"{len(universal_unique_hits)}"
)

if not universal_unique_hits:
    print("  NONE")


# ============================================================================
# [5] PREDICATE INTERSECTION SEARCH
# ============================================================================

print()
print("[5] SMALL PREDICATE INTERSECTION SEARCH")
print("-" * 78)

pair_hits = []
triple_hits = []

for a, b in combinations(property_names, 2):

    success = True

    for p, q in CASES:

        n = p * q
        K = K_value(n)

        target = true_pair_from_case(p, q)

        matching = []

        for d, e in factor_pairs(K):

            sig = structural_signature(n, d, e)

            if sig[a] and sig[b]:
                matching.append((d, e))

        if len(matching) != 1 or matching[0] != target:
            success = False
            break

    if success:
        pair_hits.append((a, b))
        print(
            f"  PAIR SELECTOR: {a} AND {b}"
        )


for a, b, c in combinations(property_names, 3):

    success = True

    for p, q in CASES:

        n = p * q
        K = K_value(n)

        target = true_pair_from_case(p, q)

        matching = []

        for d, e in factor_pairs(K):

            sig = structural_signature(n, d, e)

            if sig[a] and sig[b] and sig[c]:
                matching.append((d, e))

        if len(matching) != 1 or matching[0] != target:
            success = False
            break

    if success:
        triple_hits.append((a, b, c))
        print(
            f"  TRIPLE SELECTOR: {a} AND {b} AND {c}"
        )

print(
    f"  pair selectors   = {len(pair_hits)}"
)

print(
    f"  triple selectors = {len(triple_hits)}"
)


# ============================================================================
# [6] TARGET SUM / DIFFERENCE SIGNATURES
# ============================================================================

print()
print("[6] FACTOR-PAIR SUM / DIFFERENCE SEARCH")
print("-" * 78)

sum_hits = []
diff_hits = []

for p, q in CASES:

    n = p * q
    K = K_value(n)

    target_d, target_e = true_pair_from_case(p, q)

    target_sum = target_d + target_e
    target_diff = target_e - target_d

    sums = []
    diffs = []

    for d, e in factor_pairs(K):

        sums.append((d + e, (d, e)))
        diffs.append((e - d, (d, e)))

    sum_matches = [
        pair for value, pair in sums
        if value == target_sum
    ]

    diff_matches = [
        pair for value, pair in diffs
        if value == target_diff
    ]

    sum_hits.append(len(sum_matches) == 1)
    diff_hits.append(len(diff_matches) == 1)

    print(
        f"  ({p},{q}) "
        f"target_sum={target_sum} "
        f"sum_unique={len(sum_matches) == 1} "
        f"target_diff={target_diff} "
        f"diff_unique={len(diff_matches) == 1}"
    )


# ============================================================================
# [7] NEAREST-SQUARE SELECTION
# ============================================================================

print()
print("[7] NEAREST-SQUARE / NEAREST-SQRT SEARCH")
print("-" * 78)

for p, q in CASES:

    n = p * q
    K = K_value(n)

    target = true_pair_from_case(p, q)

    pairs = factor_pairs(K)

    ranked = []

    rootK = isqrt(K)

    for d, e in pairs:

        distance = abs(e - rootK)

        ranked.append(
            (distance, d, e)
        )

    ranked.sort()

    best = ranked[0]

    print(
        f"  ({p},{q}) "
        f"target={target} "
        f"nearest-sqrtK={best[1:]}"
        f" "
        f"TARGET_SELECTED={best[1:] == target}"
    )


# ============================================================================
# [8] PRODUCT-PAIR GAP SPECTRUM
# ============================================================================

print()
print("[8] GAP-SQUARE SPECTRUM")
print("-" * 78)

for p, q in CASES:

    n = p * q
    K = K_value(n)

    target = true_pair_from_case(p, q)
    target_gap2 = (target[1] - target[0]) ** 2

    pairs = factor_pairs(K)

    spectrum = []

    for d, e in pairs:

        gap2 = (e - d) ** 2

        spectrum.append(
            (gap2, d, e)
        )

    spectrum.sort()

    rank = None

    for i, (_, d, e) in enumerate(spectrum, start=1):
        if (d, e) == target:
            rank = i
            break

    print(
        f"  ({p},{q}) "
        f"target_gap2={target_gap2} "
        f"gap-rank={rank}/{len(spectrum)}"
    )


# ============================================================================
# [9] HISTORICAL x,y SIGNATURE
# ============================================================================

print()
print("[9] HISTORICAL x,y TARGET SELECTION")
print("-" * 78)

for p, q in CASES:

    n = p * q
    K = K_value(n)

    target = true_pair_from_case(p, q)
    tx, ty = historical_xy_from_pair(*target)

    candidates = []

    for d, e in factor_pairs(K):

        xy = historical_xy_from_pair(d, e)

        if xy is None:
            continue

        x, y = xy

        # Candidate score based ONLY on generic geometric properties.
        score = (
            abs(x),
            abs(y),
            abs(y - x),
        )

        candidates.append(
            (score, d, e, x, y)
        )

    candidates.sort()

    best = candidates[0]

    print(
        f"  ({p},{q}) "
        f"target=(x={tx},y={ty}) "
        f"min-|x| candidate="
        f"(x={best[3]},y={best[4]}) "
        f"selected={best[1:3] == target}"
    )


# ============================================================================
# [10] SAME-N CONTROL
# ============================================================================

print()
print("[10] SAME-N CONTROL")
print("-" * 78)

controls = {
    12: [(2, 6), (3, 4)],
    18: [(2, 9), (3, 6)],
    20: [(2, 10), (4, 5)],
    30: [(3, 10), (5, 6)],
}

for n, pairs in controls.items():

    K = 4 * v_value(n) - 1

    print(
        f"  N={n} "
        f"v={v_value(n)} "
        f"K={K}"
    )

    for p, q in pairs:

        target = None

        for pair in factor_pairs(K):
            normalized = normalize = None

            # For N==3 mod 4, K=N.
            # For N==1 mod 4, K=3N.
            if n % 4 == 3:
                if set(pair) == {p, q}:
                    target = pair
                    break
            else:
                candidate1 = tuple(sorted((p, 3 * q)))
                candidate2 = tuple(sorted((q, 3 * p)))

                if pair in (candidate1, candidate2):
                    target = pair
                    break

        print(
            f"      ({p},{q}) "
            f"historical pair={target}"
        )


# ============================================================================
# [11] IMPORTANT ALGEBRA
# ============================================================================

print()
print("[11] EXACT DIVISOR-SELECTION ALGEBRA")
print("-" * 78)

print(
    r"""
Let

    K = 4v - 1.

Every divisor d | K gives

    e = K/d.

The historical coordinates are

    y = (d+e)/4
    x = (e-d+6)/4.

Then

    A0 = d
    B0 = e.

The modern symmetric quantities are therefore

    S_candidate = d+e
    Delta_candidate = (e-d)^2

up to the scale-3 normalization.

This gives a precise reformulation:

    N
      -> K=4v-1
      -> divisor spectrum {d, K/d}
      -> candidate S, Delta
      -> candidate factor pair.

The remaining question is therefore not

    "Can the conic be transformed into factorization?"

It can.

The question is:

    "Does the existing KAPPA/homogeneous-layer construction
     provide a canonical rule selecting one divisor pair
     from this spectrum?"

That is the object this experiment is testing.
"""
)


# ============================================================================
# [12] FINAL STATUS
# ============================================================================

print()
print("[12] EXPERIMENT STATUS")
print("-" * 78)

print(
    f"  universal single selectors = "
    f"{len(universal_unique_hits)}"
)

print(
    f"  pair selectors             = "
    f"{len(pair_hits)}"
)

print(
    f"  triple selectors           = "
    f"{len(triple_hits)}"
)

print(
    f"  unique sum observations    = "
    f"{sum(sum_hits)}/{len(sum_hits)}"
)

print(
    f"  unique difference          = "
    f"{sum(diff_hits)}/{len(diff_hits)}"
)

print()
print(
    "NEXT RESEARCH TARGET:"
)

print(
    """
    Compare the divisor-spectrum features with the actual
    N-only homogeneous-layer observables from the previous
    experiments.

    In particular, search for an upstream quantity whose value
    equals one of:

        d+K/d
        (d-K/d)^2
        d+K/d - 4*sqrt(K)
        d+K/d - 4*b
        d-K/d
        (d+1)(K/d+1)

    for the selected historical divisor pair.

    A successful exact relation would be substantially more
    interesting than another downstream KAPPA identity because
    it would identify the missing divisor-selection interface.
    """
)

print("=" * 78)
print("EXPERIMENT 509 FINISHED")
print("=" * 78)
