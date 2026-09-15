import math
import random
import statistics
from itertools import combinations

# ======================================================================================
# CRT MULTI-MODULUS FINGERPRINT / CANDIDATE-REDUCTION EXPERIMENT
# ======================================================================================
#
# Question:
#
#   Does adding more independent prime-modulus fingerprints
#
#       n mod r1
#       n mod r2
#       ...
#
#   progressively reduce the possible prime-factor pairs?
#
# Since
#
#   M = 3*5*7*11*13*17*19*23
#
# the eight prime moduli below form the complete CRT decomposition of M.
#
# For a divisor R of M:
#
#       a*b = C (mod R)
#
# defines the candidate factor pairs compatible with the observed
# product fingerprint modulo R.
#
# We measure:
#
#   1. number of compatible unordered prime pairs
#   2. reduction factor relative to no modulus
#   3. singleton probability
#   4. average / median ambiguity
#   5. exact agreement with the previously measured C-classes at R=M
#
# We also test several non-nested modulus combinations so the result
# does not depend only on one arbitrary CRT ordering.
#
# ======================================================================================

M = 111_546_435
FACTOR_MIN = 10_000
FACTOR_MAX = 100_000
T_MIN = -25
T_MAX = 25

ANCHORS = 300
SEED = 1_511_464_998

MODULI = [3, 5, 7, 11, 13, 17, 19, 23]

assert math.prod(MODULI) == M

rng = random.Random(SEED)


# ======================================================================================
# PRIME GENERATION
# ======================================================================================

def sieve(limit):
    is_prime = bytearray(b"\x01") * (limit + 1)
    is_prime[0:2] = b"\x00\x00"

    for p in range(2, int(limit ** 0.5) + 1):
        if is_prime[p]:
            start = p * p
            is_prime[start:limit + 1:p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    return [p for p in range(limit + 1) if is_prime[p]]


ALL_PRIMES = sieve(FACTOR_MAX)
FACTOR_PRIMES = [
    p for p in ALL_PRIMES
    if FACTOR_MIN <= p <= FACTOR_MAX
]

PRIME_SET = set(FACTOR_PRIMES)

print("=" * 100)
print("CRT MULTI-MODULUS FINGERPRINT / CANDIDATE-REDUCTION EXPERIMENT")
print("=" * 100)
print(f"M                    = {M:,}")
print(f"factor range         = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
print(f"T range              = [{T_MIN}, {T_MAX}]")
print(f"actual anchors       = {ANCHORS}")
print(f"factor primes        = {len(FACTOR_PRIMES):,}")
print(f"prime CRT factors    = {MODULI}")
print(f"seed                 = {SEED:,}")


# ======================================================================================
# PRECOMPUTE INVERSES MODULO M
# ======================================================================================
#
# Every factor prime is coprime to M because all factor primes are > 23.
#
# If
#
#     a*b = C (mod R)
#
# then
#
#     b = C*a^(-1) (mod R).
#
# Since R divides M, an inverse modulo M can safely be reduced modulo R.
# ======================================================================================

INV_M = {
    p: pow(p, -1, M)
    for p in FACTOR_PRIMES
}


# ======================================================================================
# GENERATE UNIQUE ACTUAL SEMIPRIME ANCHORS
# ======================================================================================

anchors = []
seen = set()

while len(anchors) < ANCHORS:
    p, q = rng.sample(FACTOR_PRIMES, 2)
    if p > q:
        p, q = q, p

    if (p, q) in seen:
        continue

    seen.add((p, q))

    n = p * q
    C = n % M

    anchors.append({
        "p": p,
        "q": q,
        "n": n,
        "C": C,
    })


print()
print("=" * 100)
print("GENERATING ACTUAL ANCHORS")
print("=" * 100)
print(f"anchors = {len(anchors)}")


# ======================================================================================
# BUILD EXACT C-CLASSES
# ======================================================================================
#
# Because q < M, for a fixed p and C:
#
#     q = C * p^(-1) mod M
#
# produces at most one q in the factor range.
#
# This makes exact C-class construction very cheap.
# ======================================================================================

def build_c_class(C):
    result = set()

    for a in FACTOR_PRIMES:
        q = (C * INV_M[a]) % M

        if q < FACTOR_MIN or q > FACTOR_MAX:
            continue

        if q not in PRIME_SET:
            continue

        x, y = (a, q) if a <= q else (q, a)
        result.add((x, y))

    return sorted(result)


print()
print("=" * 100)
print("BUILDING EXACT C-CLASSES")
print("=" * 100)

c_classes = {}

for i, anchor in enumerate(anchors, 1):
    C = anchor["C"]

    if C not in c_classes:
        c_classes[C] = build_c_class(C)

    if i % 25 == 0 or i == len(anchors):
        print(
            f"C {i:3d}/{len(anchors):3d} "
            f"class_size={len(c_classes[C]):2d}"
        )


# ======================================================================================
# VERIFY ANCHOR MEMBERSHIP
# ======================================================================================

membership_failures = 0
class_sizes = []

for anchor in anchors:
    pair = (anchor["p"], anchor["q"])
    cls = c_classes[anchor["C"]]

    if pair not in cls:
        membership_failures += 1

    class_sizes.append(len(cls))

print()
print("=" * 100)
print("C-CLASS SANITY")
print("=" * 100)
print(f"membership failures = {membership_failures}")
print(f"unique C classes    = {len(c_classes)}")
print(f"minimum class size  = {min(class_sizes)}")
print(f"maximum class size  = {max(class_sizes)}")
print(f"mean class size     = {statistics.mean(class_sizes):.6f}")


# ======================================================================================
# GENERATE SHIFTED COLLISION EVENTS
# ======================================================================================
#
# Every other candidate in the same exact C-class gives
#
#     a*b - p*q = 2*M*t
#
# for some integer t.
#
# ======================================================================================

events = []

for anchor_id, anchor in enumerate(anchors):
    pq = anchor["n"]
    C = anchor["C"]

    for a, b in c_classes[C]:
        if (a, b) == (anchor["p"], anchor["q"]):
            continue

        diff = a * b - pq

        if diff % (2 * M) != 0:
            continue

        t = diff // (2 * M)

        if T_MIN <= t <= T_MAX:
            Nt = pq + 2 * M * t

            assert Nt == a * b
            assert Nt % M == C

            events.append({
                "anchor_id": anchor_id,
                "C": C,
                "p": anchor["p"],
                "q": anchor["q"],
                "a": a,
                "b": b,
                "t": t,
                "Nt": Nt,
            })


print()
print("=" * 100)
print("ACTUAL SHIFTED LANDSCAPE")
print("=" * 100)
print(f"collision events = {len(events)}")
print(
    f"anchors with events = "
    f"{len(set(e['anchor_id'] for e in events))}"
)


# ======================================================================================
# PRIME RESIDUE FREQUENCY TABLE
# ======================================================================================

def residue_frequency(R):
    freq = {}

    for p in FACTOR_PRIMES:
        r = p % R
        freq[r] = freq.get(r, 0) + 1

    return freq


# ======================================================================================
# EXACT NUMBER OF UNORDERED PRIME PAIRS
# ======================================================================================
#
# Count ordered solutions:
#
#     a*b = C (mod R)
#
# For each a:
#
#     b = C*a^(-1) (mod R)
#
# If D is the number of diagonal solutions a=b, then
#
#     ordered = 2*unordered_distinct + D
#
# hence
#
#     unordered = (ordered + D)/2.
#
# ======================================================================================

def candidate_pair_count(C, R, freq):
    C_R = C % R

    ordered = 0
    diagonal = 0

    for a in FACTOR_PRIMES:
        inv_a = INV_M[a] % R
        required_b = (C_R * inv_a) % R

        ordered += freq.get(required_b, 0)

        if (a * a) % R == C_R:
            diagonal += 1

    result = (ordered + diagonal) // 2

    return result


# ======================================================================================
# BUILD MODULUS SYSTEMS
# ======================================================================================
#
# First the nested CRT ladder:
#
#   3
#   3*5
#   3*5*7
#   ...
#   M
#
# Then non-nested combinations.
# ======================================================================================

systems = []

# Nested ladder
running = 1

for r in MODULI:
    running *= r
    systems.append(
        (f"LADDER [{','.join(map(str, MODULI[:len(systems)+1]))}]",
         MODULI[:len(systems) + 1])
    )

# Carefully chosen non-nested combinations
manual_systems = [
    [3, 5],
    [3, 7],
    [3, 11],
    [5, 7],
    [5, 11],
    [7, 13],
    [11, 17],
    [13, 19],
    [17, 23],

    [3, 5, 7],
    [3, 5, 11],
    [3, 7, 13],
    [5, 7, 17],
    [5, 11, 19],
    [7, 13, 23],
    [11, 17, 23],

    [3, 5, 7, 11],
    [3, 5, 13, 17],
    [3, 7, 13, 19],
    [5, 7, 17, 23],
    [5, 11, 13, 19],
    [3, 11, 17, 23],

    [3, 5, 7, 11, 13],
    [3, 5, 7, 17, 23],
    [3, 7, 11, 19, 23],
    [5, 7, 13, 17, 19],
]

for subset in manual_systems:
    if subset == MODULI:
        continue

    name = "COMBO [" + ",".join(map(str, subset)) + "]"

    if not any(existing[0] == name for existing in systems):
        systems.append((name, subset))


# ======================================================================================
# NESTED LADDER REPORT
# ======================================================================================

print()
print("=" * 100)
print("NESTED CRT FINGERPRINT LADDER")
print("=" * 100)

print(
    f"{'SYSTEM':30s} "
    f"{'R':>12s} "
    f"{'MEAN PAIRS':>14s} "
    f"{'MEDIAN':>10s} "
    f"{'MIN':>8s} "
    f"{'MAX':>10s} "
    f"{'SINGLETONS':>11s} "
    f"{'REDUCTION':>12s}"
)

print("-" * 100)

baseline_pairs = len(FACTOR_PRIMES) * (len(FACTOR_PRIMES) + 1) // 2

ladder_results = []

for idx, (name, subset) in enumerate(systems[:len(MODULI)]):
    R = math.prod(subset)

    freq = residue_frequency(R)

    counts = [
        candidate_pair_count(anchor["C"], R, freq)
        for anchor in anchors
    ]

    mean_count = statistics.mean(counts)
    median_count = statistics.median(counts)
    min_count = min(counts)
    max_count = max(counts)

    singleton_count = sum(c == 1 for c in counts)

    reduction = baseline_pairs / mean_count

    ladder_results.append({
        "name": name,
        "R": R,
        "count": counts,
        "mean": mean_count,
        "median": median_count,
        "min": min_count,
        "max": max_count,
        "singleton": singleton_count,
        "reduction": reduction,
    })

    print(
        f"{name:30s} "
        f"{R:12,} "
        f"{mean_count:14.3f} "
        f"{median_count:10.1f} "
        f"{min_count:8d} "
        f"{max_count:10d} "
        f"{singleton_count:11d} "
        f"{reduction:12.2f}x"
    )


# ======================================================================================
# MONOTONICITY CHECK
# ======================================================================================
#
# Adding congruence conditions cannot increase a candidate set.
# This checks that our implementation respects that.
# ======================================================================================

print()
print("=" * 100)
print("MONOTONIC CRT REDUCTION CHECK")
print("=" * 100)

monotonic_failures = 0

for anchor_index in range(len(anchors)):
    previous = None

    for result in ladder_results:
        current = result["count"][anchor_index]

        if previous is not None and current > previous:
            monotonic_failures += 1

        previous = current

print(f"monotonicity failures = {monotonic_failures}")


# ======================================================================================
# FULL-MATCH CHECK
# ======================================================================================
#
# At R=M, the candidate count must equal the exact C-class size from above.
# ======================================================================================

full_result = ladder_results[-1]

full_failures = 0

for anchor_index, anchor in enumerate(anchors):
    expected = len(c_classes[anchor["C"]])
    measured = full_result["count"][anchor_index]

    if expected != measured:
        full_failures += 1

print()
print("=" * 100)
print("FULL CRT / EXACT C-CLASS CHECK")
print("=" * 100)
print(f"R=M                    = {M:,}")
print(f"exact C-class failures = {full_failures}")


# ======================================================================================
# EVENT-WEIGHTED CANDIDATE REDUCTION
# ======================================================================================
#
# The same C can generate multiple t-events. We therefore also measure ambiguity
# weighted by the actual collision landscape.
# ======================================================================================

print()
print("=" * 100)
print("EVENT-WEIGHTED CRT FINGERPRINT REDUCTION")
print("=" * 100)

print(
    f"{'SYSTEM':30s} "
    f"{'MEAN PAIRS/EVENT':>18s} "
    f"{'MEDIAN':>10s} "
    f"{'SINGLETON EVENTS':>17s}"
)

print("-" * 100)

for result in ladder_results:
    counts_by_anchor = result["count"]

    event_counts = [
        counts_by_anchor[e["anchor_id"]]
        for e in events
    ]

    singleton_events = sum(c == 1 for c in event_counts)

    print(
        f"{result['name']:30s} "
        f"{statistics.mean(event_counts):18.3f} "
        f"{statistics.median(event_counts):10.1f} "
        f"{singleton_events:17d}"
    )


# ======================================================================================
# NON-NESTED MODULUS COMBINATION RESULTS
# ======================================================================================

print()
print("=" * 100)
print("NON-NESTED MODULUS COMBINATIONS")
print("=" * 100)

print(
    f"{'SYSTEM':30s} "
    f"{'R':12s} "
    f"{'MEAN PAIRS':>14s} "
    f"{'MEDIAN':>10s} "
    f"{'SINGLETONS':>11s} "
    f"{'REDUCTION':>12s}"
)

print("-" * 100)

for name, subset in systems[len(MODULI):]:
    R = math.prod(subset)
    freq = residue_frequency(R)

    counts = [
        candidate_pair_count(anchor["C"], R, freq)
        for anchor in anchors
    ]

    mean_count = statistics.mean(counts)
    median_count = statistics.median(counts)
    singleton_count = sum(c == 1 for c in counts)
    reduction = baseline_pairs / mean_count

    print(
        f"{name:30s} "
        f"{R:12,} "
        f"{mean_count:14.3f} "
        f"{median_count:10.1f} "
        f"{singleton_count:11d} "
        f"{reduction:12.2f}x"
    )


# ======================================================================================
# CANDIDATE REDUCTION FOR INDIVIDUAL ANCHORS
# ======================================================================================

print()
print("=" * 100)
print("EXAMPLE ANCHOR FINGERPRINT COLLAPSE")
print("=" * 100)

examples = sorted(
    range(len(anchors)),
    key=lambda i: len(c_classes[anchors[i]["C"]]),
    reverse=True
)[:10]

for i in examples:
    anchor = anchors[i]

    print()
    print(
        f"anchor=({anchor['p']:,},{anchor['q']:,}) "
        f"C={anchor['C']:,} "
        f"exact_class={len(c_classes[anchor['C']])}"
    )

    for result in ladder_results:
        c = result["count"][i]

        print(
            f"  R={result['R']:12,} "
            f"pairs={c:8d} "
            f"reduction={baseline_pairs / c:12.2f}x"
        )


# ======================================================================================
# FINGERPRINT ACCURACY / UNIQUENESS
# ======================================================================================
#
# Find the first point on the nested CRT ladder at which the observed product
# congruence leaves exactly one prime pair.
# ======================================================================================

print()
print("=" * 100)
print("FIRST UNIQUE-FINGERPRINT LEVEL")
print("=" * 100)

unique_levels = []

for i, anchor in enumerate(anchors):
    level = None

    for j, result in enumerate(ladder_results, 1):
        if result["count"][i] == 1:
            level = j
            break

    unique_levels.append(level)

for level in range(1, len(MODULI) + 1):
    count = sum(x == level for x in unique_levels)

    print(
        f"{level} modulus(es): "
        f"{count:3d} / {len(anchors):3d} anchors become unique here"
    )

never_unique = sum(x is None for x in unique_levels)

print(f"never unique inside ladder: {never_unique}")


# ======================================================================================
# RESIDUE IDENTITY VERIFICATION
# ======================================================================================

identity_failures = 0

for event in events:
    for r in MODULI:
        if (
            event["a"] * event["b"] % r
            != event["Nt"] % r
        ):
            identity_failures += 1

print()
print("=" * 100)
print("MODULAR IDENTITY CHECK")
print("=" * 100)
print(
    "(a*b) mod r = N_t mod r failures = "
    f"{identity_failures}"
)


# ======================================================================================
# FINAL
# ======================================================================================

print()
print("=" * 100)
print("FINAL")
print("=" * 100)

print(f"anchors                  = {len(anchors)}")
print(f"collision events         = {len(events)}")
print(f"exact C classes          = {len(c_classes)}")
print(f"membership failures      = {membership_failures}")
print(f"monotonicity failures    = {monotonic_failures}")
print(f"full CRT class failures  = {full_failures}")
print(f"identity failures        = {identity_failures}")

print()
print("The experiment measures a direct candidate-space quantity:")
print()
print("    candidate pairs satisfying a*b = N_t (mod R)")
print()
print("as more prime modulus fingerprints are combined.")
print()
print("Because every tested R divides M:")
print()
print("    N_t mod R = p*q mod R")
print()
print("so the CRT ladder gives progressively stronger information")
print("about the compatible factor-pair residue classes.")
print()
print("At R=M the result must coincide exactly with the previously")
print("measured C-conditioned prime-pair class.")
print()
print("Experiment complete.")
