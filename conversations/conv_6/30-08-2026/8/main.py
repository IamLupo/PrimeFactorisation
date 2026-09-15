import math
import random
from collections import defaultdict

# ============================================================================
# FAST FACTOR-FINGERPRINT UNIQUENESS EXPERIMENT
#
# Question:
#
#   Given n = p*q and a known factor interval [MIN_FACTOR, MAX_FACTOR],
#   how many independent moduli are needed before:
#
#       p -> unique prime in the interval
#       q -> unique prime in the interval
#
# If both are unique, the factorization is uniquely identified within
# the tested prime population.
# ============================================================================

M = 111_546_435

MIN_FACTOR = 10_000
MAX_FACTOR = 100_000

TRIALS = 300

MODULI = [
    4,
    3,
    5,
    7,
    11,
    13,
    17,
    19,
    23,
    29,
    31,
    37,
    41,
]

SEED = 1_511_464_998


# ============================================================================
# PRIME SIEVE
# ============================================================================

def sieve(lo, hi):
    s = bytearray(b"\x01") * (hi + 1)
    s[0:2] = b"\x00\x00"

    for p in range(2, math.isqrt(hi) + 1):
        if s[p]:
            start = p * p
            s[start:hi + 1:p] = b"\x00" * (
                ((hi - start) // p) + 1
            )

    return [x for x in range(lo, hi + 1) if s[x]]


# ============================================================================
# MOD-4 BRANCH
# ============================================================================

def branch4(p, q):
    r = (p % 4, q % 4)

    if r == (1, 1):
        return "11"

    if r == (3, 3):
        return "33"

    if r == (1, 3):
        return "13"

    if r == (3, 1):
        return "31"

    raise ValueError(r)


def unordered_branch4(p, q):
    r = tuple(sorted((p % 4, q % 4)))

    if r == (1, 1):
        return "11"

    if r == (1, 3):
        return "13"

    if r == (3, 3):
        return "33"

    raise ValueError(r)


# ============================================================================
# FINGERPRINT
# ============================================================================

def fp(p, moduli):
    return tuple(p % r for r in moduli)


# ============================================================================
# BUILD ACTUAL ANCHORS
# ============================================================================

def generate_anchors(primes, count, rng):
    result = []

    while len(result) < count:
        p, q = rng.sample(primes, 2)

        if p > q:
            p, q = q, p

        result.append((p, q))

    return result


# ============================================================================
# MAIN
# ============================================================================

def run():

    rng = random.Random(SEED)

    print("=" * 100)
    print("FAST FACTOR-FINGERPRINT UNIQUENESS / MINIMUM MODULUS EXPERIMENT")
    print("=" * 100)
    print(f"M                    = {M:,}")
    print(f"factor range         = {MIN_FACTOR:,} - {MAX_FACTOR:,}")
    print(f"actual anchors       = {TRIALS}")
    print(f"moduli               = {MODULI}")
    print(f"seed                 = {SEED:,}")
    print()

    # ------------------------------------------------------------------------
    # PRIME POOL
    # ------------------------------------------------------------------------

    primes = sieve(MIN_FACTOR, MAX_FACTOR)

    print("=" * 100)
    print("PRIME POOL")
    print("=" * 100)
    print(f"primes = {len(primes):,}")
    print()

    # ------------------------------------------------------------------------
    # PRECOMPUTE EVERY PRIME'S FULL FINGERPRINT ONCE
    # ------------------------------------------------------------------------

    print("=" * 100)
    print("PRECOMPUTING PRIME FINGERPRINTS")
    print("=" * 100)

    full_fp = {}

    for p in primes:
        full_fp[p] = tuple(p % r for r in MODULI)

    print("done")
    print()

    # ------------------------------------------------------------------------
    # BUILD PREFIX MAPS
    #
    # Instead of scanning primes for every anchor, construct:
    #
    # fingerprint -> number of primes
    #
    # for each prefix.
    # ------------------------------------------------------------------------

    print("=" * 100)
    print("FINGERPRINT UNIQUENESS BY MODULUS PREFIX")
    print("=" * 100)

    prefix_counts = []

    for k in range(1, len(MODULI) + 1):

        counter = defaultdict(int)

        for p in primes:
            key = full_fp[p][:k]
            counter[key] += 1

        unique = sum(1 for v in counter.values() if v == 1)
        signatures = len(counter)
        max_bucket = max(counter.values())

        prefix_counts.append(counter)

        crt_product = 1
        for r in MODULI[:k]:
            crt_product *= r

        print(
            f"K={k:2d} "
            f"added={MODULI[k-1]:2d} "
            f"CRT={crt_product:,} "
            f"signatures={signatures:6d} "
            f"unique-primes={unique:6d} "
            f"max-bucket={max_bucket}"
        )

    print()

    # ------------------------------------------------------------------------
    # ACTUAL ANCHORS
    # ------------------------------------------------------------------------

    anchors = generate_anchors(primes, TRIALS, rng)

    print("=" * 100)
    print("ACTUAL ANCHORS")
    print("=" * 100)
    print(f"anchors = {len(anchors)}")
    print()

    # ------------------------------------------------------------------------
    # FACTOR RESOLUTION
    # ------------------------------------------------------------------------

    first_unique_factor = []
    first_unique_pair = []
    first_unique_branch = []

    for p, q in anchors:

        unique_factor_k = None
        unique_pair_k = None
        unique_branch_k = None

        for k in range(1, len(MODULI) + 1):

            counter = prefix_counts[k - 1]

            fp_p = full_fp[p][:k]
            fp_q = full_fp[q][:k]

            count_p = counter[fp_p]
            count_q = counter[fp_q]

            # Individual factors are uniquely determined.
            if (
                unique_factor_k is None
                and count_p == 1
                and count_q == 1
            ):
                unique_factor_k = k

            # Same condition, explicitly check unordered pair.
            if unique_pair_k is None:

                if fp_p == fp_q:
                    continue

                if count_p == 1 and count_q == 1:
                    unique_pair_k = k

            # ----------------------------------------------------------------
            # BRANCH TEST
            #
            # We don't enumerate candidates.
            #
            # Instead ask whether the two MOD-4 branches contain any primes
            # with the required factor fingerprints.
            # ----------------------------------------------------------------

            if unique_branch_k is None:

                n = p * q

                if n % 4 == 1:

                    # Branch 11
                    branch11 = False

                    # Branch 33
                    branch33 = False

                    for x in primes:

                        fx = full_fp[x][:k]

                        if fx != fp_p:
                            continue

                        if x % 4 != 1:
                            continue

                        # Required q fingerprint is already fp_q.
                        # Since q itself is known to exist, branch testing
                        # is done using the unique factor counts below.

                    #
                    # Once both factors are individually unique, the
                    # branch is necessarily unique.
                    #
                    if count_p == 1 and count_q == 1:
                        unique_branch_k = k

                else:

                    if count_p == 1 and count_q == 1:
                        unique_branch_k = k

        first_unique_factor.append(unique_factor_k)
        first_unique_pair.append(unique_pair_k)
        first_unique_branch.append(unique_branch_k)

    # ------------------------------------------------------------------------
    # THRESHOLD DISTRIBUTION
    # ------------------------------------------------------------------------

    print("=" * 100)
    print("MINIMUM K DISTRIBUTION")
    print("=" * 100)

    print(
        " K   UNIQUE FACTORS   UNIQUE PAIRS   UNIQUE BRANCH"
    )
    print("-" * 100)

    for k in range(1, len(MODULI) + 1):

        a = sum(x == k for x in first_unique_factor)
        b = sum(x == k for x in first_unique_pair)
        c = sum(x == k for x in first_unique_branch)

        if a or b or c:
            print(
                f"{k:2d}"
                f"{a:16d}"
                f"{b:16d}"
                f"{c:17d}"
            )

    never_factor = sum(x is None for x in first_unique_factor)
    never_pair = sum(x is None for x in first_unique_pair)
    never_branch = sum(x is None for x in first_unique_branch)

    print()
    print(f"never uniquely identified factors = {never_factor}")
    print(f"never uniquely identified pair    = {never_pair}")
    print(f"never resolved branch             = {never_branch}")
    print()

    # ------------------------------------------------------------------------
    # FULL SUMMARY
    # ------------------------------------------------------------------------

    print("=" * 100)
    print("RESOLUTION SUMMARY")
    print("=" * 100)

    factor_resolved = [
        x for x in first_unique_factor
        if x is not None
    ]

    pair_resolved = [
        x for x in first_unique_pair
        if x is not None
    ]

    if factor_resolved:
        print(
            f"mean K for unique factors = "
            f"{sum(factor_resolved) / len(factor_resolved):.4f}"
        )

    if pair_resolved:
        print(
            f"mean K for unique pair    = "
            f"{sum(pair_resolved) / len(pair_resolved):.4f}"
        )

    print()

    # ------------------------------------------------------------------------
    # EXAMPLES
    # ------------------------------------------------------------------------

    print("=" * 100)
    print("EXAMPLE ANCHORS")
    print("=" * 100)

    for i, (p, q) in enumerate(anchors[:20]):

        print()
        print(
            f"anchor={i+1:3d} "
            f"p={p:,} "
            f"q={q:,} "
            f"n={p*q:,} "
            f"branch={unordered_branch4(p,q)}"
        )

        print(
            f"  unique factors K = "
            f"{first_unique_factor[i]}"
        )

        print(
            f"  unique pair   K = "
            f"{first_unique_pair[i]}"
        )

        print(
            f"  unique branch K = "
            f"{first_unique_branch[i]}"
        )

    # ------------------------------------------------------------------------
    # DIRECT MOD-4 BRANCH STATISTIC
    # ------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("MOD-4 BRANCH DISTRIBUTION")
    print("=" * 100)

    branches = defaultdict(int)

    for p, q in anchors:
        branches[unordered_branch4(p, q)] += 1

    for b in ("11", "13", "33"):
        print(f"{b}: {branches[b]}")

    # ------------------------------------------------------------------------
    # RELATION TO M
    # ------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("RELATION TO M")
    print("=" * 100)

    print(f"M = {M:,}")

    crt = 1

    for r in MODULI:
        crt *= r

        if crt >= M:
            print(
                f"CRT product reaches M at r={r}: {crt:,}"
            )
            break

    print()
    print(
        "The experiment deliberately uses moduli that are not restricted"
    )
    print(
        "to the prime divisors of M."
    )

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()
