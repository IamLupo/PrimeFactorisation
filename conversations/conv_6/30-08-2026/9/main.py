import math
import random
from collections import defaultdict, Counter

# ============================================================================
# EXACT MODULUS-SUBSET / FINGERPRINT ABLATION EXPERIMENT
#
# Goal:
#
#   The previous experiment found:
#
#       K=6 -> 7,046 signatures, max bucket = 2
#       K=7 -> 8,363 signatures, max bucket = 1
#
#   But that used one particular ordered prefix:
#
#       (4,3,5,7,11,13,17)
#
#   This experiment asks:
#
#       Is 7 really the minimum NUMBER of moduli?
#
#   We therefore examine every non-empty subset of:
#
#       {4,3,5,7,11,13,17}
#
#   There are only:
#
#       2^7 - 1 = 127 subsets.
#
#   For every subset we measure:
#
#       - number of distinct prime fingerprints
#       - maximum collision bucket
#       - number of uniquely identified primes
#
#   Then we determine, for every actual anchor:
#
#       - smallest subset size that uniquely identifies p
#       - smallest subset size that uniquely identifies q
#       - smallest subset size that uniquely identifies both
#       - all minimum resolving subsets
#
#   This directly tests whether "more moduli" is the important quantity,
#   or whether particular moduli provide disproportionate resolving power.
# ============================================================================

M = 111_546_435

MIN_FACTOR = 10_000
MAX_FACTOR = 100_000

TRIALS = 300

BASE_MODULI = (
    4,
    3,
    5,
    7,
    11,
    13,
    17,
)

SEED = 1_511_464_998


# ============================================================================
# PRIME SIEVE
# ============================================================================

def sieve(lo: int, hi: int):
    s = bytearray(b"\x01") * (hi + 1)
    s[0:2] = b"\x00\x00"

    limit = math.isqrt(hi)

    for p in range(2, limit + 1):
        if s[p]:
            start = p * p
            count = ((hi - start) // p) + 1
            s[start:hi + 1:p] = b"\x00" * count

    return [x for x in range(lo, hi + 1) if s[x]]


# ============================================================================
# MOD-4 BRANCH
# ============================================================================

def branch4(p: int, q: int) -> str:
    a = p % 4
    b = q % 4

    if a == 1 and b == 1:
        return "11"

    if a == 3 and b == 3:
        return "33"

    if {a, b} == {1, 3}:
        return "13"

    raise ValueError((p, q))


# ============================================================================
# GENERATE ACTUAL ANCHORS
# ============================================================================

def generate_anchors(primes, count, rng):
    anchors = []

    while len(anchors) < count:
        p, q = rng.sample(primes, 2)

        if p > q:
            p, q = q, p

        anchors.append((p, q))

    return anchors


# ============================================================================
# SUBSET REPRESENTATION
# ============================================================================

def subset_indices():
    """
    Yield:

        mask, tuple(indices)

    for every non-empty subset of the 7 base moduli.
    """
    n = len(BASE_MODULI)

    for mask in range(1, 1 << n):
        idx = tuple(
            i for i in range(n)
            if mask & (1 << i)
        )

        yield mask, idx


def subset_moduli(idx):
    return tuple(BASE_MODULI[i] for i in idx)


# ============================================================================
# MAIN
# ============================================================================

def run():

    rng = random.Random(SEED)

    print("=" * 100)
    print("EXACT MODULUS-SUBSET / FINGERPRINT ABLATION EXPERIMENT")
    print("=" * 100)
    print(f"M                    = {M:,}")
    print(f"factor range         = {MIN_FACTOR:,} - {MAX_FACTOR:,}")
    print(f"actual anchors       = {TRIALS}")
    print(f"base moduli          = {BASE_MODULI}")
    print(f"subsets tested       = {(1 << len(BASE_MODULI)) - 1}")
    print(f"seed                 = {SEED:,}")
    print()

    # ------------------------------------------------------------------------
    # PRIME POOL
    # ------------------------------------------------------------------------

    print("=" * 100)
    print("BUILDING PRIME POOL")
    print("=" * 100)

    primes = sieve(MIN_FACTOR, MAX_FACTOR)

    print(f"factor primes = {len(primes):,}")
    print()

    prime_set = set(primes)

    # ------------------------------------------------------------------------
    # PRECOMPUTE EACH PRIME'S 7-DIMENSIONAL RESIDUE VECTOR ONCE
    # ------------------------------------------------------------------------

    print("=" * 100)
    print("PRECOMPUTING PRIME RESIDUE VECTORS")
    print("=" * 100)

    vectors = {
        p: tuple(p % r for r in BASE_MODULI)
        for p in primes
    }

    print("done")
    print()

    # ------------------------------------------------------------------------
    # ANCHORS
    # ------------------------------------------------------------------------

    anchors = generate_anchors(primes, TRIALS, rng)

    print("=" * 100)
    print("GENERATED ACTUAL ANCHORS")
    print("=" * 100)
    print(f"anchors = {len(anchors)}")
    print()

    # ------------------------------------------------------------------------
    # ANALYZE ALL 127 SUBSETS
    # ------------------------------------------------------------------------

    subset_stats = {}
    subset_unique_masks = defaultdict(list)

    print("=" * 100)
    print("EXACT SUBSET ANALYSIS")
    print("=" * 100)

    for subset_no, (mask, idx) in enumerate(subset_indices(), 1):

        mods = subset_moduli(idx)

        buckets = defaultdict(int)

        for p in primes:
            vec = vectors[p]
            key = tuple(vec[i] for i in idx)
            buckets[key] += 1

        signatures = len(buckets)
        unique_primes = sum(
            count == 1
            for count in buckets.values()
        )
        max_bucket = max(buckets.values())

        subset_stats[mask] = {
            "idx": idx,
            "mods": mods,
            "signatures": signatures,
            "unique_primes": unique_primes,
            "max_bucket": max_bucket,
            "buckets": buckets,
        }

        if max_bucket == 1:
            subset_unique_masks[len(idx)].append(mask)

        print(
            f"{subset_no:3d}/127 "
            f"K={len(idx):1d} "
            f"mods={mods!s:<30} "
            f"signatures={signatures:6d} "
            f"unique={unique_primes:6d} "
            f"max_bucket={max_bucket}"
        )

    print()

    # ------------------------------------------------------------------------
    # BEST SUBSETS BY K
    # ------------------------------------------------------------------------

    print("=" * 100)
    print("BEST SUBSETS BY NUMBER OF MODULI")
    print("=" * 100)

    for k in range(1, len(BASE_MODULI) + 1):

        candidates = [
            (mask, data)
            for mask, data in subset_stats.items()
            if len(data["idx"]) == k
        ]

        if not candidates:
            continue

        candidates.sort(
            key=lambda x: (
                x[1]["max_bucket"],
                -x[1]["unique_primes"],
                -x[1]["signatures"],
            )
        )

        best_mask, best = candidates[0]

        print(
            f"K={k} "
            f"best={best['mods']} "
            f"signatures={best['signatures']} "
            f"unique_primes={best['unique_primes']} "
            f"max_bucket={best['max_bucket']}"
        )

    print()

    # ------------------------------------------------------------------------
    # HOW MANY SUBSETS FULLY RESOLVE THE PRIME POOL?
    # ------------------------------------------------------------------------

    print("=" * 100)
    print("FULL-POOL RESOLUTION")
    print("=" * 100)

    for k in range(1, len(BASE_MODULI) + 1):

        full = [
            data
            for data in subset_stats.values()
            if len(data["idx"]) == k
            and data["max_bucket"] == 1
        ]

        total = sum(
            1
            for data in subset_stats.values()
            if len(data["idx"]) == k
        )

        print(
            f"K={k} "
            f"fully_resolving={len(full):3d}/{total:3d}"
        )

    print()

    # ------------------------------------------------------------------------
    # PER-ANCHOR MINIMUM SUBSET SIZE
    # ------------------------------------------------------------------------

    anchor_min_k = []
    anchor_min_masks = []

    for p, q in anchors:

        minimum_k = None
        minimum_masks = []

        for mask, data in subset_stats.items():

            buckets = data["buckets"]
            idx = data["idx"]

            fp_p = tuple(vectors[p][i] for i in idx)
            fp_q = tuple(vectors[q][i] for i in idx)

            cp = buckets[fp_p]
            cq = buckets[fp_q]

            if cp == 1 and cq == 1:

                k = len(idx)

                if minimum_k is None or k < minimum_k:
                    minimum_k = k
                    minimum_masks = [mask]

                elif k == minimum_k:
                    minimum_masks.append(mask)

        anchor_min_k.append(minimum_k)
        anchor_min_masks.append(minimum_masks)

    # ------------------------------------------------------------------------
    # ANCHOR MINIMUM-K DISTRIBUTION
    # ------------------------------------------------------------------------

    print("=" * 100)
    print("MINIMUM MODULUS COUNT PER ACTUAL ANCHOR")
    print("=" * 100)

    counter = Counter(anchor_min_k)

    print(
        " K    anchors"
    )
    print("-" * 40)

    for k in sorted(counter):
        print(f"{k:2d} {counter[k]:12d}")

    print()

    unresolved = sum(
        x is None
        for x in anchor_min_k
    )

    print(f"unresolved anchors = {unresolved}")
    print()

    # ------------------------------------------------------------------------
    # HOW MANY ANCHORS CAN BE RESOLVED AT EACH K?
    # ------------------------------------------------------------------------

    print("=" * 100)
    print("CUMULATIVE ANCHOR RESOLUTION")
    print("=" * 100)

    for k in range(1, len(BASE_MODULI) + 1):

        resolved = sum(
            x is not None and x <= k
            for x in anchor_min_k
        )

        print(
            f"K <= {k}: "
            f"{resolved:3d}/{len(anchors)} "
            f"({resolved / len(anchors):.4%})"
        )

    print()

    # ------------------------------------------------------------------------
    # WHICH MODULUS MATTERS MOST?
    #
    # Count how often each modulus appears in a minimum resolving subset.
    # ------------------------------------------------------------------------

    print("=" * 100)
    print("MODULUS CONTRIBUTION TO MINIMUM RESOLUTION")
    print("=" * 100)

    contribution = Counter()

    for masks in anchor_min_masks:
        for mask in masks:
            for i, r in enumerate(BASE_MODULI):
                if mask & (1 << i):
                    contribution[r] += 1

    if contribution:

        print(
            " modulus    appearances in minimum subsets"
        )
        print("-" * 60)

        for r in BASE_MODULI:
            print(
                f"{r:7d}"
                f"{contribution[r]:35d}"
            )

    print()

    # ------------------------------------------------------------------------
    # LEAVE-ONE-OUT TEST AROUND K=7
    #
    # Test every 6-modulus subset obtained by removing one modulus.
    # ------------------------------------------------------------------------

    print("=" * 100)
    print("LEAVE-ONE-OUT FROM THE 7-MODULUS SET")
    print("=" * 100)

    for removed_index, removed_modulus in enumerate(BASE_MODULI):

        idx = tuple(
            i
            for i in range(len(BASE_MODULI))
            if i != removed_index
        )

        mask = sum(1 << i for i in idx)

        data = subset_stats[mask]

        resolved_anchors = 0

        buckets = data["buckets"]

        for p, q in anchors:

            fp_p = tuple(vectors[p][i] for i in idx)
            fp_q = tuple(vectors[q][i] for i in idx)

            if buckets[fp_p] == 1 and buckets[fp_q] == 1:
                resolved_anchors += 1

        print(
            f"without {removed_modulus:2d}: "
            f"mods={data['mods']} "
            f"prime_max_bucket={data['max_bucket']} "
            f"unique_primes={data['unique_primes']:6d} "
            f"anchors_resolved={resolved_anchors:3d}/{len(anchors)}"
        )

    print()

    # ------------------------------------------------------------------------
    # ANCHOR EXAMPLES
    # ------------------------------------------------------------------------

    print("=" * 100)
    print("ANCHOR EXAMPLES")
    print("=" * 100)

    for i, ((p, q), min_k, masks) in enumerate(
        zip(anchors, anchor_min_k, anchor_min_masks),
        1
    ):

        print()
        print(
            f"anchor={i:3d} "
            f"p={p:,} "
            f"q={q:,} "
            f"n={p*q:,} "
            f"branch={branch4(p,q)}"
        )

        print(f"minimum K = {min_k}")

        if min_k is not None:

            shown = 0

            for mask in masks:

                idx = subset_stats[mask]["idx"]
                mods = subset_stats[mask]["mods"]

                print(
                    f"  resolving subset = {mods}"
                )

                shown += 1

                if shown >= 5:
                    if len(masks) > shown:
                        print(
                            f"  ... {len(masks) - shown} more"
                        )
                    break

    # ------------------------------------------------------------------------
    # MOD-4 BRANCH
    # ------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("MOD-4 BRANCH DISTRIBUTION")
    print("=" * 100)

    branch_counts = Counter(
        branch4(p, q)
        for p, q in anchors
    )

    for branch in ("11", "13", "33"):
        print(
            f"{branch}: {branch_counts[branch]}"
        )

    # ------------------------------------------------------------------------
    # IMPORTANT INTERPRETATION
    # ------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("INTERPRETATION")
    print("=" * 100)

    print(
        "This experiment tests ALL subsets of the first seven moduli,"
    )
    print(
        "rather than assuming that the supplied order is optimal."
    )
    print()
    print(
        "A particularly important result is:"
    )
    print()
    print(
        "    whether any SIX-modulus subset uniquely identifies"
    )
    print(
        "    the entire prime pool."
    )
    print()
    print(
        "If some 6-modulus subsets have max_bucket=1, then the"
    )
    print(
        "previous K=7 threshold was order-dependent."
    )
    print()
    print(
        "If NO six-modulus subset has max_bucket=1, but the"
    )
    print(
        "seven-modulus set does, then seven moduli are genuinely"
    )
    print(
        "required for complete uniqueness over this factor interval"
    )
    print(
        "and this particular prime population."
    )
    print()
    print(
        "For individual anchors, however, the minimum K can be"
    )
    print(
        "smaller than the global threshold."
    )
    print()
    print(
        "This is the distinction between:"
    )
    print(
        "    global uniqueness of every allowed prime"
    )
    print(
        "and"
    )
    print(
        "    local uniqueness of the actual factors p and q."
    )

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()
