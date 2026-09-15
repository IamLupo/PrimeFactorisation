import math
import random
from collections import defaultdict, Counter

# =============================================================================
# MINIMAL MULTI-MODULUS RESOLVING-SET EXPERIMENT
#
# Purpose
# -------
#
# We now know that:
#
#   - some 6-modulus sets uniquely identify every prime
#   - some 5-modulus sets do so as well
#   - 13 and 17 appear particularly important in the current modulus pool
#
# This experiment asks:
#
#   1. What is the smallest resolving modulus set globally?
#   2. Which moduli occur in minimum resolving sets?
#   3. For the actual 300 factor pairs, what is the minimum K?
#   4. Does a set that uniquely identifies p and q also uniquely resolve
#      the MOD-4 branch?
#
# IMPORTANT:
#
# We work with individual prime fingerprints, not prime-pair enumeration.
#
# Therefore the expensive step is avoided completely.
# =============================================================================


# =============================================================================
# CONFIG
# =============================================================================

M = 111_546_435

MIN_FACTOR = 10_000
MAX_FACTOR = 100_000

TRIALS = 300

SEED = 1_511_464_998

# Include the original mod-4 modulus plus the divisors used previously,
# followed by additional moduli which are NOT divisors of M.
#
# We deliberately stop at 97. The resulting pool is large enough to study
# minimal resolving sets while keeping the fingerprint representation small.
#
MODULI = (
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
    43,
    47,
    53,
    59,
    61,
    67,
    71,
    73,
    79,
    83,
    89,
    97,
)

# To avoid an enormous 2^25 exhaustive search, we first perform:
#
#   - exact greedy forward selection
#   - exact leave-one-out tests
#   - randomized multi-start greedy searches
#
# This is enough to find very small resolving sets in practice and tells us
# which moduli repeatedly appear.
#
RANDOM_GREEDY_RUNS = 1000


# =============================================================================
# SIEVE
# =============================================================================

def sieve(lo: int, hi: int):
    sieve_bytes = bytearray(b"\x01") * (hi + 1)
    sieve_bytes[0:2] = b"\x00\x00"

    limit = math.isqrt(hi)

    for p in range(2, limit + 1):
        if sieve_bytes[p]:
            start = p * p
            count = ((hi - start) // p) + 1
            sieve_bytes[start:hi + 1:p] = b"\x00" * count

    return [
        x
        for x in range(lo, hi + 1)
        if sieve_bytes[x]
    ]


# =============================================================================
# MOD-4 BRANCH
# =============================================================================

def branch4(p: int, q: int) -> str:
    rp = p % 4
    rq = q % 4

    if rp == 1 and rq == 1:
        return "11"

    if rp == 3 and rq == 3:
        return "33"

    if {rp, rq} == {1, 3}:
        return "13"

    raise ValueError("Invalid odd prime pair")


# =============================================================================
# GENERATE ACTUAL ANCHORS
# =============================================================================

def generate_anchors(primes, count, rng):

    anchors = []

    while len(anchors) < count:

        p, q = rng.sample(primes, 2)

        if p > q:
            p, q = q, p

        anchors.append((p, q))

    return anchors


# =============================================================================
# FINGERPRINT UTILITIES
# =============================================================================

def signatures_for_subset(
    prime_vectors,
    indices,
):
    """
    Return the partition of primes induced by the selected moduli.

    A dictionary maps:

        fingerprint -> number of primes

    """
    buckets = defaultdict(int)

    for vec in prime_vectors:
        key = tuple(vec[i] for i in indices)
        buckets[key] += 1

    return buckets


def subset_stats(prime_vectors, indices):

    buckets = signatures_for_subset(
        prime_vectors,
        indices,
    )

    signatures = len(buckets)

    unique_primes = sum(
        1
        for count in buckets.values()
        if count == 1
    )

    max_bucket = max(buckets.values())

    collisions = sum(
        max(0, count - 1)
        for count in buckets.values()
    )

    return (
        signatures,
        unique_primes,
        max_bucket,
        collisions,
        buckets,
    )


def is_globally_resolving(
    prime_vectors,
    indices,
):
    """
    True iff every prime has a unique fingerprint.
    """
    buckets = signatures_for_subset(
        prime_vectors,
        indices,
    )

    return all(
        count == 1
        for count in buckets.values()
    )


# =============================================================================
# GREEDY RESOLUTION
# =============================================================================

def greedy_resolving_set(
    prime_vectors,
    rng,
    randomized=False,
):
    """
    Greedily select moduli.

    At every step choose the modulus that produces the greatest reduction
    in unresolved prime collisions.

    randomized=True randomly perturbs among near-best choices, allowing
    multi-start searches for different minimum-ish sets.
    """

    selected = []
    remaining = list(range(len(MODULI)))

    # Current equivalence classes.
    buckets = {
        tuple(): list(range(len(prime_vectors)))
    }

    while True:

        if all(len(v) == 1 for v in buckets.values()):
            break

        best_gain = None
        best_candidates = []

        for idx in remaining:

            new_buckets = defaultdict(int)

            # Count resulting fingerprints inside the current partition.
            for prime_index in range(len(prime_vectors)):

                vec = prime_vectors[prime_index]

                old_key = tuple(
                    vec[i]
                    for i in selected
                )

                new_key = (
                    old_key,
                    vec[idx],
                )

                new_buckets[new_key] += 1

            # Number of unresolved items after adding this modulus.
            unresolved = sum(
                count - 1
                for count in new_buckets.values()
                if count > 1
            )

            gain = len(prime_vectors) - unresolved

            if best_gain is None or gain > best_gain:
                best_gain = gain
                best_candidates = [idx]

            elif gain == best_gain:
                best_candidates.append(idx)

        if randomized and len(best_candidates) > 1:
            chosen = rng.choice(best_candidates)
        else:
            chosen = min(
                best_candidates,
                key=lambda i: MODULI[i],
            )

        selected.append(chosen)
        remaining.remove(chosen)

        # Rebuild actual buckets using the chosen modulus.
        new_buckets = defaultdict(list)

        for prime_index, vec in enumerate(prime_vectors):

            key = tuple(
                vec[i]
                for i in selected
            )

            new_buckets[key].append(prime_index)

        buckets = new_buckets

    return tuple(selected)


# =============================================================================
# ACTUAL ANCHOR RESOLUTION
# =============================================================================

def anchor_resolved(
    p,
    q,
    indices,
    prime_vectors_by_prime,
    buckets,
):
    """
    Both p and q must have unique prime fingerprints.
    """

    vp = prime_vectors_by_prime[p]
    vq = prime_vectors_by_prime[q]

    kp = tuple(vp[i] for i in indices)
    kq = tuple(vq[i] for i in indices)

    return (
        buckets.get(kp, 0) == 1
        and buckets.get(kq, 0) == 1
    )


# =============================================================================
# BRANCH RESOLUTION
#
# This is different from factor uniqueness.
#
# Given n mod 4:
#
#   n == 1 mod 4
#
# possible factor branches are:
#
#   (1,1)
#   (3,3)
#
# Given n == 3 mod 4:
#
#   (1,3)
#   (3,1)
#
# We ask whether the residue observations distinguish all primes in the
# competing branches over the chosen factor interval.
# =============================================================================

def branch_candidate_count(
    p,
    q,
    indices,
    primes,
    prime_vectors_by_prime,
):
    """
    Count how many ordered residue-compatible prime candidates exist
    for p and q under the selected modulus set.

    We only count prime-coordinate candidates inside the factor interval.

    This is a residue test, NOT a factorization algorithm.
    """

    n = p * q

    target_p_residues = tuple(
        n % MODULI[i]
        for i in indices
    )

    # For every possible first factor a, q-residue is constrained by:
    #
    #   a*b == n (mod r)
    #
    # We count compatible b values.
    #
    # The search is only over primes in the factor range and is therefore
    # manageable for 8,363 primes * 300 anchors.
    #
    compatible = set()

    for a in primes:

        ok = True

        for i in indices:

            r = MODULI[i]

            ra = a % r

            found = False

            for b in primes:

                if (ra * (b % r)) % r == n % r:
                    found = True
                    break

            if not found:
                ok = False
                break

        if ok:
            compatible.add(a)

    return len(compatible)


# =============================================================================
# FAST BRANCH TEST
#
# Instead of full pair search, use factor fingerprints:
#
# The unordered factor signature is compared against the actual factor
# signature. If another pair of primes has the same residue vector, it is
# unresolved.
#
# We construct pair buckets from the prime fingerprints only for the actual
# anchor queries.
# =============================================================================

def build_pair_fingerprint_index(
    primes,
    prime_vectors,
    indices,
):
    """
    Build an index of unordered prime-pair fingerprints for the selected
    moduli.

    This would become expensive for all 35 million prime pairs, so we DO NOT
    construct this globally.

    Function intentionally unused except as a conceptual placeholder.
    """
    raise RuntimeError(
        "Global pair index intentionally disabled."
    )


# =============================================================================
# LEAVE-ONE-OUT IMPORTANCE
# =============================================================================

def leave_one_out(
    primes,
    prime_vectors,
):
    """
    Starting from all moduli, remove one modulus at a time.

    This identifies moduli that are indispensable for the global prime
    fingerprint uniqueness.
    """

    result = []

    full = tuple(range(len(MODULI)))

    full_resolving = is_globally_resolving(
        prime_vectors,
        full,
    )

    for removed in range(len(MODULI)):

        remaining = tuple(
            i
            for i in full
            if i != removed
        )

        resolving = is_globally_resolving(
            prime_vectors,
            remaining,
        )

        buckets = signatures_for_subset(
            prime_vectors,
            remaining,
        )

        result.append(
            (
                MODULI[removed],
                resolving,
                max(buckets.values()),
                sum(
                    1
                    for x in buckets.values()
                    if x == 1
                ),
            )
        )

    return full_resolving, result


# =============================================================================
# MINIMUM-SET SEARCH
# =============================================================================

def exact_small_set_search(
    primes,
    prime_vectors,
):
    """
    Exhaustively search small K only.

    The first run established that the interesting threshold is around 4-6.
    Searching K <= 5 over 25 moduli gives:

        C(25,1)
        C(25,2)
        C(25,3)
        C(25,4)
        C(25,5)

    which is 68,405 subsets.

    Each subset is evaluated by constructing fingerprints of only 8,363
    primes. We stop immediately once a resolving set is found at a given K.

    This is deliberately implemented without pair enumeration.
    """

    from itertools import combinations

    total_primes = len(primes)

    for k in range(1, 6):

        checked = 0

        print()
        print(
            f"Searching exact K={k} subsets..."
        )

        for combo in combinations(
            range(len(MODULI)),
            k,
        ):

            checked += 1

            buckets = defaultdict(int)

            for vec in prime_vectors:

                key = tuple(
                    vec[i]
                    for i in combo
                )

                buckets[key] += 1

            if len(buckets) == total_primes:

                print(
                    f"FOUND resolving set: "
                    f"{tuple(MODULI[i] for i in combo)}"
                )

                print(
                    f"checked = {checked:,}"
                )

                return combo

        print(
            f"No resolving set at K={k}."
        )

    return None


# =============================================================================
# MAIN
# =============================================================================

def run():

    rng = random.Random(SEED)

    print("=" * 100)
    print("MINIMAL MULTI-MODULUS RESOLVING-SET EXPERIMENT")
    print("=" * 100)
    print(f"M                    = {M:,}")
    print(f"factor range         = {MIN_FACTOR:,} - {MAX_FACTOR:,}")
    print(f"actual anchors       = {TRIALS}")
    print(f"modulus count        = {len(MODULI)}")
    print(f"moduli               = {MODULI}")
    print(f"random greedy runs   = {RANDOM_GREEDY_RUNS}")
    print(f"seed                 = {SEED:,}")
    print()

    # =========================================================================
    # PRIME POOL
    # =========================================================================

    print("=" * 100)
    print("BUILDING PRIME POOL")
    print("=" * 100)

    primes = sieve(
        MIN_FACTOR,
        MAX_FACTOR,
    )

    print(
        f"factor primes = {len(primes):,}"
    )
    print()

    # =========================================================================
    # PRIME VECTORS
    # =========================================================================

    print("=" * 100)
    print("PRECOMPUTING PRIME RESIDUE VECTORS")
    print("=" * 100)

    prime_vectors_by_prime = {}

    prime_vectors = []

    for p in primes:

        vec = tuple(
            p % r
            for r in MODULI
        )

        prime_vectors_by_prime[p] = vec
        prime_vectors.append(vec)

    print("done")
    print()

    # =========================================================================
    # ACTUAL ANCHORS
    # =========================================================================

    anchors = generate_anchors(
        primes,
        TRIALS,
        rng,
    )

    print("=" * 100)
    print("GENERATING ACTUAL ANCHORS")
    print("=" * 100)

    print(
        f"anchors = {len(anchors)}"
    )
    print()

    # =========================================================================
    # SINGLE-MODULUS BASELINE
    # =========================================================================

    print("=" * 100)
    print("SINGLE-MODULUS BASELINES")
    print("=" * 100)

    for i, r in enumerate(MODULI):

        buckets = defaultdict(int)

        for vec in prime_vectors:
            buckets[vec[i]] += 1

        print(
            f"r={r:3d} "
            f"signatures={len(buckets):4d} "
            f"max_bucket={max(buckets.values()):5d}"
        )

    print()

    # =========================================================================
    # EXACT SMALL-K SEARCH
    #
    # Search until K=5 finds a globally resolving set.
    #
    # If the machine is too slow, K=5 can be skipped because randomized
    # greedy below will still find candidate sets.
    # =========================================================================

    minimum_exact = exact_small_set_search(
        primes,
        prime_vectors,
    )

    print()

    # =========================================================================
    # RANDOMIZED GREEDY
    # =========================================================================

    print("=" * 100)
    print("RANDOMIZED GREEDY RESOLVING-SET SEARCH")
    print("=" * 100)

    greedy_sets = Counter()

    best_k = None

    for run_number in range(1, RANDOM_GREEDY_RUNS + 1):

        selected = greedy_resolving_set(
            prime_vectors,
            rng,
            randomized=True,
        )

        mods = tuple(
            MODULI[i]
            for i in selected
        )

        greedy_sets[mods] += 1

        k = len(selected)

        if best_k is None or k < best_k:
            best_k = k
            print(
                f"run={run_number:4d} "
                f"NEW BEST K={k} "
                f"set={mods}"
            )

        if run_number % 100 == 0:
            print(
                f"run {run_number:4d}/{RANDOM_GREEDY_RUNS}"
            )

    print()

    print("=" * 100)
    print("MOST FREQUENT GREEDY RESOLVING SETS")
    print("=" * 100)

    for mods, count in greedy_sets.most_common(20):

        print(
            f"count={count:5d} "
            f"K={len(mods):2d} "
            f"mods={mods}"
        )

    print()

    # =========================================================================
    # MODULUS FREQUENCY IN GREEDY SETS
    # =========================================================================

    print("=" * 100)
    print("MODULUS FREQUENCY ACROSS GREEDY SETS")
    print("=" * 100)

    modulus_frequency = Counter()

    for mods, count in greedy_sets.items():

        for r in mods:
            modulus_frequency[r] += count

    for r in sorted(
        MODULI,
        key=lambda x: (
            -modulus_frequency[x],
            x,
        )
    ):

        print(
            f"r={r:3d} "
            f"frequency={modulus_frequency[r]:6d}"
        )

    print()

    # =========================================================================
    # LEAVE-ONE-OUT
    # =========================================================================

    print("=" * 100)
    print("FULL-SET LEAVE-ONE-OUT")
    print("=" * 100)

    full_resolving, loo = leave_one_out(
        primes,
        prime_vectors,
    )

    print(
        f"full modulus set globally resolving = "
        f"{full_resolving}"
    )

    print()

    for r, resolving, max_bucket, unique_count in loo:

        print(
            f"without r={r:3d}: "
            f"resolving={str(resolving):5s} "
            f"max_bucket={max_bucket:3d} "
            f"unique_primes={unique_count:5d}"
        )

    print()

    # =========================================================================
    # ACTUAL ANCHOR MINIMUM K
    #
    # We use the exact exhaustive subsets up through K=5 that were found
    # during the global search. For each anchor we determine whether its p
    # and q are individually unique.
    #
    # We do NOT claim this resolves the product-pair search beyond the prime
    # coordinate uniqueness criterion.
    # =========================================================================

    print("=" * 100)
    print("ACTUAL ANCHOR RESOLUTION")
    print("=" * 100)

    # Collect candidate resolving subsets.
    #
    # Prefer all minimum-K sets discovered by greedy search.
    # Additionally include the exact global minimum set if found.

    candidate_sets = []

    if minimum_exact is not None:
        candidate_sets.append(
            minimum_exact
        )

    for mods in greedy_sets:
        idx = tuple(
            MODULI.index(r)
            for r in mods
        )

        if idx not in candidate_sets:
            candidate_sets.append(idx)

    anchor_min_k = []
    anchor_resolving_sets = []

    for anchor_number, (p, q) in enumerate(
        anchors,
        1,
    ):

        best = None
        best_sets = []

        for idx in candidate_sets:

            buckets = signatures_for_subset(
                prime_vectors,
                idx,
            )

            if anchor_resolved(
                p,
                q,
                idx,
                prime_vectors_by_prime,
                buckets,
            ):

                k = len(idx)

                if best is None or k < best:
                    best = k
                    best_sets = [idx]

                elif k == best:
                    best_sets.append(idx)

        anchor_min_k.append(best)
        anchor_resolving_sets.append(best_sets)

    counter = Counter(anchor_min_k)

    for k in sorted(counter):

        print(
            f"K={str(k):>2s} "
            f"anchors={counter[k]:3d}"
        )

    print(
        f"unresolved anchors = "
        f"{sum(x is None for x in anchor_min_k)}"
    )

    print()

    # =========================================================================
    # MOD-4 BRANCH DISTRIBUTION
    # =========================================================================

    print("=" * 100)
    print("MOD-4 BRANCH DISTRIBUTION")
    print("=" * 100)

    branches = Counter(
        branch4(p, q)
        for p, q in anchors
    )

    for branch in (
        "11",
        "13",
        "33",
    ):

        print(
            f"{branch}: {branches[branch]}"
        )

    print()

    # =========================================================================
    # EXAMPLES
    # =========================================================================

    print("=" * 100)
    print("ANCHOR EXAMPLES")
    print("=" * 100)

    for i in range(min(20, len(anchors))):

        p, q = anchors[i]

        print()
        print(
            f"anchor={i+1:3d} "
            f"p={p:,} "
            f"q={q:,} "
            f"n={p*q:,} "
            f"branch={branch4(p,q)}"
        )

        print(
            f"minimum candidate K="
            f"{anchor_min_k[i]}"
        )

        shown = 0

        for idx in anchor_resolving_sets[i]:

            mods = tuple(
                MODULI[j]
                for j in idx
            )

            print(
                f"  resolving set={mods}"
            )

            shown += 1

            if shown >= 5:
                break

    # =========================================================================
    # INTERPRETATION
    # =========================================================================

    print()
    print("=" * 100)
    print("INTERPRETATION")
    print("=" * 100)

    print(
        "The previous experiment established that six moduli are not "
        "intrinsically required:"
    )

    print()
    print(
        "    several different 6-modulus subsets already uniquely "
        "identify every prime in the factor interval."
    )

    print()

    print(
        "The new experiment searches for the smallest resolving sets "
        "using a much larger modulus pool."
    )

    print()

    print(
        "There are three distinct notions of uniqueness:"
    )

    print()

    print(
        "1. GLOBAL PRIME UNIQUENESS"
    )
    print(
        "   Every prime in [10,000,100,000] has a unique residue fingerprint."
    )

    print()

    print(
        "2. LOCAL FACTOR UNIQUENESS"
    )
    print(
        "   The actual p and q of an anchor each have unique fingerprints."
    )

    print()

    print(
        "3. FACTOR-PAIR UNIQUENESS"
    )
    print(
        "   No other unordered prime pair reproduces all observations."
    )

    print()

    print(
        "The current experiment directly measures (1) and (2)."
    )

    print(
        "It deliberately does not confuse those with (3)."
    )

    print()

    print(
        "The most important structural output is the frequency with which"
    )
    print(
        "13, 17, and the larger non-divisors of M appear in minimum sets."
    )

    print()

    print(
        "If very small resolving sets repeatedly contain the same moduli,"
    )
    print(
        "those moduli are providing disproportionately useful information"
    )
    print(
        "for this particular factor interval."
    )

    print()

    print(
        "The mod-4 branch remains a separate question."
    )

    print(
        "For n == 1 (mod 4), the branch candidates are:"
    )

    print(
        "    (1,1) and (3,3)"
    )

    print()

    print(
        "For n == 3 (mod 4), the branch candidates are:"
    )

    print(
        "    (1,3) and (3,1)"
    )

    print()

    print(
        "A modulus set that uniquely identifies individual primes"
    )
    print(
        "necessarily resolves the branch once the candidate factor"
    )
    print(
        "domain is restricted to the tested prime interval."
    )

    print()

    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()
