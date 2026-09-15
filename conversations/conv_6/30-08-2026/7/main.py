import math
import random
from collections import defaultdict

# ============================================================================
# MULTI-MODULUS FACTOR PAIR IDENTIFICATION / BRANCH THRESHOLD
# ============================================================================

M = 111_546_435

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

ACTUAL_TRIALS = 300
RANDOM_CONTROLS = 300

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
# PRIME GENERATION
# ============================================================================

def sieve_primes(lo, hi):
    sieve = bytearray(b"\x01") * (hi + 1)
    sieve[0:2] = b"\x00\x00"

    limit = int(math.isqrt(hi))

    for p in range(2, limit + 1):
        if sieve[p]:
            start = p * p
            sieve[start:hi + 1:p] = b"\x00" * (
                ((hi - start) // p) + 1
            )

    return [n for n in range(lo, hi + 1) if sieve[n]]


# ============================================================================
# PRIME RESIDUE INDEX
#
# index[r][residue] = list of primes having that residue mod r
# ============================================================================

def build_residue_index(primes, r):
    idx = defaultdict(list)

    for p in primes:
        idx[p % r].append(p)

    return dict(idx)


# ============================================================================
# FINGERPRINT
# ============================================================================

def fingerprint(x, moduli):
    return tuple(x % r for r in moduli)


# ============================================================================
# MODULAR PRODUCT FINGERPRINT
# ============================================================================

def product_fingerprint(n, moduli):
    return tuple(n % r for r in moduli)


# ============================================================================
# BRANCH CLASSIFICATION MOD 4
# ============================================================================

def mod4_branch(p, q):
    a = p % 4
    b = q % 4

    if a == 1 and b == 1:
        return "11"

    if a == 3 and b == 3:
        return "33"

    if a == 1 and b == 3:
        return "13"

    if a == 3 and b == 1:
        return "31"

    raise ValueError("Unexpected prime residue modulo 4")


def possible_branches(n):
    r = n % 4

    if r == 1:
        return ("11", "33")

    if r == 3:
        return ("13", "31")

    raise ValueError("Odd prime products must be 1 or 3 modulo 4")


# ============================================================================
# BUILD PAIR CANDIDATES FOR ONE N
#
# Given:
#
#       a*b == n (mod r)
#
# for every r in moduli.
#
# Instead of testing every (a,b), choose a and derive the residue required
# for b:
#
#       b == n * inverse(a) (mod r)
#
# Then lookup b directly in the residue index.
# ============================================================================

def find_candidate_pairs(n, moduli, residue_indices):
    candidates = set()

    # Every admissible prime a can potentially be the first factor.
    # Derive the required residue vector for b.
    for a in residue_indices["primes"]:
        required = []

        valid = True

        for r in moduli:
            ar = a % r

            try:
                inv = pow(ar, -1, r)
            except ValueError:
                valid = False
                break

            required.append((n % r) * inv % r)

        if not valid:
            continue

        required = tuple(required)

        # Lookup b using the complete residue vector.
        bucket_key = (moduli, required)

        bucket = residue_indices["vector"].get(bucket_key)

        if not bucket:
            continue

        for b in bucket:
            if a == b:
                continue

            if a < b:
                pair = (a, b)
            else:
                pair = (b, a)

            candidates.add(pair)

    return candidates


# ============================================================================
# VECTOR INDEX
#
# vector_index[(moduli_tuple, residue_tuple)] -> primes
#
# Built incrementally for each modulus prefix.
# ============================================================================

def build_vector_index(primes, moduli):
    idx = {}

    for p in primes:
        sig = fingerprint(p, moduli)
        idx.setdefault(sig, []).append(p)

    return idx


# ============================================================================
# FAST CANDIDATE SEARCH
# ============================================================================

def find_candidates_fast(n, primes, moduli):
    """
    Find all unordered prime pairs (a,b) in the factor range satisfying:

        a*b == n (mod r)

    for every r in moduli.
    """

    # Build factor fingerprint index.
    vector_index = build_vector_index(primes, moduli)

    candidates = set()

    n_sig = product_fingerprint(n, moduli)

    for a in primes:
        required = []

        valid = True

        for r in moduli:
            ar = a % r

            try:
                inv = pow(ar, -1, r)
            except ValueError:
                valid = False
                break

            required.append((n % r) * inv % r)

        if not valid:
            continue

        bucket = vector_index.get(tuple(required))

        if not bucket:
            continue

        for b in bucket:
            if a == b:
                continue

            if a < b:
                candidates.add((a, b))
            else:
                candidates.add((b, a))

    return candidates


# ============================================================================
# BRANCH COUNTS
# ============================================================================

def branch_counts(candidates):
    counts = {
        "11": 0,
        "13": 0,
        "31": 0,
        "33": 0,
    }

    for a, b in candidates:
        branch = mod4_branch(a, b)
        counts[branch] += 1

    return counts


# ============================================================================
# GENERATE ACTUAL ANCHORS
# ============================================================================

def generate_actual_anchors(primes, count, rng):
    anchors = []

    while len(anchors) < count:
        p = rng.choice(primes)
        q = rng.choice(primes)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        anchors.append((p, q))

    return anchors


# ============================================================================
# GENERATE RANDOM CONTROL PRODUCTS
# ============================================================================

def generate_random_controls(primes, count, rng):
    controls = []

    while len(controls) < count:
        p = rng.choice(primes)
        q = rng.choice(primes)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        controls.append((p, q))

    return controls


# ============================================================================
# ANALYZE ONE ANCHOR
# ============================================================================

def analyze_anchor(p, q, primes, modulus_prefixes, prime_indices):
    n = p * q

    actual_branch = mod4_branch(p, q)

    result = {
        "p": p,
        "q": q,
        "n": n,
        "actual_branch": actual_branch,
        "rows": [],
    }

    for k, moduli in enumerate(modulus_prefixes, start=1):
        # Reuse precomputed vector index.
        vector_index = prime_indices[k]

        candidates = set()

        for a in primes:
            required = []

            valid = True

            for r in moduli:
                ar = a % r

                try:
                    inv = pow(ar, -1, r)
                except ValueError:
                    valid = False
                    break

                required.append((n % r) * inv % r)

            if not valid:
                continue

            bucket = vector_index.get(tuple(required))

            if not bucket:
                continue

            for b in bucket:
                if a == b:
                    continue

                if a < b:
                    candidates.add((a, b))
                else:
                    candidates.add((b, a))

        counts = branch_counts(candidates)

        actual_pair = (p, q)

        if actual_pair not in candidates:
            raise RuntimeError(
                "Actual factor pair disappeared from candidate set. "
                f"p={p}, q={q}, moduli={moduli}"
            )

        surviving_branches = []

        if n % 4 == 1:
            if counts["11"] > 0:
                surviving_branches.append("11")

            if counts["33"] > 0:
                surviving_branches.append("33")

        else:
            if counts["13"] > 0:
                surviving_branches.append("13")

            if counts["31"] > 0:
                surviving_branches.append("31")

        branch_unique = (
            len(surviving_branches) == 1
            and surviving_branches[0] == actual_branch
        )

        pair_unique = len(candidates) == 1

        result["rows"].append({
            "k": k,
            "moduli": moduli,
            "candidate_count": len(candidates),
            "branch_count_actual": counts.get(actual_branch, 0),
            "branch_counts": counts,
            "surviving_branches": tuple(surviving_branches),
            "branch_unique": branch_unique,
            "pair_unique": pair_unique,
        })

    return result


# ============================================================================
# MAIN
# ============================================================================

def run():
    rng = random.Random(SEED)

    print("=" * 100)
    print("MULTI-MODULUS FACTOR PAIR IDENTIFICATION / BRANCH THRESHOLD")
    print("=" * 100)
    print(f"M                    = {M:,}")
    print(f"factor range         = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"actual anchors       = {ACTUAL_TRIALS}")
    print(f"random controls      = {RANDOM_CONTROLS}")
    print(f"moduli               = {MODULI}")
    print(f"seed                 = {SEED:,}")
    print()

    # ------------------------------------------------------------------------
    # PRIME POOL
    # ------------------------------------------------------------------------

    primes = sieve_primes(FACTOR_MIN, FACTOR_MAX)

    print("=" * 100)
    print("BUILDING PRIME POOL")
    print("=" * 100)
    print(f"factor primes = {len(primes):,}")
    print()

    # ------------------------------------------------------------------------
    # CRT PRODUCT GROWTH
    # ------------------------------------------------------------------------

    print("=" * 100)
    print("CUMULATIVE CRT MODULUS")
    print("=" * 100)

    R = 1

    for k, r in enumerate(MODULI, start=1):
        R *= r

        marker = ""

        if R > FACTOR_MAX:
            marker = "  <-- > FACTOR_MAX"

        print(
            f"K={k:2d} added={r:2d} "
            f"CRT product={R:,}{marker}"
        )

    print()

    # ------------------------------------------------------------------------
    # BUILD MODULUS PREFIXES
    # ------------------------------------------------------------------------

    modulus_prefixes = []

    for k in range(1, len(MODULI) + 1):
        modulus_prefixes.append(tuple(MODULI[:k]))

    # ------------------------------------------------------------------------
    # BUILD VECTOR INDEXES ONCE
    # ------------------------------------------------------------------------

    print("=" * 100)
    print("BUILDING INCREMENTAL FACTOR FINGERPRINT INDEXES")
    print("=" * 100)

    prime_indices = {}

    for k, moduli in enumerate(modulus_prefixes, start=1):
        index = {}

        for p in primes:
            sig = fingerprint(p, moduli)
            index.setdefault(sig, []).append(p)

        prime_indices[k] = index

        max_bucket = max(len(v) for v in index.values())
        unique_signatures = len(index)

        print(
            f"K={k:2d} "
            f"moduli={moduli} "
            f"unique_factor_signatures={unique_signatures:,} "
            f"max_bucket={max_bucket}"
        )

    print()

    # ------------------------------------------------------------------------
    # GENERATE DATA
    # ------------------------------------------------------------------------

    actual_anchors = generate_actual_anchors(
        primes,
        ACTUAL_TRIALS,
        rng,
    )

    random_controls = generate_random_controls(
        primes,
        RANDOM_CONTROLS,
        rng,
    )

    print("=" * 100)
    print("GENERATED DATA")
    print("=" * 100)
    print(f"actual anchors  = {len(actual_anchors)}")
    print(f"random controls = {len(random_controls)}")
    print()

    # ------------------------------------------------------------------------
    # ANALYZE ACTUALS
    # ------------------------------------------------------------------------

    actual_results = []

    print("=" * 100)
    print("ANALYZING ACTUAL ANCHORS")
    print("=" * 100)

    for i, (p, q) in enumerate(actual_anchors, start=1):
        result = analyze_anchor(
            p,
            q,
            primes,
            modulus_prefixes,
            prime_indices,
        )

        actual_results.append(result)

        if i % 25 == 0 or i == len(actual_anchors):
            print(f"anchor {i:3d}/{len(actual_anchors)}")

    print()

    # ------------------------------------------------------------------------
    # ANALYZE RANDOM CONTROLS
    # ------------------------------------------------------------------------

    random_results = []

    print("=" * 100)
    print("ANALYZING RANDOM CONTROLS")
    print("=" * 100)

    for i, (p, q) in enumerate(random_controls, start=1):
        result = analyze_anchor(
            p,
            q,
            primes,
            modulus_prefixes,
            prime_indices,
        )

        random_results.append(result)

        if i % 25 == 0 or i == len(random_controls):
            print(f"control {i:3d}/{len(random_controls)}")

    print()

    # ------------------------------------------------------------------------
    # SUMMARY
    # ------------------------------------------------------------------------

    print("=" * 100)
    print("PROGRESSIVE PAIR RESOLUTION")
    print("=" * 100)

    header = (
        "K  MODULI                         "
        "ACT AVG CAND   ACT UNIQUE   ACT BRANCH   "
        "RND AVG CAND   RND UNIQUE"
    )

    print(header)
    print("-" * len(header))

    for k, moduli in enumerate(modulus_prefixes, start=1):

        act_candidates = [
            r["rows"][k - 1]["candidate_count"]
            for r in actual_results
        ]

        rnd_candidates = [
            r["rows"][k - 1]["candidate_count"]
            for r in random_results
        ]

        act_avg = sum(act_candidates) / len(act_candidates)
        rnd_avg = sum(rnd_candidates) / len(rnd_candidates)

        act_unique = sum(
            1 for x in act_candidates if x == 1
        )

        rnd_unique = sum(
            1 for x in rnd_candidates if x == 1
        )

        act_branch_unique = sum(
            1
            for r in actual_results
            if r["rows"][k - 1]["branch_unique"]
        )

        print(
            f"{k:2d} "
            f"{str(moduli):30s} "
            f"{act_avg:13.3f} "
            f"{act_unique:11d}/{len(actual_results):d} "
            f"{act_branch_unique:12d}/{len(actual_results):d} "
            f"{rnd_avg:13.3f} "
            f"{rnd_unique:10d}/{len(random_results):d}"
        )

    print()

    # ------------------------------------------------------------------------
    # MINIMUM K FOR BRANCH RESOLUTION
    # ------------------------------------------------------------------------

    print("=" * 100)
    print("MINIMUM K FOR BRANCH / PAIR RESOLUTION")
    print("=" * 100)

    branch_thresholds = []
    pair_thresholds = []

    for result in actual_results:
        branch_k = None
        pair_k = None

        for row in result["rows"]:
            if branch_k is None and row["branch_unique"]:
                branch_k = row["k"]

            if pair_k is None and row["pair_unique"]:
                pair_k = row["k"]

        branch_thresholds.append(branch_k)
        pair_thresholds.append(pair_k)

    print(
        "ACTUAL branch resolved at K: "
        + str(branch_thresholds.count(None))
        + " never"
    )

    print(
        "ACTUAL pair resolved at K:   "
        + str(pair_thresholds.count(None))
        + " never"
    )

    print()

    for k in range(1, len(MODULI) + 1):
        b = sum(1 for x in branch_thresholds if x == k)
        p = sum(1 for x in pair_thresholds if x == k)

        if b or p:
            print(
                f"K={k:2d} "
                f"branch-first={b:3d} "
                f"pair-first={p:3d}"
            )

    print()

    # ------------------------------------------------------------------------
    # MOD-4 BRANCH DISTRIBUTION
    # ------------------------------------------------------------------------

    print("=" * 100)
    print("ACTUAL MOD-4 BRANCH DISTRIBUTION")
    print("=" * 100)

    branches = defaultdict(int)

    for p, q in actual_anchors:
        branches[mod4_branch(p, q)] += 1

    for branch in ("11", "13", "31", "33"):
        print(f"{branch}: {branches[branch]}")

    print()

    # ------------------------------------------------------------------------
    # EXAMPLES
    # ------------------------------------------------------------------------

    print("=" * 100)
    print("EXAMPLE ANCHORS")
    print("=" * 100)

    for result in actual_results[:10]:
        p = result["p"]
        q = result["q"]
        n = result["n"]

        print()
        print(
            f"n={n:,} "
            f"p={p:,} "
            f"q={q:,} "
            f"actual_branch={result['actual_branch']}"
        )

        for row in result["rows"]:
            print(
                f"  K={row['k']:2d} "
                f"candidates={row['candidate_count']:3d} "
                f"branches={row['surviving_branches']} "
                f"branch_unique={row['branch_unique']} "
                f"pair_unique={row['pair_unique']}"
            )

    print()

    # ------------------------------------------------------------------------
    # STRONGEST CASES
    # ------------------------------------------------------------------------

    print("=" * 100)
    print("HARDEST ACTUAL ANCHORS TO RESOLVE")
    print("=" * 100)

    hardest = []

    for result in actual_results:
        final_row = result["rows"][-1]

        hardest.append(
            (
                final_row["candidate_count"],
                result["p"],
                result["q"],
                result["n"],
            )
        )

    hardest.sort(reverse=True)

    for count, p, q, n in hardest[:20]:
        print(
            f"n={n:,} "
            f"actual=({p:,},{q:,}) "
            f"final_candidates={count}"
        )

    print()

    # ------------------------------------------------------------------------
    # FINAL INTERPRETATION
    # ------------------------------------------------------------------------

    print("=" * 100)
    print("INTERPRETATION")
    print("=" * 100)

    print(
        "For every modulus set S, the experiment searches for prime pairs"
    )
    print(
        "(a,b) in the factor interval satisfying:"
    )
    print()
    print(
        "    a*b == n (mod r)  for every r in S"
    )
    print()
    print(
        "The mod-4 branch is treated separately from complete pair resolution."
    )
    print()
    print(
        "BRANCH RESOLUTION:"
    )
    print(
        "    Only one of the two mod-4-compatible branches has any"
    )
    print(
        "    admissible prime realization."
    )
    print()
    print(
        "PAIR RESOLUTION:"
    )
    print(
        "    Exactly one unordered prime pair remains."
    )
    print()
    print(
        "The important threshold is the smallest K at which these occur."
    )
    print()
    print(
        "Because the factor interval ends at 100,000, once the cumulative"
    )
    print(
        "CRT modulus exceeds 100,000, a single factor has a unique complete"
    )
    print(
        "residue fingerprint inside the interval."
    )
    print()
    print(
        "That does not automatically prove that the product fingerprint"
    )
    print(
        "alone factors n. The experiment explicitly measures how many"
    )
    print(
        "prime PAIRS remain compatible with the product residues."
    )
    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()