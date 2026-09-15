#!/usr/bin/env python3

import math
import random
from collections import defaultdict, Counter

# =============================================================================
# CONFIGURATION
# =============================================================================

M = 111_546_435
FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

TRIALS = 300
ORDER_REPETITIONS = 200

# We already know these are useful independent moduli.
# They do NOT divide M.
EXTRA_MODULI = [
    29, 31, 37, 41, 43, 47, 53, 59,
    61, 67, 71, 73, 79, 83, 89, 97
]

SEED = 1_511_464_998

# Number of examples printed in detail
PRINT_EXAMPLES = 15


# =============================================================================
# BASIC HELPERS
# =============================================================================

def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def build_prime_pool():
    return [
        p for p in range(FACTOR_MIN, FACTOR_MAX + 1)
        if is_prime(p)
    ]


def make_actual_anchors(primes, count, rng):
    anchors = set()

    while len(anchors) < count:
        p = rng.choice(primes)
        q = rng.choice(primes)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        anchors.add((p, q))

    return sorted(anchors)


def branch_from_mod4(p: int, q: int):
    rp = p % 4
    rq = q % 4

    if rp == 1 and rq == 1:
        return (1, 1)

    if rp == 3 and rq == 3:
        return (3, 3)

    if {rp, rq} == {1, 3}:
        return (1, 3)

    raise ValueError("Unexpected odd prime residues modulo 4")


def unordered_signature(pair, moduli):
    """
    Full unordered factor fingerprint:
        ({a mod r, b mod r}) for every r.

    We sort each residue pair so that (x,y) == (y,x).
    """
    a, b = pair

    return tuple(
        tuple(sorted((a % r, b % r)))
        for r in moduli
    )


def branch_signature(pair, moduli):
    """
    Fingerprint of a hypothetical pair for branch discrimination.
    Ordering of factors is ignored.
    """
    return unordered_signature(pair, moduli)


# =============================================================================
# PRECOMPUTE PRIME RESIDUES
# =============================================================================

def build_residue_sets(primes, moduli):
    result = {}

    for r in moduli:
        table = defaultdict(list)

        for p in primes:
            table[p % r].append(p)

        result[r] = table

    return result


# =============================================================================
# CANDIDATE GENERATION
# =============================================================================

def candidates_for_modulus_constraints(
    anchor_product,
    branch,
    used_moduli,
    primes,
    residue_tables,
):
    """
    Generate candidate prime pairs satisfying all product congruences and
    the specified modulo-4 branch.

    This is done incrementally rather than scanning all P^2 pairs.

    Because q is determined modulo each modulus from p, we construct the
    required residue signature and retrieve primes matching it.
    """

    # For every prime p in the pool, determine the required residue of q.
    candidates = []

    for p in primes:

        # Enforce branch first.
        rp4 = p % 4

        if branch == (1, 1):
            if rp4 != 1:
                continue

        elif branch == (3, 3):
            if rp4 != 3:
                continue

        # p cannot be zero modulo these odd prime moduli because p itself is
        # > 10k and all supplied moduli are smaller primes. Still, keep this
        # mathematically safe.
        ok = True
        required_q_residues = []

        for r in used_moduli:
            rp = p % r

            if math.gcd(rp, r) != 1:
                ok = False
                break

            inv = pow(rp, -1, r)
            rq = (anchor_product % r) * inv % r

            required_q_residues.append((r, rq))

        if not ok:
            continue

        # Intersect residue classes for q.
        q_pool = None

        for r, rq in required_q_residues:
            current = residue_tables[r].get(rq, [])

            if q_pool is None:
                q_pool = current
            else:
                # Intersect by converting the smaller side to a set.
                if len(q_pool) < len(current):
                    s = set(q_pool)
                    q_pool = [x for x in current if x in s]
                else:
                    s = set(current)
                    q_pool = [x for x in q_pool if x in s]

            if not q_pool:
                break

        if not q_pool:
            continue

        for q in q_pool:
            if q == p:
                continue

            if p > q:
                pair = (q, p)
            else:
                pair = (p, q)

            # Recheck all modular constraints exactly.
            if all(
                (pair[0] * pair[1] - anchor_product) % r == 0
                for r in used_moduli
            ):
                candidates.append(pair)

    return sorted(set(candidates))


# =============================================================================
# FASTER GLOBAL CANDIDATE INDEX
# =============================================================================

def build_fingerprint_index(primes, moduli):
    """
    For a given modulus prefix, map the entire factor fingerprint to the
    actual unordered prime pair.

    This lets us determine uniqueness without repeatedly solving every
    candidate set from scratch.
    """
    index = {}

    # Number of pairs is ~35 million, so DO NOT construct this globally.
    #
    # This function is intentionally unused. It exists to make clear that
    # we deliberately avoid the memory-heavy approach used by the failed
    # experiment.
    return index


# =============================================================================
# MAIN EXPERIMENT
# =============================================================================

def main():

    rng = random.Random(SEED)

    print("=" * 100)
    print("MINIMUM MULTI-MODULUS FINGERPRINT / BRANCH + PAIR RESOLUTION")
    print("=" * 100)
    print(f"M                    = {M:,}")
    print(f"factor range         = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"actual anchors       = {TRIALS}")
    print(f"random modulus orders = {ORDER_REPETITIONS}")
    print(f"extra moduli         = {EXTRA_MODULI}")
    print(f"seed                 = {SEED:,}")

    # -------------------------------------------------------------------------
    # PRIME POOL
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("BUILDING PRIME POOL")
    print("=" * 100)

    primes = build_prime_pool()

    print(f"factor primes = {len(primes):,}")

    # -------------------------------------------------------------------------
    # ACTUAL ANCHORS
    # -------------------------------------------------------------------------

    anchors = make_actual_anchors(primes, TRIALS, rng)

    print(f"actual anchors = {len(anchors)}")

    # -------------------------------------------------------------------------
    # RESIDUE TABLES
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("BUILDING RESIDUE INDEX")
    print("=" * 100)

    residue_tables = build_residue_sets(primes, EXTRA_MODULI)

    for r in EXTRA_MODULI:
        print(f"r={r:3d} residue classes={len(residue_tables[r]):4d}")

    # -------------------------------------------------------------------------
    # STORAGE
    # -------------------------------------------------------------------------

    branch_first_hits = Counter()
    pair_first_hits = Counter()
    exact_first_hits = Counter()

    branch_first_per_anchor = []
    pair_first_per_anchor = []

    order_survival = Counter()

    # Used to see whether one particular modulus repeatedly causes the
    # decisive reduction.
    decisive_modulus = Counter()

    # -------------------------------------------------------------------------
    # RUN RANDOM ORDERS
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("RUNNING RANDOM MODULUS ORDERS")
    print("=" * 100)

    for order_id in range(1, ORDER_REPETITIONS + 1):

        mod_order = EXTRA_MODULI[:]
        rng.shuffle(mod_order)

        if order_id % 10 == 0:
            print(
                f"order {order_id:4d}/{ORDER_REPETITIONS}"
            )

        # Each anchor is tested independently under this modulus order.
        for anchor_index, (p, q) in enumerate(anchors):

            n = p * q

            actual_branch = branch_from_mod4(p, q)

            branch_resolved_at = None
            pair_resolved_at = None

            previous_candidate_count = None

            # -----------------------------------------------------------------
            # Progressive fingerprint
            # -----------------------------------------------------------------

            for k in range(1, len(mod_order) + 1):

                prefix = mod_order[:k]

                # Test both modulo-4 same-residue branches.
                #
                # Only meaningful for n == 1 mod 4.
                # Since the purpose of this experiment is the original
                # 11-vs-33 idea, we only score those anchors.
                if n % 4 == 1:

                    cand11 = candidates_for_modulus_constraints(
                        n,
                        (1, 1),
                        prefix,
                        primes,
                        residue_tables,
                    )

                    cand33 = candidates_for_modulus_constraints(
                        n,
                        (3, 3),
                        prefix,
                        primes,
                        residue_tables,
                    )

                    branch_alive = []

                    if cand11:
                        branch_alive.append((1, 1))

                    if cand33:
                        branch_alive.append((3, 3))

                    # ---------------------------------------------------------
                    # Branch resolution
                    # ---------------------------------------------------------

                    if branch_resolved_at is None:

                        if len(branch_alive) == 1:

                            branch_resolved_at = k

                            # Determine which modulus was decisive.
                            decisive_modulus[mod_order[k - 1]] += 1

                    # ---------------------------------------------------------
                    # Actual pair uniqueness
                    # ---------------------------------------------------------

                    actual_candidates = (
                        cand11 if actual_branch == (1, 1)
                        else cand33
                    )

                    if pair_resolved_at is None:

                        # Ignore ordering of p,q.
                        if len(actual_candidates) == 1:
                            pair_resolved_at = k

                            # Check that it really is the actual pair.
                            if actual_candidates[0] == (p, q):
                                exact_first_hits[k] += 1

                else:
                    # For n == 3 mod 4 the cross branch is forced by mod 4.
                    #
                    # Pair resolution is still useful.
                    #
                    # Construct candidates with (1,3) implicitly by testing
                    # both orientations through p residues.
                    #
                    # The generic candidate routine isn't restricted to a
                    # cross branch, so handle this directly.

                    actual_candidates = []

                    for a in primes:

                        if a % 4 not in (1, 3):
                            continue

                        # q residues are forced.
                        q_pool = None
                        ok = True

                        for r in prefix:
                            ra = a % r

                            if math.gcd(ra, r) != 1:
                                ok = False
                                break

                            rq = (n % r) * pow(ra, -1, r) % r
                            current = residue_tables[r].get(rq, [])

                            if q_pool is None:
                                q_pool = current
                            else:
                                s = set(current)
                                q_pool = [
                                    x for x in q_pool
                                    if x in s
                                ]

                            if not q_pool:
                                ok = False
                                break

                        if not ok:
                            continue

                        for b in q_pool:
                            if a == b:
                                continue

                            if (a * b) % 4 != n % 4:
                                continue

                            pair = tuple(sorted((a, b)))

                            actual_candidates.append(pair)

                    actual_candidates = sorted(set(actual_candidates))

                    if pair_resolved_at is None:
                        if (
                            len(actual_candidates) == 1
                            and actual_candidates[0] == (p, q)
                        ):
                            pair_resolved_at = k
                            exact_first_hits[k] += 1

                # Track candidate count collapse.
                current_count = len(
                    actual_candidates
                    if 'actual_candidates' in locals()
                    else []
                )

                if (
                    previous_candidate_count is not None
                    and current_count < previous_candidate_count
                ):
                    order_survival[(k - 1, k)] += 1

                previous_candidate_count = current_count

                # Clean local to avoid stale use.
                if 'actual_candidates' in locals():
                    del actual_candidates

            # -----------------------------------------------------------------
            # Save anchor result
            # -----------------------------------------------------------------

            if n % 4 == 1:
                if branch_resolved_at is None:
                    branch_first_per_anchor.append(None)
                else:
                    branch_first_per_anchor.append(branch_resolved_at)

            if pair_resolved_at is None:
                pair_first_per_anchor.append(None)
            else:
                pair_first_per_anchor.append(pair_resolved_at)

    # =============================================================================
    # SUMMARY
    # =============================================================================

    print()
    print("=" * 100)
    print("BRANCH RESOLUTION DISTRIBUTION")
    print("=" * 100)

    for k in range(1, len(EXTRA_MODULI) + 1):
        v = branch_first_hits[k]
        print(
            f"K={k:2d} first branch resolution events = {v}"
        )

    print()
    print("=" * 100)
    print("PAIR RESOLUTION DISTRIBUTION")
    print("=" * 100)

    for k in range(1, len(EXTRA_MODULI) + 1):
        v = pair_first_hits[k]
        print(
            f"K={k:2d} first pair resolution events = {v}"
        )

    # -------------------------------------------------------------------------
    # Aggregate anchor statistics
    # -------------------------------------------------------------------------

    def distribution(values):
        c = Counter(v for v in values if v is not None)
        return c

    branch_dist = distribution(branch_first_per_anchor)
    pair_dist = distribution(pair_first_per_anchor)

    print()
    print("=" * 100)
    print("FIRST RESOLUTION BY ANCHOR")
    print("=" * 100)

    print(
        "branch-resolved anchors = "
        f"{sum(v is not None for v in branch_first_per_anchor):3d}"
        f"/{len(branch_first_per_anchor)}"
    )

    print(
        "pair-resolved anchors   = "
        f"{sum(v is not None for v in pair_first_per_anchor):3d}"
        f"/{len(pair_first_per_anchor)}"
    )

    print()
    print(
        f"{'K':>3} "
        f"{'BRANCH FIRST':>14} "
        f"{'PAIR FIRST':>12} "
        f"{'CUM BRANCH':>12} "
        f"{'CUM PAIR':>10}"
    )

    cumulative_branch = 0
    cumulative_pair = 0

    for k in range(1, len(EXTRA_MODULI) + 1):

        cumulative_branch += branch_dist[k]
        cumulative_pair += pair_dist[k]

        print(
            f"{k:3d} "
            f"{branch_dist[k]:14d} "
            f"{pair_dist[k]:12d} "
            f"{cumulative_branch:12d} "
            f"{cumulative_pair:10d}"
        )

    # =============================================================================
    # DECISIVE MODULI
    # =============================================================================

    print()
    print("=" * 100)
    print("MOST FREQUENTLY DECISIVE MODULI")
    print("=" * 100)

    for r, count in decisive_modulus.most_common():
        print(
            f"r={r:3d} decisive events={count}"
        )

    # =============================================================================
    # EXAMPLES
    # =============================================================================

    print()
    print("=" * 100)
    print("EXAMPLE ANCHORS")
    print("=" * 100)

    shown = 0

    for p, q in anchors:

        if shown >= PRINT_EXAMPLES:
            break

        n = p * q

        if n % 4 != 1:
            continue

        actual_branch = branch_from_mod4(p, q)

        print()
        print(
            f"n={n:,} "
            f"p={p:,} "
            f"q={q:,} "
            f"actual branch={actual_branch}"
        )

        for k in range(1, len(EXTRA_MODULI) + 1):

            prefix = EXTRA_MODULI[:k]

            cand11 = candidates_for_modulus_constraints(
                n,
                (1, 1),
                prefix,
                primes,
                residue_tables,
            )

            cand33 = candidates_for_modulus_constraints(
                n,
                (3, 3),
                prefix,
                primes,
                residue_tables,
            )

            alive = []

            if cand11:
                alive.append("11")

            if cand33:
                alive.append("33")

            actual_candidates = (
                cand11 if actual_branch == (1, 1)
                else cand33
            )

            unique_actual = (
                len(actual_candidates) == 1
                and actual_candidates[0] == (p, q)
            )

            print(
                f"  K={k:2d} "
                f"added={EXTRA_MODULI[k-1]:3d} "
                f"branches={','.join(alive):>5s} "
                f"actual_candidates={len(actual_candidates):4d} "
                f"unique_pair={unique_actual}"
            )

        shown += 1

    # =============================================================================
    # FINAL INTERPRETATION
    # =============================================================================

    print()
    print("=" * 100)
    print("INTERPRETATION")
    print("=" * 100)

    print("""
This experiment measures the information added by moduli that are independent
of M.

For each anchor:

    n = p*q

the fingerprint is

    F_S(p,q) =
        ((p mod r, q mod r) for r in S)

with factor order ignored.

For n == 1 (mod 4), the original modulo-4 ambiguity is:

    11:
        p = 1 (mod 4)
        q = 1 (mod 4)

    33:
        p = 3 (mod 4)
        q = 3 (mod 4)

For every additional modulus r:

    p*q == n (mod r)

therefore, once p is chosen,

    q == n * p^(-1) (mod r).

The experiment records:

    BRANCH FIRST
        smallest K for which only one of 11 or 33 has a prime realization
        in the tested factor interval.

    PAIR FIRST
        smallest K for which exactly one unordered prime pair remains,
        and it is the original (p,q).

The random ordering of the moduli is essential.

It tells us whether the apparent resolution threshold is a real information
requirement or merely an artifact of using 29,31,37,... in one particular
order.

The experiment does NOT claim that K moduli are intrinsically necessary for
factoring an arbitrary integer.

It measures the finite-range discrimination power of this particular
modular-fingerprint construction.
""")

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    main()
