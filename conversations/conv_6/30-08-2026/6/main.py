#!/usr/bin/env python3

import math
import random
import numpy as np

# =============================================================================
# CONFIGURATION
# =============================================================================

M = 111_546_435

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

TRIALS = 300

# We do not need 200 random orders.
# 40 orders is enough to see whether modulus order matters.
ORDER_REPETITIONS = 40

MODULI = [
    29, 31, 37, 41,
    43, 47, 53, 59,
    61, 67, 71, 73,
    79, 83, 89, 97,
]

SEED = 1_511_464_998


# =============================================================================
# PRIME TEST
# =============================================================================

def is_prime(n):
    if n < 2:
        return False

    if n == 2:
        return True

    if n % 2 == 0:
        return False

    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2

    return True


# =============================================================================
# PRIME POOL
# =============================================================================

def build_primes():
    return np.array(
        [
            p
            for p in range(FACTOR_MIN, FACTOR_MAX + 1)
            if is_prime(p)
        ],
        dtype=np.int64,
    )


# =============================================================================
# ACTUAL ANCHORS
# =============================================================================

def build_anchors(primes, count, rng):

    anchors = set()

    while len(anchors) < count:

        i = rng.randrange(len(primes))
        j = rng.randrange(len(primes))

        if i == j:
            continue

        p = int(primes[i])
        q = int(primes[j])

        if p > q:
            p, q = q, p

        anchors.add((p, q))

    return list(sorted(anchors))


# =============================================================================
# BRANCH
# =============================================================================

def branch4(p, q):

    rp = p & 3
    rq = q & 3

    if rp == 1 and rq == 1:
        return 11

    if rp == 3 and rq == 3:
        return 33

    return 13


# =============================================================================
# RESIDUE ARRAYS
# =============================================================================

def build_residue_arrays(primes):

    arrays = {}

    for r in MODULI:
        arrays[r] = primes % r

    return arrays


# =============================================================================
# PRIME LOOKUP
#
# Maps:
#
#   residue tuple -> whether a prime exists
#
# We only need existence for the branch test.
# =============================================================================

def build_signature_sets(primes):

    tables = {}

    for r in MODULI:

        residues = primes % r

        s = set(int(x) for x in residues)

        tables[r] = s

    return tables


# =============================================================================
# KEY OPTIMIZATION
#
# For a fixed anchor n and a fixed branch:
#
#     p*q == n (mod r)
#
# implies
#
#     q == n * inverse(p) (mod r)
#
# We do NOT regenerate candidate lists.
#
# Instead we calculate the required q residue vector for ALL branch primes
# simultaneously with NumPy.
#
# After enough moduli have been added, the CRT modulus exceeds 100,000.
# At that point the residue vector uniquely determines q in our factor range.
# =============================================================================

def branch_candidates(
    n,
    branch,
    prime_values,
    residue_arrays,
    signature_prime_sets,
    used_moduli,
):

    if branch == 11:

        mask = (prime_values % 4) == 1

    elif branch == 33:

        mask = (prime_values % 4) == 3

    else:

        return []

    pvals = prime_values[mask]

    if len(pvals) == 0:
        return []

    # -------------------------------------------------------------------------
    # For every possible p, determine whether a compatible q residue exists.
    # -------------------------------------------------------------------------

    possible = np.ones(len(pvals), dtype=bool)

    required_residues = []

    for r in used_moduli:

        rp = pvals % r
        nn = n % r

        inv = np.array(
            [pow(int(x), -1, r) for x in rp],
            dtype=np.int64,
        )

        rq = (nn * inv) % r

        required_residues.append(rq)

        # Since every r is prime and p is > r,
        # q residue must be non-zero.
        possible &= rq != 0

    if not np.any(possible):
        return []

    candidate_indices = np.nonzero(possible)[0]

    # -------------------------------------------------------------------------
    # We now test the required q residue tuple.
    #
    # Python sets are used only after NumPy has reduced the search space.
    # -------------------------------------------------------------------------

    out = []

    # CRT uniqueness threshold:
    crt_modulus = 1

    for r in used_moduli:
        crt_modulus *= r

    # -------------------------------------------------------------------------
    # If CRT modulus exceeds factor range, reconstruct q directly.
    # -------------------------------------------------------------------------

    if crt_modulus > FACTOR_MAX:

        # Precompute CRT reconstruction coefficients.
        #
        # q = sum(a_i * inv_i * M_i) mod CRT
        #
        # where:
        #
        #   M_i = CRT/r_i
        #
        coeffs = []

        for r in used_moduli:

            Mi = crt_modulus // r
            inv = pow(Mi, -1, r)

            coeffs.append(Mi * inv)

        coeffs = np.array(coeffs, dtype=object)

        for idx in candidate_indices:

            q = 0

            for j in range(len(used_moduli)):
                q += int(required_residues[j][idx]) * int(coeffs[j])

            q %= crt_modulus

            if q < FACTOR_MIN or q > FACTOR_MAX:
                continue

            if q == int(pvals[idx]):
                continue

            if not is_prime(q):
                continue

            a = int(pvals[idx])
            b = int(q)

            pair = tuple(sorted((a, b)))

            if (pair[0] * pair[1]) % math.prod(used_moduli) != (
                n % math.prod(used_moduli)
            ):
                continue

            out.append(pair)

        return sorted(set(out))

    # -------------------------------------------------------------------------
    # Before CRT exceeds the factor range, explicitly test tuples.
    # -------------------------------------------------------------------------

    prime_sets = []

    for r in used_moduli:

        prime_sets.append(
            set(
                int(x)
                for x in prime_values[
                    (prime_values % r) != 0
                ]
                % r
            )
        )

    for idx in candidate_indices:

        # First build required tuple.
        ok = True

        for j, r in enumerate(used_moduli):

            residue = int(required_residues[j][idx])

            if residue == 0:
                ok = False
                break

            if residue not in prime_sets[j]:
                ok = False
                break

        if not ok:
            continue

        # At this low-K stage simply check all primes of the correct
        # residue intersection.
        #
        # This stage disappears quickly because:
        #
        # 29*31*37*41 > 100,000.
        #
        possible_q = []

        for q in prime_values:

            if q == pvals[idx]:
                continue

            good = True

            for j, r in enumerate(used_moduli):

                if int(q % r) != int(required_residues[j][idx]):
                    good = False
                    break

            if good:
                possible_q.append(int(q))

        for q in possible_q:

            pair = tuple(
                sorted((int(pvals[idx]), q))
            )

            out.append(pair)

    return sorted(set(out))


# =============================================================================
# FAST BRANCH TEST
# =============================================================================

def branch_survival_fast(
    n,
    used_moduli,
    primes,
    residue_arrays,
    signature_prime_sets,
):

    c11 = branch_candidates(
        n,
        11,
        primes,
        residue_arrays,
        signature_prime_sets,
        used_moduli,
    )

    c33 = branch_candidates(
        n,
        33,
        primes,
        residue_arrays,
        signature_prime_sets,
        used_moduli,
    )

    return c11, c33


# =============================================================================
# MAIN
# =============================================================================

def main():

    rng = random.Random(SEED)

    print("=" * 100)
    print("FAST MULTI-MODULUS BRANCH DISCRIMINATION / CRT THRESHOLD")
    print("=" * 100)

    print(f"M                    = {M:,}")
    print(
        f"factor range         = "
        f"{FACTOR_MIN:,} - {FACTOR_MAX:,}"
    )
    print(f"actual anchors       = {TRIALS}")
    print(f"random orders        = {ORDER_REPETITIONS}")
    print(f"moduli               = {MODULI}")
    print(f"seed                 = {SEED:,}")

    # -------------------------------------------------------------------------
    # PRIME POOL
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("BUILDING PRIME POOL")
    print("=" * 100)

    primes = build_primes()

    print(f"factor primes = {len(primes):,}")

    # -------------------------------------------------------------------------
    # ANCHORS
    # -------------------------------------------------------------------------

    anchors = build_anchors(
        primes,
        TRIALS,
        rng,
    )

    print(f"anchors = {len(anchors)}")

    # -------------------------------------------------------------------------
    # PRECOMPUTATION
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("PRECOMPUTING RESIDUES")
    print("=" * 100)

    residue_arrays = build_residue_arrays(primes)

    signature_prime_sets = build_signature_sets(primes)

    for r in MODULI:
        print(
            f"r={r:3d} "
            f"residue classes={len(signature_prime_sets[r]):3d}"
        )

    # -------------------------------------------------------------------------
    # CRT THRESHOLD
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("CRT MODULUS GROWTH")
    print("=" * 100)

    product = 1

    for k, r in enumerate(MODULI, 1):

        product *= r

        print(
            f"K={k:2d} "
            f"added={r:3d} "
            f"CRT product={product:,}"
            f"{'  <-- > FACTOR_MAX' if product > FACTOR_MAX else ''}"
        )

    # -------------------------------------------------------------------------
    # RANDOM ORDERS
    # -------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("RANDOM MODULUS ORDER EXPERIMENT")
    print("=" * 100)

    branch_resolve_counts = np.zeros(
        len(MODULI) + 1,
        dtype=np.int64,
    )

    pair_resolve_counts = np.zeros(
        len(MODULI) + 1,
        dtype=np.int64,
    )

    never_branch = 0
    never_pair = 0

    # How often each modulus is the actual decisive modulus.
    decisive_modulus = {
        r: 0
        for r in MODULI
    }

    example_results = []

    # Only n == 1 mod 4 has the 11/33 ambiguity.
    relevant_anchors = []

    for p, q in anchors:

        n = p * q

        if n % 4 == 1:
            relevant_anchors.append((p, q))

    print(
        f"anchors with n == 1 (mod 4) = "
        f"{len(relevant_anchors)}"
    )

    # -------------------------------------------------------------------------
    # RUN
    # -------------------------------------------------------------------------

    for order_id in range(ORDER_REPETITIONS):

        order = MODULI[:]
        rng.shuffle(order)

        if (order_id + 1) % 5 == 0:
            print(
                f"order {order_id + 1:3d}/{ORDER_REPETITIONS}"
            )

        for p, q in relevant_anchors:

            n = p * q

            actual_branch = branch4(p, q)

            branch_first = None
            pair_first = None

            for k in range(1, len(order) + 1):

                used = order[:k]

                c11, c33 = branch_survival_fast(
                    n,
                    used,
                    primes,
                    residue_arrays,
                    signature_prime_sets,
                )

                alive = []

                if c11:
                    alive.append(11)

                if c33:
                    alive.append(33)

                # -------------------------------------------------------------
                # Branch discrimination
                # -------------------------------------------------------------

                if branch_first is None:

                    if len(alive) == 1:

                        branch_first = k
                        decisive_modulus[order[k - 1]] += 1

                # -------------------------------------------------------------
                # Actual branch pair uniqueness
                # -------------------------------------------------------------

                actual_candidates = (
                    c11
                    if actual_branch == 11
                    else c33
                )

                if pair_first is None:

                    if (
                        len(actual_candidates) == 1
                        and actual_candidates[0] == (p, q)
                    ):
                        pair_first = k

            # -----------------------------------------------------------------
            # Record
            # -----------------------------------------------------------------

            if branch_first is None:
                never_branch += 1
            else:
                branch_resolve_counts[branch_first] += 1

            if pair_first is None:
                never_pair += 1
            else:
                pair_resolve_counts[pair_first] += 1

            if len(example_results) < 12:

                example_results.append(
                    (
                        p,
                        q,
                        actual_branch,
                        branch_first,
                        pair_first,
                    )
                )

    # =============================================================================
    # RESULTS
    # =============================================================================

    repetitions_total = (
        len(relevant_anchors) *
        ORDER_REPETITIONS
    )

    print()
    print("=" * 100)
    print("RESULTS")
    print("=" * 100)

    print(
        f"branch observations = {repetitions_total:,}"
    )

    print(
        f"branch never resolved = "
        f"{never_branch:,}"
    )

    print(
        f"pair never resolved   = "
        f"{never_pair:,}"
    )

    print()
    print(
        f"{'K':>3} "
        f"{'MODULUS':>8} "
        f"{'BRANCH FIRST':>16} "
        f"{'PAIR FIRST':>13} "
        f"{'CUM BRANCH %':>15} "
        f"{'CUM PAIR %':>13}"
    )

    cumulative_branch = 0
    cumulative_pair = 0

    for k, r in enumerate(MODULI, 1):

        cumulative_branch += int(
            branch_resolve_counts[k]
        )

        cumulative_pair += int(
            pair_resolve_counts[k]
        )

        branch_pct = (
            100.0 *
            cumulative_branch /
            repetitions_total
        )

        pair_pct = (
            100.0 *
            cumulative_pair /
            repetitions_total
        )

        print(
            f"{k:3d} "
            f"{r:8d} "
            f"{branch_resolve_counts[k]:16d} "
            f"{pair_resolve_counts[k]:13d} "
            f"{branch_pct:14.2f}% "
            f"{pair_pct:12.2f}%"
        )

    # =============================================================================
    # DECISIVE MODULI
    # =============================================================================

    print()
    print("=" * 100)
    print("DECISIVE MODULUS FREQUENCY")
    print("=" * 100)

    total_decisive = sum(
        decisive_modulus.values()
    )

    for r in MODULI:

        count = decisive_modulus[r]

        pct = (
            100.0 * count / total_decisive
            if total_decisive
            else 0.0
        )

        print(
            f"r={r:3d} "
            f"count={count:6d} "
            f"share={pct:8.3f}%"
        )

    # =============================================================================
    # EXAMPLES
    # =============================================================================

    print()
    print("=" * 100)
    print("EXAMPLE ANCHORS")
    print("=" * 100)

    for p, q, branch, branch_first, pair_first in example_results:

        print()
        print(
            f"p={p:,} "
            f"q={q:,} "
            f"n={p*q:,} "
            f"branch={branch}"
        )

        print(
            f"  first branch resolution = "
            f"{branch_first}"
        )

        print(
            f"  first pair resolution   = "
            f"{pair_first}"
        )

    # =============================================================================
    # THEORETICAL INTERPRETATION
    # =============================================================================

    print()
    print("=" * 100)
    print("INTERPRETATION")
    print("=" * 100)

    print("""
This experiment returns to the original idea:

    n mod 4 = 1

gives two possible same-residue factor branches:

    11:
        p = 1 (mod 4)
        q = 1 (mod 4)

    33:
        p = 3 (mod 4)
        q = 3 (mod 4)

For every additional modulus r:

    p*q = n (mod r)

and therefore:

    q = n * p^(-1) (mod r).

The key distinction is between two questions.

QUESTION A:
    Does the extra modulus eliminate one modulo-4 branch?

QUESTION B:
    Does the accumulated fingerprint identify the exact
    factor pair inside the tested factor range?

These are NOT the same question.

The important structural threshold is the product of the independent
moduli.

Once:

    R = product(r)

satisfies

    R > 100,000,

a complete residue vector modulo those r uniquely identifies an individual
factor from the tested interval 10,000..100,000.

That does NOT automatically mean the factorization is solved, because we
still need compatibility between p and q.

This experiment therefore measures exactly where that compatibility causes
the wrong branch to disappear.

Randomizing the order is important because otherwise an observed threshold
could simply be caused by choosing a particularly powerful modulus early.

The final comparison to make is:

    branch resolution threshold
        versus
    factor-pair resolution threshold.

If branch resolution happens consistently much earlier than pair
resolution, then the modulo-4 idea is giving real branch discrimination but
not complete factor recovery.

If both occur at essentially the same threshold, then the factor
fingerprint itself is doing most of the work.

If neither happens reliably, then the additional moduli are not sufficient
to discriminate the branch in this factor interval.
""")

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    main()
