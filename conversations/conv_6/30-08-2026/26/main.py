#!/usr/bin/env python3

import math
import random
import time
from collections import Counter


# ================================================================================================
# CONFIG
# ================================================================================================

M = 111_546_435

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

MODULUS_MIN = 300
MODULUS_MAX = 3_000

CLOSE_RATIO = 0.20

ANCHORS = 300

SEED = 1_511_464_998


# ================================================================================================
# PRIME SIEVE
# ================================================================================================

def sieve(limit):
    a = bytearray(b"\x01") * (limit + 1)

    a[0] = 0
    a[1] = 0

    root = math.isqrt(limit)

    for p in range(2, root + 1):

        if not a[p]:
            continue

        start = p * p
        count = ((limit - start) // p) + 1

        a[start::p] = b"\x00" * count

    return [
        i
        for i in range(2, limit + 1)
        if a[i]
    ]


# ================================================================================================
# ANCHORS
# ================================================================================================

def generate_anchors(primes, count, rng):

    result = []
    seen = set()

    while len(result) < count:

        p, q = rng.sample(primes, 2)

        if p > q:
            p, q = q, p

        if (p, q) in seen:
            continue

        seen.add((p, q))

        result.append(
            (
                p,
                q,
                p * q,
            )
        )

    return result


# ================================================================================================
# SELECT MAXIMUM PRODUCT CLOSE TRIPLE
# ================================================================================================

def select_maximum_product_triple(n, modulus_primes):

    best = None
    best_R = -1

    L = len(modulus_primes)

    for i in range(L):

        r1 = modulus_primes[i]

        for j in range(i + 1, L):

            r2 = modulus_primes[j]

            if r2 > r1 * (1.0 + CLOSE_RATIO):
                break

            for k in range(j + 1, L):

                r3 = modulus_primes[k]

                if r3 > r1 * (1.0 + CLOSE_RATIO):
                    break

                R = r1 * r2 * r3

                if R < n and R > best_R:

                    best_R = R
                    best = (r1, r2, r3)

    return best


# ================================================================================================
# MODULAR INVERSE
# ================================================================================================

def inv_mod(a, p):
    a %= p

    if a == 0:
        return None

    return pow(a, p - 2, p)


# ================================================================================================
# INTEGER INTERVAL FOR QUOTIENT
# ================================================================================================

def quotient_range(residue, modulus):

    lo = math.ceil(
        (FACTOR_MIN - residue) / modulus
    )

    hi = math.floor(
        (FACTOR_MAX - residue) / modulus
    )

    return max(0, lo), hi


# ================================================================================================
# MAIN
# ================================================================================================

def run():

    t0 = time.perf_counter()

    rng = random.Random(SEED)

    print("=" * 100)
    print("THREE-CLOSE-PRIME BILINEAR RESIDUE / QUOTIENT SOLVER")
    print("=" * 100)

    print(f"M                         = {M:,}")
    print(
        f"factor range              = "
        f"{FACTOR_MIN:,} - {FACTOR_MAX:,}"
    )
    print(f"anchors                   = {ANCHORS}")
    print(
        f"modulus prime range       = "
        f"{MODULUS_MIN:,} - {MODULUS_MAX:,}"
    )
    print(f"close ratio               = {CLOSE_RATIO:.0%}")
    print(f"seed                      = {SEED:,}")
    print()

    # ============================================================================================
    # PRIME POOLS
    # ============================================================================================

    print("=" * 100)
    print("BUILDING PRIME POOLS")
    print("=" * 100)

    primes = sieve(FACTOR_MAX)

    factor_primes = [
        p
        for p in primes
        if FACTOR_MIN <= p <= FACTOR_MAX
    ]

    modulus_primes = [
        p
        for p in primes
        if MODULUS_MIN <= p <= MODULUS_MAX
    ]

    prime_set = set(factor_primes)

    print(
        f"factor primes             = "
        f"{len(factor_primes):,}"
    )

    print(
        f"modulus primes            = "
        f"{len(modulus_primes):,}"
    )

    # ============================================================================================
    # ANCHORS
    # ============================================================================================

    anchors = generate_anchors(
        factor_primes,
        ANCHORS,
        rng,
    )

    print()
    print("=" * 100)
    print("SELECTING MAXIMUM-PRODUCT CLOSE TRIPLES")
    print("=" * 100)

    records = []

    for i, (p, q, n) in enumerate(anchors, 1):

        mods = select_maximum_product_triple(
            n,
            modulus_primes,
        )

        if mods is None:
            continue

        r1, r2, r3 = mods

        R = r1 * r2 * r3
        g = n - R

        records.append(
            {
                "id": i,
                "p": p,
                "q": q,
                "n": n,
                "r1": r1,
                "r2": r2,
                "r3": r3,
                "R": R,
                "g": g,
            }
        )

        if i % 25 == 0:
            print(
                f"anchor {i:3d}/{len(anchors)}"
            )

    print()
    print(
        f"usable anchors = {len(records)}"
    )

    # ============================================================================================
    # GLOBAL COUNTERS
    # ============================================================================================

    total_a = 0
    total_ak_states = 0

    total_r2_solutions = 0
    total_r3_solutions = 0

    total_interval_candidates = 0
    total_prime_candidates = 0
    total_exact = 0

    successful = 0

    r2_candidate_distribution = Counter()
    r3_candidate_distribution = Counter()

    per_anchor = []

    # ============================================================================================
    # MAIN EXPERIMENT
    # ============================================================================================

    print()
    print("=" * 100)
    print("RUNNING BILINEAR QUOTIENT SOLVER")
    print("=" * 100)

    for idx, rec in enumerate(records, 1):

        p_actual = rec["p"]
        q_actual = rec["q"]
        n = rec["n"]

        r1 = rec["r1"]
        r2 = rec["r2"]
        r3 = rec["r3"]

        R = rec["R"]
        g = rec["g"]

        d12 = r2 - r1

        # ----------------------------------------------------------------------------------------
        # r1 residue -> q residue
        # ----------------------------------------------------------------------------------------

        anchor_a = p_actual % r1
        anchor_b = q_actual % r1

        a_count = 0
        ak_count = 0
        r2_hits = 0
        r3_hits = 0
        interval_hits = 0
        prime_hits = 0
        exact_hits = 0

        recovered = False

        r2_candidates = []
        r3_candidates = []
        exact_candidates = []

        # ----------------------------------------------------------------------------------------
        # Enumerate p residue a modulo r1.
        #
        # Since r1 < p,q, neither factor can be equal to r1.
        # We nevertheless allow all nonzero residues.
        # ----------------------------------------------------------------------------------------

        for a in range(1, r1):

            total_a += 1
            a_count += 1

            inv_a = inv_mod(a, r1)

            if inv_a is None:
                continue

            # pq = g mod r1
            b = (g * inv_a) % r1

            if b == 0:
                continue

            # p = a + k*r1
            # q = b + l*r1

            k_lo, k_hi = quotient_range(
                a,
                r1,
            )

            l_lo, l_hi = quotient_range(
                b,
                r1,
            )

            # ------------------------------------------------------------------------------------
            # For every possible k:
            #
            # p mod r2 = a - k*d12
            #
            # q mod r2 = b - l*d12
            #
            # We require:
            #
            # (a-kd)(b-ld) = g mod r2
            #
            # => (-x*d) l = g - x*b mod r2
            #
            # where x = a-kd.
            # ------------------------------------------------------------------------------------

            for k in range(k_lo, k_hi + 1):

                ak_count += 1
                total_ak_states += 1

                x = (
                    a - k * d12
                ) % r2

                rhs = (
                    g - x * b
                ) % r2

                coeff = (
                    -x * d12
                ) % r2

                # r2 is prime.
                #
                # If coeff == 0, then either every l or no l is possible.
                # For a genuine factor residue state, x == 0 would imply
                # p ≡ 0 mod r2, which cannot happen for a prime p != r2.
                if coeff == 0:

                    if rhs != 0:
                        continue

                    # Degenerate case.
                    #
                    # Fall back to bounded l scan.
                    #
                    # This should almost never happen.
                    for l in range(
                        l_lo,
                        l_hi + 1,
                    ):

                        p = a + k * r1
                        q = b + l * r1

                        if not (
                            FACTOR_MIN <= p <= FACTOR_MAX
                            and
                            FACTOR_MIN <= q <= FACTOR_MAX
                        ):
                            continue

                        if p > q:
                            continue

                        if p * q != n:
                            continue

                        r2_hits += 1
                        r3_hits += 1
                        exact_hits += 1

                        exact_candidates.append(
                            (p, q)
                        )

                    continue

                inv_coeff = inv_mod(
                    coeff,
                    r2,
                )

                if inv_coeff is None:
                    continue

                # l modulo r2.
                l0 = (
                    rhs * inv_coeff
                ) % r2

                # Since l is much smaller than r2,
                # at most one bounded l exists.
                #
                # Also allow l0 + t*r2 in case the quotient interval
                # unexpectedly exceeds one modulus.
                t_min = math.ceil(
                    (l_lo - l0) / r2
                )

                t_max = math.floor(
                    (l_hi - l0) / r2
                )

                for t in range(
                    t_min,
                    t_max + 1,
                ):

                    l = l0 + t * r2

                    if not (
                        l_lo <= l <= l_hi
                    ):
                        continue

                    r2_hits += 1

                    p = a + k * r1
                    q = b + l * r1

                    if not (
                        FACTOR_MIN <= p <= FACTOR_MAX
                        and
                        FACTOR_MIN <= q <= FACTOR_MAX
                    ):
                        continue

                    if p > q:
                        continue

                    # --------------------------------------------------------------------------------
                    # Third modulus.
                    # --------------------------------------------------------------------------------

                    if (
                        (p * q - g) % r3
                        != 0
                    ):
                        continue

                    r3_hits += 1

                    r3_candidates.append(
                        (p, q)
                    )

                    # --------------------------------------------------------------------------------
                    # Prime test.
                    # --------------------------------------------------------------------------------

                    interval_hits += 1

                    if p not in prime_set:
                        continue

                    if q not in prime_set:
                        continue

                    prime_hits += 1

                    # --------------------------------------------------------------------------------
                    # Exact product.
                    # --------------------------------------------------------------------------------

                    if p * q != n:
                        continue

                    exact_hits += 1

                    exact_candidates.append(
                        (p, q)
                    )

                    if (
                        p == p_actual
                        and
                        q == q_actual
                    ):
                        recovered = True

        # ----------------------------------------------------------------------------------------
        # Count actual state explicitly.
        # ----------------------------------------------------------------------------------------

        actual_r1 = anchor_a
        actual_r2_p = p_actual % r2
        actual_r2_q = q_actual % r2

        actual_r3 = (
            (p_actual * q_actual - g)
            % r3
        )

        if actual_r3 != 0:
            print(
                "ERROR: actual factor failed r3"
            )

        # ----------------------------------------------------------------------------------------
        # Distributions
        # ----------------------------------------------------------------------------------------

        r2_candidate_distribution[r2_hits] += 1
        r3_candidate_distribution[r3_hits] += 1

        total_r2_solutions += r2_hits
        total_r3_solutions += r3_hits

        total_interval_candidates += interval_hits
        total_prime_candidates += prime_hits
        total_exact += exact_hits

        if recovered:
            successful += 1

        per_anchor.append(
            {
                "id": rec["id"],
                "p": p_actual,
                "q": q_actual,
                "n": n,
                "r1": r1,
                "r2": r2,
                "r3": r3,
                "R": R,
                "g": g,
                "a_states": a_count,
                "ak_states": ak_count,
                "r2": r2_hits,
                "r3": r3_hits,
                "interval": interval_hits,
                "prime": prime_hits,
                "exact": exact_hits,
                "recovered": recovered,
            }
        )

        if idx % 25 == 0:
            print(
                f"anchor {idx:3d}/{len(records)}"
            )

    # ============================================================================================
    # SUMMARY
    # ============================================================================================

    N = len(records)

    print()
    print("=" * 100)
    print("SUMMARY")
    print("=" * 100)

    print(
        f"anchors analyzed                = {N}"
    )

    print(
        f"average r1 residue states       = "
        f"{total_a / N:,.3f}"
    )

    print(
        f"average (a,k) states            = "
        f"{total_ak_states / N:,.3f}"
    )

    print(
        f"average r2 solutions             = "
        f"{total_r2_solutions / N:,.3f}"
    )

    print(
        f"average r3 solutions             = "
        f"{total_r3_solutions / N:,.3f}"
    )

    print(
        f"average interval candidates      = "
        f"{total_interval_candidates / N:,.3f}"
    )

    print(
        f"average prime candidates         = "
        f"{total_prime_candidates / N:,.3f}"
    )

    print(
        f"average exact solutions          = "
        f"{total_exact / N:,.3f}"
    )

    print()

    # ============================================================================================
    # BASELINE COMPARISON
    # ============================================================================================

    baseline = len(factor_primes)

    print("=" * 100)
    print("BASELINE COMPARISON")
    print("=" * 100)

    print(
        f"ordinary prime enumeration       = "
        f"{baseline:,}"
    )

    print(
        f"average (a,k) states             = "
        f"{total_ak_states / N:,.3f}"
    )

    print(
        f"average r3 states                = "
        f"{total_r3_solutions / N:,.3f}"
    )

    print(
        f"average prime candidates         = "
        f"{total_prime_candidates / N:,.3f}"
    )

    print()

    print(
        "Relative to exhaustive factor-prime enumeration:"
    )

    print(
        f"(a,k) state ratio                = "
        f"{(total_ak_states / N) / baseline:.6f}"
    )

    print(
        f"r3 state ratio                   = "
        f"{(total_r3_solutions / N) / baseline:.6f}"
    )

    print(
        f"prime candidate ratio            = "
        f"{(total_prime_candidates / N) / baseline:.6f}"
    )

    # ============================================================================================
    # REDUCTION
    # ============================================================================================

    print()
    print("=" * 100)
    print("SEARCH REDUCTION")
    print("=" * 100)

    avg_ak = total_ak_states / N
    avg_r2 = total_r2_solutions / N
    avg_r3 = total_r3_solutions / N
    avg_prime = total_prime_candidates / N

    print(
        f"r2 bilinear reduction            = "
        f"{100.0 * (1.0 - avg_r2 / avg_ak):.6f}%"
    )

    print(
        f"r3 additional reduction          = "
        f"{100.0 * (1.0 - avg_r3 / avg_r2):.6f}%"
        if avg_r2
        else
        "r3 additional reduction          = N/A"
    )

    print(
        f"prime filtering reduction        = "
        f"{100.0 * (1.0 - avg_prime / avg_r3):.6f}%"
        if avg_r3
        else
        "prime filtering reduction        = N/A"
    )

    print()

    # ============================================================================================
    # DISTRIBUTIONS
    # ============================================================================================

    print("=" * 100)
    print("R2 SOLUTION DISTRIBUTION")
    print("=" * 100)

    for value, count in sorted(
        r2_candidate_distribution.items()
    ):

        print(
            f"r2 solutions = "
            f"{value:5d} anchors = {count:4d}"
        )

    print()
    print("=" * 100)
    print("R3 SOLUTION DISTRIBUTION")
    print("=" * 100)

    for value, count in sorted(
        r3_candidate_distribution.items()
    ):

        print(
            f"r3 solutions = "
            f"{value:5d} anchors = {count:4d}"
        )

    # ============================================================================================
    # RECOVERY
    # ============================================================================================

    print()
    print("=" * 100)
    print("RECOVERY")
    print("=" * 100)

    print(
        f"correctly recovered             = "
        f"{successful}/{N}"
    )

    print(
        f"recovery rate                    = "
        f"{100.0 * successful / N:.4f}%"
    )

    # ============================================================================================
    # STRONGEST CASES
    # ============================================================================================

    print()
    print("=" * 100)
    print("STRONGEST BILINEAR COLLAPSES")
    print("=" * 100)

    strongest = sorted(
        per_anchor,
        key=lambda x: (
            x["r3"],
            x["r2"],
            x["ak_states"],
        )
    )

    for row in strongest[:20]:

        print(
            f"n={row['n']:,} "
            f"p={row['p']:,} "
            f"q={row['q']:,} "
            f"mods=({row['r1']},{row['r2']},{row['r3']}) "
            f"R/n={row['R']/row['n']:.10f} "
            f"(a,k)={row['ak_states']} "
            f"r2={row['r2']} "
            f"r3={row['r3']} "
            f"prime={row['prime']} "
            f"exact={row['exact']}"
        )

    # ============================================================================================
    # WORST CASES
    # ============================================================================================

    print()
    print("=" * 100)
    print("WEAKEST BILINEAR COLLAPSES")
    print("=" * 100)

    weakest = sorted(
        per_anchor,
        key=lambda x: (
            -x["r3"],
            -x["r2"],
        )
    )

    for row in weakest[:20]:

        print(
            f"n={row['n']:,} "
            f"p={row['p']:,} "
            f"q={row['q']:,} "
            f"mods=({row['r1']},{row['r2']},{row['r3']}) "
            f"R/n={row['R']/row['n']:.10f} "
            f"(a,k)={row['ak_states']} "
            f"r2={row['r2']} "
            f"r3={row['r3']} "
            f"prime={row['prime']} "
            f"exact={row['exact']}"
        )

    # ============================================================================================
    # ANCHOR TABLE
    # ============================================================================================

    print()
    print("=" * 100)
    print("ANCHOR RESULTS")
    print("=" * 100)

    print(
        " ID        p        q       r1    r2    r3 "
        "      (a,k)      R2      R3    PRIME EXACT"
    )

    print("-" * 100)

    for row in per_anchor:

        print(
            f"{row['id']:3d} "
            f"{row['p']:8,d} "
            f"{row['q']:8,d} "
            f"{row['r1']:6d} "
            f"{row['r2']:6d} "
            f"{row['r3']:6d} "
            f"{row['ak_states']:10,d} "
            f"{row['r2']:8,d} "
            f"{row['r3']:7,d} "
            f"{row['prime']:7,d} "
            f"{row['exact']:5,d}"
        )

    # ============================================================================================
    # ACTUAL FACTOR STATE
    # ============================================================================================

    print()
    print("=" * 100)
    print("ACTUAL FACTOR RESIDUE STATE")
    print("=" * 100)

    failures = 0

    for row in per_anchor:

        p = row["p"]
        q = row["q"]

        r1 = row["r1"]
        r2 = row["r2"]
        r3 = row["r3"]

        g = (
            row["n"] -
            row["R"]
        )

        a = p % r1
        b = q % r1

        inv_a = inv_mod(
            a,
            r1,
        )

        if inv_a is None:
            failures += 1
            continue

        b_expected = (
            g * inv_a
        ) % r1

        if b_expected != b:
            failures += 1
            continue

        if (
            (p * q - g) % r2
            != 0
        ):
            failures += 1
            continue

        if (
            (p * q - g) % r3
            != 0
        ):
            failures += 1
            continue

    print(
        f"actual residue-state failures = "
        f"{failures}"
    )

    # ============================================================================================
    # MATHEMATICAL INTERPRETATION
    # ============================================================================================

    print()
    print("=" * 100)
    print("MATHEMATICAL INTERPRETATION")
    print("=" * 100)

    print()
    print("We start from:")
    print()
    print("    p*q = n = R + g")
    print()

    print("Therefore, for every selected modulus:")
    print()
    print("    p*q ≡ g (mod r_i)")
    print()

    print("Choose:")
    print()
    print("    a = p mod r1")
    print("    b = q mod r1")
    print()

    print("Then:")
    print()
    print("    a*b ≡ g (mod r1)")
    print()

    print("and therefore:")
    print()
    print("    b ≡ g*a^(-1) (mod r1)")
    print()

    print("Write:")
    print()
    print("    p = a + k*r1")
    print("    q = b + l*r1")
    print()

    print(
        "Because r2 is close to r1, let:"
    )

    print()
    print("    delta = r2-r1")
    print()

    print("Then:")
    print()
    print("    p mod r2 = a-k*delta (mod r2)")
    print("    q mod r2 = b-l*delta (mod r2)")
    print()

    print("The product condition becomes:")
    print()
    print(
        "    (a-k*delta)(b-l*delta) ≡ g (mod r2)"
    )
    print()

    print("For fixed a and k this is LINEAR in l:")
    print()
    print(
        "    (-x*delta) l ≡ g-x*b (mod r2)"
    )
    print()
    print("where:")
    print()
    print(
        "    x = a-k*delta"
    )
    print()

    print(
        "Thus the second close modulus does not merely act as "
        "another yes/no QR filter."
    )

    print(
        "It can algebraically solve one quotient coordinate."
    )

    print()
    print("The third modulus is then used as an independent:")
    print()
    print("    bilinear consistency check")
    print()

    print("The important statistic is consequently:")
    print()
    print("    (a,k) states")
    print("          ->")
    print("    r2 solutions")
    print("          ->")
    print("    r3 solutions")
    print("          ->")
    print("    prime candidates")
    print("          ->")
    print("    exact factor pairs")
    print()

    print(
        "Unlike the QR experiment, this construction uses the actual "
        "multiplicative relation pq ≡ g and the closeness "
        "r2-r1 directly."
    )

    print()
    print(
        "A particularly interesting outcome would be if r3 leaves "
        "only O(1) states per anchor while the number of (a,k) states "
        "remains comparable to the size of the factor interval."
    )

    print()
    print(
        "That would indicate that the close-modulus relation is doing "
        "more than ordinary independent modular filtering."
    )

    print()
    print(
        "If the r3 state count remains large, then the construction "
        "is primarily a change of coordinates rather than a new "
        "factorization mechanism."
    )

    # ============================================================================================
    # TIMING
    # ============================================================================================

    elapsed = time.perf_counter() - t0

    print()
    print("=" * 100)
    print("TIMING")
    print("=" * 100)

    print(
        f"total runtime = {elapsed:.3f} seconds"
    )

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()
