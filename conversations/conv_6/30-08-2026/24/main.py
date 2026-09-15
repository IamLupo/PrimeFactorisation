#!/usr/bin/env python3

import math
import random
import time
from collections import Counter


# ==================================================================================================
# CONFIG
# ==================================================================================================

M = 111_546_435

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

MODULUS_MIN = 300
MODULUS_MAX = 3_000

CLOSE_RATIO = 0.20

ANCHORS = 300

SEED = 1_511_464_998

# The factor interval implies:
MAX_FACTOR_DIFF = FACTOR_MAX - FACTOR_MIN

# Both factors are odd primes in this interval, so p+q and |p-q|
# are both even.
SUM_MIN = 2 * FACTOR_MIN
SUM_MAX = 2 * FACTOR_MAX


# ==================================================================================================
# SIEVE
# ==================================================================================================

def sieve(limit):
    a = bytearray(b"\x01") * (limit + 1)

    a[0] = 0
    a[1] = 0

    root = math.isqrt(limit)

    for p in range(2, root + 1):
        if a[p]:
            start = p * p
            count = ((limit - start) // p) + 1
            a[start::p] = b"\x00" * count

    return [i for i in range(2, limit + 1) if a[i]]


# ==================================================================================================
# ANCHORS
# ==================================================================================================

def generate_anchors(primes, count, rng):
    seen = set()
    result = []

    while len(result) < count:

        p, q = rng.sample(primes, 2)

        if p > q:
            p, q = q, p

        if (p, q) in seen:
            continue

        seen.add((p, q))
        result.append((p, q, p * q))

    return result


# ==================================================================================================
# CLOSE-PRIME TRIPLE
# ==================================================================================================

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


# ==================================================================================================
# QUADRATIC RESIDUE ROOT TABLE
# ==================================================================================================

def build_sqrt_table(p):
    """
    table[x] contains all y in [0,p) satisfying:

        y^2 == x (mod p)

    For prime p there are normally 0 or 2 roots, with one root
    for x=0.
    """

    table = [[] for _ in range(p)]

    for y in range(p):
        table[(y * y) % p].append(y)

    return table


# ==================================================================================================
# CRT FOR TWO MODULI
# ==================================================================================================

def crt2(a, m, b, n):
    """
    x == a (mod m)
    x == b (mod n)

    Returns x in [0,m*n).
    """

    # x = a + m*t
    #
    # a + m*t == b (mod n)
    #
    # m*t == b-a (mod n)

    t = ((b - a) * pow(m, -1, n)) % n

    return a + m * t


# ==================================================================================================
# EXACT FACTOR CHECK
# ==================================================================================================

def recover_from_sum_difference(n, s, d):
    if d < 0:
        return None

    if (s - d) & 1:
        return None

    p = (s - d) // 2
    q = (s + d) // 2

    if not (FACTOR_MIN <= p <= FACTOR_MAX):
        return None

    if not (FACTOR_MIN <= q <= FACTOR_MAX):
        return None

    if p * q != n:
        return None

    return p, q


# ==================================================================================================
# MAIN
# ==================================================================================================

def run():

    t0 = time.perf_counter()

    rng = random.Random(SEED)

    print("=" * 100)
    print("THREE-CLOSE-PRIME MODULAR-DIFFERENCE RECONSTRUCTION EXPERIMENT")
    print("=" * 100)

    print(f"M                         = {M:,}")
    print(f"factor range              = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"anchors                   = {ANCHORS}")
    print(f"modulus prime range       = {MODULUS_MIN:,} - {MODULUS_MAX:,}")
    print(f"close ratio               = {CLOSE_RATIO:.0%}")
    print(f"seed                      = {SEED:,}")
    print()

    # ==============================================================================================
    # PRIME POOLS
    # ==============================================================================================

    print("=" * 100)
    print("BUILDING PRIME POOLS")
    print("=" * 100)

    primes = sieve(FACTOR_MAX)

    factor_primes = [
        p for p in primes
        if FACTOR_MIN <= p <= FACTOR_MAX
    ]

    modulus_primes = [
        p for p in primes
        if MODULUS_MIN <= p <= MODULUS_MAX
    ]

    print(f"factor primes              = {len(factor_primes):,}")
    print(f"modulus primes             = {len(modulus_primes):,}")
    print()

    # ==============================================================================================
    # ROOT TABLE CACHE
    # ==============================================================================================

    print("=" * 100)
    print("PREPARING MODULAR SQUARE-ROOT TABLES")
    print("=" * 100)

    sqrt_cache = {}

    # We don't yet know which moduli will be selected, so tables are
    # generated lazily and cached.

    # ==============================================================================================
    # ANCHORS
    # ==============================================================================================

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

    for idx, (p, q, n) in enumerate(anchors, 1):

        mods = select_maximum_product_triple(
            n,
            modulus_primes,
        )

        if mods is None:
            continue

        r1, r2, r3 = mods

        R = r1 * r2 * r3
        gap = n - R

        records.append(
            {
                "id": idx,
                "p": p,
                "q": q,
                "n": n,
                "mods": mods,
                "R": R,
                "gap": gap,
            }
        )

        if idx % 25 == 0:
            print(f"anchor {idx:3d}/{len(anchors)}")

    print()
    print(f"usable anchors = {len(records)}")
    print()

    # ==============================================================================================
    # BUILD REQUIRED ROOT TABLES
    # ==============================================================================================

    used_moduli = sorted(
        {
            r
            for rec in records
            for r in rec["mods"]
        }
    )

    for idx, r in enumerate(used_moduli, 1):

        sqrt_cache[r] = build_sqrt_table(r)

        print(
            f"sqrt table {idx:3d}/{len(used_moduli)} "
            f"r={r}"
        )

    print()

    # ==============================================================================================
    # GLOBAL STATISTICS
    # ==============================================================================================

    total_sum_tests = 0
    total_qr1 = 0
    total_qr2 = 0
    total_two_modulus_roots = 0
    total_small_d = 0
    total_third_modulus = 0
    total_exact = 0

    recovered = 0
    unique_recovered = 0

    prediction_distances = []

    anchor_results = []

    # ==============================================================================================
    # MAIN SEARCH
    # ==============================================================================================

    print("=" * 100)
    print("SUM -> MODULAR SQUARE ROOT -> CRT DIFFERENCE SEARCH")
    print("=" * 100)

    for anchor_index, rec in enumerate(records, 1):

        n = rec["n"]
        p = rec["p"]
        q = rec["q"]
        r1, r2, r3 = rec["mods"]

        gap = rec["gap"]

        # ==========================================================================================
        # USE THE TWO LARGEST MODULI
        # ==========================================================================================

        # r2*r3 gives the largest two-modulus CRT range.
        #
        # Since p,q are in [10k,100k],
        #
        #     |p-q| <= 90,000.
        #
        # If r2*r3 > 90,000, a CRT residue below 90,000 identifies
        # the actual difference uniquely.

        m = r2 * r3

        # Sanity.
        if m <= MAX_FACTOR_DIFF:
            raise RuntimeError(
                f"Unexpectedly small two-modulus product: "
                f"{r2}*{r3}={m}"
            )

        exact_d_candidates = set()

        sum_tests = 0
        qr1_count = 0
        qr2_count = 0
        two_root_count = 0
        small_count = 0
        third_count = 0
        exact_count = 0

        # Only even sums are possible because both factors are odd.
        #
        # Also the discriminant must be nonnegative.
        #
        # This starts with the entire factor-sum interval rather than
        # enumerating prime p.
        for s in range(SUM_MIN, SUM_MAX + 1, 2):

            D = s * s - 4 * n

            if D < 0:
                continue

            sum_tests += 1

            x2 = D % r2
            x3 = D % r3

            roots2 = sqrt_cache[r2][x2]

            if not roots2:
                continue

            qr1_count += 1

            roots3 = sqrt_cache[r3][x3]

            if not roots3:
                continue

            qr2_count += 1

            # Combine all 2x2 sign possibilities.
            for a in roots2:
                for b in roots3:

                    two_root_count += 1

                    d = crt2(
                        a,
                        r2,
                        b,
                        r3,
                    )

                    # Because actual |p-q| <= 90k and m is larger,
                    # only a small CRT representative can possibly
                    # be the real difference.
                    if d > MAX_FACTOR_DIFF:
                        continue

                    small_count += 1

                    # The third modulus is now an independent
                    # consistency check.
                    x1 = D % r1

                    if d * d % r1 != x1:
                        continue

                    third_count += 1

                    pair = recover_from_sum_difference(
                        n,
                        s,
                        d,
                    )

                    if pair is None:
                        continue

                    exact_count += 1
                    exact_d_candidates.add(
                        (
                            s,
                            d,
                            pair[0],
                            pair[1],
                        )
                    )

        total_sum_tests += sum_tests
        total_qr1 += qr1_count
        total_qr2 += qr2_count
        total_two_modulus_roots += two_root_count
        total_small_d += small_count
        total_third_modulus += third_count
        total_exact += exact_count

        if exact_d_candidates:

            recovered += 1

            if len(exact_d_candidates) == 1:
                unique_recovered += 1

            # Find nearest prediction to the actual d.
            actual_d = abs(p - q)

            best_distance = min(
                abs(item[1] - actual_d)
                for item in exact_d_candidates
            )

            prediction_distances.append(
                best_distance
            )

        anchor_results.append(
            {
                "id": rec["id"],
                "p": p,
                "q": q,
                "n": n,
                "mods": rec["mods"],
                "R": rec["R"],
                "gap": gap,
                "sum_tests": sum_tests,
                "qr1": qr1_count,
                "qr2": qr2_count,
                "two_roots": two_root_count,
                "small_d": small_count,
                "third": third_count,
                "exact": exact_count,
                "predictions": sorted(exact_d_candidates),
            }
        )

        if anchor_index % 25 == 0:
            print(
                f"anchor {anchor_index:3d}/{len(records)}"
            )

    # ==============================================================================================
    # SUMMARY
    # ==============================================================================================

    print()
    print("=" * 100)
    print("SUMMARY")
    print("=" * 100)

    print(
        f"anchors analyzed                 = "
        f"{len(records)}"
    )

    print(
        f"total even sums tested           = "
        f"{total_sum_tests:,}"
    )

    print(
        f"average sums / anchor             = "
        f"{total_sum_tests / len(records):,.2f}"
    )

    print(
        f"after QR(r2)                     = "
        f"{total_qr1:,}"
    )

    print(
        f"after QR(r3)                     = "
        f"{total_qr2:,}"
    )

    print(
        f"two-modulus CRT roots             = "
        f"{total_two_modulus_roots:,}"
    )

    print(
        f"CRT d <= {MAX_FACTOR_DIFF:,}           = "
        f"{total_small_d:,}"
    )

    print(
        f"third-modulus survivors           = "
        f"{total_third_modulus:,}"
    )

    print(
        f"exact factor solutions             = "
        f"{total_exact:,}"
    )

    print()

    # ==============================================================================================
    # REDUCTION
    # ==============================================================================================

    print("=" * 100)
    print("SEARCH REDUCTION")
    print("=" * 100)

    baseline = total_sum_tests

    def reduction(before, after):
        if before == 0:
            return 0.0

        return 100.0 * (1.0 - after / before)

    print(
        f"QR(r2) reduction                  = "
        f"{reduction(baseline, total_qr1):.6f}%"
    )

    print(
        f"QR(r3) cumulative reduction       = "
        f"{reduction(baseline, total_qr2):.6f}%"
    )

    print(
        f"CRT-small-d reduction             = "
        f"{reduction(baseline, total_small_d):.6f}%"
    )

    print(
        f"third-modulus reduction           = "
        f"{reduction(baseline, total_third_modulus):.6f}%"
    )

    print(
        f"exact-solution reduction          = "
        f"{reduction(baseline, total_exact):.6f}%"
    )

    print()

    # ==============================================================================================
    # RECOVERY
    # ==============================================================================================

    print("=" * 100)
    print("RECOVERY")
    print("=" * 100)

    print(
        f"anchors with a solution            = "
        f"{recovered}/{len(records)}"
    )

    print(
        f"anchors with exactly one solution  = "
        f"{unique_recovered}/{len(records)}"
    )

    if prediction_distances:
        print(
            f"mean |predicted d - actual d|     = "
            f"{sum(prediction_distances) / len(prediction_distances):.6f}"
        )

        print(
            f"maximum prediction distance        = "
            f"{max(prediction_distances)}"
        )

    print()

    # ==============================================================================================
    # AVERAGES
    # ==============================================================================================

    print("=" * 100)
    print("AVERAGE LOAD PER ANCHOR")
    print("=" * 100)

    N = len(records)

    print(f"sums tested               = {total_sum_tests / N:,.3f}")
    print(f"QR after first modulus    = {total_qr1 / N:,.3f}")
    print(f"QR after second modulus   = {total_qr2 / N:,.3f}")
    print(f"CRT roots                 = {total_two_modulus_roots / N:,.3f}")
    print(f"small d candidates        = {total_small_d / N:,.3f}")
    print(f"third modulus survivors   = {total_third_modulus / N:,.3f}")
    print(f"exact solutions           = {total_exact / N:,.3f}")

    print()

    # ==============================================================================================
    # ANCHOR RESULTS
    # ==============================================================================================

    print("=" * 100)
    print("ANCHOR RESULTS")
    print("=" * 100)

    print(
        " ID        p        q       R/n      "
        "SUMS     QR2     CRT<90k   3rd   EXACT"
    )

    print("-" * 100)

    for row in anchor_results:

        ratio = row["R"] / row["n"]

        print(
            f"{row['id']:3d} "
            f"{row['p']:8,d} "
            f"{row['q']:8,d} "
            f"{ratio:.7f} "
            f"{row['sum_tests']:8,d} "
            f"{row['qr2']:7,d} "
            f"{row['small_d']:9,d} "
            f"{row['third']:6,d} "
            f"{row['exact']:6,d}"
        )

    # ==============================================================================================
    # STRONGEST COLLAPSES
    # ==============================================================================================

    print()
    print("=" * 100)
    print("STRONGEST DIFFERENCE-RECONSTRUCTION COLLAPSES")
    print("=" * 100)

    strongest = sorted(
        anchor_results,
        key=lambda x: (
            x["third"],
            x["small_d"],
            x["qr2"],
        )
    )

    for row in strongest[:20]:

        r1, r2, r3 = row["mods"]

        actual_d = abs(row["p"] - row["q"])

        print(
            f"n={row['n']:,} "
            f"p={row['p']:,} "
            f"q={row['q']:,} "
            f"mods=({r1},{r2},{r3}) "
            f"R/n={row['R']/row['n']:.10f} "
            f"gap={row['gap']:,} "
            f"actual_d={actual_d:,} "
            f"QR2={row['qr2']} "
            f"CRT<90k={row['small_d']} "
            f"3rd={row['third']} "
            f"exact={row['exact']}"
        )

        if row["predictions"]:
            print(
                "    solutions:",
                row["predictions"][:10]
            )

    # ==============================================================================================
    # EXPLICIT ACTUAL-FACTOR CHECK
    # ==============================================================================================

    print()
    print("=" * 100)
    print("ACTUAL FACTOR CHECK")
    print("=" * 100)

    failures = 0

    for row in anchor_results:

        actual_s = row["p"] + row["q"]
        actual_d = abs(row["p"] - row["q"])

        found = False

        for s, d, p2, q2 in row["predictions"]:

            if (
                s == actual_s
                and d == actual_d
                and p2 * q2 == row["n"]
            ):
                found = True
                break

        if not found:
            failures += 1

    print(
        f"actual (s,d) recovered            = "
        f"{len(records) - failures}/{len(records)}"
    )

    # ==============================================================================================
    # KEY ALGEBRA
    # ==============================================================================================

    print()
    print("=" * 100)
    print("KEY ALGEBRAIC EXPLOIT")
    print("=" * 100)

    print()
    print("For:")
    print()
    print("    d = |p-q|")
    print("    s = p+q")
    print()

    print("we have:")
    print()
    print("    d² = s² - 4n")
    print()

    print("Since:")
    print()
    print("    n = R + g")
    print("    R = r1*r2*r3")
    print()

    print("and therefore:")
    print()
    print("    n ≡ g (mod r_i)")
    print()

    print("we obtain:")
    print()
    print("    d² ≡ s² - 4g (mod r_i)")
    print()

    print(
        "Thus a candidate sum produces modular square roots "
        "for the unknown factor difference."
    )

    print()
    print("Using the two largest moduli:")
    print()
    print("    d mod r2")
    print("    d mod r3")
    print()
    print("can be combined into:")
    print()
    print("    d mod (r2*r3)")
    print()

    print(
        "Because r2*r3 is larger than the entire possible "
        "factor-difference interval, a CRT result below 90,000 "
        "is potentially the actual d itself."
    )

    print()
    print(
        "This is a substantially different attack direction:"
    )

    print()
    print("    candidate sum s")
    print("          |")
    print("          v")
    print("    modular discriminants")
    print("          |")
    print("          v")
    print("    square roots mod r2,r3")
    print("          |")
    print("          v")
    print("    CRT reconstruction of d")
    print("          |")
    print("          v")
    print("    require d < 90,000")
    print("          |")
    print("          v")
    print("    third-modulus verification")
    print("          |")
    print("          v")
    print("    recover p=(s-d)/2, q=(s+d)/2")
    print()

    # ==============================================================================================
    # MOST IMPORTANT TEST
    # ==============================================================================================

    print("=" * 100)
    print("MOST IMPORTANT TEST")
    print("=" * 100)

    print()
    print(
        "The critical statistic is:"
    )

    print()
    print("    CRT<90k")
    print()

    print(
        "If this number is already very small per anchor, "
        "the modular square-root CRT is doing something stronger "
        "than the previous QR filtering experiment."
    )

    print()
    print(
        "The strongest possible result would be:"
    )

    print()
    print("    many sums tested")
    print("        -> few modular roots")
    print("        -> approximately one d < 90,000")
    print("        -> third modulus confirms d")
    print("        -> exact p,q")
    print()

    print(
        "That would provide a direct route from n and the three "
        "close moduli to the factor pair without enumerating "
        "the 8,363 factor primes."
    )

    print()
    print(
        "Conversely, if thousands of CRT-small-d states survive, "
        "then this route is essentially another representation "
        "of the same search."
    )

    # ==============================================================================================
    # TIMING
    # ==============================================================================================

    elapsed = time.perf_counter() - t0

    print()
    print("=" * 100)
    print("TIMING")
    print("=" * 100)
    print(f"total runtime = {elapsed:.3f} seconds")

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()
