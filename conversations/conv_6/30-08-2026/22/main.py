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


# ==================================================================================================
# SIEVE
# ==================================================================================================

def sieve(limit):
    a = bytearray(b"\x01") * (limit + 1)

    if limit >= 0:
        a[0] = 0
    if limit >= 1:
        a[1] = 0

    root = math.isqrt(limit)

    for p in range(2, root + 1):
        if a[p]:
            start = p * p
            count = ((limit - start) // p) + 1
            a[start::p] = b"\x00" * count

    return [i for i in range(2, limit + 1) if a[i]]


# ==================================================================================================
# CLOSE MODULUS TRIPLE
# ==================================================================================================

def select_triple(n, modulus_primes):
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
# ANCHORS
# ==================================================================================================

def generate_anchors(primes, count, rng):
    seen = set()
    out = []

    while len(out) < count:
        p, q = rng.sample(primes, 2)

        if p > q:
            p, q = q, p

        pair = (p, q)

        if pair in seen:
            continue

        seen.add(pair)
        out.append((p, q, p * q))

    return out


# ==================================================================================================
# STANDARD PRIME BASELINE
# ==================================================================================================

def prime_scan(n, mods, factor_primes):
    r1, r2, r3 = mods

    hits = []

    for p in factor_primes:
        if p % r1 == 0 or p % r2 == 0 or p % r3 == 0:
            continue

        # q modulo the three moduli.
        q1 = (n % r1) * pow(p, -1, r1) % r1
        q2 = (n % r2) * pow(p, -1, r2) % r2
        q3 = (n % r3) * pow(p, -1, r3) % r3

        # CRT.
        R = r1 * r2 * r3

        x = (
            q1 * (R // r1) * pow(R // r1, -1, r1)
            + q2 * (R // r2) * pow(R // r2, -1, r2)
            + q3 * (R // r3) * pow(R // r3, -1, r3)
        ) % R

        q = x

        if FACTOR_MIN <= q <= FACTOR_MAX:
            if p * q == n:
                hits.append((p, q))

    return hits


# ==================================================================================================
# INTERVAL FOR A QUOTIENT
# ==================================================================================================

def quotient_interval(k, r):
    """
    All integers p satisfying

        k*r <= p < (k+1)*r.

    Intersected with the global factor range.
    """

    lo = max(FACTOR_MIN, k * r)
    hi = min(FACTOR_MAX, (k + 1) * r - 1)

    if lo > hi:
        return None

    return lo, hi


# ==================================================================================================
# THREE-WAY QUOTIENT CELL SEARCH FOR ONE FACTOR
# ==================================================================================================

def build_quotient_cells(mods):
    r1, r2, r3 = mods

    cells = []

    max_k1 = FACTOR_MAX // r1

    for k1 in range(max_k1 + 1):

        i1 = quotient_interval(k1, r1)

        if i1 is None:
            continue

        lo1, hi1 = i1

        # k2 values that could overlap [lo1, hi1].
        k2_min = lo1 // r2
        k2_max = hi1 // r2

        for k2 in range(k2_min, k2_max + 1):

            i2 = quotient_interval(k2, r2)

            if i2 is None:
                continue

            lo2, hi2 = i2

            lo12 = max(lo1, lo2)
            hi12 = min(hi1, hi2)

            if lo12 > hi12:
                continue

            # k3 values that can overlap the first two.
            k3_min = lo12 // r3
            k3_max = hi12 // r3

            for k3 in range(k3_min, k3_max + 1):

                i3 = quotient_interval(k3, r3)

                if i3 is None:
                    continue

                lo3, hi3 = i3

                lo = max(lo12, lo3)
                hi = min(hi12, hi3)

                if lo <= hi:
                    cells.append(
                        (
                            lo,
                            hi,
                            k1,
                            k2,
                            k3
                        )
                    )

    return cells


# ==================================================================================================
# PAIRED QUOTIENT SEARCH
# ==================================================================================================

def paired_quotient_search(n, mods, prime_set):
    """
    Search simultaneously over quotient cells for p and q.

    Each cell represents a range of possible integers that have the
    same quotient triple:

        floor(x/r1)
        floor(x/r2)
        floor(x/r3)

    We then ask whether the product of a p-cell and q-cell can contain n.

    This is deliberately integer-based: no p list is traversed initially.
    """

    p_cells = build_quotient_cells(mods)
    q_cells = p_cells

    product_states = 0
    product_interval_hits = 0
    prime_tests = 0
    exact_hits = []

    for plo, phi, pk1, pk2, pk3 in p_cells:

        # Since p*q = n and q is in [FACTOR_MIN, FACTOR_MAX],
        # narrow q range immediately.

        q_lo_bound = max(
            FACTOR_MIN,
            (n + phi - 1) // phi
        )

        q_hi_bound = min(
            FACTOR_MAX,
            n // plo
        )

        if q_lo_bound > q_hi_bound:
            continue

        for qlo, qhi, qk1, qk2, qk3 in q_cells:

            # q cell must overlap the derived q interval.
            loq = max(qlo, q_lo_bound)
            hiq = min(qhi, q_hi_bound)

            if loq > hiq:
                continue

            product_states += 1

            # Cell product interval.
            min_product = plo * loq
            max_product = phi * hiq

            if not (min_product <= n <= max_product):
                continue

            product_interval_hits += 1

            # We deliberately search only the smaller dimension.
            #
            # Given p, q = n//p is forced.
            #
            # But we do NOT enumerate the complete prime population.
            # We enumerate integers only inside the surviving cell.

            width = min(phi - plo, hiq - loq)

            if (phi - plo) <= (hiq - loq):
                for p in range(plo, phi + 1):

                    if p not in prime_set:
                        continue

                    prime_tests += 1

                    if n % p != 0:
                        continue

                    q = n // p

                    if not (loq <= q <= hiq):
                        continue

                    if q in prime_set:
                        exact_hits.append((p, q))

            else:
                for q in range(loq, hiq + 1):

                    if q not in prime_set:
                        continue

                    prime_tests += 1

                    if n % q != 0:
                        continue

                    p = n // q

                    if not (plo <= p <= phi):
                        continue

                    if p in prime_set:
                        exact_hits.append((p, q))

    # Remove duplicates.
    exact_hits = sorted(set(exact_hits))

    return {
        "cells": len(p_cells),
        "product_states": product_states,
        "product_interval_hits": product_interval_hits,
        "prime_tests": prime_tests,
        "exact_hits": exact_hits,
    }


# ==================================================================================================
# GAP-BASED TEST
# ==================================================================================================

def gap_statistics(n, mods, p, q):
    R = math.prod(mods)
    gap = n - R

    g = int(round(R ** (1.0 / 3.0)))

    # Distances of r's from cubic scale.
    rdist = tuple(r - g for r in mods)

    # Distances of factors from sqrt scale.
    s = math.isqrt(n)

    return {
        "R": R,
        "gap": gap,
        "R_over_n": R / n,
        "cube_root_R": g,
        "sqrt_n": s,
        "p_minus_sqrt": p - s,
        "q_minus_sqrt": q - s,
        "p_plus_q": p + q,
        "p_minus_q": abs(p - q),
        "rdist": rdist,
    }


# ==================================================================================================
# MAIN
# ==================================================================================================

def run():
    t0 = time.perf_counter()

    rng = random.Random(SEED)

    print("=" * 100)
    print("THREE-CLOSE-PRIME QUOTIENT-CELL / PAIR-STATE INVERSION EXPERIMENT")
    print("=" * 100)

    print(f"M                         = {M:,}")
    print(f"factor range             = {FACTOR_MIN:,} - {FACTOR_MAX:,}")
    print(f"anchors                   = {ANCHORS}")
    print(f"modulus prime range      = {MODULUS_MIN:,} - {MODULUS_MAX:,}")
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

    prime_set = set(factor_primes)

    print(f"factor primes             = {len(factor_primes):,}")
    print(f"modulus primes            = {len(modulus_primes):,}")
    print()

    # ==============================================================================================
    # ANCHORS
    # ==============================================================================================

    anchors = generate_anchors(
        factor_primes,
        ANCHORS,
        rng
    )

    # ==============================================================================================
    # SELECT TRIPLES
    # ==============================================================================================

    print("=" * 100)
    print("SELECTING MAXIMUM-PRODUCT CLOSE PRIME TRIPLES")
    print("=" * 100)

    records = []

    for i, (p, q, n) in enumerate(anchors, 1):

        mods = select_triple(
            n,
            modulus_primes
        )

        if mods is None:
            continue

        R = math.prod(mods)

        records.append(
            {
                "id": i,
                "p": p,
                "q": q,
                "n": n,
                "mods": mods,
                "R": R,
            }
        )

        if i % 25 == 0:
            print(f"anchor {i:3d}/{len(anchors)}")

    print()
    print(f"usable anchors = {len(records)}")
    print()

    # ==============================================================================================
    # MAIN ANALYSIS
    # ==============================================================================================

    baseline_tests = 0
    quotient_cells_total = 0
    quotient_states_total = 0
    interval_hits_total = 0
    quotient_prime_tests = 0

    baseline_success = 0
    quotient_success = 0

    gap_values = []
    ratio_values = []

    result_rows = []

    print("=" * 100)
    print("RUNNING QUOTIENT-CELL ANALYSIS")
    print("=" * 100)

    for idx, rec in enumerate(records, 1):

        n = rec["n"]
        p = rec["p"]
        q = rec["q"]
        mods = rec["mods"]

        # ------------------------------------------------------------------------------------------
        # Baseline
        # ------------------------------------------------------------------------------------------

        baseline_hits = prime_scan(
            n,
            mods,
            factor_primes
        )

        baseline_tests += len(factor_primes)

        if (p, q) in baseline_hits:
            baseline_success += 1

        # ------------------------------------------------------------------------------------------
        # Quotient cell search
        # ------------------------------------------------------------------------------------------

        result = paired_quotient_search(
            n,
            mods,
            prime_set
        )

        quotient_cells_total += result["cells"]
        quotient_states_total += result["product_states"]
        interval_hits_total += result["product_interval_hits"]
        quotient_prime_tests += result["prime_tests"]

        if (p, q) in result["exact_hits"]:
            quotient_success += 1

        # ------------------------------------------------------------------------------------------
        # Gap statistics
        # ------------------------------------------------------------------------------------------

        stats = gap_statistics(
            n,
            mods,
            p,
            q
        )

        gap_values.append(stats["gap"])
        ratio_values.append(stats["R_over_n"])

        result_rows.append(
            (
                rec["id"],
                p,
                q,
                mods,
                stats["R_over_n"],
                stats["gap"],
                result["cells"],
                result["product_states"],
                result["product_interval_hits"],
                result["prime_tests"],
                len(result["exact_hits"]),
            )
        )

        if idx % 25 == 0:
            print(f"anchor {idx:3d}/{len(records)}")

    # ==============================================================================================
    # SUMMARY
    # ==============================================================================================

    N = len(records)

    print()
    print("=" * 100)
    print("SUMMARY")
    print("=" * 100)

    print(f"anchors analyzed                 = {N:,}")
    print()

    print(
        f"baseline prime tests             = "
        f"{baseline_tests:,}"
    )

    print(
        f"quotient cells                   = "
        f"{quotient_cells_total:,}"
    )

    print(
        f"paired quotient states           = "
        f"{quotient_states_total:,}"
    )

    print(
        f"product-interval survivors       = "
        f"{interval_hits_total:,}"
    )

    print(
        f"final prime tests                = "
        f"{quotient_prime_tests:,}"
    )

    print()

    print(
        f"average quotient cells / anchor  = "
        f"{quotient_cells_total / N:.3f}"
    )

    print(
        f"average paired states / anchor   = "
        f"{quotient_states_total / N:.3f}"
    )

    print(
        f"average interval hits / anchor   = "
        f"{interval_hits_total / N:.3f}"
    )

    print(
        f"average final prime tests        = "
        f"{quotient_prime_tests / N:.3f}"
    )

    print()

    # ==============================================================================================
    # REDUCTION
    # ==============================================================================================

    baseline_per_anchor = len(factor_primes)

    print("=" * 100)
    print("SEARCH REDUCTION")
    print("=" * 100)

    print(
        f"prime-pool size                  = "
        f"{baseline_per_anchor:,}"
    )

    print(
        f"mean final prime tests           = "
        f"{quotient_prime_tests / N:.3f}"
    )

    if quotient_prime_tests:
        reduction = (
            1.0
            - (quotient_prime_tests / N)
            / baseline_per_anchor
        )

        print(
            f"prime-test reduction             = "
            f"{100.0 * reduction:.6f}%"
        )

    print()

    # ==============================================================================================
    # RECOVERY
    # ==============================================================================================

    print("=" * 100)
    print("RECOVERY")
    print("=" * 100)

    print(
        f"baseline correct                 = "
        f"{baseline_success}/{N}"
    )

    print(
        f"quotient-cell correct             = "
        f"{quotient_success}/{N}"
    )

    # ==============================================================================================
    # GAP ANALYSIS
    # ==============================================================================================

    print()
    print("=" * 100)
    print("R/n AND GAP ANALYSIS")
    print("=" * 100)

    mean_gap = sum(gap_values) / N
    mean_ratio = sum(ratio_values) / N

    print(
        f"mean R/n                        = "
        f"{mean_ratio:.12f}"
    )

    print(
        f"minimum R/n                     = "
        f"{min(ratio_values):.12f}"
    )

    print(
        f"maximum R/n                     = "
        f"{max(ratio_values):.12f}"
    )

    print(
        f"mean n-R                        = "
        f"{mean_gap:,.3f}"
    )

    print(
        f"minimum n-R                     = "
        f"{min(gap_values):,}"
    )

    print(
        f"maximum n-R                     = "
        f"{max(gap_values):,}"
    )

    # ==============================================================================================
    # CORRELATION
    # ==============================================================================================

    print()
    print("=" * 100)
    print("GAP VS SEARCH-STATE CORRELATION")
    print("=" * 100)

    pairs = [
        (row[5], row[7])
        for row in result_rows
    ]

    if len(pairs) > 1:

        xmean = sum(x for x, _ in pairs) / len(pairs)
        ymean = sum(y for _, y in pairs) / len(pairs)

        num = sum(
            (x - xmean) * (y - ymean)
            for x, y in pairs
        )

        den1 = math.sqrt(
            sum((x - xmean) ** 2 for x, _ in pairs)
        )

        den2 = math.sqrt(
            sum((y - ymean) ** 2 for _, y in pairs)
        )

        if den1 and den2:
            corr = num / (den1 * den2)
        else:
            corr = 0.0

        print(
            f"Pearson correlation gap vs paired states = "
            f"{corr:.6f}"
        )

    # ==============================================================================================
    # BEST CASES
    # ==============================================================================================

    print()
    print("=" * 100)
    print("MOST AGGRESSIVE QUOTIENT-CELL COLLAPSES")
    print("=" * 100)

    best = sorted(
        result_rows,
        key=lambda x: (
            x[9],
            x[7],
        )
    )

    for row in best[:20]:

        (
            idx,
            p,
            q,
            mods,
            ratio,
            gap,
            cells,
            states,
            interval_hits,
            prime_tests,
            exact,
        ) = row

        print(
            f"n={p*q:,} "
            f"p={p:,} q={q:,} "
            f"mods={mods} "
            f"R/n={ratio:.9f} "
            f"gap={gap:,} "
            f"cells={cells:,} "
            f"states={states:,} "
            f"interval={interval_hits:,} "
            f"prime_tests={prime_tests:,} "
            f"exact={exact}"
        )

    # ==============================================================================================
    # ANCHOR TABLE
    # ==============================================================================================

    print()
    print("=" * 100)
    print("ANCHOR RESULTS")
    print("=" * 100)

    print(
        " ID        p        q      r1     r2     r3       "
        "R/n      GAP    CELLS    STATES  INTERVAL  PRIME"
    )

    print("-" * 100)

    for row in result_rows[:50]:

        (
            idx,
            p,
            q,
            mods,
            ratio,
            gap,
            cells,
            states,
            interval_hits,
            prime_tests,
            exact,
        ) = row

        r1, r2, r3 = mods

        print(
            f"{idx:3d} "
            f"{p:8,d} "
            f"{q:8,d} "
            f"{r1:6d} "
            f"{r2:6d} "
            f"{r3:6d} "
            f"{ratio:0.7f} "
            f"{gap:8,d} "
            f"{cells:7,d} "
            f"{states:8,d} "
            f"{interval_hits:8,d} "
            f"{prime_tests:7,d}"
        )

    if len(result_rows) > 50:
        print()
        print(f"... {len(result_rows) - 50} additional anchors omitted")

    # ==============================================================================================
    # FINAL INTERPRETATION
    # ==============================================================================================

    print()
    print("=" * 100)
    print("INTERPRETATION")
    print("=" * 100)

    print()
    print("Previous experiment:")
    print()
    print("    enumerate p")
    print("       |")
    print("       v")
    print("    modular inversion")
    print("       |")
    print("       v")
    print("    CRT(q)")
    print("       |")
    print("       v")
    print("    exact n = p*q")
    print()

    print("This experiment removes the initial prime enumeration.")
    print()
    print("Instead:")
    print()
    print("    p belongs to a quotient cell")
    print("    q belongs to a quotient cell")
    print("             |")
    print("             v")
    print("    test whether cell products")
    print("    can contain n")
    print("             |")
    print("             v")
    print("    enumerate only surviving integers")
    print("             |")
    print("             v")
    print("    exact divisibility")
    print()

    print("The central statistic is:")
    print()
    print("    final prime tests")
    print()
    print("compared with:")
    print()
    print("    8,363 prime tests")
    print()

    print("A strong result would mean that the quotient geometry")
    print("reduces the factor search substantially BEFORE any")
    print("prime-specific modular inversion.")
    print()

    print("The second question is whether n-R matters.")
    print()
    print("If smaller:")
    print()
    print("    n - r1*r2*r3")
    print()
    print("systematically produces fewer quotient states, then")
    print("the maximum-product condition is doing more than merely")
    print("making the CRT modulus large.")
    print()

    print("If the gap has essentially no relationship with the")
    print("search state count, then R/n is probably only measuring")
    print("the size of the modular interval rather than revealing")
    print("a deeper algebraic relation between the two factors.")
    print()

    print("Most importantly, this experiment does NOT assume")
    print("that Q can be obtained without p.")
    print("It directly tests whether that missing inversion can")
    print("be replaced by a smaller quotient-state search.")

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
