#!/usr/bin/env python3

from __future__ import annotations

import math
import random
import time
from collections import Counter
from statistics import mean


# ================================================================================================
# CONFIGURATION
# ================================================================================================

M = 111_546_435

FACTOR_MIN = 10_000
FACTOR_MAX = 100_000

MOD_MIN = 300
MOD_MAX = 3_000

ANCHORS = 300

CLOSE_RATIO = 0.20

SEED = 1_511_464_998

PROGRESS = 25


# ================================================================================================
# PRIME SIEVE
# ================================================================================================

def sieve(limit: int) -> list[int]:
    if limit < 2:
        return []

    a = bytearray(b"\x01") * (limit + 1)
    a[0:2] = b"\x00\x00"

    root = math.isqrt(limit)

    for p in range(2, root + 1):
        if not a[p]:
            continue

        start = p * p
        count = ((limit - start) // p) + 1

        a[start::p] = b"\x00" * count

    return [i for i, v in enumerate(a) if v]


def is_prime(n: int) -> bool:
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


# ================================================================================================
# OUTPUT
# ================================================================================================

def banner(text: str) -> None:
    print()
    print("=" * 100)
    print(text)
    print("=" * 100)


# ================================================================================================
# CLOSE TRIPLE
# ================================================================================================

def choose_close_triple(
    n: int,
    modulus_primes: list[int],
) -> tuple[int, int, int]:

    best = None
    best_R = -1

    for i, r1 in enumerate(modulus_primes):

        max_r = int(r1 * (1.0 + CLOSE_RATIO))

        for j in range(i + 1, len(modulus_primes)):

            r2 = modulus_primes[j]

            if r2 > max_r:
                break

            for k in range(j + 1, len(modulus_primes)):

                r3 = modulus_primes[k]

                if r3 > max_r:
                    break

                R = r1 * r2 * r3

                if R < n and R > best_R:
                    best_R = R
                    best = (r1, r2, r3)

    if best is None:
        raise RuntimeError(
            f"No usable close triple for n={n}"
        )

    return best


# ================================================================================================
# ANCHORS
# ================================================================================================

def build_anchors(
    factor_primes: list[int],
    count: int,
    rng: random.Random,
) -> list[tuple[int, int]]:

    result = []

    for _ in range(count):

        p, q = rng.sample(factor_primes, 2)

        if p > q:
            p, q = q, p

        result.append((p, q))

    return result


# ================================================================================================
# INTERVAL HELPERS
# ================================================================================================

def quotient_interval(
    r: int,
) -> tuple[int, int]:

    k_min = FACTOR_MIN // r

    if FACTOR_MIN % r != 0:
        k_min += 1

    k_max = FACTOR_MAX // r

    return k_min, k_max


def cell_bounds(
    quotient: int,
    modulus: int,
) -> tuple[int, int]:

    lo = max(
        FACTOR_MIN,
        quotient * modulus,
    )

    hi = min(
        FACTOR_MAX,
        (quotient + 1) * modulus - 1,
    )

    return lo, hi


# ================================================================================================
# MIXED QUOTIENT SEARCH
# ================================================================================================

def solve_anchor(
    n: int,
    p_true: int,
    q_true: int,
    mods: tuple[int, int, int],
) -> dict:

    r1, r2, r3 = mods

    R = r1 * r2 * r3
    gap = n - R

    # --------------------------------------------------------------------------------------------
    # We deliberately use different modulus coordinates:
    #
    #     p = k*r1 + a
    #     q = l*r2 + b
    #
    # with:
    #
    #     0 <= a < r1
    #     0 <= b < r2
    #
    # Consequently:
    #
    #     pq = r1*r2*k*l + lower-order terms.
    #
    # Since pq = r1*r2*r3 + gap, the leading quotient product k*l
    # should be related to r3.
    # --------------------------------------------------------------------------------------------

    k_min, k_max = quotient_interval(r1)
    l_min, l_max = quotient_interval(r2)

    quotient_pairs = 0
    product_cell_survivors = 0
    residue_survivors = 0
    third_modulus_survivors = 0

    prime_candidates = 0
    exact_pairs: set[tuple[int, int]] = set()

    actual_pair_survived = False

    product_error_values = []

    # Product quotient distribution.
    kl_distance_distribution = Counter()

    # --------------------------------------------------------------------------------------------
    # True quotient coordinates, only for auditing.
    # --------------------------------------------------------------------------------------------

    true_k = p_true // r1
    true_l = q_true // r2

    # --------------------------------------------------------------------------------------------
    # Search the quotient lattice.
    #
    # There are roughly:
    #
    #     (100000/r1) * (100000/r2)
    #
    # states per anchor rather than r1*(100000/r1).
    #
    # This is a completely different coordinate system from the previous a,k experiments.
    # --------------------------------------------------------------------------------------------

    for k in range(k_min, k_max + 1):

        p_lo, p_hi = cell_bounds(
            k,
            r1,
        )

        if p_lo > p_hi:
            continue

        for l in range(l_min, l_max + 1):

            quotient_pairs += 1

            q_lo, q_hi = cell_bounds(
                l,
                r2,
            )

            if q_lo > q_hi:
                continue

            # ------------------------------------------------------------------------------------
            # First filter:
            #
            # Does the PRODUCT INTERVAL of the two quotient cells contain n?
            #
            # Because p,q are positive:
            #
            #     [p_lo,p_hi] * [q_lo,q_hi]
            #
            # is simply:
            #
            #     [p_lo*q_lo, p_hi*q_hi].
            # ------------------------------------------------------------------------------------

            min_product = p_lo * q_lo
            max_product = p_hi * q_hi

            if n < min_product or n > max_product:
                continue

            product_cell_survivors += 1

            # Compare KL directly with the third modulus.
            #
            # This measures whether:
            #
            #       k*l ~= r3
            #
            # contains useful information.
            #
            # We use several natural correction terms as diagnostics.
            kl = k * l

            kl_distance = abs(
                kl - r3
            )

            kl_distance_distribution[
                kl_distance
            ] += 1

            # Exact product-cell coordinates of a genuine solution must satisfy this identity.
            #
            # We don't require it yet because residues within a cell matter.

            # ------------------------------------------------------------------------------------
            # Enumerate residue coordinates only AFTER the quotient cell survives.
            #
            # This is the main search reduction mechanism.
            # ------------------------------------------------------------------------------------

            p_res_lo = p_lo - k * r1
            p_res_hi = p_hi - k * r1

            q_res_lo = q_lo - l * r2
            q_res_hi = q_hi - l * r2

            # Small correction ranges.
            for a in range(
                p_res_lo,
                p_res_hi + 1,
            ):

                if not (
                    0 <= a < r1
                ):
                    continue

                # Solve the exact product equation for b:
                #
                # (k*r1+a)(l*r2+b)=n
                #
                # => b = n/(k*r1+a) - l*r2.
                #
                # This is the crucial step: b is obtained by division,
                # not enumerated.
                #
                p = k * r1 + a

                if p <= 0:
                    continue

                if n % p != 0:
                    continue

                q = n // p

                if not (
                    q_lo <= q <= q_hi
                ):
                    continue

                b = q - l * r2

                if not (
                    0 <= b < r2
                ):
                    continue

                residue_survivors += 1

                # --------------------------------------------------------------------------------
                # Third modulus consistency.
                #
                # No inversion is used here. We directly test:
                #
                #     pq == n (mod r3)
                #
                # and also the residue representation induced by the
                # mixed quotient coordinates.
                # --------------------------------------------------------------------------------

                if (
                    p * q
                ) % r3 != n % r3:
                    continue

                third_modulus_survivors += 1

                pair = tuple(
                    sorted((p, q))
                )

                exact_pairs.add(pair)

                if (
                    is_prime(p)
                    and is_prime(q)
                ):
                    prime_candidates += 1

                    if pair == tuple(
                        sorted((p_true, q_true))
                    ):
                        actual_pair_survived = True

    # --------------------------------------------------------------------------------------------
    # Distance of the true quotient product from r3.
    # --------------------------------------------------------------------------------------------

    true_kl = true_k * true_l

    true_kl_distance = abs(
        true_kl - r3
    )

    # --------------------------------------------------------------------------------------------
    # Count states lying in progressively tighter KL shells.
    # --------------------------------------------------------------------------------------------

    shell_1 = 0
    shell_2 = 0
    shell_5 = 0
    shell_10 = 0
    shell_25 = 0
    shell_50 = 0

    for d, count in kl_distance_distribution.items():

        if d <= 1:
            shell_1 += count

        if d <= 2:
            shell_2 += count

        if d <= 5:
            shell_5 += count

        if d <= 10:
            shell_10 += count

        if d <= 25:
            shell_25 += count

        if d <= 50:
            shell_50 += count

    return {
        "r1": r1,
        "r2": r2,
        "r3": r3,
        "R": R,
        "gap": gap,
        "R_over_n": R / n,

        "k_min": k_min,
        "k_max": k_max,
        "l_min": l_min,
        "l_max": l_max,

        "quotient_pairs": quotient_pairs,
        "product_cell_survivors": product_cell_survivors,
        "residue_survivors": residue_survivors,
        "third_modulus_survivors": third_modulus_survivors,

        "prime_candidates": prime_candidates,
        "exact_pairs": sorted(exact_pairs),

        "actual_pair_survived": actual_pair_survived,

        "true_k": true_k,
        "true_l": true_l,
        "true_kl": true_kl,
        "true_kl_distance": true_kl_distance,

        "shell_1": shell_1,
        "shell_2": shell_2,
        "shell_5": shell_5,
        "shell_10": shell_10,
        "shell_25": shell_25,
        "shell_50": shell_50,
    }


# ================================================================================================
# MAIN
# ================================================================================================

def run():

    rng = random.Random(SEED)

    banner(
        "THREE-CLOSE-PRIME MIXED QUOTIENT-PRODUCT / "
        "THIRD-MODULUS EXPERIMENT"
    )

    print(
        f"M                         = {M:,}"
    )

    print(
        f"factor range              = "
        f"{FACTOR_MIN:,} - {FACTOR_MAX:,}"
    )

    print(
        f"anchors                   = {ANCHORS}"
    )

    print(
        f"modulus prime range       = "
        f"{MOD_MIN:,} - {MOD_MAX:,}"
    )

    print(
        f"close ratio               = "
        f"{CLOSE_RATIO:.0%}"
    )

    print(
        f"seed                      = "
        f"{SEED:,}"
    )

    # ============================================================================================
    # PRIME POOLS
    # ============================================================================================

    banner("BUILDING PRIME POOLS")

    all_primes = sieve(
        FACTOR_MAX
    )

    factor_primes = [
        p
        for p in all_primes
        if FACTOR_MIN <= p <= FACTOR_MAX
    ]

    modulus_primes = [
        p
        for p in all_primes
        if MOD_MIN <= p <= MOD_MAX
    ]

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

    anchors = build_anchors(
        factor_primes,
        ANCHORS,
        rng,
    )

    banner(
        "SELECTING MAXIMUM-PRODUCT CLOSE TRIPLES"
    )

    triples = []

    for i, (p, q) in enumerate(
        anchors,
        1,
    ):

        n = p * q

        mods = choose_close_triple(
            n,
            modulus_primes,
        )

        triples.append(
            (n, p, q, mods)
        )

        if i % PROGRESS == 0:
            print(
                f"anchor {i:3d}/{ANCHORS}"
            )

    print()
    print(
        f"usable anchors            = "
        f"{len(triples)}"
    )

    # ============================================================================================
    # ALGEBRA AUDIT
    # ============================================================================================

    banner(
        "VALIDATING MIXED QUOTIENT COORDINATES"
    )

    failures = 0

    for i, (
        n,
        p,
        q,
        mods,
    ) in enumerate(
        triples,
        1,
    ):

        r1, r2, r3 = mods

        a = p % r1
        b = q % r2

        k = p // r1
        l = q // r2

        if not (
            p == k * r1 + a
        ):
            failures += 1
            print(
                f"anchor {i}: p coordinate failure"
            )
            break

        if not (
            q == l * r2 + b
        ):
            failures += 1
            print(
                f"anchor {i}: q coordinate failure"
            )
            break

        if not (
            p * q == n
        ):
            failures += 1
            print(
                f"anchor {i}: product failure"
            )
            break

        if not (
            p * q % r3 == n % r3
        ):
            failures += 1
            print(
                f"anchor {i}: r3 failure"
            )
            break

        if i % PROGRESS == 0:
            print(
                f"anchor {i:3d}/{ANCHORS}"
            )

    print()
    print(
        f"coordinate failures        = "
        f"{failures}"
    )

    if failures:
        raise RuntimeError(
            "Mixed quotient coordinate validation failed."
        )

    # ============================================================================================
    # SEARCH
    # ============================================================================================

    banner(
        "RUNNING MIXED QUOTIENT-PRODUCT SEARCH"
    )

    start = time.perf_counter()

    results = []

    for i, (
        n,
        p,
        q,
        mods,
    ) in enumerate(
        triples,
        1,
    ):

        result = solve_anchor(
            n,
            p,
            q,
            mods,
        )

        results.append(
            (n, p, q, result)
        )

        if i % PROGRESS == 0:
            print(
                f"anchor {i:3d}/{ANCHORS}"
            )

    elapsed = (
        time.perf_counter()
        - start
    )

    # ============================================================================================
    # SUMMARY
    # ============================================================================================

    banner("SUMMARY")

    print(
        f"anchors analyzed              = "
        f"{len(results)}"
    )

    avg_pairs = mean(
        r["quotient_pairs"]
        for _, _, _, r in results
    )

    avg_cells = mean(
        r["product_cell_survivors"]
        for _, _, _, r in results
    )

    avg_residues = mean(
        r["residue_survivors"]
        for _, _, _, r in results
    )

    avg_r3 = mean(
        r["third_modulus_survivors"]
        for _, _, _, r in results
    )

    avg_prime = mean(
        r["prime_candidates"]
        for _, _, _, r in results
    )

    avg_exact = mean(
        len(r["exact_pairs"])
        for _, _, _, r in results
    )

    print(
        f"average quotient pairs        = "
        f"{avg_pairs:.3f}"
    )

    print(
        f"average product-cell survivors= "
        f"{avg_cells:.3f}"
    )

    print(
        f"average residue survivors     = "
        f"{avg_residues:.3f}"
    )

    print(
        f"average r3 survivors          = "
        f"{avg_r3:.3f}"
    )

    print(
        f"average prime candidates      = "
        f"{avg_prime:.3f}"
    )

    print(
        f"average exact pairs           = "
        f"{avg_exact:.3f}"
    )

    # ============================================================================================
    # BASELINE
    # ============================================================================================

    banner("SEARCH REDUCTION")

    baseline = len(
        factor_primes
    )

    print(
        f"ordinary factor-prime tests    = "
        f"{baseline:,}"
    )

    print(
        f"quotient-pair / baseline       = "
        f"{avg_pairs / baseline:.9f}"
    )

    print(
        f"product-cell / baseline       = "
        f"{avg_cells / baseline:.9f}"
    )

    print(
        f"residue / baseline             = "
        f"{avg_residues / baseline:.9f}"
    )

    print(
        f"r3 / baseline                  = "
        f"{avg_r3 / baseline:.9f}"
    )

    # ============================================================================================
    # REDUCTION PERCENTAGES
    # ============================================================================================

    banner("REDUCTION PERCENTAGES")

    print(
        f"quotient-pair reduction        = "
        f"{100.0 * (1.0 - avg_pairs / baseline):.6f}%"
    )

    print(
        f"product-cell reduction         = "
        f"{100.0 * (1.0 - avg_cells / baseline):.6f}%"
    )

    print(
        f"residue reduction              = "
        f"{100.0 * (1.0 - avg_residues / baseline):.6f}%"
    )

    print(
        f"r3 reduction                   = "
        f"{100.0 * (1.0 - avg_r3 / baseline):.6f}%"
    )

    # ============================================================================================
    # RECOVERY
    # ============================================================================================

    banner("RECOVERY")

    recovered = sum(
        1
        for _, _, _, r in results
        if r["actual_pair_survived"]
    )

    print(
        f"correctly recovered             = "
        f"{recovered}/{len(results)}"
    )

    print(
        f"recovery rate                   = "
        f"{100.0 * recovered / len(results):.4f}%"
    )

    # ============================================================================================
    # QUOTIENT PRODUCT
    # ============================================================================================

    banner(
        "QUOTIENT PRODUCT RELATION"
    )

    mean_distance = mean(
        r["true_kl_distance"]
        for _, _, _, r in results
    )

    median_distance = sorted(
        r["true_kl_distance"]
        for _, _, _, r in results
    )[len(results) // 2]

    print(
        f"mean |k*l-r3|                  = "
        f"{mean_distance:.3f}"
    )

    print(
        f"median |k*l-r3|               = "
        f"{median_distance}"
    )

    print()

    print(
        "true quotient relation:"
    )

    print(
        "    k = floor(p/r1)"
    )

    print(
        "    l = floor(q/r2)"
    )

    print(
        "    k*l compared against r3"
    )

    # ============================================================================================
    # SHELL DISTRIBUTION
    # ============================================================================================

    banner(
        "PRODUCT-SHELL SURVIVORS"
    )

    shell_counts = {
        1: sum(
            1
            for _, _, _, r in results
            if r["shell_1"] > 0
        ),
        2: sum(
            1
            for _, _, _, r in results
            if r["shell_2"] > 0
        ),
        5: sum(
            1
            for _, _, _, r in results
            if r["shell_5"] > 0
        ),
        10: sum(
            1
            for _, _, _, r in results
            if r["shell_10"] > 0
        ),
        25: sum(
            1
            for _, _, _, r in results
            if r["shell_25"] > 0
        ),
        50: sum(
            1
            for _, _, _, r in results
            if r["shell_50"] > 0
        ),
    }

    for radius, count in sorted(
        shell_counts.items()
    ):
        print(
            f"anchors with quotient pair "
            f"|k*l-r3| <= {radius:2d} = "
            f"{count:3d}/{len(results)}"
        )

    # ============================================================================================
    # ANCHOR EXAMPLES
    # ============================================================================================

    banner("ANCHOR EXAMPLES")

    for n, p, q, r in results[:20]:

        print(
            f"n={n:,} "
            f"p={p:,} "
            f"q={q:,} "
            f"mods=({r['r1']},{r['r2']},{r['r3']}) "
            f"R/n={r['R_over_n']:.10f} "
            f"gap={r['gap']:,}"
        )

        print(
            f"    k={r['true_k']} "
            f"l={r['true_l']} "
            f"k*l={r['true_kl']} "
            f"|k*l-r3|={r['true_kl_distance']}"
        )

        print(
            f"    quotient_pairs="
            f"{r['quotient_pairs']:,} "
            f"cells="
            f"{r['product_cell_survivors']:,} "
            f"residue="
            f"{r['residue_survivors']:,} "
            f"r3="
            f"{r['third_modulus_survivors']:,} "
            f"exact="
            f"{len(r['exact_pairs'])}"
        )

    # ============================================================================================
    # STRONGEST COLLAPSES
    # ============================================================================================

    banner(
        "STRONGEST QUOTIENT-PRODUCT COLLAPSES"
    )

    ranked = sorted(
        results,
        key=lambda x: (
            x[3]["product_cell_survivors"],
            x[3]["residue_survivors"],
            x[3]["third_modulus_survivors"],
        ),
    )

    for n, p, q, r in ranked[:20]:

        print(
            f"n={n:,} "
            f"p={p:,} "
            f"q={q:,} "
            f"mods=({r['r1']},{r['r2']},{r['r3']}) "
            f"R/n={r['R_over_n']:.10f} "
            f"cells={r['product_cell_survivors']} "
            f"residue={r['residue_survivors']} "
            f"r3={r['third_modulus_survivors']} "
            f"exact={len(r['exact_pairs'])}"
        )

    # ============================================================================================
    # GAP CORRELATION
    # ============================================================================================

    banner(
        "GAP / QUOTIENT-PRODUCT RELATION"
    )

    gaps = [
        r["gap"]
        for _, _, _, r in results
    ]

    cells = [
        r["product_cell_survivors"]
        for _, _, _, r in results
    ]

    def pearson(xs, ys):

        mx = mean(xs)
        my = mean(ys)

        num = sum(
            (x - mx) * (y - my)
            for x, y in zip(xs, ys)
        )

        den_x = math.sqrt(
            sum(
                (x - mx) ** 2
                for x in xs
            )
        )

        den_y = math.sqrt(
            sum(
                (y - my) ** 2
                for y in ys
            )
        )

        if den_x == 0 or den_y == 0:
            return 0.0

        return num / (den_x * den_y)

    print(
        f"corr(gap, product-cell survivors) = "
        f"{pearson(gaps, cells):.6f}"
    )

    true_dist = [
        r["true_kl_distance"]
        for _, _, _, r in results
    ]

    print(
        f"corr(gap, |k*l-r3|)               = "
        f"{pearson(gaps, true_dist):.6f}"
    )

    # ============================================================================================
    # TIMING
    # ============================================================================================

    banner("TIMING")

    print(
        f"search runtime               = "
        f"{elapsed:.3f} seconds"
    )

    # ============================================================================================
    # INTERPRETATION
    # ============================================================================================

    banner("MATHEMATICAL INTERPRETATION")

    print(
        r"""
This experiment changes coordinates again.

Instead of:

    p = a + k*r1
    q = b + l*r1

we use two different modulus scales:

    p = a + k*r1
    q = b + l*r2

where:

    0 <= a < r1
    0 <= b < r2.

The exact product is:

    n = (a+k*r1)(b+l*r2).

Expanding:

    n =
        a*b
        + a*l*r2
        + b*k*r1
        + k*l*r1*r2.

But:

    n = r1*r2*r3 + g.

Therefore the dominant quotient term is:

    k*l*r1*r2

while the target is:

    r3*r1*r2.

This suggests:

    k*l approximately r3

with the deviation caused by the residue and cross terms.

The experiment therefore asks whether:

    |k*l-r3|

is systematically small enough to collapse the quotient lattice.

The search order is:

    quotient pair (k,l)
             |
             v
    product-cell contains n
             |
             v
    solve residue b exactly
             |
             v
    third modulus consistency
             |
             v
    prime/exact verification.

This is deliberately different from the earlier (a,k)
linearization because the first enumeration contains only quotient
indices, not every possible residue a.

The critical statistics are:

    average quotient pairs
    average product-cell survivors
    average residue survivors
    average r3 survivors

and especially:

    |k*l-r3|.

If k*l stays close to r3 for the true factor pair while the number
of quotient pairs is dramatically smaller than the 8,363-prime
baseline, this gives us a new low-dimensional coordinate system.

If |k*l-r3| behaves broadly and the product-cell filter is weak,
then the mixed quotient representation is only another reparameter-
ization of the same factor search.

Every surviving candidate is checked against the exact equation:

    p*q == n

so a modular false positive cannot count as a successful factor.
"""
    )

    banner("EXPERIMENT COMPLETE")


if __name__ == "__main__":
    run()
