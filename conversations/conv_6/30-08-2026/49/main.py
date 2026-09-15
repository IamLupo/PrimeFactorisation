#!/usr/bin/env python3

"""
============================================================================================
TWO-MODULUS E/G COUPLED CARRY-INVERSION EXPERIMENT
============================================================================================

Research target
---------------

Previous experiment:

    E = c1 + c2 + c3

provided a real but incomplete pruning mechanism.

We now combine E with:

    g = n mod (r1*r2)

and test whether the remainder g can constrain the carry rectangle
BEFORE enumerating residue states.

Definitions:

    p = a + k*r1
    q = b + l*r2

    R = r1*r2

    T = floor(n/R)

    E = T - k*l

    a*l = c1*r1 + d1
    b*k = c2*r2 + d2

    E = c1 + c2 + c3

and:

    g = n mod R.

From:

    n = (a + k*r1)(b + l*r2)

we obtain:

    a(b + l*r2) + b*k*r1 = g + E*R

Therefore:

    a =
        (g + E*R - b*k*r1)
        / (b + l*r2)

This gives an exact rational function a(b).

For fixed:

    (E, k, l, c1, c2, c3)

the carry constraints give:

    a_lo <= a <= a_hi
    b_lo <= b <= b_hi.

Since a(b) is strictly decreasing, the a interval can be
inverted into a MUCH smaller b interval before any exact
divisibility testing.

The experiment measures:

    carry cells
        ->
    b states before g inversion
        ->
    b states after g inversion
        ->
    exact division hits.

The critical question:

    Does (E,g) collapse the carry rectangles strongly
    without scanning the full residue rectangle?

No true p or q is used for candidate generation.

The true factors are used only for validation/statistics.

============================================================================================
"""

import math
import random
import statistics
import time
from collections import Counter


# ==========================================================================================
# CONFIGURATION
# ==========================================================================================

N_ANCHORS = 100

P_MIN = 10_000
P_MAX = 100_000

MOD_MIN = 300
MOD_MAX = 3_000

CLOSE_RATIO = 0.20

SEED = 1_511_464_998

PROGRESS_EVERY = 10

# To keep the experiment deliberately manageable.
# True factor range remains 10k..100k.
MAX_B_SCAN_PER_CELL = None


# ==========================================================================================
# UTILITIES
# ==========================================================================================

def sieve(limit):
    is_prime = bytearray(b"\x01") * (limit + 1)
    is_prime[0:2] = b"\x00\x00"

    for p in range(2, int(limit ** 0.5) + 1):
        if is_prime[p]:
            start = p * p
            is_prime[start::p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    return [i for i in range(2, limit + 1) if is_prime[i]]


def ceil_div(a, b):
    if b <= 0:
        raise ValueError("ceil_div requires positive denominator")
    return -((-a) // b)


def floor_div(a, b):
    if b <= 0:
        raise ValueError("floor_div requires positive denominator")
    return a // b


def factor_integer(n):
    """
    Exact factorization for the small K values encountered here.

    Returns:
        [(prime, exponent), ...]
    """
    if n < 1:
        return []

    out = []

    d = 2

    while d * d <= n:
        if n % d == 0:
            e = 0
            while n % d == 0:
                n //= d
                e += 1
            out.append((d, e))

        d = 3 if d == 2 else d + 2

    if n > 1:
        out.append((n, 1))

    return out


def divisor_pairs(n):
    """
    Return all positive integer (k,l) pairs satisfying k*l=n.
    """
    pairs = []

    if n <= 0:
        return pairs

    d = 1

    while d * d <= n:
        if n % d == 0:
            q = n // d
            pairs.append((d, q))

            if d != q:
                pairs.append((q, d))

        d += 1

    return pairs


# ==========================================================================================
# PRIME POOLS
# ==========================================================================================

def build_factor_primes():
    return sieve(P_MAX)


def build_modulus_primes():
    return [
        p for p in sieve(MOD_MAX)
        if MOD_MIN <= p <= MOD_MAX
    ]


# ==========================================================================================
# ANCHORS
# ==========================================================================================

def build_anchors(factor_primes, count, seed):
    rng = random.Random(seed)

    usable = [
        p for p in factor_primes
        if P_MIN <= p <= P_MAX
    ]

    anchors = []

    seen = set()

    while len(anchors) < count:
        p, q = rng.sample(usable, 2)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        key = (p, q)

        if key in seen:
            continue

        seen.add(key)

        n = p * q

        anchors.append({
            "p": p,
            "q": q,
            "n": n,
        })

    return anchors


# ==========================================================================================
# CLOSE MODULUS PAIRS
# ==========================================================================================

def build_close_pairs(modulus_primes):
    pairs = []

    for i in range(len(modulus_primes)):
        r1 = modulus_primes[i]

        for j in range(i + 1, len(modulus_primes)):
            r2 = modulus_primes[j]

            ratio = abs(r2 - r1) / min(r1, r2)

            if ratio <= CLOSE_RATIO:
                pairs.append((r1, r2))

    return pairs


def select_pair_for_anchor(close_pairs, n, rng):
    """
    Prefer large R=r1*r2 while keeping R<n.

    This reproduces the regime used in the previous experiments:
    moduli are close and near the upper end of the allowed range.
    """
    valid = [
        pair
        for pair in close_pairs
        if pair[0] * pair[1] < n
    ]

    if not valid:
        raise RuntimeError("No valid close modulus pair for anchor.")

    valid.sort(key=lambda x: x[0] * x[1], reverse=True)

    # Pick randomly from the top few percent so the triples/pairs
    # are not identical across all anchors.
    top_count = max(1, min(len(valid), max(10, len(valid) // 20)))

    return valid[rng.randrange(top_count)]


# ==========================================================================================
# FACTOR-RANGE COORDINATE BOUNDS
# ==========================================================================================

def coordinate_bounds(base, modulus, lo, hi):
    """
    For x = residue + k*modulus with:
        lo <= x <= hi

    return possible k range.

    Actual residue bounds are handled separately.
    """
    k_min = ceil_div(lo - (modulus - 1), modulus)
    k_max = floor_div(hi, modulus)

    k_min = max(k_min, 0)

    return k_min, k_max


def residue_interval_for_cell(k, modulus, lo, hi):
    """
    Given:

        x = a + k*modulus
        0 <= a < modulus

    and:

        lo <= x <= hi

    return possible a interval.
    """
    a_lo = max(0, lo - k * modulus)
    a_hi = min(modulus - 1, hi - k * modulus)

    if a_lo > a_hi:
        return None

    return a_lo, a_hi


# ==========================================================================================
# TRUE E/CARRY DECOMPOSITION
# ==========================================================================================

def true_decomposition(p, q, r1, r2):
    R = r1 * r2

    a = p % r1
    b = q % r2

    k = (p - a) // r1
    l = (q - b) // r2

    T = p * q // R

    E = T - k * l

    z1 = a * l
    z2 = b * k

    c1 = z1 // r1
    d1 = z1 % r1

    c2 = z2 // r2
    d2 = z2 % r2

    H = d1 * r2 + d2 * r1 + a * b

    c3 = H // R
    g = H % R

    assert E == c1 + c2 + c3

    assert g == (p * q) % R

    return {
        "a": a,
        "b": b,
        "k": k,
        "l": l,
        "T": T,
        "E": E,
        "c1": c1,
        "c2": c2,
        "c3": c3,
        "d1": d1,
        "d2": d2,
        "g": g,
    }


# ==========================================================================================
# CARRY CELL GENERATION
# ==========================================================================================

def carry_states_for_pair(E, k, l, r1, r2):
    """
    Generate all possible (c1,c2,c3) states satisfying:

        E = c1+c2+c3

    with:

        0 <= c1 < l
        0 <= c2 < k
        c3 in {0,1,2}
    """
    states = []

    for c3 in (0, 1, 2):
        remaining = E - c3

        if remaining < 0:
            continue

        c1_lo = max(0, remaining - (k - 1))
        c1_hi = min(l - 1, remaining)

        for c1 in range(c1_lo, c1_hi + 1):
            c2 = remaining - c1

            if not (0 <= c2 < k):
                continue

            # Convert carries to residue intervals.

            # c1 = floor(a*l/r1)
            a_lo = ceil_div(c1 * r1, l)
            a_hi = floor_div((c1 + 1) * r1 - 1, l)

            # c2 = floor(b*k/r2)
            b_lo = ceil_div(c2 * r2, k)
            b_hi = floor_div((c2 + 1) * r2 - 1, k)

            a_lo = max(a_lo, 0)
            a_hi = min(a_hi, r1 - 1)

            b_lo = max(b_lo, 0)
            b_hi = min(b_hi, r2 - 1)

            if a_lo > a_hi or b_lo > b_hi:
                continue

            states.append({
                "c1": c1,
                "c2": c2,
                "c3": c3,
                "a_lo": a_lo,
                "a_hi": a_hi,
                "b_lo": b_lo,
                "b_hi": b_hi,
            })

    return states


# ==========================================================================================
# G-COUPLED INTERVAL INVERSION
# ==========================================================================================

def invert_g_constraint(
    n,
    g,
    E,
    k,
    l,
    r1,
    r2,
    a_lo,
    a_hi,
    b_lo,
    b_hi,
):
    """
    We have:

        a = (g + E*R - b*k*r1) / (b + l*r2)

    Define:

        C = g + E*R
        K = k*r1
        L = l*r2

    Then:

        a(b) = (C - bK)/(b+L).

    This function is strictly decreasing in b.

    We therefore invert:

        a <= a_hi

    and

        a >= a_lo

    into an integer b interval.

    No exhaustive b scan is performed here.

    Returns:
        (new_b_lo, new_b_hi) or None
    """
    R = r1 * r2

    C = g + E * R

    K = k * r1
    L = l * r2

    # From:
    #
    #     a >= a_lo
    #
    #     C - bK >= a_lo(b+L)
    #
    #     C - a_lo*L >= b(K+a_lo)
    #
    # therefore:
    #
    #     b <= (C-a_lo*L)/(K+a_lo)

    upper = floor_div(
        C - a_lo * L,
        K + a_lo,
    )

    # From:
    #
    #     a <= a_hi
    #
    #     C - bK <= a_hi(b+L)
    #
    #     C - a_hi*L <= b(K+a_hi)
    #
    # therefore:
    #
    #     b >= (C-a_hi*L)/(K+a_hi)

    lower = ceil_div(
        C - a_hi * L,
        K + a_hi,
    )

    new_lo = max(b_lo, lower)
    new_hi = min(b_hi, upper)

    if new_lo > new_hi:
        return None

    return new_lo, new_hi


# ==========================================================================================
# EXACT G-COUPLED SCAN
# ==========================================================================================

def scan_g_interval(
    g,
    E,
    k,
    l,
    r1,
    r2,
    a_lo,
    a_hi,
    b_lo,
    b_hi,
):
    """
    After interval inversion, scan only the remaining b values.

    For each b:

        q = b + l*r2

    and:

        numerator = g + E*R - b*k*r1

    We require:

        numerator % q == 0

    then:

        a = numerator / q

    Every candidate is range checked and finally checked through
    p*q == n by the caller.
    """
    R = r1 * r2
    C = g + E * R

    K = k * r1
    L = l * r2

    hits = []

    count = 0

    for b in range(b_lo, b_hi + 1):
        count += 1

        q = b + L

        numerator = C - b * K

        if numerator <= 0:
            continue

        if numerator % q != 0:
            continue

        a = numerator // q

        if not (a_lo <= a <= a_hi):
            continue

        p = a + K

        if p <= 0:
            continue

        hits.append((p, q, a, b))

    return count, hits


# ==========================================================================================
# ANCHOR SEARCH
# ==========================================================================================

def run_anchor(anchor, r1, r2):
    p_true = anchor["p"]
    q_true = anchor["q"]
    n = anchor["n"]

    R = r1 * r2

    T = n // R

    g = n % R

    decomposition = true_decomposition(
        p_true,
        q_true,
        r1,
        r2,
    )

    E_true = decomposition["E"]

    # Generic observable E upper bound.

    # For any valid factor coordinates:
    #
    #     E < k+l+1
    #
    # and k/l are bounded by the factor search interval.

    k_max = max(0, P_MAX // r1)
    l_max = max(0, P_MAX // r2)

    E_max = k_max + l_max + 1

    total_quotient_pairs = 0

    total_carry_cells = 0

    total_b_before = 0
    total_b_after = 0

    total_g_feasible_cells = 0
    total_division_hits = 0

    exact_solutions = []

    E_values_tested = 0

    true_E_retained = False
    true_carry_retained = False
    true_g_cell_retained = False

    per_cell_records = []

    for E in range(E_max + 1):

        Kprod = T - E

        if Kprod <= 0:
            continue

        E_values_tested += 1

        pair_candidates = divisor_pairs(Kprod)

        for k, l in pair_candidates:

            # Determine whether this quotient cell can intersect
            # the factor range.

            p_cell_lo = k * r1
            p_cell_hi = (k + 1) * r1 - 1

            q_cell_lo = l * r2
            q_cell_hi = (l + 1) * r2 - 1

            if p_cell_hi < P_MIN or p_cell_lo > P_MAX:
                continue

            if q_cell_hi < P_MIN or q_cell_lo > P_MAX:
                continue

            p_range_a = residue_interval_for_cell(
                k,
                r1,
                P_MIN,
                P_MAX,
            )

            q_range_b = residue_interval_for_cell(
                l,
                r2,
                P_MIN,
                P_MAX,
            )

            if p_range_a is None or q_range_b is None:
                continue

            pair_states = carry_states_for_pair(
                E,
                k,
                l,
                r1,
                r2,
            )

            total_quotient_pairs += 1

            for state in pair_states:

                total_carry_cells += 1

                # Intersect carry-derived a/b intervals with
                # the global factor interval.

                a_lo = max(
                    state["a_lo"],
                    p_range_a[0],
                )

                a_hi = min(
                    state["a_hi"],
                    p_range_a[1],
                )

                b_lo = max(
                    state["b_lo"],
                    q_range_b[0],
                )

                b_hi = min(
                    state["b_hi"],
                    q_range_b[1],
                )

                if a_lo > a_hi or b_lo > b_hi:
                    continue

                b_before = b_hi - b_lo + 1

                total_b_before += b_before

                # ------------------------------------------------------------------
                # THE NEW PART:
                #
                # Use g together with E BEFORE scanning b.
                # ------------------------------------------------------------------

                narrowed = invert_g_constraint(
                    n,
                    g,
                    E,
                    k,
                    l,
                    r1,
                    r2,
                    a_lo,
                    a_hi,
                    b_lo,
                    b_hi,
                )

                if narrowed is None:
                    continue

                total_g_feasible_cells += 1

                nb_lo, nb_hi = narrowed

                b_after = nb_hi - nb_lo + 1

                total_b_after += b_after

                scan_count, hits = scan_g_interval(
                    g,
                    E,
                    k,
                    l,
                    r1,
                    r2,
                    a_lo,
                    a_hi,
                    nb_lo,
                    nb_hi,
                )

                total_division_hits += len(hits)

                for p, q, a, b in hits:

                    if p * q != n:
                        continue

                    exact_solutions.append((p, q))

                per_cell_records.append({
                    "E": E,
                    "k": k,
                    "l": l,
                    "c1": state["c1"],
                    "c2": state["c2"],
                    "c3": state["c3"],
                    "a_lo": a_lo,
                    "a_hi": a_hi,
                    "b_lo": b_lo,
                    "b_hi": b_hi,
                    "b_before": b_before,
                    "b_after": b_after,
                    "hits": len(hits),
                })

                # True-path tracking.

                if (
                    E == E_true
                    and k == decomposition["k"]
                    and l == decomposition["l"]
                ):
                    true_carry_retained = True

                    if (
                        state["c1"] == decomposition["c1"]
                        and state["c2"] == decomposition["c2"]
                        and state["c3"] == decomposition["c3"]
                    ):
                        true_E_retained = True

                        if (
                            decomposition["b"] >= nb_lo
                            and decomposition["b"] <= nb_hi
                        ):
                            true_g_cell_retained = True

    exact_solutions = sorted(set(exact_solutions))

    recovered = (
        (p_true, q_true) in exact_solutions
        or (q_true, p_true) in exact_solutions
    )

    if E_true <= E_max:
        true_E_retained = True

    return {
        "p": p_true,
        "q": q_true,
        "n": n,
        "r1": r1,
        "r2": r2,
        "R": R,
        "T": T,
        "g": g,
        "E_true": E_true,
        "E_max": E_max,
        "true_k": decomposition["k"],
        "true_l": decomposition["l"],
        "true_c1": decomposition["c1"],
        "true_c2": decomposition["c2"],
        "true_c3": decomposition["c3"],
        "quotient_pairs": total_quotient_pairs,
        "carry_cells": total_carry_cells,
        "g_cells": total_g_feasible_cells,
        "b_before": total_b_before,
        "b_after": total_b_after,
        "division_hits": total_division_hits,
        "exact_solutions": len(exact_solutions),
        "recovered": recovered,
        "true_E_retained": true_E_retained,
        "true_carry_retained": true_carry_retained,
        "true_g_cell_retained": true_g_cell_retained,
        "records": per_cell_records,
    }


# ==========================================================================================
# MAIN
# ==========================================================================================

def run():
    rng = random.Random(SEED)

    start_total = time.perf_counter()

    print("=" * 92)
    print("TWO-MODULUS E/G COUPLED CARRY-INVERSION EXPERIMENT")
    print("=" * 92)
    print(f"N anchors                 = {N_ANCHORS:,}")
    print(f"factor range              = {P_MIN:,} - {P_MAX:,}")
    print(f"modulus range             = {MOD_MIN:,} - {MOD_MAX:,}")
    print(f"close ratio               = {CLOSE_RATIO:.0%}")
    print(f"seed                      = {SEED:,}")
    print()

    print("=" * 92)
    print("BUILDING PRIME POOLS")
    print("=" * 92)

    factor_primes = build_factor_primes()
    modulus_primes = build_modulus_primes()

    print(f"factor primes             = {len([p for p in factor_primes if P_MIN <= p <= P_MAX]):,}")
    print(f"modulus primes            = {len(modulus_primes):,}")
    print()

    print("=" * 92)
    print("BUILDING CLOSE MODULUS PAIRS")
    print("=" * 92)

    close_pairs = build_close_pairs(modulus_primes)

    print(f"close modulus pairs       = {len(close_pairs):,}")
    print()

    print("=" * 92)
    print("BUILDING ANCHORS")
    print("=" * 92)

    anchors = build_anchors(
        factor_primes,
        N_ANCHORS,
        SEED,
    )

    print(f"actual anchors            = {len(anchors):,}")
    print()

    print("=" * 92)
    print("RUNNING COUPLED E/G SEARCH")
    print("=" * 92)

    results = []

    search_start = time.perf_counter()

    for idx, anchor in enumerate(anchors, start=1):

        r1, r2 = select_pair_for_anchor(
            close_pairs,
            anchor["n"],
            rng,
        )

        result = run_anchor(
            anchor,
            r1,
            r2,
        )

        results.append(result)

        if idx % PROGRESS_EVERY == 0 or idx == len(anchors):
            print(f"anchor {idx:3d}/{len(anchors)}")

    search_runtime = time.perf_counter() - search_start

    # ======================================================================================
    # SUMMARY
    # ======================================================================================

    def avg(key):
        vals = [r[key] for r in results]
        return statistics.mean(vals) if vals else 0.0

    recovered = sum(r["recovered"] for r in results)

    true_E_retained = sum(
        r["true_E_retained"]
        for r in results
    )

    true_carry_retained = sum(
        r["true_carry_retained"]
        for r in results
    )

    true_g_retained = sum(
        r["true_g_cell_retained"]
        for r in results
    )

    avg_b_before = avg("b_before")
    avg_b_after = avg("b_after")

    if avg_b_before:
        b_reduction_ratio = avg_b_after / avg_b_before
    else:
        b_reduction_ratio = 0.0

    if avg("carry_cells"):
        g_cell_ratio = avg("g_cells") / avg("carry_cells")
    else:
        g_cell_ratio = 0.0

    print()
    print("=" * 92)
    print("SUMMARY")
    print("=" * 92)
    print(f"anchors analyzed                = {len(results):,}")
    print()
    print(f"average E                       = {avg('E_true'):.3f}")
    print(f"average E-window upper bound    = {avg('E_max'):.3f}")
    print()
    print(f"average quotient pairs          = {avg('quotient_pairs'):.3f}")
    print(f"average carry cells             = {avg('carry_cells'):.3f}")
    print(f"average g-feasible cells        = {avg('g_cells'):.3f}")
    print()
    print(f"average b states before g       = {avg_b_before:.3f}")
    print(f"average b states after g        = {avg_b_after:.3f}")
    print(f"average exact division hits     = {avg('division_hits'):.3f}")
    print(f"average exact solutions         = {avg('exact_solutions'):.3f}")
    print()

    print("=" * 92)
    print("G COUPLING REDUCTION")
    print("=" * 92)
    print(
        f"g-feasible cells / carry cells = "
        f"{g_cell_ratio:.9f}"
    )
    print(
        f"cell reduction                 = "
        f"{(1.0 - g_cell_ratio) * 100.0:.6f}%"
    )
    print()
    print(
        f"b-after / b-before             = "
        f"{b_reduction_ratio:.9f}"
    )
    print(
        f"b-state reduction              = "
        f"{(1.0 - b_reduction_ratio) * 100.0:.6f}%"
    )
    print()

    print("=" * 92)
    print("RECOVERY")
    print("=" * 92)
    print(f"correctly recovered            = {recovered}/{len(results)}")
    print(
        f"recovery rate                  = "
        f"{recovered / len(results) * 100.0:.4f}%"
    )
    print(
        f"true E retained                = "
        f"{true_E_retained}/{len(results)}"
    )
    print(
        f"true carry retained            = "
        f"{true_carry_retained}/{len(results)}"
    )
    print(
        f"true g-cell retained           = "
        f"{true_g_retained}/{len(results)}"
    )
    print()

    # ======================================================================================
    # TRUE PATH DIAGNOSTICS
    # ======================================================================================

    print("=" * 92)
    print("TRUE-PATH EXAMPLES")
    print("=" * 92)

    for r in results[:20]:
        print(
            f"n={r['n']:,} "
            f"p={r['p']:,} "
            f"q={r['q']:,} "
            f"mods=({r['r1']},{r['r2']})"
        )

        print(
            f"    E={r['E_true']} "
            f"T={r['T']} "
            f"K={r['T'] - r['E_true']} "
            f"true(k,l)=({r['true_k']},{r['true_l']}) "
            f"carry=({r['true_c1']},{r['true_c2']},{r['true_c3']})"
        )

        print(
            f"    cells={r['carry_cells']} "
            f"gCells={r['g_cells']} "
            f"bBefore={r['b_before']} "
            f"bAfter={r['b_after']} "
            f"divHits={r['division_hits']} "
            f"exact={r['exact_solutions']} "
            f"recovered={r['recovered']}"
        )

        print()

    # ======================================================================================
    # STRONGEST G CELL COLLAPSES
    # ======================================================================================

    strongest = sorted(
        results,
        key=lambda r: (
            r["carry_cells"] - r["g_cells"],
            -r["b_after"],
        ),
        reverse=True,
    )[:20]

    print("=" * 92)
    print("STRONGEST G-CELL COLLAPSES")
    print("=" * 92)

    for r in strongest:
        ratio = (
            r["g_cells"] / r["carry_cells"]
            if r["carry_cells"]
            else 0.0
        )

        print(
            f"n={r['n']:,} "
            f"cells={r['carry_cells']} "
            f"gCells={r['g_cells']} "
            f"ratio={ratio:.6f} "
            f"bBefore={r['b_before']} "
            f"bAfter={r['b_after']} "
            f"div={r['division_hits']} "
            f"exact={r['exact_solutions']}"
        )

    # ======================================================================================
    # STRONGEST B-INTERVAL COLLAPSES
    # ======================================================================================

    strongest_b = sorted(
        results,
        key=lambda r: (
            (
                r["b_before"] - r["b_after"]
            )
            if r["b_before"]
            else 0,
            -r["b_after"],
        ),
        reverse=True,
    )[:20]

    print()
    print("=" * 92)
    print("STRONGEST B-STATE COLLAPSES")
    print("=" * 92)

    for r in strongest_b:
        ratio = (
            r["b_after"] / r["b_before"]
            if r["b_before"]
            else 0.0
        )

        print(
            f"n={r['n']:,} "
            f"bBefore={r['b_before']} "
            f"bAfter={r['b_after']} "
            f"ratio={ratio:.6f} "
            f"cells={r['carry_cells']} "
            f"gCells={r['g_cells']} "
            f"div={r['division_hits']} "
            f"exact={r['exact_solutions']}"
        )

    # ======================================================================================
    # WEAKEST B-STATE COLLAPSES
    # ======================================================================================

    weakest_b = sorted(
        results,
        key=lambda r: (
            (
                r["b_after"] / r["b_before"]
            )
            if r["b_before"]
            else 1.0,
            -r["b_before"],
        ),
        reverse=True,
    )[:20]

    print()
    print("=" * 92)
    print("WEAKEST B-STATE COLLAPSES")
    print("=" * 92)

    for r in weakest_b:
        ratio = (
            r["b_after"] / r["b_before"]
            if r["b_before"]
            else 1.0
        )

        print(
            f"n={r['n']:,} "
            f"bBefore={r['b_before']} "
            f"bAfter={r['b_after']} "
            f"ratio={ratio:.6f} "
            f"cells={r['carry_cells']} "
            f"gCells={r['g_cells']} "
            f"div={r['division_hits']} "
            f"exact={r['exact_solutions']}"
        )

    # ======================================================================================
    # DISTRIBUTIONS
    # ======================================================================================

    print()
    print("=" * 92)
    print("CELL DISTRIBUTIONS")
    print("=" * 92)

    g_distribution = Counter()

    for r in results:
        g_distribution[r["g_cells"]] += 1

    for cells, count in sorted(g_distribution.items())[:30]:
        print(
            f"g-feasible cells = {cells:5d} "
            f"anchors = {count:4d}"
        )

    # ======================================================================================
    # TIMING
    # ======================================================================================

    total_runtime = time.perf_counter() - start_total

    print()
    print("=" * 92)
    print("TIMING")
    print("=" * 92)
    print(f"search runtime                  = {search_runtime:.3f} s")
    print(f"total runtime                   = {total_runtime:.3f} s")
    print()

    # ======================================================================================
    # MATHEMATICAL INTERPRETATION
    # ======================================================================================

    print("=" * 92)
    print("MATHEMATICAL INTERPRETATION")
    print("=" * 92)

    print(
        r"""
The experiment starts with:

    p = a + k*r1
    q = b + l*r2

and:

    T = floor(n/(r1*r2))
    E = T-k*l.

The carry decomposition is:

    a*l = c1*r1 + d1
    b*k = c2*r2 + d2

with:

    E = c1+c2+c3.

The previous experiment used E to construct carry rectangles
in the hidden residue plane (a,b).

This experiment additionally uses:

    g = n mod (r1*r2).

Substituting the exact product identity gives:

    a(b+l*r2) + b*k*r1 = g + E*r1*r2.

Therefore:

    a =
        (g + E*r1*r2 - b*k*r1)
        /(b+l*r2).

For fixed:

    E,k,l,c1,c2,c3

the carry equations produce:

    a_lo <= a <= a_hi
    b_lo <= b <= b_hi.

The rational function a(b) is strictly decreasing.

Therefore the a interval can be inverted directly into a
new b interval.

This means g can potentially eliminate entire portions of a
carry cell BEFORE exact divisibility testing.

The experiment measures two separate effects:

    1. CELL PRUNING

       carry cells
           ->
       g-feasible cells.

    2. WITHIN-CELL PRUNING

       all b states inside carry cells
           ->
       b states compatible with g.

The strongest possible result is:

       many carry cells
           ->
       very few g-feasible cells

AND:

       huge b ranges
           ->
       tiny g-constrained b ranges

while maintaining:

       100% exact recovery.

That would indicate that the pair:

       (E,g)

contains more predictive information than E alone.

A negative result would be:

       almost every carry cell remains feasible
       and
       b ranges barely shrink.

That would show that g is mostly confirming the same
factorization relation rather than contributing an additional
independent constraint.

The final divisibility condition is still:

       p*q == n

and only exact solutions count.

The experiment deliberately does NOT use the true p or q during
candidate generation.
"""
    )

    print("=" * 92)
    print("EXPERIMENT COMPLETE")
    print("=" * 92)


if __name__ == "__main__":
    run()
