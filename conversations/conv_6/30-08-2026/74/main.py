#!/usr/bin/env python3
"""
============================================================================================
START EXPERIMENT 74
DIRECT STATIC-K LOOKUP -> n FACTORIZATION
============================================================================================

Purpose
-------
Test whether the auxiliary n+x stage can be removed completely.

For a chosen modulus product:

    R = r1 * r2

we have:

    T = floor(n / R)
    K = k*l
    E = T-K

Instead of obtaining K from a factored n+x, this experiment:

    1. builds a static lookup table K -> divisor pairs (k,l)
    2. chooses R ~= n / K_TARGET
    3. computes T
    4. scans a bounded K-window below T
    5. gets (k,l) directly from the static table
    6. tests the quotient cells against n
    7. verifies p*q == n

No n+x is constructed.
No n+x is factored.
No factorization of n is performed by SymPy.

The only final test is exact arithmetic.
"""

from __future__ import annotations

import math
import time
import random
from collections import defaultdict
from sympy import primerange


# ============================================================================================
# CONFIGURATION
# ============================================================================================

SCALES = (
    10**9,
    10**12,
    10**16,
)

ANCHORS_PER_SCALE = 6

K_TARGET = 1000

# Test several nearby R values around n / K_TARGET.
R_OFFSETS = (
    0.90,
    0.95,
    1.00,
    1.05,
    1.10,
)

# Static table covers every K from 1..STATIC_K_MAX.
STATIC_K_MAX = 5000

# Direct K scan:
#
# Since K = T-E and E >= 0, only K <= T is possible.
#
# This is the experimental search window below T.
#
# Increase this if you want to test how sensitive the method is
# to the assumed small-defect regime.
K_WINDOW = 128

# Maximum number of (k,l) pairs tested for one K.
MAX_PAIRS_PER_K = 5000

# Balanced prime search.
PRIME_RADIUS = 20000

SEED = 1511464998


# ============================================================================================
# FORMATTING
# ============================================================================================

def fmt(n: int) -> str:
    return f"{n:,}"


def scale_label(n: int) -> str:
    return f"{n:.0e}"


# ============================================================================================
# STATIC K -> (k,l) LOOKUP
# ============================================================================================

def build_static_k_table(max_k: int):
    """
    table[K] = tuple of all positive ordered divisor pairs (k,l)
    such that k*l = K.
    """
    table = defaultdict(list)

    for k in range(1, max_k + 1):
        for l in range(k, max_k // k + 1):
            K = k * l

            table[K].append((k, l))

            if k != l:
                table[K].append((l, k))

    return {
        K: tuple(sorted(pairs))
        for K, pairs in table.items()
    }


# ============================================================================================
# PRIME / ANCHOR GENERATION
# ============================================================================================

def nearest_prime(primes, target: int) -> int | None:
    """
    Return the closest prime in a sorted prime list.
    """
    import bisect

    idx = bisect.bisect_left(primes, target)

    candidates = []

    if idx < len(primes):
        candidates.append(primes[idx])

    if idx > 0:
        candidates.append(primes[idx - 1])

    if not candidates:
        return None

    return min(candidates, key=lambda p: abs(p - target))


def generate_anchor(scale: int, rng: random.Random):
    """
    Construct n=p*q with p,q around sqrt(scale).
    """
    root = math.isqrt(scale)

    lo = max(1000, int(root * 0.72))
    hi = int(root * 1.28)

    # Ensure a sensible interval.
    lo = min(lo, hi - 100)

    prime_pool = list(primerange(lo, hi + 1))

    if len(prime_pool) < 10:
        raise RuntimeError(
            f"Not enough primes for scale={scale}"
        )

    p = rng.choice(prime_pool)
    q = rng.choice(prime_pool)

    if p == q:
        q = rng.choice(prime_pool)

    n = p * q

    return min(p, q), max(p, q), n


# ============================================================================================
# BALANCED PRIME R1,R2 SELECTION
# ============================================================================================

def choose_balanced_prime_pair(
    target_R: int,
    modulus_prime_cache: list[int],
):
    """
    Find r1,r2 with r1*r2 close to target_R and
    r1,r2 roughly sqrt(target_R).

    Unlike the previous bounded version, this function does not
    depend on a fixed 300..3000 modulus range.
    """

    target_root = math.isqrt(max(2, target_R))

    # Search around sqrt(target_R).
    low = max(2, target_root - PRIME_RADIUS)
    high = target_root + PRIME_RADIUS

    # Filter cache to the local region.
    local = [
        p for p in modulus_prime_cache
        if low <= p <= high
    ]

    if not local:
        raise RuntimeError(
            f"No modulus primes near sqrt(R)={target_root}"
        )

    best = None
    best_error = None

    for r1 in local:
        wanted_r2 = target_R // r1
        r2 = nearest_prime(local, wanted_r2)

        if r2 is None:
            continue

        R = r1 * r2
        error = abs(R - target_R)

        if best_error is None or error < best_error:
            best = (r1, r2)
            best_error = error

    if best is None:
        raise RuntimeError(
            f"Could not construct balanced prime pair near R={target_R}"
        )

    return best


# ============================================================================================
# EXACT DIRECT TEST
# ============================================================================================

def reconstruct_from_kl(n: int, r1: int, r2: int, k: int, l: int):
    """
    Given r1,r2,k,l, derive the residue rectangle:

        p = a + k*r1
        q = b + l*r2

    with:

        0 <= a < r1
        0 <= b < r2.

    We then solve:

        (a + k*r1)(b + l*r2) = n.

    Rather than enumerate all a,b, use the fact that p must divide n.

    This function performs an exact restricted test.
    """

    p_low = k * r1
    p_high = (k + 1) * r1 - 1

    q_low = l * r2
    q_high = (l + 1) * r2 - 1

    if p_low <= 0 or q_low <= 0:
        return None

    # The candidate p must lie in [p_low, p_high] and divide n.
    #
    # Test the interval via gcd-derived divisor structure:
    # We still avoid scanning the entire original factor interval.
    #
    # For the small quotient cells used here, enumerate the much smaller
    # residue interval when it is computationally reasonable.

    width_p = p_high - p_low + 1

    # Hard protection against accidentally creating a huge search.
    if width_p > 2_000_000:
        return None

    for p in range(p_low, p_high + 1):
        if n % p != 0:
            continue

        q = n // p

        if q_low <= q <= q_high:
            return p, q

    return None


# ============================================================================================
# ONE K TEST
# ============================================================================================

def test_K(
    n: int,
    r1: int,
    r2: int,
    K: int,
    static_table,
):
    """
    Test every static divisor pair (k,l) for K.
    """

    pairs = static_table.get(K, ())

    if not pairs:
        return {
            "pairs": 0,
            "tested": 0,
            "solution": None,
        }

    if len(pairs) > MAX_PAIRS_PER_K:
        pairs = pairs[:MAX_PAIRS_PER_K]

    tested = 0

    for k, l in pairs:
        tested += 1

        result = reconstruct_from_kl(
            n=n,
            r1=r1,
            r2=r2,
            k=k,
            l=l,
        )

        if result is None:
            continue

        p, q = result

        if p * q != n:
            continue

        return {
            "pairs": len(pairs),
            "tested": tested,
            "solution": (p, q, k, l),
        }

    return {
        "pairs": len(pairs),
        "tested": tested,
        "solution": None,
    }


# ============================================================================================
# ONE R TEST
# ============================================================================================

def test_R(
    n: int,
    r1: int,
    r2: int,
    static_table,
):
    R = r1 * r2
    T = n // R

    if T <= 0:
        return {
            "R": R,
            "T": T,
            "solved": False,
            "K_tests": 0,
            "pair_tests": 0,
            "K_range": (),
            "K_used": None,
            "solution": None,
        }

    K_lo = max(1, T - K_WINDOW)
    K_hi = min(T, STATIC_K_MAX)

    if K_lo > K_hi:
        return {
            "R": R,
            "T": T,
            "solved": False,
            "K_tests": 0,
            "pair_tests": 0,
            "K_range": (),
            "K_used": None,
            "solution": None,
        }

    K_tests = 0
    pair_tests = 0

    # Try the values nearest T first.
    for K in range(K_hi, K_lo - 1, -1):
        if K not in static_table:
            continue

        K_tests += 1

        result = test_K(
            n=n,
            r1=r1,
            r2=r2,
            K=K,
            static_table=static_table,
        )

        pair_tests += result["tested"]

        if result["solution"] is not None:
            return {
                "R": R,
                "T": T,
                "solved": True,
                "K_tests": K_tests,
                "pair_tests": pair_tests,
                "K_range": (K_lo, K_hi),
                "K_used": K,
                "solution": result["solution"],
            }

    return {
        "R": R,
        "T": T,
        "solved": False,
        "K_tests": K_tests,
        "pair_tests": pair_tests,
        "K_range": (K_lo, K_hi),
        "K_used": None,
        "solution": None,
    }


# ============================================================================================
# SCALE RUN
# ============================================================================================

def run_scale(
    scale: int,
    rng: random.Random,
    static_table,
    modulus_primes,
):
    print("=" * 100)
    print(f"SCALE {scale_label(scale)}")
    print("=" * 100)

    rows = []

    for anchor_idx in range(1, ANCHORS_PER_SCALE + 1):
        p_true, q_true, n = generate_anchor(scale, rng)

        print(
            f"anchor {anchor_idx:2d}/{ANCHORS_PER_SCALE} "
            f"n={fmt(n)}"
        )

        true_records = []

        for offset in R_OFFSETS:
            target_R = max(
                2,
                int(round((n / K_TARGET) * offset))
            )

            r1, r2 = choose_balanced_prime_pair(
                target_R,
                modulus_primes,
            )

            R = r1 * r2

            k_true = p_true // r1
            l_true = q_true // r2
            K_true = k_true * l_true
            T_true = n // R
            E_true = T_true - K_true

            started = time.perf_counter()

            result = test_R(
                n=n,
                r1=r1,
                r2=r2,
                static_table=static_table,
            )

            elapsed = time.perf_counter() - started

            solved = result["solved"]

            valid = False
            recovered = None

            if solved:
                rp, rq, rk, rl = result["solution"]

                recovered = (rp, rq)
                valid = (
                    rp * rq == n
                    and {rp, rq} == {p_true, q_true}
                )

            rows.append({
                "scale": scale,
                "n": n,
                "p_true": p_true,
                "q_true": q_true,
                "offset": offset,
                "r1": r1,
                "r2": r2,
                "R": R,
                "K_true": K_true,
                "T_true": T_true,
                "E_true": E_true,
                "result": result,
                "elapsed": elapsed,
                "valid": valid,
            })

            true_records.append(rows[-1])

        # Compact per-anchor result.
        for rec in true_records:
            result = rec["result"]

            print(
                f"    offset={rec['offset']:.2f} "
                f"R={fmt(rec['R'])} "
                f"(r1,r2)=({fmt(rec['r1'])},{fmt(rec['r2'])}) "
                f"T={rec['T_true']:,} "
                f"K_true={rec['K_true']:,} "
                f"E={rec['E_true']:,} "
                f"K-range={result['K_range'][0] if result['K_range'] else '-'}.."
                f"{result['K_range'][1] if result['K_range'] else '-'} "
                f"K-used={result['K_used'] if result['K_used'] is not None else '-'} "
                f"K-tests={result['K_tests']} "
                f"pair-tests={result['pair_tests']} "
                f"solved={'YES' if result['solved'] else 'NO'} "
                f"valid={'YES' if rec['valid'] else 'NO'} "
                f"time={rec['elapsed']:.4f}s"
            )

    return rows


# ============================================================================================
# SUMMARY
# ============================================================================================

def print_summary(rows):
    print()
    print("-" * 100)
    print("DIRECT STATIC-K SUMMARY")
    print("-" * 100)

    by_scale = defaultdict(list)

    for row in rows:
        by_scale[row["scale"]].append(row)

    print(
        "scale       cases   solved   valid   "
        "avg K-tests   avg pair-tests   avg time"
    )
    print("-" * 100)

    for scale in SCALES:
        records = by_scale[scale]

        solved = sum(
            1 for r in records
            if r["result"]["solved"]
        )

        valid = sum(
            1 for r in records
            if r["valid"]
        )

        avg_K_tests = (
            sum(r["result"]["K_tests"] for r in records)
            / len(records)
        )

        avg_pair_tests = (
            sum(r["result"]["pair_tests"] for r in records)
            / len(records)
        )

        avg_time = (
            sum(r["elapsed"] for r in records)
            / len(records)
        )

        print(
            f"{scale_label(scale):<10} "
            f"{len(records):>6} "
            f"{solved:>8} "
            f"{valid:>7} "
            f"{avg_K_tests:>12.2f} "
            f"{avg_pair_tests:>16.2f} "
            f"{avg_time:>10.4f}s"
        )


# ============================================================================================
# MAIN
# ============================================================================================

def main():
    total_start = time.perf_counter()

    rng = random.Random(SEED)

    print("=" * 100)
    print("START EXPERIMENT 74")
    print("DIRECT STATIC-K LOOKUP -> n FACTORIZATION")
    print("=" * 100)

    print()
    print("configuration")
    print(f"    scales                  = {[scale_label(x) for x in SCALES]}")
    print(f"    anchors / scale         = {ANCHORS_PER_SCALE}")
    print(f"    K target                = {K_TARGET}")
    print(f"    K window                = {K_WINDOW}")
    print(f"    static K table max      = {STATIC_K_MAX}")
    print(f"    R offsets               = {R_OFFSETS}")
    print(f"    seed                    = {SEED}")

    print()
    print("=" * 100)
    print("BUILDING STATIC K LOOKUP TABLE")
    print("=" * 100)

    table_start = time.perf_counter()

    static_table = build_static_k_table(STATIC_K_MAX)

    table_time = time.perf_counter() - table_start

    total_pairs = sum(
        len(v)
        for v in static_table.values()
    )

    print(f"    K entries               = {len(static_table):,}")
    print(f"    divisor pairs           = {total_pairs:,}")
    print(f"    build time              = {table_time:.4f}s")

    print()
    print("=" * 100)
    print("BUILDING MODULUS PRIME POOL")
    print("=" * 100)

    # Large enough to support sqrt(R) up to roughly 4e6 for n~1e16.
    modulus_primes = list(
        primerange(2, 5_000_000)
    )

    print(
        f"    modulus primes          = {len(modulus_primes):,}"
    )

    all_rows = []

    for scale in SCALES:
        print()

        rows = run_scale(
            scale=scale,
            rng=rng,
            static_table=static_table,
            modulus_primes=modulus_primes,
        )

        all_rows.extend(rows)

    print_summary(all_rows)

    print()
    print("=" * 100)
    print("IMPORTANT COUNTERS")
    print("=" * 100)

    print(
        "    R tests are the number of distinct modulus products tested."
    )
    print(
        "    K tests are direct static-table K candidates."
    )
    print(
        "    pair tests are the resulting (k,l) divisor-pair tests."
    )
    print(
        "    No n+x construction is used."
    )
    print(
        "    No factorization of n is used."
    )
    print(
        "    Every successful result is verified by p*q == n."
    )

    solved = sum(
        1 for r in all_rows
        if r["result"]["solved"]
    )

    valid = sum(
        1 for r in all_rows
        if r["valid"]
    )

    print()
    print(
        f"    total cases             = {len(all_rows)}"
    )
    print(
        f"    solved                  = {solved}/{len(all_rows)}"
    )
    print(
        f"    valid                   = {valid}/{len(all_rows)}"
    )

    total_time = time.perf_counter() - total_start

    print()
    print("=" * 100)
    print("TIMING")
    print("=" * 100)
    print(
        f"    static lookup build     = {table_time:.4f}s"
    )
    print(
        f"    total runtime           = {total_time:.4f}s"
    )

    print()
    print("=" * 100)
    print("FINISHED EXPERIMENT 74")
    print("=" * 100)


if __name__ == "__main__":
    main()
