#!/usr/bin/env python3

import math
import random
import statistics
import time
from collections import Counter


# ==============================================================================
# CONFIGURATION
# ==============================================================================

N_MIN = 10_000
N_MAX = 100_000

MOD_MIN = 300
MOD_MAX = 3_000

ANCHORS = 300
CLOSE_RATIO = 0.20
SEED = 1_511_464_998

M = 111_546_435


# ==============================================================================
# BASIC UTILITIES
# ==============================================================================

def sieve(limit: int) -> list[int]:
    if limit < 2:
        return []

    composite = bytearray(limit + 1)
    composite[0:2] = b"\x01\x01"

    for p in range(2, math.isqrt(limit) + 1):
        if not composite[p]:
            start = p * p
            composite[start:limit + 1:p] = b"\x01" * (
                ((limit - start) // p) + 1
            )

    return [i for i in range(2, limit + 1) if not composite[i]]


def ceil_div(a: int, b: int) -> int:
    return -((-a) // b)


def gcd3(a: int, b: int, c: int) -> int:
    return math.gcd(math.gcd(a, b), c)


# ==============================================================================
# PRIME POOLS
# ==============================================================================

def build_prime_pools():
    factor_primes = [
        p for p in sieve(N_MAX)
        if N_MIN <= p <= N_MAX
    ]

    modulus_primes = [
        p for p in sieve(MOD_MAX)
        if MOD_MIN <= p <= MOD_MAX
    ]

    return factor_primes, modulus_primes


# ==============================================================================
# ANCHORS
# ==============================================================================

def build_anchors(
    factor_primes: list[int],
    count: int,
    seed: int,
) -> list[tuple[int, int, int]]:
    """
    Return tuples:

        (n, p, q)

    with p <= q and n = p*q.
    """

    rng = random.Random(seed)

    anchors = []
    seen = set()

    # Keep products inside the intended scale.
    tries = 0
    max_tries = count * 500

    while len(anchors) < count and tries < max_tries:
        tries += 1

        p, q = rng.sample(factor_primes, 2)
        if p > q:
            p, q = q, p

        n = p * q

        # Keep the same useful experimental size range.
        if n < 100_000_000:
            continue

        # Avoid pathological duplicate products.
        if n in seen:
            continue

        seen.add(n)
        anchors.append((n, p, q))

    if len(anchors) < count:
        raise RuntimeError(
            f"Could only construct {len(anchors)} anchors out of {count}."
        )

    return anchors


# ==============================================================================
# CLOSE MODULUS TRIPLE
# ==============================================================================

def select_close_triple(
    modulus_primes: list[int],
    n: int,
    close_ratio: float,
) -> tuple[int, int, int]:
    """
    Select r1 < r2 < r3 such that:

        R = r1*r2*r3

    lies just below n whenever possible.

    We score by absolute distance |n-R|.
    """

    best = None

    # Work backwards from the cube-root scale.
    target = n ** (1.0 / 3.0)

    candidates = sorted(
        modulus_primes,
        key=lambda x: abs(x - target)
    )[:120]

    for i, r1 in enumerate(candidates):
        for r2 in candidates[i + 1:]:
            if r2 > r1 * (1.0 + close_ratio):
                continue

            for r3 in candidates:
                if r3 <= r2:
                    continue
                if r3 > r2 * (1.0 + close_ratio):
                    continue

                R = r1 * r2 * r3

                # Require the useful R < n regime.
                if R >= n:
                    continue

                score = abs(n - R)

                if best is None or score < best[0]:
                    best = (score, r1, r2, r3)

    if best is None:
        raise RuntimeError(f"No valid close triple for n={n}")

    _, r1, r2, r3 = best
    return r1, r2, r3


# ==============================================================================
# TWO-QUOTIENT THEOREM
# ==============================================================================

def two_quotient_candidates(
    n: int,
    R: int,
    p_min: int,
    p_max: int,
):
    """
    We have:

        n = p*q
        R = p*q - g
        g = n-R

    Let:

        k = floor(R/p)

    Then:

        R = k*p + s
        0 <= s < p

    and therefore:

        n = k*p + (s+g).

    If:

        0 <= g < p

    then:

        0 <= s+g < 2p

    and hence:

        q is either k or k+1.

    This function enumerates k, not p.

    It generates at most two exact divisibility candidates:

        q = k
        q = k+1

    and reconstructs:

        p = n/q.
    """

    g = n - R

    if g < 0:
        return [], {
            "valid_regime": False,
            "reason": "R >= n",
            "k_states": 0,
        }

    if g >= p_min:
        return [], {
            "valid_regime": False,
            "reason": "gap_not_less_than_p_min",
            "k_states": 0,
        }

    # Since p <= p_max:
    #
    #     floor(R/p) >= floor(R/p_max)
    #
    # and since p >= p_min:
    #
    #     floor(R/p) <= floor(R/p_min)
    #
    k_min = R // p_max
    k_max = R // p_min

    candidates = []
    k_states = 0

    for k in range(k_min, k_max + 1):
        k_states += 1

        for q in (k, k + 1):
            if q <= 0:
                continue

            if n % q != 0:
                continue

            p = n // q

            if p < p_min or p > p_max:
                continue

            if p > q:
                continue

            if p * q != n:
                continue

            if p not in candidates:
                candidates.append(p)

    return candidates, {
        "valid_regime": True,
        "reason": "g < p_min",
        "gap": g,
        "k_min": k_min,
        "k_max": k_max,
        "k_states": k_states,
    }


# ==============================================================================
# CRT-LIFT DIAGNOSTICS
# ==============================================================================

def crt_lift_diagnostics(
    n: int,
    p: int,
    q: int,
    mods: tuple[int, int, int],
):
    r1, r2, r3 = mods
    R = r1 * r2 * r3
    g = n - R

    # Since R > p,q in the intended regime:
    p_res = p % R
    q_res = q % R

    # Exact modular product.
    modular_ok = (p_res * q_res - n) % R == 0

    # The quotient index.
    k = R // p
    s = R % p

    # n = k*p + (s+g)
    identity_ok = (
        n == k * p + (s + g)
    )

    # q should be k or k+1 whenever 0 <= g < p.
    two_quotient_ok = (
        not (0 <= g < p)
        or q in (k, k + 1)
    )

    return {
        "R": R,
        "gap": g,
        "p_res": p_res,
        "q_res": q_res,
        "modular_ok": modular_ok,
        "k": k,
        "s": s,
        "identity_ok": identity_ok,
        "two_quotient_ok": two_quotient_ok,
        "R_gt_pmax": R > q,
    }


# ==============================================================================
# MAIN EXPERIMENT
# ==============================================================================

def run():
    print("=" * 100)
    print("THREE-CLOSE-PRIME CRT-LIFT / TWO-QUOTIENT GAP EXPERIMENT")
    print("=" * 100)
    print(f"M                         = {M:,}")
    print(f"factor range              = {N_MIN:,} - {N_MAX:,}")
    print(f"modulus prime range       = {MOD_MIN:,} - {MOD_MAX:,}")
    print(f"close ratio               = {CLOSE_RATIO:.0%}")
    print(f"anchors                   = {ANCHORS}")
    print(f"seed                      = {SEED}")
    print()

    t_total = time.perf_counter()

    print("=" * 100)
    print("BUILDING PRIME POOLS")
    print("=" * 100)

    factor_primes, modulus_primes = build_prime_pools()

    print(f"factor primes             = {len(factor_primes):,}")
    print(f"modulus primes            = {len(modulus_primes):,}")
    print()

    print("=" * 100)
    print("BUILDING ANCHORS")
    print("=" * 100)

    anchors = build_anchors(
        factor_primes,
        ANCHORS,
        SEED,
    )

    print(f"actual anchors            = {len(anchors)}")
    print()

    print("=" * 100)
    print("SELECTING CLOSE MODULUS TRIPLES")
    print("=" * 100)

    anchor_data = []

    for i, (n, p, q) in enumerate(anchors, 1):
        mods = select_close_triple(
            modulus_primes,
            n,
            CLOSE_RATIO,
        )

        r1, r2, r3 = mods
        R = r1 * r2 * r3

        anchor_data.append({
            "n": n,
            "p": p,
            "q": q,
            "mods": mods,
            "R": R,
        })

        if i % 25 == 0:
            print(f"anchor {i:3d}/{len(anchors)}")

    unique_triples = len({
        x["mods"] for x in anchor_data
    })

    print(f"unique modulus triples    = {unique_triples}")
    print()

    # --------------------------------------------------------------------------
    # VALIDATION
    # --------------------------------------------------------------------------

    print("=" * 100)
    print("VALIDATING CRT-LIFT IDENTITIES")
    print("=" * 100)

    validation_failures = 0
    two_quotient_failures = 0
    R_failures = 0

    for x in anchor_data:
        d = crt_lift_diagnostics(
            x["n"],
            x["p"],
            x["q"],
            x["mods"],
        )

        x["diag"] = d

        if not d["modular_ok"]:
            validation_failures += 1

        if not d["identity_ok"]:
            validation_failures += 1

        if not d["two_quotient_ok"]:
            two_quotient_failures += 1

        if not d["R_gt_pmax"]:
            R_failures += 1

    print(f"identity failures         = {validation_failures}")
    print(f"two-quotient failures     = {two_quotient_failures}")
    print(f"R <= q failures           = {R_failures}")
    print()

    # --------------------------------------------------------------------------
    # SEARCH
    # --------------------------------------------------------------------------

    print("=" * 100)
    print("RUNNING DIRECT k-SPACE SEARCH")
    print("=" * 100)

    results = []

    for i, x in enumerate(anchor_data, 1):
        n = x["n"]
        true_p = x["p"]
        true_q = x["q"]
        R = x["R"]

        t0 = time.perf_counter()

        candidates, stats = two_quotient_candidates(
            n=n,
            R=R,
            p_min=N_MIN,
            p_max=N_MAX,
        )

        elapsed = (time.perf_counter() - t0) * 1000.0

        recovered = (
            true_p in candidates
            or true_q in candidates
        )

        x["search"] = {
            **stats,
            "candidates": candidates,
            "recovered": recovered,
            "time_ms": elapsed,
        }

        results.append(x)

        if i % 25 == 0:
            print(f"anchor {i:3d}/{len(anchor_data)}")

    # --------------------------------------------------------------------------
    # SUMMARY
    # --------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("SUMMARY")
    print("=" * 100)

    active = [
        x for x in results
        if x["search"]["valid_regime"]
    ]

    inactive = len(results) - len(active)

    avg_p_states = statistics.mean(
        [x["q"] - x["p"] + 1 for x in results]
    )

    avg_prime_tests = statistics.mean(
        [
            sum(
                1
                for p in factor_primes
                if p <= math.isqrt(x["n"])
            )
            for x in results
        ]
    )

    if active:
        avg_k_states = statistics.mean(
            [x["search"]["k_states"] for x in active]
        )

        avg_candidates = statistics.mean(
            [
                len(x["search"]["candidates"])
                for x in active
            ]
        )

        recovered = sum(
            1 for x in active
            if x["search"]["recovered"]
        )

        avg_time = statistics.mean(
            [x["search"]["time_ms"] for x in active]
        )
    else:
        avg_k_states = 0.0
        avg_candidates = 0.0
        recovered = 0
        avg_time = 0.0

    print(f"anchors analyzed          = {len(results)}")
    print(f"g < p_min regime          = {len(active)}/{len(results)}")
    print(f"outside regime            = {inactive}")
    print()
    print(f"average integer p states  = {avg_p_states:,.3f}")
    print(f"average prime tests       = {avg_prime_tests:,.3f}")
    print(f"average k states          = {avg_k_states:,.3f}")
    print(f"average final candidates  = {avg_candidates:,.3f}")
    print(f"average search time       = {avg_time:.3f} ms")
    print()

    # --------------------------------------------------------------------------
    # REDUCTION
    # --------------------------------------------------------------------------

    print("=" * 100)
    print("SEARCH REDUCTION")
    print("=" * 100)

    if active:
        print(
            f"k / integer-p baseline    = "
            f"{avg_k_states / avg_p_states:.9f}"
        )

        print(
            f"k reduction                = "
            f"{100.0 * (1.0 - avg_k_states / avg_p_states):.6f}%"
        )

        print(
            f"k / prime baseline        = "
            f"{avg_k_states / avg_prime_tests:.9f}"
        )

        print(
            f"k vs prime tests          = "
            f"{avg_k_states / avg_prime_tests:.3f}x"
        )

    print()

    # --------------------------------------------------------------------------
    # GAP STATISTICS
    # --------------------------------------------------------------------------

    gaps = [x["diag"]["gap"] for x in results]

    print("=" * 100)
    print("GAP STATISTICS")
    print("=" * 100)

    print(f"mean gap                   = {statistics.mean(gaps):,.3f}")
    print(f"minimum gap                = {min(gaps):,}")
    print(f"maximum gap                = {max(gaps):,}")

    regime_count = sum(
        1 for x in results
        if x["diag"]["gap"] < x["p"]
    )

    print(
        f"gap < true p               = "
        f"{regime_count}/{len(results)}"
    )

    print()

    # --------------------------------------------------------------------------
    # DISTRIBUTION OF K STATES
    # --------------------------------------------------------------------------

    print("=" * 100)
    print("k-STATE DISTRIBUTION")
    print("=" * 100)

    if active:
        counts = Counter(
            x["search"]["k_states"]
            for x in active
        )

        for value, count in sorted(counts.items())[:25]:
            print(
                f"k states = {value:8,d} anchors = {count:4,d}"
            )

    print()

    # --------------------------------------------------------------------------
    # EXAMPLES
    # --------------------------------------------------------------------------

    print("=" * 100)
    print("EXAMPLES")
    print("=" * 100)

    for x in results[:20]:
        d = x["diag"]
        s = x["search"]

        print(
            f"n={x['n']:,} "
            f"p={x['p']:,} "
            f"q={x['q']:,} "
            f"mods={x['mods']} "
            f"R/n={x['R']/x['n']:.10f}"
        )

        print(
            f"    gap={d['gap']:,} "
            f"k={d['k']:,} "
            f"kRange="
            f"[{s.get('k_min', 0):,},{s.get('k_max', 0):,}] "
            f"kStates={s.get('k_states', 0):,} "
            f"candidates={s.get('candidates', [])} "
            f"recovered={s['recovered']}"
        )

    # --------------------------------------------------------------------------
    # STRONGEST REDUCTIONS
    # --------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("STRONGEST k-SPACE COLLAPSES")
    print("=" * 100)

    ranked = sorted(
        active,
        key=lambda x: x["search"]["k_states"]
    )

    for x in ranked[:20]:
        print(
            f"n={x['n']:,} "
            f"p={x['p']:,} "
            f"q={x['q']:,} "
            f"gap={x['diag']['gap']:,} "
            f"kStates={x['search']['k_states']:,} "
            f"candidates={len(x['search']['candidates'])} "
            f"recovered={x['search']['recovered']}"
        )

    # --------------------------------------------------------------------------
    # WEAKEST REDUCTIONS
    # --------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("WEAKEST k-SPACE COLLAPSES")
    print("=" * 100)

    for x in ranked[-20:]:
        print(
            f"n={x['n']:,} "
            f"p={x['p']:,} "
            f"q={x['q']:,} "
            f"gap={x['diag']['gap']:,} "
            f"kStates={x['search']['k_states']:,} "
            f"candidates={len(x['search']['candidates'])} "
            f"recovered={x['search']['recovered']}"
        )

    # --------------------------------------------------------------------------
    # ANCHOR TABLE
    # --------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("ANCHOR RESULTS")
    print("=" * 100)

    print(
        f"{'ID':>4} "
        f"{'p':>9} "
        f"{'q':>9} "
        f"{'gap':>9} "
        f"{'k':>8} "
        f"{'K':>8} "
        f"{'C':>4} "
        f"{'EXACT':>6}"
    )

    print("-" * 100)

    for i, x in enumerate(results, 1):
        s = x["search"]
        d = x["diag"]

        print(
            f"{i:4d} "
            f"{x['p']:9,d} "
            f"{x['q']:9,d} "
            f"{d['gap']:9,d} "
            f"{d['k']:8,d} "
            f"{s.get('k_states', 0):8,d} "
            f"{len(s['candidates']):4d} "
            f"{str(s['recovered']):>6}"
        )

    # --------------------------------------------------------------------------
    # INTERPRETATION
    # --------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("MATHEMATICAL INTERPRETATION")
    print("=" * 100)

    print(
        """
The previous experiment established that when:

    R > p_max

the modular hyperbola condition forces:

    p*q = n.

So the t-coordinate disappears.

This experiment asks whether the remaining factor search can be
reduced by using the small gap:

    g = n-R.

For a true factorization:

    n = p*q
    R = p*q-g.

Let:

    k = floor(R/p)

Then:

    R = k*p+s
    with 0 <= s < p.

Therefore:

    n = k*p + (s+g).

If:

    0 <= g < p,

then:

    0 <= s+g < 2p.

Because n = p*q, this implies:

    q is either k or k+1.

Thus the factor pair can be generated from the one-dimensional
quotient index:

    k
      |
      +--> q = k
      |
      +--> q = k+1
              |
              v
          p = n/q

The script therefore does NOT scan every integer p.

It enumerates the allowed k interval:

    floor(R/p_max) <= k <= floor(R/p_min)

and generates p only when n is exactly divisible by k or k+1.

The critical measurement is:

    k states / integer-p states

but there is a second comparison:

    k states / ordinary prime enumeration.

A useful result requires both to be small.

There is also an important limitation:

    the theorem requires g < p.

If the gap is larger than the factor p, then q is no longer
guaranteed to be only k or k+1.

Therefore the script separately reports how often the useful regime
actually occurs.

This is a direct test of whether the closeness of R to n provides
a genuine candidate-generation mechanism rather than merely a
post-enumeration filter.

Every candidate is verified with:

    p*q == n.
"""
    )

    # --------------------------------------------------------------------------
    # TIMING
    # --------------------------------------------------------------------------

    total_time = time.perf_counter() - t_total

    print()
    print("=" * 100)
    print("TIMING")
    print("=" * 100)
    print(f"total runtime              = {total_time:.3f} seconds")

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()
