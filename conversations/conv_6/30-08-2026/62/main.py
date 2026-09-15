#!/usr/bin/env python3

"""
============================================================================================
COMPACT S-SWEEP: Kx=K PROBABILITY / CROSS-CELL MATCH EXPERIMENT
============================================================================================

For

    n = p*q

and fixed moduli

    r1, r2

define

    k = floor(p/r1)
    l = floor(q/r2)
    K = k*l.

For every deterministic smooth S:

    x = (-n) mod S
    nx = n + x

so that

    nx == 0 (mod S).

We factor nx with SymPy and inspect every nontrivial factor pair:

    nx = px*qx

then

    kx = floor(px/r1)
    lx = floor(qx/r2)
    Kx = kx*lx.

We distinguish:

    EXACT-CELL:
        (kx,lx) == (k,l)

    CROSS-CELL:
        (kx,lx) != (k,l)
        but
        kx*lx == K

    NO-MATCH:
        Kx != K

The purpose is to measure whether Kx=K survives at larger n.

No search over x is performed.

x is completely determined by S.

No original p interval enumeration is performed.

Output is intentionally compact.
"""

from __future__ import annotations

import math
import multiprocessing as mp
import random
import time

import sympy as sp


# ==========================================================================================
# CONFIGURATION
# ==========================================================================================

SEED = 1_511_464_998

ANCHORS_PER_SCALE = 12

FACTOR_TIMEOUT = 10.0

MOD_MIN = 300
MOD_MAX = 3000

CLOSE_RATIO = 0.20

S_MAX = 2_000_000

SMOOTH_PRIMES = (
    2,
    3,
    5,
    7,
    11,
    13,
    17,
    19,
)


# ==========================================================================================
# SMOOTH S GENERATION
# ==========================================================================================

def build_s_values():

    values = set()

    # Squarefree smooth numbers.
    for mask in range(
        1,
        1 << len(SMOOTH_PRIMES),
    ):

        s = 1

        for i, p in enumerate(SMOOTH_PRIMES):

            if mask & (1 << i):

                s *= p

                if s > S_MAX:
                    break

        if s <= S_MAX:
            values.add(s)

    # Some smooth prime powers.
    bases = (
        2,
        3,
        5,
        6,
        10,
        30,
        42,
        70,
        210,
        330,
        390,
        570,
        870,
    )

    for base in bases:

        s = base

        while s <= S_MAX:

            values.add(s)

            s *= 2

    return sorted(values)


# ==========================================================================================
# MODULUS POOL
# ==========================================================================================

def build_modulus_pairs():

    primes = list(
        sp.primerange(
            MOD_MIN,
            MOD_MAX + 1,
        )
    )

    pairs = []

    for i, r1 in enumerate(primes):

        for r2 in primes[i + 1:]:

            if r2 - r1 <= CLOSE_RATIO * r1:

                pairs.append(
                    (r1, r2)
                )

    return pairs


# ==========================================================================================
# RANDOM PRIME GENERATION
# ==========================================================================================

def random_prime_near(
    target,
    rng,
):

    spread = max(
        1000,
        target // 8,
    )

    lo = max(
        3,
        target - spread,
    )

    hi = target + spread

    return int(
        sp.randprime(
            lo,
            hi,
        )
    )


def build_anchor(
    lo,
    hi,
    rng,
):

    center = math.isqrt(
        (lo + hi) // 2
    )

    while True:

        p = random_prime_near(
            center,
            rng,
        )

        q = random_prime_near(
            center,
            rng,
        )

        if p == q:
            continue

        n = p * q

        if lo <= n <= hi:

            if p > q:
                p, q = q, p

            return p, q


# ==========================================================================================
# FACTORIZATION
# ==========================================================================================

def factor_worker(
    n,
    queue,
):

    try:

        f = sp.factorint(
            n,
            use_trial=True,
            use_rho=True,
            use_pm1=True,
            use_ecm=True,
        )

        queue.put(
            ("ok", f)
        )

    except Exception:

        queue.put(
            ("error", None)
        )


def factor_with_timeout(
    n,
    timeout,
):

    q = mp.Queue()

    p = mp.Process(
        target=factor_worker,
        args=(n, q),
    )

    start = time.perf_counter()

    p.start()

    p.join(timeout)

    elapsed = (
        time.perf_counter()
        - start
    )

    if p.is_alive():

        p.terminate()
        p.join()

        return (
            "timeout",
            None,
            elapsed,
        )

    if q.empty():

        return (
            "error",
            None,
            elapsed,
        )

    status, factors = q.get()

    return (
        status,
        factors,
        elapsed,
    )


# ==========================================================================================
# DIVISOR PAIRS
# ==========================================================================================

def divisor_pairs(
    factors,
):

    n = 1

    for p, e in factors.items():
        n *= p ** e

    divisors = [1]

    for p, e in factors.items():

        old = list(divisors)

        power = 1

        for _ in range(e):

            power *= p

            for d in old:

                divisors.append(
                    d * power
                )

    divisors = sorted(
        set(divisors)
    )

    pairs = []

    for d in divisors:

        if d * d > n:
            break

        if n % d == 0:

            q = n // d

            pairs.append(
                (d, q)
            )

    return pairs


# ==========================================================================================
# ANCHOR TEST
# ==========================================================================================

def test_anchor(
    p,
    q,
    r1,
    r2,
    s_values,
):

    n = p * q

    k = p // r1
    l = q // r2

    K = k * l

    total_tests = 0
    successful_tests = 0
    timeouts = 0

    exact_cell_hits = 0
    cross_cell_hits = 0

    K_hits = 0

    min_distance = None

    s_good = []
    cross_examples = []

    total_aux_pairs = 0

    for S in s_values:

        x = (-n) % S

        if x == 0:
            x = S

        nx = n + x

        status, factors, elapsed = (
            factor_with_timeout(
                nx,
                FACTOR_TIMEOUT,
            )
        )

        total_tests += 1

        if status == "timeout":

            timeouts += 1
            continue

        if status != "ok":

            continue

        successful_tests += 1

        pairs = divisor_pairs(
            factors
        )

        S_has_K = False
        S_has_exact = False
        S_has_cross = False

        for px, qx in pairs:

            # Ignore very small/trivial factors.
            if px < 10_000:
                continue

            if qx < 10_000:
                continue

            total_aux_pairs += 1

            kx = px // r1
            lx = qx // r2

            Kx = kx * lx

            delta = abs(
                Kx - K
            )

            if (
                min_distance is None
                or delta < min_distance
            ):
                min_distance = delta

            if Kx != K:
                continue

            K_hits += 1
            S_has_K = True

            if (
                kx == k
                and lx == l
            ):

                exact_cell_hits += 1
                S_has_exact = True

            else:

                cross_cell_hits += 1
                S_has_cross = True

                if len(cross_examples) < 3:

                    cross_examples.append(
                        (
                            S,
                            x,
                            kx,
                            lx,
                            px,
                            qx,
                        )
                    )

        if S_has_K:

            s_good.append(
                (
                    S,
                    x,
                    S_has_exact,
                    S_has_cross,
                )
            )

    return {
        "n": n,
        "p": p,
        "q": q,
        "r1": r1,
        "r2": r2,
        "k": k,
        "l": l,
        "K": K,
        "total_tests": total_tests,
        "successful_tests": successful_tests,
        "timeouts": timeouts,
        "aux_pairs": total_aux_pairs,
        "K_hits": K_hits,
        "exact_cell_hits": exact_cell_hits,
        "cross_cell_hits": cross_cell_hits,
        "min_distance": min_distance,
        "s_good": s_good,
        "cross_examples": cross_examples,
    }


# ==========================================================================================
# SCALE
# ==========================================================================================

def run_scale(
    label,
    lo,
    hi,
    rng,
    modulus_pairs,
    s_values,
):

    print("=" * 100)
    print(
        f"SCALE {label}"
    )
    print("=" * 100)

    results = []

    for i in range(
        ANCHORS_PER_SCALE
    ):

        p, q = build_anchor(
            lo,
            hi,
            rng,
        )

        r1, r2 = rng.choice(
            modulus_pairs
        )

        result = test_anchor(
            p,
            q,
            r1,
            r2,
            s_values,
        )

        results.append(
            result
        )

        if (
            (i + 1) % max(
                1,
                ANCHORS_PER_SCALE // 4,
            )
            == 0
        ):

            print(
                f"    anchor "
                f"{i + 1:2d}/"
                f"{ANCHORS_PER_SCALE}"
            )

    # --------------------------------------------------------------------------------------
    # SUMMARY
    # --------------------------------------------------------------------------------------

    tests = sum(
        r["total_tests"]
        for r in results
    )

    successful = sum(
        r["successful_tests"]
        for r in results
    )

    timeouts = sum(
        r["timeouts"]
        for r in results
    )

    K_hits = sum(
        r["K_hits"]
        for r in results
    )

    exact_hits = sum(
        r["exact_cell_hits"]
        for r in results
    )

    cross_hits = sum(
        r["cross_cell_hits"]
        for r in results
    )

    aux_pairs = sum(
        r["aux_pairs"]
        for r in results
    )

    anchors_with_K = sum(
        bool(r["s_good"])
        for r in results
    )

    anchors_with_cross = sum(
        r["cross_cell_hits"] > 0
        for r in results
    )

    min_distances = [
        r["min_distance"]
        for r in results
        if r["min_distance"] is not None
    ]

    print()
    print(
        f"anchors                    = "
        f"{len(results)}"
    )

    print(
        f"S tests                    = "
        f"{tests:,}"
    )

    print(
        f"successful factorizations  = "
        f"{successful:,}"
    )

    print(
        f"timeouts                   = "
        f"{timeouts:,}"
    )

    print(
        f"auxiliary factor pairs     = "
        f"{aux_pairs:,}"
    )

    print()

    print(
        "Kx == K"
    )

    print(
        f"    K matches              = "
        f"{K_hits:,}"
    )

    print(
        f"    K-match / aux-pair     = "
        f"{K_hits / max(1, aux_pairs):.8f}"
    )

    print(
        f"    anchors with K match   = "
        f"{anchors_with_K}/"
        f"{len(results)}"
    )

    print()

    print(
        "CELL DECOMPOSITION"
    )

    print(
        f"    exact (kx,lx)==(k,l)  = "
        f"{exact_hits:,}"
    )

    print(
        f"    cross-cell K matches   = "
        f"{cross_hits:,}"
    )

    print(
        f"    anchors with cross     = "
        f"{anchors_with_cross}/"
        f"{len(results)}"
    )

    print()

    print(
        "MINIMUM |Kx-K|"
    )

    if min_distances:

        print(
            f"    mean                  = "
            f"{sum(min_distances) / len(min_distances):.3f}"
        )

        print(
            f"    minimum               = "
            f"{min(min_distances):,}"
        )

        print(
            f"    maximum               = "
            f"{max(min_distances):,}"
        )

    else:

        print(
            "    no valid auxiliary pairs"
        )

    print()

    print(
        "REPRESENTATIVE K-MATCHES"
    )

    shown = 0

    for r in results:

        if not r["s_good"]:
            continue

        first = r["s_good"][0]

        S, x, exact, cross = first

        print(
            f"    n={r['n']:,} "
            f"mods=({r['r1']},{r['r2']}) "
            f"K={r['K']:,} "
            f"S={S:,} "
            f"x={x:,} "
            f"exact={exact} "
            f"cross={cross}"
        )

        shown += 1

        if shown >= 5:
            break

    if shown == 0:

        print(
            "    none"
        )

    print()

    print(
        "CROSS-CELL EXAMPLES"
    )

    shown = 0

    for r in results:

        for (
            S,
            x,
            kx,
            lx,
            px,
            qx,
        ) in r["cross_examples"]:

            print(
                f"    n={r['n']:,} "
                f"S={S:,} "
                f"x={x:,} "
                f"orig=({r['k']},{r['l']}) "
                f"K={r['K']:,} "
                f"aux=({kx},{lx}) "
                f"px={px:,} "
                f"qx={qx:,}"
            )

            shown += 1

            if shown >= 5:
                break

        if shown >= 5:
            break

    if shown == 0:

        print(
            "    none"
        )

    print()

    return results


# ==========================================================================================
# MAIN
# ==========================================================================================

def run():

    total_start = time.perf_counter()

    rng = random.Random(
        SEED
    )

    print("=" * 100)
    print(
        "COMPACT S-SWEEP / "
        "Kx=K PROBABILITY EXPERIMENT"
    )
    print("=" * 100)

    print(
        f"anchors / scale           = "
        f"{ANCHORS_PER_SCALE}"
    )

    print(
        f"S maximum                 = "
        f"{S_MAX:,}"
    )

    print(
        f"factor timeout            = "
        f"{FACTOR_TIMEOUT:.1f}s"
    )

    print(
        f"seed                      = "
        f"{SEED:,}"
    )

    print()

    # --------------------------------------------------------------------------------------
    # MODULI
    # --------------------------------------------------------------------------------------

    print("=" * 100)
    print(
        "BUILDING MODULUS PAIRS"
    )
    print("=" * 100)

    modulus_pairs = build_modulus_pairs()

    print(
        f"close modulus pairs       = "
        f"{len(modulus_pairs):,}"
    )

    print()

    # --------------------------------------------------------------------------------------
    # S VALUES
    # --------------------------------------------------------------------------------------

    s_values = build_s_values()

    print("=" * 100)
    print(
        "BUILDING S SET"
    )
    print("=" * 100)

    print(
        f"S values                  = "
        f"{len(s_values):,}"
    )

    print(
        f"S range                   = "
        f"{min(s_values):,} .. "
        f"{max(s_values):,}"
    )

    print()

    # --------------------------------------------------------------------------------------
    # SCALES
    # --------------------------------------------------------------------------------------

    scales = (
        (
            "1e9",
            10**9,
            10**10 - 1,
        ),
        (
            "1e12",
            10**12,
            10**13 - 1,
        ),
        (
            "1e16",
            10**16,
            10**17 - 1,
        ),
    )

    all_results = {}

    for label, lo, hi in scales:

        all_results[label] = run_scale(
            label,
            lo,
            hi,
            rng,
            modulus_pairs,
            s_values,
        )

    # --------------------------------------------------------------------------------------
    # CROSS-SCALE COMPARISON
    # --------------------------------------------------------------------------------------

    print("=" * 100)
    print(
        "CROSS-SCALE COMPARISON"
    )
    print("=" * 100)

    print(
        f"{'scale':>8} "
        f"{'anchors':>8} "
        f"{'S-tested':>10} "
        f"{'K-match':>10} "
        f"{'K-match/S':>12} "
        f"{'exact-cell':>12} "
        f"{'cross-cell':>12}"
    )

    print("-" * 100)

    for label, results in all_results.items():

        n_anchors = len(results)

        total_s = sum(
            r["successful_tests"]
            for r in results
        )

        K_matches = sum(
            r["K_hits"]
            for r in results
        )

        exact = sum(
            r["exact_cell_hits"]
            for r in results
        )

        cross = sum(
            r["cross_cell_hits"]
            for r in results
        )

        print(
            f"{label:>8} "
            f"{n_anchors:8d} "
            f"{total_s:10d} "
            f"{K_matches:10d} "
            f"{K_matches / max(1,total_s):12.8f} "
            f"{exact:12d} "
            f"{cross:12d}"
        )

    # --------------------------------------------------------------------------------------
    # INTERPRETATION
    # --------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print(
        "INTERPRETATION"
    )
    print("=" * 100)

    print(
        """
The central observable is:

    Kx = kx*lx.

A successful auxiliary event is:

    Kx == K.

But this splits into two fundamentally different cases.

1. SAME CELL

       kx = k
       lx = l

   Then Kx=K is explained purely by quotient-cell stability.

2. CROSS CELL

       (kx,lx) != (k,l)

   but

       kx*lx = k*l.

   This is the stronger phenomenon because the auxiliary
   factorization lands in a different quotient cell while
   preserving the same multiplicative quotient invariant.

The scaling question is:

    Does the probability of Kx=K remain non-negligible
    as n grows?

An even stronger question is:

    Does the cross-cell probability remain non-negligible?

A useful recursive mechanism would require more than occasional
matches. There must eventually be a deterministic way to select
an S for which the auxiliary factorization produces a useful
K-match without testing a large fraction of all S values.

The experiment therefore treats:

    Kx == K

as an observed event, not as evidence of an algorithm.

All factor pairs are verified exactly.
"""
    )

    # --------------------------------------------------------------------------------------
    # TIMING
    # --------------------------------------------------------------------------------------

    elapsed = (
        time.perf_counter()
        - total_start
    )

    print("=" * 100)
    print("TIMING")
    print("=" * 100)

    print(
        f"total runtime              = "
        f"{elapsed:.3f}s"
    )

    print()
    print("=" * 100)
    print(
        "EXPERIMENT COMPLETE"
    )
    print("=" * 100)


if __name__ == "__main__":

    mp.freeze_support()

    run()
