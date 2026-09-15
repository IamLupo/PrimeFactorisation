#!/usr/bin/env python3

"""
============================================================================================
SCALING TEST: DETERMINISTIC n+x -> Kx -> ORIGINAL K
============================================================================================

Scales:

    1e9
    1e12
    1e16

For each anchor:

    n = p*q

Choose:

    k = floor(p/r1)
    l = floor(q/r2)
    K = k*l

Construct deterministic:

    x = (-n) mod S

so:

    n+x = S*m.

Factor n+x with SymPy.

Every divisor pair (px,qx) of n+x gives:

    kx = floor(px/r1)
    lx = floor(qx/r2)
    Kx = kx*lx

A GOOD auxiliary factorization satisfies:

    Kx == K.

The original K is also factored with SymPy.

IMPORTANT:

    No search over x.

    No enumeration of the original p interval.

    No prime pool covering the full factor range.

For 1e16-scale numbers, factor(n+x) may be difficult.
A separate process timeout prevents the entire experiment from
getting stuck on one factorization.
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

ANCHORS_PER_SCALE = 4

MOD_MIN = 300
MOD_MAX = 3000

CLOSE_RATIO = 0.20

S_VALUES = (
    30,
    210,
    2310,
    30030,
)

# Maximum time allowed for one SymPy factorint(n+x).
#
# Increase for slower machines if desired.
#
FACTOR_TIMEOUT = 10.0

SHOW_GOOD_ONLY = False


# ==========================================================================================
# FACTOR WORKER
# ==========================================================================================

def factor_worker(
    value: int,
    queue: mp.Queue,
) -> None:

    try:
        result = sp.factorint(
            value,
            use_trial=True,
            use_rho=True,
            use_pm1=True,
            use_ecm=True,
        )

        queue.put(
            ("ok", result)
        )

    except Exception as exc:
        queue.put(
            ("error", repr(exc))
        )


def factor_with_timeout(
    value: int,
    timeout: float,
) -> tuple[str, dict[int, int] | None, float]:

    queue: mp.Queue = mp.Queue()

    process = mp.Process(
        target=factor_worker,
        args=(value, queue),
    )

    start = time.perf_counter()

    process.start()
    process.join(timeout)

    elapsed = time.perf_counter() - start

    if process.is_alive():

        process.terminate()
        process.join()

        return (
            "timeout",
            None,
            elapsed,
        )

    if queue.empty():

        return (
            "error",
            None,
            elapsed,
        )

    status, value_out = queue.get()

    if status == "ok":

        return (
            "ok",
            value_out,
            elapsed,
        )

    return (
        "error",
        None,
        elapsed,
    )


# ==========================================================================================
# FACTOR PAIRS
# ==========================================================================================

def factor_pairs(
    factors: dict[int, int],
) -> list[tuple[int, int]]:

    value = 1

    for prime, exponent in factors.items():
        value *= prime ** exponent

    divisors = [1]

    for prime, exponent in factors.items():

        old = list(divisors)

        power = 1

        for _ in range(exponent):

            power *= prime

            for d in old:
                divisors.append(
                    d * power
                )

    divisors = sorted(set(divisors))

    result = []

    for d in divisors:

        if d * d > value:
            break

        if value % d == 0:
            result.append(
                (d, value // d)
            )

    return result


# ==========================================================================================
# MODULUS PAIRS
# ==========================================================================================

def build_close_pairs() -> list[tuple[int, int]]:

    primes = list(
        sp.primerange(
            MOD_MIN,
            MOD_MAX + 1,
        )
    )

    result = []

    for i, r1 in enumerate(primes):

        for r2 in primes[i + 1:]:

            if r2 - r1 <= CLOSE_RATIO * r1:
                result.append(
                    (r1, r2)
                )

    return result


# ==========================================================================================
# DIRECT LARGE ANCHOR GENERATION
# ==========================================================================================

def random_prime_near(
    center: int,
    rng: random.Random,
) -> int:

    # Random offset gives different anchors while keeping
    # both primes approximately the same size.
    spread = max(
        1000,
        center // 10,
    )

    low = max(
        3,
        center - spread,
    )

    high = center + spread

    p = sp.randprime(
        low,
        high,
    )

    return int(p)


def build_anchor(
    target_lo: int,
    target_hi: int,
    rng: random.Random,
) -> tuple[int, int]:

    center = math.isqrt(
        math.isqrt(
            target_lo
        )
    )

    # Correct approximate factor scale:
    #
    # if n ~ 1e16 then p,q ~ 1e8.
    #
    factor_center = math.isqrt(
        (target_lo + target_hi) // 2
    )

    while True:

        p = random_prime_near(
            factor_center,
            rng,
        )

        q = random_prime_near(
            factor_center,
            rng,
        )

        if p == q:
            continue

        n = p * q

        if target_lo <= n <= target_hi:

            if p > q:
                p, q = q, p

            return p, q


# ==========================================================================================
# ORIGINAL K
# ==========================================================================================

def original_state(
    p: int,
    q: int,
    r1: int,
    r2: int,
) -> tuple[int, int, int, int, int]:

    n = p * q
    R = r1 * r2

    k = p // r1
    l = q // r2

    K = k * l

    T = n // R

    E = T - K

    return (
        k,
        l,
        K,
        T,
        E,
    )


# ==========================================================================================
# DETERMINISTIC x
# ==========================================================================================

def deterministic_x(
    n: int,
    S: int,
) -> int:

    x = (-n) % S

    if x == 0:
        x = S

    return x


# ==========================================================================================
# ONE AUXILIARY TEST
# ==========================================================================================

def test_auxiliary(
    n: int,
    r1: int,
    r2: int,
    K: int,
    S: int,
) -> dict:

    x = deterministic_x(
        n,
        S,
    )

    nx = n + x

    status, factors, factor_time = factor_with_timeout(
        nx,
        FACTOR_TIMEOUT,
    )

    result = {
        "S": S,
        "x": x,
        "nx": nx,
        "status": status,
        "factor_time": factor_time,
        "aux_pairs": 0,
        "good_pairs": 0,
        "good": [],
        "min_delta": None,
        "nx_factors": factors,
    }

    if status != "ok":
        return result

    pairs = factor_pairs(
        factors
    )

    R = r1 * r2
    Tx = nx // R

    for px, qx in pairs:

        # Avoid trivial factor pairs and keep auxiliary
        # factors reasonably close to the experimental domain.
        if px < 10_000:
            continue

        if qx < 10_000:
            continue

        result["aux_pairs"] += 1

        kx = px // r1
        lx = qx // r2

        Kx = kx * lx

        delta = abs(
            Kx - K
        )

        if (
            result["min_delta"] is None
            or delta < result["min_delta"]
        ):
            result["min_delta"] = delta

        if Kx == K:

            result["good_pairs"] += 1

            result["good"].append(
                (
                    px,
                    qx,
                    kx,
                    lx,
                    Kx,
                    Tx,
                    Tx - Kx,
                )
            )

    return result


# ==========================================================================================
# SCALE
# ==========================================================================================

def run_scale(
    scale_lo: int,
    scale_hi: int,
    anchor_count: int,
    rng: random.Random,
    modulus_pairs: list[tuple[int, int]],
) -> list[dict]:

    scale_label = (
        f"{scale_lo:.0e}"
    )

    print("=" * 100)
    print(
        f"BUILDING {anchor_count} ANCHORS FOR "
        f"{scale_label}"
    )
    print("=" * 100)

    anchors = []

    for _ in range(anchor_count):

        p, q = build_anchor(
            scale_lo,
            scale_hi,
            rng,
        )

        anchors.append(
            (p, q)
        )

    print(
        f"anchors                   = "
        f"{len(anchors)}"
    )

    print()

    results = []

    for index, (p, q) in enumerate(
        anchors,
        start=1,
    ):

        n = p * q

        r1, r2 = rng.choice(
            modulus_pairs
        )

        k, l, K, T, E = original_state(
            p,
            q,
            r1,
            r2,
        )

        K_status, K_factors, K_time = factor_with_timeout(
            K,
            FACTOR_TIMEOUT,
        )

        print(
            f"[{scale_label}] "
            f"anchor {index}/{anchor_count} "
            f"n={n:,} "
            f"mods=({r1},{r2}) "
            f"K={K:,}"
        )

        print(
            f"    original "
            f"(p,q)=({p:,},{q:,}) "
            f"(k,l)=({k},{l}) "
            f"T={T:,} "
            f"E={E:,}"
        )

        print(
            f"    factor(K): "
            f"status={K_status} "
            f"time={K_time:.4f}s "
            f"factors={K_factors}"
        )

        for S in S_VALUES:

            result = test_auxiliary(
                n=n,
                r1=r1,
                r2=r2,
                K=K,
                S=S,
            )

            record = {
                "scale": scale_lo,
                "p": p,
                "q": q,
                "n": n,
                "r1": r1,
                "r2": r2,
                "k": k,
                "l": l,
                "K": K,
                "T": T,
                "E": E,
                "K_status": K_status,
                "K_factors": K_factors,
                "K_time": K_time,
                **result,
            }

            results.append(record)

            if (
                SHOW_GOOD_ONLY
                and result["good_pairs"] == 0
            ):
                continue

            min_delta = result["min_delta"]

            if min_delta is None:
                min_delta_text = "n/a"
            else:
                min_delta_text = str(
                    min_delta
                )

            print(
                f"    S={S:5d} "
                f"x={result['x']:6d} "
                f"status={result['status']:<7} "
                f"auxPairs={result['aux_pairs']:4d} "
                f"good={result['good_pairs']:3d} "
                f"min|Kx-K|={min_delta_text:<8} "
                f"time={result['factor_time']:.4f}s"
            )

            if result["good_pairs"]:

                for (
                    px,
                    qx,
                    kx,
                    lx,
                    Kx,
                    Tx,
                    Ex,
                ) in result["good"]:

                    print(
                        f"        GOOD: "
                        f"(px,qx)=({px:,},{qx:,}) "
                        f"(kx,lx)=({kx},{lx}) "
                        f"Kx={Kx}"
                    )

        print()

    return results


# ==========================================================================================
# MAIN
# ==========================================================================================

def run() -> None:

    total_start = time.perf_counter()

    rng = random.Random(
        SEED
    )

    print("=" * 100)
    print(
        "DETERMINISTIC n+x -> Kx SCALING EXPERIMENT"
    )
    print("=" * 100)

    print(
        f"anchors / scale           = "
        f"{ANCHORS_PER_SCALE}"
    )

    print(
        f"S values                  = "
        f"{S_VALUES}"
    )

    print(
        f"factor timeout             = "
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
    print("BUILDING MODULUS POOL")
    print("=" * 100)

    modulus_pairs = build_close_pairs()

    print(
        f"close modulus pairs       = "
        f"{len(modulus_pairs):,}"
    )

    print()

    # --------------------------------------------------------------------------------------
    # SCALES
    # --------------------------------------------------------------------------------------

    scales = [
        (10**9, 10**10 - 1),
        (10**12, 10**13 - 1),
        (10**16, 10**17 - 1),
    ]

    all_results = []

    for scale_lo, scale_hi in scales:

        results = run_scale(
            scale_lo,
            scale_hi,
            ANCHORS_PER_SCALE,
            rng,
            modulus_pairs,
        )

        all_results.extend(
            results
        )

    # --------------------------------------------------------------------------------------
    # COMPACT FINAL RECORD
    # --------------------------------------------------------------------------------------

    print("=" * 100)
    print("FINAL RECORD")
    print("=" * 100)

    for scale in sorted(
        set(
            r["scale"]
            for r in all_results
        )
    ):

        records = [
            r
            for r in all_results
            if r["scale"] == scale
        ]

        print(
            f"scale={scale:.0e} "
            f"cases={len(records)}"
        )

        for r in records:

            print(
                f"    n={r['n']:,} "
                f"S={r['S']} "
                f"x={r['x']} "
                f"status={r['status']} "
                f"K={r['K']:,} "
                f"auxPairs={r['aux_pairs']} "
                f"good={r['good_pairs']} "
                f"minDelta={r['min_delta']}"
            )

    # --------------------------------------------------------------------------------------
    # TIMING
    # --------------------------------------------------------------------------------------

    elapsed = (
        time.perf_counter()
        - total_start
    )

    print()
    print("=" * 100)
    print("TIMING")
    print("=" * 100)

    print(
        f"total runtime              = "
        f"{elapsed:.3f}s"
    )

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    mp.freeze_support()
    run()
