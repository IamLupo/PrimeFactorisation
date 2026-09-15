#!/usr/bin/env python3

"""
============================================================================================
DETERMINISTIC S-SWEEP / n+x -> Kx SCALING EXPERIMENT
============================================================================================

Research target
---------------

Given

    n = p*q

and fixed close moduli

    r1 < r2,

define

    k = floor(p/r1)
    l = floor(q/r2)
    K = k*l.

Instead of searching over x, choose x deterministically from S:

    x_S = (-n) mod S

so that

    n + x_S = 0 mod S.

The auxiliary integer is then factored:

    n+x = product(prime powers)

and every factor pair (px,qx) is converted into

    kx = floor(px/r1)
    lx = floor(qx/r2)

    Kx = kx*lx.

The important event is:

    Kx == K.

This experiment sweeps MANY deterministic smooth S values to test
whether the occurrence of Kx == K depends systematically on S.

No search over x is performed.

No enumeration of the original p interval is performed.

SymPy is used for factorization.

The script is deliberately compact in output.
"""

from __future__ import annotations

import math
import multiprocessing as mp
import random
import time
from dataclasses import dataclass

import sympy as sp


# ==========================================================================================
# CONFIGURATION
# ==========================================================================================

SEED = 1_511_464_998

ANCHORS_PER_SCALE = 4

FACTOR_TIMEOUT = 10.0

MOD_MIN = 300
MOD_MAX = 3000
CLOSE_RATIO = 0.20

# ------------------------------------------------------------------------------------------
# Smooth S values
#
# They are generated as products of small primes.
#
# We deliberately include multiple sizes and families rather than
# only 30, 210, 2310, 30030.
#
# ------------------------------------------------------------------------------------------

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

S_MAX = 2_000_000

SHOW_GOOD_ONLY = False


# ==========================================================================================
# DATA
# ==========================================================================================

@dataclass
class Case:
    scale: int
    n: int

    p: int
    q: int

    r1: int
    r2: int

    k: int
    l: int
    K: int

    S: int
    x: int
    nx: int

    factor_status: str
    factor_time: float

    aux_pairs: int
    good_pairs: int

    min_delta: int | None

    good: list[tuple[int, int, int, int]]


# ==========================================================================================
# SMOOTH NUMBERS
# ==========================================================================================

def build_s_values() -> list[int]:

    values = set()

    # Products of distinct small primes.
    #
    # We deliberately include a broad collection of smooth S values.
    #
    values.add(2)
    values.add(6)
    values.add(10)
    values.add(14)
    values.add(30)
    values.add(42)
    values.add(70)
    values.add(210)

    for mask in range(
        1,
        1 << len(SMOOTH_PRIMES),
    ):

        value = 1

        for i, p in enumerate(SMOOTH_PRIMES):

            if mask & (1 << i):
                value *= p

                if value > S_MAX:
                    break

        if value <= S_MAX:
            values.add(value)

    # Additional smooth powers.
    bases = [
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
    ]

    for base in bases:

        value = base

        while value <= S_MAX:
            values.add(value)
            value *= 2

    return sorted(values)


# ==========================================================================================
# MODULI
# ==========================================================================================

def build_close_modulus_pairs() -> list[tuple[int, int]]:

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
# LARGE ANCHOR GENERATION
# ==========================================================================================

def random_prime_near(
    center: int,
    rng: random.Random,
) -> int:

    spread = max(
        1000,
        center // 10,
    )

    lo = max(
        3,
        center - spread,
    )

    hi = center + spread

    return int(
        sp.randprime(
            lo,
            hi,
        )
    )


def build_anchor(
    lo: int,
    hi: int,
    rng: random.Random,
) -> tuple[int, int]:

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
# ORIGINAL STATE
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
# FACTORIZATION WORKER
# ==========================================================================================

def factor_worker(
    n: int,
    queue: mp.Queue,
) -> None:

    try:

        factors = sp.factorint(
            n,
            use_trial=True,
            use_rho=True,
            use_pm1=True,
            use_ecm=True,
        )

        queue.put(
            ("ok", factors)
        )

    except Exception as exc:

        queue.put(
            ("error", repr(exc))
        )


def factor_with_timeout(
    n: int,
    timeout: float,
) -> tuple[str, dict[int, int] | None, float]:

    queue: mp.Queue = mp.Queue()

    process = mp.Process(
        target=factor_worker,
        args=(n, queue),
    )

    start = time.perf_counter()

    process.start()

    process.join(
        timeout
    )

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

    status, value = queue.get()

    if status != "ok":

        return (
            "error",
            None,
            elapsed,
        )

    reconstructed = 1

    for p, e in value.items():
        reconstructed *= p ** e

    if reconstructed != n:

        return (
            "error",
            None,
            elapsed,
        )

    return (
        "ok",
        value,
        elapsed,
    )


# ==========================================================================================
# DIVISOR PAIRS
# ==========================================================================================

def divisor_pairs(
    factors: dict[int, int],
) -> list[tuple[int, int]]:

    value = 1

    for p, e in factors.items():
        value *= p ** e

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
# SINGLE S VALUE
# ==========================================================================================

def test_S(
    scale: int,
    p: int,
    q: int,
    r1: int,
    r2: int,
    K: int,
    S: int,
) -> Case:

    n = p * q

    x = (-n) % S

    if x == 0:
        x = S

    nx = n + x

    if nx % S != 0:

        raise AssertionError(
            "Deterministic S construction failed."
        )

    status, factors, factor_time = (
        factor_with_timeout(
            nx,
            FACTOR_TIMEOUT,
        )
    )

    if status != "ok":

        return Case(
            scale=scale,
            n=n,
            p=p,
            q=q,
            r1=r1,
            r2=r2,
            k=p // r1,
            l=q // r2,
            K=K,
            S=S,
            x=x,
            nx=nx,
            factor_status=status,
            factor_time=factor_time,
            aux_pairs=0,
            good_pairs=0,
            min_delta=None,
            good=[],
        )

    pairs = divisor_pairs(
        factors
    )

    aux_pairs = 0
    good_pairs = 0

    min_delta = None

    good = []

    for px, qx in pairs:

        # Ignore trivial factors.
        if px < 10_000:
            continue

        if qx < 10_000:
            continue

        aux_pairs += 1

        kx = px // r1
        lx = qx // r2

        Kx = kx * lx

        delta = abs(
            Kx - K
        )

        if (
            min_delta is None
            or delta < min_delta
        ):
            min_delta = delta

        if Kx == K:

            good_pairs += 1

            good.append(
                (
                    px,
                    qx,
                    kx,
                    lx,
                )
            )

    return Case(
        scale=scale,
        n=n,
        p=p,
        q=q,
        r1=r1,
        r2=r2,
        k=p // r1,
        l=q // r2,
        K=K,
        S=S,
        x=x,
        nx=nx,
        factor_status=status,
        factor_time=factor_time,
        aux_pairs=aux_pairs,
        good_pairs=good_pairs,
        min_delta=min_delta,
        good=good,
    )


# ==========================================================================================
# SCALE TEST
# ==========================================================================================

def run_scale(
    scale_lo: int,
    scale_hi: int,
    rng: random.Random,
    modulus_pairs: list[tuple[int, int]],
    S_values: list[int],
) -> list[Case]:

    label = f"{scale_lo:.0e}"

    print("=" * 100)
    print(
        f"BUILDING ANCHORS FOR {label}"
    )
    print("=" * 100)

    anchors = []

    for _ in range(
        ANCHORS_PER_SCALE
    ):

        anchors.append(
            build_anchor(
                scale_lo,
                scale_hi,
                rng,
            )
        )

    print(
        f"anchors                   = "
        f"{len(anchors)}"
    )

    print()

    all_cases = []

    for anchor_id, (p, q) in enumerate(
        anchors,
        start=1,
    ):

        r1, r2 = rng.choice(
            modulus_pairs
        )

        k, l, K, T, E = (
            original_state(
                p,
                q,
                r1,
                r2,
            )
        )

        print(
            f"[{label}] anchor "
            f"{anchor_id}/{len(anchors)} "
            f"n={p*q:,} "
            f"mods=({r1},{r2}) "
            f"(k,l)=({k},{l}) "
            f"K={K:,}"
        )

        for S in S_values:

            case = test_S(
                scale=scale_lo,
                p=p,
                q=q,
                r1=r1,
                r2=r2,
                K=K,
                S=S,
            )

            all_cases.append(
                case
            )

            if (
                SHOW_GOOD_ONLY
                and case.good_pairs == 0
            ):
                continue

            delta = (
                "n/a"
                if case.min_delta is None
                else str(case.min_delta)
            )

            print(
                f"    S={S:7d} "
                f"x={case.x:8d} "
                f"status={case.factor_status:<7} "
                f"pairs={case.aux_pairs:4d} "
                f"good={case.good_pairs:2d} "
                f"minDelta={delta:<12} "
                f"time={case.factor_time:.4f}s"
            )

            if case.good_pairs:

                for (
                    px,
                    qx,
                    kx,
                    lx,
                ) in case.good:

                    print(
                        f"        GOOD "
                        f"(px,qx)=("
                        f"{px:,},{qx:,}) "
                        f"(kx,lx)=("
                        f"{kx},{lx})"
                    )

        print()

    return all_cases


# ==========================================================================================
# MAIN
# ==========================================================================================

def run() -> None:

    total_start = time.perf_counter()

    rng = random.Random(
        SEED
    )

    S_values = build_s_values()

    print("=" * 100)
    print(
        "DETERMINISTIC S-SWEEP / n+x -> Kx SCALING EXPERIMENT"
    )
    print("=" * 100)

    print(
        f"anchors / scale           = "
        f"{ANCHORS_PER_SCALE}"
    )

    print(
        f"S count                   = "
        f"{len(S_values)}"
    )

    print(
        f"S maximum                 = "
        f"{max(S_values):,}"
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
    # MODULUS POOL
    # --------------------------------------------------------------------------------------

    print("=" * 100)
    print(
        "BUILDING MODULUS POOL"
    )
    print("=" * 100)

    modulus_pairs = build_close_modulus_pairs()

    print(
        f"close modulus pairs       = "
        f"{len(modulus_pairs):,}"
    )

    print()

    # --------------------------------------------------------------------------------------
    # S VALUES
    # --------------------------------------------------------------------------------------

    print("=" * 100)
    print(
        "DETERMINISTIC S VALUES"
    )
    print("=" * 100)

    print(
        "S values are generated from small prime products."
    )

    print(
        f"first S values            = "
        f"{S_values[:20]}"
    )

    print(
        f"last S values             = "
        f"{S_values[-20:]}"
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

    all_cases = []

    for lo, hi in scales:

        cases = run_scale(
            lo,
            hi,
            rng,
            modulus_pairs,
            S_values,
        )

        all_cases.extend(
            cases
        )

    # --------------------------------------------------------------------------------------
    # GOOD EVENTS
    # --------------------------------------------------------------------------------------

    print("=" * 100)
    print(
        "GOOD Kx EVENTS"
    )
    print("=" * 100)

    good = [
        c
        for c in all_cases
        if c.good_pairs > 0
    ]

    if not good:

        print(
            "no Kx == K event was encountered."
        )

    else:

        for c in good:

            print(
                f"scale={c.scale:.0e} "
                f"n={c.n:,} "
                f"S={c.S:,} "
                f"x={c.x:,}"
            )

            print(
                f"    original "
                f"(k,l)=({c.k},{c.l}) "
                f"K={c.K:,}"
            )

            print(
                f"    auxiliary n+x={c.nx:,}"
            )

            print(
                f"    good pairs={c.good_pairs}"
            )

            for (
                px,
                qx,
                kx,
                lx,
            ) in c.good:

                print(
                    f"    "
                    f"(px,qx)=("
                    f"{px:,},{qx:,}) "
                    f"(kx,lx)=("
                    f"{kx},{lx})"
                )

    # --------------------------------------------------------------------------------------
    # ALGEBRAIC VALIDATION
    # --------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print(
        "ALGEBRAIC VALIDATION"
    )
    print("=" * 100)

    failures = 0

    for c in all_cases:

        if c.factor_status != "ok":
            continue

        R = c.r1 * c.r2

        T = c.n // R
        Tx = c.nx // R

        E = T - c.K

        factors = sp.factorint(
            c.nx,
            use_trial=True,
            use_rho=True,
        )

        for px, qx in divisor_pairs(
            factors
        ):

            kx = px // c.r1
            lx = qx // c.r2

            Kx = kx * lx
            Ex = Tx - Kx

            if Kx - c.K != (
                (Tx - T)
                - (Ex - E)
            ):

                failures += 1
                break

    print(
        f"trajectory identity failures = "
        f"{failures}"
    )

    # --------------------------------------------------------------------------------------
    # FINAL METRICS
    # --------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print(
        "FINAL COUNTS"
    )
    print("=" * 100)

    print(
        f"total auxiliary tests       = "
        f"{len(all_cases):,}"
    )

    print(
        f"factorization timeouts     = "
        f"{sum(c.factor_status == 'timeout' for c in all_cases):,}"
    )

    print(
        f"successful factorizations  = "
        f"{sum(c.factor_status == 'ok' for c in all_cases):,}"
    )

    print(
        f"Kx == K cases              = "
        f"{len(good):,}"
    )

    print()

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
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":

    mp.freeze_support()

    run()
