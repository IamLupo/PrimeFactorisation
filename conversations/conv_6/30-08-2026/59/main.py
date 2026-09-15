#!/usr/bin/env python3

"""
====================================================================================================
SCALING TEST: DETERMINISTIC n+x -> Kx -> ORIGINAL K
====================================================================================================

Tests the recursive hypothesis at:

    ~10^9
    ~10^12
    ~10^16

For each anchor:

    n = p*q

Choose close primes:

    r1 < r2

Then:

    k = floor(p/r1)
    l = floor(q/r2)
    K = k*l

For several smooth S:

    x = (-n) mod S

so:

    n+x = S*m.

SymPy factors n+x.

Every divisor pair (px,qx) of n+x is converted to:

    kx = floor(px/r1)
    lx = floor(qx/r2)

    Kx = kx*lx

A GOOD Kx occurs when:

    Kx == K.

The experiment also factors K itself.

IMPORTANT:

    The original p,q are used only to score the experiment.

Candidate generation for n+x uses only:

    n
    S
    r1
    r2

No search over x is performed.

The output reports the probability of obtaining Kx == K among
the valid auxiliary factorizations.
"""

from __future__ import annotations

import math
import random
import time
from dataclasses import dataclass

import sympy as sp


# ==========================================================================================
# CONFIGURATION
# ==========================================================================================

SEED = 1_511_464_998

# Number of anchors per scale.
ANCHORS_PER_SCALE = 4

# Smooth deterministic multipliers.
#
# They all contain factor 2.
# Since generated n is odd, x will be odd.
#
SMOOTH_S = (
    30,
    210,
    2310,
    30030,
)

MOD_MIN = 300
MOD_MAX = 3000

CLOSE_RATIO = 0.20

# For the large scale, keep factor construction near sqrt(n).
BALANCE_RATIO = 0.80


# ==========================================================================================
# DATA
# ==========================================================================================

@dataclass
class Anchor:
    p: int
    q: int

    @property
    def n(self) -> int:
        return self.p * self.q


@dataclass
class Case:
    scale: int
    anchor_id: int

    n: int
    p: int
    q: int

    r1: int
    r2: int

    k: int
    l: int
    K: int

    K_factors: dict[int, int]
    K_factor_time: float

    S: int
    x: int
    nx: int

    nx_factors: dict[int, int]
    nx_factor_time: float

    auxiliary_pairs: int
    good_pairs: int

    good_k_pairs: list[tuple[int, int]]

    min_abs_K_delta: int | None


# ==========================================================================================
# PRIME UTILITIES
# ==========================================================================================

def prime_pool(lo: int, hi: int) -> list[int]:
    return list(
        sp.primerange(lo, hi + 1)
    )


def build_close_pairs(
    primes: list[int],
    ratio: float,
) -> list[tuple[int, int]]:

    out = []

    for i, r1 in enumerate(primes):

        for r2 in primes[i + 1:]:

            if r2 - r1 <= ratio * r1:
                out.append((r1, r2))

    return out


# ==========================================================================================
# ANCHOR GENERATION
# ==========================================================================================

def generate_balanced_anchor(
    lo: int,
    hi: int,
    rng: random.Random,
) -> Anchor:

    """
    Generate two distinct primes whose product lies inside
    [lo,hi].

    The factors are kept reasonably close so the test resembles
    the previous experiments.
    """

    sqrt_lo = math.isqrt(lo)
    sqrt_hi = math.isqrt(hi)

    # Stay in the interior of the square-root interval.
    lower = max(
        3,
        int(sqrt_lo * BALANCE_RATIO),
    )

    upper = max(
        lower + 100,
        int(sqrt_hi / BALANCE_RATIO),
    )

    # Clamp to a sensible range.
    lower = max(3, lower)
    upper = max(lower + 100, upper)

    prime_list = list(
        sp.primerange(
            lower,
            upper + 1,
        )
    )

    while True:

        p = rng.choice(prime_list)
        q = rng.choice(prime_list)

        if p == q:
            continue

        n = p * q

        if lo <= n <= hi:
            if p > q:
                p, q = q, p

            return Anchor(p, q)


def build_anchors_for_scale(
    lo: int,
    hi: int,
    count: int,
    rng: random.Random,
) -> list[Anchor]:

    result = []
    seen = set()

    while len(result) < count:

        anchor = generate_balanced_anchor(
            lo,
            hi,
            rng,
        )

        if anchor.n in seen:
            continue

        seen.add(anchor.n)
        result.append(anchor)

    return result


# ==========================================================================================
# ORIGINAL K
# ==========================================================================================

def original_state(
    anchor: Anchor,
    r1: int,
    r2: int,
) -> tuple[int, int, int, int, int]:

    p = anchor.p
    q = anchor.q
    n = anchor.n

    k = p // r1
    l = q // r2

    K = k * l

    R = r1 * r2

    T = n // R
    E = T - K

    return k, l, K, T, E


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
# DIVISOR PAIRS FROM SYMPY FACTORIZATION
# ==========================================================================================

def all_divisors(
    factors: dict[int, int],
) -> list[int]:

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

    return sorted(
        set(divisors)
    )


def factor_pairs_from_factorization(
    factors: dict[int, int],
) -> list[tuple[int, int]]:

    value = 1

    for prime, exponent in factors.items():
        value *= prime ** exponent

    divisors = all_divisors(
        factors
    )

    pairs = []

    for d in divisors:

        if d * d > value:
            break

        if value % d == 0:
            pairs.append(
                (d, value // d)
            )

    return pairs


# ==========================================================================================
# FACTOR ONE NUMBER
# ==========================================================================================

def factor_with_timing(
    n: int,
) -> tuple[dict[int, int], float]:

    start = time.perf_counter()

    factors = sp.factorint(
        n,
        use_trial=True,
        use_rho=True,
        use_pm1=True,
        use_ecm=True,
    )

    elapsed = time.perf_counter() - start

    reconstructed = 1

    for p, e in factors.items():
        reconstructed *= p ** e

    if reconstructed != n:
        raise AssertionError(
            f"Factorization reconstruction failed for {n}"
        )

    return factors, elapsed


# ==========================================================================================
# AUXILIARY CASE
# ==========================================================================================

def run_auxiliary_case(
    scale: int,
    anchor_id: int,
    anchor: Anchor,
    r1: int,
    r2: int,
    K: int,
    S: int,
) -> Case:

    n = anchor.n

    # ------------------------------------------------------------------
    # Original K factorization
    # ------------------------------------------------------------------

    K_factors, K_factor_time = factor_with_timing(
        K
    )

    # ------------------------------------------------------------------
    # Deterministic n+x
    # ------------------------------------------------------------------

    x = deterministic_x(
        n,
        S,
    )

    nx = n + x

    if nx % S != 0:
        raise AssertionError(
            "n+x is not divisible by S."
        )

    # ------------------------------------------------------------------
    # Factor n+x
    # ------------------------------------------------------------------

    nx_factors, nx_factor_time = factor_with_timing(
        nx
    )

    # ------------------------------------------------------------------
    # Generate factor pairs
    # ------------------------------------------------------------------

    pairs = factor_pairs_from_factorization(
        nx_factors
    )

    auxiliary_pairs = 0
    good_pairs = 0

    good_k_pairs = []

    min_abs_K_delta = None

    for px, qx in pairs:

        # Keep both factors in the requested experimental range.
        #
        # This also prevents trivial factor pairs such as
        # (1,n+x) from dominating the count.
        if px < 10_000:
            continue

        if px > 100_000_000:
            continue

        if qx < 10_000:
            continue

        if qx > 100_000_000:
            continue

        auxiliary_pairs += 1

        kx = px // r1
        lx = qx // r2

        Kx = kx * lx

        delta = abs(
            Kx - K
        )

        if (
            min_abs_K_delta is None
            or delta < min_abs_K_delta
        ):
            min_abs_K_delta = delta

        if Kx == K:

            good_pairs += 1

            good_k_pairs.append(
                (kx, lx)
            )

    return Case(
        scale=scale,
        anchor_id=anchor_id,
        n=n,
        p=anchor.p,
        q=anchor.q,
        r1=r1,
        r2=r2,
        k=anchor.p // r1,
        l=anchor.q // r2,
        K=K,
        K_factors=K_factors,
        K_factor_time=K_factor_time,
        S=S,
        x=x,
        nx=nx,
        nx_factors=nx_factors,
        nx_factor_time=nx_factor_time,
        auxiliary_pairs=auxiliary_pairs,
        good_pairs=good_pairs,
        good_k_pairs=good_k_pairs,
        min_abs_K_delta=min_abs_K_delta,
    )


# ==========================================================================================
# PRINT FACTORIZATION COMPACTLY
# ==========================================================================================

def fmt_factors(
    factors: dict[int, int],
) -> str:

    parts = []

    for p in sorted(factors):

        e = factors[p]

        if e == 1:
            parts.append(
                str(p)
            )
        else:
            parts.append(
                f"{p}^{e}"
            )

    return "*".join(parts)


# ==========================================================================================
# MAIN
# ==========================================================================================

def run() -> None:

    total_start = time.perf_counter()

    rng = random.Random(
        SEED
    )

    scales = [
        (10**9, 10**10 - 1),
        (10**12, 10**13 - 1),
        (10**16, 10**17 - 1),
    ]

    print("=" * 100)
    print(
        "SCALING TEST: DETERMINISTIC n+x -> Kx -> ORIGINAL K"
    )
    print("=" * 100)

    print(
        f"anchors / scale           = {ANCHORS_PER_SCALE}"
    )

    print(
        f"S values                  = {SMOOTH_S}"
    )

    print(
        f"modulus range             = "
        f"{MOD_MIN} - {MOD_MAX}"
    )

    print(
        f"close ratio               = "
        f"{CLOSE_RATIO:.2%}"
    )

    print(
        f"seed                      = {SEED:,}"
    )

    print()

    # --------------------------------------------------------------------------------------
    # MODULUS POOL
    # --------------------------------------------------------------------------------------

    print("=" * 100)
    print("BUILDING MODULUS POOL")
    print("=" * 100)

    modulus_primes = prime_pool(
        MOD_MIN,
        MOD_MAX,
    )

    close_pairs = build_close_pairs(
        modulus_primes,
        CLOSE_RATIO,
    )

    print(
        f"modulus primes            = "
        f"{len(modulus_primes):,}"
    )

    print(
        f"close modulus pairs       = "
        f"{len(close_pairs):,}"
    )

    print()

    # --------------------------------------------------------------------------------------
    # MAIN SCALING TEST
    # --------------------------------------------------------------------------------------

    all_cases: list[Case] = []

    for scale_lo, scale_hi in scales:

        scale_label = f"{scale_lo:.0e}"

        print("=" * 100)
        print(
            f"BUILDING ANCHORS FOR {scale_label}"
        )
        print("=" * 100)

        anchors = build_anchors_for_scale(
            scale_lo,
            scale_hi,
            ANCHORS_PER_SCALE,
            rng,
        )

        print(
            f"anchors                   = "
            f"{len(anchors)}"
        )

        print()

        for anchor_id, anchor in enumerate(
            anchors,
            start=1,
        ):

            r1, r2 = rng.choice(
                close_pairs
            )

            k, l, K, T, E = original_state(
                anchor,
                r1,
                r2,
            )

            print(
                f"[{scale_label}] "
                f"anchor {anchor_id}/{len(anchors)} "
                f"n={anchor.n:,} "
                f"mods=({r1},{r2})"
            )

            for S in SMOOTH_S:

                case_start = time.perf_counter()

                case = run_auxiliary_case(
                    scale=scale_lo,
                    anchor_id=anchor_id,
                    anchor=anchor,
                    r1=r1,
                    r2=r2,
                    K=K,
                    S=S,
                )

                case_elapsed = (
                    time.perf_counter()
                    - case_start
                )

                all_cases.append(
                    case
                )

                good = (
                    "YES"
                    if case.good_pairs > 0
                    else "NO"
                )

                delta_text = (
                    str(case.min_abs_K_delta)
                    if case.min_abs_K_delta is not None
                    else "n/a"
                )

                print(
                    f"    S={S:5d} "
                    f"x={case.x:6d} "
                    f"auxPairs={case.auxiliary_pairs:3d} "
                    f"K={case.K} "
                    f"Kx=K={good:<3} "
                    f"good={case.good_pairs} "
                    f"min|Kx-K|={delta_text:<10} "
                    f"factor(n+x)={case.nx_factor_time:.4f}s "
                    f"factor(K)={case.K_factor_time:.4f}s "
                    f"case={case_elapsed:.4f}s"
                )

            print()

    # --------------------------------------------------------------------------------------
    # SUMMARY
    # --------------------------------------------------------------------------------------

    print("=" * 100)
    print("SUMMARY")
    print("=" * 100)

    by_scale: dict[int, list[Case]] = {}

    for case in all_cases:
        by_scale.setdefault(
            case.scale,
            []
        ).append(case)

    for scale in sorted(by_scale):

        cases = by_scale[scale]

        good_cases = sum(
            1
            for c in cases
            if c.good_pairs > 0
        )

        total_aux = sum(
            c.auxiliary_pairs
            for c in cases
        )

        good_pairs = sum(
            c.good_pairs
            for c in cases
        )

        factor_nx_time = sum(
            c.nx_factor_time
            for c in cases
        )

        factor_K_time = sum(
            c.K_factor_time
            for c in cases
        )

        print(
            f"scale={scale:.0e}"
        )

        print(
            f"    cases                    = "
            f"{len(cases)}"
        )

        print(
            f"    good Kx cases            = "
            f"{good_cases}/{len(cases)}"
        )

        print(
            f"    case probability         = "
            f"{good_cases / len(cases):.6f}"
        )

        print(
            f"    auxiliary factor pairs  = "
            f"{total_aux}"
        )

        print(
            f"    exact Kx==K pairs        = "
            f"{good_pairs}"
        )

        print(
            f"    total factor(n+x) time  = "
            f"{factor_nx_time:.4f}s"
        )

        print(
            f"    total factor(K) time    = "
            f"{factor_K_time:.4f}s"
        )

        print()

    # --------------------------------------------------------------------------------------
    # PER-S VALUE SUMMARY
    # --------------------------------------------------------------------------------------

    print("=" * 100)
    print("BY SMOOTH MULTIPLIER S")
    print("=" * 100)

    for S in SMOOTH_S:

        cases = [
            c
            for c in all_cases
            if c.S == S
        ]

        good_cases = sum(
            1
            for c in cases
            if c.good_pairs > 0
        )

        good_pairs = sum(
            c.good_pairs
            for c in cases
        )

        print(
            f"S={S:5d} "
            f"cases={len(cases):3d} "
            f"good-cases={good_cases:3d} "
            f"probability="
            f"{good_cases / len(cases):.6f} "
            f"good-pairs={good_pairs:4d}"
        )

    # --------------------------------------------------------------------------------------
    # GOOD CASES
    # --------------------------------------------------------------------------------------

    print()
    print("=" * 100)
    print("GOOD Kx CASES")
    print("=" * 100)

    good_cases = [
        c
        for c in all_cases
        if c.good_pairs > 0
    ]

    if not good_cases:

        print(
            "No exact Kx == K cases were found."
        )

    else:

        for case in good_cases:

            print(
                f"scale={case.scale:.0e} "
                f"n={case.n:,} "
                f"S={case.S} "
                f"x={case.x} "
                f"(p,q)=({case.p:,},{case.q:,})"
            )

            print(
                f"    original "
                f"(k,l)=({case.k},{case.l}) "
                f"K={case.K}"
            )

            print(
                f"    auxiliary n+x={case.nx:,}"
            )

            print(
                f"    K factorization="
                f"{fmt_factors(case.K_factors)}"
            )

            print(
                f"    good Kx pairs="
                f"{case.good_k_pairs}"
            )

            print(
                f"    n+x factorization="
                f"{fmt_factors(case.nx_factors)}"
            )

            print()

    # --------------------------------------------------------------------------------------
    # ALGEBRAIC CHECKS
    # --------------------------------------------------------------------------------------

    print("=" * 100)
    print("ALGEBRAIC CHECKS")
    print("=" * 100)

    failures = 0

    for case in all_cases:

        R = case.r1 * case.r2

        T = case.n // R
        Tx = case.nx // R

        E = T - case.K

        for pair in factor_pairs_from_factorization(
            case.nx_factors
        ):

            px, qx = pair

            kx = px // case.r1
            lx = qx // case.r2
            Kx = kx * lx

            Ex = Tx - Kx

            if Kx - case.K != (
                (Tx - T)
                - (Ex - E)
            ):
                failures += 1
                break

    print(
        f"trajectory identity failures = "
        f"{failures}"
    )

    if failures == 0:
        print(
            "all tested trajectories satisfy"
        )
        print(
            "    Kx-K=(Tx-T)-(Ex-E)"
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
        f"{elapsed:.3f} s"
    )

    print()
    print("=" * 100)
    print("EXPERIMENT COMPLETE")
    print("=" * 100)


if __name__ == "__main__":
    run()
