#!/usr/bin/env python3

"""
====================================================================================================
START EXPERIMENT 75
AUXILIARY RESIDUE TRAJECTORY / (ax,bx) -> (a,b) EXPERIMENT
====================================================================================================

Purpose
-------
Investigate whether the exact auxiliary residue coordinates

    ax = px - kx*r1
    bx = qx - lx*r2

contain a useful relationship to the original

    a = p - k*r1
    b = q - l*r2.

The experiment keeps the complete auxiliary tuple:

    (px,qx,kx,lx,ax,bx,Kx,Tx,Ex)

and compares it against the original:

    (p,q,k,l,a,b,K,T,E).

Important:
------------
Kx, Tx and Ex are NOT assumed to determine a,b.

This experiment directly investigates the residue layer.

No data analysis is performed beyond the requested experiment.
====================================================================================================
"""

from __future__ import annotations

import math
import random
import statistics
import time
from dataclasses import dataclass

from sympy import factorint, integer_nthroot, isprime


# ================================================================================================
# CONFIGURATION
# ================================================================================================

SEED = 1_511_464_998

SCALES = (
    10**9,
    10**12,
    10**16,
)

ANCHORS_PER_SCALE = 4

K_TARGET = 1000

R_OFFSETS = (
    0.90,
    0.95,
    1.00,
    1.05,
    1.10,
)

S_VALUES = (
    30,
    210,
    2310,
    30030,
)

# Global modulus-prime pool.
# This must cover sqrt(n / K_TARGET) at 1e16:
#
# sqrt(1e16 / 1000) ~= 3.16e6
#
# We therefore go comfortably beyond that.
MODULUS_MIN = 2
MODULUS_MAX = 6_000_000


# ================================================================================================
# DATA TYPES
# ================================================================================================

@dataclass(frozen=True)
class Anchor:
    n: int
    p: int
    q: int


@dataclass(frozen=True)
class AuxiliaryRecord:
    S: int
    x: int
    n_x: int

    px: int
    qx: int

    kx: int
    lx: int

    ax: int
    bx: int

    Kx: int
    Tx: int
    Ex: int

    same_K: bool
    same_cell: bool


# ================================================================================================
# PRIME HELPERS
# ================================================================================================

def next_prime(x: int) -> int:
    x = max(2, int(x))

    if x == 2:
        return 2

    if x % 2 == 0:
        x += 1

    while not bool(isprime(x)):
        x += 2

    return x


def generate_prime_pool(lo: int, hi: int) -> list[int]:
    primes = []

    if lo <= 2 <= hi:
        primes.append(2)

    start = max(3, lo)
    if start % 2 == 0:
        start += 1

    for x in range(start, hi + 1, 2):
        if bool(isprime(x)):
            primes.append(x)

    return primes


def closest_prime(primes: list[int], target: float) -> int:
    """
    Binary-search nearest prime instead of scanning the whole pool.
    """

    import bisect

    i = bisect.bisect_left(primes, target)

    if i <= 0:
        return primes[0]

    if i >= len(primes):
        return primes[-1]

    a = primes[i - 1]
    b = primes[i]

    if abs(a - target) <= abs(b - target):
        return a

    return b


def choose_balanced_prime_pair(
    target_R: int,
    prime_pool: list[int],
) -> tuple[int, int]:
    """
    Find a balanced prime pair r1,r2 such that r1*r2 is close to target_R.

    The previous experiment created a narrow local pool around sqrt(R).
    That could incorrectly produce:

        "Insufficient modulus primes near sqrt(R)"

    even though valid global pairs existed.

    This routine searches the GLOBAL prime pool and therefore avoids
    that failure mode.
    """

    if not prime_pool:
        raise RuntimeError("empty modulus prime pool")

    sqrt_R = math.sqrt(target_R)

    # First candidates are around sqrt(R), which favors balanced pairs.
    # We take a reasonably large neighborhood from the GLOBAL pool.
    import bisect

    center = bisect.bisect_left(prime_pool, sqrt_R)

    candidate_count = 250

    lo_index = max(0, center - candidate_count)
    hi_index = min(len(prime_pool), center + candidate_count)

    candidates = prime_pool[lo_index:hi_index]

    best = None

    for r1 in candidates:

        desired_r2 = target_R / r1

        r2 = closest_prime(prime_pool, desired_r2)

        if r1 <= 0 or r2 <= 0:
            continue

        R = r1 * r2

        relative_error = abs(R - target_R) / target_R
        balance_penalty = abs(math.log(r1 / r2))

        # Strong preference for product accuracy, weak preference for balance.
        score = relative_error + 0.01 * balance_penalty

        candidate = (
            score,
            abs(R - target_R),
            balance_penalty,
            r1,
            r2,
        )

        if best is None or candidate < best:
            best = candidate

    if best is None:
        raise RuntimeError(
            f"Could not construct modulus pair near R={target_R}"
        )

    _, _, _, r1, r2 = best

    return int(r1), int(r2)


# ================================================================================================
# ANCHOR GENERATION
# ================================================================================================

def random_prime_in_range(
    rng: random.Random,
    lo: int,
    hi: int,
) -> int:

    if lo > hi:
        raise ValueError(
            f"invalid prime range {lo}..{hi}"
        )

    candidate = rng.randint(lo, hi)
    return next_prime(candidate)


def generate_anchor(
    scale: int,
    rng: random.Random,
) -> Anchor:

    root, _ = integer_nthroot(scale, 2)
    root = int(root)

    lo = max(10_000, int(root * 0.70))
    hi = int(root * 1.30)

    while True:

        p = random_prime_in_range(rng, lo, hi)
        q = random_prime_in_range(rng, lo, hi)

        if p == q:
            continue

        p, q = sorted((p, q))
        n = p * q

        if 0.45 * scale <= n <= 1.70 * scale:
            return Anchor(
                n=n,
                p=p,
                q=q,
            )


# ================================================================================================
# QUOTIENT / RESIDUE DECOMPOSITION
# ================================================================================================

def decompose(
    value: int,
    modulus: int,
) -> tuple[int, int]:

    quotient, remainder = divmod(value, modulus)
    return quotient, remainder


def original_structure(
    n: int,
    p: int,
    q: int,
    r1: int,
    r2: int,
) -> dict[str, int]:

    k, a = decompose(p, r1)
    l, b = decompose(q, r2)

    R = r1 * r2
    K = k * l
    T = n // R
    E = T - K

    return {
        "n": n,
        "p": p,
        "q": q,

        "r1": r1,
        "r2": r2,
        "R": R,

        "k": k,
        "l": l,

        "a": a,
        "b": b,

        "K": K,
        "T": T,
        "E": E,
    }


# ================================================================================================
# FACTORIZATION / DIVISORS
# ================================================================================================

def factor_exact(n: int) -> dict[int, int]:
    return {
        int(p): int(e)
        for p, e in factorint(n).items()
    }


def divisor_pairs_from_factorization(
    n: int,
    factors: dict[int, int],
) -> list[tuple[int, int]]:

    divisors = [1]

    for prime, exponent in factors.items():

        old = list(divisors)
        multiplier = 1

        for _ in range(exponent):

            multiplier *= prime

            for d in old:
                divisors.append(d * multiplier)

    pairs = []

    root = math.isqrt(n)

    for d in sorted(set(divisors)):

        if d > root:
            break

        if n % d == 0:
            pairs.append(
                (d, n // d)
            )

    return pairs


# ================================================================================================
# AUXILIARY TRAJECTORY
# ================================================================================================

def build_auxiliary_records(
    original: dict[str, int],
) -> list[AuxiliaryRecord]:

    n = original["n"]

    r1 = original["r1"]
    r2 = original["r2"]

    k = original["k"]
    l = original["l"]

    K = original["K"]
    R = original["R"]

    records: list[AuxiliaryRecord] = []

    for S in S_VALUES:

        # Deterministic construction.
        x = (-n) % S

        # x=0 would reproduce n and is not an auxiliary displacement.
        if x == 0:
            continue

        n_x = n + x

        factors = factor_exact(n_x)
        pairs = divisor_pairs_from_factorization(
            n_x,
            factors,
        )

        Tx = n_x // R

        for px, qx in pairs:

            kx, ax = decompose(px, r1)
            lx, bx = decompose(qx, r2)

            Kx = kx * lx
            Ex = Tx - Kx

            records.append(
                AuxiliaryRecord(
                    S=S,
                    x=x,
                    n_x=n_x,

                    px=px,
                    qx=qx,

                    kx=kx,
                    lx=lx,

                    ax=ax,
                    bx=bx,

                    Kx=Kx,
                    Tx=Tx,
                    Ex=Ex,

                    same_K=(Kx == K),
                    same_cell=(
                        kx == k and
                        lx == l
                    ),
                )
            )

    return records


# ================================================================================================
# RESIDUE ANALYSIS
# ================================================================================================

def residue_metrics(
    original: dict[str, int],
    records: list[AuxiliaryRecord],
) -> dict[str, float | int]:

    a = original["a"]
    b = original["b"]

    if not records:
        return {
            "records": 0,
            "same_a": 0,
            "same_b": 0,
            "same_ab": 0,
            "near_a": 0,
            "near_b": 0,
            "near_ab": 0,
            "mean_da": 0.0,
            "mean_db": 0.0,
            "mean_sum": 0.0,
            "mean_diff": 0.0,
        }

    da_values = [
        abs(rec.ax - a)
        for rec in records
    ]

    db_values = [
        abs(rec.bx - b)
        for rec in records
    ]

    sum_values = [
        abs(
            (rec.ax + rec.bx) -
            (a + b)
        )
        for rec in records
    ]

    diff_values = [
        abs(
            (rec.ax - rec.bx) -
            (a - b)
        )
        for rec in records
    ]

    r1 = original["r1"]
    r2 = original["r2"]

    a_threshold = max(1, r1 // 100)
    b_threshold = max(1, r2 // 100)

    return {
        "records": len(records),

        "same_a": sum(
            rec.ax == a
            for rec in records
        ),

        "same_b": sum(
            rec.bx == b
            for rec in records
        ),

        "same_ab": sum(
            rec.ax == a and
            rec.bx == b
            for rec in records
        ),

        "near_a": sum(
            d <= a_threshold
            for d in da_values
        ),

        "near_b": sum(
            d <= b_threshold
            for d in db_values
        ),

        "near_ab": sum(
            abs(rec.ax - a) <= a_threshold and
            abs(rec.bx - b) <= b_threshold
            for rec in records
        ),

        "mean_da": statistics.mean(da_values),
        "mean_db": statistics.mean(db_values),

        "mean_sum": statistics.mean(sum_values),
        "mean_diff": statistics.mean(diff_values),
    }


# ================================================================================================
# IDENTITY CHECK
# ================================================================================================

def verify_difference_identity(
    original: dict[str, int],
    rec: AuxiliaryRecord,
) -> bool:

    r1 = original["r1"]
    r2 = original["r2"]

    k = original["k"]
    l = original["l"]

    a = original["a"]
    b = original["b"]

    K = original["K"]

    # Exact expansion:
    #
    # n+x - n
    #
    # = R(Kx-K)
    # + r1*(kx*bx - k*b)
    # + r2*(lx*ax - l*a)
    # + (ax*bx - a*b)

    rhs = (
        r1 * r2 * (rec.Kx - K)
        + r1 * (
            rec.kx * rec.bx -
            k * b
        )
        + r2 * (
            rec.lx * rec.ax -
            l * a
        )
        + (
            rec.ax * rec.bx -
            a * b
        )
    )

    return rhs == rec.x


# ================================================================================================
# BEST MATCHES
# ================================================================================================

def closest_residue_records(
    original: dict[str, int],
    records: list[AuxiliaryRecord],
    count: int = 8,
) -> list[AuxiliaryRecord]:

    a = original["a"]
    b = original["b"]

    ranked = sorted(
        records,
        key=lambda rec: (
            abs(rec.ax - a) +
            abs(rec.bx - b),
            abs(rec.ax - a),
            abs(rec.bx - b),
        ),
    )

    return ranked[:count]


# ================================================================================================
# SCALE RUNNER
# ================================================================================================

def run_scale(
    scale: int,
    rng: random.Random,
    prime_pool: list[int],
) -> dict[str, object]:

    print("=" * 100)
    print(f"SCALE {scale:.0e}")
    print("=" * 100)

    anchors = []

    for i in range(ANCHORS_PER_SCALE):

        anchor = generate_anchor(
            scale,
            rng,
        )

        anchors.append(anchor)

        print(
            f"anchor {i+1}/{ANCHORS_PER_SCALE} "
            f"n={anchor.n:,}"
        )

    print()

    total_records = 0
    total_K_matches = 0
    total_same_a = 0
    total_same_b = 0
    total_same_ab = 0

    per_scale = []

    for ai, anchor in enumerate(anchors, 1):

        print("-" * 100)
        print(
            f"anchor {ai}/{len(anchors)} "
            f"n={anchor.n:,}"
        )

        # Target R:
        #
        #     R ~= n/K_TARGET
        #
        R_target_base = max(
            1,
            anchor.n // K_TARGET,
        )

        for offset in R_OFFSETS:

            target_R = int(
                R_target_base * offset
            )

            r1, r2 = choose_balanced_prime_pair(
                target_R,
                prime_pool,
            )

            original = original_structure(
                anchor.n,
                anchor.p,
                anchor.q,
                r1,
                r2,
            )

            start_case = time.perf_counter()

            records = build_auxiliary_records(
                original,
            )

            case_time = time.perf_counter() - start_case

            metrics = residue_metrics(
                original,
                records,
            )

            identity_failures = sum(
                not verify_difference_identity(
                    original,
                    rec,
                )
                for rec in records
            )

            K_matches = sum(
                rec.same_K
                for rec in records
            )

            total_records += len(records)
            total_K_matches += K_matches
            total_same_a += int(metrics["same_a"])
            total_same_b += int(metrics["same_b"])
            total_same_ab += int(metrics["same_ab"])

            per_scale.append(
                {
                    "n": anchor.n,
                    "R": original["R"],
                    "r1": r1,
                    "r2": r2,
                    "K": original["K"],
                    "T": original["T"],
                    "E": original["E"],
                    "a": original["a"],
                    "b": original["b"],
                    "records": len(records),
                    "K_matches": K_matches,
                    "same_a": metrics["same_a"],
                    "same_b": metrics["same_b"],
                    "same_ab": metrics["same_ab"],
                }
            )

            print(
                f"    R={original['R']:,} "
                f"(r1,r2)=({r1:,},{r2:,})"
            )

            print(
                f"        "
                f"TRUE "
                f"(k,l)=({original['k']},{original['l']}) "
                f"K={original['K']:,} "
                f"T={original['T']:,} "
                f"E={original['E']:,}"
            )

            print(
                f"        "
                f"TRUE residues "
                f"a={original['a']:,} "
                f"b={original['b']:,}"
            )

            print(
                f"        "
                f"auxiliary records={len(records):,} "
                f"Kx=K={K_matches:,}"
            )

            print(
                f"        "
                f"exact ax=a={metrics['same_a']:,} "
                f"bx=b={metrics['same_b']:,} "
                f"(ax,bx)=(a,b)={metrics['same_ab']:,}"
            )

            print(
                f"        "
                f"within 1%: "
                f"ax/a={metrics['near_a']:,} "
                f"bx/b={metrics['near_b']:,} "
                f"both={metrics['near_ab']:,}"
            )

            print(
                f"        "
                f"mean |ax-a|={metrics['mean_da']:.3f} "
                f"mean |bx-b|={metrics['mean_db']:.3f}"
            )

            print(
                f"        "
                f"mean sum-error={metrics['mean_sum']:.3f} "
                f"mean diff-error={metrics['mean_diff']:.3f}"
            )

            print(
                f"        "
                f"difference identity failures="
                f"{identity_failures}"
            )

            best = closest_residue_records(
                original,
                records,
                count=6,
            )

            if best:

                print(
                    "        closest residue trajectories:"
                )

                for rec in best:

                    print(
                        f"            "
                        f"S={rec.S:,} "
                        f"x={rec.x:,} "
                        f"(px,qx)=({rec.px:,},{rec.qx:,}) "
                        f"(kx,lx)=({rec.kx},{rec.lx}) "
                        f"(ax,bx)=({rec.ax},{rec.bx}) "
                        f"Kx={rec.Kx:,} "
                        f"Tx={rec.Tx:,} "
                        f"Ex={rec.Ex:,} "
                        f"|ax-a|={abs(rec.ax-original['a']):,} "
                        f"|bx-b|={abs(rec.bx-original['b']):,} "
                        f"Kx=K={rec.same_K}"
                    )

            print(
                f"        time={case_time:.4f}s"
            )

    return {
        "scale": scale,
        "records": total_records,
        "K_matches": total_K_matches,
        "same_a": total_same_a,
        "same_b": total_same_b,
        "same_ab": total_same_ab,
        "cases": per_scale,
    }


# ================================================================================================
# MAIN
# ================================================================================================

def main() -> None:

    print("=" * 100)
    print("START EXPERIMENT 75")
    print("AUXILIARY RESIDUE TRAJECTORY / (ax,bx) -> (a,b) EXPERIMENT")
    print("=" * 100)
    print()

    print("configuration")
    print(
        f"    scales                    = "
        f"{[f'{s:.0e}' for s in SCALES]}"
    )
    print(
        f"    anchors / scale           = "
        f"{ANCHORS_PER_SCALE}"
    )
    print(
        f"    K target                  = "
        f"{K_TARGET}"
    )
    print(
        f"    R offsets                 = "
        f"{R_OFFSETS}"
    )
    print(
        f"    S values                  = "
        f"{S_VALUES}"
    )
    print(
        f"    modulus range             = "
        f"{MODULUS_MIN:,} .. {MODULUS_MAX:,}"
    )
    print(
        f"    seed                      = "
        f"{SEED:,}"
    )

    print()

    rng = random.Random(SEED)

    print("=" * 100)
    print("BUILDING GLOBAL MODULUS PRIME POOL")
    print("=" * 100)

    pool_start = time.perf_counter()

    prime_pool = generate_prime_pool(
        MODULUS_MIN,
        MODULUS_MAX,
    )

    pool_time = time.perf_counter() - pool_start

    print(
        f"    modulus primes            = "
        f"{len(prime_pool):,}"
    )

    print(
        f"    build time                = "
        f"{pool_time:.4f}s"
    )

    print()

    total_start = time.perf_counter()

    scale_results = []

    for scale in SCALES:

        result = run_scale(
            scale,
            rng,
            prime_pool,
        )

        scale_results.append(result)

        print()

    total_time = time.perf_counter() - total_start

    # ============================================================================================
    # CROSS-SCALE SUMMARY
    # ============================================================================================

    print("=" * 100)
    print("CROSS-SCALE SUMMARY")
    print("=" * 100)

    print(
        "scale       auxRecords      Kx=K       ax=a       bx=b    "
        "(ax,bx)=(a,b)"
    )

    print("-" * 100)

    for result in scale_results:

        print(
            f"{result['scale']:.0e}"
            f"{result['records']:>16,}"
            f"{result['K_matches']:>11,}"
            f"{result['same_a']:>11,}"
            f"{result['same_b']:>11,}"
            f"{result['same_ab']:>18,}"
        )

    print()

    print("=" * 100)
    print("INTERPRETATION TARGET")
    print("=" * 100)

    print(
        "The primary question is whether the auxiliary residue"
    )
    print(
        "coordinates (ax,bx) show a deterministic relation to (a,b)."
    )

    print()

    print(
        "The experiment distinguishes:"
    )

    print(
        "    quotient layer:"
    )
    print(
        "        Kx, Tx, Ex, kx, lx"
    )

    print()

    print(
        "    residue layer:"
    )
    print(
        "        ax, bx"
    )

    print()

    print(
        "The exact identity checked is:"
    )

    print(
        "    x = R*(Kx-K)"
    )

    print(
        "        + r1*(kx*bx-k*b)"
    )

    print(
        "        + r2*(lx*ax-l*a)"
    )

    print(
        "        + (ax*bx-a*b)"
    )

    print()

    print(
        "This lets us explicitly examine whether information from"
    )

    print(
        "the auxiliary residue layer can constrain the original"
    )

    print(
        "(a,b), rather than merely reproducing K."
    )

    print()

    print("=" * 100)
    print("TIMING")
    print("=" * 100)

    print(
        f"    prime pool build          = "
        f"{pool_time:.4f}s"
    )

    print(
        f"    total experiment runtime  = "
        f"{total_time:.4f}s"
    )

    print()

    print("=" * 100)
    print("FINISHED EXPERIMENT 75")
    print("=" * 100)


if __name__ == "__main__":
    main()