#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 54
HOLDOUT QUADRATIC-RESIDUE STRUCTURE / CYCLOTOMIC VS RANDOM
NO CSV OUTPUT
==============================================================================

PURPOSE
-------
Test whether unseen cyclotomic prime factors behave differently from
unrelated random primes when applied AFTER the existing base cyclotomic
quadratic-residue sieve.

Main questions:

1. Do unseen cyclotomic primes reject base survivors at a rate different
   from the generic ~1/2 quadratic-residue expectation?

2. Do they behave differently from random-prime holdouts?

3. Are holdout tests approximately independent, or correlated?

4. Does the holdout family reduce the survivor set substantially while
   preserving the true factorization?

No candidate-pair enumeration.
No CSV files.
==============================================================================

"""

from __future__ import annotations

import math
import random
import time
from typing import Iterable

import numpy as np


# ============================================================================
# 0. CONFIGURATION
# ============================================================================

SEED = 20260814

PRIME_LO = 2_000_000
PRIME_HI = 4_200_000

S_MIN = 4_000_000
S_MAX = 8_399_998

BASE_R = [
    2, 3, 5, 7, 11,
    13, 17, 19, 23, 29,
    31, 37, 41, 43, 47,
]

HOLDOUT_R_POOL = [
    53, 59, 61, 67, 71, 73, 79, 83,
    89, 97, 101, 103, 107, 109, 113,
    127, 131, 137, 139, 149, 151,
]

HOLDOUT_COUNT = 12

RANDOM_HOLDOUT_LO = 2_000
RANDOM_HOLDOUT_HI = 50_000

TARGETS = [
    (3318013, 4042603),
    (2129167, 3402323),
    (2224517, 3978749),
    (3685051, 4020281),
    (2399627, 2452649),
    (2593039, 2996527),
    (2149859, 2772097),
    (2060543, 2514401),
    (2675423, 2883973),
    (2828887, 3960137),
    (3497381, 3793241),
    (2193509, 4011353),
]


# ============================================================================
# 1. PRIME / FACTORIZATION HELPERS
# ============================================================================

def sieve_primes(limit: int) -> list[int]:
    """Return all primes <= limit."""
    if limit < 2:
        return []

    sieve = np.ones(limit + 1, dtype=np.bool_)
    sieve[:2] = False

    root = int(math.isqrt(limit))

    for p in range(2, root + 1):
        if sieve[p]:
            sieve[p * p : limit + 1 : p] = False

    return np.flatnonzero(sieve).tolist()


def factor_integer(n: int) -> list[int]:
    """Return distinct prime factors of a small integer."""
    factors: list[int] = []
    x = n

    if x % 2 == 0:
        factors.append(2)
        while x % 2 == 0:
            x //= 2

    d = 3
    while d * d <= x:
        if x % d == 0:
            factors.append(d)
            while x % d == 0:
                x //= d
        d += 2

    if x > 1:
        factors.append(x)

    return factors


def F(r: int) -> int:
    return r * r + r + 1


# ============================================================================
# 2. TARGET PREPARATION
# ============================================================================

def prepare_targets() -> list[dict[str, int]]:
    result = []

    for i, (p, q) in enumerate(TARGETS, start=1):
        n = p * q
        s = p + q
        d = abs(q - p)

        if s % 2 != 0:
            raise ValueError(f"Target {i}: p+q must be even.")

        if d * d != s * s - 4 * n:
            raise ValueError(f"Target {i}: discriminant identity failed.")

        result.append(
            {
                "id": i,
                "p": p,
                "q": q,
                "n": n,
                "s": s,
                "d": d,
            }
        )

    return result


# ============================================================================
# 3. CYCLOTOMIC FAMILIES
# ============================================================================

def build_base_cyclotomic_family() -> tuple[list[int], dict[int, list[int]]]:
    sources: dict[int, list[int]] = {}

    for r in BASE_R:
        for ell in factor_integer(F(r)):
            sources.setdefault(ell, []).append(r)

    return sorted(sources), sources


def build_cyclotomic_holdouts(
    base_primes: set[int],
) -> tuple[list[int], dict[int, list[int]], dict[int, int]]:
    """
    Find NEW prime factors from F(r)=r^2+r+1 for larger r values.
    """
    sources: dict[int, list[int]] = {}
    source_values: dict[int, int] = {}

    for r in HOLDOUT_R_POOL:
        value = F(r)
        source_values[r] = value

        for ell in factor_integer(value):
            if ell in base_primes:
                continue

            sources.setdefault(ell, []).append(r)

    holdouts = sorted(sources)[:HOLDOUT_COUNT]

    if len(holdouts) < HOLDOUT_COUNT:
        raise RuntimeError(
            f"Only discovered {len(holdouts)} new cyclotomic primes."
        )

    return holdouts, sources, source_values


# ============================================================================
# 4. RANDOM HOLDOUTS
# ============================================================================

def choose_random_holdouts(
    all_primes: list[int],
    excluded: set[int],
    count: int,
) -> list[int]:
    """
    Corrected parameter name.

    Choose 'count' deterministic random primes from the requested range,
    excluding all known cyclotomic factors.
    """
    candidates = [
        p for p in all_primes
        if RANDOM_HOLDOUT_LO <= p <= RANDOM_HOLDOUT_HI
        and p not in excluded
    ]

    if len(candidates) < count:
        raise RuntimeError(
            f"Only {len(candidates)} random holdout primes available."
        )

    rng = random.Random(SEED)
    return sorted(rng.sample(candidates, count))


# ============================================================================
# 5. QR TABLES
# ============================================================================

def build_qr_table(modulus: int) -> np.ndarray:
    """
    qr[a] = True iff a is a quadratic residue modulo modulus.
    Zero is included.
    """
    qr = np.zeros(modulus, dtype=np.bool_)

    values = np.arange(modulus, dtype=np.int64)
    residues = (values * values) % modulus
    qr[residues] = True

    return qr


# ============================================================================
# 6. EVEN-SUM DOMAIN
# ============================================================================

def all_even_sums() -> np.ndarray:
    count = (S_MAX - S_MIN) // 2 + 1

    k = np.arange(count, dtype=np.int64)
    return S_MIN + 2 * k


# ============================================================================
# 7. BASE QR FILTER
# ============================================================================

def qr_filter_sums(
    n: int,
    moduli: Iterable[int],
    qr_tables: dict[int, np.ndarray],
    sums: np.ndarray,
) -> np.ndarray:
    """
    Apply all supplied QR filters simultaneously.
    """
    disc = sums * sums - 4 * np.int64(n)

    mask = np.ones(sums.shape, dtype=np.bool_)

    for ell in moduli:
        residues = disc % ell
        mask &= qr_tables[ell][residues]

        if not np.any(mask):
            return np.empty(0, dtype=np.int64)

    return sums[mask]


def apply_holdout(
    survivors: np.ndarray,
    n: int,
    ell: int,
    qr: np.ndarray,
) -> np.ndarray:
    if len(survivors) == 0:
        return survivors

    disc = survivors * survivors - 4 * np.int64(n)
    residues = disc % ell

    return survivors[qr[residues]]


# ============================================================================
# 8. HOLDOUT STATISTICS
# ============================================================================

def independent_holdout_statistics(
    survivors: np.ndarray,
    n: int,
    moduli: list[int],
    qr_tables: dict[int, np.ndarray],
) -> list[dict[str, float]]:
    """
    Evaluate each holdout independently against the SAME survivor set.
    """
    before = len(survivors)
    result = []

    if before == 0:
        return [
            {
                "ell": ell,
                "before": 0,
                "retained": 0,
                "fraction": 0.0,
                "bits": float("inf"),
            }
            for ell in moduli
        ]

    disc = survivors * survivors - 4 * np.int64(n)

    for ell in moduli:
        residues = disc % ell
        retained = int(
            np.count_nonzero(qr_tables[ell][residues])
        )

        fraction = retained / before

        bits = (
            -math.log2(fraction)
            if fraction > 0
            else float("inf")
        )

        result.append(
            {
                "ell": ell,
                "before": before,
                "retained": retained,
                "fraction": fraction,
                "bits": bits,
            }
        )

    return result


# ============================================================================
# 9. PAIRWISE DEPENDENCE
# ============================================================================

def pairwise_holdout_correlation(
    survivors: np.ndarray,
    n: int,
    moduli: list[int],
    qr_tables: dict[int, np.ndarray],
) -> list[tuple[int, int, float, float]]:
    """
    Compare observed joint retention to independent expectation.

    ratio = observed_joint / (marginal_a * marginal_b)

    ratio ~= 1  -> approximately independent
    ratio < 1   -> negative dependence
    ratio > 1   -> positive dependence
    """
    if len(survivors) == 0:
        return []

    disc = survivors * survivors - 4 * np.int64(n)

    masks: dict[int, np.ndarray] = {}

    for ell in moduli:
        residues = disc % ell
        masks[ell] = qr_tables[ell][residues]

    result = []

    for i, a in enumerate(moduli):
        ma = masks[a]
        pa = float(np.mean(ma))

        for b in moduli[i + 1 :]:
            mb = masks[b]
            pb = float(np.mean(mb))

            observed = float(np.mean(ma & mb))
            expected = pa * pb

            ratio = (
                observed / expected
                if expected > 0
                else 0.0
            )

            result.append((a, b, observed, ratio))

    return result


# ============================================================================
# 10. EXACT RECOVERY
# ============================================================================

def exact_factor_from_sum(
    n: int,
    s: int,
) -> tuple[int, int] | None:
    disc = s * s - 4 * n

    if disc < 0:
        return None

    d = math.isqrt(disc)

    if d * d != disc:
        return None

    if (s - d) % 2 != 0:
        return None

    p = (s - d) // 2
    q = (s + d) // 2

    if p <= 0 or q <= 0:
        return None

    if p * q != n:
        return None

    return min(p, q), max(p, q)


def recover_from_survivors(
    n: int,
    survivors: np.ndarray,
) -> tuple[int, int] | None:
    for s in survivors:
        pair = exact_factor_from_sum(n, int(s))
        if pair is not None:
            return pair

    return None


# ============================================================================
# 11. MAIN
# ============================================================================

def main() -> None:
    wall_start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 54")
    print("HOLDOUT QUADRATIC-RESIDUE STRUCTURE / CYCLOTOMIC VS RANDOM")
    print("NO CSV OUTPUT")
    print("=" * 78)
    print()

    # ----------------------------------------------------------------------
    # Prime population
    # ----------------------------------------------------------------------

    t0 = time.perf_counter()

    prime_population_full = sieve_primes(PRIME_HI - 1)

    prime_population = [
        p for p in prime_population_full
        if PRIME_LO <= p < PRIME_HI
    ]

    prime_time = time.perf_counter() - t0

    print("-" * 78)
    print("1. CONFIGURATION")
    print("-" * 78)
    print(f"random seed            = {SEED}")
    print(f"targets                = {len(TARGETS)}")
    print(
        f"prime interval         = "
        f"[{PRIME_LO:,}, {PRIME_HI:,})"
    )
    print(f"base r count           = {len(BASE_R)}")
    print(f"holdout count          = {HOLDOUT_COUNT}")
    print()

    print("-" * 78)
    print("2. PRIME POPULATION")
    print("-" * 78)
    print(f"prime population       = {len(prime_population):,}")
    print(f"generation time        = {prime_time:.6f}s")
    print()

    # ----------------------------------------------------------------------
    # Targets
    # ----------------------------------------------------------------------

    targets = prepare_targets()

    print("-" * 78)
    print("3. TARGETS")
    print("-" * 78)

    for t in targets:
        print(
            f"target {t['id']:2d}: "
            f"p={t['p']} q={t['q']} "
            f"n={t['n']} s={t['s']}"
        )

    print()

    # ----------------------------------------------------------------------
    # Base family
    # ----------------------------------------------------------------------

    base_cyclo, base_sources = build_base_cyclotomic_family()

    print("-" * 78)
    print("4. BASE CYCLOTOMIC FAMILY")
    print("-" * 78)

    for ell in base_cyclo:
        print(
            f"ell={ell:5d} "
            f"source_r={base_sources[ell]}"
        )

    print()
    print(
        f"distinct base cyclotomic primes = "
        f"{len(base_cyclo)}"
    )
    print()

    # ----------------------------------------------------------------------
    # Cyclotomic holdouts
    # ----------------------------------------------------------------------

    (
        holdout_cyclo,
        holdout_sources,
        holdout_values,
    ) = build_cyclotomic_holdouts(set(base_cyclo))

    print("-" * 78)
    print("5. UNSEEN CYCLOTOMIC HOLDOUT FAMILY")
    print("-" * 78)

    for ell in holdout_cyclo:
        print(
            f"ell={ell:5d} "
            f"source_r={holdout_sources[ell]}"
        )

    print()
    print(
        f"holdout cyclotomic primes = "
        f"{len(holdout_cyclo)}"
    )
    print()

    print("holdout source values:")

    printed_r = sorted(
        {
            r
            for ell in holdout_cyclo
            for r in holdout_sources[ell]
        }
    )

    for r in printed_r:
        print(
            f"  r={r:3d} "
            f"F(r)={F(r):6d}"
        )

    print()

    # ----------------------------------------------------------------------
    # Random holdouts
    # ----------------------------------------------------------------------

    excluded = (
        set(base_cyclo)
        | set(holdout_cyclo)
    )

    # IMPORTANT:
    # The bug in the previous version came from passing
    # 'population=' instead of the function's 'all_primes=' argument.
    #
    # Here we pass the complete prime list explicitly.

    random_holdout_pool = [
        p for p in sieve_primes(RANDOM_HOLDOUT_HI)
        if RANDOM_HOLDOUT_LO <= p <= RANDOM_HOLDOUT_HI
    ]

    holdout_random = choose_random_holdouts(
        all_primes=random_holdout_pool,
        excluded=excluded,
        count=HOLDOUT_COUNT,
    )

    print("-" * 78)
    print("6. RANDOM HOLDOUT FAMILY")
    print("-" * 78)
    print(holdout_random)
    print()

    # ----------------------------------------------------------------------
    # QR tables
    # ----------------------------------------------------------------------

    all_moduli = sorted(
        set(base_cyclo)
        | set(holdout_cyclo)
        | set(holdout_random)
    )

    t0 = time.perf_counter()

    qr_tables = {
        ell: build_qr_table(ell)
        for ell in all_moduli
    }

    qr_time = time.perf_counter() - t0

    print("-" * 78)
    print("7. QR TABLE PREPARATION")
    print("-" * 78)
    print(f"total QR tables        = {len(qr_tables)}")
    print(f"preparation time       = {qr_time:.6f}s")
    print()

    # ----------------------------------------------------------------------
    # Complete even-sum domain
    # ----------------------------------------------------------------------

    sums = all_even_sums()
    full_domain_size = len(sums)

    print("-" * 78)
    print("8. BASE SIEVE")
    print("-" * 78)

    all_results = []

    for t in targets:
        t0 = time.perf_counter()

        survivors = qr_filter_sums(
            n=t["n"],
            moduli=base_cyclo,
            qr_tables=qr_tables,
            sums=sums,
        )

        elapsed = time.perf_counter() - t0

        true_survives = (
            t["s"] in set(survivors.tolist())
        )

        information = (
            math.log2(full_domain_size / len(survivors))
            if len(survivors) > 0
            else float("inf")
        )

        print(
            f"target {t['id']:2d}: "
            f"base_survivors={len(survivors):6d} "
            f"information={information:8.3f} bits "
            f"time={elapsed:.6f}s "
            f"true_survives={true_survives}"
        )

        all_results.append(
            {
                "target": t,
                "base_survivors": survivors,
            }
        )

    print()

    # ----------------------------------------------------------------------
    # Independent holdout analysis
    # ----------------------------------------------------------------------

    print("-" * 78)
    print("9. HOLDOUT MARGINAL RETENTION")
    print("-" * 78)

    aggregate_cyclo = []
    aggregate_random = []

    for result in all_results:
        t = result["target"]
        survivors = result["base_survivors"]
        n = t["n"]

        print(
            f"TARGET {t['id']:2d} "
            f"base_survivors={len(survivors)}"
        )

        print("  CYCLOTOMIC HOLDOUT")

        cstats = independent_holdout_statistics(
            survivors,
            n,
            holdout_cyclo,
            qr_tables,
        )

        for st in cstats:
            aggregate_cyclo.append(st["fraction"])

            print(
                f"    ell={st['ell']:5d} "
                f"retained={st['retained']:5d} "
                f"fraction={st['fraction']:.6f} "
                f"bits={st['bits']:.4f}"
            )

        print("  RANDOM HOLDOUT")

        rstats = independent_holdout_statistics(
            survivors,
            n,
            holdout_random,
            qr_tables,
        )

        for st in rstats:
            aggregate_random.append(st["fraction"])

            print(
                f"    ell={st['ell']:5d} "
                f"retained={st['retained']:5d} "
                f"fraction={st['fraction']:.6f} "
                f"bits={st['bits']:.4f}"
            )

        print()

    # ----------------------------------------------------------------------
    # Marginal summary
    # ----------------------------------------------------------------------

    print("-" * 78)
    print("10. HOLDOUT MARGINAL SUMMARY")
    print("-" * 78)

    cyclo_median = float(np.median(aggregate_cyclo))
    random_median = float(np.median(aggregate_random)

    )

    cyclo_mean = float(np.mean(aggregate_cyclo))
    random_mean = float(np.mean(aggregate_random))

    print(
        f"cyclotomic holdout median retention = "
        f"{cyclo_median:.6f}"
    )
    print(
        f"random holdout median retention     = "
        f"{random_median:.6f}"
    )
    print(
        f"cyclotomic holdout mean retention   = "
        f"{cyclo_mean:.6f}"
    )
    print(
        f"random holdout mean retention       = "
        f"{random_mean:.6f}"
    )
    print()
    print("generic QR expectation ≈ 0.5")
    print()

    # ----------------------------------------------------------------------
    # Joint holdout analysis
    # ----------------------------------------------------------------------

    print("-" * 78)
    print("11. JOINT HOLDOUT SIEVE")
    print("-" * 78)

    joint_summary = []

    for result in all_results:
        t = result["target"]
        n = t["n"]
        base = result["base_survivors"]

        print(f"TARGET {t['id']:2d}")

        current_c = base.copy()

        for ell in holdout_cyclo:
            before = len(current_c)

            current_c = apply_holdout(
                current_c,
                n,
                ell,
                qr_tables[ell],
            )

            after = len(current_c)

            retention = (
                after / before
                if before else 0.0
            )

            print(
                f"  C ell={ell:5d} "
                f"before={before:5d} "
                f"after={after:5d} "
                f"retained={retention:.6f}"
            )

            if after == 0:
                break

        cyclo_final = current_c

        current_r = base.copy()

        for ell in holdout_random:
            before = len(current_r)

            current_r = apply_holdout(
                current_r,
                n,
                ell,
                qr_tables[ell],
            )

            after = len(current_r)

            retention = (
                after / before
                if before else 0.0
            )

            print(
                f"  R ell={ell:5d} "
                f"before={before:5d} "
                f"after={after:5d} "
                f"retained={retention:.6f}"
            )

            if after == 0:
                break

        random_final = current_r

        true_s = t["s"]

        cyclo_correct_survival = (
            true_s in set(cyclo_final.tolist())
        )

        random_correct_survival = (
            true_s in set(random_final.tolist())
        )

        print(
            f"  cyclotomic final={len(cyclo_final):5d} "
            f"true_survives={cyclo_correct_survival}"
        )

        print(
            f"  random final    ={len(random_final):5d} "
            f"true_survives={random_correct_survival}"
        )

        joint_summary.append(
            {
                "target": t["id"],
                "base": len(base),
                "cyclo_final": len(cyclo_final),
                "random_final": len(random_final),
                "cyclo_survivor_array": cyclo_final,
                "random_survivor_array": random_final,
            }
        )

        print()

    # ----------------------------------------------------------------------
    # Pairwise dependence
    # ----------------------------------------------------------------------

    print("-" * 78)
    print("12. HOLDOUT FILTER CORRELATION")
    print("-" * 78)

    for result in all_results:
        t = result["target"]
        n = t["n"]
        survivors = result["base_survivors"]

        print(f"TARGET {t['id']:2d}")

        cyclo_corr = pairwise_holdout_correlation(
            survivors,
            n,
            holdout_cyclo,
            qr_tables,
        )

        random_corr = pairwise_holdout_correlation(
            survivors,
            n,
            holdout_random,
            qr_tables,
        )

        if cyclo_corr:
            cyclo_ratios = [x[3] for x in cyclo_corr]

            print(
                f"  cyclotomic median ratio = "
                f"{float(np.median(cyclo_ratios)):.6f}"
            )
            print(
                f"  cyclotomic mean ratio   = "
                f"{float(np.mean(cyclo_ratios)):.6f}"
            )

        if random_corr:
            random_ratios = [x[3] for x in random_corr]

            print(
                f"  random median ratio     = "
                f"{float(np.median(random_ratios)):.6f}"
            )
            print(
                f"  random mean ratio       = "
                f"{float(np.mean(random_ratios)):.6f}"
            )

        print("  independence benchmark = 1.0")
        print()

    # ----------------------------------------------------------------------
    # Exact recovery
    # ----------------------------------------------------------------------

    print("-" * 78)
    print("13. EXACT RECOVERY AFTER HOLDOUT")
    print("-" * 78)

    cyclo_recovery = 0
    random_recovery = 0

    for js in joint_summary:
        target_id = js["target"]
        t = targets[target_id - 1]

        true_pair = tuple(
            sorted((t["p"], t["q"]))
        )

        cyclo_pair = recover_from_survivors(
            t["n"],
            js["cyclo_survivor_array"],
        )

        random_pair = recover_from_survivors(
            t["n"],
            js["random_survivor_array"],
        )

        c_ok = cyclo_pair == true_pair
        r_ok = random_pair == true_pair

        cyclo_recovery += int(c_ok)
        random_recovery += int(r_ok)

        print(
            f"target {target_id:2d}: "
            f"cyclo_candidates="
            f"{len(js['cyclo_survivor_array']):5d} "
            f"pair={cyclo_pair} "
            f"correct={c_ok}"
        )

        print(
            f"             "
            f"random_candidates="
            f"{len(js['random_survivor_array']):5d} "
            f"pair={random_pair} "
            f"correct={r_ok}"
        )

    print()
    print(
        f"cyclotomic recovery = "
        f"{cyclo_recovery}/{len(targets)}"
    )
    print(
        f"random recovery     = "
        f"{random_recovery}/{len(targets)}"
    )
    print()

    # ----------------------------------------------------------------------
    # Distribution diagnostics
    # ----------------------------------------------------------------------

    print("-" * 78)
    print("14. BASE SURVIVOR DISTRIBUTION")
    print("-" * 78)

    for result in all_results:
        t = result["target"]
        survivors = result["base_survivors"]

        if len(survivors) > 1:
            gaps = np.diff(np.sort(survivors))

            print(
                f"target {t['id']:2d}: "
                f"count={len(survivors):5d} "
                f"gap_median={float(np.median(gaps)):6.1f} "
                f"gap_mean={float(np.mean(gaps)):6.1f} "
                f"gap_max={int(np.max(gaps)):7d}"
            )
        else:
            print(
                f"target {t['id']:2d}: "
                f"count={len(survivors):5d}"
            )

    print()

    # ----------------------------------------------------------------------
    # Final summary
    # ----------------------------------------------------------------------

    print("-" * 78)
    print("15. FINAL SUMMARY")
    print("-" * 78)

    base_counts = [
        len(r["base_survivors"])
        for r in all_results
    ]

    cyclo_counts = [
        x["cyclo_final"]
        for x in joint_summary
    ]

    random_counts = [
        x["random_final"]
        for x in joint_summary
    ]

    print(
        f"median base survivors          = "
        f"{float(np.median(base_counts)):.1f}"
    )

    print(
        f"median cyclotomic holdout end  = "
        f"{float(np.median(cyclo_counts)):.1f}"
    )

    print(
        f"median random holdout end      = "
        f"{float(np.median(random_counts)):.1f}"
    )

    print()
    print(
        f"cyclotomic median retention    = "
        f"{cyclo_median:.6f}"
    )

    print(
        f"random median retention        = "
        f"{random_median:.6f}"
    )

    print()
    print(
        "Interpretation:"
    )
    print(
        "  - retention ~0.5 suggests generic QR behaviour."
    )
    print(
        "  - systematic cyclotomic deviation from random suggests"
    )
    print(
        "    additional structure tied to F(r)=r^2+r+1."
    )
    print(
        "  - holdout correlations far from 1 suggest redundancy"
    )
    print(
        "    or hidden relations among the residue constraints."
    )
    print(
        "  - a large extra reduction with both holdout families"
    )
    print(
        "    means the base survivors were not yet fully characterized."
    )
    print()

    total_runtime = time.perf_counter() - wall_start

    print("=" * 78)
    print("EXPERIMENT 54 COMPLETE")
    print("=" * 78)
    print(f"total runtime = {total_runtime:.6f}s")
    print()


if __name__ == "__main__":
    main()