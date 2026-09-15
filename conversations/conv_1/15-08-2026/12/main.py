#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 70
KAPPA-GUIDED QUADRATIC-RELATION GENERATION
DIVISIBILITY / SMOOTHNESS ENRICHMENT TEST
UNRESTRICTED VS EQUAL-DENSITY RANDOM VS KAPPA FILTERS
NO CSV OUTPUT
NO SKLEARN
==============================================================================

Purpose
-------
Test whether the accumulated Kappa/cyclotomic information can improve
quadratic-sieve-style relation generation.

For each semiprime n=p*q define

    Q(x) = x^2 - 4n.

At the unknown true factor sum s=p+q,

    Q(s) = (p-q)^2,

so Q(s) is a square.

The experiment does NOT use s, p, q, or p-q to select candidates.
They are used only for post-hoc diagnostics.

Compared methods
----------------
1. SEQUENTIAL_BASELINE
   First EVAL_BUDGET candidates from the raw x interval.

2. RANDOM_CONTROL
   EVAL_BUDGET uniformly sampled candidates from the same raw interval.

3. KAPPA_L7
   Accept x when Q(x) is a quadratic residue OR zero modulo ell=7.

4. KAPPA_C3
   Accept x when Q(x) is QR-or-zero modulo 7,13,19.

5. KAPPA_C5
   Accept x when Q(x) is QR-or-zero modulo 7,13,19,31,37.

6. KAPPA_C7
   Accept x when Q(x) is QR-or-zero modulo
   7,13,19,31,37,61,67.

For every Kappa filter, the number of selected candidates is exactly the
same EVAL_BUDGET as the controls. Thus the key comparison is:

    relation yield among equally many Q(x) evaluations.

The Kappa scan overhead is also measured separately.

A "usable relation" here means either:

    FULL_SMOOTH:
        Q(x) is completely B-smooth,

or

    ONE_LARGE_PRIME:
        after removing all factors <= B, the remaining cofactor is
        a prime <= LARGE_PRIME.

This is a relation-generation diagnostic, not a complete QS implementation.
It counts relation candidates but does not attempt linear-algebra
independence.
==============================================================================
"""

from __future__ import annotations

import argparse
import math
import random
import statistics
import time
from dataclasses import dataclass
from typing import Callable, Iterable


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

DEFAULT_TARGETS = 120
DEFAULT_EVALS = 750
DEFAULT_POOL = 120_000
DEFAULT_FACTOR_BASE = 1000
DEFAULT_LARGE_PRIME = 1_000_000
DEFAULT_SEED = 70070

PRIME_LOW = 2_000_000
PRIME_HIGH = 4_200_000

CYCLOTOMIC = [7, 13, 19, 31, 37, 61, 67, 79, 127, 307, 331, 631, 1723]

FILTERS = {
    "KAPPA_L7": [7],
    "KAPPA_C3": [7, 13, 19],
    "KAPPA_C5": [7, 13, 19, 31, 37],
    "KAPPA_C7": [7, 13, 19, 31, 37, 61, 67],
}


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Target:
    index: int
    p: int
    q: int
    n: int
    s: int


@dataclass
class RelationStats:
    evaluations: int = 0
    full_smooth: int = 0
    one_large_prime: int = 0
    usable: int = 0
    total_cofactor: int = 0
    max_cofactor: int = 0
    min_cofactor: int | None = None

    @property
    def full_rate(self) -> float:
        return self.full_smooth / self.evaluations if self.evaluations else 0.0

    @property
    def one_lp_rate(self) -> float:
        return self.one_large_prime / self.evaluations if self.evaluations else 0.0

    @property
    def usable_rate(self) -> float:
        return self.usable / self.evaluations if self.evaluations else 0.0

    @property
    def mean_cofactor(self) -> float:
        return self.total_cofactor / self.evaluations if self.evaluations else 0.0


@dataclass
class MethodResult:
    stats: RelationStats
    raw_candidates_examined: int
    scan_seconds: float
    relation_seconds: float
    total_seconds: float
    selected_x: list[int]

    @property
    def total_candidate_work(self) -> int:
        return self.raw_candidates_examined

    @property
    def usable_per_second(self) -> float:
        return self.stats.usable / self.total_seconds if self.total_seconds else 0.0

    @property
    def usable_per_raw_candidate(self) -> float:
        return (
            self.stats.usable / self.raw_candidates_examined
            if self.raw_candidates_examined
            else 0.0
        )


# ---------------------------------------------------------------------------
# Prime generation
# ---------------------------------------------------------------------------

def sieve_primes(lo: int, hi: int) -> list[int]:
    """Return primes in [lo, hi]."""
    if hi < 2 or hi < lo:
        return []

    sieve = bytearray(b"\x01") * (hi + 1)
    sieve[:2] = b"\x00\x00"

    limit = math.isqrt(hi)
    for p in range(2, limit + 1):
        if sieve[p]:
            start = p * p
            sieve[start : hi + 1 : p] = b"\x00" * (
                ((hi - start) // p) + 1
            )

    return [x for x in range(lo, hi + 1) if sieve[x]]


# ---------------------------------------------------------------------------
# Legendre / QR tables
# ---------------------------------------------------------------------------

def qr_table(ell: int) -> bytearray:
    """
    table[r] = 1 iff r is a quadratic residue modulo ell, including 0.

    We only need QR-or-zero, so this is deliberately binary.
    """
    table = bytearray(ell)
    for x in range(ell):
        table[(x * x) % ell] = 1
    return table


def build_qr_tables(moduli: Iterable[int]) -> dict[int, bytearray]:
    return {ell: qr_table(ell) for ell in moduli}


def kappa_accept(
    q_value: int,
    moduli: list[int],
    tables: dict[int, bytearray],
) -> bool:
    """
    Accept iff Q(x) is QR-or-zero for every selected modulus.
    """
    for ell in moduli:
        if not tables[ell][q_value % ell]:
            return False
    return True


# ---------------------------------------------------------------------------
# Deterministic primality for the small cofactors encountered here
# ---------------------------------------------------------------------------

def is_prime_64(n: int) -> bool:
    """Deterministic Miller-Rabin for unsigned 64-bit integers."""
    if n < 2:
        return False

    small_primes = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for p in small_primes:
        if n == p:
            return True
        if n % p == 0:
            return False

    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1

    # Deterministic for 64-bit integers.
    bases = (
        2,
        325,
        9375,
        28178,
        450775,
        9780504,
        1795265022,
    )

    for a in bases:
        if a % n == 0:
            continue

        x = pow(a, d, n)
        if x in (1, n - 1):
            continue

        for _ in range(s - 1):
            x = (x * x) % n
            if x == n - 1:
                break
        else:
            return False

    return True


# ---------------------------------------------------------------------------
# Factor-base / smoothness test
# ---------------------------------------------------------------------------

def relation_status(
    q_value: int,
    factor_base: list[int],
    large_prime_bound: int,
) -> tuple[bool, bool, int]:
    """
    Return:

        full_smooth
        one_large_prime
        remaining_cofactor

    The cofactor is obtained after stripping all factor-base primes.
    """
    rem = q_value

    for p in factor_base:
        if p * p > rem:
            break

        while rem % p == 0:
            rem //= p

        if rem == 1:
            break

    if rem == 1:
        return True, False, 1

    if rem <= large_prime_bound and is_prime_64(rem):
        return False, True, rem

    return False, False, rem


def evaluate_candidates(
    n: int,
    xs: list[int],
    factor_base: list[int],
    large_prime_bound: int,
) -> RelationStats:
    stats = RelationStats()

    for x in xs:
        q_value = x * x - 4 * n

        if q_value <= 0:
            continue

        stats.evaluations += 1

        full, one_lp, rem = relation_status(
            q_value,
            factor_base,
            large_prime_bound,
        )

        if full:
            stats.full_smooth += 1
            stats.usable += 1
        elif one_lp:
            stats.one_large_prime += 1
            stats.usable += 1

        stats.total_cofactor += rem
        stats.max_cofactor = max(stats.max_cofactor, rem)

        if stats.min_cofactor is None:
            stats.min_cofactor = rem
        else:
            stats.min_cofactor = min(stats.min_cofactor, rem)

    return stats


# ---------------------------------------------------------------------------
# Target generation
# ---------------------------------------------------------------------------

def generate_targets(
    primes: list[int],
    count: int,
    rng: random.Random,
) -> list[Target]:
    targets: list[Target] = []
    seen: set[tuple[int, int]] = set()

    while len(targets) < count:
        p = rng.choice(primes)
        q = rng.choice(primes)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        key = (p, q)
        if key in seen:
            continue

        seen.add(key)

        targets.append(
            Target(
                index=len(targets) + 1,
                p=p,
                q=q,
                n=p * q,
                s=p + q,
            )
        )

    return targets


# ---------------------------------------------------------------------------
# Candidate interval
# ---------------------------------------------------------------------------

def candidate_interval(n: int, pool_size: int) -> tuple[int, int]:
    """
    Start just above sqrt(4n). The interval contains pool_size candidates.
    """
    x0 = math.isqrt(4 * n)
    while x0 * x0 <= 4 * n:
        x0 += 1

    return x0, x0 + pool_size - 1


# ---------------------------------------------------------------------------
# Baseline / random / Kappa selection
# ---------------------------------------------------------------------------

def sequential_candidates(
    x0: int,
    x1: int,
    budget: int,
) -> tuple[list[int], int]:
    xs = []
    examined = 0

    for x in range(x0, x1 + 1):
        examined += 1
        xs.append(x)

        if len(xs) >= budget:
            break

    return xs, examined


def random_candidates(
    x0: int,
    x1: int,
    budget: int,
    rng: random.Random,
) -> tuple[list[int], int]:
    population_size = x1 - x0 + 1

    if budget > population_size:
        raise RuntimeError("Random budget exceeds candidate-pool size.")

    offsets = rng.sample(range(population_size), budget)
    xs = [x0 + off for off in offsets]

    return xs, budget


def kappa_candidates(
    n: int,
    x0: int,
    x1: int,
    budget: int,
    moduli: list[int],
    tables: dict[int, bytearray],
) -> tuple[list[int], int, float]:
    """
    Scan the raw interval until 'budget' accepted candidates are found.

    Returns:
        selected x values
        number of raw x values examined
        acceptance rate in examined prefix
    """
    selected: list[int] = []
    examined = 0

    for x in range(x0, x1 + 1):
        examined += 1
        q_value = x * x - 4 * n

        if kappa_accept(q_value, moduli, tables):
            selected.append(x)

            if len(selected) >= budget:
                break

    acceptance = len(selected) / examined if examined else 0.0

    return selected, examined, acceptance


# ---------------------------------------------------------------------------
# Method execution
# ---------------------------------------------------------------------------

def run_method(
    name: str,
    n: int,
    xs: list[int],
    raw_examined: int,
    scan_seconds: float,
    factor_base: list[int],
    large_prime_bound: int,
) -> MethodResult:
    t0 = time.perf_counter()

    stats = evaluate_candidates(
        n,
        xs,
        factor_base,
        large_prime_bound,
    )

    relation_seconds = time.perf_counter() - t0

    return MethodResult(
        stats=stats,
        raw_candidates_examined=raw_examined,
        scan_seconds=scan_seconds,
        relation_seconds=relation_seconds,
        total_seconds=scan_seconds + relation_seconds,
        selected_x=xs,
    )


# ---------------------------------------------------------------------------
# Reporting helpers
# ---------------------------------------------------------------------------

def pct(x: float) -> str:
    return f"{100.0 * x:.4f}%"


def safe_mean(values: list[float]) -> float:
    return statistics.fmean(values) if values else float("nan")


def print_method_result(result: MethodResult) -> None:
    s = result.stats

    print(
        f"      evals={s.evaluations:5d}"
        f" full={s.full_smooth:4d}"
        f" oneLP={s.one_large_prime:4d}"
        f" usable={s.usable:4d}"
        f" usable_rate={pct(s.usable_rate)}"
        f" raw_examined={result.raw_candidates_examined:7d}"
        f" time={result.total_seconds:.4f}s"
        f" usable/s={result.usable_per_second:.3f}"
    )


# ---------------------------------------------------------------------------
# Aggregate experiment
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Kappa 70 relation-generation experiment"
    )
    parser.add_argument("--targets", type=int, default=DEFAULT_TARGETS)
    parser.add_argument("--evals", type=int, default=DEFAULT_EVALS)
    parser.add_argument("--pool", type=int, default=DEFAULT_POOL)
    parser.add_argument("--B", type=int, default=DEFAULT_FACTOR_BASE)
    parser.add_argument(
        "--large-prime",
        type=int,
        default=DEFAULT_LARGE_PRIME,
    )
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)

    args = parser.parse_args()

    if args.targets <= 0:
        raise ValueError("--targets must be positive")
    if args.evals <= 0:
        raise ValueError("--evals must be positive")
    if args.pool < args.evals:
        raise ValueError("--pool must be >= --evals")
    if args.B < 2:
        raise ValueError("--B must be >= 2")
    if args.large_prime < args.B:
        raise ValueError("--large-prime must be >= --B")

    rng = random.Random(args.seed)

    print("=" * 78)
    print("KAPPA EXPERIMENT 70")
    print("KAPPA-GUIDED QUADRATIC-RELATION GENERATION")
    print("DIVISIBILITY / SMOOTHNESS ENRICHMENT TEST")
    print("UNRESTRICTED VS EQUAL-DENSITY RANDOM VS KAPPA")
    print("NO CSV OUTPUT")
    print("NO SKLEARN")
    print("=" * 78)
    print()

    # -----------------------------------------------------------------------
    # 1. Prime population
    # -----------------------------------------------------------------------
    t0 = time.perf_counter()
    primes = sieve_primes(PRIME_LOW, PRIME_HIGH)
    generation_time = time.perf_counter() - t0

    print("1. PRIME POPULATION")
    print("-" * 78)
    print(f"prime population = {len(primes)}")
    print(f"generation time  = {generation_time:.6f}s")
    print()

    # -----------------------------------------------------------------------
    # 2. Targets
    # -----------------------------------------------------------------------
    targets = generate_targets(primes, args.targets, rng)

    print("2. TARGET SUMMARY")
    print("-" * 78)
    print(f"targets = {len(targets)}")
    for t in targets[:24]:
        print(
            f"target {t.index:4d}: "
            f"p={t.p} q={t.q} n={t.n} s={t.s}"
        )
    if len(targets) > 24:
        print("... remaining generated targets omitted")
    print()

    # -----------------------------------------------------------------------
    # 3. Kappa moduli
    # -----------------------------------------------------------------------
    print("3. KAPPA MODULUS FAMILIES")
    print("-" * 78)
    print(f"full cyclotomic family = {CYCLOTOMIC}")
    for name, moduli in FILTERS.items():
        print(f"{name:12s} = {moduli}")
    print()

    # -----------------------------------------------------------------------
    # 4. QR table construction
    # -----------------------------------------------------------------------
    used_moduli = sorted({ell for xs in FILTERS.values() for ell in xs})

    t0 = time.perf_counter()
    tables = build_qr_tables(used_moduli)
    qr_time = time.perf_counter() - t0

    print("4. QR TABLE PREPARATION")
    print("-" * 78)
    print(f"tables = {len(tables)}")
    print(f"time   = {qr_time:.6f}s")
    print()

    # -----------------------------------------------------------------------
    # 5. Factor base
    # -----------------------------------------------------------------------
    factor_base = sieve_primes(2, args.B)

    print("5. RELATION PARAMETERS")
    print("-" * 78)
    print(f"B factor base      = {args.B}")
    print(f"factor-base primes = {len(factor_base)}")
    print(f"large-prime bound  = {args.large_prime}")
    print(f"candidate pool     = {args.pool}")
    print(f"evaluation budget  = {args.evals}")
    print()

    # -----------------------------------------------------------------------
    # Aggregate storage
    # -----------------------------------------------------------------------
    method_names = [
        "SEQUENTIAL_BASELINE",
        "RANDOM_CONTROL",
        *FILTERS.keys(),
    ]

    aggregate: dict[str, dict[str, list[float]]] = {
        name: {
            "full_rate": [],
            "one_lp_rate": [],
            "usable_rate": [],
            "usable_per_second": [],
            "usable_per_raw_candidate": [],
            "raw_examined": [],
            "scan_time": [],
            "relation_time": [],
            "total_time": [],
            "acceptance": [],
            "true_in_pool": [],
            "true_passes": [],
        }
        for name in method_names
    }

    # We'll also keep per-target Kappa multipliers.
    multiplier_rows: dict[str, list[float]] = {
        name: [] for name in FILTERS
    }

    # -----------------------------------------------------------------------
    # 6. Run targets
    # -----------------------------------------------------------------------
    experiment_start = time.perf_counter()

    for target in targets:
        x0, x1 = candidate_interval(target.n, args.pool)

        true_in_pool = x0 <= target.s <= x1

        # ---------------------------------------------------------------
        # Sequential baseline
        # ---------------------------------------------------------------
        t_scan = time.perf_counter()
        seq_xs, seq_examined = sequential_candidates(
            x0, x1, args.evals
        )
        seq_scan = time.perf_counter() - t_scan

        seq_result = run_method(
            "SEQUENTIAL_BASELINE",
            target.n,
            seq_xs,
            seq_examined,
            seq_scan,
            factor_base,
            args.large_prime,
        )

        # ---------------------------------------------------------------
        # Random matched control
        # ---------------------------------------------------------------
        t_scan = time.perf_counter()
        rnd_xs, rnd_examined = random_candidates(
            x0, x1, args.evals, rng
        )
        rnd_scan = time.perf_counter() - t_scan

        rnd_result = run_method(
            "RANDOM_CONTROL",
            target.n,
            rnd_xs,
            rnd_examined,
            rnd_scan,
            factor_base,
            args.large_prime,
        )

        results: dict[str, MethodResult] = {
            "SEQUENTIAL_BASELINE": seq_result,
            "RANDOM_CONTROL": rnd_result,
        }

        # ---------------------------------------------------------------
        # Kappa filters
        # ---------------------------------------------------------------
        for name, moduli in FILTERS.items():
            t_scan = time.perf_counter()

            kappa_xs, raw_examined, acceptance = kappa_candidates(
                target.n,
                x0,
                x1,
                args.evals,
                moduli,
                tables,
            )

            scan_time = time.perf_counter() - t_scan

            # Pool exhausted: not enough accepted candidates.
            # We do not silently weaken the filter.
            if len(kappa_xs) < args.evals:
                print(
                    f"WARNING target {target.index}: "
                    f"{name} found only {len(kappa_xs)} "
                    f"of {args.evals} candidates "
                    f"within pool={args.pool}"
                )

            result = run_method(
                name,
                target.n,
                kappa_xs,
                raw_examined,
                scan_time,
                factor_base,
                args.large_prime,
            )

            result_acceptance = acceptance
            results[name] = result

            aggregate[name]["acceptance"].append(result_acceptance)

            passes_true = (
                true_in_pool
                and kappa_accept(
                    target.s * target.s - 4 * target.n,
                    moduli,
                    tables,
                )
            )

            aggregate[name]["true_passes"].append(
                1.0 if passes_true else 0.0
            )

        # ---------------------------------------------------------------
        # Aggregate all methods
        # ---------------------------------------------------------------
        for name, result in results.items():
            stats = result.stats

            aggregate[name]["full_rate"].append(stats.full_rate)
            aggregate[name]["one_lp_rate"].append(stats.one_lp_rate)
            aggregate[name]["usable_rate"].append(stats.usable_rate)
            aggregate[name]["usable_per_second"].append(
                result.usable_per_second
            )
            aggregate[name]["usable_per_raw_candidate"].append(
                result.usable_per_raw_candidate
            )
            aggregate[name]["raw_examined"].append(
                float(result.raw_candidates_examined)
            )
            aggregate[name]["scan_time"].append(result.scan_seconds)
            aggregate[name]["relation_time"].append(
                result.relation_seconds
            )
            aggregate[name]["total_time"].append(
                result.total_seconds
            )
            aggregate[name]["true_in_pool"].append(
                1.0 if true_in_pool else 0.0
            )

        # Compare each Kappa method directly to random control.
        for name in FILTERS:
            base_rate = rnd_result.stats.usable_rate
            k_rate = results[name].stats.usable_rate

            if base_rate > 0:
                multiplier_rows[name].append(
                    k_rate / base_rate
                )

        # Progress only every 10 targets.
        if target.index % 10 == 0 or target.index == len(targets):
            print(
                f"processed {target.index:4d}/{len(targets)} targets"
            )

    total_experiment_time = time.perf_counter() - experiment_start

    # -----------------------------------------------------------------------
    # 7. Aggregate results
    # -----------------------------------------------------------------------
    print()
    print("6. AGGREGATE RELATION-GENERATION RESULTS")
    print("-" * 78)

    for name in method_names:
        a = aggregate[name]

        print()
        print(name)
        print(
            f"  mean full-smooth rate = "
            f"{safe_mean(a['full_rate']):.8f}"
        )
        print(
            f"  mean one-large-prime = "
            f"{safe_mean(a['one_lp_rate']):.8f}"
        )
        print(
            f"  mean usable rate     = "
            f"{safe_mean(a['usable_rate']):.8f}"
        )
        print(
            f"  mean usable/raw-x    = "
            f"{safe_mean(a['usable_per_raw_candidate']):.8f}"
        )
        print(
            f"  mean usable/sec      = "
            f"{safe_mean(a['usable_per_second']):.4f}"
        )
        print(
            f"  mean scan time       = "
            f"{safe_mean(a['scan_time']):.6f}s"
        )
        print(
            f"  mean relation time   = "
            f"{safe_mean(a['relation_time']):.6f}s"
        )
        print(
            f"  mean total time      = "
            f"{safe_mean(a['total_time']):.6f}s"
        )

    # -----------------------------------------------------------------------
    # 8. Kappa vs random-control enrichment
    # -----------------------------------------------------------------------
    print()
    print("7. KAPPA ENRICHMENT RELATIVE TO RANDOM CONTROL")
    print("-" * 78)

    random_usable = safe_mean(
        aggregate["RANDOM_CONTROL"]["usable_rate"]
    )
    random_full = safe_mean(
        aggregate["RANDOM_CONTROL"]["full_rate"]
    )

    for name in FILTERS:
        k_usable = safe_mean(aggregate[name]["usable_rate"])
        k_full = safe_mean(aggregate[name]["full_rate"])

        usable_ratio = (
            k_usable / random_usable
            if random_usable > 0
            else float("nan")
        )

        full_ratio = (
            k_full / random_full
            if random_full > 0
            else float("nan")
        )

        raw_acceptance = safe_mean(
            aggregate[name]["acceptance"]
        )

        print(
            f"{name:12s} "
            f"acceptance={pct(raw_acceptance)} "
            f"usable_ratio={usable_ratio:.6f} "
            f"full_smooth_ratio={full_ratio:.6f}"
        )

    # -----------------------------------------------------------------------
    # 9. Target-level multipliers
    # -----------------------------------------------------------------------
    print()
    print("8. TARGET-LEVEL ENRICHMENT DISTRIBUTION")
    print("-" * 78)

    for name in FILTERS:
        values = multiplier_rows[name]

        if not values:
            print(f"{name:12s} no valid multiplier data")
            continue

        positives = sum(1 for x in values if x > 1.0)
        zeros = sum(1 for x in values if x == 0.0)

        print(
            f"{name:12s} "
            f"mean={statistics.fmean(values):.6f} "
            f"median={statistics.median(values):.6f} "
            f"min={min(values):.6f} "
            f"max={max(values):.6f} "
            f">1.0={positives}/{len(values)} "
            f"zero={zeros}"
        )

    # -----------------------------------------------------------------------
    # 10. True-sum inclusion diagnostics
    # -----------------------------------------------------------------------
    print()
    print("9. TRUE-SUM COVERAGE / KAPPA ADMISSIBILITY")
    print("-" * 78)

    true_pool_rate = safe_mean(
        aggregate["SEQUENTIAL_BASELINE"]["true_in_pool"]
    )

    print(
        f"true s inside candidate pool = "
        f"{pct(true_pool_rate)}"
    )

    for name in FILTERS:
        a = aggregate[name]
        print(
            f"{name:12s} "
            f"true_s_passes_filter="
            f"{int(sum(a['true_passes']))}/{len(a['true_passes'])} "
            f"({pct(safe_mean(a['true_passes']))})"
        )

    # -----------------------------------------------------------------------
    # 11. Divisibility sanity test
    # -----------------------------------------------------------------------
    print()
    print("10. LOCAL DIVISIBILITY SANITY")
    print("-" * 78)
    print(
        "This section asks whether Kappa filters merely force more"
        " small-prime divisibility, rather than genuinely improving"
    )
    print(
        "smoothness. The smoothness comparison above should therefore"
        " be interpreted together with these acceptance rates."
    )

    for name, moduli in FILTERS.items():
        expected = 1.0
        for ell in moduli:
            expected *= (ell + 1) / (2.0 * ell)

        observed = safe_mean(
            aggregate[name]["acceptance"]
        )

        ratio = observed / expected if expected else float("nan")

        print(
            f"{name:12s} "
            f"theory~={pct(expected)} "
            f"observed={pct(observed)} "
            f"observed/theory={ratio:.6f}"
        )

    # -----------------------------------------------------------------------
    # 12. Interpretation
    # -----------------------------------------------------------------------
    print()
    print("11. FINAL DIAGNOSTIC")
    print("-" * 78)

    random_usable_rate = safe_mean(
        aggregate["RANDOM_CONTROL"]["usable_rate"]
    )
    random_usable_sec = safe_mean(
        aggregate["RANDOM_CONTROL"]["usable_per_second"]
    )

    print(
        f"random-control usable rate = "
        f"{random_usable_rate:.8f}"
    )
    print(
        f"random-control usable/sec  = "
        f"{random_usable_sec:.6f}"
    )
    print()

    for name in FILTERS:
        usable_rate = safe_mean(
            aggregate[name]["usable_rate"]
        )
        usable_sec = safe_mean(
            aggregate[name]["usable_per_second"]
        )
        ratio = (
            usable_rate / random_usable_rate
            if random_usable_rate > 0
            else float("nan")
        )

        print(
            f"{name:12s}: "
            f"usable_rate={usable_rate:.8f} "
            f"ratio_vs_random={ratio:.6f} "
            f"usable/sec={usable_sec:.6f}"
        )

    print()
    print("Interpretation rules")
    print("--------------------")
    print(
        "NO EFFECT:"
        "  Kappa usable-rate ratios remain near 1."
    )
    print(
        "LOCAL EFFECT:"
        "  A Kappa filter raises usable relation yield,"
        " but only because it strongly enriches expected"
        " small-prime divisibility."
    )
    print(
        "INTERESTING EFFECT:"
        "  Kappa raises full/usable smoothness substantially"
        " beyond the matched random control."
    )
    print(
        "ALGORITHMIC EFFECT:"
        "  The gain survives when total CPU cost includes"
        " the Kappa scan, and usable relations/sec improves."
    )

    print()
    print(
        "IMPORTANT: this experiment does not claim that a"
        " Kappa-selected relation is independent in a QS matrix."
    )
    print(
        "It measures relation-generation enrichment only."
    )

    print()
    print("=" * 78)
    print("EXPERIMENT 70 COMPLETE")
    print(f"total runtime = {total_experiment_time:.6f}s")
    print("=" * 78)


if __name__ == "__main__":
    main()

