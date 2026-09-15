#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 71R
CACHED SIZE-MATCHED KAPPA RELATION ENRICHMENT
WITHIN-Q-SIZE / WITHIN-X CONTROL
L7 PRIMARY + C3 SECONDARY
STRICT OUT-OF-SAMPLE TARGET TEST
NO CSV OUTPUT
NO SKLEARN
==============================================================================

Purpose
-------
Experiment 70 showed apparent enrichment for Kappa-selected candidates, but
the Kappa and random populations may have had different Q(x) / x-size
distributions.

Experiment 71R asks a cleaner question:

    Among candidates of comparable Q-size, does Kappa selection itself
    increase the probability of obtaining a usable QS relation?

Primary family:
    KAPPA_L7 = [7]

Secondary family:
    KAPPA_C3 = [7, 13, 19]

The expensive factorization of Q(x) is performed exactly once per candidate
and cached. Matched controls and permutation tests operate only on cached
Boolean labels, making them cheap.

NO CSV FILES ARE PRODUCED.
"""

from __future__ import annotations

import math
import random
import statistics
import time
from dataclasses import dataclass
from typing import Iterable


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

SEED = 71071

NUM_TARGETS = 40
TRAIN_TARGETS = 30
TEST_TARGETS = 10

PRIME_MIN = 2_000_000
PRIME_MAX = 4_200_000

B = 1000
LARGE_PRIME_BOUND = 1_000_000

CANDIDATE_POOL = 30_000

# Small enough to preserve runtime while still controlling size.
SIZE_BINS = 30

# Maximum number of matched control observations per Kappa observation.
MAX_MATCHED_PER_TARGET = 1_500

# Number of cheap permutations.
PERMUTATIONS = 100

# Only these two families are tested.
KAPPA_FAMILIES = {
    "KAPPA_L7": [7],
    "KAPPA_C3": [7, 13, 19],
}


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Target:
    idx: int
    p: int
    q: int
    n: int
    s: int


@dataclass
class Candidate:
    x: int
    q_value: int
    log_size: float
    x_bin: int
    q_bin: int

    full_smooth: bool
    one_large_prime: bool
    usable: bool

    kappa_pass: dict[str, bool]


# ---------------------------------------------------------------------------
# Prime generation
# ---------------------------------------------------------------------------

def sieve_primes(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"

    root = math.isqrt(limit)
    for p in range(2, root + 1):
        if sieve[p]:
            start = p * p
            sieve[start : limit + 1 : p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    return [i for i, flag in enumerate(sieve) if flag]


# ---------------------------------------------------------------------------
# Targets
# ---------------------------------------------------------------------------

def generate_targets(
    primes: list[int],
    count: int,
    rng: random.Random,
) -> list[Target]:
    prime_set = set(primes)

    out: list[Target] = []
    seen_n: set[int] = set()

    while len(out) < count:
        p = rng.choice(primes)
        q = rng.choice(primes)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        # Keep the same rough population used by the prior experiments.
        if not (PRIME_MIN <= p <= PRIME_MAX):
            continue
        if not (PRIME_MIN <= q <= PRIME_MAX):
            continue

        n = p * q
        if n in seen_n:
            continue

        s = p + q

        # Sanity.
        if p not in prime_set or q not in prime_set:
            continue

        seen_n.add(n)
        out.append(
            Target(
                idx=len(out) + 1,
                p=p,
                q=q,
                n=n,
                s=s,
            )
        )

    return out


# ---------------------------------------------------------------------------
# Legendre / QR
# ---------------------------------------------------------------------------

def is_quadratic_residue(a: int, p: int) -> bool:
    """
    QR including zero.

    For prime p:
        a is QR iff a^(p-1)/2 == 1 mod p,
        with zero also accepted.
    """
    a %= p
    if a == 0:
        return True
    return pow(a, (p - 1) // 2, p) == 1


def kappa_pass(q_value: int, family: list[int]) -> bool:
    return all(is_quadratic_residue(q_value, ell) for ell in family)


# ---------------------------------------------------------------------------
# Factorization classification
# ---------------------------------------------------------------------------

def classify_q(
    q_value: int,
    factor_base: list[int],
) -> tuple[bool, bool, bool]:
    """
    Return:
        full_smooth
        one_large_prime
        usable

    A usable relation is either:
        - completely B-smooth
        - B-smooth times exactly one prime <= LARGE_PRIME_BOUND
    """
    if q_value == 0:
        # Degenerate relation. Do not count it.
        return False, False, False

    if q_value < 0:
        q_value = -q_value

    remainder = q_value

    for p in factor_base:
        while remainder % p == 0:
            remainder //= p

        if remainder == 1:
            return True, False, True

    if remainder == 1:
        return True, False, True

    # Exactly one remaining prime, bounded by the LP limit.
    if remainder <= LARGE_PRIME_BOUND and is_prime(remainder):
        return False, True, True

    return False, False, False


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    if n % 3 == 0:
        return n == 3

    r = math.isqrt(n)
    f = 5
    step = 2

    while f <= r:
        if n % f == 0:
            return False
        f += step
        step = 6 - step

    return True


# ---------------------------------------------------------------------------
# Candidate generation
# ---------------------------------------------------------------------------

def ceil_sqrt(n: int) -> int:
    r = math.isqrt(n)
    return r if r * r == n else r + 1


def build_candidates(
    target: Target,
    rng: random.Random,
    factor_base: list[int],
    x_bin_count: int,
    q_bin_count: int,
) -> list[Candidate]:

    x0 = 2 * ceil_sqrt(target.n)

    raw: list[tuple[int, int, float]] = []

    # We generate the exact candidate pool once.
    # Even x only, matching p+q parity for odd p,q.
    x_values = [
        x0 + 2 * k
        for k in range(CANDIDATE_POOL)
    ]

    # Shuffle only the order in which we evaluate, not the candidate set.
    rng.shuffle(x_values)

    for x in x_values:
        q_value = x * x - 4 * target.n

        if q_value <= 0:
            continue

        raw.append(
            (
                x,
                q_value,
                math.log2(q_value),
            )
        )

    if not raw:
        return []

    xs = [r[0] for r in raw]
    qs = [r[1] for r in raw]

    x_min = min(xs)
    x_max = max(xs)
    q_logs = [r[2] for r in raw]
    q_min = min(q_logs)
    q_max = max(q_logs)

    def get_bin(v: float, lo: float, hi: float, count: int) -> int:
        if hi <= lo:
            return 0
        t = (v - lo) / (hi - lo)
        b = int(t * count)
        return min(count - 1, max(0, b))

    out: list[Candidate] = []

    for x, q_value, log_size in raw:
        full_smooth, one_lp, usable = classify_q(
            q_value,
            factor_base,
        )

        kappa = {
            name: kappa_pass(q_value, family)
            for name, family in KAPPA_FAMILIES.items()
        }

        out.append(
            Candidate(
                x=x,
                q_value=q_value,
                log_size=log_size,
                x_bin=get_bin(x, x_min, x_max, x_bin_count),
                q_bin=get_bin(log_size, q_min, q_max, q_bin_count),
                full_smooth=full_smooth,
                one_large_prime=one_lp,
                usable=usable,
                kappa_pass=kappa,
            )
        )

    return out


# ---------------------------------------------------------------------------
# Matched sampling
# ---------------------------------------------------------------------------

def grouped_candidates(
    candidates: list[Candidate],
) -> dict[tuple[int, int], list[Candidate]]:
    groups: dict[tuple[int, int], list[Candidate]] = {}

    for c in candidates:
        key = (c.x_bin, c.q_bin)
        groups.setdefault(key, []).append(c)

    return groups


def make_matched_pairs(
    candidates: list[Candidate],
    family: str,
    rng: random.Random,
    max_pairs: int,
) -> list[tuple[Candidate, Candidate]]:
    """
    Return (kappa_candidate, matched_non_kappa_candidate).

    Matching is done inside exact (x_bin, q_bin) strata.

    Each Kappa observation gets at most one control so the experiment does not
    manufacture effective sample size by repeatedly reusing the same null.
    """
    groups = grouped_candidates(candidates)

    pairs: list[tuple[Candidate, Candidate]] = []

    for members in groups.values():
        kappa = [
            c for c in members
            if c.kappa_pass[family]
        ]

        null = [
            c for c in members
            if not c.kappa_pass[family]
        ]

        if not kappa or not null:
            continue

        rng.shuffle(kappa)
        rng.shuffle(null)

        limit = min(len(kappa), len(null), max_pairs - len(pairs))
        for i in range(limit):
            pairs.append((kappa[i], null[i]))

        if len(pairs) >= max_pairs:
            break

    return pairs


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def rate(values: Iterable[bool]) -> float:
    vals = list(values)
    return sum(vals) / len(vals) if vals else float("nan")


def matched_stats(
    pairs: list[tuple[Candidate, Candidate]],
) -> dict[str, float]:
    if not pairs:
        return {
            "pairs": 0,
            "kappa_usable": float("nan"),
            "null_usable": float("nan"),
            "usable_delta": float("nan"),
            "kappa_full": float("nan"),
            "null_full": float("nan"),
            "full_delta": float("nan"),
            "kappa_lp": float("nan"),
            "null_lp": float("nan"),
            "lp_delta": float("nan"),
        }

    kappa = [a for a, _ in pairs]
    null = [b for _, b in pairs]

    ku = rate(c.usable for c in kappa)
    nu = rate(c.usable for c in null)

    kf = rate(c.full_smooth for c in kappa)
    nf = rate(c.full_smooth for c in null)

    kl = rate(c.one_large_prime for c in kappa)
    nl = rate(c.one_large_prime for c in null)

    return {
        "pairs": len(pairs),
        "kappa_usable": ku,
        "null_usable": nu,
        "usable_delta": ku - nu,
        "kappa_full": kf,
        "null_full": nf,
        "full_delta": kf - nf,
        "kappa_lp": kl,
        "null_lp": nl,
        "lp_delta": kl - nl,
    }


def permutation_test(
    pairs: list[tuple[Candidate, Candidate]],
    rng: random.Random,
    permutations: int,
) -> tuple[float, float]:
    """
    Swap Kappa/null labels independently within each matched pair.

    Statistic:
        mean(usable[Kappa]) - mean(usable[Null])

    Returns:
        observed_delta
        empirical_p
    """
    if not pairs:
        return float("nan"), float("nan")

    observed = statistics.fmean(
        int(a.usable) - int(b.usable)
        for a, b in pairs
    )

    exceed = 0

    for _ in range(permutations):
        total = 0

        for a, b in pairs:
            if rng.random() < 0.5:
                total += int(a.usable) - int(b.usable)
            else:
                total += int(b.usable) - int(a.usable)

        delta = total / len(pairs)

        if delta >= observed:
            exceed += 1

    p_value = (exceed + 1) / (permutations + 1)

    return observed, p_value


def permutation_test_full_smooth(
    pairs: list[tuple[Candidate, Candidate]],
    rng: random.Random,
    permutations: int,
) -> tuple[float, float]:

    if not pairs:
        return float("nan"), float("nan")

    observed = statistics.fmean(
        int(a.full_smooth) - int(b.full_smooth)
        for a, b in pairs
    )

    exceed = 0

    for _ in range(permutations):
        total = 0

        for a, b in pairs:
            if rng.random() < 0.5:
                total += int(a.full_smooth) - int(b.full_smooth)
            else:
                total += int(b.full_smooth) - int(a.full_smooth)

        delta = total / len(pairs)

        if delta >= observed:
            exceed += 1

    return observed, (exceed + 1) / (permutations + 1)


# ---------------------------------------------------------------------------
# Per-target experiment
# ---------------------------------------------------------------------------

def run_target(
    target: Target,
    factor_base: list[int],
    rng: random.Random,
    train: bool,
) -> dict:
    started = time.perf_counter()

    candidates = build_candidates(
        target=target,
        rng=rng,
        factor_base=factor_base,
        x_bin_count=SIZE_BINS,
        q_bin_count=SIZE_BINS,
    )

    result: dict = {
        "target": target,
        "train": train,
        "candidate_count": len(candidates),
        "true_s_in_pool": target.s in {c.x for c in candidates},
    }

    for family in KAPPA_FAMILIES:
        matched_rng = random.Random(
            SEED + target.idx * 1009 + sum(KAPPA_FAMILIES[family])
        )

        pairs = make_matched_pairs(
            candidates=candidates,
            family=family,
            rng=matched_rng,
            max_pairs=MAX_MATCHED_PER_TARGET,
        )

        stats = matched_stats(pairs)

        perm_rng = random.Random(
            SEED + 500_000 + target.idx * 7919
            + sum(KAPPA_FAMILIES[family])
        )

        usable_delta, usable_p = permutation_test(
            pairs,
            perm_rng,
            PERMUTATIONS,
        )

        full_delta, full_p = permutation_test_full_smooth(
            pairs,
            perm_rng,
            PERMUTATIONS,
        )

        result[family] = {
            "pairs": pairs,
            "stats": stats,
            "usable_delta": usable_delta,
            "usable_p": usable_p,
            "full_delta": full_delta,
            "full_p": full_p,
        }

    result["time"] = time.perf_counter() - started

    return result


# ---------------------------------------------------------------------------
# Aggregate reporting
# ---------------------------------------------------------------------------

def aggregate(results: list[dict], family: str) -> dict[str, float]:
    usable_deltas = []
    full_deltas = []
    kappa_usable = []
    null_usable = []
    kappa_full = []
    null_full = []
    pair_counts = []
    p_values = []

    for r in results:
        st = r[family]["stats"]

        if not math.isnan(st["usable_delta"]):
            usable_deltas.append(st["usable_delta"])
            full_deltas.append(st["full_delta"])
            kappa_usable.append(st["kappa_usable"])
            null_usable.append(st["null_usable"])
            kappa_full.append(st["kappa_full"])
            null_full.append(st["null_full"])
            pair_counts.append(st["pairs"])
            p_values.append(r[family]["usable_p"])

    if not usable_deltas:
        return {
            "targets": 0,
            "mean_delta": float("nan"),
            "median_delta": float("nan"),
            "mean_full_delta": float("nan"),
            "mean_kappa_usable": float("nan"),
            "mean_null_usable": float("nan"),
            "mean_kappa_full": float("nan"),
            "mean_null_full": float("nan"),
            "mean_pairs": 0.0,
            "median_p": float("nan"),
        }

    return {
        "targets": len(usable_deltas),
        "mean_delta": statistics.fmean(usable_deltas),
        "median_delta": statistics.median(usable_deltas),
        "mean_full_delta": statistics.fmean(full_deltas),
        "mean_kappa_usable": statistics.fmean(kappa_usable),
        "mean_null_usable": statistics.fmean(null_usable),
        "mean_kappa_full": statistics.fmean(kappa_full),
        "mean_null_full": statistics.fmean(null_full),
        "mean_pairs": statistics.fmean(pair_counts),
        "median_p": statistics.median(p_values),
    }


def print_aggregate(
    name: str,
    stats: dict[str, float],
) -> None:
    print(f"{name}")
    print(
        f"  targets with matches      = {stats['targets']}"
    )
    print(
        f"  mean matched pairs        = {stats['mean_pairs']:.2f}"
    )
    print(
        f"  mean Kappa usable         = {stats['mean_kappa_usable']:.8f}"
    )
    print(
        f"  mean matched null usable = {stats['mean_null_usable']:.8f}"
    )
    print(
        f"  mean usable delta         = {stats['mean_delta']:+.8f}"
    )
    print(
        f"  median usable delta       = "
        f"{stats['median_delta']:+.8f}"
    )
    print(
        f"  mean Kappa full-smooth   = {stats['mean_kappa_full']:.8f}"
    )
    print(
        f"  mean matched null smooth = {stats['mean_null_full']:.8f}"
    )
    print(
        f"  mean full-smooth delta   = "
        f"{stats['mean_full_delta']:+.8f}"
    )
    print(
        f"  median pair p-value      = {stats['median_p']:.6f}"
    )
    print()


# ---------------------------------------------------------------------------
# Family acceptance diagnostic
# ---------------------------------------------------------------------------

def acceptance_rates(
    results: list[dict],
    family: str,
) -> tuple[float, float]:
    k = []
    all_candidates = []

    for r in results:
        for c in reconstruct_candidates_from_result(r):
            all_candidates.append(c)
            k.append(c.kappa_pass[family])

    return rate(k), rate(not x for x in k)


def reconstruct_candidates_from_result(_result: dict) -> list[Candidate]:
    # Candidate objects remain stored through matched-pair references.
    # De-duplicate object identities.
    seen: set[int] = set()
    out: list[Candidate] = []

    for family_data in _result.values():
        if not isinstance(family_data, dict):
            continue

        pairs = family_data.get("pairs")
        if not pairs:
            continue

        for a, b in pairs:
            for c in (a, b):
                key = id(c)
                if key not in seen:
                    seen.add(key)
                    out.append(c)

    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    total_start = time.perf_counter()

    rng = random.Random(SEED)

    print("=" * 78)
    print("KAPPA EXPERIMENT 71R")
    print("CACHED SIZE-MATCHED KAPPA RELATION ENRICHMENT")
    print("WITHIN-Q-SIZE / WITHIN-X CONTROL")
    print("SEQUENTIAL CLASSIFICATION ONCE; MATCHING AFTER CACHE")
    print("L7 PRIMARY + C3 SECONDARY")
    print("GLOBAL MATCHED-PERMUTATION NULL")
    print("NO CSV OUTPUT")
    print("NO SKLEARN")
    print("=" * 78)
    print()

    # Prime population.
    pop_start = time.perf_counter()
    all_primes = sieve_primes(PRIME_MAX)
    population = [
        p for p in all_primes
        if PRIME_MIN <= p <= PRIME_MAX
    ]
    pop_time = time.perf_counter() - pop_start

    print("1. PRIME POPULATION")
    print("-" * 78)
    print(f"prime population = {len(population)}")
    print(f"generation time  = {pop_time:.6f}s")
    print()

    # Targets.
    targets = generate_targets(
        population,
        NUM_TARGETS,
        rng,
    )

    print("2. TARGET SUMMARY")
    print("-" * 78)
    print(f"total targets = {len(targets)}")

    for t in targets[:24]:
        print(
            f"target {t.idx:4d}: "
            f"p={t.p} q={t.q} n={t.n} s={t.s}"
        )

    if len(targets) > 24:
        print("... remaining targets omitted")
    print()

    print("3. KAPPA FAMILIES")
    print("-" * 78)
    for name, family in KAPPA_FAMILIES.items():
        print(f"{name:<12} = {family}")
    print()

    # Factor base.
    factor_base = [
        p for p in all_primes
        if p <= B
    ]

    print("4. RELATION PARAMETERS")
    print("-" * 78)
    print(f"B factor base        = {B}")
    print(f"factor-base primes   = {len(factor_base)}")
    print(f"large-prime bound    = {LARGE_PRIME_BOUND}")
    print(f"candidate pool       = {CANDIDATE_POOL}")
    print(f"size bins            = {SIZE_BINS} x-bins and {SIZE_BINS} q-bins")
    print(f"matched-pair cap     = {MAX_MATCHED_PER_TARGET}")
    print(f"permutations         = {PERMUTATIONS}")
    print()

    # Split by target, not candidate.
    indices = list(range(len(targets)))
    rng.shuffle(indices)

    train_idx = set(indices[:TRAIN_TARGETS])
    test_idx = set(indices[TRAIN_TARGETS:TRAIN_TARGETS + TEST_TARGETS])

    print("5. TARGET HOLDOUT")
    print("-" * 78)
    print(f"training targets = {len(train_idx)}")
    print(f"test targets     = {len(test_idx)}")
    print()

    all_results: list[dict] = []

    for pos, target in enumerate(targets, start=1):
        started = time.perf_counter()

        r = run_target(
            target=target,
            factor_base=factor_base,
            rng=random.Random(SEED + target.idx * 37),
            train=(target.idx in train_idx),
        )

        all_results.append(r)

        if pos % 5 == 0 or pos == len(targets):
            elapsed = time.perf_counter() - started
            print(
                f"processed {pos:3d}/{len(targets)} "
                f"target={target.idx:3d} "
                f"candidates={r['candidate_count']:5d} "
                f"time={elapsed:.3f}s"
            )

    print()
    print("6. DATASET / CACHE SUMMARY")
    print("-" * 78)

    total_candidates = sum(
        r["candidate_count"] for r in all_results
    )

    true_in_pool = sum(
        bool(r["true_s_in_pool"])
        for r in all_results
    )

    print(f"cached candidate classifications = {total_candidates:,}")
    print(
        f"true s inside candidate pool    = "
        f"{true_in_pool}/{len(all_results)} "
        f"({true_in_pool / len(all_results):.6%})"
    )
    print()

    # Family summaries.
    print("7. SIZE-MATCHED RESULTS")
    print("-" * 78)
    print()

    for family in KAPPA_FAMILIES:
        print(f"[{family}]")

        train_results = [
            r for r in all_results
            if r["train"]
        ]

        test_results = [
            r for r in all_results
            if not r["train"]
        ]

        train_stats = aggregate(
            train_results,
            family,
        )

        test_stats = aggregate(
            test_results,
            family,
        )

        print_aggregate(
            "TRAIN",
            train_stats,
        )

        print_aggregate(
            "TEST",
            test_stats,
        )

    # Target-level test enrichment distribution.
    print("8. TEST TARGET ENRICHMENT")
    print("-" * 78)

    for family in KAPPA_FAMILIES:
        deltas = []

        for r in all_results:
            if r["train"]:
                continue

            d = r[family]["stats"]["usable_delta"]
            if not math.isnan(d):
                deltas.append(d)

        if deltas:
            print(
                f"{family:<12} "
                f"mean={statistics.fmean(deltas):+.8f} "
                f"median={statistics.median(deltas):+.8f} "
                f"min={min(deltas):+.8f} "
                f"max={max(deltas):+.8f}"
            )
        else:
            print(f"{family:<12} no matched test data")

    print()

    # Overall paired permutation summary.
    print("9. PERMUTATION NULL")
    print("-" * 78)

    for family in KAPPA_FAMILIES:
        test_p = [
            r[family]["usable_p"]
            for r in all_results
            if not r["train"]
            and not math.isnan(r[family]["usable_p"])
        ]

        if test_p:
            print(
                f"{family:<12} "
                f"median target p={statistics.median(test_p):.6f} "
                f"fraction p<=0.05="
                f"{sum(p <= 0.05 for p in test_p)/len(test_p):.6f}"
            )
        else:
            print(
                f"{family:<12} no permutation results"
            )

    print()

    # Aggregate candidate-level test data, using equal target weighting.
    print("10. TARGET-WEIGHTED TEST SUMMARY")
    print("-" * 78)

    for family in KAPPA_FAMILIES:
        test_stats = aggregate(
            [r for r in all_results if not r["train"]],
            family,
        )

        print(
            f"{family:<12} "
            f"usable_delta={test_stats['mean_delta']:+.8f} "
            f"full_delta={test_stats['mean_full_delta']:+.8f}"
        )

    print()

    # Final diagnostic.
    print("=" * 78)
    print("11. FINAL DIAGNOSTIC")
    print("=" * 78)

    print()
    print(
        "The central comparison is now size matched:"
    )
    print(
        "    Kappa candidate  vs  non-Kappa candidate"
    )
    print(
        "    within the same x-bin and log2(|Q|)-bin."
    )
    print()

    for family in KAPPA_FAMILIES:
        test_stats = aggregate(
            [r for r in all_results if not r["train"]],
            family,
        )

        delta = test_stats["mean_delta"]
        full_delta = test_stats["mean_full_delta"]

        print(f"{family}")
        print(
            f"  test usable delta      = {delta:+.8f}"
        )
        print(
            f"  test full-smooth delta = {full_delta:+.8f}"
        )

        if math.isnan(delta):
            print("  interpretation         = NO MATCHED DATA")
        elif delta <= 0:
            print(
                "  interpretation         = NO ENRICHMENT "
                "AFTER SIZE MATCHING"
            )
        elif full_delta <= 0:
            print(
                "  interpretation         = USABLE GAIN "
                "WITHOUT FULL-SMOOTH GAIN"
            )
        else:
            print(
                "  interpretation         = "
                "PROMISING SIZE-MATCHED ENRICHMENT"
            )

        print()

    print(
        "Important:"
    )
    print(
        "The result is relation-generation evidence only."
    )
    print(
        "It does not establish a faster quadratic-sieve factorization."
    )
    print(
        "A positive result should next be tested by inserting the"
    )
    print(
        "selected relations into an actual QS matrix and measuring"
    )
    print(
        "usable independent relations per CPU second."
    )
    print()

    print(
        f"total runtime = {time.perf_counter() - total_start:.6f}s"
    )
    print()
    print("=" * 78)
    print("EXPERIMENT 71R COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()