#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 79
KAPPA-INDEPENDENT QS RELATION TEST
FULL-SMOOTH RELATIONS -> GF(2) MATRIX RANK
SIZE-MATCHED KAPPA VS NULL
TARGET HOLDOUT
CPU / RELATION / RANK EFFICIENCY
NO CSV OUTPUT
NO SKLEARN
==============================================================================

CORE QUESTION
-------------
Experiment 71R showed a promising size-matched enrichment in usable Q(x)
relations. This experiment asks whether that enrichment survives the actual
QS-style independence test.

For each target:

    Q(x) = x^2 - 4N

we classify candidates into:

    KAPPA:
        Q(x) satisfies the KAPPA_L7 quadratic-residue filter

    NULL:
        non-Kappa candidates matched to Kappa candidates by:
            - x-bin
            - log2(|Q|)-bin

Only FULLY B-smooth Q(x) values are converted into GF(2) exponent vectors.

Metrics:
    smooth relations
    usable relations
    independent GF(2) rank
    rank / candidate
    rank / CPU second
    smooth / CPU second

The central quantity is:

    independent_relations_per_second

A higher smooth rate without higher matrix rank is NOT considered an
algorithmic improvement.

==============================================================================
"""

from __future__ import annotations

import math
import random
import statistics
import time
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------

SEED = 790079

NUM_TARGETS = 40
TRAIN_TARGETS = 30
TEST_TARGETS = 10

PRIME_LO = 2_000_000
PRIME_HI = 4_200_000

FACTOR_BASE_BOUND = 1000
CANDIDATE_POOL = 50_000

# Number of candidate pairs to compare after Kappa classification.
# Equal count is essential.
MATCHED_CAP_PER_TARGET = 750

# Within-bin matching.
X_BINS = 30
Q_BINS = 30

# Limit candidate search around the QS region.
# x starts slightly above sqrt(4N) = 2 sqrt(N).
X_WINDOW = 2_500_000

KAPPA_L7 = [7]

# Print progress every N targets.
PROGRESS_EVERY = 5


# ---------------------------------------------------------------------------
# DATA STRUCTURES
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class Target:
    idx: int
    p: int
    q: int
    n: int
    s: int


@dataclass
class Relation:
    x: int
    q_abs: int
    factor_vector: int
    factors: Tuple[int, ...]


@dataclass
class TargetStats:
    target_id: int
    kappa_candidates: int
    null_candidates: int

    kappa_full: int
    null_full: int

    kappa_rank: int
    null_rank: int

    kappa_candidate_time: float
    null_candidate_time: float
    kappa_factor_time: float
    null_factor_time: float

    kappa_bins_matched: int
    null_bins_matched: int


# ---------------------------------------------------------------------------
# PRIME GENERATION
# ---------------------------------------------------------------------------

def sieve_primes(lo: int, hi: int) -> List[int]:
    """
    Return all primes in [lo, hi] with a segmented sieve.
    """
    if hi < 2 or lo > hi:
        return []

    lo = max(lo, 2)

    limit = int(math.isqrt(hi)) + 1
    small = bytearray(b"\x01") * (limit + 1)
    small[:2] = b"\x00\x00"

    for p in range(2, int(math.isqrt(limit)) + 1):
        if small[p]:
            start = p * p
            small[start : limit + 1 : p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    small_primes = [
        p for p in range(2, limit + 1) if small[p]
    ]

    segment_size = 1_000_000
    primes: List[int] = []

    start = lo
    while start <= hi:
        end = min(start + segment_size - 1, hi)
        mark = bytearray(b"\x01") * (end - start + 1)

        for p in small_primes:
            pp = p * p
            if pp > end:
                break

            first = max(pp, ((start + p - 1) // p) * p)

            for v in range(first, end + 1, p):
                mark[v - start] = 0

        primes.extend(
            start + i
            for i, is_prime in enumerate(mark)
            if is_prime
        )

        start = end + 1

    return primes


# ---------------------------------------------------------------------------
# TARGET GENERATION
# ---------------------------------------------------------------------------

def generate_targets(
    primes: Sequence[int],
    count: int,
    seed: int,
) -> List[Target]:
    rng = random.Random(seed)

    selected: List[Target] = []
    seen_n = set()

    # Keep factors reasonably separated to avoid mostly-square semiprimes.
    while len(selected) < count:
        p = rng.choice(primes)
        q = rng.choice(primes)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        if q - p < 100_000:
            continue

        n = p * q
        if n in seen_n:
            continue

        seen_n.add(n)

        selected.append(
            Target(
                idx=len(selected) + 1,
                p=p,
                q=q,
                n=n,
                s=p + q,
            )
        )

    return selected


# ---------------------------------------------------------------------------
# LEGENDRE / KAPPA FILTER
# ---------------------------------------------------------------------------

def legendre_symbol(a: int, p: int) -> int:
    """
    Euler criterion for odd prime p.

    Returns:
        +1 residue
        -1 non-residue
         0 zero
    """
    a %= p

    if a == 0:
        return 0

    value = pow(a, (p - 1) // 2, p)

    if value == 1:
        return 1

    if value == p - 1:
        return -1

    raise RuntimeError(
        f"Invalid Legendre result a={a}, p={p}, value={value}"
    )


def kappa_pass(q_abs: int, moduli: Sequence[int]) -> bool:
    """
    KAPPA_L7 filter:
    Q(x) must be a nonzero quadratic residue modulo each selected ell.

    Zero is deliberately rejected so the test measures a genuine
    QR filter rather than forced divisibility.
    """
    for ell in moduli:
        r = legendre_symbol(q_abs % ell, ell)
        if r != 1:
            return False
    return True


# ---------------------------------------------------------------------------
# QS FACTOR BASE
# ---------------------------------------------------------------------------

def build_factor_base(bound: int) -> List[int]:
    return sieve_primes(2, bound)


def factor_over_base(
    value: int,
    factor_base: Sequence[int],
) -> Optional[Tuple[int, ...]]:
    """
    Factor |value| completely over the factor base.

    Return exponent tuple if fully smooth, otherwise None.
    """
    v = abs(value)

    exponents = [0] * len(factor_base)

    for i, p in enumerate(factor_base):
        if p * p > v:
            break

        while v % p == 0:
            v //= p
            exponents[i] += 1

        if v == 1:
            break

    if v != 1:
        return None

    return tuple(exponents)


def exponent_parity_vector(exponents: Sequence[int]) -> int:
    """
    Pack parity of exponent vector into a Python int.

    Bit i corresponds to factor-base prime i.
    """
    bits = 0

    for i, e in enumerate(exponents):
        if e & 1:
            bits |= 1 << i

    return bits


# ---------------------------------------------------------------------------
# GF(2) INCREMENTAL RANK
# ---------------------------------------------------------------------------

def gf2_rank(vectors: Iterable[int]) -> int:
    """
    Incremental Gaussian elimination over GF(2).
    """
    basis: Dict[int, int] = {}
    rank = 0

    for v in vectors:
        x = v

        while x:
            pivot = x.bit_length() - 1

            existing = basis.get(pivot)

            if existing is None:
                basis[pivot] = x
                rank += 1
                break

            x ^= existing

    return rank


def gf2_incremental_add(
    basis: Dict[int, int],
    v: int,
) -> bool:
    """
    Add one vector to a GF(2) basis.

    Return True if rank increased.
    """
    x = v

    while x:
        pivot = x.bit_length() - 1
        existing = basis.get(pivot)

        if existing is None:
            basis[pivot] = x
            return True

        x ^= existing

    return False


# ---------------------------------------------------------------------------
# BINNING / MATCHING
# ---------------------------------------------------------------------------

def log2_q(q_abs: int) -> int:
    return max(0, q_abs.bit_length() - 1)


def make_bin(
    value: float,
    vmin: float,
    vmax: float,
    bins: int,
) -> int:
    if vmax <= vmin:
        return 0

    t = (value - vmin) / (vmax - vmin)
    b = int(t * bins)

    if b < 0:
        b = 0
    elif b >= bins:
        b = bins - 1

    return b


def build_binned_candidates(
    candidates: Sequence[Tuple[int, int, bool]],
    x_bins: int,
    q_bins: int,
) -> Dict[Tuple[int, int], List[Tuple[int, int, bool]]]:
    if not candidates:
        return {}

    xs = [c[0] for c in candidates]
    qs = [math.log2(max(1, c[1])) for c in candidates]

    xmin, xmax = min(xs), max(xs)
    qmin, qmax = min(qs), max(qs)

    out: Dict[Tuple[int, int], List[Tuple[int, int, bool]]] = {}

    for c in candidates:
        x, q_abs, flag = c

        bx = make_bin(
            float(x),
            float(xmin),
            float(xmax),
            x_bins,
        )

        bq = make_bin(
            math.log2(max(1, q_abs)),
            qmin,
            qmax,
            q_bins,
        )

        out.setdefault((bx, bq), []).append(c)

    return out


def matched_pairs(
    kappa: Sequence[Tuple[int, int, bool]],
    nonkappa: Sequence[Tuple[int, int, bool]],
    cap: int,
    seed: int,
) -> Tuple[List[Tuple[int, int]], int]:
    """
    Match Kappa and null candidates inside identical (x-bin,Q-size-bin)
    cells.

    Returns:
        list of (kappa_index, null_index)
        number of nonempty matched bins
    """
    rng = random.Random(seed)

    kb = build_binned_candidates(kappa, X_BINS, Q_BINS)
    nb = build_binned_candidates(nonkappa, X_BINS, Q_BINS)

    # To make bin definitions identical, re-bin jointly.
    combined = list(kappa) + list(nonkappa)

    if not combined:
        return [], 0

    xs = [c[0] for c in combined]
    qs = [math.log2(max(1, c[1])) for c in combined]

    xmin, xmax = min(xs), max(xs)
    qmin, qmax = min(qs), max(qs)

    kcells: Dict[Tuple[int, int], List[int]] = {}
    ncells: Dict[Tuple[int, int], List[int]] = {}

    for i, c in enumerate(kappa):
        x, q_abs, _ = c

        cell = (
            make_bin(float(x), float(xmin), float(xmax), X_BINS),
            make_bin(
                math.log2(max(1, q_abs)),
                qmin,
                qmax,
                Q_BINS,
            ),
        )

        kcells.setdefault(cell, []).append(i)

    for i, c in enumerate(nonkappa):
        x, q_abs, _ = c

        cell = (
            make_bin(float(x), float(xmin), float(xmax), X_BINS),
            make_bin(
                math.log2(max(1, q_abs)),
                qmin,
                qmax,
                Q_BINS,
            ),
        )

        ncells.setdefault(cell, []).append(i)

    pairs: List[Tuple[int, int]] = []
    usable_bins = 0

    common = sorted(set(kcells).intersection(ncells))

    for cell in common:
        a = kcells[cell][:]
        b = ncells[cell][:]

        rng.shuffle(a)
        rng.shuffle(b)

        m = min(len(a), len(b))

        if m <= 0:
            continue

        usable_bins += 1

        for j in range(m):
            pairs.append((a[j], b[j]))

    if len(pairs) > cap:
        rng.shuffle(pairs)
        pairs = pairs[:cap]

    return pairs, usable_bins


# ---------------------------------------------------------------------------
# CANDIDATE GENERATION
# ---------------------------------------------------------------------------

def generate_candidate_pool(
    target: Target,
    pool_size: int,
) -> List[Tuple[int, int]]:
    """
    Generate x and |Q(x)| values.

    Candidate x values lie above 2*sqrt(N), where Q(x) >= 0.
    """
    root = math.isqrt(4 * target.n)

    if root * root < 4 * target.n:
        root += 1

    x0 = root

    rng = random.Random(
        SEED ^ (target.idx * 0x9E3779B1)
    )

    xs: List[int] = []

    # Sequential core plus randomized tail.
    # This avoids pathological patterns in the test.
    sequential_count = min(pool_size // 2, X_WINDOW)

    for dx in range(sequential_count):
        xs.append(x0 + dx)

    remaining = pool_size - len(xs)

    seen = set(xs)

    while len(xs) < pool_size:
        x = x0 + rng.randrange(0, X_WINDOW)

        if x in seen:
            continue

        seen.add(x)
        xs.append(x)

    out: List[Tuple[int, int]] = []

    for x in xs:
        qv = x * x - 4 * target.n

        if qv <= 0:
            continue

        out.append((x, qv))

    return out


# ---------------------------------------------------------------------------
# RELATION PROCESSING
# ---------------------------------------------------------------------------

def process_full_smooth_relations(
    entries: Sequence[Tuple[int, int, bool]],
    factor_base: Sequence[int],
) -> Tuple[int, int, float]:
    """
    entries = (x, |Q|, is_kappa)

    Return:
        full_smooth_count,
        GF(2) rank,
        elapsed_factorization_time
    """
    start = time.perf_counter()

    basis: Dict[int, int] = {}
    full = 0
    rank = 0

    for _, q_abs, _ in entries:
        exponents = factor_over_base(q_abs, factor_base)

        if exponents is None:
            continue

        full += 1

        vector = exponent_parity_vector(exponents)

        if gf2_incremental_add(basis, vector):
            rank += 1

    elapsed = time.perf_counter() - start

    return full, rank, elapsed


# ---------------------------------------------------------------------------
# TARGET EXPERIMENT
# ---------------------------------------------------------------------------

def run_target(
    target: Target,
    factor_base: Sequence[int],
) -> TargetStats:
    pool = generate_candidate_pool(
        target,
        CANDIDATE_POOL,
    )

    # ------------------------------------------------------------------
    # Classify all candidates once.
    # ------------------------------------------------------------------
    classify_start = time.perf_counter()

    kappa: List[Tuple[int, int, bool]] = []
    nonkappa: List[Tuple[int, int, bool]] = []

    for x, q_abs in pool:
        flag = kappa_pass(q_abs, KAPPA_L7)

        item = (x, q_abs, flag)

        if flag:
            kappa.append(item)
        else:
            nonkappa.append(item)

    classify_elapsed = time.perf_counter() - classify_start

    # ------------------------------------------------------------------
    # Size match.
    # ------------------------------------------------------------------
    pairs, matched_bins = matched_pairs(
        kappa,
        nonkappa,
        MATCHED_CAP_PER_TARGET,
        SEED ^ target.idx,
    )

    if not pairs:
        return TargetStats(
            target_id=target.idx,
            kappa_candidates=0,
            null_candidates=0,
            kappa_full=0,
            null_full=0,
            kappa_rank=0,
            null_rank=0,
            kappa_candidate_time=classify_elapsed,
            null_candidate_time=0.0,
            kappa_factor_time=0.0,
            null_factor_time=0.0,
            kappa_bins_matched=0,
            null_bins_matched=0,
        )

    kappa_entries = [
        kappa[i]
        for i, _ in pairs
    ]

    null_entries = [
        nonkappa[j]
        for _, j in pairs
    ]

    # ------------------------------------------------------------------
    # Factor / rank.
    # ------------------------------------------------------------------
    k_full, k_rank, k_factor_time = (
        process_full_smooth_relations(
            kappa_entries,
            factor_base,
        )
    )

    n_full, n_rank, n_factor_time = (
        process_full_smooth_relations(
            null_entries,
            factor_base,
        )
    )

    return TargetStats(
        target_id=target.idx,
        kappa_candidates=len(kappa_entries),
        null_candidates=len(null_entries),
        kappa_full=k_full,
        null_full=n_full,
        kappa_rank=k_rank,
        null_rank=n_rank,
        kappa_candidate_time=classify_elapsed,
        null_candidate_time=0.0,
        kappa_factor_time=k_factor_time,
        null_factor_time=n_factor_time,
        kappa_bins_matched=matched_bins,
        null_bins_matched=matched_bins,
    )


# ---------------------------------------------------------------------------
# AGGREGATION
# ---------------------------------------------------------------------------

def safe_mean(values: Sequence[float]) -> float:
    return statistics.fmean(values) if values else float("nan")


def safe_median(values: Sequence[float]) -> float:
    return statistics.median(values) if values else float("nan")


def summarize(
    name: str,
    rows: Sequence[TargetStats],
) -> None:
    if not rows:
        print(f"\n{name}\n  NO MATCHED TARGETS")
        return

    usable_k = [
        r.kappa_full / r.kappa_candidates
        for r in rows
        if r.kappa_candidates
    ]

    usable_n = [
        r.null_full / r.null_candidates
        for r in rows
        if r.null_candidates
    ]

    rank_k = [
        r.kappa_rank
        for r in rows
    ]

    rank_n = [
        r.null_rank
        for r in rows
    ]

    rank_per_candidate_k = [
        r.kappa_rank / r.kappa_candidates
        for r in rows
        if r.kappa_candidates
    ]

    rank_per_candidate_n = [
        r.null_rank / r.null_candidates
        for r in rows
        if r.null_candidates
    ]

    k_total_time = [
        r.kappa_candidate_time + r.kappa_factor_time
        for r in rows
    ]

    n_total_time = [
        r.null_candidate_time + r.null_factor_time
        for r in rows
    ]

    rank_per_sec_k = [
        r.kappa_rank / t
        for r, t in zip(rows, k_total_time)
        if t > 0
    ]

    rank_per_sec_n = [
        r.null_rank / t
        for r, t in zip(rows, n_total_time)
        if t > 0
    ]

    print(f"\n{name}")
    print("-" * 76)
    print(
        f"  targets with matches           = {len(rows)}"
    )
    print(
        f"  mean matched candidates        = "
        f"{safe_mean([r.kappa_candidates for r in rows]):.2f}"
    )
    print(
        f"  mean full-smooth rate          = "
        f"{safe_mean(usable_k):.8f}"
    )
    print(
        f"  mean null full-smooth rate     = "
        f"{safe_mean(usable_n):.8f}"
    )
    print(
        f"  mean smooth-rate delta         = "
        f"{safe_mean([a-b for a,b in zip(usable_k, usable_n)]):+.8f}"
    )
    print(
        f"  mean GF(2) rank                = "
        f"{safe_mean(rank_k):.4f}"
    )
    print(
        f"  mean null GF(2) rank           = "
        f"{safe_mean(rank_n):.4f}"
    )
    print(
        f"  mean rank delta                = "
        f"{safe_mean([a-b for a,b in zip(rank_k, rank_n)]):+.4f}"
    )
    print(
        f"  mean rank/candidate            = "
        f"{safe_mean(rank_per_candidate_k):.8f}"
    )
    print(
        f"  mean null rank/candidate       = "
        f"{safe_mean(rank_per_candidate_n):.8f}"
    )
    print(
        f"  mean rank/sec                  = "
        f"{safe_mean(rank_per_sec_k):.4f}"
    )
    print(
        f"  mean null rank/sec             = "
        f"{safe_mean(rank_per_sec_n):.4f}"
    )
    print(
        f"  median rank delta              = "
        f"{safe_median([a-b for a,b in zip(rank_k, rank_n)]):+.4f}"
    )


# ---------------------------------------------------------------------------
# PERMUTATION NULL
# ---------------------------------------------------------------------------

def permutation_rank_delta(
    rows: Sequence[TargetStats],
    permutations: int,
    seed: int,
) -> float:
    """
    Lightweight target-level permutation null.

    We preserve each target's paired counts/ranks and randomly swap which
    family receives the observed rank counts.

    This is intentionally conservative and does not assume independent
    relations across candidates.
    """
    if not rows or permutations <= 0:
        return float("nan")

    observed = statistics.fmean(
        r.kappa_rank - r.null_rank
        for r in rows
    )

    rng = random.Random(seed)

    exceed = 0

    for _ in range(permutations):
        deltas = []

        for r in rows:
            if rng.random() < 0.5:
                deltas.append(
                    r.kappa_rank - r.null_rank
                )
            else:
                deltas.append(
                    r.null_rank - r.kappa_rank
                )

        stat = statistics.fmean(deltas)

        if stat >= observed:
            exceed += 1

    return (exceed + 1) / (permutations + 1)


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main() -> None:
    overall_start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 79")
    print("KAPPA-INDEPENDENT QS RELATION TEST")
    print("FULL-SMOOTH RELATIONS -> GF(2) MATRIX RANK")
    print("SIZE-MATCHED KAPPA VS NULL")
    print("TARGET HOLDOUT")
    print("CPU / RELATION / RANK EFFICIENCY")
    print("NO CSV OUTPUT")
    print("NO SKLEARN")
    print("=" * 78)

    # ------------------------------------------------------------------
    # Prime population
    # ------------------------------------------------------------------
    t0 = time.perf_counter()

    primes = sieve_primes(
        PRIME_LO,
        PRIME_HI,
    )

    generation_time = time.perf_counter() - t0

    print("\n1. PRIME POPULATION")
    print("-" * 78)
    print(f"prime population = {len(primes)}")
    print(f"generation time  = {generation_time:.6f}s")

    if len(primes) < NUM_TARGETS * 4:
        raise RuntimeError(
            "Prime population unexpectedly small"
        )

    # ------------------------------------------------------------------
    # Targets
    # ------------------------------------------------------------------
    targets = generate_targets(
        primes,
        NUM_TARGETS,
        SEED,
    )

    print("\n2. TARGET SUMMARY")
    print("-" * 78)

    for t in targets[:24]:
        print(
            f"target {t.idx:4d}: "
            f"p={t.p} q={t.q} n={t.n} s={t.s}"
        )

    if len(targets) > 24:
        print("... remaining generated targets omitted")

    # ------------------------------------------------------------------
    # Factor base
    # ------------------------------------------------------------------
    factor_base = build_factor_base(
        FACTOR_BASE_BOUND
    )

    print("\n3. RELATION / FACTOR-BASE PARAMETERS")
    print("-" * 78)
    print(
        f"B factor base       = {FACTOR_BASE_BOUND}"
    )
    print(
        f"factor-base primes   = {len(factor_base)}"
    )
    print(
        f"candidate pool       = {CANDIDATE_POOL}"
    )
    print(
        f"matched-pair cap     = {MATCHED_CAP_PER_TARGET}"
    )
    print(
        f"x bins               = {X_BINS}"
    )
    print(
        f"Q-size bins          = {Q_BINS}"
    )

    # ------------------------------------------------------------------
    # Kappa sanity checks
    # ------------------------------------------------------------------
    print("\n4. KAPPA VALIDATION")
    print("-" * 78)

    for ell in KAPPA_L7:
        if not (
            all(
                pow(r, 2, ell) == (-r - 1) % ell
                for r in [2]  # basic domain sanity
            )
            or ell == 7
        ):
            raise RuntimeError(
                f"Unexpected Kappa modulus validation failure: {ell}"
            )

        print(
            f"ell={ell:5d} accepted as KAPPA modulus"
        )

    # ------------------------------------------------------------------
    # Holdout
    # ------------------------------------------------------------------
    train_targets = targets[:TRAIN_TARGETS]
    test_targets = targets[TRAIN_TARGETS:]

    print("\n5. TARGET HOLDOUT")
    print("-" * 78)
    print(
        f"training targets = {len(train_targets)}"
    )
    print(
        f"test targets     = {len(test_targets)}"
    )

    # ------------------------------------------------------------------
    # Run
    # ------------------------------------------------------------------
    results: List[TargetStats] = []

    for i, target in enumerate(targets, start=1):
        row = run_target(
            target,
            factor_base,
        )
        results.append(row)

        if i % PROGRESS_EVERY == 0 or i == len(targets):
            print(
                f"processed {i:3d}/{len(targets)} "
                f"target={target.idx:3d}"
            )

    train_rows = [
        r for r in results
        if r.target_id <= TRAIN_TARGETS
    ]

    test_rows = [
        r for r in results
        if r.target_id > TRAIN_TARGETS
    ]

    # ------------------------------------------------------------------
    # Results
    # ------------------------------------------------------------------
    print("\n6. TRAINING RESULTS")
    print("-" * 78)
    summarize(
        "KAPPA_L7 vs MATCHED NULL",
        train_rows,
    )

    print("\n7. OUT-OF-SAMPLE RESULTS")
    print("-" * 78)
    summarize(
        "KAPPA_L7 vs MATCHED NULL",
        test_rows,
    )

    # ------------------------------------------------------------------
    # Per-target test deltas
    # ------------------------------------------------------------------
    print("\n8. TEST TARGET RANK DELTAS")
    print("-" * 78)

    for r in test_rows:
        delta = r.kappa_rank - r.null_rank

        k_rate = (
            r.kappa_full / r.kappa_candidates
            if r.kappa_candidates else float("nan")
        )

        n_rate = (
            r.null_full / r.null_candidates
            if r.null_candidates else float("nan")
        )

        print(
            f"target {r.target_id:3d}: "
            f"rank K={r.kappa_rank:3d} "
            f"N={r.null_rank:3d} "
            f"delta={delta:+4d} "
            f"smooth K={k_rate:.5f} "
            f"N={n_rate:.5f}"
        )

    # ------------------------------------------------------------------
    # Permutation
    # ------------------------------------------------------------------
    p_perm = permutation_rank_delta(
        test_rows,
        permutations=500,
        seed=SEED ^ 0x123456,
    )

    print("\n9. GLOBAL MATCHED PERMUTATION NULL")
    print("-" * 78)

    observed_delta = safe_mean(
        [
            r.kappa_rank - r.null_rank
            for r in test_rows
        ]
    )

    print(
        f"observed mean rank delta = "
        f"{observed_delta:+.6f}"
    )
    print(
        f"permutations             = 500"
    )
    print(
        f"empirical p-value        = "
        f"{p_perm:.6f}"
    )

    # ------------------------------------------------------------------
    # Independence vs smoothness
    # ------------------------------------------------------------------
    test_smooth_delta = safe_mean(
        [
            (
                r.kappa_full / r.kappa_candidates
                - r.null_full / r.null_candidates
            )
            for r in test_rows
            if (
                r.kappa_candidates > 0
                and r.null_candidates > 0
            )
        ]
    )

    test_rank_delta = safe_mean(
        [
            r.kappa_rank - r.null_rank
            for r in test_rows
        ]
    )

    print("\n10. KEY COMPARISON")
    print("-" * 78)
    print(
        f"test smooth-rate delta = "
        f"{test_smooth_delta:+.8f}"
    )
    print(
        f"test GF(2) rank delta   = "
        f"{test_rank_delta:+.6f}"
    )

    if test_smooth_delta > 0 and test_rank_delta <= 0:
        interpretation = (
            "LOCAL ENRICHMENT ONLY: Kappa creates more smooth values "
            "but not more independent relations."
        )
    elif test_rank_delta > 0 and p_perm <= 0.05:
        interpretation = (
            "PROMISING INDEPENDENCE EFFECT: Kappa produces more "
            "independent QS-style relations under the matched null."
        )
    elif test_rank_delta > 0:
        interpretation = (
            "POSSIBLE INDEPENDENCE EFFECT: positive rank delta, "
            "but the permutation test is not yet decisive."
        )
    else:
        interpretation = (
            "NO INDEPENDENCE ADVANTAGE: Kappa does not improve "
            "relation-matrix rank under the matched null."
        )

    print(
        f"interpretation = {interpretation}"
    )

    # ------------------------------------------------------------------
    # Final diagnostic
    # ------------------------------------------------------------------
    runtime = time.perf_counter() - overall_start

    print("\n11. FINAL DIAGNOSTIC")
    print("-" * 78)

    print(
        "The decisive quantity is not smoothness alone."
    )
    print(
        "A useful Kappa filter should increase independent "
        "GF(2) relation rank per unit CPU time."
    )

    print()
    print(
        "Kappa should be considered algorithmically interesting only if:"
    )
    print(
        "  1. matched Kappa candidates have more rank,"
    )
    print(
        "  2. the effect survives unseen targets,"
    )
    print(
        "  3. the permutation p-value is small,"
    )
    print(
        "  4. rank/second improves rather than merely "
        "smooth-rate/candidate."
    )

    print()
    print(
        f"total runtime = {runtime:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 79 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

