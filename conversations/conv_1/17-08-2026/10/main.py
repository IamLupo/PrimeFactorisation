#!/usr/bin/env python3
"""
==============================================================================
KAPPA EXPERIMENT 80
LARGE-PRIME QS RELATION PAIRING
KAPPA VS SIZE-MATCHED NULL
ONE-LARGE-PRIME + FULL-SMOOTH RELATIONS
GF(2) RANK AFTER LARGE-PRIME PAIRING
TARGET HOLDOUT
PERMUTATION NULL
NO CSV OUTPUT
NO SKLEARN
==============================================================================

CORE QUESTION
-------------
Experiment 79 only counted fully B-smooth relations.

Experiment 70 showed substantial one-large-prime enrichment.

This experiment asks:

    Does Kappa produce more USABLE and INDEPENDENT QS relations
    after one-large-prime relations are paired?

For each matched Kappa/null candidate:

    Q(x) = x^2 - 4N

Classify Q(x) as:

    FULL:
        Q is completely B-smooth

    ONE-LP:
        Q = smooth_part * r
        where r is prime and r <= LARGE_PRIME_BOUND

Then:

    FULL relations -> one exponent-parity vector

    ONE-LP relations sharing the same r are paired:
        (smooth_a * r) * (smooth_b * r)
        = smooth_a * smooth_b * r^2

    Their GF(2) relation vector is:
        parity(a) XOR parity(b)

For each target/family we measure:

    full smooth count
    one-large-prime count
    paired relation count
    total usable relation count
    GF(2) rank
    rank/candidate
    rank/second

The Kappa and null candidates are size matched by:

    x-bin
    log2(|Q|)-bin

Only equal numbers of candidates enter the comparison.

The decisive metric is:

    held-out GF(2) rank
    and GF(2) rank per CPU second

A smoothness advantage without rank advantage is NOT sufficient.

==============================================================================
"""

from __future__ import annotations

import math
import random
import statistics
import time
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


# ============================================================================
# CONFIG
# ============================================================================

SEED = 800080

NUM_TARGETS = 40
TRAIN_TARGETS = 30
TEST_TARGETS = 10

PRIME_LO = 2_000_000
PRIME_HI = 4_200_000

FACTOR_BASE_BOUND = 1000
LARGE_PRIME_BOUND = 1_000_000

CANDIDATE_POOL = 50_000
MATCHED_PAIR_CAP = 750

X_BINS = 30
Q_BINS = 30

PERMUTATIONS = 500
PROGRESS_EVERY = 5

KAPPA_MODULI = [7]


# ============================================================================
# DATA STRUCTURES
# ============================================================================

@dataclass(frozen=True)
class Target:
    idx: int
    p: int
    q: int
    n: int
    s: int


@dataclass(frozen=True)
class Candidate:
    x: int
    q_abs: int
    kappa: bool


@dataclass
class RelationData:
    full_count: int
    one_lp_count: int
    paired_count: int
    usable_count: int
    rank: int
    factor_time: float
    lp_groups: int


@dataclass
class TargetResult:
    target_id: int
    matched_candidates: int

    kappa: RelationData
    null: RelationData

    classification_time: float
    matching_time: float


# ============================================================================
# PRIME GENERATION
# ============================================================================

def sieve_primes(lo: int, hi: int) -> List[int]:
    if hi < 2 or lo > hi:
        return []

    lo = max(lo, 2)

    limit = int(math.isqrt(hi)) + 1

    small = bytearray(b"\x01") * (limit + 1)
    small[:2] = b"\x00\x00"

    for p in range(2, int(math.isqrt(limit)) + 1):
        if small[p]:
            start = p * p
            small[start: limit + 1: p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    small_primes = [
        p for p in range(2, limit + 1)
        if small[p]
    ]

    result: List[int] = []

    segment_size = 1_000_000
    start = lo

    while start <= hi:
        end = min(start + segment_size - 1, hi)

        mark = bytearray(b"\x01") * (end - start + 1)

        for p in small_primes:
            pp = p * p
            if pp > end:
                break

            first = max(
                pp,
                ((start + p - 1) // p) * p,
            )

            for v in range(first, end + 1, p):
                mark[v - start] = 0

        result.extend(
            start + i
            for i, flag in enumerate(mark)
            if flag
        )

        start = end + 1

    return result


# ============================================================================
# TARGET GENERATION
# ============================================================================

def generate_targets(
    primes: Sequence[int],
    count: int,
    seed: int,
) -> List[Target]:
    rng = random.Random(seed)

    out: List[Target] = []
    seen = set()

    while len(out) < count:
        p = rng.choice(primes)
        q = rng.choice(primes)

        if p == q:
            continue

        if p > q:
            p, q = q, p

        if q - p < 100_000:
            continue

        n = p * q

        if n in seen:
            continue

        seen.add(n)

        out.append(
            Target(
                idx=len(out) + 1,
                p=p,
                q=q,
                n=n,
                s=p + q,
            )
        )

    return out


# ============================================================================
# LEGENDRE / KAPPA
# ============================================================================

def legendre_symbol(a: int, p: int) -> int:
    a %= p

    if a == 0:
        return 0

    v = pow(a, (p - 1) // 2, p)

    if v == 1:
        return 1

    if v == p - 1:
        return -1

    raise RuntimeError(
        f"Invalid Legendre result: a={a}, p={p}, result={v}"
    )


def kappa_pass(q_abs: int) -> bool:
    for ell in KAPPA_MODULI:
        if legendre_symbol(q_abs, ell) != 1:
            return False

    return True


# ============================================================================
# DETERMINISTIC MILLER-RABIN
# ============================================================================

def is_probable_prime(n: int) -> bool:
    """
    Deterministic for the range used by this experiment.
    """

    if n < 2:
        return False

    small_primes = (
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37
    )

    for p in small_primes:
        if n == p:
            return True
        if n % p == 0:
            return False

    d = n - 1
    s = 0

    while d % 2 == 0:
        s += 1
        d //= 2

    for a in (2, 3, 5, 7, 11, 13, 17):
        if a >= n:
            continue

        x = pow(a, d, n)

        if x in (1, n - 1):
            continue

        composite = True

        for _ in range(s - 1):
            x = (x * x) % n

            if x == n - 1:
                composite = False
                break

        if composite:
            return False

    return True


# ============================================================================
# FACTOR BASE
# ============================================================================

def build_factor_base(bound: int) -> List[int]:
    return sieve_primes(2, bound)


def factor_smooth_part(
    value: int,
    factor_base: Sequence[int],
) -> Tuple[Optional[Tuple[int, ...]], int]:
    """
    Divide value completely by the factor base.

    Returns:
        exponent tuple or None if a base factorization is impossible,
        residual cofactor.
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

    return tuple(exponents), v


# ============================================================================
# PARITY VECTOR
# ============================================================================

def parity_vector(exponents: Sequence[int]) -> int:
    bits = 0

    for i, e in enumerate(exponents):
        if e & 1:
            bits |= 1 << i

    return bits


# ============================================================================
# GF(2) BASIS
# ============================================================================

def add_to_basis(
    basis: Dict[int, int],
    vector: int,
) -> bool:
    """
    Incremental GF(2) elimination.

    Returns True iff rank increases.
    """
    x = vector

    while x:
        pivot = x.bit_length() - 1

        existing = basis.get(pivot)

        if existing is None:
            basis[pivot] = x
            return True

        x ^= existing

    return False


# ============================================================================
# CANDIDATE GENERATION
# ============================================================================

def generate_candidate_pool(
    target: Target,
    pool_size: int,
) -> List[Tuple[int, int]]:
    """
    Generate Q(x) values with Q(x) > 0.
    """
    four_n = 4 * target.n

    root = math.isqrt(four_n)

    if root * root < four_n:
        root += 1

    x0 = root

    rng = random.Random(
        SEED ^ (target.idx * 0x9E3779B1)
    )

    xs: List[int] = []
    seen = set()

    sequential = min(pool_size // 2, 2_500_000)

    for dx in range(sequential):
        x = x0 + dx
        xs.append(x)
        seen.add(x)

    while len(xs) < pool_size:
        x = x0 + rng.randrange(0, 2_500_000)

        if x in seen:
            continue

        seen.add(x)
        xs.append(x)

    result = []

    for x in xs:
        q_abs = x * x - four_n

        if q_abs > 0:
            result.append((x, q_abs))

    return result


# ============================================================================
# BIN MATCHING
# ============================================================================

def bin_index(
    value: float,
    low: float,
    high: float,
    bins: int,
) -> int:
    if high <= low:
        return 0

    t = (value - low) / (high - low)

    b = int(t * bins)

    if b < 0:
        return 0

    if b >= bins:
        return bins - 1

    return b


def match_kappa_null(
    kappa: Sequence[Candidate],
    nulls: Sequence[Candidate],
    cap: int,
    seed: int,
) -> Tuple[List[Candidate], List[Candidate], int]:
    """
    Match by common x-bin and common log2(Q)-bin.
    """

    all_candidates = list(kappa) + list(nulls)

    if not all_candidates:
        return [], [], 0

    xs = [c.x for c in all_candidates]
    qs = [
        math.log2(max(1, c.q_abs))
        for c in all_candidates
    ]

    xmin = min(xs)
    xmax = max(xs)
    qmin = min(qs)
    qmax = max(qs)

    kcells: Dict[Tuple[int, int], List[Candidate]] = {}
    ncells: Dict[Tuple[int, int], List[Candidate]] = {}

    for c in kappa:
        cell = (
            bin_index(
                c.x, xmin, xmax, X_BINS
            ),
            bin_index(
                math.log2(max(1, c.q_abs)),
                qmin,
                qmax,
                Q_BINS,
            ),
        )

        kcells.setdefault(cell, []).append(c)

    for c in nulls:
        cell = (
            bin_index(
                c.x, xmin, xmax, X_BINS
            ),
            bin_index(
                math.log2(max(1, c.q_abs)),
                qmin,
                qmax,
                Q_BINS,
            ),
        )

        ncells.setdefault(cell, []).append(c)

    rng = random.Random(seed)

    kselected: List[Candidate] = []
    nselected: List[Candidate] = []

    matched_bins = 0

    for cell in sorted(set(kcells).intersection(ncells)):
        ka = kcells[cell][:]
        na = ncells[cell][:]

        rng.shuffle(ka)
        rng.shuffle(na)

        m = min(len(ka), len(na))

        if m <= 0:
            continue

        matched_bins += 1

        kselected.extend(ka[:m])
        nselected.extend(na[:m])

    if len(kselected) > cap:
        idx = list(range(len(kselected)))
        rng.shuffle(idx)
        idx = idx[:cap]

        kselected = [kselected[i] for i in idx]
        nselected = [nselected[i] for i in idx]

    return kselected, nselected, matched_bins


# ============================================================================
# RELATION EXTRACTION
# ============================================================================

def extract_relations(
    candidates: Sequence[Candidate],
    factor_base: Sequence[int],
) -> RelationData:
    """
    Convert selected candidates into:

        full smooth relations

    and:

        one-large-prime relations.

    Same-large-prime partial relations are paired using a spanning tree.
    """

    start = time.perf_counter()

    full_count = 0
    one_lp_count = 0

    full_vectors: List[int] = []

    lp_groups: Dict[int, List[int]] = {}

    for c in candidates:
        exponents, residual = factor_smooth_part(
            c.q_abs,
            factor_base,
        )

        # Fully smooth.
        if residual == 1:
            full_count += 1

            v = parity_vector(exponents)
            full_vectors.append(v)

            continue

        # Single large prime.
        if (
            residual > 1
            and residual <= LARGE_PRIME_BOUND
            and is_probable_prime(residual)
        ):
            one_lp_count += 1

            v = parity_vector(exponents)

            lp_groups.setdefault(residual, []).append(v)

    # ------------------------------------------------------------------
    # Build GF(2) basis.
    # ------------------------------------------------------------------
    basis: Dict[int, int] = {}
    rank = 0

    # Full smooth rows.
    for v in full_vectors:
        if add_to_basis(basis, v):
            rank += 1

    # ------------------------------------------------------------------
    # Pair one-large-prime relations.
    #
    # For m relations sharing r:
    #
    #   v1 xor v2
    #   v2 xor v3
    #   ...
    #
    # This creates a spanning tree on the LP component.
    # ------------------------------------------------------------------
    paired_count = 0

    for vectors in lp_groups.values():
        if len(vectors) < 2:
            continue

        for i in range(len(vectors) - 1):
            pair_vector = vectors[i] ^ vectors[i + 1]

            paired_count += 1

            if add_to_basis(basis, pair_vector):
                rank += 1

    factor_time = time.perf_counter() - start

    usable_count = full_count + paired_count

    return RelationData(
        full_count=full_count,
        one_lp_count=one_lp_count,
        paired_count=paired_count,
        usable_count=usable_count,
        rank=rank,
        factor_time=factor_time,
        lp_groups=sum(
            1
            for group in lp_groups.values()
            if len(group) >= 2
        ),
    )


# ============================================================================
# TARGET RUN
# ============================================================================

def run_target(
    target: Target,
    factor_base: Sequence[int],
) -> TargetResult:

    # ---------------------------------------------------------------
    # Generate and classify.
    # ---------------------------------------------------------------
    classification_start = time.perf_counter()

    pool = generate_candidate_pool(
        target,
        CANDIDATE_POOL,
    )

    kappa: List[Candidate] = []
    nulls: List[Candidate] = []

    for x, q_abs in pool:
        flag = kappa_pass(q_abs)

        c = Candidate(
            x=x,
            q_abs=q_abs,
            kappa=flag,
        )

        if flag:
            kappa.append(c)
        else:
            nulls.append(c)

    classification_time = (
        time.perf_counter()
        - classification_start
    )

    # ---------------------------------------------------------------
    # Match.
    # ---------------------------------------------------------------
    matching_start = time.perf_counter()

    k_selected, n_selected, matched_bins = (
        match_kappa_null(
            kappa,
            nulls,
            MATCHED_PAIR_CAP,
            SEED ^ target.idx,
        )
    )

    matching_time = (
        time.perf_counter()
        - matching_start
    )

    if not k_selected:
        empty = RelationData(
            full_count=0,
            one_lp_count=0,
            paired_count=0,
            usable_count=0,
            rank=0,
            factor_time=0.0,
            lp_groups=0,
        )

        return TargetResult(
            target_id=target.idx,
            matched_candidates=0,
            kappa=empty,
            null=empty,
            classification_time=classification_time,
            matching_time=matching_time,
        )

    # ---------------------------------------------------------------
    # Relation construction.
    # ---------------------------------------------------------------
    k_result = extract_relations(
        k_selected,
        factor_base,
    )

    n_result = extract_relations(
        n_selected,
        factor_base,
    )

    return TargetResult(
        target_id=target.idx,
        matched_candidates=len(k_selected),
        kappa=k_result,
        null=n_result,
        classification_time=classification_time,
        matching_time=matching_time,
    )


# ============================================================================
# STATISTICS
# ============================================================================

def safe_mean(values: Sequence[float]) -> float:
    return (
        statistics.fmean(values)
        if values
        else float("nan")
    )


def safe_median(values: Sequence[float]) -> float:
    return (
        statistics.median(values)
        if values
        else float("nan")
    )


def summarize_family(
    name: str,
    rows: Sequence[TargetResult],
    family: str,
) -> None:

    if not rows:
        print(f"{name}: no rows")
        return

    rels = [
        getattr(r, family)
        for r in rows
    ]

    cand = [
        r.matched_candidates
        for r in rows
        if r.matched_candidates > 0
    ]

    full_rates = [
        x.full_count / r.matched_candidates
        for r, x in zip(rows, rels)
        if r.matched_candidates > 0
    ]

    lp_rates = [
        x.one_lp_count / r.matched_candidates
        for r, x in zip(rows, rels)
        if r.matched_candidates > 0
    ]

    pair_rates = [
        x.paired_count / r.matched_candidates
        for r, x in zip(rows, rels)
        if r.matched_candidates > 0
    ]

    usable_rates = [
        x.usable_count / r.matched_candidates
        for r, x in zip(rows, rels)
        if r.matched_candidates > 0
    ]

    ranks = [
        x.rank
        for x in rels
    ]

    rank_per_candidate = [
        x.rank / r.matched_candidates
        for r, x in zip(rows, rels)
        if r.matched_candidates > 0
    ]

    total_times = [
        r.classification_time
        + r.matching_time
        + x.factor_time
        for r, x in zip(rows, rels)
    ]

    rank_per_second = [
        x.rank / t
        for x, t in zip(rels, total_times)
        if t > 0
    ]

    print(f"\n{name}")
    print("-" * 76)
    print(
        f"  targets                 = {len(rows)}"
    )
    print(
        f"  mean matched candidates = "
        f"{safe_mean(cand):.2f}"
    )
    print(
        f"  mean full smooth        = "
        f"{safe_mean(full_rates):.8f}"
    )
    print(
        f"  mean one-large-prime    = "
        f"{safe_mean(lp_rates):.8f}"
    )
    print(
        f"  mean paired relations   = "
        f"{safe_mean(pair_rates):.8f}"
    )
    print(
        f"  mean usable rate        = "
        f"{safe_mean(usable_rates):.8f}"
    )
    print(
        f"  mean usable relations   = "
        f"{safe_mean([x.usable_count for x in rels]):.4f}"
    )
    print(
        f"  mean GF(2) rank         = "
        f"{safe_mean(ranks):.4f}"
    )
    print(
        f"  mean rank/candidate     = "
        f"{safe_mean(rank_per_candidate):.8f}"
    )
    print(
        f"  mean rank/second        = "
        f"{safe_mean(rank_per_second):.6f}"
    )
    print(
        f"  median rank             = "
        f"{safe_median(ranks):.4f}"
    )


# ============================================================================
# PERMUTATION NULL
# ============================================================================

def rank_delta_permutation_p(
    rows: Sequence[TargetResult],
    permutations: int,
    seed: int,
) -> Tuple[float, float]:
    """
    Target-level sign permutation of observed Kappa-vs-null rank deltas.

    This is intentionally conservative: it does not pretend individual
    relations are independent observations.
    """

    deltas = [
        r.kappa.rank - r.null.rank
        for r in rows
    ]

    if not deltas:
        return float("nan"), float("nan")

    observed = statistics.fmean(deltas)

    rng = random.Random(seed)

    exceed = 0

    for _ in range(permutations):
        shuffled = [
            d if rng.random() < 0.5 else -d
            for d in deltas
        ]

        stat = statistics.fmean(shuffled)

        if stat >= observed:
            exceed += 1

    p = (exceed + 1) / (permutations + 1)

    return observed, p


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    overall_start = time.perf_counter()

    print("=" * 78)
    print("KAPPA EXPERIMENT 80")
    print("LARGE-PRIME QS RELATION PAIRING")
    print("KAPPA VS SIZE-MATCHED NULL")
    print("FULL-SMOOTH + ONE-LARGE-PRIME RELATIONS")
    print("GF(2) RANK AFTER LARGE-PRIME PAIRING")
    print("TARGET HOLDOUT")
    print("PERMUTATION NULL")
    print("NO CSV OUTPUT")
    print("NO SKLEARN")
    print("=" * 78)

    # ------------------------------------------------------------------
    # Prime population.
    # ------------------------------------------------------------------
    t0 = time.perf_counter()

    primes = sieve_primes(
        PRIME_LO,
        PRIME_HI,
    )

    prime_time = time.perf_counter() - t0

    print("\n1. PRIME POPULATION")
    print("-" * 78)
    print(f"prime population = {len(primes)}")
    print(f"generation time  = {prime_time:.6f}s")

    # ------------------------------------------------------------------
    # Targets.
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
    # Factor base.
    # ------------------------------------------------------------------
    factor_base = build_factor_base(
        FACTOR_BASE_BOUND
    )

    print("\n3. RELATION PARAMETERS")
    print("-" * 78)
    print(
        f"B factor base       = "
        f"{FACTOR_BASE_BOUND}"
    )
    print(
        f"factor-base primes   = "
        f"{len(factor_base)}"
    )
    print(
        f"large-prime bound    = "
        f"{LARGE_PRIME_BOUND}"
    )
    print(
        f"candidate pool       = "
        f"{CANDIDATE_POOL}"
    )
    print(
        f"matched-pair cap     = "
        f"{MATCHED_PAIR_CAP}"
    )
    print(
        f"x bins               = "
        f"{X_BINS}"
    )
    print(
        f"Q-size bins          = "
        f"{Q_BINS}"
    )

    # ------------------------------------------------------------------
    # Kappa validation.
    # ------------------------------------------------------------------
    print("\n4. KAPPA VALIDATION")
    print("-" * 78)

    for ell in KAPPA_MODULI:
        # x^2+x+1 has two roots mod ell iff ell == 1 mod 3.
        if ell % 3 != 1:
            raise RuntimeError(
                f"Kappa modulus {ell} is not 1 mod 3"
            )

        roots = [
            x
            for x in range(1, ell)
            if (x * x + x + 1) % ell == 0
        ]

        if len(roots) != 2:
            raise RuntimeError(
                f"Expected two nontrivial roots for ell={ell}"
            )

        print(
            f"ell={ell:5d} roots={tuple(roots)}"
        )

    # ------------------------------------------------------------------
    # Target split.
    # ------------------------------------------------------------------
    train_targets = targets[:TRAIN_TARGETS]
    test_targets = targets[TRAIN_TARGETS:]

    print("\n5. TARGET HOLDOUT")
    print("-" * 78)
    print(
        f"training targets = "
        f"{len(train_targets)}"
    )
    print(
        f"test targets     = "
        f"{len(test_targets)}"
    )

    # ------------------------------------------------------------------
    # Run.
    # ------------------------------------------------------------------
    results: List[TargetResult] = []

    for i, target in enumerate(targets, 1):
        row = run_target(
            target,
            factor_base,
        )

        results.append(row)

        if (
            i % PROGRESS_EVERY == 0
            or i == len(targets)
        ):
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
    # Training.
    # ------------------------------------------------------------------
    print("\n6. TRAINING RESULTS")
    print("-" * 78)

    summarize_family(
        "KAPPA_L7",
        train_rows,
        "kappa",
    )

    summarize_family(
        "MATCHED_NULL",
        train_rows,
        "null",
    )

    # ------------------------------------------------------------------
    # Test.
    # ------------------------------------------------------------------
    print("\n7. OUT-OF-SAMPLE RESULTS")
    print("-" * 78)

    summarize_family(
        "KAPPA_L7",
        test_rows,
        "kappa",
    )

    summarize_family(
        "MATCHED_NULL",
        test_rows,
        "null",
    )

    # ------------------------------------------------------------------
    # Deltas.
    # ------------------------------------------------------------------
    usable_deltas = [
        (
            r.kappa.usable_count / r.matched_candidates
            - r.null.usable_count / r.matched_candidates
        )
        for r in test_rows
        if r.matched_candidates > 0
    ]

    rank_deltas = [
        r.kappa.rank - r.null.rank
        for r in test_rows
    ]

    lp_deltas = [
        (
            r.kappa.one_lp_count
            - r.null.one_lp_count
        )
        for r in test_rows
    ]

    pair_deltas = [
        (
            r.kappa.paired_count
            - r.null.paired_count
        )
        for r in test_rows
    ]

    print("\n8. TEST TARGET DELTAS")
    print("-" * 78)

    for r in test_rows:
        print(
            f"target {r.target_id:3d}: "
            f"LP K={r.kappa.one_lp_count:3d} "
            f"N={r.null.one_lp_count:3d} "
            f"dLP={r.kappa.one_lp_count-r.null.one_lp_count:+4d} "
            f"pair K={r.kappa.paired_count:3d} "
            f"N={r.null.paired_count:3d} "
            f"drank={r.kappa.rank-r.null.rank:+3d}"
        )

    # ------------------------------------------------------------------
    # Permutation.
    # ------------------------------------------------------------------
    observed_rank_delta, p_rank = (
        rank_delta_permutation_p(
            test_rows,
            PERMUTATIONS,
            SEED ^ 0x5A5A5A,
        )
    )

    print("\n9. TARGET-LEVEL PERMUTATION NULL")
    print("-" * 78)
    print(
        f"observed mean rank delta = "
        f"{observed_rank_delta:+.6f}"
    )
    print(
        f"permutations             = "
        f"{PERMUTATIONS}"
    )
    print(
        f"empirical p-value        = "
        f"{p_rank:.6f}"
    )

    # ------------------------------------------------------------------
    # Core metrics.
    # ------------------------------------------------------------------
    k_test = [r.kappa for r in test_rows]
    n_test = [r.null for r in test_rows]

    k_time = [
        r.classification_time
        + r.matching_time
        + r.kappa.factor_time
        for r in test_rows
    ]

    n_time = [
        r.classification_time
        + r.matching_time
        + r.null.factor_time
        for r in test_rows
    ]

    k_rank_sec = [
        a.rank / t
        for a, t in zip(k_test, k_time)
        if t > 0
    ]

    n_rank_sec = [
        a.rank / t
        for a, t in zip(n_test, n_time)
        if t > 0
    ]

    print("\n10. DECISIVE COMPARISON")
    print("-" * 78)
    print(
        f"test usable-rate delta      = "
        f"{safe_mean(usable_deltas):+.8f}"
    )
    print(
        f"test one-LP count delta     = "
        f"{safe_mean(lp_deltas):+.4f}"
    )
    print(
        f"test paired-relation delta  = "
        f"{safe_mean(pair_deltas):+.4f}"
    )
    print(
        f"test GF(2) rank delta       = "
        f"{safe_mean(rank_deltas):+.4f}"
    )
    print(
        f"Kappa rank/sec              = "
        f"{safe_mean(k_rank_sec):.6f}"
    )
    print(
        f"Null rank/sec               = "
        f"{safe_mean(n_rank_sec):.6f}"
    )

    if (
        observed_rank_delta > 0
        and p_rank <= 0.05
        and safe_mean(k_rank_sec) > safe_mean(n_rank_sec)
    ):
        interpretation = (
            "ALGORITHMICALLY INTERESTING: Kappa creates more "
            "independent QS-style relations after large-prime "
            "pairing and improves rank/sec."
        )
    elif (
        observed_rank_delta > 0
        and p_rank <= 0.05
    ):
        interpretation = (
            "POSSIBLE RELATION EFFECT: Kappa improves held-out "
            "GF(2) rank, but CPU efficiency does not yet improve."
        )
    elif (
        safe_mean(usable_deltas) > 0
        and observed_rank_delta <= 0
    ):
        interpretation = (
            "SMOOTH/LARGE-PRIME ENRICHMENT ONLY: Kappa may create "
            "more partial relations but they do not translate "
            "into greater independent rank."
        )
    else:
        interpretation = (
            "NO DETECTABLE KAPPA ADVANTAGE after one-large-prime "
            "pairing."
        )

    print(
        f"interpretation = {interpretation}"
    )

    # ------------------------------------------------------------------
    # Final.
    # ------------------------------------------------------------------
    runtime = time.perf_counter() - overall_start

    print("\n11. FINAL DIAGNOSTIC")
    print("-" * 78)

    print(
        "This experiment is the direct follow-up to Experiments 70, "
        "71R, and 79."
    )
    print()
    print(
        "Experiment 70 showed relation-generation enrichment."
    )
    print(
        "Experiment 71R showed size-matched enrichment."
    )
    print(
        "Experiment 79 found no advantage when only fully smooth "
        "relations were converted into GF(2) rank."
    )
    print()
    print(
        "Experiment 80 restores the one-large-prime mechanism and "
        "tests whether shared large primes create additional "
        "independent QS relations."
    )
    print()
    print(
        "The strongest positive outcome would be:"
    )
    print(
        "  Kappa > null in paired relations"
    )
    print(
        "  Kappa > null in GF(2) rank"
    )
    print(
        "  positive permutation result"
    )
    print(
        "  Kappa > null in rank/second"
    )
    print()
    print(
        "A positive smoothness effect without a rank effect should "
        "be classified as enrichment rather than factorization progress."
    )
    print()
    print(
        f"total runtime = {runtime:.6f}s"
    )

    print("=" * 78)
    print("EXPERIMENT 80 COMPLETE")
    print("=" * 78)


if __name__ == "__main__":
    main()

