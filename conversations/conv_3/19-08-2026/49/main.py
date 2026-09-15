from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Tuple, Dict, Optional


# =============================================================================
# EXPERIMENT 426
# STATIC PRIME-RANKING / GREEDY MODULAR GAP-SQUARE SIEVE
# =============================================================================
#
# Experiment 424:
#   fixed prime order
#
# Experiment 425:
#   dynamically chooses the best next prime, but evaluates every unused
#   prime against the current candidate set at every step.
#
# Experiment 426:
#   measure every prime ONCE against the ORIGINAL candidate set,
#   rank the primes by rejection power, then use that fixed ranking.
#
# QUESTION
#
#   Can a one-time static ranking recover most of the benefit of
#   adaptive prime selection without paying the repeated adaptive
#   selection cost?
#
# RULES
#   exact integer arithmetic only
#   no resultants
#   no Groebner basis
#   no symbolic factorization
#   modular conditions are necessary only
#   exact integer square reconstruction remains the final authority
#
# IMPORTANT
#
#   A modular survivor is NOT accepted as a factorization.
#   It must still pass the exact integer gap-square reconstruction.
#
# =============================================================================


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

UNKNOWN_BITS = [
    20,
    24,
    28,
    32,
    36,
    40,
    44,
    48,
]

MODULAR_PRIMES = [
    3,
    5,
    7,
    11,
    13,
    17,
    19,
    23,
    29,
    31,
    37,
    41,
    43,
    47,
]

ENUMERATION_LIMIT = 250_000

# Only show detailed traces for cases where:
#
#   - the f-candidate set contains multiple candidates, and
#   - the ranking has something meaningful to show.
#
# This keeps output compact.
SHOW_DETAILS = True

# Maximum number of detailed cases shown.
DETAIL_LIMIT = 20


# -----------------------------------------------------------------------------
# Instance representation
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class Instance:
    p: int
    q: int
    N: int
    S: int
    d: int
    gap: int
    K: int


def make_instance(p: int, q: int) -> Instance:
    N = p * q
    S = p + q

    # d = N + 1 - 2S
    d = N + 1 - 2 * S

    gap = q - p

    # d^2 = -4K - 3N^2 + 6N + 1
    K_num = -d * d - 3 * N * N + 6 * N + 1

    assert K_num % 4 == 0

    K = K_num // 4

    # Cross-check the exact identities.
    assert -4 * K - 3 * N * N + 6 * N + 1 == d * d
    assert N * N - N + K == S * (N + 1 - S)

    return Instance(
        p=p,
        q=q,
        N=N,
        S=S,
        d=d,
        gap=gap,
        K=K,
    )


# -----------------------------------------------------------------------------
# Partial K
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class PartialK:
    u: int
    K0: int
    k: int
    C0: int
    d0: int
    r: int
    x_true: int


def partial_k_data(inst: Instance, u: int) -> PartialK:
    scale = 1 << u

    # Floor-aligned decomposition:
    #
    #   K = K0 + k
    #   0 <= k < 2^u
    #
    K0 = (inst.K // scale) * scale
    k = inst.K - K0

    assert 0 <= k < scale
    assert K0 + k == inst.K

    C0 = (
        -4 * K0
        - 3 * inst.N * inst.N
        + 6 * inst.N
        + 1
    )

    d0 = math.isqrt(C0)
    r = C0 - d0 * d0

    x_true = inst.d - d0

    residual = (
        x_true * x_true
        + 2 * d0 * x_true
        + 4 * k
        - r
    )

    assert residual == 0

    return PartialK(
        u=u,
        K0=K0,
        k=k,
        C0=C0,
        d0=d0,
        r=r,
        x_true=x_true,
    )


# -----------------------------------------------------------------------------
# Exact interval
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class IntervalData:
    d_low: int
    d_high: int
    x_low: int
    x_high: int
    candidate_count: int


def exact_d_interval(pk: PartialK) -> IntervalData:
    lo_sq = pk.C0 - 4 * ((1 << pk.u) - 1)
    hi_sq = pk.C0

    if lo_sq <= 0:
        d_low = 0
    else:
        d_low = math.isqrt(lo_sq)

        if d_low * d_low < lo_sq:
            d_low += 1

    d_high = math.isqrt(hi_sq)

    x_low = d_low - pk.d0
    x_high = d_high - pk.d0

    count = max(0, d_high - d_low + 1)

    return IntervalData(
        d_low=d_low,
        d_high=d_high,
        x_low=x_low,
        x_high=x_high,
        candidate_count=count,
    )


# -----------------------------------------------------------------------------
# Exact f candidates
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class Candidate:
    d: int
    x: int
    k: int


def exact_f_candidates(
    pk: PartialK,
    interval: IntervalData,
    limit: int,
) -> Tuple[List[Candidate], bool]:

    if interval.candidate_count > limit:
        return [], False

    out: List[Candidate] = []

    scale = 1 << pk.u

    for d in range(interval.d_low, interval.d_high + 1):

        d2 = d * d

        if d2 > pk.C0:
            continue

        delta = pk.C0 - d2

        # Since delta = 4k.
        if delta % 4 != 0:
            continue

        k = delta // 4

        if not (0 <= k < scale):
            continue

        x = d - pk.d0

        # Explicit exact f check.
        residual = (
            x * x
            + 2 * pk.d0 * x
            + 4 * k
            - pk.r
        )

        if residual != 0:
            continue

        out.append(
            Candidate(
                d=d,
                x=x,
                k=k,
            )
        )

    return out, True


# -----------------------------------------------------------------------------
# Gap-square arithmetic
# -----------------------------------------------------------------------------

def gap_square_value(inst: Instance, d: int) -> int:
    """
    A(d) = (N+1-d)^2 - 16N

    Since

        4g^2 = A(d),

    A(d) must be a perfect square multiple of 4.
    """
    return (
        (inst.N + 1 - d) ** 2
        - 16 * inst.N
    )


def modular_gap_passes(inst: Instance, d: int, p: int) -> bool:
    """
    Necessary quadratic-residue test for

        A(d) = 4g^2.

    Since 4 is itself a square modulo every prime,
    A(d) must be a quadratic residue modulo p.

    Exact zero modulo p is accepted.
    """
    a = gap_square_value(inst, d) % p
    return pow(a, (p - 1) // 2, p) in (0, 1)


def exact_gap_reconstruct(
    inst: Instance,
    d: int,
) -> Optional[Tuple[int, int, int]]:
    """
    Exact reconstruction:

        S = (N+1-d)/2
        g^2 = S^2 - 4N
        p = (S-g)/2
        q = (S+g)/2

    Returns (p, q, g) if everything is exact.
    """
    numerator = inst.N + 1 - d

    if numerator % 2 != 0:
        return None

    S = numerator // 2

    disc = S * S - 4 * inst.N

    if disc < 0:
        return None

    g = math.isqrt(disc)

    if g * g != disc:
        return None

    if (S - g) % 2 != 0:
        return None

    if (S + g) % 2 != 0:
        return None

    p = (S - g) // 2
    q = (S + g) // 2

    if p <= 0 or q <= 0:
        return None

    if p * q != inst.N:
        return None

    if p > q:
        p, q = q, p

    if q - p != g:
        return None

    return p, q, g


# -----------------------------------------------------------------------------
# Candidate-set filtering
# -----------------------------------------------------------------------------

def apply_prime(
    inst: Instance,
    candidates: List[Candidate],
    p: int,
) -> List[Candidate]:

    return [
        c
        for c in candidates
        if modular_gap_passes(inst, c.d, p)
    ]


# -----------------------------------------------------------------------------
# Static ranking
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class PrimeScore:
    prime: int
    survivors: int
    rejected: int
    rejection_ratio_num: int
    rejection_ratio_den: int


def score_prime(
    inst: Instance,
    candidates: List[Candidate],
    p: int,
) -> PrimeScore:

    survivors = 0

    for candidate in candidates:
        if modular_gap_passes(inst, candidate.d, p):
            survivors += 1

    rejected = len(candidates) - survivors

    return PrimeScore(
        prime=p,
        survivors=survivors,
        rejected=rejected,
        rejection_ratio_num=rejected,
        rejection_ratio_den=max(1, len(candidates)),
    )


def rank_primes_once(
    inst: Instance,
    candidates: List[Candidate],
) -> List[PrimeScore]:

    scores = [
        score_prime(inst, candidates, p)
        for p in MODULAR_PRIMES
    ]

    # Best rejection power first.
    #
    # Ties are resolved deterministically:
    #   1. fewer survivors
    #   2. larger rejection count
    #   3. smaller prime
    scores.sort(
        key=lambda s: (
            s.survivors,
            -s.rejected,
            s.prime,
        )
    )

    return scores


# -----------------------------------------------------------------------------
# Fixed modular sieve
# -----------------------------------------------------------------------------

@dataclass
class SieveResult:
    survivors: List[Candidate]
    evaluations: int
    selected_primes: List[int]
    trace: List[Tuple[int, int, int]]


def fixed_sieve(
    inst: Instance,
    candidates: List[Candidate],
) -> SieveResult:

    current = list(candidates)
    evaluations = 0
    trace = []

    for p in MODULAR_PRIMES:

        before = len(current)

        current = apply_prime(
            inst,
            current,
            p,
        )

        evaluations += before

        after = len(current)

        trace.append(
            (p, before, after)
        )

    return SieveResult(
        survivors=current,
        evaluations=evaluations,
        selected_primes=list(MODULAR_PRIMES),
        trace=trace,
    )


# -----------------------------------------------------------------------------
# Static ranking sieve
# -----------------------------------------------------------------------------

@dataclass
class StaticRankedResult:
    survivors: List[Candidate]
    scoring_evaluations: int
    filtering_evaluations: int
    total_evaluations: int
    selected_primes: List[int]
    trace: List[Tuple[int, int, int]]
    scores: List[PrimeScore]


def static_ranked_sieve(
    inst: Instance,
    candidates: List[Candidate],
) -> StaticRankedResult:

    n = len(candidates)

    # ---------------------------------------------------------
    # Phase 1:
    # score every prime once against the ORIGINAL candidate set.
    # ---------------------------------------------------------

    scores = rank_primes_once(
        inst,
        candidates,
    )

    scoring_evaluations = len(MODULAR_PRIMES) * n

    # ---------------------------------------------------------
    # Phase 2:
    # apply the resulting static order.
    # ---------------------------------------------------------

    current = list(candidates)

    filtering_evaluations = 0
    selected_primes: List[int] = []
    trace: List[Tuple[int, int, int]] = []

    for score in scores:

        if len(current) <= 1:
            break

        p = score.prime

        before = len(current)

        current = apply_prime(
            inst,
            current,
            p,
        )

        filtering_evaluations += before

        after = len(current)

        selected_primes.append(p)

        trace.append(
            (p, before, after)
        )

    total = (
        scoring_evaluations
        + filtering_evaluations
    )

    return StaticRankedResult(
        survivors=current,
        scoring_evaluations=scoring_evaluations,
        filtering_evaluations=filtering_evaluations,
        total_evaluations=total,
        selected_primes=selected_primes,
        trace=trace,
        scores=scores,
    )


# -----------------------------------------------------------------------------
# Exact final check
# -----------------------------------------------------------------------------

@dataclass
class ExactResult:
    reconstructed: List[Tuple[Candidate, Tuple[int, int, int]]]
    exact_tests: int


def exact_reconstruction(
    inst: Instance,
    candidates: List[Candidate],
) -> ExactResult:

    reconstructed = []

    for candidate in candidates:

        result = exact_gap_reconstruct(
            inst,
            candidate.d,
        )

        if result is not None:
            reconstructed.append(
                (
                    candidate,
                    result,
                )
            )

    return ExactResult(
        reconstructed=reconstructed,
        exact_tests=len(candidates),
    )


# -----------------------------------------------------------------------------
# Per-test result
# -----------------------------------------------------------------------------

@dataclass
class TestResult:
    instance_index: int
    u: int
    gap: int
    x_true: int

    f_candidate_count: int
    f_candidates: List[Candidate]

    fixed_survivors: List[Candidate]
    static_survivors: List[Candidate]

    fixed_exact: ExactResult
    static_exact: ExactResult

    static_scores: List[PrimeScore]
    static_selected_primes: List[int]

    fixed_modular_evaluations: int
    static_scoring_evaluations: int
    static_filtering_evaluations: int
    static_total_evaluations: int

    fixed_trace: List[Tuple[int, int, int]]
    static_trace: List[Tuple[int, int, int]]

    ranking_changed_result: bool

    enumerated: bool


# -----------------------------------------------------------------------------
# One complete test
# -----------------------------------------------------------------------------

def run_test(
    inst: Instance,
    instance_index: int,
    u: int,
) -> TestResult:

    pk = partial_k_data(
        inst,
        u,
    )

    interval = exact_d_interval(pk)

    candidates, enumerated = exact_f_candidates(
        pk,
        interval,
        ENUMERATION_LIMIT,
    )

    if not enumerated:
        return TestResult(
            instance_index=instance_index,
            u=u,
            gap=inst.gap,
            x_true=pk.x_true,
            f_candidate_count=interval.candidate_count,
            f_candidates=[],
            fixed_survivors=[],
            static_survivors=[],
            fixed_exact=ExactResult([], 0),
            static_exact=ExactResult([], 0),
            static_scores=[],
            static_selected_primes=[],
            fixed_modular_evaluations=0,
            static_scoring_evaluations=0,
            static_filtering_evaluations=0,
            static_total_evaluations=0,
            fixed_trace=[],
            static_trace=[],
            ranking_changed_result=False,
            enumerated=False,
        )

    fixed = fixed_sieve(
        inst,
        candidates,
    )

    ranked = static_ranked_sieve(
        inst,
        candidates,
    )

    fixed_exact = exact_reconstruction(
        inst,
        fixed.survivors,
    )

    static_exact = exact_reconstruction(
        inst,
        ranked.survivors,
    )

    fixed_set = {
        c.d
        for c in fixed.survivors
    }

    static_set = {
        c.d
        for c in ranked.survivors
    }

    ranking_changed_result = (
        fixed_set != static_set
    )

    return TestResult(
        instance_index=instance_index,
        u=u,
        gap=inst.gap,
        x_true=pk.x_true,
        f_candidate_count=len(candidates),
        f_candidates=candidates,
        fixed_survivors=fixed.survivors,
        static_survivors=ranked.survivors,
        fixed_exact=fixed_exact,
        static_exact=static_exact,
        static_scores=ranked.scores,
        static_selected_primes=ranked.selected_primes,
        fixed_modular_evaluations=fixed.evaluations,
        static_scoring_evaluations=ranked.scoring_evaluations,
        static_filtering_evaluations=ranked.filtering_evaluations,
        static_total_evaluations=ranked.total_evaluations,
        fixed_trace=fixed.trace,
        static_trace=ranked.trace,
        ranking_changed_result=ranking_changed_result,
        enumerated=True,
    )


# -----------------------------------------------------------------------------
# Formatting helpers
# -----------------------------------------------------------------------------

def yes_no(value: bool) -> str:
    return "YES" if value else "NO"


def fmt_ratio(n: int, d: int) -> str:
    if d == 0:
        return "0.00000000"
    return f"{n / d:.8f}"


def unique_recovery(
    inst: Instance,
    exact: ExactResult,
) -> bool:

    if len(exact.reconstructed) != 1:
        return False

    candidate, factors = exact.reconstructed[0]

    p, q, g = factors

    return (
        candidate.d == inst.d
        and p * q == inst.N
        and q - p == inst.gap
        and g == inst.gap
    )


def true_survives(
    inst: Instance,
    candidates: List[Candidate],
) -> bool:
    return any(
        c.d == inst.d
        for c in candidates
    )


def survivor_ratio(
    before: int,
    after: int,
) -> float:
    if before == 0:
        return 1.0
    return after / before


# -----------------------------------------------------------------------------
# Main experiment
# -----------------------------------------------------------------------------

def run_experiment() -> None:

    print("=" * 120)
    print("EXPERIMENT 426 START")
    print("=" * 120)
    print()
    print("STATIC PRIME-RANKING / GREEDY MODULAR GAP-SQUARE SIEVE")
    print()
    print("QUESTION")
    print(
        "  Can a one-time prime ranking retain most of the"
        " filtering benefit without repeated adaptive selection?"
    )
    print()
    print("RULES")
    print("  exact integer arithmetic only")
    print("  no resultants")
    print("  no Groebner basis")
    print("  no symbolic factorization")
    print("  modular conditions are necessary only")
    print("  exact gap-square reconstruction remains authoritative")
    print()

    print("=" * 120)
    print("CONFIGURATION")
    print("=" * 120)
    print(f"  instances            = 8")
    print(f"  K bits               = {UNKNOWN_BITS}")
    print(f"  modular primes       = {MODULAR_PRIMES}")
    print(f"  enumeration max      = {ENUMERATION_LIMIT}")
    print("  ranking strategy     = one-time scoring on original f-candidates")

    # -------------------------------------------------------------------------
    # Fresh independent instances
    # -------------------------------------------------------------------------

    instance_specs = [
        (50411, 282599),
        (1013, 10009),
        (10009, 1000033),
        (10009, 10037),
        (50023, 50051),
        (100019, 100043),
        (200009, 200017),
        (300017, 900007),
    ]

    instances = [
        make_instance(p, q)
        for p, q in instance_specs
    ]

    # -------------------------------------------------------------------------
    # Result collection
    # -------------------------------------------------------------------------

    results: List[TestResult] = []

    for idx, inst in enumerate(instances, start=1):

        for u in UNKNOWN_BITS:

            result = run_test(
                inst,
                idx,
                u,
            )

            results.append(result)

    # -------------------------------------------------------------------------
    # Compact summary
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("COMPACT STATIC-RANKING SUMMARY")
    print("=" * 120)

    print(
        " i   u       gap        x    f-cand  fixed  static  "
        "funiq  funiq  ranking primes evals"
    )
    print("-" * 120)

    for result in results:

        if not result.enumerated:
            print(
                f"{result.instance_index:2d} "
                f"{result.u:3d} "
                f"{result.gap:10d} "
                f"{result.x_true:9d} "
                f"{result.f_candidate_count:9d} "
                f"SKIPPED"
            )
            continue

        fixed_count = len(result.fixed_survivors)
        static_count = len(result.static_survivors)

        f_unique = unique_recovery(
            instances[result.instance_index - 1],
            result.fixed_exact,
        )

        static_unique = unique_recovery(
            instances[result.instance_index - 1],
            result.static_exact,
        )

        print(
            f"{result.instance_index:2d} "
            f"{result.u:3d} "
            f"{result.gap:10d} "
            f"{result.x_true:9d} "
            f"{result.f_candidate_count:9d} "
            f"{fixed_count:6d} "
            f"{static_count:6d} "
            f"{yes_no(f_unique):>5s} "
            f"{yes_no(static_unique):>6s} "
            f"{len(result.static_selected_primes):7d} "
            f"{result.static_total_evaluations:10d}"
        )

    # -------------------------------------------------------------------------
    # Global statistics
    # -------------------------------------------------------------------------

    enumerated = [
        r for r in results
        if r.enumerated
    ]

    skipped = [
        r for r in results
        if not r.enumerated
    ]

    fixed_unique_count = 0
    static_unique_count = 0

    fixed_true_survive = 0
    static_true_survive = 0

    total_f_candidates = 0

    total_fixed_evals = 0
    total_static_score_evals = 0
    total_static_filter_evals = 0
    total_static_evals = 0

    total_fixed_exact = 0
    total_static_exact = 0

    ranking_changed = 0

    for r in enumerated:

        inst = instances[r.instance_index - 1]

        total_f_candidates += len(r.f_candidates)

        total_fixed_evals += r.fixed_modular_evaluations
        total_static_score_evals += r.static_scoring_evaluations
        total_static_filter_evals += r.static_filtering_evaluations
        total_static_evals += r.static_total_evaluations

        total_fixed_exact += len(r.fixed_survivors)
        total_static_exact += len(r.static_survivors)

        if unique_recovery(inst, r.fixed_exact):
            fixed_unique_count += 1

        if unique_recovery(inst, r.static_exact):
            static_unique_count += 1

        if true_survives(inst, r.fixed_survivors):
            fixed_true_survive += 1

        if true_survives(inst, r.static_survivors):
            static_true_survive += 1

        if r.ranking_changed_result:
            ranking_changed += 1

    print()
    print("=" * 120)
    print("GLOBAL STATIC-RANKING SUMMARY")
    print("=" * 120)

    print(f"  total tests                         = {len(results)}")
    print(f"  enumerated tests                    = {len(enumerated)}")
    print(f"  skipped tests                       = {len(skipped)}")
    print(
        f"  unique after fixed modular sieve   = "
        f"{fixed_unique_count}/{len(enumerated)}"
    )
    print(
        f"  unique after static ranked sieve   = "
        f"{static_unique_count}/{len(enumerated)}"
    )
    print(
        f"  true survives fixed sieve          = "
        f"{fixed_true_survive}/{len(enumerated)}"
    )
    print(
        f"  true survives static ranking       = "
        f"{static_true_survive}/{len(enumerated)}"
    )

    print()
    print(f"  total f-candidates                 = {total_f_candidates}")
    print(
        f"  fixed modular evaluations          = "
        f"{total_fixed_evals}"
    )
    print(
        f"  static ranking evaluations         = "
        f"{total_static_score_evals}"
    )
    print(
        f"  static filtering evaluations       = "
        f"{total_static_filter_evals}"
    )
    print(
        f"  static total modular evaluations   = "
        f"{total_static_evals}"
    )
    print(
        f"  fixed exact survivor tests         = "
        f"{total_fixed_exact}"
    )
    print(
        f"  static exact survivor tests        = "
        f"{total_static_exact}"
    )

    if total_fixed_evals:
        print(
            "  static/fixed modular evaluation ratio = "
            f"{total_static_evals / total_fixed_evals:.8f}"
        )

    if total_static_score_evals:
        print(
            "  ranking-only fraction of static work = "
            f"{total_static_score_evals / total_static_evals:.8f}"
        )

    print(
        "  cases where fixed/static survivor sets differ = "
        f"{ranking_changed}/{len(enumerated)}"
    )

    # -------------------------------------------------------------------------
    # Aggregate by u
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("AGGREGATE BY UNKNOWN K BITS")
    print("=" * 120)

    print(
        " u   tests enum f-uniq static-uniq "
        "avg-f  avg-fixed  avg-static  avg-primes"
    )
    print("-" * 120)

    for u in UNKNOWN_BITS:

        group = [
            r
            for r in enumerated
            if r.u == u
        ]

        if not group:
            continue

        f_unique = 0
        static_unique = 0

        f_sum = 0
        fixed_sum = 0
        static_sum = 0
        prime_sum = 0

        for r in group:

            inst = instances[r.instance_index - 1]

            f_sum += len(r.f_candidates)
            fixed_sum += len(r.fixed_survivors)
            static_sum += len(r.static_survivors)
            prime_sum += len(r.static_selected_primes)

            if unique_recovery(inst, r.fixed_exact):
                f_unique += 1

            if unique_recovery(inst, r.static_exact):
                static_unique += 1

        n = len(group)

        print(
            f"{u:3d} "
            f"{len(results) // len(UNKNOWN_BITS):5d} "
            f"{n:4d} "
            f"{f_unique:6d} "
            f"{static_unique:11d} "
            f"{f_sum / n:7.2f} "
            f"{fixed_sum / n:10.2f} "
            f"{static_sum / n:10.2f} "
            f"{prime_sum / n:10.2f}"
        )

    # -------------------------------------------------------------------------
    # Ranking quality
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("STATIC PRIME RANKING QUALITY")
    print("=" * 120)

    print(
        "  For each test, primes are ranked once against the original"
        " f-candidate population."
    )
    print(
        "  The table below shows the strongest one-shot primes by"
        " aggregate normalized rejection."
    )

    prime_aggregate: Dict[int, Dict[str, int]] = {
        p: {
            "tests": 0,
            "candidates": 0,
            "survivors": 0,
            "rejected": 0,
        }
        for p in MODULAR_PRIMES
    }

    for r in enumerated:

        for score in r.static_scores:

            a = prime_aggregate[score.prime]

            a["tests"] += 1
            a["candidates"] += score.survivors + score.rejected
            a["survivors"] += score.survivors
            a["rejected"] += score.rejected

    print()
    print(
        " prime    tests     candidates    survivors    "
        "rejected    survival-ratio"
    )
    print("-" * 90)

    ranking_rows = []

    for p in MODULAR_PRIMES:

        a = prime_aggregate[p]

        if a["candidates"] == 0:
            ratio = 1.0
        else:
            ratio = a["survivors"] / a["candidates"]

        ranking_rows.append(
            (
                ratio,
                p,
                a,
            )
        )

    ranking_rows.sort(
        key=lambda item: (
            item[0],
            item[1],
        )
    )

    for ratio, p, a in ranking_rows:

        print(
            f"{p:6d} "
            f"{a['tests']:8d} "
            f"{a['candidates']:13d} "
            f"{a['survivors']:11d} "
            f"{a['rejected']:10d} "
            f"{ratio:.8f}"
        )

    # -------------------------------------------------------------------------
    # Strongest static isolations
    # -------------------------------------------------------------------------

    candidates_for_ranking = []

    for r in enumerated:

        fixed_count = len(r.fixed_survivors)
        static_count = len(r.static_survivors)

        if len(r.f_candidates) == 0:
            continue

        ratio = static_count / len(r.f_candidates)

        candidates_for_ranking.append(
            (
                ratio,
                len(r.f_candidates),
                r,
            )
        )

    candidates_for_ranking.sort(
        key=lambda item: (
            item[0],
            -item[1],
        )
    )

    print()
    print("=" * 120)
    print("STRONGEST STATIC-RANKING ISOLATIONS")
    print("=" * 120)

    for ratio, f_count, r in candidates_for_ranking[:20]:

        print(
            f"  instance={r.instance_index:2d} "
            f"u={r.u:2d} "
            f"f={f_count:8d} "
            f"fixed={len(r.fixed_survivors):6d} "
            f"static={len(r.static_survivors):6d} "
            f"primes={len(r.static_selected_primes):2d} "
            f"ratio={ratio:.8f}"
        )

    # -------------------------------------------------------------------------
    # Detailed cases
    # -------------------------------------------------------------------------

    interesting = []

    for r in enumerated:

        # Show cases where:
        #   * static ranking produces >1 survivor,
        #   * or ranking differs from fixed,
        #   * or candidate population is large.
        #
        # This gives useful diagnostics without printing every test.

        if (
            len(r.static_survivors) > 1
            or r.ranking_changed_result
            or len(r.f_candidates) >= 1000
        ):
            interesting.append(r)

    print()
    print("=" * 120)
    print("INTERESTING STATIC-RANKING CASES")
    print("=" * 120)

    shown = 0

    for r in interesting:

        if shown >= DETAIL_LIMIT:
            break

        inst = instances[r.instance_index - 1]

        print()
        print(
            f"DETAIL: INSTANCE {r.instance_index}, "
            f"UNKNOWN K BITS = {r.u}"
        )
        print("-" * 120)

        print(
            f"p={inst.p} q={inst.q} "
            f"gap={inst.gap} N={inst.N}"
        )
        print(
            f"d={inst.d} K={inst.K}"
        )

        print()
        print("CANDIDATES")
        print(
            f"  f-candidates             = "
            f"{len(r.f_candidates)}"
        )
        print(
            f"  fixed survivors          = "
            f"{len(r.fixed_survivors)}"
        )
        print(
            f"  static ranked survivors  = "
            f"{len(r.static_survivors)}"
        )

        print()
        print("STATIC RANKING")

        for rank, score in enumerate(
            r.static_scores,
            start=1,
        ):
            print(
                f"  rank={rank:2d} "
                f"prime={score.prime:2d} "
                f"survivors={score.survivors:8d} "
                f"rejected={score.rejected:8d} "
                f"survival={score.survivors / max(1, len(r.f_candidates)):.8f}"
            )

        print()
        print("STATIC APPLICATION TRACE")

        for step, (p, before, after) in enumerate(
            r.static_trace,
            start=1,
        ):
            print(
                f"  step={step:2d} "
                f"mod={p:2d} "
                f"{before:8d} -> {after:8d} "
                f"rejected={before-after:8d}"
            )

        print()
        print("WORK")

        print(
            f"  fixed modular evaluations = "
            f"{r.fixed_modular_evaluations}"
        )
        print(
            f"  static scoring evaluations = "
            f"{r.static_scoring_evaluations}"
        )
        print(
            f"  static filtering evaluations = "
            f"{r.static_filtering_evaluations}"
        )
        print(
            f"  static total evaluations = "
            f"{r.static_total_evaluations}"
        )

        print()
        print("EXACT")

        print(
            f"  fixed exact tests = "
            f"{len(r.fixed_survivors)}"
        )
        print(
            f"  static exact tests = "
            f"{len(r.static_survivors)}"
        )

        print(
            "  fixed unique recovery = "
            f"{unique_recovery(inst, r.fixed_exact)}"
        )
        print(
            "  static unique recovery = "
            f"{unique_recovery(inst, r.static_exact)}"
        )

        print(
            "  true survives static = "
            f"{true_survives(inst, r.static_survivors)}"
        )

        shown += 1

    if not interesting:
        print("  none")

    # -------------------------------------------------------------------------
    # Final exact consistency
    # -------------------------------------------------------------------------

    all_fixed_true = all(
        true_survives(
            instances[r.instance_index - 1],
            r.fixed_survivors,
        )
        for r in enumerated
    )

    all_static_true = all(
        true_survives(
            instances[r.instance_index - 1],
            r.static_survivors,
        )
        for r in enumerated
    )

    all_fixed_exact = all(
        unique_recovery(
            instances[r.instance_index - 1],
            r.fixed_exact,
        )
        or len(r.fixed_survivors) != 1
        for r in enumerated
    )

    all_static_exact_consistent = all(
        all(
            factors[0] * factors[1]
            == instances[r.instance_index - 1].N
            for _, factors
            in r.static_exact.reconstructed
        )
        for r in enumerated
    )

    print()
    print("=" * 120)
    print("EXACT INTERNAL CHECKS")
    print("=" * 120)

    print(
        "  all enumerated cases preserve true d through fixed sieve = "
        f"{all_fixed_true}"
    )
    print(
        "  all enumerated cases preserve true d through static ranking = "
        f"{all_static_true}"
    )
    print(
        "  all static reconstructed pairs satisfy pq=N = "
        f"{all_static_exact_consistent}"
    )
    print(
        "  static ranking remains necessary-condition only = True"
    )

    print()
    print("=" * 120)
    print("INTERPRETATION")
    print("=" * 120)

    print(
        "  Experiment 424 used a fixed prime order."
    )
    print(
        "  Experiment 425 dynamically selected the best next prime,"
        " but repeatedly scored every unused prime."
    )
    print(
        "  Experiment 426 scores every prime exactly once on the"
        " original candidate set, then applies the resulting order."
    )
    print()
    print(
        "  This isolates the cost of prime-order discovery from the"
        " cost of repeated adaptive re-ranking."
    )
    print()
    print(
        "  A static ranking can use fewer or equally many selected"
        " primes than the fixed order, but the one-time scoring phase"
        " is itself computational work."
    )
    print()
    print(
        "  Therefore the correct comparison is TOTAL modular"
        " evaluations, not merely the number of selected primes."
    )
    print()
    print(
        "  The modular sieve remains a necessary-condition filter."
    )
    print(
        "  Exact integer gap-square reconstruction remains the"
        " final acceptance condition."
    )
    print()
    print(
        "  This is still an exact candidate-isolation experiment,"
        " not a factoring theorem."
    )

    print()
    print("=" * 120)
    print("EXPERIMENT 426 FINAL STATUS")
    print("=" * 120)

    print(
        f"  ENUMERATED TESTS = {len(enumerated)}"
    )
    print(
        f"  FIXED UNIQUE RECOVERIES = "
        f"{fixed_unique_count}/{len(enumerated)}"
    )
    print(
        f"  STATIC-RANKED UNIQUE RECOVERIES = "
        f"{static_unique_count}/{len(enumerated)}"
    )
    print(
        f"  TOTAL F-CANDIDATES = {total_f_candidates}"
    )
    print(
        f"  FIXED MODULAR EVALUATIONS = {total_fixed_evals}"
    )
    print(
        f"  STATIC RANKING EVALUATIONS = "
        f"{total_static_score_evals}"
    )
    print(
        f"  STATIC FILTERING EVALUATIONS = "
        f"{total_static_filter_evals}"
    )
    print(
        f"  STATIC TOTAL MODULAR EVALUATIONS = "
        f"{total_static_evals}"
    )

    if total_fixed_evals:
        print(
            f"  STATIC/FIXED MODULAR RATIO = "
            f"{total_static_evals / total_fixed_evals:.8f}"
        )

    print(
        "  ALL ENUMERATED TRUE-D-SURVIVAL CHECKS = "
        f"{all_fixed_true and all_static_true}"
    )
    print(
        "  ALL STATIC RECONSTRUCTION CHECKS = "
        f"{all_static_exact_consistent}"
    )

    print("=" * 120)
    print("EXPERIMENT 426 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    run_experiment()
