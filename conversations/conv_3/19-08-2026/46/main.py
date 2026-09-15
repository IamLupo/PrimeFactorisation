from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Tuple


# =============================================================================
# EXPERIMENT 423
# EXACT GAP-SQUARE CANDIDATE FILTER AUDIT
# =============================================================================
#
# Purpose
# -------
#
# Experiment 422 established that the g-free relations produced by the
# augmented (x,k,g) lattice were still exact multiples of
#
#     f(x,k) = x^2 + 2*d0*x + 4*k-r.
#
# Therefore the next question is not "can LLL produce another polynomial?"
# but:
#
#     Can the independent gap identity directly eliminate the surviving
#     d-candidates?
#
# For every candidate d coming from the exact f-equation, define
#
#     S = (N + 1 - d) / 2
#
# and test
#
#     g^2 = S^2 - 4N.
#
# If the RHS is a perfect square and parity is correct, reconstruct
#
#     p = (S-g)/2
#     q = (S+g)/2
#
# and verify
#
#     p*q = N
#     p+q = S
#     q-p = g.
#
# No floating point.
# No resultants.
# No Groebner basis.
# No symbolic factorization.
#
# This is an exact candidate-filter experiment.
# =============================================================================


# -----------------------------------------------------------------------------
# Number-theoretic instance
# -----------------------------------------------------------------------------

@dataclass
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
    d = N + 1 - 2 * S
    gap = q - p

    K = (
        -d * d
        - 3 * N * N
        + 6 * N
        + 1
    ) // 4

    # Exact identities.
    assert (
        -4 * K
        - 3 * N * N
        + 6 * N
        + 1
        == d * d
    )

    assert (
        N * N
        - N
        + K
        == S * (N + 1 - S)
    )

    assert (
        gap * gap
        == S * S - 4 * N
    )

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
# Partial-K decomposition
# -----------------------------------------------------------------------------

@dataclass
class PartialK:
    u: int
    K0: int
    k: int
    C0: int
    d0: int
    r: int
    x_true: int


def partial_k_data(
    inst: Instance,
    u: int,
) -> PartialK:

    modulus = 1 << u

    K0 = (
        inst.K // modulus
    ) * modulus

    k = inst.K - K0

    assert 0 <= k < modulus
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

@dataclass
class IntervalData:
    d_low: int
    d_high: int
    x_low: int
    x_high: int
    candidate_count: int


def exact_interval(
    pk: PartialK,
) -> IntervalData:

    lo_sq = (
        pk.C0
        - 4 * ((1 << pk.u) - 1)
    )

    hi_sq = pk.C0

    if lo_sq <= 0:
        d_low = 0
    else:
        d_low = math.isqrt(lo_sq)

        while d_low * d_low < lo_sq:
            d_low += 1

    d_high = math.isqrt(hi_sq)

    x_low = d_low - pk.d0
    x_high = d_high - pk.d0

    candidate_count = max(
        0,
        d_high - d_low + 1,
    )

    return IntervalData(
        d_low=d_low,
        d_high=d_high,
        x_low=x_low,
        x_high=x_high,
        candidate_count=candidate_count,
    )


# -----------------------------------------------------------------------------
# Candidate from f
# -----------------------------------------------------------------------------

@dataclass
class FCandidate:
    d: int
    x: int
    k: int


def enumerate_f_candidates(
    inst: Instance,
    pk: PartialK,
    interval: IntervalData,
    limit: int,
) -> Tuple[List[FCandidate], bool]:

    if interval.candidate_count > limit:
        return [], False

    candidates: List[FCandidate] = []

    for d in range(
        interval.d_low,
        interval.d_high + 1,
    ):

        x = d - pk.d0

        # Solve the exact polynomial for k:
        #
        #   x^2 + 2*d0*x + 4k-r = 0
        #
        numerator = (
            pk.r
            - x * x
            - 2 * pk.d0 * x
        )

        if numerator % 4 != 0:
            continue

        k = numerator // 4

        if not (
            0 <= k < (1 << pk.u)
        ):
            continue

        residual = (
            x * x
            + 2 * pk.d0 * x
            + 4 * k
            - pk.r
        )

        if residual != 0:
            continue

        candidates.append(
            FCandidate(
                d=d,
                x=x,
                k=k,
            )
        )

    return candidates, True


# -----------------------------------------------------------------------------
# Exact gap-square filter
# -----------------------------------------------------------------------------

@dataclass
class GapFilterResult:
    d: int
    valid_parity: bool
    square: bool
    valid_gap_parity: bool
    reconstructed: bool
    p: int | None
    q: int | None
    g: int | None
    exact_product: bool
    exact_sum: bool
    exact_gap: bool


def check_gap_condition(
    inst: Instance,
    d: int,
) -> GapFilterResult:

    # From
    #
    #   d = N + 1 - 2S
    #
    # obtain
    #
    #   S = (N+1-d)/2.
    #

    numerator = (
        inst.N + 1 - d
    )

    valid_parity = (
        numerator % 2 == 0
    )

    if not valid_parity:
        return GapFilterResult(
            d=d,
            valid_parity=False,
            square=False,
            valid_gap_parity=False,
            reconstructed=False,
            p=None,
            q=None,
            g=None,
            exact_product=False,
            exact_sum=False,
            exact_gap=False,
        )

    S = numerator // 2

    if S < 0:
        return GapFilterResult(
            d=d,
            valid_parity=True,
            square=False,
            valid_gap_parity=False,
            reconstructed=False,
            p=None,
            q=None,
            g=None,
            exact_product=False,
            exact_sum=False,
            exact_gap=False,
        )

    gap_sq = (
        S * S
        - 4 * inst.N
    )

    if gap_sq < 0:
        return GapFilterResult(
            d=d,
            valid_parity=True,
            square=False,
            valid_gap_parity=False,
            reconstructed=False,
            p=None,
            q=None,
            g=None,
            exact_product=False,
            exact_sum=False,
            exact_gap=False,
        )

    g = math.isqrt(gap_sq)

    square = (
        g * g == gap_sq
    )

    if not square:
        return GapFilterResult(
            d=d,
            valid_parity=True,
            square=False,
            valid_gap_parity=False,
            reconstructed=False,
            p=None,
            q=None,
            g=g,
            exact_product=False,
            exact_sum=False,
            exact_gap=False,
        )

    # Since
    #
    #   p = (S-g)/2
    #   q = (S+g)/2
    #
    # p and q must have the same parity as S.
    valid_gap_parity = (
        (S - g) % 2 == 0
        and
        (S + g) % 2 == 0
    )

    if not valid_gap_parity:
        return GapFilterResult(
            d=d,
            valid_parity=True,
            square=True,
            valid_gap_parity=False,
            reconstructed=False,
            p=None,
            q=None,
            g=g,
            exact_product=False,
            exact_sum=False,
            exact_gap=False,
        )

    p = (S - g) // 2
    q = (S + g) // 2

    if p <= 0 or q <= 0:
        return GapFilterResult(
            d=d,
            valid_parity=True,
            square=True,
            valid_gap_parity=True,
            reconstructed=False,
            p=p,
            q=q,
            g=g,
            exact_product=False,
            exact_sum=False,
            exact_gap=False,
        )

    exact_product = (
        p * q == inst.N
    )

    exact_sum = (
        p + q == S
    )

    exact_gap = (
        q - p == g
    )

    reconstructed = (
        exact_product
        and exact_sum
        and exact_gap
    )

    return GapFilterResult(
        d=d,
        valid_parity=True,
        square=True,
        valid_gap_parity=True,
        reconstructed=reconstructed,
        p=p,
        q=q,
        g=g,
        exact_product=exact_product,
        exact_sum=exact_sum,
        exact_gap=exact_gap,
    )


# -----------------------------------------------------------------------------
# One test
# -----------------------------------------------------------------------------

@dataclass
class TestResult:
    instance: int
    u: int

    gap: int
    x_true: int

    interval_count: int
    f_candidates: int

    enumerated: bool

    parity_pass: int
    square_pass: int
    reconstructed_count: int

    true_d_in_f: bool
    true_d_in_gap: bool

    unique_f: bool
    unique_gap: bool

    reduction_ratio_num: int
    reduction_ratio_den: int


def run_test(
    index: int,
    inst: Instance,
    u: int,
    enumeration_limit: int,
) -> TestResult:

    pk = partial_k_data(
        inst,
        u,
    )

    interval = exact_interval(
        pk
    )

    candidates, enumerated = (
        enumerate_f_candidates(
            inst,
            pk,
            interval,
            enumeration_limit,
        )
    )

    if not enumerated:
        return TestResult(
            instance=index,
            u=u,
            gap=inst.gap,
            x_true=pk.x_true,
            interval_count=interval.candidate_count,
            f_candidates=0,
            enumerated=False,
            parity_pass=0,
            square_pass=0,
            reconstructed_count=0,
            true_d_in_f=False,
            true_d_in_gap=False,
            unique_f=False,
            unique_gap=False,
            reduction_ratio_num=0,
            reduction_ratio_den=0,
        )

    parity_pass = 0
    square_pass = 0
    reconstructed_count = 0

    gap_valid_d: List[int] = []

    for candidate in candidates:

        check = check_gap_condition(
            inst,
            candidate.d,
        )

        if check.valid_parity:
            parity_pass += 1

        if check.square:
            square_pass += 1

        if check.reconstructed:
            reconstructed_count += 1
            gap_valid_d.append(
                candidate.d
            )

    true_d_in_f = (
        any(
            c.d == inst.d
            for c in candidates
        )
    )

    true_d_in_gap = (
        inst.d in gap_valid_d
    )

    unique_f = (
        len(candidates) == 1
        and candidates[0].d == inst.d
    )

    unique_gap = (
        len(gap_valid_d) == 1
        and gap_valid_d[0] == inst.d
    )

    return TestResult(
        instance=index,
        u=u,
        gap=inst.gap,
        x_true=pk.x_true,
        interval_count=interval.candidate_count,
        f_candidates=len(candidates),
        enumerated=True,
        parity_pass=parity_pass,
        square_pass=square_pass,
        reconstructed_count=reconstructed_count,
        true_d_in_f=true_d_in_f,
        true_d_in_gap=true_d_in_gap,
        unique_f=unique_f,
        unique_gap=unique_gap,
        reduction_ratio_num=reconstructed_count,
        reduction_ratio_den=len(candidates),
    )


# -----------------------------------------------------------------------------
# Main experiment
# -----------------------------------------------------------------------------

def run_experiment() -> None:

    print("=" * 120)
    print("EXPERIMENT 423 START")
    print("=" * 120)

    print()
    print("EXACT GAP-SQUARE CANDIDATE FILTER AUDIT")

    print()
    print("CORE f RELATION")
    print(
        "  f(x,k) = x^2 + 2*d0*x + 4*k-r"
    )

    print()
    print("GAP RELATION")
    print(
        "  S = (N+1-d)/2"
    )
    print(
        "  g^2 = S^2 - 4N"
    )
    print(
        "  p = (S-g)/2"
    )
    print(
        "  q = (S+g)/2"
    )

    print()
    print("QUESTION")
    print(
        "  How many exact f-candidates survive the independent"
    )
    print(
        "  integer gap/factor consistency conditions?"
    )

    print()
    print("RULES")
    print(
        "  exact integer arithmetic only"
    )
    print(
        "  no resultants"
    )
    print(
        "  no Groebner basis"
    )
    print(
        "  no symbolic factorization"
    )
    print(
        "  no floating point"
    )

    unknown_bits = [
        20,
        24,
        28,
        32,
        36,
        40,
        44,
        48,
    ]

    instances = [
        (50411, 282599),
        (1013, 10009),
        (10009, 1000033),
        (10009, 10037),
        (50023, 50051),
        (100019, 100043),
        (200009, 200017),
        (300017, 900007),
    ]

    enumeration_limit = 200000

    print()
    print("CONFIGURATION")
    print(
        f"  instances       = {len(instances)}"
    )
    print(
        f"  K bits          = {unknown_bits}"
    )
    print(
        f"  enumeration max = {enumeration_limit}"
    )

    generated = [
        make_instance(p, q)
        for p, q in instances
    ]

    results: List[TestResult] = []

    # -------------------------------------------------------------------------
    # Compact summary
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("COMPACT GAP-FILTER SUMMARY")
    print("=" * 120)

    print(
        " i   u       gap        x        interval    f-cand "
        "parity square recon   F-unique G-unique"
    )

    print("-" * 120)

    for index, inst in enumerate(
        generated,
        start=1,
    ):

        for u in unknown_bits:

            result = run_test(
                index,
                inst,
                u,
                enumeration_limit,
            )

            results.append(
                result
            )

            if not result.enumerated:
                print(
                    f"{index:2d} "
                    f"{u:3d} "
                    f"{inst.gap:10d} "
                    f"{result.x_true:9d} "
                    f"{result.interval_count:11d} "
                    f"{'SKIPPED':>7}"
                )
                continue

            print(
                f"{index:2d} "
                f"{u:3d} "
                f"{inst.gap:10d} "
                f"{result.x_true:9d} "
                f"{result.interval_count:11d} "
                f"{result.f_candidates:7d} "
                f"{result.parity_pass:6d} "
                f"{result.square_pass:6d} "
                f"{result.reconstructed_count:5d} "
                f"{'YES' if result.unique_f else 'NO ':>8} "
                f"{'YES' if result.unique_gap else 'NO ':>8}"
            )

    # -------------------------------------------------------------------------
    # Global summary
    # -------------------------------------------------------------------------

    enumerated_results = [
        r
        for r in results
        if r.enumerated
    ]

    print()
    print("=" * 120)
    print("GLOBAL GAP-FILTER SUMMARY")
    print("=" * 120)

    print(
        f"  total tests                    = {len(results)}"
    )

    print(
        f"  enumerated tests               = {len(enumerated_results)}"
    )

    print(
        f"  skipped tests                  = "
        f"{len(results) - len(enumerated_results)}"
    )

    print(
        "  unique from f alone            = "
        f"{sum(r.unique_f for r in enumerated_results)}"
        f"/{len(enumerated_results)}"
    )

    print(
        "  unique after gap filtering     = "
        f"{sum(r.unique_gap for r in enumerated_results)}"
        f"/{len(enumerated_results)}"
    )

    print(
        "  exact reconstructed factor pairs = "
        f"{sum(r.reconstructed_count for r in enumerated_results)}"
    )

    print(
        "  tests where true d survives gap = "
        f"{sum(r.true_d_in_gap for r in enumerated_results)}"
        f"/{len(enumerated_results)}"
    )

    # -------------------------------------------------------------------------
    # Aggregate by u
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("AGGREGATE BY UNKNOWN K BITS")
    print("=" * 120)

    print(
        " u   tests  f-unique gap-unique "
        "avg-f-cand avg-gap-cand total-reconstructed"
    )

    print("-" * 120)

    for u in unknown_bits:

        subset = [
            r
            for r in enumerated_results
            if r.u == u
        ]

        if not subset:
            continue

        total_f = sum(
            r.f_candidates
            for r in subset
        )

        total_gap = sum(
            r.reconstructed_count
            for r in subset
        )

        print(
            f"{u:2d} "
            f"{len(subset):6d} "
            f"{sum(r.unique_f for r in subset):9d} "
            f"{sum(r.unique_gap for r in subset):10d} "
            f"{total_f / len(subset):11.2f} "
            f"{total_gap / len(subset):12.2f} "
            f"{total_gap:18d}"
        )

    # -------------------------------------------------------------------------
    # Interesting cases
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("INTERESTING GAP-FILTER CASES")
    print("=" * 120)

    interesting = False

    for result in enumerated_results:

        # Show cases where gap filtering changes the candidate count,
        # or where f was non-unique but gap becomes unique.
        if (
            result.reconstructed_count
            != result.f_candidates
            or (
                not result.unique_f
                and result.unique_gap
            )
        ):

            interesting = True

            print()
            print(
                f"INSTANCE {result.instance}, "
                f"UNKNOWN K BITS = {result.u}"
            )

            print(
                f"  gap                  = {result.gap}"
            )

            print(
                f"  true x              = {result.x_true}"
            )

            print(
                f"  f interval candidates = {result.interval_count}"
            )

            print(
                f"  f-valid candidates    = {result.f_candidates}"
            )

            print(
                f"  parity-valid         = {result.parity_pass}"
            )

            print(
                f"  square-gap-valid     = {result.square_pass}"
            )

            print(
                f"  reconstructed pairs  = {result.reconstructed_count}"
            )

            print(
                f"  f unique             = {result.unique_f}"
            )

            print(
                f"  gap unique           = {result.unique_gap}"
            )

            print(
                f"  true d survives gap  = {result.true_d_in_gap}"
            )

    if not interesting:
        print(
            "  No candidate count changed under exact gap filtering."
        )

    # -------------------------------------------------------------------------
    # Strongest cases
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("STRONGEST EXACT FILTER CASES")
    print("=" * 120)

    reductions = [
        r
        for r in enumerated_results
        if r.f_candidates > 0
    ]

    reductions.sort(
        key=lambda r: (
            r.reconstructed_count
            / r.f_candidates,
            -r.f_candidates,
        )
    )

    shown = 0

    for result in reductions:

        if shown >= 15:
            break

        if (
            result.reconstructed_count
            == result.f_candidates
        ):
            continue

        fraction = (
            result.reconstructed_count
            / result.f_candidates
        )

        print(
            f"  instance={result.instance:2d} "
            f"u={result.u:2d} "
            f"f={result.f_candidates:8d} "
            f"gap={result.reconstructed_count:8d} "
            f"survival={fraction:.6f}"
        )

        shown += 1

    if shown == 0:
        print(
            "  No strict reductions observed."
        )

    # -------------------------------------------------------------------------
    # Final interpretation
    # -------------------------------------------------------------------------

    total_f_unique = sum(
        r.unique_f
        for r in enumerated_results
    )

    total_gap_unique = sum(
        r.unique_gap
        for r in enumerated_results
    )

    total_reconstructed = sum(
        r.reconstructed_count
        for r in enumerated_results
    )

    total_f_candidates = sum(
        r.f_candidates
        for r in enumerated_results
    )

    print()
    print("=" * 120)
    print("INTERPRETATION")
    print("=" * 120)

    print(
        "  The f-equation produces a finite interval of exact candidates."
    )

    print(
        "  Experiment 423 applies the independent identity"
    )

    print(
        "      g^2 = S^2 - 4N"
    )

    print(
        "  directly to those candidates."
    )

    print()

    print(
        "  A candidate survives only when the implied S gives"
    )

    print(
        "  an exact square discriminant and an exact integer"
    )

    print(
        "  reconstruction of p and q."
    )

    print()

    print(
        f"  Total f-candidates examined = {total_f_candidates}"
    )

    print(
        f"  Total reconstructed factor-pairs = {total_reconstructed}"
    )

    print(
        f"  Unique f recoveries = {total_f_unique}"
    )

    print(
        f"  Unique gap-filter recoveries = {total_gap_unique}"
    )

    print()

    if total_gap_unique > total_f_unique:
        print(
            "  The exact gap condition improves candidate isolation."
        )
    elif total_gap_unique == total_f_unique:
        print(
            "  The exact gap condition does not improve the"
        )
        print(
            "  tested unique-recovery count."
        )
    else:
        print(
            "  Gap filtering exposes cases where the f-only"
        )
        print(
            "  candidate set contains additional false candidates."
        )

    print()
    print(
        "  IMPORTANT:"
    )

    print(
        "  This is an exact candidate filter, not a new factoring theorem."
    )

    print(
        "  The experiment explicitly checks the original N=pq structure"
    )

    print(
        "  through the discriminant identity."
    )

    # -------------------------------------------------------------------------
    # Final status
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("EXPERIMENT 423 FINAL STATUS")
    print("=" * 120)

    print(
        "  ALL EXACT INTERNAL CHECKS = TRUE"
    )

    print(
        f"  F-ONLY UNIQUE RECOVERIES = "
        f"{total_f_unique}/{len(enumerated_results)}"
    )

    print(
        f"  GAP-FILTER UNIQUE RECOVERIES = "
        f"{total_gap_unique}/{len(enumerated_results)}"
    )

    print(
        f"  TOTAL EXACT RECONSTRUCTED PAIRS = "
        f"{total_reconstructed}"
    )

    print("=" * 120)
    print("EXPERIMENT 423 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    run_experiment()
