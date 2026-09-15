from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Tuple


# =============================================================================
# EXPERIMENT 424
# MODULAR GAP-SQUARE SIEVE / EXACT CANDIDATE-ISOLATION AUDIT
# =============================================================================
#
# Motivation
# ----------
#
# Experiment 423 showed:
#
#   f(x,k)=0
#
# can leave many exact d-candidates, while the independent identity
#
#   S = (N+1-d)/2
#   g^2 = S^2 - 4N
#
# reduces the candidate set to exactly one valid factor pair in every
# enumerated test.
#
# Experiment 424 asks a more structural question:
#
#   How much of this exact gap-square filtering can be reproduced using
#   only modular square-residue constraints?
#
# For each candidate d:
#
#   D = N + 1 - d
#
# and
#
#   4*g^2 = D^2 - 16N.
#
# Therefore, for every odd prime l,
#
#   D^2 - 16N
#
# must be a quadratic residue modulo l.
#
# We use several small primes as an exact modular sieve.
#
# No floating point.
# No resultants.
# No Groebner basis.
# No symbolic factorization.
#
# The experiment compares:
#
#   1. f-candidates
#   2. candidates surviving modular gap tests
#   3. candidates surviving exact square testing
#   4. exact reconstruction of p,q
#
# The modular stage is only a necessary-condition sieve.
# A modular survivor is NOT automatically a valid factorization.
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

# Primes used for the modular square sieve.
#
# 2 is handled separately through exact parity.
#
# These primes are deliberately small so the sieve is cheap.
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

# We keep the exact candidate enumeration bounded.
ENUMERATION_LIMIT = 250_000

# Only print detailed information for cases where the sieve actually removes
# candidates.
MAX_INTERESTING_CASES = 40


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
# Partial K
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

    K0 = (inst.K // modulus) * modulus
    k = inst.K - K0

    assert 0 <= k < modulus

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

    return IntervalData(
        d_low=d_low,
        d_high=d_high,
        candidate_count=max(
            0,
            d_high - d_low + 1,
        ),
    )


# -----------------------------------------------------------------------------
# Exact f-candidate enumeration
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

    # The relation implies d has the same parity pattern as required by
    # the exact K construction. We nevertheless test every d here.
    for d in range(
        interval.d_low,
        interval.d_high + 1,
    ):

        x = d - pk.d0

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
# Modular square tables
# -----------------------------------------------------------------------------

@dataclass
class PrimeResidueTable:
    prime: int
    square_residues: Tuple[bool, ...]


def build_square_table(p: int) -> PrimeResidueTable:
    table = [False] * p

    for x in range(p):
        table[(x * x) % p] = True

    return PrimeResidueTable(
        prime=p,
        square_residues=tuple(table),
    )


def build_residue_tables(
    primes: List[int],
) -> List[PrimeResidueTable]:

    return [
        build_square_table(p)
        for p in primes
    ]


# -----------------------------------------------------------------------------
# Candidate modular test
# -----------------------------------------------------------------------------
#
# We use:
#
#   S = (N+1-d)/2
#
# but modular division by 2 is avoided completely.
#
# Instead:
#
#   4*g^2 = (N+1-d)^2 - 16N.
#
# So for every odd prime p, the RHS must be a quadratic residue mod p.
#
# This avoids any modular inverse and works directly on integers.
# -----------------------------------------------------------------------------

def modular_gap_passes(
    N: int,
    d: int,
    tables: List[PrimeResidueTable],
) -> bool:

    D = N + 1 - d

    rhs = (
        D * D
        - 16 * N
    )

    # Exact parity condition.
    #
    # For integer S=(N+1-d)/2 we need:
    #
    #   N+1-d == 0 mod 2.
    #
    if D & 1:
        return False

    for table in tables:

        p = table.prime

        residue = rhs % p

        if not table.square_residues[residue]:
            return False

    return True


# -----------------------------------------------------------------------------
# Progressive modular sieve
# -----------------------------------------------------------------------------

@dataclass
class ModularStage:
    prime: int
    survivors: int


def progressive_modular_sieve(
    inst: Instance,
    candidates: List[FCandidate],
    tables: List[PrimeResidueTable],
) -> Tuple[List[FCandidate], List[ModularStage]]:

    current = candidates[:]
    stages: List[ModularStage] = []

    # First impose the exact parity condition.
    parity_survivors = []

    for candidate in current:
        D = inst.N + 1 - candidate.d

        if (D & 1) == 0:
            parity_survivors.append(candidate)

    current = parity_survivors

    stages.append(
        ModularStage(
            prime=2,
            survivors=len(current),
        )
    )

    for table in tables:

        p = table.prime

        Deltas = []

        for candidate in current:

            D = inst.N + 1 - candidate.d

            rhs = (
                D * D
                - 16 * inst.N
            )

            if table.square_residues[rhs % p]:
                Deltas.append(candidate)

        current = Deltas

        stages.append(
            ModularStage(
                prime=p,
                survivors=len(current),
            )
        )

        if not current:
            break

    return current, stages


# -----------------------------------------------------------------------------
# Exact gap test
# -----------------------------------------------------------------------------

@dataclass
class ExactGapResult:
    d: int
    valid: bool
    S: int | None
    g: int | None
    p: int | None
    q: int | None


def exact_gap_test(
    inst: Instance,
    d: int,
) -> ExactGapResult:

    D = inst.N + 1 - d

    if D & 1:
        return ExactGapResult(
            d=d,
            valid=False,
            S=None,
            g=None,
            p=None,
            q=None,
        )

    S = D // 2

    discriminant = (
        S * S
        - 4 * inst.N
    )

    if discriminant < 0:
        return ExactGapResult(
            d=d,
            valid=False,
            S=S,
            g=None,
            p=None,
            q=None,
        )

    g = math.isqrt(discriminant)

    if g * g != discriminant:
        return ExactGapResult(
            d=d,
            valid=False,
            S=S,
            g=g,
            p=None,
            q=None,
        )

    if (S - g) & 1:
        return ExactGapResult(
            d=d,
            valid=False,
            S=S,
            g=g,
            p=None,
            q=None,
        )

    p = (S - g) // 2
    q = (S + g) // 2

    if p <= 0 or q <= 0:
        return ExactGapResult(
            d=d,
            valid=False,
            S=S,
            g=g,
            p=p,
            q=q,
        )

    if p * q != inst.N:
        return ExactGapResult(
            d=d,
            valid=False,
            S=S,
            g=g,
            p=p,
            q=q,
        )

    if p + q != S:
        return ExactGapResult(
            d=d,
            valid=False,
            S=S,
            g=g,
            p=p,
            q=q,
        )

    if q - p != g:
        return ExactGapResult(
            d=d,
            valid=False,
            S=S,
            g=g,
            p=p,
            q=q,
        )

    return ExactGapResult(
        d=d,
        valid=True,
        S=S,
        g=g,
        p=p,
        q=q,
    )


# -----------------------------------------------------------------------------
# Test result
# -----------------------------------------------------------------------------

@dataclass
class TestResult:
    instance: int
    u: int
    gap: int
    x_true: int

    interval_count: int
    f_candidate_count: int

    enumerated: bool

    stages: List[ModularStage]

    final_modular_count: int
    exact_gap_count: int

    f_unique: bool
    modular_unique: bool
    exact_unique: bool

    true_d_in_f: bool
    true_d_in_modular: bool
    true_d_in_exact: bool


# -----------------------------------------------------------------------------
# One experiment point
# -----------------------------------------------------------------------------

def run_test(
    index: int,
    inst: Instance,
    u: int,
    tables: List[PrimeResidueTable],
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
            f_candidate_count=0,
            enumerated=False,
            stages=[],
            final_modular_count=0,
            exact_gap_count=0,
            f_unique=False,
            modular_unique=False,
            exact_unique=False,
            true_d_in_f=False,
            true_d_in_modular=False,
            true_d_in_exact=False,
        )

    modular_candidates, stages = (
        progressive_modular_sieve(
            inst,
            candidates,
            tables,
        )
    )

    exact_candidates: List[FCandidate] = []

    for candidate in modular_candidates:

        result = exact_gap_test(
            inst,
            candidate.d,
        )

        if result.valid:
            exact_candidates.append(
                candidate
            )

    f_unique = (
        len(candidates) == 1
        and candidates[0].d == inst.d
    )

    modular_unique = (
        len(modular_candidates) == 1
        and modular_candidates[0].d == inst.d
    )

    exact_unique = (
        len(exact_candidates) == 1
        and exact_candidates[0].d == inst.d
    )

    return TestResult(
        instance=index,
        u=u,
        gap=inst.gap,
        x_true=pk.x_true,
        interval_count=interval.candidate_count,
        f_candidate_count=len(candidates),
        enumerated=True,
        stages=stages,
        final_modular_count=len(modular_candidates),
        exact_gap_count=len(exact_candidates),
        f_unique=f_unique,
        modular_unique=modular_unique,
        exact_unique=exact_unique,
        true_d_in_f=any(
            c.d == inst.d
            for c in candidates
        ),
        true_d_in_modular=any(
            c.d == inst.d
            for c in modular_candidates
        ),
        true_d_in_exact=any(
            c.d == inst.d
            for c in exact_candidates
        ),
    )


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def run_experiment() -> None:

    print("=" * 120)
    print("EXPERIMENT 424 START")
    print("=" * 120)

    print()
    print("MODULAR GAP-SQUARE SIEVE / EXACT CANDIDATE-ISOLATION AUDIT")

    print()
    print("CORE f RELATION")
    print(
        "  f(x,k) = x^2 + 2*d0*x + 4*k-r"
    )

    print()
    print("GAP CONDITION")
    print(
        "  4*g^2 = (N+1-d)^2 - 16N"
    )

    print()
    print("QUESTION")
    print(
        "  How much candidate reduction can be obtained from"
    )
    print(
        "  cheap quadratic-residue conditions before the exact"
    )
    print(
        "  integer square test?"
    )

    print()
    print("RULES")
    print("  exact integer arithmetic only")
    print("  no resultants")
    print("  no Groebner basis")
    print("  no symbolic factorization")
    print("  no floating point")
    print("  modular conditions are necessary only")

    print()
    print("CONFIGURATION")
    print(
        f"  instances          = 8"
    )
    print(
        f"  K bits             = {UNKNOWN_BITS}"
    )
    print(
        f"  modular primes     = {MODULAR_PRIMES}"
    )
    print(
        f"  enumeration max    = {ENUMERATION_LIMIT}"
    )

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

    generated = [
        make_instance(p, q)
        for p, q in instances
    ]

    tables = build_residue_tables(
        MODULAR_PRIMES
    )

    results: List[TestResult] = []

    # -------------------------------------------------------------------------
    # Compact summary
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("COMPACT MODULAR-SIEVE SUMMARY")
    print("=" * 120)

    print(
        " i   u       gap        x   interval   f-cand   mod-final "
        " exact  F-uniq M-uniq E-uniq"
    )

    print("-" * 120)

    for index, inst in enumerate(
        generated,
        start=1,
    ):

        for u in UNKNOWN_BITS:

            result = run_test(
                index,
                inst,
                u,
                tables,
                ENUMERATION_LIMIT,
            )

            results.append(result)

            if not result.enumerated:

                print(
                    f"{index:2d} "
                    f"{u:3d} "
                    f"{inst.gap:10d} "
                    f"{result.x_true:9d} "
                    f"{result.interval_count:10d} "
                    f"{'SKIPPED':>8}"
                )

                continue

            print(
                f"{index:2d} "
                f"{u:3d} "
                f"{inst.gap:10d} "
                f"{result.x_true:9d} "
                f"{result.interval_count:10d} "
                f"{result.f_candidate_count:8d} "
                f"{result.final_modular_count:9d} "
                f"{result.exact_gap_count:7d} "
                f"{'YES' if result.f_unique else 'NO ':>6} "
                f"{'YES' if result.modular_unique else 'NO ':>6} "
                f"{'YES' if result.exact_unique else 'NO ':>6}"
            )

    enumerated = [
        r
        for r in results
        if r.enumerated
    ]

    # -------------------------------------------------------------------------
    # Global summary
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("GLOBAL MODULAR-SIEVE SUMMARY")
    print("=" * 120)

    print(
        f"  total tests                    = {len(results)}"
    )

    print(
        f"  enumerated tests               = {len(enumerated)}"
    )

    print(
        f"  skipped tests                  = "
        f"{len(results) - len(enumerated)}"
    )

    print(
        f"  unique from f alone            = "
        f"{sum(r.f_unique for r in enumerated)}/{len(enumerated)}"
    )

    print(
        f"  unique after modular sieve    = "
        f"{sum(r.modular_unique for r in enumerated)}/{len(enumerated)}"
    )

    print(
        f"  unique after exact gap test   = "
        f"{sum(r.exact_unique for r in enumerated)}/{len(enumerated)}"
    )

    print(
        f"  true d survives modular sieve = "
        f"{sum(r.true_d_in_modular for r in enumerated)}/{len(enumerated)}"
    )

    print(
        f"  true d survives exact gap      = "
        f"{sum(r.true_d_in_exact for r in enumerated)}/{len(enumerated)}"
    )

    # -------------------------------------------------------------------------
    # Aggregate by bit size
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("AGGREGATE BY UNKNOWN K BITS")
    print("=" * 120)

    print(
        " u   tests f-uniq mod-uniq exact-uniq "
        "avg-f avg-mod avg-exact"
    )

    print("-" * 120)

    for u in UNKNOWN_BITS:

        subset = [
            r
            for r in enumerated
            if r.u == u
        ]

        if not subset:
            continue

        avg_f = (
            sum(r.f_candidate_count for r in subset)
            / len(subset)
        )

        avg_mod = (
            sum(r.final_modular_count for r in subset)
            / len(subset)
        )

        avg_exact = (
            sum(r.exact_gap_count for r in subset)
            / len(subset)
        )

        print(
            f"{u:2d} "
            f"{len(subset):6d} "
            f"{sum(r.f_unique for r in subset):6d} "
            f"{sum(r.modular_unique for r in subset):8d} "
            f"{sum(r.exact_unique for r in subset):9d} "
            f"{avg_f:10.2f} "
            f"{avg_mod:10.2f} "
            f"{avg_exact:11.2f}"
        )

    # -------------------------------------------------------------------------
    # Progressive sieve examples
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("PROGRESSIVE MODULAR FILTER EXAMPLES")
    print("=" * 120)

    interesting_count = 0

    for result in enumerated:

        if (
            result.f_candidate_count <= 1
            or result.final_modular_count
            == result.f_candidate_count
        ):
            continue

        print()
        print(
            f"INSTANCE {result.instance}, "
            f"u={result.u}"
        )

        print(
            f"  initial f-candidates = "
            f"{result.f_candidate_count}"
        )

        for stage in result.stages:

            print(
                f"  after modulus {stage.prime:2d}"
                f" = {stage.survivors}"
            )

        print(
            f"  exact-square survivors = "
            f"{result.exact_gap_count}"
        )

        interesting_count += 1

        if (
            interesting_count
            >= MAX_INTERESTING_CASES
        ):
            break

    if interesting_count == 0:
        print(
            "  No strict modular reductions observed."
        )

    # -------------------------------------------------------------------------
    # Best modular reductions
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("STRONGEST MODULAR REDUCTIONS")
    print("=" * 120)

    scored = [
        r
        for r in enumerated
        if r.f_candidate_count > 0
    ]

    scored.sort(
        key=lambda r: (
            r.final_modular_count
            / r.f_candidate_count,
            -r.f_candidate_count,
        )
    )

    shown = 0

    for result in scored:

        if shown >= 15:
            break

        if (
            result.final_modular_count
            >= result.f_candidate_count
        ):
            continue

        ratio = (
            result.final_modular_count
            / result.f_candidate_count
        )

        print(
            f"  instance={result.instance:2d} "
            f"u={result.u:2d} "
            f"f={result.f_candidate_count:8d} "
            f"mod={result.final_modular_count:8d} "
            f"ratio={ratio:.8f}"
        )

        shown += 1

    if shown == 0:
        print(
            "  No strict reductions observed."
        )

    # -------------------------------------------------------------------------
    # Final interpretation
    # -------------------------------------------------------------------------

    total_f = sum(
        r.f_candidate_count
        for r in enumerated
    )

    total_mod = sum(
        r.final_modular_count
        for r in enumerated
    )

    total_exact = sum(
        r.exact_gap_count
        for r in enumerated
    )

    print()
    print("=" * 120)
    print("INTERPRETATION")
    print("=" * 120)

    print(
        "  Experiment 423 used the exact discriminant test"
    )

    print(
        "  directly on every f-candidate."
    )

    print()

    print(
        "  Experiment 424 replaces that expensive exact-square"
    )

    print(
        "  condition with a chain of necessary quadratic-residue"
    )

    print(
        "  conditions modulo small primes."
    )

    print()

    print(
        f"  Total f-candidates = {total_f}"
    )

    print(
        f"  Total modular survivors = {total_mod}"
    )

    print(
        f"  Total exact survivors = {total_exact}"
    )

    if total_f:
        print(
            f"  Overall modular survival = "
            f"{total_mod / total_f:.8f}"
        )

        print(
            f"  Overall exact survival = "
            f"{total_exact / total_f:.8f}"
        )

    print()

    print(
        "  The modular sieve can only reject impossible candidates."
    )

    print(
        "  A modular survivor still requires the exact integer"
    )

    print(
        "  square test before it can be accepted."
    )

    print()

    print(
        "  This experiment therefore measures whether the independent"
    )

    print(
        "  gap constraint has a useful cheap arithmetic footprint"
    )

    print(
        "  before exact square-root evaluation."
    )

    # -------------------------------------------------------------------------
    # Final status
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("EXPERIMENT 424 FINAL STATUS")
    print("=" * 120)

    print(
        "  ALL EXACT INTERNAL CHECKS = TRUE"
    )

    print(
        f"  F-CANDIDATES = {total_f}"
    )

    print(
        f"  MODULAR SURVIVORS = {total_mod}"
    )

    print(
        f"  EXACT GAP SURVIVORS = {total_exact}"
    )

    print(
        f"  MODULAR UNIQUE RECOVERIES = "
        f"{sum(r.modular_unique for r in enumerated)}"
        f"/{len(enumerated)}"
    )

    print(
        f"  EXACT UNIQUE RECOVERIES = "
        f"{sum(r.exact_unique for r in enumerated)}"
        f"/{len(enumerated)}"
    )

    print("=" * 120)
    print("EXPERIMENT 424 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    run_experiment()
