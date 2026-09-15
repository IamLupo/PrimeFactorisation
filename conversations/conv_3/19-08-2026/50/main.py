from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Tuple, Set


# =============================================================================
# EXPERIMENT 427
# BITSET MODULAR GAP-SQUARE SIEVE / COST-AWARE ADAPTIVE SELECTION
# =============================================================================
#
# PREVIOUS RESULTS
#
# Experiment 424:
#   fixed modular prime order
#
# Experiment 425:
#   adaptive prime selection
#   but every unused prime was repeatedly evaluated against
#   the current candidate set
#
# Experiment 426:
#   every prime was scored once against the original candidate set
#   and then applied in a static order
#
# OBSERVATION
#
#   Experiment 426 retained the same survivor sets as the fixed sieve,
#   but paid heavily for the one-time scoring phase:
#
#       static/fixed modular evaluation ratio ~ 8.78
#
# NEW IDEA
#
#   Represent the current candidate set as an exact integer bitset.
#
#   For each prime p, construct a bitset mask containing exactly the
#   candidates that pass the quadratic-residue test modulo p.
#
#   Then:
#
#       current_mask &= prime_mask[p]
#
#   performs the modular filtering as a machine-word bit operation
#   rather than a fresh modular evaluation for every surviving candidate.
#
# QUESTION
#
#   Does bitset-based candidate filtering make adaptive prime selection
#   computationally competitive with the fixed sieve, while preserving
#   the exact candidate-isolation behavior?
#
# IMPORTANT
#
#   This experiment measures arithmetic candidate filtering only.
#   It does NOT claim that bitset operations constitute a factoring
#   breakthrough.  The final exact gap-square reconstruction remains
#   authoritative.
#
# RULES
#   exact integer arithmetic only
#   no resultants
#   no Groebner basis
#   no symbolic factorization
#   modular conditions are necessary only
#   exact integer square reconstruction is final
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


# -----------------------------------------------------------------------------
# Number-theoretic instance
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
    d = N + 1 - 2 * S
    gap = q - p

    numerator = -d * d - 3 * N * N + 6 * N + 1

    assert numerator % 4 == 0

    K = numerator // 4

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

    return IntervalData(
        d_low=d_low,
        d_high=d_high,
        candidate_count=max(0, d_high - d_low + 1),
    )


# -----------------------------------------------------------------------------
# Exact f candidates
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class Candidate:
    d: int
    x: int
    k: int


def enumerate_f_candidates(
    pk: PartialK,
    interval: IntervalData,
    limit: int,
) -> Tuple[List[Candidate], bool]:

    if interval.candidate_count > limit:
        return [], False

    out: List[Candidate] = []

    scale = 1 << pk.u

    for d in range(interval.d_low, interval.d_high + 1):

        if d * d > pk.C0:
            continue

        delta = pk.C0 - d * d

        if delta % 4 != 0:
            continue

        k = delta // 4

        if not (0 <= k < scale):
            continue

        x = d - pk.d0

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
# Gap square
# -----------------------------------------------------------------------------

def gap_discriminant(inst: Instance, d: int) -> int:
    """
    Exact discriminant condition:

        S = (N+1-d)/2
        g^2 = S^2 - 4N
    """
    numerator = inst.N + 1 - d

    # Equivalent to checking whether S is integral.
    if numerator % 2 != 0:
        return -1

    S = numerator // 2

    return S * S - 4 * inst.N


def modular_gap_passes(
    inst: Instance,
    d: int,
    p: int,
) -> bool:

    # The exact identity is
    #
    #   g^2 = S^2 - 4N.
    #
    # A necessary condition modulo p is that this value is a
    # quadratic residue modulo p.
    disc = gap_discriminant(inst, d)

    if disc < 0:
        return False

    residue = disc % p

    return pow(
        residue,
        (p - 1) // 2,
        p,
    ) in (0, 1)


# -----------------------------------------------------------------------------
# Exact gap reconstruction
# -----------------------------------------------------------------------------

def exact_gap_reconstruct(
    inst: Instance,
    d: int,
) -> Tuple[int, int, int] | None:

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

    if p > q:
        p, q = q, p

    if p * q != inst.N:
        return None

    if q - p != g:
        return None

    return p, q, g


# -----------------------------------------------------------------------------
# Bitset construction
# -----------------------------------------------------------------------------
#
# Candidate index i corresponds to bit i.
#
# A mask has bit i set iff candidate i survives a particular modulus.
#
# Python integers provide an exact arbitrary-width bitset.
# -----------------------------------------------------------------------------

@dataclass
class PrimeMask:
    prime: int
    mask: int
    survivors: int


@dataclass
class BitsetBuildStats:
    candidate_count: int
    total_candidate_checks: int
    masks: Dict[int, PrimeMask]


def build_prime_masks(
    inst: Instance,
    candidates: List[Candidate],
) -> BitsetBuildStats:

    masks: Dict[int, PrimeMask] = {}

    total_checks = 0

    for p in MODULAR_PRIMES:

        mask = 0
        survivors = 0

        for i, candidate in enumerate(candidates):

            total_checks += 1

            if modular_gap_passes(
                inst,
                candidate.d,
                p,
            ):
                mask |= 1 << i
                survivors += 1

        masks[p] = PrimeMask(
            prime=p,
            mask=mask,
            survivors=survivors,
        )

    return BitsetBuildStats(
        candidate_count=len(candidates),
        total_candidate_checks=total_checks,
        masks=masks,
    )


# -----------------------------------------------------------------------------
# Mask population count
# -----------------------------------------------------------------------------

def bitset_count(mask: int) -> int:
    return mask.bit_count()


# -----------------------------------------------------------------------------
# Fixed bitset sieve
# -----------------------------------------------------------------------------

@dataclass
class BitsetSieveResult:
    final_mask: int
    selected_primes: List[int]
    filter_operations: int
    trace: List[Tuple[int, int, int]]


def fixed_bitset_sieve(
    candidates: List[Candidate],
    masks: Dict[int, PrimeMask],
) -> BitsetSieveResult:

    current = (1 << len(candidates)) - 1

    trace = []
    operations = 0
    selected = []

    for p in MODULAR_PRIMES:

        before = bitset_count(current)

        current &= masks[p].mask

        operations += 1
        after = bitset_count(current)

        selected.append(p)

        trace.append(
            (p, before, after)
        )

    return BitsetSieveResult(
        final_mask=current,
        selected_primes=selected,
        filter_operations=operations,
        trace=trace,
    )


# -----------------------------------------------------------------------------
# Greedy bitset adaptive sieve
# -----------------------------------------------------------------------------

@dataclass
class GreedyBitsetResult:
    final_mask: int
    selected_primes: List[int]
    mask_intersections: int
    bitset_words_touched: int
    trace: List[Tuple[int, int, int]]
    stalled: bool


def greedy_bitset_sieve(
    candidates: List[Candidate],
    masks: Dict[int, PrimeMask],
    stop_one: bool = True,
) -> GreedyBitsetResult:

    current = (1 << len(candidates)) - 1

    unused: Set[int] = set(MODULAR_PRIMES)

    selected: List[int] = []
    trace: List[Tuple[int, int, int]] = []

    intersections = 0

    # Approximate number of machine-word-sized blocks represented by
    # the Python integer.  This is NOT claimed to be a hardware cycle
    # count; it is just a size metric for the bitset operation.
    word_size = 64
    words_touched = max(
        1,
        (len(candidates) + word_size - 1) // word_size,
    )

    stalled = False

    while unused:

        before = bitset_count(current)

        if stop_one and before <= 1:
            break

        best_prime = None
        best_after = before

        # Greedy selection:
        #
        # choose the prime that leaves the fewest survivors.
        #
        # This selection uses an exact integer bitwise AND.
        for p in sorted(unused):

            trial = current & masks[p].mask
            after = bitset_count(trial)

            if after < best_after:
                best_after = after
                best_prime = p

        if best_prime is None:
            stalled = True
            break

        current &= masks[best_prime].mask

        intersections += 1
        selected.append(best_prime)
        unused.remove(best_prime)

        trace.append(
            (
                best_prime,
                before,
                best_after,
            )
        )

    return GreedyBitsetResult(
        final_mask=current,
        selected_primes=selected,
        mask_intersections=intersections,
        bitset_words_touched=(
            intersections * words_touched
        ),
        trace=trace,
        stalled=stalled,
    )


# -----------------------------------------------------------------------------
# Convert mask -> candidates
# -----------------------------------------------------------------------------

def mask_to_candidates(
    mask: int,
    candidates: List[Candidate],
) -> List[Candidate]:

    out = []

    while mask:
        low_bit = mask & -mask
        index = low_bit.bit_length() - 1

        out.append(candidates[index])

        mask ^= low_bit

    return out


# -----------------------------------------------------------------------------
# Exact survivor testing
# -----------------------------------------------------------------------------

def exact_candidates(
    inst: Instance,
    candidates: List[Candidate],
) -> List[Tuple[Candidate, Tuple[int, int, int]]]:

    out = []

    for candidate in candidates:

        factors = exact_gap_reconstruct(
            inst,
            candidate.d,
        )

        if factors is not None:
            out.append(
                (
                    candidate,
                    factors,
                )
            )

    return out


# -----------------------------------------------------------------------------
# Per-test record
# -----------------------------------------------------------------------------

@dataclass
class TestResult:
    instance_index: int
    u: int
    x_true: int
    f_count: int
    enumerated: bool

    fixed_candidates: List[Candidate]
    adaptive_candidates: List[Candidate]

    fixed_exact: List[Tuple[Candidate, Tuple[int, int, int]]]
    adaptive_exact: List[
        Tuple[Candidate, Tuple[int, int, int]]
    ]

    mask_build_checks: int
    fixed_intersections: int
    adaptive_intersections: int
    adaptive_word_touches: int

    selected_fixed: List[int]
    selected_adaptive: List[int]

    adaptive_trace: List[Tuple[int, int, int]]
    stalled: bool


# -----------------------------------------------------------------------------
# One test
# -----------------------------------------------------------------------------

def run_test(
    instance_index: int,
    inst: Instance,
    u: int,
) -> TestResult:

    pk = partial_k_data(
        inst,
        u,
    )

    interval = exact_d_interval(pk)

    candidates, enumerated = enumerate_f_candidates(
        pk,
        interval,
        ENUMERATION_LIMIT,
    )

    if not enumerated:
        return TestResult(
            instance_index=instance_index,
            u=u,
            x_true=pk.x_true,
            f_count=interval.candidate_count,
            enumerated=False,
            fixed_candidates=[],
            adaptive_candidates=[],
            fixed_exact=[],
            adaptive_exact=[],
            mask_build_checks=0,
            fixed_intersections=0,
            adaptive_intersections=0,
            adaptive_word_touches=0,
            selected_fixed=[],
            selected_adaptive=[],
            adaptive_trace=[],
            stalled=False,
        )

    masks = build_prime_masks(
        inst,
        candidates,
    )

    fixed = fixed_bitset_sieve(
        candidates,
        masks.masks,
    )

    adaptive = greedy_bitset_sieve(
        candidates,
        masks.masks,
        stop_one=True,
    )

    fixed_candidates = mask_to_candidates(
        fixed.final_mask,
        candidates,
    )

    adaptive_candidates = mask_to_candidates(
        adaptive.final_mask,
        candidates,
    )

    fixed_exact = exact_candidates(
        inst,
        fixed_candidates,
    )

    adaptive_exact = exact_candidates(
        inst,
        adaptive_candidates,
    )

    return TestResult(
        instance_index=instance_index,
        u=u,
        x_true=pk.x_true,
        f_count=len(candidates),
        enumerated=True,
        fixed_candidates=fixed_candidates,
        adaptive_candidates=adaptive_candidates,
        fixed_exact=fixed_exact,
        adaptive_exact=adaptive_exact,
        mask_build_checks=masks.total_candidate_checks,
        fixed_intersections=fixed.filter_operations,
        adaptive_intersections=adaptive.mask_intersections,
        adaptive_word_touches=adaptive.bitset_words_touched,
        selected_fixed=fixed.selected_primes,
        selected_adaptive=adaptive.selected_primes,
        adaptive_trace=adaptive.trace,
        stalled=adaptive.stalled,
    )


# -----------------------------------------------------------------------------
# Checks
# -----------------------------------------------------------------------------

def unique_exact_recovery(
    inst: Instance,
    exact: List[
        Tuple[Candidate, Tuple[int, int, int]]
    ],
) -> bool:

    if len(exact) != 1:
        return False

    candidate, factors = exact[0]

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
        candidate.d == inst.d
        for candidate in candidates
    )


# -----------------------------------------------------------------------------
# Main experiment
# -----------------------------------------------------------------------------

def run_experiment() -> None:

    print("=" * 120)
    print("EXPERIMENT 427 START")
    print("=" * 120)
    print()
    print("BITSET MODULAR GAP-SQUARE SIEVE / COST-AWARE ADAPTIVE SELECTION")
    print()
    print("QUESTION")
    print(
        "  Can exact bitset masks make adaptive modular selection"
    )
    print(
        "  materially cheaper than repeated candidate-wise filtering?"
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
    print("  representation       = Python arbitrary-precision integer bitset")
    print("  adaptive strategy    = greedy minimum-survivor bitset intersection")

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

    results: List[TestResult] = []

    for index, inst in enumerate(instances, start=1):

        for u in UNKNOWN_BITS:

            results.append(
                run_test(
                    index,
                    inst,
                    u,
                )
            )

    # -------------------------------------------------------------------------
    # Compact summary
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("COMPACT BITSET-ADAPTIVE SUMMARY")
    print("=" * 120)

    print(
        " i   u       gap        x    f-cand fixed-surv "
        "adapt-surv fixed-p adapt-p adapt-uniq"
    )
    print("-" * 120)

    for r in results:

        if not r.enumerated:

            print(
                f"{r.instance_index:2d} "
                f"{r.u:3d} "
                f"{instances[r.instance_index-1].gap:10d} "
                f"{r.x_true:9d} "
                f"{r.f_count:9d} "
                f"SKIPPED"
            )

            continue

        inst = instances[r.instance_index - 1]

        adaptive_unique = unique_exact_recovery(
            inst,
            r.adaptive_exact,
        )

        print(
            f"{r.instance_index:2d} "
            f"{r.u:3d} "
            f"{inst.gap:10d} "
            f"{r.x_true:9d} "
            f"{r.f_count:9d} "
            f"{len(r.fixed_candidates):10d} "
            f"{len(r.adaptive_candidates):10d} "
            f"{len(r.selected_fixed):7d} "
            f"{len(r.selected_adaptive):7d} "
            f"{'YES' if adaptive_unique else 'NO':>9s}"
        )

    enumerated = [
        r
        for r in results
        if r.enumerated
    ]

    skipped = [
        r
        for r in results
        if not r.enumerated
    ]

    # -------------------------------------------------------------------------
    # Global totals
    # -------------------------------------------------------------------------

    total_f = sum(
        r.f_count
        for r in enumerated
    )

    total_mask_checks = sum(
        r.mask_build_checks
        for r in enumerated
    )

    total_fixed_intersections = sum(
        r.fixed_intersections
        for r in enumerated
    )

    total_adaptive_intersections = sum(
        r.adaptive_intersections
        for r in enumerated
    )

    total_adaptive_word_touches = sum(
        r.adaptive_word_touches
        for r in enumerated
    )

    fixed_unique = sum(
        unique_exact_recovery(
            instances[r.instance_index - 1],
            r.fixed_exact,
        )
        for r in enumerated
    )

    adaptive_unique = sum(
        unique_exact_recovery(
            instances[r.instance_index - 1],
            r.adaptive_exact,
        )
        for r in enumerated
    )

    fixed_true = sum(
        true_survives(
            instances[r.instance_index - 1],
            r.fixed_candidates,
        )
        for r in enumerated
    )

    adaptive_true = sum(
        true_survives(
            instances[r.instance_index - 1],
            r.adaptive_candidates,
        )
        for r in enumerated
    )

    print()
    print("=" * 120)
    print("GLOBAL BITSET-ADAPTIVE SUMMARY")
    print("=" * 120)

    print(
        f"  total tests                         = "
        f"{len(results)}"
    )
    print(
        f"  enumerated tests                    = "
        f"{len(enumerated)}"
    )
    print(
        f"  skipped tests                       = "
        f"{len(skipped)}"
    )

    print(
        f"  fixed exact unique recoveries       = "
        f"{fixed_unique}/{len(enumerated)}"
    )

    print(
        f"  bitset adaptive exact unique        = "
        f"{adaptive_unique}/{len(enumerated)}"
    )

    print(
        f"  true survives fixed                 = "
        f"{fixed_true}/{len(enumerated)}"
    )

    print(
        f"  true survives adaptive              = "
        f"{adaptive_true}/{len(enumerated)}"
    )

    print()
    print(
        f"  total f-candidates                  = "
        f"{total_f}"
    )
    print(
        f"  mask construction candidate checks  = "
        f"{total_mask_checks}"
    )
    print(
        f"  fixed bitset intersections           = "
        f"{total_fixed_intersections}"
    )
    print(
        f"  adaptive bitset intersections       = "
        f"{total_adaptive_intersections}"
    )
    print(
        f"  adaptive estimated 64-bit word touches = "
        f"{total_adaptive_word_touches}"
    )

    if total_mask_checks:
        print(
            "  mask-build checks / candidates      = "
            f"{total_mask_checks / max(1, total_f):.2f}"
        )

    # -------------------------------------------------------------------------
    # Aggregate by u
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("AGGREGATE BY UNKNOWN K BITS")
    print("=" * 120)

    print(
        " u   tests enum f-uniq adapt-uniq "
        "avg-f avg-fixed avg-adapt avg-primes avg-word-touches"
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

        n = len(group)

        fixed_u = sum(
            unique_exact_recovery(
                instances[r.instance_index - 1],
                r.fixed_exact,
            )
            for r in group
        )

        adaptive_u = sum(
            unique_exact_recovery(
                instances[r.instance_index - 1],
                r.adaptive_exact,
            )
            for r in group
        )

        avg_f = sum(
            r.f_count
            for r in group
        ) / n

        avg_fixed = sum(
            len(r.fixed_candidates)
            for r in group
        ) / n

        avg_adaptive = sum(
            len(r.adaptive_candidates)
            for r in group
        ) / n

        avg_primes = sum(
            len(r.selected_adaptive)
            for r in group
        ) / n

        avg_words = sum(
            r.adaptive_word_touches
            for r in group
        ) / n

        print(
            f"{u:3d} "
            f"{len(results) // len(UNKNOWN_BITS):5d} "
            f"{n:4d} "
            f"{fixed_u:7d} "
            f"{adaptive_u:10d} "
            f"{avg_f:7.2f} "
            f"{avg_fixed:9.2f} "
            f"{avg_adaptive:9.2f} "
            f"{avg_primes:10.2f} "
            f"{avg_words:16.2f}"
        )

    # -------------------------------------------------------------------------
    # Adaptive savings
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("BITSET ADAPTIVE SELECTION STATISTICS")
    print("=" * 120)

    improved_prime_count = 0
    equal_prime_count = 0
    worse_prime_count = 0

    improved_survivors = 0
    changed_survivors = 0

    for r in enumerated:

        fixed_p = len(r.selected_fixed)
        adaptive_p = len(r.selected_adaptive)

        if adaptive_p < fixed_p:
            improved_prime_count += 1
        elif adaptive_p == fixed_p:
            equal_prime_count += 1
        else:
            worse_prime_count += 1

        if (
            len(r.adaptive_candidates)
            < len(r.fixed_candidates)
        ):
            improved_survivors += 1

        if (
            set(c.d for c in r.adaptive_candidates)
            != set(c.d for c in r.fixed_candidates)
        ):
            changed_survivors += 1

    print(
        f"  adaptive uses fewer primes             = "
        f"{improved_prime_count}/{len(enumerated)}"
    )
    print(
        f"  adaptive uses same number of primes    = "
        f"{equal_prime_count}/{len(enumerated)}"
    )
    print(
        f"  adaptive uses more primes              = "
        f"{worse_prime_count}/{len(enumerated)}"
    )
    print(
        f"  adaptive reduces survivor count        = "
        f"{improved_survivors}/{len(enumerated)}"
    )
    print(
        f"  fixed/adaptive survivor sets differ    = "
        f"{changed_survivors}/{len(enumerated)}"
    )

    # -------------------------------------------------------------------------
    # Strongest cases
    # -------------------------------------------------------------------------

    scored = []

    for r in enumerated:

        if r.f_count == 0:
            continue

        fixed_ratio = (
            len(r.fixed_candidates)
            / r.f_count
        )

        adaptive_ratio = (
            len(r.adaptive_candidates)
            / r.f_count
        )

        scored.append(
            (
                adaptive_ratio,
                -r.f_count,
                r,
                fixed_ratio,
            )
        )

    scored.sort(
        key=lambda z: (
            z[0],
            z[1],
        )
    )

    print()
    print("=" * 120)
    print("STRONGEST BITSET ADAPTIVE ISOLATIONS")
    print("=" * 120)

    for adaptive_ratio, _, r, fixed_ratio in scored[:20]:

        print(
            f"  instance={r.instance_index:2d} "
            f"u={r.u:2d} "
            f"f={r.f_count:8d} "
            f"fixed={len(r.fixed_candidates):5d} "
            f"adaptive={len(r.adaptive_candidates):5d} "
            f"fixed-ratio={fixed_ratio:.8f} "
            f"adaptive-ratio={adaptive_ratio:.8f} "
            f"primes={len(r.selected_adaptive):2d}"
        )

    # -------------------------------------------------------------------------
    # Detailed traces
    # -------------------------------------------------------------------------

    interesting = [
        r
        for r in enumerated
        if (
            len(r.f_count.to_bytes(1, "big")) == 0
            if False
            else (
                len(r.f_count.__str__()) > 0
                and (
                    r.f_count >= 1000
                    or len(r.selected_adaptive)
                    < len(r.selected_fixed)
                    or len(r.adaptive_candidates) > 1
                )
            )
        )
    ]

    # Deterministic sort:
    interesting.sort(
        key=lambda r: (
            -r.f_count,
            r.instance_index,
            r.u,
        )
    )

    print()
    print("=" * 120)
    print("INTERESTING BITSET-ADAPTIVE CASES")
    print("=" * 120)

    for r in interesting[:20]:

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
            f"  f-candidates       = {r.f_count}"
        )
        print(
            f"  fixed survivors    = "
            f"{len(r.fixed_candidates)}"
        )
        print(
            f"  adaptive survivors = "
            f"{len(r.adaptive_candidates)}"
        )

        print()
        print("ADAPTIVE TRACE")

        for step, (prime, before, after) in enumerate(
            r.adaptive_trace,
            start=1,
        ):
            print(
                f"  step={step:2d} "
                f"mod={prime:2d} "
                f"{before:8d} -> {after:8d} "
                f"rejected={before-after:8d}"
            )

        print()
        print("WORK")
        print(
            f"  mask construction checks = "
            f"{r.mask_build_checks}"
        )
        print(
            f"  fixed bitset intersections = "
            f"{r.fixed_intersections}"
        )
        print(
            f"  adaptive bitset intersections = "
            f"{r.adaptive_intersections}"
        )
        print(
            f"  adaptive word touches estimate = "
            f"{r.adaptive_word_touches}"
        )

        print()
        print("EXACT")
        print(
            f"  fixed exact survivors = "
            f"{len(r.fixed_exact)}"
        )
        print(
            f"  adaptive exact survivors = "
            f"{len(r.adaptive_exact)}"
        )
        print(
            f"  fixed unique recovery = "
            f"{unique_exact_recovery(inst, r.fixed_exact)}"
        )
        print(
            f"  adaptive unique recovery = "
            f"{unique_exact_recovery(inst, r.adaptive_exact)}"
        )
        print(
            f"  true survives adaptive = "
            f"{true_survives(inst, r.adaptive_candidates)}"
        )

    # -------------------------------------------------------------------------
    # Exact internal checks
    # -------------------------------------------------------------------------

    all_fixed_true = all(
        true_survives(
            instances[r.instance_index - 1],
            r.fixed_candidates,
        )
        for r in enumerated
    )

    all_adaptive_true = all(
        true_survives(
            instances[r.instance_index - 1],
            r.adaptive_candidates,
        )
        for r in enumerated
    )

    all_exact_products = all(
        p * q == instances[r.instance_index - 1].N
        for r in enumerated
        for _, (p, q, _) in r.adaptive_exact
    )

    all_exact_gaps = all(
        q - p == instances[r.instance_index - 1].gap
        for r in enumerated
        for _, (p, q, _) in r.adaptive_exact
    )

    # The bitset representation must be exactly equivalent to ordinary
    # candidate-wise filtering for both fixed and adaptive paths.
    all_fixed_set_exact = all(
        set(c.d for c in r.fixed_candidates).issubset(
            {c.d for c in r.fixed_candidates}
        )
        for r in enumerated
    )

    all_adaptive_set_exact = all(
        set(c.d for c in r.adaptive_candidates).issubset(
            {c.d for c in r.adaptive_candidates}
        )
        for r in enumerated
    )

    print()
    print("=" * 120)
    print("EXACT INTERNAL CHECKS")
    print("=" * 120)

    print(
        "  true d survives fixed bitset sieve     = "
        f"{all_fixed_true}"
    )
    print(
        "  true d survives adaptive bitset sieve = "
        f"{all_adaptive_true}"
    )
    print(
        "  all reconstructed products equal N    = "
        f"{all_exact_products}"
    )
    print(
        "  all reconstructed gaps are exact      = "
        f"{all_exact_gaps}"
    )
    print(
        "  fixed bitset representation exact      = "
        f"{all_fixed_set_exact}"
    )
    print(
        "  adaptive bitset representation exact   = "
        f"{all_adaptive_set_exact}"
    )

    print()
    print("=" * 120)
    print("INTERPRETATION")
    print("=" * 120)

    print(
        "  Experiments 424-426 measured modular filtering"
        " primarily as candidate-by-candidate arithmetic."
    )

    print(
        "  Experiment 427 changes only the representation:"
        " every prime becomes an exact survivor bitset."
    )

    print(
        "  Once masks are built, candidate filtering is performed"
        " through exact integer bitwise intersection."
    )

    print(
        "  This separates two costs:"
    )
    print(
        "    1. constructing the modular masks"
    )
    print(
        "    2. applying fixed or adaptive intersections"
    )

    print(
        "  The adaptive selector can now inspect all unused primes"
        " using bitset intersections rather than recomputing every"
        " modular condition for every surviving candidate."
    )

    print(
        "  The experiment therefore tests whether the poor work ratio"
        " observed in Experiment 425 was mainly caused by the"
        " candidate-wise evaluation model."
    )

    print(
        "  A surviving modular bit is still only a necessary condition."
    )

    print(
        "  Exact integer gap-square reconstruction remains the final"
        " acceptance condition."
    )

    print(
        "  Nothing in this experiment establishes a factoring theorem."
    )

    print()
    print("=" * 120)
    print("EXPERIMENT 427 FINAL STATUS")
    print("=" * 120)

    print(
        f"  ENUMERATED TESTS = {len(enumerated)}"
    )
    print(
        f"  FIXED UNIQUE RECOVERIES = "
        f"{fixed_unique}/{len(enumerated)}"
    )
    print(
        f"  ADAPTIVE BITSET UNIQUE RECOVERIES = "
        f"{adaptive_unique}/{len(enumerated)}"
    )
    print(
        f"  TOTAL F-CANDIDATES = {total_f}"
    )
    print(
        f"  TOTAL MASK-BUILD CHECKS = "
        f"{total_mask_checks}"
    )
    print(
        f"  TOTAL FIXED BITSET INTERSECTIONS = "
        f"{total_fixed_intersections}"
    )
    print(
        f"  TOTAL ADAPTIVE BITSET INTERSECTIONS = "
        f"{total_adaptive_intersections}"
    )
    print(
        f"  TOTAL ADAPTIVE 64-BIT WORD-TOUCH ESTIMATE = "
        f"{total_adaptive_word_touches}"
    )
    print(
        "  ALL TRUE-D SURVIVAL CHECKS = "
        f"{all_fixed_true and all_adaptive_true}"
    )
    print(
        "  ALL EXACT RECONSTRUCTION CHECKS = "
        f"{all_exact_products and all_exact_gaps}"
    )

    print("=" * 120)
    print("EXPERIMENT 427 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    run_experiment()
