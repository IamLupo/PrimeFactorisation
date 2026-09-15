from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple


# =============================================================================
# EXPERIMENT 428
# LAZY BITSET PRIME-MASK CACHE / GREEDY GAP-SQUARE SIEVE
# =============================================================================
#
# Goal:
#
#   Compare eager bitset mask construction against lazy, cached mask
#   construction while preserving the exact same modular survivor sets.
#
# Rules:
#
#   exact integer arithmetic only
#   no resultants
#   no Groebner basis
#   no symbolic factorization
#   modular conditions are necessary only
#   exact gap-square reconstruction remains authoritative
#
# =============================================================================


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

WORD_SIZE = 64


# =============================================================================
# NUMBER-THEORETIC INSTANCE
# =============================================================================

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

    value = (
        -d * d
        - 3 * N * N
        + 6 * N
        + 1
    )

    assert value % 4 == 0

    K = value // 4

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

    return Instance(
        p=p,
        q=q,
        N=N,
        S=S,
        d=d,
        gap=gap,
        K=K,
    )


# =============================================================================
# PARTIAL-K DECOMPOSITION
# =============================================================================

@dataclass(frozen=True)
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

    scale = 1 << u

    K0 = (
        inst.K // scale
    ) * scale

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


# =============================================================================
# EXACT D INTERVAL
# =============================================================================

@dataclass(frozen=True)
class IntervalData:
    d_low: int
    d_high: int
    candidate_count: int


def exact_d_interval(
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

        if d_low * d_low < lo_sq:
            d_low += 1

    d_high = math.isqrt(
        hi_sq
    )

    return IntervalData(
        d_low=d_low,
        d_high=d_high,
        candidate_count=max(
            0,
            d_high - d_low + 1,
        ),
    )


# =============================================================================
# EXACT F CANDIDATES
# =============================================================================

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

    for d in range(
        interval.d_low,
        interval.d_high + 1,
    ):

        if d * d > pk.C0:
            continue

        delta = pk.C0 - d * d

        if delta % 4 != 0:
            continue

        k = delta // 4

        if not (
            0 <= k < scale
        ):
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


# =============================================================================
# GAP DISCRIMINANT
# =============================================================================

def gap_discriminant(
    inst: Instance,
    d: int,
) -> int:

    numerator = (
        inst.N
        + 1
        - d
    )

    if numerator % 2 != 0:
        return -1

    S = numerator // 2

    return (
        S * S
        - 4 * inst.N
    )


# =============================================================================
# MODULAR GAP CONDITION
# =============================================================================

def modular_gap_passes(
    inst: Instance,
    d: int,
    prime: int,
) -> bool:

    disc = gap_discriminant(
        inst,
        d,
    )

    if disc < 0:
        return False

    residue = disc % prime

    return pow(
        residue,
        (prime - 1) // 2,
        prime,
    ) in (0, 1)


# =============================================================================
# EXACT GAP RECONSTRUCTION
# =============================================================================

def exact_gap_reconstruct(
    inst: Instance,
    d: int,
) -> Optional[
    Tuple[int, int, int]
]:

    numerator = (
        inst.N
        + 1
        - d
    )

    if numerator % 2 != 0:
        return None

    S = numerator // 2

    disc = (
        S * S
        - 4 * inst.N
    )

    if disc < 0:
        return None

    g = math.isqrt(disc)

    if g * g != disc:
        return None

    if (
        (S - g) % 2 != 0
        or
        (S + g) % 2 != 0
    ):
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

    return (
        p,
        q,
        g,
    )


# =============================================================================
# BITSET HELPERS
# =============================================================================

def bit_count(
    mask: int,
) -> int:

    return mask.bit_count()


def full_mask(
    count: int,
) -> int:

    if count <= 0:
        return 0

    return (
        1 << count
    ) - 1


def word_count(
    candidate_count: int,
) -> int:

    if candidate_count <= 0:
        return 0

    return (
        candidate_count
        + WORD_SIZE
        - 1
    ) // WORD_SIZE


# =============================================================================
# LAZY PRIME MASK
# =============================================================================

@dataclass
class LazyPrimeMask:
    prime: int
    mask: int
    survivor_count: int
    candidate_checks: int


class LazyMaskCache:

    def __init__(
        self,
        inst: Instance,
        candidates: List[Candidate],
    ) -> None:

        self.inst = inst
        self.candidates = candidates

        self.cache: Dict[
            int,
            LazyPrimeMask,
        ] = {}

        self.total_mask_build_checks = 0
        self.total_masks_built = 0

    def get(
        self,
        prime: int,
    ) -> LazyPrimeMask:

        cached = self.cache.get(
            prime
        )

        if cached is not None:
            return cached

        mask = 0
        survivors = 0

        for index, candidate in enumerate(
            self.candidates
        ):

            self.total_mask_build_checks += 1

            if modular_gap_passes(
                self.inst,
                candidate.d,
                prime,
            ):

                mask |= (
                    1 << index
                )

                survivors += 1

        result = LazyPrimeMask(
            prime=prime,
            mask=mask,
            survivor_count=survivors,
            candidate_checks=len(
                self.candidates
            ),
        )

        self.cache[prime] = result
        self.total_masks_built += 1

        return result

    def built_primes(
        self,
    ) -> List[int]:

        return sorted(
            self.cache
        )

    def built_count(
        self,
    ) -> int:

        return self.total_masks_built


# =============================================================================
# FIXED LAZY SIEVE
# =============================================================================

@dataclass
class FixedLazyResult:
    final_mask: int
    selected_primes: List[int]
    intersections: int
    word_touches: int
    cache: LazyMaskCache


def fixed_lazy_sieve(
    inst: Instance,
    candidates: List[Candidate],
) -> FixedLazyResult:

    cache = LazyMaskCache(
        inst,
        candidates,
    )

    current = full_mask(
        len(candidates)
    )

    selected: List[int] = []

    intersections = 0

    words = word_count(
        len(candidates)
    )

    for prime in MODULAR_PRIMES:

        pmask = cache.get(
            prime
        )

        current &= pmask.mask

        selected.append(
            prime
        )

        intersections += 1

    return FixedLazyResult(
        final_mask=current,
        selected_primes=selected,
        intersections=intersections,
        word_touches=(
            intersections * words
        ),
        cache=cache,
    )


# =============================================================================
# ADAPTIVE LAZY SIEVE
# =============================================================================

@dataclass
class AdaptiveLazyResult:
    final_mask: int
    selected_primes: List[int]
    intersections: int
    word_touches: int
    selector_mask_probes: int
    selector_bitset_ands: int
    cache: LazyMaskCache
    trace: List[
        Tuple[
            int,
            int,
            int,
        ]
    ]
    stalled: bool


def adaptive_lazy_sieve(
    inst: Instance,
    candidates: List[Candidate],
) -> AdaptiveLazyResult:

    cache = LazyMaskCache(
        inst,
        candidates,
    )

    current = full_mask(
        len(candidates)
    )

    unused: Set[int] = set(
        MODULAR_PRIMES
    )

    selected: List[int] = []

    trace: List[
        Tuple[int, int, int]
    ] = []

    intersections = 0
    selector_probes = 0
    selector_ands = 0

    words = word_count(
        len(candidates)
    )

    word_touches = 0
    stalled = False

    while unused:

        before = bit_count(
            current
        )

        if before <= 1:
            break

        best_prime: Optional[int] = None
        best_after = before

        for prime in sorted(
            unused
        ):

            selector_probes += 1

            pmask = cache.get(
                prime
            )

            trial = (
                current
                & pmask.mask
            )

            selector_ands += 1

            after = bit_count(
                trial
            )

            if after < best_after:
                best_after = after
                best_prime = prime

        if best_prime is None:

            stalled = True

            break

        current &= cache.get(
            best_prime
        ).mask

        selected.append(
            best_prime
        )

        intersections += 1

        word_touches += words

        unused.remove(
            best_prime
        )

        trace.append(
            (
                best_prime,
                before,
                best_after,
            )
        )

    return AdaptiveLazyResult(
        final_mask=current,
        selected_primes=selected,
        intersections=intersections,
        word_touches=word_touches,
        selector_mask_probes=selector_probes,
        selector_bitset_ands=selector_ands,
        cache=cache,
        trace=trace,
        stalled=stalled,
    )


# =============================================================================
# BITSET -> CANDIDATES
# =============================================================================

def mask_to_candidates(
    mask: int,
    candidates: List[Candidate],
) -> List[Candidate]:

    result: List[Candidate] = []

    while mask:

        low_bit = (
            mask & -mask
        )

        index = (
            low_bit.bit_length()
            - 1
        )

        result.append(
            candidates[index]
        )

        mask ^= low_bit

    return result


# =============================================================================
# EXACT SURVIVORS
# =============================================================================

def exact_survivors(
    inst: Instance,
    candidates: List[Candidate],
) -> List[
    Tuple[
        Candidate,
        Tuple[int, int, int],
    ]
]:

    result = []

    for candidate in candidates:

        factors = exact_gap_reconstruct(
            inst,
            candidate.d,
        )

        if factors is not None:

            result.append(
                (
                    candidate,
                    factors,
                )
            )

    return result


# =============================================================================
# RECOVERY CHECKS
# =============================================================================

def unique_recovery(
    inst: Instance,
    exact: List[
        Tuple[
            Candidate,
            Tuple[int, int, int],
        ]
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
        c.d == inst.d
        for c in candidates
    )


# =============================================================================
# TEST RESULT
# =============================================================================

@dataclass
class TestResult:
    instance_index: int
    u: int
    x_true: int
    f_count: int
    enumerated: bool

    fixed_candidates: List[Candidate]
    adaptive_candidates: List[Candidate]

    fixed_exact: List[
        Tuple[
            Candidate,
            Tuple[int, int, int],
        ]
    ]

    adaptive_exact: List[
        Tuple[
            Candidate,
            Tuple[int, int, int],
        ]
    ]

    fixed_mask_build_checks: int
    fixed_masks_built: int

    adaptive_mask_build_checks: int
    adaptive_masks_built: int

    fixed_intersections: int
    adaptive_intersections: int

    adaptive_selector_probes: int
    adaptive_selector_bitset_ands: int

    fixed_word_touches: int
    adaptive_word_touches: int

    fixed_primes: List[int]
    adaptive_primes: List[int]

    adaptive_trace: List[
        Tuple[int, int, int]
    ]

    stalled: bool


# =============================================================================
# ONE TEST
# =============================================================================

def run_test(
    instance_index: int,
    inst: Instance,
    u: int,
) -> TestResult:

    pk = partial_k_data(
        inst,
        u,
    )

    interval = exact_d_interval(
        pk
    )

    candidates, enumerated = (
        enumerate_f_candidates(
            pk,
            interval,
            ENUMERATION_LIMIT,
        )
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
            fixed_mask_build_checks=0,
            fixed_masks_built=0,
            adaptive_mask_build_checks=0,
            adaptive_masks_built=0,
            fixed_intersections=0,
            adaptive_intersections=0,
            adaptive_selector_probes=0,
            adaptive_selector_bitset_ands=0,
            fixed_word_touches=0,
            adaptive_word_touches=0,
            fixed_primes=[],
            adaptive_primes=[],
            adaptive_trace=[],
            stalled=False,
        )

    fixed = fixed_lazy_sieve(
        inst,
        candidates,
    )

    adaptive = adaptive_lazy_sieve(
        inst,
        candidates,
    )

    fixed_candidates = (
        mask_to_candidates(
            fixed.final_mask,
            candidates,
        )
    )

    adaptive_candidates = (
        mask_to_candidates(
            adaptive.final_mask,
            candidates,
        )
    )

    fixed_exact = exact_survivors(
        inst,
        fixed_candidates,
    )

    adaptive_exact = exact_survivors(
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
        fixed_mask_build_checks=(
            fixed.cache.total_mask_build_checks
        ),
        fixed_masks_built=(
            fixed.cache.built_count()
        ),
        adaptive_mask_build_checks=(
            adaptive.cache.total_mask_build_checks
        ),
        adaptive_masks_built=(
            adaptive.cache.built_count()
        ),
        fixed_intersections=fixed.intersections,
        adaptive_intersections=(
            adaptive.intersections
        ),
        adaptive_selector_probes=(
            adaptive.selector_mask_probes
        ),
        adaptive_selector_bitset_ands=(
            adaptive.selector_bitset_ands
        ),
        fixed_word_touches=(
            fixed.word_touches
        ),
        adaptive_word_touches=(
            adaptive.word_touches
        ),
        fixed_primes=(
            fixed.selected_primes
        ),
        adaptive_primes=(
            adaptive.selected_primes
        ),
        adaptive_trace=(
            adaptive.trace
        ),
        stalled=adaptive.stalled,
    )


# =============================================================================
# MAIN EXPERIMENT
# =============================================================================

def run_experiment() -> None:

    print("=" * 120)
    print("EXPERIMENT 428 START")
    print("=" * 120)
    print()

    print(
        "LAZY BITSET PRIME-MASK CACHE / "
        "GREEDY GAP-SQUARE SIEVE"
    )

    print()

    print("QUESTION")

    print(
        "  Can lazy mask construction reduce "
        "the dominant mask-building cost"
    )

    print(
        "  observed in Experiment 427 while "
        "preserving the same exact"
    )

    print(
        "  modular survivor sets and exact "
        "recoveries?"
    )

    print()

    print("=" * 120)
    print("CONFIGURATION")
    print("=" * 120)

    print("  instances             = 8")
    print(
        f"  K bits                = "
        f"{UNKNOWN_BITS}"
    )
    print(
        f"  modular primes        = "
        f"{MODULAR_PRIMES}"
    )
    print(
        f"  enumeration max       = "
        f"{ENUMERATION_LIMIT}"
    )
    print(
        "  representation        = "
        "lazy exact integer bitset masks"
    )
    print(
        "  adaptive selector     = "
        "greedy minimum-survivor"
    )
    print(
        "  mask construction     = "
        "lazy + cached"
    )
    print(
        f"  estimated word size   = "
        f"{WORD_SIZE} bits"
    )

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

    for index, inst in enumerate(
        instances,
        start=1,
    ):

        for u in UNKNOWN_BITS:

            results.append(
                run_test(
                    index,
                    inst,
                    u,
                )
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
    # COMPACT SUMMARY
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("COMPACT LAZY-BITSET SUMMARY")
    print("=" * 120)

    print(
        " i   u       gap        x    f-cand "
        "fixed adapt fixed-masks adapt-masks "
        "adapt-p unique"
    )

    print("-" * 120)

    for r in results:

        inst = instances[
            r.instance_index - 1
        ]

        if not r.enumerated:

            print(
                f"{r.instance_index:2d} "
                f"{r.u:3d} "
                f"{inst.gap:10d} "
                f"{r.x_true:9d} "
                f"{r.f_count:9d} "
                "SKIPPED"
            )

            continue

        adaptive_unique = (
            unique_recovery(
                inst,
                r.adaptive_exact,
            )
        )

        print(
            f"{r.instance_index:2d} "
            f"{r.u:3d} "
            f"{inst.gap:10d} "
            f"{r.x_true:9d} "
            f"{r.f_count:9d} "
            f"{len(r.fixed_candidates):6d} "
            f"{len(r.adaptive_candidates):5d} "
            f"{r.fixed_masks_built:11d} "
            f"{r.adaptive_masks_built:11d} "
            f"{len(r.adaptive_primes):8d} "
            f"{'YES' if adaptive_unique else 'NO':>6s}"
        )

    # -------------------------------------------------------------------------
    # GLOBAL TOTALS
    # -------------------------------------------------------------------------

    total_f = sum(
        r.f_count
        for r in enumerated
    )

    eager_equivalent_checks = (
        total_f
        * len(MODULAR_PRIMES)
    )

    fixed_lazy_checks = sum(
        r.fixed_mask_build_checks
        for r in enumerated
    )

    adaptive_lazy_checks = sum(
        r.adaptive_mask_build_checks
        for r in enumerated
    )

    fixed_masks = sum(
        r.fixed_masks_built
        for r in enumerated
    )

    adaptive_masks = sum(
        r.adaptive_masks_built
        for r in enumerated
    )

    fixed_intersections = sum(
        r.fixed_intersections
        for r in enumerated
    )

    adaptive_intersections = sum(
        r.adaptive_intersections
        for r in enumerated
    )

    selector_probes = sum(
        r.adaptive_selector_probes
        for r in enumerated
    )

    selector_ands = sum(
        r.adaptive_selector_bitset_ands
        for r in enumerated
    )

    fixed_words = sum(
        r.fixed_word_touches
        for r in enumerated
    )

    adaptive_words = sum(
        r.adaptive_word_touches
        for r in enumerated
    )

    fixed_unique = sum(
        unique_recovery(
            instances[
                r.instance_index - 1
            ],
            r.fixed_exact,
        )
        for r in enumerated
    )

    adaptive_unique = sum(
        unique_recovery(
            instances[
                r.instance_index - 1
            ],
            r.adaptive_exact,
        )
        for r in enumerated
    )

    fixed_true = all(
        true_survives(
            instances[
                r.instance_index - 1
            ],
            r.fixed_candidates,
        )
        for r in enumerated
    )

    adaptive_true = all(
        true_survives(
            instances[
                r.instance_index - 1
            ],
            r.adaptive_candidates,
        )
        for r in enumerated
    )

    fixed_adaptive_same = all(
        {
            c.d
            for c in r.fixed_candidates
        }
        ==
        {
            c.d
            for c in r.adaptive_candidates
        }
        for r in enumerated
    )

    all_adaptive_products = all(
        p * q
        ==
        instances[
            r.instance_index - 1
        ].N
        for r in enumerated
        for _, (p, q, _) in r.adaptive_exact
    )

    all_adaptive_gaps = all(
        q - p
        ==
        instances[
            r.instance_index - 1
        ].gap
        for r in enumerated
        for _, (p, q, _) in r.adaptive_exact
    )

    all_selected_valid = all(
        set(
            r.adaptive_primes
        ).issubset(
            set(MODULAR_PRIMES)
        )
        for r in enumerated
    )

    print()
    print("=" * 120)
    print("GLOBAL LAZY-MASK SUMMARY")
    print("=" * 120)

    print(
        f"  total tests                     = "
        f"{len(results)}"
    )

    print(
        f"  enumerated tests                = "
        f"{len(enumerated)}"
    )

    print(
        f"  skipped tests                   = "
        f"{len(skipped)}"
    )

    print(
        f"  fixed unique recoveries         = "
        f"{fixed_unique}/{len(enumerated)}"
    )

    print(
        f"  adaptive unique recoveries      = "
        f"{adaptive_unique}/{len(enumerated)}"
    )

    print(
        f"  true survives fixed             = "
        f"{fixed_true}"
    )

    print(
        f"  true survives adaptive          = "
        f"{adaptive_true}"
    )

    print()

    print(
        f"  total f-candidates              = "
        f"{total_f}"
    )

    print(
        f"  eager mask-build checks         = "
        f"{eager_equivalent_checks}"
    )

    print(
        f"  fixed lazy mask-build checks    = "
        f"{fixed_lazy_checks}"
    )

    print(
        f"  adaptive lazy mask-build checks = "
        f"{adaptive_lazy_checks}"
    )

    print(
        f"  fixed masks actually built      = "
        f"{fixed_masks}"
    )

    print(
        f"  adaptive masks actually built   = "
        f"{adaptive_masks}"
    )

    print(
        f"  fixed bitset intersections      = "
        f"{fixed_intersections}"
    )

    print(
        f"  adaptive bitset intersections   = "
        f"{adaptive_intersections}"
    )

    print(
        f"  adaptive selector mask probes   = "
        f"{selector_probes}"
    )

    print(
        f"  adaptive selector bitset ANDs   = "
        f"{selector_ands}"
    )

    print(
        f"  fixed 64-bit word touches       = "
        f"{fixed_words}"
    )

    print(
        f"  adaptive 64-bit word touches    = "
        f"{adaptive_words}"
    )

    if eager_equivalent_checks:

        print(
            "  fixed lazy / eager mask-check ratio = "
            f"{fixed_lazy_checks / eager_equivalent_checks:.8f}"
        )

        print(
            "  adaptive lazy / eager mask-check ratio = "
            f"{adaptive_lazy_checks / eager_equivalent_checks:.8f}"
        )

    if fixed_lazy_checks:

        print(
            "  adaptive / fixed lazy mask-check ratio = "
            f"{adaptive_lazy_checks / fixed_lazy_checks:.8f}"
        )

    # -------------------------------------------------------------------------
    # AGGREGATE BY U
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("AGGREGATE BY UNKNOWN K BITS")
    print("=" * 120)

    print(
        " u   tests enum avg-f fixed-masks adapt-masks "
        "fixed-checks adapt-checks fixed-uniq adapt-uniq"
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

        avg_f = (
            sum(
                r.f_count
                for r in group
            )
            / n
        )

        avg_fixed_masks = (
            sum(
                r.fixed_masks_built
                for r in group
            )
            / n
        )

        avg_adaptive_masks = (
            sum(
                r.adaptive_masks_built
                for r in group
            )
            / n
        )

        avg_fixed_checks = (
            sum(
                r.fixed_mask_build_checks
                for r in group
            )
            / n
        )

        avg_adaptive_checks = (
            sum(
                r.adaptive_mask_build_checks
                for r in group
            )
            / n
        )

        fixed_u = sum(
            unique_recovery(
                instances[
                    r.instance_index - 1
                ],
                r.fixed_exact,
            )
            for r in group
        )

        adaptive_u = sum(
            unique_recovery(
                instances[
                    r.instance_index - 1
                ],
                r.adaptive_exact,
            )
            for r in group
        )

        print(
            f"{u:3d} "
            f"{len(group):7d} "
            f"{len(group):4d} "
            f"{avg_f:9.2f} "
            f"{avg_fixed_masks:11.2f} "
            f"{avg_adaptive_masks:11.2f} "
            f"{avg_fixed_checks:12.2f} "
            f"{avg_adaptive_checks:13.2f} "
            f"{fixed_u:10d} "
            f"{adaptive_u:10d}"
        )

    # -------------------------------------------------------------------------
    # STRONGEST LAZY SAVINGS
    # -------------------------------------------------------------------------

    savings = []

    for r in enumerated:

        eager = (
            r.f_count
            * len(MODULAR_PRIMES)
        )

        adaptive = (
            r.adaptive_mask_build_checks
        )

        ratio = (
            adaptive / eager
            if eager
            else 0.0
        )

        saved = (
            eager - adaptive
        )

        savings.append(
            (
                ratio,
                -saved,
                r,
                eager,
                adaptive,
            )
        )

    savings.sort(
        key=lambda z: (
            z[0],
            z[1],
        )
    )

    print()
    print("=" * 120)
    print("STRONGEST LAZY MASK SAVINGS")
    print("=" * 120)

    for (
        ratio,
        _,
        r,
        eager,
        adaptive,
    ) in savings[:20]:

        print(
            f"  instance={r.instance_index:2d} "
            f"u={r.u:2d} "
            f"f={r.f_count:8d} "
            f"eager={eager:9d} "
            f"lazy={adaptive:9d} "
            f"ratio={ratio:.8f} "
            f"masks={r.adaptive_masks_built:2d}"
        )

    # -------------------------------------------------------------------------
    # INTERESTING CASES
    # -------------------------------------------------------------------------

    interesting = [
        r
        for r in enumerated
        if (
            r.f_count >= 1000
            or
            r.adaptive_masks_built
            < len(MODULAR_PRIMES)
        )
    ]

    interesting.sort(
        key=lambda r: (
            -(
                r.f_count
                * len(MODULAR_PRIMES)
                - r.adaptive_mask_build_checks
            ),
            -r.f_count,
        )
    )

    print()
    print("=" * 120)
    print("INTERESTING LAZY-CACHE CASES")
    print("=" * 120)

    for r in interesting[:20]:

        inst = instances[
            r.instance_index - 1
        ]

        eager = (
            r.f_count
            * len(MODULAR_PRIMES)
        )

        print()

        print(
            f"DETAIL: INSTANCE "
            f"{r.instance_index}, "
            f"UNKNOWN K BITS = {r.u}"
        )

        print("-" * 120)

        print(
            f"p={inst.p} "
            f"q={inst.q} "
            f"gap={inst.gap} "
            f"N={inst.N}"
        )

        print(
            f"d={inst.d} "
            f"K={inst.K}"
        )

        print()
        print("CANDIDATES")

        print(
            f"  f-candidates       = "
            f"{r.f_count}"
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
        print("MASK CONSTRUCTION")

        print(
            f"  eager equivalent checks = "
            f"{eager}"
        )

        print(
            f"  adaptive lazy checks    = "
            f"{r.adaptive_mask_build_checks}"
        )

        print(
            f"  masks actually built    = "
            f"{r.adaptive_masks_built}"
        )

        print(
            f"  masks not built         = "
            f"{len(MODULAR_PRIMES) - r.adaptive_masks_built}"
        )

        print()
        print("BUILT PRIME MASKS")

        print(
            f"  {r.adaptive_primes}"
        )

        print()
        print("ADAPTIVE PRIME TRACE")

        for step, (
            prime,
            before,
            after,
        ) in enumerate(
            r.adaptive_trace,
            start=1,
        ):

            print(
                f"  step={step:2d} "
                f"mod={prime:2d} "
                f"{before:8d} -> "
                f"{after:8d} "
                f"rejected="
                f"{before-after:8d}"
            )

        print()
        print("SELECTOR WORK")

        print(
            f"  selector mask probes     = "
            f"{r.adaptive_selector_probes}"
        )

        print(
            f"  selector bitset ANDs     = "
            f"{r.adaptive_selector_bitset_ands}"
        )

        print(
            f"  final adaptive ANDs      = "
            f"{r.adaptive_intersections}"
        )

        print(
            f"  adaptive word touches    = "
            f"{r.adaptive_word_touches}"
        )

        print(
            f"  stalled                  = "
            f"{r.stalled}"
        )

        print()
        print("EXACT")

        print(
            f"  fixed exact survivors    = "
            f"{len(r.fixed_exact)}"
        )

        print(
            f"  adaptive exact survivors = "
            f"{len(r.adaptive_exact)}"
        )

        print(
            f"  fixed unique recovery    = "
            f"{unique_recovery(inst, r.fixed_exact)}"
        )

        print(
            f"  adaptive unique recovery = "
            f"{unique_recovery(inst, r.adaptive_exact)}"
        )

        print(
            f"  true survives adaptive   = "
            f"{true_survives(inst, r.adaptive_candidates)}"
        )

    # -------------------------------------------------------------------------
    # EXACT INTERNAL CHECKS
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("EXACT INTERNAL CHECKS")
    print("=" * 120)

    print(
        "  true d survives fixed lazy sieve    = "
        f"{fixed_true}"
    )

    print(
        "  true d survives adaptive lazy sieve = "
        f"{adaptive_true}"
    )

    print(
        "  fixed/adaptive survivor sets equal  = "
        f"{fixed_adaptive_same}"
    )

    print(
        "  all reconstructed products equal N = "
        f"{all_adaptive_products}"
    )

    print(
        "  all reconstructed gaps exact        = "
        f"{all_adaptive_gaps}"
    )

    print(
        "  all selected primes valid           = "
        f"{all_selected_valid}"
    )

    # -------------------------------------------------------------------------
    # INTERPRETATION
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("INTERPRETATION")
    print("=" * 120)

    print(
        "  Experiment 427 showed that exact integer bitsets make adaptive"
    )
    print(
        "  modular intersection substantially cheaper than repeated"
    )
    print(
        "  candidate-wise modular evaluation."
    )

    print()

    print(
        "  Experiment 427 still built every prime mask eagerly, producing"
    )
    print(
        "  one modular mask-build check for every candidate-prime pair."
    )

    print()

    print(
        "  Experiment 428 makes those masks lazy and cached."
    )

    print(
        "  A prime mask is constructed only when the adaptive selector"
    )

    print(
        "  actually inspects that prime. Once constructed, it is reused."
    )

    print()

    print(
        "  The decisive implementation metric is adaptive lazy mask-build"
    )

    print(
        "  checks versus the eager baseline:"
    )

    print(
        "    candidates × number_of_primes"
    )

    print()

    print(
        "  The exact survivor set is unchanged because laziness changes"
    )

    print(
        "  only when a necessary-condition mask is constructed, not the"
    )

    print(
        "  mathematical condition represented by the mask."
    )

    print()

    print(
        "  This remains a candidate-isolation experiment, not a factoring"
    )

    print(
        "  theorem. The modular sieve is necessary-only, and exact integer"
    )

    print(
        "  gap-square reconstruction remains authoritative."
    )

    # -------------------------------------------------------------------------
    # FINAL STATUS
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("EXPERIMENT 428 FINAL STATUS")
    print("=" * 120)

    print(
        f"  ENUMERATED TESTS = "
        f"{len(enumerated)}"
    )

    print(
        f"  FIXED UNIQUE RECOVERIES = "
        f"{fixed_unique}/{len(enumerated)}"
    )

    print(
        f"  ADAPTIVE UNIQUE RECOVERIES = "
        f"{adaptive_unique}/{len(enumerated)}"
    )

    print(
        f"  TOTAL F-CANDIDATES = "
        f"{total_f}"
    )

    print(
        f"  EAGER EQUIVALENT MASK CHECKS = "
        f"{eager_equivalent_checks}"
    )

    print(
        f"  FIXED LAZY MASK CHECKS = "
        f"{fixed_lazy_checks}"
    )

    print(
        f"  ADAPTIVE LAZY MASK CHECKS = "
        f"{adaptive_lazy_checks}"
    )

    if eager_equivalent_checks:

        print(
            f"  ADAPTIVE/EAGER MASK-CHECK RATIO = "
            f"{adaptive_lazy_checks / eager_equivalent_checks:.8f}"
        )

    print(
        f"  ADAPTIVE MASKS BUILT = "
        f"{adaptive_masks}"
    )

    print(
        f"  ADAPTIVE BITSET INTERSECTIONS = "
        f"{adaptive_intersections}"
    )

    print(
        f"  ALL TRUE-D SURVIVAL CHECKS = "
        f"{fixed_true and adaptive_true}"
    )

    print(
        f"  ALL EXACT RECONSTRUCTION CHECKS = "
        f"{all_adaptive_products and all_adaptive_gaps}"
    )

    print("=" * 120)
    print("EXPERIMENT 428 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    run_experiment()