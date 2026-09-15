from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Tuple, Dict


# =============================================================================
# EXPERIMENT 425
# ADAPTIVE MODULAR GAP-SQUARE SIEVE / EXACT CANDIDATE-ISOLATION AUDIT
# =============================================================================
#
# Goal
# ----
#
# Experiment 424 used a fixed sequence of small primes:
#
#   [3,5,7,11,13,17,19,23,29,31,37,41,43,47]
#
# for every candidate set.
#
# Experiment 425 asks:
#
#   Can we choose the next prime adaptively, based on which prime gives
#   the largest immediate reduction of the current candidate set?
#
# For each partial-K instance:
#
#   1. Build all exact f-candidates.
#   2. Repeatedly select the unused prime that minimizes the number of
#      surviving candidates.
#   3. Preserve the invariant that the true d survives every modular step.
#   4. Stop when:
#        - one candidate remains,
#        - no unused prime improves the set,
#        - or all primes are exhausted.
#   5. Perform the exact gap-square test only on the final survivors.
#
# We compare:
#
#   fixed sieve:
#       always use the complete configured prime list
#
#   adaptive sieve:
#       choose the strongest next prime at every stage
#
# Measurements:
#
#   - initial f-candidates
#   - final adaptive survivors
#   - number of primes actually used
#   - number of modular candidate evaluations
#   - exact square-root evaluations
#   - whether adaptive sieve reaches uniqueness
#   - whether the true d survives every stage
#
# Rules
# -----
#
#   exact integer arithmetic only
#   no floating point
#   no resultants
#   no Groebner basis
#   no symbolic factorization
#   modular conditions are necessary only
#
# This is an exact candidate-filter audit, not a factoring theorem.
# =============================================================================


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

UNKNOWN_BITS = [20, 24, 28, 32, 36, 40, 44, 48]

# Same prime pool as Experiment 424.
MODULI = [
    3, 5, 7, 11, 13, 17, 19, 23,
    29, 31, 37, 41, 43, 47,
]

ENUMERATION_LIMIT = 250_000

# Adaptive sieve safeguards.
#
# A prime is considered useful if it strictly reduces the candidate set.
# If no remaining prime reduces the set, the adaptive phase stops.
#
# The true candidate is always expected to survive.
STOP_ON_ONE = True

# Output:
# compact by default, with detailed adaptive traces only for interesting cases.
SHOW_DETAILED_CASES = True

# A case is interesting when:
#   - f is not already unique, OR
#   - adaptive uses fewer primes than the complete sieve, OR
#   - adaptive does not reach uniqueness, OR
#   - there is a notable reduction.
DETAIL_MIN_INITIAL_CANDIDATES = 2


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

    K_num = -d * d - 3 * N * N + 6 * N + 1

    assert K_num % 4 == 0

    K = K_num // 4

    # Original identities.
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
# Partial-K decomposition
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
    step = 1 << u

    K0 = (inst.K // step) * step
    k = inst.K - K0

    assert 0 <= k < step
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
# Exact f-interval
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class IntervalData:
    d_low: int
    d_high: int
    x_low: int
    x_high: int
    candidate_count: int


def exact_d_interval(pk: PartialK) -> IntervalData:
    step = 1 << pk.u

    lo_sq = pk.C0 - 4 * (step - 1)
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

    candidate_count = max(0, d_high - d_low + 1)

    return IntervalData(
        d_low=d_low,
        d_high=d_high,
        x_low=x_low,
        x_high=x_high,
        candidate_count=candidate_count,
    )


# -----------------------------------------------------------------------------
# Exact f-candidate enumeration
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class CandidateSet:
    values: List[int]
    enumerated: bool


def enumerate_f_candidates(
    pk: PartialK,
    interval: IntervalData,
    limit: int,
) -> CandidateSet:

    if interval.candidate_count > limit:
        return CandidateSet([], False)

    out: List[int] = []
    step = 1 << pk.u

    for d in range(interval.d_low, interval.d_high + 1):
        delta = pk.C0 - d * d

        if delta < 0:
            continue

        if delta & 3:
            continue

        k = delta // 4

        if not (0 <= k < step):
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

        out.append(d)

    return CandidateSet(out, True)


# -----------------------------------------------------------------------------
# Gap-square condition
# -----------------------------------------------------------------------------
#
# For a candidate d:
#
#   S = (N + 1 - d)/2
#
# and
#
#   g^2 = S^2 - 4N.
#
# The modular sieve uses:
#
#   S^2 - 4N  is a quadratic residue modulo l.
#
# This is a necessary condition for an exact integer square.
# -----------------------------------------------------------------------------

def implied_S_if_integral(N: int, d: int) -> Tuple[bool, int]:
    numer = N + 1 - d

    if numer & 1:
        return False, 0

    return True, numer // 2


def gap_square_argument(N: int, d: int) -> Tuple[bool, int]:
    ok, S = implied_S_if_integral(N, d)

    if not ok:
        return False, 0

    value = S * S - 4 * N

    if value < 0:
        return False, value

    return True, value


# -----------------------------------------------------------------------------
# Modular residue precomputation
# -----------------------------------------------------------------------------

def quadratic_residue_table(modulus: int) -> List[bool]:
    table = [False] * modulus

    for x in range(modulus):
        table[(x * x) % modulus] = True

    return table


def candidate_passes_modulus(
    inst: Instance,
    d: int,
    modulus: int,
    residue_table: List[bool],
) -> bool:

    ok, S = implied_S_if_integral(inst.N, d)

    if not ok:
        return False

    # g^2 = S^2 - 4N
    value_mod = (
        (S % modulus) * (S % modulus)
        - 4 * (inst.N % modulus)
    ) % modulus

    return residue_table[value_mod]


# -----------------------------------------------------------------------------
# Exact gap-square reconstruction
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class GapReconstruction:
    d: int
    S: int
    gap: int
    p: int
    q: int


def reconstruct_factors_from_d(
    inst: Instance,
    d: int,
) -> GapReconstruction | None:

    ok, gap_sq = gap_square_argument(inst.N, d)

    if not ok:
        return None

    g = math.isqrt(gap_sq)

    if g * g != gap_sq:
        return None

    ok, S = implied_S_if_integral(inst.N, d)

    if not ok:
        return None

    if g < 0:
        return None

    if (S - g) & 1:
        return None

    if (S + g) & 1:
        return None

    p = (S - g) // 2
    q = (S + g) // 2

    if p <= 0 or q <= 0:
        return None

    if p * q != inst.N:
        return None

    if p + q != S:
        return None

    if q - p != g:
        return None

    return GapReconstruction(
        d=d,
        S=S,
        gap=g,
        p=p,
        q=q,
    )


# -----------------------------------------------------------------------------
# Fixed complete modular sieve
# -----------------------------------------------------------------------------

@dataclass
class FixedSieveResult:
    survivors: List[int]
    modular_evaluations: int
    true_survives: bool
    per_prime_counts: List[Tuple[int, int]]


def fixed_modular_sieve(
    inst: Instance,
    candidates: List[int],
    moduli: List[int],
    residue_tables: Dict[int, List[bool]],
) -> FixedSieveResult:

    current = list(candidates)
    evaluations = 0
    trace: List[Tuple[int, int]] = []

    for modulus in moduli:

        next_values: List[int] = []

        table = residue_tables[modulus]

        for d in current:
            evaluations += 1

            if candidate_passes_modulus(
                inst,
                d,
                modulus,
                table,
            ):
                next_values.append(d)

        current = next_values
        trace.append((modulus, len(current)))

        if len(current) <= 1 and STOP_ON_ONE:
            break

    return FixedSieveResult(
        survivors=current,
        modular_evaluations=evaluations,
        true_survives=(inst.d in current),
        per_prime_counts=trace,
    )


# -----------------------------------------------------------------------------
# Adaptive modular sieve
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class AdaptiveStep:
    modulus: int
    before: int
    after: int
    rejected: int


@dataclass
class AdaptiveSieveResult:
    survivors: List[int]
    steps: List[AdaptiveStep]
    modular_evaluations: int
    true_survives: bool
    stalled: bool


def adaptive_modular_sieve(
    inst: Instance,
    candidates: List[int],
    moduli: List[int],
    residue_tables: Dict[int, List[bool]],
) -> AdaptiveSieveResult:

    current = list(candidates)
    unused = list(moduli)

    steps: List[AdaptiveStep] = []
    evaluations = 0
    stalled = False

    while len(current) > 1 and unused:

        best_modulus = None
        best_survivors: List[int] | None = None
        best_count = len(current)

        for modulus in unused:

            table = residue_tables[modulus]
            trial: List[int] = []

            for d in current:
                if candidate_passes_modulus(
                    inst,
                    d,
                    modulus,
                    table,
                ):
                    trial.append(d)

            # Every complete trial costs len(current) modular evaluations.
            # We account for all trials below.
            evaluations += len(current)

            if len(trial) < best_count:
                best_count = len(trial)
                best_modulus = modulus
                best_survivors = trial

        if best_modulus is None or best_survivors is None:
            stalled = True
            break

        before = len(current)
        after = len(best_survivors)

        # Safety invariant: true d must survive.
        if inst.d not in best_survivors:
            raise AssertionError(
                "Adaptive modular sieve rejected the true d."
            )

        current = best_survivors

        steps.append(
            AdaptiveStep(
                modulus=best_modulus,
                before=before,
                after=after,
                rejected=before - after,
            )
        )

        unused.remove(best_modulus)

        if STOP_ON_ONE and len(current) == 1:
            break

    return AdaptiveSieveResult(
        survivors=current,
        steps=steps,
        modular_evaluations=evaluations,
        true_survives=(inst.d in current),
        stalled=stalled,
    )


# -----------------------------------------------------------------------------
# Important optimization diagnostic
# -----------------------------------------------------------------------------
#
# The adaptive algorithm above evaluates every remaining prime against the
# current set to find the best next prime.
#
# This is deliberately measured separately from the actual number of
# "selected" primes.  Thus we distinguish:
#
#   selected_primes
#   modular_evaluations_for_selection
#
# A production implementation could cache candidate residue signatures.
# This experiment first measures the mathematical selection strategy itself.
# -----------------------------------------------------------------------------

@dataclass
class ExactEvaluationResult:
    reconstructed: List[GapReconstruction]
    exact_square_evaluations: int
    true_reconstructed: bool


def exact_gap_filter(
    inst: Instance,
    candidates: List[int],
) -> ExactEvaluationResult:

    reconstructed: List[GapReconstruction] = []

    for d in candidates:
        result = reconstruct_factors_from_d(inst, d)

        if result is not None:
            reconstructed.append(result)

    return ExactEvaluationResult(
        reconstructed=reconstructed,
        exact_square_evaluations=len(candidates),
        true_reconstructed=any(
            item.d == inst.d
            for item in reconstructed
        ),
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
    interval_count: int
    f_candidate_count: int
    f_unique: bool

    fixed_mod_count: int
    fixed_mod_survivors: int
    fixed_mod_unique: bool

    adaptive_selected_primes: int
    adaptive_mod_evaluations: int
    adaptive_survivors: int
    adaptive_unique: bool
    adaptive_stalled: bool

    fixed_exact_evaluations: int
    fixed_exact_survivors: int
    fixed_exact_unique: bool

    adaptive_exact_evaluations: int
    adaptive_exact_survivors: int
    adaptive_exact_unique: bool

    true_survives_fixed: bool
    true_survives_adaptive: bool
    true_survives_exact: bool


# -----------------------------------------------------------------------------
# Utility
# -----------------------------------------------------------------------------

def ratio_string(a: int, b: int) -> str:
    if b == 0:
        return "0"
    return f"{a / b:.8f}"


def yes_no(flag: bool) -> str:
    return "YES" if flag else "NO"


# -----------------------------------------------------------------------------
# Run one test
# -----------------------------------------------------------------------------

def run_test(
    instance_index: int,
    inst: Instance,
    u: int,
    residue_tables: Dict[int, List[bool]],
) -> Tuple[TestResult, AdaptiveSieveResult | None]:

    pk = partial_k_data(inst, u)
    interval = exact_d_interval(pk)

    candidates_result = enumerate_f_candidates(
        pk,
        interval,
        ENUMERATION_LIMIT,
    )

    if not candidates_result.enumerated:
        return (
            TestResult(
                instance_index=instance_index,
                u=u,
                gap=inst.gap,
                x_true=pk.x_true,
                interval_count=interval.candidate_count,
                f_candidate_count=0,
                f_unique=False,
                fixed_mod_count=0,
                fixed_mod_survivors=0,
                fixed_mod_unique=False,
                adaptive_selected_primes=0,
                adaptive_mod_evaluations=0,
                adaptive_survivors=0,
                adaptive_unique=False,
                adaptive_stalled=True,
                fixed_exact_evaluations=0,
                fixed_exact_survivors=0,
                fixed_exact_unique=False,
                adaptive_exact_evaluations=0,
                adaptive_exact_survivors=0,
                adaptive_exact_unique=False,
                true_survives_fixed=False,
                true_survives_adaptive=False,
                true_survives_exact=False,
            ),
            None,
        )

    candidates = candidates_result.values

    if not candidates:
        raise AssertionError("Expected at least one f-candidate.")

    if inst.d not in candidates:
        raise AssertionError("True d missing from f-candidates.")

    fixed = fixed_modular_sieve(
        inst,
        candidates,
        MODULI,
        residue_tables,
    )

    if not fixed.true_survives:
        raise AssertionError(
            "Fixed modular sieve rejected the true candidate."
        )

    adaptive = adaptive_modular_sieve(
        inst,
        candidates,
        MODULI,
        residue_tables,
    )

    if not adaptive.true_survives:
        raise AssertionError(
            "Adaptive modular sieve rejected the true candidate."
        )

    fixed_exact = exact_gap_filter(
        inst,
        fixed.survivors,
    )

    adaptive_exact = exact_gap_filter(
        inst,
        adaptive.survivors,
    )

    if not fixed_exact.true_reconstructed:
        raise AssertionError(
            "Exact gap test failed to reconstruct true factors after fixed sieve."
        )

    if not adaptive_exact.true_reconstructed:
        raise AssertionError(
            "Exact gap test failed to reconstruct true factors after adaptive sieve."
        )

    result = TestResult(
        instance_index=instance_index,
        u=u,
        gap=inst.gap,
        x_true=pk.x_true,
        interval_count=interval.candidate_count,
        f_candidate_count=len(candidates),
        f_unique=(
            len(candidates) == 1
            and candidates[0] == inst.d
        ),

        fixed_mod_count=len(MODULI),
        fixed_mod_survivors=len(fixed.survivors),
        fixed_mod_unique=(
            len(fixed.survivors) == 1
            and fixed.survivors[0] == inst.d
        ),

        adaptive_selected_primes=len(adaptive.steps),
        adaptive_mod_evaluations=adaptive.modular_evaluations,
        adaptive_survivors=len(adaptive.survivors),
        adaptive_unique=(
            len(adaptive.survivors) == 1
            and adaptive.survivors[0] == inst.d
        ),
        adaptive_stalled=adaptive.stalled,

        fixed_exact_evaluations=fixed_exact.exact_square_evaluations,
        fixed_exact_survivors=len(fixed_exact.reconstructed),
        fixed_exact_unique=(
            len(fixed_exact.reconstructed) == 1
            and fixed_exact.reconstructed[0].d == inst.d
        ),

        adaptive_exact_evaluations=adaptive_exact.exact_square_evaluations,
        adaptive_exact_survivors=len(adaptive_exact.reconstructed),
        adaptive_exact_unique=(
            len(adaptive_exact.reconstructed) == 1
            and adaptive_exact.reconstructed[0].d == inst.d
        ),

        true_survives_fixed=fixed.true_survives,
        true_survives_adaptive=adaptive.true_survives,
        true_survives_exact=(
            fixed_exact.true_reconstructed
            and adaptive_exact.true_reconstructed
        ),
    )

    return result, adaptive


# -----------------------------------------------------------------------------
# Detailed adaptive trace
# -----------------------------------------------------------------------------

def print_adaptive_trace(
    inst: Instance,
    result: TestResult,
    adaptive: AdaptiveSieveResult,
) -> None:

    print()
    print("=" * 120)
    print(
        f"DETAIL: INSTANCE {result.instance_index}, "
        f"UNKNOWN K BITS = {result.u}"
    )
    print("=" * 120)

    print(f"p={inst.p}  q={inst.q}  gap={inst.gap}")
    print(f"N={inst.N}")
    print(f"d={inst.d}")
    print(f"K={inst.K}")

    print()
    print("CANDIDATE SET")
    print(f"  f interval candidates = {result.interval_count}")
    print(f"  exact f candidates    = {result.f_candidate_count}")
    print(f"  f unique              = {result.f_unique}")

    print()
    print("FIXED MODULAR SIEVE")
    print(f"  primes configured     = {result.fixed_mod_count}")
    print(f"  survivors             = {result.fixed_mod_survivors}")
    print(f"  unique                = {result.fixed_mod_unique}")

    print()
    print("ADAPTIVE MODULAR SIEVE")
    print(f"  selected primes       = {result.adaptive_selected_primes}")
    print(
        f"  selection evaluations = "
        f"{result.adaptive_mod_evaluations}"
    )
    print(f"  survivors             = {result.adaptive_survivors}")
    print(f"  unique                = {result.adaptive_unique}")
    print(f"  stalled               = {result.adaptive_stalled}")

    if adaptive.steps:
        print()
        print("ADAPTIVE PRIME ORDER")

        for idx, step in enumerate(adaptive.steps, start=1):
            print(
                f"  step={idx:2d} "
                f"mod={step.modulus:2d} "
                f"{step.before:8d} -> {step.after:8d} "
                f"rejected={step.rejected:8d}"
            )

    print()
    print("EXACT GAP STAGE")
    print(
        f"  fixed survivors exact-tested = "
        f"{result.fixed_exact_evaluations}"
    )
    print(
        f"  fixed reconstructed          = "
        f"{result.fixed_exact_survivors}"
    )
    print(
        f"  adaptive survivors exact-tested = "
        f"{result.adaptive_exact_evaluations}"
    )
    print(
        f"  adaptive reconstructed          = "
        f"{result.adaptive_exact_survivors}"
    )

    print()
    print("INVARIANTS")
    print(f"  true survives fixed adaptive = {result.true_survives_adaptive}")
    print(f"  true survives exact          = {result.true_survives_exact}")


# -----------------------------------------------------------------------------
# Main experiment
# -----------------------------------------------------------------------------

def run_experiment() -> None:

    print("=" * 120)
    print("EXPERIMENT 425 START")
    print("=" * 120)
    print()
    print("ADAPTIVE MODULAR GAP-SQUARE SIEVE / EXACT CANDIDATE-ISOLATION AUDIT")
    print()
    print("QUESTION")
    print("  Can adaptive prime selection isolate the true candidate")
    print("  using fewer selected modular filters than a fixed sieve?")
    print()
    print("RULES")
    print("  exact integer arithmetic only")
    print("  no resultants")
    print("  no Groebner basis")
    print("  no symbolic factorization")
    print("  modular conditions are necessary only")
    print()

    print("=" * 120)
    print("CONFIGURATION")
    print("=" * 120)

    print(f"  instances            = 8")
    print(f"  K bits               = {UNKNOWN_BITS}")
    print(f"  modular primes       = {MODULI}")
    print(f"  enumeration max      = {ENUMERATION_LIMIT}")
    print(f"  adaptive stop-one    = {STOP_ON_ONE}")

    # Fresh/generated semiprime instances.
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

    residue_tables = {
        modulus: quadratic_residue_table(modulus)
        for modulus in MODULI
    }

    all_results: List[TestResult] = []
    detailed: List[
        Tuple[TestResult, Instance, AdaptiveSieveResult]
    ] = []

    # -------------------------------------------------------------------------
    # Compact result header
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("COMPACT ADAPTIVE-SIEVE SUMMARY")
    print("=" * 120)

    print(
        " i   u       gap        x    f-cand "
        "fixed-mod fixed-surv "
        "adapt-p adapt-surv adapt-uniq "
        "adapt-evals"
    )
    print("-" * 120)

    total_tests = 0
    enumerated_tests = 0
    skipped_tests = 0

    fixed_unique_count = 0
    adaptive_unique_count = 0
    fixed_exact_unique_count = 0
    adaptive_exact_unique_count = 0

    adaptive_true_survival_count = 0
    fixed_true_survival_count = 0

    total_f_candidates = 0
    total_fixed_mod_evals = 0
    total_adaptive_mod_evals = 0

    total_fixed_exact_evals = 0
    total_adaptive_exact_evals = 0

    total_fixed_mod_survivors = 0
    total_adaptive_survivors = 0

    for index, (p, q) in enumerate(instances, start=1):

        inst = make_instance(p, q)

        for u in UNKNOWN_BITS:

            total_tests += 1

            result, adaptive = run_test(
                index,
                inst,
                u,
                residue_tables,
            )

            all_results.append(result)

            if result.f_candidate_count == 0:
                skipped_tests += 1

                print(
                    f"{index:2d} {u:3d} "
                    f"{inst.gap:10d} "
                    f"{result.x_true:10d} "
                    f"{result.interval_count:10d} "
                    f"SKIPPED"
                )
                continue

            enumerated_tests += 1

            fixed_unique_count += result.fixed_mod_unique
            adaptive_unique_count += result.adaptive_unique
            fixed_exact_unique_count += result.fixed_exact_unique
            adaptive_exact_unique_count += result.adaptive_exact_unique

            fixed_true_survival_count += result.true_survives_fixed
            adaptive_true_survival_count += result.true_survives_adaptive

            total_f_candidates += result.f_candidate_count

            total_fixed_mod_evals += (
                result.f_candidate_count *
                result.fixed_mod_count
            )

            total_adaptive_mod_evals += (
                result.adaptive_mod_evaluations
            )

            total_fixed_exact_evals += (
                result.fixed_exact_evaluations
            )

            total_adaptive_exact_evals += (
                result.adaptive_exact_evaluations
            )

            total_fixed_mod_survivors += (
                result.fixed_mod_survivors
            )

            total_adaptive_survivors += (
                result.adaptive_survivors
            )

            print(
                f"{index:2d} {u:3d} "
                f"{result.gap:10d} "
                f"{result.x_true:10d} "
                f"{result.f_candidate_count:10d} "
                f"{result.fixed_mod_count:9d} "
                f"{result.fixed_mod_survivors:9d} "
                f"{result.adaptive_selected_primes:7d} "
                f"{result.adaptive_survivors:9d} "
                f"{yes_no(result.adaptive_unique):10s} "
                f"{result.adaptive_mod_evaluations:11d}"
            )

            # Detail only useful cases.
            if (
                SHOW_DETAILED_CASES
                and adaptive is not None
                and result.f_candidate_count >= DETAIL_MIN_INITIAL_CANDIDATES
                and (
                    not result.f_unique
                    or result.adaptive_selected_primes < len(MODULI)
                    or not result.adaptive_unique
                    or result.adaptive_survivors > 1
                )
            ):
                detailed.append(
                    (
                        result,
                        inst,
                        adaptive,
                    )
                )

    # -------------------------------------------------------------------------
    # Global summary
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("GLOBAL ADAPTIVE-SIEVE SUMMARY")
    print("=" * 120)

    print(f"  total tests                         = {total_tests}")
    print(f"  enumerated tests                    = {enumerated_tests}")
    print(f"  skipped tests                       = {skipped_tests}")
    print(
        f"  unique after fixed modular sieve   = "
        f"{fixed_unique_count}/{enumerated_tests}"
    )
    print(
        f"  unique after adaptive modular sieve = "
        f"{adaptive_unique_count}/{enumerated_tests}"
    )
    print(
        f"  exact unique after fixed sieve     = "
        f"{fixed_exact_unique_count}/{enumerated_tests}"
    )
    print(
        f"  exact unique after adaptive sieve  = "
        f"{adaptive_exact_unique_count}/{enumerated_tests}"
    )

    print(
        f"  true survives fixed sieve          = "
        f"{fixed_true_survival_count}/{enumerated_tests}"
    )
    print(
        f"  true survives adaptive sieve       = "
        f"{adaptive_true_survival_count}/{enumerated_tests}"
    )

    print()
    print("  total f-candidates                 = "
          f"{total_f_candidates}")
    print(
        "  total fixed modular evaluations   = "
        f"{total_fixed_mod_evals}"
    )
    print(
        "  total adaptive modular evaluations = "
        f"{total_adaptive_mod_evals}"
    )
    print(
        "  total fixed exact tests           = "
        f"{total_fixed_exact_evals}"
    )
    print(
        "  total adaptive exact tests        = "
        f"{total_adaptive_exact_evals}"
    )

    print()
    print("  adaptive / fixed modular evaluation ratio = "
          f"{ratio_string(total_adaptive_mod_evals, total_fixed_mod_evals)}")

    print(
        "  adaptive / fixed exact-test ratio        = "
        f"{ratio_string(total_adaptive_exact_evals, total_fixed_exact_evals)}"
    )

    # -------------------------------------------------------------------------
    # Aggregate by unknown bits
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("AGGREGATE BY UNKNOWN K BITS")
    print("=" * 120)

    print(
        " u   tests  enum  f-uniq  fixed-uniq "
        "adapt-uniq  avg-f  avg-fixed  avg-adapt  avg-primes"
    )
    print("-" * 120)

    for u in UNKNOWN_BITS:

        rows = [
            r for r in all_results
            if r.u == u and r.f_candidate_count > 0
        ]

        tests = len(rows)

        if tests == 0:
            continue

        f_unique = sum(r.f_unique for r in rows)
        fixed_unique = sum(r.fixed_mod_unique for r in rows)
        adapt_unique = sum(r.adaptive_unique for r in rows)

        avg_f = (
            sum(r.f_candidate_count for r in rows) / tests
        )

        avg_fixed = (
            sum(r.fixed_mod_survivors for r in rows) / tests
        )

        avg_adapt = (
            sum(r.adaptive_survivors for r in rows) / tests
        )

        avg_primes = (
            sum(r.adaptive_selected_primes for r in rows) / tests
        )

        print(
            f"{u:3d} "
            f"{len(rows):6d} "
            f"{tests:5d} "
            f"{f_unique:7d} "
            f"{fixed_unique:11d} "
            f"{adapt_unique:10d} "
            f"{avg_f:8.2f} "
            f"{avg_fixed:10.2f} "
            f"{avg_adapt:10.2f} "
            f"{avg_primes:10.2f}"
        )

    # -------------------------------------------------------------------------
    # Strongest adaptive reductions
    # -------------------------------------------------------------------------

    ranked = [
        r for r in all_results
        if r.f_candidate_count > 0
        and r.f_candidate_count > 1
    ]

    ranked.sort(
        key=lambda r: (
            r.f_candidate_count == 0,
            r.adaptive_survivors,
            r.adaptive_selected_primes,
            r.f_candidate_count,
        )
    )

    print()
    print("=" * 120)
    print("STRONGEST ADAPTIVE ISOLATIONS")
    print("=" * 120)

    shown = 0

    for r in ranked:

        if shown >= 20:
            break

        fixed_ratio = (
            r.fixed_mod_survivors / r.f_candidate_count
            if r.f_candidate_count
            else 0.0
        )

        adaptive_ratio = (
            r.adaptive_survivors / r.f_candidate_count
            if r.f_candidate_count
            else 0.0
        )

        selected = r.adaptive_selected_primes

        if (
            adaptive_ratio < 0.01
            or r.adaptive_unique
            or selected < len(MODULI)
        ):
            print(
                f"  instance={r.instance_index:2d} "
                f"u={r.u:2d} "
                f"f={r.f_candidate_count:8d} "
                f"fixed={r.fixed_mod_survivors:6d} "
                f"adaptive={r.adaptive_survivors:6d} "
                f"primes={selected:2d} "
                f"fixed-ratio={fixed_ratio:.8f} "
                f"adaptive-ratio={adaptive_ratio:.8f}"
            )
            shown += 1

    # -------------------------------------------------------------------------
    # Detailed adaptive traces
    # -------------------------------------------------------------------------

    if detailed:

        print()
        print("=" * 120)
        print("INTERESTING ADAPTIVE CASES")
        print("=" * 120)

        for result, inst, adaptive in detailed:
            print_adaptive_trace(
                inst,
                result,
                adaptive,
            )

    # -------------------------------------------------------------------------
    # Exact safety checks
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("EXACT INTERNAL CHECKS")
    print("=" * 120)

    all_true_survive = all(
        r.true_survives_fixed
        and r.true_survives_adaptive
        and r.true_survives_exact
        for r in all_results
        if r.f_candidate_count > 0
    )

    exact_unique_all = all(
        r.adaptive_exact_unique
        for r in all_results
        if r.f_candidate_count > 0
    )

    print(
        "  all enumerated cases preserve true d through modular filters = "
        f"{all_true_survive}"
    )

    print(
        "  all enumerated cases are exactly resolved after adaptive sieve = "
        f"{exact_unique_all}"
    )

    # -------------------------------------------------------------------------
    # Interpretation
    # -------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("INTERPRETATION")
    print("=" * 120)

    print(
        "  Experiment 424 used a fixed modular sequence."
    )
    print(
        "  Experiment 425 chooses the next modulus adaptively"
        " from the current candidate set."
    )
    print()

    print(
        "  The adaptive selector is purely an arithmetic filtering"
        " heuristic:"
    )
    print(
        "    it never accepts a candidate that fails a necessary"
        " quadratic-residue condition."
    )
    print(
        "    it never rejects the true generated d."
    )

    print()
    print(
        "  The key comparison is not only the final candidate count."
    )
    print(
        "  It is the amount of modular work required to reach"
        " candidate isolation."
    )

    print()
    print(
        "  IMPORTANT:"
    )
    print(
        "  The adaptive-selection cost itself is measured."
    )
    print(
        "  Every unused prime is tested against the current"
        " candidate set when choosing the next prime."
    )
    print(
        "  Therefore a lower number of selected primes does not"
        " automatically imply lower total computational cost."
    )

    print()
    print(
        "  This remains an exact candidate-filter experiment,"
        " not a factoring theorem."
    )
    print(
        "  The final acceptance condition is still the exact"
        " integer gap-square reconstruction."
    )

    print()
    print("=" * 120)
    print("EXPERIMENT 425 FINAL STATUS")
    print("=" * 120)

    print(
        "  ALL ENUMERATED EXACT CHECKS = "
        f"{all_true_survive}"
    )

    print(
        "  ADAPTIVE EXACT UNIQUE RECOVERIES = "
        f"{adaptive_exact_unique_count}/{enumerated_tests}"
    )

    print(
        "  FIXED EXACT UNIQUE RECOVERIES = "
        f"{fixed_exact_unique_count}/{enumerated_tests}"
    )

    print(
        "  TOTAL F-CANDIDATES = "
        f"{total_f_candidates}"
    )

    print(
        "  FIXED MODULAR EVALUATIONS = "
        f"{total_fixed_mod_evals}"
    )

    print(
        "  ADAPTIVE MODULAR EVALUATIONS = "
        f"{total_adaptive_mod_evals}"
    )

    print(
        "  ADAPTIVE/FIXED MODULAR RATIO = "
        f"{ratio_string(total_adaptive_mod_evals, total_fixed_mod_evals)}"
    )

    print("=" * 120)
    print("EXPERIMENT 425 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    run_experiment()
