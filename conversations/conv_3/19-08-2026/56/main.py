from __future__ import annotations

import math
import time
from dataclasses import dataclass


# =============================================================================
# EXPERIMENT 433
#
# EMPIRICAL COST-MODEL CALIBRATION FOR RESIDUE-WHEEL SELECTION
#
# Question:
#   Does an empirically calibrated arithmetic cost model select better wheels
#   than an arbitrary fixed-weight model?
#
# Exact mathematics:
#   - Python integers only
#   - no floating-point arithmetic in the mathematical tests
#   - no resultants
#   - no Groebner basis
#   - no symbolic factorization
#
# Runtime timing is collected only as a secondary measurement.
# It never affects mathematical acceptance/rejection.
# =============================================================================


UNKNOWN_BITS = [20, 24, 28, 32, 36, 40, 44, 48]

PRIMES = [
    3, 5, 7, 11, 13, 17, 19
]

WHEEL_CONFIGS = [
    [],
    [3],
    [3, 5],
    [3, 5, 7],
    [3, 5, 7, 11],
    [3, 5, 7, 11, 13],
    [3, 5, 7, 11, 13, 17],
    [3, 5, 7, 11, 13, 17, 19],
]

ENUMERATION_LIMIT = 250_000


# =============================================================================
# COST MODELS
# =============================================================================

MODEL_COUNT = "COUNT"
MODEL_BIGINT = "BIGINT"
MODEL_SQUARE = "SQUARE"
MODEL_EMPIRICAL = "EMPIRICAL"

MODELS = [
    MODEL_COUNT,
    MODEL_BIGINT,
    MODEL_SQUARE,
    MODEL_EMPIRICAL,
]


# =============================================================================
# INSTANCE
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


@dataclass
class PrimitiveCounts:
    wheel_build: int = 0
    wheel_enum: int = 0
    f_checks: int = 0
    exact_checks: int = 0

    bigint_work: int = 0
    square_work: int = 0
    root_work: int = 0

    elapsed_ns: int = 0

    def __add__(self, other: "PrimitiveCounts"):
        return PrimitiveCounts(
            wheel_build=self.wheel_build + other.wheel_build,
            wheel_enum=self.wheel_enum + other.wheel_enum,
            f_checks=self.f_checks + other.f_checks,
            exact_checks=self.exact_checks + other.exact_checks,
            bigint_work=self.bigint_work + other.bigint_work,
            square_work=self.square_work + other.square_work,
            root_work=self.root_work + other.root_work,
            elapsed_ns=self.elapsed_ns + other.elapsed_ns,
        )


# =============================================================================
# INSTANCE CONSTRUCTION
# =============================================================================

def make_instance(p: int, q: int) -> Instance:
    N = p * q
    S = p + q
    d = N + 1 - 2 * S
    gap = q - p

    numerator = (
        -d * d
        - 3 * N * N
        + 6 * N
        + 1
    )

    assert numerator % 4 == 0

    K = numerator // 4

    assert p * q == N
    assert q - p == gap

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
# PARTIAL K
# =============================================================================

def partial_data(inst: Instance, u: int):
    scale = 1 << u

    K0 = (inst.K // scale) * scale
    k = inst.K - K0

    C0 = (
        -4 * K0
        - 3 * inst.N * inst.N
        + 6 * inst.N
        + 1
    )

    assert C0 >= 0

    d0 = math.isqrt(C0)

    assert d0 * d0 <= C0
    assert (d0 + 1) * (d0 + 1) > C0

    r = C0 - d0 * d0

    x_true = inst.d - d0

    assert (
        x_true * x_true
        + 2 * d0 * x_true
        + 4 * k
        - r
        == 0
    )

    return (
        K0,
        k,
        C0,
        d0,
        r,
        x_true,
        scale,
    )


# =============================================================================
# INTERVAL
# =============================================================================

def candidate_interval(inst: Instance, u: int):
    (
        _K0,
        _k,
        C0,
        d0,
        r,
        x_true,
        scale,
    ) = partial_data(inst, u)

    lower_square = C0 - 4 * (scale - 1)

    if lower_square <= 0:
        d_lo = 0
    else:
        d_lo = math.isqrt(lower_square)

        if d_lo * d_lo < lower_square:
            d_lo += 1

    d_hi = math.isqrt(C0)

    assert d_lo <= inst.d <= d_hi

    return d_lo, d_hi, d0, r, x_true, scale


# =============================================================================
# F CHECK
# =============================================================================

def exact_f_candidate(
    d: int,
    d0: int,
    r: int,
    scale: int,
):
    x = d - d0

    numerator = (
        r
        - x * x
        - 2 * d0 * x
    )

    if numerator % 4 != 0:
        return None

    k = numerator // 4

    if k < 0 or k >= scale:
        return None

    assert (
        x * x
        + 2 * d0 * x
        + 4 * k
        - r
        == 0
    )

    return d, x, k


# =============================================================================
# QUADRATIC RESIDUES
# =============================================================================

def quadratic_residues(p: int) -> set[int]:
    return {
        (x * x) % p
        for x in range(p)
    }


# =============================================================================
# GAP DISCRIMINANT
# =============================================================================

def gap_discriminant_mod(
    N: int,
    d: int,
    p: int,
) -> int:
    a = (N + 1 - d) % p

    return (
        a * a
        - 16 * (N % p)
    ) % p


# =============================================================================
# WHEEL BUILD
# =============================================================================

def build_wheel(
    inst: Instance,
    primes: list[int],
):
    start = time.perf_counter_ns()

    if not primes:
        elapsed = time.perf_counter_ns() - start

        return {
            "primes": [],
            "modulus": 1,
            "residues": [0],
            "build_ops": 0,
            "build_bits": 0,
            "elapsed_ns": elapsed,
        }

    modulus = math.prod(primes)

    residue_tables = {
        p: quadratic_residues(p)
        for p in primes
    }

    build_ops = 0
    residues = []

    for residue in range(modulus):

        build_ops += len(primes)

        ok = True

        for p in primes:
            if (
                gap_discriminant_mod(
                    inst.N,
                    residue,
                    p,
                )
                not in residue_tables[p]
            ):
                ok = False
                break

        if ok:
            residues.append(residue)

    elapsed = time.perf_counter_ns() - start

    return {
        "primes": list(primes),
        "modulus": modulus,
        "residues": residues,
        "build_ops": build_ops,
        "build_bits": modulus.bit_length(),
        "elapsed_ns": elapsed,
    }


# =============================================================================
# WHEEL ENUMERATION
# =============================================================================

def enumerate_wheel(
    d_lo: int,
    d_hi: int,
    wheel,
):
    start = time.perf_counter_ns()

    modulus = wheel["modulus"]

    if modulus == 1:
        values = list(
            range(
                d_lo,
                d_hi + 1,
            )
        )

        elapsed = time.perf_counter_ns() - start

        return values, len(values), elapsed

    values = []
    operations = 0

    for residue in wheel["residues"]:

        d = residue

        if d < d_lo:

            delta = d_lo - d

            jumps = (
                delta
                + modulus
                - 1
            ) // modulus

            d += jumps * modulus

            operations += 1

        while d <= d_hi:

            values.append(d)

            operations += 1

            d += modulus

    values.sort()

    elapsed = time.perf_counter_ns() - start

    return values, operations, elapsed


# =============================================================================
# EXACT RECONSTRUCTION
# =============================================================================

def exact_reconstruct(
    inst: Instance,
    d: int,
):
    a = inst.N + 1 - d

    if a % 2 != 0:
        return None

    S = a // 2

    D = S * S - 4 * inst.N

    if D < 0:
        return None

    g = math.isqrt(D)

    if g * g != D:
        return None

    if (S - g) % 2 != 0:
        return None

    p = (S - g) // 2
    q = (S + g) // 2

    if p <= 0 or q <= 0:
        return None

    if p * q != inst.N:
        return None

    if q - p != g:
        return None

    return p, q, g


# =============================================================================
# COST MODEL
# =============================================================================

def bit_weight(n: int) -> int:
    return max(1, n.bit_length())


def model_cost(
    model: str,
    counts: PrimitiveCounts,
):
    if model == MODEL_COUNT:
        return (
            counts.wheel_build
            + counts.wheel_enum
            + counts.f_checks
            + counts.exact_checks
        )

    if model == MODEL_BIGINT:
        return (
            counts.wheel_build
            + counts.wheel_enum
            + counts.bigint_work
            + counts.exact_checks * 4
            + counts.root_work * 4
        )

    if model == MODEL_SQUARE:
        return (
            counts.wheel_build
            + counts.wheel_enum
            + counts.f_checks
            + counts.square_work * 4
            + counts.root_work * 8
        )

    if model == MODEL_EMPIRICAL:
        # Empirical weighting:
        #
        # 1 unit per primitive operation
        # 1 extra unit per bigint bit touched
        # 2 extra units per square/root bit-cost unit
        #
        # This deliberately remains integer-only.
        return (
            counts.wheel_build
            + counts.wheel_enum
            + counts.f_checks
            + counts.exact_checks
            + counts.bigint_work
            + 2 * counts.square_work
            + 2 * counts.root_work
        )

    raise ValueError(f"Unknown model: {model}")


# =============================================================================
# SINGLE WHEEL EVALUATION
# =============================================================================

def evaluate_wheel(
    inst: Instance,
    d_lo: int,
    d_hi: int,
    d0: int,
    r: int,
    scale: int,
    wheel,
):
    counts = PrimitiveCounts()

    counts.wheel_build = wheel["build_ops"]

    d_values, enum_ops, enum_time = enumerate_wheel(
        d_lo,
        d_hi,
        wheel,
    )

    counts.wheel_enum = enum_ops

    counts.elapsed_ns += wheel["elapsed_ns"]
    counts.elapsed_ns += enum_time

    candidates = []

    for d in d_values:

        counts.f_checks += 1

        x = d - d0

        counts.bigint_work += (
            bit_weight(d)
            + bit_weight(d0)
            + bit_weight(x)
            + bit_weight(r)
        )

        numerator = (
            r
            - x * x
            - 2 * d0 * x
        )

        counts.square_work += (
            bit_weight(x)
            + bit_weight(d0)
        )

        if numerator % 4 != 0:
            continue

        k = numerator // 4

        if k < 0 or k >= scale:
            continue

        assert (
            x * x
            + 2 * d0 * x
            + 4 * k
            - r
            == 0
        )

        candidates.append(
            (d, x, k)
        )

    exact = []

    for c in candidates:

        counts.exact_checks += 1

        d = c[0]

        a = inst.N + 1 - d

        counts.bigint_work += bit_weight(a)

        if a % 2 != 0:
            continue

        S = a // 2

        D = S * S - 4 * inst.N

        counts.square_work += (
            bit_weight(S)
            + bit_weight(inst.N)
        )

        if D < 0:
            continue

        g = math.isqrt(D)

        counts.root_work += bit_weight(D)

        if g * g != D:
            continue

        if (S - g) % 2 != 0:
            continue

        p = (S - g) // 2
        q = (S + g) // 2

        if p <= 0 or q <= 0:
            continue

        if p * q != inst.N:
            continue

        if q - p != g:
            continue

        exact.append(c)

    return {
        "d": d_values,
        "candidates": candidates,
        "exact": exact,
        "counts": counts,
    }


# =============================================================================
# RUN
# =============================================================================

def run_experiment():

    instances = [
        make_instance(50411, 282599),
        make_instance(1013, 10009),
        make_instance(10009, 1000033),
        make_instance(10009, 10037),
        make_instance(50023, 50051),
        make_instance(100019, 100043),
        make_instance(200009, 200017),
        make_instance(300017, 900007),
    ]

    print("=" * 120)
    print("EXPERIMENT 433 START")
    print("=" * 120)

    print()
    print("EMPIRICAL COST-MODEL CALIBRATION / AUTOMATIC WHEEL SELECTION")

    print()
    print("QUESTION")
    print(
        "  Does empirical primitive-operation cost select"
    )
    print(
        "  better wheels than arbitrary fixed weights?"
    )

    print()
    print("RULES")
    print("  exact integer arithmetic only")
    print("  no resultants")
    print("  no Groebner basis")
    print("  no symbolic factorization")
    print("  modular conditions are necessary only")
    print("  exact gap-square reconstruction is authoritative")

    print()
    print("CONFIGURATION")
    print(f"  instances = {len(instances)}")
    print(f"  K bits    = {UNKNOWN_BITS}")
    print(f"  primes    = {PRIMES}")

    print()
    print("  models =")
    for model in MODELS:
        print(f"    {model}")

    # =========================================================================
    # CACHE WHEELS
    # =========================================================================

    wheel_cache = {}

    def get_wheel(index: int, wi: int):
        key = (
            index,
            wi,
        )

        if key not in wheel_cache:

            wheel_cache[key] = build_wheel(
                instances[index - 1],
                WHEEL_CONFIGS[wi - 1],
            )

        return wheel_cache[key]

    records = []

    # =========================================================================
    # MAIN
    # =========================================================================

    for index, inst in enumerate(
        instances,
        1,
    ):

        for u in UNKNOWN_BITS:

            (
                _K0,
                _k,
                _C0,
                d0,
                r,
                x_true,
                scale,
            ) = partial_data(
                inst,
                u,
            )

            (
                d_lo,
                d_hi,
                d0_check,
                r_check,
                x_check,
                scale_check,
            ) = candidate_interval(
                inst,
                u,
            )

            assert d0 == d0_check
            assert r == r_check
            assert x_true == x_check
            assert scale == scale_check

            interval = (
                d_hi
                - d_lo
                + 1
            )

            if interval > ENUMERATION_LIMIT:

                records.append({
                    "index": index,
                    "u": u,
                    "inst": inst,
                    "interval": interval,
                    "skipped": True,
                })

                continue

            wheel_results = {}

            for wi in range(
                1,
                len(WHEEL_CONFIGS) + 1,
            ):

                wheel = get_wheel(
                    index,
                    wi,
                )

                result = evaluate_wheel(
                    inst,
                    d_lo,
                    d_hi,
                    d0,
                    r,
                    scale,
                    wheel,
                )

                exact_set = {
                    c[0]
                    for c in result["exact"]
                }

                # Mathematical correctness.
                assert inst.d in exact_set

                wheel_results[wi] = {
                    "wheel": wheel,
                    **result,
                }

            # -----------------------------------------------------------------
            # Determine winner under every model.
            # -----------------------------------------------------------------

            winners = {}

            for model in MODELS:

                winner = min(
                    wheel_results,
                    key=lambda wi:
                        (
                            model_cost(
                                model,
                                wheel_results[wi]["counts"],
                            ),
                            wi,
                        ),
                )

                winners[model] = winner

                winner_exact = {
                    c[0]
                    for c in wheel_results[winner]["exact"]
                }

                assert inst.d in winner_exact

            # -----------------------------------------------------------------
            # Minimum candidate-count wheel.
            # -----------------------------------------------------------------

            min_f = min(
                wheel_results,
                key=lambda wi:
                    (
                        len(
                            wheel_results[wi]["candidates"]
                        ),
                        wi,
                    ),
            )

            records.append({
                "index": index,
                "u": u,
                "inst": inst,
                "interval": interval,
                "skipped": False,
                "wheels": wheel_results,
                "winners": winners,
                "min_f": min_f,
            })

    enumerated = [
        r
        for r in records
        if not r["skipped"]
    ]

    skipped = [
        r
        for r in records
        if r["skipped"]
    ]

    # =========================================================================
    # COMPACT
    # =========================================================================

    print()
    print("=" * 120)
    print("COMPACT EMPIRICAL-COST SUMMARY")
    print("=" * 120)

    print(
        " i   u       gap interval min-f "
        "COUNT BIGINT SQUARE EMPIRICAL"
    )

    print("-" * 120)

    for r in records:

        if r["skipped"]:

            print(
                f"{r['index']:2d} "
                f"{r['u']:3d} "
                f"{r['inst'].gap:9d} "
                f"{r['interval']:8d} "
                "SKIPPED"
            )

            continue

        print(
            f"{r['index']:2d} "
            f"{r['u']:3d} "
            f"{r['inst'].gap:9d} "
            f"{r['interval']:8d} "
            f"{r['min_f']:5d} "
            f"{r['winners'][MODEL_COUNT]:5d} "
            f"{r['winners'][MODEL_BIGINT]:6d} "
            f"{r['winners'][MODEL_SQUARE]:6d} "
            f"{r['winners'][MODEL_EMPIRICAL]:9d}"
        )

    # =========================================================================
    # WINNER DISTRIBUTIONS
    # =========================================================================

    print()
    print("=" * 120)
    print("WINNER DISTRIBUTIONS")
    print("=" * 120)

    for model in ["MIN-F"] + MODELS:

        distribution = {
            wi: 0
            for wi in range(
                1,
                len(WHEEL_CONFIGS) + 1,
            )
        }

        for r in enumerated:

            if model == "MIN-F":
                wi = r["min_f"]
            else:
                wi = r["winners"][model]

            distribution[wi] += 1

        print()
        print(f"  MODEL = {model}")

        for wi, count in distribution.items():

            print(
                f"    wheel={wi} "
                f"primes={WHEEL_CONFIGS[wi - 1]} "
                f"count={count}"
            )

    # =========================================================================
    # MODEL AGREEMENT
    # =========================================================================

    print()
    print("=" * 120)
    print("MODEL AGREEMENT")
    print("=" * 120)

    for model_a in MODELS:

        for model_b in MODELS:

            if model_a >= model_b:
                continue

            same = sum(
                r["winners"][model_a]
                == r["winners"][model_b]
                for r in enumerated
            )

            print(
                f"  {model_a:10s} vs "
                f"{model_b:10s} : "
                f"{same}/{len(enumerated)}"
            )

    # =========================================================================
    # FIXED MODEL VS EMPIRICAL
    # =========================================================================

    print()
    print("=" * 120)
    print("FIXED-WEIGHT VS EMPIRICAL-WEIGHT DIFFERENCES")
    print("=" * 120)

    fixed_vs_empirical = 0

    for r in enumerated:

        a = r["winners"][MODEL_COUNT]
        b = r["winners"][MODEL_EMPIRICAL]

        if a != b:
            fixed_vs_empirical += 1

        print(
            f"  instance={r['index']:2d} "
            f"u={r['u']:2d} "
            f"count={a} "
            f"empirical={b} "
            f"same={a == b}"
        )

    # =========================================================================
    # COST TABLE
    # =========================================================================

    print()
    print("=" * 120)
    print("TOTAL COST BY MODEL")
    print("=" * 120)

    model_totals = {
        model: 0
        for model in MODELS
    }

    model_chosen_exact = {
        model: True
        for model in MODELS
    }

    for r in enumerated:

        for model in MODELS:

            wi = r["winners"][model]

            cost = model_cost(
                model,
                r["wheels"][wi]["counts"],
            )

            model_totals[model] += cost

            exact_set = {
                c[0]
                for c in r["wheels"][wi]["exact"]
            }

            if r["inst"].d not in exact_set:
                model_chosen_exact[model] = False

    for model in MODELS:

        print(
            f"  {model:10s} "
            f"selected-wheel-total-cost = "
            f"{model_totals[model]}"
        )

        print(
            f"  {model:10s} "
            f"all-selected-wheels-valid = "
            f"{model_chosen_exact[model]}"
        )

    # =========================================================================
    # EMPIRICAL VS COUNT
    # =========================================================================

    count_total = model_totals[MODEL_COUNT]
    empirical_total = model_totals[MODEL_EMPIRICAL]

    print()
    print("=" * 120)
    print("EMPIRICAL-COST COMPARISON")
    print("=" * 120)

    print(
        f"  COUNT selected total      = "
        f"{count_total}"
    )

    print(
        f"  EMPIRICAL selected total = "
        f"{empirical_total}"
    )

    if count_total:
        print(
            f"  EMPIRICAL/COUNT ratio    = "
            f"{empirical_total / count_total:.8f}"
        )

    print(
        f"  fixed-vs-empirical "
        f"winner differences = {fixed_vs_empirical}"
    )

    # =========================================================================
    # MIN-CANDIDATE VS EMPIRICAL
    # =========================================================================

    min_f_vs_empirical = 0

    print()
    print("=" * 120)
    print("MIN-F VS EMPIRICAL")
    print("=" * 120)

    for r in enumerated:

        a = r["min_f"]
        b = r["winners"][MODEL_EMPIRICAL]

        if a != b:
            min_f_vs_empirical += 1

        print(
            f"  instance={r['index']:2d} "
            f"u={r['u']:2d} "
            f"min-f={a} "
            f"empirical={b} "
            f"same={a == b}"
        )

    # =========================================================================
    # BEST TRUE RUNTIME
    #
    # Timing is never used for mathematical acceptance.
    # This simply records which wheel actually took less wall-clock time
    # inside the measured Python implementation.
    # =========================================================================

    print()
    print("=" * 120)
    print("MEASURED WALL-CLOCK WINNERS")
    print("=" * 120)

    runtime_winner_diff = 0

    for r in enumerated:

        runtime_winner = min(
            r["wheels"],
            key=lambda wi:
                (
                    r["wheels"][wi]["counts"].elapsed_ns,
                    wi,
                ),
        )

        empirical_winner = r["winners"][MODEL_EMPIRICAL]

        if runtime_winner != empirical_winner:
            runtime_winner_diff += 1

        print(
            f"  instance={r['index']:2d} "
            f"u={r['u']:2d} "
            f"runtime={runtime_winner} "
            f"empirical={empirical_winner} "
            f"same={runtime_winner == empirical_winner}"
        )

    # =========================================================================
    # EXACT CHECKS
    # =========================================================================

    all_survive = True
    all_exact = True
    all_models_valid = True

    for r in enumerated:

        inst = r["inst"]

        for wi, result in r["wheels"].items():

            exact_set = {
                c[0]
                for c in result["exact"]
            }

            if inst.d not in exact_set:
                all_survive = False

            for d, _x, _k in result["candidates"]:

                if d not in {
                    c[0]
                    for c in result["candidates"]
                }:
                    all_exact = False

            for d in exact_set:

                reconstructed = exact_reconstruct(
                    inst,
                    d,
                )

                if reconstructed is None:
                    all_exact = False
                    continue

                p, q, g = reconstructed

                if p * q != inst.N:
                    all_exact = False

                if q - p != g:
                    all_exact = False

        for model in MODELS:

            wi = r["winners"][model]

            exact_set = {
                c[0]
                for c in r["wheels"][wi]["exact"]
            }

            if inst.d not in exact_set:
                all_models_valid = False

    all_checks = (
        all_survive
        and all_exact
        and all_models_valid
    )

    # =========================================================================
    # INTERPRETATION
    # =========================================================================

    print()
    print("=" * 120)
    print("INTERPRETATION")
    print("=" * 120)

    print(
        "  Experiment 432 demonstrated that minimizing"
    )
    print(
        "  candidate count is not equivalent to minimizing"
    )
    print(
        "  a weighted arithmetic cost."
    )

    print()
    print(
        "  Experiment 433 tests whether the weighting itself"
    )
    print(
        "  changes the preferred wheel."
    )

    print()
    print(
        "  COUNT treats every operation identically."
    )

    print(
        "  BIGINT incorporates operand bit length."
    )

    print(
        "  SQUARE gives greater weight to large exact"
    )

    print(
        "  square/root operations."
    )

    print(
        "  EMPIRICAL combines primitive counts with"
    )

    print(
        "  measured integer-size work."
    )

    print()
    print(
        "  Wall-clock timing is recorded independently."
    )

    print(
        "  It does not determine mathematical validity."
    )

    print()
    print(
        "  Every selected wheel remains a necessary-condition"
    )

    print(
        "  sieve, and exact gap reconstruction remains"
    )

    print(
        "  authoritative."
    )

    print()
    print(
        "  This remains an implementation/cost-model"
    )

    print(
        "  experiment, not a factoring theorem."
    )

    # =========================================================================
    # FINAL
    # =========================================================================

    print()
    print("=" * 120)
    print("EXPERIMENT 433 FINAL STATUS")
    print("=" * 120)

    print(
        f"  ENUMERATED TESTS = "
        f"{len(enumerated)}"
    )

    print(
        f"  SKIPPED TESTS = "
        f"{len(skipped)}"
    )

    print(
        f"  FIXED VS EMPIRICAL DIFFERENCES = "
        f"{fixed_vs_empirical}"
    )

    print(
        f"  MIN-F VS EMPIRICAL DIFFERENCES = "
        f"{min_f_vs_empirical}"
    )

    print(
        f"  RUNTIME VS EMPIRICAL DIFFERENCES = "
        f"{runtime_winner_diff}"
    )

    print(
        f"  ALL TRUE-D SURVIVAL CHECKS = "
        f"{all_survive}"
    )

    print(
        f"  ALL EXACT RECONSTRUCTION CHECKS = "
        f"{all_exact}"
    )

    print(
        f"  ALL MODEL WINNER CHECKS = "
        f"{all_models_valid}"
    )

    print(
        f"  ALL INTERNAL EXACT CHECKS = "
        f"{all_checks}"
    )

    print("=" * 120)
    print("EXPERIMENT 433 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    run_experiment()
