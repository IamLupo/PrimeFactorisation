from __future__ import annotations

import math
from dataclasses import dataclass


# =============================================================================
# EXPERIMENT 432
#
# COST-AWARE AUTOMATIC RESIDUE-WHEEL SELECTION
#
# Goal:
#   Automatically choose the wheel size using an explicit cost model rather
#   than minimizing only the number of surviving f-candidates.
#
# No floating point is used in the mathematical computations.
# =============================================================================


UNKNOWN_BITS = [20, 24, 28, 32, 36, 40, 44, 48]

PRIMES = [
    3, 5, 7, 11, 13, 17, 19
]

# Wheel prefixes.
#
# Keep the largest modulus bounded.  The 7-prime wheel has modulus:
#
#   3*5*7*11*13*17*19 = 4,849,845
#
# which is still practical for an exact residue table.
#
# More primes can be added later once the cost model itself is validated.
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
# DATA
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
class WheelInfo:
    primes: list[int]
    modulus: int
    residues: list[int]
    density_num: int
    density_den: int
    build_operations: int
    residue_tables: dict[int, set[int]]


# =============================================================================
# INSTANCE CREATION
# =============================================================================

def make_instance(p: int, q: int) -> Instance:
    N = p * q
    S = p + q

    # d = N + 1 - 2S
    d = N + 1 - 2 * S

    gap = q - p

    # Exact K relation used by the previous experiments:
    #
    #   d^2 = -4K - 3N^2 + 6N + 1
    #
    # Hence:
    #
    #   K = (-d^2 - 3N^2 + 6N + 1)/4
    #
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
# PARTIAL K TRANSFORMATION
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
# EXACT INTERVAL
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

    return (
        d_lo,
        d_hi,
        d0,
        r,
        x_true,
        scale,
    )


# =============================================================================
# EXACT F RELATION
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

    # Each valid candidate must make k integral.
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
# STANDARD SEARCH
# =============================================================================

def standard_candidates(
    d_lo: int,
    d_hi: int,
    d0: int,
    r: int,
    scale: int,
):
    interval = d_hi - d_lo + 1

    if interval > ENUMERATION_LIMIT:
        return None

    result = []

    for d in range(d_lo, d_hi + 1):
        candidate = exact_f_candidate(
            d,
            d0,
            r,
            scale,
        )

        if candidate is not None:
            result.append(candidate)

    return result


# =============================================================================
# QUADRATIC RESIDUES
# =============================================================================

def quadratic_residues(p: int) -> set[int]:
    return {
        (x * x) % p
        for x in range(p)
    }


# =============================================================================
# GAP DISCRIMINANT MOD P
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


def wheel_prime_valid(
    inst: Instance,
    d: int,
    p: int,
    residues: set[int],
) -> bool:
    return (
        gap_discriminant_mod(
            inst.N,
            d,
            p,
        )
        in residues
    )


# =============================================================================
# WHEEL CONSTRUCTION
#
# Cost proxy:
#   - generating quadratic residues: p operations
#   - testing each residue of the wheel against every prime
#
# This is deliberately explicit because Experiment 431 showed that survivor
# count alone is not sufficient for choosing the optimal wheel.
# =============================================================================

def build_wheel(
    inst: Instance,
    primes: list[int],
) -> WheelInfo:

    if not primes:
        return WheelInfo(
            primes=[],
            modulus=1,
            residues=[0],
            density_num=1,
            density_den=1,
            build_operations=0,
            residue_tables={},
        )

    modulus = math.prod(primes)

    residue_tables = {
        p: quadratic_residues(p)
        for p in primes
    }

    build_operations = sum(primes)

    admissible = []

    for residue in range(modulus):

        build_operations += len(primes)

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
            admissible.append(residue)

    return WheelInfo(
        primes=list(primes),
        modulus=modulus,
        residues=admissible,
        density_num=len(admissible),
        density_den=modulus,
        build_operations=build_operations,
        residue_tables=residue_tables,
    )


# =============================================================================
# WHEEL INTERVAL ENUMERATION
#
# Every produced d satisfies the wheel exactly by construction.
#
# Cost model:
#   one residue hit + one modulus jump is one wheel enumeration operation.
# =============================================================================

def enumerate_wheel_d(
    d_lo: int,
    d_hi: int,
    wheel: WheelInfo,
):
    if wheel.modulus == 1:
        values = list(range(d_lo, d_hi + 1))
        return values, len(values)

    values = []
    operations = 0

    for residue in wheel.residues:

        d = residue

        if d < d_lo:
            delta = d_lo - d
            jumps = (
                delta + wheel.modulus - 1
            ) // wheel.modulus

            d += jumps * wheel.modulus
            operations += 1

        while d <= d_hi:
            values.append(d)
            operations += 1
            d += wheel.modulus

    values.sort()

    return values, operations


# =============================================================================
# EXACT GAP RECONSTRUCTION
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
#
# This is intentionally an integer-weighted arithmetic proxy.
#
# The purpose is not to claim machine-cycle accuracy. It provides a stable
# comparison between wheel sizes using the actual mathematical operations
# performed by the experiment.
#
# Costs:
#   wheel construction        = 1
#   wheel enumeration         = 1
#   f evaluation              = 4
#   exact reconstruction      = 12
#
# Changing these weights is supported and does not affect correctness.
# =============================================================================

COST_WHEEL_BUILD = 1
COST_WHEEL_ENUM = 1
COST_F_CHECK = 4
COST_EXACT_CHECK = 12


def score_cost(
    wheel: WheelInfo,
    wheel_enum_ops: int,
    f_checks: int,
    exact_checks: int,
):
    return (
        wheel.build_operations * COST_WHEEL_BUILD
        + wheel_enum_ops * COST_WHEEL_ENUM
        + f_checks * COST_F_CHECK
        + exact_checks * COST_EXACT_CHECK
    )


# =============================================================================
# RUN EXPERIMENT
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
    print("EXPERIMENT 432 START")
    print("=" * 120)

    print()
    print("COST-AWARE AUTOMATIC RESIDUE-WHEEL SELECTION")

    print()
    print("QUESTION")
    print("  Which wheel size minimizes total exact arithmetic work,")
    print("  rather than merely minimizing the final candidate count?")

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
    print("  wheel configurations =")

    for i, config in enumerate(
        WHEEL_CONFIGS,
        1,
    ):
        print(f"    {i}: {config}")

    print()
    print("  COST MODEL")
    print(f"    wheel-build operation = {COST_WHEEL_BUILD}")
    print(f"    wheel-enumeration op  = {COST_WHEEL_ENUM}")
    print(f"    f-check               = {COST_F_CHECK}")
    print(f"    exact reconstruction   = {COST_EXACT_CHECK}")

    # Cache wheels per (instance, wheel-index).
    wheel_cache = {}

    def get_wheel(index: int, wheel_index: int):
        key = (index, wheel_index)

        if key not in wheel_cache:
            wheel_cache[key] = build_wheel(
                instances[index - 1],
                WHEEL_CONFIGS[
                    wheel_index - 1
                ],
            )

        return wheel_cache[key]

    records = []

    # =========================================================================
    # MAIN TEST LOOP
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

            interval_count = (
                d_hi - d_lo + 1
            )

            if interval_count > ENUMERATION_LIMIT:

                records.append({
                    "index": index,
                    "u": u,
                    "inst": inst,
                    "interval": interval_count,
                    "skipped": True,
                    "wheels": {},
                })

                continue

            standard = standard_candidates(
                d_lo,
                d_hi,
                d0,
                r,
                scale,
            )

            assert standard is not None

            standard_set = {
                c[0]
                for c in standard
            }

            assert inst.d in standard_set

            wheel_results = {}

            for wheel_index in range(
                1,
                len(WHEEL_CONFIGS) + 1,
            ):

                wheel = get_wheel(
                    index,
                    wheel_index,
                )

                wheel_d, enum_ops = enumerate_wheel_d(
                    d_lo,
                    d_hi,
                    wheel,
                )

                wheel_d_set = set(wheel_d)

                assert inst.d in wheel_d_set

                wheel_candidates = []

                for d in wheel_d:

                    candidate = exact_f_candidate(
                        d,
                        d0,
                        r,
                        scale,
                    )

                    if candidate is not None:
                        wheel_candidates.append(
                            candidate
                        )

                wheel_f_set = {
                    c[0]
                    for c in wheel_candidates
                }

                expected = (
                    standard_set
                    & wheel_d_set
                )

                assert wheel_f_set == expected

                exact_candidates = []

                for candidate in wheel_candidates:

                    if (
                        exact_reconstruct(
                            inst,
                            candidate[0],
                        )
                        is not None
                    ):
                        exact_candidates.append(
                            candidate
                        )

                exact_set = {
                    c[0]
                    for c in exact_candidates
                }

                assert inst.d in exact_set

                cost = score_cost(
                    wheel,
                    enum_ops,
                    len(wheel_d),
                    len(exact_candidates),
                )

                wheel_results[wheel_index] = {
                    "wheel": wheel,
                    "d": wheel_d,
                    "enum_ops": enum_ops,
                    "f_count": len(wheel_d),
                    "candidate_count": len(wheel_candidates),
                    "exact_count": len(
                        exact_candidates
                    ),
                    "cost": cost,
                    "exact": exact_candidates,
                }

            best_index = min(
                wheel_results,
                key=lambda wi:
                    (
                        wheel_results[wi]["cost"],
                        wi,
                    ),
            )

            best = wheel_results[best_index]

            # Final correctness checks.
            assert inst.d in {
                c[0]
                for c in best["exact"]
            }

            records.append({
                "index": index,
                "u": u,
                "inst": inst,
                "interval": interval_count,
                "skipped": False,
                "standard": standard,
                "standard_count": len(standard),
                "wheels": wheel_results,
                "best_index": best_index,
                "best": best,
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
    # COMPACT SUMMARY
    # =========================================================================

    print()
    print("=" * 120)
    print("COMPACT COST-AWARE SUMMARY")
    print("=" * 120)

    print(
        " i   u       gap   interval f-cand "
        "best-wheel best-f exact cost"
    )

    print("-" * 120)

    for r in records:

        if r["skipped"]:
            print(
                f"{r['index']:2d} "
                f"{r['u']:3d} "
                f"{r['inst'].gap:10d} "
                f"{r['interval']:9d} "
                "SKIPPED"
            )
            continue

        best = r["best"]

        print(
            f"{r['index']:2d} "
            f"{r['u']:3d} "
            f"{r['inst'].gap:10d} "
            f"{r['interval']:9d} "
            f"{r['standard_count']:6d} "
            f"{r['best_index']:10d} "
            f"{best['f_count']:8d} "
            f"{best['exact_count']:5d} "
            f"{best['cost']:12d}"
        )

    # =========================================================================
    # PER-WHEEL AGGREGATES
    # =========================================================================

    print()
    print("=" * 120)
    print("PER-WHEEL COST AGGREGATES")
    print("=" * 120)

    print(
        " wheel primes                                "
        "build      enum       f      exact       total"
    )

    print("-" * 120)

    aggregate = {}

    for wheel_index, primes in enumerate(
        WHEEL_CONFIGS,
        1,
    ):

        build = 0
        enum_ops = 0
        f_count = 0
        exact_count = 0
        total_cost = 0

        for r in enumerated:

            w = r["wheels"][wheel_index]

            build += (
                w["wheel"].build_operations
            )

            enum_ops += w["enum_ops"]
            f_count += w["f_count"]
            exact_count += w["exact_count"]
            total_cost += w["cost"]

        aggregate[wheel_index] = {
            "build": build,
            "enum": enum_ops,
            "f": f_count,
            "exact": exact_count,
            "cost": total_cost,
        }

        print(
            f"{wheel_index:2d} "
            f"{str(primes):42s} "
            f"{build:10d} "
            f"{enum_ops:10d} "
            f"{f_count:10d} "
            f"{exact_count:10d} "
            f"{total_cost:14d}"
        )

    # =========================================================================
    # BEST-WHEEL DISTRIBUTION
    # =========================================================================

    distribution = {
        i: 0
        for i in range(
            1,
            len(WHEEL_CONFIGS) + 1,
        )
    }

    for r in enumerated:
        distribution[r["best_index"]] += 1

    print()
    print("=" * 120)
    print("BEST-WHEEL DISTRIBUTION")
    print("=" * 120)

    for wheel_index, primes in enumerate(
        WHEEL_CONFIGS,
        1,
    ):

        print(
            f"  wheel={wheel_index:2d} "
            f"primes={primes!s:36s} "
            f"best on "
            f"{distribution[wheel_index]} tests"
        )

    # =========================================================================
    # COST SAVINGS
    # =========================================================================

    adaptive_total = 0
    no_wheel_total = 0

    for r in enumerated:

        best = r["best"]

        adaptive_total += best["cost"]

        interval = r["interval"]
        standard_count = r["standard_count"]

        no_wheel_total += (
            standard_count * COST_F_CHECK
            + standard_count * COST_EXACT_CHECK
            + interval
        )

    print()
    print("=" * 120)
    print("AUTOMATIC-SELECTION COST")
    print("=" * 120)

    print(
        f"  no-wheel cost proxy       = "
        f"{no_wheel_total}"
    )

    print(
        f"  selected-wheel cost proxy = "
        f"{adaptive_total}"
    )

    if no_wheel_total > 0:
        print(
            f"  selected/no-wheel ratio   = "
            f"{adaptive_total / no_wheel_total:.8f}"
        )

        print(
            f"  estimated reduction       = "
            f"{1.0 - adaptive_total / no_wheel_total:.8f}"
        )

    # =========================================================================
    # STRONGEST COST SAVINGS
    # =========================================================================

    print()
    print("=" * 120)
    print("STRONGEST COST-AWARE SELECTIONS")
    print("=" * 120)

    ranked = sorted(
        enumerated,
        key=lambda r:
            r["best"]["cost"],
    )

    for r in ranked[:20]:

        best = r["best"]

        print(
            f"  instance={r['index']:2d} "
            f"u={r['u']:2d} "
            f"interval={r['interval']:8d} "
            f"selected-wheel={r['best_index']} "
            f"primes={WHEEL_CONFIGS[r['best_index'] - 1]} "
            f"f={best['f_count']:8d} "
            f"exact={best['exact_count']:3d} "
            f"cost={best['cost']:12d}"
        )

    # =========================================================================
    # IMPORTANT COUNTER-COMPARISON
    #
    # Compare:
    #   A) wheel minimizing f-candidates
    #   B) wheel minimizing total measured cost
    #
    # This is the main novelty of Experiment 432.
    # =========================================================================

    print()
    print("=" * 120)
    print("MIN-F VS MIN-COST")
    print("=" * 120)

    min_f_diff = 0
    min_cost_diff = 0

    for r in enumerated:

        min_f_index = min(
            r["wheels"],
            key=lambda wi:
                (
                    r["wheels"][wi]["candidate_count"],
                    wi,
                ),
        )

        min_cost_index = r["best_index"]

        if min_f_index != min_cost_index:
            min_f_diff += 1

        # Cost winner should always be valid.
        if r["inst"].d not in {
            c[0]
            for c in r["wheels"][min_cost_index]["exact"]
        }:
            min_cost_diff += 1

        print(
            f"  instance={r['index']:2d} "
            f"u={r['u']:2d} "
            f"min-f={min_f_index} "
            f"min-cost={min_cost_index} "
            f"same={min_f_index == min_cost_index}"
        )

    # =========================================================================
    # EXACT INTERNAL CHECKS
    # =========================================================================

    all_true_survival = True
    all_exact_valid = True
    all_cost_winners_valid = True

    for r in enumerated:

        inst = r["inst"]

        for wheel_index in r["wheels"]:

            w = r["wheels"][wheel_index]

            d_set = set(w["d"])

            exact_set = {
                c[0]
                for c in w["exact"]
            }

            if inst.d not in d_set:
                all_true_survival = False

            if inst.d not in exact_set:
                all_true_survival = False

            for d in exact_set:

                reconstructed = exact_reconstruct(
                    inst,
                    d,
                )

                if reconstructed is None:
                    all_exact_valid = False
                    continue

                p, q, g = reconstructed

                if p * q != inst.N:
                    all_exact_valid = False

                if q - p != g:
                    all_exact_valid = False

        best = r["best"]

        best_exact = {
            c[0]
            for c in best["exact"]
        }

        if inst.d not in best_exact:
            all_cost_winners_valid = False

    all_checks = (
        all_true_survival
        and all_exact_valid
        and all_cost_winners_valid
    )

    # =========================================================================
    # INTERPRETATION
    # =========================================================================

    print()
    print("=" * 120)
    print("INTERPRETATION")
    print("=" * 120)

    print(
        "  Experiment 431 optimized only the number of d-values"
    )

    print(
        "  reaching the f equation."
    )

    print()
    print(
        "  Experiment 432 assigns an explicit arithmetic cost"
    )

    print(
        "  to wheel construction, wheel enumeration, f evaluation,"
    )

    print(
        "  and exact reconstruction."
    )

    print()
    print(
        "  This allows the experiment to distinguish:"
    )

    print(
        "    minimum survivor count"
    )

    print(
        "  from:"
    )

    print(
        "    minimum total arithmetic work."
    )

    print()
    print(
        "  A larger wheel is not automatically better."
    )

    print(
        "  Its lower candidate count can be outweighed by"
    )

    print(
        "  the cost of constructing and enumerating the larger"
    )

    print(
        "  residue table."
    )

    print()
    print(
        "  The selected wheel is still only a necessary-condition"
    )

    print(
        "  filter. Exact integer gap reconstruction remains the"
    )

    print(
        "  authoritative acceptance test."
    )

    print()
    print(
        "  This is an implementation/candidate-isolation experiment,"
    )

    print(
        "  not a factorization theorem."
    )

    # =========================================================================
    # FINAL STATUS
    # =========================================================================

    print()
    print("=" * 120)
    print("EXPERIMENT 432 FINAL STATUS")
    print("=" * 120)

    print(
        f"  ENUMERATED TESTS = {len(enumerated)}"
    )

    print(
        f"  SKIPPED TESTS = {len(skipped)}"
    )

    print(
        f"  MIN-F / MIN-COST DIFFERENCES = "
        f"{min_f_diff}"
    )

    print(
        f"  INVALID COST WINNERS = "
        f"{min_cost_diff}"
    )

    print(
        f"  ALL TRUE-D SURVIVAL CHECKS = "
        f"{all_true_survival}"
    )

    print(
        f"  ALL EXACT RECONSTRUCTION CHECKS = "
        f"{all_exact_valid}"
    )

    print(
        f"  ALL COST-WINNER CHECKS = "
        f"{all_cost_winners_valid}"
    )

    print(
        f"  ALL INTERNAL EXACT CHECKS = "
        f"{all_checks}"
    )

    print("=" * 120)
    print("EXPERIMENT 432 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    run_experiment()
