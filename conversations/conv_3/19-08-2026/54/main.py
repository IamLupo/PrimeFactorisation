from __future__ import annotations

import math
from dataclasses import dataclass


# =============================================================================
# EXPERIMENT 431
# PROGRESSIVE RESIDUE-WHEEL SIZE / COST-TRADEOFF AUDIT
# =============================================================================

UNKNOWN_BITS = [20, 24, 28, 32, 36, 40, 44, 48]

PRIMES = [
    3, 5, 7, 11, 13, 17, 19,
    23, 29, 31, 37, 41, 43, 47,
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


# =============================================================================
# INSTANCE CREATION
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

    check = (
        -4 * K
        - 3 * N * N
        + 6 * N
        + 1
    )

    assert check == d * d

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
# PARTIAL-K TRANSFORMATION
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
    r = C0 - d0 * d0

    x_true = inst.d - d0

    assert (
        x_true * x_true
        + 2 * d0 * x_true
        + 4 * k
        - r
        == 0
    )

    return K0, k, C0, d0, r, x_true, scale


# =============================================================================
# EXACT D INTERVAL
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

    return d_lo, d_hi, d0, r, x_true, scale


# =============================================================================
# EXACT F CANDIDATE
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
# STANDARD ENUMERATION
# =============================================================================

def enumerate_standard(
    d_lo: int,
    d_hi: int,
    d0: int,
    r: int,
    scale: int,
):
    interval_count = max(0, d_hi - d_lo + 1)

    if interval_count > ENUMERATION_LIMIT:
        return {
            "skipped": True,
            "interval_count": interval_count,
            "candidates": [],
        }

    candidates = []

    for d in range(d_lo, d_hi + 1):
        candidate = exact_f_candidate(
            d,
            d0,
            r,
            scale,
        )

        if candidate is not None:
            candidates.append(candidate)

    return {
        "skipped": False,
        "interval_count": interval_count,
        "candidates": candidates,
    }


# =============================================================================
# QUADRATIC RESIDUES
# =============================================================================

def quadratic_residues(p: int):
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
):
    a = (N + 1 - d) % p
    return (a * a - 16 * (N % p)) % p


def gap_mod_valid(
    N: int,
    d: int,
    p: int,
    residues: set[int],
):
    return (
        gap_discriminant_mod(N, d, p)
        in residues
    )


# =============================================================================
# BUILD EXACT WHEEL
# =============================================================================

def build_wheel(
    inst: Instance,
    wheel_primes: list[int],
):
    if not wheel_primes:
        return {
            "modulus": 1,
            "residues": [0],
            "density": 1.0,
            "tables": {},
        }

    modulus = math.prod(wheel_primes)

    tables = {
        p: quadratic_residues(p)
        for p in wheel_primes
    }

    admissible = []

    for residue in range(modulus):

        ok = True

        for p in wheel_primes:
            if not gap_mod_valid(
                inst.N,
                residue,
                p,
                tables[p],
            ):
                ok = False
                break

        if ok:
            admissible.append(residue)

    return {
        "modulus": modulus,
        "residues": admissible,
        "density": len(admissible) / modulus,
        "tables": tables,
    }


# =============================================================================
# WHEEL D ENUMERATION
# =============================================================================

def enumerate_wheel_d_values(
    d_lo: int,
    d_hi: int,
    wheel_info,
):
    modulus = wheel_info["modulus"]
    residues = wheel_info["residues"]

    d_values = []

    for residue in residues:

        d = residue

        if d < d_lo:
            delta = d_lo - d
            jumps = (
                delta + modulus - 1
            ) // modulus
            d += jumps * modulus

        while d <= d_hi:
            d_values.append(d)
            d += modulus

    d_values.sort()

    return d_values


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


def exact_filter(
    inst: Instance,
    candidates,
):
    return [
        candidate
        for candidate in candidates
        if exact_reconstruct(
            inst,
            candidate[0],
        ) is not None
    ]


# =============================================================================
# MAIN
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
    print("EXPERIMENT 431 START")
    print("=" * 120)

    print()
    print("PROGRESSIVE RESIDUE-WHEEL SIZE / COST-TRADEOFF AUDIT")

    print()
    print("QUESTION")
    print("  How many wheel primes should be pushed before")
    print("  exact f-candidate computation?")

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
    print("  instances =", len(instances))
    print("  K bits    =", UNKNOWN_BITS)
    print("  primes    =", PRIMES)

    print("  wheel configurations =")
    for i, config in enumerate(WHEEL_CONFIGS, 1):
        print(f"    {i}: {config}")

    records = []

    # =========================================================================
    # COMPUTATION
    # =========================================================================

    for index, inst in enumerate(instances, 1):

        for u in UNKNOWN_BITS:

            (
                _K0,
                _k,
                _C0,
                d0,
                r,
                x_true,
                scale,
            ) = partial_data(inst, u)

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

            standard = enumerate_standard(
                d_lo,
                d_hi,
                d0,
                r,
                scale,
            )

            record = {
                "index": index,
                "u": u,
                "inst": inst,
                "d_lo": d_lo,
                "d_hi": d_hi,
                "d0": d0,
                "r": r,
                "x_true": x_true,
                "scale": scale,
                "standard": standard,
                "wheels": {},
            }

            if standard["skipped"]:
                records.append(record)
                continue

            standard_candidates = standard["candidates"]

            standard_set = {
                c[0]
                for c in standard_candidates
            }

            assert inst.d in standard_set

            for wheel_index, wheel_primes in enumerate(
                WHEEL_CONFIGS,
                1,
            ):

                wheel_info = build_wheel(
                    inst,
                    wheel_primes,
                )

                wheel_d_values = enumerate_wheel_d_values(
                    d_lo,
                    d_hi,
                    wheel_info,
                )

                wheel_d_set = set(wheel_d_values)

                # -------------------------------------------------------------
                # Independent exact wheel validation.
                #
                # Expected set is ALL d values satisfying the modular wheel,
                # not merely f-candidates.
                # -------------------------------------------------------------

                expected_d_set = set()

                for d in range(d_lo, d_hi + 1):

                    ok = True

                    for p in wheel_primes:

                        residues = wheel_info["tables"][p]

                        if not gap_mod_valid(
                            inst.N,
                            d,
                            p,
                            residues,
                        ):
                            ok = False
                            break

                    if ok:
                        expected_d_set.add(d)

                assert wheel_d_set == expected_d_set

                # True d must always survive the necessary wheel.
                assert inst.d in wheel_d_set

                # -------------------------------------------------------------
                # Now apply the exact f equation only to wheel survivors.
                # -------------------------------------------------------------

                wheel_candidates = []

                for d in wheel_d_values:

                    candidate = exact_f_candidate(
                        d,
                        d0,
                        r,
                        scale,
                    )

                    if candidate is not None:
                        wheel_candidates.append(candidate)

                wheel_candidate_set = {
                    c[0]
                    for c in wheel_candidates
                }

                # Every wheel f-candidate must also be an ordinary f-candidate.
                assert wheel_candidate_set <= standard_set

                # The wheel candidate set must equal:
                # standard f candidates ∩ wheel-admissible d values.
                expected_candidate_set = (
                    standard_set
                    & wheel_d_set
                )

                assert (
                    wheel_candidate_set
                    == expected_candidate_set
                )

                assert inst.d in wheel_candidate_set

                # -------------------------------------------------------------
                # Exact gap reconstruction.
                # -------------------------------------------------------------

                exact_candidates = exact_filter(
                    inst,
                    wheel_candidates,
                )

                exact_set = {
                    c[0]
                    for c in exact_candidates
                }

                assert inst.d in exact_set

                record["wheels"][wheel_index] = {
                    "primes": wheel_primes,
                    "wheel": wheel_info,
                    "d_values": wheel_d_values,
                    "wheel_hits": len(wheel_d_values),
                    "f_checks": len(wheel_d_values),
                    "candidates": wheel_candidates,
                    "exact": exact_candidates,
                }

            records.append(record)

    # =========================================================================
    # BASIC COUNTS
    # =========================================================================

    enumerated = [
        r for r in records
        if not r["standard"]["skipped"]
    ]

    skipped = [
        r for r in records
        if r["standard"]["skipped"]
    ]

    total_interval = sum(
        r["standard"]["interval_count"]
        for r in enumerated
    )

    total_standard_f = sum(
        len(r["standard"]["candidates"])
        for r in enumerated
    )

    # =========================================================================
    # COMPACT SUMMARY
    # =========================================================================

    print()
    print("=" * 120)
    print("COMPACT PROGRESSIVE-WHEEL SUMMARY")
    print("=" * 120)

    print(
        " i   u   gap        interval f-cand "
        + " ".join(
            f"w{i}"
            for i in range(1, len(WHEEL_CONFIGS) + 1)
        )
    )

    print("-" * 120)

    for r in records:

        if r["standard"]["skipped"]:
            print(
                f"{r['index']:2d} "
                f"{r['u']:3d} "
                f"{r['inst'].gap:10d} "
                f"{r['standard']['interval_count']:9d} "
                f"SKIPPED"
            )
            continue

        values = []

        for wheel_index in range(
            1,
            len(WHEEL_CONFIGS) + 1,
        ):
            values.append(
                len(
                    r["wheels"][
                        wheel_index
                    ]["candidates"]
                )
            )

        print(
            f"{r['index']:2d} "
            f"{r['u']:3d} "
            f"{r['inst'].gap:10d} "
            f"{r['standard']['interval_count']:9d} "
            f"{len(r['standard']['candidates']):6d} "
            + " ".join(
                f"{x:4d}"
                for x in values
            )
        )

    # =========================================================================
    # GLOBAL WHEEL COST
    # =========================================================================

    print()
    print("=" * 120)
    print("GLOBAL WHEEL COST SUMMARY")
    print("=" * 120)

    print(
        " wheel primes                                modulus "
        "density      wheel-hits   f-ratio    reduction"
    )
    print("-" * 120)

    global_metrics = {}

    for wheel_index, primes in enumerate(
        WHEEL_CONFIGS,
        1,
    ):

        total_hits = sum(
            r["wheels"][wheel_index]["wheel_hits"]
            for r in enumerated
        )

        first = enumerated[0]["wheels"][wheel_index]["wheel"]

        modulus = first["modulus"]
        density = first["density"]

        ratio = (
            total_hits
            / max(1, total_interval)
        )

        reduction = 1.0 - ratio

        global_metrics[wheel_index] = {
            "hits": total_hits,
            "modulus": modulus,
            "density": density,
            "ratio": ratio,
            "reduction": reduction,
        }

        print(
            f"{wheel_index:5d} "
            f"{str(primes):40s} "
            f"{modulus:8d} "
            f"{density:10.8f} "
            f"{total_hits:12d} "
            f"{ratio:10.8f} "
            f"{reduction:10.8f}"
        )

    # =========================================================================
    # WHEEL QUALITY
    # =========================================================================

    print()
    print("=" * 120)
    print("WHEEL QUALITY")
    print("=" * 120)

    for wheel_index, primes in enumerate(
        WHEEL_CONFIGS,
        1,
    ):

        info = enumerated[0]["wheels"][wheel_index]["wheel"]

        print(
            f"  wheel={wheel_index:2d} "
            f"primes={primes!s:32s} "
            f"modulus={info['modulus']:8d} "
            f"admissible-residues={len(info['residues']):8d} "
            f"density={info['density']:.8f}"
        )

    # =========================================================================
    # PER-TEST BEST WHEEL
    # =========================================================================

    print()
    print("=" * 120)
    print("PER-TEST BEST WHEEL")
    print("=" * 120)

    print(
        " i   u   gap        interval f-cand "
        "best-wheel best-f reduction exact"
    )

    print("-" * 120)

    best_distribution = {}

    for r in enumerated:

        interval = r["standard"]["interval_count"]
        standard_f = len(r["standard"]["candidates"])

        best_index = min(
            range(1, len(WHEEL_CONFIGS) + 1),
            key=lambda wi:
                r["wheels"][wi]["f_checks"],
        )

        best = r["wheels"][best_index]

        ratio = (
            best["f_checks"]
            / max(1, interval)
        )

        reduction = 1.0 - ratio

        best_distribution[best_index] = (
            best_distribution.get(best_index, 0)
            + 1
        )

        print(
            f"{r['index']:2d} "
            f"{r['u']:3d} "
            f"{r['inst'].gap:10d} "
            f"{interval:9d} "
            f"{standard_f:6d} "
            f"{best_index:10d} "
            f"{best['f_checks']:8d} "
            f"{reduction:9.6f} "
            f"{len(best['exact']):5d}"
        )

    # =========================================================================
    # BEST WHEEL DISTRIBUTION
    # =========================================================================

    print()
    print("=" * 120)
    print("BEST-WHEEL DISTRIBUTION")
    print("=" * 120)

    for wheel_index, primes in enumerate(
        WHEEL_CONFIGS,
        1,
    ):
        print(
            f"  wheel {wheel_index}: "
            f"{str(primes):40s} "
            f"best on "
            f"{best_distribution.get(wheel_index, 0)} tests"
        )

    # =========================================================================
    # STRONGEST REDUCTIONS
    # =========================================================================

    ranked = sorted(
        enumerated,
        key=lambda r: min(
            r["wheels"][wi]["f_checks"]
            / max(1, r["standard"]["interval_count"])
            for wi in range(
                1,
                len(WHEEL_CONFIGS) + 1,
            )
        ),
    )

    print()
    print("=" * 120)
    print("STRONGEST WHEEL REDUCTIONS")
    print("=" * 120)

    for r in ranked[:20]:

        best_index = min(
            range(1, len(WHEEL_CONFIGS) + 1),
            key=lambda wi:
                r["wheels"][wi]["f_checks"],
        )

        best = r["wheels"][best_index]

        ratio = (
            best["f_checks"]
            / max(
                1,
                r["standard"]["interval_count"],
            )
        )

        print(
            f"  instance={r['index']:2d} "
            f"u={r['u']:2d} "
            f"interval={r['standard']['interval_count']:8d} "
            f"wheel={best_index} "
            f"primes={WHEEL_CONFIGS[best_index - 1]} "
            f"wheel-f={best['f_checks']:8d} "
            f"ratio={ratio:.8f}"
        )

    # =========================================================================
    # FINAL EXACT CHECKS
    # =========================================================================

    all_true_survival = True
    all_exact_valid = True
    all_set_relations = True

    for r in enumerated:

        inst = r["inst"]

        standard_f_set = {
            c[0]
            for c in r["standard"]["candidates"]
        }

        standard_exact_set = {
            c[0]
            for c in exact_filter(
                inst,
                r["standard"]["candidates"],
            )
        }

        for wheel_index in range(
            1,
            len(WHEEL_CONFIGS) + 1,
        ):

            wheel = r["wheels"][wheel_index]

            wheel_d_set = set(
                wheel["d_values"]
            )

            wheel_f_set = {
                c[0]
                for c in wheel["candidates"]
            }

            wheel_exact_set = {
                c[0]
                for c in wheel["exact"]
            }

            if inst.d not in wheel_d_set:
                all_true_survival = False

            if inst.d not in wheel_f_set:
                all_true_survival = False

            if inst.d not in wheel_exact_set:
                all_true_survival = False

            if not wheel_f_set <= standard_f_set:
                all_set_relations = False

            if (
                wheel_f_set
                != standard_f_set & wheel_d_set
            ):
                all_set_relations = False

            # Exact reconstruction can only reduce to actual factor pairs.
            for d in wheel_exact_set:

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

            # Since the wheel is a necessary condition,
            # every exact standard reconstruction must survive it.
            if not standard_exact_set <= wheel_d_set:
                all_set_relations = False

    all_checks = (
        all_true_survival
        and all_exact_valid
        and all_set_relations
    )

    # =========================================================================
    # INTERPRETATION
    # =========================================================================

    print()
    print("=" * 120)
    print("INTERPRETATION")
    print("=" * 120)

    print(
        "  Experiment 430 showed that a small residue wheel"
    )
    print(
        "  can reject most d-values before evaluating f."
    )

    print()
    print(
        "  Experiment 431 varies the number of wheel primes"
    )
    print(
        "  to measure the trade-off between wheel density"
    )
    print(
        "  and the number of d-values reaching f."
    )

    print()
    print(
        "  For every wheel:"
    )

    print(
        "    wheel_d_set = all interval d-values satisfying"
    )
    print(
        "                   the modular necessary conditions."
    )

    print(
        "    wheel_f_set = wheel_d_set intersected with"
    )
    print(
        "                  the exact f relation."
    )

    print()
    print(
        "  The wheel cannot reject the true d because every"
    )

    print(
        "  wheel predicate is a necessary consequence of"
    )

    print(
        "  the exact gap-square identity."
    )

    print()
    print(
        "  Increasing the wheel should lower the f-evaluation"
    )

    print(
        "  population, but larger wheels also require larger"
    )

    print(
        "  residue tables and wheel moduli."
    )

    print()
    print(
        "  This remains a candidate-generation optimization."
    )

    print(
        "  It does not establish a factorization theorem."
    )

    # =========================================================================
    # FINAL STATUS
    # =========================================================================

    print()
    print("=" * 120)
    print("EXPERIMENT 431 FINAL STATUS")
    print("=" * 120)

    print(
        f"  ENUMERATED TESTS = {len(enumerated)}"
    )

    print(
        f"  SKIPPED TESTS = {len(skipped)}"
    )

    print(
        f"  TOTAL INTERVAL D-VALUES = {total_interval}"
    )

    print(
        f"  TOTAL STANDARD F-CANDIDATES = "
        f"{total_standard_f}"
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
        f"  ALL WHEEL/F SET RELATION CHECKS = "
        f"{all_set_relations}"
    )

    print(
        f"  ALL INTERNAL EXACT CHECKS = "
        f"{all_checks}"
    )

    print("=" * 120)
    print("EXPERIMENT 431 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    run_experiment()