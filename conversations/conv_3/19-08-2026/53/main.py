from __future__ import annotations

import math
from dataclasses import dataclass


# =============================================================================
# EXPERIMENT 430
# RESIDUE-WHEEL PRE-ENUMERATION / GAP-SQUARE CANDIDATE ISOLATION
# =============================================================================

UNKNOWN_BITS = [20, 24, 28, 32, 36, 40, 44, 48]

PRIMES = [
    3, 5, 7, 11, 13, 17, 19,
    23, 29, 31, 37, 41, 43, 47,
]

WHEEL_PRIMES = [3, 5, 7, 11]

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
# INSTANCE GENERATION
# =============================================================================

def make_instance(p: int, q: int) -> Instance:
    if p <= 1 or q <= 1:
        raise ValueError("p and q must be > 1")

    if p > q:
        p, q = q, p

    N = p * q
    S = p + q
    d = N + 1 - 2 * S
    gap = q - p

    # Exact identity:
    #
    #   d^2 = -4K - 3N^2 + 6N + 1
    #
    # hence:
    #
    #   K = (-d^2 - 3N^2 + 6N + 1) / 4

    numerator = (
        -d * d
        - 3 * N * N
        + 6 * N
        + 1
    )

    if numerator % 4 != 0:
        raise ValueError("K is not integral")

    K = numerator // 4

    assert N == p * q
    assert S == p + q
    assert d == N + 1 - 2 * S
    assert gap == q - p

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

def partial_data(
    inst: Instance,
    u: int,
):
    scale = 1 << u

    K0 = (inst.K // scale) * scale
    k = inst.K - K0

    C0 = (
        -4 * K0
        - 3 * inst.N * inst.N
        + 6 * inst.N
        + 1
    )

    if C0 < 0:
        raise ValueError("C0 < 0")

    d0 = math.isqrt(C0)
    r = C0 - d0 * d0

    x_true = inst.d - d0

    f_value = (
        x_true * x_true
        + 2 * d0 * x_true
        + 4 * k
        - r
    )

    assert f_value == 0

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
# CANDIDATE INTERVAL
# =============================================================================

def candidate_interval(
    inst: Instance,
    u: int,
):
    (
        _K0,
        _k,
        C0,
        d0,
        r,
        x_true,
        scale,
    ) = partial_data(inst, u)

    lower_sq = C0 - 4 * (scale - 1)

    if lower_sq <= 0:
        d_lo = 0
    else:
        d_lo = math.isqrt(lower_sq)

        if d_lo * d_lo < lower_sq:
            d_lo += 1

    d_hi = math.isqrt(C0)

    return (
        d_lo,
        d_hi,
        d0,
        r,
        x_true,
        scale,
    )


# =============================================================================
# EXACT F CANDIDATE
# =============================================================================

def f_candidate(
    d: int,
    d0: int,
    r: int,
    scale: int,
):
    x = d - d0

    # f(x,k) = x^2 + 2*d0*x + 4*k-r = 0
    #
    # therefore:
    #
    # 4k = r - x^2 - 2*d0*x

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

    return (d, x, k)


# =============================================================================
# ORDINARY ENUMERATION
# =============================================================================

def standard_enumeration(
    d_lo: int,
    d_hi: int,
    d0: int,
    r: int,
    scale: int,
    limit: int,
):
    interval_count = max(
        0,
        d_hi - d_lo + 1,
    )

    if interval_count > limit:
        return {
            "enumerated": False,
            "interval_count": interval_count,
            "scan_count": 0,
            "candidates": [],
        }

    candidates = []

    for d in range(d_lo, d_hi + 1):
        item = f_candidate(
            d,
            d0,
            r,
            scale,
        )

        if item is not None:
            candidates.append(item)

    return {
        "enumerated": True,
        "interval_count": interval_count,
        "scan_count": interval_count,
        "candidates": candidates,
    }


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

def gap_discriminant(
    N: int,
    d: int,
) -> int:
    a = N + 1 - d
    return a * a - 16 * N


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


def modular_gap_valid(
    N: int,
    d: int,
    p: int,
    residues: set[int],
) -> bool:
    D = gap_discriminant_mod(
        N,
        d,
        p,
    )

    return D in residues


# =============================================================================
# WHEEL
# =============================================================================

def build_wheel(
    inst: Instance,
):
    modulus = math.prod(WHEEL_PRIMES)

    tables = {
        p: quadratic_residues(p)
        for p in WHEEL_PRIMES
    }

    admissible = []

    for residue in range(modulus):
        good = True

        for p in WHEEL_PRIMES:
            if not modular_gap_valid(
                inst.N,
                residue,
                p,
                tables[p],
            ):
                good = False
                break

        if good:
            admissible.append(residue)

    return (
        modulus,
        admissible,
        tables,
    )


# =============================================================================
# WHEEL ENUMERATION
# =============================================================================

def wheel_enumeration(
    inst: Instance,
    d_lo: int,
    d_hi: int,
    d0: int,
    r: int,
    scale: int,
    limit: int,
):
    interval_count = max(
        0,
        d_hi - d_lo + 1,
    )

    if interval_count > limit:
        return {
            "enumerated": False,
            "interval_count": interval_count,
            "wheel_hits": 0,
            "f_checks": 0,
            "candidates": [],
            "modulus": 0,
            "residues": [],
        }

    (
        modulus,
        admissible,
        _tables,
    ) = build_wheel(inst)

    candidates = []
    wheel_hits = 0
    f_checks = 0

    for residue in admissible:

        d = residue

        if d < d_lo:
            distance = d_lo - d
            jumps = (
                distance + modulus - 1
            ) // modulus
            d += jumps * modulus

        while d <= d_hi:

            wheel_hits += 1

            item = f_candidate(
                d,
                d0,
                r,
                scale,
            )

            f_checks += 1

            if item is not None:
                candidates.append(item)

            d += modulus

    return {
        "enumerated": True,
        "interval_count": interval_count,
        "wheel_hits": wheel_hits,
        "f_checks": f_checks,
        "candidates": candidates,
        "modulus": modulus,
        "residues": admissible,
    }


# =============================================================================
# FIXED MODULAR FILTER
# =============================================================================

def fixed_filter(
    inst: Instance,
    candidates,
):
    current = list(candidates)
    evaluations = 0
    stages = []

    tables = {
        p: quadratic_residues(p)
        for p in PRIMES
    }

    for p in PRIMES:

        residues = tables[p]
        next_candidates = []

        for item in current:

            evaluations += 1

            if modular_gap_valid(
                inst.N,
                item[0],
                p,
                residues,
            ):
                next_candidates.append(item)

        current = next_candidates
        stages.append((p, len(current)))

        if len(current) <= 1:
            break

    return (
        current,
        evaluations,
        stages,
    )


# =============================================================================
# EXACT FACTOR RECONSTRUCTION
# =============================================================================

def reconstruct(
    inst: Instance,
    d: int,
):
    A = inst.N + 1 - d

    if A % 2 != 0:
        return None

    S = A // 2

    if S * S < 4 * inst.N:
        return None

    D = S * S - 4 * inst.N

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

    return (p, q, g)


def exact_filter(
    inst: Instance,
    candidates,
):
    result = []

    for item in candidates:
        if reconstruct(
            inst,
            item[0],
        ) is not None:
            result.append(item)

    return result


# =============================================================================
# TEST INSTANCES
# =============================================================================

def build_instances():
    return [
        make_instance(50411, 282599),
        make_instance(1013, 10009),
        make_instance(10009, 1000033),
        make_instance(10009, 10037),
        make_instance(50023, 50051),
        make_instance(100019, 100043),
        make_instance(200009, 200017),
        make_instance(300017, 900007),
    ]


# =============================================================================
# MAIN EXPERIMENT
# =============================================================================

def run_experiment():

    instances = build_instances()

    print("=" * 120)
    print("EXPERIMENT 430 START")
    print("=" * 120)

    print()
    print("RESIDUE-WHEEL PRE-ENUMERATION / GAP-SQUARE CANDIDATE ISOLATION")

    print()
    print("QUESTION")
    print("  Can cheap necessary modular conditions be moved")
    print("  before the exact f-candidate computation?")

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
    print(f"  instances          = {len(instances)}")
    print(f"  K bits             = {UNKNOWN_BITS}")
    print(f"  modular primes     = {PRIMES}")
    print(f"  wheel primes       = {WHEEL_PRIMES}")
    print(f"  enumeration max    = {ENUMERATION_LIMIT}")

    records = []

    for index, inst in enumerate(
        instances,
        start=1,
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
                x_true_check,
                scale_check,
            ) = candidate_interval(
                inst,
                u,
            )

            assert d0 == d0_check
            assert r == r_check
            assert x_true == x_true_check
            assert scale == scale_check

            standard = standard_enumeration(
                d_lo,
                d_hi,
                d0,
                r,
                scale,
                ENUMERATION_LIMIT,
            )

            wheel = wheel_enumeration(
                inst,
                d_lo,
                d_hi,
                d0,
                r,
                scale,
                ENUMERATION_LIMIT,
            )

            record = {
                "index": index,
                "u": u,
                "inst": inst,
                "standard": standard,
                "wheel": wheel,
                "d0": d0,
                "r": r,
                "x_true": x_true,
                "scale": scale,
                "skipped": not standard["enumerated"],
            }

            if record["skipped"]:
                records.append(record)
                continue

            standard_d = {
                item[0]
                for item in standard["candidates"]
            }

            wheel_d = {
                item[0]
                for item in wheel["candidates"]
            }

            # The true d must survive the wheel.
            assert inst.d in standard_d
            assert inst.d in wheel_d

            # Verify that the wheel condition is exact as a
            # necessary-condition pre-filter.
            expected_wheel_d = set()

            wheel_tables = {
                p: quadratic_residues(p)
                for p in WHEEL_PRIMES
            }

            for d in standard_d:

                valid = True

                for p in WHEEL_PRIMES:

                    if not modular_gap_valid(
                        inst.N,
                        d,
                        p,
                        wheel_tables[p],
                    ):
                        valid = False
                        break

                if valid:
                    expected_wheel_d.add(d)

            assert wheel_d == expected_wheel_d

            (
                standard_modular,
                standard_modular_evals,
                standard_stages,
            ) = fixed_filter(
                inst,
                standard["candidates"],
            )

            # Remaining primes after the wheel.
            remaining_primes = [
                p
                for p in PRIMES
                if p not in WHEEL_PRIMES
            ]

            remaining_tables = {
                p: quadratic_residues(p)
                for p in remaining_primes
            }

            current = list(
                wheel["candidates"]
            )

            wheel_remaining_evals = 0
            wheel_stages = []

            for p in remaining_primes:

                residues = remaining_tables[p]
                next_candidates = []

                for item in current:

                    wheel_remaining_evals += 1

                    if modular_gap_valid(
                        inst.N,
                        item[0],
                        p,
                        residues,
                    ):
                        next_candidates.append(item)

                current = next_candidates

                wheel_stages.append(
                    (p, len(current))
                )

                if len(current) <= 1:
                    break

            standard_modular_d = {
                item[0]
                for item in standard_modular
            }

            wheel_modular_d = {
                item[0]
                for item in current
            }

            assert (
                standard_modular_d
                == wheel_modular_d
            )

            standard_exact = exact_filter(
                inst,
                standard_modular,
            )

            wheel_exact = exact_filter(
                inst,
                current,
            )

            standard_exact_d = {
                item[0]
                for item in standard_exact
            }

            wheel_exact_d = {
                item[0]
                for item in wheel_exact
            }

            assert (
                standard_exact_d
                == wheel_exact_d
            )

            assert inst.d in wheel_exact_d

            record["standard_modular"] = standard_modular
            record["standard_modular_evals"] = (
                standard_modular_evals
            )
            record["standard_stages"] = (
                standard_stages
            )

            record["wheel_remaining"] = current
            record["wheel_remaining_evals"] = (
                wheel_remaining_evals
            )
            record["wheel_stages"] = (
                wheel_stages
            )

            record["standard_exact"] = standard_exact
            record["wheel_exact"] = wheel_exact

            records.append(record)

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
    print("COMPACT RESIDUE-WHEEL SUMMARY")
    print("=" * 120)

    print(
        " i   u       gap        x    interval "
        "f-cand wheel-hit wheel-f final exact"
    )

    print("-" * 120)

    for r in records:

        if r["skipped"]:
            print(
                f"{r['index']:2d} "
                f"{r['u']:3d} "
                f"{r['inst'].gap:10d} "
                f"{r['x_true']:9d} "
                f"{r['standard']['interval_count']:9d} "
                "SKIPPED"
            )
            continue

        print(
            f"{r['index']:2d} "
            f"{r['u']:3d} "
            f"{r['inst'].gap:10d} "
            f"{r['x_true']:9d} "
            f"{r['standard']['interval_count']:9d} "
            f"{len(r['standard']['candidates']):6d} "
            f"{r['wheel']['wheel_hits']:9d} "
            f"{r['wheel']['f_checks']:7d} "
            f"{len(r['wheel_remaining']):5d} "
            f"{len(r['wheel_exact']):5d}"
        )

    # =========================================================================
    # TOTALS
    # =========================================================================

    total_interval = sum(
        r["standard"]["interval_count"]
        for r in enumerated
    )

    total_f = sum(
        len(r["standard"]["candidates"])
        for r in enumerated
    )

    total_wheel_hits = sum(
        r["wheel"]["wheel_hits"]
        for r in enumerated
    )

    total_standard_f_checks = sum(
        r["standard"]["scan_count"]
        for r in enumerated
    )

    total_wheel_f_checks = sum(
        r["wheel"]["f_checks"]
        for r in enumerated
    )

    total_standard_modular = sum(
        r["standard_modular_evals"]
        for r in enumerated
    )

    total_wheel_remaining = sum(
        r["wheel_remaining_evals"]
        for r in enumerated
    )

    total_standard_exact = sum(
        len(r["standard_exact"])
        for r in enumerated
    )

    total_wheel_exact = sum(
        len(r["wheel_exact"])
        for r in enumerated
    )

    # =========================================================================
    # GLOBAL CHECKS
    # =========================================================================

    true_survives_wheel = all(
        r["inst"].d in {
            item[0]
            for item in r["wheel"]["candidates"]
        }
        for r in enumerated
    )

    modular_equal = all(
        {
            item[0]
            for item in r["standard_modular"]
        }
        == {
            item[0]
            for item in r["wheel_remaining"]
        }
        for r in enumerated
    )

    exact_equal = all(
        {
            item[0]
            for item in r["standard_exact"]
        }
        == {
            item[0]
            for item in r["wheel_exact"]
        }
        for r in enumerated
    )

    reconstruction_valid = all(
        all(
            reconstruct(
                r["inst"],
                item[0],
            ) is not None
            for item in r["wheel_exact"]
        )
        for r in enumerated
    )

    true_d_reconstructs = all(
        reconstruct(
            r["inst"],
            r["inst"].d,
        ) is not None
        for r in enumerated
    )

    all_checks = (
        true_survives_wheel
        and modular_equal
        and exact_equal
        and reconstruction_valid
        and true_d_reconstructs
    )

    # =========================================================================
    # GLOBAL SUMMARY
    # =========================================================================

    print()
    print("=" * 120)
    print("GLOBAL RESIDUE-WHEEL SUMMARY")
    print("=" * 120)

    print(f"  total tests                  = {len(records)}")
    print(f"  enumerated tests             = {len(enumerated)}")
    print(f"  skipped tests                = {len(skipped)}")
    print(f"  total interval d-values      = {total_interval}")
    print(f"  total f-candidates           = {total_f}")
    print(f"  wheel-admissible d-hits      = {total_wheel_hits}")
    print(f"  standard f checks            = {total_standard_f_checks}")
    print(f"  wheel f checks               = {total_wheel_f_checks}")

    check_ratio = (
        total_wheel_f_checks
        / max(1, total_standard_f_checks)
    )

    check_reduction = 1.0 - check_ratio

    print(
        f"  wheel/standard f-check ratio = "
        f"{check_ratio:.8f}"
    )

    print(
        f"  f-check reduction            = "
        f"{check_reduction:.8f}"
    )

    print(
        f"  standard remaining modular   = "
        f"{total_standard_modular}"
    )

    print(
        f"  wheel remaining modular      = "
        f"{total_wheel_remaining}"
    )

    print(
        f"  standard exact survivors     = "
        f"{total_standard_exact}"
    )

    print(
        f"  wheel exact survivors        = "
        f"{total_wheel_exact}"
    )

    # =========================================================================
    # AGGREGATE
    # =========================================================================

    print()
    print("=" * 120)
    print("AGGREGATE BY UNKNOWN K BITS")
    print("=" * 120)

    print(
        " u   tests avg-int avg-f avg-wheel "
        "avg-f-check avg-wheel-f reduction"
    )

    print("-" * 120)

    for u in UNKNOWN_BITS:

        group = [
            r
            for r in enumerated
            if r["u"] == u
        ]

        if not group:
            continue

        avg_interval = (
            sum(
                r["standard"]["interval_count"]
                for r in group
            )
            / len(group)
        )

        avg_f = (
            sum(
                len(r["standard"]["candidates"])
                for r in group
            )
            / len(group)
        )

        avg_wheel = (
            sum(
                r["wheel"]["wheel_hits"]
                for r in group
            )
            / len(group)
        )

        avg_standard_checks = (
            sum(
                r["standard"]["scan_count"]
                for r in group
            )
            / len(group)
        )

        avg_wheel_checks = (
            sum(
                r["wheel"]["f_checks"]
                for r in group
            )
            / len(group)
        )

        reduction = (
            1.0
            - avg_wheel_checks / avg_standard_checks
            if avg_standard_checks
            else 0.0
        )

        print(
            f"{u:3d} "
            f"{len(group):7d} "
            f"{avg_interval:9.2f} "
            f"{avg_f:8.2f} "
            f"{avg_wheel:10.2f} "
            f"{avg_standard_checks:13.2f} "
            f"{avg_wheel_checks:13.2f} "
            f"{reduction:10.6f}"
        )

    # =========================================================================
    # WHEEL QUALITY
    # =========================================================================

    modulus, admissible, _ = build_wheel(
        instances[0]
    )

    density = (
        len(admissible)
        / modulus
    )

    print()
    print("=" * 120)
    print("WHEEL QUALITY")
    print("=" * 120)

    print(f"  wheel modulus                 = {modulus}")
    print(f"  admissible residues           = {len(admissible)}")
    print(f"  residue density               = {density:.8f}")
    print(
        f"  theoretical rejection         = "
        f"{1.0 - density:.8f}"
    )

    # =========================================================================
    # STRONGEST CASES
    # =========================================================================

    ranked = sorted(
        enumerated,
        key=lambda r: (
            r["wheel"]["f_checks"]
            / max(
                1,
                r["standard"]["scan_count"],
            )
        ),
    )

    print()
    print("=" * 120)
    print("STRONGEST PRE-ENUMERATION REDUCTIONS")
    print("=" * 120)

    for r in ranked[:20]:

        ratio = (
            r["wheel"]["f_checks"]
            / max(
                1,
                r["standard"]["scan_count"],
            )
        )

        print(
            f"  instance={r['index']:2d} "
            f"u={r['u']:2d} "
            f"interval={r['standard']['interval_count']:8d} "
            f"f={len(r['standard']['candidates']):7d} "
            f"wheel-hit={r['wheel']['wheel_hits']:7d} "
            f"wheel-f={r['wheel']['f_checks']:7d} "
            f"ratio={ratio:.8f}"
        )

    # =========================================================================
    # INTERESTING DETAILS
    # =========================================================================

    print()
    print("=" * 120)
    print("INTERESTING WHEEL CASES")
    print("=" * 120)

    for r in ranked[:10]:

        inst = r["inst"]

        wheel_d = {
            item[0]
            for item in r["wheel"]["candidates"]
        }

        standard_mod_d = {
            item[0]
            for item in r["standard_modular"]
        }

        wheel_mod_d = {
            item[0]
            for item in r["wheel_remaining"]
        }

        standard_exact_d = {
            item[0]
            for item in r["standard_exact"]
        }

        wheel_exact_d = {
            item[0]
            for item in r["wheel_exact"]
        }

        all_recon_valid = all(
            reconstruct(
                inst,
                item[0],
            ) is not None
            for item in r["wheel_exact"]
        )

        print()
        print(
            f"DETAIL: INSTANCE {r['index']}, "
            f"UNKNOWN K BITS = {r['u']}"
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
            f"  d interval values      = "
            f"{r['standard']['interval_count']}"
        )

        print(
            f"  ordinary d scans       = "
            f"{r['standard']['scan_count']}"
        )

        print(
            f"  wheel hits             = "
            f"{r['wheel']['wheel_hits']}"
        )

        print(
            f"  ordinary f candidates  = "
            f"{len(r['standard']['candidates'])}"
        )

        print(
            f"  wheel f candidates     = "
            f"{len(r['wheel']['candidates'])}"
        )

        print()
        print("MODULAR")

        print(
            f"  standard survivors     = "
            f"{len(r['standard_modular'])}"
        )

        print(
            f"  wheel survivors        = "
            f"{len(r['wheel_remaining'])}"
        )

        print(
            f"  exact survivors        = "
            f"{len(r['wheel_exact'])}"
        )

        print()
        print("INVARIANTS")

        print(
            f"  true d survives wheel  = "
            f"{inst.d in wheel_d}"
        )

        print(
            f"  modular sets equal     = "
            f"{standard_mod_d == wheel_mod_d}"
        )

        print(
            f"  exact sets equal       = "
            f"{standard_exact_d == wheel_exact_d}"
        )

        print(
            f"  all exact recon valid  = "
            f"{all_recon_valid}"
        )

    # =========================================================================
    # FINAL
    # =========================================================================

    print()
    print("=" * 120)
    print("EXACT INTERNAL CHECKS")
    print("=" * 120)

    print(
        f"  true d survives wheel pre-filter = "
        f"{true_survives_wheel}"
    )

    print(
        f"  standard/wheel modular sets equal = "
        f"{modular_equal}"
    )

    print(
        f"  standard/wheel exact sets equal   = "
        f"{exact_equal}"
    )

    print(
        f"  all exact reconstructions valid   = "
        f"{reconstruction_valid}"
    )

    print(
        f"  all true d reconstruct exactly   = "
        f"{true_d_reconstructs}"
    )

    print()
    print("=" * 120)
    print("INTERPRETATION")
    print("=" * 120)

    print(
        "  Experiment 430 moves several necessary quadratic-"
        "residue conditions before exact f evaluation."
    )

    print(
        "  The wheel consists of residue classes modulo the "
        "product of the selected small primes."
    )

    print(
        "  A d outside the wheel cannot satisfy the exact "
        "gap-square condition."
    )

    print(
        "  Therefore such d values can be rejected before "
        "computing the exact k/f relation."
    )

    print(
        "  The experiment measures candidate-generation work, "
        "not the mathematical strength of factorization."
    )

    print(
        "  Exact gap-square reconstruction remains authoritative."
    )

    print(
        "  This is an implementation optimization experiment, "
        "not a factorization theorem."
    )

    print()
    print("=" * 120)
    print("EXPERIMENT 430 FINAL STATUS")
    print("=" * 120)

    print(
        f"  ENUMERATED TESTS = "
        f"{len(enumerated)}"
    )

    print(
        f"  TOTAL INTERVAL D-VALUES = "
        f"{total_interval}"
    )

    print(
        f"  TOTAL STANDARD F-CANDIDATES = "
        f"{total_f}"
    )

    print(
        f"  TOTAL WHEEL HITS = "
        f"{total_wheel_hits}"
    )

    print(
        f"  STANDARD F CHECKS = "
        f"{total_standard_f_checks}"
    )

    print(
        f"  WHEEL F CHECKS = "
        f"{total_wheel_f_checks}"
    )

    print(
        f"  WHEEL/STANDARD F-CHECK RATIO = "
        f"{check_ratio:.8f}"
    )

    print(
        f"  F-CHECK REDUCTION = "
        f"{check_reduction:.8f}"
    )

    print(
        f"  ALL INTERNAL EXACT CHECKS = "
        f"{all_checks}"
    )

    print("=" * 120)
    print("EXPERIMENT 430 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    run_experiment()