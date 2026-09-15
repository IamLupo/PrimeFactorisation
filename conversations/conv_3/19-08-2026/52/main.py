from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Tuple


# =============================================================================
# EXPERIMENT 429
# ROLLING DISCRIMINANT / INCREMENTAL MODULAR GAP-SQUARE SIEVE
# =============================================================================

UNKNOWN_BITS = [20, 24, 28, 32, 36, 40, 44, 48]

PRIMES = [
    3, 5, 7, 11, 13, 17, 19,
    23, 29, 31, 37, 41, 43, 47,
]

ENUMERATION_LIMIT = 250_000


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

    K0 = (
        inst.K // scale
    ) * scale

    k = inst.K - K0

    C0 = (
        -4 * K0
        - 3 * inst.N * inst.N
        + 6 * inst.N
        + 1
    )

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

    return (
        K0,
        k,
        C0,
        d0,
        r,
        x_true,
    )


# =============================================================================
# F-CANDIDATES
# =============================================================================

def candidate_list(
    inst: Instance,
    u: int,
    limit: int,
):

    (
        _,
        _,
        C0,
        d0,
        r,
        x_true,
    ) = partial_data(inst, u)

    lo_sq = (
        C0
        - 4 * ((1 << u) - 1)
    )

    if lo_sq <= 0:
        d_lo = 0
    else:
        d_lo = math.isqrt(lo_sq)

        if d_lo * d_lo < lo_sq:
            d_lo += 1

    d_hi = math.isqrt(C0)

    interval_count = max(
        0,
        d_hi - d_lo + 1,
    )

    if interval_count > limit:
        return (
            x_true,
            d0,
            r,
            [],
            False,
            interval_count,
        )

    candidates = []

    scale = 1 << u

    for d in range(
        d_lo,
        d_hi + 1,
    ):

        delta = C0 - d * d

        if delta & 3:
            continue

        k = delta // 4

        if not (
            0 <= k < scale
        ):
            continue

        x = d - d0

        if (
            x * x
            + 2 * d0 * x
            + 4 * k
            - r
            != 0
        ):
            continue

        candidates.append(
            (
                d,
                x,
                k,
            )
        )

    return (
        x_true,
        d0,
        r,
        candidates,
        True,
        interval_count,
    )


# =============================================================================
# DIRECT DISCRIMINANT
# =============================================================================

def discriminant4(
    inst: Instance,
    d: int,
) -> int:

    a = (
        inst.N
        + 1
        - d
    )

    return (
        a * a
        - 16 * inst.N
    )


# =============================================================================
# EXACT GAP RECONSTRUCTION
# =============================================================================

def exact_reconstruct(
    inst: Instance,
    d: int,
):

    a = (
        inst.N
        + 1
        - d
    )

    if a & 1:
        return None

    S = a // 2

    D = (
        S * S
        - 4 * inst.N
    )

    if D < 0:
        return None

    g = math.isqrt(D)

    if g * g != D:
        return None

    if (S - g) & 1:
        return None

    p = (S - g) // 2
    q = (S + g) // 2

    if (
        p <= 0
        or q <= 0
        or p * q != inst.N
        or q - p != g
    ):
        return None

    return (
        p,
        q,
        g,
    )


def exact_filter(
    inst: Instance,
    candidates,
):

    return [
        c
        for c in candidates
        if exact_reconstruct(
            inst,
            c[0],
        ) is not None
    ]


# =============================================================================
# QUADRATIC RESIDUE TABLE
# =============================================================================

def residue_table(p: int):

    return {
        (x * x) % p
        for x in range(p)
    }


# =============================================================================
# DIRECT BASELINE
# =============================================================================

def baseline_masks(
    inst: Instance,
    candidates,
):

    masks: Dict[int, int] = {}

    checks = 0
    square_evaluations = 0

    for p in PRIMES:

        residues = residue_table(p)

        mask = 0

        for index, (
            d,
            _,
            _,
        ) in enumerate(candidates):

            D4 = discriminant4(
                inst,
                d,
            )

            checks += 1
            square_evaluations += 1

            if (
                D4 % p
                in residues
            ):
                mask |= (
                    1 << index
                )

        masks[p] = mask

    return (
        masks,
        checks,
        square_evaluations,
    )


# =============================================================================
# ROLLING DISCRIMINANT
# =============================================================================

def rolling_masks(
    inst: Instance,
    candidates,
):

    masks: Dict[int, int] = {}

    checks = 0
    initial_square_evaluations = 0
    updates = 0

    if not candidates:
        return (
            masks,
            checks,
            initial_square_evaluations,
            updates,
        )

    d_values = [
        c[0]
        for c in candidates
    ]

    N_mod_cache = {
        p: inst.N % p
        for p in PRIMES
    }

    for p in PRIMES:

        residues = residue_table(p)

        N16 = (
            16
            * N_mod_cache[p]
        ) % p

        d0 = d_values[0]

        a = (
            inst.N
            + 1
            - d0
        ) % p

        D4 = (
            a * a
            - N16
        ) % p

        mask = 0

        checks += 1
        initial_square_evaluations += 1

        if D4 in residues:
            mask |= 1

        for index in range(
            1,
            len(d_values),
        ):

            delta = (
                d_values[index]
                - d_values[index - 1]
            )

            delta_mod = delta % p

            # Exact recurrence:
            #
            # D4(a-delta)
            #   = (a-delta)^2 - 16N
            #   = D4 - 2*a*delta + delta^2
            #

            D4 = (
                D4
                - 2 * a * delta_mod
                + delta_mod * delta_mod
            ) % p

            a = (
                a
                - delta_mod
            ) % p

            updates += 1
            checks += 1

            if D4 in residues:
                mask |= (
                    1 << index
                )

        masks[p] = mask

    return (
        masks,
        checks,
        initial_square_evaluations,
        updates,
    )


# =============================================================================
# BITSET GREEDY
# =============================================================================

def greedy_bitset(
    masks: Dict[int, int],
    candidate_count: int,
):

    if candidate_count == 0:
        return (
            0,
            [],
            0,
        )

    current = (
        (1 << candidate_count)
        - 1
    )

    unused = set(PRIMES)

    selected = []
    intersections = 0

    while (
        unused
        and current.bit_count() > 1
    ):

        best_prime = None
        best_count = current.bit_count()

        for p in sorted(unused):

            trial = (
                current
                & masks[p]
            )

            count = trial.bit_count()

            if count < best_count:
                best_prime = p
                best_count = count

        if best_prime is None:
            break

        current &= masks[
            best_prime
        ]

        selected.append(
            best_prime
        )

        intersections += 1

        unused.remove(
            best_prime
        )

    return (
        current,
        selected,
        intersections,
    )


# =============================================================================
# MASK -> CANDIDATES
# =============================================================================

def mask_indices(
    mask: int,
    candidates,
):

    result = []

    while mask:

        bit = (
            mask
            & -mask
        )

        index = (
            bit.bit_length()
            - 1
        )

        result.append(
            candidates[index]
        )

        mask ^= bit

    return result


# =============================================================================
# MAIN
# =============================================================================

def run():

    instances = [
        make_instance(
            50411,
            282599,
        ),
        make_instance(
            1013,
            10009,
        ),
        make_instance(
            10009,
            1000033,
        ),
        make_instance(
            10009,
            10037,
        ),
        make_instance(
            50023,
            50051,
        ),
        make_instance(
            100019,
            100043,
        ),
        make_instance(
            200009,
            200017,
        ),
        make_instance(
            300017,
            900007,
        ),
    ]

    results = []

    print("=" * 120)
    print("EXPERIMENT 429 START")
    print("=" * 120)

    print()
    print(
        "ROLLING DISCRIMINANT / "
        "INCREMENTAL MODULAR GAP-SQUARE SIEVE"
    )

    print()
    print("QUESTION")
    print(
        "  Can the modular gap sieve avoid recomputing"
    )
    print(
        "  (N+1-d)^2-16N from scratch for every candidate?"
    )

    print()
    print("RULES")
    print("  exact integer arithmetic only")
    print("  no resultants")
    print("  no Groebner basis")
    print("  no symbolic factorization")
    print("  modular conditions are necessary only")
    print(
        "  exact gap-square reconstruction remains authoritative"
    )

    print()
    print("CONFIGURATION")
    print(
        f"  instances       = {len(instances)}"
    )
    print(
        f"  K bits          = {UNKNOWN_BITS}"
    )
    print(
        f"  modular primes  = {PRIMES}"
    )
    print(
        f"  enumeration max = {ENUMERATION_LIMIT}"
    )
    print(
        "  rolling update  = "
        "D4' = D4 - 2*a*delta + delta^2 (mod p)"
    )

    for index, inst in enumerate(
        instances,
        start=1,
    ):

        for u in UNKNOWN_BITS:

            (
                x_true,
                d0,
                r,
                candidates,
                enumerated,
                interval_count,
            ) = candidate_list(
                inst,
                u,
                ENUMERATION_LIMIT,
            )

            if not enumerated:

                results.append(
                    {
                        "index": index,
                        "u": u,
                        "inst": inst,
                        "x_true": x_true,
                        "count": interval_count,
                        "skipped": True,
                    }
                )

                continue

            (
                baseline,
                baseline_checks,
                baseline_squares,
            ) = baseline_masks(
                inst,
                candidates,
            )

            (
                rolling,
                rolling_checks,
                rolling_initial,
                rolling_updates,
            ) = rolling_masks(
                inst,
                candidates,
            )

            # Fundamental exact implementation check.
            assert baseline == rolling

            (
                baseline_final,
                baseline_primes,
                baseline_intersections,
            ) = greedy_bitset(
                baseline,
                len(candidates),
            )

            (
                rolling_final,
                rolling_primes,
                rolling_intersections,
            ) = greedy_bitset(
                rolling,
                len(candidates),
            )

            assert (
                baseline_final
                == rolling_final
            )

            baseline_survivors = (
                mask_indices(
                    baseline_final,
                    candidates,
                )
            )

            rolling_survivors = (
                mask_indices(
                    rolling_final,
                    candidates,
                )
            )

            baseline_exact = exact_filter(
                inst,
                baseline_survivors,
            )

            rolling_exact = exact_filter(
                inst,
                rolling_survivors,
            )

            assert (
                {
                    c[0]
                    for c in baseline_survivors
                }
                ==
                {
                    c[0]
                    for c in rolling_survivors
                }
            )

            assert (
                {
                    c[0]
                    for c in baseline_exact
                }
                ==
                {
                    c[0]
                    for c in rolling_exact
                }
            )

            assert any(
                c[0] == inst.d
                for c in rolling_survivors
            )

            results.append(
                {
                    "index": index,
                    "u": u,
                    "inst": inst,
                    "x_true": x_true,
                    "count": len(candidates),
                    "candidates": candidates,
                    "baseline_checks": baseline_checks,
                    "rolling_checks": rolling_checks,
                    "baseline_squares": baseline_squares,
                    "rolling_initial": rolling_initial,
                    "rolling_updates": rolling_updates,
                    "baseline_primes": baseline_primes,
                    "rolling_primes": rolling_primes,
                    "baseline_intersections": (
                        baseline_intersections
                    ),
                    "rolling_intersections": (
                        rolling_intersections
                    ),
                    "baseline_survivors": (
                        baseline_survivors
                    ),
                    "rolling_survivors": (
                        rolling_survivors
                    ),
                    "baseline_exact": baseline_exact,
                    "rolling_exact": rolling_exact,
                    "skipped": False,
                }
            )

    enumerated = [
        r
        for r in results
        if not r["skipped"]
    ]

    skipped = [
        r
        for r in results
        if r["skipped"]
    ]

    # =========================================================================
    # COMPACT
    # =========================================================================

    print()
    print("=" * 120)
    print("COMPACT ROLLING-SIEVE SUMMARY")
    print("=" * 120)

    print(
        " i   u       gap        x    f-cand "
        "baseline-checks rolling-checks updates exact-uniq"
    )

    print("-" * 120)

    for r in results:

        if r["skipped"]:

            print(
                f"{r['index']:2d} "
                f"{r['u']:3d} "
                f"{r['inst'].gap:10d} "
                f"{r['x_true']:9d} "
                f"{r['count']:9d} "
                f"SKIPPED"
            )

            continue

        unique = (
            len(r["rolling_exact"])
            == 1
        )

        print(
            f"{r['index']:2d} "
            f"{r['u']:3d} "
            f"{r['inst'].gap:10d} "
            f"{r['x_true']:9d} "
            f"{r['count']:9d} "
            f"{r['baseline_checks']:14d} "
            f"{r['rolling_checks']:14d} "
            f"{r['rolling_updates']:7d} "
            f"{'YES' if unique else 'NO':>10s}"
        )

    # =========================================================================
    # GLOBAL
    # =========================================================================

    total_candidates = sum(
        r["count"]
        for r in enumerated
    )

    baseline_checks = sum(
        r["baseline_checks"]
        for r in enumerated
    )

    rolling_checks = sum(
        r["rolling_checks"]
        for r in enumerated
    )

    baseline_square_evaluations = sum(
        r["baseline_squares"]
        for r in enumerated
    )

    rolling_initial_evaluations = sum(
        r["rolling_initial"]
        for r in enumerated
    )

    rolling_updates = sum(
        r["rolling_updates"]
        for r in enumerated
    )

    exact_unique = sum(
        len(r["rolling_exact"]) == 1
        for r in enumerated
    )

    survivor_sets_equal = all(
        {
            c[0]
            for c in r["baseline_survivors"]
        }
        ==
        {
            c[0]
            for c in r["rolling_survivors"]
        }
        for r in enumerated
    )

    exact_sets_equal = all(
        {
            c[0]
            for c in r["baseline_exact"]
        }
        ==
        {
            c[0]
            for c in r["rolling_exact"]
        }
        for r in enumerated
    )

    all_true_survive = all(
        any(
            c[0] == r["inst"].d
            for c in r["rolling_survivors"]
        )
        for r in enumerated
    )

    all_products = all(
        (
            exact_reconstruct(
                r["inst"],
                c[0],
            )[0]
            *
            exact_reconstruct(
                r["inst"],
                c[0],
            )[1]
            ==
            r["inst"].N
        )
        for r in enumerated
        for c in r["rolling_exact"]
    )

    print()
    print("=" * 120)
    print("GLOBAL ROLLING-SIEVE SUMMARY")
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
        f"  total f-candidates                  = "
        f"{total_candidates}"
    )

    print(
        f"  baseline discriminant evaluations   = "
        f"{baseline_square_evaluations}"
    )

    print(
        f"  rolling initial discriminants       = "
        f"{rolling_initial_evaluations}"
    )

    print(
        f"  rolling discriminant updates        = "
        f"{rolling_updates}"
    )

    print(
        f"  baseline modular candidate checks   = "
        f"{baseline_checks}"
    )

    print(
        f"  rolling modular candidate checks    = "
        f"{rolling_checks}"
    )

    print(
        f"  baseline/rolling check ratio        = "
        f"{baseline_checks / rolling_checks:.8f}"
    )

    print(
        f"  exact unique recoveries             = "
        f"{exact_unique}/{len(enumerated)}"
    )

    print(
        f"  baseline/rolling survivor sets same = "
        f"{survivor_sets_equal}"
    )

    print(
        f"  baseline/rolling exact sets same    = "
        f"{exact_sets_equal}"
    )

    # =========================================================================
    # AGGREGATE
    # =========================================================================

    print()
    print("=" * 120)
    print("AGGREGATE BY UNKNOWN K BITS")
    print("=" * 120)

    print(
        " u   tests avg-f baseline-checks "
        "rolling-checks updates reduction"
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

        avg_f = (
            sum(
                r["count"]
                for r in group
            )
            / len(group)
        )

        bc = sum(
            r["baseline_checks"]
            for r in group
        )

        rc = sum(
            r["rolling_checks"]
            for r in group
        )

        up = sum(
            r["rolling_updates"]
            for r in group
        )

        reduction = (
            1.0
            - rc / bc
            if bc
            else 0.0
        )

        print(
            f"{u:3d} "
            f"{len(group):7d} "
            f"{avg_f:9.2f} "
            f"{bc:14d} "
            f"{rc:13d} "
            f"{up:7d} "
            f"{reduction:9.6f}"
        )

    # =========================================================================
    # LARGEST CASES
    # =========================================================================

    print()
    print("=" * 120)
    print("ROLLING-UPDATE DETAILS")
    print("=" * 120)

    largest = sorted(
        enumerated,
        key=lambda r: r["rolling_updates"],
        reverse=True,
    )[:20]

    for r in largest:

        print(
            f"  instance={r['index']:2d} "
            f"u={r['u']:2d} "
            f"f={r['count']:7d} "
            f"baseline={r['baseline_checks']:8d} "
            f"rolling={r['rolling_checks']:8d} "
            f"updates={r['rolling_updates']:8d} "
            f"exact={len(r['rolling_exact'])}"
        )

    # =========================================================================
    # IMPORTANT DETAIL
    # =========================================================================

    print()
    print("=" * 120)
    print("ARITHMETIC KERNEL COMPARISON")
    print("=" * 120)

    print(
        "  Direct implementation:"
    )

    print(
        "    D4 = (N+1-d)^2 - 16N"
    )

    print(
        "    one full square is evaluated for every"
    )

    print(
        "    candidate/modulus pair."
    )

    print()

    print(
        "  Rolling implementation:"
    )

    print(
        "    D4' = D4 - 2*a*delta + delta^2 (mod p)"
    )

    print(
        "    only the first candidate performs a square"
    )

    print(
        "    initialization for each modulus."
    )

    print(
        "    Every later candidate uses the recurrence."
    )

    print()

    print(
        f"  baseline square evaluations = "
        f"{baseline_square_evaluations}"
    )

    print(
        f"  rolling initial squares     = "
        f"{rolling_initial_evaluations}"
    )

    print(
        f"  rolling updates             = "
        f"{rolling_updates}"
    )

    print(
        "  theoretical square-operation "
        "reduction = "
        f"{1.0 - rolling_initial_evaluations / baseline_square_evaluations:.8f}"
    )

    # =========================================================================
    # CHECKS
    # =========================================================================

    print()
    print("=" * 120)
    print("EXACT INTERNAL CHECKS")
    print("=" * 120)

    print(
        f"  all rolling masks equal baseline masks = "
        f"{survivor_sets_equal}"
    )

    print(
        f"  all exact survivor sets equal          = "
        f"{exact_sets_equal}"
    )

    print(
        f"  all enumerated true d values survive  = "
        f"{all_true_survive}"
    )

    print(
        f"  all reconstructed products equal N    = "
        f"{all_products}"
    )

    # =========================================================================
    # INTERPRETATION
    # =========================================================================

    print()
    print("=" * 120)
    print("INTERPRETATION")
    print("=" * 120)

    print(
        "  Experiments 424-428 optimized the representation"
    )

    print(
        "  and selection strategy of the modular sieve."
    )

    print()

    print(
        "  Experiment 429 changes the arithmetic kernel."
    )

    print(
        "  The discriminant is a quadratic polynomial in d:"
    )

    print(
        "    D4(d) = (N+1-d)^2 - 16N."
    )

    print()

    print(
        "  Therefore consecutive candidate values can be"
    )

    print(
        "  evaluated through an exact modular recurrence"
    )

    print(
        "  instead of repeatedly recomputing the square."
    )

    print()

    print(
        "  This does not alter the mathematical sieve."
    )

    print(
        "  It only changes how the same necessary condition"
    )

    print(
        "  is evaluated."
    )

    print()

    print(
        "  The exact factor reconstruction remains the"
    )

    print(
        "  authoritative acceptance condition."
    )

    print()

    print(
        "  This remains an implementation/arithmetic"
    )

    print(
        "  optimization experiment, not a factoring theorem."
    )

    # =========================================================================
    # FINAL
    # =========================================================================

    print()
    print("=" * 120)
    print("EXPERIMENT 429 FINAL STATUS")
    print("=" * 120)

    print(
        f"  ENUMERATED TESTS = "
        f"{len(enumerated)}"
    )

    print(
        f"  EXACT UNIQUE RECOVERIES = "
        f"{exact_unique}/{len(enumerated)}"
    )

    print(
        f"  TOTAL F-CANDIDATES = "
        f"{total_candidates}"
    )

    print(
        f"  BASELINE MODULAR CHECKS = "
        f"{baseline_checks}"
    )

    print(
        f"  ROLLING MODULAR CHECKS = "
        f"{rolling_checks}"
    )

    print(
        f"  BASELINE/ROLLING CHECK RATIO = "
        f"{baseline_checks / rolling_checks:.8f}"
    )

    print(
        f"  ROLLING DISCRIMINANT UPDATES = "
        f"{rolling_updates}"
    )

    print(
        f"  ALL INTERNAL EXACT CHECKS = "
        f"{survivor_sets_equal and exact_sets_equal and all_true_survive and all_products}"
    )

    print("=" * 120)
    print("EXPERIMENT 429 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    run()
