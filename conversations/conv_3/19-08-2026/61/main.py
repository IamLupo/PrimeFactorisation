#!/usr/bin/env python3

from __future__ import annotations

from dataclasses import dataclass
from math import isqrt
from typing import List, Tuple, Optional


# =============================================================================
# EXPERIMENT 438R2
# =============================================================================
#
# MEMORY-SAFE BOUNDED-K RECONSTRUCTION
#
# Goal:
#   Determine whether the reachable 2-adic K-image becomes sufficient to
#   identify the actual integer K inside a finite integer range.
#
# Important:
#   We do NOT assume that |K| is bounded by N*S.
#
#   Instead we explicitly test ranges:
#
#       -2^B <= K < 2^B
#
#   for several B.
#
# Memory policy:
#   - never materialize the complete bounded K population;
#   - count residue lifts analytically;
#   - materialize only when the total count is small.
#
# Exact arithmetic only.
# =============================================================================


EXPERIMENT_NAME = "EXPERIMENT 438R2"

ENUMERATION_LIMIT = 2_000_000

# Only materialize actual integer K values when the population is small.
MATERIALIZE_MAX = 50_000

# 2-adic levels already validated by Experiment 437.
U_VALUES = [8, 10, 12, 14, 16, 18, 20]

# Explicit K-bit bounds.
#
# B=64 means:
#       -2^64 <= K < 2^64
#
# Larger bounds are measured analytically without constructing all K.
K_BOUND_BITS = [
    32,
    40,
    48,
    56,
    64,
    72,
    80,
]


# =============================================================================
# BASE INSTANCES
# =============================================================================

BASE_FACTORS = [
    (50411, 282599),
    (1013, 10009),
    (10009, 1000033),
    (10009, 10037),
    (50023, 50051),
    (100019, 100043),
    (200009, 200017),
    (300017, 900007),
]


@dataclass(frozen=True)
class Instance:
    index: int
    p: int
    q: int
    N: int
    S: int
    d: int
    d0: int
    x_true: int
    K_true: int
    r: int


def build_instances() -> List[Instance]:
    out: List[Instance] = []

    for idx, (p, q) in enumerate(BASE_FACTORS, start=1):

        N = p * q
        S = p + q

        # Same d construction used in the previous experiments.
        d = N + 1 - 2 * S

        # Deterministic d0 shift.
        d0 = d + (7 * idx + 3)

        x_true = d - d0

        # Deterministic hidden K.
        #
        # This is only used internally to validate the generated test.
        # The reconstruction measurements do NOT use K_true.
        K_true = -(
            (N + 101 * idx)
            * (S + 37 * idx)
            * (idx + 11)
        )

        r = (
            x_true * x_true
            + 2 * d0 * x_true
            + 4 * K_true
        )

        assert (
            x_true * x_true
            + 2 * d0 * x_true
            + 4 * K_true
            - r
            == 0
        )

        out.append(
            Instance(
                index=idx,
                p=p,
                q=q,
                N=N,
                S=S,
                d=d,
                d0=d0,
                x_true=x_true,
                K_true=K_true,
                r=r,
            )
        )

    return out


INSTANCES = build_instances()


# =============================================================================
# CORE
# =============================================================================

def f_value(inst: Instance, x: int, K: int) -> int:
    return (
        x * x
        + 2 * inst.d0 * x
        + 4 * K
        - inst.r
    )


def K_from_x(inst: Instance, x: int) -> int:
    numerator = (
        inst.r
        - x * x
        - 2 * inst.d0 * x
    )

    if numerator % 4 != 0:
        raise ValueError("K is not integral for this x")

    return numerator // 4


# =============================================================================
# 2-ADIC K IMAGE
# =============================================================================

def compute_K_image(
    inst: Instance,
    u: int,
) -> Tuple[int, List[int]]:
    """
    Compute reachable K residues modulo 2^(u-2).

    We enumerate x modulo 2^u only when that is below the configured
    enumeration limit.
    """

    x_modulus = 1 << u

    if x_modulus > ENUMERATION_LIMIT:
        raise RuntimeError(
            f"2^{u}={x_modulus} exceeds enumeration limit"
        )

    k_modulus = 1 << (u - 2)

    image = set()

    for x in range(x_modulus):

        numerator = (
            inst.r
            - x * x
            - 2 * inst.d0 * x
        )

        if numerator & 3:
            continue

        K = numerator // 4

        image.add(K % k_modulus)

    return k_modulus, sorted(image)


# =============================================================================
# BOUNDED K RANGE
# =============================================================================

def bounded_range(bits: int) -> Tuple[int, int]:
    """
    Inclusive lower, exclusive upper:

        -2^bits <= K < 2^bits
    """

    radius = 1 << bits

    return -radius, radius


# =============================================================================
# ANALYTIC RESIDUE-LIFT COUNT
# =============================================================================

def count_residue_lifts(
    residue: int,
    modulus: int,
    lower: int,
    upper: int,
) -> int:
    """
    Count K in [lower, upper) such that

        K == residue (mod modulus)

    without enumerating K.
    """

    if lower >= upper:
        return 0

    first = lower + ((residue - lower) % modulus)

    if first >= upper:
        return 0

    return 1 + (upper - 1 - first) // modulus


def count_bounded_K(
    image: List[int],
    modulus: int,
    lower: int,
    upper: int,
) -> int:
    """
    Total number of bounded integer K values represented by the
    reachable 2-adic image.
    """

    total = 0

    for residue in image:
        total += count_residue_lifts(
            residue,
            modulus,
            lower,
            upper,
        )

    return total


# =============================================================================
# LAZY K GENERATOR
# =============================================================================

def iter_bounded_K(
    image: List[int],
    modulus: int,
    lower: int,
    upper: int,
):
    """
    Lazily generate bounded K values.

    No global K list/set is constructed.
    """

    for residue in image:

        first = lower + ((residue - lower) % modulus)

        if first >= upper:
            continue

        K = first

        while K < upper:
            yield K
            K += modulus


# =============================================================================
# EXACT QUADRATIC ROOTS
# =============================================================================

def roots_from_K(
    inst: Instance,
    K: int,
) -> List[int]:
    """
    Solve

        x^2 + 2*d0*x + 4K-r = 0

    exactly.

    Equivalent form:

        (x+d0)^2 = d0^2-r-4K
    """

    discriminant = (
        inst.d0 * inst.d0
        - inst.r
        - 4 * K
    )

    if discriminant < 0:
        return []

    root = isqrt(discriminant)

    if root * root != discriminant:
        return []

    return sorted({
        -inst.d0 + root,
        -inst.d0 - root,
    })


# =============================================================================
# FACTOR RECONSTRUCTION
# =============================================================================

def reconstruct_factors(
    N: int,
    d: int,
) -> Optional[Tuple[int, int]]:
    """
    Exact reconstruction from

        d = N+1-(p+q).
    """

    numerator = N + 1 - d

    if numerator & 1:
        return None

    S = numerator // 2

    if S <= 0:
        return None

    discriminant = S * S - 4 * N

    if discriminant < 0:
        return None

    gap = isqrt(discriminant)

    if gap * gap != discriminant:
        return None

    if (S - gap) & 1:
        return None

    p = (S - gap) // 2
    q = (S + gap) // 2

    if p <= 1 or q <= 1:
        return None

    if p * q != N:
        return None

    return p, q


# =============================================================================
# CASE
# =============================================================================

def run_case(
    inst: Instance,
    u: int,
    bound_bits: int,
    materialize: bool = True,
):

    k_modulus, image = compute_K_image(
        inst,
        u,
    )

    lower, upper = bounded_range(
        bound_bits
    )

    # -------------------------------------------------------------------------
    # This is only a validation invariant.
    # We do NOT require every arbitrary bound to contain K_true.
    # -------------------------------------------------------------------------

    true_in_range = (
        lower <= inst.K_true < upper
    )

    true_residue = inst.K_true % k_modulus

    true_residue_reachable = (
        true_residue in image
    )

    # Experiment 437 tells us the latter should always hold.
    assert true_residue_reachable

    bounded_count = count_bounded_K(
        image,
        k_modulus,
        lower,
        upper,
    )

    materialized_count = 0
    exact_root_candidates = 0
    exact_factor_candidates = 0

    unique_K = False
    unique_root = False
    unique_factor = False

    factor_records = []

    # -------------------------------------------------------------------------
    # Only materialize when small.
    # -------------------------------------------------------------------------

    if materialize and bounded_count <= MATERIALIZE_MAX:

        materialized_count = bounded_count

        seen_K = []

        for K in iter_bounded_K(
            image,
            k_modulus,
            lower,
            upper,
        ):
            seen_K.append(K)

        # Exact consistency.
        assert len(seen_K) == bounded_count

        unique_K = (
            len(seen_K) == 1
        )

        seen_roots = set()

        for K in seen_K:

            roots = roots_from_K(
                inst,
                K,
            )

            for x in roots:

                assert f_value(
                    inst,
                    x,
                    K,
                ) == 0

                seen_roots.add(
                    (K, x)
                )

                d = inst.d0 + x

                pair = reconstruct_factors(
                    inst.N,
                    d,
                )

                if pair is None:
                    continue

                factor_records.append(
                    (
                        K,
                        x,
                        d,
                        pair,
                    )
                )

        exact_root_candidates = len(
            seen_roots
        )

        unique_root = (
            exact_root_candidates == 1
        )

        exact_factor_candidates = len(
            factor_records
        )

        unique_factor = (
            exact_factor_candidates == 1
        )

    return {
        "u": u,
        "bound_bits": bound_bits,
        "Kmod": k_modulus,
        "image_size": len(image),
        "K_lower": lower,
        "K_upper": upper,
        "bounded_count": bounded_count,
        "true_in_range": true_in_range,
        "true_residue_reachable": true_residue_reachable,
        "materialized_count": materialized_count,
        "exact_root_candidates": exact_root_candidates,
        "exact_factor_candidates": exact_factor_candidates,
        "unique_K": unique_K,
        "unique_root": unique_root,
        "unique_factor": unique_factor,
        "factor_records": factor_records,
    }


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 120)
    print(EXPERIMENT_NAME)
    print("=" * 120)
    print()
    print("MEMORY-SAFE BOUNDED-K RECONSTRUCTION")
    print()
    print("K RANGE MODEL")
    print("  -2^B <= K < 2^B")
    print()
    print(
        f"Materialization threshold = "
        f"{MATERIALIZE_MAX:,} K values"
    )
    print()

    results = []

    for inst in INSTANCES:

        print("-" * 120)

        print(
            f"INSTANCE {inst.index}: "
            f"N={inst.N} "
            f"d={inst.d} "
            f"d0={inst.d0} "
            f"x_true={inst.x_true}"
        )

        # This is deliberately printed only as validation information.
        # The reconstruction itself does not use K_true.
        true_abs_bits = max(
            1,
            abs(inst.K_true).bit_length(),
        )

        print(
            f"  hidden-test K magnitude bits = "
            f"{true_abs_bits}"
        )

        for u in U_VALUES:

            if (1 << u) > ENUMERATION_LIMIT:
                print(
                    f"  u={u:2d} "
                    f"SKIPPED "
                    f"(2^{u} > {ENUMERATION_LIMIT})"
                )
                continue

            print()
            print(
                f"  u={u:2d}"
            )

            for B in K_BOUND_BITS:

                result = run_case(
                    inst,
                    u,
                    B,
                )

                results.append(
                    (inst, result)
                )

                status = (
                    "IN"
                    if result["true_in_range"]
                    else "OUT"
                )

                if result["bounded_count"] <= MATERIALIZE_MAX:
                    materialized = "YES"
                else:
                    materialized = "NO"

                print(
                    f"    B={B:2d} "
                    f"Kmod={result['Kmod']:>8} "
                    f"image={result['image_size']:>7} "
                    f"bounded={result['bounded_count']:>12} "
                    f"true={status:>3} "
                    f"materialized={materialized:>3} "
                    f"roots={result['exact_root_candidates']:>4} "
                    f"factors={result['exact_factor_candidates']:>4} "
                    f"uniqueK={str(result['unique_K']):>5} "
                    f"uniqueFactor={str(result['unique_factor']):>5}"
                )

    # =========================================================================
    # GLOBAL SUMMARY
    # =========================================================================

    print()
    print("=" * 120)
    print("GLOBAL BOUNDED-K RECONSTRUCTION SUMMARY")
    print("=" * 120)

    evaluated = len(results)

    true_in_range_cases = sum(
        r["true_in_range"]
        for _, r in results
    )

    reachable_cases = sum(
        r["true_residue_reachable"]
        for _, r in results
    )

    unique_K_cases = sum(
        r["unique_K"]
        for _, r in results
    )

    unique_factor_cases = sum(
        r["unique_factor"]
        for _, r in results
    )

    materialized_cases = sum(
        r["materialized_count"] > 0
        or r["bounded_count"] == 0
        for _, r in results
    )

    print(
        f"  evaluated cases                 = {evaluated}"
    )

    print(
        f"  true K residue reachable        = "
        f"{reachable_cases}/{evaluated}"
    )

    print(
        f"  cases where true K is in range  = "
        f"{true_in_range_cases}/{evaluated}"
    )

    print(
        f"  unique bounded K reconstructions= "
        f"{unique_K_cases}/{evaluated}"
    )

    print(
        f"  unique factor reconstructions   = "
        f"{unique_factor_cases}/{evaluated}"
    )

    print()
    print("IMPORTANT")

    print(
        "  A bound with true=OUT is not a failed mathematical test."
    )

    print(
        "  It simply means that particular artificial K range"
    )

    print(
        "  does not contain the generated hidden K."
    )

    print(
        "  Such a range must not be interpreted as a reconstruction failure."
    )

    print()
    print("SMALLEST VALID BOUNDS")

    # For each instance/u find the first tested B containing K_true.
    for inst in INSTANCES:

        print()
        print(
            f"  INSTANCE {inst.index}"
        )

        for u in U_VALUES:

            candidates = [
                r
                for i, r in results
                if i.index == inst.index
                and r["u"] == u
                and r["true_in_range"]
            ]

            if not candidates:
                continue

            best = min(
                candidates,
                key=lambda r: r["bound_bits"]
            )

            print(
                f"    u={u:2d} "
                f"first-containing-B={best['bound_bits']:2d} "
                f"bounded-K={best['bounded_count']:>12} "
                f"image={best['image_size']:>7} "
                f"uniqueK={best['unique_K']} "
                f"uniqueFactor={best['unique_factor']}"
            )

    print()
    print("SMALLEST BOUNDED POPULATIONS")

    ranked = sorted(
        results,
        key=lambda z: (
            z[1]["bounded_count"],
            z[1]["bound_bits"],
            z[1]["u"],
            z[0].index,
        ),
    )

    for inst, r in ranked[:50]:

        print(
            f"  inst={inst.index} "
            f"u={r['u']:2d} "
            f"B={r['bound_bits']:2d} "
            f"image={r['image_size']:>7} "
            f"bounded={r['bounded_count']:>12} "
            f"true={r['true_in_range']} "
            f"uniqueK={r['unique_K']} "
            f"factor={r['unique_factor']}"
        )

    print()
    print("INTERPRETATION")
    print()
    print(
        "  Experiment 437 showed that the 2-adic K-image has"
    )
    print(
        "  approximately one reachable residue for every six"
    )
    print(
        "  possible residues asymptotically."
    )
    print()
    print(
        "  Experiment 438R2 now asks whether a finite integer"
    )
    print(
        "  K interval converts that residue information into"
    )
    print(
        "  an actual K determination."
    )
    print()
    print(
        "  The critical quantity is bounded_count:"
    )
    print(
        "      number of integer K values in the range whose"
    )
    print(
        "      residue modulo 2^(u-2) belongs to the 2-adic image."
    )
    print()
    print(
        "  If bounded_count = 1, the bounded K is uniquely"
    )
    print(
        "  determined by the 2-adic image."
    )
    print()
    print(
        "  If bounded_count > 1, the 2-adic information alone"
    )
    print(
        "  does not identify the integer K inside that range."
    )
    print()
    print(
        "  Exact quadratic-root and factor reconstruction are"
    )
    print(
        "  performed only for small candidate populations."
    )
    print()
    print(
        "  No giant K set is constructed."
    )

    print()
    print("=" * 120)
    print("EXPERIMENT 438R2 FINAL STATUS")
    print("=" * 120)

    print(
        f"  EVALUATED CASES = {evaluated}"
    )

    print(
        f"  TRUE RESIDUE REACHABILITY = "
        f"{reachable_cases}/{evaluated}"
    )

    print(
        f"  UNIQUE BOUNDED-K CASES = "
        f"{unique_K_cases}/{evaluated}"
    )

    print(
        f"  UNIQUE FACTOR CASES = "
        f"{unique_factor_cases}/{evaluated}"
    )

    print(
        "  GIANT-K-MATERIALIZATION = False"
    )

    print(
        "  ALL INTERNAL EXACT CHECKS = True"
    )

    print("=" * 120)
    print("EXPERIMENT 438R2 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    main()