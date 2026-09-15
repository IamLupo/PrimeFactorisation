#!/usr/bin/env python3

from __future__ import annotations

from dataclasses import dataclass
from math import isqrt


# ======================================================================================================================
# EXPERIMENT 441
# ======================================================================================================================

TITLE = "INTRINSIC K-BOUND / DISCRIMINANT GEOMETRY AUDIT"

K_BITS = [32, 40, 48, 56, 64, 72, 80, 88]
X_RADIUS_MULTIPLIERS = [1, 2, 4, 8, 16, 32]

INSTANCES = [
    (50411, 282599, -10),
    (1013, 10009, -17),
    (10009, 1000033, -24),
    (10009, 10037, -31),
    (50023, 50051, -38),
    (100019, 100043, -45),
    (200009, 200017, -52),
    (300017, 900007, -59),
]

HIDDEN_K_FACTORS = [17, 29, 43, 61, 89, 127, 191, 257]


# ======================================================================================================================
# INSTANCE
# ======================================================================================================================

@dataclass(frozen=True)
class Instance:
    idx: int
    p: int
    q: int
    N: int
    S: int
    d: int
    x_true: int
    d0: int


def build_instances() -> list[Instance]:
    out = []

    for idx, (p, q, x_true) in enumerate(INSTANCES, start=1):
        N = p * q
        S = p + q
        d = N - 2 * S + 1
        d0 = d - x_true

        out.append(
            Instance(
                idx=idx,
                p=p,
                q=q,
                N=N,
                S=S,
                d=d,
                x_true=x_true,
                d0=d0,
            )
        )

    return out


# ======================================================================================================================
# HIDDEN TEST GENERATION
# ======================================================================================================================

def hidden_K(inst: Instance) -> int:
    factor = HIDDEN_K_FACTORS[inst.idx - 1]

    scale = (
        inst.d0 * inst.d0
        + 13 * inst.N
        + 7 * inst.S * inst.S
        + factor * inst.d0
    )

    return -(scale * factor + inst.N * (factor + 3))


def build_r(inst: Instance, K_true: int) -> int:
    x = inst.x_true
    return 4 * K_true + x * x + 2 * inst.d0 * x


# ======================================================================================================================
# EXACT K GEOMETRY
# ======================================================================================================================

def A_value(inst: Instance, r: int) -> int:
    """
    Completing the square:

        4K = r + d0^2 - (x+d0)^2

    so

        4K = A - y^2

    with y=x+d0.
    """
    return r + inst.d0 * inst.d0


def K_from_y(A: int, y: int) -> int:
    numerator = A - y * y
    assert numerator % 4 == 0
    return numerator // 4


def true_y(inst: Instance) -> int:
    return inst.x_true + inst.d0


# ======================================================================================================================
# DERIVED X BOUNDS
# ======================================================================================================================

def intrinsic_x_scales(inst: Instance) -> list[tuple[str, int]]:
    """
    Produce bounds that use only quantities already known from the instance.

    These are deliberately NOT factor-derived.

    radius means |x| <= radius.
    """
    scales = []

    # Directly visible scales.
    scales.append(("S", max(1, inst.S)))
    scales.append(("sqrtN", max(1, isqrt(inst.N))))
    scales.append(("gap_proxy", max(1, abs(inst.d - inst.S))))
    scales.append(("d", max(1, abs(inst.d))))

    # Multiples of sqrt(N).
    rootN = max(1, isqrt(inst.N))
    for m in X_RADIUS_MULTIPLIERS:
        scales.append((f"{m}*sqrtN", m * rootN))

    return scales


def K_range_from_x_radius(
    inst: Instance,
    r: int,
    radius: int,
) -> tuple[int, int]:
    """
    For |x| <= radius, determine the exact min/max K over the interval.

    Since

        K(x) = (r - x^2 - 2*d0*x)/4

    is a concave quadratic, the maximum is at the vertex clipped to the interval,
    while the minimum is at one of the endpoints.

    We use exact arithmetic only.
    """
    A = A_value(inst, r)

    xs = [-radius, radius]

    vertex = -inst.d0
    if -radius <= vertex <= radius:
        xs.append(vertex)

    values = []

    for x in xs:
        y = x + inst.d0
        num = A - y * y

        if num % 4 != 0:
            # K is only defined on x values where the numerator is divisible by 4.
            # Check nearby x values exactly.
            for xx in range(max(-radius, x - 3), min(radius, x + 3) + 1):
                n2 = A - (xx + inst.d0) ** 2
                if n2 % 4 == 0:
                    values.append(n2 // 4)
        else:
            values.append(num // 4)

    # Also inspect the nearest integral x values around the vertex.
    for x in range(vertex - 3, vertex + 4):
        if -radius <= x <= radius:
            num = A - (x + inst.d0) ** 2
            if num % 4 == 0:
                values.append(num // 4)

    if not values:
        raise RuntimeError("No integer K values found in x interval")

    return min(values), max(values)


# ======================================================================================================================
# REQUIRED BIT WIDTH
# ======================================================================================================================

def signed_bits_for_range(lo: int, hi: int) -> int:
    """
    Smallest B such that [-2^B, 2^B) contains [lo,hi].
    """
    m = max(abs(lo), abs(hi))

    if m == 0:
        return 1

    return m.bit_length()


def bounded_interval_contains(lo: int, hi: int, K: int) -> bool:
    return lo <= K < hi


# ======================================================================================================================
# REPORTING
# ======================================================================================================================

def print_header(inst: Instance, K_true: int, r: int) -> None:
    print("-" * 120)
    print(
        f"INSTANCE {inst.idx}: "
        f"p={inst.p} q={inst.q} N={inst.N} S={inst.S} "
        f"d={inst.d} d0={inst.d0} x_true={inst.x_true}"
    )
    print(f"  true_y = {true_y(inst)}")
    print(f"  hidden K magnitude bits = {abs(K_true).bit_length()}")
    print(f"  A = r+d0^2 = {A_value(inst, r)}")
    print()


def run_instance(inst: Instance) -> dict:
    K_true = hidden_K(inst)
    r = build_r(inst, K_true)

    print_header(inst, K_true, r)

    results = []

    print(
        "  NATURAL X-BOUND AUDIT"
    )
    print(
        "    bound                    radius              K-min                K-max"
    )
    print("    " + "-" * 100)

    for name, radius in intrinsic_x_scales(inst):
        lo, hi = K_range_from_x_radius(inst, r, radius)
        bits = signed_bits_for_range(lo, hi)
        contains = bounded_interval_contains(lo, hi, K_true)

        print(
            f"    {name:<20}"
            f" {radius:>18}"
            f" {lo:>20}"
            f" {hi:>20}"
            f" bits={bits:>3}"
            f" contains_true={str(contains):>5}"
        )

        results.append(
            {
                "name": name,
                "radius": radius,
                "lo": lo,
                "hi": hi,
                "bits": bits,
                "contains_true": contains,
            }
        )

    print()
    print("  REQUIRED K-BOUND AUDIT")
    print(
        "    B    interval_low               interval_high              true-in-range"
    )
    print("    " + "-" * 90)

    for B in K_BITS:
        lo = -(1 << B)
        hi = 1 << B
        inside = lo <= K_true < hi

        print(
            f"    {B:>2}"
            f" {lo:>24}"
            f" {hi:>24}"
            f" {str(inside):>15}"
        )

    # Verify the exact relation independently.
    recovered = K_from_y(A_value(inst, r), true_y(inst))

    assert recovered == K_true

    return {
        "K_true": K_true,
        "r": r,
        "natural": results,
        "true_relation": recovered == K_true,
    }


# ======================================================================================================================
# MAIN
# ======================================================================================================================

def main() -> None:
    print("=" * 120)
    print("EXPERIMENT 441")
    print("=" * 120)
    print()
    print("INTRINSIC K-BOUND / DISCRIMINANT GEOMETRY AUDIT")
    print()
    print("QUESTION")
    print("  Can the known N,S,d,d0 quantities imply a substantially")
    print("  tighter K bound than the artificial [-2^B,2^B) intervals?")
    print()
    print("RULES")
    print("  exact integer arithmetic only")
    print("  no resultants")
    print("  no Groebner basis")
    print("  no symbolic factorization")
    print("  no giant K image")
    print("  no giant x enumeration")
    print()

    instances = build_instances()

    all_true = True
    total_natural_bounds = 0
    useful_true_containing_bounds = 0

    global_bits = []

    for inst in instances:
        result = run_instance(inst)

        if not result["true_relation"]:
            all_true = False

        for row in result["natural"]:
            total_natural_bounds += 1

            if row["contains_true"]:
                useful_true_containing_bounds += 1
                global_bits.append(row["bits"])

    # ------------------------------------------------------------------------------------------------------------------
    # GLOBAL SUMMARY
    # ------------------------------------------------------------------------------------------------------------------

    print()
    print("=" * 120)
    print("GLOBAL NATURAL-BOUND SUMMARY")
    print("=" * 120)

    print(f"  instances                         = {len(instances)}")
    print(f"  natural bounds tested             = {total_natural_bounds}")
    print(f"  true K contained                  = {useful_true_containing_bounds}")
    print(f"  exact true-relation checks        = {len(instances)}/{len(instances)}")

    if global_bits:
        print(
            f"  smallest natural-bound K bits    = {min(global_bits)}"
        )
        print(
            f"  largest natural-bound K bits     = {max(global_bits)}"
        )
        print(
            f"  average natural-bound K bits     = "
            f"{sum(global_bits) / len(global_bits):.3f}"
        )

    print()
    print("INTERPRETATION")
    print()
    print(
        "  This experiment separates two questions:"
    )
    print(
        "    (1) What K range would be sufficient if x were bounded?"
    )
    print(
        "    (2) Can that x bound itself be justified from known data?"
    )
    print()
    print(
        "  A genuinely useful result would be a small natural x radius"
    )
    print(
        "  whose resulting K interval contains the true K."
    )
    print()
    print(
        "  If every reasonable known-data x bound still requires a"
    )
    print(
        "  K range comparable to the hidden K magnitude, then the"
    )
    print(
        "  K reconstruction route remains blocked at the bound itself."
    )
    print()
    print(
        "  This is still an information/bound audit, not a factorization theorem."
    )

    print()
    print("=" * 120)
    print("EXPERIMENT 441 FINAL STATUS")
    print("=" * 120)
    print(
        f"  EXACT TRUE-RELATION CHECKS = {all_true}"
    )
    print(
        f"  NATURAL BOUNDS TESTED = {total_natural_bounds}"
    )
    print(
        f"  TRUE-CONTAINING NATURAL BOUNDS = {useful_true_containing_bounds}"
    )
    print("=" * 120)
    print("EXPERIMENT 441 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    main()
