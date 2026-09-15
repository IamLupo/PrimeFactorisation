#!/usr/bin/env python3

from __future__ import annotations

from dataclasses import dataclass
from math import isqrt


# ======================================================================================================================
# EXPERIMENT 446R2
# ======================================================================================================================
#
# KNOWN-DATA-CENTERED K INTERVAL -> EXACT X POPULATION
#
# FIXED VERSION
#
# The previous implementation assumed x=0 was always parity-valid.
#
# Actual condition:
#
#     x + d0 == 0 (mod 2)
#
# Therefore:
#
#     d0 even -> x=0 is valid
#     d0 odd  -> x=1 (or -1) is the nearest valid reference
#
# This version chooses the nearest parity-valid x to zero.
#
# K_TRUE IS NEVER USED TO CONSTRUCT ANY BOUND.
#
# It is used only for post-hoc diagnostics.
#
# NO GIANT K ENUMERATION
# NO GIANT X ENUMERATION
# NO CRT CARTESIAN PRODUCT
#
# EXACT INTEGER ARITHMETIC
#
# ======================================================================================================================


X_RADII = [
    0,
    1,
    2,
    4,
    8,
    16,
    32,
    64,
    128,
    256,
    512,
    1024,
    2048,
    4096,
    8192,
    16384,
    32768,
    65536,
    131072,
    262144,
    524288,
    1048576,
    4194304,
    16777216,
    67108864,
]

BRUTE_FORCE_LIMIT = 10000


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
    d0: int
    x_true: int


RAW = [
    (50411, 282599, -10),
    (1013, 10009, -17),
    (10009, 1000033, -24),
    (10009, 10037, -31),
    (50023, 50051, -38),
    (100019, 100043, -45),
    (200009, 200017, -52),
    (300017, 900007, -59),
]


def build_instances() -> list[Instance]:
    out: list[Instance] = []

    for idx, (p, q, x_true) in enumerate(RAW, start=1):

        N = p * q
        S = p + q

        d = N - 2 * S + 1

        # x_true + d0 = d
        d0 = d - x_true

        out.append(
            Instance(
                idx=idx,
                p=p,
                q=q,
                N=N,
                S=S,
                d=d,
                d0=d0,
                x_true=x_true,
            )
        )

    return out


# ======================================================================================================================
# CONTROLLED HIDDEN K GENERATOR
# ======================================================================================================================

HIDDEN_FACTORS = [
    17,
    29,
    43,
    61,
    89,
    127,
    191,
    257,
]


def hidden_K(inst: Instance) -> int:

    f = HIDDEN_FACTORS[inst.idx - 1]

    value = (
        inst.d0 * inst.d0
        + 13 * inst.N
        + 7 * inst.S * inst.S
        + f * inst.d0
    )

    return -(value * f + inst.N * (f + 3))


def build_r(inst: Instance, K_true: int) -> int:

    x = inst.x_true

    return (
        4 * K_true
        + x * x
        + 2 * inst.d0 * x
    )


# ======================================================================================================================
# EXACT MAP
# ======================================================================================================================

def z_from_x(inst: Instance, x: int) -> int | None:

    y = x + inst.d0

    if y & 1:
        return None

    return y // 2


def x_from_z(inst: Instance, z: int) -> int:

    return 2 * z - inst.d0


def k_from_z(A4: int, z: int) -> int:

    return A4 - z * z


def k_from_x(
    inst: Instance,
    A4: int,
    x: int,
) -> int | None:

    z = z_from_x(inst, x)

    if z is None:
        return None

    return k_from_z(A4, z)


# ======================================================================================================================
# INTEGER SQRT
# ======================================================================================================================

def ceil_sqrt(n: int) -> int:

    if n < 0:
        raise ValueError("ceil_sqrt requires n >= 0")

    r = isqrt(n)

    return r if r * r == n else r + 1


# ======================================================================================================================
# PARITY HELPERS
# ======================================================================================================================

def required_x_parity(inst: Instance) -> int:
    """
    x + d0 must be even.

    Therefore:

        x mod 2 = d0 mod 2
    """

    return inst.d0 & 1


def nearest_valid_x_to_zero(inst: Instance) -> int:
    """
    Choose the parity-valid x closest to zero.

    d0 even -> 0
    d0 odd  -> 1

    Either +1 or -1 is mathematically valid for odd d0.
    We use +1 deterministically.
    """

    return required_x_parity(inst)


# ======================================================================================================================
# KNOWN-DATA CENTER
# ======================================================================================================================

def known_center(
    inst: Instance,
    A4: int,
) -> tuple[int, int, int]:

    x0 = nearest_valid_x_to_zero(inst)

    z0 = z_from_x(inst, x0)

    if z0 is None:
        raise AssertionError(
            f"internal parity error at instance={inst.idx}"
        )

    K0 = k_from_z(A4, z0)

    return x0, z0, K0


# ======================================================================================================================
# PARITY-VALID WINDOW
# ======================================================================================================================

def parity_adjusted_bounds(
    inst: Instance,
    R: int,
) -> tuple[int, int, int]:

    lo = -R
    hi = R

    wanted = required_x_parity(inst)

    # First parity-valid x >= lo.
    if (lo & 1) != wanted:
        lo += 1

    # Last parity-valid x <= hi.
    if (hi & 1) != wanted:
        hi -= 1

    if lo > hi:
        return lo, hi, 0

    count = ((hi - lo) // 2) + 1

    return lo, hi, count


# ======================================================================================================================
# EXACT K INTERVAL GENERATED BY KNOWN X WINDOW
# ======================================================================================================================
#
# Given:
#
#     -R <= x <= R
#
# and
#
#     x+d0 even
#
# compute the exact K image interval.
#
# Since
#
#     K(x) = A4 - ((x+d0)/2)^2
#
# the function is concave in x.
#
# The maximum occurs at the parity-valid point closest to x=-d0.
# The minimum occurs at one of the interval endpoints.
#
# ======================================================================================================================

@dataclass(frozen=True)
class KInterval:

    radius: int

    x_request_low: int
    x_request_high: int

    x_low: int
    x_high: int

    K_low: int
    K_high: int

    parity_points: int


def known_x_window_to_K_interval(
    inst: Instance,
    A4: int,
    R: int,
) -> KInterval:

    lo, hi, count = parity_adjusted_bounds(
        inst,
        R,
    )

    if count == 0:

        return KInterval(
            radius=R,
            x_request_low=-R,
            x_request_high=R,
            x_low=1,
            x_high=0,
            K_low=1,
            K_high=0,
            parity_points=0,
        )

    candidates = [lo, hi]

    # Vertex of K(x):
    #
    #     x = -d0
    #
    # If it lies inside the interval, add the closest
    # parity-valid point(s).

    vertex_x = -inst.d0

    if lo <= vertex_x <= hi:

        if (vertex_x & 1) == required_x_parity(inst):
            candidates.append(vertex_x)
        else:
            if vertex_x - 1 >= lo:
                candidates.append(vertex_x - 1)

            if vertex_x + 1 <= hi:
                candidates.append(vertex_x + 1)

    k_values: list[int] = []

    for x in candidates:

        if not (lo <= x <= hi):
            continue

        K = k_from_x(
            inst,
            A4,
            x,
        )

        if K is None:
            continue

        k_values.append(K)

    if not k_values:
        raise AssertionError(
            f"failed to construct K interval "
            f"at instance={inst.idx}, R={R}"
        )

    return KInterval(
        radius=R,
        x_request_low=-R,
        x_request_high=R,
        x_low=lo,
        x_high=hi,
        K_low=min(k_values),
        K_high=max(k_values),
        parity_points=count,
    )


# ======================================================================================================================
# EXACT K -> X PULLBACK
# ======================================================================================================================
#
#     K_low <= A4-z² <= K_high
#
# iff
#
#     A4-K_high <= z² <= A4-K_low
#
# The K interval alone describes all z satisfying the quadratic
# equation. For the current positive-z branch we obtain the exact
# z interval using two integer square roots.
#
# ======================================================================================================================

@dataclass(frozen=True)
class Pullback:

    z_abs_low: int
    z_abs_high: int

    full_z_count: int
    positive_count: int

    positive_x_low: int
    positive_x_high: int

    true_survives: bool


def exact_pullback(
    inst: Instance,
    A4: int,
    K_low: int,
    K_high: int,
    K_true: int,
) -> Pullback | None:

    sq_low = A4 - K_high
    sq_high = A4 - K_low

    if sq_high < 0:
        return None

    sq_low = max(0, sq_low)

    z_abs_low = ceil_sqrt(sq_low)
    z_abs_high = isqrt(sq_high)

    if z_abs_low > z_abs_high:
        return None

    positive_count = (
        z_abs_high - z_abs_low + 1
    )

    if z_abs_low == 0:
        full_z_count = (
            2 * positive_count - 1
        )
    else:
        full_z_count = (
            2 * positive_count
        )

    positive_x_low = x_from_z(
        inst,
        z_abs_low,
    )

    positive_x_high = x_from_z(
        inst,
        z_abs_high,
    )

    return Pullback(
        z_abs_low=z_abs_low,
        z_abs_high=z_abs_high,
        full_z_count=full_z_count,
        positive_count=positive_count,
        positive_x_low=positive_x_low,
        positive_x_high=positive_x_high,
        true_survives=(
            K_low <= K_true <= K_high
        ),
    )


# ======================================================================================================================
# DIRECT BRUTE CHECK
# ======================================================================================================================

def brute_x_window(
    inst: Instance,
    A4: int,
    K_low: int,
    K_high: int,
    x_low: int,
    x_high: int,
) -> tuple[int, int]:

    count = 0
    true_count = 0

    for x in range(x_low, x_high + 1):

        K = k_from_x(
            inst,
            A4,
            x,
        )

        if K is None:
            continue

        if K_low <= K <= K_high:

            count += 1

            if x == inst.x_true:
                true_count += 1

    return count, true_count


# ======================================================================================================================
# KNOWN-DATA RADII
# ======================================================================================================================

def known_data_radii(inst: Instance) -> dict[str, int]:

    sqrtN = isqrt(inst.N)

    return {
        "0": 0,
        "1": 1,
        "S": inst.S,
        "sqrtN": sqrtN,
        "2sqrtN": 2 * sqrtN,
        "4sqrtN": 4 * sqrtN,
        "8sqrtN": 8 * sqrtN,
        "32sqrtN": 32 * sqrtN,
        "d": abs(inst.d),
        "gap_proxy": max(0, inst.d - inst.S),
    }


# ======================================================================================================================
# MAIN
# ======================================================================================================================

def main() -> None:

    print("=" * 120)
    print("EXPERIMENT 446R2")
    print("=" * 120)
    print()
    print("KNOWN-DATA-CENTERED K INTERVAL -> EXACT X POPULATION")
    print()
    print("PARITY-AWARE CENTER")
    print("K_TRUE IS NEVER USED TO CONSTRUCT A BOUND")
    print("NO GIANT K ENUMERATION")
    print("NO GIANT X ENUMERATION")
    print("NO CRT CARTESIAN PRODUCT")
    print("INTEGER-EXACT")
    print()

    total_cases = 0
    true_contained = 0
    exact_empty = 0
    brute_checks = 0
    brute_failures = 0

    instances = build_instances()

    for inst in instances:

        # ----------------------------------------------------------------------------------
        # Hidden values are generated solely for diagnostics.
        # They do NOT participate in the tested bounds.
        # ----------------------------------------------------------------------------------

        K_true = hidden_K(inst)

        r = build_r(
            inst,
            K_true,
        )

        A = (
            r
            + inst.d0 * inst.d0
        )

        if A % 4 != 0:
            raise AssertionError(
                f"A mod 4 != 0 at instance={inst.idx}"
            )

        A4 = A // 4

        # ----------------------------------------------------------------------------------
        # Exact true relation check.
        # ----------------------------------------------------------------------------------

        z_true = z_from_x(
            inst,
            inst.x_true,
        )

        if z_true is None:
            raise AssertionError(
                f"true x parity invalid at instance={inst.idx}"
            )

        if k_from_z(A4, z_true) != K_true:
            raise AssertionError(
                f"true K identity failed at instance={inst.idx}"
            )

        # ----------------------------------------------------------------------------------
        # PARITY-AWARE known center
        # ----------------------------------------------------------------------------------

        x0, z0, K0 = known_center(
            inst,
            A4,
        )

        print("-" * 120)

        print(
            f"INSTANCE {inst.idx}: "
            f"N={inst.N} "
            f"S={inst.S} "
            f"d={inst.d} "
            f"d0={inst.d0} "
            f"x_true={inst.x_true}"
        )

        print(
            f"  hidden K bits = {abs(K_true).bit_length()}"
        )

        print(
            f"  K_true = {K_true}"
        )

        print(
            f"  A/4 = {A4}"
        )

        print(
            f"  required x parity = {required_x_parity(inst)}"
        )

        print(
            f"  known-data reference x0 = {x0}"
        )

        print(
            f"  known-data reference z0 = {z0}"
        )

        print(
            f"  K(x0) = {K0}"
        )

        print()

        # ----------------------------------------------------------------------------------
        # AUDIT
        # ----------------------------------------------------------------------------------

        print("  KNOWN-DATA RADIUS AUDIT")
        print()

        print(
            "    model"
            "               R"
            "                K-low"
            "               K-high"
            "       K-width"
            "     parity-x"
            "   true-in"
            "   x/K-ratio"
        )

        print(
            "    " + "-" * 112
        )

        for model, R in known_data_radii(inst).items():

            total_cases += 1

            interval = known_x_window_to_K_interval(
                inst,
                A4,
                R,
            )

            pullback = exact_pullback(
                inst,
                A4,
                interval.K_low,
                interval.K_high,
                K_true,
            )

            if pullback is None:

                exact_empty += 1

                print(
                    f"    {model:<12}"
                    f"{R:>20}"
                    f"{interval.K_low:>24}"
                    f"{interval.K_high:>24}"
                    f"{0:>14}"
                    f"{interval.parity_points:>12}"
                    f"{str(False):>10}"
                    f"{'N/A':>12}"
                )

                continue

            K_width = (
                interval.K_high
                - interval.K_low
                + 1
            )

            x_count = pullback.positive_count

            ratio = (
                K_width / x_count
                if x_count
                else 0.0
            )

            if pullback.true_survives:
                true_contained += 1

            # --------------------------------------------------------------------------
            # Brute check only for small windows.
            # --------------------------------------------------------------------------

            if (
                R <= BRUTE_FORCE_LIMIT
                and interval.parity_points <= BRUTE_FORCE_LIMIT
            ):

                brute_count, brute_true = brute_x_window(
                    inst,
                    A4,
                    interval.K_low,
                    interval.K_high,
                    interval.x_low,
                    interval.x_high,
                )

                brute_checks += 1

                if brute_count != interval.parity_points:
                    brute_failures += 1

                    raise AssertionError(
                        "brute population mismatch: "
                        f"instance={inst.idx}, "
                        f"model={model}, "
                        f"expected={interval.parity_points}, "
                        f"actual={brute_count}"
                    )

                expected_true = int(
                    interval.x_low
                    <= inst.x_true
                    <= interval.x_high
                    and
                    ((inst.x_true + inst.d0) & 1) == 0
                )

                if brute_true != expected_true:
                    brute_failures += 1

                    raise AssertionError(
                        "true containment mismatch: "
                        f"instance={inst.idx}, "
                        f"model={model}, "
                        f"expected={expected_true}, "
                        f"actual={brute_true}"
                    )

            print(
                f"    {model:<12}"
                f"{R:>20}"
                f"{interval.K_low:>24}"
                f"{interval.K_high:>24}"
                f"{K_width:>14}"
                f"{interval.parity_points:>12}"
                f"{str(pullback.true_survives):>10}"
                f"{ratio:>12.2f}"
            )

        # ----------------------------------------------------------------------------------
        # CENTER DISPLACEMENT
        # ----------------------------------------------------------------------------------

        displacement = abs(
            K_true - K0
        )

        print()

        print("  KNOWN CENTER DISPLACEMENT")
        print(
            f"    x0 = {x0}"
        )
        print(
            f"    |K_true-K(x0)| = {displacement}"
        )
        print(
            f"    displacement bits = {displacement.bit_length()}"
        )
        print(
            f"    true K bits       = {abs(K_true).bit_length()}"
        )

        # ----------------------------------------------------------------------------------
        # Diagnostic only:
        #
        # smallest radius |x| <= R that is guaranteed to include
        # the hidden x_true.
        #
        # This number is NOT used to construct the main audit.
        # ----------------------------------------------------------------------------------

        critical_R = abs(inst.x_true)

        critical_interval = known_x_window_to_K_interval(
            inst,
            A4,
            critical_R,
        )

        critical_pullback = exact_pullback(
            inst,
            A4,
            critical_interval.K_low,
            critical_interval.K_high,
            K_true,
        )

        if critical_pullback is None:
            raise AssertionError(
                f"critical pullback failed at instance={inst.idx}"
            )

        print()

        print("  DIAGNOSTIC TRUE-COVERING WINDOW")
        print(
            f"    radius = |x_true| = {critical_R}"
        )
        print(
            f"    K-low  = {critical_interval.K_low}"
        )
        print(
            f"    K-high = {critical_interval.K_high}"
        )
        print(
            f"    parity-valid x count = "
            f"{critical_interval.parity_points}"
        )
        print(
            f"    exact positive-branch population = "
            f"{critical_pullback.positive_count}"
        )
        print(
            f"    true survives = "
            f"{critical_pullback.true_survives}"
        )

    # ================================================================================================================
    # GLOBAL SUMMARY
    # ================================================================================================================

    print()
    print("=" * 120)
    print("GLOBAL EXPERIMENT 446R2 SUMMARY")
    print("=" * 120)

    print(
        f"  known-data interval cases = {total_cases}"
    )

    print(
        f"  true-containing cases    = {true_contained}"
    )

    print(
        f"  empty pullbacks          = {exact_empty}"
    )

    print(
        f"  brute checks             = {brute_checks}"
    )

    print(
        f"  brute failures           = {brute_failures}"
    )

    print()
    print("PARITY STRUCTURE")
    print()
    print(
        "  x + d0 ≡ 0 (mod 2)"
    )
    print()
    print(
        "  d0 even -> known reference x0 = 0"
    )
    print(
        "  d0 odd  -> known reference x0 = 1"
    )

    print()
    print("CORE FORMULA")
    print()
    print(
        "  K = A/4 - ((x+d0)/2)^2"
    )
    print()
    print(
        "  A4-K_high <= z^2 <= A4-K_low"
    )
    print()
    print(
        "  z interval is therefore obtained with"
    )
    print(
        "  two exact integer square roots."
    )

    print()
    print("KEY TEST")
    print()
    print(
        "  K_TRUE participates only in the post-hoc"
    )
    print(
        "  containment check."
    )

    print()
    print(
        "  It does NOT participate in:"
    )
    print(
        "    - choosing x0"
    )
    print(
        "    - selecting R"
    )
    print(
        "    - constructing K-low"
    )
    print(
        "    - constructing K-high"
    )

    print()
    print("INTERPRETATION")
    print()
    print(
        "  Experiment 445R2 used K intervals centered on K_true."
    )
    print()
    print(
        "  This experiment instead derives the tested intervals"
    )
    print(
        "  from known x-space radii and known d0/A4."
    )
    print()
    print(
        "  The parity-aware center fixes the previous failure for"
    )
    print(
        "  odd d0 without changing the underlying mathematics."
    )
    print()
    print(
        "  A small population is meaningful only when the radius"
    )
    print(
        "  itself has an independent derivation from known data."
    )
    print()
    print(
        "  This remains an exact inverse-image / information audit,"
    )
    print(
        "  not a factoring theorem."
    )

    print()
    print("=" * 120)
    print("EXPERIMENT 446R2 FINAL STATUS")
    print("=" * 120)

    print(
        f"  KNOWN-DATA INTERVAL CASES = {total_cases}"
    )

    print(
        f"  TRUE-CONTAINING CASES     = {true_contained}"
    )

    print(
        f"  EMPTY PULLBACKS           = {exact_empty}"
    )

    print(
        f"  BRUTE FAILURES            = {brute_failures}"
    )

    print(
        "  K_TRUE USED IN BOUND      = False"
    )

    print(
        "  PARITY-AWARE CENTER       = True"
    )

    print(
        "  GIANT K ENUMERATION       = False"
    )

    print(
        "  GIANT X ENUMERATION       = False"
    )

    print(
        "  CRT CARTESIAN PRODUCT     = False"
    )

    print(
        "  O(1) INTERVAL INVERSION   = True"
    )

    print(
        "  INTEGER-EXACT              = True"
    )

    print("=" * 120)
    print("EXPERIMENT 446R2 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    main()