#!/usr/bin/env python3

from __future__ import annotations

from dataclasses import dataclass
from math import isqrt


# ======================================================================================================================
# EXPERIMENT 445
# ======================================================================================================================
#
# EXACT K-INTERVAL -> X-COUNT / MONOTONE PULLBACK AUDIT
#
# QUESTION
#   Once the CRT branch is known to be tautological after substituting
#   K = A/4 - ((x+d0)/2)^2, how much does a finite K interval actually
#   reduce the integer x population?
#
# OBJECTIVE
#   Replace K enumeration entirely by an exact inverse-image calculation.
#
#   For
#
#       z = (x+d0)/2
#       K = A4 - z^2
#
#   and an exact integer interval
#
#       K_low <= K <= K_high,
#
#   determine the exact admissible z and x interval.
#
# KEY INEQUALITY
#
#       K_low <= A4 - z^2 <= K_high
#
#   is equivalent to
#
#       A4 - K_high <= z^2 <= A4 - K_low.
#
#   This is an exact interval of squares.
#
# RULES
#   exact integer arithmetic only
#   no K enumeration
#   no giant x enumeration
#   no CRT Cartesian product
#   no factor search
#   no symbolic algebra systems
#
# TESTS
#   1. exact inverse K interval -> z interval
#   2. exact z interval -> x interval
#   3. parity consistency
#   4. monotonicity / uniqueness on each branch
#   5. direct point checks around interval boundaries
#   6. compare exact candidate counts with brute force only for tiny
#      intervals as an implementation cross-check
#
# IMPORTANT
#   The K intervals themselves are measurement models.
#   This experiment does NOT claim that any tested K radius is justified.
#
# ======================================================================================================================


K_RADII = [
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

# These are deliberate diagnostic radii in x-space.
X_CHECK_RADII = [
    1,
    2,
    4,
    8,
    16,
    32,
    64,
]

# Small intervals can be brute-force checked without risk.
BRUTE_FORCE_LIMIT = 200000


# ======================================================================================================================
# INSTANCE DATA
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


INSTANCES_RAW = [
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
    result: list[Instance] = []

    for idx, (p, q, x_true) in enumerate(
        INSTANCES_RAW,
        start=1,
    ):
        N = p * q
        S = p + q
        d = N - 2 * S + 1
        d0 = d - x_true

        result.append(
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

    return result


# ======================================================================================================================
# SAME CONTROL GENERATOR USED BY THE PREVIOUS EXPERIMENTS
# ======================================================================================================================

HIDDEN_FACTORS = [17, 29, 43, 61, 89, 127, 191, 257]


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

    r = (
        4 * K_true
        + x * x
        + 2 * inst.d0 * x
    )

    numerator = r - x * x - 2 * inst.d0 * x

    assert numerator % 4 == 0
    assert numerator // 4 == K_true

    return r


# ======================================================================================================================
# K / X RELATION
# ======================================================================================================================

def z_from_x(inst: Instance, x: int) -> int | None:
    y = x + inst.d0

    if y & 1:
        return None

    return y // 2


def k_from_x(
    inst: Instance,
    A4: int,
    x: int,
) -> int | None:

    z = z_from_x(inst, x)

    if z is None:
        return None

    return A4 - z * z


def x_from_z(
    inst: Instance,
    z: int,
) -> int:

    return 2 * z - inst.d0


# ======================================================================================================================
# EXACT INTEGER SQUARE ROOT HELPERS
# ======================================================================================================================

def ceil_sqrt(n: int) -> int:
    """
    Smallest t >= 0 such that t^2 >= n.

    Requires n >= 0.
    """

    if n < 0:
        raise ValueError("ceil_sqrt requires n >= 0")

    r = isqrt(n)

    if r * r == n:
        return r

    return r + 1


def floor_sqrt(n: int) -> int:
    """
    Largest t >= 0 such that t^2 <= n.

    Requires n >= 0.
    """

    if n < 0:
        raise ValueError("floor_sqrt requires n >= 0")

    return isqrt(n)


# ======================================================================================================================
# EXACT K-INTERVAL PULLBACK
# ======================================================================================================================

@dataclass(frozen=True)
class Pullback:
    k_low: int
    k_high: int

    q_low: int
    q_high: int

    z_abs_low: int
    z_abs_high: int

    z_low: int
    z_high: int

    x_low: int
    x_high: int

    x_count: int
    parity_valid_count: int

    contains_true: bool
    true_z: int
    true_k: int


def pullback_k_interval(
    inst: Instance,
    A4: int,
    k_low: int,
    k_high: int,
) -> Pullback | None:
    """
    Pull back the exact K interval through

        K = A4 - z^2.

    We require

        k_low <= k_high.

    Then

        A4-k_high <= z^2 <= A4-k_low.
    """

    if k_low > k_high:
        raise ValueError("invalid K interval")

    q_low = A4 - k_high
    q_high = A4 - k_low

    # z^2 is never negative.
    if q_high < 0:
        return None

    square_low = max(0, q_low)
    square_high = q_high

    z_abs_low = ceil_sqrt(square_low)
    z_abs_high = floor_sqrt(square_high)

    if z_abs_low > z_abs_high:
        return None

    # There are two sign branches:
    #
    #   -z_abs_high <= z <= -z_abs_low
    #       OR
    #    z_abs_low <= z <= z_abs_high
    #
    # This gives the complete z preimage.
    #
    # The corresponding x values are obtained by x=2z-d0.
    #
    # Because the two branches are disjoint except when zero occurs,
    # count them exactly.

    negative_count = z_abs_high - z_abs_low + 1
    positive_count = z_abs_high - z_abs_low + 1

    if z_abs_low == 0:
        # Zero is shared by the two sign branches.
        z_count = (
            negative_count
            + positive_count
            - 1
        )
    else:
        z_count = negative_count + positive_count

    z_low = -z_abs_high
    z_high = z_abs_high

    x_low = x_from_z(inst, z_low)
    x_high = x_from_z(inst, z_high)

    true_z = z_from_x(inst, inst.x_true)

    if true_z is None:
        raise AssertionError(
            f"true x has invalid parity at instance={inst.idx}"
        )

    true_k = A4 - true_z * true_z

    contains_true = (
        k_low <= true_k <= k_high
    )

    return Pullback(
        k_low=k_low,
        k_high=k_high,
        q_low=q_low,
        q_high=q_high,
        z_abs_low=z_abs_low,
        z_abs_high=z_abs_high,
        z_low=z_low,
        z_high=z_high,
        x_low=x_low,
        x_high=x_high,
        x_count=z_count,
        parity_valid_count=z_count,
        contains_true=contains_true,
        true_z=true_z,
        true_k=true_k,
    )


# ======================================================================================================================
# DIRECT BRUTE FORCE CROSS-CHECK
# ======================================================================================================================

def brute_force_x_count(
    inst: Instance,
    A4: int,
    k_low: int,
    k_high: int,
) -> tuple[int, bool]:
    """
    Only used when the exact pullback predicts a small population.

    Returns:
        count,
        true_survives
    """

    values: set[int] = set()

    radius_limit = 0

    # We do not know a priori how large the x interval is.
    # Instead obtain the exact pullback first outside this function
    # and let callers only invoke this when the interval is small.
    #
    # This function is intentionally implemented as a direct scan
    # over the exact x interval supplied separately elsewhere.
    raise RuntimeError("use brute_force_exact_interval")


def brute_force_exact_interval(
    inst: Instance,
    A4: int,
    pullback: Pullback,
) -> tuple[int, bool]:

    population = pullback.x_count

    if population > BRUTE_FORCE_LIMIT:
        raise MemoryError(
            f"brute-force cross-check too large: "
            f"{population} > {BRUTE_FORCE_LIMIT}"
        )

    survivors = 0
    true_survives = False

    # x = 2z-d0, so enumerate z directly.
    for z in range(
        pullback.z_low,
        pullback.z_high + 1,
    ):
        square = z * z

        if not (
            pullback.q_low
            <= square
            <= pullback.q_high
        ):
            continue

        x = x_from_z(inst, z)

        # Exact definition of K.
        K = A4 - square

        if not (
            pullback.k_low
            <= K
            <= pullback.k_high
        ):
            raise AssertionError(
                "inverse interval mismatch"
            )

        survivors += 1

        if x == inst.x_true:
            true_survives = True

    return survivors, true_survives


# ======================================================================================================================
# WIDTH / LOCAL SENSITIVITY
# ======================================================================================================================

def x_branch_distance(
    inst: Instance,
    x: int,
) -> int:
    return abs(x - inst.x_true)


def k_slope_near_true(
    inst: Instance,
    A4: int,
    x: int,
) -> int:
    """
    Exact one-step change for an allowed parity step x -> x+2.

        z -> z+1

        K(x+2)-K(x)
          = -(z+1)^2 + z^2
          = -2z-1
          = -(x+d0)-1.
    """

    z = z_from_x(inst, x)

    if z is None:
        raise ValueError(
            "slope requested at parity-invalid x"
        )

    return -(2 * z + 1)


def exact_k_difference(
    inst: Instance,
    A4: int,
    x: int,
) -> int:
    k0 = k_from_x(inst, A4, x)
    k1 = k_from_x(inst, A4, x + 2)

    if k0 is None or k1 is None:
        raise AssertionError(
            "unexpected parity mismatch"
        )

    return k1 - k0


# ======================================================================================================================
# INTERVAL MODELS
# ======================================================================================================================

def make_k_interval(
    k_true: int,
    radius: int,
) -> tuple[int, int]:

    return (
        k_true - radius,
        k_true + radius,
    )


def x_interval_width(
    pullback: Pullback,
) -> int:
    if pullback.x_count == 0:
        return 0

    return pullback.x_high - pullback.x_low + 1


# ======================================================================================================================
# MAIN
# ======================================================================================================================

def main() -> None:

    print("=" * 120)
    print("EXPERIMENT 445")
    print("=" * 120)
    print()
    print("EXACT K-INTERVAL -> X-COUNT / MONOTONE PULLBACK AUDIT")
    print()
    print("QUESTION")
    print("  How does a finite K interval translate into an exact")
    print("  integer x population?")
    print()
    print("K_RADII =", K_RADII)
    print("BRUTE_FORCE_LIMIT =", BRUTE_FORCE_LIMIT)
    print()

    total_interval_tests = 0
    total_empty_pullbacks = 0
    total_true_containment = 0
    total_monotonicity_failures = 0
    total_bruteforce_mismatches = 0

    instances = build_instances()

    for inst in instances:

        K_true = hidden_K(inst)
        r = build_r(inst, K_true)

        A = r + inst.d0 * inst.d0

        if A % 4 != 0:
            raise AssertionError(
                f"A % 4 != 0 at instance={inst.idx}"
            )

        A4 = A // 4

        true_z = z_from_x(inst, inst.x_true)

        if true_z is None:
            raise AssertionError(
                f"true x is not parity-valid at instance={inst.idx}"
            )

        reconstructed_K = A4 - true_z * true_z

        if reconstructed_K != K_true:
            raise AssertionError(
                f"K identity mismatch at instance={inst.idx}"
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
            f"  true z = {true_z}"
        )
        print()

        # ==================================================================
        # TRUE-BRANCH LOCAL SLOPE
        # ==================================================================

        delta_k = exact_k_difference(
            inst,
            A4,
            inst.x_true,
        )

        predicted_delta = k_slope_near_true(
            inst,
            A4,
            inst.x_true,
        )

        if delta_k != predicted_delta:
            raise AssertionError(
                f"slope mismatch at instance={inst.idx}"
            )

        print("  LOCAL K-SPACING AROUND TRUE x")
        print(
            "    K(x+2)-K(x) exact       = "
            f"{delta_k}"
        )
        print(
            "    -(x+d0)-1 predicted     = "
            f"{predicted_delta}"
        )
        print(
            "    exact spacing identity  = True"
        )
        print()

        # ==================================================================
        # FINITE K INTERVAL AUDIT
        # ==================================================================

        print("  BOUNDED-K INVERSE IMAGE")
        print()
        print(
            "    radius"
            "       K-width"
            "     pullback-x-count"
            "       x-span"
            " true-in-bound"
            "   brute-check"
        )
        print("    " + "-" * 92)

        for radius in K_RADII:

            k_low, k_high = make_k_interval(
                K_true,
                radius,
            )

            pullback = pullback_k_interval(
                inst,
                A4,
                k_low,
                k_high,
            )

            total_interval_tests += 1

            if pullback is None:
                total_empty_pullbacks += 1

                print(
                    f"    {radius:>8}"
                    f"{2 * radius + 1:>15}"
                    f"{0:>20}"
                    f"{0:>14}"
                    f"{str(False):>15}"
                    f"{'NONE':>14}"
                )

                continue

            if pullback.contains_true:
                total_true_containment += 1

            brute_status = "SKIP"

            if pullback.x_count <= BRUTE_FORCE_LIMIT:

                brute_count, brute_true = (
                    brute_force_exact_interval(
                        inst,
                        A4,
                        pullback,
                    )
                )

                brute_status = (
                    f"{brute_count}/{str(brute_true)}"
                )

                if brute_count != pullback.x_count:
                    total_bruteforce_mismatches += 1

                    raise AssertionError(
                        f"brute-force count mismatch "
                        f"instance={inst.idx}, "
                        f"radius={radius}"
                    )

                if brute_true != pullback.contains_true:
                    total_bruteforce_mismatches += 1

                    raise AssertionError(
                        f"brute-force true mismatch "
                        f"instance={inst.idx}, "
                        f"radius={radius}"
                    )

            print(
                f"    {radius:>8}"
                f"{2 * radius + 1:>15}"
                f"{pullback.x_count:>20}"
                f"{x_interval_width(pullback):>14}"
                f"{str(pullback.contains_true):>15}"
                f"{brute_status:>14}"
            )

        print()

        # ==================================================================
        # TRUE-CENTERED X WINDOWS -> EXACT K RANGE
        # ==================================================================

        print("  TRUE-CENTERED X-WINDOW AUDIT")
        print()
        print(
            "    x-radius"
            "    admissible-x"
            "     K-min"
            "     K-max"
            "       K-width"
            "      monotone"
        )
        print("    " + "-" * 95)

        for xr in X_CHECK_RADII:

            xs = []

            # Only parity-valid x belong to the exact K domain.
            for x in range(
                inst.x_true - xr,
                inst.x_true + xr + 1,
            ):
                if z_from_x(inst, x) is not None:
                    xs.append(x)

            if not xs:
                print(
                    f"    {xr:>8}"
                    f"{0:>16}"
                    f"{'N/A':>16}"
                    f"{'N/A':>16}"
                    f"{0:>16}"
                    f"{False:>12}"
                )
                continue

            k_values = [
                k_from_x(inst, A4, x)
                for x in xs
            ]

            if any(k is None for k in k_values):
                raise AssertionError(
                    "unexpected None in x-window"
                )

            # Narrow the type for static reasoning.
            k_values_int = [int(k) for k in k_values]

            k_min = min(k_values_int)
            k_max = max(k_values_int)

            # Since x_true is far to the right of the parabola vertex,
            # K must decrease as x increases.
            monotone_decreasing = all(
                k_values_int[i + 1] < k_values_int[i]
                for i in range(len(k_values_int) - 1)
            )

            if not monotone_decreasing and len(k_values_int) > 1:
                total_monotonicity_failures += 1

                raise AssertionError(
                    f"monotonicity failure "
                    f"instance={inst.idx}, "
                    f"x-radius={xr}"
                )

            print(
                f"    {xr:>8}"
                f"{len(xs):>16}"
                f"{k_min:>16}"
                f"{k_max:>16}"
                f"{k_max - k_min:>16}"
                f"{monotone_decreasing:>12}"
            )

        print()

        # ==================================================================
        # EXACT BOUNDARY CHECKS
        # ==================================================================

        print("  BOUNDARY INVERSION CHECK")

        boundary_radii = [
            0,
            1,
            2,
            4,
            8,
            16,
            32,
            64,
        ]

        boundary_failures = 0

        for xr in boundary_radii:

            for direction in (-1, 1):

                x = inst.x_true + direction * xr

                z = z_from_x(inst, x)

                if z is None:
                    continue

                Kx = k_from_x(
                    inst,
                    A4,
                    x,
                )

                if Kx is None:
                    raise AssertionError(
                        "Kx unexpectedly None"
                    )

                reverse_z_sq = A4 - Kx

                if reverse_z_sq != z * z:
                    boundary_failures += 1

                    raise AssertionError(
                        f"boundary inversion mismatch "
                        f"instance={inst.idx}, "
                        f"x={x}"
                    )

        print(
            f"    tested x-boundaries        = "
            f"{len(boundary_radii) * 2}"
        )
        print(
            f"    exact inversion failures   = "
            f"{boundary_failures}"
        )
        print(
            "    boundary inversion = TRUE"
        )
        print()

    # ======================================================================
    # GLOBAL SUMMARY
    # ======================================================================

    print("=" * 120)
    print("GLOBAL EXACT K->X PULLBACK SUMMARY")
    print("=" * 120)

    print(
        f"  total K interval tests          = "
        f"{total_interval_tests}"
    )

    print(
        f"  empty inverse images            = "
        f"{total_empty_pullbacks}"
    )

    print(
        f"  true-containing K intervals     = "
        f"{total_true_containment}"
    )

    print(
        f"  monotonicity failures           = "
        f"{total_monotonicity_failures}"
    )

    print(
        f"  brute-force cross-check errors  = "
        f"{total_bruteforce_mismatches}"
    )

    print()
    print("EXACT STRUCTURE")
    print()
    print(
        "  K = A4 - z^2"
    )
    print(
        "  x = 2z - d0"
    )
    print()
    print(
        "  Therefore a finite K interval is pulled back by"
    )
    print(
        "  solving an exact interval-of-squares problem."
    )
    print()
    print(
        "  No K enumeration is required."
    )
    print()
    print(
        "  On the tested true branch, the exact step relation is:"
    )
    print()
    print(
        "      K(x+2)-K(x) = -(x+d0)-1."
    )
    print()
    print(
        "  Hence sufficiently close to the true x values, where"
    )
    print(
        "  x+d0 > 0, K is strictly decreasing as x increases."
    )
    print()
    print(
        "  The remaining candidate count is therefore controlled"
    )
    print(
        "  directly by the width and placement of the independent K bound."
    )
    print()
    print("INTERPRETATION")
    print()
    print(
        "  Experiment 444 showed that the local quadratic-residue"
    )
    print(
        "  conditions become identities after substituting K(x)."
    )
    print()
    print(
        "  Experiment 445 removes the remaining computational"
    )
    print(
        "  artifact: there is no need to enumerate K candidates"
    )
    print(
        "  to determine how restrictive a K interval is in x-space."
    )
    print()
    print(
        "  A K interval can be inverted exactly into one or two"
    )
    print(
        "  z branches, hence into an exact x population."
    )
    print()
    print(
        "  This separates the problem cleanly:"
    )
    print()
    print(
        "      independent K information"
    )
    print(
        "              -> exact x interval"
    )
    print(
        "              -> exact candidate count."
    )
    print()
    print(
        "  The hard question is consequently no longer the CRT"
    )
    print(
        "  pullback. It is whether an independently justified"
    )
    print(
        "  K bound can be made narrow enough to force a tiny"
    )
    print(
        "  x population."
    )
    print()
    print(
        "  The radii tested here are diagnostic intervals centered"
    )
    print(
        "  at the hidden K and are not presented as independently"
    )
    print(
        "  justified factorization bounds."
    )
    print()
    print(
        "  This is an exact inverse-image / information-width audit,"
    )
    print(
        "  not a factoring theorem."
    )

    print()
    print("=" * 120)
    print("EXPERIMENT 445 FINAL STATUS")
    print("=" * 120)
    print(
        f"  K INTERVAL TESTS          = "
        f"{total_interval_tests}"
    )
    print(
        f"  EMPTY PULLBACKS           = "
        f"{total_empty_pullbacks}"
    )
    print(
        f"  TRUE-CONTAINING INTERVALS = "
        f"{total_true_containment}"
    )
    print(
        f"  MONOTONICITY FAILURES     = "
        f"{total_monotonicity_failures}"
    )
    print(
        f"  BRUTE-FORCE MISMATCHES    = "
        f"{total_bruteforce_mismatches}"
    )
    print("  GIANT K ENUMERATION = False")
    print("  GIANT X ENUMERATION = False")
    print("  CRT CARTESIAN PRODUCT = False")
    print("  EXACT INTEGER ARITHMETIC = True")
    print("=" * 120)
    print("EXPERIMENT 445 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    main()

