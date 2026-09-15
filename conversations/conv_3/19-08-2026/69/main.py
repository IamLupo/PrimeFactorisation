#!/usr/bin/env python3

from __future__ import annotations

from dataclasses import dataclass
from math import isqrt


# ======================================================================================================================
# EXPERIMENT 445R2
# ======================================================================================================================
#
# O(1) EXACT K-INTERVAL -> X-COUNT PULLBACK
#
# Core identity:
#
#   y = x + d0
#   A = r + d0^2
#   4K = A - y^2
#
# Since A is divisible by 4:
#
#   y = 2z
#   K = A4 - z^2
#
# where
#
#   A4 = A/4
#   z  = (x+d0)/2
#
# For a finite K interval
#
#   K_low <= K <= K_high
#
# we require
#
#   A4-K_high <= z^2 <= A4-K_low.
#
# The exact integer population is therefore determined entirely
# by square-root boundaries.
#
# NO ENUMERATION OF K
# NO ENUMERATION OF LARGE X RANGES
# NO CRT PRODUCT
#
# Everything below the analytical boundary calculation is only
# a tiny brute-force verification on very small populations.
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
    268435456,
]

# Only use brute force when the exact answer is this small.
BRUTE_FORCE_LIMIT = 10000


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
    out = []

    for idx, (p, q, x_true) in enumerate(RAW, start=1):
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
                d0=d0,
                x_true=x_true,
            )
        )

    return out


# ======================================================================================================================
# CONTROLLED HIDDEN K GENERATOR
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

    return r


# ======================================================================================================================
# CORE EXACT MAP
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


def k_from_x(inst: Instance, A4: int, x: int) -> int | None:
    z = z_from_x(inst, x)

    if z is None:
        return None

    return k_from_z(A4, z)


# ======================================================================================================================
# EXACT SQRT HELPERS
# ======================================================================================================================

def ceil_sqrt(n: int) -> int:
    if n < 0:
        raise ValueError("ceil_sqrt requires n >= 0")

    r = isqrt(n)

    if r * r == n:
        return r

    return r + 1


# ======================================================================================================================
# EXACT PREIMAGE COUNTS
# ======================================================================================================================

@dataclass(frozen=True)
class PullbackResult:
    k_low: int
    k_high: int

    sq_low: int
    sq_high: int

    z_abs_low: int
    z_abs_high: int

    z_count: int
    x_count: int

    z_low: int
    z_high: int

    x_low: int
    x_high: int

    true_in: bool


def exact_pullback(
    inst: Instance,
    A4: int,
    k_low: int,
    k_high: int,
    true_k: int,
) -> PullbackResult | None:

    if k_low > k_high:
        raise ValueError("K interval invalid")

    # K_low <= A4-z² <= K_high
    #
    # A4-K_high <= z² <= A4-K_low

    sq_low = A4 - k_high
    sq_high = A4 - k_low

    # z² >= 0
    if sq_high < 0:
        return None

    sq_low = max(0, sq_low)

    z_abs_low = ceil_sqrt(sq_low)
    z_abs_high = isqrt(sq_high)

    if z_abs_low > z_abs_high:
        return None

    # Number of NONNEGATIVE z values.
    nonnegative_count = z_abs_high - z_abs_low + 1

    # Reflect negative branch.
    #
    # If 0 is included, it occurs only once globally.
    if z_abs_low == 0:
        z_count = 2 * nonnegative_count - 1
    else:
        z_count = 2 * nonnegative_count

    z_low = -z_abs_high
    z_high = z_abs_high

    x_low = x_from_z(inst, z_low)
    x_high = x_from_z(inst, z_high)

    return PullbackResult(
        k_low=k_low,
        k_high=k_high,
        sq_low=sq_low,
        sq_high=sq_high,
        z_abs_low=z_abs_low,
        z_abs_high=z_abs_high,
        z_count=z_count,
        x_count=z_count,
        z_low=z_low,
        z_high=z_high,
        x_low=x_low,
        x_high=x_high,
        true_in=(k_low <= true_k <= k_high),
    )


# ======================================================================================================================
# LOCAL / TRUE-BRANCH PULLBACK
# ======================================================================================================================
#
# The full inverse image contains +z and -z.
#
# The factorization construction normally has the relevant x near x_true.
# Since z_true is huge and positive in these test instances, the local
# branch is the positive branch:
#
#     z_abs_low <= z <= z_abs_high
#
# Its count is simply:
#
#     z_abs_high-z_abs_low+1
#
# and therefore the exact number of x candidates on the positive branch
# is identical.
#
# ======================================================================================================================

def exact_positive_branch_count(
    inst: Instance,
    A4: int,
    k_low: int,
    k_high: int,
) -> tuple[int, int, int] | None:

    sq_low = A4 - k_high
    sq_high = A4 - k_low

    if sq_high < 0:
        return None

    sq_low = max(0, sq_low)

    lo = ceil_sqrt(sq_low)
    hi = isqrt(sq_high)

    if lo > hi:
        return None

    count = hi - lo + 1

    x_lo = x_from_z(inst, lo)
    x_hi = x_from_z(inst, hi)

    return count, x_lo, x_hi


# ======================================================================================================================
# BRUTE CROSS-CHECK
# ======================================================================================================================

def brute_positive_branch(
    inst: Instance,
    A4: int,
    k_low: int,
    k_high: int,
    lo_z: int,
    hi_z: int,
) -> tuple[int, bool]:

    count = 0
    true_survives = False

    for z in range(lo_z, hi_z + 1):

        K = k_from_z(A4, z)

        if k_low <= K <= k_high:

            count += 1

            x = x_from_z(inst, z)

            if x == inst.x_true:
                true_survives = True

    return count, true_survives


# ======================================================================================================================
# LOCAL K SPACING
# ======================================================================================================================

def exact_k_step(
    inst: Instance,
    A4: int,
) -> tuple[int, int]:

    x = inst.x_true
    z = z_from_x(inst, x)

    if z is None:
        raise AssertionError("true x parity invalid")

    k0 = k_from_z(A4, z)
    k1 = k_from_z(A4, z + 1)

    exact = k1 - k0

    predicted = -(2 * z + 1)

    return exact, predicted


# ======================================================================================================================
# MAIN
# ======================================================================================================================

def main() -> None:

    print("=" * 120)
    print("EXPERIMENT 445R2")
    print("=" * 120)
    print()
    print("O(1) EXACT K-INTERVAL -> X-COUNT / MONOTONE PULLBACK AUDIT")
    print()
    print("NO GIANT K ENUMERATION")
    print("NO GIANT X ENUMERATION")
    print("NO CRT CARTESIAN PRODUCT")
    print()

    total_tests = 0
    total_true = 0
    total_empty = 0
    total_brute = 0
    total_brute_fail = 0

    instances = build_instances()

    for inst in instances:

        K_true = hidden_K(inst)
        r = build_r(inst, K_true)

        A = r + inst.d0 * inst.d0

        if A % 4 != 0:
            raise AssertionError(
                f"A mod 4 != 0 at instance={inst.idx}"
            )

        A4 = A // 4

        z_true = z_from_x(
            inst,
            inst.x_true,
        )

        if z_true is None:
            raise AssertionError(
                f"true x parity invalid at instance={inst.idx}"
            )

        reconstructed_K = k_from_z(
            A4,
            z_true,
        )

        if reconstructed_K != K_true:
            raise AssertionError(
                f"true K identity failed at instance={inst.idx}"
            )

        exact_step, predicted_step = exact_k_step(
            inst,
            A4,
        )

        if exact_step != predicted_step:
            raise AssertionError(
                f"K step identity failed at instance={inst.idx}"
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
            f"  true z = {z_true}"
        )
        print()
        print("  LOCAL K-SPACING")
        print(
            f"    exact K(z+1)-K(z) = {exact_step}"
        )
        print(
            f"    predicted         = {predicted_step}"
        )
        print(
            "    identity = True"
        )
        print()

        print("  EXACT BOUNDED-K -> POSITIVE-X-BRANCH")
        print()
        print(
            "    radius"
            "        K-width"
            "   x-candidates"
            "          x-low"
            "         x-high"
            "   true-in"
            " brute"
        )
        print(
            "    " + "-" * 90
        )

        for radius in K_RADII:

            k_low = K_true - radius
            k_high = K_true + radius

            total_tests += 1

            full = exact_pullback(
                inst,
                A4,
                k_low,
                k_high,
                K_true,
            )

            positive = exact_positive_branch_count(
                inst,
                A4,
                k_low,
                k_high,
            )

            if full is None or positive is None:

                total_empty += 1

                print(
                    f"    {radius:>10}"
                    f"{2 * radius + 1:>16}"
                    f"{0:>16}"
                    f"{'N/A':>16}"
                    f"{'N/A':>16}"
                    f"{False:>10}"
                    f"{'NONE':>8}"
                )

                continue

            if full.true_in:
                total_true += 1

            count, x_low, x_high = positive

            brute_text = "SKIP"

            # This is the ONLY place where actual enumeration occurs.
            # It is bounded by BRUTE_FORCE_LIMIT.
            if count <= BRUTE_FORCE_LIMIT:

                lo_z = exact_positive_branch_count(
                    inst,
                    A4,
                    k_low,
                    k_high,
                )

                if lo_z is None:
                    raise AssertionError("unexpected empty positive branch")

                # Recover the z bounds exactly.
                sq_low = max(
                    0,
                    A4 - k_high,
                )
                sq_high = A4 - k_low

                brute_lo = ceil_sqrt(sq_low)
                brute_hi = isqrt(sq_high)

                brute_count, brute_true = brute_positive_branch(
                    inst,
                    A4,
                    k_low,
                    k_high,
                    brute_lo,
                    brute_hi,
                )

                total_brute += 1

                expected_true = (
                    k_low <= K_true <= k_high
                )

                if brute_count != count:
                    total_brute_fail += 1

                    raise AssertionError(
                        f"brute count mismatch "
                        f"instance={inst.idx}, "
                        f"radius={radius}"
                    )

                if brute_true != expected_true:
                    total_brute_fail += 1

                    raise AssertionError(
                        f"brute true mismatch "
                        f"instance={inst.idx}, "
                        f"radius={radius}"
                    )

                brute_text = f"{brute_count}"

            print(
                f"    {radius:>10}"
                f"{2 * radius + 1:>16}"
                f"{count:>16}"
                f"{x_low:>16}"
                f"{x_high:>16}"
                f"{str(full.true_in):>10}"
                f"{brute_text:>8}"
            )

        print()

        # ==================================================================
        # DIRECT FORMULA FOR SMALL K RADIUS
        # ==================================================================

        print("  LOCAL WIDTH CHECK")

        for radius in [0, 1, 2, 4, 8, 16, 32]:

            k_low = K_true - radius
            k_high = K_true + radius

            positive = exact_positive_branch_count(
                inst,
                A4,
                k_low,
                k_high,
            )

            if positive is None:
                count = 0
            else:
                count = positive[0]

            # A centered K interval around the true K must contain
            # the true point, so the positive branch should be nonempty.
            if count == 0:
                raise AssertionError(
                    f"unexpected local empty branch "
                    f"instance={inst.idx}, radius={radius}"
                )

            print(
                f"    radius={radius:>3}"
                f"  positive-x-count={count:>3}"
            )

        print()

    # ======================================================================
    # SUMMARY
    # ======================================================================

    print("=" * 120)
    print("GLOBAL EXPERIMENT 445R2 SUMMARY")
    print("=" * 120)

    print(
        f"  exact K interval tests       = {total_tests}"
    )
    print(
        f"  empty pullbacks              = {total_empty}"
    )
    print(
        f"  true-containing intervals    = {total_true}"
    )
    print(
        f"  brute cross-checks           = {total_brute}"
    )
    print(
        f"  brute cross-check failures   = {total_brute_fail}"
    )
    print()

    print("CORE FORMULA")
    print()
    print("  K = A4 - z^2")
    print("  x = 2z - d0")
    print()
    print("  K_low <= K <= K_high")
    print("      <=>")
    print("  A4-K_high <= z^2 <= A4-K_low")
    print()
    print("  Therefore the exact positive-branch x population")
    print("  is obtained from TWO integer square roots.")
    print()
    print("  No candidate enumeration is needed.")
    print()

    print("LOCAL MONOTONICITY")
    print()
    print("  K(z+1)-K(z) = -(2z+1)")
    print()
    print("  Equivalently, for parity-valid x:")
    print()
    print("  K(x+2)-K(x) = -(x+d0)-1")
    print()
    print("  Hence on the positive z branch, K decreases strictly.")
    print()

    print("INTERPRETATION")
    print()
    print("  Experiment 444 showed that the CRT conditions become")
    print("  identities after substituting K(x).")
    print()
    print("  Experiment 445R2 now removes the computational issue")
    print("  encountered by the first 445 implementation.")
    print()
    print("  A bounded K interval is inverted analytically.")
    print("  The population is computed in O(1) integer operations")
    print("  per interval rather than by scanning candidate x values.")
    print()
    print("  Therefore very large K radii are safe to test.")
    print()
    print("  The remaining question is purely informational:")
    print()
    print("      how narrow must an independently justified")
    print("      K bound be before only a tiny x population")
    print("      remains?")
    print()
    print("  The tested centered K intervals are diagnostic.")
    print("  They do not constitute independently derived")
    print("  factorization bounds.")
    print()
    print("  This is an exact inverse-image / width experiment,")
    print("  not a factoring theorem.")

    print()
    print("=" * 120)
    print("EXPERIMENT 445R2 FINAL STATUS")
    print("=" * 120)
    print(
        f"  EXACT INTERVAL TESTS       = {total_tests}"
    )
    print(
        f"  TRUE-CONTAINING INTERVALS  = {total_true}"
    )
    print(
        f"  BRUTE CHECKS               = {total_brute}"
    )
    print(
        f"  BRUTE FAILURES             = {total_brute_fail}"
    )
    print("  GIANT K ENUMERATION = False")
    print("  GIANT X ENUMERATION = False")
    print("  CRT CARTESIAN PRODUCT = False")
    print("  O(1) INTERVAL INVERSION = True")
    print("  INTEGER-EXACT = True")
    print("=" * 120)
    print("EXPERIMENT 445R2 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    main()