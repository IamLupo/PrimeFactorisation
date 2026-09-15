#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 450R3
==============================================================================

KNOWN-DATA K-PRECISION -> EXACT d0 POPULATION AUDIT

CORRECTIONS
  1. d0 parity is derived from the exact mod-4 integrality condition.
  2. No global monotonicity assertion is made for an arbitrary model list.
  3. Nested-interval monotonicity is checked explicitly.
  4. Every population is obtained by exact square-root inversion.

CORE IDENTITY

    4K = r + d0^2 - d^2

therefore

    d0^2 = d^2 + 4K - r

RULES

    K_true is NEVER used to construct a bound.
    d0_true is NEVER used to construct a bound.
    x_true is NEVER used to construct a bound.

    no giant K enumeration
    no giant d0 enumeration
    no CRT Cartesian product
    integer-exact
==============================================================================
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isqrt


# ============================================================================
# DATA
# ============================================================================

@dataclass(frozen=True)
class Instance:
    idx: int
    N: int
    S: int
    d: int
    r: int

    # POST-HOC ONLY
    true_d0: int
    true_x: int


INSTANCES = [
    Instance(
        1,
        14246098189,
        333010,
        14245432170,
        -13799399066930810668576,
        14245432180,
        -10,
    ),
    Instance(
        2,
        10139117,
        11022,
        10117074,
        -11873391125935945,
        10117091,
        -17,
    ),
    Instance(
        3,
        10009330297,
        1010042,
        10007310214,
        -17225157737347240276460,
        10007310238,
        -24,
    ),
    Instance(
        4,
        100460333,
        20046,
        100420242,
        -2460554951578020025,
        100420273,
        -31,
    ),
    Instance(
        5,
        2503701173,
        100074,
        2503501026,
        -2231236374548471123872,
        2503501064,
        -38,
    ),
    Instance(
        6,
        10006200817,
        200062,
        10005800694,
        -50858953461762196260849,
        10005800739,
        -45,
    ),
    Instance(
        7,
        40005200153,
        400026,
        40004400102,
        -1222668959330484861858204,
        40004400154,
        -52,
    ),
    Instance(
        8,
        270017400119,
        1200024,
        270015000072,
        -74949527189645489927322133,
        270015000131,
        -59,
    ),
]


# ============================================================================
# EXACT CORE
# ============================================================================

def exact_k(inst: Instance, d0: int) -> int:
    num = inst.r + d0 * d0 - inst.d * inst.d

    if num % 4 != 0:
        raise AssertionError(
            f"non-integral K(d0): instance={inst.idx}, d0={d0}"
        )

    return num // 4


def valid_d0_parities(inst: Instance) -> list[int]:
    """
    Solve

        r + d0^2 - d^2 == 0 (mod 4)

    for d0 parity.
    """
    out = []

    for p in (0, 1):
        if (inst.r + p * p - inst.d * inst.d) % 4 == 0:
            out.append(p)

    return out


def required_d0_parity(inst: Instance) -> int:
    vals = valid_d0_parities(inst)

    if len(vals) != 1:
        raise AssertionError(
            f"unexpected d0 parity classes: "
            f"instance={inst.idx}, values={vals}"
        )

    return vals[0]


def known_reference(inst: Instance) -> tuple[int, int, int]:
    """
    Choose x0 in {0,1} so that

        d0_ref = d - x0

    lies in the actual admissible parity class.
    """
    parity = required_d0_parity(inst)

    for x0 in (0, 1):
        d0_ref = inst.d - x0

        if (d0_ref & 1) == parity:
            return x0, d0_ref, exact_k(inst, d0_ref)

    raise AssertionError(
        f"unable to construct parity-compatible reference: "
        f"instance={inst.idx}"
    )


# ============================================================================
# EXACT K -> d0 INVERSION
# ============================================================================

def interval_d0_population(
    inst: Instance,
    k_low: int,
    k_high: int,
) -> tuple[int, int, int]:
    """
    Count all d0 satisfying

        k_low <= K(d0) <= k_high

    exactly.

    Returns:
        count, smallest nonnegative d0, largest nonnegative d0

    No candidate enumeration.
    """
    if k_low > k_high:
        raise ValueError("k_low > k_high")

    T_low = inst.d * inst.d + 4 * k_low - inst.r
    T_high = inst.d * inst.d + 4 * k_high - inst.r

    if T_high < 0:
        return 0, 0, -1

    T_low = max(T_low, 0)

    lo = isqrt(T_low)
    if lo * lo < T_low:
        lo += 1

    hi = isqrt(T_high)

    if lo > hi:
        return 0, lo, hi

    parity = required_d0_parity(inst)

    if (lo & 1) != parity:
        lo += 1

    if (hi & 1) != parity:
        hi -= 1

    if lo > hi:
        return 0, lo, hi

    count = (hi - lo) // 2 + 1

    return count, lo, hi


def interval_contains(
    k_low: int,
    k_high: int,
    k: int,
) -> bool:
    return k_low <= k <= k_high


def classify(count: int) -> str:
    if count == 0:
        return "EMPTY"
    if count == 1:
        return "UNIQUE"
    if count <= 10:
        return "SMALL"
    return "LARGE"


# ============================================================================
# KNOWN-DATA WIDTHS
# ============================================================================

def width_models(inst: Instance, d0_ref: int):
    """
    All widths are constructed from known data only.

    The reference K itself is also known-data-derived.
    """
    nroot = isqrt(inst.N)

    models: list[tuple[str, int]] = [
        ("0", 0),
        ("1", 1),
        ("2", 2),
        ("4", 4),
        ("8", 8),
        ("16", 16),
        ("32", 32),
        ("64", 64),
        ("128", 128),
        ("256", 256),
        ("512", 512),
        ("1024", 1024),
        ("2048", 2048),
        ("4096", 4096),
        ("8192", 8192),
        ("16384", 16384),
        ("32768", 32768),
        ("65536", 65536),
        ("131072", 131072),
        ("262144", 262144),
        ("524288", 524288),
        ("1048576", 1048576),
        ("4194304", 4194304),
        ("16777216", 16777216),
        ("67108864", 67108864),
        ("268435456", 268435456),
        ("1073741824", 1073741824),
        ("4294967296", 4294967296),

        ("S", inst.S),
        ("sqrtN", nroot),
        ("2sqrtN", 2 * nroot),
        ("4sqrtN", 4 * nroot),
        ("8sqrtN", 8 * nroot),
        ("32sqrtN", 32 * nroot),
        ("d", inst.d),
        ("2d", 2 * inst.d),
        ("N", inst.N),
        ("N*S", inst.N * inst.S),
        ("S^2", inst.S * inst.S),
    ]

    return models


# ============================================================================
# NESTED WIDTH MONOTONICITY
# ============================================================================

def check_nested_widths(
    inst: Instance,
    k_ref: int,
    widths: list[tuple[str, int]],
) -> tuple[int, list[tuple[str, int, int]]]:
    """
    Sort distinct widths numerically and verify:

        W2 >= W1
        =>
        I(W1) subset I(W2)
        =>
        count(W2) >= count(W1)

    This is the mathematically valid monotonicity test.

    We do NOT assume the original model ordering is sorted.
    """
    normalized = sorted(
        {(name, W) for name, W in widths},
        key=lambda x: (x[1], x[0]),
    )

    previous_count = None
    failures = []

    rows = []

    for name, W in normalized:
        k_low = k_ref - W
        k_high = k_ref + W

        count, _, _ = interval_d0_population(
            inst,
            k_low,
            k_high,
        )

        rows.append((name, W, count))

        if previous_count is not None and count < previous_count:
            failures.append(
                (name, W, count)
            )

        previous_count = count

    return len(failures), failures


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    print("=" * 120)
    print("EXPERIMENT 450R3")
    print("=" * 120)
    print()
    print("KNOWN-DATA K-PRECISION -> EXACT d0 POPULATION AUDIT")
    print()
    print("CORRECTIONS")
    print("  1. d0 parity is derived from the exact mod-4 condition.")
    print("  2. Width models are not assumed to be ordered.")
    print("  3. Monotonicity is checked only on nested intervals.")
    print("  4. Every d0 population is obtained by O(1) exact arithmetic.")
    print()
    print("RULES")
    print("  K_true is NEVER used to construct a bound")
    print("  d0_true is NEVER used to construct a bound")
    print("  x_true is NEVER used to construct a bound")
    print("  no giant K enumeration")
    print("  no giant d0 enumeration")
    print("  no CRT Cartesian product")
    print("  integer-exact")
    print()

    total_tests = 0
    true_containing = 0
    empty_cases = 0
    unique_cases = 0
    small_cases = 0
    nested_monotonicity_failures = 0

    # ========================================================================
    # INSTANCES
    # ========================================================================

    for inst in INSTANCES:
        print("-" * 120)
        print(
            f"INSTANCE {inst.idx}: "
            f"N={inst.N} "
            f"S={inst.S} "
            f"d={inst.d} "
            f"r={inst.r}"
        )
        print()

        # --------------------------------------------------------------------
        # PARITY
        # --------------------------------------------------------------------

        valid_parities = valid_d0_parities(inst)
        parity = required_d0_parity(inst)

        print("  PARITY DERIVATION")
        print(f"    valid d0 parities = {valid_parities}")
        print(f"    required parity   = {parity}")
        print()

        # --------------------------------------------------------------------
        # REFERENCE
        # --------------------------------------------------------------------

        x0, d0_ref, k_ref = known_reference(inst)

        print("  KNOWN-DATA REFERENCE")
        print(f"    x0     = {x0}")
        print(f"    d0_ref = {d0_ref}")
        print(f"    K_ref  = {k_ref}")
        print()

        # --------------------------------------------------------------------
        # POST-HOC TRUE CHECK
        # --------------------------------------------------------------------

        k_true = exact_k(inst, inst.true_d0)

        print("  POST-HOC TRUE VALUES")
        print(f"    true d0 = {inst.true_d0}")
        print(f"    true x  = {inst.true_x}")
        print(f"    true K  = {k_true}")
        print()

        if inst.true_d0 % 2 != parity:
            raise AssertionError(
                f"true d0 parity failure: instance={inst.idx}"
            )

        # --------------------------------------------------------------------
        # REFERENCE IDENTITIES
        # --------------------------------------------------------------------

        if d0_ref != inst.d - x0:
            raise AssertionError(
                f"reference x/d0 identity failure: instance={inst.idx}"
            )

        if (inst.r + d0_ref * d0_ref - inst.d * inst.d) % 4 != 0:
            raise AssertionError(
                f"reference K integrality failure: instance={inst.idx}"
            )

        print("  REFERENCE IDENTITIES")
        print(
            f"    d - d0_ref = {inst.d - d0_ref}"
            f"  == x0 = {x0}"
        )
        print(
            "    K integrality = True"
        )
        print()

        # --------------------------------------------------------------------
        # LOCAL K SPACING
        # --------------------------------------------------------------------

        k_next = exact_k(inst, d0_ref + 2)
        predicted = d0_ref + 1

        if k_next - k_ref != predicted:
            raise AssertionError(
                f"local spacing failure: instance={inst.idx}"
            )

        print("  LOCAL SPACING")
        print(f"    K(d0+2)-K(d0) = {k_next-k_ref}")
        print(f"    predicted      = {predicted}")
        print("    identity       = True")
        print()

        # --------------------------------------------------------------------
        # MODEL AUDIT
        # --------------------------------------------------------------------

        models = width_models(inst, d0_ref)

        print("  K-WIDTH -> d0 POPULATION")
        print()
        print(
            "    model"
            "                     W"
            "                    K-low"
            "                    K-high"
            "       d0-count"
            "   true-in"
            "   class"
        )
        print("    " + "-" * 112)

        for name, W in models:
            k_low = k_ref - W
            k_high = k_ref + W

            count, d0_lo, d0_hi = interval_d0_population(
                inst,
                k_low,
                k_high,
            )

            contains = interval_contains(
                k_low,
                k_high,
                k_true,
            )

            total_tests += 1

            if contains:
                true_containing += 1

            if count == 0:
                empty_cases += 1

            if count == 1:
                unique_cases += 1

            if 2 <= count <= 10:
                small_cases += 1

            print(
                f"    {name:<24}"
                f"{W:>22}"
                f"{k_low:>24}"
                f"{k_high:>24}"
                f"{count:>14}"
                f"{str(contains):>10}"
                f"{classify(count):>10}"
            )

        print()

        # --------------------------------------------------------------------
        # NESTED MONOTONICITY CHECK
        # --------------------------------------------------------------------

        failures, failure_rows = check_nested_widths(
            inst,
            k_ref,
            models,
        )

        nested_monotonicity_failures += failures

        print("  NESTED INTERVAL MONOTONICITY")
        print(
            "    count(W) must be nondecreasing after widths are "
            "sorted numerically."
        )
        print(f"    failures = {failures}")

        if failures:
            for name, W, count in failure_rows[:10]:
                print(
                    f"      width={W} model={name} count={count}"
                )

        print()

        # --------------------------------------------------------------------
        # DIAGNOSTIC TRUE-COVERING WINDOW
        # --------------------------------------------------------------------

        true_displacement = abs(k_true - k_ref)

        true_count, true_lo, true_hi = interval_d0_population(
            inst,
            k_ref - true_displacement,
            k_ref + true_displacement,
        )

        if not (
            k_ref - true_displacement
            <= k_true
            <= k_ref + true_displacement
        ):
            raise AssertionError(
                f"true covering window construction failed: "
                f"instance={inst.idx}"
            )

        if not (true_lo <= inst.true_d0 <= true_hi):
            raise AssertionError(
                f"true d0 not inside diagnostic pullback: "
                f"instance={inst.idx}"
            )

        print("  DIAGNOSTIC TRUE-COVERING WINDOW")
        print(f"    |K_true-K_ref| = {true_displacement}")
        print(f"    d0-count      = {true_count}")
        print(f"    d0-range      = [{true_lo}, {true_hi}]")
        print("    true survives  = True")
        print()

        # --------------------------------------------------------------------
        # EXACT x FORM
        # --------------------------------------------------------------------

        x_form = (
            inst.r
            - 2 * inst.d * inst.true_x
            + inst.true_x * inst.true_x
        )

        if 4 * k_true != x_form:
            raise AssertionError(
                f"x-form identity failed: instance={inst.idx}"
            )

    # ========================================================================
    # GLOBAL SUMMARY
    # ========================================================================

    print("=" * 120)
    print("GLOBAL EXPERIMENT 450R3 SUMMARY")
    print("=" * 120)

    print(f"  K-width interval tests        = {total_tests}")
    print(f"  true-containing tests        = {true_containing}")
    print(f"  empty d0 pullbacks            = {empty_cases}")
    print(f"  unique d0 cases               = {unique_cases}")
    print(f"  small d0 cases (2..10)        = {small_cases}")
    print(
        f"  nested monotonicity failures = "
        f"{nested_monotonicity_failures}"
    )
    print()

    print("EXACT PULLBACK")
    print()
    print("  K_low <= K <= K_high")
    print()
    print("      <=>")
    print()
    print("  d^2 + 4K_low - r")
    print("      <= d0^2 <=")
    print("  d^2 + 4K_high - r")
    print()
    print("  The candidate count therefore requires only:")
    print("    - two integer square roots")
    print("    - one parity adjustment")
    print("    - one integer count")
    print()

    print("PARITY")
    print()
    print("  The admissible d0 parity is obtained from")
    print()
    print("      r + d0^2 - d^2 == 0 (mod 4).")
    print()

    print("MONOTONICITY")
    print()
    print("  Arbitrary model ordering is NOT used.")
    print()
    print("  Only genuinely nested intervals are compared.")
    print()
    print("  If W2 >= W1 then")
    print()
    print("      [Kref-W1, Kref+W1]")
    print("          is contained in")
    print("      [Kref-W2, Kref+W2],")
    print()
    print("  so the exact d0 population cannot decrease.")
    print()

    print("INTERPRETATION")
    print()
    print("  This isolates the amount of d0 ambiguity left after")
    print("  specifying a finite K interval.")
    print()
    print("  The computation itself is O(1) per interval.")
    print()
    print("  The important unresolved question remains whether")
    print("  the K interval can be independently justified from")
    print("  the available public data.")
    print()

    # ========================================================================
    # FINAL STRUCTURAL CHECKS
    # ========================================================================

    if nested_monotonicity_failures != 0:
        raise AssertionError(
            "nested interval population monotonicity failed"
        )

    for inst in INSTANCES:
        vals = valid_d0_parities(inst)

        if len(vals) != 1:
            raise AssertionError(
                f"invalid parity structure: instance={inst.idx}"
            )

        k = exact_k(inst, inst.true_d0)

        if 4 * k != (
            inst.r
            + inst.true_d0 * inst.true_d0
            - inst.d * inst.d
        ):
            raise AssertionError(
                f"quadratic identity failure: instance={inst.idx}"
            )

    print("=" * 120)
    print("EXPERIMENT 450R3 FINAL STATUS")
    print("=" * 120)
    print(f"  K-WIDTH INTERVAL TESTS       = {total_tests}")
    print(f"  TRUE-CONTAINING TESTS       = {true_containing}")
    print(f"  EMPTY PULLBACKS             = {empty_cases}")
    print(f"  UNIQUE d0 CASES             = {unique_cases}")
    print(f"  SMALL d0 CASES              = {small_cases}")
    print(
        f"  NESTED MONOTONICITY FAILS   = "
        f"{nested_monotonicity_failures}"
    )
    print("  CORRECT PARITY DERIVATION   = True")
    print("  K_TRUE USED IN BOUNDS       = False")
    print("  d0_TRUE USED IN BOUNDS      = False")
    print("  X_TRUE USED IN BOUNDS       = False")
    print("  GIANT K ENUMERATION         = False")
    print("  GIANT d0 ENUMERATION        = False")
    print("  CRT CARTESIAN PRODUCT       = False")
    print("  O(1) INTERVAL INVERSION     = True")
    print("  INTEGER-EXACT               = True")
    print()
    print("  CONCLUSION:")
    print("    The previous failure was a test-design error:")
    print("    arbitrary model order is not a monotonicity invariant.")
    print("    This version checks only nested K intervals.")
    print("=" * 120)
    print("EXPERIMENT 450R3 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    main()