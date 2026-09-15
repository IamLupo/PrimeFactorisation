#!/usr/bin/env python3
"""
EXPERIMENT 451

EXACT K-PRECISION THRESHOLD / d0 POPULATION AUDIT

QUESTION
  How much independent K precision is required before the true d0
  can enter the exact K -> d0 pullback, and how many d0 candidates
  remain at that threshold?

IMPORTANT
  K_true is NEVER used to construct an experimental/model bound.
  d0_true is NEVER used to construct an experimental/model bound.
  x_true is NEVER used to construct an experimental/model bound.

  True values are used ONLY for post-hoc threshold measurement.

NO GIANT K ENUMERATION
NO GIANT d0 ENUMERATION
NO CRT CARTESIAN PRODUCT
INTEGER-EXACT

CORE RELATION
  4K = r + d0^2 - d^2

  d0^2 = d^2 + 4K - r

KNOWN-DATA REFERENCE
  parity is derived from the exact mod-4 condition.

  x0 = 0 if d is even
  x0 = 1 if d is odd

  d0_ref = d - x0

  K_ref = (r + d0_ref^2 - d^2) / 4

THRESHOLD IDEA

  Let

      DeltaK = |K_true - K_ref|.

  For a symmetric interval

      [K_ref - W, K_ref + W]

  the true solution first becomes reachable at

      W = DeltaK.

  This threshold is measured post-hoc only.

  We then compare:

      count(W = DeltaK - 1)
      count(W = DeltaK)
      count(W = DeltaK + 1)

  and report the exact d0 population at the first reachable radius.

The model widths themselves remain independent diagnostic bounds.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isqrt
from typing import Iterable


# ======================================================================
# DATA
# ======================================================================

@dataclass(frozen=True)
class Instance:
    idx: int
    N: int
    S: int
    d: int
    d0_true: int
    K_true: int


INSTANCES = (
    Instance(
        1,
        14246098189,
        333010,
        14245432170,
        14245432180,
        -3449849766661475506269,
    ),
    Instance(
        2,
        10139117,
        11022,
        10117074,
        10117091,
        -2968347695488785,
    ),
    Instance(
        3,
        10009330297,
        1010042,
        10007310214,
        10007310238,
        -4306289434216722346403,
    ),
    Instance(
        4,
        100460333,
        20046,
        100420242,
        100420273,
        -615138736337991015,
    ),
    Instance(
        5,
        2503701173,
        100074,
        2503501026,
        2503501064,
        -557809093589551261113,
    ),
    Instance(
        6,
        10006200817,
        200062,
        10005800694,
        10005800739,
        -12714738365215418549091,
    ),
    Instance(
        7,
        40005200153,
        400026,
        40004400102,
        40004400154,
        -305667239831581101061223,
    ),
    Instance(
        8,
        270017400119,
        1200024,
        270015000072,
        270015000131,
        -18737381797403407039327539,
    ),
)


# Independent diagnostic width models.
# These are NOT selected from K_true.
MODEL_WIDTHS = (
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
    1073741824,
    4294967296,
)


# ======================================================================
# EXACT ARITHMETIC
# ======================================================================

def mod4_valid_parities(r: int, d: int) -> list[int]:
    """
    Return all d0 parities p in {0,1} satisfying

        r + d0^2 - d^2 == 0 (mod 4)

    Since p^2 == p (mod 2), but d0^2 mod 4 depends on parity,
    direct exact testing is simplest and O(1).
    """
    out = []
    for p in (0, 1):
        if (r + p * p - d * d) % 4 == 0:
            out.append(p)
    return out


def reference_from_known_data(d: int, r: int) -> tuple[int, int, int]:
    """
    Construct a parity-valid known-data reference.

    x0 = 0 for even d
    x0 = 1 for odd d

    d0_ref = d - x0

    Returns:
        x0, d0_ref, K_ref
    """
    valid = mod4_valid_parities(r, d)
    if not valid:
        raise AssertionError("no admissible d0 parity")

    # For this construction the admissible parity is unique.
    if len(valid) != 1:
        raise AssertionError(
            f"expected unique parity class, got {valid}"
        )

    parity = valid[0]

    # Choose the smallest nonnegative x0 with the required parity.
    x0 = parity
    d0_ref = d - x0

    numerator = r + d0_ref * d0_ref - d * d
    if numerator % 4 != 0:
        raise AssertionError("reference K is not integral")

    K_ref = numerator // 4

    return x0, d0_ref, K_ref


def k_from_d0(r: int, d: int, d0: int) -> int:
    numerator = r + d0 * d0 - d * d
    if numerator % 4 != 0:
        raise ValueError("K is not integral for this d0")
    return numerator // 4


def exact_d0_population(
    *,
    d: int,
    r: int,
    k_low: int,
    k_high: int,
) -> tuple[int, int | None, int | None, int]:
    """
    Count nonnegative d0 satisfying:

        K_low <= K(d0) <= K_high

    and the exact mod-4 integrality/parity condition.

    Using

        d0^2 = d^2 + 4K - r,

    we obtain

        d^2 + 4K_low  - r <= d0^2
        <= d^2 + 4K_high - r.

    Returns:
        count,
        d0_low,
        d0_high,
        parity
    """
    if k_low > k_high:
        raise ValueError("invalid K interval")

    valid_parities = mod4_valid_parities(r, d)
    if len(valid_parities) != 1:
        raise AssertionError(
            f"expected one d0 parity class, got {valid_parities}"
        )

    parity = valid_parities[0]

    q_low = d * d + 4 * k_low - r
    q_high = d * d + 4 * k_high - r

    if q_high < 0:
        return 0, None, None, parity

    # d0 >= 0, so negative square roots are excluded.
    # Clamp lower radicand to zero.
    q_low = max(q_low, 0)

    if q_low > q_high:
        return 0, None, None, parity

    lo = isqrt(q_low)
    if lo * lo < q_low:
        lo += 1

    hi = isqrt(q_high)

    if lo > hi:
        return 0, None, None, parity

    # Adjust lower endpoint to required parity.
    if (lo & 1) != parity:
        lo += 1

    # Adjust upper endpoint to required parity.
    if (hi & 1) != parity:
        hi -= 1

    if lo > hi:
        return 0, None, None, parity

    count = ((hi - lo) // 2) + 1
    return count, lo, hi, parity


def interval_for_reference(k_ref: int, width: int) -> tuple[int, int]:
    return k_ref - width, k_ref + width


def count_at_width(
    *,
    d: int,
    r: int,
    k_ref: int,
    width: int,
) -> tuple[int, int | None, int | None]:
    k_low, k_high = interval_for_reference(k_ref, width)
    count, lo, hi, _ = exact_d0_population(
        d=d,
        r=r,
        k_low=k_low,
        k_high=k_high,
    )
    return count, lo, hi


def classify_population(
    count: int,
    true_in: bool,
) -> str:
    if count == 0:
        return "EMPTY"
    if count == 1:
        return "UNIQUE"
    if count <= 10:
        return "SMALL"
    if count <= 100:
        return "MEDIUM"
    return "LARGE"


def bit_length_abs(x: int) -> int:
    return abs(x).bit_length()


def ceil_log2_positive(x: int) -> int:
    if x <= 0:
        raise ValueError("x must be positive")
    return (x - 1).bit_length()


def nearest_lower_model_width(delta: int) -> int | None:
    lower = [w for w in MODEL_WIDTHS if w < delta]
    return max(lower) if lower else None


def nearest_upper_model_width(delta: int) -> int | None:
    upper = [w for w in MODEL_WIDTHS if w >= delta]
    return min(upper) if upper else None


# ======================================================================
# BRUTE CROSS-CHECK
# ======================================================================

def brute_population_small(
    *,
    d: int,
    r: int,
    k_ref: int,
    width: int,
) -> int:
    """
    Small diagnostic only.

    This is deliberately limited so it cannot turn into giant enumeration.
    """
    if width > 200_000:
        raise ValueError("brute check width too large")

    k_low, k_high = interval_for_reference(k_ref, width)

    valid_parity = mod4_valid_parities(r, d)
    if len(valid_parity) != 1:
        raise AssertionError("unexpected parity multiplicity")

    parity = valid_parity[0]

    count = 0

    # Enumerate K only for tiny diagnostic widths.
    for k in range(k_low, k_high + 1):
        value = d * d + 4 * k - r
        if value < 0:
            continue

        z = isqrt(value)
        if z * z != value:
            continue

        if z >= 0 and (z & 1) == parity:
            count += 1

    return count


# ======================================================================
# EXPERIMENT
# ======================================================================

def run_instance(inst: Instance) -> dict:
    r = 4 * inst.K_true - inst.d0_true * inst.d0_true + inst.d * inst.d

    # --------------------------------------------------------------
    # Known-data-only reference construction
    # --------------------------------------------------------------
    x0, d0_ref, k_ref = reference_from_known_data(inst.d, r)

    # --------------------------------------------------------------
    # Post-hoc consistency checks
    # --------------------------------------------------------------
    k_true_recovered = k_from_d0(r, inst.d, inst.d0_true)

    if k_true_recovered != inst.K_true:
        raise AssertionError(
            f"true K reconstruction failure: instance={inst.idx}"
        )

    if inst.d - d0_ref != x0:
        raise AssertionError(
            f"reference x identity failure: instance={inst.idx}"
        )

    # IMPORTANT:
    # This is post-hoc measurement only.
    delta_k = abs(inst.K_true - k_ref)

    if delta_k == 0:
        raise AssertionError(
            f"degenerate threshold delta=0: instance={inst.idx}"
        )

    # --------------------------------------------------------------
    # Exact threshold diagnostics
    # --------------------------------------------------------------
    below_width = delta_k - 1
    exact_width = delta_k
    above_width = delta_k + 1

    below_count, below_lo, below_hi = count_at_width(
        d=inst.d,
        r=r,
        k_ref=k_ref,
        width=below_width,
    )

    exact_count, exact_lo, exact_hi = count_at_width(
        d=inst.d,
        r=r,
        k_ref=k_ref,
        width=exact_width,
    )

    above_count, above_lo, above_hi = count_at_width(
        d=inst.d,
        r=r,
        k_ref=k_ref,
        width=above_width,
    )

    # Exact threshold must contain the true d0.
    if not (exact_lo <= inst.d0_true <= exact_hi):
        raise AssertionError(
            f"threshold interval does not contain true d0: instance={inst.idx}"
        )

    if exact_count <= 0:
        raise AssertionError(
            f"threshold population unexpectedly empty: instance={inst.idx}"
        )

    # The interval at delta-1 must not contain the true K.
    k_low_below, k_high_below = interval_for_reference(
        k_ref, below_width
    )
    if k_low_below <= inst.K_true <= k_high_below:
        raise AssertionError(
            f"below-threshold interval contains true K: instance={inst.idx}"
        )

    # --------------------------------------------------------------
    # Model-scale diagnostic bounds
    # These bounds do NOT use K_true.
    # --------------------------------------------------------------
    model_rows = []

    for width in MODEL_WIDTHS:
        count, lo, hi = count_at_width(
            d=inst.d,
            r=r,
            k_ref=k_ref,
            width=width,
        )

        true_in = (
            lo is not None
            and hi is not None
            and lo <= inst.d0_true <= hi
        )

        model_rows.append(
            {
                "width": width,
                "count": count,
                "true_in": true_in,
                "class": classify_population(count, true_in),
                "lo": lo,
                "hi": hi,
            }
        )

    # --------------------------------------------------------------
    # First independent model that reaches the true solution,
    # when one exists in our fixed diagnostic model list.
    # --------------------------------------------------------------
    containing_models = [
        row for row in model_rows if row["true_in"]
    ]

    first_model = (
        min(containing_models, key=lambda row: row["width"])
        if containing_models
        else None
    )

    lower_model = nearest_lower_model_width(delta_k)
    upper_model = nearest_upper_model_width(delta_k)

    # --------------------------------------------------------------
    # Near-threshold local diagnostic counts.
    # --------------------------------------------------------------
    local_widths = []

    candidate_widths = {
        0,
        1,
        delta_k // 4,
        delta_k // 2,
        max(0, delta_k - 1),
        delta_k,
        delta_k + 1,
        delta_k * 2,
    }

    for w in sorted(candidate_widths):
        count, lo, hi = count_at_width(
            d=inst.d,
            r=r,
            k_ref=k_ref,
            width=w,
        )

        true_in = (
            lo is not None
            and hi is not None
            and lo <= inst.d0_true <= hi
        )

        local_widths.append(
            {
                "width": w,
                "count": count,
                "true_in": true_in,
                "class": classify_population(count, true_in),
                "lo": lo,
                "hi": hi,
            }
        )

    # --------------------------------------------------------------
    # Brute cross-check for tiny widths only.
    # --------------------------------------------------------------
    brute_checks = []
    brute_failures = 0

    brute_widths = sorted(
        {w for w in MODEL_WIDTHS if w <= 200_000}
    )

    for width in brute_widths:
        exact_count, _, _ = count_at_width(
            d=inst.d,
            r=r,
            k_ref=k_ref,
            width=width,
        )

        brute_count = brute_population_small(
            d=inst.d,
            r=r,
            k_ref=k_ref,
            width=width,
        )

        ok = exact_count == brute_count
        if not ok:
            brute_failures += 1

        brute_checks.append(
            {
                "width": width,
                "exact": exact_count,
                "brute": brute_count,
                "ok": ok,
            }
        )

    return {
        "instance": inst,
        "r": r,
        "x0": x0,
        "d0_ref": d0_ref,
        "k_ref": k_ref,
        "delta_k": delta_k,
        "below": {
            "width": below_width,
            "count": below_count,
            "lo": below_lo,
            "hi": below_hi,
        },
        "exact": {
            "width": exact_width,
            "count": exact_count,
            "lo": exact_lo,
            "hi": exact_hi,
        },
        "above": {
            "width": above_width,
            "count": above_count,
            "lo": above_lo,
            "hi": above_hi,
        },
        "model_rows": model_rows,
        "first_model": first_model,
        "lower_model": lower_model,
        "upper_model": upper_model,
        "local_widths": local_widths,
        "brute_checks": brute_checks,
        "brute_failures": brute_failures,
    }


# ======================================================================
# OUTPUT
# ======================================================================

def print_instance(result: dict) -> None:
    inst: Instance = result["instance"]

    print("-" * 120)
    print(
        f"INSTANCE {inst.idx}: "
        f"N={inst.N} S={inst.S} d={inst.d} "
        f"r={result['r']}"
    )
    print()

    print("KNOWN-DATA REFERENCE")
    print(f"  x0     = {result['x0']}")
    print(f"  d0_ref = {result['d0_ref']}")
    print(f"  K_ref  = {result['k_ref']}")
    print()

    print("POST-HOC TRUE VALUES")
    print(f"  true d0 = {inst.d0_true}")
    print(f"  true x  = {inst.d - inst.d0_true}")
    print(f"  true K  = {inst.K_true}")
    print()

    print("EXACT CONSISTENCY")
    print(
        "  K(d0_true) == K_true : "
        f"{k_from_d0(result['r'], inst.d, inst.d0_true) == inst.K_true}"
    )
    print()

    print("K-PRECISION THRESHOLD")
    print(f"  |K_true-K_ref| = {result['delta_k']}")
    print(
        f"  threshold bits = "
        f"{bit_length_abs(result['delta_k'])}"
    )
    print(
        f"  ceil(log2(delta+1)) = "
        f"{ceil_log2_positive(result['delta_k'] + 1)}"
    )
    print()

    print("THRESHOLD CROSSING")
    print(
        "  width        count        d0-low        d0-high      true-in"
    )
    print("  " + "-" * 82)

    for label, row in (
        ("delta-1", result["below"]),
        ("delta", result["exact"]),
        ("delta+1", result["above"]),
    ):
        true_in = (
            row["lo"] is not None
            and row["hi"] is not None
            and row["lo"] <= inst.d0_true <= row["hi"]
        )

        print(
            f"  {label:<9} "
            f"{row['width']:>14} "
            f"{row['count']:>12} "
            f"{str(row['lo']):>14} "
            f"{str(row['hi']):>14} "
            f"{str(true_in):>10}"
        )

    print()
    print("INTERPRETATION OF THRESHOLD")
    print(
        "  delta-1 excludes K_true = ",
        result["instance"].K_true not in range(
            result["k_ref"] - result["below"]["width"],
            result["k_ref"] + result["below"]["width"] + 1,
        ),
    )
    print(
        "  delta contains K_true  = ",
        result["k_ref"] - result["exact"]["width"]
        <= inst.K_true
        <= result["k_ref"] + result["exact"]["width"],
    )
    print(
        f"  exact d0 population at first reachable width = "
        f"{result['exact']['count']}"
    )
    print()

    print("NEAREST DIAGNOSTIC MODEL WIDTHS")
    print(
        f"  lower model width < delta : "
        f"{result['lower_model']}"
    )
    print(
        f"  upper model width >= delta: "
        f"{result['upper_model']}"
    )

    first_model = result["first_model"]
    if first_model is None:
        print("  first model containing true d0 = NONE")
    else:
        print(
            "  first model containing true d0 = "
            f"{first_model['width']}"
        )
        print(
            "  population at first model     = "
            f"{first_model['count']}"
        )
        print(
            "  d0 range                      = "
            f"[{first_model['lo']}, {first_model['hi']}]"
        )

    print()
    print("NEAR-THRESHOLD POPULATION")
    print(
        "  width                         d0-count   true-in   class"
    )
    print("  " + "-" * 82)

    for row in result["local_widths"]:
        print(
            f"  {row['width']:>28} "
            f"{row['count']:>12} "
            f"{str(row['true_in']):>10} "
            f"{row['class']:>10}"
        )

    print()
    print("BRUTE CROSS-CHECK")
    print(
        f"  tiny-width checks = "
        f"{len(result['brute_checks'])}"
    )
    print(
        f"  failures          = "
        f"{result['brute_failures']}"
    )

    if result["brute_failures"] != 0:
        raise AssertionError(
            f"brute mismatch: instance={inst.idx}"
        )

    print()
    print("THRESHOLD STRUCTURE")
    print(
        "  A symmetric interval around K_ref reaches the true"
    )
    print(
        "  solution for the first time at W = |K_true-K_ref|."
    )
    print(
        "  The threshold is a post-hoc measurement, not a bound."
    )
    print()


def main() -> None:
    print("=" * 120)
    print("EXPERIMENT 451")
    print("=" * 120)
    print()
    print("EXACT K-PRECISION THRESHOLD / d0 POPULATION AUDIT")
    print()
    print("PURPOSE")
    print(
        "  Measure the independent K precision required to reach "
        "the true d0 branch."
    )
    print()
    print("RULES")
    print("  K_true is NEVER used to construct an experimental bound")
    print("  d0_true is NEVER used to construct an experimental bound")
    print("  x_true is NEVER used to construct an experimental bound")
    print("  no giant K enumeration")
    print("  no giant d0 enumeration")
    print("  no CRT Cartesian product")
    print("  integer-exact")
    print()

    results = []

    for inst in INSTANCES:
        result = run_instance(inst)
        results.append(result)
        print_instance(result)

    # ==================================================================
    # GLOBAL SUMMARY
    # ==================================================================

    total = len(results)
    threshold_populations = [r["exact"]["count"] for r in results]

    model_reaches_true = sum(
        r["first_model"] is not None
        for r in results
    )

    unique_at_threshold = sum(
        r["exact"]["count"] == 1
        for r in results
    )

    small_at_threshold = sum(
        2 <= r["exact"]["count"] <= 10
        for r in results
    )

    medium_at_threshold = sum(
        11 <= r["exact"]["count"] <= 100
        for r in results
    )

    large_at_threshold = sum(
        r["exact"]["count"] > 100
        for r in results
    )

    brute_failures = sum(
        r["brute_failures"]
        for r in results
    )

    print("=" * 120)
    print("GLOBAL EXPERIMENT 451 SUMMARY")
    print("=" * 120)

    print(f"  instances                         = {total}")
    print(
        f"  diagnostic model reaches true    = "
        f"{model_reaches_true}/{total}"
    )
    print(
        f"  unique d0 at exact threshold     = "
        f"{unique_at_threshold}"
    )
    print(
        f"  small d0 at exact threshold      = "
        f"{small_at_threshold}"
    )
    print(
        f"  medium d0 at exact threshold     = "
        f"{medium_at_threshold}"
    )
    print(
        f"  large d0 at exact threshold      = "
        f"{large_at_threshold}"
    )
    print(f"  brute failures                   = {brute_failures}")
    print()

    print("EXACT THRESHOLD FORMULA")
    print()
    print("  K interval:")
    print()
    print("      [K_ref-W, K_ref+W]")
    print()
    print("  true solution first enters when:")
    print()
    print("      W = |K_true-K_ref|")
    print()
    print("  This is a post-hoc information requirement,")
    print("  not an independently justified bound.")
    print()

    print("EXACT d0 PULLBACK")
    print()
    print("      d^2 + 4(K_ref-W) - r")
    print("        <= d0^2 <=")
    print("      d^2 + 4(K_ref+W) - r")
    print()
    print("  The population is computed from integer square roots")
    print("  and the unique admissible parity class.")
    print()

    print("THRESHOLD TABLE")
    print()
    print(
        "  inst   DeltaK                  bits   threshold-count   "
        "true d0 offset"
    )
    print("  " + "-" * 100)

    for r in results:
        inst = r["instance"]
        delta = r["delta_k"]
        offset = abs(inst.d0_true - r["d0_ref"])

        print(
            f"  {inst.idx:>4} "
            f"{delta:>24} "
            f"{bit_length_abs(delta):>8} "
            f"{r['exact']['count']:>17} "
            f"{offset:>18}"
        )

    print()
    print("INTERPRETATION")
    print()
    print(
        "  This experiment separates two quantities:"
    )
    print()
    print(
        "    1. K precision required to reach the correct d0 branch."
    )
    print(
        "    2. Remaining d0 ambiguity once that precision is available."
    )
    print()
    print(
        "  A unique/very small threshold population is NOT by itself"
    )
    print(
        "  evidence of a factoring method, because the threshold is"
    )
    print(
        "  measured using the hidden solution after the fact."
    )
    print()
    print(
        "  The meaningful next question is whether an independently"
    )
    print(
        "  derivable bound on K can achieve comparable precision."
    )
    print()

    if brute_failures != 0:
        raise AssertionError("global brute cross-check failure")

    print("=" * 120)
    print("EXPERIMENT 451 FINAL STATUS")
    print("=" * 120)
    print(f"  INSTANCES                     = {total}")
    print(f"  BRUTE FAILURES                = {brute_failures}")
    print("  O(1) THRESHOLD INVERSION      = True")
    print("  O(1) d0 POPULATION            = True")
    print("  GIANT K ENUMERATION           = False")
    print("  GIANT d0 ENUMERATION          = False")
    print("  CRT CARTESIAN PRODUCT         = False")
    print("  INTEGER-EXACT                 = True")
    print("  TRUE VALUES USED FOR BOUNDS   = False")
    print()
    print(
        "  CONCLUSION:"
    )
    print(
        "    The first reachable K radius can be measured exactly,"
    )
    print(
        "    and the remaining d0 population at that precision can"
    )
    print(
        "    be computed in O(1) arithmetic."
    )
    print("=" * 120)
    print("EXPERIMENT 451 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    main()
