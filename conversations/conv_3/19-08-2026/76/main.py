#!/usr/bin/env python3
"""
========================================================================================================================
EXPERIMENT 452
========================================================================================================================

KNOWN-DATA K-UNIQUENESS RADIUS / TRUE-ENTRY GAP AUDIT

QUESTION

  Around the known-data reference K_ref, what is the largest symmetric
  K interval that contains exactly one admissible d0?

  How does that known-data uniqueness radius compare with the
  post-hoc K precision required to reach the true d0?

CORE IDENTITY

  4K = r + d0^2 - d^2

  therefore

  d0^2 = d^2 + 4K - r

KNOWN-DATA REFERENCE

  Valid d0 parity is derived directly from

      r + d0^2 - d^2 == 0 (mod 4).

  We choose the canonical known-data reference

      x0 = parity
      d0_ref = d - x0

  and therefore

      K_ref = (r + d0_ref^2 - d^2) / 4.

UNIQUENESS QUESTION

  Consider

      I(W) = [K_ref-W, K_ref+W].

  We seek the largest W for which I(W) contains exactly one
  admissible nonnegative d0.

  Because K(d0) is exact and monotone on each positive d0 branch,
  this threshold can be obtained directly from the nearest other
  admissible d0.

IMPORTANT

  All experimental bounds and uniqueness radii are constructed
  exclusively from known data.

  K_true, d0_true and x_true are used only for POST-HOC comparison.

NO GIANT K ENUMERATION
NO GIANT d0 ENUMERATION
NO CRT CARTESIAN PRODUCT
INTEGER-EXACT
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isqrt


# ======================================================================
# DATA
# ======================================================================

@dataclass(frozen=True)
class Instance:
    idx: int
    N: int
    S: int
    d: int

    # Hidden values: used ONLY post-hoc.
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


# ======================================================================
# EXACT CORE
# ======================================================================

def valid_d0_parities(r: int, d: int) -> list[int]:
    """
    Exact admissible parities from

        r + d0^2 - d^2 == 0 (mod 4).

    Since only parity is relevant, test 0 and 1 directly.
    """
    out = []

    for p in (0, 1):
        if (r + p * p - d * d) % 4 == 0:
            out.append(p)

    return out


def derive_r(inst: Instance) -> int:
    """
    Recover the observed r from the hidden construction.

    This is only how the synthetic instance is reconstructed.
    The experiment itself never uses K_true to construct a bound.
    """
    return (
        4 * inst.K_true
        - inst.d0_true * inst.d0_true
        + inst.d * inst.d
    )


def known_reference(r: int, d: int) -> tuple[int, int, int]:
    """
    Construct the canonical known-data reference.

        x0 = required parity
        d0_ref = d - x0

    Returns:
        x0, d0_ref, K_ref
    """
    parities = valid_d0_parities(r, d)

    if len(parities) != 1:
        raise AssertionError(
            f"expected one admissible parity, got {parities}"
        )

    parity = parities[0]
    x0 = parity
    d0_ref = d - x0

    numerator = r + d0_ref * d0_ref - d * d

    if numerator % 4 != 0:
        raise AssertionError("reference K not integral")

    k_ref = numerator // 4

    return x0, d0_ref, k_ref


def k_of_d0(r: int, d: int, d0: int) -> int:
    numerator = r + d0 * d0 - d * d

    if numerator % 4 != 0:
        raise AssertionError(
            f"non-integral K for d0={d0}"
        )

    return numerator // 4


def interval_count(
    r: int,
    d: int,
    k_low: int,
    k_high: int,
) -> tuple[int, int | None, int | None]:
    """
    Exact O(1) pullback.

        K_low <= K <= K_high

    iff

        d^2 + 4K_low - r
          <= d0^2 <=
        d^2 + 4K_high - r.

    Only the unique admissible parity class is retained.

    Returns:

        count, d0_low, d0_high
    """
    if k_low > k_high:
        raise ValueError("invalid interval")

    parities = valid_d0_parities(r, d)

    if len(parities) != 1:
        raise AssertionError(
            f"expected unique parity, got {parities}"
        )

    parity = parities[0]

    q_low = d * d + 4 * k_low - r
    q_high = d * d + 4 * k_high - r

    if q_high < 0:
        return 0, None, None

    q_low = max(q_low, 0)

    if q_low > q_high:
        return 0, None, None

    lo = isqrt(q_low)

    if lo * lo < q_low:
        lo += 1

    hi = isqrt(q_high)

    if lo > hi:
        return 0, None, None

    if (lo & 1) != parity:
        lo += 1

    if (hi & 1) != parity:
        hi -= 1

    if lo > hi:
        return 0, None, None

    count = ((hi - lo) // 2) + 1

    return count, lo, hi


def interval_at(k_ref: int, width: int) -> tuple[int, int]:
    return k_ref - width, k_ref + width


def population_at(
    r: int,
    d: int,
    k_ref: int,
    width: int,
) -> tuple[int, int | None, int | None]:
    k_low, k_high = interval_at(k_ref, width)

    return interval_count(
        r,
        d,
        k_low,
        k_high,
    )


# ======================================================================
# UNIQUENESS RADIUS
# ======================================================================

def exact_uniqueness_radius(
    r: int,
    d: int,
    d0_ref: int,
    k_ref: int,
) -> tuple[int, int | None, int | None]:
    """
    Determine the largest W such that the symmetric interval

        [K_ref-W, K_ref+W]

    contains exactly the reference d0.

    The nearest competing admissible d0 is one step of 2 away.

    We inspect both directions where d0 >= 0.

    If K_ref itself is represented by d0_ref, the distance to the
    nearest other represented K is

        min(|K(d0_ref-2)-K_ref|,
            |K(d0_ref+2)-K_ref|).

    The largest integer width containing only the reference is

        nearest_distance - 1.

    Returns:

        uniqueness_width,
        nearest_competitor_d0,
        nearest_competitor_K
    """
    competitors: list[tuple[int, int, int]] = []

    for candidate in (d0_ref - 2, d0_ref + 2):
        if candidate < 0:
            continue

        try:
            k_candidate = k_of_d0(
                r,
                d,
                candidate,
            )
        except AssertionError:
            continue

        distance = abs(k_candidate - k_ref)

        if distance == 0:
            raise AssertionError(
                "distinct d0 values mapped to identical K"
            )

        competitors.append(
            (
                distance,
                candidate,
                k_candidate,
            )
        )

    if not competitors:
        raise AssertionError(
            "reference d0 has no admissible neighboring competitor"
        )

    competitors.sort(key=lambda row: row[0])

    nearest_distance, nearest_d0, nearest_k = competitors[0]

    uniqueness_width = nearest_distance - 1

    return (
        uniqueness_width,
        nearest_d0,
        nearest_k,
    )


# ======================================================================
# PRECISION GAP
# ======================================================================

def signed_gap(delta_true: int, uniqueness_width: int) -> int:
    """
    Positive means the true-entry radius lies outside the
    known-data uniqueness interval.

        gap = Delta_true - U

    Negative would mean the true solution already lies inside
    the uniqueness interval.
    """
    return delta_true - uniqueness_width


def ratio_floor(delta_true: int, uniqueness_width: int) -> int:
    if uniqueness_width <= 0:
        return 0
    return delta_true // uniqueness_width


def exact_ratio_text(delta_true: int, uniqueness_width: int) -> str:
    if uniqueness_width == 0:
        return "INF"

    q = delta_true // uniqueness_width
    rem = delta_true % uniqueness_width

    return f"{q} + {rem}/{uniqueness_width}"


# ======================================================================
# POSTSCRIPT DIAGNOSTICS
# ======================================================================

def classify(count: int) -> str:
    if count == 0:
        return "EMPTY"
    if count == 1:
        return "UNIQUE"
    if count <= 10:
        return "SMALL"
    if count <= 100:
        return "MEDIUM"
    return "LARGE"


def test_known_uniqueness_threshold(
    r: int,
    d: int,
    d0_ref: int,
    k_ref: int,
    uniqueness_width: int,
) -> tuple[int, int, int]:
    """
    Confirm the exact boundary:

      W = U       -> exactly one d0
      W = U+1     -> at least two d0, provided the nearest
                     competitor is at an integer K distance.
    """
    count_u, lo_u, hi_u = population_at(
        r,
        d,
        k_ref,
        uniqueness_width,
    )

    count_next, lo_next, hi_next = population_at(
        r,
        d,
        k_ref,
        uniqueness_width + 1,
    )

    if count_u != 1:
        raise AssertionError(
            "uniqueness width did not produce population 1"
        )

    if lo_u != d0_ref or hi_u != d0_ref:
        raise AssertionError(
            "uniqueness interval does not isolate d0_ref"
        )

    if count_next < 2:
        raise AssertionError(
            "uniqueness boundary failed to expose competitor"
        )

    return (
        count_u,
        count_next,
        hi_next - lo_next,
    )


# ======================================================================
# MAIN EXPERIMENT
# ======================================================================

def run_instance(inst: Instance) -> dict:
    r = derive_r(inst)

    x0, d0_ref, k_ref = known_reference(
        r,
        inst.d,
    )

    # Exact construction checks.
    if inst.d - d0_ref != x0:
        raise AssertionError(
            f"reference x mismatch: instance={inst.idx}"
        )

    k_true_check = k_of_d0(
        r,
        inst.d,
        inst.d0_true,
    )

    if k_true_check != inst.K_true:
        raise AssertionError(
            f"true K reconstruction mismatch: instance={inst.idx}"
        )

    # ------------------------------------------------------------------
    # KNOWN-DATA UNIQUENESS RADIUS
    # ------------------------------------------------------------------
    U, competitor_d0, competitor_k = exact_uniqueness_radius(
        r,
        inst.d,
        d0_ref,
        k_ref,
    )

    # ------------------------------------------------------------------
    # Boundary verification
    # ------------------------------------------------------------------
    count_u, count_next, next_span = test_known_uniqueness_threshold(
        r,
        inst.d,
        d0_ref,
        k_ref,
        U,
    )

    # ------------------------------------------------------------------
    # Post-hoc true-entry radius
    # ------------------------------------------------------------------
    delta_true = abs(inst.K_true - k_ref)

    gap = signed_gap(
        delta_true,
        U,
    )

    true_inside_uniqueness = (
        abs(inst.K_true - k_ref) <= U
    )

    # ------------------------------------------------------------------
    # Population at:
    #
    #   U
    #   U+1
    #   Delta-1
    #   Delta
    #   Delta+1
    #
    # These are all O(1).
    # ------------------------------------------------------------------
    diagnostic_widths = sorted(
        {
            0,
            1,
            U,
            U + 1,
            max(0, delta_true - 1),
            delta_true,
            delta_true + 1,
        }
    )

    diagnostic_rows = []

    for width in diagnostic_widths:
        count, lo, hi = population_at(
            r,
            inst.d,
            k_ref,
            width,
        )

        true_in = (
            lo is not None
            and hi is not None
            and lo <= inst.d0_true <= hi
        )

        ref_in = (
            lo is not None
            and hi is not None
            and lo <= d0_ref <= hi
        )

        diagnostic_rows.append(
            {
                "width": width,
                "count": count,
                "lo": lo,
                "hi": hi,
                "true_in": true_in,
                "ref_in": ref_in,
                "class": classify(count),
            }
        )

    return {
        "inst": inst,
        "r": r,
        "x0": x0,
        "d0_ref": d0_ref,
        "k_ref": k_ref,
        "competitor_d0": competitor_d0,
        "competitor_k": competitor_k,
        "uniqueness_width": U,
        "uniqueness_count": count_u,
        "next_count": count_next,
        "next_span": next_span,
        "delta_true": delta_true,
        "gap": gap,
        "true_inside_uniqueness": true_inside_uniqueness,
        "diagnostic_rows": diagnostic_rows,
    }


# ======================================================================
# PRINT
# ======================================================================

def print_instance(result: dict) -> None:
    inst: Instance = result["inst"]

    print("-" * 120)
    print(
        f"INSTANCE {inst.idx}: "
        f"N={inst.N} S={inst.S} d={inst.d} "
        f"r={result['r']}"
    )
    print()

    print("PARITY / KNOWN-DATA REFERENCE")
    print(
        f"  valid d0 parities = "
        f"{valid_d0_parities(result['r'], inst.d)}"
    )
    print(
        f"  x0              = {result['x0']}"
    )
    print(
        f"  d0_ref          = {result['d0_ref']}"
    )
    print(
        f"  K_ref           = {result['k_ref']}"
    )
    print()

    print("POST-HOC TRUE VALUES")
    print(
        f"  true d0 = {inst.d0_true}"
    )
    print(
        f"  true x  = {inst.d - inst.d0_true}"
    )
    print(
        f"  true K  = {inst.K_true}"
    )
    print()

    print("NEAREST KNOWN-DATA COMPETITOR")
    print(
        f"  competitor d0 = {result['competitor_d0']}"
    )
    print(
        f"  competitor K  = {result['competitor_k']}"
    )
    print(
        f"  |K_comp-K_ref| = "
        f"{abs(result['competitor_k'] - result['k_ref'])}"
    )
    print()

    print("EXACT KNOWN-DATA UNIQUENESS RADIUS")
    print(
        f"  U = {result['uniqueness_width']}"
    )
    print(
        "  [K_ref-U, K_ref+U] contains exactly one admissible d0"
    )
    print(
        f"  population at U   = "
        f"{result['uniqueness_count']}"
    )
    print(
        f"  population at U+1 = "
        f"{result['next_count']}"
    )
    print()

    print("TRUE-ENTRY RADIUS")
    print(
        f"  Delta_true = |K_true-K_ref| = "
        f"{result['delta_true']}"
    )
    print(
        f"  Delta bits = "
        f"{result['delta_true'].bit_length()}"
    )
    print()

    print("PRECISION GAP")
    print(
        f"  Delta_true - U = "
        f"{result['gap']}"
    )
    print(
        f"  Delta_true / U = "
        f"{result['delta_true']} / {result['uniqueness_width']}"
    )
    print(
        f"  integer ratio = "
        f"{ratio_floor(result['delta_true'], result['uniqueness_width'])}"
    )
    print(
        f"  exact ratio   = "
        f"{exact_ratio_text(result['delta_true'], result['uniqueness_width'])}"
    )
    print(
        f"  true inside known-data uniqueness interval = "
        f"{result['true_inside_uniqueness']}"
    )
    print()

    print("PRECISION-BOUNDARY POPULATION")
    print(
        "  width                 d0-count   d0-low        d0-high   "
        "ref-in   true-in   class"
    )
    print("  " + "-" * 105)

    for row in result["diagnostic_rows"]:
        print(
            f"  {row['width']:>18} "
            f"{row['count']:>11} "
            f"{str(row['lo']):>14} "
            f"{str(row['hi']):>14} "
            f"{str(row['ref_in']):>8} "
            f"{str(row['true_in']):>9} "
            f"{row['class']:>8}"
        )

    print()

    if result["true_inside_uniqueness"]:
        print("STRUCTURAL OUTCOME")
        print(
            "  The true K lies inside the known-data interval that"
        )
        print(
            "  uniquely identifies the reference d0."
        )
        print(
            "  This would be a significant consistency result."
        )
    else:
        print("STRUCTURAL OUTCOME")
        print(
            "  The true K lies OUTSIDE the known-data uniqueness"
        )
        print(
            "  interval around K_ref."
        )
        print(
            "  Therefore the known-data reference precision is"
        )
        print(
            "  insufficient to reach the true branch."
        )

    print()


# ======================================================================
# GLOBAL SUMMARY
# ======================================================================

def main() -> None:
    print("=" * 120)
    print("EXPERIMENT 452")
    print("=" * 120)
    print()
    print("KNOWN-DATA K-UNIQUENESS RADIUS / TRUE-ENTRY GAP AUDIT")
    print()
    print("QUESTION")
    print(
        "  What is the largest known-data K interval around K_ref"
    )
    print(
        "  that still identifies exactly one admissible d0?"
    )
    print()
    print("RULES")
    print(
        "  K_true is NEVER used to construct a known-data bound"
    )
    print(
        "  d0_true is NEVER used to construct a known-data bound"
    )
    print(
        "  x_true is NEVER used to construct a known-data bound"
    )
    print(
        "  all interval inversion is O(1)"
    )
    print(
        "  integer-exact arithmetic only"
    )
    print()

    results = []

    for inst in INSTANCES:
        result = run_instance(inst)
        results.append(result)
        print_instance(result)

    # ==================================================================
    # GLOBAL METRICS
    # ==================================================================

    true_inside = sum(
        r["true_inside_uniqueness"]
        for r in results
    )

    positive_gap = sum(
        r["gap"] > 0
        for r in results
    )

    equal_gap = sum(
        r["gap"] == 0
        for r in results
    )

    negative_gap = sum(
        r["gap"] < 0
        for r in results
    )

    unique_thresholds = sum(
        r["uniqueness_count"] == 1
        for r in results
    )

    boundary_expands = sum(
        r["next_count"] >= 2
        for r in results
    )

    print("=" * 120)
    print("GLOBAL EXPERIMENT 452 SUMMARY")
    print("=" * 120)
    print(
        f"  instances                           = {len(results)}"
    )
    print(
        f"  exact uniqueness-radius checks      = {unique_thresholds}"
    )
    print(
        f"  U+1 exposes another d0             = {boundary_expands}"
    )
    print(
        f"  true K inside uniqueness region     = {true_inside}"
    )
    print(
        f"  positive precision gaps             = {positive_gap}"
    )
    print(
        f"  zero precision gaps                 = {equal_gap}"
    )
    print(
        f"  negative precision gaps             = {negative_gap}"
    )
    print()

    print("CORE KNOWN-DATA RESULT")
    print()
    print(
        "  Let U be the largest W such that"
    )
    print()
    print(
        "      [K_ref-W, K_ref+W]"
    )
    print()
    print(
        "  contains exactly one admissible d0."
    )
    print()
    print(
        "  U is determined entirely from the nearest alternative"
    )
    print(
        "  admissible d0 around the known-data reference."
    )
    print()
    print(
        "  No hidden value is required to calculate U."
    )
    print()

    print("TRUE-ENTRY COMPARISON")
    print()
    print(
        "  Delta_true = |K_true-K_ref|"
    )
    print()
    print(
        "  Delta_true > U"
    )
    print(
        "  means the true branch requires more K precision than"
    )
    print(
        "  the known-data uniqueness interval provides."
    )
    print()

    print("GLOBAL THRESHOLD TABLE")
    print()
    print(
        "  inst    U (unique)        Delta_true          gap"
    )
    print("  " + "-" * 78)

    for result in results:
        print(
            f"  {result['inst'].idx:>4} "
            f"{result['uniqueness_width']:>15} "
            f"{result['delta_true']:>20} "
            f"{result['gap']:>18}"
        )

    print()
    print("IMPORTANT INTERPRETATION")
    print()
    print(
        "  U is a genuinely known-data quantity."
    )
    print(
        "  Delta_true is a post-hoc diagnostic quantity."
    )
    print()
    print(
        "  Therefore a small U does NOT imply that the true d0"
    )
    print(
        "  can be recovered: the independently derived K interval"
    )
    print(
        "  must actually be narrower than the corresponding"
    )
    print(
        "  true-entry distance."
    )
    print()

    print("RESEARCH CONSEQUENCE")
    print()
    print(
        "  The next mathematically useful object is no longer"
    )
    print(
        "  another local CRT test of K."
    )
    print()
    print(
        "  It is an independently derivable upper/lower bound on K"
    )
    print(
        "  whose width can be compared directly against U and"
    )
    print(
        "  Delta_true."
    )
    print()

    print("=" * 120)
    print("EXPERIMENT 452 FINAL STATUS")
    print("=" * 120)
    print(
        f"  KNOWN-DATA UNIQUENESS RADIUS = True"
    )
    print(
        f"  O(1) d0 POPULATION            = True"
    )
    print(
        f"  INTEGER-EXACT                 = True"
    )
    print(
        f"  GIANT K ENUMERATION           = False"
    )
    print(
        f"  GIANT d0 ENUMERATION          = False"
    )
    print(
        f"  CRT CARTESIAN PRODUCT         = False"
    )
    print(
        f"  K_TRUE USED IN BOUNDS         = False"
    )
    print(
        f"  d0_TRUE USED IN BOUNDS        = False"
    )
    print(
        f"  X_TRUE USED IN BOUNDS         = False"
    )
    print()
    print(
        "  CONCLUSION:"
    )
    print(
        "    The exact known-data K precision window that preserves"
    )
    print(
        "    uniqueness can be calculated independently, then"
    )
    print(
        "    compared against the post-hoc precision required by"
    )
    print(
        "    the true branch."
    )
    print("=" * 120)
    print("EXPERIMENT 452 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    main()
