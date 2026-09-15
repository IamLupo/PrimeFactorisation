#!/usr/bin/env python3
"""
========================================================================================================================
EXPERIMENT 453
========================================================================================================================

EXACT K-DISPLACEMENT / x-OFFSET STRUCTURE AUDIT

QUESTION

  Experiment 452 measured:

      U = largest symmetric K-width around K_ref
          containing exactly one admissible d0.

  The observed values strongly suggest the exact identity

      U = d0_ref - 2

  and the exact K displacement should be expressible directly
  from the d0/x displacement.

GOALS

  1. Prove/verify the exact nearest-neighbour spacing:

       K(d0+2) - K(d0) = d0 + 1
       K(d0) - K(d0-2) = d0 - 1

  2. Verify the uniqueness-radius identity:

       U = d0_ref - 2

     for every supplied instance.

  3. Let

       j = (d0_true - d0_ref) / 2

     because admissible d0 values differ by 2.

     Derive and verify

       K_true - K_ref
         = j*d0_ref + j^2

     exactly.

  4. Express the same displacement through x.

     Since

       x = d-d0

     and

       x0 = d-d0_ref,

     we have

       d0_true - d0_ref = x0 - x_true.

     Therefore

       j = (x0-x_true)/2.

  5. Determine whether Delta_true can be normalized by the
     x displacement and reference d0 without enumerating anything.

  6. Compare the exact true-entry radius with the exact
     known-data uniqueness radius.

RULES

  K_true is NEVER used to construct a bound.
  d0_true is NEVER used to construct a bound.
  x_true is NEVER used to construct a bound.

  Hidden values are used ONLY for post-hoc verification.

  No giant K enumeration.
  No giant d0 enumeration.
  No CRT Cartesian product.
  Integer-exact arithmetic only.
"""


from __future__ import annotations

from dataclasses import dataclass


# ======================================================================
# DATA
# ======================================================================

@dataclass(frozen=True)
class Instance:
    idx: int
    N: int
    S: int
    d: int

    # Hidden/post-hoc values.
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
# CORE ALGEBRA
# ======================================================================

def derive_r(inst: Instance) -> int:
    """
    Recover the synthetic observed r.

    This is instance-generation consistency only.
    It is not used to construct hidden-data bounds.
    """
    return (
        4 * inst.K_true
        - inst.d0_true * inst.d0_true
        + inst.d * inst.d
    )


def valid_parities(r: int, d: int) -> list[int]:
    out: list[int] = []

    for p in (0, 1):
        if (r + p * p - d * d) % 4 == 0:
            out.append(p)

    return out


def reference_data(
    r: int,
    d: int,
) -> tuple[int, int, int]:
    """
    Construct reference entirely from known N,S,d,r.

    x0 is the unique admissible parity.
    d0_ref = d-x0.
    """
    parities = valid_parities(r, d)

    if len(parities) != 1:
        raise AssertionError(
            f"expected exactly one parity class, got {parities}"
        )

    x0 = parities[0]
    d0_ref = d - x0

    numerator = r + d0_ref * d0_ref - d * d

    if numerator % 4 != 0:
        raise AssertionError("K_ref not integral")

    k_ref = numerator // 4

    return x0, d0_ref, k_ref


def k_of_d0(
    r: int,
    d: int,
    d0: int,
) -> int:
    num = r + d0 * d0 - d * d

    if num % 4 != 0:
        raise AssertionError(
            f"non-integral K for d0={d0}"
        )

    return num // 4


# ======================================================================
# EXACT SPACING IDENTITIES
# ======================================================================

def predicted_up_spacing(d0: int) -> int:
    """
    K(d0+2)-K(d0) = d0+1.
    """
    return d0 + 1


def predicted_down_spacing(d0: int) -> int:
    """
    K(d0)-K(d0-2) = d0-1.
    """
    return d0 - 1


def verify_spacing(
    r: int,
    d: int,
    d0: int,
) -> tuple[int, int, bool, int, int, bool]:
    k0 = k_of_d0(r, d, d0)

    k_up = k_of_d0(r, d, d0 + 2)
    up_exact = k_up - k0
    up_pred = predicted_up_spacing(d0)

    if d0 >= 2:
        k_down = k_of_d0(r, d, d0 - 2)
        down_exact = k0 - k_down
        down_pred = predicted_down_spacing(d0)
    else:
        down_exact = 0
        down_pred = 0

    return (
        up_exact,
        up_pred,
        up_exact == up_pred,
        down_exact,
        down_pred,
        down_exact == down_pred,
    )


# ======================================================================
# UNIQUENESS RADIUS
# ======================================================================

def exact_uniqueness_radius_from_reference(
    d0_ref: int,
) -> int:
    """
    The nearest lower admissible d0 is d0_ref-2.

    Its K distance from K_ref is

        d0_ref - 1.

    Therefore the largest integer width which still excludes
    that competitor is

        U = d0_ref - 2.
    """
    if d0_ref < 2:
        raise AssertionError(
            "reference d0 too small for this construction"
        )

    return d0_ref - 2


# ======================================================================
# K DISPLACEMENT IN TERMS OF d0 OFFSET
# ======================================================================

def displacement_parameter(
    d0_true: int,
    d0_ref: int,
) -> int:
    diff = d0_true - d0_ref

    if diff % 2 != 0:
        raise AssertionError(
            "true d0 and reference d0 are not in same parity class"
        )

    return diff // 2


def predicted_k_delta(
    d0_ref: int,
    j: int,
) -> int:
    """
    Exact identity:

        K(d0_ref + 2j) - K(d0_ref)
          = j*d0_ref + j^2
    """
    return j * d0_ref + j * j


def alternative_factored_delta(
    d0_ref: int,
    j: int,
) -> int:
    """
    Same quantity:

        j(d0_ref+j)
    """
    return j * (d0_ref + j)


# ======================================================================
# x-SPACE FORM
# ======================================================================

def x_offset(
    d: int,
    d0_ref: int,
    d0_true: int,
) -> tuple[int, int, int]:
    x0 = d - d0_ref
    x_true = d - d0_true

    diff = x0 - x_true

    if diff % 2 != 0:
        raise AssertionError(
            "x displacement parity mismatch"
        )

    j = diff // 2

    return x0, x_true, j


def predicted_k_delta_from_x(
    d0_ref: int,
    x0: int,
    x_true: int,
) -> int:
    """
    Since

        j = (x0-x_true)/2,

    substitute directly.
    """
    dx = x0 - x_true

    if dx % 2 != 0:
        raise AssertionError(
            "x displacement must be even"
        )

    j = dx // 2

    return j * d0_ref + j * j


# ======================================================================
# TRUE ENTRY RADIUS
# ======================================================================

def true_entry_radius(
    k_true: int,
    k_ref: int,
) -> int:
    return abs(k_true - k_ref)


# ======================================================================
# NORMALIZED STRUCTURE
# ======================================================================

def derive_normalized_relation(
    delta_k: int,
    j: int,
    d0_ref: int,
) -> tuple[bool, int, int]:
    """
    Verify

        DeltaK = j(d0_ref+j)

    and expose

        DeltaK / j = d0_ref+j

    when j != 0.
    """
    pred = j * (d0_ref + j)

    if j == 0:
        return delta_k == 0, 0, d0_ref

    divisible = (delta_k % j == 0)

    if not divisible:
        return False, 0, 0

    quotient = delta_k // j

    return (
        delta_k == pred,
        quotient,
        d0_ref + j,
    )


# ======================================================================
# MAIN INSTANCE
# ======================================================================

def run_instance(inst: Instance) -> dict:
    r = derive_r(inst)

    x0, d0_ref, k_ref = reference_data(
        r,
        inst.d,
    )

    # Exact reference relation.
    if inst.d - d0_ref != x0:
        raise AssertionError(
            f"x0 relation failed: instance={inst.idx}"
        )

    # Verify hidden K reconstruction.
    k_true_reconstructed = k_of_d0(
        r,
        inst.d,
        inst.d0_true,
    )

    if k_true_reconstructed != inst.K_true:
        raise AssertionError(
            f"hidden K reconstruction failed: instance={inst.idx}"
        )

    # --------------------------------------------------------------
    # LOCAL SPACING
    # --------------------------------------------------------------
    (
        up_exact,
        up_pred,
        up_ok,
        down_exact,
        down_pred,
        down_ok,
    ) = verify_spacing(
        r,
        inst.d,
        d0_ref,
    )

    # --------------------------------------------------------------
    # UNIQUENESS RADIUS
    # --------------------------------------------------------------
    U = exact_uniqueness_radius_from_reference(
        d0_ref
    )

    nearest_lower_d0 = d0_ref - 2
    nearest_lower_k = k_of_d0(
        r,
        inst.d,
        nearest_lower_d0,
    )

    lower_distance = abs(
        nearest_lower_k - k_ref
    )

    expected_lower_distance = d0_ref - 1

    if lower_distance != expected_lower_distance:
        raise AssertionError(
            f"lower spacing identity failed: instance={inst.idx}"
        )

    # --------------------------------------------------------------
    # TRUE DISPLACEMENT
    # --------------------------------------------------------------
    j = displacement_parameter(
        inst.d0_true,
        d0_ref,
    )

    delta_k_signed = inst.K_true - k_ref
    delta_k_abs = abs(delta_k_signed)

    predicted_delta = predicted_k_delta(
        d0_ref,
        j,
    )

    predicted_delta_factored = alternative_factored_delta(
        d0_ref,
        j,
    )

    delta_formula_ok = (
        delta_k_signed == predicted_delta
        and predicted_delta == predicted_delta_factored
    )

    # --------------------------------------------------------------
    # x form
    # --------------------------------------------------------------
    x0_check, x_true, j_from_x = x_offset(
        inst.d,
        d0_ref,
        inst.d0_true,
    )

    if x0_check != x0:
        raise AssertionError(
            f"x offset reference mismatch: instance={inst.idx}"
        )

    x_formula_delta = predicted_k_delta_from_x(
        d0_ref,
        x0,
        x_true,
    )

    x_formula_ok = (
        delta_k_signed == x_formula_delta
        and j == j_from_x
    )

    # --------------------------------------------------------------
    # NORMALIZED RELATION
    # --------------------------------------------------------------
    (
        normalized_ok,
        delta_over_j,
        expected_delta_over_j,
    ) = derive_normalized_relation(
        delta_k_signed,
        j,
        d0_ref,
    )

    # --------------------------------------------------------------
    # GAP
    # --------------------------------------------------------------
    gap = delta_k_abs - U

    if gap <= 0:
        # It should not occur in the supplied construction, but keep
        # this as an explicit diagnostic.
        relative_status = "TRUE_WITHIN_UNIQUENESS"
    else:
        relative_status = "TRUE_OUTSIDE_UNIQUENESS"

    return {
        "inst": inst,
        "r": r,
        "x0": x0,
        "x_true": x_true,
        "d0_ref": d0_ref,
        "k_ref": k_ref,
        "up_exact": up_exact,
        "up_pred": up_pred,
        "up_ok": up_ok,
        "down_exact": down_exact,
        "down_pred": down_pred,
        "down_ok": down_ok,
        "U": U,
        "nearest_lower_d0": nearest_lower_d0,
        "nearest_lower_k": nearest_lower_k,
        "lower_distance": lower_distance,
        "expected_lower_distance": expected_lower_distance,
        "j": j,
        "delta_k_signed": delta_k_signed,
        "delta_k_abs": delta_k_abs,
        "predicted_delta": predicted_delta,
        "predicted_delta_factored": predicted_delta_factored,
        "delta_formula_ok": delta_formula_ok,
        "j_from_x": j_from_x,
        "x_formula_delta": x_formula_delta,
        "x_formula_ok": x_formula_ok,
        "normalized_ok": normalized_ok,
        "delta_over_j": delta_over_j,
        "expected_delta_over_j": expected_delta_over_j,
        "gap": gap,
        "relative_status": relative_status,
    }


# ======================================================================
# PRINT
# ======================================================================

def print_instance(res: dict) -> None:
    inst: Instance = res["inst"]

    print("-" * 120)
    print(
        f"INSTANCE {inst.idx}: "
        f"N={inst.N} S={inst.S} d={inst.d} r={res['r']}"
    )
    print()

    print("KNOWN-DATA REFERENCE")
    print(f"  x0     = {res['x0']}")
    print(f"  d0_ref = {res['d0_ref']}")
    print(f"  K_ref  = {res['k_ref']}")
    print()

    print("POST-HOC TRUE VALUES")
    print(f"  true d0 = {inst.d0_true}")
    print(f"  true x  = {res['x_true']}")
    print(f"  true K  = {inst.K_true}")
    print()

    print("LOCAL SPACING")
    print(
        f"  K(d0+2)-K(d0) exact      = {res['up_exact']}"
    )
    print(
        f"  K(d0+2)-K(d0) predicted  = {res['up_pred']}"
    )
    print(
        f"  upward identity          = {res['up_ok']}"
    )
    print(
        f"  K(d0)-K(d0-2) exact      = {res['down_exact']}"
    )
    print(
        f"  K(d0)-K(d0-2) predicted  = {res['down_pred']}"
    )
    print(
        f"  downward identity        = {res['down_ok']}"
    )
    print()

    print("UNIQUENESS-RADIUS IDENTITY")
    print(
        f"  nearest lower d0         = {res['nearest_lower_d0']}"
    )
    print(
        f"  nearest lower K          = {res['nearest_lower_k']}"
    )
    print(
        f"  K-distance               = {res['lower_distance']}"
    )
    print(
        f"  predicted d0_ref-1       = "
        f"{res['expected_lower_distance']}"
    )
    print(
        f"  spacing identity         = "
        f"{res['lower_distance'] == res['expected_lower_distance']}"
    )
    print(
        f"  uniqueness radius U      = {res['U']}"
    )
    print(
        f"  predicted U=d0_ref-2    = "
        f"{res['U'] == res['d0_ref'] - 2}"
    )
    print()

    print("EXACT d0 DISPLACEMENT")
    print(
        f"  j = (d0_true-d0_ref)/2   = {res['j']}"
    )
    print(
        f"  x0-x_true               = "
        f"{res['x0'] - res['x_true']}"
    )
    print(
        f"  j from x                 = {res['j_from_x']}"
    )
    print(
        f"  j consistency            = "
        f"{res['j'] == res['j_from_x']}"
    )
    print()

    print("K DISPLACEMENT IDENTITY")
    print(
        f"  DeltaK signed             = "
        f"{res['delta_k_signed']}"
    )
    print(
        f"  |DeltaK|                 = "
        f"{res['delta_k_abs']}"
    )
    print(
        f"  j*d0_ref+j^2             = "
        f"{res['predicted_delta']}"
    )
    print(
        f"  j*(d0_ref+j)             = "
        f"{res['predicted_delta_factored']}"
    )
    print(
        f"  exact d0-offset identity = "
        f"{res['delta_formula_ok']}"
    )
    print()

    print("x-SPACE PULLBACK OF DELTA K")
    print(
        "  DeltaK = ((x0-x_true)/2) * "
        "(d0_ref + (x0-x_true)/2)"
    )
    print(
        f"  x-derived DeltaK         = "
        f"{res['x_formula_delta']}"
    )
    print(
        f"  x-derived identity       = "
        f"{res['x_formula_ok']}"
    )
    print()

    print("NORMALIZED FORM")
    print(
        f"  DeltaK / j               = "
        f"{res['delta_over_j']}"
    )
    print(
        f"  d0_ref + j               = "
        f"{res['expected_delta_over_j']}"
    )
    print(
        f"  normalized identity      = "
        f"{res['normalized_ok']}"
    )
    print()

    print("TRUE-ENTRY VS UNIQUENESS")
    print(
        f"  uniqueness radius U      = {res['U']}"
    )
    print(
        f"  true-entry radius        = {res['delta_k_abs']}"
    )
    print(
        f"  precision gap            = "
        f"{res['gap']}"
    )
    print(
        f"  status                   = "
        f"{res['relative_status']}"
    )
    print()


# ======================================================================
# MAIN
# ======================================================================

def main() -> None:
    print("=" * 120)
    print("EXPERIMENT 453")
    print("=" * 120)
    print()
    print("EXACT K-DISPLACEMENT / x-OFFSET STRUCTURE AUDIT")
    print()
    print("QUESTIONS")
    print("  Is U exactly d0_ref-2?")
    print("  Is DeltaK exactly j*d0_ref+j^2?")
    print("  Can DeltaK be expressed directly through x0-x_true?")
    print()

    results = []

    for inst in INSTANCES:
        result = run_instance(inst)
        results.append(result)
        print_instance(result)

    # ==================================================================
    # GLOBAL CHECKS
    # ==================================================================

    spacing_up_failures = sum(
        not r["up_ok"]
        for r in results
    )

    spacing_down_failures = sum(
        not r["down_ok"]
        for r in results
    )

    uniqueness_failures = sum(
        not (r["U"] == r["d0_ref"] - 2)
        for r in results
    )

    delta_formula_failures = sum(
        not r["delta_formula_ok"]
        for r in results
    )

    x_formula_failures = sum(
        not r["x_formula_ok"]
        for r in results
    )

    normalized_failures = sum(
        not r["normalized_ok"]
        for r in results
    )

    true_outside = sum(
        r["gap"] > 0
        for r in results
    )

    print("=" * 120)
    print("GLOBAL EXPERIMENT 453 SUMMARY")
    print("=" * 120)

    print(
        f"  instances                         = {len(results)}"
    )
    print(
        f"  upward spacing failures           = "
        f"{spacing_up_failures}"
    )
    print(
        f"  downward spacing failures         = "
        f"{spacing_down_failures}"
    )
    print(
        f"  U=d0_ref-2 failures               = "
        f"{uniqueness_failures}"
    )
    print(
        f"  exact d0-offset delta failures    = "
        f"{delta_formula_failures}"
    )
    print(
        f"  exact x-offset delta failures     = "
        f"{x_formula_failures}"
    )
    print(
        f"  normalized identity failures      = "
        f"{normalized_failures}"
    )
    print(
        f"  true branch outside U              = "
        f"{true_outside}/{len(results)}"
    )
    print()

    print("EXACT STRUCTURAL IDENTITIES")
    print()
    print(
        "  For every parity-compatible d0:"
    )
    print()
    print(
        "      K(d0+2)-K(d0) = d0+1"
    )
    print(
        "      K(d0)-K(d0-2) = d0-1"
    )
    print()
    print(
        "  Therefore, for the canonical reference:"
    )
    print()
    print(
        "      U = d0_ref - 2"
    )
    print()
    print(
        "  This is an exact algebraic identity, not an empirical"
    )
    print(
        "  property of the tested instances."
    )
    print()

    print("EXACT GLOBAL DISPLACEMENT FORM")
    print()
    print(
        "  Let"
    )
    print(
        "      j = (d0_true-d0_ref)/2."
    )
    print()
    print(
        "  Then"
    )
    print(
        "      K_true-K_ref = j*d0_ref + j^2"
    )
    print(
        "                    = j(d0_ref+j)."
    )
    print()
    print(
        "  Since"
    )
    print(
        "      j = (x0-x_true)/2,"
    )
    print()
    print(
        "  the same displacement is"
    )
    print(
        "      K_true-K_ref"
    )
    print(
        "        = ((x0-x_true)/2)"
    )
    print(
        "          * (d0_ref + (x0-x_true)/2)."
    )
    print()

    print("GLOBAL TABLE")
    print()
    print(
        "  inst   d0_ref          U          |DeltaK|       "
        "j        x0-x_true"
    )
    print("  " + "-" * 100)

    for r in results:
        print(
            f"  {r['inst'].idx:>4} "
            f"{r['d0_ref']:>14} "
            f"{r['U']:>12} "
            f"{r['delta_k_abs']:>18} "
            f"{r['j']:>9} "
            f"{r['x0'] - r['x_true']:>11}"
        )

    print()
    print("INTERPRETATION")
    print()
    print(
        "  Experiment 452's uniqueness radius is structurally fixed"
    )
    print(
        "  by the local spacing at the known-data reference."
    )
    print()
    print(
        "  It does not provide an independent estimate of the hidden"
    )
    print(
        "  d0. It only describes how much K uncertainty can be tolerated"
    )
    print(
        "  before a neighboring admissible d0 enters."
    )
    print()
    print(
        "  The true-entry distance, in contrast, depends on the hidden"
    )
    print(
        "  displacement j and is therefore a post-hoc diagnostic."
    )
    print()
    print(
        "  The next useful question is whether some independently"
    )
    print(
        "  observable quantity can determine or bound j itself."
    )
    print()

    print("=" * 120)
    print("EXPERIMENT 453 FINAL STATUS")
    print("=" * 120)
    print(
        f"  EXACT SPACING IDENTITIES        = "
        f"{spacing_up_failures == 0 and spacing_down_failures == 0}"
    )
    print(
        f"  U = d0_ref - 2                  = "
        f"{uniqueness_failures == 0}"
    )
    print(
        f"  EXACT d0-OFFSET DELTA           = "
        f"{delta_formula_failures == 0}"
    )
    print(
        f"  EXACT x-OFFSET DELTA            = "
        f"{x_formula_failures == 0}"
    )
    print(
        f"  NORMALIZED IDENTITY              = "
        f"{normalized_failures == 0}"
    )
    print(
        "  O(1) ARITHMETIC                 = True"
    )
    print(
        "  GIANT K ENUMERATION             = False"
    )
    print(
        "  GIANT d0 ENUMERATION            = False"
    )
    print(
        "  CRT CARTESIAN PRODUCT           = False"
    )
    print(
        "  INTEGER-EXACT                    = True"
    )
    print(
        "  HIDDEN VALUES USED FOR BOUNDS   = False"
    )
    print()
    print(
        "  CONCLUSION:"
    )
    print(
        "    The K-precision gap measured in Experiment 452"
    )
    print(
        "    has an exact closed-form displacement structure."
    )
    print("=" * 120)
    print("EXPERIMENT 453 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    main()
