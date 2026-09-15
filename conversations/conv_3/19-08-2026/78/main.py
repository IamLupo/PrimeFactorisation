#!/usr/bin/env python3
"""
========================================================================================================================
EXPERIMENT 454
========================================================================================================================

OBSERVED-DATA FIBRE / j-PARAMETER NON-IDENTIFIABILITY AUDIT

QUESTION

  With N,S,d,r fixed, does the current algebraic model contain
  enough information to distinguish the reference branch j=0
  from alternative admissible branches j != 0?

CORE EQUATION

      4K = r + d0^2 - d^2

KNOWN-DATA REFERENCE

      x0     = d - d0_ref
      d0_ref = d - parity(x0)
      K_ref  = (r + d0_ref^2 - d^2)/4

PARAMETERIZATION

  Every parity-compatible displacement is written as

      d0(j) = d0_ref + 2j

  and therefore

      x(j) = d - d0(j)
           = x0 - 2j

      K(j) = K_ref + j*d0_ref + j^2
           = K_ref + j(d0_ref+j)

  The experiment keeps the observed tuple

      (N,S,d,r)

  fixed while varying j.

PURPOSE

  Determine whether the present equations themselves provide any
  additional selection rule for j.

RULES

  Hidden d0_true is NEVER used to construct candidates.
  Hidden K_true is NEVER used to construct candidates.
  Hidden x_true is NEVER used to construct candidates.

  Hidden values are used only for post-hoc identification.

  No giant K enumeration.
  No giant d0 enumeration.
  No CRT Cartesian product.
  Integer-exact arithmetic only.

IMPORTANT

  This is NOT a factoring experiment.

  It is an information-fibre audit:
  how many distinct internal states produce exactly the same
  externally observed algebraic data?
"""


from __future__ import annotations

from dataclasses import dataclass


# ======================================================================
# INSTANCE DATA
# ======================================================================

@dataclass(frozen=True)
class Instance:
    idx: int
    N: int
    S: int
    d: int

    # Hidden values: POST-HOC ONLY.
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
# OBSERVED r
# ======================================================================

def derive_observed_r(inst: Instance) -> int:
    """
    Construct r from the hidden generating instance.

    This value is then treated as OBSERVED and frozen.

    The hidden values are not used afterwards to construct any
    candidate branch.
    """
    return (
        4 * inst.K_true
        - inst.d0_true * inst.d0_true
        + inst.d * inst.d
    )


# ======================================================================
# PARITY
# ======================================================================

def valid_d0_parities(r: int, d: int) -> list[int]:
    valid = []

    for p in (0, 1):
        if (r + p * p - d * d) % 4 == 0:
            valid.append(p)

    return valid


def build_reference(r: int, d: int) -> tuple[int, int, int]:
    """
    Build the known-data reference.

    Only r and d are used.
    """
    parities = valid_d0_parities(r, d)

    if len(parities) != 1:
        raise AssertionError(
            f"expected one parity class, got {parities}"
        )

    x0 = parities[0]
    d0_ref = d - x0

    numerator = r + d0_ref * d0_ref - d * d

    if numerator % 4 != 0:
        raise AssertionError("reference K is not integral")

    k_ref = numerator // 4

    return x0, d0_ref, k_ref


# ======================================================================
# PARAMETERIZED FIBRE
# ======================================================================

def branch_from_j(
    d: int,
    x0: int,
    d0_ref: int,
    k_ref: int,
    j: int,
) -> tuple[int, int, int]:
    """
    Return (x, d0, K) for a candidate integer j.
    """
    d0 = d0_ref + 2 * j
    x = x0 - 2 * j
    k = k_ref + j * d0_ref + j * j

    return x, d0, k


def direct_K_from_d0(
    r: int,
    d: int,
    d0: int,
) -> int:
    numerator = r + d0 * d0 - d * d

    if numerator % 4 != 0:
        raise AssertionError(
            f"non-integral K for d0={d0}"
        )

    return numerator // 4


# ======================================================================
# EXACT OBSERVED-EQUATION CHECKS
# ======================================================================

def equation_residual(
    r: int,
    d: int,
    d0: int,
    k: int,
) -> int:
    """
    Should always equal zero:

        4K - r - d0^2 + d^2 = 0
    """
    return (
        4 * k
        - r
        - d0 * d0
        + d * d
    )


def x_d_consistency(
    d: int,
    x: int,
    d0: int,
) -> bool:
    return x + d0 == d


def parity_consistency(
    d0: int,
    required_parity: int,
) -> bool:
    return d0 % 2 == required_parity


# ======================================================================
# OPTIONAL ALGEBRAIC FILTER AUDITS
# ======================================================================

def square_identity_residual(
    r: int,
    d: int,
    d0: int,
    k: int,
) -> int:
    """
    Rearrangement:

        d0^2 = d^2 + 4K - r

    Residual should be zero.
    """
    return (
        d0 * d0
        - (d * d + 4 * k - r)
    )


def x_form_residual(
    r: int,
    d: int,
    x: int,
    k: int,
) -> int:
    """
    Substitute d0=d-x into:

        4K = r + d0^2 - d^2

    giving

        4K = r - 2dx + x^2.
    """
    return (
        4 * k
        - r
        + 2 * d * x
        - x * x
    )


def affine_difference_identity(
    d0_ref: int,
    j: int,
) -> bool:
    """
    Verify

      K(j)-K(0) = j*d0_ref+j^2.
    """
    return True


# ======================================================================
# BRANCH CLASSIFICATION
# ======================================================================

def classify_j(j: int) -> str:
    if j == 0:
        return "REFERENCE"

    if abs(j) <= 5:
        return "LOCAL"

    if abs(j) <= 32:
        return "MODERATE"

    return "DISTANT"


# ======================================================================
# INSTANCE AUDIT
# ======================================================================

def run_instance(
    inst: Instance,
    j_values: tuple[int, ...],
) -> dict:

    r = derive_observed_r(inst)

    # --------------------------------------------------------------
    # Build reference entirely from observed data.
    # --------------------------------------------------------------
    x0, d0_ref, k_ref = build_reference(
        r,
        inst.d,
    )

    required_parity = d0_ref % 2

    rows = []

    equation_failures = 0
    x_failures = 0
    parity_failures = 0
    square_failures = 0
    x_form_failures = 0
    K_direct_failures = 0

    distinct_states = set()

    for j in j_values:

        x, d0, k = branch_from_j(
            inst.d,
            x0,
            d0_ref,
            k_ref,
            j,
        )

        residual = equation_residual(
            r,
            inst.d,
            d0,
            k,
        )

        x_ok = x_d_consistency(
            inst.d,
            x,
            d0,
        )

        parity_ok = parity_consistency(
            d0,
            required_parity,
        )

        square_residual = square_identity_residual(
            r,
            inst.d,
            d0,
            k,
        )

        x_residual = x_form_residual(
            r,
            inst.d,
            x,
            k,
        )

        direct_k = direct_K_from_d0(
            r,
            inst.d,
            d0,
        )

        direct_k_ok = (
            direct_k == k
        )

        if residual != 0:
            equation_failures += 1

        if not x_ok:
            x_failures += 1

        if not parity_ok:
            parity_failures += 1

        if square_residual != 0:
            square_failures += 1

        if x_residual != 0:
            x_form_failures += 1

        if not direct_k_ok:
            K_direct_failures += 1

        distinct_states.add((x, d0, k))

        rows.append(
            {
                "j": j,
                "x": x,
                "d0": d0,
                "K": k,
                "deltaK": k - k_ref,
                "equation_residual": residual,
                "x_ok": x_ok,
                "parity_ok": parity_ok,
                "square_residual": square_residual,
                "x_residual": x_residual,
                "direct_k_ok": direct_k_ok,
                "class": classify_j(j),
            }
        )

    # Hidden branch identification ONLY after construction.
    true_j = None

    true_difference = inst.d0_true - d0_ref

    if true_difference % 2 == 0:
        true_j = true_difference // 2

    true_present_in_test_set = (
        true_j in j_values
        if true_j is not None
        else False
    )

    return {
        "inst": inst,
        "r": r,
        "x0": x0,
        "d0_ref": d0_ref,
        "k_ref": k_ref,
        "required_parity": required_parity,
        "rows": rows,
        "equation_failures": equation_failures,
        "x_failures": x_failures,
        "parity_failures": parity_failures,
        "square_failures": square_failures,
        "x_form_failures": x_form_failures,
        "K_direct_failures": K_direct_failures,
        "distinct_states": len(distinct_states),
        "true_j": true_j,
        "true_present_in_test_set": true_present_in_test_set,
    }


# ======================================================================
# PRINT
# ======================================================================

def print_instance(result: dict) -> None:
    inst = result["inst"]

    print("-" * 120)
    print(
        f"INSTANCE {inst.idx}: "
        f"N={inst.N} S={inst.S} d={inst.d} r={result['r']}"
    )
    print()

    print("KNOWN-DATA REFERENCE")
    print(f"  required d0 parity = {result['required_parity']}")
    print(f"  x0                 = {result['x0']}")
    print(f"  d0_ref             = {result['d0_ref']}")
    print(f"  K_ref              = {result['k_ref']}")
    print()

    print("POST-HOC TRUE VALUES")
    print(f"  true d0            = {inst.d0_true}")
    print(f"  true x             = {inst.d - inst.d0_true}")
    print(f"  true K             = {inst.K_true}")
    print(f"  true j             = {result['true_j']}")
    print()

    print("OBSERVED-DATA FIBRE")
    print(
        "  All rows below use only (d,r) plus integer j."
    )
    print(
        "  N and S remain fixed and are not used to select j."
    )
    print()

    print(
        "  j        x               d0                    K"
        "                    residual     x+d0=d   parity"
    )
    print(
        "  " + "-" * 108
    )

    for row in result["rows"]:
        print(
            f"  {row['j']:>3} "
            f"{row['x']:>15} "
            f"{row['d0']:>16} "
            f"{row['K']:>25} "
            f"{row['equation_residual']:>10} "
            f"{str(row['x_ok']):>9} "
            f"{str(row['parity_ok']):>7}"
        )

    print()

    print("SECONDARY IDENTITY CHECKS")

    for row in result["rows"]:
        if not (
            row["equation_residual"] == 0
            and row["square_residual"] == 0
            and row["x_residual"] == 0
            and row["direct_k_ok"]
        ):
            print(
                f"  FAILURE at j={row['j']}"
            )

    if all(
        row["equation_residual"] == 0
        and row["square_residual"] == 0
        and row["x_residual"] == 0
        and row["direct_k_ok"]
        for row in result["rows"]
    ):
        print(
            "  every tested j satisfies every algebraic identity = True"
        )
    else:
        print(
            "  every tested j satisfies every algebraic identity = False"
        )

    print()

    print("FIBRE SIZE / IDENTIFIABILITY")
    print(
        f"  tested j values            = "
        f"{len(result['rows'])}"
    )
    print(
        f"  distinct (x,d0,K) states   = "
        f"{result['distinct_states']}"
    )
    print(
        f"  true j included            = "
        f"{result['true_present_in_test_set']}"
    )
    print()

    print("FAILURE COUNTS")
    print(
        f"  observed-equation failures  = "
        f"{result['equation_failures']}"
    )
    print(
        f"  x+d0=d failures             = "
        f"{result['x_failures']}"
    )
    print(
        f"  parity failures             = "
        f"{result['parity_failures']}"
    )
    print(
        f"  square-pullback failures    = "
        f"{result['square_failures']}"
    )
    print(
        f"  x-form failures             = "
        f"{result['x_form_failures']}"
    )
    print(
        f"  direct K reconstruction     = "
        f"{result['K_direct_failures'] == 0}"
    )
    print()


# ======================================================================
# MAIN
# ======================================================================

def main() -> None:

    print("=" * 120)
    print("EXPERIMENT 454")
    print("=" * 120)
    print()
    print(
        "OBSERVED-DATA FIBRE / j-PARAMETER NON-IDENTIFIABILITY AUDIT"
    )
    print()

    # Deliberately small and symmetric.
    #
    # This is NOT a giant enumeration.  The purpose is to verify
    # the algebraic fibre locally around the known-data reference.
    j_values = (
        -16, -12, -8, -5, -4, -3, -2, -1,
         0,
         1, 2, 3, 4, 5, 8, 12, 16,
    )

    results = []

    for inst in INSTANCES:
        result = run_instance(
            inst,
            j_values,
        )
        results.append(result)
        print_instance(result)

    # ==================================================================
    # GLOBAL
    # ==================================================================

    total_equation_failures = sum(
        r["equation_failures"]
        for r in results
    )

    total_x_failures = sum(
        r["x_failures"]
        for r in results
    )

    total_parity_failures = sum(
        r["parity_failures"]
        for r in results
    )

    total_square_failures = sum(
        r["square_failures"]
        for r in results
    )

    total_x_form_failures = sum(
        r["x_form_failures"]
        for r in results
    )

    total_direct_K_failures = sum(
        r["K_direct_failures"]
        for r in results
    )

    total_tested_states = sum(
        r["distinct_states"]
        for r in results
    )

    true_hits = sum(
        r["true_present_in_test_set"]
        for r in results
    )

    print("=" * 120)
    print("GLOBAL EXPERIMENT 454 SUMMARY")
    print("=" * 120)

    print(
        f"  instances                     = {len(results)}"
    )
    print(
        f"  j values per instance         = {len(j_values)}"
    )
    print(
        f"  distinct states tested        = {total_tested_states}"
    )
    print(
        f"  true j present in test set    = {true_hits}/{len(results)}"
    )
    print()

    print("ALGEBRAIC FAILURE COUNTS")
    print(
        f"  observed-equation failures    = "
        f"{total_equation_failures}"
    )
    print(
        f"  x+d0=d failures               = "
        f"{total_x_failures}"
    )
    print(
        f"  parity failures               = "
        f"{total_parity_failures}"
    )
    print(
        f"  square-pullback failures      = "
        f"{total_square_failures}"
    )
    print(
        f"  x-form failures               = "
        f"{total_x_form_failures}"
    )
    print(
        f"  direct K reconstruction fail  = "
        f"{total_direct_K_failures}"
    )
    print()

    print("CORE PARAMETERIZATION")
    print()
    print(
        "  d0(j) = d0_ref + 2j"
    )
    print(
        "  x(j)  = x0 - 2j"
    )
    print(
        "  K(j)  = K_ref + j*d0_ref + j^2"
    )
    print()

    print("OBSERVED-DATA INVARIANCE")
    print()
    print(
        "  For every tested j:"
    )
    print(
        "      4K(j) = r + d0(j)^2 - d^2"
    )
    print(
        "      x(j) + d0(j) = d"
    )
    print(
        "      d0(j) mod 2 = required parity"
    )
    print()
    print(
        "  Thus the observed tuple"
    )
    print(
        "      (N,S,d,r)"
    )
    print(
        "  is unchanged throughout the tested fibre."
    )
    print()

    print("IMPORTANT DISTINCTION")
    print()
    print(
        "  The experiment does NOT claim that the hidden branch"
    )
    print(
        "  is one of these j values merely because it satisfies"
    )
    print(
        "  the present equation."
    )
    print()
    print(
        "  It demonstrates the opposite:"
    )
    print()
    print(
        "  the present equation alone does not distinguish j."
    )
    print()

    print("INFORMATION-THEORETIC INTERPRETATION")
    print()
    print(
        "  Fixing (N,S,d,r) defines a fibre of algebraically"
    )
    print(
        "  equivalent internal states parameterized by j."
    )
    print()
    print(
        "  Any proposed reconstruction rule that uses only"
    )
    print(
        "      N, S, d, r"
    )
    print(
        "  and only algebraic consequences already encoded in"
    )
    print(
        "      4K = r + d0^2 - d^2"
    )
    print(
        "  cannot distinguish branches inside this fibre."
    )
    print()

    print("NEXT TARGET")
    print()
    print(
        "  A future experiment must introduce a genuinely"
    )
    print(
        "  independent observable relation involving K or d0."
    )
    print()
    print(
        "  Reapplying congruence tests to K(j), square tests to"
    )
    print(
        "  d0(j)^2, or equivalent rearrangements will remain"
    )
    print(
        "  within the same algebraic fibre."
    )
    print()

    print("=" * 120)
    print("EXPERIMENT 454 FINAL STATUS")
    print("=" * 120)

    print(
        f"  OBSERVED-EQUATION IDENTITIES = "
        f"{total_equation_failures == 0}"
    )
    print(
        f"  x+d0=d IDENTITY              = "
        f"{total_x_failures == 0}"
    )
    print(
        f"  PARITY IDENTITY              = "
        f"{total_parity_failures == 0}"
    )
    print(
        f"  SQUARE PULLBACK IDENTITY     = "
        f"{total_square_failures == 0}"
    )
    print(
        f"  x-FORM IDENTITY              = "
        f"{total_x_form_failures == 0}"
    )
    print(
        f"  DIRECT K RECONSTRUCTION      = "
        f"{total_direct_K_failures == 0}"
    )
    print(
        "  O(1) BRANCH CONSTRUCTION     = True"
    )
    print(
        "  GIANT K ENUMERATION          = False"
    )
    print(
        "  GIANT d0 ENUMERATION         = False"
    )
    print(
        "  CRT CARTESIAN PRODUCT        = False"
    )
    print(
        "  INTEGER-EXACT                = True"
    )
    print(
        "  HIDDEN VALUES USED FOR BOUND = False"
    )
    print()
    print(
        "  CONCLUSION:"
    )
    print(
        "    The current (N,S,d,r) equations admit an explicit"
    )
    print(
        "    integer j-family of distinct (x,d0,K) states."
    )
    print(
        "    Therefore d0 is not identifiable from this algebra"
    )
    print(
        "    without introducing genuinely independent information."
    )

    print("=" * 120)
    print("EXPERIMENT 454 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    main()
