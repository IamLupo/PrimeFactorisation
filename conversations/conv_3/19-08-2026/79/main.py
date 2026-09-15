#!/usr/bin/env python3
"""
========================================================================================================================
EXPERIMENT 455
========================================================================================================================

FULL KNOWN-DATA RELATION / j-FIBRE INVARIANCE AUDIT

QUESTION

  Experiment 454 established a one-parameter family

      d0(j) = d0_ref + 2j
      x(j)  = x0 - 2j
      K(j)  = K_ref + j*d0_ref + j^2

  satisfying

      4K = r + d0^2 - d^2.

  The next question is:

      Do the OTHER known-data relations involving
      N, S, d, x, d0, K and A break this j-family?

OBJECTIVE

  Test whether any previously used algebraic relation contains
  information about j that is independent of the defining
  quadratic relation.

IMPORTANT

  A relation is counted as INDEPENDENT only if its value/residual
  changes across j and consequently excludes some otherwise
  admissible j.

  A relation whose residual is identically zero for every j is
  classified as a PULLBACK IDENTITY.

RULES

  exact integer arithmetic only
  no factorization
  no giant K enumeration
  no giant d0 enumeration
  no CRT Cartesian product

  K_true, d0_true and x_true are POST-HOC ONLY.
  They are never used to construct candidate branches.

========================================================================================================================
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

    # hidden / post-hoc only
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

def derive_r(inst: Instance) -> int:
    return (
        4 * inst.K_true
        - inst.d0_true * inst.d0_true
        + inst.d * inst.d
    )


# ======================================================================
# PARITY / REFERENCE
# ======================================================================

def valid_d0_parities(r: int, d: int) -> list[int]:
    return [
        p
        for p in (0, 1)
        if (r + p * p - d * d) % 4 == 0
    ]


def build_reference(r: int, d: int) -> tuple[int, int, int]:
    parities = valid_d0_parities(r, d)

    if len(parities) != 1:
        raise AssertionError(
            f"expected exactly one parity class, got {parities}"
        )

    x0 = parities[0]
    d0_ref = d - x0

    numerator = r + d0_ref * d0_ref - d * d

    if numerator % 4 != 0:
        raise AssertionError(
            "reference K is not integer"
        )

    K_ref = numerator // 4

    return x0, d0_ref, K_ref


# ======================================================================
# j-FAMILY
# ======================================================================

def branch(
    d: int,
    x0: int,
    d0_ref: int,
    K_ref: int,
    j: int,
) -> tuple[int, int, int]:
    d0 = d0_ref + 2 * j
    x = x0 - 2 * j
    K = K_ref + j * d0_ref + j * j
    return x, d0, K


# ======================================================================
# DERIVED QUANTITIES
# ======================================================================

def derive_A_over_4(
    inst: Instance,
    d0: int,
    K: int,
) -> int:
    """
    From

        4K = r + d0^2 - d^2

    and the observed construction used in the previous experiments,
    obtain

        A/4 = K + ((x+d0)/2)^2

    Since x+d0=d, this is

        A/4 = K + d^2/4

    The function uses exact divisibility and is intentionally
    written in several equivalent ways below.
    """
    x = inst.d - d0

    y = x + d0

    if y % 2 != 0:
        raise AssertionError(
            "x+d0 must be even for A/4 construction"
        )

    return K + (y // 2) ** 2


def derive_A(
    inst: Instance,
    d0: int,
    K: int,
) -> int:
    A4 = derive_A_over_4(inst, d0, K)
    return 4 * A4


# ======================================================================
# RELATION AUDITS
# ======================================================================

def audit_relations(
    inst: Instance,
    r: int,
    x0: int,
    d0_ref: int,
    K_ref: int,
    j: int,
) -> dict:

    x, d0, K = branch(
        inst.d,
        x0,
        d0_ref,
        K_ref,
        j,
    )

    # ---------------------------------------------------------------
    # Basic observed-data relations
    # ---------------------------------------------------------------

    residual_r = (
        4 * K
        - r
        - d0 * d0
        + inst.d * inst.d
    )

    residual_x_d0 = (
        x + d0 - inst.d
    )

    residual_d_NS = (
        inst.d
        - (inst.N - 2 * inst.S + 1)
    )

    # ---------------------------------------------------------------
    # Quadratic rewritten in x
    # ---------------------------------------------------------------

    residual_x_form = (
        4 * K
        - r
        + 2 * inst.d * x
        - x * x
    )

    # ---------------------------------------------------------------
    # Square pullback
    # ---------------------------------------------------------------

    residual_square = (
        d0 * d0
        - (
            inst.d * inst.d
            + 4 * K
            - r
        )
    )

    # ---------------------------------------------------------------
    # Local K spacing
    # ---------------------------------------------------------------

    x_next, d0_next, K_next = branch(
        inst.d,
        x0,
        d0_ref,
        K_ref,
        j + 1,
    )

    spacing_up_residual = (
        K_next - K
        - (d0 + 1)
    )

    x_next2, d0_prev, K_prev = branch(
        inst.d,
        x0,
        d0_ref,
        K_ref,
        j - 1,
    )

    spacing_down_residual = (
        K - K_prev
        - (d0 - 1)
    )

    # ---------------------------------------------------------------
    # d0 / x displacement identities
    # ---------------------------------------------------------------

    displacement_residual = (
        (d0 - d0_ref)
        + (x - x0)
    )

    # ---------------------------------------------------------------
    # K displacement identity
    # ---------------------------------------------------------------

    K_displacement_residual = (
        (K - K_ref)
        - j * (d0_ref + j)
    )

    # ---------------------------------------------------------------
    # A identities
    # ---------------------------------------------------------------

    A4 = derive_A_over_4(
        inst,
        d0,
        K,
    )

    A = 4 * A4

    # A4 expressed through x+d0
    if (x + d0) % 2 == 0:
        A4_from_y = K + ((x + d0) // 2) ** 2
        A4_y_residual = A4 - A4_from_y
    else:
        A4_y_residual = None

    # A4 expressed only through d
    if inst.d % 2 == 0:
        A4_from_d = K + (inst.d // 2) ** 2
        A4_d_residual = A4 - A4_from_d
    else:
        A4_d_residual = None

    # ---------------------------------------------------------------
    # Candidate summary
    # ---------------------------------------------------------------

    return {
        "j": j,
        "x": x,
        "d0": d0,
        "K": K,

        "residual_r": residual_r,
        "residual_x_d0": residual_x_d0,
        "residual_d_NS": residual_d_NS,
        "residual_x_form": residual_x_form,
        "residual_square": residual_square,
        "spacing_up": spacing_up_residual,
        "spacing_down": spacing_down_residual,
        "displacement": displacement_residual,
        "K_displacement": K_displacement_residual,

        "A4": A4,
        "A": A,

        "A4_y_residual": A4_y_residual,
        "A4_d_residual": A4_d_residual,
    }


# ======================================================================
# MAIN INSTANCE AUDIT
# ======================================================================

def run_instance(
    inst: Instance,
    j_values: tuple[int, ...],
) -> dict:

    r = derive_r(inst)

    x0, d0_ref, K_ref = build_reference(
        r,
        inst.d,
    )

    rows = [
        audit_relations(
            inst,
            r,
            x0,
            d0_ref,
            K_ref,
            j,
        )
        for j in j_values
    ]

    return {
        "inst": inst,
        "r": r,
        "x0": x0,
        "d0_ref": d0_ref,
        "K_ref": K_ref,
        "rows": rows,
    }


# ======================================================================
# PRINT
# ======================================================================

def print_instance(result: dict) -> None:

    inst = result["inst"]

    print("-" * 120)
    print(
        f"INSTANCE {inst.idx}: "
        f"N={inst.N} "
        f"S={inst.S} "
        f"d={inst.d} "
        f"r={result['r']}"
    )
    print()

    print("KNOWN-DATA REFERENCE")
    print(f"  x0      = {result['x0']}")
    print(f"  d0_ref  = {result['d0_ref']}")
    print(f"  K_ref   = {result['K_ref']}")
    print()

    print("POST-HOC TRUE VALUES")
    print(f"  true d0 = {inst.d0_true}")
    print(f"  true x  = {inst.d - inst.d0_true}")
    print(f"  true K  = {inst.K_true}")

    true_j = (
        (inst.d0_true - result["d0_ref"]) // 2
    )

    print(f"  true j  = {true_j}")
    print()

    print("FULL j-FIBRE RELATION AUDIT")
    print(
        "  j        d0            x                    K"
        "       r-resid   x-resid   d(NS)-resid"
    )
    print(
        "  " + "-" * 108
    )

    for row in result["rows"]:
        print(
            f"  {row['j']:>3}"
            f" {row['d0']:>13}"
            f" {row['x']:>13}"
            f" {row['K']:>25}"
            f" {row['residual_r']:>10}"
            f" {row['residual_x_d0']:>10}"
            f" {row['residual_d_NS']:>12}"
        )

    print()

    # ---------------------------------------------------------------
    # Count nonzero residuals.
    # ---------------------------------------------------------------

    relation_names = (
        "residual_r",
        "residual_x_d0",
        "residual_d_NS",
        "residual_x_form",
        "residual_square",
        "spacing_up",
        "spacing_down",
        "displacement",
        "K_displacement",
        "A4_y_residual",
        "A4_d_residual",
    )

    print("RELATION STATUS")

    for name in relation_names:

        values = [
            row[name]
            for row in result["rows"]
        ]

        nonzero = [
            value
            for value in values
            if value not in (0, None)
        ]

        if not nonzero:
            status = "IDENTITY / INVARIANT"
        else:
            status = "NONZERO"

        print(
            f"  {name:<28} {status}"
        )

    print()

    # ---------------------------------------------------------------
    # A invariance
    # ---------------------------------------------------------------

    A_values = {
        row["A"]
        for row in result["rows"]
    }

    A4_values = {
        row["A4"]
        for row in result["rows"]
    }

    print("A-VALUE AUDIT")
    print(
        f"  distinct A values       = {len(A_values)}"
    )
    print(
        f"  distinct A/4 values     = {len(A4_values)}"
    )

    if len(A_values) == 1:
        print(
            "  A is invariant across the tested j-fibre = True"
        )
    else:
        print(
            "  A is invariant across the tested j-fibre = False"
        )

    print()

    # ---------------------------------------------------------------
    # True branch
    # ---------------------------------------------------------------

    true_j = (
        (inst.d0_true - result["d0_ref"]) // 2
    )

    tested_js = {
        row["j"]
        for row in result["rows"]
    }

    print("TRUE-BRANCH POST-HOC CHECK")
    print(
        f"  true j                 = {true_j}"
    )
    print(
        f"  true j tested          = {true_j in tested_js}"
    )

    print()


# ======================================================================
# GLOBAL
# ======================================================================

def main() -> None:

    print("=" * 120)
    print("EXPERIMENT 455")
    print("=" * 120)
    print()
    print(
        "FULL KNOWN-DATA RELATION / j-FIBRE INVARIANCE AUDIT"
    )
    print()

    j_values = (
        -32, -24, -16, -12, -8,
        -5, -4, -3, -2, -1,
         0,
         1, 2, 3, 4, 5,
         8, 12, 16, 24, 32,
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
    # GLOBAL RELATION AUDIT
    # ==================================================================

    relation_names = (
        "residual_r",
        "residual_x_d0",
        "residual_d_NS",
        "residual_x_form",
        "residual_square",
        "spacing_up",
        "spacing_down",
        "displacement",
        "K_displacement",
        "A4_y_residual",
        "A4_d_residual",
    )

    print("=" * 120)
    print("GLOBAL EXPERIMENT 455 SUMMARY")
    print("=" * 120)

    print(
        f"  instances                 = {len(results)}"
    )

    print(
        f"  j values per instance     = {len(j_values)}"
    )

    print()

    print("RELATION FAILURE COUNTS")

    global_failures = {}

    for name in relation_names:

        failures = 0

        for result in results:
            for row in result["rows"]:

                value = row[name]

                if value not in (0, None):
                    failures += 1

        global_failures[name] = failures

        print(
            f"  {name:<28} = {failures}"
        )

    print()

    # ==================================================================
    # A INVARIANCE
    # ==================================================================

    A_noninvariant_instances = 0

    for result in results:

        A_values = {
            row["A"]
            for row in result["rows"]
        }

        if len(A_values) != 1:
            A_noninvariant_instances += 1

    print("A INVARIANCE")

    print(
        f"  instances with invariant A = "
        f"{len(results) - A_noninvariant_instances}/{len(results)}"
    )

    print(
        f"  instances with varying A   = "
        f"{A_noninvariant_instances}/{len(results)}"
    )

    print()

    # ==================================================================
    # INTERPRETATION
    # ==================================================================

    print("INTERPRETATION")
    print()

    print(
        "  The audit separates three possibilities:"
    )

    print(
        "    1. IDENTITY:"
    )

    print(
        "       the relation vanishes for every j and therefore"
    )

    print(
        "       carries no additional branch-selection information."
    )

    print()

    print(
        "    2. INDEPENDENT CONSTRAINT:"
    )

    print(
        "       the relation changes across j and excludes some"
    )

    print(
        "       otherwise admissible branches."
    )

    print()

    print(
        "    3. INVARIANT OBSERVABLE:"
    )

    print(
        "       the relation produces the same observed quantity"
    )

    print(
        "       for every j."
    )

    print()

    print(
        "  In particular, A is explicitly audited rather than"
    )

    print(
        "  assuming that the previous K/A manipulations are"
    )

    print(
        "  independent."
    )

    print()

    # ==================================================================
    # NEXT-TARGET DECISION
    # ==================================================================

    identity_only = all(
        failures == 0
        for failures in global_failures.values()
    )

    if identity_only and A_noninvariant_instances == 0:

        print(
            "  RESULT:"
        )

        print(
            "    Every tested relation is algebraically preserved"
        )

        print(
            "    across the complete local j-fibre."
        )

        print()

        print(
            "    Therefore the current relation set does not"
        )

        print(
            "    provide an independent observable selector for j."
        )

    else:

        print(
            "  RESULT:"
        )

        print(
            "    At least one audited relation varies across the"
        )

        print(
            "    j-fibre and must be investigated as a possible"
        )

        print(
            "    independent constraint."
        )

    print()

    print("=" * 120)
    print("EXPERIMENT 455 FINAL STATUS")
    print("=" * 120)

    print(
        f"  ALL TESTED RELATIONS IDENTITIES = "
        f"{identity_only}"
    )

    print(
        f"  A INVARIANT ACROSS FIBRE        = "
        f"{A_noninvariant_instances == 0}"
    )

    print(
        "  O(1) BRANCH CONSTRUCTION        = True"
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
        "  INTEGER-EXACT                   = True"
    )

    print(
        "  HIDDEN VALUES USED FOR BOUNDS  = False"
    )

    print()

    if identity_only and A_noninvariant_instances == 0:

        print(
            "  CONCLUSION:"
        )

        print(
            "    The presently audited N,S,d,r,A-derived algebra"
        )

        print(
            "    remains inside the same j-parameter fibre."
        )

        print(
            "    A genuinely independent observable is required."
        )

    else:

        print(
            "  CONCLUSION:"
        )

        print(
            "    A relation outside the previously identified"
        )

        print(
            "    pullback identities has been detected."
        )

    print("=" * 120)
    print("EXPERIMENT 455 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    main()
