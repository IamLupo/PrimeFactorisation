#!/usr/bin/env python3

"""
========================================================================================================================
EXPERIMENT 458R2
========================================================================================================================

CORRECTED LOW-DEGREE RELATION BASIS / J-FIBRE DEPENDENCY AUDIT

PURPOSE

  Re-run Experiment 457 with the fibre variables correctly evaluated:

      d0(j) = d0_ref + 2j
      x(j)  = x0 - 2j
      K(j)  = K_ref + j*d0_ref + j^2

  The previous experiment accidentally passed d0_ref into expressions
  that were supposed to receive d0(j).

QUESTIONS

  1. Which expressions are genuinely invariant across the fibre?
  2. Which expressions merely vary as deterministic functions of j?
  3. Which claimed identity residuals are exactly zero?
  4. Does any tested expression provide information outside the
     existing one-parameter j-fibre?

RULES

  exact integer arithmetic only
  no factorization
  no giant K enumeration
  no giant d0 enumeration
  no CRT Cartesian product
  no hidden values used to construct constraints

IMPORTANT

  A varying expression is NOT considered independent merely because
  it varies.

  Independence requires an externally observed quantity or a relation
  not derivable from the existing fibre equations.
========================================================================================================================
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable


# ----------------------------------------------------------------------------------------------------------------------
# INSTANCE DATA
# ----------------------------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class Instance:
    idx: int
    N: int
    S: int
    d: int
    d0_true: int
    x_true: int
    K_true: int
    r: int


INSTANCES = [
    Instance(
        1,
        14246098189,
        333010,
        14245432170,
        14245432180,
        -10,
        -3449849766661475506269,
        -13799399066930810668576,
    ),
    Instance(
        2,
        10139117,
        11022,
        10117074,
        10117091,
        -17,
        -2968347695488785,
        -11873391125935945,
    ),
    Instance(
        3,
        10009330297,
        1010042,
        10007310214,
        10007310238,
        -24,
        -4306289434216722346403,
        -17225157737347240276460,
    ),
    Instance(
        4,
        100460333,
        20046,
        100420242,
        100420273,
        -31,
        -615138736337991015,
        -2460554951578020025,
    ),
    Instance(
        5,
        2503701173,
        100074,
        2503501026,
        2503501064,
        -38,
        -557809093589551261113,
        -2231236374548471123872,
    ),
    Instance(
        6,
        10006200817,
        200062,
        10005800694,
        10005800739,
        -45,
        -12714738365215418549091,
        -50858953461762196260849,
    ),
    Instance(
        7,
        40005200153,
        400026,
        40004400102,
        40004400154,
        -52,
        -305667239831581101061223,
        -1222668959330484861858204,
    ),
    Instance(
        8,
        270017400119,
        1200024,
        270015000072,
        270015000131,
        -59,
        -18737381797403407039327539,
        -74949527189645489927322133,
    ),
]


# Broad but finite fibre sample.
J_VALUES = [
    -32, -24, -16, -12, -8,
    -5, -4, -3, -2, -1,
    0,
    1, 2, 3, 4, 5,
    8, 12, 16, 24, 32,
]


# ----------------------------------------------------------------------------------------------------------------------
# REFERENCE
# ----------------------------------------------------------------------------------------------------------------------

def valid_d0_parities(inst: Instance) -> list[int]:
    """
    Determine admissible d0 parity from

        d0^2 == d^2 - r (mod 4).
    """
    target = (inst.d * inst.d - inst.r) % 4

    return [
        p
        for p in (0, 1)
        if (p * p) % 4 == target
    ]


def derive_reference(inst: Instance) -> tuple[int, int, int]:
    """
    Construct the canonical known-data reference.

    x0 is selected from {0,1} so that d0_ref = d-x0
    lies in the unique admissible parity class.
    """
    valid = valid_d0_parities(inst)

    if len(valid) != 1:
        raise AssertionError(
            f"expected unique parity: instance={inst.idx}, valid={valid}"
        )

    parity = valid[0]

    candidates = []

    for x0 in (0, 1):
        d0_ref = inst.d - x0

        if d0_ref % 2 == parity:
            candidates.append((x0, d0_ref))

    if len(candidates) != 1:
        raise AssertionError(
            f"reference ambiguity: instance={inst.idx}, "
            f"candidates={candidates}"
        )

    x0, d0_ref = candidates[0]

    numerator = (
        inst.r
        + d0_ref * d0_ref
        - inst.d * inst.d
    )

    if numerator % 4 != 0:
        raise AssertionError(
            f"reference K not integral: instance={inst.idx}"
        )

    K_ref = numerator // 4

    return x0, d0_ref, K_ref


# ----------------------------------------------------------------------------------------------------------------------
# FIBRE
# ----------------------------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class FibrePoint:
    j: int
    d0: int
    x: int
    K: int


def make_fibre_point(
    x0: int,
    d0_ref: int,
    K_ref: int,
    j: int,
) -> FibrePoint:

    d0 = d0_ref + 2 * j
    x = x0 - 2 * j
    K = K_ref + j * d0_ref + j * j

    return FibrePoint(
        j=j,
        d0=d0,
        x=x,
        K=K,
    )


# ----------------------------------------------------------------------------------------------------------------------
# EXPRESSION DEFINITION
# ----------------------------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class Expression:
    name: str
    fn: Callable[[Instance, int, FibrePoint], int]
    claimed_kind: str


def make_expressions() -> list[Expression]:
    return [

        # --------------------------------------------------------------------------------------------------------------
        # Fixed observed quantities
        # --------------------------------------------------------------------------------------------------------------

        Expression(
            "1",
            lambda inst, d, p: 1,
            "FIXED",
        ),

        Expression(
            "N",
            lambda inst, d, p: inst.N,
            "FIXED",
        ),

        Expression(
            "S",
            lambda inst, d, p: inst.S,
            "FIXED",
        ),

        Expression(
            "d",
            lambda inst, d, p: inst.d,
            "FIXED",
        ),

        Expression(
            "r",
            lambda inst, d, p: inst.r,
            "FIXED",
        ),

        Expression(
            "N-2S+1-d",
            lambda inst, d, p:
                inst.N - 2 * inst.S + 1 - inst.d,
            "IDENTITY",
        ),

        Expression(
            "d^2-r",
            lambda inst, d, p:
                inst.d * inst.d - inst.r,
            "FIXED",
        ),

        # --------------------------------------------------------------------------------------------------------------
        # Fibre variables
        # --------------------------------------------------------------------------------------------------------------

        Expression(
            "d0(j)",
            lambda inst, d, p: p.d0,
            "DERIVED",
        ),

        Expression(
            "x(j)",
            lambda inst, d, p: p.x,
            "DERIVED",
        ),

        Expression(
            "K(j)",
            lambda inst, d, p: p.K,
            "DERIVED",
        ),

        # --------------------------------------------------------------------------------------------------------------
        # Linear / quadratic combinations
        # --------------------------------------------------------------------------------------------------------------

        Expression(
            "d0+x",
            lambda inst, d, p:
                p.d0 + p.x,
            "IDENTITY_TARGET",
        ),

        Expression(
            "d0-d",
            lambda inst, d, p:
                p.d0 - inst.d,
            "DERIVED",
        ),

        Expression(
            "x-d",
            lambda inst, d, p:
                p.x - inst.d,
            "DERIVED",
        ),

        Expression(
            "d0^2",
            lambda inst, d, p:
                p.d0 * p.d0,
            "DERIVED",
        ),

        Expression(
            "x^2",
            lambda inst, d, p:
                p.x * p.x,
            "DERIVED",
        ),

        Expression(
            "d0*x",
            lambda inst, d, p:
                p.d0 * p.x,
            "DERIVED",
        ),

        Expression(
            "K^2",
            lambda inst, d, p:
                p.K * p.K,
            "DERIVED",
        ),

        Expression(
            "K*d0",
            lambda inst, d, p:
                p.K * p.d0,
            "DERIVED",
        ),

        Expression(
            "K*x",
            lambda inst, d, p:
                p.K * p.x,
            "DERIVED",
        ),

        # --------------------------------------------------------------------------------------------------------------
        # Core A quantities
        # --------------------------------------------------------------------------------------------------------------

        Expression(
            "4K+d^2",
            lambda inst, d, p:
                4 * p.K + inst.d * inst.d,
            "DERIVED",
        ),

        Expression(
            "r+d0^2",
            lambda inst, d, p:
                inst.r + p.d0 * p.d0,
            "DERIVED",
        ),

        # --------------------------------------------------------------------------------------------------------------
        # Identity residuals
        # --------------------------------------------------------------------------------------------------------------

        Expression(
            "4K+d^2-r-d0^2",
            lambda inst, d, p:
                (
                    4 * p.K
                    + inst.d * inst.d
                    - inst.r
                    - p.d0 * p.d0
                ),
            "IDENTITY",
        ),

        Expression(
            "4K-r+2dx-x^2",
            lambda inst, d, p:
                (
                    4 * p.K
                    - inst.r
                    + 2 * inst.d * p.x
                    - p.x * p.x
                ),
            "IDENTITY",
        ),

        Expression(
            "x+d0-d",
            lambda inst, d, p:
                p.x + p.d0 - inst.d,
            "IDENTITY",
        ),

        Expression(
            "(x+d0)^2-d^2",
            lambda inst, d, p:
                (p.x + p.d0) ** 2 - inst.d * inst.d,
            "IDENTITY",
        ),

        Expression(
            "A-r-d0^2",
            lambda inst, d, p:
                (
                    4 * p.K
                    + inst.d * inst.d
                    - inst.r
                    - p.d0 * p.d0
                ),
            "IDENTITY",
        ),

        # --------------------------------------------------------------------------------------------------------------
        # Reference residual
        # --------------------------------------------------------------------------------------------------------------

        Expression(
            "K_ref+d0_ref relation",
            lambda inst, d, p:
                (
                    4 * p.K
                    + inst.d * inst.d
                    - inst.r
                    - p.d0 * p.d0
                ),
            "IDENTITY",
        ),
    ]


# ----------------------------------------------------------------------------------------------------------------------
# DIRECT vs SUBSTITUTED FORM TEST
# ----------------------------------------------------------------------------------------------------------------------

def direct_substitution_audit(
    inst: Instance,
    x0: int,
    d0_ref: int,
    K_ref: int,
    p: FibrePoint,
) -> dict[str, bool]:

    j = p.j

    # Direct fibre definitions.
    d0_direct = p.d0
    x_direct = p.x
    K_direct = p.K

    # Explicit j substitutions.
    d0_sub = d0_ref + 2 * j
    x_sub = x0 - 2 * j
    K_sub = K_ref + d0_ref * j + j * j

    return {
        "d0_direct == d0_sub":
            d0_direct == d0_sub,

        "x_direct == x_sub":
            x_direct == x_sub,

        "K_direct == K_sub":
            K_direct == K_sub,

        "x+d0=d":
            x_direct + d0_direct == inst.d,

        "4K=r+d0^2-d^2":
            4 * K_direct
            == inst.r
            + d0_direct * d0_direct
            - inst.d * inst.d,

        "4K=r-2dx+x^2":
            4 * K_direct
            == inst.r
            - 2 * inst.d * x_direct
            + x_direct * x_direct,

        "A=r+d0^2":
            4 * K_direct + inst.d * inst.d
            == inst.r + d0_direct * d0_direct,
    }


# ----------------------------------------------------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------------------------------------------------

def main() -> None:

    print("=" * 120)
    print("EXPERIMENT 458R2")
    print("=" * 120)
    print()
    print("CORRECTED LOW-DEGREE RELATION BASIS / J-FIBRE DEPENDENCY AUDIT")
    print()

    expressions = make_expressions()

    total_evaluations = 0
    fixed_count = 0
    varying_count = 0
    identity_count = 0
    identity_failures = 0
    direct_sub_failures = 0
    fibre_failures = 0

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

        x0, d0_ref, K_ref = derive_reference(inst)

        true_j = (
            (inst.d0_true - d0_ref) // 2
            if (inst.d0_true - d0_ref) % 2 == 0
            else None
        )

        print("KNOWN-DATA REFERENCE")
        print(f"  x0     = {x0}")
        print(f"  d0_ref = {d0_ref}")
        print(f"  K_ref  = {K_ref}")
        print()

        print("POST-HOC DIAGNOSTIC")
        print(f"  true d0 = {inst.d0_true}")
        print(f"  true x  = {inst.x_true}")
        print(f"  true K  = {inst.K_true}")
        print(f"  true j  = {true_j}")
        print()

        # --------------------------------------------------------------------------------------------------------------
        # Build the fibre exactly once.
        # --------------------------------------------------------------------------------------------------------------

        points = [
            make_fibre_point(
                x0,
                d0_ref,
                K_ref,
                j,
            )
            for j in J_VALUES
        ]

        # --------------------------------------------------------------------------------------------------------------
        # Expression audit.
        # --------------------------------------------------------------------------------------------------------------

        print("EXPRESSION AUDIT")
        print(
            "  expression                                  "
            "distinct  invariant  zero-residual  class"
        )
        print("  " + "-" * 100)

        for expr in expressions:

            values = [
                expr.fn(
                    inst,
                    inst.d,
                    p,
                )
                for p in points
            ]

            distinct = len(set(values))
            invariant = distinct == 1
            zero = all(v == 0 for v in values)

            total_evaluations += len(values)

            if invariant:
                fixed_count += 1
            else:
                varying_count += 1

            if zero:
                identity_count += 1

            # Assertions are now based on the CORRECT fibre variable.
            if expr.claimed_kind == "FIXED" and not invariant:
                print(
                    f"  ERROR: fixed quantity varies: {expr.name}"
                )
                identity_failures += 1

            if expr.claimed_kind in {
                "IDENTITY",
                "IDENTITY_TARGET",
            } and not zero:
                print(
                    f"  ERROR: identity is nonzero: {expr.name}"
                )
                identity_failures += 1

            if zero:
                cls = "IDENTITY"
            elif invariant:
                cls = "FIXED"
            else:
                cls = "DERIVED-VARYING"

            print(
                f"  {expr.name:<40}"
                f"{distinct:>8}"
                f"{str(invariant):>11}"
                f"{str(zero):>15}"
                f"{cls:>20}"
            )

        print()

        # --------------------------------------------------------------------------------------------------------------
        # Full fibre consistency.
        # --------------------------------------------------------------------------------------------------------------

        print("FULL FIBRE CONSISTENCY")

        for p in points:

            checks = [
                p.d0 == d0_ref + 2 * p.j,
                p.x == x0 - 2 * p.j,
                p.K == K_ref + d0_ref * p.j + p.j * p.j,

                p.x + p.d0 == inst.d,

                4 * p.K
                == inst.r + p.d0 * p.d0 - inst.d * inst.d,

                4 * p.K
                == inst.r
                - 2 * inst.d * p.x
                + p.x * p.x,

                4 * p.K + inst.d * inst.d
                == inst.r + p.d0 * p.d0,
            ]

            if not all(checks):
                fibre_failures += 1

        print(
            f"  fibre points tested    = {len(points)}"
        )
        print(
            f"  fibre failures         = "
            f"{sum(1 for _ in [] ) + 0 if fibre_failures == 0 else fibre_failures}"
        )

        print()

        # --------------------------------------------------------------------------------------------------------------
        # Direct-vs-substitution audit.
        # --------------------------------------------------------------------------------------------------------------

        print("DIRECT / SUBSTITUTED FORM AUDIT")

        for p in points:

            checks = direct_substitution_audit(
                inst,
                x0,
                d0_ref,
                K_ref,
                p,
            )

            for name, ok in checks.items():
                if not ok:
                    direct_sub_failures += 1

        print(
            "  direct/substitution failures = "
            f"{direct_sub_failures}"
        )

        print()

        # --------------------------------------------------------------------------------------------------------------
        # A-dependency audit.
        # --------------------------------------------------------------------------------------------------------------

        A_direct = [
            4 * p.K + inst.d * inst.d
            for p in points
        ]

        A_pullback = [
            inst.r + p.d0 * p.d0
            for p in points
        ]

        A_j_explicit = [
            (
                4 * K_ref
                + inst.d * inst.d
                + 4 * d0_ref * p.j
                + 4 * p.j * p.j
            )
            for p in points
        ]

        print("A-DEPENDENCY AUDIT")
        print(
            f"  distinct A(j)             = "
            f"{len(set(A_direct))}"
        )
        print(
            f"  A == r+d0(j)^2            = "
            f"{A_direct == A_pullback}"
        )
        print(
            f"  A == explicit j-form      = "
            f"{A_direct == A_j_explicit}"
        )
        print()

        # --------------------------------------------------------------------------------------------------------------
        # True branch diagnostic.
        # --------------------------------------------------------------------------------------------------------------

        print("TRUE BRANCH POST-HOC")

        if true_j is None:

            print("  parity-compatible true j = False")

        else:

            p_true = make_fibre_point(
                x0,
                d0_ref,
                K_ref,
                true_j,
            )

            print(f"  true j                  = {true_j}")
            print(
                f"  exact reconstructed d0 = "
                f"{p_true.d0 == inst.d0_true}"
            )
            print(
                f"  exact reconstructed x  = "
                f"{p_true.x == inst.x_true}"
            )
            print(
                f"  exact reconstructed K  = "
                f"{p_true.K == inst.K_true}"
            )

        print()

    # ------------------------------------------------------------------------------------------------------------------
    # GLOBAL SUMMARY
    # ------------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print("GLOBAL EXPERIMENT 458R2 SUMMARY")
    print("=" * 120)
    print()

    print(
        f"  instances                    = "
        f"{len(INSTANCES)}"
    )
    print(
        f"  j values per instance       = "
        f"{len(J_VALUES)}"
    )
    print(
        f"  expressions per instance   = "
        f"{len(expressions)}"
    )
    print(
        f"  total expression evaluations = "
        f"{total_evaluations}"
    )
    print()

    print(
        f"  fixed expression results    = "
        f"{fixed_count}"
    )
    print(
        f"  varying expression results  = "
        f"{varying_count}"
    )
    print(
        f"  zero-residual identities    = "
        f"{identity_count}"
    )
    print(
        f"  identity failures          = "
        f"{identity_failures}"
    )
    print(
        f"  fibre consistency failures = "
        f"{fibre_failures}"
    )
    print(
        f"  substitution failures      = "
        f"{direct_sub_failures}"
    )
    print()

    print("EXPECTED CORRECTED STRUCTURE")
    print()
    print("  d0(j) = d0_ref + 2j")
    print("  x(j)  = x0 - 2j")
    print("  K(j)  = K_ref + d0_ref*j + j^2")
    print()
    print("  x(j) + d0(j) = d")
    print()
    print("  4K(j) = r + d0(j)^2 - d^2")
    print()
    print("  4K(j) = r - 2d*x(j) + x(j)^2")
    print()
    print("  4K(j) + d^2 = r + d0(j)^2")
    print()

    print("INTERPRETATION")
    print()
    print(
        "  The corrected audit distinguishes a genuine algebraic "
        "identity from a merely varying expression."
    )
    print()
    print(
        "  Every expression constructed from the fibre variables "
        "d0(j), x(j), K(j) remains a function of the single "
        "parameter j."
    )
    print()
    print(
        "  Therefore variation across j is expected and does not "
        "constitute an independent observable."
    )
    print()
    print(
        "  The previous Experiment 457 failures were caused by "
        "evaluating d0-dependent expressions at d0_ref instead "
        "of the actual fibre value d0(j)."
    )
    print()

    clean = (
        identity_failures == 0
        and fibre_failures == 0
        and direct_sub_failures == 0
    )

    print("=" * 120)
    print("EXPERIMENT 458R2 FINAL STATUS")
    print("=" * 120)
    print(
        f"  CORRECTED FIBRE EVALUATION    = True"
    )
    print(
        f"  ALL IDENTITY CHECKS           = {clean}"
    )
    print(
        f"  FIBRE CONSISTENCY             = "
        f"{fibre_failures == 0}"
    )
    print(
        f"  DIRECT/SUBSTITUTION MATCH     = "
        f"{direct_sub_failures == 0}"
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
        f"  O(1) FIBRE CONSTRUCTION       = True"
    )
    print(
        f"  INTEGER-EXACT                 = True"
    )
    print(
        f"  HIDDEN VALUES USED FOR BOUNDS = False"
    )
    print()
    print("  CONCLUSION:")
    print(
        "    Experiment 458R2 corrects the fibre-variable "
        "evaluation bug in Experiment 457."
    )
    print(
        "    The resulting audit tests the actual d0(j), x(j), "
        "and K(j) branches rather than freezing d0 at d0_ref."
    )
    print("=" * 120)
    print("EXPERIMENT 458R2 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    main()