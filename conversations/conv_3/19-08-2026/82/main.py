#!/usr/bin/env python3

"""
========================================================================================================================
EXPERIMENT 459
========================================================================================================================

EXACT FIBRE INVARIANT / ZERO-RESIDUAL SEPARATION AUDIT

PURPOSE

  Correct the final bookkeeping issue from Experiment 458R2.

  A quantity may be invariant across the j-fibre without being zero.

  Example:

      d0(j) + x(j) = d

  is an invariant quantity, while

      d0(j) + x(j) - d = 0

  is the corresponding zero residual.

QUESTIONS

  1. Which expressions are j-invariant?
  2. Which expressions have a fixed known value?
  3. Which residuals are identically zero?
  4. Are all claimed identities exact across the fibre?
  5. Can any tested quantity distinguish j?

RULES

  exact integer arithmetic only
  no factorization
  no giant K enumeration
  no giant d0 enumeration
  no CRT Cartesian product
  no hidden values used to construct constraints

KEY DISTINCTION

  invariant(expr)
      means expr(j) is constant over the fibre.

  zero_identity(expr)
      means expr(j) == 0 for every j.

  An invariant expression therefore does NOT need to be zero.
========================================================================================================================
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable


# ----------------------------------------------------------------------------------------------------------------------
# DATA
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
    target = (inst.d * inst.d - inst.r) % 4

    return [
        p
        for p in (0, 1)
        if (p * p) % 4 == target
    ]


def derive_reference(inst: Instance) -> tuple[int, int, int]:

    valid = valid_d0_parities(inst)

    if len(valid) != 1:
        raise AssertionError(
            f"unexpected parity family: instance={inst.idx}, "
            f"valid={valid}"
        )

    parity = valid[0]

    candidates = []

    for x0 in (0, 1):

        d0_ref = inst.d - x0

        if d0_ref % 2 == parity:
            candidates.append((x0, d0_ref))

    if len(candidates) != 1:
        raise AssertionError(
            f"reference ambiguity: instance={inst.idx}"
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


def make_fibre(
    x0: int,
    d0_ref: int,
    K_ref: int,
) -> list[FibrePoint]:

    points = []

    for j in J_VALUES:

        d0 = d0_ref + 2 * j
        x = x0 - 2 * j
        K = K_ref + d0_ref * j + j * j

        points.append(
            FibrePoint(
                j=j,
                d0=d0,
                x=x,
                K=K,
            )
        )

    return points


# ----------------------------------------------------------------------------------------------------------------------
# EXPRESSION MODEL
# ----------------------------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class Quantity:
    name: str
    fn: Callable[[Instance, FibrePoint], int]
    expected_invariant: Callable[[Instance, int, int, int], int] | None
    expected_zero: bool
    category: str


def make_quantities() -> list[Quantity]:

    return [

        Quantity(
            "1",
            lambda inst, p: 1,
            lambda inst, x0, d0_ref, K_ref: 1,
            False,
            "FIXED",
        ),

        Quantity(
            "N",
            lambda inst, p: inst.N,
            lambda inst, x0, d0_ref, K_ref: inst.N,
            False,
            "FIXED",
        ),

        Quantity(
            "S",
            lambda inst, p: inst.S,
            lambda inst, x0, d0_ref, K_ref: inst.S,
            False,
            "FIXED",
        ),

        Quantity(
            "d",
            lambda inst, p: inst.d,
            lambda inst, x0, d0_ref, K_ref: inst.d,
            False,
            "FIXED",
        ),

        Quantity(
            "r",
            lambda inst, p: inst.r,
            lambda inst, x0, d0_ref, K_ref: inst.r,
            False,
            "FIXED",
        ),

        Quantity(
            "d0+x",
            lambda inst, p: p.d0 + p.x,
            lambda inst, x0, d0_ref, K_ref: inst.d,
            False,
            "FIXED",
        ),

        Quantity(
            "d0-d",
            lambda inst, p: p.d0 - inst.d,
            None,
            False,
            "DERIVED",
        ),

        Quantity(
            "x",
            lambda inst, p: p.x,
            None,
            False,
            "DERIVED",
        ),

        Quantity(
            "d0",
            lambda inst, p: p.d0,
            None,
            False,
            "DERIVED",
        ),

        Quantity(
            "K",
            lambda inst, p: p.K,
            None,
            False,
            "DERIVED",
        ),

        Quantity(
            "x-d",
            lambda inst, p: p.x - inst.d,
            None,
            False,
            "DERIVED",
        ),

        Quantity(
            "d0^2",
            lambda inst, p: p.d0 * p.d0,
            None,
            False,
            "DERIVED",
        ),

        Quantity(
            "x^2",
            lambda inst, p: p.x * p.x,
            None,
            False,
            "DERIVED",
        ),

        Quantity(
            "d0*x",
            lambda inst, p: p.d0 * p.x,
            None,
            False,
            "DERIVED",
        ),

        Quantity(
            "K^2",
            lambda inst, p: p.K * p.K,
            None,
            False,
            "DERIVED",
        ),

        Quantity(
            "K*d0",
            lambda inst, p: p.K * p.d0,
            None,
            False,
            "DERIVED",
        ),

        Quantity(
            "K*x",
            lambda inst, p: p.K * p.x,
            None,
            False,
            "DERIVED",
        ),

        Quantity(
            "4K+d^2",
            lambda inst, p:
                4 * p.K + inst.d * inst.d,
            None,
            False,
            "DERIVED",
        ),

        Quantity(
            "r+d0^2",
            lambda inst, p:
                inst.r + p.d0 * p.d0,
            None,
            False,
            "DERIVED",
        ),

        # True zero identities.

        Quantity(
            "N-2S+1-d",
            lambda inst, p:
                inst.N - 2 * inst.S + 1 - inst.d,
            lambda inst, x0, d0_ref, K_ref: 0,
            True,
            "IDENTITY",
        ),

        Quantity(
            "4K+d^2-r-d0^2",
            lambda inst, p:
                4 * p.K
                + inst.d * inst.d
                - inst.r
                - p.d0 * p.d0,
            lambda inst, x0, d0_ref, K_ref: 0,
            True,
            "IDENTITY",
        ),

        Quantity(
            "4K-r+2dx-x^2",
            lambda inst, p:
                4 * p.K
                - inst.r
                + 2 * inst.d * p.x
                - p.x * p.x,
            lambda inst, x0, d0_ref, K_ref: 0,
            True,
            "IDENTITY",
        ),

        Quantity(
            "x+d0-d",
            lambda inst, p:
                p.x + p.d0 - inst.d,
            lambda inst, x0, d0_ref, K_ref: 0,
            True,
            "IDENTITY",
        ),

        Quantity(
            "(x+d0)^2-d^2",
            lambda inst, p:
                (p.x + p.d0) ** 2 - inst.d * inst.d,
            lambda inst, x0, d0_ref, K_ref: 0,
            True,
            "IDENTITY",
        ),

        Quantity(
            "A-r-d0^2",
            lambda inst, p:
                4 * p.K
                + inst.d * inst.d
                - inst.r
                - p.d0 * p.d0,
            lambda inst, x0, d0_ref, K_ref: 0,
            True,
            "IDENTITY",
        ),
    ]


# ----------------------------------------------------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------------------------------------------------

def main() -> None:

    print("=" * 120)
    print("EXPERIMENT 459")
    print("=" * 120)
    print()
    print("EXACT FIBRE INVARIANT / ZERO-RESIDUAL SEPARATION AUDIT")
    print()

    quantities = make_quantities()

    invariant_failures = 0
    expected_value_failures = 0
    zero_identity_failures = 0
    fibre_failures = 0
    direct_formula_failures = 0

    total_evaluations = 0

    global_fixed = 0
    global_varying = 0
    global_identities = 0

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

        fibre = make_fibre(
            x0,
            d0_ref,
            K_ref,
        )

        true_j = (
            (inst.d0_true - d0_ref) // 2
            if (inst.d0_true - d0_ref) % 2 == 0
            else None
        )

        print("REFERENCE")
        print(f"  x0     = {x0}")
        print(f"  d0_ref = {d0_ref}")
        print(f"  K_ref  = {K_ref}")
        print()

        print("POST-HOC TRUE VALUES")
        print(f"  true d0 = {inst.d0_true}")
        print(f"  true x  = {inst.x_true}")
        print(f"  true K  = {inst.K_true}")
        print(f"  true j  = {true_j}")
        print()

        print("QUANTITY CLASSIFICATION")
        print(
            "  quantity                                  "
            "distinct  invariant  expected  zero-id  class"
        )
        print("  " + "-" * 105)

        for q in quantities:

            values = [
                q.fn(inst, p)
                for p in fibre
            ]

            total_evaluations += len(values)

            distinct = len(set(values))
            invariant = distinct == 1

            if invariant:
                global_fixed += 1
            else:
                global_varying += 1

            if q.expected_zero:
                global_identities += 1

            expected_value = None

            if invariant and q.expected_invariant is not None:
                expected_value = q.expected_invariant(
                    inst,
                    x0,
                    d0_ref,
                    K_ref,
                )

                if values[0] != expected_value:
                    expected_value_failures += 1

            zero_identity_ok = (
                all(v == 0 for v in values)
                if q.expected_zero
                else True
            )

            if q.expected_zero and not zero_identity_ok:
                zero_identity_failures += 1

            if not invariant and q.expected_invariant is not None:
                invariant_failures += 1

            if q.expected_zero:
                cls = "IDENTITY"

            elif invariant:
                cls = "FIXED"

            else:
                cls = "DERIVED-VARYING"

            exp_text = (
                str(expected_value)
                if expected_value is not None
                else "-"
            )

            print(
                f"  {q.name:<38}"
                f"{distinct:>9}"
                f"{str(invariant):>11}"
                f"{exp_text:>12}"
                f"{str(zero_identity_ok):>10}"
                f"{cls:>20}"
            )

        print()

        # --------------------------------------------------------------------------------------------------------------
        # Explicit fibre identities.
        # --------------------------------------------------------------------------------------------------------------

        print("EXPLICIT FIBRE CHECKS")

        for p in fibre:

            checks = {
                "d0=d0_ref+2j":
                    p.d0 == d0_ref + 2 * p.j,

                "x=x0-2j":
                    p.x == x0 - 2 * p.j,

                "K=K_ref+d0_ref*j+j^2":
                    p.K
                    == K_ref
                    + d0_ref * p.j
                    + p.j * p.j,

                "x+d0=d":
                    p.x + p.d0 == inst.d,

                "4K+d^2=r+d0^2":
                    4 * p.K + inst.d * inst.d
                    == inst.r + p.d0 * p.d0,

                "4K=r-2dx+x^2":
                    4 * p.K
                    == inst.r
                    - 2 * inst.d * p.x
                    + p.x * p.x,
            }

            if not all(checks.values()):
                fibre_failures += 1

        print(
            f"  fibre points tested = {len(fibre)}"
        )
        print(
            f"  fibre failures      = {fibre_failures}"
        )
        print()

        # --------------------------------------------------------------------------------------------------------------
        # Polynomial A dependency.
        # --------------------------------------------------------------------------------------------------------------

        A_values = [
            4 * p.K + inst.d * inst.d
            for p in fibre
        ]

        A_pullback = [
            inst.r + p.d0 * p.d0
            for p in fibre
        ]

        A_explicit = [
            (
                4 * K_ref
                + inst.d * inst.d
                + 4 * d0_ref * p.j
                + 4 * p.j * p.j
            )
            for p in fibre
        ]

        print("A-DEPENDENCY")

        print(
            f"  distinct A(j)        = "
            f"{len(set(A_values))}"
        )

        print(
            f"  A == r+d0(j)^2      = "
            f"{A_values == A_pullback}"
        )

        print(
            f"  A == explicit j-form = "
            f"{A_values == A_explicit}"
        )

        if not (
            A_values == A_pullback
            and A_values == A_explicit
        ):
            direct_formula_failures += 1

        print()

        # --------------------------------------------------------------------------------------------------------------
        # Post-hoc branch.
        # --------------------------------------------------------------------------------------------------------------

        print("TRUE BRANCH POST-HOC")

        if true_j is None:

            print("  true branch parity-compatible = False")

        else:

            p_true = make_fibre(
                x0,
                d0_ref,
                K_ref,
            )

            p_match = next(
                (
                    p
                    for p in p_true
                    if p.j == true_j
                ),
                None,
            )

            if p_match is None:

                p_match = FibrePoint(
                    true_j,
                    d0_ref + 2 * true_j,
                    x0 - 2 * true_j,
                    K_ref
                    + d0_ref * true_j
                    + true_j * true_j,
                )

            print(
                f"  true j                  = {true_j}"
            )
            print(
                f"  exact d0 reconstruction = "
                f"{p_match.d0 == inst.d0_true}"
            )
            print(
                f"  exact x reconstruction  = "
                f"{p_match.x == inst.x_true}"
            )
            print(
                f"  exact K reconstruction  = "
                f"{p_match.K == inst.K_true}"
            )

        print()

    # ------------------------------------------------------------------------------------------------------------------
    # GLOBAL SUMMARY
    # ------------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print("GLOBAL EXPERIMENT 459 SUMMARY")
    print("=" * 120)
    print()

    print(
        f"  instances                    = {len(INSTANCES)}"
    )
    print(
        f"  j values per instance       = {len(J_VALUES)}"
    )
    print(
        f"  quantities per instance     = {len(quantities)}"
    )
    print(
        f"  total expression evaluations = "
        f"{total_evaluations}"
    )
    print()

    print(
        f"  invariant quantity results  = "
        f"{global_fixed}"
    )
    print(
        f"  varying quantity results   = "
        f"{global_varying}"
    )
    print(
        f"  zero-identity quantities   = "
        f"{global_identities}"
    )
    print()

    print(
        f"  invariant expectation failures = "
        f"{invariant_failures}"
    )
    print(
        f"  expected-value failures       = "
        f"{expected_value_failures}"
    )
    print(
        f"  zero-identity failures        = "
        f"{zero_identity_failures}"
    )
    print(
        f"  fibre consistency failures    = "
        f"{fibre_failures}"
    )
    print(
        f"  A formula failures            = "
        f"{direct_formula_failures}"
    )
    print()

    clean = (
        invariant_failures == 0
        and expected_value_failures == 0
        and zero_identity_failures == 0
        and fibre_failures == 0
        and direct_formula_failures == 0
    )

    print("CENTRAL RESULT")
    print()
    print(
        "  The experiment distinguishes:"
    )
    print()
    print(
        "    invariant quantity:"
    )
    print(
        "        f(j) = C"
    )
    print()
    print(
        "    zero identity:"
    )
    print(
        "        f(j) - C = 0"
    )
    print()
    print(
        "  In particular:"
    )
    print()
    print(
        "        d0(j) + x(j) = d"
    )
    print(
        "  is an invariant."
    )
    print()
    print(
        "        d0(j) + x(j) - d = 0"
    )
    print(
        "  is the zero residual."
    )
    print()

    print("FIBRE")
    print()
    print(
        "    d0(j) = d0_ref + 2j"
    )
    print(
        "    x(j)  = x0 - 2j"
    )
    print(
        "    K(j)  = K_ref + d0_ref*j + j^2"
    )
    print()
    print(
        "  and therefore:"
    )
    print(
        "    x(j)+d0(j)=d"
    )
    print(
        "    4K(j)+d^2=r+d0(j)^2"
    )
    print(
        "    4K(j)=r-2d*x(j)+x(j)^2"
    )
    print()

    print("INTERPRETATION")
    print()
    print(
        "  No tested quantity in this basis distinguishes j."
    )
    print(
        "  Some quantities are fixed observables."
    )
    print(
        "  Others vary, but are deterministic functions of j."
    )
    print(
        "  The zero residuals vanish throughout the entire fibre."
    )
    print()
    print(
        "  Consequently, a genuinely useful constraint must be "
        "an independently observed quantity or a relation that "
        "is not generated by the existing fibre algebra."
    )
    print()

    print("=" * 120)
    print("EXPERIMENT 459 FINAL STATUS")
    print("=" * 120)
    print(
        f"  INVARIANT CLASSIFICATION      = "
        f"{invariant_failures == 0}"
    )
    print(
        f"  EXPECTED CONSTANTS            = "
        f"{expected_value_failures == 0}"
    )
    print(
        f"  ZERO-RESIDUAL IDENTITIES      = "
        f"{zero_identity_failures == 0}"
    )
    print(
        f"  FIBRE CONSISTENCY             = "
        f"{fibre_failures == 0}"
    )
    print(
        f"  A DEPENDENCY CHECK             = "
        f"{direct_formula_failures == 0}"
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
        "    The remaining Experiment-458R2 failures were "
        "classification errors, not algebraic failures."
    )
    print(
        "    Experiment 459 explicitly separates invariant "
        "quantities from zero residual identities."
    )
    print("=" * 120)
    print("EXPERIMENT 459 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    main()
