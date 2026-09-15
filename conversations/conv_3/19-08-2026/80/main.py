#!/usr/bin/env python3

"""
========================================================================================================================
EXPERIMENT 456R2
========================================================================================================================

OBSERVED-vs-DERIVED QUANTITY / j-FIBRE SEPARATION AUDIT

CORRECTIONS FROM 456
  1. d0 parity is derived from
         r + d0^2 - d^2 == 0 (mod 4)
     rather than inferred from d parity.
  2. x0 is then derived from
         x0 = d - d0_ref.
  3. Both possible d0 parity classes are tested explicitly.
  4. The canonical reference branch uses the unique admissible
     parity class together with the smallest nonnegative x0
     satisfying x0 + d0_ref = d.
  5. No hidden value is used to construct the reference.

QUESTION

  Experiment 455 observed that an A-like quantity varied across
  the j-fibre.

  This experiment asks whether that variation is:

      (a) an independently observed quantity, or
      (b) merely a derived quantity that changes because j changes.

RULES

  exact integer arithmetic only
  no factorization
  no giant K enumeration
  no giant d0 enumeration
  no CRT Cartesian product

  Hidden true values are post-hoc diagnostics only.
========================================================================================================================
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable


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
    -32, -24, -16, -12, -8, -5, -4, -3, -2, -1,
    0, 1, 2, 3, 4, 5, 8, 12, 16, 24, 32,
]


# ----------------------------------------------------------------------------------------------------------------------
# PARITY / REFERENCE CONSTRUCTION
# ----------------------------------------------------------------------------------------------------------------------

def valid_d0_parities(inst: Instance) -> list[int]:
    """
    From

        4K = r + d0^2 - d^2

    integrality requires

        d0^2 == d^2 - r (mod 4).

    Since d0^2 mod 4 is either 0 or 1, this determines the admissible
    parity class.
    """
    target = (inst.d * inst.d - inst.r) % 4

    valid: list[int] = []

    for parity in (0, 1):
        if (parity * parity) % 4 == target:
            valid.append(parity)

    return valid


def canonical_reference(inst: Instance) -> tuple[int, int, int, list[int]]:
    """
    Construct the canonical known-data reference.

    Let p be the admissible parity of d0.

    We need
        x0 + d0_ref = d
    with
        d0_ref == p (mod 2).

    Choose x0 in {0,1} such that
        d0_ref = d - x0
    has the required parity.

    This produces the canonical nearby reference:
        x0 = 0 or 1.
    """
    valid = valid_d0_parities(inst)

    if len(valid) != 1:
        raise AssertionError(
            f"expected unique admissible d0 parity: "
            f"instance={inst.idx}, valid={valid}"
        )

    required_d0_parity = valid[0]

    candidates: list[tuple[int, int]] = []

    for x0 in (0, 1):
        d0_ref = inst.d - x0

        if d0_ref % 2 == required_d0_parity:
            candidates.append((x0, d0_ref))

    if len(candidates) != 1:
        raise AssertionError(
            f"reference branch construction failed: "
            f"instance={inst.idx}, candidates={candidates}, "
            f"required_parity={required_d0_parity}"
        )

    x0, d0_ref = candidates[0]

    numerator = (
        inst.r
        + d0_ref * d0_ref
        - inst.d * inst.d
    )

    if numerator % 4 != 0:
        raise AssertionError(
            f"reference K not integral after parity derivation: "
            f"instance={inst.idx}, "
            f"d0_ref={d0_ref}, "
            f"numerator_mod4={numerator % 4}"
        )

    K_ref = numerator // 4

    return x0, d0_ref, K_ref, valid


# ----------------------------------------------------------------------------------------------------------------------
# CORE FIBRE
# ----------------------------------------------------------------------------------------------------------------------

def fibre_state(
    x0: int,
    d0_ref: int,
    K_ref: int,
    j: int,
) -> tuple[int, int, int]:
    """
    Exact fibre parameterization:

        d0(j) = d0_ref + 2j
        x(j)  = x0 - 2j
        K(j)  = K_ref + j*d0_ref + j^2
    """
    d0 = d0_ref + 2 * j
    x = x0 - 2 * j
    K = K_ref + j * d0_ref + j * j
    return d0, x, K


# ----------------------------------------------------------------------------------------------------------------------
# QUANTITY DEFINITIONS
# ----------------------------------------------------------------------------------------------------------------------

def q_fixed_N(inst: Instance, *_):
    return inst.N


def q_fixed_S(inst: Instance, *_):
    return inst.S


def q_fixed_d(inst: Instance, *_):
    return inst.d


def q_fixed_r(inst: Instance, *_):
    return inst.r


def q_fixed_NS_relation(inst: Instance, *_):
    return inst.N - 2 * inst.S + 1


def q_fixed_d2_minus_r(inst: Instance, *_):
    return inst.d * inst.d - inst.r


def q_derived_d0(inst, x0, d0_ref, K_ref, j):
    return fibre_state(x0, d0_ref, K_ref, j)[0]


def q_derived_x(inst, x0, d0_ref, K_ref, j):
    return fibre_state(x0, d0_ref, K_ref, j)[1]


def q_derived_K(inst, x0, d0_ref, K_ref, j):
    return fibre_state(x0, d0_ref, K_ref, j)[2]


def q_derived_A_455(inst, x0, d0_ref, K_ref, j):
    """
    Experiment-455 A-like quantity:

        A(j) = 4K(j) + d^2.
    """
    _, _, K = fibre_state(
        x0,
        d0_ref,
        K_ref,
        j,
    )
    return 4 * K + inst.d * inst.d


def q_derived_A_pullback(inst, x0, d0_ref, K_ref, j):
    d0, _, _ = fibre_state(
        x0,
        d0_ref,
        K_ref,
        j,
    )
    return inst.r + d0 * d0


def q_derived_square(inst, x0, d0_ref, K_ref, j):
    d0, x, _ = fibre_state(
        x0,
        d0_ref,
        K_ref,
        j,
    )
    return (x + d0) ** 2


def q_derived_K_from_x(inst, x0, d0_ref, K_ref, j):
    _, x, _ = fibre_state(
        x0,
        d0_ref,
        K_ref,
        j,
    )

    numerator = (
        inst.r
        - 2 * inst.d * x
        + x * x
    )

    if numerator % 4 != 0:
        raise AssertionError(
            f"K_from_x nonintegral: "
            f"instance={inst.idx}, j={j}"
        )

    return numerator // 4


QUANTITIES: list[
    tuple[str, str, Callable[..., int]]
] = [
    ("N", "OBSERVED/FIXED", q_fixed_N),
    ("S", "OBSERVED/FIXED", q_fixed_S),
    ("d", "OBSERVED/FIXED", q_fixed_d),
    ("r", "OBSERVED/FIXED", q_fixed_r),
    ("N-2S+1", "OBSERVED/FIXED", q_fixed_NS_relation),
    ("d^2-r", "OBSERVED/FIXED", q_fixed_d2_minus_r),
    ("d0(j)", "DERIVED/FIBRE", q_derived_d0),
    ("x(j)", "DERIVED/FIBRE", q_derived_x),
    ("K(j)", "DERIVED/FIBRE", q_derived_K),
    ("A_455(j)=4K+d^2", "DERIVED/FIBRE", q_derived_A_455),
    ("A_pullback(j)=r+d0^2", "DERIVED/FIBRE", q_derived_A_pullback),
    ("(x+d0)^2", "DERIVED/FIBRE", q_derived_square),
    ("K_from_x(j)", "DERIVED/FIBRE", q_derived_K_from_x),
]


# ----------------------------------------------------------------------------------------------------------------------
# IDENTITY AUDITS
# ----------------------------------------------------------------------------------------------------------------------

def audit_identities(
    inst: Instance,
    x0: int,
    d0_ref: int,
    K_ref: int,
) -> dict[str, bool]:

    failures = {
        "parity_equation": False,
        "r_relation": False,
        "x_d0_relation": False,
        "d_NS_relation": False,
        "x_form": False,
        "A_455_equals_pullback": False,
        "A_variation": False,
        "K_displacement": False,
        "K_from_x": False,
    }

    A0 = 4 * K_ref + inst.d * inst.d

    for j in J_VALUES:
        d0, x, K = fibre_state(
            x0,
            d0_ref,
            K_ref,
            j,
        )

        # Parity/integrality condition.
        if (
            inst.r
            + d0 * d0
            - inst.d * inst.d
        ) % 4 != 0:
            failures["parity_equation"] = True

        # Original relation.
        if (
            4 * K
            != inst.r
            + d0 * d0
            - inst.d * inst.d
        ):
            failures["r_relation"] = True

        # x + d0 = d.
        if x + d0 != inst.d:
            failures["x_d0_relation"] = True

        # d = N - 2S + 1.
        if inst.d != inst.N - 2 * inst.S + 1:
            failures["d_NS_relation"] = True

        # 4K = r - 2dx + x^2.
        if (
            4 * K
            != inst.r
            - 2 * inst.d * x
            + x * x
        ):
            failures["x_form"] = True

        # A_455 = pullback.
        A1 = 4 * K + inst.d * inst.d
        A2 = inst.r + d0 * d0

        if A1 != A2:
            failures["A_455_equals_pullback"] = True

        # A(j)-A(0).
        predicted_A_delta = (
            4 * j * d0_ref
            + 4 * j * j
        )

        if A1 - A0 != predicted_A_delta:
            failures["A_variation"] = True

        # K displacement.
        predicted_K_delta = (
            j * d0_ref
            + j * j
        )

        if K - K_ref != predicted_K_delta:
            failures["K_displacement"] = True

        # K reconstructed from x.
        numerator = (
            inst.r
            - 2 * inst.d * x
            + x * x
        )

        if (
            numerator % 4 != 0
            or numerator // 4 != K
        ):
            failures["K_from_x"] = True

    return failures


# ----------------------------------------------------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------------------------------------------------

def main() -> None:
    print("=" * 120)
    print("EXPERIMENT 456R2")
    print("=" * 120)
    print()
    print("OBSERVED-vs-DERIVED QUANTITY / j-FIBRE SEPARATION AUDIT")
    print()
    print("PARITY CORRECTION")
    print("  d0 parity is derived from:")
    print("      d0^2 == d^2-r (mod 4)")
    print()
    print("RULES")
    print("  exact integer arithmetic only")
    print("  no factorization")
    print("  no giant K enumeration")
    print("  no giant d0 enumeration")
    print("  no CRT Cartesian product")
    print()

    global_failures = {
        "parity_equation": 0,
        "r_relation": 0,
        "x_d0_relation": 0,
        "d_NS_relation": 0,
        "x_form": 0,
        "A_455_equals_pullback": 0,
        "A_variation": 0,
        "K_displacement": 0,
        "K_from_x": 0,
    }

    fixed_quantity_count = 0
    derived_invariant_count = 0
    derived_varying_count = 0

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

        x0, d0_ref, K_ref, valid_parities = (
            canonical_reference(inst)
        )

        print("PARITY DERIVATION")
        print(f"  valid d0 parities = {valid_parities}")
        print(f"  required parity   = {valid_parities[0]}")
        print()

        print("KNOWN-DATA REFERENCE")
        print(f"  x0     = {x0}")
        print(f"  d0_ref = {d0_ref}")
        print(f"  K_ref  = {K_ref}")
        print()

        print("REFERENCE CONSISTENCY")
        print(
            f"  x0 + d0_ref = d = "
            f"{x0 + d0_ref == inst.d}"
        )
        print(
            f"  K_ref integer = "
            f"{(
                inst.r
                + d0_ref * d0_ref
                - inst.d * inst.d
            ) % 4 == 0}"
        )
        print()

        print("POST-HOC TRUE VALUES")
        print(f"  true d0 = {inst.d0_true}")
        print(f"  true x  = {inst.x_true}")
        print(f"  true K  = {inst.K_true}")

        if (
            inst.d0_true - d0_ref
        ) % 2 == 0:
            true_j = (
                inst.d0_true - d0_ref
            ) // 2
        else:
            true_j = None

        print(f"  true j  = {true_j}")
        print()

        # --------------------------------------------------------------------------------------------------------------
        # Quantity classification
        # --------------------------------------------------------------------------------------------------------------

        print("QUANTITY CLASSIFICATION")
        print(
            "  quantity                              class                 "
            "distinct-values   j-invariant"
        )
        print("  " + "-" * 105)

        for name, category, fn in QUANTITIES:
            values = [
                fn(
                    inst,
                    x0,
                    d0_ref,
                    K_ref,
                    j,
                )
                for j in J_VALUES
            ]

            distinct = set(values)
            invariant = len(distinct) == 1

            if category == "OBSERVED/FIXED":
                if not invariant:
                    raise AssertionError(
                        f"observed quantity changed across fibre: "
                        f"instance={inst.idx}, quantity={name}"
                    )

                fixed_quantity_count += 1

            else:
                if invariant:
                    derived_invariant_count += 1
                else:
                    derived_varying_count += 1

            print(
                f"  {name:<38}"
                f"{category:<21}"
                f"{len(distinct):>8}"
                f"{str(invariant):>16}"
            )

        print()

        # --------------------------------------------------------------------------------------------------------------
        # A audit
        # --------------------------------------------------------------------------------------------------------------

        print("A-VARIATION AUDIT")

        A_values = [
            q_derived_A_455(
                inst,
                x0,
                d0_ref,
                K_ref,
                j,
            )
            for j in J_VALUES
        ]

        A0 = A_values[J_VALUES.index(0)]

        print(f"  A(0) = {A0}")
        print(
            f"  distinct A(j) values = "
            f"{len(set(A_values))}"
        )
        print(
            f"  A varies across fibre = "
            f"{len(set(A_values)) != 1}"
        )
        print()

        print("  EXACT PULLBACK")
        print("    A(j)")
        print("      = 4K(j) + d^2")
        print("      = r + d0(j)^2")
        print()

        print("  EXACT VARIATION")
        print(
            "    A(j)-A(0)"
            " = 4*j*d0_ref + 4*j^2"
        )
        print()

        print("  SELECTED BRANCHES")

        for j in [-16, -8, -1, 0, 1, 8, 16]:
            d0, x, K = fibre_state(
                x0,
                d0_ref,
                K_ref,
                j,
            )

            A = (
                4 * K
                + inst.d * inst.d
            )

            predicted = (
                A0
                + 4 * j * d0_ref
                + 4 * j * j
            )

            print(
                f"    j={j:>3} "
                f"d0={d0} "
                f"x={x:>5} "
                f"A={A} "
                f"predicted={predicted} "
                f"match={A == predicted}"
            )

        print()

        # --------------------------------------------------------------------------------------------------------------
        # Identity checks
        # --------------------------------------------------------------------------------------------------------------

        failures = audit_identities(
            inst,
            x0,
            d0_ref,
            K_ref,
        )

        print("IDENTITY FAILURE AUDIT")

        for key, failed in failures.items():
            if failed:
                global_failures[key] += 1

            print(
                f"  {key:<30} = "
                f"{'FAIL' if failed else 'IDENTITY'}"
            )

        print()

        # --------------------------------------------------------------------------------------------------------------
        # True branch post-hoc
        # --------------------------------------------------------------------------------------------------------------

        print("TRUE BRANCH POST-HOC CHECK")

        if true_j is None:
            print(
                "  true branch has incompatible reference "
                "parity displacement"
            )
        else:
            recovered = fibre_state(
                x0,
                d0_ref,
                K_ref,
                true_j,
            )

            print(f"  true j = {true_j}")
            print(
                f"  reconstructed d0 = "
                f"{recovered[0]}"
            )
            print(
                f"  reconstructed x  = "
                f"{recovered[1]}"
            )
            print(
                f"  reconstructed K  = "
                f"{recovered[2]}"
            )
            print(
                "  exact true branch "
                f"reconstruction = "
                f"{recovered == (inst.d0_true, inst.x_true, inst.K_true)}"
            )

        print()

    # ------------------------------------------------------------------------------------------------------------------
    # Global summary
    # ------------------------------------------------------------------------------------------------------------------

    print("=" * 120)
    print("GLOBAL EXPERIMENT 456R2 SUMMARY")
    print("=" * 120)
    print()

    print(f"  instances                     = {len(INSTANCES)}")
    print(
        f"  j values per instance        = "
        f"{len(J_VALUES)}"
    )
    print(
        f"  observed/fixed quantities    = "
        f"{fixed_quantity_count}"
    )
    print(
        f"  derived invariant quantities = "
        f"{derived_invariant_count}"
    )
    print(
        f"  derived varying quantities   = "
        f"{derived_varying_count}"
    )
    print()

    print("GLOBAL FAILURE COUNTS")

    for key, count in global_failures.items():
        print(
            f"  {key:<32} = {count}"
        )

    print()

    print("CENTRAL ALGEBRAIC RESULT")
    print()
    print("  For every tested j:")
    print()
    print("      d0(j) = d0_ref + 2j")
    print("      x(j)  = x0 - 2j")
    print("      K(j)  = K_ref + j*d0_ref + j^2")
    print()
    print("  and:")
    print()
    print("      4K(j) = r + d0(j)^2 - d^2")
    print()
    print("  therefore:")
    print()
    print("      4K(j) + d^2 = r + d0(j)^2")
    print()
    print("  Hence the varying A quantity from Experiment 455")
    print("  is exactly a derived function of the existing fibre.")
    print()

    print("PARITY RESULT")
    print()
    print("  The admissible d0 parity is NOT inferred from d alone.")
    print()
    print("  It is derived from:")
    print()
    print("      d0^2 == d^2-r (mod 4).")
    print()
    print("  This is why instance 2 has")
    print()
    print("      d  = 10117074  (even)")
    print("      d0 = 10117073  (odd)")
    print()
    print("  and nevertheless admits an integral K.")
    print()

    print("INDEPENDENT-CONSTRAINT TEST")
    print()
    print("  A derived quantity may vary across j.")
    print("  That variation alone does NOT identify the true branch.")
    print()
    print("  To break the fibre, a quantity must be independently")
    print("  observed/fixed and cannot be reconstructed from")
    print("  d0(j), x(j), K(j) through the existing identities.")
    print()

    all_clean = all(
        value == 0
        for value in global_failures.values()
    )

    print("=" * 120)
    print("EXPERIMENT 456R2 FINAL STATUS")
    print("=" * 120)
    print(
        f"  ALL ALGEBRAIC IDENTITY CHECKS = "
        f"{all_clean}"
    )
    print("  OBSERVED QUANTITIES J-INVARIANT = True")
    print("  A_455 VARIES ACROSS FIBRE      = True")
    print("  A_455 IS INDEPENDENT           = False")
    print("  CORRECT MOD-4 PARITY DERIVATION = True")
    print("  O(1) BRANCH CONSTRUCTION       = True")
    print("  GIANT K ENUMERATION            = False")
    print("  GIANT d0 ENUMERATION           = False")
    print("  CRT CARTESIAN PRODUCT          = False")
    print("  INTEGER-EXACT                  = True")
    print()
    print("  CONCLUSION:")
    print("    The Experiment-455 A variation is a derived fibre")
    print("    effect. The reference construction is parity-correct")
    print("    because d0 parity is derived directly from the mod-4")
    print("    integrality condition.")
    print("=" * 120)
    print("EXPERIMENT 456R2 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    main()