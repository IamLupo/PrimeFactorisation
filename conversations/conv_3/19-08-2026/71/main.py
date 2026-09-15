#!/usr/bin/env python3

from __future__ import annotations

from dataclasses import dataclass
from math import isqrt


# ======================================================================================================================
# EXPERIMENT 447
# ======================================================================================================================
#
# D-D0 / X INFORMATION-FLOW AUDIT
#
# QUESTION
#
#   Is x actually an independent reconstruction variable in the present
#   formulation, or is it already determined by the supplied d and d0?
#
# CORE IDENTITY
#
#   x + d0 = d
#
# therefore
#
#   x = d - d0
#
# The experiment deliberately separates:
#
#   A) KNOWN N,S,d ONLY
#   B) KNOWN N,S,d,d0
#   C) K/x quadratic identity
#
# The goal is to determine exactly where x becomes known.
#
# RULES
#   exact integer arithmetic only
#   no factorization
#   no giant K enumeration
#   no giant X enumeration
#   no CRT Cartesian product
#
# The hidden p,q values are used only to construct controlled test data.
# Reconstruction tests do not use them.
#
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

    result = []

    for idx, (p, q, x_true) in enumerate(RAW, start=1):

        N = p * q
        S = p + q

        # Same relation used throughout the preceding experiments.
        d = N - 2 * S + 1

        # Controlled construction.
        d0 = d - x_true

        result.append(
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

    return result


# ======================================================================================================================
# BASIC EXACT IDENTITIES
# ======================================================================================================================

def recover_x_from_d_d0(inst: Instance) -> int:
    return inst.d - inst.d0


def recover_d0_from_d_x(inst: Instance) -> int:
    return inst.d - inst.x_true


def check_d_relation(inst: Instance) -> bool:
    return inst.d == inst.N - 2 * inst.S + 1


def check_x_relation(inst: Instance) -> bool:
    return inst.x_true + inst.d0 == inst.d


# ======================================================================================================================
# KNOWN-N,S-ONLY RADII
# ======================================================================================================================

def candidate_radii_NS(inst: Instance) -> dict[str, int]:

    sqrt_n = isqrt(inst.N)

    candidates = {
        "0": 0,
        "1": 1,
        "sqrtN": sqrt_n,
        "2sqrtN": 2 * sqrt_n,
        "4sqrtN": 4 * sqrt_n,
        "8sqrtN": 8 * sqrt_n,
        "S": inst.S,
        "2S": 2 * inst.S,
        "N/S_floor": inst.N // inst.S,
    }

    return candidates


# ======================================================================================================================
# INFORMATION-CHANNEL AUDIT
# ======================================================================================================================

def audit_instance(inst: Instance) -> dict[str, bool | int]:

    x_from_d0 = recover_x_from_d_d0(inst)

    return {
        "d_relation": check_d_relation(inst),
        "x_relation": check_x_relation(inst),
        "recovered_x_equals_true": (
            x_from_d0 == inst.x_true
        ),
        "x_recovered_exactly": x_from_d0,
    }


# ======================================================================================================================
# NEGATIVE CONTROL
# ======================================================================================================================
#
# Keep N,S,d fixed.
#
# Replace d0 with another value.
#
# Then x = d-d0 changes immediately.
#
# This demonstrates that N,S,d alone do NOT determine the specific x
# used by the construction.
#
# No p,q factorization is needed for this argument.
#
# ======================================================================================================================

def negative_control(
    inst: Instance,
    alternate_x: int,
) -> tuple[int, int]:

    alternate_d0 = inst.d - alternate_x

    recovered_x = inst.d - alternate_d0

    return alternate_d0, recovered_x


# ======================================================================================================================
# REQUIRED-RADIUS TEST
# ======================================================================================================================

def smallest_ns_radius(inst: Instance) -> int:

    return abs(inst.x_true)


def ns_radius_contains_true(
    inst: Instance,
    radius: int,
) -> bool:

    return abs(inst.x_true) <= radius


# ======================================================================================================================
# MAIN
# ======================================================================================================================

def main() -> None:

    instances = build_instances()

    print("=" * 120)
    print("EXPERIMENT 447")
    print("=" * 120)
    print()
    print("D-D0 / X INFORMATION-FLOW AUDIT")
    print()
    print("QUESTION")
    print(
        "  Is x an independent reconstruction variable, or is"
    )
    print(
        "  x already determined exactly by d and d0?"
    )
    print()
    print("CORE IDENTITY")
    print("  x + d0 = d")
    print("  x = d - d0")
    print()
    print("RULES")
    print("  exact integer arithmetic only")
    print("  no factorization")
    print("  no giant K enumeration")
    print("  no giant X enumeration")
    print("  no CRT Cartesian product")
    print()

    total = 0
    exact_recoveries = 0
    ns_only_failures = 0
    d_identity_failures = 0
    negative_control_failures = 0

    for inst in instances:

        total += 1

        audit = audit_instance(inst)

        if not audit["d_relation"]:
            d_identity_failures += 1

        if audit["recovered_x_equals_true"]:
            exact_recoveries += 1

        print("-" * 120)
        print(
            f"INSTANCE {inst.idx}: "
            f"N={inst.N} "
            f"S={inst.S} "
            f"d={inst.d} "
            f"d0={inst.d0} "
            f"x_true={inst.x_true}"
        )

        print()
        print("  EXACT DATA RELATIONS")

        print(
            f"    d = N - 2S + 1          : "
            f"{audit['d_relation']}"
        )

        print(
            f"    x_true + d0 = d         : "
            f"{audit['x_relation']}"
        )

        print(
            f"    d - d0                 = "
            f"{audit['x_recovered_exactly']}"
        )

        print(
            f"    recovered x == x_true : "
            f"{audit['recovered_x_equals_true']}"
        )

        print()

        # --------------------------------------------------------------------------------------------------------------
        # N,S-only audit
        # --------------------------------------------------------------------------------------------------------------

        print("  N,S-ONLY RADIUS AUDIT")
        print()
        print(
            "    radius model"
            "                radius"
            "       contains x"
            "       excess"
        )
        print(
            "    " + "-" * 78
        )

        for name, radius in candidate_radii_NS(inst).items():

            contains = ns_radius_contains_true(
                inst,
                radius,
            )

            excess = radius - abs(inst.x_true)

            print(
                f"    {name:<16}"
                f"{radius:>20}"
                f"{str(contains):>18}"
                f"{excess:>18}"
            )

            if not contains:
                ns_only_failures += 1

        print()

        # --------------------------------------------------------------------------------------------------------------
        # Exact minimum radius after d0 is available
        # --------------------------------------------------------------------------------------------------------------

        exact_radius = smallest_ns_radius(inst)

        print("  EXACT RADIUS AFTER d0 IS KNOWN")
        print(
            f"    |x_true| = |d-d0| = {exact_radius}"
        )

        print(
            "    This is an exact value, not an upper bound."
        )

        print()

        # --------------------------------------------------------------------------------------------------------------
        # Negative controls
        # --------------------------------------------------------------------------------------------------------------

        print("  NEGATIVE CONTROL: SAME N,S,d, DIFFERENT d0")

        alternate_x = inst.x_true + 12345

        alternate_d0, recovered_x = negative_control(
            inst,
            alternate_x,
        )

        print(
            f"    original d0       = {inst.d0}"
        )

        print(
            f"    alternate d0      = {alternate_d0}"
        )

        print(
            f"    alternate x       = {alternate_x}"
        )

        print(
            f"    d-alternate_d0    = {recovered_x}"
        )

        print(
            f"    recovery correct  = "
            f"{recovered_x == alternate_x}"
        )

        if recovered_x != alternate_x:
            negative_control_failures += 1

        print()

        # --------------------------------------------------------------------------------------------------------------
        # Information separation
        # --------------------------------------------------------------------------------------------------------------

        print("  INFORMATION SEPARATION")

        print(
            "    N,S,d only determine:"
        )

        print(
            "      d = N - 2S + 1"
        )

        print(
            "    N,S,d,d0 determine:"
        )

        print(
            "      x = d - d0"
        )

        print(
            "    Therefore x is not an independent unknown once d0"
        )

        print(
            "    is included among the known reconstruction inputs."
        )

        print()

    # ==================================================================================================================
    # GLOBAL SUMMARY
    # ==================================================================================================================

    print("=" * 120)
    print("GLOBAL INFORMATION-FLOW SUMMARY")
    print("=" * 120)

    print(
        f"  instances                         = {total}"
    )

    print(
        f"  exact d relation failures        = {d_identity_failures}"
    )

    print(
        f"  exact x recoveries from d-d0      = "
        f"{exact_recoveries}/{total}"
    )

    print(
        f"  N,S-only radius containment fails = "
        f"{ns_only_failures}"
    )

    print(
        f"  negative-control failures         = "
        f"{negative_control_failures}"
    )

    print()
    print("CORE STRUCTURAL RESULT")
    print()
    print(
        "  Given d and d0:"
    )

    print(
        "      x = d - d0"
    )

    print()
    print(
        "  Therefore any experiment that treats x as hidden while"
    )

    print(
        "  simultaneously supplying d0 as known is not performing"
    )

    print(
        "  an x-recovery experiment."
    )

    print()
    print(
        "  The actual unresolved question is whether d0 itself can"
    )

    print(
        "  be obtained independently from the public/known data."
    )

    print()
    print("N,S,d-ONLY CONTROL")
    print()
    print(
        "  Holding N,S,d fixed while changing d0 changes x exactly."
    )

    print(
        "  Thus N,S,d alone do not identify the constructed x."
    )

    print()
    print("CONNECTION TO EXPERIMENT 446R2")
    print()
    print(
        "  Experiment 446R2 used d0 inside"
    )

    print(
        "      K = A/4 - ((x+d0)/2)^2"
    )

    print(
        "  and also treated x as hidden."
    )

    print(
        "  But with d0 known, the separate identity"
    )

    print(
        "      x = d-d0"
    )

    print(
        "  must be accounted for before interpreting the K pullback"
    )

    print(
        "  as an x-reconstruction mechanism."
    )

    print()
    print("NEXT INFORMATION QUESTION")
    print()
    print(
        "  The mathematically meaningful next audit is therefore:"
    )

    print(
        "      Can d0 be derived from N,S and the observed r/K relation"
    )

    print(
        "      without already knowing x?"
    )

    print()
    print(
        "  That is a different problem from inverting a bounded K interval."
    )

    print()
    print("=" * 120)
    print("EXPERIMENT 447 FINAL STATUS")
    print("=" * 120)

    print(
        f"  EXACT x=d-d0 RECOVERY = "
        f"{exact_recoveries}/{total}"
    )

    print(
        f"  N,S-ONLY TEST FAILURES = "
        f"{ns_only_failures}"
    )

    print(
        f"  NEGATIVE CONTROL FAILURES = "
        f"{negative_control_failures}"
    )

    print(
        "  GIANT K ENUMERATION = False"
    )

    print(
        "  GIANT X ENUMERATION = False"
    )

    print(
        "  CRT CARTESIAN PRODUCT = False"
    )

    print(
        "  INTEGER-EXACT = True"
    )

    print()
    print(
        "  CONCLUSION:"
    )

    print(
        "    With d0 supplied, x is exactly known as d-d0."
    )

    print(
        "    The remaining research target is recovery of d0 itself."
    )

    print("=" * 120)
    print("EXPERIMENT 447 FINISHED")
    print("=" * 120)


if __name__ == "__main__":
    main()
