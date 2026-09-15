#!/usr/bin/env python3

"""
======================================================================
2^z RECURSIVE PARAMETERIZATION EXPERIMENT
======================================================================

For a factor pair

    n = p*q

and

    k_z = 2^(z-1),

define

    p = k_z (y_z - x_z) + a_z
    q = k_z (y_z + x_z) - a_z

Equivalently:

    a_z = p mod 2^(z-1)

    y_z = (p+q) / 2^z

    x_z = (q-p+2*a_z) / 2^z

whenever x_z,y_z are integers.

This experiment investigates transitions between levels:

    z -> z+1

and tests candidate recurrences for

    a
    x
    y

including:

    x' = floor(x/2)
    x' = ceil(x/2)
    x' = floor(x/2)+1
    x' = (x+1)/2

and analogous rules for y.

It also checks the exact algebraic transition.

======================================================================
"""

from __future__ import annotations

from dataclasses import dataclass


# ======================================================================
# TEST CASES
# ======================================================================

CASES = [
    # p, q
    (59, 101),
    (53, 61),
    (101, 103),
    (97, 193),
    (223, 449),
]


# ======================================================================
# DATA
# ======================================================================

@dataclass
class Level:
    z: int
    modulus: int
    a: int
    x: int
    y: int
    p: int
    q: int


# ======================================================================
# EXACT LEVEL
# ======================================================================

def get_level(
    p: int,
    q: int,
    z: int,
) -> Level | None:

    modulus = 1 << z
    a_modulus = 1 << (z - 1)

    # Low residue of p.
    a = p % a_modulus

    # Reconstruction.
    y_num = p + q
    x_num = q - p + 2 * a

    if y_num % modulus != 0:
        return None

    if x_num % modulus != 0:
        return None

    x = x_num // modulus
    y = y_num // modulus

    # Verify original equations.
    k = 1 << (z - 1)

    p_check = k * (y - x) + a
    q_check = k * (y + x) - a

    if p_check != p or q_check != q:
        raise AssertionError(
            f"Reconstruction failed at z={z}"
        )

    return Level(
        z=z,
        modulus=modulus,
        a=a,
        x=x,
        y=y,
        p=p,
        q=q,
    )


# ======================================================================
# COLLECT VALID LEVELS
# ======================================================================

def collect_levels(
    p: int,
    q: int,
    max_z: int = 32,
) -> list[Level]:

    levels = []

    for z in range(2, max_z + 1):

        level = get_level(
            p,
            q,
            z,
        )

        if level is not None:
            levels.append(level)

    return levels


# ======================================================================
# PRINT LEVEL TABLE
# ======================================================================

def print_levels(
    p: int,
    q: int,
    levels: list[Level],
):

    n = p * q

    print()
    print("=" * 110)
    print(f"N = {n} = {p} × {q}")
    print("=" * 110)

    print()
    print(
        f"{'z':>3} "
        f"{'2^z':>10} "
        f"{'a_z':>12} "
        f"{'x_z':>12} "
        f"{'y_z':>12} "
        f"{'p mod 2^(z-1)':>18}"
    )

    print("-" * 85)

    for L in levels:

        print(
            f"{L.z:3d} "
            f"{L.modulus:10d} "
            f"{L.a:12d} "
            f"{L.x:12d} "
            f"{L.y:12d} "
            f"{L.a:18d}"
        )


# ======================================================================
# TEST X RECURRENCES
# ======================================================================

def test_x_recursions(
    levels: list[Level],
):

    print()
    print("=" * 110)
    print("X-RECURRENCE TEST")
    print("=" * 110)

    rules = {
        "floor(x/2)": lambda x: x // 2,

        "ceil(x/2)": lambda x: (x + 1) // 2,

        "floor(x/2)+1": lambda x: x // 2 + 1,

        "(x+1)/2": lambda x: (
            (x + 1) // 2
        ),

        "floor((x+1)/2)": lambda x: (
            (x + 1) // 2
        ),
    }

    for name, rule in rules.items():

        successes = 0
        failures = 0

        print()
        print(f"RULE: x' = {name}")
        print("-" * 80)

        for i in range(len(levels) - 1):

            L = levels[i]
            N = levels[i + 1]

            predicted = rule(L.x)
            actual = N.x

            ok = predicted == actual

            if ok:
                successes += 1
            else:
                failures += 1

            print(
                f"z={L.z:2d}->{N.z:2d} "
                f"x={L.x:<8d} "
                f"pred={predicted:<8d} "
                f"actual={actual:<8d} "
                f"{'PASS' if ok else 'FAIL'}"
            )

        print(
            f"RESULT: "
            f"{successes} PASS, "
            f"{failures} FAIL"
        )


# ======================================================================
# TEST Y RECURRENCES
# ======================================================================

def test_y_recursions(
    levels: list[Level],
):

    print()
    print("=" * 110)
    print("Y-RECURRENCE TEST")
    print("=" * 110)

    rules = {
        "floor(y/2)": lambda y: y // 2,

        "ceil(y/2)": lambda y: (y + 1) // 2,

        "y/2": lambda y: (
            y // 2
        ),

        "floor(y/2)+1": lambda y: (
            y // 2 + 1
        ),
    }

    for name, rule in rules.items():

        successes = 0
        failures = 0

        print()
        print(f"RULE: y' = {name}")
        print("-" * 80)

        for i in range(len(levels) - 1):

            L = levels[i]
            N = levels[i + 1]

            predicted = rule(L.y)
            actual = N.y

            ok = predicted == actual

            if ok:
                successes += 1
            else:
                failures += 1

            print(
                f"z={L.z:2d}->{N.z:2d} "
                f"y={L.y:<8d} "
                f"pred={predicted:<8d} "
                f"actual={actual:<8d} "
                f"{'PASS' if ok else 'FAIL'}"
            )

        print(
            f"RESULT: "
            f"{successes} PASS, "
            f"{failures} FAIL"
        )


# ======================================================================
# TEST a TRANSITION
# ======================================================================

def test_a_recursion(
    levels: list[Level],
):

    print()
    print("=" * 110)
    print("a-RECURRENCE TEST")
    print("=" * 110)

    for i in range(len(levels) - 1):

        L = levels[i]
        N = levels[i + 1]

        # Since:
        #
        # a_z = p mod 2^(z-1)
        #
        # a_{z+1} is either:
        #
        # a_z
        #
        # or
        #
        # a_z + 2^(z-1)
        #
        # depending on the next bit of p.

        low_bit = (
            p_bit(
                L.p,
                L.z - 1,
            )
        )

        predicted = (
            L.a
            + low_bit * L.modulus // 2
        )

        actual = N.a

        print(
            f"z={L.z:2d}->{N.z:2d} "
            f"a={L.a:<8d} "
            f"next_bit={low_bit} "
            f"pred={predicted:<8d} "
            f"actual={actual:<8d} "
            f"{'PASS' if predicted == actual else 'FAIL'}"
        )


# ======================================================================
# p BIT
# ======================================================================

def p_bit(
    p: int,
    bit: int,
) -> int:

    return (
        p >> bit
    ) & 1


# ======================================================================
# EXACT TRANSITION DERIVATION
# ======================================================================

def exact_transition(
    L: Level,
    N: Level,
):

    print()
    print(
        f"z={L.z} -> z={N.z}"
    )

    print(
        f"  a: {L.a} -> {N.a}"
    )

    print(
        f"  x: {L.x} -> {N.x}"
    )

    print(
        f"  y: {L.y} -> {N.y}"
    )

    # ----------------------------------------------------------
    # Derive from equations.
    # ----------------------------------------------------------

    #
    # a_{z+1} = a_z + b * 2^(z-1)
    #

    b = (
        N.a - L.a
    ) // L.modulus * 2

    # More directly:
    #
    # N.a - L.a is either 0 or 2^(z-1).
    #

    delta_a = N.a - L.a
    half_modulus = L.modulus // 2

    next_bit = (
        1
        if delta_a == half_modulus
        else 0
    )

    print(
        f"  Δa = {delta_a}"
    )

    print(
        f"  next factor bit = {next_bit}"
    )

    # ----------------------------------------------------------
    # Exact x transition.
    # ----------------------------------------------------------

    #
    # x_z = (q-p+2a_z)/2^z
    #
    # Therefore:
    #
    # x_{z+1}
    # =
    # (x_z + next_bit)/2
    #

    numerator = (
        L.x + next_bit
    )

    print(
        f"  exact x formula:"
    )

    print(
        f"      x' = ({L.x} + {next_bit}) / 2"
    )

    if numerator % 2 == 0:

        predicted_x = numerator // 2

        print(
            f"          = {predicted_x}"
        )

    else:

        print(
            f"          = {numerator}/2 "
            f"(non-integer)"
        )

    # ----------------------------------------------------------
    # Exact y transition.
    # ----------------------------------------------------------

    #
    # y_z = (p+q)/2^z
    #
    # Therefore:
    #
    # y_{z+1}=y_z/2
    #

    print(
        f"  exact y formula:"
    )

    print(
        f"      y' = {L.y}/2"
    )

    if L.y % 2 == 0:

        print(
            f"         = {L.y // 2}"
        )

    else:

        print(
            f"         = non-integer"
        )


# ======================================================================
# ALL TRANSITIONS
# ======================================================================

def show_transitions(
    levels: list[Level],
):

    print()
    print("=" * 110)
    print("EXACT LEVEL TRANSITIONS")
    print("=" * 110)

    for i in range(
        len(levels) - 1
    ):

        exact_transition(
            levels[i],
            levels[i + 1],
        )


# ======================================================================
# VERIFY DIRECTLY FROM p,q
# ======================================================================

def verify_all_levels(
    p: int,
    q: int,
    levels: list[Level],
):

    print()
    print("=" * 110)
    print("DIRECT VERIFICATION")
    print("=" * 110)

    failures = 0

    for L in levels:

        k = 1 << (L.z - 1)

        p_check = (
            k * (L.y - L.x)
            + L.a
        )

        q_check = (
            k * (L.y + L.x)
            - L.a
        )

        product = (
            p_check * q_check
        )

        ok = (
            p_check == p
            and q_check == q
            and product == p * q
        )

        print(
            f"z={L.z:2d} "
            f"p_check={p_check:<10d} "
            f"q_check={q_check:<10d} "
            f"{'PASS' if ok else 'FAIL'}"
        )

        if not ok:
            failures += 1

    print()
    print(
        f"TOTAL FAILURES = {failures}"
    )


# ======================================================================
# ONE CASE
# ======================================================================

def experiment_case(
    p: int,
    q: int,
):

    levels = collect_levels(
        p,
        q,
        max_z=32,
    )

    if not levels:

        print(
            f"No levels for {p} × {q}"
        )

        return

    print_levels(
        p,
        q,
        levels,
    )

    test_a_recursion(
        levels,
    )

    test_x_recursions(
        levels,
    )

    test_y_recursions(
        levels,
    )

    show_transitions(
        levels,
    )

    verify_all_levels(
        p,
        q,
        levels,
    )


# ======================================================================
# MAIN
# ======================================================================

def main():

    print("=" * 110)
    print(
        "2^z RECURSIVE COORDINATE EXPERIMENT"
    )
    print("=" * 110)

    for p, q in CASES:

        experiment_case(
            p,
            q,
        )


if __name__ == "__main__":
    main()
