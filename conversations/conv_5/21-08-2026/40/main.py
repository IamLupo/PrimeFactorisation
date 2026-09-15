#!/usr/bin/env python3

# ==============================================================================
# EXPERIMENT 599
# ==============================================================================
# EXACT 2-ADIC CHILD-RESIDUE LAW
#
# Goal
# ----
#
# Stop searching for arbitrary residue recurrences.
#
# Derive the child residue directly from the A/B factor functions.
#
# Let:
#
#     k = 2^(z-1)
#
# FRAME A:
#
#     p = k(y-x)+3
#     q = k(y+x)-3
#
# FRAME B:
#
#     p = (k(x+y)-3)/3
#     q = k(y-x)+3
#
# Then:
#
#     n_A
#       = k^2(y^2-x^2) + 6kx - 9
#
#     n_B
#       = [k^2(y^2-x^2) + 6kx - 9] / 3
#
# modulo:
#
#     2^(z+1) = 4k
#
# the child residue depends only on x mod 2.
#
# Predicted:
#
#     A:
#         child_residue =
#             -9 + 2^z*(x mod 2)
#             mod 2^(z+1)
#
#     B:
#         child_residue =
#             -3 + 2^z*(x mod 2)
#             mod 2^(z+1)
#
# For a valid coordinate renormalization:
#
#     x' = x/2
#     y' = y/2
#
# therefore x must be even.
#
# The experiment tests all of this against the generated semiprimes.
#
# NO FILES
# NO EXTERNAL SOURCES
# ==============================================================================

from __future__ import annotations

from dataclasses import dataclass
from collections import Counter, defaultdict
from math import isqrt


# ==============================================================================
# CONFIGURATION
# ==============================================================================

MAX_N = 2_000_000

MIN_Z = 2
MAX_Z = 12

SHOW_MISMATCHES = 10
SHOW_EXAMPLES = 8


# ==============================================================================
# STATE
# ==============================================================================

@dataclass(frozen=True)
class State:

    n: int
    p: int
    q: int

    z: int
    residue: int

    frame: str

    x: int
    y: int


# ==============================================================================
# PRIME SIEVE
# ==============================================================================

def prime_sieve(limit: int):

    sieve = bytearray(
        b"\x01"
    ) * (limit + 1)

    sieve[0] = 0
    sieve[1] = 0

    root = isqrt(limit)

    for p in range(
        2,
        root + 1,
    ):

        if not sieve[p]:
            continue

        start = p * p

        count = (
            (limit - start) // p
        ) + 1

        sieve[
            start:
            limit + 1:
            p
        ] = b"\x00" * count

    return sieve


# ==============================================================================
# SEMIPRIMES
# ==============================================================================

def generate_semiprimes(
    limit,
    sieve,
):

    primes = [
        p
        for p in range(
            3,
            limit + 1,
            2,
        )
        if sieve[p]
    ]

    result = []

    for i, p in enumerate(primes):

        if p * p > limit:
            break

        max_q = limit // p

        for q in primes[i:]:

            if q > max_q:
                break

            result.append(
                (
                    p * q,
                    p,
                    q,
                )
            )

    return result


# ==============================================================================
# FRAME A
# ==============================================================================

def frame_A(
    p,
    q,
    z,
):

    k = 1 << (z - 1)
    d = 2 * k

    x_num = q - p + 6
    y_num = p + q

    if x_num % d != 0:
        return None

    if y_num % d != 0:
        return None

    return (
        x_num // d,
        y_num // d,
    )


# ==============================================================================
# FRAME B
# ==============================================================================

def frame_B(
    p,
    q,
    z,
):

    k = 1 << (z - 1)
    d = 2 * k

    # Correct B frame:
    #
    #     p = (k(x+y)-3)/3
    #     q = k(y-x)+3
    #
    # Therefore:
    #
    #     x = (3p-q+6)/(2k)
    #     y = (3p+q)/(2k)

    x_num = (
        3 * p
        - q
        + 6
    )

    y_num = (
        3 * p
        + q
    )

    if x_num % d != 0:
        return None

    if y_num % d != 0:
        return None

    return (
        x_num // d,
        y_num // d,
    )


# ==============================================================================
# BUILD REPRESENTATIONS
# ==============================================================================

def build_levels(
    semiprimes,
):

    levels = {
        z: {}
        for z in range(
            MIN_Z,
            MAX_Z + 1,
        )
    }

    for n, p, q in semiprimes:

        for z in range(
            MIN_Z,
            MAX_Z + 1,
        ):

            modulus = 1 << z

            residue = n % modulus

            a = frame_A(
                p,
                q,
                z,
            )

            b = frame_B(
                p,
                q,
                z,
            )

            if a is not None:

                levels[z][n] = State(
                    n=n,
                    p=p,
                    q=q,
                    z=z,
                    residue=residue,
                    frame="A",
                    x=a[0],
                    y=a[1],
                )

            elif b is not None:

                levels[z][n] = State(
                    n=n,
                    p=p,
                    q=q,
                    z=z,
                    residue=residue,
                    frame="B",
                    x=b[0],
                    y=b[1],
                )

    return levels


# ==============================================================================
# SYMBOLIC NUMERATOR
# ==============================================================================

def symbolic_n_A(
    x,
    y,
    z,
):

    k = 1 << (z - 1)

    p = (
        k * (y - x)
        + 3
    )

    q = (
        k * (y + x)
        - 3
    )

    return p * q


def symbolic_n_B(
    x,
    y,
    z,
):

    k = 1 << (z - 1)

    numerator = (
        k * (x + y)
        - 3
    )

    q = (
        k * (y - x)
        + 3
    )

    if numerator % 3 != 0:
        return None

    p = numerator // 3

    return p * q


# ==============================================================================
# SYMBOLIC EXPANSION CHECK
# ==============================================================================

def expansion_A(
    x,
    y,
    z,
):

    k = 1 << (z - 1)

    return (
        k * k * (y * y - x * x)
        + 6 * k * x
        - 9
    )


def expansion_B(
    x,
    y,
    z,
):

    numerator = expansion_A(
        x,
        y,
        z,
    )

    if numerator % 3 != 0:
        return None

    return numerator // 3


# ==============================================================================
# PREDICTED CHILD RESIDUE
# ==============================================================================

def predicted_child_residue(
    state,
):

    z = state.z

    child_modulus = 1 << (z + 1)

    parity = state.x & 1

    if state.frame == "A":

        return (
            -9
            + (1 << z) * parity
        ) % child_modulus

    if state.frame == "B":

        return (
            -3
            + (1 << z) * parity
        ) % child_modulus

    raise ValueError(
        state.frame
    )


# ==============================================================================
# PREDICTED CHILD FRAME
# ==============================================================================

def predicted_child_frame(
    state,
):

    #
    # If x is even:
    #
    #     x' = x/2
    #     y' = y/2
    #
    # and the same frame survives.
    #
    # If x is odd, this particular coordinate frame cannot
    # perform an integral half-coordinate transformation.
    #

    if state.x % 2 != 0:
        return None

    if state.y % 2 != 0:
        return None

    return state.frame


# ==============================================================================
# SAME-N CHILD CHECK
# ==============================================================================

def child_state(
    levels,
    state,
):

    child = levels.get(
        state.z + 1,
        {}
    )

    return child.get(
        state.n
    )


# ==============================================================================
# [1] SYMBOLIC EXPANSION AUDIT
# ==============================================================================

def symbolic_expansion_audit():

    print()
    print("=" * 90)
    print(
        "SYMBOLIC EXPANSION AUDIT"
    )
    print("=" * 90)

    failures_A = 0
    failures_B = 0

    # Test many small integer coordinates.
    for z in range(
        2,
        8,
    ):

        for x in range(
            -20,
            21,
        ):

            for y in range(
                -20,
                21,
            ):

                a1 = symbolic_n_A(
                    x,
                    y,
                    z,
                )

                a2 = expansion_A(
                    x,
                    y,
                    z,
                )

                if a1 != a2:
                    failures_A += 1

                b1 = symbolic_n_B(
                    x,
                    y,
                    z,
                )

                b2 = expansion_B(
                    x,
                    y,
                    z,
                )

                if b1 != b2:
                    failures_B += 1

    print(
        f"    A expansion failures = "
        f"{failures_A}"
    )

    print(
        f"    B expansion failures = "
        f"{failures_B}"
    )


# ==============================================================================
# [2] CHILD RESIDUE LAW AUDIT
# ==============================================================================

def child_residue_audit(
    levels,
):

    print()
    print("=" * 90)
    print(
        "EXACT CHILD RESIDUE LAW"
    )
    print("=" * 90)

    total = 0
    failures = 0

    per_level = {}

    examples = []

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        checked = 0
        failed = 0

        parent = levels[z]
        child = levels[z + 1]

        for n, state in parent.items():

            actual_child = child.get(
                n
            )

            predicted = (
                predicted_child_residue(
                    state
                )
            )

            actual = (
                actual_child.residue
                if actual_child is not None
                else n % (1 << (z + 1))
            )

            checked += 1
            total += 1

            if predicted != actual:

                failed += 1
                failures += 1

                if (
                    len(examples)
                    < SHOW_MISMATCHES
                ):

                    examples.append(
                        (
                            z,
                            state,
                            predicted,
                            actual,
                            actual_child,
                        )
                    )

        per_level[z] = (
            checked,
            failed,
        )

    for z in sorted(per_level):

        checked, failed = (
            per_level[z]
        )

        print(
            f"z={z} -> z={z+1}"
            f"  checked={checked}"
            f"  failures={failed}"
        )

    print()
    print(
        f"GLOBAL checked={total}"
    )

    print(
        f"GLOBAL failures={failures}"
    )

    if examples:

        print()
        print(
            "FIRST MISMATCHES"
        )

        for (
            z,
            state,
            predicted,
            actual,
            actual_child,
        ) in examples:

            print()
            print(
                f"    z={z}"
            )

            print(
                f"    n={state.n}"
            )

            print(
                f"    frame={state.frame}"
            )

            print(
                f"    residue={state.residue}"
            )

            print(
                f"    x={state.x}"
            )

            print(
                f"    y={state.y}"
            )

            print(
                f"    xbit={state.x & 1}"
            )

            print(
                f"    predicted={predicted}"
            )

            print(
                f"    actual={actual}"
            )

            if actual_child is not None:

                print(
                    f"    child="
                    f"{actual_child.frame}"
                    f"{actual_child.residue}"
                )

            else:

                print(
                    "    child representation=None"
                )


# ==============================================================================
# [3] SURVIVAL CONDITION
# ==============================================================================

def survival_audit(
    levels,
):

    print()
    print("=" * 90)
    print(
        "SURVIVAL CONDITION: x,y EVEN"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        parent = levels[z]
        child = levels[z + 1]

        total = 0
        even_xy = 0
        odd_x = 0
        odd_y = 0

        even_and_child = 0
        odd_and_child = 0

        mismatches = 0

        for n, state in parent.items():

            actual_child = child.get(
                n
            )

            total += 1

            x_even = (
                state.x % 2 == 0
            )

            y_even = (
                state.y % 2 == 0
            )

            if x_even and y_even:

                even_xy += 1

                if actual_child is not None:
                    even_and_child += 1

            else:

                if not x_even:
                    odd_x += 1

                if not y_even:
                    odd_y += 1

                if actual_child is not None:
                    odd_and_child += 1

            expected_child = (
                x_even and y_even
            )

            actual_exists = (
                actual_child is not None
            )

            if expected_child != actual_exists:

                mismatches += 1

        print()
        print(
            f"z={z} -> z={z+1}"
        )

        print(
            f"    total={total}"
        )

        print(
            f"    x,y even={even_xy}"
        )

        print(
            f"    odd x={odd_x}"
        )

        print(
            f"    odd y={odd_y}"
        )

        print(
            f"    even(x,y) with child="
            f"{even_and_child}"
        )

        print(
            f"    odd(x/y) with child="
            f"{odd_and_child}"
        )

        print(
            f"    survival mismatches="
            f"{mismatches}"
        )


# ==============================================================================
# [4] SAME-FRAME HALVING
# ==============================================================================

def halving_audit(
    levels,
):

    print()
    print("=" * 90)
    print(
        "SAME-FRAME HALVING"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        parent = levels[z]
        child = levels[z + 1]

        common = set(parent).intersection(
            child
        )

        x_fail = 0
        y_fail = 0
        frame_fail = 0

        for n in common:

            s0 = parent[n]
            s1 = child[n]

            if s1.x * 2 != s0.x:
                x_fail += 1

            if s1.y * 2 != s0.y:
                y_fail += 1

            if s1.frame != s0.frame:
                frame_fail += 1

        print(
            f"z={z} -> z={z+1}"
            f"  common={len(common)}"
            f"  x_fail={x_fail}"
            f"  y_fail={y_fail}"
            f"  frame_fail={frame_fail}"
        )


# ==============================================================================
# [5] RESIDUE FORMULAS IN CLOSED FORM
# ==============================================================================

def closed_form_residue(
    frame,
    z,
):

    modulus = 1 << z

    if frame == "A":
        return (-9) % modulus

    if frame == "B":
        return (-3) % modulus

    raise ValueError(
        frame
    )


def residue_formula_audit(
    levels,
):

    print()
    print("=" * 90)
    print(
        "CLOSED-FORM ACTIVE RESIDUES"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        by_frame = {
            "A": Counter(),
            "B": Counter(),
        }

        for state in levels[z].values():

            by_frame[
                state.frame
            ][
                state.residue
            ] += 1

        print()
        print(
            f"z={z}"
        )

        for frame in (
            "A",
            "B",
        ):

            counts = by_frame[frame]

            predicted = closed_form_residue(
                frame,
                z,
            )

            print(
                f"    frame={frame}"
                f"  predicted={predicted}"
                f"  observed={dict(counts)}"
            )


# ==============================================================================
# [6] CHILD RESIDUE SYMBOLIC TABLE
# ==============================================================================

def symbolic_child_table():

    print()
    print("=" * 90)
    print(
        "SYMBOLIC CHILD RESIDUE TABLE"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        modulus = 1 << (
            z + 1
        )

        shift = 1 << z

        print()
        print(
            f"z={z}"
        )

        for frame in (
            "A",
            "B",
        ):

            even = (
                (
                    -9
                    if frame == "A"
                    else -3
                )
                % modulus
            )

            odd = (
                (
                    -9
                    if frame == "A"
                    else -3
                )
                + shift
            ) % modulus

            print(
                f"    "
                f"{frame}: "
                f"x even -> {even:<5} "
                f"x odd -> {odd:<5}"
            )


# ==============================================================================
# [7] FRAME-PRESERVING CHILD PREDICTION
# ==============================================================================

def frame_preservation_audit(
    levels,
):

    print()
    print("=" * 90)
    print(
        "FRAME-PRESERVING CHILD PREDICTION"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        parent = levels[z]
        child = levels[z + 1]

        total = 0
        failures = 0

        for n, state in parent.items():

            if (
                state.x % 2
                != 0
                or
                state.y % 2
                != 0
            ):
                continue

            total += 1

            actual = child.get(
                n
            )

            if actual is None:

                failures += 1
                continue

            if actual.frame != state.frame:

                failures += 1

        print(
            f"z={z} -> z={z+1}"
            f"  eligible={total}"
            f"  failures={failures}"
        )


# ==============================================================================
# [8] EXAMPLES
# ==============================================================================

def examples(
    levels,
):

    print()
    print("=" * 90)
    print(
        "EXAMPLES"
    )
    print("=" * 90)

    for z in range(
        MIN_Z,
        MAX_Z,
    ):

        parent = levels[z]
        child = levels[z + 1]

        shown = 0

        print()
        print(
            f"z={z} -> z={z+1}"
        )

        for n in sorted(parent):

            if shown >= SHOW_EXAMPLES:
                break

            state = parent[n]

            actual = child.get(
                n
            )

            predicted = (
                predicted_child_residue(
                    state
                )
            )

            print(
                f"    "
                f"n={n} "
                f"{state.frame}{state.residue}"
                f" "
                f"x={state.x}"
                f" y={state.y}"
                f" "
                f"xbit={state.x & 1}"
                f" "
                f"pred={predicted}"
                f" "
                f"actual="
                f"{actual.residue if actual else None}"
            )

            shown += 1


# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

def final_summary():

    print()
    print("=" * 90)
    print(
        "FINAL STRUCTURAL SUMMARY"
    )
    print("=" * 90)

    print(
        r"""
For k = 2^(z-1):

FRAME A

    n
      =
    (k(y-x)+3)(k(y+x)-3)

      =
    k^2(y^2-x^2)
      + 6kx
      - 9.

Therefore modulo 2^(z+1):

    n
      =
    -9 + 2^z (x mod 2)
      mod 2^(z+1).


FRAME B

    n
      =
    [(k(x+y)-3)(k(y-x)+3)] / 3

      =
    [k^2(y^2-x^2)
      + 6kx
      - 9] / 3.

Since 3 is invertible modulo powers of 2:

    n
      =
    -3 + 2^z (x mod 2)
      mod 2^(z+1).


Thus:

    A:
        x even -> child residue = -9
        x odd  -> child residue = -9 + 2^z

    B:
        x even -> child residue = -3
        x odd  -> child residue = -3 + 2^z.

The coordinate renormalization requires:

    x' = x/2
    y' = y/2,

so the surviving child has:

    x even
    y even

and therefore:

    A -> A
    B -> B

with the deterministic active residues:

    A_z = -9 mod 2^z
    B_z = -3 mod 2^z.

This is the exact finite-state hierarchy we are testing.
"""
    )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print(
        "EXPERIMENT 599 START"
    )
    print("=" * 90)

    print()
    print(
        "EXACT 2-ADIC CHILD-RESIDUE LAW"
    )

    # --------------------------------------------------------------------------
    # Prime sieve
    # --------------------------------------------------------------------------

    print()
    print(
        "[1] PRIME SIEVE"
    )

    sieve = prime_sieve(
        MAX_N
    )

    # --------------------------------------------------------------------------
    # Semiprimes
    # --------------------------------------------------------------------------

    print()
    print(
        "[2] SEMIPRIME GENERATION"
    )

    semiprimes = generate_semiprimes(
        MAX_N,
        sieve,
    )

    print(
        f"    semiprimes="
        f"{len(semiprimes)}"
    )

    # --------------------------------------------------------------------------
    # Levels
    # --------------------------------------------------------------------------

    print()
    print(
        "[3] BUILD LEVELS"
    )

    levels = build_levels(
        semiprimes
    )

    for z in range(
        MIN_Z,
        MAX_Z + 1,
    ):

        print(
            f"    z={z}"
            f" states={len(levels[z])}"
        )

    # --------------------------------------------------------------------------
    # Symbolic algebra
    # --------------------------------------------------------------------------

    symbolic_expansion_audit()

    # --------------------------------------------------------------------------
    # Main child law
    # --------------------------------------------------------------------------

    child_residue_audit(
        levels
    )

    # --------------------------------------------------------------------------
    # Survival
    # --------------------------------------------------------------------------

    survival_audit(
        levels
    )

    # --------------------------------------------------------------------------
    # Same-frame coordinate recursion
    # --------------------------------------------------------------------------

    halving_audit(
        levels
    )

    # --------------------------------------------------------------------------
    # Closed-form residue sequence
    # --------------------------------------------------------------------------

    residue_formula_audit(
        levels
    )

    # --------------------------------------------------------------------------
    # Symbolic child table
    # --------------------------------------------------------------------------

    symbolic_child_table()

    # --------------------------------------------------------------------------
    # Frame preservation
    # --------------------------------------------------------------------------

    frame_preservation_audit(
        levels
    )

    # --------------------------------------------------------------------------
    # Examples
    # --------------------------------------------------------------------------

    if SHOW_EXAMPLES:

        examples(
            levels
        )

    # --------------------------------------------------------------------------
    # Summary
    # --------------------------------------------------------------------------

    final_summary()

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 599 FINISHED"
    )
    print("=" * 90)


if __name__ == "__main__":
    main()
