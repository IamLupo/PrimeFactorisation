#!/usr/bin/env python3

"""
======================================================================
2^z RECURSIVE STATE-MACHINE EXPERIMENT
======================================================================

Hypothesis:

    State at level z:

        (a_z, x_z, y_z)

    with

        p = 2^(z-1)(y_z-x_z) + a_z
        q = 2^(z-1)(y_z+x_z) - a_z

Hypothesized transition:

        b_z = x_z mod 2

        a_(z+1) = a_z + b_z * 2^(z-1)

        x_(z+1) = (x_z + b_z) / 2

        y_(z+1) = y_z / 2

The experiment:

    1. Computes the true state from p,q.
    2. Generates the next state recursively WITHOUT p,q.
    3. Compares recursive state against true state.
    4. Checks reconstruction of p and q.
    5. Detects the first level where recursion fails.
    6. Tests all possible initial branches from low-bit states.

======================================================================
"""

from __future__ import annotations

from dataclasses import dataclass


# ======================================================================
# KNOWN FACTOR CASES
# ======================================================================

CASES = [
    (59, 101),
    (101, 103),
    (223, 449),
]


# ======================================================================
# DATA
# ======================================================================

@dataclass
class State:
    z: int
    a: int
    x: int
    y: int


# ======================================================================
# TRUE STATE FROM p,q
# ======================================================================

def true_state(
    p: int,
    q: int,
    z: int,
) -> State | None:

    modulus = 1 << z
    k = 1 << (z - 1)

    a = p % k

    y_num = p + q
    x_num = q - p + 2 * a

    if y_num % modulus != 0:
        return None

    if x_num % modulus != 0:
        return None

    x = x_num // modulus
    y = y_num // modulus

    # Direct verification.
    p_check = k * (y - x) + a
    q_check = k * (y + x) - a

    if p_check != p or q_check != q:
        raise AssertionError(
            "True-state reconstruction failure"
        )

    return State(
        z=z,
        a=a,
        x=x,
        y=y,
    )


# ======================================================================
# RECURSIVE TRANSITION
# ======================================================================

def recursive_transition(
    state: State,
) -> State | None:

    z = state.z

    # --------------------------------------------------------------
    # Proposed next bit.
    # --------------------------------------------------------------

    b = state.x & 1

    # --------------------------------------------------------------
    # New a.
    # --------------------------------------------------------------

    new_a = (
        state.a
        + b * (1 << (z - 1))
    )

    # --------------------------------------------------------------
    # New x.
    # --------------------------------------------------------------

    x_num = state.x + b

    if x_num % 2 != 0:
        return None

    new_x = x_num // 2

    # --------------------------------------------------------------
    # New y.
    # --------------------------------------------------------------

    if state.y % 2 != 0:
        return None

    new_y = state.y // 2

    return State(
        z=z + 1,
        a=new_a,
        x=new_x,
        y=new_y,
    )


# ======================================================================
# VERIFY STATE
# ======================================================================

def verify_state(
    state: State,
    p: int,
    q: int,
) -> bool:

    k = 1 << (state.z - 1)

    p_check = (
        k * (state.y - state.x)
        + state.a
    )

    q_check = (
        k * (state.y + state.x)
        - state.a
    )

    return (
        p_check == p
        and q_check == q
    )


# ======================================================================
# PRINT STATE
# ======================================================================

def print_state(
    state: State,
    p: int,
    q: int,
    prefix: str = "",
):

    valid = verify_state(
        state,
        p,
        q,
    )

    print(
        f"{prefix}"
        f"z={state.z:2d} "
        f"a={state.a:<8d} "
        f"x={state.x:<8d} "
        f"y={state.y:<8d} "
        f"x%2={state.x & 1} "
        f"y%2={state.y & 1} "
        f"valid={valid}"
    )


# ======================================================================
# COMPARE RECURSION
# ======================================================================

def compare_case(
    p: int,
    q: int,
    max_z: int = 32,
):

    print()
    print("=" * 100)
    print(f"N = {p*q} = {p} × {q}")
    print("=" * 100)

    # --------------------------------------------------------------
    # Find starting state.
    # --------------------------------------------------------------

    initial = None

    for z in range(2, max_z + 1):

        state = true_state(
            p,
            q,
            z,
        )

        if state is not None:

            initial = state
            break

    if initial is None:

        print(
            "No valid starting state."
        )

        return

    print()
    print("INITIAL STATE")
    print("-" * 100)

    print_state(
        initial,
        p,
        q,
    )

    # --------------------------------------------------------------
    # Recursive comparison.
    # --------------------------------------------------------------

    recursive = initial

    print()
    print("RECURSIVE TRAJECTORY")
    print("-" * 100)

    for next_z in range(
        initial.z + 1,
        max_z + 1,
    ):

        predicted = recursive_transition(
            recursive
        )

        actual = true_state(
            p,
            q,
            next_z,
        )

        print()

        print(
            f"TRANSITION "
            f"{recursive.z} -> {next_z}"
        )

        # ----------------------------------------------------------
        # Recursive state ceased to exist.
        # ----------------------------------------------------------

        if predicted is None:

            print(
                "  RECURSIVE TRANSITION: INVALID"
            )

            if actual is None:

                print(
                    "  ACTUAL STATE: INVALID"
                )

            else:

                print(
                    "  ACTUAL STATE STILL EXISTS:"
                )

                print_state(
                    actual,
                    p,
                    q,
                    prefix="    ",
                )

            break

        # ----------------------------------------------------------
        # Print predicted.
        # ----------------------------------------------------------

        print(
            "  PREDICTED:"
        )

        print_state(
            predicted,
            p,
            q,
            prefix="    ",
        )

        # ----------------------------------------------------------
        # Print actual.
        # ----------------------------------------------------------

        if actual is None:

            print(
                "  ACTUAL: INVALID"
            )

            print(
                "  >>> RECURSION OUTLIVES "
                "THE VALID PARAMETERIZATION"
            )

            break

        print(
            "  ACTUAL:"
        )

        print_state(
            actual,
            p,
            q,
            prefix="    ",
        )

        # ----------------------------------------------------------
        # Compare.
        # ----------------------------------------------------------

        same = (
            predicted.a == actual.a
            and predicted.x == actual.x
            and predicted.y == actual.y
        )

        print(
            f"  MATCH = {same}"
        )

        if not same:

            print(
                "  >>> RECURSIVE MODEL FAILED"
            )

            break

        recursive = predicted

    # --------------------------------------------------------------
    # Final summary.
    # --------------------------------------------------------------

    print()
    print("FINAL SUMMARY")
    print("-" * 100)

    valid_actual = []

    for z in range(
        2,
        max_z + 1,
    ):

        actual = true_state(
            p,
            q,
            z,
        )

        if actual is None:
            break

        valid_actual.append(actual)

    print(
        f"Number of valid levels = "
        f"{len(valid_actual)}"
    )

    if valid_actual:

        print(
            f"First z = "
            f"{valid_actual[0].z}"
        )

        print(
            f"Last z = "
            f"{valid_actual[-1].z}"
        )


# ======================================================================
# BRANCH EXPLORER
# ======================================================================

def explore_branches(
    p: int,
    q: int,
    max_z: int = 16,
):

    print()
    print("=" * 100)
    print(
        f"BRANCH EXPLORATION: {p} × {q}"
    )
    print("=" * 100)

    start = true_state(
        p,
        q,
        2,
    )

    if start is None:

        print(
            "No z=2 state."
        )

        return

    print()
    print(
        "Starting state:"
    )

    print_state(
        start,
        p,
        q,
    )

    # --------------------------------------------------------------
    # At each level show whether each possible b would be legal.
    # --------------------------------------------------------------

    current = start

    for z in range(
        current.z,
        max_z,
    ):

        print()
        print(
            f"z={current.z}"
        )

        print(
            f"  current x = {current.x}"
        )

        print(
            f"  current y = {current.y}"
        )

        print()

        for b in [0, 1]:

            x_num = current.x + b
            y = current.y

            x_valid = (
                x_num % 2 == 0
            )

            y_valid = (
                y % 2 == 0
            )

            if x_valid and y_valid:

                next_x = x_num // 2
                next_y = y // 2

                next_a = (
                    current.a
                    + b * (1 << (current.z - 1))
                )

                print(
                    f"  b={b}: VALID "
                    f"-> "
                    f"a'={next_a} "
                    f"x'={next_x} "
                    f"y'={next_y}"
                )

            else:

                print(
                    f"  b={b}: INVALID "
                    f"(x parity={x_num % 2}, "
                    f"y parity={y % 2})"
                )

        predicted = recursive_transition(
            current
        )

        if predicted is None:
            break

        current = predicted


# ======================================================================
# INDEPENDENT RECURRENCE TEST
# ======================================================================

def test_parity_rule(
    p: int,
    q: int,
    max_z: int = 32,
):

    print()
    print("=" * 100)
    print(
        f"PARITY RULE TEST: {p} × {q}"
    )
    print("=" * 100)

    failures = 0
    transitions = 0

    for z in range(
        2,
        max_z,
    ):

        current = true_state(
            p,
            q,
            z,
        )

        nxt = true_state(
            p,
            q,
            z + 1,
        )

        if current is None or nxt is None:
            continue

        transitions += 1

        # ----------------------------------------------------------
        # Actual a difference.
        # ----------------------------------------------------------

        delta_a = (
            nxt.a - current.a
        )

        bit = (
            delta_a
            // (1 << (z - 1))
        )

        # ----------------------------------------------------------
        # Proposed bit from x parity.
        # ----------------------------------------------------------

        predicted_bit = (
            current.x & 1
        )

        ok_bit = (
            bit == predicted_bit
        )

        # ----------------------------------------------------------
        # x recurrence.
        # ----------------------------------------------------------

        predicted_x = (
            current.x
            + predicted_bit
        ) // 2

        ok_x = (
            predicted_x == nxt.x
        )

        # ----------------------------------------------------------
        # y recurrence.
        # ----------------------------------------------------------

        ok_y = (
            current.y % 2 == 0
            and current.y // 2 == nxt.y
        )

        ok = (
            ok_bit
            and ok_x
            and ok_y
        )

        print(
            f"z={z:2d}->{z+1:2d} "
            f"x={current.x:<8d} "
            f"y={current.y:<8d} "
            f"b={bit} "
            f"pred_b={predicted_bit} "
            f"x'={nxt.x:<8d} "
            f"y'={nxt.y:<8d} "
            f"{'PASS' if ok else 'FAIL'}"
        )

        if not ok:
            failures += 1

    print()
    print(
        f"Transitions = {transitions}"
    )

    print(
        f"Failures    = {failures}"
    )

    print(
        f"PASS        = {failures == 0}"
    )


# ======================================================================
# MAIN
# ======================================================================

def main():

    print("=" * 100)
    print(
        "2^z RECURSIVE STATE-MACHINE EXPERIMENT"
    )
    print("=" * 100)

    for p, q in CASES:

        compare_case(
            p,
            q,
            max_z=32,
        )

        test_parity_rule(
            p,
            q,
            max_z=32,
        )

        explore_branches(
            p,
            q,
            max_z=16,
        )


if __name__ == "__main__":
    main()
