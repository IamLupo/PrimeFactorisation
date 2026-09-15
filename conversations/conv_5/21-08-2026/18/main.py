#!/usr/bin/env python3

"""
======================================================================
2^z CLOSED-FORM / v2 THEOREM EXPERIMENT
======================================================================

For odd p,q define

    N = p*q

At z=2:

    a_2 = 1

    y_2 = (p+q)/4

    x_2 = (q-p+2)/4

when these are integers.

Observed recursive system:

    y_(z+1) = y_z / 2

    x_(z+1) = ceil(x_z / 2)

    a_(z+1) = a_z + (x_z mod 2)*2^(z-1)

Observed termination:

    last valid z = v2(p+q)

Equivalent:

    number of transitions = v2(y_2)

This experiment tests:

    H1:
        z=2 exists iff p+q ≡ 0 (mod 4)

    H2:
        last valid z = v2(p+q)

    H3:
        y_z = y_2 / 2^(z-2)

    H4:
        x_z = ceil(x_2 / 2^(z-2))

    H5:
        recursive state == direct state

    H6:
        a_z recursion remains exact

Unlike the previous experiment, this one tests ALL odd pairs,
not only primes.

======================================================================
"""

from __future__ import annotations

from dataclasses import dataclass


# ======================================================================
# CONFIGURATION
# ======================================================================

MAX_VALUE = 1000

MAX_Z = 32

TEST_COPRIME_ONLY = False


# ======================================================================
# DATA
# ======================================================================

@dataclass
class State:
    z: int
    a: int
    x: int
    y: int


@dataclass
class CaseResult:
    p: int
    q: int
    z2_exists: bool
    v2_sum: int | None
    predicted_last_z: int | None
    actual_last_z: int | None

    x_closed_form: bool
    y_closed_form: bool
    recursive_match: bool
    a_recurrence: bool

    passed: bool


# ======================================================================
# BASIC FUNCTIONS
# ======================================================================

def v2(n: int) -> int:

    if n == 0:
        return 10**9

    n = abs(n)

    result = 0

    while n % 2 == 0:
        n //= 2
        result += 1

    return result


def ceil_div_pow2(
    x: int,
    r: int,
) -> int:

    if r == 0:
        return x

    denominator = 1 << r

    return (
        x + denominator - 1
    ) // denominator


# ======================================================================
# DIRECT STATE
# ======================================================================

def direct_state(
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

    # Verify.

    p_check = (
        k * (y - x)
        + a
    )

    q_check = (
        k * (y + x)
        - a
    )

    if p_check != p:
        raise AssertionError()

    if q_check != q:
        raise AssertionError()

    return State(
        z=z,
        a=a,
        x=x,
        y=y,
    )


# ======================================================================
# RECURSIVE STATE
# ======================================================================

def recursive_next(
    state: State,
) -> State | None:

    # Next bit of p is forced by x parity.

    b = state.x & 1

    # y must be even.

    if state.y % 2 != 0:
        return None

    next_x = (
        state.x + b
    ) // 2

    next_y = (
        state.y // 2
    )

    next_a = (
        state.a
        + b * (1 << (state.z - 1))
    )

    return State(
        z=state.z + 1,
        a=next_a,
        x=next_x,
        y=next_y,
    )


# ======================================================================
# TRUE LEVELS
# ======================================================================

def true_levels(
    p: int,
    q: int,
) -> list[State]:

    levels = []

    for z in range(
        2,
        MAX_Z + 1,
    ):

        state = direct_state(
            p,
            q,
            z,
        )

        if state is None:
            break

        levels.append(state)

    return levels


# ======================================================================
# RECURSIVE LEVELS
# ======================================================================

def recursive_levels(
    initial: State,
) -> list[State]:

    levels = [
        initial
    ]

    current = initial

    while current.z < MAX_Z:

        nxt = recursive_next(
            current
        )

        if nxt is None:
            break

        levels.append(nxt)

        current = nxt

    return levels


# ======================================================================
# CASE TEST
# ======================================================================

def test_case(
    p: int,
    q: int,
) -> CaseResult:

    sum_pq = p + q

    # --------------------------------------------------------------
    # z=2 existence
    # --------------------------------------------------------------

    z2 = direct_state(
        p,
        q,
        2,
    )

    z2_exists = (
        z2 is not None
    )

    # --------------------------------------------------------------
    # Cases that do not enter the system.
    # --------------------------------------------------------------

    if not z2_exists:

        return CaseResult(
            p=p,
            q=q,
            z2_exists=False,
            v2_sum=None,
            predicted_last_z=None,
            actual_last_z=None,
            x_closed_form=True,
            y_closed_form=True,
            recursive_match=True,
            a_recurrence=True,
            passed=True,
        )

    # --------------------------------------------------------------
    # v2 prediction
    # --------------------------------------------------------------

    valuation = v2(
        sum_pq
    )

    predicted_last_z = valuation

    # --------------------------------------------------------------
    # True levels
    # --------------------------------------------------------------

    levels = true_levels(
        p,
        q,
    )

    actual_last_z = (
        levels[-1].z
        if levels
        else None
    )

    # --------------------------------------------------------------
    # H2: last z
    # --------------------------------------------------------------

    last_z_ok = (
        actual_last_z
        == predicted_last_z
    )

    # --------------------------------------------------------------
    # H3: y closed form
    # --------------------------------------------------------------

    y_closed = True

    y2 = levels[0].y

    for state in levels:

        r = state.z - 2

        numerator = y2

        denominator = 1 << r

        if numerator % denominator != 0:
            y_closed = False
            break

        predicted = (
            numerator // denominator
        )

        if predicted != state.y:
            y_closed = False
            break

    # --------------------------------------------------------------
    # H4: x closed form
    # --------------------------------------------------------------

    x_closed = True

    x2 = levels[0].x

    for state in levels:

        r = state.z - 2

        predicted = ceil_div_pow2(
            x2,
            r,
        )

        if predicted != state.x:
            x_closed = False
            break

    # --------------------------------------------------------------
    # H5: recursive match
    # --------------------------------------------------------------

    recursive = recursive_levels(
        levels[0]
    )

    recursive_match = (
        len(recursive)
        == len(levels)
        and all(
            a.z == b.z
            and a.a == b.a
            and a.x == b.x
            and a.y == b.y
            for a, b in zip(
                recursive,
                levels,
            )
        )
    )

    # --------------------------------------------------------------
    # H6: a recurrence
    # --------------------------------------------------------------

    a_recurrence = True

    for i in range(
        len(levels) - 1
    ):

        current = levels[i]
        nxt = levels[i + 1]

        b = current.x & 1

        expected = (
            current.a
            + b * (1 << (current.z - 1))
        )

        if nxt.a != expected:

            a_recurrence = False
            break

    # --------------------------------------------------------------
    # Final result
    # --------------------------------------------------------------

    passed = (
        last_z_ok
        and x_closed
        and y_closed
        and recursive_match
        and a_recurrence
    )

    return CaseResult(
        p=p,
        q=q,
        z2_exists=True,
        v2_sum=valuation,
        predicted_last_z=predicted_last_z,
        actual_last_z=actual_last_z,
        x_closed_form=x_closed,
        y_closed_form=y_closed,
        recursive_match=recursive_match,
        a_recurrence=a_recurrence,
        passed=passed,
    )


# ======================================================================
# EXHAUSTIVE ODD PAIR TEST
# ======================================================================

def exhaustive_test():

    print()
    print("=" * 100)
    print("EXHAUSTIVE ODD-PAIR TEST")
    print("=" * 100)

    total = 0
    entered = 0
    failures = []

    for p in range(
        3,
        MAX_VALUE + 1,
        2,
    ):

        for q in range(
            3,
            MAX_VALUE + 1,
            2,
        ):

            if TEST_COPRIME_ONLY:

                import math

                if math.gcd(p, q) != 1:
                    continue

            total += 1

            result = test_case(
                p,
                q,
            )

            if not result.z2_exists:
                continue

            entered += 1

            if not result.passed:

                failures.append(
                    result
                )

    print()
    print(
        f"odd ordered pairs tested = {total}"
    )

    print(
        f"pairs entering z=2       = {entered}"
    )

    print(
        f"failures                 = {len(failures)}"
    )

    print(
        f"PASS                     = "
        f"{len(failures) == 0}"
    )

    if failures:

        print()
        print("FIRST FAILURES")
        print("-" * 100)

        for result in failures[:20]:

            print(
                f"p={result.p} "
                f"q={result.q} "
                f"v2(p+q)={result.v2_sum} "
                f"pred={result.predicted_last_z} "
                f"actual={result.actual_last_z} "
                f"x={result.x_closed_form} "
                f"y={result.y_closed_form} "
                f"rec={result.recursive_match} "
                f"a={result.a_recurrence}"
            )


# ======================================================================
# DETAILED EXAMPLES
# ======================================================================

def detailed(
    p: int,
    q: int,
):

    result = test_case(
        p,
        q,
    )

    print()
    print("=" * 100)
    print(
        f"DETAILED: {p} × {q}"
    )
    print("=" * 100)

    print()
    print(
        f"p+q = {p+q}"
    )

    print(
        f"v2(p+q) = {v2(p+q)}"
    )

    print(
        f"predicted last z = "
        f"{result.predicted_last_z}"
    )

    print(
        f"actual last z = "
        f"{result.actual_last_z}"
    )

    if not result.z2_exists:

        print()
        print(
            "No z=2 parameterization."
        )

        return

    levels = true_levels(
        p,
        q,
    )

    print()
    print(
        f"{'z':>3} "
        f"{'a':>10} "
        f"{'x':>10} "
        f"{'y':>10} "
        f"{'closed x':>12} "
        f"{'closed y':>12}"
    )

    print("-" * 75)

    x2 = levels[0].x
    y2 = levels[0].y

    for state in levels:

        r = state.z - 2

        closed_x = ceil_div_pow2(
            x2,
            r,
        )

        closed_y = (
            y2 // (1 << r)
        )

        print(
            f"{state.z:3d} "
            f"{state.a:10d} "
            f"{state.x:10d} "
            f"{state.y:10d} "
            f"{closed_x:12d} "
            f"{closed_y:12d}"
        )

    print()
    print(
        f"x closed form = "
        f"{result.x_closed_form}"
    )

    print(
        f"y closed form = "
        f"{result.y_closed_form}"
    )

    print(
        f"last-z theorem = "
        f"{result.actual_last_z == result.predicted_last_z}"
    )

    print(
        f"recursive match = "
        f"{result.recursive_match}"
    )

    print(
        f"a recurrence = "
        f"{result.a_recurrence}"
    )

    print()
    print(
        f"CASE PASS = "
        f"{result.passed}"
    )


# ======================================================================
# MAIN
# ======================================================================

def main():

    print("=" * 100)
    print(
        "2^z CLOSED-FORM / v2 THEOREM EXPERIMENT"
    )
    print("=" * 100)

    # --------------------------------------------------------------
    # Known examples.
    # --------------------------------------------------------------

    detailed(
        59,
        101,
    )

    detailed(
        223,
        449,
    )

    detailed(
        101,
        103,
    )

    # --------------------------------------------------------------
    # Exhaustive odd-pair population.
    # --------------------------------------------------------------

    exhaustive_test()


if __name__ == "__main__":
    main()
