#!/usr/bin/env python3

"""
======================================================================
DIFFERENCE-OF-SQUARES vs 2^z RECURSION -- CORRECTED
======================================================================

At z=2:

    p = 2y - 2x + 1
    q = 2y + 2x - 1

Therefore:

    N = p*q

      = (2y)^2 - (2x-1)^2

Let:

    A = 2y
    B = 2x - 1

Then:

    N = A^2 - B^2

and

    p = A - B
    q = A + B

The recursive refinement is:

    b = x mod 2

    a' = a + b*2^(z-1)

    x' = ceil(x/2)

    y' = y/2

This experiment compares the direct difference-of-squares
representation against the recursive 2^z representation.
======================================================================
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isqrt

import sympy as sp


# ======================================================================
# CONFIGURATION
# ======================================================================

N_VALUES = [
    5959,
    3233,
    10403,
    18721,
    100127,
    100127 * 100129,
]

MAX_Z = 20

# Maximum number of A increments in the direct search.
MAX_FERMAT_STEPS = 1_000_000

MAX_RESULTS = 50


# ======================================================================
# DATA
# ======================================================================

@dataclass
class State:
    z: int
    a: int
    x: int
    y: int
    p: int
    q: int


# ======================================================================
# z=2 STATE
# ======================================================================

def state_from_xy(
    x: int,
    y: int,
) -> State:

    p = 2 * y - 2 * x + 1
    q = 2 * y + 2 * x - 1

    return State(
        z=2,
        a=1,
        x=x,
        y=y,
        p=p,
        q=q,
    )


# ======================================================================
# DIRECT DIFFERENCE-OF-SQUARES SEARCH
# ======================================================================

def search_z2_states(
    n: int,
) -> list[State]:

    """
    Solve:

        N = A^2 - B^2

    where:

        A = 2y
        B = 2x-1

    Therefore:

        A must be even
        B must be odd

    and:

        p = A-B
        q = A+B
    """

    target = n

    # Start at ceil(sqrt(N)).
    A = isqrt(target)

    if A * A < target:
        A += 1

    states: list[State] = []

    for _ in range(MAX_FERMAT_STEPS):

        B2 = A * A - target

        if B2 >= 0:

            B = isqrt(B2)

            if B * B == B2:

                # Our parameterization requires:
                #
                # A = 2y
                # B = 2x-1

                if A % 2 == 0 and B % 2 == 1:

                    y = A // 2
                    x = (B + 1) // 2

                    state = state_from_xy(
                        x,
                        y,
                    )

                    if (
                        state.p > 1
                        and state.q > 1
                        and state.p * state.q == n
                    ):
                        states.append(state)

        # Once A is larger than N, the remaining
        # factors are trivial.
        if A > n + 2:
            break

        A += 1

    return states


# ======================================================================
# DIRECT GENERAL z STATE
# ======================================================================

def direct_state(
    p: int,
    q: int,
    z: int,
) -> State | None:

    modulus = 1 << z
    k = 1 << (z - 1)

    # Low z-1 bits of p.
    a = p % k

    y_num = p + q
    x_num = q - p + 2 * a

    if y_num % modulus != 0:
        return None

    if x_num % modulus != 0:
        return None

    x = x_num // modulus
    y = y_num // modulus

    p_check = (
        k * (y - x)
        + a
    )

    q_check = (
        k * (y + x)
        - a
    )

    if p_check != p or q_check != q:
        return None

    return State(
        z=z,
        a=a,
        x=x,
        y=y,
        p=p,
        q=q,
    )


# ======================================================================
# RECURSIVE TRANSITION
# ======================================================================

def next_state(
    state: State,
) -> State | None:

    # y must be even.
    if state.y % 2 != 0:
        return None

    # Next factor bit.
    b = state.x & 1

    next_z = state.z + 1

    next_a = (
        state.a
        + b * (1 << (state.z - 1))
    )

    next_x = (
        state.x + b
    ) // 2

    next_y = (
        state.y // 2
    )

    next_k = 1 << (next_z - 1)

    next_p = (
        next_k * (next_y - next_x)
        + next_a
    )

    next_q = (
        next_k * (next_y + next_x)
        - next_a
    )

    return State(
        z=next_z,
        a=next_a,
        x=next_x,
        y=next_y,
        p=next_p,
        q=next_q,
    )


# ======================================================================
# BUILD RECURSIVE CHAIN
# ======================================================================

def build_chain(
    initial: State,
) -> list[State]:

    chain = [initial]

    current = initial

    while current.z < MAX_Z:

        nxt = next_state(current)

        if nxt is None:
            break

        if nxt.p * nxt.q != initial.p * initial.q:
            raise AssertionError(
                "Recursive transition changed N"
            )

        chain.append(nxt)
        current = nxt

    return chain


# ======================================================================
# PRINT CHAIN
# ======================================================================

def print_chain(
    chain: list[State],
) -> None:

    for state in chain:

        pp = bool(sp.isprime(state.p))
        qq = bool(sp.isprime(state.q))

        print(
            f"  z={state.z:2d} "
            f"a={state.a:<8d} "
            f"x={state.x:<8d} "
            f"y={state.y:<8d} "
            f"p={state.p:<12d} "
            f"q={state.q:<12d} "
            f"prime=({pp},{qq})"
        )


# ======================================================================
# DIFFERENCE-OF-SQUARES VERIFICATION
# ======================================================================

def verify_difference_of_squares(
    n: int,
    state: State,
) -> bool:

    A = 2 * state.y
    B = 2 * state.x - 1

    lhs = (
        A * A
        - B * B
    )

    return lhs == n


# ======================================================================
# CLOSED FORM CHECK
# ======================================================================

def check_closed_form(
    initial: State,
    chain: list[State],
) -> None:

    x2 = initial.x
    y2 = initial.y

    print()
    print("CLOSED-FORM CHECK")
    print("-" * 90)

    for state in chain:

        r = state.z - 2
        divisor = 1 << r

        predicted_x = (
            x2 + divisor - 1
        ) // divisor

        if y2 % divisor == 0:
            predicted_y = (
                y2 // divisor
            )
        else:
            predicted_y = None

        ok_x = (
            state.x == predicted_x
        )

        ok_y = (
            predicted_y is not None
            and state.y == predicted_y
        )

        ok = ok_x and ok_y

        print(
            f"  z={state.z:2d} "
            f"x={state.x:<8d} "
            f"x_closed={predicted_x:<8d} "
            f"y={state.y:<8d} "
            f"y_closed={str(predicted_y):<8s} "
            f"PASS={ok}"
        )


# ======================================================================
# METHOD EQUIVALENCE
# ======================================================================

def show_difference_of_squares(
    n: int,
    state: State,
) -> None:

    A = 2 * state.y
    B = 2 * state.x - 1

    p_from_dos = A - B
    q_from_dos = A + B

    print()
    print(
        f"  A = 2y   = {A}"
    )

    print(
        f"  B = 2x-1 = {B}"
    )

    print(
        f"  A²-B²   = {A*A - B*B}"
    )

    print(
        f"  N       = {n}"
    )

    print(
        f"  p=A-B   = {p_from_dos}"
    )

    print(
        f"  q=A+B   = {q_from_dos}"
    )

    print(
        f"  MATCH    = "
        f"{A * A - B * B == n}"
    )


# ======================================================================
# COMPARE ONE N
# ======================================================================

def compare_n(
    n: int,
) -> None:

    print()
    print("=" * 110)
    print(
        f"N = {n}"
    )
    print("=" * 110)

    # --------------------------------------------------------------
    # Reference factorization
    # --------------------------------------------------------------

    try:

        factors = sp.factorint(n)

        factor_text = " × ".join(
            (
                f"{prime}^{exponent}"
                if exponent != 1
                else str(prime)
            )
            for prime, exponent in factors.items()
        )

        print(
            f"REFERENCE FACTORIZATION: "
            f"{factor_text}"
        )

    except Exception as exc:

        print(
            f"REFERENCE FACTORIZATION FAILED: "
            f"{exc}"
        )

    # --------------------------------------------------------------
    # METHOD A
    # --------------------------------------------------------------

    print()
    print(
        "METHOD A: DIFFERENCE OF SQUARES"
    )
    print("-" * 110)

    states = search_z2_states(n)

    if not states:

        print(
            "No z=2 difference-of-squares state found."
        )

        return

    for state in states[:MAX_RESULTS]:

        pp = bool(sp.isprime(state.p))
        qq = bool(sp.isprime(state.q))

        dos_ok = verify_difference_of_squares(
            n,
            state,
        )

        print(
            f"  x={state.x:<10d} "
            f"y={state.y:<10d} "
            f"p={state.p:<12d} "
            f"q={state.q:<12d} "
            f"prime=({pp},{qq}) "
            f"DOS={dos_ok}"
        )

    # --------------------------------------------------------------
    # METHOD B
    # --------------------------------------------------------------

    print()
    print(
        "METHOD B: RECURSIVE 2^z REFINEMENT"
    )
    print("-" * 110)

    for state in states[:MAX_RESULTS]:

        chain = build_chain(
            state
        )

        print()
        print(
            f"START: "
            f"x={state.x}, "
            f"y={state.y}, "
            f"p={state.p}, "
            f"q={state.q}"
        )

        print_chain(
            chain
        )

        check_closed_form(
            state,
            chain,
        )

    # --------------------------------------------------------------
    # Difference-of-squares details
    # --------------------------------------------------------------

    print()
    print(
        "DIFFERENCE-OF-SQUARES DETAILS"
    )
    print("-" * 110)

    for state in states[:MAX_RESULTS]:

        print()
        print(
            f"Factor pair: "
            f"{state.p} × {state.q}"
        )

        show_difference_of_squares(
            n,
            state,
        )

    # --------------------------------------------------------------
    # Direct z comparison
    # --------------------------------------------------------------

    print()
    print(
        "DIRECT z-LEVEL CHECK"
    )
    print("-" * 110)

    for state in states[:MAX_RESULTS]:

        print()
        print(
            f"Initial factor pair: "
            f"{state.p} × {state.q}"
        )

        for z in range(
            2,
            MAX_Z + 1,
        ):

            direct = direct_state(
                state.p,
                state.q,
                z,
            )

            if direct is None:
                break

            print(
                f"  z={z:2d} "
                f"a={direct.a:<8d} "
                f"x={direct.x:<8d} "
                f"y={direct.y:<8d}"
            )


# ======================================================================
# MAIN
# ======================================================================

def main():

    print("=" * 110)
    print(
        "CORRECTED DIFFERENCE-OF-SQUARES "
        "vs 2^z RECURSION"
    )
    print("=" * 110)

    for n in N_VALUES:

        compare_n(n)


if __name__ == "__main__":
    main()