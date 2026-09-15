#!/usr/bin/env python3

"""
======================================================================
N -> 2^z RECURSIVE STATE SEARCH EXPERIMENT
======================================================================

GOAL

Start with N only.

Do NOT know p and q.

Use the discovered parameterization:

    k = 2^(z-1)

    p = k(y-x) + a
    q = k(y+x) - a

Therefore:

    N = p*q

         = k^2(y^2-x^2) + 2*a*k*x - a^2

and:

    N ≡ -a^2 (mod 2^z)

DISCOVERED RECURSION

    b = x mod 2

    a' = a + b*2^(z-1)

    x' = ceil(x/2)

    y' = y/2

A transition requires:

    y even

This experiment:

    1. Starts from N.
    2. Finds admissible a values.
    3. Searches possible (x,y) states.
    4. Checks N exactly.
    5. Recursively propagates each state.
    6. Reconstructs p,q.
    7. Checks primality.
    8. Compares against known factors when available.

The key question:

    Can the recursive state itself constrain the factorization
    strongly enough that we can recover p,q from N?

======================================================================
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isqrt
from sympy import isprime


# ======================================================================
# CONFIGURATION
# ======================================================================

N_VALUES = [
    5959,
    3233,
    10403,
    18721,
    100127,
    10025616383,
]

# Maximum z to search.
MAX_Z = 16

# Search window around candidate y.
#
# For an unknown factorization this is the expensive part.
# Keep small initially.
Y_SEARCH = 100_000

# Only display this many states per N.
MAX_RESULTS = 50

# Known factors, only used for comparison.
KNOWN = {
    5959: (59, 101),
    3233: (53, 61),
    10403: (101, 103),
    18721: (97, 193),
    100127: (223, 449),
    10025616383: (223, 449, 100129),
}


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
class Solution:
    z: int
    a: int
    x: int
    y: int
    p: int
    q: int
    p_prime: bool
    q_prime: bool


# ======================================================================
# MODULAR HELPERS
# ======================================================================

def admissible_a(
    n: int,
    z: int,
) -> list[int]:

    modulus = 1 << z
    a_modulus = 1 << (z - 1)

    target = n % modulus

    result = []

    for a in range(a_modulus):

        if (
            -a * a
        ) % modulus == target:

            result.append(a)

    return result


# ======================================================================
# EXACT N FROM STATE
# ======================================================================

def state_to_factors(
    state: State,
):
    """
    p = 2^(z-1)(y-x) + a
    q = 2^(z-1)(y+x) - a
    """

    k = 1 << (state.z - 1)

    p = (
        k * (state.y - state.x)
        + state.a
    )

    q = (
        k * (state.y + state.x)
        - state.a
    )

    return p, q


def state_to_n(
    state: State,
) -> int:

    p, q = state_to_factors(
        state
    )

    return p * q


# ======================================================================
# RECURSION
# ======================================================================

def next_state(
    state: State,
) -> State | None:

    # --------------------------------------------------------------
    # Transition only possible when y is even.
    # --------------------------------------------------------------

    if state.y % 2 != 0:
        return None

    # --------------------------------------------------------------
    # Next bit comes from x parity.
    # --------------------------------------------------------------

    b = state.x & 1

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

    return State(
        z=state.z + 1,
        a=next_a,
        x=next_x,
        y=next_y,
    )


# ======================================================================
# STATE SEARCH FROM N
# ======================================================================

def search_states(
    n: int,
    z: int,
):
    """
    Search integer x,y for the chosen z.

    We use:

        N = k^2(y^2-x^2) + 2akx - a^2

    Rearranged:

        N + a^2 - 2akx = k^2(y^2-x^2)

    For each x we solve for y^2.

    This is much cheaper than blindly scanning p,q.
    """

    k = 1 << (z - 1)

    candidates_a = admissible_a(
        n,
        z,
    )

    states = []

    for a in candidates_a:

        # ----------------------------------------------------------
        # We need:
        #
        # N = k^2(y^2-x^2) + 2akx - a^2
        #
        # Therefore:
        #
        # y^2 =
        # x^2
        # + (N + a^2 - 2akx)/k^2
        # ----------------------------------------------------------

        #
        # Crude x bound.
        #
        # The factors must be positive, so begin around zero.
        #

        # We don't want an enormous brute-force x range.
        #
        # Use the fact that:
        #
        # p = k(y-x)+a > 1
        # q = k(y+x)-a > 1
        #
        # so x is generally smaller than y.
        #

        max_x = Y_SEARCH

        for x in range(
            0,
            max_x + 1,
        ):

            numerator = (
                n
                + a * a
                - 2 * a * k * x
            )

            denominator = k * k

            if numerator <= 0:
                continue

            if numerator % denominator != 0:
                continue

            y2 = (
                x * x
                + numerator // denominator
            )

            if y2 <= 0:
                continue

            y = isqrt(y2)

            if y * y != y2:
                continue

            if y <= x:
                continue

            state = State(
                z=z,
                a=a,
                x=x,
                y=y,
            )

            p, q = state_to_factors(
                state
            )

            if p <= 1 or q <= 1:
                continue

            if p * q != n:
                continue

            states.append(
                state
            )

    return states


# ======================================================================
# VERIFY RECURSIVE STATE
# ======================================================================

def verify_state_chain(
    n: int,
    initial: State,
):
    """
    Starting from one state, recursively walk upward.
    """

    chain = []

    current = initial

    while True:

        p, q = state_to_factors(
            current
        )

        if p * q != n:
            break

        chain.append(
            (
                current,
                p,
                q,
            )
        )

        nxt = next_state(
            current
        )

        if nxt is None:
            break

        current = nxt

        if current.z > MAX_Z:
            break

    return chain


# ======================================================================
# PRINT CHAIN
# ======================================================================

def print_chain(
    n: int,
    chain,
):

    print()

    for state, p, q in chain:

        print(
            f"z={state.z:2d} "
            f"a={state.a:<8d} "
            f"x={state.x:<8d} "
            f"y={state.y:<8d} "
            f"p={p:<12d} "
            f"q={q:<12d} "
            f"prime=({isprime(p)},{isprime(q)})"
        )


# ======================================================================
# ONE N
# ======================================================================

def experiment_n(
    n: int,
):

    print()
    print("=" * 110)
    print(f"N = {n}")
    print("=" * 110)

    if n in KNOWN:

        print(
            "KNOWN FACTORIZATION:",
            KNOWN[n]
        )

    total_states = 0
    solutions = []

    # --------------------------------------------------------------
    # Search each z independently.
    # --------------------------------------------------------------

    for z in range(
        2,
        MAX_Z + 1,
    ):

        states = search_states(
            n,
            z,
        )

        if not states:
            continue

        print()
        print(
            f"z={z} "
            f"admissible states={len(states)}"
        )

        for state in states:

            total_states += 1

            p, q = state_to_factors(
                state
            )

            pp = isprime(p)
            qq = isprime(q)

            solution = Solution(
                z=z,
                a=state.a,
                x=state.x,
                y=state.y,
                p=p,
                q=q,
                p_prime=pp,
                q_prime=qq,
            )

            solutions.append(
                solution
            )

            print(
                f"  a={state.a:<8d} "
                f"x={state.x:<8d} "
                f"y={state.y:<8d} "
                f"p={p:<12d} "
                f"q={q:<12d} "
                f"prime=({pp},{qq})"
            )

            # Stop massive output.
            if len(solutions) >= MAX_RESULTS:
                break

        if len(solutions) >= MAX_RESULTS:
            break

    # --------------------------------------------------------------
    # Recursive chains.
    # --------------------------------------------------------------

    print()
    print("RECURSIVE CHAINS")
    print("-" * 110)

    for solution in solutions:

        initial = State(
            z=solution.z,
            a=solution.a,
            x=solution.x,
            y=solution.y,
        )

        chain = verify_state_chain(
            n,
            initial,
        )

        print()
        print(
            f"START z={initial.z} "
            f"a={initial.a} "
            f"x={initial.x} "
            f"y={initial.y}"
        )

        print_chain(
            n,
            chain,
        )

    # --------------------------------------------------------------
    # Summary.
    # --------------------------------------------------------------

    print()
    print("-" * 110)

    prime_pairs = [
        s
        for s in solutions
        if s.p_prime and s.q_prime
    ]

    print(
        f"states found     = {len(solutions)}"
    )

    print(
        f"prime x prime    = {len(prime_pairs)}"
    )


# ======================================================================
# SPECIAL SEARCH: BUILD FACTORS FROM STATE
# ======================================================================

def special_state_search(
    n: int,
):
    """
    Instead of searching x up to Y_SEARCH for every z,
    examine only states which can recursively descend from
    lower levels.

    This is closer to the SAT-state idea.
    """

    print()
    print("=" * 110)
    print("RECURSIVE STATE SEARCH")
    print("=" * 110)

    # --------------------------------------------------------------
    # Start at z=2.
    # --------------------------------------------------------------

    initial_states = search_states(
        n,
        2,
    )

    if not initial_states:

        print(
            "No z=2 states."
        )

        return

    print(
        f"z=2 initial states = "
        f"{len(initial_states)}"
    )

    # --------------------------------------------------------------
    # Propagate every state recursively.
    # --------------------------------------------------------------

    for initial in initial_states:

        print()
        print(
            "INITIAL:",
            initial
        )

        chain = verify_state_chain(
            n,
            initial,
        )

        print_chain(
            n,
            chain,
        )


# ======================================================================
# MAIN
# ======================================================================

def main():

    print("=" * 110)
    print(
        "N-ONLY RECURSIVE 2^z PARAMETER SEARCH"
    )
    print("=" * 110)

    for n in N_VALUES:

        experiment_n(
            n
        )

        special_state_search(
            n
        )


if __name__ == "__main__":
    main()
