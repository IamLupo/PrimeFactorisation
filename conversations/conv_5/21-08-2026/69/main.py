#!/usr/bin/env python3

from __future__ import annotations

from dataclasses import dataclass
from collections import Counter
from math import gcd, isqrt


# ==============================================================================
# EXPERIMENT 631
# ==============================================================================
#
# DEGENERATE GCD BRANCH + EXACT INTEGER GCD TRANSFORMATION
#
# Experiment 630 found:
#
#     depth = 1 + v2(gcd(X,Y))
#
# exactly for every state.
#
# The attempted stronger identity
#
#     gcd(X,Y) = gcd(A,B)
#                 or gcd(A,B)/2
#
# failed only for:
#
#     n=9
#     p=q=3
#     frame=B
#     A=4
#     B=0
#
# because:
#
#     X=(3A-B)/2
#     Y=(3A+B)/2
#
# gives:
#
#     X=Y=6
#
# and therefore:
#
#     gcd(X,Y)=6
#
# whereas:
#
#     gcd(A,B)=4.
#
# This experiment isolates all degenerate branches and searches
# for the exact integer gcd law.
#
# ==============================================================================


PRIME_LIMIT = 6250
INF = 10**9


# ==============================================================================
# STATE
# ==============================================================================

@dataclass(frozen=True)
class State:
    n: int
    p: int
    q: int


# ==============================================================================
# v2
# ==============================================================================

def v2(x: int) -> int:

    if x == 0:
        return INF

    x = abs(x)

    return (x & -x).bit_length() - 1


# ==============================================================================
# PRIME SIEVE
# ==============================================================================

def prime_sieve(
    limit: int,
) -> list[int]:

    sieve = bytearray(
        b"\x01"
    ) * (limit + 1)

    sieve[0] = 0
    sieve[1] = 0

    for p in range(
        2,
        isqrt(limit) + 1,
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

    return [
        p
        for p in range(
            3,
            limit + 1,
            2,
        )
        if sieve[p]
    ]


# ==============================================================================
# SEMIPRIME GENERATION
# ==============================================================================

def generate_states(
    primes: list[int],
) -> list[State]:

    states: list[State] = []

    for i, p in enumerate(primes):

        for q in primes[i:]:

            states.append(
                State(
                    n=p * q,
                    p=p,
                    q=q,
                )
            )

    return states


# ==============================================================================
# FRAME
# ==============================================================================

def frame_of(
    n: int,
) -> str:

    r = n & 3

    if r == 3:
        return "A"

    if r == 1:
        return "B"

    raise ValueError(
        f"unexpected odd n mod 4={r}"
    )


# ==============================================================================
# RAW RESIDUALS
# ==============================================================================

def raw_residuals(
    frame: str,
    p: int,
    q: int,
) -> tuple[int, int]:

    if frame == "A":

        return (
            p - 3,
            q + 3,
        )

    return (
        p + 1,
        q - 3,
    )


# ==============================================================================
# GLOBAL X,Y
# ==============================================================================

def global_xy(
    frame: str,
    p: int,
    q: int,
) -> tuple[int, int]:

    if frame == "A":

        X = (
            q - p + 6
        ) // 2

        Y = (
            p + q
        ) // 2

    else:

        X = (
            3 * p - q + 6
        ) // 2

        Y = (
            3 * p + q
        ) // 2

    return X, Y


# ==============================================================================
# DEPTH
# ==============================================================================

def depth_from_xy(
    X: int,
    Y: int,
) -> int:

    return (
        1
        + min(
            v2(X),
            v2(Y),
        )
    )


# ==============================================================================
# SIGN / ZERO CLASSIFICATION
# ==============================================================================

def residual_class(
    A: int,
    B: int,
) -> str:

    if A == 0 and B == 0:
        return "A=0,B=0"

    if A == 0:
        return "A=0"

    if B == 0:
        return "B=0"

    if A == B:
        return "A=B"

    if A == -B:
        return "A=-B"

    return "generic"


# ==============================================================================
# TEST 0
# ==============================================================================

def test_zero_branches(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 0: DEGENERATE RESIDUAL BRANCHES")
    print("=" * 90)

    counts = Counter()

    examples: dict[str, list[State]] = {}

    for state in states:

        frame = frame_of(
            state.n
        )

        A, B = raw_residuals(
            frame,
            state.p,
            state.q,
        )

        cls = residual_class(
            A,
            B,
        )

        counts[cls] += 1

        examples.setdefault(
            cls,
            [],
        )

        if len(
            examples[cls]
        ) < 10:

            examples[cls].append(
                state
            )

    for cls in [
        "A=0,B=0",
        "A=0",
        "B=0",
        "A=B",
        "A=-B",
        "generic",
    ]:

        print(
            f"    {cls:10s}: "
            f"{counts[cls]}"
        )

    print()

    for cls in [
        "A=0,B=0",
        "A=0",
        "B=0",
        "A=B",
        "A=-B",
    ]:

        if not examples.get(cls):
            continue

        print(
            f"    examples [{cls}]"
        )

        for state in examples[cls]:

            frame = frame_of(
                state.n
            )

            A, B = raw_residuals(
                frame,
                state.p,
                state.q,
            )

            X, Y = global_xy(
                frame,
                state.p,
                state.q,
            )

            print(
                f"        n={state.n}"
                f" p={state.p}"
                f" q={state.q}"
                f" frame={frame}"
                f" A={A}"
                f" B={B}"
                f" X={X}"
                f" Y={Y}"
            )


# ==============================================================================
# TEST 1
# ==============================================================================

def test_generic_gcd_formula(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 1: GENERIC GCD FORMULA")
    print("=" * 90)

    failures = 0
    checked = 0

    ratio_counts = Counter()

    for state in states:

        frame = frame_of(
            state.n
        )

        A, B = raw_residuals(
            frame,
            state.p,
            state.q,
        )

        if A == 0 or B == 0:
            continue

        checked += 1

        X, Y = global_xy(
            frame,
            state.p,
            state.q,
        )

        d = gcd(
            abs(A),
            abs(B),
        )

        g = gcd(
            abs(X),
            abs(Y),
        )

        ratio = (
            d // g
            if g and d % g == 0
            else None
        )

        ratio_counts[ratio] += 1

        alpha = v2(A)
        beta = v2(B)

        if alpha == beta:

            expected = d

        else:

            expected = d // 2

        if g != expected:

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch"
                    f" n={state.n}"
                    f" frame={frame}"
                    f" A={A}"
                    f" B={B}"
                    f" gcdAB={d}"
                    f" gcdXY={g}"
                    f" v2A={alpha}"
                    f" v2B={beta}"
                    f" expected={expected}"
                )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )

    print(
        "    ratios:"
    )

    for key in sorted(
        ratio_counts,
        key=lambda x: (
            x is None,
            x if x is not None else 0,
        ),
    ):

        print(
            f"        {key}: "
            f"{ratio_counts[key]}"
        )


# ==============================================================================
# TEST 2
# ==============================================================================

def test_zero_A_branch(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 2: A=0 BRANCH")
    print("=" * 90)

    failures = 0
    checked = 0

    for state in states:

        frame = frame_of(
            state.n
        )

        A, B = raw_residuals(
            frame,
            state.p,
            state.q,
        )

        if A != 0:
            continue

        checked += 1

        X, Y = global_xy(
            frame,
            state.p,
            state.q,
        )

        d = abs(B)

        g = gcd(
            abs(X),
            abs(Y),
        )

        if frame == "A":

            # X=B/2, Y=B/2
            expected = abs(B) // 2

        else:

            # X=-B/2, Y=B/2
            expected = abs(B) // 2

        if g != expected:

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch"
                    f" n={state.n}"
                    f" frame={frame}"
                    f" B={B}"
                    f" gcdXY={g}"
                    f" expected={expected}"
                )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 3
# ==============================================================================

def test_zero_B_branch(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 3: B=0 BRANCH")
    print("=" * 90)

    failures = 0
    checked = 0

    for state in states:

        frame = frame_of(
            state.n
        )

        A, B = raw_residuals(
            frame,
            state.p,
            state.q,
        )

        if B != 0:
            continue

        checked += 1

        X, Y = global_xy(
            frame,
            state.p,
            state.q,
        )

        d = abs(A)

        g = gcd(
            abs(X),
            abs(Y),
        )

        if frame == "A":

            # X=-A/2, Y=A/2
            expected = abs(A) // 2

        else:

            # X=3A/2, Y=3A/2
            expected = (
                3 * abs(A)
            ) // 2

        if g != expected:

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch"
                    f" n={state.n}"
                    f" frame={frame}"
                    f" A={A}"
                    f" gcdXY={g}"
                    f" expected={expected}"
                )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 4
# ==============================================================================

def test_exact_piecewise_gcd(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 4: EXACT PIECEWISE gcd(X,Y) LAW")
    print("=" * 90)

    failures = 0

    for state in states:

        frame = frame_of(
            state.n
        )

        A, B = raw_residuals(
            frame,
            state.p,
            state.q,
        )

        X, Y = global_xy(
            frame,
            state.p,
            state.q,
        )

        g = gcd(
            abs(X),
            abs(Y),
        )

        alpha = v2(A)
        beta = v2(B)

        # Explicit degenerate zero handling.

        if A == 0:

            if B == 0:

                # This does not occur in the tested
                # odd semiprime domain, but keep exact.

                expected = 0

            elif frame == "A":

                expected = abs(B) // 2

            else:

                expected = abs(B) // 2

        elif B == 0:

            if frame == "A":

                expected = abs(A) // 2

            else:

                expected = (
                    3 * abs(A)
                ) // 2

        else:

            d = gcd(
                abs(A),
                abs(B),
            )

            if alpha == beta:

                expected = d

            else:

                expected = d // 2

        if g != expected:

            failures += 1

            if failures <= 30:

                print(
                    "    mismatch"
                    f" n={state.n}"
                    f" frame={frame}"
                    f" A={A}"
                    f" B={B}"
                    f" X={X}"
                    f" Y={Y}"
                    f" v2A={alpha}"
                    f" v2B={beta}"
                    f" gcdXY={g}"
                    f" expected={expected}"
                )

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 5
# ==============================================================================

def test_odd_part_generic(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 5: ODD-PART GCD RELATION")
    print("=" * 90)

    failures = 0

    for state in states:

        frame = frame_of(
            state.n
        )

        A, B = raw_residuals(
            frame,
            state.p,
            state.q,
        )

        if A == 0 or B == 0:
            continue

        X, Y = global_xy(
            frame,
            state.p,
            state.q,
        )

        d = gcd(
            abs(A),
            abs(B),
        )

        g = gcd(
            abs(X),
            abs(Y),
        )

        vd = v2(d)
        vg = v2(g)

        odd_d = (
            d >> vd
        )

        odd_g = (
            g >> vg
        )

        if odd_d != odd_g:

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch"
                    f" n={state.n}"
                    f" frame={frame}"
                    f" A={A}"
                    f" B={B}"
                    f" gcdAB={d}"
                    f" gcdXY={g}"
                    f" oddAB={odd_d}"
                    f" oddXY={odd_g}"
                )

    print(
        f"generic checked "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 6
# ==============================================================================

def test_frame_B_prime_factor_three(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 6: FRAME-B FACTOR-3 EFFECT")
    print("=" * 90)

    print(
        """
FRAME B:

    X = (3A-B)/2
    Y = (3A+B)/2

Hence:

    gcd(X,Y)
        is controlled by
    gcd(3A-B, 3A+B).

And:

    gcd(3A-B, 3A+B)
        divides
    6*gcd(A,B).

The only possible extra odd factor introduced
by the transformation is therefore 3.

This test checks whether that factor 3 occurs
only in the B=0 branch.
"""
    )

    failures = 0
    checked = 0

    extra_three_states = []

    for state in states:

        frame = frame_of(
            state.n
        )

        if frame != "B":
            continue

        A, B = raw_residuals(
            frame,
            state.p,
            state.q,
        )

        if A == 0 or B == 0:
            continue

        checked += 1

        X, Y = global_xy(
            frame,
            state.p,
            state.q,
        )

        d = gcd(
            abs(A),
            abs(B),
        )

        g = gcd(
            abs(X),
            abs(Y),
        )

        # Normalize both by their powers of two.

        vd = v2(d)
        vg = v2(g)

        odd_d = d >> vd
        odd_g = g >> vg

        if odd_g % 3 == 0:

            extra_three_states.append(
                state
            )

    if extra_three_states:

        failures = len(
            extra_three_states
        )

        for state in (
            extra_three_states[:20]
        ):

            frame = frame_of(
                state.n
            )

            A, B = raw_residuals(
                frame,
                state.p,
                state.q,
            )

            print(
                "    unexpected 3-factor"
                f" n={state.n}"
                f" p={state.p}"
                f" q={state.q}"
                f" A={A}"
                f" B={B}"
            )

    print(
        f"generic frame-B checked={checked}"
    )

    print(
        f"extra-three failures={failures}"
    )


# ==============================================================================
# TEST 7
# ==============================================================================

def test_depth_from_factor_gcd(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 7: DEPTH FROM FACTOR-SIDE INTEGER DATA")
    print("=" * 90)

    failures = 0

    for state in states:

        frame = frame_of(
            state.n
        )

        A, B = raw_residuals(
            frame,
            state.p,
            state.q,
        )

        X, Y = global_xy(
            frame,
            state.p,
            state.q,
        )

        d = gcd(
            abs(A),
            abs(B),
        )

        alpha = v2(A)
        beta = v2(B)

        # Exact depth law.

        if alpha == beta:

            predicted = (
                v2(d) + 1
            )

        else:

            predicted = (
                v2(d)
            )

        actual = depth_from_xy(
            X,
            Y,
        )

        if predicted != actual:

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch"
                    f" n={state.n}"
                    f" frame={frame}"
                    f" A={A}"
                    f" B={B}"
                    f" predicted={predicted}"
                    f" actual={actual}"
                )

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 8
# ==============================================================================

def test_degenerate_depth_branches(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 8: DEGENERATE DEPTH BRANCHES")
    print("=" * 90)

    failures = 0
    counts = Counter()

    for state in states:

        frame = frame_of(
            state.n
        )

        A, B = raw_residuals(
            frame,
            state.p,
            state.q,
        )

        X, Y = global_xy(
            frame,
            state.p,
            state.q,
        )

        cls = residual_class(
            A,
            B,
        )

        counts[cls] += 1

        actual = depth_from_xy(
            X,
            Y,
        )

        if cls == "B=0":

            # Only possible exceptional branch expected:
            # frame B, p=q=3.

            expected = (
                1
                + v2(
                    3 * abs(A) // 2
                )
            )

        elif cls == "A=0":

            expected = (
                1
                + v2(
                    abs(B) // 2
                )
            )

        elif cls == "A=B":

            expected = (
                1
                + v2(
                    abs(A)
                )
            )

        elif cls == "A=-B":

            expected = (
                1
                + v2(
                    abs(A)
                )
            )

        else:

            continue

        if expected != actual:

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch"
                    f" n={state.n}"
                    f" frame={frame}"
                    f" class={cls}"
                    f" A={A}"
                    f" B={B}"
                    f" expected={expected}"
                    f" actual={actual}"
                )

    print(
        "    branch counts:"
    )

    for key, value in sorted(
        counts.items()
    ):

        print(
            f"        {key}: {value}"
        )

    print(
        f"failures={failures}"
    )


# ==============================================================================
# EXAMPLES
# ==============================================================================

def print_examples(
    states: list[State],
) -> None:

    targets = [
        9,
        15,
        21,
        33,
        39,
        51,
        57,
        69,
        87,
        93,
        111,
        141,
        183,
        213,
    ]

    by_n = {
        state.n: state
        for state in states
    }

    print()
    print("=" * 90)
    print("EXAMPLES")
    print("=" * 90)

    for n in targets:

        state = by_n.get(n)

        if state is None:
            continue

        frame = frame_of(n)

        A, B = raw_residuals(
            frame,
            state.p,
            state.q,
        )

        X, Y = global_xy(
            frame,
            state.p,
            state.q,
        )

        d = gcd(
            abs(A),
            abs(B),
        )

        g = gcd(
            abs(X),
            abs(Y),
        )

        alpha = v2(A)
        beta = v2(B)

        print()
        print(
            f"n={n}"
            f" p={state.p}"
            f" q={state.q}"
        )

        print(
            f"    frame={frame}"
        )

        print(
            f"    A={A}"
            f" B={B}"
        )

        print(
            f"    X={X}"
            f" Y={Y}"
        )

        print(
            f"    gcd(A,B)={d}"
        )

        print(
            f"    gcd(X,Y)={g}"
        )

        print(
            f"    v2(A)={alpha}"
            f" v2(B)={beta}"
        )

        print(
            f"    class="
            f"{residual_class(A,B)}"
        )

        print(
            f"    depth="
            f"{depth_from_xy(X,Y)}"
        )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print("EXPERIMENT 631 START")
    print("=" * 90)

    print()
    print("[1] PRIME SIEVE")

    primes = prime_sieve(
        PRIME_LIMIT
    )

    print(
        f"    odd primes={len(primes)}"
    )

    print()
    print("[2] SEMIPRIME GENERATION")

    states = generate_states(
        primes
    )

    print(
        f"    semiprimes={len(states)}"
    )

    print()
    print("[3] GLOBAL STATES")

    print(
        f"    states={len(states)}"
    )

    test_zero_branches(
        states
    )

    test_generic_gcd_formula(
        states
    )

    test_zero_A_branch(
        states
    )

    test_zero_B_branch(
        states
    )

    test_exact_piecewise_gcd(
        states
    )

    test_odd_part_generic(
        states
    )

    test_frame_B_prime_factor_three(
        states
    )

    test_depth_from_factor_gcd(
        states
    )

    test_degenerate_depth_branches(
        states
    )

    print_examples(
        states
    )

    print()
    print("=" * 90)
    print("FINAL SUMMARY")
    print("=" * 90)

    print(
        """
The global law is already exact:

    depth = 1 + v2(gcd(X,Y)).

The remaining question is the exact INTEGER gcd bridge.

For nonzero residuals:

    alpha = v2(A)
    beta  = v2(B)

    alpha = beta
        ->
    gcd(X,Y) = gcd(A,B)

    alpha != beta
        ->
    gcd(X,Y) = gcd(A,B)/2.

The exceptional state from Experiment 630 is:

    frame B
    B = 0

where:

    X = 3A/2
    Y = 3A/2

and therefore:

    gcd(X,Y) = 3A/2.

The goal of Experiment 631 is to determine whether
this is the ONLY exceptional branch, and whether the
factor 3 in frame B can occur anywhere else.

If yes, the integer hierarchy becomes:

    GENERIC:

        gcd(X,Y)
          =
        gcd(A,B) / 2^[alpha != beta]

    SPECIAL B=0:

        gcd(X,Y)
          =
        3|A|/2.

Together with

    depth = 1 + v2(gcd(X,Y))

this gives a complete integer-level description of
the observed 2-adic hierarchy.
"""
    )

    print()
    print("=" * 90)
    print("EXPERIMENT 631 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
