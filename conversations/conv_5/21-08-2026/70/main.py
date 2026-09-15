#!/usr/bin/env python3

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import gcd, isqrt


# ==============================================================================
# EXPERIMENT 632
# ==============================================================================
#
# EXACT GCD TRANSFORMATION THEOREM
#
# Goal:
#
#   Turn the empirical Experiment 631 result into an exact algebraic
#   classification of gcd(X,Y) from the raw factor residuals A,B.
#
#
# RAW RESIDUALS
#
# FRAME A:
#
#     A = p - 3
#     B = q + 3
#
#     X = (B-A)/2
#     Y = (A+B)/2
#
#
# FRAME B:
#
#     A = p + 1
#     B = q - 3
#
#     X = (3A-B)/2
#     Y = (3A+B)/2
#
#
# EMPIRICAL LAW
#
# For A,B != 0:
#
#     v2(A) = v2(B)
#         => gcd(X,Y) = gcd(A,B)
#
#     v2(A) != v2(B)
#         => gcd(X,Y) = gcd(A,B)/2
#
# Special branch:
#
#     FRAME B, B=0:
#
#         X=Y=3A/2
#
#         gcd(X,Y)=3|A|/2.
#
# This experiment:
#
#   1. derives the exact gcd identities symbolically;
#   2. isolates the powers of two;
#   3. isolates the odd part;
#   4. proves the factor-3 issue is impossible generically;
#   5. verifies every state;
#   6. verifies the depth theorem;
#   7. tests for any additional exceptional states.
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

def prime_sieve(limit: int) -> list[int]:

    sieve = bytearray(
        b"\x01"
    ) * (limit + 1)

    sieve[0] = 0
    sieve[1] = 0

    for p in range(
        2,
        isqrt(limit) + 1,
    ):

        if sieve[p]:

            start = p * p

            sieve[
                start:
                limit + 1:
                p
            ] = b"\x00" * (
                ((limit - start) // p) + 1
            )

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
# SEMIPRIMES
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

def frame_of(n: int) -> str:

    r = n & 3

    if r == 3:
        return "A"

    if r == 1:
        return "B"

    raise ValueError(
        f"odd semiprime expected, got n mod 4={r}"
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

        return (
            (q - p + 6) // 2,
            (p + q) // 2,
        )

    return (
        (3 * p - q + 6) // 2,
        (3 * p + q) // 2,
    )


# ==============================================================================
# DEPTH
# ==============================================================================

def depth_xy(
    X: int,
    Y: int,
) -> int:

    return 1 + min(
        v2(X),
        v2(Y),
    )


# ==============================================================================
# ODD PART
# ==============================================================================

def odd_part(x: int) -> int:

    if x == 0:
        return 0

    x = abs(x)

    return x >> v2(x)


# ==============================================================================
# CLASSIFICATION
# ==============================================================================

def branch(
    frame: str,
    A: int,
    B: int,
) -> str:

    if A == 0 and B == 0:
        return "A=B=0"

    if A == 0:
        return "A=0"

    if B == 0:
        return "B=0"

    if v2(A) == v2(B):
        return "equal-v2"

    return "unequal-v2"


# ==============================================================================
# TEST 0
# ==============================================================================

def test_algebraic_coordinate_maps(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 0: EXACT LINEAR COORDINATE MAPS")
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

        if frame == "A":

            ex = (B - A) // 2
            ey = (A + B) // 2

        else:

            ex = (3 * A - B) // 2
            ey = (3 * A + B) // 2

        if (X, Y) != (ex, ey):

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch"
                    f" n={state.n}"
                    f" frame={frame}"
                    f" A={A}"
                    f" B={B}"
                    f" X={X}"
                    f" Y={Y}"
                    f" expected=({ex},{ey})"
                )

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 1
# ==============================================================================

def test_frame_A_gcd_identity(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 1: FRAME-A gcd IDENTITY")
    print("=" * 90)

    failures = 0
    checked = 0

    for state in states:

        if frame_of(state.n) != "A":
            continue

        A, B = raw_residuals(
            "A",
            state.p,
            state.q,
        )

        if A == 0 or B == 0:
            continue

        checked += 1

        X, Y = global_xy(
            "A",
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
                    f" A={A}"
                    f" B={B}"
                    f" gcdAB={d}"
                    f" gcdXY={g}"
                    f" alpha={alpha}"
                    f" beta={beta}"
                    f" expected={expected}"
                )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 2
# ==============================================================================

def test_frame_B_gcd_identity(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 2: FRAME-B gcd IDENTITY")
    print("=" * 90)

    failures = 0
    checked = 0
    odd_three_cases = 0

    for state in states:

        if frame_of(state.n) != "B":
            continue

        A, B = raw_residuals(
            "B",
            state.p,
            state.q,
        )

        if A == 0 or B == 0:
            continue

        checked += 1

        X, Y = global_xy(
            "B",
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

        # Algebra:
        #
        #   2X = 3A-B
        #   2Y = 3A+B
        #
        # Let d=gcd(A,B), A=d*a, B=d*b.
        #
        # Then:
        #
        #   gcd(3A-B, 3A+B)
        #     = d * gcd(3a-b,3a+b)
        #
        # with gcd(a,b)=1.
        #
        # Since a,b are coprime and A,B are both nonzero,
        # the only odd common divisor that could appear is 3.
        #
        # But if 3 divides both 3a-b and 3a+b, then:
        #
        #   3 | b
        #
        # and gcd(a,b)=1 implies 3 ∤ a.
        #
        # That would force B itself to carry the 3 factor.
        #
        # The tested domain will determine whether such a state
        # is possible under the particular frame-B factor structure.

        odd_g = odd_part(g)
        odd_d = odd_part(d)

        if (
            odd_g != odd_d
        ):

            odd_three_cases += 1

            if odd_g % 3 != 0:

                failures += 1

                if failures <= 20:

                    print(
                        "    unexpected odd gcd change"
                        f" n={state.n}"
                        f" A={A}"
                        f" B={B}"
                        f" d={d}"
                        f" g={g}"
                    )

        if alpha == beta:
            expected = d
        else:
            expected = d // 2

        if g != expected:

            failures += 1

            if failures <= 20:

                print(
                    "    gcd mismatch"
                    f" n={state.n}"
                    f" A={A}"
                    f" B={B}"
                    f" d={d}"
                    f" g={g}"
                    f" alpha={alpha}"
                    f" beta={beta}"
                    f" expected={expected}"
                )

    print(
        f"checked={checked}"
    )

    print(
        f"odd-part changes={odd_three_cases}"
    )

    print(
        f"failures={failures}"
    )


# ==============================================================================
# TEST 3
# ==============================================================================

def test_exact_common_divisor(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 3: gcd(TRANSFORMED NUMERATORS)")
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

        if A == 0 or B == 0:
            continue

        checked += 1

        d = gcd(
            abs(A),
            abs(B),
        )

        a = A // d
        b = B // d

        if frame == "A":

            u = B - A
            v = A + B

        else:

            u = 3 * A - B
            v = 3 * A + B

        h = gcd(
            abs(u),
            abs(v),
        )

        normalized = gcd(
            abs(
                u // d
            ),
            abs(
                v // d
            ),
        )

        # For the generic branches the transformed numerators
        # must have gcd 2 or 4.
        #
        # Then dividing both by 2 gives gcd 1 or 2.

        if frame == "A":

            expected_normalized = (
                2
                if v2(A) != v2(B)
                else 4
            )

        else:

            expected_normalized = (
                2
                if v2(A) != v2(B)
                else 4
            )

        if normalized != expected_normalized:

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch"
                    f" n={state.n}"
                    f" frame={frame}"
                    f" A={A}"
                    f" B={B}"
                    f" a={a}"
                    f" b={b}"
                    f" transformed=({u},{v})"
                    f" normalized_gcd={normalized}"
                    f" expected={expected_normalized}"
                )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 4
# ==============================================================================

def test_degenerate_branches(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 4: COMPLETE DEGENERATE BRANCH CLASSIFICATION")
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

        cls = branch(
            frame,
            A,
            B,
        )

        counts[
            (frame, cls)
        ] += 1

        X, Y = global_xy(
            frame,
            state.p,
            state.q,
        )

        g = gcd(
            abs(X),
            abs(Y),
        )

        if cls == "A=0":

            expected = abs(B) // 2

        elif cls == "B=0":

            expected = 3 * abs(A) // 2

        elif cls == "A=B":

            expected = abs(A)

        elif cls == "A=-B":

            expected = abs(A)

        elif cls == "equal-v2":

            d = gcd(
                abs(A),
                abs(B),
            )

            expected = d

        elif cls == "unequal-v2":

            d = gcd(
                abs(A),
                abs(B),
            )

            expected = d // 2

        else:

            continue

        if g != expected:

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch"
                    f" n={state.n}"
                    f" frame={frame}"
                    f" class={cls}"
                    f" A={A}"
                    f" B={B}"
                    f" g={g}"
                    f" expected={expected}"
                )

    for key in sorted(
        counts
    ):

        print(
            f"    {key}: "
            f"{counts[key]}"
        )

    print(
        f"failures={failures}"
    )


# ==============================================================================
# TEST 5
# ==============================================================================

def test_depth_theorem(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 5: EXACT DEPTH THEOREM")
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

        actual = depth_xy(
            X,
            Y,
        )

        alpha = v2(A)
        beta = v2(B)

        if A == 0:

            predicted = (
                beta
            )

        elif B == 0:

            # This special branch is only valid for frame B
            # in the present odd-prime domain.

            predicted = (
                alpha - 1
                + 1
                + v2(3)
            )

            # Simplifies to alpha.

        elif alpha == beta:

            predicted = alpha + 1

        else:

            predicted = min(
                alpha,
                beta,
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
                    f" alpha={alpha}"
                    f" beta={beta}"
                    f" predicted={predicted}"
                    f" actual={actual}"
                )

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 6
# ==============================================================================

def test_single_exception(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 6: B=0 EXCEPTION UNIQUENESS")
    print("=" * 90)

    b_zero = []

    for state in states:

        frame = frame_of(
            state.n
        )

        A, B = raw_residuals(
            frame,
            state.p,
            state.q,
        )

        if B == 0:

            b_zero.append(
                (
                    state.n,
                    state.p,
                    state.q,
                    frame,
                    A,
                )
            )

    print(
        f"B=0 states={len(b_zero)}"
    )

    for record in b_zero:

        print(
            "    B=0:"
            f" n={record[0]}"
            f" p={record[1]}"
            f" q={record[2]}"
            f" frame={record[3]}"
            f" A={record[4]}"
        )

    failures = 0

    if len(b_zero) != 1:

        failures = 1

    elif b_zero[0][:4] != (
        9,
        3,
        3,
        "B",
    ):

        failures = 1

    print(
        f"uniqueness failures={failures}"
    )


# ==============================================================================
# TEST 7
# ==============================================================================

def test_unified_depth_from_gcd(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 7: UNIFIED gcd -> DEPTH")
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

        gxy = gcd(
            abs(X),
            abs(Y),
        )

        actual = (
            1 + v2(gxy)
        )

        if A == 0:

            expected = (
                1
                + v2(
                    abs(B) // 2
                )
            )

        elif B == 0:

            expected = (
                1
                + v2(
                    3 * abs(A) // 2
                )
            )

        else:

            gab = gcd(
                abs(A),
                abs(B),
            )

            correction = (
                0
                if v2(A) == v2(B)
                else -1
            )

            expected = (
                1
                + v2(gab)
                + correction
            )

        if actual != expected:

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch"
                    f" n={state.n}"
                    f" frame={frame}"
                    f" A={A}"
                    f" B={B}"
                    f" gcdAB={gcd(abs(A),abs(B))}"
                    f" gcdXY={gxy}"
                    f" expected_depth={expected}"
                    f" actual_depth={actual}"
                )

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 8
# ==============================================================================

def test_factor_three_symbolically(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 8: SYMBOLIC FACTOR-3 CONDITION")
    print("=" * 90)

    failures = 0
    checked = 0
    unexpected = []

    for state in states:

        if frame_of(state.n) != "B":
            continue

        A, B = raw_residuals(
            "B",
            state.p,
            state.q,
        )

        if B == 0:
            continue

        checked += 1

        d = gcd(
            abs(A),
            abs(B),
        )

        a = A // d
        b = B // d

        # gcd(a,b)=1.
        #
        # If 3 divides both transformed normalized
        # coordinates:
        #
        #     3a-b
        #     3a+b
        #
        # then adding/subtracting gives:
        #
        #     3 | 2b
        #
        # hence:
        #
        #     3 | b.
        #
        # But if b is divisible by 3 and gcd(a,b)=1,
        # then 3 does not divide a.
        #
        # This branch is mathematically possible in arbitrary
        # coprime (a,b), but the present residual construction
        # may constrain it further.
        #
        # Check exact occurrences.

        u = 3 * a - b
        v = 3 * a + b

        h = gcd(
            abs(u),
            abs(v),
        )

        if h % 3 == 0:

            unexpected.append(
                (
                    state,
                    A,
                    B,
                    a,
                    b,
                    h,
                )
            )

    print(
        f"checked={checked}"
    )

    print(
        f"unexpected 3-factor states="
        f"{len(unexpected)}"
    )

    for (
        state,
        A,
        B,
        a,
        b,
        h,
    ) in unexpected[:20]:

        print(
            "    unexpected:"
            f" n={state.n}"
            f" p={state.p}"
            f" q={state.q}"
            f" A={A}"
            f" B={B}"
            f" a={a}"
            f" b={b}"
            f" transformed_gcd={h}"
        )

    if unexpected:
        failures = len(
            unexpected
        )

    print(
        f"failures={failures}"
    )


# ==============================================================================
# TEST 9
# ==============================================================================

def test_depth_distribution(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 9: DEPTH DISTRIBUTION")
    print("=" * 90)

    counts = Counter()

    for state in states:

        frame = frame_of(
            state.n
        )

        X, Y = global_xy(
            frame,
            state.p,
            state.q,
        )

        counts[
            depth_xy(X, Y)
        ] += 1

    for depth in sorted(
        counts
    ):

        print(
            f"    depth={depth}: "
            f"{counts[depth]}"
        )


# ==============================================================================
# EXAMPLES
# ==============================================================================

def print_examples(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("EXAMPLES")
    print("=" * 90)

    targets = [
        9,
        15,
        21,
        33,
        39,
        51,
        57,
        69,
        77,
        87,
        93,
        141,
        183,
        213,
    ]

    lookup = {
        s.n: s
        for s in states
    }

    for n in targets:

        state = lookup.get(n)

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

        gab = gcd(
            abs(A),
            abs(B),
        )

        gxy = gcd(
            abs(X),
            abs(Y),
        )

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
            f"    gcdAB={gab}"
        )

        print(
            f"    gcdXY={gxy}"
        )

        print(
            f"    v2A={v2(A)}"
            f" v2B={v2(B)}"
        )

        print(
            f"    branch="
            f"{branch(frame,A,B)}"
        )

        print(
            f"    depth="
            f"{depth_xy(X,Y)}"
        )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print("EXPERIMENT 632 START")
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

    test_algebraic_coordinate_maps(
        states
    )

    test_frame_A_gcd_identity(
        states
    )

    test_frame_B_gcd_identity(
        states
    )

    test_exact_common_divisor(
        states
    )

    test_degenerate_branches(
        states
    )

    test_depth_theorem(
        states
    )

    test_single_exception(
        states
    )

    test_unified_depth_from_gcd(
        states
    )

    test_factor_three_symbolically(
        states
    )

    test_depth_distribution(
        states
    )

    print_examples(
        states
    )

    print()
    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)

    print(
r"""
For the tested odd semiprime domain:

    depth = 1 + v2(gcd(X,Y))

EXACT GENERIC GCD TRANSFORMATION:

    let

        A,B = raw factor residuals

    and

        alpha = v2(A)
        beta  = v2(B).

    For A,B != 0:

        alpha = beta
            =>
        gcd(X,Y) = gcd(A,B)

        alpha != beta
            =>
        gcd(X,Y) = gcd(A,B)/2.

FRAME A:

    X = (B-A)/2
    Y = (A+B)/2

FRAME B:

    X = (3A-B)/2
    Y = (3A+B)/2

The only observed B=0 state is:

    p=q=3
    n=9
    frame=B
    A=4
    B=0

where:

    X=Y=6

and therefore:

    gcd(X,Y)=6=3A/2.

If every test passes, the complete integer theorem is:

    gcd(X,Y)
      =
    gcd(A,B) / 2^[v2(A) != v2(B)]

for nonzero A,B,

plus exactly one degenerate B=0 branch.

Consequently:

    depth
      =
    1 + v2(gcd(A,B))
      - [v2(A) != v2(B)]

for A,B != 0,

and the zero branches are handled explicitly.

This is the integer-level counterpart of the earlier
2-adic residual automaton.
"""
    )

    print()
    print("=" * 90)
    print("EXPERIMENT 632 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
