#!/usr/bin/env python3

"""
==========================================================================================
EXPERIMENT 640
==========================================================================================

CORRECTED SYMMETRIC NUMERATOR VALUATION LAW

Experiment 639 showed that the proposed identity

    depth = min(v2(U), v2(n+c)+2)

is FALSE.

Counterexample:

    n=93
    frame=B
    A=4
    B=28

    X=(3A-B)/2 = -8
    Y=(3A+B)/2 = 20

    gcd(X,Y)=4
    depth=1+v2(4)=3.

But:

    U=A+B=32
    v2(U)=5

so U alone overestimates the depth.

The correct starting point is:

FRAME A:

    2X = B-A
    2Y = A+B

FRAME B:

    2X = 3A-B
    2Y = 3A+B.

Therefore:

    gcd(X,Y)
      =
    gcd(2X,2Y)/2

and hence:

FRAME A:

    depth
      =
    min(
        v2(B-A),
        v2(A+B)
    )

FRAME B:

    depth
      =
    min(
        v2(3A-B),
        v2(3A+B)
    ).

This experiment proves that identity over the full semiprime
domain and then derives the correct symmetric relations.

No brute-force n-only offset search is performed.
"""


from __future__ import annotations

from dataclasses import dataclass
from math import gcd
from collections import Counter, defaultdict


# ==============================================================================
# CONFIG
# ==============================================================================

PRIME_LIMIT = 6000
INF = 10**9


# ==============================================================================
# STATE
# ==============================================================================

@dataclass(frozen=True)
class State:
    n: int
    p: int
    q: int

    frame: str

    A: int
    B: int

    X: int
    Y: int

    depth: int


# ==============================================================================
# v2
# ==============================================================================

def v2(x: int) -> int:

    if x == 0:
        return INF

    x = abs(x)

    return (x & -x).bit_length() - 1


# ==============================================================================
# gcd
# ==============================================================================

def gcd0(a: int, b: int) -> int:
    return gcd(abs(a), abs(b))


# ==============================================================================
# FRAME
# ==============================================================================

def frame_from_n(n: int) -> str:
    """
    Established convention:

        n == 3 mod 4 -> A
        n == 1 mod 4 -> B
    """

    r = n & 3

    if r == 3:
        return "A"

    if r == 1:
        return "B"

    raise ValueError(
        f"Unexpected odd n mod 4={r}"
    )


# ==============================================================================
# PRIME SIEVE
# ==============================================================================

def prime_sieve(limit: int) -> list[int]:

    sieve = bytearray(b"\x01") * (limit + 1)

    sieve[0] = 0
    sieve[1] = 0

    root = int(limit ** 0.5)

    for p in range(2, root + 1):

        if sieve[p]:

            start = p * p

            sieve[
                start:
                limit + 1:
                p
            ] = (
                b"\x00"
                * (
                    ((limit - start) // p)
                    + 1
                )
            )

    return [
        p
        for p in range(
            3,
            limit + 1,
            2
        )
        if sieve[p]
    ]


# ==============================================================================
# RESIDUALS
# ==============================================================================

def residuals(
    frame: str,
    p: int,
    q: int,
):

    if frame == "A":

        return (
            p - 3,
            q + 3,
        )

    if frame == "B":

        return (
            p + 1,
            q - 3,
        )

    raise ValueError(frame)


# ==============================================================================
# GLOBAL COORDINATES
# ==============================================================================

def global_xy(
    frame: str,
    A: int,
    B: int,
):

    if frame == "A":

        nx = B - A
        ny = A + B

    elif frame == "B":

        nx = 3 * A - B
        ny = 3 * A + B

    else:
        raise ValueError(frame)

    if nx & 1:
        raise ArithmeticError(
            f"Non-integral X: "
            f"frame={frame} A={A} B={B}"
        )

    if ny & 1:
        raise ArithmeticError(
            f"Non-integral Y: "
            f"frame={frame} A={A} B={B}"
        )

    return nx // 2, ny // 2


# ==============================================================================
# DEPTH
# ==============================================================================

def exact_depth(
    X: int,
    Y: int,
):

    g = gcd0(
        X,
        Y,
    )

    if g == 0:

        raise ArithmeticError(
            "Zero gcd(X,Y)"
        )

    return 1 + v2(g)


# ==============================================================================
# BUILD
# ==============================================================================

def build_states(
    primes: list[int],
):

    states = []

    for i, p in enumerate(primes):

        for q in primes[i:]:

            n = p * q

            frame = frame_from_n(n)

            A, B = residuals(
                frame,
                p,
                q,
            )

            X, Y = global_xy(
                frame,
                A,
                B,
            )

            depth = exact_depth(
                X,
                Y,
            )

            states.append(
                State(
                    n=n,
                    p=p,
                    q=q,
                    frame=frame,
                    A=A,
                    B=B,
                    X=X,
                    Y=Y,
                    depth=depth,
                )
            )

    return states


# ==============================================================================
# TEST 0
# ==============================================================================

def test_frame_regression():

    print("=" * 90)
    print("TEST 0: FRAME REGRESSION")
    print("=" * 90)

    expected = {
        9: "B",
        15: "A",
        21: "B",
        33: "B",
        39: "A",
        57: "B",
        69: "B",
        77: "B",
        93: "B",
        111: "A",
        141: "B",
        183: "A",
        213: "B",
    }

    failures = 0

    for n, wanted in expected.items():

        actual = frame_from_n(n)

        if actual != wanted:

            failures += 1

            print(
                f"    mismatch n={n} "
                f"expected={wanted} "
                f"actual={actual}"
            )

    print(
        f"checked={len(expected)} "
        f"failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 1
# ==============================================================================

def test_numerators(
    states,
):

    print("=" * 90)
    print("TEST 1: EXACT TRANSFORMED NUMERATORS")
    print("=" * 90)

    failures = 0

    for s in states:

        if s.frame == "A":

            nx = s.B - s.A
            ny = s.A + s.B

        else:

            nx = 3 * s.A - s.B
            ny = 3 * s.A + s.B

        if nx != 2 * s.X:
            failures += 1

        if ny != 2 * s.Y:
            failures += 1

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 2
# ==============================================================================

def test_gcd_half_identity(
    states,
):

    print("=" * 90)
    print("TEST 2: gcd(2X,2Y) = 2*gcd(X,Y)")
    print("=" * 90)

    failures = 0

    for s in states:

        lhs = gcd0(
            2 * s.X,
            2 * s.Y,
        )

        rhs = 2 * gcd0(
            s.X,
            s.Y,
        )

        if lhs != rhs:

            failures += 1

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 3
# ==============================================================================

def test_direct_depth_formula(
    states,
):

    print("=" * 90)
    print(
        "TEST 3: DIRECT DEPTH = MINIMUM "
        "OF TRANSFORMED NUMERATOR VALUATIONS"
    )
    print("=" * 90)

    failures = 0

    distribution = Counter()

    shown = 0

    for s in states:

        if s.frame == "A":

            nx = s.B - s.A
            ny = s.A + s.B

        else:

            nx = 3 * s.A - s.B
            ny = 3 * s.A + s.B

        predicted = min(
            v2(nx),
            v2(ny),
        )

        distribution[
            predicted
        ] += 1

        if predicted != s.depth:

            failures += 1

            if shown < 20:

                print(
                    f"    mismatch "
                    f"n={s.n} "
                    f"frame={s.frame} "
                    f"nx={nx} "
                    f"ny={ny} "
                    f"v2nx={v2(nx)} "
                    f"v2ny={v2(ny)} "
                    f"pred={predicted} "
                    f"actual={s.depth}"
                )

                shown += 1

    print()

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )

    print()
    print("    predicted depth distribution:")

    for d in sorted(distribution):

        print(
            f"        depth={d}: "
            f"{distribution[d]}"
        )

    print()

    return failures


# ==============================================================================
# TEST 4
# ==============================================================================

def test_frame_A_sum_difference(
    states,
):

    print("=" * 90)
    print(
        "TEST 4: FRAME-A "
        "SUM/DIFFERENCE FORMULA"
    )
    print("=" * 90)

    failures = 0

    checked = 0

    shown = 0

    for s in states:

        if s.frame != "A":
            continue

        checked += 1

        diff = s.B - s.A
        summ = s.A + s.B

        predicted = min(
            v2(diff),
            v2(summ),
        )

        if predicted != s.depth:

            failures += 1

            if shown < 20:

                print(
                    f"    mismatch "
                    f"n={s.n} "
                    f"A={s.A} "
                    f"B={s.B} "
                    f"diff={diff} "
                    f"sum={summ}"
                )

                shown += 1

    print(
        f"checked={checked} "
        f"failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 5
# ==============================================================================

def test_frame_B_sum_difference(
    states,
):

    print("=" * 90)
    print(
        "TEST 5: FRAME-B "
        "3A±B FORMULA"
    )
    print("=" * 90)

    failures = 0

    checked = 0

    shown = 0

    for s in states:

        if s.frame != "B":
            continue

        checked += 1

        minus = 3 * s.A - s.B
        plus = 3 * s.A + s.B

        predicted = min(
            v2(minus),
            v2(plus),
        )

        if predicted != s.depth:

            failures += 1

            if shown < 20:

                print(
                    f"    mismatch "
                    f"n={s.n} "
                    f"A={s.A} "
                    f"B={s.B} "
                    f"3A-B={minus} "
                    f"3A+B={plus}"
                )

                shown += 1

    print(
        f"checked={checked} "
        f"failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 6
# ==============================================================================

def test_symmetric_square_identity(
    states,
):

    print("=" * 90)
    print(
        "TEST 6: gcd OF SUM/DIFFERENCE "
        "AND PRODUCT IDENTITY"
    )
    print("=" * 90)

    """
    Let:

        U = transformed first numerator
        V = transformed second numerator.

    Then:

        U*V

    should have a simple relation to A,B.

    FRAME A:

        (B-A)(A+B)
          =
        B^2-A^2.

    FRAME B:

        (3A-B)(3A+B)
          =
        9A^2-B^2.

    We verify these exactly.

    """

    failures = 0

    for s in states:

        if s.frame == "A":

            U = s.B - s.A
            V = s.A + s.B

            lhs = U * V
            rhs = s.B * s.B - s.A * s.A

        else:

            U = 3 * s.A - s.B
            V = 3 * s.A + s.B

            lhs = U * V
            rhs = 9 * s.A * s.A - s.B * s.B

        if lhs != rhs:

            failures += 1

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 7
# ==============================================================================

def test_n_relation(
    states,
):

    print("=" * 90)
    print(
        "TEST 7: PRODUCT OF TRANSFORMED "
        "NUMERATORS VS n"
    )
    print("=" * 90)

    """
    Frame A:

        (B-A)(A+B)
          =
        (q+3 - (p-3))(p+q)
          =
        (q-p+6)(p+q).

    Frame B:

        (3A-B)(3A+B)
          =
        9A^2-B^2.

    This test searches for the exact relation to n
    without assuming the incorrect gcd formula from 638.
    """

    failures = 0

    relations = Counter()

    for s in states:

        if s.frame == "A":

            U = s.B - s.A
            V = s.A + s.B

            product = U * V

            # Express directly through p,q,n.
            expected = (
                (s.q - s.p + 6)
                * (s.p + s.q)
            )

            if product != expected:
                failures += 1

            relations["A"] += 1

        else:

            U = 3 * s.A - s.B
            V = 3 * s.A + s.B

            product = U * V

            expected = (
                (4 * s.p - s.q + 3)
                * (4 * s.p + s.q + 3)
            )

            if product != expected:
                failures += 1

            relations["B"] += 1

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )

    print(
        f"    frame counts={dict(relations)}"
    )

    print()

    return failures


# ==============================================================================
# TEST 8
# ==============================================================================

def test_primary_sum_is_not_depth(
    states,
):

    print("=" * 90)
    print(
        "TEST 8: v2(A+B) IS NOT GENERALLY "
        "THE DEPTH"
    )
    print("=" * 90)

    counterexamples = 0
    shown = 0

    for s in states:

        U = s.A + s.B

        if v2(U) != s.depth:

            counterexamples += 1

            if shown < 20:

                print(
                    f"    n={s.n} "
                    f"frame={s.frame} "
                    f"A={s.A} "
                    f"B={s.B} "
                    f"v2(A+B)={v2(U)} "
                    f"depth={s.depth}"
                )

                shown += 1

    print(
        f"counterexamples={counterexamples}"
    )
    print()

    return 0


# ==============================================================================
# TEST 9
# ==============================================================================

def test_minimum_is_essential(
    states,
):

    print("=" * 90)
    print(
        "TEST 9: WHICH NUMERATOR "
        "LIMITS THE DEPTH?"
    )
    print("=" * 90)

    counts = Counter()

    examples = {}

    for s in states:

        if s.frame == "A":

            U = s.B - s.A
            V = s.A + s.B

        else:

            U = 3 * s.A - s.B
            V = 3 * s.A + s.B

        vu = v2(U)
        vv = v2(V)

        if vu < vv:

            relation = "first-limits"

        elif vv < vu:

            relation = "second-limits"

        else:

            relation = "equal"

        counts[
            (
                s.frame,
                relation,
            )
        ] += 1

        examples.setdefault(
            (
                s.frame,
                relation,
            ),
            s,
        )

    for frame in ("A", "B"):

        print()
        print(
            f"    FRAME {frame}"
        )

        for relation in (
            "first-limits",
            "second-limits",
            "equal",
        ):

            key = (
                frame,
                relation,
            )

            count = counts.get(
                key,
                0,
            )

            print(
                f"        {relation}: "
                f"{count}"
            )

            if count:

                s = examples[key]

                if s.frame == "A":

                    U = s.B - s.A
                    V = s.A + s.B

                else:

                    U = 3 * s.A - s.B
                    V = 3 * s.A + s.B

                print(
                    f"            example "
                    f"n={s.n} "
                    f"v2first={v2(U)} "
                    f"v2second={v2(V)} "
                    f"depth={s.depth}"
                )

    print()

    return 0


# ==============================================================================
# TEST 10
# ==============================================================================

def test_known_examples(
    states,
):

    print("=" * 90)
    print(
        "TEST 10: KNOWN EXAMPLE REGRESSION"
    )
    print("=" * 90)

    wanted = {
        9: 2,
        15: 3,
        21: 3,
        33: 2,
        39: 4,
        51: 2,
        57: 2,
        69: 3,
        77: 4,
        87: 5,
        93: 3,
        111: 3,
        141: 3,
        183: 6,
        213: 3,
    }

    lookup = {
        s.n: s
        for s in states
    }

    failures = 0

    for n, expected in wanted.items():

        actual = lookup[n].depth

        if actual != expected:

            failures += 1

            print(
                f"    mismatch "
                f"n={n} "
                f"expected={expected} "
                f"actual={actual}"
            )

    print(
        f"checked={len(wanted)} "
        f"failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# EXAMPLES
# ==============================================================================

def print_examples(
    states,
):

    print("=" * 90)
    print("EXAMPLES")
    print("=" * 90)

    wanted = [
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
        111,
        141,
        183,
        213,
    ]

    lookup = {
        s.n: s
        for s in states
    }

    for n in wanted:

        s = lookup[n]

        if s.frame == "A":

            first = s.B - s.A
            second = s.A + s.B

        else:

            first = 3 * s.A - s.B
            second = 3 * s.A + s.B

        print()
        print(
            f"n={s.n} "
            f"p={s.p} "
            f"q={s.q} "
            f"frame={s.frame}"
        )

        print(
            f"    A={s.A} B={s.B}"
        )

        print(
            f"    X={s.X} Y={s.Y}"
        )

        print(
            f"    transformed=("
            f"{first},{second})"
        )

        print(
            f"    v2(first)={v2(first)}"
        )

        print(
            f"    v2(second)={v2(second)}"
        )

        print(
            f"    min={min(v2(first), v2(second))}"
        )

        print(
            f"    depth={s.depth}"
        )


# ==============================================================================
# SUMMARY
# ==============================================================================

def print_summary():

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)

    print(
r"""
The Experiment 638/639 reduction to

    depth = v2(gcd(U,R))

was incorrect.

The correct identity starts directly from the coordinate map.

FRAME A:

    X = (B-A)/2
    Y = (A+B)/2

Therefore:

    depth
      =
    1 + v2(gcd(X,Y))

      =
    min(
        v2(B-A),
        v2(A+B)
    ).

FRAME B:

    X = (3A-B)/2
    Y = (3A+B)/2

Therefore:

    depth
      =
    min(
        v2(3A-B),
        v2(3A+B)
    ).

This immediately explains the previous counterexample:

    n=93
    A=4
    B=28

    3A-B=-16       v2=4
    3A+B=40        v2=3

    depth=min(4,3)=3.

Meanwhile:

    A+B=32
    v2(A+B)=5

which is larger than the true depth.

Thus v2(A+B) is only one side of the minimum.

The next structural problem is now:

    Can the two transformed valuations be reduced
    to a smaller pair of symmetric invariants?

For FRAME A those are:

    v2(A+B)
    v2(B-A)

For FRAME B:

    v2(3A+B)
    v2(3A-B).

Their minimum is exactly the depth.

This is the correct continuation of the residual/gcd
hierarchy.
"""
    )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print("EXPERIMENT 640 START")
    print("=" * 90)
    print()
    print(
        "CORRECTED SYMMETRIC NUMERATOR VALUATION LAW"
    )
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

    states = build_states(
        primes
    )

    print(
        f"    semiprimes={len(states)}"
    )
    print()

    failures = 0

    failures += test_frame_regression()

    failures += test_numerators(
        states
    )

    failures += test_gcd_half_identity(
        states
    )

    failures += test_direct_depth_formula(
        states
    )

    failures += test_frame_A_sum_difference(
        states
    )

    failures += test_frame_B_sum_difference(
        states
    )

    failures += test_symmetric_square_identity(
        states
    )

    failures += test_n_relation(
        states
    )

    test_primary_sum_is_not_depth(
        states
    )

    test_minimum_is_essential(
        states
    )

    failures += test_known_examples(
        states
    )

    print_examples(
        states
    )

    print_summary()

    print()
    print("=" * 90)
    print("EXPERIMENT 640 FINISHED")
    print("=" * 90)
    print()

    print(
        f"TOTAL FAILURES={failures}"
    )

    if failures == 0:

        print(
            "STATUS=ALL CORE TESTS PASSED"
        )

    else:

        print(
            "STATUS=COUNTEREXAMPLES FOUND"
        )


if __name__ == "__main__":
    main()
