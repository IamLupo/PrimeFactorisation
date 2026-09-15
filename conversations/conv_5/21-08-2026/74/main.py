#!/usr/bin/env python3

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import gcd, isqrt


# ==============================================================================
# EXPERIMENT 636
# ==============================================================================
#
# SUM / DISCRIMINANT gcd COLLAPSE
#
# Experiment 633/634 established:
#
#     d = gcd(A,B)
#
#     depth = 1 + v2(gcd(X,Y))
#
# and:
#
# FRAME A:
#     A = p-3
#     B = q+3
#
# FRAME B:
#     A = p+1
#     B = q-3
#
# The next structural question is whether d can be expressed through
# the symmetric root data:
#
#     S = p+q
#     D = p-q
#
# because:
#
# FRAME A:
#
#     A+B = p+q       = S
#     B-A = q-p+6     = -D+6
#
#     d = gcd(S, D-6)
#
# FRAME B:
#
#     A+B = p+q-2     = S-2
#     B-A = q-p-4     = -D-4
#
#     d = gcd(S-2, D+4)
#
# But D itself satisfies:
#
#     D^2 = S^2 - 4n.
#
# Therefore we can eliminate the signed D from a gcd using:
#
#     (D-c)(D+c) = D^2-c^2.
#
# The experiment tests whether this produces an exact integer-level
# description of d up to only tiny, explicitly characterized factors.
#
# It also tests whether:
#
#     d = gcd(linear sum form, discriminant expression)
#
# and whether the remaining factor is completely explained by 2 and 3.
#
# No large offset search.
# No level traversal.
#
# ==============================================================================


PRIME_LIMIT = 6250
INF = 10**9


@dataclass(frozen=True)
class State:
    n: int
    p: int
    q: int


# ==============================================================================
# NUMBER THEORY
# ==============================================================================

def v2(x: int) -> int:

    if x == 0:
        return INF

    x = abs(x)

    return (x & -x).bit_length() - 1


def sieve(limit: int) -> list[int]:

    a = bytearray(b"\x01") * (limit + 1)

    a[0] = 0
    a[1] = 0

    for p in range(
        2,
        isqrt(limit) + 1,
    ):

        if a[p]:

            start = p * p

            a[
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
        if a[p]
    ]


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
# FRAME / RESIDUALS
# ==============================================================================

def frame_of(n: int) -> str:

    r = n & 3

    if r == 3:
        return "A"

    if r == 1:
        return "B"

    raise ValueError(
        f"Unexpected odd semiprime residue: "
        f"n={n}, n mod 4={r}"
    )


def residuals(
    frame: str,
    p: int,
    q: int,
) -> tuple[int, int]:

    if frame == "A":
        return p - 3, q + 3

    return p + 1, q - 3


# ==============================================================================
# GLOBAL COORDINATES
# ==============================================================================

def xy(
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


def observed_depth(
    frame: str,
    p: int,
    q: int,
) -> int:

    X, Y = xy(
        frame,
        p,
        q,
    )

    return 1 + v2(
        gcd(
            abs(X),
            abs(Y),
        )
    )


# ==============================================================================
# TEST 0
# ==============================================================================

def test_exact_linear_gcd(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 0: d = gcd(SHIFTED SUM, SHIFTED DIFFERENCE)")
    print("=" * 90)

    failures = 0
    checked = 0

    for s in states:

        frame = frame_of(s.n)

        A, B = residuals(
            frame,
            s.p,
            s.q,
        )

        # Skip the degenerate zero branch here.
        if A == 0 or B == 0:
            continue

        d = gcd(
            abs(A),
            abs(B),
        )

        S = s.p + s.q
        D = s.p - s.q

        if frame == "A":

            candidate = gcd(
                abs(S),
                abs(D - 6),
            )

        else:

            candidate = gcd(
                abs(S - 2),
                abs(D + 4),
            )

        checked += 1

        if candidate != d:

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch:"
                    f" n={s.n}"
                    f" frame={frame}"
                    f" A={A}"
                    f" B={B}"
                    f" d={d}"
                    f" candidate={candidate}"
                )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 1
# ==============================================================================

def test_discriminant_gcd(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 1: gcd(S, D-c) -> "
        "DISCRIMINANT-ONLY FORM"
    )
    print("=" * 90)

    #
    # Frame A:
    #
    #     d = gcd(S, D-6)
    #
    # Since:
    #
    #     (D-6)(D+6)
    #       = D^2 - 36
    #
    # and:
    #
    #     D^2 = S^2 - 4n
    #
    # define:
    #
    #     R_A = S^2 - 4n - 36.
    #
    # Any d divides both S and R_A.
    #
    # Frame B:
    #
    #     d = gcd(S-2, D+4)
    #
    # and:
    #
    #     (D+4)(D-4)
    #       = D^2 - 16
    #
    # so:
    #
    #     R_B = S^2 - 4n - 16.
    #
    # We test the resulting gcd.
    #

    failures = 0
    checked = 0

    ratio_distribution = Counter()

    for s in states:

        frame = frame_of(s.n)

        A, B = residuals(
            frame,
            s.p,
            s.q,
        )

        if A == 0 or B == 0:
            continue

        d = gcd(
            abs(A),
            abs(B),
        )

        S = s.p + s.q

        if frame == "A":

            U = S
            R = (
                S * S
                - 4 * s.n
                - 36
            )

        else:

            U = S - 2
            R = (
                S * S
                - 4 * s.n
                - 16
            )

        g = gcd(
            abs(U),
            abs(R),
        )

        if d == 0:
            continue

        checked += 1

        if g % d != 0:

            failures += 1

            if failures <= 20:

                print(
                    "    d does not divide "
                    "discriminant gcd:"
                    f" n={s.n}"
                    f" frame={frame}"
                    f" d={d}"
                    f" g={g}"
                )

            continue

        ratio = g // d

        ratio_distribution[
            ratio
        ] += 1

    print(
        f"checked={checked} "
        f"failures={failures}"
    )

    print()
    print("    gcd/discriminant ratio distribution:")

    for ratio, count in sorted(
        ratio_distribution.items()
    ):

        print(
            f"        {ratio}: {count}"
        )


# ==============================================================================
# TEST 2
# ==============================================================================

def test_ratio_prime_content(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 2: IS THE DISCRIMINANT RATIO "
        "ONLY A SMALL CONSTANT?"
    )
    print("=" * 90)

    ratios = Counter()

    for s in states:

        frame = frame_of(s.n)

        A, B = residuals(
            frame,
            s.p,
            s.q,
        )

        if A == 0 or B == 0:
            continue

        d = gcd(
            abs(A),
            abs(B),
        )

        S = s.p + s.q

        if frame == "A":

            U = S

            R = (
                S * S
                - 4 * s.n
                - 36
            )

        else:

            U = S - 2

            R = (
                S * S
                - 4 * s.n
                - 16
            )

        g = gcd(
            abs(U),
            abs(R),
        )

        ratios[
            g // d
        ] += 1

    print(
        f"    distinct ratios={len(ratios)}"
    )

    print()

    for ratio, count in sorted(ratios.items()):

        print(
            f"    ratio={ratio}"
            f" count={count}"
        )


# ==============================================================================
# TEST 3
# ==============================================================================

def test_2_3_normalization(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 3: REMOVE POSSIBLE 2/3 CONTENT"
    )
    print("=" * 90)

    failures = 0
    checked = 0

    def strip_small(x: int) -> int:

        x = abs(x)

        while x and x % 2 == 0:
            x //= 2

        while x and x % 3 == 0:
            x //= 3

        return x

    for s in states:

        frame = frame_of(s.n)

        A, B = residuals(
            frame,
            s.p,
            s.q,
        )

        if A == 0 or B == 0:
            continue

        d = gcd(
            abs(A),
            abs(B),
        )

        S = s.p + s.q

        if frame == "A":

            U = S

            R = (
                S * S
                - 4 * s.n
                - 36
            )

        else:

            U = S - 2

            R = (
                S * S
                - 4 * s.n
                - 16
            )

        g = gcd(
            abs(U),
            abs(R),
        )

        checked += 1

        if strip_small(g) != strip_small(d):

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch:"
                    f" n={s.n}"
                    f" frame={frame}"
                    f" d={d}"
                    f" g={g}"
                    f" reduced_d={strip_small(d)}"
                    f" reduced_g={strip_small(g)}"
                )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 4
# ==============================================================================

def test_sum_only_depth(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 4: DEPTH FROM SYMMETRIC ROOT DATA"
    )
    print("=" * 90)

    failures = 0
    checked = 0

    for s in states:

        frame = frame_of(s.n)

        A, B = residuals(
            frame,
            s.p,
            s.q,
        )

        if A == 0 or B == 0:
            continue

        d = gcd(
            abs(A),
            abs(B),
        )

        S = s.p + s.q
        D = s.p - s.q

        if frame == "A":

            g = gcd(
                abs(S),
                abs(D - 6),
            )

        else:

            g = gcd(
                abs(S - 2),
                abs(D + 4),
            )

        predicted = 1 + v2(g)

        actual = observed_depth(
            frame,
            s.p,
            s.q,
        )

        checked += 1

        if predicted != actual:

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch:"
                    f" n={s.n}"
                    f" frame={frame}"
                    f" predicted={predicted}"
                    f" actual={actual}"
                    f" d={d}"
                    f" g={g}"
                )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 5
# ==============================================================================

def test_discriminant_only_depth(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 5: DEPTH FROM SUM + n "
        "WITHOUT SIGNED DIFFERENCE"
    )
    print("=" * 90)

    failures = 0
    checked = 0

    ratio = Counter()

    for s in states:

        frame = frame_of(s.n)

        A, B = residuals(
            frame,
            s.p,
            s.q,
        )

        if A == 0 or B == 0:
            continue

        S = s.p + s.q

        if frame == "A":

            U = S

            R = (
                S * S
                - 4 * s.n
                - 36
            )

        else:

            U = S - 2

            R = (
                S * S
                - 4 * s.n
                - 16
            )

        g = gcd(
            abs(U),
            abs(R),
        )

        d = gcd(
            abs(A),
            abs(B),
        )

        ratio[
            g // d
        ] += 1

        checked += 1

        predicted = 1 + v2(d)

        # The test here deliberately asks whether the extra gcd
        # does at least preserve the relevant 2-adic valuation.
        #
        # In other words:
        #
        #     v2(g) == v2(d)
        #
        # If true, then the discriminant-only gcd already preserves
        # the exact depth even if its odd part is enlarged.

        if v2(g) != v2(d):

            failures += 1

            if failures <= 20:

                print(
                    "    v2 mismatch:"
                    f" n={s.n}"
                    f" frame={frame}"
                    f" d={d}"
                    f" g={g}"
                    f" v2d={v2(d)}"
                    f" v2g={v2(g)}"
                    f" depth={predicted}"
                )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )

    print()
    print("    g/d ratios:")

    for k, count in sorted(ratio.items()):

        print(
            f"        {k}: {count}"
        )


# ==============================================================================
# TEST 6
# ==============================================================================

def test_degenerate_branches(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 6: DEGENERATE BRANCHES")
    print("=" * 90)

    failures = 0
    checked = 0

    for s in states:

        frame = frame_of(s.n)

        A, B = residuals(
            frame,
            s.p,
            s.q,
        )

        if A != 0 and B != 0:
            continue

        X, Y = xy(
            frame,
            s.p,
            s.q,
        )

        actual = 1 + v2(
            gcd(
                abs(X),
                abs(Y),
            )
        )

        if frame == "A" and A == 0:

            predicted = 1 + v2(
                abs(B) // 2
            )

        elif frame == "B" and B == 0:

            predicted = 1 + v2(
                3 * abs(A) // 2
            )

        else:

            print(
                "    unexpected degenerate:"
                f" n={s.n}"
                f" frame={frame}"
                f" A={A}"
                f" B={B}"
            )

            failures += 1
            continue

        checked += 1

        if predicted != actual:

            failures += 1

            print(
                "    mismatch:"
                f" n={s.n}"
                f" frame={frame}"
                f" A={A}"
                f" B={B}"
                f" predicted={predicted}"
                f" actual={actual}"
            )

    print(
        f"checked={checked} "
        f"failures={failures}"
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

    wanted = [
        9,
        15,
        21,
        33,
        39,
        57,
        69,
        77,
        93,
        141,
        183,
        213,
    ]

    by_n = {
        s.n: s
        for s in states
    }

    for n in wanted:

        s = by_n[n]

        frame = frame_of(n)

        A, B = residuals(
            frame,
            s.p,
            s.q,
        )

        X, Y = xy(
            frame,
            s.p,
            s.q,
        )

        actual_d = gcd(
            abs(A),
            abs(B),
        )

        S = s.p + s.q
        D = s.p - s.q

        if frame == "A":

            U = S
            R = (
                S * S
                - 4 * n
                - 36
            )

        else:

            U = S - 2
            R = (
                S * S
                - 4 * n
                - 16
            )

        discriminant_gcd = gcd(
            abs(U),
            abs(R),
        )

        print()
        print(
            f"n={n}"
            f" p={s.p}"
            f" q={s.q}"
            f" frame={frame}"
        )

        print(
            f"    A={A}"
            f" B={B}"
            f" d={actual_d}"
        )

        print(
            f"    S=p+q={S}"
            f" D=p-q={D}"
        )

        print(
            f"    X={X}"
            f" Y={Y}"
        )

        print(
            f"    shifted-sum gcd={gcd(abs(U),abs(D-6)) if frame=='A' else gcd(abs(U),abs(D+4))}"
        )

        print(
            f"    U={U}"
            f" R={R}"
        )

        print(
            f"    gcd(U,R)={discriminant_gcd}"
        )

        print(
            f"    v2(d)={v2(actual_d)}"
            f" v2(gcd(U,R))={v2(discriminant_gcd)}"
        )

        print(
            f"    depth={observed_depth(frame,s.p,s.q)}"
        )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print("EXPERIMENT 636 START")
    print("=" * 90)

    print()
    print("SUM / DISCRIMINANT gcd COLLAPSE")

    print()
    print("[1] PRIME SIEVE")

    primes = sieve(
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

    test_exact_linear_gcd(
        states
    )

    test_discriminant_gcd(
        states
    )

    test_ratio_prime_content(
        states
    )

    test_2_3_normalization(
        states
    )

    test_sum_only_depth(
        states
    )

    test_discriminant_only_depth(
        states
    )

    test_degenerate_branches(
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
The residual gcd is:

    d = gcd(A,B).

Experiment 633 gave:

    FRAME A:

        A=p-3
        B=q+3

        X=(B-A)/2
        Y=(A+B)/2

    FRAME B:

        A=p+1
        B=q-3

        X=(3A-B)/2
        Y=(3A+B)/2.

A symmetric reformulation is:

    FRAME A:

        d = gcd(
                p+q,
                p-q-6
            )

    FRAME B:

        d = gcd(
                p+q-2,
                p-q+4
            ).

Now define the discriminant:

    Delta = (p-q)^2
           = (p+q)^2 - 4n.

Therefore:

    FRAME A:

        (p-q-6)(p-q+6)
          =
        (p+q)^2 - 4n - 36.

    FRAME B:

        (p-q+4)(p-q-4)
          =
        (p+q)^2 - 4n - 16.

So the experiment tests whether the gcd can be recovered,
up to a harmless small factor, from:

    sum + n + discriminant.

The particularly important result is:

    v2(gcd(U,R)) == v2(d)

which would mean that the complete depth is already
recoverable from the symmetric integer data:

    S=p+q
    n
    frame.

This is materially different from the previous arbitrary
n-only valuation searches.

The possible structural endpoint is:

    depth
      =
    1 + v2(
        gcd(
            S-shift,
            S^2 - 4n - constant
        )
      )

with the small constants determined solely by FRAME.

That would provide a clean quadratic/discriminant
description of the hierarchy.
"""
    )

    print()
    print("=" * 90)
    print("EXPERIMENT 636 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
