#!/usr/bin/env python3

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import gcd, isqrt


# ==============================================================================
# EXPERIMENT 634
# ==============================================================================
#
# QUOTIENT-LADDER / gcd-RESIDUAL DECOMPOSITION
#
# Experiment 633 established:
#
#     depth = 1 + v2(gcd(X,Y))
#
# and, for nonzero A,B:
#
#     d = gcd(A,B)
#
#     gcd(X,Y) =
#         d       if v2(A)=v2(B)
#         d / 2   otherwise.
#
# This experiment now focuses on:
#
#     n + c
#
# where:
#
# FRAME A:
#
#     A = p-3
#     B = q+3
#
#     n+9 = AB + 3B - 3A
#
# FRAME B:
#
#     A = p+1
#     B = q-3
#
#     n+3 = AB + 3A - B
#
# Writing:
#
#     A=d*a
#     B=d*b
#     gcd(a,b)=1
#
# gives:
#
# FRAME A:
#
#     (n+9)/d = d*a*b + 3*(b-a)
#
# FRAME B:
#
#     (n+3)/d = d*a*b + 3*a - b
#
# The key question:
#
#     Is the entire depth mechanism visible in the parity / valuation
#     of this quotient?
#
# Expected:
#
#     unequal v2(A),v2(B):
#         quotient is odd
#
#     equal v2(A),v2(B):
#         quotient is even
#
# and therefore:
#
#     v2(n+c) =
#         v2(d)             unequal
#         > v2(d)           equal.
#
# This is exactly the equality-branch mechanism discovered earlier,
# now expressed as an integer quotient identity.
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
# BASIC NUMBER THEORY
# ==============================================================================

def v2(x: int) -> int:

    if x == 0:
        return INF

    return (abs(x) & -abs(x)).bit_length() - 1


def odd_part(x: int) -> int:

    if x == 0:
        return 0

    return abs(x) >> v2(x)


def sieve(limit: int) -> list[int]:

    s = bytearray(b"\x01") * (limit + 1)

    s[0] = 0
    s[1] = 0

    for p in range(
        2,
        isqrt(limit) + 1,
    ):

        if s[p]:

            start = p * p

            s[
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
            2
        )
        if s[p]
    ]


def generate_states(
    primes: list[int],
) -> list[State]:

    out: list[State] = []

    for i, p in enumerate(primes):

        for q in primes[i:]:

            out.append(
                State(
                    n=p * q,
                    p=p,
                    q=q,
                )
            )

    return out


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
        f"Expected odd semiprime: "
        f"n={n}, n mod 4={r}"
    )


# ==============================================================================
# RAW FACTOR RESIDUALS
# ==============================================================================

def residuals(
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
# GLOBAL COORDINATES
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

def global_depth(
    frame: str,
    p: int,
    q: int,
) -> int:

    X, Y = global_xy(
        frame,
        p,
        q,
    )

    return (
        1
        + v2(
            gcd(
                abs(X),
                abs(Y),
            )
        )
    )


# ==============================================================================
# TEST 0
# ==============================================================================

def test_exact_n_shift_identity(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 0: EXACT n-SHIFT IDENTITIES")
    print("=" * 90)

    failures = 0

    for s in states:

        frame = frame_of(s.n)

        A, B = residuals(
            frame,
            s.p,
            s.q,
        )

        if frame == "A":

            expected = (
                A * B
                + 3 * B
                - 3 * A
            )

            actual = s.n + 9

        else:

            expected = (
                A * B
                + 3 * A
                - B
            )

            actual = s.n + 3

        if actual != expected:

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch:"
                    f" n={s.n}"
                    f" frame={frame}"
                    f" A={A}"
                    f" B={B}"
                    f" actual={actual}"
                    f" expected={expected}"
                )

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 1
# ==============================================================================

def test_d_divides_shift(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 1: gcd(A,B) DIVIDES n+c")
    print("=" * 90)

    failures = 0

    for s in states:

        frame = frame_of(s.n)

        A, B = residuals(
            frame,
            s.p,
            s.q,
        )

        d = gcd(
            abs(A),
            abs(B),
        )

        if frame == "A":

            shifted = s.n + 9

        else:

            shifted = s.n + 3

        if d == 0 or shifted % d != 0:

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch:"
                    f" n={s.n}"
                    f" frame={frame}"
                    f" d={d}"
                    f" shifted={shifted}"
                )

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 2
# ==============================================================================

def test_quotient_formula(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 2: EXACT QUOTIENT FORMULA")
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

        d = gcd(
            abs(A),
            abs(B),
        )

        if d == 0:
            continue

        a = A // d
        b = B // d

        if frame == "A":

            shifted = s.n + 9

            predicted = (
                d * a * b
                + 3 * (b - a)
            )

        else:

            shifted = s.n + 3

            predicted = (
                d * a * b
                + 3 * a
                - b
            )

        quotient = shifted // d

        checked += 1

        if quotient != predicted:

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch:"
                    f" n={s.n}"
                    f" frame={frame}"
                    f" A={A}"
                    f" B={B}"
                    f" d={d}"
                    f" a={a}"
                    f" b={b}"
                    f" quotient={quotient}"
                    f" predicted={predicted}"
                )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 3
# ==============================================================================

def test_quotient_parity(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 3: QUOTIENT PARITY == v2 EQUALITY")
    print("=" * 90)

    failures = 0
    checked = 0

    distribution = Counter()

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

        a = A // d
        b = B // d

        if frame == "A":

            shifted = s.n + 9

        else:

            shifted = s.n + 3

        quotient = shifted // d

        equal_v2 = (
            v2(A) == v2(B)
        )

        quotient_even = (
            quotient & 1
        ) == 0

        checked += 1

        distribution[
            (
                frame,
                equal_v2,
                quotient_even,
            )
        ] += 1

        if equal_v2 != quotient_even:

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch:"
                    f" n={s.n}"
                    f" frame={frame}"
                    f" A={A}"
                    f" B={B}"
                    f" quotient={quotient}"
                    f" equal_v2={equal_v2}"
                )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )

    print()
    print("    distribution:")

    for key in sorted(
        distribution
    ):

        print(
            f"        {key}: "
            f"{distribution[key]}"
        )


# ==============================================================================
# TEST 4
# ==============================================================================

def test_exact_shift_valuation(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 4: v2(n+c) = v2(d) + v2(QUOTIENT)")
    print("=" * 90)

    failures = 0
    checked = 0

    distribution = Counter()

    for s in states:

        frame = frame_of(s.n)

        A, B = residuals(
            frame,
            s.p,
            s.q,
        )

        d = gcd(
            abs(A),
            abs(B),
        )

        if d == 0:
            continue

        if frame == "A":

            shifted = s.n + 9

        else:

            shifted = s.n + 3

        if shifted == 0:

            shift_v = INF

        else:

            shift_v = v2(
                abs(shifted)
            )

        quotient = shifted // d

        d_v = v2(d)
        q_v = v2(quotient)

        checked += 1

        if shift_v != d_v + q_v:

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch:"
                    f" n={s.n}"
                    f" frame={frame}"
                    f" d={d}"
                    f" quotient={quotient}"
                    f" v2shift={shift_v}"
                    f" v2d={d_v}"
                    f" v2q={q_v}"
                )

        distribution[
            (
                frame,
                d_v,
                q_v,
                shift_v,
            )
        ] += 1

    print(
        f"checked={checked} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 5
# ==============================================================================

def test_depth_from_quotient(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 5: DEPTH FROM d + QUOTIENT PARITY")
    print("=" * 90)

    failures = 0

    for s in states:

        frame = frame_of(s.n)

        A, B = residuals(
            frame,
            s.p,
            s.q,
        )

        actual = global_depth(
            frame,
            s.p,
            s.q,
        )

        # Degenerate A=0 branch.
        if A == 0:

            d = abs(B)

            predicted = (
                1 + v2(d // 2)
            )

        # Degenerate B=0 branch.
        elif B == 0:

            d = abs(A)

            predicted = (
                1 + v2(
                    (3 * d) // 2
                )
            )

        else:

            d = gcd(
                abs(A),
                abs(B),
            )

            if frame == "A":

                shifted = s.n + 9

            else:

                shifted = s.n + 3

            quotient = shifted // d

            # The quotient parity is precisely the equality flag.
            #
            # quotient odd:
            #     unequal v2
            #
            # quotient even:
            #     equal v2

            predicted = (
                v2(d)
                + 1
                + (
                    1
                    if (
                        quotient & 1
                    ) == 0
                    else 0
                )
                - 1
            )

            # Simplifies to:
            #
            #   v2(d)+1, if quotient even
            #   v2(d),   if quotient odd.
            #
            # Write it explicitly for readability.

            if quotient & 1:

                predicted = v2(d)

            else:

                predicted = (
                    v2(d) + 1
                )

        if predicted != actual:

            failures += 1

            if failures <= 20:

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
        f"checked={len(states)} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 6
# ==============================================================================

def test_odd_part_relation(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 6: ODD gcd PART VS n+c QUOTIENT")
    print("=" * 90)

    # We do not expect the quotient to uniquely determine
    # the odd part of d. Instead we test exact divisibility
    # identities.
    #
    # Since:
    #
    #     n+c = d * Q
    #
    # every odd prime of d is also an odd prime divisor
    # of n+c.
    #
    # We record the gcd between odd(d) and odd(Q) to see
    # how much odd overlap survives.

    failures = 0
    checked = 0

    overlap = Counter()

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

        if frame == "A":

            shifted = s.n + 9

        else:

            shifted = s.n + 3

        quotient = shifted // d

        od = odd_part(d)
        oq = odd_part(quotient)

        g = gcd(
            od,
            oq,
        )

        checked += 1

        overlap[
            (
                frame,
                g == 1,
                v2(d),
            )
        ] += 1

        # Fundamental divisibility check.
        if shifted != d * quotient:

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch:"
                    f" n={s.n}"
                    f" shifted={shifted}"
                    f" d={d}"
                    f" quotient={quotient}"
                )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )

    print()
    print(
        "    odd-part overlap categories="
    )

    for key in sorted(
        overlap
    ):

        print(
            f"        {key}: "
            f"{overlap[key]}"
        )


# ==============================================================================
# TEST 7
# ==============================================================================

def test_primary_valuation_statement(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 7: PRIMARY VALUATION THEOREM")
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

        if frame == "A":

            shifted = s.n + 9

        else:

            shifted = s.n + 3

        m = v2(d)
        primary = v2(shifted)

        equal = (
            v2(A) == v2(B)
        )

        if equal:

            expected_relation = (
                primary > m
            )

        else:

            expected_relation = (
                primary == m
            )

        checked += 1

        if not expected_relation:

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch:"
                    f" n={s.n}"
                    f" frame={frame}"
                    f" d={d}"
                    f" m={m}"
                    f" primary={primary}"
                    f" equal={equal}"
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

    target_n = [
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
        213,
    ]

    by_n = {
        s.n: s
        for s in states
    }

    for n in target_n:

        s = by_n[n]

        frame = frame_of(n)

        A, B = residuals(
            frame,
            s.p,
            s.q,
        )

        d = gcd(
            abs(A),
            abs(B),
        )

        if frame == "A":

            c = 9

        else:

            c = 3

        shifted = n + c

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
            f" d={d}"
        )

        print(
            f"    n+c={shifted}"
        )

        if d:

            quotient = shifted // d

            print(
                f"    quotient=(n+c)/d="
                f"{quotient}"
            )

            print(
                f"    v2(d)={v2(d)}"
                f" v2(quotient)={v2(quotient)}"
            )

            print(
                f"    v2(n+c)={v2(shifted)}"
            )

        if A != 0 and B != 0:

            print(
                f"    v2(A)={v2(A)}"
                f" v2(B)={v2(B)}"
            )

            print(
                f"    equal-v2="
                f"{v2(A) == v2(B)}"
            )

        print(
            f"    depth="
            f"{global_depth(frame,s.p,s.q)}"
        )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print("EXPERIMENT 634 START")
    print("=" * 90)

    print()
    print(
        "QUOTIENT-LADDER / gcd-RESIDUAL DECOMPOSITION"
    )

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

    test_exact_n_shift_identity(
        states
    )

    test_d_divides_shift(
        states
    )

    test_quotient_formula(
        states
    )

    test_quotient_parity(
        states
    )

    test_exact_shift_valuation(
        states
    )

    test_depth_from_quotient(
        states
    )

    test_odd_part_relation(
        states
    )

    test_primary_valuation_statement(
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
For the raw factor residuals:

    FRAME A:

        A = p-3
        B = q+3

        n+9
          =
        AB + 3B - 3A

          =
        d * (
            d*a*b
            + 3(b-a)
        )

    where:

        d=gcd(A,B)
        A=d*a
        B=d*b
        gcd(a,b)=1.

    FRAME B:

        A = p+1
        B = q-3

        n+3
          =
        AB + 3A - B

          =
        d * (
            d*a*b
            + 3a-b
        ).

Therefore:

    d | n+9       FRAME A

    d | n+3       FRAME B.

The quotient is the exact carrier of the extra
2-adic layer.

For A,B != 0:

    unequal v2(A),v2(B)
        ->
    quotient is odd

        ->
    v2(n+c)=v2(d)

        ->
    depth=v2(d).

For equal v2(A),v2(B):

    quotient is even

        ->
    v2(n+c)>v2(d)

        ->
    depth=v2(d)+1.

Thus:

    depth =
        v2(d)
            if quotient odd,

        v2(d)+1
            if quotient even,

    where:

        d=gcd(A,B).

This is equivalent to:

    depth =
        v2(gcd(A,B))
        + [v2(A)=v2(B)].

The important new integer identity is therefore:

    n+c = gcd(A,B) * Q

where the parity of Q itself is exactly the equality
branch.

The remaining information hidden from n-only searches is
not the equality flag anymore: it is the 2-adic valuation
of the divisor

    gcd(A,B).

That is the remaining factor-dependent quantity.
"""
    )

    print()
    print("=" * 90)
    print("EXPERIMENT 634 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
