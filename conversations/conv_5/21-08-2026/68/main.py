#!/usr/bin/env python3

from __future__ import annotations

from dataclasses import dataclass
from math import gcd, isqrt
from collections import Counter


# ==============================================================================
# EXPERIMENT 630
# ==============================================================================
#
# EXACT INTEGER GCD BRIDGE BETWEEN FACTOR RESIDUALS AND GLOBAL (X,Y)
#
# Experiment 629 established:
#
#     depth = 1 + min(v2(X), v2(Y))
#
# exactly.
#
# We now test the stronger integer identity:
#
#     d = gcd(A,B)
#
#     gcd(X,Y) =
#
#         d       if v2(A) == v2(B)
#         d / 2   if v2(A) != v2(B)
#
# for the two factor frames.
#
# This is stronger than checking only v2().
#
# If exact, it gives:
#
#     depth
#       = 1 + v2(gcd(X,Y))
#
#       = v2(gcd(A,B))
#         + [v2(A)=v2(B)].
#
# No n-only search is performed.
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

        A = p - 3
        B = q + 3

    else:

        A = p + 1
        B = q - 3

    return A, B


# ==============================================================================
# GLOBAL X,Y
# ==============================================================================

def global_xy(
    frame: str,
    p: int,
    q: int,
) -> tuple[int, int]:

    if frame == "A":

        X = (q - p + 6) // 2
        Y = (p + q) // 2

    else:

        X = (3 * p - q + 6) // 2
        Y = (3 * p + q) // 2

    return X, Y


# ==============================================================================
# DEPTH
# ==============================================================================

def observed_depth(
    X: int,
    Y: int,
) -> int:

    return (
        1 + min(
            v2(X),
            v2(Y),
        )
    )


# ==============================================================================
# TEST 0
# ==============================================================================

def test_domain(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 0: ODD SEMIPRIME DOMAIN")
    print("=" * 90)

    failures = 0

    for state in states:

        if state.n & 1 == 0:

            failures += 1

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 1
# ==============================================================================

def test_global_gcd_depth(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 1: DEPTH = 1 + v2(gcd(X,Y))")
    print("=" * 90)

    failures = 0

    for state in states:

        frame = frame_of(state.n)

        X, Y = global_xy(
            frame,
            state.p,
            state.q,
        )

        gxy = gcd(
            abs(X),
            abs(Y),
        )

        predicted = (
            1 + v2(gxy)
        )

        actual = observed_depth(
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
                    f" X={X}"
                    f" Y={Y}"
                    f" gcd={gxy}"
                    f" predicted={predicted}"
                    f" actual={actual}"
                )

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 2
# ==============================================================================

def test_exact_gcd_identity(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 2: EXACT INTEGER GCD IDENTITY")
    print("=" * 90)

    failures = 0

    equal_count = 0
    unequal_count = 0

    for state in states:

        frame = frame_of(state.n)

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

        gxy = gcd(
            abs(X),
            abs(Y),
        )

        alpha = v2(A)
        beta = v2(B)

        if alpha == beta:

            expected = d
            equal_count += 1

        else:

            expected = d // 2
            unequal_count += 1

        if gxy != expected:

            failures += 1

            if failures <= 30:

                print(
                    "    mismatch"
                    f" n={state.n}"
                    f" frame={frame}"
                    f" A={A}"
                    f" B={B}"
                    f" gcdAB={d}"
                    f" gcdXY={gxy}"
                    f" v2A={alpha}"
                    f" v2B={beta}"
                    f" expected={expected}"
                )

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )

    print(
        f"    equal valuation states="
        f"{equal_count}"
    )

    print(
        f"    unequal valuation states="
        f"{unequal_count}"
    )


# ==============================================================================
# TEST 3
# ==============================================================================

def test_gcd_valuation_bridge(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 3: v2(gcdAB) -> v2(gcdXY)")
    print("=" * 90)

    failures = 0

    for state in states:

        frame = frame_of(state.n)

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

        gxy = gcd(
            abs(X),
            abs(Y),
        )

        alpha = v2(A)
        beta = v2(B)

        predicted = (
            v2(d)
            - int(alpha != beta)
        )

        actual = v2(gxy)

        if predicted != actual:

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch"
                    f" n={state.n}"
                    f" frame={frame}"
                    f" v2A={alpha}"
                    f" v2B={beta}"
                    f" v2gcdAB={v2(d)}"
                    f" v2gcdXY={actual}"
                    f" predicted={predicted}"
                )

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 4
# ==============================================================================

def test_depth_raw_formula(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 4: DEPTH FROM RAW GCD + EQUALITY")
    print("=" * 90)

    failures = 0

    for state in states:

        frame = frame_of(state.n)

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

        predicted = (
            v2(d)
            + int(alpha == beta)
        )

        actual = observed_depth(
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
                    f" gcdAB={d}"
                    f" predicted={predicted}"
                    f" actual={actual}"
                )

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 5
# ==============================================================================

def test_same_gcd_odd_part(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 5: ODD PART OF GCD(X,Y) == ODD PART OF GCD(A,B)")
    print("=" * 90)

    failures = 0

    for state in states:

        frame = frame_of(state.n)

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

        gxy = gcd(
            abs(X),
            abs(Y),
        )

        vd = v2(d)
        vg = v2(gxy)

        odd_d = (
            d >> vd
            if d != 0
            else 0
        )

        odd_g = (
            gxy >> vg
            if gxy != 0
            else 0
        )

        if odd_d != odd_g:

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch"
                    f" n={state.n}"
                    f" frame={frame}"
                    f" gcdAB={d}"
                    f" gcdXY={gxy}"
                    f" oddAB={odd_d}"
                    f" oddXY={odd_g}"
                )

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 6
# ==============================================================================

def test_full_gcd_ratio_distribution(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 6: gcd(A,B) / gcd(X,Y) RATIO")
    print("=" * 90)

    counts = Counter()

    unexpected = 0

    for state in states:

        frame = frame_of(state.n)

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

        gxy = gcd(
            abs(X),
            abs(Y),
        )

        if gxy == 0:

            ratio = "undefined"

        else:

            if d % gxy != 0:

                ratio = "noninteger"

            else:

                ratio = d // gxy

        counts[ratio] += 1

        if ratio not in (1, 2):

            unexpected += 1

            if unexpected <= 20:

                print(
                    "    unexpected ratio"
                    f" n={state.n}"
                    f" frame={frame}"
                    f" gcdAB={d}"
                    f" gcdXY={gxy}"
                    f" ratio={ratio}"
                )

    print(
        "    ratio distribution:"
    )

    for key in sorted(
        counts,
        key=lambda x: str(x),
    ):

        print(
            f"        {key}: "
            f"{counts[key]}"
        )

    print(
        f"    unexpected={unexpected}"
    )


# ==============================================================================
# TEST 7
# ==============================================================================

def test_frame_specific_linear_map(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 7: SYMBOLIC LINEAR-MAP GCD CHECK")
    print("=" * 90)

    failures = 0

    for state in states:

        frame = frame_of(state.n)

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

            Xa = (
                B - A
            ) // 2

            Ya = (
                A + B
            ) // 2

        else:

            # For frame B:
            #
            # X = (3A-B)/2
            # Y = (3A+B)/2

            Xa = (
                3 * A - B
            ) // 2

            Ya = (
                3 * A + B
            ) // 2

        if (
            Xa != X
            or Ya != Y
        ):

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch"
                    f" n={state.n}"
                    f" frame={frame}"
                    f" XY=({X},{Y})"
                    f" map=({Xa},{Ya})"
                )

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )


# ==============================================================================
# EXAMPLES
# ==============================================================================

def print_examples(
    states: list[State],
) -> None:

    examples = [
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

    for n in examples:

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

        gxy = gcd(
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
            f"    gcd(X,Y)={gxy}"
        )

        print(
            f"    v2(A)={alpha}"
            f" v2(B)={beta}"
        )

        print(
            f"    equal={alpha == beta}"
        )

        print(
            f"    ratio="
            f"{d // gxy}"
        )

        print(
            f"    depth="
            f"{1 + v2(gxy)}"
        )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print("EXPERIMENT 630 START")
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

    test_domain(
        states
    )

    test_global_gcd_depth(
        states
    )

    test_exact_gcd_identity(
        states
    )

    test_gcd_valuation_bridge(
        states
    )

    test_depth_raw_formula(
        states
    )

    test_same_gcd_odd_part(
        states
    )

    test_full_gcd_ratio_distribution(
        states
    )

    test_frame_specific_linear_map(
        states
    )

    print_examples(
        states
    )

    print()
    print("=" * 90)
    print("FINAL SUMMARY")
    print("=" * 90)

    print()
    print(
        "Target identities:"
    )

    print()
    print(
        "    gcd(X,Y)"
        " ="
    )

    print(
        "        gcd(A,B)"
        "      if v2(A)=v2(B)"
    )

    print(
        "        gcd(A,B)/2"
        "  if v2(A)!=v2(B)"
    )

    print()
    print(
        "and therefore:"
    )

    print()
    print(
        "    depth"
        " = 1 + v2(gcd(X,Y))"
    )

    print()
    print(
        "    depth"
        " = v2(gcd(A,B))"
        " + [v2(A)=v2(B)]"
    )

    print()
    print(
        "If all tests pass, the equality branch is not"
        " an arbitrary extra state."
    )

    print(
        "It is exactly the choice between:"
    )

    print(
        "    gcd(X,Y) = gcd(A,B)"
    )

    print(
        "and:"
    )

    print(
        "    gcd(X,Y) = gcd(A,B)/2."
    )

    print()
    print(
        "This gives an integer-level explanation for the"
        " 2-adic depth collapse."
    )

    print()
    print("=" * 90)
    print("EXPERIMENT 630 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
