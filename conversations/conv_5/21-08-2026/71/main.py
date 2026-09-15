#!/usr/bin/env python3

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import gcd, isqrt


# ==============================================================================
# EXPERIMENT 633
# ==============================================================================
#
# EXACT COPRIME-RESIDUAL GCD THEOREM
#
# Experiment 632 established the final depth law empirically, but Test 3
# contained an incorrect expectation for the normalized transformed gcd.
#
# This experiment removes the valuation language and derives the result
# directly from:
#
#       A = d*a
#       B = d*b
#       gcd(a,b) = 1.
#
# FRAME A:
#
#       X = (B-A)/2
#       Y = (A+B)/2
#
# FRAME B:
#
#       X = (3A-B)/2
#       Y = (3A+B)/2
#
# The target theorem is:
#
#       gcd(X,Y) = d
#           if a,b are both odd
#
#       gcd(X,Y) = d/2
#           if a,b have opposite parity.
#
# Because gcd(a,b)=1, the only possible parity classes are:
#
#       (odd,odd)
#       (odd,even)
#       (even,odd).
#
# (even,even) is impossible.
#
# The Frame-B factor 3 must also be shown to be absent from the
# nondegenerate branch B != 0.
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

    sieve = bytearray(b"\x01") * (limit + 1)

    sieve[0] = 0
    sieve[1] = 0

    for p in range(2, isqrt(limit) + 1):

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
            2
        )
        if sieve[p]
    ]


# ==============================================================================
# SEMIPRIME GENERATION
# ==============================================================================

def generate_states(
    primes: list[int],
) -> list[State]:

    result: list[State] = []

    for i, p in enumerate(primes):

        for q in primes[i:]:

            result.append(
                State(
                    n=p * q,
                    p=p,
                    q=q,
                )
            )

    return result


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
        f"odd semiprime required: n={n}, n mod 4={r}"
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
# ODD PART
# ==============================================================================

def odd_part(x: int) -> int:

    if x == 0:
        return 0

    return abs(x) >> v2(x)


# ==============================================================================
# TEST 0
# ==============================================================================

def test_coprime_normalization(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 0: COPRIME RESIDUAL NORMALIZATION")
    print("=" * 90)

    failures = 0
    checked = 0

    for s in states:

        frame = frame_of(s.n)

        A, B = raw_residuals(
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

        checked += 1

        if gcd(
            abs(a),
            abs(b),
        ) != 1:

            failures += 1

            if failures <= 20:
                print(
                    "    gcd failure:"
                    f" n={s.n}"
                    f" A={A}"
                    f" B={B}"
                    f" d={d}"
                    f" a={a}"
                    f" b={b}"
                )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 1
# ==============================================================================

def test_parity_classes(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 1: COPRIME PARITY CLASSIFICATION")
    print("=" * 90)

    counts = Counter()
    failures = 0
    checked = 0

    for s in states:

        frame = frame_of(s.n)

        A, B = raw_residuals(
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

        checked += 1

        key = (
            a & 1,
            b & 1,
        )

        counts[key] += 1

        # Coprimality makes (0,0) impossible.
        if key == (0, 0):

            failures += 1

            if failures <= 20:
                print(
                    "    impossible even/even pair:"
                    f" n={s.n}"
                    f" a={a}"
                    f" b={b}"
                )

    print(
        f"checked={checked}"
    )

    for key in sorted(counts):

        print(
            f"    parity={key}: "
            f"{counts[key]}"
        )

    print(
        f"failures={failures}"
    )


# ==============================================================================
# TEST 2
# ==============================================================================

def test_frame_A_theorem(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 2: FRAME-A COPRIME gcd THEOREM")
    print("=" * 90)

    failures = 0
    checked = 0

    for s in states:

        if frame_of(s.n) != "A":
            continue

        A, B = raw_residuals(
            "A",
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

        X, Y = global_xy(
            "A",
            s.p,
            s.q,
        )

        actual = gcd(
            abs(X),
            abs(Y),
        )

        if (a & 1) and (b & 1):

            expected = d

        else:

            expected = d // 2

        checked += 1

        if actual != expected:

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch:"
                    f" n={s.n}"
                    f" A={A}"
                    f" B={B}"
                    f" a={a}"
                    f" b={b}"
                    f" d={d}"
                    f" actual={actual}"
                    f" expected={expected}"
                )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 3
# ==============================================================================

def test_frame_B_theorem(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 3: FRAME-B COPRIME gcd THEOREM")
    print("=" * 90)

    failures = 0
    checked = 0

    for s in states:

        if frame_of(s.n) != "B":
            continue

        A, B = raw_residuals(
            "B",
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

        X, Y = global_xy(
            "B",
            s.p,
            s.q,
        )

        actual = gcd(
            abs(X),
            abs(Y),
        )

        if (a & 1) and (b & 1):

            expected = d

        else:

            expected = d // 2

        checked += 1

        if actual != expected:

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch:"
                    f" n={s.n}"
                    f" A={A}"
                    f" B={B}"
                    f" a={a}"
                    f" b={b}"
                    f" d={d}"
                    f" actual={actual}"
                    f" expected={expected}"
                )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 4
# ==============================================================================

def test_normalized_transformed_gcd(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 4: NORMALIZED TRANSFORMED gcd")
    print("=" * 90)

    failures = 0
    checked = 0

    distribution = Counter()

    for s in states:

        frame = frame_of(s.n)

        A, B = raw_residuals(
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

            u = b - a
            v = a + b

        else:

            u = 3 * a - b
            v = 3 * a + b

        h = gcd(
            abs(u),
            abs(v),
        )

        distribution[
            (
                frame,
                a & 1,
                b & 1,
                h,
            )
        ] += 1

        checked += 1

        if (a & 1) and (b & 1):

            expected_h = 2

        else:

            expected_h = 1

        if h != expected_h:

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch:"
                    f" n={s.n}"
                    f" frame={frame}"
                    f" a={a}"
                    f" b={b}"
                    f" h={h}"
                    f" expected={expected_h}"
                )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )

    print()
    print("    normalized gcd distribution:")

    for key in sorted(
        distribution
    ):

        print(
            f"        {key}: "
            f"{distribution[key]}"
        )


# ==============================================================================
# TEST 5
# ==============================================================================

def test_factor_three_absence(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 5: FRAME-B FACTOR-3 ABSENCE")
    print("=" * 90)

    failures = 0
    checked = 0
    three_states = []

    for s in states:

        if frame_of(s.n) != "B":
            continue

        A, B = raw_residuals(
            "B",
            s.p,
            s.q,
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

        h = gcd(
            abs(3 * a - b),
            abs(3 * a + b),
        )

        if h % 3 == 0:

            three_states.append(
                (
                    s.n,
                    s.p,
                    s.q,
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
        f"factor-3 states="
        f"{len(three_states)}"
    )

    for rec in three_states[:20]:

        print(
            "    unexpected:"
            f" n={rec[0]}"
            f" p={rec[1]}"
            f" q={rec[2]}"
            f" A={rec[3]}"
            f" B={rec[4]}"
            f" a={rec[5]}"
            f" b={rec[6]}"
            f" h={rec[7]}"
        )

    if three_states:
        failures = len(
            three_states
        )

    print(
        f"failures={failures}"
    )

    print()
    print(
        "    Algebraic reason:"
    )
    print(
        "        3 | (3a-b) and (3a+b)"
    )
    print(
        "        => 3 | 2b"
    )
    print(
        "        => 3 | b"
    )
    print(
        "        but B=d*b=q-3."
    )
    print(
        "        For B!=0, q!=3."
    )
    print(
        "        Since q is prime, 3 cannot divide q-3."
    )
    print(
        "        Therefore 3 cannot divide b."
    )


# ==============================================================================
# TEST 6
# ==============================================================================

def test_depth_from_coprime_parity(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 6: DEPTH FROM COPRIME RESIDUAL PARITY")
    print("=" * 90)

    failures = 0

    for s in states:

        frame = frame_of(s.n)

        A, B = raw_residuals(
            frame,
            s.p,
            s.q,
        )

        X, Y = global_xy(
            frame,
            s.p,
            s.q,
        )

        actual = (
            1 + v2(
                gcd(
                    abs(X),
                    abs(Y),
                )
            )
        )

        # Degenerate branches first.

        if A == 0:

            predicted = v2(
                abs(B) // 2
            ) + 1

        elif B == 0:

            predicted = v2(
                3 * abs(A) // 2
            ) + 1

        else:

            d = gcd(
                abs(A),
                abs(B),
            )

            a = A // d
            b = B // d

            if (a & 1) and (b & 1):

                predicted = (
                    v2(d) + 1
                )

            else:

                predicted = (
                    v2(d)
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
# TEST 7
# ==============================================================================

def test_valuation_equivalence(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 7: PARITY <-> v2 EQUALITY")
    print("=" * 90)

    failures = 0
    checked = 0

    for s in states:

        A, B = raw_residuals(
            frame_of(s.n),
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

        parity_equal = (
            (a & 1) == (b & 1)
        )

        # Since gcd(a,b)=1, equality of parity
        # can only mean both are odd.

        parity_equal = (
            (a & 1) == 1
            and
            (b & 1) == 1
        )

        valuation_equal = (
            v2(A) == v2(B)
        )

        checked += 1

        if parity_equal != valuation_equal:

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch:"
                    f" n={s.n}"
                    f" A={A}"
                    f" B={B}"
                    f" a={a}"
                    f" b={b}"
                    f" parity_equal={parity_equal}"
                    f" valuation_equal={valuation_equal}"
                )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 8
# ==============================================================================

def test_exact_integer_formula(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 8: EXACT INTEGER gcd FORMULA")
    print("=" * 90)

    failures = 0

    for s in states:

        frame = frame_of(s.n)

        A, B = raw_residuals(
            frame,
            s.p,
            s.q,
        )

        X, Y = global_xy(
            frame,
            s.p,
            s.q,
        )

        actual = gcd(
            abs(X),
            abs(Y),
        )

        if A == 0:

            predicted = (
                abs(B) // 2
            )

        elif B == 0:

            predicted = (
                3 * abs(A) // 2
            )

        else:

            d = gcd(
                abs(A),
                abs(B),
            )

            if v2(A) == v2(B):

                predicted = d

            else:

                predicted = d // 2

        if actual != predicted:

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch:"
                    f" n={s.n}"
                    f" frame={frame}"
                    f" A={A}"
                    f" B={B}"
                    f" actual={actual}"
                    f" predicted={predicted}"
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
        57,
        69,
        77,
        93,
        141,
        213,
    ]

    lookup = {
        s.n: s
        for s in states
    }

    for n in targets:

        s = lookup.get(n)

        if s is None:
            continue

        frame = frame_of(n)

        A, B = raw_residuals(
            frame,
            s.p,
            s.q,
        )

        X, Y = global_xy(
            frame,
            s.p,
            s.q,
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
        )

        print(
            f"    X={X}"
            f" Y={Y}"
        )

        if A != 0 and B != 0:

            d = gcd(
                abs(A),
                abs(B),
            )

            a = A // d
            b = B // d

            print(
                f"    d=gcd(A,B)={d}"
            )

            print(
                f"    normalized=(a,b)=({a},{b})"
            )

            print(
                f"    normalized parity="
                f"({a & 1},{b & 1})"
            )

            print(
                f"    equal-v2="
                f"{v2(A) == v2(B)}"
            )

        else:

            print(
                "    degenerate residual branch"
            )

        print(
            f"    gcdXY="
            f"{gcd(abs(X),abs(Y))}"
        )

        print(
            f"    depth="
            f"{1 + v2(gcd(abs(X),abs(Y)))}"
        )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print("EXPERIMENT 633 START")
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
    print("[3] BUILD GLOBAL STATES")

    print(
        f"    states={len(states)}"
    )

    test_coprime_normalization(
        states
    )

    test_parity_classes(
        states
    )

    test_frame_A_theorem(
        states
    )

    test_frame_B_theorem(
        states
    )

    test_normalized_transformed_gcd(
        states
    )

    test_factor_three_absence(
        states
    )

    test_depth_from_coprime_parity(
        states
    )

    test_valuation_equivalence(
        states
    )

    test_exact_integer_formula(
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
Let:

    d = gcd(A,B)

and for A,B != 0:

    A=d*a
    B=d*b
    gcd(a,b)=1.

Therefore the normalized pair has exactly three
possible parity classes:

    (a,b)=(odd,odd)
    (a,b)=(odd,even)
    (a,b)=(even,odd).

No (even,even) state is possible.

FRAME A:

    X = d(b-a)/2
    Y = d(a+b)/2.

Hence:

    gcd(X,Y)
      =
    d
      when a,b are both odd,

    d/2
      when a,b have opposite parity.

FRAME B:

    X = d(3a-b)/2
    Y = d(3a+b)/2.

For B != 0, the possible factor 3 in the transformed gcd
is impossible because:

    3 | (3a-b)
    3 | (3a+b)

    => 3 | b

but:

    B=d*b=q-3

and q != 3 in the nondegenerate branch, so 3 ∤ B.

Therefore again:

    gcd(X,Y)
      =
    d
      when a,b are both odd,

    d/2
      when a,b have opposite parity.

Because gcd(a,b)=1:

    a,b both odd
        <=> v2(A)=v2(B).

Thus:

    gcd(X,Y)
      =
    gcd(A,B)
        if v2(A)=v2(B),

    gcd(A,B)/2
        otherwise.

Finally:

    depth
      =
    1 + v2(gcd(X,Y)).

Therefore, for A,B != 0:

    depth =
        v2(gcd(A,B))
        + 1,
        if v2(A)=v2(B),

        v2(gcd(A,B)),
        otherwise.

The remaining degenerate branches are:

    Frame A, A=0:
        gcd(X,Y)=|B|/2.

    Frame B, B=0:
        gcd(X,Y)=3|A|/2.

The observed B=0 branch is uniquely:

    p=q=3,
    n=9.

This gives an integer-gcd formulation of the complete
2-adic depth mechanism without any level-by-level traversal.
"""
    )

    print()
    print("=" * 90)
    print("EXPERIMENT 633 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
