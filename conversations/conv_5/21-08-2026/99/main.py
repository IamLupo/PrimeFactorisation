#!/usr/bin/env python3

"""
==========================================================================================
EXPERIMENT 662

EXACT INTEGER gcd FACTORIZATION THROUGH d = gcd(A,B)

GOAL
----

Experiment 661 established empirically:

    generic:
        gcd(2A,n+c) = gcd(2B,n+c)
                      = gcd(2A,2B,n+c)

    A=0:
        gcd(2A,n+9) = 3*gcd(2B,n+9)

    B=0:
        gcd(2B,n+3) = 3*gcd(2A,n+3)

and:

    depth = v2(gcd(2A,2B,n+c)).

This experiment tests the stronger exact integer decomposition.

Define:

    d = gcd(A,B)

and, since d | (n+c),

    Q = (n+c)/d.

Write:

    A = d*a
    B = d*b
    gcd(a,b)=1.

Then:

    gcd(2A,2B,n+c)
      =
    gcd(2da,2db,dQ)

      =
    d * gcd(2a,2b,Q).

Because gcd(a,b)=1:

    gcd(2a,2b) = 2.

Therefore the expected exact theorem is:

    G = gcd(2A,2B,n+c)
      = d * gcd(2,Q).

Hence:

    Q odd:
        G = d

    Q even:
        G = 2d.

The experiment also tests:

    Q even <=> v2(A)=v2(B)

and therefore:

    G = d * 2^[v2(A)=v2(B)].

This should give a complete integer decomposition of
the symmetric gcd.

The one-sided gcd relation from Experiment 661 is then
tested separately:

    generic A!=0,B!=0:
        gcd(2A,n+c) = G = gcd(2B,n+c)

    Frame A, A=0:
        gcd(2A,n+9) = 3G
        gcd(2B,n+9) = G

    Frame B, B=0:
        gcd(2A,n+3) = G
        gcd(2B,n+3) = 3G.

==========================================================================================
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import gcd, isqrt


INF = 10**9


# =============================================================================
# BASIC NUMBER THEORY
# =============================================================================

def v2(x: int) -> int:
    x = abs(x)

    if x == 0:
        return INF

    return (x & -x).bit_length() - 1


def odd_part(x: int) -> int:
    x = abs(x)

    if x == 0:
        return 0

    return x >> v2(x)


def sieve_primes(limit: int) -> list[int]:

    sieve = bytearray(b"\x01") * (limit + 1)

    if limit >= 0:
        sieve[0] = 0

    if limit >= 1:
        sieve[1] = 0

    for p in range(2, isqrt(limit) + 1):

        if sieve[p]:

            start = p * p

            sieve[start:limit + 1:p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    return [
        p
        for p in range(3, limit + 1, 2)
        if sieve[p]
    ]


# =============================================================================
# STATE
# =============================================================================

@dataclass(frozen=True)
class State:

    n: int
    p: int
    q: int

    frame: str

    A: int
    B: int

    c: int

    depth: int


def make_state(p: int, q: int) -> State:

    n = p * q

    if n % 4 == 3:

        # FRAME A
        #
        # A = p - 3
        # B = q + 3
        # n + 9 = AB + 3(B-A)

        frame = "A"

        A = p - 3
        B = q + 3
        c = 9

    else:

        # FRAME B
        #
        # A = p + 1
        # B = q - 3
        # n + 3 = AB + 3A - B

        frame = "B"

        A = p + 1
        B = q - 3
        c = 3

    depth = min(
        v2(A) + 1,
        v2(n + c),
    )

    return State(
        n=n,
        p=p,
        q=q,
        frame=frame,
        A=A,
        B=B,
        c=c,
        depth=depth,
    )


def build_states(limit: int) -> list[State]:

    primes = sieve_primes(limit)

    states: list[State] = []

    for i, p in enumerate(primes):

        for q in primes[i:]:

            states.append(
                make_state(p, q)
            )

    return states


# =============================================================================
# TEST 0
# =============================================================================

def test_baseline(states: list[State]) -> int:

    print("=" * 90)
    print("TEST 0: BASELINE DEPTH LAW")
    print("=" * 90)

    failures = 0

    for s in states:

        predicted = min(
            v2(s.A) + 1,
            v2(s.n + s.c),
        )

        if predicted != s.depth:
            failures += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 1
# =============================================================================

def test_d_divides_shift(states: list[State]) -> int:

    print("=" * 90)
    print("TEST 1: d = gcd(A,B) DIVIDES n+c")
    print("=" * 90)

    failures = 0

    for s in states:

        d = gcd(
            abs(s.A),
            abs(s.B),
        )

        shift = abs(s.n + s.c)

        if d == 0:

            # Only possible if A=B=0.
            # This does not occur in the tested domain.

            failures += 1
            continue

        if shift % d != 0:
            failures += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 2
# =============================================================================

def test_exact_quotient(states: list[State]) -> int:

    print("=" * 90)
    print("TEST 2: EXACT QUOTIENT Q = (n+c)/d")
    print("=" * 90)

    failures = 0

    quotient_parity = Counter()

    examples = []

    for s in states:

        d = gcd(
            abs(s.A),
            abs(s.B),
        )

        shift = abs(s.n + s.c)

        if d == 0 or shift % d != 0:

            failures += 1
            continue

        Q = shift // d

        quotient_parity[
            (s.frame, Q % 2)
        ] += 1

        if Q <= 0:

            failures += 1

            if len(examples) < 20:
                examples.append(
                    (s, d, Q)
                )

    print(
        "quotient parity distribution:"
    )

    for key, count in sorted(
        quotient_parity.items()
    ):
        print(
            f"    {key}: {count}"
        )

    if examples:

        print()
        print("bad quotient examples:")

        for s, d, Q in examples:

            print(
                f"    n={s.n} "
                f"p={s.p} "
                f"q={s.q} "
                f"frame={s.frame} "
                f"d={d} "
                f"Q={Q}"
            )

    print()
    print(
        f"checked={len(states)}"
    )

    print(
        f"failures={failures}"
    )

    print()

    return failures


# =============================================================================
# TEST 3
# =============================================================================

def test_exact_gcd_factorization(
    states: list[State],
) -> int:

    print("=" * 90)
    print("TEST 3: EXACT INTEGER gcd FACTORIZATION")
    print("=" * 90)

    failures = 0

    factor_counter = Counter()

    for s in states:

        d = gcd(
            abs(s.A),
            abs(s.B),
        )

        shift = abs(s.n + s.c)

        Q = shift // d

        G = gcd(
            gcd(
                2 * abs(s.A),
                2 * abs(s.B),
            ),
            shift,
        )

        predicted = d * gcd(
            2,
            Q,
        )

        factor_counter[
            gcd(2, Q)
        ] += 1

        if G != predicted:

            failures += 1

    print(
        "gcd(2,Q) distribution:"
    )

    for key, count in sorted(
        factor_counter.items()
    ):
        print(
            f"    {key}: {count}"
        )

    print()
    print(
        f"checked={len(states)}"
    )

    print(
        f"failures={failures}"
    )

    print()

    return failures


# =============================================================================
# TEST 4
# =============================================================================

def test_parity_equivalence(
    states: list[State],
) -> int:

    print("=" * 90)
    print(
        "TEST 4: Q PARITY <-> EQUAL 2-ADIC RESIDUAL VALUATIONS"
    )
    print("=" * 90)

    failures = 0

    distribution = Counter()

    for s in states:

        d = gcd(
            abs(s.A),
            abs(s.B),
        )

        shift = abs(s.n + s.c)

        Q = shift // d

        equal_v2 = (
            v2(s.A) == v2(s.B)
        )

        predicted_even = equal_v2

        actual_even = (
            Q % 2 == 0
        )

        distribution[
            (s.frame, actual_even, equal_v2)
        ] += 1

        if actual_even != predicted_even:

            failures += 1

    print(
        "distribution:"
    )

    for key, count in sorted(
        distribution.items()
    ):
        print(
            f"    {key}: {count}"
        )

    print()
    print(
        f"checked={len(states)}"
    )

    print(
        f"failures={failures}"
    )

    print()

    return failures


# =============================================================================
# TEST 5
# =============================================================================

def test_d_times_equality_bit(
    states: list[State],
) -> int:

    print("=" * 90)
    print(
        "TEST 5: G = d * 2^[v2(A)=v2(B)]"
    )
    print("=" * 90)

    failures = 0

    distribution = Counter()

    for s in states:

        d = gcd(
            abs(s.A),
            abs(s.B),
        )

        G = gcd(
            gcd(
                2 * abs(s.A),
                2 * abs(s.B),
            ),
            abs(s.n + s.c),
        )

        equality = (
            v2(s.A) == v2(s.B)
        )

        predicted = d * (
            2 if equality else 1
        )

        distribution[
            (s.frame, equality)
        ] += 1

        if G != predicted:

            failures += 1

    print(
        "branch distribution:"
    )

    for key, count in sorted(
        distribution.items()
    ):
        print(
            f"    {key}: {count}"
        )

    print()
    print(
        f"checked={len(states)}"
    )

    print(
        f"failures={failures}"
    )

    print()

    return failures


# =============================================================================
# TEST 6
# =============================================================================

def test_one_sided_against_d(
    states: list[State],
) -> int:

    print("=" * 90)
    print(
        "TEST 6: ONE-SIDED gcd VS d"
    )
    print("=" * 90)

    failures = 0

    ratio_counter = Counter()

    examples = []

    for s in states:

        d = gcd(
            abs(s.A),
            abs(s.B),
        )

        shift = abs(s.n + s.c)

        gA = gcd(
            2 * abs(s.A),
            shift,
        )

        gB = gcd(
            2 * abs(s.B),
            shift,
        )

        G = d * (
            2
            if v2(s.A) == v2(s.B)
            else 1
        )

        if s.A != 0 and s.B != 0:

            if gA != G or gB != G:

                failures += 1

                if len(examples) < 20:
                    examples.append(
                        (
                            s,
                            d,
                            gA,
                            gB,
                            G,
                        )
                    )

        elif s.frame == "A" and s.A == 0:

            if gB != G:
                failures += 1

            if gA != 3 * G:
                failures += 1

        elif s.frame == "B" and s.B == 0:

            if gA != G:
                failures += 1

            if gB != 3 * G:
                failures += 1

        else:

            failures += 1

        ratio_counter[
            (
                s.frame,
                "equal" if gA == gB else (
                    "A=3B"
                    if gA == 3 * gB
                    else
                    "B=3A"
                    if gB == 3 * gA
                    else "other"
                ),
            )
        ] += 1

    print(
        "one-sided relation:"
    )

    for key, count in sorted(
        ratio_counter.items()
    ):
        print(
            f"    {key}: {count}"
        )

    if examples:

        print()
        print("generic counterexamples:")

        for s, d, gA, gB, G in examples:

            print(
                f"    n={s.n} "
                f"p={s.p} "
                f"q={s.q} "
                f"frame={s.frame} "
                f"d={d} "
                f"gA={gA} "
                f"gB={gB} "
                f"G={G}"
            )

    print()
    print(
        f"checked={len(states)}"
    )

    print(
        f"failures={failures}"
    )

    print()

    return failures


# =============================================================================
# TEST 7
# =============================================================================

def test_depth_from_factorization(
    states: list[State],
) -> int:

    print("=" * 90)
    print(
        "TEST 7: DEPTH FROM d AND QUOTIENT PARITY"
    )
    print("=" * 90)

    failures = 0

    for s in states:

        d = gcd(
            abs(s.A),
            abs(s.B),
        )

        shift = abs(s.n + s.c)

        Q = shift // d

        predicted = (
            v2(d)
            + (1 if Q % 2 == 0 else 0)
        )

        if predicted != s.depth:
            failures += 1

    print(
        f"checked={len(states)}"
    )

    print(
        f"failures={failures}"
    )

    print()

    return failures


# =============================================================================
# TEST 8
# =============================================================================

def test_residual_normalization(
    states: list[State],
) -> int:

    print("=" * 90)
    print(
        "TEST 8: NORMALIZED RESIDUAL PARITY"
    )
    print("=" * 90)

    failures = 0

    parity = Counter()

    for s in states:

        d = gcd(
            abs(s.A),
            abs(s.B),
        )

        if d == 0:
            failures += 1
            continue

        a = abs(s.A) // d
        b = abs(s.B) // d

        if gcd(a, b) != 1:

            failures += 1
            continue

        pa = a & 1
        pb = b & 1

        parity[
            (pa, pb)
        ] += 1

        if (pa, pb) == (0, 0):

            failures += 1

        if (pa, pb) not in {
            (0, 1),
            (1, 0),
            (1, 1),
        }:

            failures += 1

    print(
        "normalized parity distribution:"
    )

    for key, count in sorted(
        parity.items()
    ):
        print(
            f"    {key}: {count}"
        )

    print()
    print(
        f"checked={len(states)}"
    )

    print(
        f"failures={failures}"
    )

    print()

    return failures


# =============================================================================
# EXAMPLES
# =============================================================================

def print_examples(
    states: list[State],
) -> None:

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
        87,
        93,
        111,
        141,
        183,
        213,
        485879,
        5579767,
    ]

    by_n = {
        s.n: s
        for s in states
    }

    for n in wanted:

        s = by_n.get(n)

        if s is None:
            continue

        d = gcd(
            abs(s.A),
            abs(s.B),
        )

        shift = abs(s.n + s.c)

        Q = shift // d

        G = gcd(
            gcd(
                2 * abs(s.A),
                2 * abs(s.B),
            ),
            shift,
        )

        gA = gcd(
            2 * abs(s.A),
            shift,
        )

        gB = gcd(
            2 * abs(s.B),
            shift,
        )

        equality = (
            v2(s.A) == v2(s.B)
        )

        print()
        print(
            f"n={s.n} "
            f"p={s.p} "
            f"q={s.q} "
            f"frame={s.frame}"
        )

        print(
            f"    A={s.A} "
            f"B={s.B}"
        )

        print(
            f"    d=gcd(A,B)={d}"
        )

        print(
            f"    Q=(n+c)/d={Q}"
        )

        print(
            f"    Q parity={'even' if Q % 2 == 0 else 'odd'}"
        )

        print(
            f"    v2(A)={v2(s.A)} "
            f"v2(B)={v2(s.B)}"
        )

        print(
            f"    equal-v2={equality}"
        )

        print(
            f"    gcd(2,Q)={gcd(2,Q)}"
        )

        print(
            f"    G=d*gcd(2,Q)={d*gcd(2,Q)}"
        )

        print(
            f"    symmetric G={G}"
        )

        print(
            f"    gcd(2A,n+c)={gA}"
        )

        print(
            f"    gcd(2B,n+c)={gB}"
        )

        print(
            f"    depth=v2(G)={v2(G)}"
        )


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:

    PRIME_LIMIT = 6000

    print("=" * 90)
    print("EXPERIMENT 662 START")
    print("=" * 90)

    print()
    print(
        f"prime limit={PRIME_LIMIT}"
    )

    primes = sieve_primes(
        PRIME_LIMIT
    )

    print(
        f"odd primes={len(primes)}"
    )

    states = build_states(
        PRIME_LIMIT
    )

    print(
        f"semiprimes={len(states)}"
    )

    print()

    total_failures = 0

    total_failures += (
        test_baseline(states)
    )

    total_failures += (
        test_d_divides_shift(states)
    )

    total_failures += (
        test_exact_quotient(states)
    )

    total_failures += (
        test_exact_gcd_factorization(states)
    )

    total_failures += (
        test_parity_equivalence(states)
    )

    total_failures += (
        test_d_times_equality_bit(states)
    )

    total_failures += (
        test_one_sided_against_d(states)
    )

    total_failures += (
        test_depth_from_factorization(states)
    )

    total_failures += (
        test_residual_normalization(states)
    )

    print_examples(
        states
    )

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)

    print()

    print(
        "Let:"
    )

    print(
        "    d = gcd(A,B)"
    )

    print(
        "    Q = (n+c)/d"
    )

    print()

    print(
        "The tested exact decomposition is:"
    )

    print()

    print(
        "    gcd(2A,2B,n+c)"
    )

    print(
        "      ="
    )

    print(
        "    d * gcd(2,Q)"
    )

    print()

    print(
        "Since gcd(A/d,B/d)=1:"
    )

    print()

    print(
        "    Q odd"
    )

    print(
        "        => G=d"
    )

    print()

    print(
        "    Q even"
    )

    print(
        "        => G=2d"
    )

    print()

    print(
        "The quotient parity is expected to satisfy:"
    )

    print()

    print(
        "    Q even"
    )

    print(
        "        <=>"
    )

    print(
        "    v2(A)=v2(B)"
    )

    print()

    print(
        "Therefore:"
    )

    print()

    print(
        "    G"
    )

    print(
        "      ="
    )

    print(
        "    d * 2^[v2(A)=v2(B)]"
    )

    print()

    print(
        "and consequently:"
    )

    print()

    print(
        "    depth"
    )

    print(
        "      ="
    )

    print(
        "    v2(d) + [v2(A)=v2(B)]"
    )

    print()

    print(
        "This would connect the two previously separate"
    )

    print(
        "descriptions:"
    )

    print(
        "    residual gcd description"
    )

    print(
        "        <=>"
    )

    print(
        "    canonical integer gcd description"
    )

    print()

    print(
        "The one-sided gcd is then:"
    )

    print()

    print(
        "    generic A!=0,B!=0:"
    )

    print(
        "        gcd(2A,n+c)"
    )

    print(
        "          = G"
    )

    print(
        "        gcd(2B,n+c)"
    )

    print(
        "          = G"
    )

    print()

    print(
        "    Frame A, A=0:"
    )

    print(
        "        gcd(2A,n+9)=3G"
    )

    print(
        "        gcd(2B,n+9)=G"
    )

    print()

    print(
        "    Frame B, B=0:"
    )

    print(
        "        gcd(2A,n+3)=G"
    )

    print(
        "        gcd(2B,n+3)=3G"
    )

    print()

    print(
        f"TOTAL FAILURES={total_failures}"
    )

    if total_failures == 0:

        print(
            "STATUS=ALL INTEGER DECOMPOSITION TESTS PASSED"
        )

    else:

        print(
            "STATUS=COUNTEREXAMPLES FOUND"
        )

    print()
    print("=" * 90)
    print("EXPERIMENT 662 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
