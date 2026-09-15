#!/usr/bin/env python3

"""
==========================================================================================
EXPERIMENT 664

CANONICAL gcd WITHOUT THE EXPLICIT FACTOR 2

Experiment 663 established exactly:

    FRAME A:
        C = q + 3
        depth = v2(gcd(2C, n+9))

    FRAME B:
        C = p + 1
        depth = v2(gcd(2C, n+3))

This experiment asks whether the leading factor 2 is merely
a compact way of encoding one additional parity bit.

Define:

    H = gcd(C, n+c)

    R = (n+c) / H

Candidate theorem:

    depth
      =
    v2(H) + [R is even]

Equivalently:

    depth
      =
    v2(gcd(C,n+c))
      +
    [2 | ((n+c)/gcd(C,n+c))].

This is potentially a cleaner canonical integer formulation.

Additional questions:

1. Is the correction bit exactly the old equality-v2 branch?
2. Does it equal the parity of Q=(n+c)/gcd(A,B)?
3. Is H's odd part irrelevant to depth?
4. Can the correction bit be determined directly from
   v2(C) and v2(n+c)?
5. Is there any smaller formula than:

       v2(H) + parity(R)

==========================================================================================
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import gcd, isqrt


INF = 10**9


# =============================================================================
# NUMBER THEORY
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

    sieve[0] = 0
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
        # A = p-3
        # B = q+3
        # c = 9

        frame = "A"

        A = p - 3
        B = q + 3
        c = 9

    else:

        # FRAME B
        #
        # A = p+1
        # B = q-3
        # c = 3

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

    states = []

    for i, p in enumerate(primes):

        for q in primes[i:]:

            states.append(
                make_state(p, q)
            )

    return states


# =============================================================================
# CANONICAL RESIDUAL
# =============================================================================

def canonical_residual(s: State) -> int:

    if s.frame == "A":

        return s.B

    return s.A


# =============================================================================
# TEST 0
# =============================================================================

def test_baseline(states: list[State]) -> int:

    print("=" * 90)
    print("TEST 0: BASELINE CANONICAL gcd LAW")
    print("=" * 90)

    failures = 0

    for s in states:

        C = canonical_residual(s)

        predicted = v2(
            gcd(
                2 * abs(C),
                abs(s.n + s.c),
            )
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
# TEST 1
# =============================================================================

def test_half_gcd_formula(states: list[State]) -> int:

    print("=" * 90)
    print(
        "TEST 1: depth = v2(H) + [R even]"
    )
    print("=" * 90)

    failures = 0

    distribution = Counter()

    for s in states:

        C = canonical_residual(s)

        shift = abs(
            s.n + s.c
        )

        H = gcd(
            abs(C),
            shift,
        )

        if H == 0:

            failures += 1
            continue

        R = shift // H

        correction = int(
            R % 2 == 0
        )

        predicted = (
            v2(H) + correction
        )

        distribution[
            (
                correction,
                s.depth - v2(H),
            )
        ] += 1

        if predicted != s.depth:

            failures += 1

    print(
        "correction distribution:"
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
# TEST 2
# =============================================================================

def test_correction_equals_equal_v2(
    states: list[State],
) -> int:

    print("=" * 90)
    print(
        "TEST 2: CORRECTION BIT <-> EQUAL RESIDUAL v2"
    )
    print("=" * 90)

    failures = 0

    distribution = Counter()

    for s in states:

        C = canonical_residual(s)

        shift = abs(
            s.n + s.c
        )

        H = gcd(
            abs(C),
            shift,
        )

        R = shift // H

        correction = int(
            R % 2 == 0
        )

        equal_v2 = int(
            v2(s.A) == v2(s.B)
        )

        distribution[
            (
                correction,
                equal_v2,
            )
        ] += 1

        if correction != equal_v2:

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
# TEST 3
# =============================================================================

def test_correction_equals_Q_parity(
    states: list[State],
) -> int:

    print("=" * 90)
    print(
        "TEST 3: CORRECTION BIT <-> PARITY OF Q"
    )
    print("=" * 90)

    failures = 0

    distribution = Counter()

    for s in states:

        C = canonical_residual(s)

        shift = abs(
            s.n + s.c
        )

        H = gcd(
            abs(C),
            shift,
        )

        d = gcd(
            abs(s.A),
            abs(s.B),
        )

        Q = shift // d
        R = shift // H

        correction = int(
            R % 2 == 0
        )

        q_even = int(
            Q % 2 == 0
        )

        distribution[
            (
                correction,
                q_even,
            )
        ] += 1

        if correction != q_even:

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
# TEST 4
# =============================================================================

def test_h_v2_relation(
    states: list[State],
) -> int:

    print("=" * 90)
    print(
        "TEST 4: v2(H) = MIN(v2(C), v2(n+c))"
    )
    print("=" * 90)

    failures = 0

    for s in states:

        C = canonical_residual(s)

        shift = abs(
            s.n + s.c
        )

        H = gcd(
            abs(C),
            shift,
        )

        predicted = min(
            v2(C),
            v2(shift),
        )

        if v2(H) != predicted:

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
# TEST 5
# =============================================================================

def test_odd_part_irrelevance(
    states: list[State],
) -> int:

    print("=" * 90)
    print(
        "TEST 5: ODD PART OF H IS IRRELEVANT"
    )
    print("=" * 90)

    failures = 0

    depth_odd_parts = {}

    for s in states:

        C = canonical_residual(s)

        shift = abs(
            s.n + s.c
        )

        H = gcd(
            abs(C),
            shift,
        )

        key = (
            s.frame,
            s.depth,
        )

        depth_odd_parts.setdefault(
            key,
            set(),
        ).add(
            odd_part(H)
        )

    print(
        "distinct odd parts:"
    )

    for key in sorted(
        depth_odd_parts
    ):

        values = depth_odd_parts[key]

        print(
            f"    {key}: "
            f"count={len(values)} "
            f"sample={sorted(values)[:10]}"
        )

    # This test intentionally does not expect
    # the odd part to be unique. It only verifies
    # that changing it does not change v2(H)
    # or the final depth.
    #
    # Therefore the actual correctness condition
    # is checked directly below.

    for s in states:

        C = canonical_residual(s)

        shift = abs(
            s.n + s.c
        )

        H = gcd(
            abs(C),
            shift,
        )

        reconstructed = (
            v2(H)
            +
            int(
                (shift // H) % 2 == 0
            )
        )

        if reconstructed != s.depth:

            failures += 1

    print()

    print(
        f"failures={failures}"
    )

    print()

    return failures


# =============================================================================
# TEST 6
# =============================================================================

def test_no_explicit_two_needed(
    states: list[State],
) -> int:

    print("=" * 90)
    print(
        "TEST 6: NO EXPLICIT FACTOR-2 NEEDED"
    )
    print("=" * 90)

    failures = 0

    for s in states:

        C = canonical_residual(s)

        shift = abs(
            s.n + s.c
        )

        H = gcd(
            abs(C),
            shift,
        )

        R = shift // H

        # Main candidate:
        #
        # depth = v2(H) + parity(R)

        predicted = (
            v2(H)
            +
            int(R % 2 == 0)
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
# TEST 7
# =============================================================================

def test_min_formula(
    states: list[State],
) -> int:

    print("=" * 90)
    print(
        "TEST 7: EQUIVALENT MIN FORMULA"
    )
    print("=" * 90)

    failures = 0

    for s in states:

        C = canonical_residual(s)

        cval = v2(C)
        w = v2(
            s.n + s.c
        )

        predicted = min(
            cval + 1,
            w,
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

def test_canonical_signature(
    states: list[State],
) -> int:

    print("=" * 90)
    print(
        "TEST 8: CANONICAL (v2(C), v2(n+c)) SIGNATURE"
    )
    print("=" * 90)

    buckets = {}

    for s in states:

        C = canonical_residual(s)

        signature = (
            s.frame,
            v2(C),
            v2(s.n + s.c),
        )

        buckets.setdefault(
            signature,
            set(),
        ).add(
            s.depth
        )

    ambiguous = {
        key: values
        for key, values in buckets.items()
        if len(values) > 1
    }

    print(
        f"signatures={len(buckets)}"
    )

    print(
        f"ambiguous={len(ambiguous)}"
    )

    if ambiguous:

        print(
            "counterexamples:"
        )

        for key, values in list(
            ambiguous.items()
        )[:20]:

            print(
                f"    {key} -> "
                f"{sorted(values)}"
            )

    print()

    return len(ambiguous)


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

        C = canonical_residual(s)

        shift = abs(
            s.n + s.c
        )

        H = gcd(
            abs(C),
            shift,
        )

        R = shift // H

        d = gcd(
            abs(s.A),
            abs(s.B),
        )

        Q = shift // d

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
            f"    canonical C={C}"
        )

        print(
            f"    n+c={shift}"
        )

        print(
            f"    H=gcd(C,n+c)={H}"
        )

        print(
            f"    v2(H)={v2(H)}"
        )

        print(
            f"    R=(n+c)/H={R}"
        )

        print(
            f"    R parity="
            f"{'even' if R % 2 == 0 else 'odd'}"
        )

        print(
            f"    d=gcd(A,B)={d}"
        )

        print(
            f"    Q=(n+c)/d={Q}"
        )

        print(
            f"    Q parity="
            f"{'even' if Q % 2 == 0 else 'odd'}"
        )

        print(
            f"    depth candidate="
            f"v2(H)+[R even]="
            f"{v2(H) + int(R % 2 == 0)}"
        )

        print(
            f"    depth={s.depth}"
        )


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:

    PRIME_LIMIT = 6000

    print("=" * 90)
    print("EXPERIMENT 664 START")
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
        test_half_gcd_formula(states)
    )

    total_failures += (
        test_correction_equals_equal_v2(states)
    )

    total_failures += (
        test_correction_equals_Q_parity(states)
    )

    total_failures += (
        test_h_v2_relation(states)
    )

    total_failures += (
        test_odd_part_irrelevance(states)
    )

    total_failures += (
        test_no_explicit_two_needed(states)
    )

    total_failures += (
        test_min_formula(states)
    )

    total_failures += (
        test_canonical_signature(states)
    )

    print_examples(
        states
    )

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)

    print()

    print(
        "Experiment 663 established the canonical formula:"
    )

    print()

    print(
        "    depth = v2(gcd(2C,n+c))"
    )

    print()

    print(
        "where:"
    )

    print(
        "    FRAME A: C=q+3, c=9"
    )

    print(
        "    FRAME B: C=p+1, c=3"
    )

    print()

    print(
        "This experiment removes the explicit factor 2."
    )

    print()

    print(
        "Define:"
    )

    print(
        "    H = gcd(C,n+c)"
    )

    print(
        "    R = (n+c)/H"
    )

    print()

    print(
        "Then the candidate identity is:"
    )

    print()

    print(
        "    depth = v2(H) + [R is even]"
    )

    print()

    print(
        "Equivalently:"
    )

    print()

    print(
        "    depth = min(v2(C)+1, v2(n+c))"
    )

    print()

    print(
        "The correction bit should coincide with:"
    )

    print()

    print(
        "    [v2(A)=v2(B)]"
    )

    print(
        "and with:"
    )

    print()

    print(
        "    [Q=(n+c)/gcd(A,B) is even]."
    )

    print()

    print(
        "If all tests pass, the canonical depth theorem"
    )

    print(
        "has an especially compact form:"
    )

    print()

    print(
        "    depth"
    )

    print(
        "      = v2(gcd(C,n+c))"
    )

    print(
        "        + parity((n+c)/gcd(C,n+c))."
    )

    print()

    print(
        "This is an integer-gcd theorem with only one"
    )

    print(
        "binary correction bit."
    )

    print()

    print(
        f"TOTAL FAILURES={total_failures}"
    )

    if total_failures == 0:

        print(
            "STATUS=ALL TESTS PASSED"
        )

    else:

        print(
            "STATUS=COUNTEREXAMPLES FOUND"
        )

    print()

    print("=" * 90)
    print("EXPERIMENT 664 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
