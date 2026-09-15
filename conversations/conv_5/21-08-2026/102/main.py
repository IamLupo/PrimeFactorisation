#!/usr/bin/env python3

"""
==========================================================================================
EXPERIMENT 665

CANONICAL gcd = RESIDUAL gcd

Experiment 664 established:

    depth = v2(gcd(C, n+c)) + [ (n+c)/gcd(C,n+c) is even ]

with:

    FRAME A:
        C = q + 3
        c = 9

    FRAME B:
        C = p + 1
        c = 3

The examples strongly suggest the stronger integer identity:

    gcd(C, n+c) = gcd(A,B)

for the canonical residual C.

If this is exact, then the entire hierarchy collapses to:

    d = gcd(A,B)

    depth =
        v2(d)
        + [v2(A)=v2(B)]

and simultaneously:

    d = gcd(C,n+c).

This experiment explicitly tests that bridge.

It also checks the quotient identity:

    (n+c)/gcd(C,n+c)
      =
    (n+c)/gcd(A,B)

and therefore proves that the correction parity in
Experiment 664 is exactly the parity of the residual gcd
quotient.

The intended endpoint is the exact chain:

    factor residuals
        ->
    d = gcd(A,B)
        =
    gcd(C,n+c)
        ->
    depth.

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


# =============================================================================
# PRIME SIEVE
# =============================================================================

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

    C: int
    c: int

    depth: int


def make_state(p: int, q: int) -> State:

    n = p * q

    if n % 4 == 3:

        # FRAME A

        frame = "A"

        A = p - 3
        B = q + 3

        C = B
        c = 9

    else:

        # FRAME B

        frame = "B"

        A = p + 1
        B = q - 3

        C = A
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
        C=C,
        c=c,
        depth=depth,
    )


# =============================================================================
# STATE GENERATION
# =============================================================================

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

def test_canonical_gcd_bridge(
    states: list[State],
) -> int:

    print("=" * 90)
    print(
        "TEST 1: gcd(C,n+c) = gcd(A,B)"
    )
    print("=" * 90)

    failures = 0

    ratio_distribution = Counter()
    frame_distribution = Counter()

    for s in states:

        d = gcd(
            abs(s.A),
            abs(s.B),
        )

        H = gcd(
            abs(s.C),
            abs(s.n + s.c),
        )

        frame_distribution[
            (
                s.frame,
                H == d,
            )
        ] += 1

        if d != 0:

            ratio_distribution[
                H // d
                if H % d == 0
                else "noninteger"
            ] += 1

        if H != d:

            failures += 1

    print(
        "ratio H/d:"
    )

    for key, count in sorted(
        ratio_distribution.items(),
        key=lambda item: str(item[0]),
    ):

        print(
            f"    {key}: {count}"
        )

    print()

    print(
        "frame/equality distribution:"
    )

    for key, count in sorted(
        frame_distribution.items()
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

def test_quotient_identity(
    states: list[State],
) -> int:

    print("=" * 90)
    print(
        "TEST 2: EXACT QUOTIENT IDENTITY"
    )
    print("=" * 90)

    failures = 0

    for s in states:

        d = gcd(
            abs(s.A),
            abs(s.B),
        )

        H = gcd(
            abs(s.C),
            abs(s.n + s.c),
        )

        if d == 0 or H == 0:

            failures += 1
            continue

        Qd = (
            abs(s.n + s.c) // d
        )

        QH = (
            abs(s.n + s.c) // H
        )

        if Qd != QH:

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
# TEST 3
# =============================================================================

def test_symbolic_frame_relations(
    states: list[State],
) -> int:

    print("=" * 90)
    print(
        "TEST 3: FRAME-SPECIFIC gcd REDUCTION"
    )
    print("=" * 90)

    failures = 0

    for s in states:

        d = gcd(
            abs(s.A),
            abs(s.B),
        )

        H = gcd(
            abs(s.C),
            abs(s.n + s.c),
        )

        if s.frame == "A":

            # C = B
            #
            # n+9 = p(q+3) - 3(p-3)
            #
            # Therefore:
            #
            # gcd(B,n+9)
            #   = gcd(B,3A).

            symbolic = gcd(
                abs(s.B),
                abs(3 * s.A),
            )

        else:

            # C = A
            #
            # n+3 = q(p+1) - (q-3)
            #
            # Therefore:
            #
            # gcd(A,n+3)
            #   = gcd(A,B).

            symbolic = gcd(
                abs(s.A),
                abs(s.B),
            )

        if H != symbolic:

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
# TEST 4
# =============================================================================

def test_frame_a_factor_three(
    states: list[State],
) -> int:

    print("=" * 90)
    print(
        "TEST 4: FRAME-A FACTOR-3 CANCELLATION"
    )
    print("=" * 90)

    failures = 0

    categories = Counter()

    for s in states:

        if s.frame != "A":
            continue

        d = gcd(
            abs(s.A),
            abs(s.B),
        )

        H = gcd(
            abs(s.B),
            abs(s.n + 9),
        )

        if d == 0:

            failures += 1
            continue

        ratio = (
            H // d
            if H % d == 0
            else None
        )

        categories[
            (
                "A=0" if s.A == 0 else "generic",
                ratio,
            )
        ] += 1

        if H != d:

            failures += 1

    for key, count in sorted(
        categories.items(),
        key=lambda item: str(item[0]),
    ):

        print(
            f"    {key}: {count}"
        )

    print()

    print(
        "Algebraic reduction:"
    )

    print(
        "    gcd(B,n+9)"
    )

    print(
        "      = gcd(B,3A)"
    )

    print(
        "      = d * gcd(b,3a)"
    )

    print(
        "where gcd(a,b)=1."
    )

    print()

    print(
        "The only possible extra factor is 3."
    )

    print(
        "The prime-domain constraints must eliminate"
    )

    print(
        "that factor in the Frame-A states."
    )

    print()

    print(
        f"checked={sum("
        "1 for s in states if s.frame == 'A'"
        ")}"
    )

    print(
        f"failures={failures}"
    )

    print()

    return failures


# =============================================================================
# TEST 5
# =============================================================================

def test_frame_b_reduction(
    states: list[State],
) -> int:

    print("=" * 90)
    print(
        "TEST 5: FRAME-B gcd REDUCTION"
    )
    print("=" * 90)

    failures = 0

    for s in states:

        if s.frame != "B":
            continue

        d = gcd(
            abs(s.A),
            abs(s.B),
        )

        H = gcd(
            abs(s.A),
            abs(s.n + 3),
        )

        if H != d:

            failures += 1

    print(
        f"checked={sum("
        "1 for s in states if s.frame == 'B'"
        ")}"
    )

    print(
        f"failures={failures}"
    )

    print()

    return failures


# =============================================================================
# TEST 6
# =============================================================================

def test_depth_from_canonical_gcd(
    states: list[State],
) -> int:

    print("=" * 90)
    print(
        "TEST 6: DEPTH FROM gcd(C,n+c)"
    )
    print("=" * 90)

    failures = 0

    distribution = Counter()

    for s in states:

        d = gcd(
            abs(s.A),
            abs(s.B),
        )

        H = gcd(
            abs(s.C),
            abs(s.n + s.c),
        )

        Q = (
            abs(s.n + s.c) // H
        )

        predicted = (
            v2(H)
            +
            int(Q % 2 == 0)
        )

        equal_v2 = (
            v2(s.A) == v2(s.B)
        )

        distribution[
            (
                s.frame,
                "equal-v2" if equal_v2 else "unequal-v2",
                Q % 2,
            )
        ] += 1

        if predicted != s.depth:

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
# TEST 7
# =============================================================================

def test_odd_part_identity(
    states: list[State],
) -> int:

    print("=" * 90)
    print(
        "TEST 7: ODD PART gcd(C,n+c) = ODD PART gcd(A,B)"
    )
    print("=" * 90)

    failures = 0

    for s in states:

        d = gcd(
            abs(s.A),
            abs(s.B),
        )

        H = gcd(
            abs(s.C),
            abs(s.n + s.c),
        )

        if odd_part(d) != odd_part(H):

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

    buckets: dict[
        tuple[str, int, int],
        set[int],
    ] = {}

    for s in states:

        signature = (
            s.frame,
            v2(s.C),
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

        d = gcd(
            abs(s.A),
            abs(s.B),
        )

        H = gcd(
            abs(s.C),
            abs(s.n + s.c),
        )

        Q = (
            abs(s.n + s.c) // H
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
            f"    C={s.C}"
        )

        print(
            f"    n+c={abs(s.n+s.c)}"
        )

        print(
            f"    d=gcd(A,B)={d}"
        )

        print(
            f"    H=gcd(C,n+c)={H}"
        )

        print(
            f"    H/d="
            f"{H // d if d else 'undefined'}"
        )

        print(
            f"    odd(d)={odd_part(d)}"
        )

        print(
            f"    odd(H)={odd_part(H)}"
        )

        print(
            f"    Q=(n+c)/H={Q}"
        )

        print(
            f"    Q parity="
            f"{'even' if Q % 2 == 0 else 'odd'}"
        )

        print(
            f"    equal-v2="
            f"{v2(s.A) == v2(s.B)}"
        )

        print(
            f"    depth="
            f"v2(H)+[Q even]="
            f"{v2(H) + int(Q % 2 == 0)}"
        )

        print(
            f"    actual depth={s.depth}"
        )


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:

    PRIME_LIMIT = 6000

    print("=" * 90)
    print("EXPERIMENT 665 START")
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

    total_failures += test_baseline(
        states
    )

    total_failures += test_canonical_gcd_bridge(
        states
    )

    total_failures += test_quotient_identity(
        states
    )

    total_failures += test_symbolic_frame_relations(
        states
    )

    total_failures += test_frame_a_factor_three(
        states
    )

    total_failures += test_frame_b_reduction(
        states
    )

    total_failures += test_depth_from_canonical_gcd(
        states
    )

    total_failures += test_odd_part_identity(
        states
    )

    total_failures += test_canonical_signature(
        states
    )

    print_examples(
        states
    )

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)

    print()

    print(
        "Experiment 664 established:"
    )

    print(
        "    depth = v2(gcd(C,n+c))"
        " + [(n+c)/gcd(C,n+c) is even]."
    )

    print()

    print(
        "The key new candidate identity is:"
    )

    print()

    print(
        "    gcd(C,n+c) = gcd(A,B)."
    )

    print()

    print(
        "Thus the two descriptions become the same"
        " integer object:"
    )

    print()

    print(
        "    d = gcd(A,B)"
    )

    print(
        "      = gcd(C,n+c)."
    )

    print()

    print(
        "Then:"
    )

    print()

    print(
        "    Q = (n+c)/d"
    )

    print()

    print(
        "and:"
    )

    print()

    print(
        "    Q even"
    )

    print(
        "        <=>"
    )

    print(
        "    v2(A)=v2(B)."
    )

    print()

    print(
        "Therefore:"
    )

    print()

    print(
        "    depth"
    )

    print(
        "      = v2(d)"
    )

    print(
        "        + [Q even]"
    )

    print()

    print(
        "or equivalently:"
    )

    print()

    print(
        "    depth"
    )

    print(
        "      = v2(gcd(C,n+c))"
    )

    print(
        "        + [(n+c)/gcd(C,n+c) is even]."
    )

    print()

    print(
        "If the bridge passes, the canonical gcd is not"
    )

    print(
        "merely a convenient representation."
    )

    print(
        "It is exactly the original residual gcd d."
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
    print("EXPERIMENT 665 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
