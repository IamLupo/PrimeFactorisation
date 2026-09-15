#!/usr/bin/env python3

"""
==========================================================================================
EXPERIMENT 661

CANONICAL ONE-SIDED INTEGER gcd THEOREM

Established:

    G = gcd(2A, 2B, n+c)

    depth = v2(G)

Experiment 660 found empirically:

    gA = gcd(2A, n+c)
    gB = gcd(2B, n+c)

and almost always:

    gA = gB.

The only exceptions were the degenerate residual branches:

    FRAME A, A=0:
        gA = 3*gB

    FRAME B, B=0:
        gB = 3*gA.

This experiment tries to establish the exact integer theorem:

    GENERIC:
        A != 0 and B != 0
        =>
        gcd(2A,n+c) = gcd(2B,n+c)

    FRAME A, A=0:
        gcd(2A,n+c) = 3*gcd(2B,n+c)

    FRAME B, B=0:
        gcd(2B,n+c) = 3*gcd(2A,n+c).

If this passes, then the complete 2-adic hierarchy has a
canonical one-sided integer gcd representation:

    generic:
        G = gcd(2A,n+c) = gcd(2B,n+c)

    Frame A, A=0:
        G = gcd(2B,n+c)

    Frame B, B=0:
        G = gcd(2A,n+c).

The experiment also checks the algebraic source of the
factor 3 and verifies that it cannot occur elsewhere.
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
        # c = 9

        frame = "A"
        A = p - 3
        B = q + 3
        c = 9

    else:

        # FRAME B
        #
        # A = p + 1
        # B = q - 3
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

def test_one_sided_gcds(states: list[State]) -> int:

    print("=" * 90)
    print("TEST 1: ONE-SIDED gcd VALUES")
    print("=" * 90)

    failures = 0

    equality = 0
    A_zero = 0
    B_zero = 0
    generic = 0

    ratio_counter: Counter[
        tuple[str, int | str]
    ] = Counter()

    examples: list[
        tuple[
            State,
            int,
            int,
            int | str,
        ]
    ] = []

    for s in states:

        shift = abs(s.n + s.c)

        gA = gcd(
            2 * abs(s.A),
            shift,
        )

        gB = gcd(
            2 * abs(s.B),
            shift,
        )

        if s.A == 0:
            A_zero += 1

            if gB == 0:
                ratio = "undefined"
            else:
                ratio = (
                    gA // gB
                    if gA % gB == 0
                    else "noninteger"
                )

            if gA != 3 * gB:
                failures += 1

            ratio_counter[
                (s.frame, ratio)
            ] += 1

            continue

        if s.B == 0:
            B_zero += 1

            if gA == 0:
                ratio = "undefined"
            else:
                ratio = (
                    gB // gA
                    if gB % gA == 0
                    else "noninteger"
                )

            if gB != 3 * gA:
                failures += 1

            ratio_counter[
                (s.frame, ratio)
            ] += 1

            continue

        generic += 1

        if gA != gB:
            failures += 1

            if len(examples) < 20:
                examples.append(
                    (
                        s,
                        gA,
                        gB,
                        "generic mismatch",
                    )
                )

        equality += int(gA == gB)

        ratio_counter[
            (s.frame, 1)
        ] += 1

    print(
        f"generic states={generic}"
    )

    print(
        f"A=0 states={A_zero}"
    )

    print(
        f"B=0 states={B_zero}"
    )

    print(
        f"generic exact equality={equality}"
    )

    print()

    print("ratio classification:")

    for key, count in sorted(
        ratio_counter.items(),
        key=lambda x: (x[0][0], str(x[0][1])),
    ):
        print(
            f"    {key}: {count}"
        )

    if examples:

        print()
        print("generic counterexamples:")

        for s, gA, gB, reason in examples:

            print(
                f"    n={s.n} "
                f"p={s.p} "
                f"q={s.q} "
                f"frame={s.frame} "
                f"A={s.A} "
                f"B={s.B} "
                f"gA={gA} "
                f"gB={gB} "
                f"{reason}"
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

def test_canonical_gcd_equals_symmetric(
    states: list[State],
) -> int:

    print("=" * 90)
    print("TEST 2: CANONICAL ONE-SIDED gcd = SYMMETRIC gcd")
    print("=" * 90)

    failures = 0

    for s in states:

        shift = abs(s.n + s.c)

        gA = gcd(
            2 * abs(s.A),
            shift,
        )

        gB = gcd(
            2 * abs(s.B),
            shift,
        )

        G = gcd(
            gcd(
                2 * abs(s.A),
                2 * abs(s.B),
            ),
            shift,
        )

        if s.A != 0 and s.B != 0:

            if gA != G:
                failures += 1

            if gB != G:
                failures += 1

        elif s.frame == "A" and s.A == 0:

            if gB != G:
                failures += 1

        elif s.frame == "B" and s.B == 0:

            if gA != G:
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

def test_factor_three_uniqueness(
    states: list[State],
) -> int:

    print("=" * 90)
    print("TEST 3: FACTOR-3 UNIQUENESS")
    print("=" * 90)

    failures = 0

    unexpected = 0

    examples: list[State] = []

    for s in states:

        shift = abs(s.n + s.c)

        gA = gcd(
            2 * abs(s.A),
            shift,
        )

        gB = gcd(
            2 * abs(s.B),
            shift,
        )

        if gA == gB:
            ratio = 1

        elif gA % gB == 0:
            ratio = gA // gB

        elif gB % gA == 0:
            ratio = -(gB // gA)

        else:
            ratio = None

        if s.A != 0 and s.B != 0:

            if ratio != 1:

                unexpected += 1
                failures += 1

                if len(examples) < 20:
                    examples.append(s)

        elif s.frame == "A" and s.A == 0:

            if ratio != 3:

                unexpected += 1
                failures += 1

                if len(examples) < 20:
                    examples.append(s)

        elif s.frame == "B" and s.B == 0:

            if ratio != -3:

                unexpected += 1
                failures += 1

                if len(examples) < 20:
                    examples.append(s)

        else:

            failures += 1

    print(
        f"unexpected ratio states={unexpected}"
    )

    if examples:

        print()
        print("counterexamples:")

        for s in examples:

            shift = abs(s.n + s.c)

            gA = gcd(
                2 * abs(s.A),
                shift,
            )

            gB = gcd(
                2 * abs(s.B),
                shift,
            )

            print(
                f"    n={s.n} "
                f"p={s.p} "
                f"q={s.q} "
                f"frame={s.frame} "
                f"A={s.A} "
                f"B={s.B} "
                f"gA={gA} "
                f"gB={gB}"
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

def test_v2_and_depth(
    states: list[State],
) -> int:

    print("=" * 90)
    print("TEST 4: ONE-SIDED gcd v2 = DEPTH")
    print("=" * 90)

    failures = 0

    for s in states:

        shift = abs(s.n + s.c)

        gA = gcd(
            2 * abs(s.A),
            shift,
        )

        gB = gcd(
            2 * abs(s.B),
            shift,
        )

        if v2(gA) != s.depth:
            failures += 1

        if v2(gB) != s.depth:
            failures += 1

    print(
        f"checked={2 * len(states)}"
    )

    print(
        f"failures={failures}"
    )

    print()

    return failures


# =============================================================================
# TEST 5
# =============================================================================

def test_odd_parts_match(
    states: list[State],
) -> int:

    print("=" * 90)
    print("TEST 5: GENERIC ODD-PART EQUALITY")
    print("=" * 90)

    failures = 0

    mismatches = Counter()

    for s in states:

        shift = abs(s.n + s.c)

        gA = gcd(
            2 * abs(s.A),
            shift,
        )

        gB = gcd(
            2 * abs(s.B),
            shift,
        )

        if s.A == 0 or s.B == 0:
            continue

        if odd_part(gA) != odd_part(gB):

            failures += 1

            mismatches[
                (
                    s.frame,
                    odd_part(gA),
                    odd_part(gB),
                )
            ] += 1

    print(
        f"generic states="
        f"{sum(1 for s in states if s.A != 0 and s.B != 0)}"
    )

    print(
        f"distinct odd-part mismatches="
        f"{len(mismatches)}"
    )

    if mismatches:

        print()
        print(
            "first mismatches:"
        )

        for key, count in list(
            mismatches.items()
        )[:20]:

            print(
                f"    {key}: {count}"
            )

    print()
    print(
        f"failures={failures}"
    )

    print()

    return failures


# =============================================================================
# TEST 6
# =============================================================================

def test_algebraic_relations(
    states: list[State],
) -> int:

    print("=" * 90)
    print("TEST 6: ALGEBRAIC ONE-SIDED gcd RELATIONS")
    print("=" * 90)

    failures = 0

    generic_checked = 0
    A_zero_checked = 0
    B_zero_checked = 0

    for s in states:

        shift = abs(s.n + s.c)

        gA = gcd(
            2 * abs(s.A),
            shift,
        )

        gB = gcd(
            2 * abs(s.B),
            shift,
        )

        # -------------------------------------------------------------
        # Generic
        # -------------------------------------------------------------

        if s.A != 0 and s.B != 0:

            generic_checked += 1

            # The difference between 2A and 2B must divide
            # every common divisor after the appropriate
            # frame identity is applied.

            if s.frame == "A":

                # n+9 = AB + 3(B-A)
                #
                # Hence every common divisor of:
                #
                #     2A and n+9
                #
                # and every common divisor of:
                #
                #     2B and n+9
                #
                # has the same odd and 2-adic content.

                lhs = gcd(
                    gA,
                    2 * abs(s.A),
                )

                rhs = gcd(
                    gB,
                    2 * abs(s.B),
                )

            else:

                # Frame B:
                #
                # n+3 = AB + 3A - B
                #
                # The same common-gcd mechanism applies.

                lhs = gcd(
                    gA,
                    2 * abs(s.A),
                )

                rhs = gcd(
                    gB,
                    2 * abs(s.B),
                )

            if lhs != abs(gA):
                failures += 1

            if rhs != abs(gB):
                failures += 1

            if gA != gB:
                failures += 1

        # -------------------------------------------------------------
        # Frame A zero branch
        # -------------------------------------------------------------

        elif s.frame == "A" and s.A == 0:

            A_zero_checked += 1

            # Here:
            #
            #   p=3
            #   A=0
            #   B=q+3
            #
            # n+9 = 3(q+3) = 3B.
            #
            # Therefore:
            #
            #   gcd(2A,n+9)=n+9=3B
            #
            # while:
            #
            #   gcd(2B,n+9)=B.

            expectedA = 3 * abs(s.B)
            expectedB = abs(s.B)

            if gA != expectedA:
                failures += 1

            if gB != expectedB:
                failures += 1

        # -------------------------------------------------------------
        # Frame B zero branch
        # -------------------------------------------------------------

        elif s.frame == "B" and s.B == 0:

            B_zero_checked += 1

            # Here:
            #
            #   q=3
            #   B=0
            #   A=p+1
            #
            # n+3 = 3(p+1) = 3A.

            expectedA = abs(s.A)
            expectedB = 3 * abs(s.A)

            if gA != expectedA:
                failures += 1

            if gB != expectedB:
                failures += 1

        else:

            failures += 1

    print(
        f"generic checked={generic_checked}"
    )

    print(
        f"A=0 checked={A_zero_checked}"
    )

    print(
        f"B=0 checked={B_zero_checked}"
    )

    print()
    print(
        f"failures={failures}"
    )

    print()

    return failures


# =============================================================================
# TEST 7
# =============================================================================

def test_canonical_depth_formula(
    states: list[State],
) -> int:

    print("=" * 90)
    print("TEST 7: CANONICAL gcd -> DEPTH")
    print("=" * 90)

    failures = 0

    for s in states:

        shift = abs(s.n + s.c)

        if s.A != 0:
            G = gcd(
                2 * abs(s.A),
                shift,
            )

        else:
            G = gcd(
                2 * abs(s.B),
                shift,
            )

        if v2(G) != s.depth:
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

        shift = abs(s.n + s.c)

        gA = gcd(
            2 * abs(s.A),
            shift,
        )

        gB = gcd(
            2 * abs(s.B),
            shift,
        )

        G = gcd(
            gcd(
                2 * abs(s.A),
                2 * abs(s.B),
            ),
            shift,
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
            f"B={s.B} "
            f"n+c={s.n+s.c}"
        )

        print(
            f"    gcd(2A,n+c)={gA}"
        )

        print(
            f"    gcd(2B,n+c)={gB}"
        )

        print(
            f"    symmetric G={G}"
        )

        print(
            f"    gA/gB="
            f"{gA}/{gB}"
        )

        print(
            f"    odd(gA)={odd_part(gA)}"
        )

        print(
            f"    odd(gB)={odd_part(gB)}"
        )

        print(
            f"    v2(gA)={v2(gA)}"
        )

        print(
            f"    v2(gB)={v2(gB)}"
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
    print("EXPERIMENT 661 START")
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
        test_one_sided_gcds(states)
    )

    total_failures += (
        test_canonical_gcd_equals_symmetric(states)
    )

    total_failures += (
        test_factor_three_uniqueness(states)
    )

    total_failures += (
        test_v2_and_depth(states)
    )

    total_failures += (
        test_odd_parts_match(states)
    )

    total_failures += (
        test_algebraic_relations(states)
    )

    total_failures += (
        test_canonical_depth_formula(states)
    )

    print_examples(
        states
    )

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)

    print()

    print(
        "The exact depth identity remains:"
    )

    print()

    print(
        "    depth = v2(gcd(2A,2B,n+c))"
    )

    print()

    print(
        "The new candidate integer theorem is:"
    )

    print()

    print(
        "    GENERIC A!=0,B!=0:"
    )

    print(
        "        gcd(2A,n+c)"
    )

    print(
        "          ="
    )

    print(
        "        gcd(2B,n+c)"
    )

    print(
        "          ="
    )

    print(
        "        gcd(2A,2B,n+c)"
    )

    print()

    print(
        "The only expected exceptions are:"
    )

    print()

    print(
        "    FRAME A, A=0:"
    )

    print(
        "        gcd(2A,n+9)"
    )

    print(
        "          ="
    )

    print(
        "        3*gcd(2B,n+9)"
    )

    print()

    print(
        "    FRAME B, B=0:"
    )

    print(
        "        gcd(2B,n+3)"
    )

    print(
        "          ="
    )

    print(
        "        3*gcd(2A,n+3)"
    )

    print()

    print(
        "If all tests pass, the canonical depth computation"
    )

    print(
        "can use exactly ONE residual coordinate:"
    )

    print()

    print(
        "    generic:"
    )

    print(
        "        depth = v2(gcd(2A,n+c))"
    )

    print()

    print(
        "    A=0:"
    )

    print(
        "        depth = v2(gcd(2B,n+9))"
    )

    print()

    print(
        "    B=0:"
    )

    print(
        "        depth = v2(gcd(2A,n+3))"
    )

    print()

    print(
        "This would collapse the symmetric gcd"
    )

    print(
        "to a canonical one-sided integer object."
    )

    print()

    print(
        f"TOTAL FAILURES={total_failures}"
    )

    if total_failures == 0:

        print(
            "STATUS=ALL CORE THEOREMS PASSED"
        )

    else:

        print(
            "STATUS=COUNTEREXAMPLES FOUND"
        )

    print()
    print("=" * 90)
    print("EXPERIMENT 661 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
