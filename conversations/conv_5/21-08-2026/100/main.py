#!/usr/bin/env python3

"""
==========================================================================================
EXPERIMENT 663

CANONICAL ONE-SIDED INTEGER gcd THEOREM

Goal
----

Experiment 662 established:

    d = gcd(A,B)
    Q = (n+c)/d

    G = gcd(2A,2B,n+c)
      = d*gcd(2,Q)

and:

    depth = v2(G).

Experiment 661 showed that, generically,

    gcd(2A,n+c) = G
    gcd(2B,n+c) = G,

with only the degenerate zero-residual branches producing
a factor 3 on the non-canonical side.

This experiment asks whether a CANONICAL residual coordinate
can eliminate all exceptions.

Define:

    FRAME A:
        A = p-3
        B = q+3
        c = 9

        canonical C = B

    FRAME B:
        A = p+1
        B = q-3
        c = 3

        canonical C = A

Then test:

    Gc = gcd(2C,n+c)

Does:

    Gc = gcd(2A,2B,n+c)

hold for every state?

If yes:

    depth = v2(gcd(2C,n+c))

is a universal exact integer formula with no
A=0/B=0 exception.

Additional tests:

1. Exact canonical gcd == symmetric gcd.
2. Canonical gcd == d*gcd(2,Q).
3. Canonical gcd / d is always 1 or 2.
4. Q even iff canonical gcd = 2d.
5. depth from canonical gcd.
6. Canonical formula using only one residual coordinate.
7. Verify that the NON-canonical side is the only place
   where the factor 3 appears.
8. Search for a frame-independent selector rule for C.

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

    states: list[State] = []

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

        # Frame A uses B=q+3.
        return s.B

    if s.frame == "B":

        # Frame B uses A=p+1.
        return s.A

    raise ValueError(
        f"Unknown frame: {s.frame}"
    )


def noncanonical_residual(s: State) -> int:

    if s.frame == "A":

        return s.A

    if s.frame == "B":

        return s.B

    raise ValueError(
        f"Unknown frame: {s.frame}"
    )


# =============================================================================
# COMMON QUANTITIES
# =============================================================================

def residual_gcd(s: State) -> int:

    return gcd(
        abs(s.A),
        abs(s.B),
    )


def symmetric_gcd(s: State) -> int:

    return gcd(
        gcd(
            2 * abs(s.A),
            2 * abs(s.B),
        ),
        abs(s.n + s.c),
    )


def canonical_gcd(s: State) -> int:

    C = canonical_residual(s)

    return gcd(
        2 * abs(C),
        abs(s.n + s.c),
    )


def other_gcd(s: State) -> int:

    C = noncanonical_residual(s)

    return gcd(
        2 * abs(C),
        abs(s.n + s.c),
    )


# =============================================================================
# TEST 0
# =============================================================================

def test_baseline(states: list[State]) -> int:

    print("=" * 90)
    print("TEST 0: BASELINE")
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

def test_canonical_equals_symmetric(
    states: list[State],
) -> int:

    print("=" * 90)
    print("TEST 1: CANONICAL gcd == SYMMETRIC gcd")
    print("=" * 90)

    failures = 0

    first_examples = []

    for s in states:

        gc = canonical_gcd(s)
        gs = symmetric_gcd(s)

        if gc != gs:

            failures += 1

            if len(first_examples) < 20:

                first_examples.append(
                    (
                        s,
                        gc,
                        gs,
                    )
                )

    if first_examples:

        print(
            "counterexamples:"
        )

        for s, gc, gs in first_examples:

            print(
                f"    n={s.n} "
                f"p={s.p} "
                f"q={s.q} "
                f"frame={s.frame} "
                f"A={s.A} "
                f"B={s.B} "
                f"gc={gc} "
                f"gs={gs}"
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

def test_canonical_factorization(
    states: list[State],
) -> int:

    print("=" * 90)
    print("TEST 2: CANONICAL gcd = d * gcd(2,Q)")
    print("=" * 90)

    failures = 0

    distribution = Counter()

    for s in states:

        d = residual_gcd(s)

        shift = abs(
            s.n + s.c
        )

        if d == 0:

            failures += 1
            continue

        if shift % d != 0:

            failures += 1
            continue

        Q = shift // d

        gc = canonical_gcd(s)

        predicted = (
            d * gcd(2, Q)
        )

        distribution[
            gcd(2, Q)
        ] += 1

        if gc != predicted:

            failures += 1

    print(
        "canonical factor distribution:"
    )

    for key, count in sorted(
        distribution.items()
    ):

        print(
            f"    gcd(2,Q)={key}: {count}"
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

def test_canonical_ratio(
    states: list[State],
) -> int:

    print("=" * 90)
    print("TEST 3: CANONICAL gcd / d")
    print("=" * 90)

    failures = 0

    distribution = Counter()

    examples = []

    for s in states:

        d = residual_gcd(s)
        gc = canonical_gcd(s)

        if d == 0 or gc % d != 0:

            failures += 1

            if len(examples) < 20:

                examples.append(
                    (
                        s,
                        d,
                        gc,
                    )
                )

            continue

        ratio = gc // d

        distribution[
            ratio
        ] += 1

        if ratio not in (1, 2):

            failures += 1

            if len(examples) < 20:

                examples.append(
                    (
                        s,
                        d,
                        gc,
                    )
                )

    print(
        "ratio distribution:"
    )

    for key, count in sorted(
        distribution.items()
    ):

        print(
            f"    {key}: {count}"
        )

    if examples:

        print()
        print(
            "counterexamples:"
        )

        for s, d, gc in examples:

            print(
                f"    n={s.n} "
                f"p={s.p} "
                f"q={s.q} "
                f"frame={s.frame} "
                f"d={d} "
                f"Gc={gc}"
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

def test_q_parity(
    states: list[State],
) -> int:

    print("=" * 90)
    print(
        "TEST 4: Q PARITY <-> CANONICAL gcd / d"
    )
    print("=" * 90)

    failures = 0

    distribution = Counter()

    for s in states:

        d = residual_gcd(s)

        shift = abs(
            s.n + s.c
        )

        Q = shift // d

        gc = canonical_gcd(s)

        ratio = gc // d

        expected = 2 if (
            Q % 2 == 0
        ) else 1

        distribution[
            (
                "even"
                if Q % 2 == 0
                else "odd",
                ratio,
            )
        ] += 1

        if ratio != expected:

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

def test_canonical_depth(
    states: list[State],
) -> int:

    print("=" * 90)
    print(
        "TEST 5: DEPTH = v2(gcd(2C,n+c))"
    )
    print("=" * 90)

    failures = 0

    for s in states:

        predicted = v2(
            canonical_gcd(s)
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
# TEST 6
# =============================================================================

def test_noncanonical_factor3(
    states: list[State],
) -> int:

    print("=" * 90)
    print(
        "TEST 6: NON-CANONICAL FACTOR-3 EXCEPTION"
    )
    print("=" * 90)

    failures = 0

    distribution = Counter()

    examples = []

    for s in states:

        gc = canonical_gcd(s)
        go = other_gcd(s)

        if gc == 0:

            failures += 1
            continue

        if go == gc:

            relation = "equal"

        elif go == 3 * gc:

            relation = "other=3*canonical"

        elif gc == 3 * go:

            relation = "canonical=3*other"

        else:

            relation = "other"

            failures += 1

            if len(examples) < 20:

                examples.append(
                    (
                        s,
                        gc,
                        go,
                    )
                )

        distribution[
            (
                s.frame,
                relation,
            )
        ] += 1

    print(
        "relation distribution:"
    )

    for key, count in sorted(
        distribution.items()
    ):

        print(
            f"    {key}: {count}"
        )

    if examples:

        print()
        print(
            "unexpected relations:"
        )

        for s, gc, go in examples:

            print(
                f"    n={s.n} "
                f"p={s.p} "
                f"q={s.q} "
                f"frame={s.frame} "
                f"canonical={gc} "
                f"other={go}"
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

def test_canonical_residual_formula(
    states: list[State],
) -> int:

    print("=" * 90)
    print(
        "TEST 7: CANONICAL RESIDUAL FORMULA"
    )
    print("=" * 90)

    failures = 0

    for s in states:

        C = canonical_residual(s)

        if s.frame == "A":

            expected = s.B

        elif s.frame == "B":

            expected = s.A

        else:

            failures += 1
            continue

        if C != expected:

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

def test_canonical_vs_first_difference(
    states: list[State],
) -> int:

    print("=" * 90)
    print(
        "TEST 8: CANONICAL gcd VS FIRST-DIFFERING-BIT LAW"
    )
    print("=" * 90)

    failures = 0

    for s in states:

        C = canonical_residual(s)

        if s.frame == "A":

            anchor = -3

        else:

            anchor = 1  # placeholder replaced below

            # Canonical C = A = p+1.
            # Therefore C itself is the distance
            # from p=-1, represented by A.

        gc = canonical_gcd(s)

        if s.frame == "A":

            branch_v = v2(
                s.q + 3
            )

        else:

            branch_v = v2(
                s.p + 1
            )

        predicted = min(
            branch_v + 1,
            v2(s.n + s.c),
        )

        if predicted != v2(gc):

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
# TEST 9
# =============================================================================

def test_canonical_state_signature(
    states: list[State],
) -> int:

    print("=" * 90)
    print(
        "TEST 9: CANONICAL SIGNATURE"
    )
    print("=" * 90)

    buckets: dict[
        tuple[str, int, int],
        set[int],
    ] = {}

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
            "examples:"
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

        d = residual_gcd(s)

        shift = abs(
            s.n + s.c
        )

        Q = shift // d

        C = canonical_residual(s)

        gc = canonical_gcd(s)
        gs = symmetric_gcd(s)
        go = other_gcd(s)

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
            f"    symmetric G={gs}"
        )

        print(
            f"    canonical Gc={gc}"
        )

        print(
            f"    noncanonical Go={go}"
        )

        print(
            f"    v2(Gc)={v2(gc)}"
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
    print("EXPERIMENT 663 START")
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
        test_canonical_equals_symmetric(states)
    )

    total_failures += (
        test_canonical_factorization(states)
    )

    total_failures += (
        test_canonical_ratio(states)
    )

    total_failures += (
        test_q_parity(states)
    )

    total_failures += (
        test_canonical_depth(states)
    )

    total_failures += (
        test_noncanonical_factor3(states)
    )

    total_failures += (
        test_canonical_residual_formula(states)
    )

    total_failures += (
        test_canonical_vs_first_difference(states)
    )

    # Test 9 is an ambiguity count rather than a direct
    # failure count.
    total_failures += (
        test_canonical_state_signature(states)
    )

    print_examples(
        states
    )

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)

    print()

    print(
        "Define the canonical residual:"
    )

    print()

    print(
        "    FRAME A:"
    )

    print(
        "        C = B = q+3"
    )

    print()

    print(
        "    FRAME B:"
    )

    print(
        "        C = A = p+1"
    )

    print()

    print(
        "The central candidate theorem is:"
    )

    print()

    print(
        "    gcd(2C,n+c)"
    )

    print(
        "      ="
    )

    print(
        "    gcd(2A,2B,n+c)"
    )

    print()

    print(
        "Therefore:"
    )

    print()

    print(
        "    depth = v2(gcd(2C,n+c))"
    )

    print()

    print(
        "With:"
    )

    print(
        "    d = gcd(A,B)"
    )

    print(
        "    Q = (n+c)/d"
    )

    print()

    print(
        "the expected integer factorization is:"
    )

    print()

    print(
        "    gcd(2C,n+c)"
    )

    print(
        "      ="
    )

    print(
        "    d * gcd(2,Q)"
    )

    print()

    print(
        "Hence:"
    )

    print(
        "    Q odd  -> G=d"
    )

    print(
        "    Q even -> G=2d"
    )

    print()

    print(
        "The non-canonical residual is expected to differ"
    )

    print(
        "only on the two degenerate factor-3 branches."
    )

    print()

    print(
        "If all tests pass, the complete depth mechanism"
    )

    print(
        "has a canonical one-sided integer gcd form:"
    )

    print()

    print(
        "    FRAME A:"
    )

    print(
        "        depth = v2(gcd(2(q+3), n+9))"
    )

    print()

    print(
        "    FRAME B:"
    )

    print(
        "        depth = v2(gcd(2(p+1), n+3))"
    )

    print()

    print(
        f"TOTAL FAILURES={total_failures}"
    )

    if total_failures == 0:

        print(
            "STATUS=ALL CANONICAL gcd TESTS PASSED"
        )

    else:

        print(
            "STATUS=COUNTEREXAMPLES FOUND"
        )

    print()
    print("=" * 90)
    print("EXPERIMENT 663 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
