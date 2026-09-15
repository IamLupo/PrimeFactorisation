#!/usr/bin/env python3

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from math import gcd, isqrt


INF = 10**9
PRIME_LIMIT = 6000


# =============================================================================
# BASIC HELPERS
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
            count = ((limit - start) // p) + 1
            sieve[start:limit + 1:p] = b"\x00" * count

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


def build_state(p: int, q: int) -> State:
    n = p * q

    if n % 4 == 3:
        # FRAME A
        #
        # A = p-3
        # B = q+3
        # C = B = q+3
        # c = 9

        frame = "A"

        A = p - 3
        B = q + 3
        C = B
        c = 9

    else:
        # FRAME B
        #
        # A = p+1
        # B = q-3
        # C = A = p+1
        # c = 3

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


def build_states(primes: list[int]) -> list[State]:
    states: list[State] = []

    for i, p in enumerate(primes):
        for q in primes[i:]:
            states.append(build_state(p, q))

    return states


# =============================================================================
# DIVISOR GENERATION
# =============================================================================

def divisors(n: int) -> list[int]:
    n = abs(n)

    if n == 0:
        return []

    small: list[int] = []
    large: list[int] = []

    r = isqrt(n)

    for d in range(1, r + 1):
        if n % d != 0:
            continue

        small.append(d)

        other = n // d

        if other != d:
            large.append(other)

    return small + list(reversed(large))


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

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 1
# =============================================================================

def test_exact_canonical_gcd(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 1: CANONICAL gcd IDENTITY")
    print("=" * 90)

    failures = 0

    for s in states:
        d = gcd(abs(s.A), abs(s.B))
        H = gcd(abs(s.C), abs(s.n + s.c))

        if H != d:
            failures += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 2
#
# Enumerate divisors of n+c and inspect where H sits.
# =============================================================================

def test_divisor_position(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 2: DIVISOR-LATTICE POSITION OF H")
    print("=" * 90)

    failures = 0

    distributions = Counter()

    for s in states:
        value = abs(s.n + s.c)

        if value == 0:
            failures += 1
            continue

        H = gcd(abs(s.C), value)
        h_v2 = v2(H)

        ds = divisors(value)

        same_v2 = [
            d
            for d in ds
            if v2(d) == h_v2
        ]

        try:
            rank = same_v2.index(H) + 1
        except ValueError:
            failures += 1
            continue

        distributions[
            (
                s.frame,
                h_v2,
                len(same_v2),
                rank,
            )
        ] += 1

    print("Number of divisors sharing v2(H):")

    summary: Counter[tuple[str, int]] = Counter()

    for key, count in distributions.items():
        frame, h_v2, width, rank = key
        summary[(frame, width)] += count

    for key, count in sorted(summary.items()):
        print(f"    {key}: {count}")

    print()
    print("Question:")
    print("    Is H uniquely determined by its 2-adic level?")
    print()

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 3
#
# Compare H with divisors d having the same 2-adic valuation
# and simple quotient parity.
# =============================================================================

def test_quotient_signature(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 3: DIVISOR + QUOTIENT SIGNATURE")
    print("=" * 90)

    failures = 0

    ambiguous = 0
    exact = 0

    examples: list[tuple] = []

    for s in states:
        value = abs(s.n + s.c)

        if value == 0:
            failures += 1
            continue

        H = gcd(abs(s.C), value)
        target_v2 = v2(H)

        ds = divisors(value)

        candidates = []

        for d in ds:
            if v2(d) != target_v2:
                continue

            quotient = value // d

            signature = (
                v2(d),
                quotient % 2,
                quotient % 3,
                quotient % 4,
            )

            candidates.append((d, signature))

        target_q = value // H

        target_signature = (
            target_v2,
            target_q % 2,
            target_q % 3,
            target_q % 4,
        )

        matching = [
            d
            for d, sig in candidates
            if sig == target_signature
        ]

        if H not in matching:
            failures += 1
            continue

        if len(matching) == 1:
            exact += 1
        else:
            ambiguous += 1

            if len(examples) < 10:
                examples.append(
                    (
                        s.n,
                        s.frame,
                        H,
                        target_v2,
                        target_q,
                        matching[:10],
                    )
                )

    print(f"exact={exact}")
    print(f"ambiguous={ambiguous}")

    if examples:
        print()
        print("first ambiguous examples:")

        for item in examples:
            print(
                "    "
                f"n={item[0]} "
                f"frame={item[1]} "
                f"H={item[2]} "
                f"v2={item[3]} "
                f"Q={item[4]} "
                f"candidates={item[5]}"
            )

    print()
    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 4
#
# Does H have a maximality property?
# =============================================================================

def test_maximal_divisor_property(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 4: MAXIMAL DIVISOR CHARACTERIZATION")
    print("=" * 90)

    failures = 0
    exact = 0
    counterexamples = []

    for s in states:
        value = abs(s.n + s.c)

        H = gcd(abs(s.C), value)
        target_v2 = v2(H)

        ds = divisors(value)

        candidates = [
            d
            for d in ds
            if v2(d) == target_v2
        ]

        larger_same_v2 = [
            d
            for d in candidates
            if d > H
        ]

        if not larger_same_v2:
            exact += 1
        else:
            if len(counterexamples) < 20:
                counterexamples.append(
                    (
                        s.n,
                        s.p,
                        s.q,
                        s.frame,
                        H,
                        target_v2,
                        larger_same_v2[:8],
                    )
                )

    print(f"H maximal among same-v2 divisors={exact}")
    print(
        f"H not maximal="
        f"{len(states) - exact}"
    )

    if counterexamples:
        print()
        print("counterexamples:")

        for item in counterexamples:
            print(
                "    "
                f"n={item[0]} "
                f"p={item[1]} "
                f"q={item[2]} "
                f"frame={item[3]} "
                f"H={item[4]} "
                f"v2={item[5]} "
                f"larger={item[6]}"
            )

    print()
    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 5
#
# Try quotient congruences.
# =============================================================================

def test_quotient_congruence(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 5: QUOTIENT CONGRUENCE FINGERPRINT")
    print("=" * 90)

    failures = 0

    moduli = [2, 4, 8, 3, 5, 7, 9, 11, 13]

    signatures: dict[
        tuple[str, int, tuple[int, ...]],
        set[int],
    ] = {}

    for s in states:
        value = abs(s.n + s.c)

        H = gcd(abs(s.C), value)
        Q = value // H

        signature = (
            s.frame,
            v2(H),
            tuple(Q % m for m in moduli),
        )

        if signature not in signatures:
            signatures[signature] = set()

        signatures[signature].add(H)

    ambiguous = {
        key: values
        for key, values in signatures.items()
        if len(values) > 1
    }

    print(f"signatures={len(signatures)}")
    print(f"ambiguous={len(ambiguous)}")

    if ambiguous:
        print()
        print("first ambiguous signatures:")

        for key, values in list(ambiguous.items())[:20]:
            print(
                f"    {key} -> "
                f"{sorted(values)[:10]}"
            )

    print()
    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 6
#
# Can the full divisor H be distinguished solely from n+c?
# This intentionally compares all states with identical n+c.
# =============================================================================

def test_shifted_value_collisions(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 6: SAME SHIFTED VALUE")
    print("=" * 90)

    buckets: dict[
        tuple[str, int],
        set[int],
    ] = defaultdict(set)

    for s in states:
        value = s.n + s.c

        H = gcd(abs(s.C), abs(value))

        buckets[
            (s.frame, value)
        ].add(H)

    ambiguous = {
        key: values
        for key, values in buckets.items()
        if len(values) > 1
    }

    print(f"shifted-value buckets={len(buckets)}")
    print(f"ambiguous shifted-value buckets={len(ambiguous)}")

    for key, values in list(ambiguous.items())[:20]:
        print(f"    {key} -> {sorted(values)}")

    print()

    return len(ambiguous)


# =============================================================================
# EXAMPLES
# =============================================================================

def print_examples(states: list[State]) -> None:
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

    by_n = {s.n: s for s in states}

    shown = 0

    for n in wanted:
        s = by_n.get(n)

        if s is None:
            continue

        shown += 1

        value = abs(s.n + s.c)

        d = gcd(abs(s.A), abs(s.B))
        H = gcd(abs(s.C), value)
        Q = value // H

        ds = divisors(value)

        same_v2 = [
            x
            for x in ds
            if v2(x) == v2(H)
        ]

        print()
        print(
            f"n={s.n} p={s.p} q={s.q} "
            f"frame={s.frame}"
        )
        print(f"    A={s.A} B={s.B}")
        print(f"    C={s.C}")
        print(f"    n+c={value}")
        print(f"    d=gcd(A,B)={d}")
        print(f"    H=gcd(C,n+c)={H}")
        print(f"    v2(H)={v2(H)}")
        print(f"    odd(H)={odd_part(H)}")
        print(f"    Q=(n+c)/H={Q}")
        print(
            "    Q parity="
            + ("even" if Q % 2 == 0 else "odd")
        )
        print(
            "    divisors with same v2(H)="
            f"{len(same_v2)}"
        )
        print(
            "    same-v2 divisors="
            f"{same_v2[:20]}"
        )
        print(f"    depth={s.depth}")

    print()
    print(f"examples shown={shown}")
    print()


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    print("=" * 90)
    print("EXPERIMENT 666 START")
    print("=" * 90)
    print()

    primes = sieve_primes(PRIME_LIMIT)

    print(f"prime limit={PRIME_LIMIT}")
    print(f"odd primes={len(primes)}")

    states = build_states(primes)

    print(f"semiprimes={len(states)}")
    print()

    total_failures = 0

    total_failures += test_baseline(states)
    total_failures += test_exact_canonical_gcd(states)
    total_failures += test_divisor_position(states)
    total_failures += test_quotient_signature(states)
    total_failures += test_maximal_divisor_property(states)
    total_failures += test_quotient_congruence(states)

    # This is diagnostic, not a theorem-failure count.
    shifted_collisions = test_shifted_value_collisions(states)

    print_examples(states)

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)
    print()

    print("Experiment 665 established the exact identity:")
    print()
    print("    H = gcd(C,n+c) = gcd(A,B).")
    print()

    print("Experiment 666 asks a different question:")
    print()
    print(
        "    Can H be recognized as a distinguished divisor"
    )
    print(
        "    of the shifted integer n+c?"
    )
    print()

    print("The candidate information layers are:")
    print()
    print("    1. v2(H)")
    print()
    print("    2. quotient Q = (n+c)/H")
    print()
    print("    3. parity of Q")
    print()
    print("    4. Q modulo several small integers")
    print()
    print(
        "    5. maximality of H among divisors"
    )
    print(
        "       sharing the same 2-adic valuation"
    )
    print()

    print("The strongest possible result would be:")
    print()
    print(
        "    H = distinguished_divisor(n+c, frame)"
    )
    print()
    print(
        "which would turn the canonical gcd into a"
    )
    print(
        "factor-free divisor selection rule."
    )
    print()

    print(
        "A negative result is equally useful: it would"
    )
    print(
        "show exactly which divisor information is"
    )
    print(
        "missing from n+c alone."
    )
    print()

    print(f"TOTAL TEST FAILURES={total_failures}")
    print(
        f"SHIFTED-VALUE AMBIGUITIES={shifted_collisions}"
    )

    if total_failures == 0:
        print("STATUS=ALL TESTS PASSED")
    else:
        print("STATUS=COUNTEREXAMPLES FOUND")

    print("=" * 90)
    print("EXPERIMENT 666 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
