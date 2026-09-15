#!/usr/bin/env python3

"""
==========================================================================================
EXPERIMENT 660

EXACT INTEGER gcd COMPRESSION AND ODD-PART RESIDUAL ANALYSIS

Established:

    depth = v2(gcd(2A, 2B, n+c))

The remaining question is whether the FULL integer gcd can be
collapsed to gcd(A,B) plus the equality bit.

Candidate theorem:

    G = gcd(2A,2B,n+c)

    d = gcd(A,B)

    G = d
        if v2(A) != v2(B)

    G = 2d
        if v2(A) == v2(B)

Therefore:

    depth = v2(d) + [v2(A)=v2(B)].

This experiment tests the integer identity directly and then
separates:

    1. 2-adic information
    2. odd information

The odd part is NOT used to predict depth.

The goal is to determine exactly what information survives in
the full gcd after the 2-adic depth has been extracted.
==========================================================================================
"""

from __future__ import annotations

from collections import Counter, defaultdict
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

    if limit >= 0:
        sieve[0] = 0

    if limit >= 1:
        sieve[1] = 0

    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            start = p * p
            sieve[start : limit + 1 : p] = b"\x00" * (
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


# =============================================================================
# STATE GENERATION
# =============================================================================

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

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 1
# =============================================================================

def test_core_gcd_identity(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 1: CORE gcd(A,B) DIVISIBILITY")
    print("=" * 90)

    failures = 0

    for s in states:
        d = gcd(
            abs(s.A),
            abs(s.B),
        )

        shifted = abs(s.n + s.c)

        if shifted % d != 0:
            failures += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 2
# =============================================================================

def test_exact_full_gcd_formula(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 2: EXACT FULL INTEGER gcd FORMULA")
    print("=" * 90)

    failures = 0

    ratio_counter: Counter[int | str] = Counter()

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

        equal = (
            v2(s.A) == v2(s.B)
        )

        predicted = (
            2 * d
            if equal
            else d
        )

        if G != predicted:
            failures += 1

        if d != 0:
            if G % d == 0:
                ratio_counter[G // d] += 1
            else:
                ratio_counter["noninteger"] += 1
        else:
            ratio_counter["zero"] += 1

    print("ratio G / gcd(A,B):")

    for key in sorted(
        ratio_counter,
        key=lambda x: (isinstance(x, str), x),
    ):
        print(
            f"    {key}: {ratio_counter[key]}"
        )

    print()
    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 3
# =============================================================================

def test_one_sided_gcd_exact(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 3: ONE-SIDED gcd vs gcd(A,B)")
    print("=" * 90)

    failures = 0

    branch_counter: Counter[
        tuple[str, int]
    ] = Counter()

    for s in states:

        d = gcd(
            abs(s.A),
            abs(s.B),
        )

        gA = gcd(
            2 * abs(s.A),
            abs(s.n + s.c),
        )

        gB = gcd(
            2 * abs(s.B),
            abs(s.n + s.c),
        )

        if v2(gA) != s.depth:
            failures += 1

        if v2(gB) != s.depth:
            failures += 1

        # Determine how the one-sided gcd differs
        # from the common d.

        if d != 0 and gA % d == 0:
            ratioA = gA // d
        else:
            ratioA = -1

        if d != 0 and gB % d == 0:
            ratioB = gB // d
        else:
            ratioB = -1

        branch_counter[
            (s.frame, ratioA)
        ] += 1

        branch_counter[
            (s.frame, ratioB)
        ] += 1

    for frame in ("A", "B"):
        print(f"FRAME {frame}")

        vals = sorted(
            (
                ratio,
                count,
            )
            for (fr, ratio), count
            in branch_counter.items()
            if fr == frame
        )

        print(
            f"    ratios={vals}"
        )

    print()
    print(f"checked={2 * len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 4
# =============================================================================

def test_odd_part_exact(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 4: ODD-PART gcd IDENTITY")
    print("=" * 90)

    failures = 0

    examples: list[
        tuple[
            int,
            int,
            int,
            int,
            int,
        ]
    ] = []

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

        od = odd_part(d)
        oG = odd_part(G)

        if od != oG:
            failures += 1

            if len(examples) < 20:
                examples.append(
                    (
                        s.n,
                        d,
                        G,
                        od,
                        oG,
                    )
                )

    print(f"checked={len(states)}")
    print(f"failures={failures}")

    if examples:
        print()
        print("counterexamples:")

        for row in examples:
            print(
                f"    n={row[0]:<10d} "
                f"d={row[1]:<10d} "
                f"G={row[2]:<10d} "
                f"odd(d)={row[3]:<8d} "
                f"odd(G)={row[4]:<8d}"
            )

    print()

    return failures


# =============================================================================
# TEST 5
# =============================================================================

def test_shifted_quotient(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 5: SHIFTED n QUOTIENT AFTER EXTRACTING DEPTH")
    print("=" * 90)

    failures = 0

    parity_counter: Counter[
        tuple[str, int, int]
    ] = Counter()

    for s in states:
        shifted = abs(s.n + s.c)

        if shifted == 0:
            failures += 1
            continue

        if shifted % (1 << s.depth) != 0:
            failures += 1
            continue

        quotient = shifted >> s.depth

        if quotient % 2 == 0:
            failures += 1

        parity_counter[
            (
                s.frame,
                s.depth,
                quotient & 1,
            )
        ] += 1

    print(
        "quotient parity states="
        f"{len(parity_counter)}"
    )

    for key, count in sorted(
        parity_counter.items()
    ):
        print(
            f"    {key}: {count}"
        )

    print()
    print(
        "NOTE:"
    )
    print(
        "    This test intentionally expects quotient parity = 1."
    )
    print(
        "    Any failure would contradict"
    )
    print(
        "        depth = v2(gcd(...))."
    )

    print()
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 6
# =============================================================================

def test_odd_part_vs_shifted_odd_part(
    states: list[State],
) -> int:

    print("=" * 90)
    print("TEST 6: ODD gcd PART VS ODD PART OF n+c")
    print("=" * 90)

    failures = 0

    distribution: dict[
        str,
        Counter[int]
    ] = {
        "A": Counter(),
        "B": Counter(),
    }

    for s in states:
        d = gcd(
            abs(s.A),
            abs(s.B),
        )

        shifted_odd = odd_part(
            s.n + s.c
        )

        d_odd = odd_part(d)

        if d_odd == 0:
            failures += 1
            continue

        if shifted_odd % d_odd != 0:
            failures += 1
            continue

        quotient = shifted_odd // d_odd

        distribution[s.frame][
            quotient
        ] += 1

    for frame in ("A", "B"):

        vals = distribution[frame]

        print(
            f"FRAME {frame}"
        )

        print(
            f"    distinct quotient values="
            f"{len(vals)}"
        )

        print(
            "    most common="
            f"{vals.most_common(15)}"
        )

    print()
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 7
# =============================================================================

def test_depth_from_d_and_equality(
    states: list[State],
) -> int:

    print("=" * 90)
    print("TEST 7: DEPTH FROM gcd(A,B) + EQUALITY BIT")
    print("=" * 90)

    failures = 0

    branch_counter: Counter[
        tuple[str, str]
    ] = Counter()

    for s in states:
        d = gcd(
            abs(s.A),
            abs(s.B),
        )

        equal = (
            v2(s.A) == v2(s.B)
        )

        predicted = (
            v2(d)
            + int(equal)
        )

        if predicted != s.depth:
            failures += 1

        branch = (
            "equal"
            if equal
            else "unequal"
        )

        branch_counter[
            (
                s.frame,
                branch,
            )
        ] += 1

    for key, count in sorted(
        branch_counter.items()
    ):
        print(
            f"    {key}: {count}"
        )

    print()
    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 8
# =============================================================================

def test_n_only_shifted_value(
    states: list[State],
) -> int:

    print("=" * 90)
    print("TEST 8: SAME n+c CAN HAVE DIFFERENT DEPTH")
    print("=" * 90)

    buckets: dict[
        tuple[str, int],
        dict[int, list[tuple[int, int, int]]]
    ] = defaultdict(
        lambda: defaultdict(list)
    )

    for s in states:
        key = (
            s.frame,
            s.n + s.c,
        )

        buckets[key][
            s.depth
        ].append(
            (
                s.n,
                s.p,
                s.q,
            )
        )

    ambiguous = 0
    shown = 0

    for key, depth_map in sorted(
        buckets.items()
    ):

        if len(depth_map) <= 1:
            continue

        ambiguous += 1

        if shown < 10:
            depths = sorted(
                depth_map
            )

            print(
                f"    frame={key[0]} "
                f"n+c={key[1]} "
                f"depths={depths}"
            )

            for depth in depths:
                row = depth_map[depth][0]

                print(
                    f"        depth={depth} "
                    f"n={row[0]} "
                    f"p={row[1]} "
                    f"q={row[2]}"
                )

            shown += 1

    print()
    print(
        f"shifted-value buckets={len(buckets)}"
    )
    print(
        f"ambiguous shifted-value buckets={ambiguous}"
    )
    print()

    # This test is diagnostic only.
    return 0


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

        G = gcd(
            gcd(
                2 * abs(s.A),
                2 * abs(s.B),
            ),
            abs(s.n + s.c),
        )

        equal = (
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
            f"B={s.B} "
            f"c={s.c}"
        )

        print(
            f"    d=gcd(A,B)={d}"
        )

        print(
            f"    v2(d)={v2(d)}"
        )

        print(
            f"    equal-v2={equal}"
        )

        print(
            f"    2^depth={1 << s.depth}"
        )

        print(
            f"    G={G}"
        )

        print(
            f"    G / d={G // d}"
        )

        print(
            f"    odd(d)={odd_part(d)}"
        )

        print(
            f"    odd(G)={odd_part(G)}"
        )

        print(
            f"    n+c={s.n+s.c}"
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
    print("EXPERIMENT 660 START")
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
        test_core_gcd_identity(states)
    )

    total_failures += (
        test_exact_full_gcd_formula(states)
    )

    total_failures += (
        test_one_sided_gcd_exact(states)
    )

    total_failures += (
        test_odd_part_exact(states)
    )

    total_failures += (
        test_shifted_quotient(states)
    )

    total_failures += (
        test_odd_part_vs_shifted_odd_part(states)
    )

    total_failures += (
        test_depth_from_d_and_equality(states)
    )

    test_n_only_shifted_value(states)

    print_examples(states)

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)

    print()

    print("The currently established exact identity is:")

    print()
    print(
        "    depth = v2(gcd(2A, 2B, n+c))"
    )

    print()

    print("The proposed stronger integer decomposition is:")

    print()
    print(
        "    d = gcd(A,B)"
    )

    print()

    print(
        "    G = d"
    )

    print(
        "        when v2(A) != v2(B)"
    )

    print()

    print(
        "    G = 2d"
    )

    print(
        "        when v2(A) = v2(B)"
    )

    print()

    print(
        "which implies:"
    )

    print()
    print(
        "    depth = v2(d) + [v2(A)=v2(B)]."
    )

    print()

    print(
        "The odd-part test asks whether:"
    )

    print()
    print(
        "    odd(G) = odd(gcd(A,B))"
    )

    print()

    print(
        "If this passes, then the complete integer gcd"
    )
    print(
        "contains no new odd information after the"
    )
    print(
        "factor-side gcd has been identified."
    )

    print()

    print(
        "The remaining unresolved quantity is therefore"
    )
    print()
    print(
        "    v2(gcd(A,B))"
    )
    print()
    print(
        "or equivalently the common 2-adic depth of the"
    )
    print(
        "two residual branches."
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
            "STATUS=CORE COUNTEREXAMPLES FOUND"
        )

    print()
    print("=" * 90)
    print("EXPERIMENT 660 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
