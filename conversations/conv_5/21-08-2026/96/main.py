#!/usr/bin/env python3

"""
==========================================================================================
EXPERIMENT 659
==========================================================================================

2-ADIC gcd BRIDGE / DISTINGUISHED DIVISOR ANALYSIS

Previous exact theorem:

    depth =
        v2(gcd(2A, 2B, n+c))

equivalently:

    depth =
        v2(gcd(2A, n+c))
      =
        v2(gcd(2B, n+c))

where:

FRAME A:
    A = p - 3
    B = q + 3
    c = 9

FRAME B:
    A = p + 1
    B = q - 3
    c = 3

This experiment asks:

    1. Does the symmetric gcd have a canonical 2-adic form?
    2. Is its odd part irrelevant to depth?
    3. Can the relevant gcd be reduced to:
           2^depth
       times an odd divisor of n+c?
    4. Is there a simpler exact relation between:
           gcd(2A, n+c)
       and:
           gcd(2B, n+c)?
    5. Can the depth be recovered from a purely divisor-level
       object once the branch-dependent gcd is factored?

The important distinction:

    FULL GCD:
        gcd(2A, n+c)

    2-ADIC PART:
        2^depth

The goal is to determine whether anything besides the
power of two survives in a structurally meaningful way.
==========================================================================================
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from math import gcd, isqrt


INF = 10**9


# =============================================================================
# BASIC NUMBER THEORY
# =============================================================================

def v2(x: int) -> int:
    """Return v2(|x|). Treat 0 as infinity."""
    x = abs(x)
    if x == 0:
        return INF
    return (x & -x).bit_length() - 1


def odd_part(x: int) -> int:
    """Return odd part of |x|."""
    x = abs(x)
    if x == 0:
        return 0
    return x >> v2(x)


def sieve_primes(limit: int) -> list[int]:
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"

    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            sieve[p * p : limit + 1 : p] = b"\x00" * (
                ((limit - p * p) // p) + 1
            )

    return [p for p in range(3, limit + 1, 2) if sieve[p]]


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
# FRAME / RESIDUAL CONSTRUCTION
# =============================================================================

def make_state(p: int, q: int) -> State:
    n = p * q

    # The established frame convention.
    if n % 4 == 3:
        # FRAME A
        A = p - 3
        B = q + 3
        c = 9
        frame = "A"
    else:
        # FRAME B
        A = p + 1
        B = q - 3
        c = 3
        frame = "B"

    tp = v2(A)
    w = v2(n + c)

    depth = min(tp + 1, w)

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


def build_states(prime_limit: int) -> list[State]:
    primes = sieve_primes(prime_limit)

    states: list[State] = []

    for i, p in enumerate(primes):
        for q in primes[i:]:
            states.append(make_state(p, q))

    return states


# =============================================================================
# TEST 0
# =============================================================================

def test_baseline(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 0: PREVIOUS EXACT DEPTH LAW")
    print("=" * 90)

    failures = 0

    for s in states:
        tp = v2(s.A)
        w = v2(s.n + s.c)

        predicted = min(tp + 1, w)

        if predicted != s.depth:
            failures += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 1
# =============================================================================

def test_one_sided_gcd(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 1: ONE-SIDED INTEGER gcd BRIDGE")
    print("=" * 90)

    failures = 0

    for s in states:
        gA = gcd(2 * abs(s.A), abs(s.n + s.c))
        gB = gcd(2 * abs(s.B), abs(s.n + s.c))

        if v2(gA) != s.depth:
            failures += 1

        if v2(gB) != s.depth:
            failures += 1

    print(f"checked={2 * len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 2
# =============================================================================

def test_symmetric_gcd(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 2: SYMMETRIC INTEGER gcd")
    print("=" * 90)

    failures = 0

    for s in states:
        g = gcd(
            gcd(2 * abs(s.A), 2 * abs(s.B)),
            abs(s.n + s.c),
        )

        if v2(g) != s.depth:
            failures += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 3
# =============================================================================

def test_odd_part_irrelevance(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 3: ODD PART IRRELEVANCE")
    print("=" * 90)

    failures = 0

    distributions: dict[str, Counter[int]] = {
        "A": Counter(),
        "B": Counter(),
    }

    for s in states:
        g = gcd(
            gcd(2 * abs(s.A), 2 * abs(s.B)),
            abs(s.n + s.c),
        )

        power = 1 << s.depth

        if g % power != 0:
            failures += 1
            continue

        odd = g // power

        if odd % 2 == 0:
            failures += 1
            continue

        distributions[s.frame][odd] += 1

    for frame in ("A", "B"):
        vals = distributions[frame]

        print(f"FRAME {frame}")
        print(f"    distinct odd gcd parts={len(vals)}")

        if vals:
            print(
                "    most common:",
                vals.most_common(15),
            )

    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 4
# =============================================================================

def test_gcd_ratio_structure(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 4: RELATION BETWEEN gcd(2A,n+c) AND gcd(2B,n+c)")
    print("=" * 90)

    failures = 0

    ratio_counter: Counter[str] = Counter()

    examples: dict[str, list[tuple[int, int, int, int, int]]] = defaultdict(list)

    for s in states:
        gA = gcd(2 * abs(s.A), abs(s.n + s.c))
        gB = gcd(2 * abs(s.B), abs(s.n + s.c))

        if gA == 0 or gB == 0:
            failures += 1
            continue

        if v2(gA) != v2(gB):
            failures += 1
            continue

        if gA == gB:
            key = "equal"
        elif gA % gB == 0:
            key = "gA/gB-integer"
        elif gB % gA == 0:
            key = "gB/gA-integer"
        else:
            key = "nondividing"

        ratio_counter[key] += 1

        if len(examples[key]) < 10:
            examples[key].append(
                (s.n, s.p, s.q, gA, gB)
            )

    for key, count in ratio_counter.items():
        print(f"    {key}: {count}")

    for key, rows in examples.items():
        if key != "equal":
            print()
            print(f"    examples [{key}]")
            for row in rows:
                print(
                    "        n=%-10d p=%-6d q=%-6d gcdA=%-8d gcdB=%-8d"
                    % row
                )

    print()
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 5
# =============================================================================

def test_normalized_gcd(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 5: NORMALIZED gcd BY ITS 2-ADIC DEPTH")
    print("=" * 90)

    failures = 0

    counter: Counter[tuple[str, int, int]] = Counter()

    for s in states:
        g = gcd(
            gcd(2 * abs(s.A), 2 * abs(s.B)),
            abs(s.n + s.c),
        )

        d = 1 << s.depth

        if g % d != 0:
            failures += 1
            continue

        normalized = g // d

        if normalized % 2 == 0:
            failures += 1
            continue

        counter[(s.frame, s.depth, normalized)] += 1

    for frame in ("A", "B"):
        print(f"FRAME {frame}")

        frame_data = [
            (depth, odd)
            for (fr, depth, odd), count in counter.items()
            if fr == frame
        ]

        depths = sorted(set(depth for depth, _ in frame_data))

        for depth in depths:
            odds = sorted(
                odd
                for d, odd in frame_data
                if d == depth
            )

            print(
                f"    depth={depth:<2} "
                f"distinct_normalized_gcd={len(odds):<4} "
                f"sample={odds[:12]}"
            )

    print()
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 6
# =============================================================================

def test_divisor_lattice_relation(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 6: 2-ADIC DIVISOR-LATTICE POSITION")
    print("=" * 90)

    failures = 0

    relation_counter: Counter[tuple[str, int]] = Counter()

    for s in states:
        shifted = abs(s.n + s.c)
        g = gcd(
            gcd(2 * abs(s.A), 2 * abs(s.B)),
            shifted,
        )

        quotient = shifted // (1 << s.depth)

        if shifted % (1 << s.depth) != 0:
            failures += 1
            continue

        if quotient % 2 == 0:
            failures += 1
            continue

        relation_counter[(s.frame, v2(quotient))] += 1

    for frame in ("A", "B"):
        print(f"FRAME {frame}")

        vals = sorted(
            (k[1], count)
            for k, count in relation_counter.items()
            if k[0] == frame
        )

        print(f"    quotient valuations={vals}")

    print()
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 7
# =============================================================================

def test_factor_free_candidate(states: list[State]) -> int:
    """
    Search a SMALL deterministic family of n-only gcd expressions.

    Important:
        gcd(n+c, n+r)
        = gcd(n+c, r-c),

    so this family can only see fixed constants and therefore
    cannot encode an unbounded branch valuation.

    We test it explicitly rather than assuming this.
    """

    print("=" * 90)
    print("TEST 7: SMALL n-ONLY gcd FAMILY")
    print("=" * 90)

    candidates = range(-32, 33)

    exact: list[int] = []

    for r in candidates:
        ok = True

        for s in states:
            value = v2(
                gcd(
                    abs(s.n + s.c),
                    abs(s.n + r),
                )
            )

            if value != s.depth:
                ok = False
                break

        if ok:
            exact.append(r)

    print(f"tested shifts={len(list(candidates))}")
    print(f"exact shifts={exact}")

    # This test is diagnostic; failure is expected.
    print()
    return 0


# =============================================================================
# TEST 8
# =============================================================================

def test_branch_information_needed(states: list[State]) -> int:
    """
    Information-loss test.

    Group by:
        (frame, w=v2(n+c))

    and inspect whether multiple depths occur.

    If multiple depths occur, w alone cannot reconstruct depth.
    """

    print("=" * 90)
    print("TEST 8: BRANCH INFORMATION IS STILL REQUIRED")
    print("=" * 90)

    buckets: dict[tuple[str, int], set[int]] = defaultdict(set)

    examples: dict[tuple[str, int], list[tuple[int, int, int, int]]] = defaultdict(list)

    for s in states:
        w = v2(s.n + s.c)
        key = (s.frame, w)

        buckets[key].add(s.depth)

        if len(examples[key]) < 4:
            examples[key].append(
                (s.n, s.p, s.q, s.depth)
            )

    ambiguous = 0

    for key, depths in sorted(buckets.items()):
        if len(depths) > 1:
            ambiguous += 1

    print(f"w-buckets={len(buckets)}")
    print(f"ambiguous w-buckets={ambiguous}")

    shown = 0

    for key, depths in sorted(buckets.items()):
        if len(depths) <= 1:
            continue

        print(
            f"    frame={key[0]} w={key[1]} depths={sorted(depths)}"
        )

        for row in examples[key]:
            print(
                f"        n={row[0]:<10d} "
                f"p={row[1]:<6d} "
                f"q={row[2]:<6d} "
                f"depth={row[3]}"
            )

        shown += 1

        if shown >= 10:
            break

    print()

    # Diagnostic only. It is expected that ambiguity exists.
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

    by_n = {s.n: s for s in states}

    shown = 0

    for n in wanted:
        if n not in by_n:
            continue

        s = by_n[n]

        gA = gcd(2 * abs(s.A), abs(s.n + s.c))
        gB = gcd(2 * abs(s.B), abs(s.n + s.c))
        gS = gcd(
            gcd(2 * abs(s.A), 2 * abs(s.B)),
            abs(s.n + s.c),
        )

        print(
            f"\nn={s.n} p={s.p} q={s.q} frame={s.frame}"
        )
        print(
            f"    A={s.A} B={s.B} c={s.c}"
        )
        print(
            f"    n+c={s.n+s.c}"
        )
        print(
            f"    gcd(2A,n+c)={gA} v2={v2(gA)}"
        )
        print(
            f"    gcd(2B,n+c)={gB} v2={v2(gB)}"
        )
        print(
            f"    gcd(2A,2B,n+c)={gS} v2={v2(gS)}"
        )
        print(
            f"    odd(gcd-sym)={odd_part(gS)}"
        )
        print(
            f"    depth={s.depth}"
        )

        shown += 1

    print(f"\nexamples shown={shown}")
    print()


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    PRIME_LIMIT = 6000

    print("=" * 90)
    print("EXPERIMENT 659 START")
    print("=" * 90)
    print()
    print(f"prime limit={PRIME_LIMIT}")

    primes = sieve_primes(PRIME_LIMIT)
    print(f"odd primes={len(primes)}")

    states = build_states(PRIME_LIMIT)
    print(f"semiprimes={len(states)}")
    print()

    total_failures = 0

    total_failures += test_baseline(states)
    total_failures += test_one_sided_gcd(states)
    total_failures += test_symmetric_gcd(states)
    total_failures += test_odd_part_irrelevance(states)
    total_failures += test_gcd_ratio_structure(states)
    total_failures += test_normalized_gcd(states)
    total_failures += test_divisor_lattice_relation(states)

    # Diagnostic searches / information tests.
    test_factor_free_candidate(states)
    test_branch_information_needed(states)

    print_examples(states)

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)
    print()
    print("Exact theorem:")
    print()
    print("    depth = v2(gcd(2A, 2B, n+c))")
    print()
    print("Equivalently:")
    print()
    print("    depth = v2(gcd(2A, n+c))")
    print("          = v2(gcd(2B, n+c))")
    print()
    print("The symmetric gcd therefore has the form:")
    print()
    print("    G = 2^depth * odd_part(G)")
    print()
    print("This experiment asks whether the odd part of G")
    print("contains any structurally relevant information, or")
    print("whether the complete depth mechanism is genuinely")
    print("the 2-adic valuation alone.")
    print()
    print("The factor-free n-only gcd family is also tested.")
    print("Because gcd(n+c,n+r)=gcd(n+c,r-c), such a fixed")
    print("shift family can only produce a bounded truncation")
    print("of v2(n+c).")
    print()
    print(f"TOTAL FAILURES={total_failures}")
    print(
        "STATUS="
        + (
            "ALL CORE TESTS PASSED"
            if total_failures == 0
            else "CORE COUNTEREXAMPLES FOUND"
        )
    )
    print()
    print("=" * 90)
    print("EXPERIMENT 659 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
