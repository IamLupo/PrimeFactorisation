#!/usr/bin/env python3

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from math import gcd, isqrt


# =============================================================================
# CONFIGURATION
# =============================================================================

PRIME_LIMIT = 6000
INF = 10**9


# =============================================================================
# BASIC ARITHMETIC
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
# SIEVE
# =============================================================================

def sieve(limit: int) -> list[int]:
    a = bytearray(b"\x01") * (limit + 1)

    a[0] = 0
    a[1] = 0

    for p in range(2, isqrt(limit) + 1):
        if a[p]:
            start = p * p
            count = ((limit - start) // p) + 1
            a[start:limit + 1:p] = b"\x00" * count

    return [
        p
        for p in range(3, limit + 1, 2)
        if a[p]
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

    H: int
    m: int
    N: int


def make_state(p: int, q: int) -> State:
    n = p * q

    if n % 4 == 3:
        # FRAME A
        #
        # A = p-3
        # B = q+3
        # C = q+3
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
        # C = p+1
        # c = 3

        frame = "B"

        A = p + 1
        B = q - 3
        C = A
        c = 3

    N = n + c

    H = gcd(abs(C), abs(N))

    m = v2(H)

    depth = m

    if v2(A) == v2(B):
        depth += 1

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
        H=H,
        m=m,
        N=N,
    )


def build_states(primes: list[int]) -> list[State]:
    out = []

    for i, p in enumerate(primes):
        for q in primes[i:]:
            out.append(make_state(p, q))

    return out


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
            v2(s.N),
        )

        if predicted != s.depth:
            failures += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 1
#
# H = gcd(C,N)
# m = v2(H)
# Compare m with the 2-adic threshold structure of N.
# =============================================================================

def test_2adic_threshold(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 1: 2-ADIC H THRESHOLD")
    print("=" * 90)

    failures = 0

    relation = Counter()

    for s in states:
        w = v2(s.N)
        m = s.m

        relation[
            (
                s.frame,
                m,
                w - m,
            )
        ] += 1

        if m > w:
            failures += 1

    print("distribution of w-m:")

    for key, count in sorted(relation.items()):
        print(f"    {key}: {count}")

    print()
    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 2
#
# Strip the odd part from N and ask whether the canonical residual
# C has exactly the same 2-adic intersection.
# =============================================================================

def test_odd_part_removed(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 2: ODD-PART REMOVAL")
    print("=" * 90)

    failures = 0

    for s in states:
        N2 = 1 << v2(s.N)
        C2 = 1 << v2(s.C)

        expected = 1 << min(
            v2(s.C),
            v2(s.N),
        )

        actual = gcd(C2, N2)

        if actual != expected:
            failures += 1

    print(
        "Testing:"
    )
    print(
        "    gcd(2^v2(C), 2^v2(N))"
        " = 2^v2(H)"
    )
    print()

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 3
#
# The genuinely new question:
#
# For every t <= m, the canonical branch survives:
#
#      C == anchor (mod 2^t)
#      N == 0     (mod 2^t)
#
# At t=m+1, one of those should fail.
#
# We locate which channel is responsible.
# =============================================================================

def test_2adic_terminal_channel(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 3: 2-ADIC TERMINAL CHANNEL")
    print("=" * 90)

    failures = 0
    distribution = Counter()

    for s in states:
        w = v2(s.N)
        tc = v2(s.C)

        m = min(tc, w)

        if m != s.m:
            failures += 1
            continue

        if tc < w:
            kind = "C-first"
        elif w < tc:
            kind = "N-first"
        else:
            kind = "simultaneous"

        distribution[
            (s.frame, kind, m)
        ] += 1

    for key, count in sorted(distribution.items()):
        print(f"    {key}: {count}")

    print()
    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 4
#
# Does the odd part of H affect the terminal channel?
# =============================================================================

def test_odd_part_independence(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 4: ODD PART INDEPENDENCE")
    print("=" * 90)

    failures = 0

    buckets: dict[
        tuple[str, int],
        set[tuple[int, int, int]],
    ] = defaultdict(set)

    for s in states:
        h_odd = odd_part(s.H)

        buckets[
            (s.frame, s.m)
        ].add(
            (
                h_odd,
                v2(s.C),
                v2(s.N),
            )
        )

    # The actual theorem we can establish is that depth depends
    # only on the two valuations. The odd part must therefore
    # not alter the depth for fixed valuation pair.

    for s in states:
        predicted = min(
            v2(s.C),
            v2(s.N),
        )

        if predicted != s.m:
            failures += 1

    print("fixed valuation pair -> fixed m")

    for frame in ("A", "B"):
        pairs = Counter(
            (
                s.frame,
                v2(s.C),
                v2(s.N),
            )
            for s in states
            if s.frame == frame
        )

        print(
            f"    FRAME {frame}: "
            f"valuation-pairs={len(pairs)}"
        )

    print()

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 5
#
# Can m be described as the largest t for which BOTH:
#
#     C == 0 (after subtracting the frame anchor)
#     N == 0
#
# survive modulo 2^t?
#
# This is essentially a valuation theorem, but encoded as
# a pure threshold set.
# =============================================================================

def test_nested_threshold_sets(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 5: NESTED 2-ADIC THRESHOLD SETS")
    print("=" * 90)

    failures = 0

    for s in states:
        m = s.m

        for t in range(1, m + 1):
            if s.C % (1 << t) != 0:
                failures += 1
                break

            if s.N % (1 << t) != 0:
                failures += 1
                break

        if m + 1 <= 20:
            c_ok = s.C % (1 << (m + 1)) == 0
            n_ok = s.N % (1 << (m + 1)) == 0

            # At least one channel must fail.
            if c_ok and n_ok:
                failures += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 6
#
# NEW STRUCTURAL TEST:
#
# Does the actual H correspond to the intersection of the two
# principal 2-adic divisor chains?
#
#       1 | 2 | 4 | 8 | ...
#
# C contributes one chain.
# N contributes one chain.
#
# H_2 is their last common element.
# =============================================================================

def test_2adic_divisor_chain(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 6: 2-ADIC DIVISOR CHAIN INTERSECTION")
    print("=" * 90)

    failures = 0

    for s in states:
        chain_c = {
            1 << t
            for t in range(v2(s.C) + 1)
        }

        chain_n = {
            1 << t
            for t in range(v2(s.N) + 1)
        }

        intersection = chain_c & chain_n

        if not intersection:
            failures += 1
            continue

        largest = max(intersection)

        expected = 1 << s.m

        if largest != expected:
            failures += 1

    print("Testing:")
    print("    largest common power-of-two divisor")
    print("    of C and n+c")
    print()

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 7
#
# Compare the canonical integer H with its 2-adic projection.
# =============================================================================

def test_projection(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 7: INTEGER gcd -> 2-ADIC PROJECTION")
    print("=" * 90)

    failures = 0

    for s in states:
        H2 = 1 << v2(s.H)

        if H2 != gcd(
            1 << v2(s.C),
            1 << v2(s.N),
        ):
            failures += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 8
#
# Signature:
#
#     (frame, v2(C), v2(N))
#
# Does it determine both m and depth?
# =============================================================================

def test_signature(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 8: VALUATION-ONLY SIGNATURE")
    print("=" * 90)

    buckets: dict[
        tuple[str, int, int],
        set[tuple[int, int]],
    ] = defaultdict(set)

    for s in states:
        key = (
            s.frame,
            v2(s.C),
            v2(s.N),
        )

        buckets[key].add(
            (
                s.m,
                s.depth,
            )
        )

    ambiguous = {
        key: values
        for key, values in buckets.items()
        if len(values) > 1
    }

    print(f"signatures={len(buckets)}")
    print(f"ambiguous={len(ambiguous)}")

    if ambiguous:
        print("first ambiguous signatures:")

        for key, values in list(ambiguous.items())[:20]:
            print(
                f"    {key} -> "
                f"{sorted(values)}"
            )

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

        vc = v2(s.C)
        vn = v2(s.N)
        common = min(vc, vn)

        print()
        print(
            f"n={s.n} "
            f"p={s.p} "
            f"q={s.q} "
            f"frame={s.frame}"
        )

        print(
            f"    C={s.C}"
        )

        print(
            f"    n+c={s.N}"
        )

        print(
            f"    H=gcd(C,n+c)={s.H}"
        )

        print(
            f"    v2(C)={vc}"
        )

        print(
            f"    v2(n+c)={vn}"
        )

        print(
            f"    m=v2(H)={s.m}"
        )

        print(
            f"    common 2-adic level={common}"
        )

        print(
            f"    odd(H)={odd_part(s.H)}"
        )

        print(
            f"    depth={s.depth}"
        )

        print(
            "    terminal channel="
            + (
                "C"
                if vc < vn
                else "N"
                if vn < vc
                else "simultaneous"
            )
        )

    print()
    print(f"examples shown={shown}")
    print()


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    print("=" * 90)
    print("EXPERIMENT 667 START")
    print("=" * 90)
    print()

    primes = sieve(PRIME_LIMIT)

    print(f"prime limit={PRIME_LIMIT}")
    print(f"odd primes={len(primes)}")

    states = build_states(primes)

    print(f"semiprimes={len(states)}")
    print()

    total_failures = 0

    total_failures += test_baseline(states)
    total_failures += test_2adic_threshold(states)
    total_failures += test_odd_part_removed(states)
    total_failures += test_2adic_terminal_channel(states)
    total_failures += test_odd_part_independence(states)
    total_failures += test_nested_threshold_sets(states)
    total_failures += test_2adic_divisor_chain(states)
    total_failures += test_projection(states)

    ambiguity = test_signature(states)

    print_examples(states)

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)
    print()

    print("Experiment 665 established:")
    print()
    print("    H = gcd(C,n+c) = gcd(A,B).")
    print()

    print("Experiment 666 showed:")
    print()
    print("    H is generally NOT uniquely characterized")
    print("    by ordinary divisor-lattice properties of n+c.")
    print()

    print("Experiment 667 therefore isolates the 2-adic part:")
    print()
    print("    m = v2(H)")
    print()
    print("Since:")
    print()
    print("    H = gcd(C,n+c)")
    print()
    print("we have:")
    print()
    print("    m = min(v2(C), v2(n+c)).")
    print()
    print("The experiment tests whether this is best understood")
    print("as the intersection of two nested 2-adic divisor chains:")
    print()
    print("    C:")
    print("        1 | 2 | 4 | 8 | ... | 2^v2(C)")
    print()
    print("    n+c:")
    print("        1 | 2 | 4 | 8 | ... | 2^v2(n+c)")
    print()
    print("The largest common power of two is:")
    print()
    print("    2^m")
    print()
    print("and therefore:")
    print()
    print("    m = min(v2(C), v2(n+c)).")
    print()
    print("The remaining question is whether this 2-adic")
    print("intersection has a deeper interpretation in terms")
    print("of the original factor branches.")
    print()

    print(f"TOTAL TEST FAILURES={total_failures}")
    print(f"SIGNATURE AMBIGUITIES={ambiguity}")

    if total_failures == 0 and ambiguity == 0:
        print("STATUS=ALL TESTS PASSED")
    else:
        print("STATUS=COUNTEREXAMPLES FOUND")

    print("=" * 90)
    print("EXPERIMENT 667 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
