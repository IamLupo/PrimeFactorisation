#!/usr/bin/env python3

from __future__ import annotations

from collections import defaultdict, Counter
from dataclasses import dataclass
from math import gcd, isqrt


# =============================================================================
# CONFIG
# =============================================================================

PRIME_LIMIT = 6000
INF = 10**9


# =============================================================================
# ARITHMETIC
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

def sieve(limit: int) -> list[int]:
    a = bytearray(b"\x01") * (limit + 1)
    a[0] = a[1] = 0

    for p in range(2, isqrt(limit) + 1):
        if a[p]:
            start = p * p
            count = ((limit - start) // p) + 1
            a[start:limit + 1:p] = b"\x00" * count

    return [p for p in range(3, limit + 1, 2) if a[p]]


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

    S: int
    D: int

    c: int

    w: int
    vc: int
    vd: int
    depth: int


def make_state(p: int, q: int) -> State:
    n = p * q
    S = p + q
    D = p - q

    if n % 4 == 3:
        # FRAME A
        #
        # A = p-3
        # B = q+3
        # C = B
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
        # C = A
        # c = 3

        frame = "B"
        A = p + 1
        B = q - 3
        C = A
        c = 3

    N = n + c

    vc = v2(C)
    w = v2(N)
    vd = v2(D)

    depth = min(vc + 1, w)

    return State(
        n=n,
        p=p,
        q=q,
        frame=frame,
        A=A,
        B=B,
        C=C,
        S=S,
        D=D,
        c=c,
        w=w,
        vc=vc,
        vd=vd,
        depth=depth,
    )


def build_states(primes: list[int]) -> list[State]:
    states = []

    for i, p in enumerate(primes):
        for q in primes[i:]:
            states.append(make_state(p, q))

    return states


# =============================================================================
# TEST 0
# =============================================================================

def test_baseline(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 0: BASELINE CANONICAL LAW")
    print("=" * 90)

    failures = 0

    for s in states:
        predicted = min(s.vc + 1, s.w)

        if predicted != s.depth:
            failures += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# =============================================================================
# TEST 1
#
# Can v2(C) be recovered from v2(S)?
# =============================================================================

def test_sum_valuation(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 1: C VALUATION VS SUM VALUATION")
    print("=" * 90)

    buckets: dict[
        tuple[str, int, int],
        set[int],
    ] = defaultdict(set)

    for s in states:
        buckets[
            (s.frame, v2(s.S), s.w)
        ].add(s.vc)

    ambiguous = {
        k: v
        for k, v in buckets.items()
        if len(v) > 1
    }

    print(f"signatures={len(buckets)}")
    print(f"ambiguous={len(ambiguous)}")

    if ambiguous:
        print("first ambiguous signatures:")

        for key, values in list(ambiguous.items())[:20]:
            print(f"    {key} -> {sorted(values)}")

    print()

    return len(ambiguous)


# =============================================================================
# TEST 2
#
# Can v2(C) be recovered from v2(p-q)?
# =============================================================================

def test_difference_valuation(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 2: C VALUATION VS DIFFERENCE VALUATION")
    print("=" * 90)

    buckets: dict[
        tuple[str, int, int],
        set[int],
    ] = defaultdict(set)

    for s in states:
        buckets[
            (s.frame, s.vd, s.w)
        ].add(s.vc)

    ambiguous = {
        k: v
        for k, v in buckets.items()
        if len(v) > 1
    }

    print(f"signatures={len(buckets)}")
    print(f"ambiguous={len(ambiguous)}")

    if ambiguous:
        print("first ambiguous signatures:")

        for key, values in list(ambiguous.items())[:20]:
            print(f"    {key} -> {sorted(values)}")

    print()

    return len(ambiguous)


# =============================================================================
# TEST 3
#
# Can the pair (v2(S), v2(D)) recover v2(C)?
# =============================================================================

def test_sum_difference_pair(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 3: (v2(S), v2(D)) -> v2(C)")
    print("=" * 90)

    buckets: dict[
        tuple[str, int, int],
        set[int],
    ] = defaultdict(set)

    for s in states:
        buckets[
            (s.frame, v2(s.S), v2(s.D))
        ].add(s.vc)

    ambiguous = {
        k: v
        for k, v in buckets.items()
        if len(v) > 1
    }

    print(f"signatures={len(buckets)}")
    print(f"ambiguous={len(ambiguous)}")

    if ambiguous:
        print("first ambiguous signatures:")

        for key, values in list(ambiguous.items())[:20]:
            print(f"    {key} -> {sorted(values)}")

    print()

    return len(ambiguous)


# =============================================================================
# TEST 4
#
# Can only the minimum/max relation recover C valuation?
#
# For odd p,q:
#
#     v2(S) and v2(D)
#
# encode the equal/unequal branch.
# =============================================================================

def test_minmax_signature(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 4: MIN/MAX 2-ADIC SIGNATURE")
    print("=" * 90)

    buckets: dict[
        tuple[str, int, int, int],
        set[int],
    ] = defaultdict(set)

    for s in states:
        a = v2(s.p - (3 if s.frame == "A" else -1))

        # second residual valuation
        b = (
            v2(s.q + 3)
            if s.frame == "A"
            else
            v2(s.q - 3)
        )

        lo = min(a, b)
        hi = max(a, b)

        buckets[
            (s.frame, lo, hi, s.w)
        ].add(s.vc)

    ambiguous = {
        k: v
        for k, v in buckets.items()
        if len(v) > 1
    }

    print(f"signatures={len(buckets)}")
    print(f"ambiguous={len(ambiguous)}")

    if ambiguous:
        print("first ambiguous signatures:")

        for key, values in list(ambiguous.items())[:20]:
            print(f"    {key} -> {sorted(values)}")

    print()

    return len(ambiguous)


# =============================================================================
# TEST 5
#
# Pure n-side signatures versus C valuation.
#
# This should reproduce the information loss:
#
#     (frame, v2(n+c))
#
# does NOT determine v2(C).
# =============================================================================

def test_n_only(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 5: N-ONLY C-VALUATION LOSS")
    print("=" * 90)

    buckets: dict[
        tuple[str, int],
        set[int],
    ] = defaultdict(set)

    for s in states:
        buckets[
            (s.frame, s.w)
        ].add(s.vc)

    ambiguous = {
        k: v
        for k, v in buckets.items()
        if len(v) > 1
    }

    print(f"signatures={len(buckets)}")
    print(f"ambiguous={len(ambiguous)}")

    if ambiguous:
        print("first collisions:")

        for key, values in list(ambiguous.items())[:20]:
            print(f"    {key} -> {sorted(values)}")

    print()

    return len(ambiguous)


# =============================================================================
# TEST 6
#
# Minimal integer invariant search.
#
# Candidates:
#
#   S
#   D
#   S mod powers of two
#   D mod powers of two
#   n+c
#   combinations of the above
#
# Search signatures built from small 2-adic observables.
# =============================================================================

def test_small_signature_family(states: list[State]) -> None:
    print("=" * 90)
    print("TEST 6: SMALL INTEGER SIGNATURE FAMILY")
    print("=" * 90)

    candidates = {
        "v2S": lambda s: v2(s.S),
        "v2D": lambda s: v2(s.D),
        "v2N": lambda s: s.w,
        "v2S+v2N": lambda s: (v2(s.S), s.w),
        "v2D+v2N": lambda s: (v2(s.D), s.w),
        "v2S+v2D": lambda s: (v2(s.S), v2(s.D)),
        "v2S+v2D+v2N": lambda s: (
            v2(s.S),
            v2(s.D),
            s.w,
        ),
    }

    for name, fn in candidates.items():
        buckets: dict[tuple, set[int]] = defaultdict(set)

        for s in states:
            key = (s.frame, fn(s))
            buckets[key].add(s.vc)

        ambiguous = sum(
            1
            for values in buckets.values()
            if len(values) > 1
        )

        print(
            f"    {name:<18}"
            f" signatures={len(buckets):5d}"
            f" ambiguous={ambiguous:5d}"
        )

    print()


# =============================================================================
# TEST 7
#
# Exact relationship between C valuation and residual valuation.
# =============================================================================

def test_c_is_one_residual(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 7: CANONICAL C IS ONE RESIDUAL")
    print("=" * 90)

    failures = 0

    for s in states:
        expected = (
            v2(s.B)
            if s.frame == "A"
            else
            v2(s.A)
        )

        if s.vc != expected:
            failures += 1

    print(
        "FRAME A: C=B=q+3"
    )
    print(
        "FRAME B: C=A=p+1"
    )
    print()

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


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
        93,
        141,
        183,
        213,
        485879,
        5579767,
    ]

    by_n = {s.n: s for s in states}

    for n in wanted:
        s = by_n.get(n)

        if s is None:
            continue

        print()
        print(
            f"n={s.n} p={s.p} q={s.q} frame={s.frame}"
        )

        print(
            f"    S=p+q={s.S}"
        )

        print(
            f"    D=p-q={s.D}"
        )

        print(
            f"    C={s.C}"
        )

        print(
            f"    v2(S)={v2(s.S)}"
        )

        print(
            f"    v2(D)={v2(s.D)}"
        )

        print(
            f"    v2(C)={s.vc}"
        )

        print(
            f"    v2(n+c)={s.w}"
        )

        print(
            f"    depth={s.depth}"
        )


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:
    print("=" * 90)
    print("EXPERIMENT 668 START")
    print("=" * 90)

    primes = sieve(PRIME_LIMIT)

    print()
    print(f"prime limit={PRIME_LIMIT}")
    print(f"odd primes={len(primes)}")

    states = build_states(primes)

    print(f"semiprimes={len(states)}")
    print()

    failures = 0

    failures += test_baseline(states)
    failures += test_sum_valuation(states)
    failures += test_difference_valuation(states)
    failures += test_sum_difference_pair(states)

    test_minmax_signature(states)

    n_only_ambiguous = test_n_only(states)

    test_small_signature_family(states)

    failures += test_c_is_one_residual(states)

    print_examples(states)

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)
    print()

    print("Experiment 667 established:")
    print()
    print("    H = gcd(C,n+c)")
    print()
    print("and therefore:")
    print()
    print("    v2(H) = min(v2(C), v2(n+c)).")
    print()
    print("Experiment 668 asks what information is actually")
    print("required to determine v2(C).")
    print()
    print("The key candidates are:")
    print()
    print("    v2(p+q)")
    print("    v2(p-q)")
    print("    (v2(p+q), v2(p-q))")
    print()
    print("If one of these becomes exact, then the canonical")
    print("residual valuation can be reconstructed from symmetric")
    print("factor invariants instead of from the selected factor.")
    print()
    print("If all of them remain ambiguous, that is a much sharper")
    print("information-loss result:")
    print()
    print("    n+c -> min(v2(C),v2(n+c))")
    print()
    print("but:")
    print()
    print("    v2(C)")
    print()
    print("still contains genuinely branch-specific information.")
    print()
    print(f"TOTAL TEST FAILURES={failures}")
    print(f"N-ONLY AMBIGUITY BUCKETS={n_only_ambiguous}")
    print()

    if failures == 0:
        print("STATUS=CORE TESTS PASSED")
    else:
        print("STATUS=COUNTEREXAMPLES FOUND")

    print("=" * 90)
    print("EXPERIMENT 668 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
