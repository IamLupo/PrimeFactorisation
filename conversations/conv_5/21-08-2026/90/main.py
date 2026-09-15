#!/usr/bin/env python3

"""
==========================================================================================
EXPERIMENT 652
==========================================================================================

BRANCH-BREAK LEVEL FROM n AND FACTOR ORDER

EXACT LAW FROM EXPERIMENT 651
--------------------------------

    depth = min(tp + 1, w)

where:

    tp = v2(p - p0)
    w  = v2(n + c)

FRAME A:
    p0 =  3
    c  =  9

FRAME B:
    p0 = -1
    c  =  3


Experiment 651 also established the equivalent threshold law:

    depth >= d
        iff
    p == p0 (mod 2^(d-1))
    AND
    n == -c (mod 2^d).

The remaining question is whether the branch-break level

    tp + 1

can be reconstructed from simple INTEGER information about
the factor pair without explicitly evaluating v2(p-p0).

This experiment tests increasingly compressed descriptions:

    1. p mod 2^k
    2. q mod 2^k
    3. p+q mod 2^k
    4. p-q mod 2^k
    5. n and one factor modulo 2^k
    6. the first differing bit of p from p0
    7. a direct threshold predicate

The most important test is:

    depth >= d

    <=> threshold_d(frame, n, p)

and whether threshold_d can be rewritten using:

    n,
    p,
    q,
    p+q,
    p-q

without an explicit valuation call.

No arbitrary offset search is performed.
"""


from __future__ import annotations

import math
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass


INF = 10**9


# ==========================================================================================
# BASIC UTILITIES
# ==========================================================================================

def v2(x: int) -> int:
    if x == 0:
        return INF

    x = abs(x)
    return (x & -x).bit_length() - 1


def sieve(limit: int) -> list[int]:
    if limit < 2:
        return []

    s = bytearray(b"\x01") * (limit + 1)
    s[:2] = b"\x00\x00"

    for p in range(2, math.isqrt(limit) + 1):
        if s[p]:
            start = p * p
            s[start::p] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    return [p for p in range(3, limit + 1, 2) if s[p]]


# ==========================================================================================
# STATE
# ==========================================================================================

@dataclass(frozen=True)
class State:
    p: int
    q: int
    n: int
    frame: str

    p0: int
    c: int

    tp: int
    tq: int
    w: int

    depth: int


# ==========================================================================================
# FRAME
# ==========================================================================================

def frame_parameters(p: int, q: int):
    n = p * q

    if n % 4 == 3:
        return "A", 3, 9

    return "B", -1, 3


def exact_depth(p: int, q: int, frame: str) -> int:
    if frame == "A":
        u = q - p + 6
        v = p + q
    else:
        u = 3 * p - q + 6
        v = 3 * p + q

    return min(v2(u), v2(v))


def make_states(primes: list[int]) -> list[State]:
    states: list[State] = []

    for i, p in enumerate(primes):
        for q in primes[i:]:
            n = p * q

            frame, p0, c = frame_parameters(p, q)

            if frame == "A":
                q0 = -3
            else:
                q0 = 3

            tp = v2(p - p0)
            tq = v2(q - q0)
            w = v2(n + c)

            depth = exact_depth(p, q, frame)

            states.append(
                State(
                    p=p,
                    q=q,
                    n=n,
                    frame=frame,
                    p0=p0,
                    c=c,
                    tp=tp,
                    tq=tq,
                    w=w,
                    depth=depth,
                )
            )

    return states


# ==========================================================================================
# TEST 0
# ==========================================================================================

def test_baseline(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 0: EXPERIMENT 651 BASELINE")
    print("=" * 90)

    failures = 0

    for s in states:
        predicted = min(s.tp + 1, s.w)

        if predicted != s.depth:
            if failures < 20:
                print(
                    f"mismatch n={s.n} p={s.p} q={s.q} "
                    f"frame={s.frame} tp={s.tp} w={s.w} "
                    f"actual={s.depth} predicted={predicted}"
                )

            failures += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# ==========================================================================================
# THRESHOLD
# ==========================================================================================

def threshold(state: State, d: int) -> bool:
    if d < 2:
        return True

    return (
        (state.p - state.p0) % (1 << (d - 1)) == 0
        and
        (state.n + state.c) % (1 << d) == 0
    )


# ==========================================================================================
# TEST 1
# ==========================================================================================

def test_threshold_theorem(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 1: THRESHOLD THEOREM")
    print("=" * 90)

    failures = 0

    for s in states:
        for d in range(2, s.depth + 2):
            lhs = s.depth >= d
            rhs = threshold(s, d)

            if lhs != rhs:
                if failures < 20:
                    print(
                        f"mismatch n={s.n} frame={s.frame} "
                        f"d={d} depth={s.depth} "
                        f"lhs={lhs} rhs={rhs}"
                    )

                failures += 1

    checks = sum(s.depth for s in states)

    print(f"checks={checks}")
    print(f"failures={failures}")
    print()

    return failures


# ==========================================================================================
# TEST 2
# ==========================================================================================

def test_first_break(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 2: FIRST BRANCH-BREAK LEVEL")
    print("=" * 90)

    failures = 0
    counts = Counter()

    for s in states:
        if s.tp >= INF:
            branch_level = INF
        else:
            branch_level = s.tp + 1

        if branch_level < s.w:
            kind = "branch-first"
        elif branch_level > s.w:
            kind = "n-first"
        else:
            kind = "simultaneous"

        counts[(s.frame, kind)] += 1

        # Direct verification.
        expected = min(branch_level, s.w)

        if expected != s.depth:
            failures += 1

    for key in sorted(counts):
        print(f"    {key}: {counts[key]}")

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# ==========================================================================================
# TEST 3
# ==========================================================================================

def test_modular_factor_signatures(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 3: MODULAR FACTOR SIGNATURES")
    print("=" * 90)

    failures = 0

    # A depth d is determined by a factor prefix of d-1 bits.
    #
    # Test whether:
    #
    #     (frame, n mod 2^d, p mod 2^(d-1))
    #
    # determines depth exactly.

    for d in range(2, 13):
        groups: dict[tuple, set[int]] = defaultdict(set)

        for s in states:
            if s.depth < d:
                continue

            sig = (
                s.frame,
                s.n % (1 << d),
                s.p % (1 << (d - 1)),
            )

            groups[sig].add(s.depth)

        ambiguous = sum(
            1
            for depths in groups.values()
            if len(depths) > 1
        )

        if ambiguous:
            print(
                f"    d={d:2d} "
                f"signatures={len(groups):6d} "
                f"ambiguous={ambiguous:6d}"
            )
            failures += ambiguous
        else:
            print(
                f"    d={d:2d} "
                f"signatures={len(groups):6d} "
                f"ambiguous=0"
            )

    print()
    return failures


# ==========================================================================================
# TEST 4
# ==========================================================================================

def test_sum_difference_signatures(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 4: SUM/DIFFERENCE PREFIX SIGNATURE")
    print("=" * 90)

    failures = 0

    for d in range(2, 13):
        mod_n = 1 << d
        mod_f = 1 << (d - 1)

        groups_sum: dict[tuple, set[int]] = defaultdict(set)
        groups_diff: dict[tuple, set[int]] = defaultdict(set)

        for s in states:
            if s.depth < d:
                continue

            sig_sum = (
                s.frame,
                s.n % mod_n,
                (s.p + s.q) % mod_f,
            )

            sig_diff = (
                s.frame,
                s.n % mod_n,
                (s.p - s.q) % mod_f,
            )

            groups_sum[sig_sum].add(s.depth)
            groups_diff[sig_diff].add(s.depth)

        amb_sum = sum(
            len(v) > 1 for v in groups_sum.values()
        )

        amb_diff = sum(
            len(v) > 1 for v in groups_diff.values()
        )

        print(
            f"    d={d:2d} "
            f"sum_ambiguous={amb_sum:6d} "
            f"diff_ambiguous={amb_diff:6d}"
        )

        failures += amb_sum + amb_diff

    print()
    return failures


# ==========================================================================================
# TEST 5
# ==========================================================================================

def test_threshold_without_v2(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 5: DEPTH AS FIRST FAILED CONGRUENCE")
    print("=" * 90)

    failures = 0

    for s in states:
        reconstructed = None

        # Search only until the actual w-bound.
        for d in range(2, s.w + 2):
            if not threshold(s, d):
                reconstructed = d - 1
                break

        if reconstructed is None:
            reconstructed = s.w

        if reconstructed != s.depth:
            if failures < 20:
                print(
                    f"mismatch n={s.n} p={s.p} q={s.q} "
                    f"actual={s.depth} reconstructed={reconstructed}"
                )

            failures += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()

    return failures


# ==========================================================================================
# TEST 6
# ==========================================================================================

def test_branch_bit(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 6: FIRST-DIFFERING-BIT ENCODING")
    print("=" * 90)

    failures = 0
    distribution = Counter()

    for s in states:
        if s.tp >= INF:
            kind = "never-breaks-before-n-bound"
        else:
            # First differing bit is bit tp.
            #
            # p-p0 has:
            #   bits 0..tp-1 = 0
            #   bit tp       = 1

            mask_low = (1 << s.tp) - 1

            low_zero = (
                (s.p - s.p0) % (mask_low + 1) == 0
            )

            next_bit_one = (
                ((s.p - s.p0) & (1 << s.tp)) != 0
            )

            if not low_zero or not next_bit_one:
                if failures < 20:
                    print(
                        f"branch-bit failure n={s.n} "
                        f"p={s.p} p0={s.p0} tp={s.tp}"
                    )

                failures += 1

            kind = f"break@{s.tp}"

        distribution[(s.frame, kind)] += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")

    print()
    print("    branch distribution:")

    for key, count in sorted(distribution.items()):
        print(f"        {key}: {count}")

    print()

    return failures


# ==========================================================================================
# TEST 7
# ==========================================================================================

def test_compressed_branch_encoding(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 7: COMPRESSED BRANCH ENCODING")
    print("=" * 90)

    failures = 0

    groups: dict[tuple, set[int]] = defaultdict(set)

    for s in states:
        if s.tp < s.w:
            # Branch breaks before the n-bound.
            #
            # The only needed branch information is tp.
            #
            # Everything else should be irrelevant.
            branch_code = ("break", s.tp)
        else:
            # The n-bound wins.
            branch_code = ("bound",)

        groups[(s.frame, s.w, branch_code)].add(s.depth)

    ambiguous = 0

    for sig, depths in groups.items():
        if len(depths) > 1:
            ambiguous += 1

            if ambiguous <= 20:
                print(
                    f"    ambiguous={sig} depths={sorted(depths)}"
                )

    print(f"signatures={len(groups)}")
    print(f"ambiguous={ambiguous}")

    failures += ambiguous

    print()
    return failures


# ==========================================================================================
# EXAMPLES
# ==========================================================================================

def print_examples(states: list[State]) -> None:
    print("=" * 90)
    print("EXAMPLES")
    print("=" * 90)

    wanted = {
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
        62922749,
    }

    seen = set()

    for s in states:
        if s.n not in wanted:
            continue

        seen.add(s.n)

        print()
        print(
            f"n={s.n} p={s.p} q={s.q} frame={s.frame}"
        )
        print(
            f"    p0={s.p0} c={s.c} "
            f"tp={s.tp} w={s.w}"
        )
        print(
            f"    depth={s.depth} "
            f"min(tp+1,w)={min(s.tp + 1, s.w)}"
        )

        if s.tp < INF:
            print(
                f"    branch break level={s.tp + 1}"
            )
        else:
            print(
                "    branch never breaks before n-bound"
            )

        print("    threshold prefix:")

        for d in range(2, min(s.depth + 1, 6)):
            print(
                f"        d={d}: "
                f"p-prefix="
                f"{(s.p - s.p0) % (1 << (d - 1)) == 0} "
                f"n-prefix="
                f"{(s.n + s.c) % (1 << d) == 0}"
            )

    print()
    print(
        f"examples shown={len(seen)}"
    )
    print()


# ==========================================================================================
# MAIN
# ==========================================================================================

def main() -> None:
    prime_limit = 6000

    if len(sys.argv) > 1:
        prime_limit = int(sys.argv[1])

    print("=" * 90)
    print("EXPERIMENT 652 START")
    print("=" * 90)
    print()
    print(f"prime limit={prime_limit}")

    primes = sieve(prime_limit)

    print()
    print(f"odd primes={len(primes)}")

    states = make_states(primes)

    print(f"semiprimes={len(states)}")
    print()

    failures = 0

    failures += test_baseline(states)
    failures += test_threshold_theorem(states)
    failures += test_first_break(states)
    failures += test_modular_factor_signatures(states)
    failures += test_sum_difference_signatures(states)
    failures += test_threshold_without_v2(states)
    failures += test_branch_bit(states)
    failures += test_compressed_branch_encoding(states)

    print_examples(states)

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)
    print(
        """
Experiment 651 established:

    depth = min(tp+1,w)

with:

    tp=v2(p-p0)
    w =v2(n+c).

Experiment 652 reformulates that as a prefix theorem:

    depth >= d

        iff

    p == p0 (mod 2^(d-1))
    and
    n == -c (mod 2^d).

Thus the branch valuation itself is equivalent to a
sequence of binary survival decisions:

    level 1:
        p matches anchor

    level 2:
        p matches anchor

    ...

    level tp:
        p still matches anchor

    level tp+1:
        first mismatch.

The hierarchy therefore has two independent stopping
mechanisms:

    BRANCH STOP:
        first differing bit of p from p0

    N STOP:
        first missing bit in n from -c.

The exact depth is the earlier stopping event.

The important question after this experiment is whether
the first branch mismatch itself has a simpler algebraic
description for these prime branches.

Candidates include:

    p mod 2^d
    p+q mod 2^d
    p-q mod 2^d
    gcd(p-p0, 2^d)
    the first set bit of p-p0.

If those collapse further, the remaining factor-side
information may reduce from an integer valuation to a
single discrete branch event.
"""
    )

    print()
    print("=" * 90)
    print("EXPERIMENT 652 FINISHED")
    print("=" * 90)
    print(f"TOTAL FAILURES={failures}")

    if failures == 0:
        print("STATUS=ALL TESTS PASSED")
    else:
        print("STATUS=COUNTEREXAMPLES FOUND")


if __name__ == "__main__":
    main()
