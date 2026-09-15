#!/usr/bin/env python3
"""
==========================================================================================
EXPERIMENT 651
==========================================================================================

2-ADIC SURVIVAL PREFIX / THRESHOLD DEPTH THEOREM

Goal
----

Experiment 650 showed that a fixed cap on

    tp = v2(p-p0)

is not universal.

The exact law was:

    depth = min(tp + 1, w)

where

    w = v2(n+c)

and:

    FRAME A:
        p0 =  3
        c  =  9

    FRAME B:
        p0 = -1
        c  =  3

This experiment removes the explicit valuation formula.

For d >= 2 we test the equivalent threshold statement:

    depth >= d

        iff

    p == p0 (mod 2^(d-1))
    AND
    n == -c (mod 2^d).

Therefore depth is the largest d for which the branch survives
both simultaneous 2-adic constraints.

The experiment also identifies the exact terminal event:

    BRANCH BREAK:
        p != p0 (mod 2^depth)

    or

    N-BOUND:
        n != -c (mod 2^(depth+1))

The first case means the prime branch terminates the hierarchy.
The second means the n-only valuation terminates it.

This gives a finite-prefix interpretation of the hierarchy
without introducing any artificial valuation cap.
"""

from __future__ import annotations

import math
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Iterable


INF = 10**9


# ==========================================================================================
# BASIC 2-ADIC UTILITIES
# ==========================================================================================

def v2(x: int) -> int:
    """Return v2(|x|), with v2(0)=INF."""
    if x == 0:
        return INF
    x = abs(x)
    return (x & -x).bit_length() - 1


def is_odd_prime(n: int) -> bool:
    if n < 3 or n % 2 == 0:
        return False
    if n % 3 == 0:
        return n == 3

    r = math.isqrt(n)
    f = 5
    step = 2

    while f <= r:
        if n % f == 0:
            return False
        f += step
        step = 6 - step

    return True


def odd_primes_upto(limit: int) -> list[int]:
    """Fast sieve for odd primes."""
    if limit < 3:
        return []

    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"

    for p in range(2, math.isqrt(limit) + 1):
        if sieve[p]:
            start = p * p
            sieve[start::p] = b"\x00" * (((limit - start) // p) + 1)

    return [p for p in range(3, limit + 1, 2) if sieve[p]]


# ==========================================================================================
# STATE
# ==========================================================================================

@dataclass(frozen=True)
class State:
    n: int
    p: int
    q: int
    frame: str

    # Raw residuals
    A: int
    B: int

    # Exact direct depth
    depth: int

    # Branch quantities
    p0: int
    c: int

    tp: int
    tq: int
    w: int


# ==========================================================================================
# FRAME / COORDINATE MAP
# ==========================================================================================

def frame_data(p: int, q: int) -> tuple[str, int, int, int, int]:
    """
    Return:

        frame, A, B, p0, c

    FRAME A:
        n = pq == 3 (mod 4)
        A = p - 3
        B = q + 3
        p0 = 3
        c  = 9

    FRAME B:
        n = pq == 1 (mod 4)
        A = p + 1
        B = q - 3
        p0 = -1
        c  = 3
    """
    n = p * q

    if n % 4 == 3:
        return "A", p - 3, q + 3, 3, 9

    return "B", p + 1, q - 3, -1, 3


def direct_depth(p: int, q: int, frame: str) -> int:
    """
    Exact depth from the transformed numerators.

    FRAME A:
        2X = q - p + 6
        2Y = p + q

    FRAME B:
        2X = 3p - q + 6
        2Y = 3p + q

    Therefore:

        depth = min(v2(2X), v2(2Y)).
    """
    if frame == "A":
        first = q - p + 6
        second = p + q
    else:
        first = 3 * p - q + 6
        second = 3 * p + q

    return min(v2(first), v2(second))


def build_states(primes: list[int]) -> list[State]:
    states: list[State] = []

    for i, p in enumerate(primes):
        for q in primes[i:]:
            n = p * q

            frame, A, B, p0, c = frame_data(p, q)

            d = direct_depth(p, q, frame)

            tp = v2(p - p0)
            tq = v2(q - (3 if frame == "B" else -3))
            w = v2(n + c)

            states.append(
                State(
                    n=n,
                    p=p,
                    q=q,
                    frame=frame,
                    A=A,
                    B=B,
                    depth=d,
                    p0=p0,
                    c=c,
                    tp=tp,
                    tq=tq,
                    w=w,
                )
            )

    return states


# ==========================================================================================
# EXPECTED STRUCTURAL LAW
# ==========================================================================================

def expected_depth_from_branch(state: State) -> int:
    return min(state.tp + 1, state.w)


# ==========================================================================================
# THRESHOLD CONDITION
# ==========================================================================================

def survives_level(state: State, level: int) -> bool:
    """
    Test the proposed survival-prefix condition.

    For level d >= 2:

        depth >= d

        iff

        p == p0 (mod 2^(d-1))
        and
        n == -c (mod 2^d)
    """
    if level < 2:
        return True

    p_mod = 1 << (level - 1)
    n_mod = 1 << level

    return (
        (state.p - state.p0) % p_mod == 0
        and
        (state.n + state.c) % n_mod == 0
    )


def reconstructed_depth(state: State) -> int:
    """
    Reconstruct depth by finding the first level where the
    simultaneous threshold condition fails.

    This does not use tp directly.
    """
    # Depth is never smaller than 2 in this domain.
    d = 2

    # The exact theorem gives a safe finite upper bound.
    upper = max(state.w, state.tp + 1 if state.tp < INF else state.w) + 2

    while d <= upper and survives_level(state, d):
        d += 1

    return d - 1


# ==========================================================================================
# TEST 0
# ==========================================================================================

def test_domain(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 0: ODD SEMIPRIME DOMAIN")
    print("=" * 90)

    failures = 0

    for s in states:
        if s.p % 2 == 0 or s.q % 2 == 0:
            failures += 1
        if not is_odd_prime(s.p) or not is_odd_prime(s.q):
            failures += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()
    return failures


# ==========================================================================================
# TEST 1
# ==========================================================================================

def test_known_depth_law(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 1: EXPERIMENT 648 BASELINE")
    print("=" * 90)

    failures = 0

    for s in states:
        predicted = expected_depth_from_branch(s)

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
# TEST 2
# ==========================================================================================

def test_threshold_equivalence(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 2: EXACT SURVIVAL THRESHOLD THEOREM")
    print("=" * 90)

    failures = 0

    for s in states:
        max_level = s.depth + 1

        # Exact statement:
        #
        # depth >= d  <=>  threshold condition holds at d
        #
        # Test all relevant d, including the first failing level.
        for d in range(2, max_level + 1):
            lhs = s.depth >= d
            rhs = survives_level(s, d)

            if lhs != rhs:
                if failures < 20:
                    print(
                        f"mismatch n={s.n} p={s.p} q={s.q} "
                        f"frame={s.frame} d={d} "
                        f"depth={s.depth} lhs={lhs} rhs={rhs}"
                    )
                failures += 1

    total_checks = sum(s.depth for s in states) - len(states)

    print(f"states={len(states)}")
    print(f"threshold_checks={total_checks}")
    print(f"failures={failures}")
    print()
    return failures


# ==========================================================================================
# TEST 3
# ==========================================================================================

def test_reconstruction(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 3: DEPTH RECONSTRUCTION WITHOUT tp")
    print("=" * 90)

    failures = 0

    for s in states:
        predicted = reconstructed_depth(s)

        if predicted != s.depth:
            if failures < 20:
                print(
                    f"mismatch n={s.n} p={s.p} q={s.q} "
                    f"frame={s.frame} actual={s.depth} "
                    f"reconstructed={predicted}"
                )
            failures += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()
    return failures


# ==========================================================================================
# TEST 4
# ==========================================================================================

def test_terminal_event(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 4: TERMINAL EVENT CLASSIFICATION")
    print("=" * 90)

    failures = 0
    counts = Counter()

    examples: dict[str, State] = {}

    for s in states:
        d = s.depth

        # At level d, survival must hold.
        if not survives_level(s, d):
            if failures < 20:
                print(
                    f"depth-level failure n={s.n} "
                    f"frame={s.frame} depth={d}"
                )
            failures += 1
            continue

        # At d+1, at least one of the two constraints must fail.
        branch_ok = (
            (s.p - s.p0) % (1 << d) == 0
        )
        n_ok = (
            (s.n + s.c) % (1 << (d + 1)) == 0
        )

        if branch_ok and n_ok:
            if failures < 20:
                print(
                    f"no terminal event n={s.n} p={s.p} q={s.q} "
                    f"frame={s.frame} depth={d}"
                )
            failures += 1
            continue

        if not branch_ok and not n_ok:
            kind = "both-break"
        elif not branch_ok:
            kind = "branch-break"
        else:
            kind = "n-bound"

        counts[(s.frame, kind)] += 1
        examples.setdefault((s.frame, kind), s)

    for key in sorted(counts):
        print(f"    {key}: {counts[key]}")

    print()
    print("    examples:")

    for key in sorted(examples):
        s = examples[key]
        print(
            f"        {key}: "
            f"n={s.n} p={s.p} q={s.q} "
            f"depth={s.depth} tp={s.tp} w={s.w}"
        )

    print()
    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()
    return failures


# ==========================================================================================
# TEST 5
# ==========================================================================================

def test_piecewise_equivalence(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 5: PIECEWISE SURVIVAL LAW")
    print("=" * 90)

    failures = 0

    for s in states:
        # Threshold formulation
        a = s.tp + 1
        b = s.w
        expected = min(a, b)

        if expected != s.depth:
            if failures < 20:
                print(
                    f"mismatch n={s.n} frame={s.frame} "
                    f"tp+1={a} w={b} depth={s.depth}"
                )
            failures += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()
    return failures


# ==========================================================================================
# TEST 6
# ==========================================================================================

def test_unbounded_branch_structure(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 6: UNBOUNDED BRANCH / FIRST-DIFFERING-BIT")
    print("=" * 90)

    finite_tp = [s.tp for s in states if s.tp < INF]

    max_tp = max(finite_tp) if finite_tp else 0
    max_w = max(s.w for s in states)

    distribution = Counter()

    for s in states:
        if s.tp >= s.w:
            distribution["n-bound"] += 1
        elif s.tp < s.w:
            distribution["branch-break"] += 1
        else:
            distribution["degenerate"] += 1

    print(f"max finite tp={max_tp}")
    print(f"max w={max_w}")

    for key in sorted(distribution):
        print(f"    {key}: {distribution[key]}")

    # The mathematical branch claim:
    #
    # tp = r
    # =>
    # p == p0 mod 2^r
    # but
    # p != p0 mod 2^(r+1)
    #
    # Check it directly.
    failures = 0

    for s in states:
        if s.tp >= INF:
            continue

        if (s.p - s.p0) % (1 << s.tp) != 0:
            failures += 1

        if (s.p - s.p0) % (1 << (s.tp + 1)) == 0:
            failures += 1

    print(f"checked={len(states)}")
    print(f"failures={failures}")
    print()
    return failures


# ==========================================================================================
# TEST 7
# ==========================================================================================

def test_information_minimality(states: list[State]) -> int:
    print("=" * 90)
    print("TEST 7: DOES (frame, w, RELATIVE BRANCH POSITION) DETERMINE DEPTH?")
    print("=" * 90)

    groups: dict[tuple, set[int]] = defaultdict(set)

    for s in states:
        # This is intentionally much smaller than storing full tp.
        #
        # We only record:
        #
        #   tp < w
        #   tp >= w
        #
        # and, in the branch-breaking case, the exact first
        # differing level tp.
        #
        # The signature below should reconstruct depth exactly.
        if s.tp < s.w:
            branch = ("break", s.tp)
        else:
            branch = ("bound",)

        groups[(s.frame, s.w, branch)].add(s.depth)

    ambiguous = 0

    for signature, depths in groups.items():
        if len(depths) > 1:
            ambiguous += 1
            if ambiguous <= 20:
                print(
                    f"    ambiguous signature={signature} depths={sorted(depths)}"
                )

    failures = ambiguous

    print(f"signatures={len(groups)}")
    print(f"ambiguous={ambiguous}")
    print()

    return failures


# ==========================================================================================
# EXAMPLES
# ==========================================================================================

def print_examples(states: list[State]) -> None:
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
        1047,
        4135,
        40999,
        485879,
        5579767,
        62922749,
    }

    print("=" * 90)
    print("EXAMPLES")
    print("=" * 90)

    found = 0

    for s in states:
        if s.n not in wanted:
            continue

        found += 1

        print()
        print(f"n={s.n} p={s.p} q={s.q} frame={s.frame}")
        print(f"    p0={s.p0} c={s.c}")
        print(f"    tp=v2(p-p0)={s.tp}")
        print(f"    tq=v2(q-q0)={s.tq}")
        print(f"    w=v2(n+c)={s.w}")
        print(f"    actual depth={s.depth}")
        print(f"    predicted depth=min(tp+1,w)={min(s.tp + 1, s.w)}")

        for d in range(2, min(s.depth + 1, 5)):
            print(
                f"    d={d} "
                f"p-prefix={survives_level(s, d) and ((s.p - s.p0) % (1 << (d - 1)) == 0)} "
                f"n-bound={((s.n + s.c) % (1 << d) == 0)} "
                f"survives={survives_level(s, d)}"
            )

        branch_break = (
            (s.p - s.p0) % (1 << s.depth) != 0
        )

        n_bound_break = (
            (s.n + s.c) % (1 << (s.depth + 1)) != 0
        )

        print(
            f"    terminal branch-break={branch_break} "
            f"n-bound={n_bound_break}"
        )

    if found == 0:
        print("no selected examples found")

    print()


# ==========================================================================================
# MAIN
# ==========================================================================================

def run(prime_limit: int) -> int:
    print("=" * 90)
    print("EXPERIMENT 651 START")
    print("=" * 90)
    print()
    print(f"prime limit={prime_limit}")

    primes = odd_primes_upto(prime_limit)

    print()
    print(f"odd primes={len(primes)}")

    states = build_states(primes)

    print(f"semiprimes={len(states)}")
    print()

    failures = 0

    failures += test_domain(states)
    failures += test_known_depth_law(states)
    failures += test_threshold_equivalence(states)
    failures += test_reconstruction(states)
    failures += test_terminal_event(states)
    failures += test_piecewise_equivalence(states)
    failures += test_unbounded_branch_structure(states)
    failures += test_information_minimality(states)

    print_examples(states)

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)
    print(
        """
The tested exact law is:

    depth = min(tp + 1, w)

where:

    tp = v2(p-p0)
    w  = v2(n+c)

with:

    FRAME A:
        p0 = 3
        c  = 9

    FRAME B:
        p0 = -1
        c  = 3.

The equivalent threshold theorem is:

    depth >= d

    iff

    p == p0 (mod 2^(d-1))
    and
    n == -c (mod 2^d).

Thus the hierarchy is a simultaneous survival-prefix:

    prime branch:
        survives through d
        iff p matches the frame anchor through d-1 bits

    n-bound:
        survives through d
        iff n matches -c through d bits.

The depth is the first level at which one of these two
constraints terminates the branch.

This removes the need for an artificial finite cap.

The remaining factor-side information is exactly the
first 2-adic divergence level of p from its frame anchor.
"""
    )

    print()
    print("=" * 90)
    print("EXPERIMENT 651 FINISHED")
    print("=" * 90)
    print(f"TOTAL FAILURES={failures}")

    if failures == 0:
        print("STATUS=ALL TESTS PASSED")
    else:
        print("STATUS=COUNTEREXAMPLES FOUND")

    return failures


def main() -> None:
    if len(sys.argv) > 1:
        try:
            prime_limit = int(sys.argv[1])
        except ValueError:
            raise SystemExit("Usage: python main.py [prime_limit]")
    else:
        prime_limit = 6000

    if prime_limit < 3:
        raise SystemExit("prime_limit must be >= 3")

    raise SystemExit(run(prime_limit))


if __name__ == "__main__":
    main()
