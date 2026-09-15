#!/usr/bin/env python3

from __future__ import annotations

from collections import defaultdict, Counter
from math import isqrt


INF = 10**9


# =============================================================================
# EXPERIMENT 649
# BRANCH-VALUATION TRUNCATION / FIRST-DIFFERING-BIT TEST
#
# Experiment 648 established exactly:
#
#     m = min(tp, w)
#
#     depth = min(tp + 1, w)
#
# where
#
# FRAME A:
#     tp = v2(p-3)
#     w  = v2(n+9)
#
# FRAME B:
#     tp = v2(p+1)
#     w  = v2(n+3)
#
# This experiment asks whether the full integer tp is necessary.
#
# Candidate compressed observables:
#
#     min(tp, w)
#     min(tp, w-1)
#     tp capped at small thresholds
#     tp parity
#     tp relative to w
#     first differing bit of p from p0
#
# The central test is:
#
#     (frame, w, tp_cap)
#         -> depth
#
# for progressively small caps.
#
# If a small cap is exact, the infinite-looking valuation
# collapses to a finite-state branch rule.
# =============================================================================


def v2(x: int) -> int:
    if x == 0:
        return INF
    return (abs(x) & -abs(x)).bit_length() - 1


# =============================================================================
# PRIME SIEVE
# =============================================================================

def sieve(limit: int) -> list[int]:
    flags = bytearray(b"\x01") * (limit + 1)
    flags[0] = 0
    flags[1] = 0

    for p in range(2, isqrt(limit) + 1):
        if flags[p]:
            start = p * p
            flags[start::p] = b"\x00" * (
                (limit - start) // p + 1
            )

    return [
        p
        for p in range(3, limit + 1, 2)
        if flags[p]
    ]


# =============================================================================
# FRAME
# =============================================================================

def frame_from_n(n: int) -> str:
    r = n & 3

    if r == 3:
        return "A"

    if r == 1:
        return "B"

    raise ValueError(
        f"unexpected odd n={n}, n mod 4={r}"
    )


def frame_params(frame: str):
    if frame == "A":
        # A = p - 3
        # B = q + 3
        # n+c = n+9
        return 3, -3, 9

    # FRAME B
    # A = p + 1
    # B = q - 3
    # n+c = n+3
    return -1, 3, 3


# =============================================================================
# STATES
# =============================================================================

def build_states(primes: list[int]):
    states = []

    for i, p in enumerate(primes):
        for q in primes[i:]:

            n = p * q
            frame = frame_from_n(n)

            p0, q0, c = frame_params(frame)

            A = p - p0
            B = q - q0

            tp = v2(A)
            tq = v2(B)
            w = v2(n + c)

            m = min(tp, tq)
            depth = min(tp + 1, w)

            # Baseline consistency with Experiment 648.
            assert m == min(tp, w)
            assert depth == min(tp + 1, w)

            states.append(
                {
                    "n": n,
                    "p": p,
                    "q": q,
                    "frame": frame,
                    "p0": p0,
                    "q0": q0,
                    "A": A,
                    "B": B,
                    "tp": tp,
                    "tq": tq,
                    "w": w,
                    "m": m,
                    "depth": depth,
                }
            )

    return states


# =============================================================================
# TEST 0
# BASELINE
# =============================================================================

def test_baseline(states) -> int:

    print("=" * 90)
    print("TEST 0: EXPERIMENT 648 BASELINE")
    print("=" * 90)

    failures = 0

    for s in states:

        if s["m"] != min(s["tp"], s["w"]):
            failures += 1

        if s["depth"] != min(
            s["tp"] + 1,
            s["w"],
        ):
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
# CAP FUNCTIONS
# =============================================================================

def cap_value(x: int, cap: int) -> int:
    return min(x, cap)


def branch_state_cap(
    tp: int,
    w: int,
    cap: int,
):
    """
    Finite-state candidate.

    The obvious safe cap is w, but we deliberately
    test fixed caps independent of w.
    """
    return min(tp, cap)


# =============================================================================
# TEST 1
# FIXED tp CAP
#
#     (frame, w, min(tp,C))
#         -> depth
#
# =============================================================================

def test_fixed_caps(states) -> int:

    print("=" * 90)
    print("TEST 1: FIXED BRANCH-VALUATION CAPS")
    print("=" * 90)

    total_failures = 0

    for cap in range(0, 14):

        buckets = defaultdict(set)

        for s in states:

            sig = (
                s["frame"],
                s["w"],
                branch_state_cap(
                    s["tp"],
                    s["w"],
                    cap,
                ),
            )

            buckets[sig].add(
                s["depth"]
            )

        ambiguous = {
            k: v
            for k, v in buckets.items()
            if len(v) > 1
        }

        print(
            f"    cap={cap:2d}"
            f" signatures={len(buckets):6d}"
            f" ambiguous={len(ambiguous):6d}"
        )

        total_failures += len(ambiguous)

    print()

    return total_failures


# =============================================================================
# TEST 2
# DOES ONLY tp < w MATTER?
#
# The exact formula is:
#
#     depth = tp+1   when tp < w
#     depth = w      when tp >= w
#
# Therefore candidate state:
#
#     relation = compare(tp,w)
#
# plus min(tp,w).
# =============================================================================

def test_relative_branch_state(states) -> int:

    print("=" * 90)
    print(
        "TEST 2: RELATIVE tp-vs-w STATE"
    )
    print("=" * 90)

    buckets = defaultdict(set)

    for s in states:

        tp = s["tp"]
        w = s["w"]

        if tp < w:
            relation = "tp<w"
        elif tp == w:
            relation = "tp=w"
        else:
            relation = "tp>w"

        sig = (
            s["frame"],
            relation,
            min(tp, w),
        )

        buckets[sig].add(
            s["depth"]
        )

    ambiguous = {
        k: v
        for k, v in buckets.items()
        if len(v) > 1
    }

    print(
        f"signatures={len(buckets)}"
    )
    print(
        f"ambiguous={len(ambiguous)}"
    )

    if ambiguous:

        print(
            "    first collisions:"
        )

        for key, vals in list(
            ambiguous.items()
        )[:20]:

            print(
                f"        {key}"
                f" -> {sorted(vals)}"
            )

    print()

    return len(ambiguous)


# =============================================================================
# TEST 3
# FIRST-DIFFERING-BIT CLASSIFICATION
#
# tp = v2(p-p0)
#
# means:
#
#     p == p0 mod 2^tp
#
# and:
#
#     p != p0 mod 2^(tp+1)
#
# This test asks whether only the first differing bit
# relative to the anchor is enough once w is known.
#
# We encode:
#
#     bit_position = min(tp, w-1)
#
# plus whether the branch survives all w bits.
# =============================================================================

def test_first_differing_bit(states) -> int:

    print("=" * 90)
    print(
        "TEST 3: FIRST-DIFFERING-BIT SIGNATURE"
    )
    print("=" * 90)

    buckets = defaultdict(set)

    for s in states:

        tp = s["tp"]
        w = s["w"]

        if tp >= w:
            bit = w
            branch = "survives_to_w"
        else:
            bit = tp
            branch = "dies_at_bit"

        sig = (
            s["frame"],
            w,
            bit,
            branch,
        )

        buckets[sig].add(
            s["depth"]
        )

    ambiguous = {
        k: v
        for k, v in buckets.items()
        if len(v) > 1
    }

    print(
        f"signatures={len(buckets)}"
    )
    print(
        f"ambiguous={len(ambiguous)}"
    )

    if ambiguous:

        print(
            "    first collisions:"
        )

        for key, vals in list(
            ambiguous.items()
        )[:20]:

            print(
                f"        {key}"
                f" -> {sorted(vals)}"
            )

    print()

    return len(ambiguous)


# =============================================================================
# TEST 4
# CAN DEPTH BE WRITTEN WITHOUT tp?
#
# Since:
#
#     depth = min(tp+1,w)
#
# a tempting finite-state form is:
#
#     depth = w                  if tp >= w
#             branch_position+1  otherwise
#
# Test the branch-position representation directly.
# =============================================================================

def test_branch_position(states) -> int:

    print("=" * 90)
    print(
        "TEST 4: BRANCH POSITION -> DEPTH"
    )
    print("=" * 90)

    failures = 0

    for s in states:

        tp = s["tp"]
        w = s["w"]

        if tp >= w:
            position = w
            predicted = w
        else:
            position = tp
            predicted = position + 1

        if predicted != s["depth"]:

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch"
                    f" n={s['n']}"
                    f" frame={s['frame']}"
                    f" tp={tp}"
                    f" w={w}"
                    f" position={position}"
                    f" actual={s['depth']}"
                    f" predicted={predicted}"
                )

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
# TEST WHETHER tp IS NEEDED EXACTLY OR ONLY tp+1
#
# Candidate:
#
#     branch_level = min(tp+1, w)
#
# This is numerically equal to depth by definition.
#
# We instead test whether the raw branch valuation can
# be replaced by the valuation of a SINGLE affine bit:
#
#     v2((p-p0) / 2^r - 1)
#
# at small r.
#
# This is mainly diagnostic.
# =============================================================================

def test_affine_branch_layers(states) -> int:

    print("=" * 90)
    print(
        "TEST 5: BRANCH LAYER STATISTICS"
    )
    print("=" * 90)

    counts = Counter()

    for s in states:

        tp = s["tp"]
        w = s["w"]

        # Which level terminates the branch?
        if tp < w:
            level = tp
            kind = "branch-break"
        else:
            level = w
            kind = "n-bound"

        counts[
            (
                s["frame"],
                level,
                kind,
            )
        ] += 1

    print(
        "    first branch layers:"
    )

    for key in sorted(
        counts
    ):

        frame, level, kind = key

        if level <= 10:

            print(
                f"        frame={frame}"
                f" level={level}"
                f" {kind}"
                f" count={counts[key]}"
            )

    print()

    return 0


# =============================================================================
# TEST 6
# KNOWN EXAMPLES
# =============================================================================

def print_examples(states):

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
        1047,
        4135,
        40999,
    }

    for s in states:

        if s["n"] not in wanted:
            continue

        print()

        print(
            f"n={s['n']}"
            f" p={s['p']}"
            f" q={s['q']}"
            f" frame={s['frame']}"
        )

        print(
            f"    p0={s['p0']}"
            f" q0={s['q0']}"
        )

        print(
            f"    tp=v2(p-p0)={s['tp']}"
        )

        print(
            f"    tq=v2(q-q0)={s['tq']}"
        )

        print(
            f"    w=v2(n+c)={s['w']}"
        )

        print(
            f"    m={s['m']}"
        )

        print(
            f"    depth={s['depth']}"
        )

        print(
            f"    tp<w = {s['tp'] < s['w']}"
        )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 90)
    print("EXPERIMENT 649 START")
    print("=" * 90)
    print()

    primes = sieve(6000)

    print(
        "[1] PRIME SIEVE"
    )

    print(
        f"    odd primes={len(primes)}"
    )

    print()

    states = build_states(
        primes
    )

    print(
        "[2] SEMIPRIME GENERATION"
    )

    print(
        f"    states={len(states)}"
    )

    print()

    failures = 0

    failures += test_baseline(
        states
    )

    failures += test_fixed_caps(
        states
    )

    failures += test_relative_branch_state(
        states
    )

    failures += test_first_differing_bit(
        states
    )

    failures += test_branch_position(
        states
    )

    test_affine_branch_layers(
        states
    )

    print_examples(
        states
    )

    print()
    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)

    print(
r"""
Experiment 648 established the exact compression:

    m = min(
        v2(p-p0),
        v2(n+c)
    )

and therefore:

    depth =
        min(
            v2(p-p0)+1,
            v2(n+c)
        ).

Experiment 649 asks whether even the integer valuation

    v2(p-p0)

is necessary as a stored quantity.

The valuation has a direct bit interpretation:

    tp = v2(p-p0)

means:

    p == p0 (mod 2^tp)

but:

    p != p0 (mod 2^(tp+1))

when tp is finite.

Thus tp identifies the FIRST 2-adic bit at which the
prime branch leaves the frame anchor.

The exact depth law can therefore be viewed as:

    depth = first terminal level

where termination occurs either because:

    1. the prime branch differs from the anchor, or
    2. the n-only bound w is reached.

The key structural representation is:

    branch survives through levels < tp
    branch breaks at level tp.

Consequently:

    if tp < w:

        depth = tp + 1

    if tp >= w:

        depth = w.

The experiment tests whether this first-differing-bit
description can replace the full valuation and whether
small fixed truncations of tp already determine depth.

A successful small-cap result would mean:

    finite branch state
        +
    v2(n+c)
        ->
    exact depth.

A failure means the branch position must remain
unbounded, although it is still only ONE factor-side
valuation rather than the original two-factor state.
"""
    )

    print()

    print("=" * 90)
    print("EXPERIMENT 649 FINISHED")
    print("=" * 90)

    print(
        f"TOTAL FAILURES={failures}"
    )

    if failures == 0:
        print(
            "STATUS=ALL TESTS PASSED"
        )
    else:
        print(
            "STATUS=COUNTEREXAMPLES FOUND"
        )


if __name__ == "__main__":
    main()
