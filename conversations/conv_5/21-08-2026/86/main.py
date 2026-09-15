#!/usr/bin/env python3

from __future__ import annotations

from collections import defaultdict, Counter
from math import isqrt


# =============================================================================
# EXPERIMENT 648
# SINGLE-VALUATION BRANCH COLLAPSE
#
# Target:
#
#     m = min(
#         v2(branch_displacement),
#         v2(n+c)
#     )
#
# and therefore:
#
#     depth =
#         min(
#             v2(branch_displacement) + 1,
#             v2(n+c)
#         )
#
# FRAME A:
#     branch displacement = p - 3
#     n+c = n+9
#
# FRAME B:
#     branch displacement = p + 1
#     n+c = n+3
#
# This tests whether the full p mod 2^k signature from
# Experiment 647 collapses to ONE integer valuation.
# =============================================================================


INF = 10**9


# =============================================================================
# 2-ADIC VALUATION
# =============================================================================

def v2(x: int) -> int:
    if x == 0:
        return INF

    x = abs(x)
    return (x & -x).bit_length() - 1


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
        f"unexpected odd semiprime n={n}, n mod 4={r}"
    )


def frame_constants(frame: str):
    if frame == "A":
        # A = p-3
        # B = q+3
        # n+c = n+9
        return 3, -3, 9

    # FRAME B:
    # A = p+1
    # B = q-3
    # n+c = n+3
    return -1, 3, 3


# =============================================================================
# SIEVE
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
# BUILD SEMIPRIME STATES
# =============================================================================

def build_states(primes: list[int]):
    states = []

    for i, p in enumerate(primes):

        for q in primes[i:]:

            n = p * q
            frame = frame_from_n(n)

            if frame == "A":

                A = p - 3
                B = q + 3

            else:

                A = p + 1
                B = q - 3

            m = min(
                v2(A),
                v2(B),
            )

            equality = (
                v2(A) == v2(B)
            )

            depth = (
                m + int(equality)
            )

            states.append(
                (
                    n,
                    p,
                    q,
                    frame,
                    A,
                    B,
                    m,
                    depth,
                )
            )

    return states


# =============================================================================
# TEST 0
# BASELINE DOMAIN
# =============================================================================

def test_domain(states) -> int:

    print("=" * 90)
    print("TEST 0: ODD SEMIPRIME DOMAIN")
    print("=" * 90)

    failures = 0

    for (
        n,
        p,
        q,
        frame,
        A,
        B,
        m,
        depth,
    ) in states:

        if n % 2 != 1:
            failures += 1

        if p * q != n:
            failures += 1

        if frame_from_n(n) != frame:
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
# RAW m == min(branch valuation, n+c valuation)
# =============================================================================

def test_m_from_p_branch(states) -> int:

    print("=" * 90)
    print(
        "TEST 1: m = min(v2(P-BRANCH), v2(n+c))"
    )
    print("=" * 90)

    failures = 0

    for (
        n,
        p,
        q,
        frame,
        A,
        B,
        m,
        depth,
    ) in states:

        p0, q0, c = frame_constants(frame)

        tp = v2(
            p - p0
        )

        w = v2(
            n + c
        )

        predicted = min(
            tp,
            w,
        )

        if predicted != m:

            failures += 1

            if failures <= 25:

                print(
                    "    mismatch"
                    f" n={n}"
                    f" p={p}"
                    f" q={q}"
                    f" frame={frame}"
                    f" tp={tp}"
                    f" w={w}"
                    f" actual_m={m}"
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
# TEST 2
# SAME RESULT FROM q BRANCH
# =============================================================================

def test_m_from_q_branch(states) -> int:

    print("=" * 90)
    print(
        "TEST 2: m = min(v2(Q-BRANCH), v2(n+c))"
    )
    print("=" * 90)

    failures = 0

    for (
        n,
        p,
        q,
        frame,
        A,
        B,
        m,
        depth,
    ) in states:

        p0, q0, c = frame_constants(frame)

        tq = v2(
            q - q0
        )

        w = v2(
            n + c
        )

        predicted = min(
            tq,
            w,
        )

        if predicted != m:

            failures += 1

            if failures <= 25:

                print(
                    "    mismatch"
                    f" n={n}"
                    f" p={p}"
                    f" q={q}"
                    f" frame={frame}"
                    f" tq={tq}"
                    f" w={w}"
                    f" actual_m={m}"
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
# TEST 3
# DIRECT DEPTH COLLAPSE
#
#     depth = min(branch valuation + 1, v2(n+c))
#
# =============================================================================

def test_depth_from_p_branch(states) -> int:

    print("=" * 90)
    print(
        "TEST 3: DEPTH = min(v2(P-BRANCH)+1, v2(n+c))"
    )
    print("=" * 90)

    failures = 0

    for (
        n,
        p,
        q,
        frame,
        A,
        B,
        m,
        depth,
    ) in states:

        p0, q0, c = frame_constants(frame)

        tp = v2(
            p - p0
        )

        w = v2(
            n + c
        )

        predicted = min(
            tp + 1,
            w,
        )

        if predicted != depth:

            failures += 1

            if failures <= 25:

                print(
                    "    mismatch"
                    f" n={n}"
                    f" p={p}"
                    f" q={q}"
                    f" frame={frame}"
                    f" tp={tp}"
                    f" w={w}"
                    f" actual_depth={depth}"
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
# TEST 4
# COLLAPSED SIGNATURE EXACTNESS
#
# Signature:
#
#     (frame, v2(n+c), v2(p-p0))
#
# must determine both m and depth.
# =============================================================================

def test_collapsed_signature(states) -> int:

    print("=" * 90)
    print(
        "TEST 4: COLLAPSED TWO-VALUATION SIGNATURE"
    )
    print("=" * 90)

    signatures = defaultdict(
        set
    )

    depth_signatures = defaultdict(
        set
    )

    for (
        n,
        p,
        q,
        frame,
        A,
        B,
        m,
        depth,
    ) in states:

        p0, q0, c = frame_constants(frame)

        tp = v2(
            p - p0
        )

        w = v2(
            n + c
        )

        signature = (
            frame,
            tp,
            w,
        )

        signatures[
            signature
        ].add(
            m
        )

        depth_signatures[
            signature
        ].add(
            depth
        )

    ambiguous_m = {
        k: v
        for k, v in signatures.items()
        if len(v) > 1
    }

    ambiguous_depth = {
        k: v
        for k, v in depth_signatures.items()
        if len(v) > 1
    }

    print(
        f"signatures={len(signatures)}"
    )
    print(
        f"ambiguous_m={len(ambiguous_m)}"
    )
    print(
        f"ambiguous_depth={len(ambiguous_depth)}"
    )

    if ambiguous_m:

        print(
            "    first m collisions:"
        )

        for key, values in list(
            ambiguous_m.items()
        )[:20]:

            print(
                f"        {key}"
                f" -> {sorted(values)}"
            )

    print()

    return (
        len(ambiguous_m)
        + len(ambiguous_depth)
    )


# =============================================================================
# TEST 5
# WHICH SIDE LIMITS?
# =============================================================================

def test_limiting_cases(states):

    print("=" * 90)
    print(
        "TEST 5: LIMITING CASE CLASSIFICATION"
    )
    print("=" * 90)

    counts = Counter()

    for (
        n,
        p,
        q,
        frame,
        A,
        B,
        m,
        depth,
    ) in states:

        p0, q0, c = frame_constants(frame)

        tp = v2(
            p - p0
        )

        tq = v2(
            q - q0
        )

        w = v2(
            n + c
        )

        if tp < tq:
            branch = "p-limits"

        elif tq < tp:
            branch = "q-limits"

        else:
            branch = "equal"

        # Relation between w and m:
        #
        # unequal residual valuations:
        #     w=m
        #
        # equal residual valuations:
        #     w>m

        if w == m:
            n_relation = "w=m"
        elif w > m:
            n_relation = "w>m"
        else:
            n_relation = "INVALID"

        counts[
            (
                frame,
                branch,
                n_relation,
            )
        ] += 1

    for key in sorted(
        counts
    ):

        print(
            f"    {key}:"
            f" {counts[key]}"
        )

    print()

    return 0


# =============================================================================
# TEST 6
# STRONGER FORM:
#
# depth =
#     if w <= tp:
#         w
#     else:
#         tp + 1
#
# equivalently:
#
#     depth = min(w, tp+1)
# =============================================================================

def test_piecewise_depth(states) -> int:

    print("=" * 90)
    print(
        "TEST 6: PIECEWISE DEPTH LAW"
    )
    print("=" * 90)

    failures = 0

    for (
        n,
        p,
        q,
        frame,
        A,
        B,
        m,
        depth,
    ) in states:

        p0, q0, c = frame_constants(frame)

        tp = v2(
            p - p0
        )

        w = v2(
            n + c
        )

        if w <= tp:
            predicted = w
        else:
            predicted = tp + 1

        if predicted != depth:

            failures += 1

            if failures <= 25:

                print(
                    "    mismatch"
                    f" n={n}"
                    f" frame={frame}"
                    f" tp={tp}"
                    f" w={w}"
                    f" actual={depth}"
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
# TEST 7
# KNOWN EXAMPLES
# =============================================================================

def print_examples(states):

    print("=" * 90)
    print(
        "EXAMPLES"
    )
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

    for state in states:

        (
            n,
            p,
            q,
            frame,
            A,
            B,
            m,
            depth,
        ) = state

        if n not in wanted:
            continue

        p0, q0, c = frame_constants(
            frame
        )

        tp = v2(
            p - p0
        )

        tq = v2(
            q - q0
        )

        w = v2(
            n + c
        )

        print()

        print(
            f"n={n}"
            f" p={p}"
            f" q={q}"
            f" frame={frame}"
        )

        print(
            f"    A={A}"
            f" B={B}"
        )

        print(
            f"    v2(p-p0)={tp}"
            f" p0={p0}"
        )

        print(
            f"    v2(q-q0)={tq}"
            f" q0={q0}"
        )

        print(
            f"    v2(n+c)={w}"
            f" c={c}"
        )

        print(
            f"    actual m={m}"
        )

        print(
            f"    min(tp,w)={min(tp,w)}"
        )

        print(
            f"    actual depth={depth}"
        )

        print(
            f"    min(tp+1,w)="
            f"{min(tp + 1, w)}"
        )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 90)
    print(
        "EXPERIMENT 648 START"
    )
    print("=" * 90)
    print()

    primes = sieve(
        6000
    )

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

    failures += test_domain(
        states
    )

    failures += test_m_from_p_branch(
        states
    )

    failures += test_m_from_q_branch(
        states
    )

    failures += test_depth_from_p_branch(
        states
    )

    failures += test_collapsed_signature(
        states
    )

    test_limiting_cases(
        states
    )

    failures += test_piecewise_depth(
        states
    )

    print_examples(
        states
    )

    print()
    print("=" * 90)
    print(
        "FINAL STRUCTURAL SUMMARY"
    )
    print("=" * 90)

    print(
r"""
Experiment 647 established:

    (frame, n mod 2^k, p mod 2^k)
        ->
    m_k

exactly.

Experiment 648 tests whether the full p-residue can be
collapsed to only its 2-adic displacement from the
frame anchor.

Define:

FRAME A:

    p0 = 3
    q0 = -3
    c  = 9

    A = p-p0
    B = q-q0.

FRAME B:

    p0 = -1
    q0 = 3
    c  = 3

    A = p-p0
    B = q-q0.

Let:

    tp = v2(p-p0)
    tq = v2(q-q0)
    w  = v2(n+c).

The exact residual result suggests:

    m = min(tp,tq).

Experiment 648 tests the stronger identity:

    m = min(tp,w).

Since:

    w = m
        when tp != tq,

and:

    w > m
        when tp = tq,

this would imply:

    m = min(tp,w).

Then the complete depth becomes:

    depth
      =
    m + [tp=tq]

      =
    min(
        tp+1,
        w
    ).

Therefore the whole hierarchy could collapse to:

    FRAME A:

        depth =
        min(
            v2(p-3)+1,
            v2(n+9)
        ).

    FRAME B:

        depth =
        min(
            v2(p+1)+1,
            v2(n+3)
        ).

This would be a major reduction:

    full factor pair (p,q)
        ->
    one factor-side valuation
        +
    one n-only valuation.

The remaining information would no longer be
the full 2-adic factor branch.

It would be only the distance of ONE prime factor
from its frame anchor:

    p = 3  (FRAME A)

or:

    p = -1 (FRAME B).

If this experiment passes, the next question becomes
whether even that one valuation can be replaced by a
single factor bit, a prime-residue class, or another
smaller branch invariant.
"""
    )

    print()

    print("=" * 90)
    print(
        "EXPERIMENT 648 FINISHED"
    )
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
