#!/usr/bin/env python3

from __future__ import annotations

from collections import defaultdict
from math import isqrt


INF = 10**9


# =============================================================================
# EXPERIMENT 650
#
# SCALING TEST FOR THE FINITE BRANCH-VALUATION COLLAPSE
#
# Experiment 649 found:
#
#     cap=9 -> exact
#
# on the tested prime population.
#
# This experiment asks whether that remains true as the prime
# search range increases.
#
# Exact law:
#
#     FRAME A:
#         p0 = 3
#         c  = 9
#
#     FRAME B:
#         p0 = -1
#         c  = 3
#
#     tp = v2(p-p0)
#     w  = v2(n+c)
#
#     depth = min(tp+1, w)
#
# We test whether:
#
#     (frame, w, min(tp, CAP))
#
# determines depth exactly.
# =============================================================================


def v2(x: int) -> int:
    if x == 0:
        return INF

    x = abs(x)
    return (x & -x).bit_length() - 1


# =============================================================================
# PRIME SIEVE
# =============================================================================

def sieve(limit: int) -> list[int]:
    if limit < 3:
        return []

    flags = bytearray(b"\x01") * (limit + 1)
    flags[0] = 0
    flags[1] = 0

    for p in range(2, isqrt(limit) + 1):
        if flags[p]:
            start = p * p
            count = ((limit - start) // p) + 1
            flags[start::p] = b"\x00" * count

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
        f"Odd semiprime has invalid n mod 4={r}"
    )


def frame_params(frame: str) -> tuple[int, int, int]:
    if frame == "A":
        # A = p - 3
        # B = q + 3
        # n + 9
        return 3, -3, 9

    # FRAME B
    # A = p + 1
    # B = q - 3
    # n + 3
    return -1, 3, 3


# =============================================================================
# BASIC STATE
# =============================================================================

def state_for(
    p: int,
    q: int,
) -> tuple[str, int, int, int, int, int]:
    n = p * q
    frame = frame_from_n(n)

    p0, q0, c = frame_params(frame)

    tp = v2(p - p0)
    tq = v2(q - q0)
    w = v2(n + c)

    m = min(tp, tq)
    depth = min(tp + 1, w)

    # Exact Experiment-648 identity.
    expected_depth = m + int(tp == tq)

    if depth != expected_depth:
        raise AssertionError(
            "Experiment 648 baseline failed:\n"
            f"p={p} q={q} n={n} frame={frame}\n"
            f"tp={tp} tq={tq} w={w}\n"
            f"depth={depth} expected={expected_depth}"
        )

    return (
        frame,
        tp,
        tq,
        w,
        m,
        depth,
    )


# =============================================================================
# CAP SCAN
# =============================================================================

def scan_caps(
    limit: int,
    caps: list[int],
) -> dict[int, tuple[int | None, int, int]]:

    primes = sieve(limit)

    print()
    print("=" * 90)
    print(f"PRIME LIMIT = {limit}")
    print("=" * 90)
    print(
        f"odd primes={len(primes)}"
    )

    buckets: dict[int, dict[tuple, set[int]]] = {
        cap: defaultdict(set)
        for cap in caps
    }

    first_collision: dict[int, tuple | None] = {
        cap: None
        for cap in caps
    }

    states = 0
    max_tp = 0
    max_w = 0

    for i, p in enumerate(primes):

        for q in primes[i:]:

            (
                frame,
                tp,
                tq,
                w,
                m,
                depth,
            ) = state_for(p, q)

            states += 1
            max_tp = max(max_tp, tp)
            max_w = max(max_w, w)

            for cap in caps:

                signature = (
                    frame,
                    w,
                    min(tp, cap),
                )

                bucket = buckets[cap][signature]
                bucket.add(depth)

                if (
                    len(bucket) > 1
                    and first_collision[cap] is None
                ):
                    first_collision[cap] = (
                        p,
                        q,
                        p * q,
                        frame,
                        tp,
                        tq,
                        w,
                        min(tp, cap),
                        tuple(sorted(bucket)),
                    )

    print(
        f"states={states}"
    )
    print(
        f"max observed tp={max_tp}"
    )
    print(
        f"max observed w={max_w}"
    )

    print()
    print("CAP RESULTS")

    results = {}

    for cap in caps:

        ambiguous = sum(
            1
            for depths in buckets[cap].values()
            if len(depths) > 1
        )

        collision = first_collision[cap]

        print(
            f"    cap={cap:2d}"
            f" signatures={len(buckets[cap]):6d}"
            f" ambiguous={ambiguous:6d}"
        )

        if collision is not None:
            (
                p,
                q,
                n,
                frame,
                tp,
                tq,
                w,
                truncated_tp,
                depths,
            ) = collision

            print(
                "        first collision:"
            )
            print(
                f"            n={n}"
                f" p={p}"
                f" q={q}"
                f" frame={frame}"
            )
            print(
                f"            tp={tp}"
                f" tq={tq}"
                f" w={w}"
            )
            print(
                f"            truncated_tp={truncated_tp}"
            )
            print(
                f"            depths={depths}"
            )

        results[cap] = (
            collision,
            len(buckets[cap]),
            ambiguous,
        )

    return results


# =============================================================================
# TEST 1
#
# Search specifically for the first counterexample to CAP=9.
# =============================================================================

def test_cap9_scaling() -> None:

    print()
    print("=" * 90)
    print("TEST 1: DOES CAP=9 SURVIVE PRIME-RANGE SCALING?")
    print("=" * 90)

    limits = [
        100,
        300,
        1000,
        3000,
        6000,
        10000,
        20000,
        30000,
    ]

    found = False

    for limit in limits:

        results = scan_caps(
            limit,
            [7, 8, 9, 10, 11, 12],
        )

        collision = results[9][0]

        if collision is not None:

            found = True

            (
                p,
                q,
                n,
                frame,
                tp,
                tq,
                w,
                truncated_tp,
                depths,
            ) = collision

            print()
            print(
                "CAP=9 FIRST COUNTEREXAMPLE"
            )
            print(
                f"    prime_limit={limit}"
            )
            print(
                f"    n={n}"
                f" p={p}"
                f" q={q}"
                f" frame={frame}"
            )
            print(
                f"    tp={tp}"
                f" tq={tq}"
                f" w={w}"
            )
            print(
                f"    truncated_tp={truncated_tp}"
            )
            print(
                f"    depths={depths}"
            )

            break

    if not found:
        print()
        print(
            "CAP=9 remained exact for every tested prime limit."
        )


# =============================================================================
# TEST 2
#
# Smallest exact CAP at each prime limit.
#
# We stop at CAP=16.
# =============================================================================

def test_minimal_cap() -> None:

    print()
    print("=" * 90)
    print("TEST 2: MINIMAL EXACT CAP BY PRIME LIMIT")
    print("=" * 90)

    limits = [
        100,
        300,
        1000,
        3000,
        6000,
        10000,
    ]

    caps = list(range(0, 17))

    for limit in limits:

        results = scan_caps(
            limit,
            caps,
        )

        minimal = None

        for cap in caps:

            if results[cap][0] is None:
                minimal = cap
                break

        print(
            f"    limit={limit:6d}"
            f" minimal_exact_cap={minimal}"
        )


# =============================================================================
# TEST 3
#
# Relative branch state without the actual tp value.
#
# This is intentionally expected to fail:
#
#     tp < w
#
# tells us that depth=tp+1, but does not tell us tp itself.
# =============================================================================

def test_relative_state(limit: int) -> int:

    print()
    print("=" * 90)
    print("TEST 3: RELATIVE BRANCH STATE")
    print("=" * 90)

    primes = sieve(limit)

    buckets: dict[tuple, set[int]] = defaultdict(set)

    states = 0

    for i, p in enumerate(primes):

        for q in primes[i:]:

            (
                frame,
                tp,
                tq,
                w,
                m,
                depth,
            ) = state_for(p, q)

            relation = (
                "tp<w"
                if tp < w
                else "tp>=w"
            )

            signature = (
                frame,
                w,
                relation,
            )

            buckets[signature].add(depth)
            states += 1

    ambiguous = {
        signature: depths
        for signature, depths in buckets.items()
        if len(depths) > 1
    }

    print(
        f"states={states}"
    )
    print(
        f"signatures={len(buckets)}"
    )
    print(
        f"ambiguous={len(ambiguous)}"
    )

    for signature, depths in list(
        ambiguous.items()
    )[:10]:

        print(
            f"    {signature} -> "
            f"{sorted(depths)}"
        )

    return len(ambiguous)


# =============================================================================
# TEST 4
#
# Explicit first-differing-bit interpretation.
# =============================================================================

def test_first_difference(limit: int) -> int:

    print()
    print("=" * 90)
    print("TEST 4: FIRST-DIFFERING-BIT INTERPRETATION")
    print("=" * 90)

    primes = sieve(limit)

    failures = 0
    checked = 0

    for i, p in enumerate(primes):

        for q in primes[i:]:

            (
                frame,
                tp,
                tq,
                w,
                m,
                depth,
            ) = state_for(p, q)

            p0, _, c = frame_params(frame)

            direct_tp = v2(p - p0)

            predicted = min(
                direct_tp + 1,
                v2(p * q + c),
            )

            checked += 1

            if predicted != depth:

                failures += 1

                if failures <= 10:
                    print(
                        "    mismatch:"
                        f" n={p*q}"
                        f" p={p}"
                        f" q={q}"
                        f" frame={frame}"
                        f" tp={direct_tp}"
                        f" w={v2(p*q+c)}"
                        f" predicted={predicted}"
                        f" depth={depth}"
                    )

    print(
        f"checked={checked}"
    )
    print(
        f"failures={failures}"
    )

    return failures


# =============================================================================
# EXAMPLES
# =============================================================================

def print_examples() -> None:

    print()
    print("=" * 90)
    print("EXAMPLES")
    print("=" * 90)

    examples = [
        (9, 3, 3),
        (15, 3, 5),
        (21, 3, 7),
        (33, 3, 11),
        (39, 3, 13),
        (57, 3, 19),
        (69, 3, 23),
        (77, 7, 11),
        (87, 3, 29),
        (93, 3, 31),
        (111, 3, 37),
        (141, 3, 47),
        (183, 3, 61),
        (213, 3, 71),
        (1047, 3, 349),
        (4135, 5, 827),
        (40999, 7, 5857),
    ]

    for n, p, q in examples:

        (
            frame,
            tp,
            tq,
            w,
            m,
            depth,
        ) = state_for(p, q)

        p0, q0, c = frame_params(frame)

        print()
        print(
            f"n={n}"
            f" p={p}"
            f" q={q}"
            f" frame={frame}"
        )

        print(
            f"    p0={p0}"
            f" q0={q0}"
            f" c={c}"
        )

        print(
            f"    tp=v2(p-p0)={tp}"
        )
        print(
            f"    tq=v2(q-q0)={tq}"
        )
        print(
            f"    w=v2(n+c)={w}"
        )
        print(
            f"    m={m}"
        )
        print(
            f"    depth={depth}"
        )
        print(
            f"    tp clipped at 9={min(tp, 9)}"
        )


# =============================================================================
# MAIN
# =============================================================================

def main() -> None:

    print("=" * 90)
    print("EXPERIMENT 650 START")
    print("=" * 90)

    test_cap9_scaling()

    test_minimal_cap()

    relative_failures = test_relative_state(
        6000
    )

    bit_failures = test_first_difference(
        6000
    )

    print_examples()

    print()
    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)

    print(
        """
Experiment 649 established the exact compression:

    depth =
        min(
            v2(p-p0)+1,
            v2(n+c)
        )

with:

    FRAME A:
        p0=3
        c=9

    FRAME B:
        p0=-1
        c=3.

Experiment 649 also found:

    cap=9
        ambiguous=0

for the then-current prime population.

That result is not automatically a theorem.

Experiment 650 therefore tests whether a fixed finite cap
survives when the prime search range grows.

The tested signature is:

    (
        frame,
        w=v2(n+c),
        min(v2(p-p0), cap)
    )

and the target is the exact depth.

There are two possible outcomes.

FINITE-STATE RESULT:

    A bounded cap remains exact while the prime range
    grows.

Then the factor branch can genuinely be compressed to a
finite 2-adic state.

SCALING RESULT:

    The minimum exact cap increases with the prime range.

Then cap=9 was only a finite-domain artifact and the
branch valuation is genuinely unbounded.

The exact uncompressed law remains:

    depth =
        min(
            v2(p-p0)+1,
            v2(n+c)
        ).

The most important datum is the FIRST counterexample to
cap=9.
"""
    )

    print()
    print("=" * 90)
    print("EXPERIMENT 650 FINISHED")
    print("=" * 90)

    total_failures = (
        relative_failures
        + bit_failures
    )

    print(
        f"TOTAL DIAGNOSTIC FAILURES={total_failures}"
    )


if __name__ == "__main__":
    main()