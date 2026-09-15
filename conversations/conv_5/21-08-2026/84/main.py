#!/usr/bin/env python3

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from math import isqrt


# =============================================================================
# EXPERIMENT 646
# =============================================================================
#
# EXACT 2-ADIC FACTOR-LIFT FIBER
#
# We now distinguish two different statements:
#
#   IMAGE STATEMENT:
#
#       n == -9 (mod 2^t)   FRAME A
#       n == -3 (mod 2^t)   FRAME B
#
#   means that there EXISTS a compatible factor lift with m >= t.
#
# It does NOT mean that the particular factorization (p,q) of n
# has m >= t.
#
# The question is now:
#
#   For fixed n mod 2^k,
#
#       F(n,k) = { possible m }
#
# over ALL compatible odd factor residues
#
#       p*q == n (mod 2^k).
#
# Candidate theorem:
#
#       F(n,k)
#         =
#       {1,2,...,min(v2(n+c),k)}
#
# with
#
#       c=9  FRAME A
#       c=3  FRAME B.
#
# If this is exact, then we have completely characterized the
# information lost by the projection
#
#       (p,q) -> n mod 2^k.
#
# ==============================================================================


INF = 10**9


# =============================================================================
# 2-ADIC UTILITIES
# =============================================================================

def v2(x: int) -> int:
    if x == 0:
        return INF

    x = abs(x)
    return (x & -x).bit_length() - 1


def v2_capped(x: int, k: int) -> int:
    """
    v2(x), capped at k.

    Thus:
        x == 0 mod 2^k -> k
    """
    if x == 0:
        return k

    x = abs(x)

    if x == 0:
        return k

    z = (x & -x).bit_length() - 1
    return min(z, k)


# =============================================================================
# FRAME / TARGET RESIDUES
# =============================================================================

def frame_from_residue(n: int) -> str:
    r = n & 3

    if r == 3:
        return "A"

    if r == 1:
        return "B"

    raise ValueError(
        f"Expected odd residue, got n mod 4 = {r}"
    )


def targets(frame: str) -> tuple[int, int, int]:
    """
    Return:

        p target residue
        q target residue
        c in n+c
    """

    if frame == "A":
        return 3, -3, 9

    return -1, 3, 3


# =============================================================================
# EXACT FIBER FOR ONE n-RESIDUE
# =============================================================================

def exact_fiber(
    n: int,
    k: int,
) -> tuple[str, set[int]]:

    M = 1 << k

    frame = frame_from_residue(n)

    rp, rq, c = targets(frame)

    rp %= M
    rq %= M

    fiber: set[int] = set()

    # Every odd residue p is a unit modulo 2^k.
    for p in range(1, M, 2):

        q = (
            n
            * pow(p, -1, M)
        ) % M

        ma = v2_capped(
            (p - rp) % M,
            k,
        )

        mb = v2_capped(
            (q - rq) % M,
            k,
        )

        m = min(
            ma,
            mb,
        )

        fiber.add(m)

    return frame, fiber


# =============================================================================
# EXPECTED FIBER
# =============================================================================

def expected_fiber(
    n: int,
    k: int,
) -> set[int]:

    frame = frame_from_residue(n)

    if frame == "A":
        c = 9
    else:
        c = 3

    t = v2_capped(
        n + c,
        k,
    )

    return set(
        range(
            1,
            t + 1,
        )
    )


# =============================================================================
# TEST 0: EXACT FIBER THEOREM
# =============================================================================

def test_exact_fibers(k: int) -> int:

    print("=" * 90)
    print(
        f"TEST 0: EXACT FIBER THEOREM k={k}"
    )
    print("=" * 90)

    M = 1 << k

    failures = 0
    checked = 0

    widths = defaultdict(int)

    for n in range(
        1,
        M,
        2,
    ):

        checked += 1

        frame, actual = exact_fiber(
            n,
            k,
        )

        expected = expected_fiber(
            n,
            k,
        )

        widths[
            (
                frame,
                tuple(sorted(actual)),
            )
        ] += 1

        if actual != expected:

            failures += 1

            if failures <= 30:

                print(
                    f"    mismatch "
                    f"frame={frame} "
                    f"n={n} "
                    f"actual={sorted(actual)} "
                    f"expected={sorted(expected)}"
                )

    print(
        f"checked={checked}"
    )

    print(
        f"failures={failures}"
    )

    print()

    return failures


# =============================================================================
# TEST 1: FIBER IS ALWAYS AN INTERVAL
# =============================================================================

def test_interval_structure(k: int) -> int:

    print("=" * 90)
    print(
        f"TEST 1: FIBER INTERVAL STRUCTURE k={k}"
    )
    print("=" * 90)

    M = 1 << k

    failures = 0

    for n in range(
        1,
        M,
        2,
    ):

        _, fiber = exact_fiber(
            n,
            k,
        )

        if not fiber:
            failures += 1
            continue

        lo = min(fiber)
        hi = max(fiber)

        expected = set(
            range(
                lo,
                hi + 1,
            )
        )

        if fiber != expected:

            failures += 1

            if failures <= 20:

                print(
                    f"    non-interval "
                    f"n={n} "
                    f"fiber={sorted(fiber)}"
                )

    print(
        f"failures={failures}"
    )

    print()

    return failures


# =============================================================================
# TEST 2: MAXIMUM POSSIBLE m
# =============================================================================

def test_maximum_fiber(k: int) -> int:

    print("=" * 90)
    print(
        f"TEST 2: MAXIMUM FIBER VALUE k={k}"
    )
    print("=" * 90)

    M = 1 << k

    failures = 0

    for n in range(
        1,
        M,
        2,
    ):

        frame, fiber = exact_fiber(
            n,
            k,
        )

        if frame == "A":
            c = 9
        else:
            c = 3

        expected_max = v2_capped(
            n + c,
            k,
        )

        actual_max = max(fiber)

        if actual_max != expected_max:

            failures += 1

            if failures <= 30:

                print(
                    f"    mismatch "
                    f"frame={frame} "
                    f"n={n} "
                    f"actual_max={actual_max} "
                    f"expected_max={expected_max}"
                )

    print(
        f"failures={failures}"
    )

    print()

    return failures


# =============================================================================
# TEST 3: EVERY m BELOW THE MAXIMUM EXISTS
# =============================================================================

def test_every_level_exists(k: int) -> int:

    print("=" * 90)
    print(
        f"TEST 3: EVERY m <= MAX EXISTS k={k}"
    )
    print("=" * 90)

    M = 1 << k

    failures = 0

    for n in range(
        1,
        M,
        2,
    ):

        frame, fiber = exact_fiber(
            n,
            k,
        )

        if frame == "A":
            c = 9
        else:
            c = 3

        vmax = v2_capped(
            n + c,
            k,
        )

        for m in range(
            1,
            vmax + 1,
        ):

            if m not in fiber:

                failures += 1

                if failures <= 30:

                    print(
                        f"    missing "
                        f"frame={frame} "
                        f"n={n} "
                        f"m={m} "
                        f"vmax={vmax}"
                    )

    print(
        f"failures={failures}"
    )

    print()

    return failures


# =============================================================================
# TEST 4: CONSTRUCTIVE LIFT
# =============================================================================

def constructive_lift(
    n: int,
    frame: str,
    m: int,
) -> tuple[int, int]:

    M = 1 << m

    if frame == "A":

        # Need:
        #
        #   p == 3  mod 2^m
        #   q == -3 mod 2^m
        #
        p = 3 % M

        q = (
            n
            * pow(
                p,
                -1,
                M,
            )
        ) % M

    else:

        # Need:
        #
        #   p == -1 mod 2^m
        #   q == 3 mod 2^m
        #
        p = (-1) % M

        q = (
            n
            * pow(
                p,
                -1,
                M,
            )
        ) % M

    return p, q


def test_constructive_lifts(
    k: int,
) -> int:

    print("=" * 90)
    print(
        f"TEST 4: CONSTRUCTIVE LIFT k={k}"
    )
    print("=" * 90)

    M = 1 << k

    failures = 0

    for n in range(
        1,
        M,
        2,
    ):

        frame = frame_from_residue(
            n
        )

        if frame == "A":
            c = 9
        else:
            c = 3

        vmax = v2_capped(
            n + c,
            k,
        )

        for m in range(
            1,
            vmax + 1,
        ):

            p, q = constructive_lift(
                n,
                frame,
                m,
            )

            if (
                p * q
            ) % (1 << m) != (
                n % (1 << m)
            ):

                failures += 1

                if failures <= 20:

                    print(
                        f"    product mismatch "
                        f"frame={frame} "
                        f"n={n} "
                        f"m={m} "
                        f"p={p} "
                        f"q={q}"
                    )

                continue

            rp, rq, _ = targets(
                frame
            )

            if (
                p - rp
            ) % (1 << m) != 0:

                failures += 1

                if failures <= 20:

                    print(
                        f"    p-lift mismatch "
                        f"frame={frame} "
                        f"n={n} "
                        f"m={m}"
                    )

            if (
                q - rq
            ) % (1 << m) != 0:

                failures += 1

                if failures <= 20:

                    print(
                        f"    q-lift mismatch "
                        f"frame={frame} "
                        f"n={n} "
                        f"m={m}"
                    )

    print(
        f"failures={failures}"
    )

    print()

    return failures


# =============================================================================
# TEST 5: ACTUAL SEMIPRIME FIBER VS THEORETICAL FIBER
# =============================================================================

def sieve(
    limit: int,
) -> list[int]:

    a = bytearray(
        b"\x01"
    ) * (
        limit + 1
    )

    a[0] = a[1] = 0

    for p in range(
        2,
        isqrt(limit) + 1,
    ):

        if a[p]:

            a[
                p * p:
                limit + 1:
                p
            ] = b"\x00" * (
                (
                    limit
                    - p * p
                )
                // p
                + 1
            )

    return [
        p
        for p in range(
            3,
            limit + 1,
            2
        )
        if a[p]
    ]


def build_semiprimes(
    primes: list[int],
) -> list[tuple[int, int, int]]:

    out = []

    for i, p in enumerate(
        primes
    ):

        for q in primes[i:]:

            n = p * q

            if n % 4 == 3:
                frame = "A"
                A = p - 3
                B = q + 3
                c = 9
            else:
                frame = "B"
                A = p + 1
                B = q - 3
                c = 3

            m = min(
                v2(A),
                v2(B),
            )

            out.append(
                (
                    n,
                    frame,
                    m,
                )
            )

    return out


def test_actual_projection(
    limit_prime: int = 6000,
    k: int = 10,
) -> int:

    print("=" * 90)
    print(
        "TEST 5: ACTUAL SEMIPRIME PROJECTION"
    )
    print("=" * 90)

    primes = sieve(
        limit_prime
    )

    states = build_semiprimes(
        primes
    )

    M = 1 << k

    buckets = defaultdict(
        set
    )

    actual_state_count = 0

    for n, frame, m in states:

        buckets[
            (
                frame,
                n % M,
            )
        ].add(m)

        actual_state_count += 1

    failures = 0

    ambiguous = 0

    for (
        frame,
        residue
    ), observed in buckets.items():

        vmax = v2_capped(
            residue
            +
            (
                9
                if frame == "A"
                else 3
            ),
            k,
        )

        theoretical = set(
            range(
                1,
                vmax + 1,
            )
        )

        # Actual semiprimes do not have to realize every
        # p-adic lift, so observed is expected to be a subset.
        #
        # But it may never contain an m > vmax.
        if not observed.issubset(
            theoretical
        ):

            failures += 1

            print(
                f"    impossible observed m "
                f"frame={frame} "
                f"residue={residue} "
                f"observed={sorted(observed)} "
                f"theoretical={sorted(theoretical)}"
            )

        if len(observed) > 1:
            ambiguous += 1

    print(
        f"states={actual_state_count}"
    )

    print(
        f"residue buckets={len(buckets)}"
    )

    print(
        f"ambiguous actual buckets={ambiguous}"
    )

    print(
        f"failures={failures}"
    )

    print()

    return failures


# =============================================================================
# EXAMPLES
# =============================================================================

def print_examples() -> None:

    print("=" * 90)
    print(
        "EXAMPLES"
    )
    print("=" * 90)
    print()

    examples = [
        (9,  "B"),
        (15, "A"),
        (21, "B"),
        (33, "B"),
        (39, "A"),
        (57, "B"),
        (69, "B"),
        (77, "B"),
        (87, "A"),
        (93, "B"),
        (111, "A"),
        (141, "B"),
        (183, "A"),
        (213, "B"),
    ]

    for n, frame in examples:

        c = (
            9
            if frame == "A"
            else 3
        )

        vmax = v2(
            n + c
        )

        theoretical = set(
            range(
                1,
                vmax + 1,
            )
        )

        print(
            f"n={n:<5} "
            f"frame={frame} "
            f"v2(n+c)={vmax}"
        )

        print(
            f"    theoretical fiber="
            f"{sorted(theoretical)}"
        )

        for m in range(
            1,
            min(vmax, 3) + 1,
        ):

            p, q = constructive_lift(
                n,
                frame,
                m,
            )

            print(
                f"    m={m:<2} "
                f"lift=("
                f"{p} mod 2^{m}, "
                f"{q} mod 2^{m})"
            )

        print()


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 90)
    print(
        "EXPERIMENT 646 START"
    )
    print("=" * 90)
    print()

    total_failures = 0

    # Small exact full-residue tests.
    for k in (
        6,
        8,
        10,
        12,
    ):

        total_failures += (
            test_exact_fibers(k)
        )

        total_failures += (
            test_interval_structure(k)
        )

        total_failures += (
            test_maximum_fiber(k)
        )

        total_failures += (
            test_every_level_exists(k)
        )

    # Constructive proof check.
    for k in (
        8,
        10,
        12,
    ):

        total_failures += (
            test_constructive_lifts(k)
        )

    # Actual prime-product projection.
    total_failures += (
        test_actual_projection(
            limit_prime=6000,
            k=10,
        )
    )

    print_examples()

    print("=" * 90)
    print(
        "FINAL STRUCTURAL SUMMARY"
    )
    print("=" * 90)

    print(
r"""
The previous experiment failed because it attempted to
identify the actual factorization with the existence of
a compatible 2-adic factor lift.

The correct distinction is:

    ACTUAL FACTORIZATION:

        m = min(v2(A),v2(B))

    2-ADIC FIBER:

        F_k(n)
          =
        { m :
          some odd factor residues p,q satisfy
          p*q == n (mod 2^k)
          and produce this m
        }.

The candidate exact fiber theorem is:

    F_k(n)
      =
    {1,2,...,T}

where

    T = min(
        v2(n+9),
        k
    )

for FRAME A,

and

    T = min(
        v2(n+3),
        k
    )

for FRAME B.

Thus:

    max F_k(n)
        =
    min(v2(n+c), k),

but the actual factorization may realize
a strictly smaller m.

This explains the observed examples:

    n=93, FRAME B

        v2(n+3)=5

    but the actual factorization

        93 = 3*31

    has

        m=2.

The residue fiber can nevertheless contain

        m=1,2,3,4,5.

Therefore n determines the MAXIMUM compatible
2-adic residual depth, not necessarily the depth of
the actual prime factorization.

If TEST 0/1/2/3 all pass, the information-loss
problem is completely characterized:

    n mod 2^k
        ->
    maximum compatible m

but

    actual factorization
        ->
    chooses one element of the fiber.

The remaining question would then be whether the
prime restriction itself introduces a useful
selection rule inside that fiber.
"""
    )

    print()

    print("=" * 90)
    print(
        "EXPERIMENT 646 FINISHED"
    )
    print("=" * 90)

    print(
        f"TOTAL FAILURES={total_failures}"
    )

    if total_failures == 0:

        print(
            "STATUS=ALL FIBER TESTS PASSED"
        )

    else:

        print(
            "STATUS=FIBER COUNTEREXAMPLES FOUND"
        )


if __name__ == "__main__":
    main()
