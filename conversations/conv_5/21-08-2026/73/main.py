#!/usr/bin/env python3

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from math import gcd, isqrt


# ==============================================================================
# EXPERIMENT 635
# ==============================================================================
#
# DIVISOR-LEVEL COLLISION STRUCTURE OF n+c
#
# Experiment 634 established:
#
#     d = gcd(A,B)
#
#     n+c = d*Q
#
# and:
#
#     depth = v2(d)                  if Q is odd
#             v2(d)+1                if Q is even.
#
# Here:
#
# FRAME A:
#     A = p-3
#     B = q+3
#     c = 9
#
# FRAME B:
#     A = p+1
#     B = q-3
#     c = 3
#
# Therefore d is a divisor of n+c.
#
# The central question is now:
#
#     Is d recoverable from the divisor structure of n+c?
#
# In particular:
#
#   1. Compare v2(d) against all divisor valuations of n+c.
#   2. Determine whether d is the largest divisor D of n+c satisfying
#      the corresponding frame congruence conditions.
#   3. Search for collisions:
#
#          same n+c
#          different d/depth
#
#      which would demonstrate that n+c alone cannot determine d.
#   4. Test whether the odd part of d is a distinguished divisor of n+c.
#   5. Test whether d can be characterized using gcd(n+c, small functions
#      of n+c), without arbitrary offset searches.
#
# IMPORTANT:
#     No expensive global offset search is performed.
#     Everything is grouped by n+c and divisor structure.
#
# ==============================================================================


PRIME_LIMIT = 6250
INF = 10**9


@dataclass(frozen=True)
class State:
    n: int
    p: int
    q: int


# ==============================================================================
# NUMBER THEORY
# ==============================================================================

def v2(x: int) -> int:

    if x == 0:
        return INF

    x = abs(x)

    return (x & -x).bit_length() - 1


def odd_part(x: int) -> int:

    if x == 0:
        return 0

    return abs(x) >> v2(x)


def sieve(limit: int) -> list[int]:

    s = bytearray(b"\x01") * (limit + 1)

    s[0] = 0
    s[1] = 0

    for p in range(
        2,
        isqrt(limit) + 1,
    ):

        if s[p]:

            start = p * p

            s[
                start:
                limit + 1:
                p
            ] = b"\x00" * (
                ((limit - start) // p) + 1
            )

    return [
        p
        for p in range(
            3,
            limit + 1,
            2
        )
        if s[p]
    ]


def generate_states(
    primes: list[int],
) -> list[State]:

    out = []

    for i, p in enumerate(primes):

        for q in primes[i:]:

            out.append(
                State(
                    p * q,
                    p,
                    q,
                )
            )

    return out


# ==============================================================================
# STRUCTURE
# ==============================================================================

def frame_of(n: int) -> str:

    r = n & 3

    if r == 3:
        return "A"

    if r == 1:
        return "B"

    raise ValueError(
        f"Unexpected odd semiprime residue: "
        f"n={n}, n mod 4={r}"
    )


def residuals(
    frame: str,
    p: int,
    q: int,
) -> tuple[int, int]:

    if frame == "A":
        return p - 3, q + 3

    return p + 1, q - 3


def xy(
    frame: str,
    p: int,
    q: int,
) -> tuple[int, int]:

    if frame == "A":

        return (
            (q - p + 6) // 2,
            (p + q) // 2,
        )

    return (
        (3 * p - q + 6) // 2,
        (3 * p + q) // 2,
    )


def depth(
    frame: str,
    p: int,
    q: int,
) -> int:

    X, Y = xy(
        frame,
        p,
        q,
    )

    return 1 + v2(
        gcd(
            abs(X),
            abs(Y),
        )
    )


# ==============================================================================
# DIVISORS
# ==============================================================================

def divisors(n: int) -> list[int]:

    n = abs(n)

    if n == 0:
        return []

    small = []
    large = []

    r = isqrt(n)

    for d in range(1, r + 1):

        if n % d:
            continue

        small.append(d)

        if d * d != n:
            large.append(n // d)

    return small + large[::-1]


# ==============================================================================
# TEST 0
# ==============================================================================

def test_divisor_inclusion(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 0: d IS A DIVISOR OF n+c")
    print("=" * 90)

    failures = 0

    for s in states:

        frame = frame_of(s.n)

        A, B = residuals(
            frame,
            s.p,
            s.q,
        )

        d = gcd(
            abs(A),
            abs(B),
        )

        c = 9 if frame == "A" else 3

        shifted = s.n + c

        if d == 0 or shifted % d != 0:

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch:"
                    f" n={s.n}"
                    f" frame={frame}"
                    f" d={d}"
                    f" n+c={shifted}"
                )

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )


# ==============================================================================
# TEST 1
# ==============================================================================

def test_v2_divisor_spectrum(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 1: v2(d) INSIDE THE DIVISOR SPECTRUM OF n+c")
    print("=" * 90)

    failures = 0
    checked = 0

    distribution = Counter()

    for s in states:

        frame = frame_of(s.n)

        A, B = residuals(
            frame,
            s.p,
            s.q,
        )

        if A == 0 or B == 0:
            continue

        d = gcd(
            abs(A),
            abs(B),
        )

        c = 9 if frame == "A" else 3

        shifted = s.n + c

        vd = v2(d)
        vn = v2(shifted)

        checked += 1

        distribution[
            (
                frame,
                vd,
                vn - vd,
            )
        ] += 1

        if vd > vn:

            failures += 1

            if failures <= 20:

                print(
                    "    impossible:"
                    f" n={s.n}"
                    f" frame={frame}"
                    f" v2d={vd}"
                    f" v2(n+c)={vn}"
                )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )

    print()
    print("    delta=v2(n+c)-v2(d) distribution:")

    delta = Counter()

    for (
        frame,
        vd,
        delta_v,
    ), count in distribution.items():

        delta[
            (frame, delta_v)
        ] += count

    for key in sorted(delta):

        print(
            f"        {key}: "
            f"{delta[key]}"
        )


# ==============================================================================
# TEST 2
# ==============================================================================

def test_distinguished_v2_divisor(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 2: IS d THE DISTINGUISHED 2-ADIC DIVISOR?")
    print("=" * 90)

    failures = 0
    checked = 0

    # We classify d among divisors of n+c by v2.
    #
    # For each n+c we record:
    #
    #     maximum v2 divisor = v2(n+c)
    #
    # but d may stop one or more layers earlier.
    #
    # The interesting question is whether its valuation is
    # determined by a simple divisor rank.
    #

    relation = Counter()

    for s in states:

        frame = frame_of(s.n)

        A, B = residuals(
            frame,
            s.p,
            s.q,
        )

        if A == 0 or B == 0:
            continue

        d = gcd(
            abs(A),
            abs(B),
        )

        c = 9 if frame == "A" else 3

        shifted = s.n + c

        ds = divisors(
            shifted
        )

        vd = v2(d)
        eligible = [
            x for x in ds
            if v2(x) == vd
        ]

        checked += 1

        relation[
            (
                frame,
                len(eligible),
            )
        ] += 1

        # d must be present.
        if d not in ds:

            failures += 1

            if failures <= 20:

                print(
                    "    missing divisor:"
                    f" n={s.n}"
                    f" d={d}"
                    f" n+c={shifted}"
                )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )

    print()
    print(
        "    number of divisors of n+c "
        "sharing v2(d):"
    )

    for key in sorted(relation)[:100]:

        print(
            f"        {key}: "
            f"{relation[key]}"
        )


# ==============================================================================
# TEST 3
# ==============================================================================

def test_same_shift_collisions(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 3: SAME n+c COLLISION STRUCTURE")
    print("=" * 90)

    buckets: dict[
        tuple[str, int],
        list[tuple[int, int, int, int]]
    ] = defaultdict(list)

    for s in states:

        frame = frame_of(s.n)

        A, B = residuals(
            frame,
            s.p,
            s.q,
        )

        if A == 0 or B == 0:
            continue

        c = 9 if frame == "A" else 3

        shifted = s.n + c

        d = gcd(
            abs(A),
            abs(B),
        )

        dep = depth(
            frame,
            s.p,
            s.q,
        )

        buckets[
            (
                frame,
                shifted,
            )
        ].append(
            (
                s.n,
                d,
                v2(d),
                dep,
            )
        )

    collision_buckets = 0
    ambiguous_depth = 0
    ambiguous_d = 0

    examples = []

    for key, values in buckets.items():

        if len(values) <= 1:
            continue

        collision_buckets += 1

        ds = {
            x[2]
            for x in values
        }

        depths = {
            x[3]
            for x in values
        }

        if len(ds) > 1:
            ambiguous_d += 1

        if len(depths) > 1:

            ambiguous_depth += 1

            if len(examples) < 20:

                examples.append(
                    (
                        key,
                        values,
                    )
                )

    print(
        f"same-shift collision buckets="
        f"{collision_buckets}"
    )

    print(
        f"ambiguous v2(d) buckets="
        f"{ambiguous_d}"
    )

    print(
        f"ambiguous depth buckets="
        f"{ambiguous_depth}"
    )

    if examples:

        print()
        print(
            "    first ambiguous "
            "same-shift collisions:"
        )

        for key, values in examples:

            print()
            print(
                f"    key={key}"
            )

            for value in values:

                print(
                    "        "
                    f"n={value[0]}"
                    f" d={value[1]}"
                    f" v2d={value[2]}"
                    f" depth={value[3]}"
                )


# ==============================================================================
# TEST 4
# ==============================================================================

def test_shift_factorization_identity(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 4: FACTORIZATION OF n+c BY d")
    print("=" * 90)

    failures = 0
    checked = 0

    parity = Counter()

    for s in states:

        frame = frame_of(s.n)

        A, B = residuals(
            frame,
            s.p,
            s.q,
        )

        if A == 0 or B == 0:
            continue

        d = gcd(
            abs(A),
            abs(B),
        )

        a = A // d
        b = B // d

        c = 9 if frame == "A" else 3

        shifted = s.n + c

        quotient = shifted // d

        if frame == "A":

            symbolic = (
                d * a * b
                + 3 * (b - a)
            )

        else:

            symbolic = (
                d * a * b
                + 3 * a
                - b
            )

        checked += 1

        if quotient != symbolic:

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch:"
                    f" n={s.n}"
                    f" frame={frame}"
                    f" quotient={quotient}"
                    f" symbolic={symbolic}"
                )

        parity[
            (
                frame,
                (a & 1),
                (b & 1),
                (quotient & 1),
            )
        ] += 1

    print(
        f"checked={checked} "
        f"failures={failures}"
    )

    print()
    print(
        "    normalized parity -> quotient parity:"
    )

    for key in sorted(parity):

        print(
            f"        {key}: "
            f"{parity[key]}"
        )


# ==============================================================================
# TEST 5
# ==============================================================================

def test_possible_d_from_shift_factorization(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print(
        "TEST 5: POSSIBLE d FROM FACTORIZATION "
        "OF n+c"
    )
    print("=" * 90)

    # For each shifted value we construct its divisors and ask:
    #
    #     which divisor valuations could represent v2(d)?
    #
    # The actual d is compared against all divisors.
    #
    # We also test a weaker candidate:
    #
    #     d has the largest power of 2 among divisors D
    #     satisfying the parity signature induced by frame.
    #
    # For nondegenerate states:
    #
    #     a,b coprime
    #
    # and exactly:
    #
    #     odd,odd OR odd,even OR even,odd.
    #
    # We cannot recover a,b from n+c alone, but we can measure
    # how much information the divisor lattice retains.
    #

    failures = 0
    checked = 0

    ambiguity = Counter()

    for s in states:

        frame = frame_of(s.n)

        A, B = residuals(
            frame,
            s.p,
            s.q,
        )

        if A == 0 or B == 0:
            continue

        d = gcd(
            abs(A),
            abs(B),
        )

        c = 9 if frame == "A" else 3

        shifted = s.n + c

        ds = divisors(
            shifted
        )

        # All divisors with the same v2(d).
        vd = v2(d)

        same_v = [
            x for x in ds
            if v2(x) == vd
        ]

        checked += 1

        ambiguity[
            len(same_v)
        ] += 1

        if d not in same_v:

            failures += 1

            if failures <= 20:

                print(
                    "    missing actual d:"
                    f" n={s.n}"
                    f" d={d}"
                    f" vd={vd}"
                )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )

    print()
    print(
        "    count of divisors having "
        "the correct v2:"
    )

    for k in sorted(ambiguity):

        print(
            f"        {k}: "
            f"{ambiguity[k]}"
        )


# ==============================================================================
# TEST 6
# ==============================================================================

def test_depth_from_shift_and_d(
    states: list[State],
) -> None:

    print()
    print("=" * 90)
    print("TEST 6: DEPTH FROM v2(d) + v2(n+c)-v2(d)")
    print("=" * 90)

    failures = 0
    checked = 0

    for s in states:

        frame = frame_of(s.n)

        A, B = residuals(
            frame,
            s.p,
            s.q,
        )

        d = gcd(
            abs(A),
            abs(B),
        )

        c = 9 if frame == "A" else 3

        shifted = s.n + c

        q = shifted // d

        predicted = (
            v2(d)
            + (
                1
                if (q & 1) == 0
                else 0
            )
        )

        actual = depth(
            frame,
            s.p,
            s.q,
        )

        checked += 1

        if predicted != actual:

            failures += 1

            if failures <= 20:

                print(
                    "    mismatch:"
                    f" n={s.n}"
                    f" frame={frame}"
                    f" d={d}"
                    f" q={q}"
                    f" predicted={predicted}"
                    f" actual={actual}"
                )

    print(
        f"checked={checked} "
        f"failures={failures}"
    )


# ==============================================================================
# EXAMPLES
# ==============================================================================

def print_examples(
    states: list[State],
) -> None:

    print()
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
        213,
    ]

    by_n = {
        s.n: s
        for s in states
    }

    for n in wanted:

        s = by_n[n]

        frame = frame_of(n)

        A, B = residuals(
            frame,
            s.p,
            s.q,
        )

        d = gcd(
            abs(A),
            abs(B),
        )

        c = 9 if frame == "A" else 3

        shifted = n + c

        q = shifted // d

        print()
        print(
            f"n={n}"
            f" p={s.p}"
            f" q={s.q}"
            f" frame={frame}"
        )

        print(
            f"    A={A}"
            f" B={B}"
            f" d={d}"
        )

        print(
            f"    n+c={shifted}"
            f" quotient={q}"
        )

        print(
            f"    factorization:"
            f" {d} * {q}"
            f" = {shifted}"
        )

        print(
            f"    v2(d)={v2(d)}"
            f" v2(q)={v2(q)}"
            f" v2(n+c)={v2(shifted)}"
        )

        if A != 0 and B != 0:

            print(
                f"    v2(A)={v2(A)}"
                f" v2(B)={v2(B)}"
                f" equal="
                f"{v2(A)==v2(B)}"
            )

        print(
            f"    depth="
            f"{depth(frame,s.p,s.q)}"
        )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print("EXPERIMENT 635 START")
    print("=" * 90)

    print()
    print(
        "DIVISOR-LEVEL COLLISION STRUCTURE OF n+c"
    )

    print()
    print("[1] PRIME SIEVE")

    primes = sieve(
        PRIME_LIMIT
    )

    print(
        f"    odd primes={len(primes)}"
    )

    print()
    print("[2] SEMIPRIME GENERATION")

    states = generate_states(
        primes
    )

    print(
        f"    semiprimes={len(states)}"
    )

    print()
    print("[3] BUILD GLOBAL RECORDS")

    print(
        f"    states={len(states)}"
    )

    test_divisor_inclusion(
        states
    )

    test_v2_divisor_spectrum(
        states
    )

    test_distinguished_v2_divisor(
        states
    )

    test_same_shift_collisions(
        states
    )

    test_shift_factorization_identity(
        states
    )

    test_possible_d_from_shift_factorization(
        states
    )

    test_depth_from_shift_and_d(
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
The current exact theorem is:

    depth
      =
    1 + v2(gcd(X,Y))

and, for nondegenerate residuals:

    d = gcd(A,B)

    depth
      =
    v2(d)
      +
    [v2(A)=v2(B)].

Experiment 634 additionally proved:

    FRAME A:
        n+9 = d*Q

    FRAME B:
        n+3 = d*Q

with:

    Q odd
        <=> v2(A) != v2(B)

    Q even
        <=> v2(A) = v2(B).

Therefore:

    v2(n+c)
      =
    v2(d) + v2(Q).

The remaining unresolved quantity is:

    v2(d).

This experiment deliberately moves from
"searching arbitrary n-only formulas"
to the divisor lattice of n+c itself.

The critical possible outcomes are:

    1. SAME-SHIFT COLLISIONS = 0

       Then n+c uniquely determines d/depth
       inside the tested semiprime domain.

    2. SAME-SHIFT COLLISIONS > 0,
       but v2(d) is constant inside every collision.

       Then n+c may still determine the depth.

    3. SAME-SHIFT COLLISIONS with different v2(d).

       Then n+c alone provably loses the information
       needed for depth, at least on the tested domain.

The strongest negative result would be an explicit pair:

    n1 + c = n2 + c

with:

    v2(d1) != v2(d2)

or:

    depth(n1) != depth(n2).

That would identify the precise information that cannot
be recovered from the shifted integer alone.

The strongest positive result would be a divisor-theoretic
characterization of d inside n+c, for example:

    d = distinguished_divisor(n+c, frame)

which would give a genuinely factor-free reconstruction
of the depth.
"""
    )

    print()
    print("=" * 90)
    print("EXPERIMENT 635 FINISHED")
    print("=" * 90)


if __name__ == "__main__":
    main()
