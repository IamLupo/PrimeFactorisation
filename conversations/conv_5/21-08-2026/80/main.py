#!/usr/bin/env python3

# ==============================================================================
# EXPERIMENT 642
# ==============================================================================
#
# NORMALIZED RESIDUAL PARITY THEOREM
#
# Exact target:
#
#   d = gcd(A,B)
#   A = d*a
#   B = d*b
#   gcd(a,b) = 1
#
# For nondegenerate states:
#
#   depth = v2(d) + [a,b both odd]
#
# Equivalently:
#
#   depth = v2(gcd(A,B)) + [v2(A)=v2(B)]
#
# This experiment independently verifies:
#
#   FRAME A:
#       2X = B-A
#       2Y = A+B
#
#   FRAME B:
#       2X = 3A-B
#       2Y = 3A+B
#
# and proves that after extracting gcd(A,B), only the
# normalized parity class matters for the extra valuation.
#
# No guessed n+c product identities are used.
#
# ==============================================================================

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from math import gcd


# ==============================================================================
# CONFIG
# ==============================================================================

PRIME_LIMIT = 6000
INF = 10**9


# ==============================================================================
# STATE
# ==============================================================================

@dataclass(frozen=True)
class State:
    n: int
    p: int
    q: int
    frame: str
    A: int
    B: int
    X: int
    Y: int
    depth: int


# ==============================================================================
# BASIC HELPERS
# ==============================================================================

def v2(x: int) -> int:
    if x == 0:
        return INF

    x = abs(x)
    return (x & -x).bit_length() - 1


def ggcd(a: int, b: int) -> int:
    return gcd(abs(a), abs(b))


def frame_from_n(n: int) -> str:
    r = n & 3

    if r == 3:
        return "A"

    if r == 1:
        return "B"

    raise ValueError(
        f"Unexpected odd n mod 4={r}"
    )


# ==============================================================================
# PRIME SIEVE
# ==============================================================================

def prime_sieve(limit: int) -> list[int]:

    sieve = bytearray(b"\x01") * (limit + 1)

    sieve[0] = 0
    sieve[1] = 0

    root = int(limit ** 0.5)

    for p in range(2, root + 1):

        if not sieve[p]:
            continue

        start = p * p

        sieve[start:limit + 1:p] = (
            b"\x00"
            * (
                (limit - start) // p + 1
            )
        )

    return [
        p
        for p in range(3, limit + 1, 2)
        if sieve[p]
    ]


# ==============================================================================
# FACTOR RESIDUALS
# ==============================================================================

def residuals(
    frame: str,
    p: int,
    q: int,
) -> tuple[int, int]:

    if frame == "A":
        return (
            p - 3,
            q + 3,
        )

    return (
        p + 1,
        q - 3,
    )


# ==============================================================================
# GLOBAL X,Y
# ==============================================================================

def global_xy(
    frame: str,
    A: int,
    B: int,
) -> tuple[int, int]:

    if frame == "A":

        nx = B - A
        ny = A + B

    else:

        nx = 3 * A - B
        ny = 3 * A + B

    if (nx & 1) or (ny & 1):

        raise ArithmeticError(
            f"Non-integral X,Y: "
            f"frame={frame} A={A} B={B}"
        )

    return (
        nx // 2,
        ny // 2,
    )


# ==============================================================================
# DEPTH
# ==============================================================================

def depth_from_xy(
    X: int,
    Y: int,
) -> int:

    return 1 + v2(
        ggcd(X, Y)
    )


# ==============================================================================
# BUILD STATES
# ==============================================================================

def build_states(
    primes: list[int],
) -> list[State]:

    states: list[State] = []

    for i, p in enumerate(primes):

        for q in primes[i:]:

            n = p * q

            frame = frame_from_n(n)

            A, B = residuals(
                frame,
                p,
                q,
            )

            X, Y = global_xy(
                frame,
                A,
                B,
            )

            depth = depth_from_xy(
                X,
                Y,
            )

            states.append(
                State(
                    n=n,
                    p=p,
                    q=q,
                    frame=frame,
                    A=A,
                    B=B,
                    X=X,
                    Y=Y,
                    depth=depth,
                )
            )

    return states


# ==============================================================================
# TEST 0
# ==============================================================================

def test_domain(
    states: list[State],
) -> int:

    print("=" * 90)
    print("TEST 0: ODD SEMIPRIME DOMAIN")
    print("=" * 90)

    failures = 0

    for s in states:

        if not (s.p & 1):
            failures += 1

        if not (s.q & 1):
            failures += 1

        if s.n != s.p * s.q:
            failures += 1

        if frame_from_n(s.n) != s.frame:
            failures += 1

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 1
# ==============================================================================

def test_transformed_numerators(
    states: list[State],
) -> int:

    print("=" * 90)
    print("TEST 1: EXACT TRANSFORMED NUMERATORS")
    print("=" * 90)

    failures = 0

    for s in states:

        if s.frame == "A":

            u = s.B - s.A
            v = s.A + s.B

        else:

            u = 3 * s.A - s.B
            v = 3 * s.A + s.B

        if u != 2 * s.X:
            failures += 1

        if v != 2 * s.Y:
            failures += 1

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 2
# ==============================================================================

def test_direct_depth(
    states: list[State],
) -> int:

    print("=" * 90)
    print("TEST 2: DEPTH = MINIMUM TRANSFORMED v2")
    print("=" * 90)

    failures = 0

    for s in states:

        if s.frame == "A":

            u = s.B - s.A
            v = s.A + s.B

        else:

            u = 3 * s.A - s.B
            v = 3 * s.A + s.B

        predicted = min(
            v2(u),
            v2(v),
        )

        if predicted != s.depth:
            failures += 1

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 3
# ==============================================================================

def test_gcd_normalization(
    states: list[State],
) -> int:

    print("=" * 90)
    print("TEST 3: COPRIME RESIDUAL NORMALIZATION")
    print("=" * 90)

    failures = 0
    checked = 0

    parity_counts = Counter()

    for s in states:

        A = s.A
        B = s.B

        if A == 0 or B == 0:
            continue

        checked += 1

        d = ggcd(A, B)

        if d == 0:
            failures += 1
            continue

        a = A // d
        b = B // d

        if ggcd(a, b) != 1:
            failures += 1

        pa = a & 1
        pb = b & 1

        parity_counts[(pa, pb)] += 1

        if pa == 0 and pb == 0:
            failures += 1

    print(
        f"checked={checked} "
        f"failures={failures}"
    )

    print("    parity distribution:")

    for key in (
        (0, 1),
        (1, 0),
        (1, 1),
        (0, 0),
    ):

        print(
            f"        {key}: "
            f"{parity_counts[key]}"
        )

    print()

    return failures


# ==============================================================================
# TEST 4
# ==============================================================================

def test_frame_a_normalized(
    states: list[State],
) -> int:

    print("=" * 90)
    print("TEST 4: FRAME-A NORMALIZED THEOREM")
    print("=" * 90)

    failures = 0
    checked = 0

    branch_counts = Counter()

    for s in states:

        if s.frame != "A":
            continue

        A = s.A
        B = s.B

        if A == 0 or B == 0:
            continue

        checked += 1

        d = ggcd(A, B)
        a = A // d
        b = B // d

        u = b - a
        v = a + b

        increment = min(
            v2(u),
            v2(v),
        )

        predicted = v2(d) + increment

        if predicted != s.depth:
            failures += 1

        parity = (
            a & 1,
            b & 1,
        )

        if parity in (
            (0, 1),
            (1, 0),
        ):

            branch = "opposite"
            expected_increment = 0

        elif parity == (1, 1):

            branch = "both-odd"
            expected_increment = 1

        else:

            branch = "invalid"
            expected_increment = None

        branch_counts[branch] += 1

        if (
            expected_increment is None
            or increment != expected_increment
        ):

            failures += 1

    print(
        f"checked={checked} "
        f"failures={failures}"
    )

    for branch in (
        "opposite",
        "both-odd",
        "invalid",
    ):

        print(
            f"    {branch}: "
            f"{branch_counts[branch]}"
        )

    print()

    return failures


# ==============================================================================
# TEST 5
# ==============================================================================

def test_frame_b_normalized(
    states: list[State],
) -> int:

    print("=" * 90)
    print("TEST 5: FRAME-B NORMALIZED THEOREM")
    print("=" * 90)

    failures = 0
    checked = 0

    branch_counts = Counter()

    for s in states:

        if s.frame != "B":
            continue

        A = s.A
        B = s.B

        if A == 0 or B == 0:
            continue

        checked += 1

        d = ggcd(A, B)
        a = A // d
        b = B // d

        u = 3 * a - b
        v = 3 * a + b

        increment = min(
            v2(u),
            v2(v),
        )

        predicted = v2(d) + increment

        if predicted != s.depth:
            failures += 1

        parity = (
            a & 1,
            b & 1,
        )

        if parity in (
            (0, 1),
            (1, 0),
        ):

            branch = "opposite"
            expected_increment = 0

        elif parity == (1, 1):

            branch = "both-odd"
            expected_increment = 1

        else:

            branch = "invalid"
            expected_increment = None

        branch_counts[branch] += 1

        if (
            expected_increment is None
            or increment != expected_increment
        ):

            failures += 1

    print(
        f"checked={checked} "
        f"failures={failures}"
    )

    for branch in (
        "opposite",
        "both-odd",
        "invalid",
    ):

        print(
            f"    {branch}: "
            f"{branch_counts[branch]}"
        )

    print()

    return failures


# ==============================================================================
# TEST 6
# ==============================================================================

def test_exact_gcd_formula(
    states: list[State],
) -> int:

    print("=" * 90)
    print("TEST 6: EXACT gcd(X,Y) NORMALIZED FORMULA")
    print("=" * 90)

    failures = 0

    for s in states:

        A = s.A
        B = s.B

        gx = ggcd(
            s.X,
            s.Y,
        )

        if A == 0 and B == 0:

            failures += 1
            continue

        if A == 0:

            expected = abs(B) // 2

            if gx != expected:
                failures += 1

            continue

        if B == 0:

            expected = 3 * abs(A) // 2

            if gx != expected:
                failures += 1

            continue

        d = ggcd(A, B)

        a = A // d
        b = B // d

        if (a & 1) and (b & 1):
            expected = d
        else:
            expected = d // 2

        if gx != expected:
            failures += 1

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 7
# ==============================================================================

def test_equal_v2_equivalence(
    states: list[State],
) -> int:

    print("=" * 90)
    print("TEST 7: BOTH-ODD <-> EQUAL v2")
    print("=" * 90)

    failures = 0
    checked = 0

    for s in states:

        A = s.A
        B = s.B

        if A == 0 or B == 0:
            continue

        checked += 1

        d = ggcd(A, B)

        a = A // d
        b = B // d

        both_odd = (
            (a & 1) == 1
            and
            (b & 1) == 1
        )

        equal_v2 = (
            v2(A)
            ==
            v2(B)
        )

        if both_odd != equal_v2:
            failures += 1

    print(
        f"checked={checked} "
        f"failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 8
# ==============================================================================

def test_universal_depth_formula(
    states: list[State],
) -> int:

    print("=" * 90)
    print("TEST 8: UNIVERSAL FACTOR-SIDE DEPTH FORMULA")
    print("=" * 90)

    failures = 0

    branch_counts = Counter()

    for s in states:

        A = s.A
        B = s.B

        if A == 0:

            predicted = v2(abs(B))
            branch = "A=0"

        elif B == 0:

            predicted = v2(abs(A))
            branch = "B=0"

        else:

            d = ggcd(A, B)

            if v2(A) == v2(B):

                predicted = v2(d) + 1
                branch = "equal-v2"

            else:

                predicted = v2(d)
                branch = "unequal-v2"

        branch_counts[
            (s.frame, branch)
        ] += 1

        if predicted != s.depth:
            failures += 1

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )

    for key in sorted(branch_counts):
        print(
            f"    {key}: "
            f"{branch_counts[key]}"
        )

    print()

    return failures


# ==============================================================================
# TEST 9
# ==============================================================================

def test_no_higher_bits_needed(
    states: list[State],
) -> int:

    print("=" * 90)
    print(
        "TEST 9: NORMALIZED PARITY COMPLETELY "
        "DETERMINES THE EXTRA VALUATION"
    )
    print("=" * 90)

    failures = 0
    checked = 0

    # IMPORTANT:
    #
    # This is intentionally a defaultdict(set), not a normal
    # dictionary of integers. The previous script attempted
    # ".add()" on an integer and crashed.
    #
    increment_values = defaultdict(set)

    # Also keep a plain Counter for statistics.
    parity_counts = Counter()

    for s in states:

        A = s.A
        B = s.B

        if A == 0 or B == 0:
            continue

        checked += 1

        d = ggcd(A, B)

        a = A // d
        b = B // d

        if s.frame == "A":

            u = b - a
            v = a + b

        else:

            u = 3 * a - b
            v = 3 * a + b

        increment = min(
            v2(u),
            v2(v),
        )

        parity = (
            a & 1,
            b & 1,
        )

        increment_values[
            parity
        ].add(
            increment
        )

        parity_counts[
            parity
        ] += 1

    expected_sets = {
        (0, 1): {0},
        (1, 0): {0},
        (1, 1): {1},
    }

    for parity, values in increment_values.items():

        expected = expected_sets.get(
            parity
        )

        if expected is None:

            failures += 1
            continue

        if values != expected:

            failures += 1

    print(
        f"checked={checked} "
        f"failures={failures}"
    )

    for parity in (
        (0, 1),
        (1, 0),
        (1, 1),
        (0, 0),
    ):

        values = increment_values.get(
            parity,
            set(),
        )

        print(
            f"    parity={parity} "
            f"count={parity_counts[parity]} "
            f"increment_values="
            f"{sorted(values)}"
        )

    print()

    return failures


# ==============================================================================
# TEST 10
# ==============================================================================

def test_examples(
    states: list[State],
) -> int:

    print("=" * 90)
    print("TEST 10: EXAMPLES")
    print("=" * 90)

    lookup = {
        s.n: s
        for s in states
    }

    examples = [
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
    ]

    failures = 0

    for n in examples:

        if n not in lookup:
            print(
                f"    missing example n={n}"
            )
            failures += 1
            continue

        s = lookup[n]

        A = s.A
        B = s.B

        print(
            f"n={n} "
            f"p={s.p} "
            f"q={s.q} "
            f"frame={s.frame}"
        )

        print(
            f"    A={A} "
            f"B={B} "
            f"X={s.X} "
            f"Y={s.Y}"
        )

        if A == 0:

            predicted = v2(
                abs(B)
            )

            print(
                "    branch=A=0"
            )

        elif B == 0:

            predicted = v2(
                abs(A)
            )

            print(
                "    branch=B=0"
            )

        else:

            d = ggcd(A, B)

            a = A // d
            b = B // d

            print(
                f"    d={d}"
            )

            print(
                f"    normalized="
                f"({a},{b})"
            )

            print(
                f"    parity="
                f"({a & 1},{b & 1})"
            )

            print(
                f"    v2(A)={v2(A)} "
                f"v2(B)={v2(B)}"
            )

            if (a & 1) and (b & 1):

                predicted = v2(d) + 1

                print(
                    "    branch=both-odd"
                )

            else:

                predicted = v2(d)

                print(
                    "    branch=opposite-parity"
                )

        print(
            f"    predicted_depth="
            f"{predicted} "
            f"actual_depth="
            f"{s.depth}"
        )

        print()

        if predicted != s.depth:
            failures += 1

    print(
        f"checked={len(examples)} "
        f"failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# SUMMARY
# ==============================================================================

def summary():

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)

    print(
r"""
The direct coordinate theorem is:

FRAME A:

    2X = B-A
    2Y = A+B

FRAME B:

    2X = 3A-B
    2Y = 3A+B.

Therefore:

    depth = min(v2(2X), v2(2Y)).

Let:

    d = gcd(A,B)
    A = d*a
    B = d*b
    gcd(a,b)=1.

Only three normalized parity states can occur:

    (a,b)=(0,1)
    (a,b)=(1,0)
    (a,b)=(1,1).

FRAME A:

    2X = d(b-a)
    2Y = d(a+b).

FRAME B:

    2X = d(3a-b)
    2Y = d(3a+b).

Because 3 is odd, both frames have the same normalized
2-adic behavior.

Opposite parity:

    (a,b)=(0,1) or (1,0)

    transformed numerators are odd,

    so:

        depth = v2(d).

Both odd:

    (a,b)=(1,1)

    transformed numerators are even.

    Since a,b are odd, one of the two transformed
    normalized numerators is exactly 2 modulo 4.

    Therefore:

        min transformed valuation = 1,

    and:

        depth = v2(d)+1.

Thus for every nondegenerate state:

    depth =
        v2(gcd(A,B))
        + [v2(A)=v2(B)].

This is an exact finite-state collapse.

The previous crash was only an implementation error:
the normalized increment-state container is now explicitly
a defaultdict(set), so multiple observed increments per
parity class can be collected safely.

The remaining unresolved quantity is:

    v2(gcd(A,B)).

That is now the only factor-side information needed.
"""
    )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print("EXPERIMENT 642 START")
    print("=" * 90)
    print()
    print(
        "NORMALIZED RESIDUAL PARITY THEOREM"
    )
    print()

    print("[1] PRIME SIEVE")

    primes = prime_sieve(
        PRIME_LIMIT
    )

    print(
        f"    odd primes={len(primes)}"
    )
    print()

    print("[2] SEMIPRIME GENERATION")

    states = build_states(
        primes
    )

    print(
        f"    semiprimes={len(states)}"
    )
    print()

    total_failures = 0

    total_failures += test_domain(
        states
    )

    total_failures += test_transformed_numerators(
        states
    )

    total_failures += test_direct_depth(
        states
    )

    total_failures += test_gcd_normalization(
        states
    )

    total_failures += test_frame_a_normalized(
        states
    )

    total_failures += test_frame_b_normalized(
        states
    )

    total_failures += test_exact_gcd_formula(
        states
    )

    total_failures += test_equal_v2_equivalence(
        states
    )

    total_failures += test_universal_depth_formula(
        states
    )

    total_failures += test_no_higher_bits_needed(
        states
    )

    total_failures += test_examples(
        states
    )

    summary()

    print()
    print("=" * 90)
    print("EXPERIMENT 642 FINISHED")
    print("=" * 90)
    print()

    print(
        f"TOTAL FAILURES={total_failures}"
    )

    if total_failures == 0:
        print(
            "STATUS=ALL TESTS PASSED"
        )
    else:
        print(
            "STATUS=COUNTEREXAMPLES FOUND"
        )


if __name__ == "__main__":
    main()