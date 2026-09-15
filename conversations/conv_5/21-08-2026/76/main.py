#!/usr/bin/env python3

"""
==========================================================================================
EXPERIMENT 638
==========================================================================================

SYMMETRIC gcd -> MINIMUM OF TWO 2-ADIC VALUATIONS

Established candidate:

    depth = v2(gcd(U,R))

with

    FRAME A:
        U = p + q
        R = (p-q)^2 - 36

    FRAME B:
        U = p + q - 2
        R = (p-q)^2 - 16

This experiment proves/test:

    gcd(U,R)
        =
    gcd(U, 4(n+c))

where

    c=9 for FRAME A
    c=3 for FRAME B.

Therefore:

    depth
        =
    min(v2(U), v2(n+c)+2).

Then investigate:

    1. Does this reproduce every observed depth?
    2. Is U exactly A+B in both frames?
    3. Can v2(U) be expressed through A,B,d?
    4. What is the relation between v2(U), v2(A), v2(B)?
    5. Can the whole depth be written without the discriminant?
    6. Can the remaining U valuation be reduced further?
"""


from __future__ import annotations

from dataclasses import dataclass
from math import gcd
from collections import Counter, defaultdict


# ==============================================================================
# CONSTANTS
# ==============================================================================

PRIME_LIMIT = 6000
INF = 10**9


# ==============================================================================
# DATA
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


def gcd0(a: int, b: int) -> int:
    return gcd(abs(a), abs(b))


def frame_from_n(n: int) -> str:
    r = n & 3

    if r == 1:
        return "A"

    if r == 3:
        return "B"

    raise ValueError(
        f"Expected odd n, got n mod 4={r}"
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

        if sieve[p]:

            start = p * p

            sieve[start:limit + 1:p] = (
                b"\x00"
                * (((limit - start) // p) + 1)
            )

    return [
        p
        for p in range(3, limit + 1, 2)
        if sieve[p]
    ]


# ==============================================================================
# RESIDUALS
# ==============================================================================

def residuals(
    frame: str,
    p: int,
    q: int,
):
    if frame == "A":
        return p - 3, q + 3

    if frame == "B":
        return p + 1, q - 3

    raise ValueError(frame)


# ==============================================================================
# GLOBAL X,Y
# ==============================================================================

def global_xy(
    frame: str,
    A: int,
    B: int,
):
    if frame == "A":

        nx = B - A
        ny = A + B

    elif frame == "B":

        nx = 3 * A - B
        ny = 3 * A + B

    else:
        raise ValueError(frame)

    if nx & 1:
        raise ArithmeticError(
            f"Non-integral X: frame={frame}, A={A}, B={B}"
        )

    if ny & 1:
        raise ArithmeticError(
            f"Non-integral Y: frame={frame}, A={A}, B={B}"
        )

    return nx // 2, ny // 2


# ==============================================================================
# GLOBAL DEPTH
# ==============================================================================

def global_depth(X: int, Y: int) -> int:

    g = gcd0(X, Y)

    if g == 0:
        raise ArithmeticError(
            "gcd(X,Y)=0"
        )

    return 1 + v2(g)


# ==============================================================================
# STATES
# ==============================================================================

def build_states(
    primes: list[int],
) -> list[State]:

    states = []

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

            depth = global_depth(
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
# SYMMETRIC VARIABLES
# ==============================================================================

def symmetric_values(
    s: State,
):
    S = s.p + s.q
    D = s.p - s.q

    if s.frame == "A":

        U = S
        c = 9
        R = D * D - 36

    else:

        U = S - 2
        c = 3
        R = D * D - 16

    return S, D, U, R, c


# ==============================================================================
# TEST 0
# ==============================================================================

def test_U_equals_A_plus_B(
    states,
):
    print("=" * 90)
    print("TEST 0: U == A+B")
    print("=" * 90)

    failures = 0

    for s in states:

        _, _, U, _, _ = symmetric_values(s)

        if U != s.A + s.B:
            failures += 1

    print(
        f"checked={len(states)} failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 1
# ==============================================================================

def test_exact_congruence(
    states,
):
    print("=" * 90)
    print("TEST 1: EXACT R CONGRUENCE")
    print("=" * 90)

    failures = 0

    for s in states:

        _, _, U, R, c = symmetric_values(s)

        # Candidate:
        #
        #     R + 4(n+c)
        #
        # should be divisible by U.

        delta = R + 4 * (s.n + c)

        if delta % U != 0:
            failures += 1

    print(
        f"checked={len(states)} failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 2
# ==============================================================================

def test_exact_gcd_reduction(
    states,
):
    print("=" * 90)
    print("TEST 2: gcd(U,R) == gcd(U,4(n+c))")
    print("=" * 90)

    failures = 0
    shown = 0

    for s in states:

        _, _, U, R, c = symmetric_values(s)

        g1 = gcd0(
            U,
            R,
        )

        g2 = gcd0(
            U,
            4 * (s.n + c),
        )

        if g1 != g2:

            failures += 1

            if shown < 20:

                print(
                    "    mismatch:",
                    f"n={s.n}",
                    f"frame={s.frame}",
                    f"U={U}",
                    f"R={R}",
                    f"n+c={s.n+c}",
                    f"g1={g1}",
                    f"g2={g2}",
                )

                shown += 1

    print(
        f"checked={len(states)} failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 3
# ==============================================================================

def test_depth_min_formula(
    states,
):
    print("=" * 90)
    print("TEST 3: DEPTH = min(v2(U), v2(n+c)+2)")
    print("=" * 90)

    failures = 0
    shown = 0

    distribution = Counter()

    for s in states:

        _, _, U, _, c = symmetric_values(s)

        lhs = min(
            v2(U),
            v2(s.n + c) + 2,
        )

        distribution[lhs] += 1

        if lhs != s.depth:

            failures += 1

            if shown < 20:

                print(
                    "    mismatch:",
                    f"n={s.n}",
                    f"frame={s.frame}",
                    f"U={U}",
                    f"v2U={v2(U)}",
                    f"n+c={s.n+c}",
                    f"v2(n+c)={v2(s.n+c)}",
                    f"pred={lhs}",
                    f"actual={s.depth}",
                )

                shown += 1

    print(
        f"checked={len(states)} failures={failures}"
    )

    print()
    print("    predicted depth distribution:")

    for d in sorted(distribution):
        print(
            f"        depth={d}: {distribution[d]}"
        )

    print()

    return failures


# ==============================================================================
# TEST 4
# ==============================================================================

def test_depth_against_symmetric_gcd(
    states,
):
    print("=" * 90)
    print("TEST 4: SYMMETRIC gcd DEPTH CROSSCHECK")
    print("=" * 90)

    failures = 0

    for s in states:

        _, _, U, R, _ = symmetric_values(s)

        g = gcd0(
            U,
            R,
        )

        if v2(g) != s.depth:
            failures += 1

    print(
        f"checked={len(states)} failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 5
# ==============================================================================

def test_U_valuation_vs_raw_residuals(
    states,
):
    print("=" * 90)
    print("TEST 5: v2(U) VS v2(A),v2(B)")
    print("=" * 90)

    failures = 0

    relation_counts = Counter()

    examples = {}

    for s in states:

        if s.A == 0 or s.B == 0:
            continue

        U = s.A + s.B

        alpha = v2(s.A)
        beta = v2(s.B)
        gamma = v2(U)

        m = min(
            alpha,
            beta,
        )

        if alpha != beta:

            # Standard valuation identity:
            #
            # v2(A+B)=min(v2(A),v2(B))

            expected = m
            relation = "unequal"

        else:

            # Equal valuation means the sum gains
            # at least one extra factor of two.

            expected = gamma
            relation = "equal"

            if gamma <= m:
                failures += 1

        if alpha != beta and gamma != expected:
            failures += 1

        relation_counts[
            (
                s.frame,
                relation,
                gamma - m,
            )
        ] += 1

        examples.setdefault(
            (
                s.frame,
                relation,
                gamma - m,
            ),
            s,
        )

    print(
        f"checked={sum(relation_counts.values())} failures={failures}"
    )

    print()
    print("    valuation relations:")

    for key in sorted(
        relation_counts,
        key=str,
    ):

        print(
            f"        {key}: {relation_counts[key]}"
        )

    print()

    return failures


# ==============================================================================
# TEST 6
# ==============================================================================

def test_depth_as_min_raw_sum(
    states,
):
    print("=" * 90)
    print("TEST 6: DEPTH = min(v2(A+B), v2(n+c)+2)")
    print("=" * 90)

    failures = 0
    shown = 0

    for s in states:

        U = s.A + s.B

        if s.frame == "A":
            c = 9
        else:
            c = 3

        predicted = min(
            v2(U),
            v2(s.n + c) + 2,
        )

        if predicted != s.depth:

            failures += 1

            if shown < 20:

                print(
                    "    mismatch:",
                    f"n={s.n}",
                    f"frame={s.frame}",
                    f"A={s.A}",
                    f"B={s.B}",
                    f"A+B={U}",
                    f"pred={predicted}",
                    f"depth={s.depth}",
                )

                shown += 1

    print(
        f"checked={len(states)} failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 7
# ==============================================================================

def test_primary_cannot_exceed_U(
    states,
):
    print("=" * 90)
    print("TEST 7: PRIMARY-vs-SUM VALUATION")
    print("=" * 90)

    failures = 0

    counts = Counter()

    for s in states:

        U = s.A + s.B

        if s.frame == "A":
            c = 9
        else:
            c = 3

        vu = v2(U)
        vn = v2(s.n + c) + 2

        depth = s.depth

        branch = (
            "U-limited"
            if vu < vn
            else "n-limited"
            if vn < vu
            else "tie"
        )

        counts[
            (s.frame, branch)
        ] += 1

        if depth != min(vu, vn):
            failures += 1

    print(
        f"checked={len(states)} failures={failures}"
    )

    print()
    print("    limiting branches:")

    for key in sorted(
        counts,
        key=str,
    ):
        print(
            f"        {key}: {counts[key]}"
        )

    print()

    return failures


# ==============================================================================
# TEST 8
# ==============================================================================

def test_frame_specific_closed_forms(
    states,
):
    print("=" * 90)
    print("TEST 8: FRAME-SPECIFIC CLOSED FORMS")
    print("=" * 90)

    failures = 0

    for s in states:

        if s.frame == "A":

            # U = p+q
            # R = (p-q)^2-36
            #
            # gcd reduction:
            #
            # gcd(
            #     p+q,
            #     (p-q)^2-36
            # )
            #
            # =
            #
            # gcd(
            #     p+q,
            #     4(n+9)
            # ).

            U = s.p + s.q
            target = 4 * (s.n + 9)

        else:

            # U = p+q-2
            #
            # gcd(
            #     p+q-2,
            #     (p-q)^2-16
            # )
            #
            # =
            #
            # gcd(
            #     p+q-2,
            #     4(n+3)
            # ).

            U = s.p + s.q - 2
            target = 4 * (s.n + 3)

        info = symmetric_values(s)

        R = info[3]

        if gcd0(U, R) != gcd0(U, target):
            failures += 1

    print(
        f"checked={len(states)} failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 9
# ==============================================================================

def test_discriminant_elimination(
    states,
):
    print("=" * 90)
    print("TEST 9: DISCRIMINANT ELIMINATION")
    print("=" * 90)

    failures = 0

    for s in states:

        S = s.p + s.q

        if s.frame == "A":

            U = S
            c = 9

        else:

            U = S - 2
            c = 3

        # No p-q is used here.

        R_reconstructed = (
            S * S
            - 4 * s.n
            - (
                36
                if s.frame == "A"
                else 16
            )
        )

        # Same gcd.
        g = gcd0(
            U,
            R_reconstructed,
        )

        predicted = min(
            v2(U),
            v2(s.n + c) + 2,
        )

        if v2(g) != predicted:
            failures += 1

    print(
        f"checked={len(states)} failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 10
# ==============================================================================

def test_only_remaining_information(
    states,
):
    print("=" * 90)
    print("TEST 10: REMAINING INFORMATION = v2(A+B)")
    print("=" * 90)

    """
    At this point:

        depth =
            min(
                v2(A+B),
                v2(n+c)+2
            ).

    This test checks that replacing U by A+B loses nothing.

    """

    failures = 0

    for s in states:

        if s.frame == "A":
            c = 9
        else:
            c = 3

        predicted = min(
            v2(s.A + s.B),
            v2(s.n + c) + 2,
        )

        if predicted != s.depth:
            failures += 1

    print(
        f"checked={len(states)} failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# EXAMPLES
# ==============================================================================

def print_examples(
    states,
):
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
        111,
        141,
        183,
        213,
    ]

    lookup = {
        s.n: s
        for s in states
    }

    for n in wanted:

        s = lookup.get(n)

        if s is None:
            continue

        S = s.p + s.q

        D = s.p - s.q

        if s.frame == "A":

            U = S
            c = 9
            R = D * D - 36

        else:

            U = S - 2
            c = 3
            R = D * D - 16

        g = gcd0(
            U,
            R,
        )

        print()
        print(
            f"n={n} p={s.p} q={s.q} frame={s.frame}"
        )

        print(
            f"    A={s.A} B={s.B}"
        )

        print(
            f"    A+B={s.A+s.B}"
        )

        print(
            f"    S=p+q={S}"
        )

        print(
            f"    D=p-q={D}"
        )

        print(
            f"    U={U}"
        )

        print(
            f"    R={R}"
        )

        print(
            f"    gcd(U,R)={g}"
        )

        print(
            f"    v2(U)={v2(U)}"
        )

        print(
            f"    n+c={s.n+c}"
        )

        print(
            f"    v2(n+c)+2={v2(s.n+c)+2}"
        )

        print(
            f"    min={min(v2(U),v2(s.n+c)+2)}"
        )

        print(
            f"    depth={s.depth}"
        )


# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

def print_summary():

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)

    print(
r"""
The next reduction is:

    FRAME A:

        U = p+q

        R = (p-q)^2-36

        R
          =
        (p+q)^2 - 4n - 36

        R
          =
        U^2 - 4(n+9).

    Therefore:

        R ≡ -4(n+9) (mod U)

    and hence:

        gcd(U,R)
          =
        gcd(U,4(n+9)).


    FRAME B:

        U = p+q-2

        R = (p-q)^2-16

        R
          =
        (p+q)^2 - 4n - 16

        since p+q=U+2,

        R
          =
        U^2 + 4U - 4(n+3).

    Therefore:

        R ≡ -4(n+3) (mod U)

    and:

        gcd(U,R)
          =
        gcd(U,4(n+3)).


Thus the depth candidate collapses to:

    FRAME A:

        depth =
        min(
            v2(p+q),
            v2(n+9)+2
        ).

    FRAME B:

        depth =
        min(
            v2(p+q-2),
            v2(n+3)+2
        ).

Because:

    U = A+B,

the remaining factor-side information is simply:

    v2(A+B).

This is substantially smaller than the original state:

    (p,q)
        ->
    (A,B)
        ->
    (X,Y)
        ->
    level traversal.

The reduced chain becomes:

    (p,q)
        ->
    U=A+B

and

    n+c

with

    depth=min(
        v2(U),
        v2(n+c)+2
    ).

The next question after this experiment is therefore very
specific:

    Can v2(A+B) be characterized directly from n and frame,
    or is v2(A+B) genuinely additional information?

Unlike the previous brute-force n-only signature searches,
this is now testing one exact structural quantity.
"""
    )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print("EXPERIMENT 638 START")
    print("=" * 90)

    print()
    print(
        "SYMMETRIC gcd -> MINIMUM OF TWO 2-ADIC VALUATIONS"
    )
    print()

    # --------------------------------------------------------------------------
    # PRIME SIEVE
    # --------------------------------------------------------------------------

    print("[1] PRIME SIEVE")

    primes = prime_sieve(
        PRIME_LIMIT
    )

    print(
        f"    odd primes={len(primes)}"
    )
    print()

    # --------------------------------------------------------------------------
    # SEMIPRIMES
    # --------------------------------------------------------------------------

    print("[2] SEMIPRIME GENERATION")

    states = build_states(
        primes
    )

    print(
        f"    semiprimes={len(states)}"
    )
    print()

    # --------------------------------------------------------------------------
    # TESTS
    # --------------------------------------------------------------------------

    failures = 0

    failures += test_U_equals_A_plus_B(
        states
    )

    failures += test_exact_congruence(
        states
    )

    failures += test_exact_gcd_reduction(
        states
    )

    failures += test_depth_min_formula(
        states
    )

    failures += test_depth_against_symmetric_gcd(
        states
    )

    failures += test_U_valuation_vs_raw_residuals(
        states
    )

    failures += test_depth_as_min_raw_sum(
        states
    )

    failures += test_primary_cannot_exceed_U(
        states
    )

    failures += test_frame_specific_closed_forms(
        states
    )

    failures += test_discriminant_elimination(
        states
    )

    failures += test_only_remaining_information(
        states
    )

    # --------------------------------------------------------------------------
    # EXAMPLES
    # --------------------------------------------------------------------------

    print_examples(
        states
    )

    # --------------------------------------------------------------------------
    # SUMMARY
    # --------------------------------------------------------------------------

    print_summary()

    print()
    print("=" * 90)
    print("EXPERIMENT 638 FINISHED")
    print("=" * 90)
    print()

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
