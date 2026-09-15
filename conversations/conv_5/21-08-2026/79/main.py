#!/usr/bin/env python3

"""
==========================================================================================
EXPERIMENT 641 START
==========================================================================================

EXACT FRAME-SYMMETRIC NUMERATOR PRODUCT LAW

PURPOSE
-------

Experiment 640 established exactly:

FRAME A:

    2X = B-A
    2Y = A+B

    depth = min(v2(B-A), v2(A+B))


FRAME B:

    2X = 3A-B
    2Y = 3A+B

    depth = min(v2(3A-B), v2(3A+B))

Experiment 640 Test 7 was wrong because it attempted to express
the Frame-B product using an incorrect formula.

This experiment derives the exact product/sum identities and
then tests whether the transformed numerator valuations can
be reduced to symmetric combinations of:

    p+q
    p-q
    n+c

without doing a large brute-force offset search.

The experiment is deliberately algebraic and exact.
"""


from __future__ import annotations

from dataclasses import dataclass
from collections import Counter
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
# v2
# ==============================================================================

def v2(x: int) -> int:

    if x == 0:
        return INF

    x = abs(x)

    return (x & -x).bit_length() - 1


# ==============================================================================
# gcd
# ==============================================================================

def g(a: int, b: int) -> int:
    return gcd(abs(a), abs(b))


# ==============================================================================
# FRAME
# ==============================================================================

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
# SIEVE
# ==============================================================================

def prime_sieve(limit: int) -> list[int]:

    sieve = bytearray(
        b"\x01"
    ) * (limit + 1)

    sieve[0] = 0
    sieve[1] = 0

    root = int(limit ** 0.5)

    for p in range(
        2,
        root + 1,
    ):

        if sieve[p]:

            start = p * p

            sieve[
                start:
                limit + 1:
                p
            ] = (
                b"\x00"
                * (
                    (limit - start) // p
                    + 1
                )
            )

    return [
        p
        for p in range(
            3,
            limit + 1,
            2,
        )
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

        return (
            p - 3,
            q + 3,
        )

    if frame == "B":

        return (
            p + 1,
            q - 3,
        )

    raise ValueError(frame)


# ==============================================================================
# GLOBAL XY
# ==============================================================================

def global_xy(
    frame: str,
    A: int,
    B: int,
):

    if frame == "A":

        nx = B - A
        ny = A + B

    else:

        nx = 3 * A - B
        ny = 3 * A + B

    if nx & 1 or ny & 1:

        raise ArithmeticError(
            "Non-integral global coordinates"
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
):

    return (
        1
        + v2(
            g(X, Y)
        )
    )


# ==============================================================================
# BUILD STATES
# ==============================================================================

def build_states(
    primes: list[int],
):

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

def test_direct_depth(
    states,
):

    print("=" * 90)
    print(
        "TEST 0: DIRECT TRANSFORMED-VALUATION DEPTH"
    )
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
# TEST 1
# ==============================================================================

def test_exact_frame_products(
    states,
):

    print("=" * 90)
    print(
        "TEST 1: EXACT FRAME PRODUCT IDENTITIES"
    )
    print("=" * 90)

    failures = 0
    shown = 0

    for s in states:

        if s.frame == "A":

            u = s.B - s.A
            v = s.A + s.B

            lhs = u * v

            rhs = (
                (s.q - s.p + 6)
                * (s.p + s.q)
            )

            if lhs != rhs:

                failures += 1

                if shown < 10:
                    print(
                        f"    A mismatch "
                        f"n={s.n} "
                        f"lhs={lhs} "
                        f"rhs={rhs}"
                    )
                    shown += 1

        else:

            u = 3 * s.A - s.B
            v = 3 * s.A + s.B

            lhs = u * v

            # Derive directly:
            #
            # A=p+1
            # B=q-3
            #
            # 3A-B = 3(p+1)-(q-3)
            #       = 3p-q+6
            #
            # 3A+B = 3(p+1)+(q-3)
            #       = 3p+q
            #
            # Therefore:
            #
            # product=(3p-q+6)(3p+q)

            rhs = (
                (3 * s.p - s.q + 6)
                * (3 * s.p + s.q)
            )

            if lhs != rhs:

                failures += 1

                if shown < 10:
                    print(
                        f"    B mismatch "
                        f"n={s.n} "
                        f"lhs={lhs} "
                        f"rhs={rhs}"
                    )
                    shown += 1

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 2
# ==============================================================================

def test_exact_frame_sums(
    states,
):

    print("=" * 90)
    print(
        "TEST 2: EXACT SUM/DIFferences IN p,q"
    )
    print("=" * 90)

    failures = 0

    for s in states:

        if s.frame == "A":

            u = s.B - s.A
            v = s.A + s.B

            expected_u = (
                s.q - s.p + 6
            )

            expected_v = (
                s.p + s.q
            )

        else:

            u = 3 * s.A - s.B
            v = 3 * s.A + s.B

            expected_u = (
                3 * s.p - s.q + 6
            )

            expected_v = (
                3 * s.p + s.q
            )

        if u != expected_u:
            failures += 1

        if v != expected_v:
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

def test_frame_A_n_identity(
    states,
):

    print("=" * 90)
    print(
        "TEST 3: FRAME-A PRODUCT VS n"
    )
    print("=" * 90)

    failures = 0

    shown = 0

    for s in states:

        if s.frame != "A":
            continue

        u = s.B - s.A
        v = s.A + s.B

        product = u * v

        # u*v
        # = (q-p+6)(p+q)
        #
        # = q^2-p^2+6p+6q
        #
        # Use:
        #
        # p^2+q^2 = (p+q)^2 - 2pq.
        #
        # This is tested both directly and symbolically
        # against n.

        S = s.p + s.q

        expected = (
            S * S
            - 4 * s.n
            - 36
        )

        if product != expected:

            failures += 1

            if shown < 10:

                print(
                    f"    mismatch "
                    f"n={s.n} "
                    f"product={product} "
                    f"expected={expected}"
                )

                shown += 1

    print(
        f"checked="
        f"{sum(1 for s in states if s.frame == 'A')} "
        f"failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 4
# ==============================================================================

def test_frame_B_n_identity(
    states,
):

    print("=" * 90)
    print(
        "TEST 4: FRAME-B PRODUCT VS n"
    )
    print("=" * 90)

    failures = 0

    shown = 0

    for s in states:

        if s.frame != "B":
            continue

        u = 3 * s.A - s.B
        v = 3 * s.A + s.B

        product = u * v

        # Correct Frame-B product:
        #
        # (3p-q+6)(3p+q)
        #
        # Expand:
        #
        # = 9p^2 + 18p - q^2 - 6q.
        #
        # We can eliminate q^2 with
        #
        # q^2 = (p+q)^2 - 2pq - p^2
        #
        # but retain a direct symbolic candidate.

        expected = (
            9 * s.p * s.p
            + 18 * s.p
            - s.q * s.q
            - 6 * s.q
        )

        if product != expected:

            failures += 1

            if shown < 10:

                print(
                    f"    direct mismatch "
                    f"n={s.n} "
                    f"product={product} "
                    f"expected={expected}"
                )

                shown += 1

    print(
        f"checked="
        f"{sum(1 for s in states if s.frame == 'B')} "
        f"failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 5
# ==============================================================================

def test_frame_B_symmetric_variables(
    states,
):

    print("=" * 90)
    print(
        "TEST 5: FRAME-B PRODUCT IN S=p+q, D=p-q"
    )
    print("=" * 90)

    failures = 0
    shown = 0

    for s in states:

        if s.frame != "B":
            continue

        S = s.p + s.q
        D = s.p - s.q

        u = 3 * s.A - s.B
        v = 3 * s.A + s.B

        product = u * v

        # p=(S+D)/2, q=(S-D)/2.
        #
        # Product:
        #
        # (3p-q+6)(3p+q)
        #
        # = ((S+2D+6)/2)
        #   * ((2S+D)/2)
        #
        # = (S+2D+6)(2S+D)/4.

        numerator = (
            (S + 2 * D + 6)
            * (2 * S + D)
        )

        if numerator % 4 != 0:

            failures += 1

            if shown < 10:

                print(
                    f"    non-integral "
                    f"n={s.n}"
                )

                shown += 1

            continue

        expected = numerator // 4

        if expected != product:

            failures += 1

            if shown < 10:

                print(
                    f"    mismatch "
                    f"n={s.n} "
                    f"product={product} "
                    f"expected={expected}"
                )

                shown += 1

    print(
        f"checked="
        f"{sum(1 for s in states if s.frame == 'B')} "
        f"failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 6
# ==============================================================================

def test_depth_limiting_valuation(
    states,
):

    print("=" * 90)
    print(
        "TEST 6: DEPTH LIMITING VALUATION"
    )
    print("=" * 90)

    distribution = Counter()
    failures = 0

    for s in states:

        if s.frame == "A":

            u = s.B - s.A
            v = s.A + s.B

        else:

            u = 3 * s.A - s.B
            v = 3 * s.A + s.B

        vu = v2(u)
        vv = v2(v)

        if vu < vv:

            branch = "first"

        elif vv < vu:

            branch = "second"

        else:

            branch = "equal"

        distribution[
            (
                s.frame,
                branch,
            )
        ] += 1

        if min(vu, vv) != s.depth:

            failures += 1

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )

    print()

    for frame in (
        "A",
        "B",
    ):

        print(
            f"    FRAME {frame}"
        )

        for branch in (
            "first",
            "second",
            "equal",
        ):

            print(
                f"        {branch}: "
                f"{distribution[(frame, branch)]}"
            )

        print()

    return failures


# ==============================================================================
# TEST 7
# ==============================================================================

def test_sum_difference_valuation_theorems(
    states,
):

    print("=" * 90)
    print(
        "TEST 7: 2-ADIC SUM/DIFFERENCE THEOREMS"
    )
    print("=" * 90)

    failures = 0

    # For odd coprime-normalized components:
    #
    # If v2(A) != v2(B):
    #
    #     v2(A+B)=v2(A-B)=min(v2(A),v2(B)).
    #
    # If v2(A)=v2(B):
    #
    #     both are strictly larger.
    #
    # For B frame the same phenomenon applies to
    #
    #     3A+B
    #     3A-B.
    #
    # because 3 is odd.

    checked = 0

    shown = 0

    for s in states:

        # Skip zero-degenerate cases.
        if s.A == 0 or s.B == 0:
            continue

        checked += 1

        alpha = v2(s.A)
        beta = v2(s.B)

        m = min(
            alpha,
            beta,
        )

        if s.frame == "A":

            plus = s.A + s.B
            minus = s.B - s.A

        else:

            plus = 3 * s.A + s.B
            minus = 3 * s.A - s.B

        vp = v2(plus)
        vm = v2(minus)

        if alpha != beta:

            if vp != m or vm != m:

                failures += 1

                if shown < 20:

                    print(
                        f"    mismatch "
                        f"n={s.n} "
                        f"frame={s.frame} "
                        f"alpha={alpha} "
                        f"beta={beta} "
                        f"vp={vp} "
                        f"vm={vm} "
                        f"m={m}"
                    )

                    shown += 1

        else:

            if not (
                vp > m
                and vm > m
            ):

                failures += 1

                if shown < 20:

                    print(
                        f"    mismatch equal "
                        f"n={s.n} "
                        f"frame={s.frame} "
                        f"alpha={alpha} "
                        f"beta={beta} "
                        f"vp={vp} "
                        f"vm={vm}"
                    )

                    shown += 1

    print(
        f"checked={checked} "
        f"failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 8
# ==============================================================================

def test_correct_n_c_relation(
    states,
):

    print("=" * 90)
    print(
        "TEST 8: CORRECT n+c RELATIONS"
    )
    print("=" * 90)

    """
    We do NOT assume that the product of transformed
    numerators is simply U^2 - 4(n+c) in Frame B.

    Instead derive the exact relations.

    FRAME A:

        U = p+q

        (B-A)(A+B)
          = (q-p+6)(p+q)
          = U^2 - 4(n+9).

    FRAME B:

        Let

            S = p+q
            D = p-q.

        Then:

            3A-B = 3p-q+6
            3A+B = 3p+q

        Their product is:

            (3p-q+6)(3p+q)

        Expanding in S,D gives:

            (S+2D+6)(2S+D)/4.

    We verify these exact forms.
    """

    failures = 0
    counts = Counter()

    for s in states:

        if s.frame == "A":

            U = s.p + s.q

            actual = (
                (3 * 0 + s.q - s.p + 6)
                * U
            )

            expected = (
                U * U
                - 4 * (s.n + 9)
            )

        else:

            S = s.p + s.q
            D = s.p - s.q

            actual = (
                (3 * s.p - s.q + 6)
                * (3 * s.p + s.q)
            )

            numerator = (
                (S + 2 * D + 6)
                * (2 * S + D)
            )

            if numerator % 4 != 0:

                failures += 1
                continue

            expected = numerator // 4

        if actual != expected:

            failures += 1

        counts[
            s.frame
        ] += 1

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )

    print(
        f"    frame counts={dict(counts)}"
    )

    print()

    return failures


# ==============================================================================
# TEST 9
# ==============================================================================

def test_small_exact_symmetric_signatures(
    states,
):

    print("=" * 90)
    print(
        "TEST 9: SMALL SYMMETRIC DEPTH SIGNATURES"
    )
    print("=" * 90)

    candidates = [
        (
            "A: min(v2(p+q),v2(p-q+6))",
            lambda s: min(
                v2(s.p + s.q),
                v2(s.p - s.q + 6),
            ),
        ),
        (
            "A: min(v2(p+q),v2(q-p+6))",
            lambda s: min(
                v2(s.p + s.q),
                v2(s.q - s.p + 6),
            ),
        ),
        (
            "B: min(v2(3p-q+6),v2(3p+q))",
            lambda s: min(
                v2(3 * s.p - s.q + 6),
                v2(3 * s.p + s.q),
            ),
        ),
    ]

    for name, fn in candidates:

        failures = 0

        checked = 0

        for s in states:

            if name.startswith("A:") and s.frame != "A":
                continue

            if name.startswith("B:") and s.frame != "B":
                continue

            checked += 1

            predicted = fn(s)

            if predicted != s.depth:

                failures += 1

        print(
            f"    {name}"
        )

        print(
            f"        checked={checked} "
            f"failures={failures}"
        )

    print()


# ==============================================================================
# TEST 10
# ==============================================================================

def test_known_examples(
    states,
):

    print("=" * 90)
    print(
        "TEST 10: KNOWN COUNTEREXAMPLE REGRESSION"
    )
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

        s = lookup[n]

        if s.frame == "A":

            first = s.B - s.A
            second = s.A + s.B

        else:

            first = 3 * s.A - s.B
            second = 3 * s.A + s.B

        predicted = min(
            v2(first),
            v2(second),
        )

        print(
            f"    n={n:<7}"
            f"frame={s.frame} "
            f"A={s.A:<5} "
            f"B={s.B:<5} "
            f"v2u={v2(first):<3} "
            f"v2v={v2(second):<3} "
            f"depth={s.depth}"
        )

        if predicted != s.depth:

            failures += 1

    print()

    print(
        f"checked={len(examples)} "
        f"failures={failures}"
    )

    print()

    return failures


# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

def final_summary():

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)

    print(
r"""
Experiment 640 established the exact coordinate-level theorem:

    FRAME A:

        2X = B-A
        2Y = A+B

        depth =
            min(
                v2(B-A),
                v2(A+B)
            ).

    FRAME B:

        2X = 3A-B
        2Y = 3A+B

        depth =
            min(
                v2(3A-B),
                v2(3A+B)
            ).

The failure in Experiment 640 Test 7 was only a bad
Frame-B product identity.

The correct Frame-B numerator product is:

    (3A-B)(3A+B)
      =
    9A^2-B^2

and, with:

    A=p+1
    B=q-3,

it is:

    (3p-q+6)(3p+q).

The important structural observation is now:

    DEPTH IS THE MINIMUM OF TWO 2-ADIC NUMERATORS.

For both frames the two numerators have the same
2-adic behavior determined by the valuations of A,B.

If:

    v2(A) != v2(B),

then for either frame:

    v2(first)
      =
    v2(second)
      =
    min(v2(A),v2(B)).

If:

    v2(A) = v2(B),

then both transformed numerators gain
an additional factor of 2.

Therefore the entire depth is again controlled by
the common normalized residual parity, but now in an
exact symmetric form:

    depth =
        min(v2(transformed_1),
            v2(transformed_2)).

The remaining research question is:

    Can this minimum be expressed through a smaller
    symmetric invariant without recovering p and q?

The natural candidates are now:

    FRAME A:
        p+q
        q-p+6

    FRAME B:
        3p+q
        3p-q+6

rather than the incorrect single quantity v2(A+B).
"""
    )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print("EXPERIMENT 641 START")
    print("=" * 90)
    print()
    print(
        "EXACT FRAME-SYMMETRIC NUMERATOR PRODUCT LAW"
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

    failures = 0

    failures += test_direct_depth(
        states
    )

    failures += test_exact_frame_products(
        states
    )

    failures += test_exact_frame_sums(
        states
    )

    failures += test_frame_A_n_identity(
        states
    )

    failures += test_frame_B_n_identity(
        states
    )

    failures += test_frame_B_symmetric_variables(
        states
    )

    failures += test_depth_limiting_valuation(
        states
    )

    failures += test_sum_difference_valuation_theorems(
        states
    )

    failures += test_correct_n_c_relation(
        states
    )

    test_small_exact_symmetric_signatures(
        states
    )

    failures += test_known_examples(
        states
    )

    final_summary()

    print()
    print("=" * 90)
    print("EXPERIMENT 641 FINISHED")
    print("=" * 90)
    print()

    print(
        f"TOTAL FAILURES={failures}"
    )

    if failures == 0:
        print(
            "STATUS=ALL CORE TESTS PASSED"
        )
    else:
        print(
            "STATUS=COUNTEREXAMPLES FOUND"
        )


if __name__ == "__main__":
    main()
