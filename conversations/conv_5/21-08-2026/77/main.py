#!/usr/bin/env python3

"""
==========================================================================================
EXPERIMENT 639
==========================================================================================

CORRECTED SYMMETRIC SUM VALUATION REDUCTION

IMPORTANT FRAME CONVENTION
---------------------------

Established by the previous experiments:

    FRAME A:
        n == 3 mod 4

    FRAME B:
        n == 1 mod 4

Examples:

    n=15 -> A
    n=21 -> B
    n=9  -> B

The previous Experiment 638 accidentally reversed this mapping.
Its internal tests therefore passed against a self-consistent,
but incorrect, frame assignment.

This experiment fixes that first.

CORE VARIABLES
--------------

FRAME A:

    A = p - 3
    B = q + 3

    X = (B-A)/2
    Y = (A+B)/2

    U = A+B
      = p+q

FRAME B:

    A = p + 1
    B = q - 3

    X = (3A-B)/2
    Y = (3A+B)/2

    U = A+B
      = p+q-2

GLOBAL DEPTH:

    depth = 1 + v2(gcd(X,Y))

SYMMETRIC gcd:

FRAME A:

    R = (p-q)^2 - 36

FRAME B:

    R = (p-q)^2 - 16

and:

    gcd(U,R)
      =
    gcd(U,4(n+c))

where:

    c=9  for FRAME A
    c=3  for FRAME B.

Therefore:

    depth
      =
    min(
        v2(U),
        v2(n+c)+2
    ).

The remaining structural quantity is:

    v2(U)=v2(A+B).

This experiment determines how that valuation relates to:

    v2(A)
    v2(B)
    v2(A-B)
    gcd(A,B)
    p+q
    p-q
    n
    frame.

No brute-force offset search is performed.
"""


from __future__ import annotations

from dataclasses import dataclass
from math import gcd
from collections import Counter, defaultdict


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
# HELPERS
# ==============================================================================

def v2(x: int) -> int:
    if x == 0:
        return INF

    x = abs(x)
    return (x & -x).bit_length() - 1


def gcd0(a: int, b: int) -> int:
    return gcd(abs(a), abs(b))


# ==============================================================================
# CORRECT FRAME DECODER
# ==============================================================================

def frame_from_n(n: int) -> str:
    """
    Established convention:

        n % 4 == 3 -> FRAME A
        n % 4 == 1 -> FRAME B
    """

    r = n & 3

    if r == 3:
        return "A"

    if r == 1:
        return "B"

    raise ValueError(
        f"Expected odd semiprime, got n mod 4={r}"
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
# RESIDUAL MAP
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
            f"Non-integral X: frame={frame} A={A} B={B}"
        )

    if ny & 1:
        raise ArithmeticError(
            f"Non-integral Y: frame={frame} A={A} B={B}"
        )

    return nx // 2, ny // 2


# ==============================================================================
# GLOBAL DEPTH
# ==============================================================================

def global_depth(
    X: int,
    Y: int,
) -> int:

    g = gcd0(
        X,
        Y,
    )

    if g == 0:
        raise ArithmeticError(
            "gcd(X,Y)=0"
        )

    return 1 + v2(g)


# ==============================================================================
# BUILD STATES
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
# TEST 0
# ==============================================================================

def test_frame_regression():
    print("=" * 90)
    print("TEST 0: FRAME REGRESSION CHECK")
    print("=" * 90)

    checks = {
        9: "B",
        15: "A",
        21: "B",
        33: "B",
        39: "A",
        57: "B",
        69: "B",
        77: "B",
        93: "B",
        111: "A",
        141: "B",
        183: "A",
        213: "B",
    }

    failures = 0

    for n, expected in checks.items():

        actual = frame_from_n(n)

        if actual != expected:

            failures += 1

            print(
                f"    mismatch n={n} "
                f"expected={expected} "
                f"actual={actual}"
            )

    print(
        f"checked={len(checks)} failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 1
# ==============================================================================

def test_state_consistency(
    states,
):
    print("=" * 90)
    print("TEST 1: GLOBAL STATE CONSISTENCY")
    print("=" * 90)

    failures = 0

    for s in states:

        if frame_from_n(s.n) != s.frame:
            failures += 1
            continue

        A, B = residuals(
            s.frame,
            s.p,
            s.q,
        )

        if (A, B) != (s.A, s.B):
            failures += 1
            continue

        X, Y = global_xy(
            s.frame,
            s.A,
            s.B,
        )

        if (X, Y) != (s.X, s.Y):
            failures += 1
            continue

        if global_depth(
            s.X,
            s.Y,
        ) != s.depth:
            failures += 1

    print(
        f"checked={len(states)} failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 2
# ==============================================================================

def test_U_equals_A_plus_B(
    states,
):
    print("=" * 90)
    print("TEST 2: U == A+B")
    print("=" * 90)

    failures = 0

    for s in states:

        if s.frame == "A":

            U = s.p + s.q

        else:

            U = s.p + s.q - 2

        if U != s.A + s.B:
            failures += 1

    print(
        f"checked={len(states)} failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 3
# ==============================================================================

def test_frame_specific_U(
    states,
):
    print("=" * 90)
    print("TEST 3: FRAME-SPECIFIC U")
    print("=" * 90)

    failures = 0

    for s in states:

        if s.frame == "A":

            expected = s.p + s.q

        else:

            expected = s.p + s.q - 2

        actual = s.A + s.B

        if actual != expected:
            failures += 1

    print(
        f"checked={len(states)} failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 4
# ==============================================================================

def test_symmetric_congruence(
    states,
):
    print("=" * 90)
    print("TEST 4: R + 4(n+c) IS DIVISIBLE BY U")
    print("=" * 90)

    failures = 0
    shown = 0

    for s in states:

        D = s.p - s.q

        if s.frame == "A":

            U = s.p + s.q
            c = 9
            R = D * D - 36

        else:

            U = s.p + s.q - 2
            c = 3
            R = D * D - 16

        delta = R + 4 * (s.n + c)

        if delta % U != 0:

            failures += 1

            if shown < 20:

                print(
                    f"    mismatch n={s.n} "
                    f"frame={s.frame} "
                    f"U={U} "
                    f"R={R} "
                    f"n+c={s.n+c}"
                )

                shown += 1

    print(
        f"checked={len(states)} failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 5
# ==============================================================================

def test_gcd_reduction(
    states,
):
    print("=" * 90)
    print("TEST 5: gcd(U,R) == gcd(U,4(n+c))")
    print("=" * 90)

    failures = 0

    for s in states:

        D = s.p - s.q

        if s.frame == "A":

            U = s.p + s.q
            c = 9
            R = D * D - 36

        else:

            U = s.p + s.q - 2
            c = 3
            R = D * D - 16

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

    print(
        f"checked={len(states)} failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 6
# ==============================================================================

def test_depth_min(
    states,
):
    print("=" * 90)
    print("TEST 6: DEPTH = min(v2(U), v2(n+c)+2)")
    print("=" * 90)

    failures = 0
    distribution = Counter()

    for s in states:

        U = s.A + s.B

        c = 9 if s.frame == "A" else 3

        predicted = min(
            v2(U),
            v2(s.n + c) + 2,
        )

        distribution[predicted] += 1

        if predicted != s.depth:
            failures += 1

    print(
        f"checked={len(states)} failures={failures}"
    )

    print()
    print("    predicted depth distribution:")

    for depth in sorted(distribution):

        print(
            f"        depth={depth}: "
            f"{distribution[depth]}"
        )

    print()

    return failures


# ==============================================================================
# TEST 7
# ==============================================================================

def test_v2U_residual_structure(
    states,
):
    print("=" * 90)
    print("TEST 7: v2(A+B) RESIDUAL STRUCTURE")
    print("=" * 90)

    failures = 0
    relation_counts = Counter()

    for s in states:

        if s.A == 0 or s.B == 0:
            continue

        alpha = v2(s.A)
        beta = v2(s.B)
        sigma = v2(s.A + s.B)

        m = min(
            alpha,
            beta,
        )

        if alpha != beta:

            if sigma != m:
                failures += 1

            relation = "unequal"

        else:

            if sigma <= m:
                failures += 1

            relation = "equal"

        relation_counts[
            (
                s.frame,
                relation,
                sigma - m,
            )
        ] += 1

    print(
        f"checked={sum(relation_counts.values())} "
        f"failures={failures}"
    )

    print()
    print("    relation distribution:")

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
# TEST 8
# ==============================================================================

def test_U_vs_pq(
    states,
):
    print("=" * 90)
    print("TEST 8: U VALUATION VS p+q")
    print("=" * 90)

    failures = 0

    counts = Counter()

    for s in states:

        if s.frame == "A":

            U = s.p + s.q

        else:

            U = s.p + s.q - 2

        counts[
            (
                s.frame,
                v2(U),
            )
        ] += 1

        # direct structural identity
        if U != s.A + s.B:
            failures += 1

    print(
        f"checked={len(states)} failures={failures}"
    )

    print()
    print("    v2(U) distribution by frame:")

    for frame in ("A", "B"):

        print()
        print(f"    FRAME {frame}")

        rows = [
            (key[1], count)
            for key, count in counts.items()
            if key[0] == frame
        ]

        for valuation, count in sorted(rows):

            print(
                f"        v2(U)={valuation}: {count}"
            )

    print()

    return failures


# ==============================================================================
# TEST 9
# ==============================================================================

def test_possible_n_only_reduction(
    states,
):
    print("=" * 90)
    print("TEST 9: n-ONLY COLLISION ANALYSIS")
    print("=" * 90)

    """
    This deliberately does NOT search hundreds of offsets.

    Instead it measures the strongest direct candidate:

        (frame, v2(n+c))

    -> v2(U)

    If this is not exact, we print the first collision.

    """

    buckets = defaultdict(set)

    examples = {}

    for s in states:

        c = 9 if s.frame == "A" else 3

        signature = (
            s.frame,
            v2(s.n + c),
        )

        U = s.A + s.B

        target = v2(U)

        buckets[signature].add(target)

        examples.setdefault(
            (signature, target),
            s,
        )

    ambiguous = 0
    shown = 0

    for signature in sorted(
        buckets,
        key=str,
    ):

        values = buckets[signature]

        if len(values) > 1:

            ambiguous += 1

            if shown < 20:

                print(
                    f"    collision signature={signature} "
                    f"v2U={sorted(values)}"
                )

                for target in sorted(values)[:4]:

                    s = examples[
                        (signature, target)
                    ]

                    print(
                        f"        "
                        f"n={s.n} "
                        f"p={s.p} "
                        f"q={s.q} "
                        f"U={s.A+s.B} "
                        f"depth={s.depth}"
                    )

                shown += 1

    print()
    print(
        f"signature buckets={len(buckets)}"
    )

    print(
        f"ambiguous buckets={ambiguous}"
    )

    print()

    return 0


# ==============================================================================
# TEST 10
# ==============================================================================

def test_exact_final_formula(
    states,
):
    print("=" * 90)
    print("TEST 10: FINAL MINIMUM FORMULA")
    print("=" * 90)

    failures = 0

    for s in states:

        U = s.A + s.B

        c = 9 if s.frame == "A" else 3

        predicted = min(
            v2(U),
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
        51,
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

    lookup = {
        s.n: s
        for s in states
    }

    for n in wanted:

        s = lookup.get(n)

        if s is None:
            continue

        U = s.A + s.B

        c = 9 if s.frame == "A" else 3

        print()
        print(
            f"n={n} p={s.p} q={s.q} "
            f"frame={s.frame}"
        )

        print(
            f"    A={s.A} B={s.B}"
        )

        print(
            f"    U=A+B={U}"
        )

        if s.frame == "A":

            expected_U = s.p + s.q

        else:

            expected_U = s.p + s.q - 2

        print(
            f"    U formula={expected_U}"
        )

        print(
            f"    v2(A)={v2(s.A)}"
        )

        print(
            f"    v2(B)={v2(s.B)}"
        )

        print(
            f"    v2(A+B)={v2(U)}"
        )

        print(
            f"    n+c={s.n+c}"
        )

        print(
            f"    v2(n+c)+2={v2(s.n+c)+2}"
        )

        print(
            f"    predicted depth="
            f"{min(v2(U), v2(s.n+c)+2)}"
        )

        print(
            f"    actual depth={s.depth}"
        )


# ==============================================================================
# SUMMARY
# ==============================================================================

def print_summary():

    print("=" * 90)
    print("FINAL STRUCTURAL SUMMARY")
    print("=" * 90)

    print(
r"""
The frame convention is:

    FRAME A:
        n ≡ 3 (mod 4)

    FRAME B:
        n ≡ 1 (mod 4)

and the residual coordinates are:

    FRAME A:

        A = p-3
        B = q+3

        U=A+B=p+q

    FRAME B:

        A = p+1
        B = q-3

        U=A+B=p+q-2.

The corrected symmetric identities are:

    FRAME A:

        R=(p-q)^2-36
         = U^2 - 4(n+9)

    FRAME B:

        R=(p-q)^2-16
         = U^2 + 4U - 4(n+3).

Hence:

    gcd(U,R)
      =
    gcd(U,4(n+c))

with:

    c=9  in FRAME A
    c=3  in FRAME B.

Because:

    depth=v2(gcd(U,R)),

we obtain:

    depth
      =
    min(
        v2(U),
        v2(n+c)+2
    ).

Equivalently:

    FRAME A:

        depth
          =
        min(
            v2(p+q),
            v2(n+9)+2
        )

    FRAME B:

        depth
          =
        min(
            v2(p+q-2),
            v2(n+3)+2
        ).

The remaining structural quantity is therefore exactly:

    v2(A+B).

The residual valuation law is:

    if v2(A) != v2(B):

        v2(A+B)
          =
        min(v2(A),v2(B))

    if v2(A) == v2(B):

        v2(A+B)
          >
        v2(A).

So the sum valuation is precisely the quantity carrying
the additional common 2-adic structure of the two residuals.

The next logical question is no longer the generic
'n-only offset search'.

It is:

    Can v2(A+B) be reconstructed from the known
    2-adic structure of n and the frame?

The first direct candidate is:

    (frame, v2(n+c))
        -> v2(A+B).

If this is not exact, the collision examples printed by
TEST 9 identify exactly where the missing information begins.
"""
    )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print("EXPERIMENT 639 START")
    print("=" * 90)

    print()
    print(
        "CORRECTED SYMMETRIC SUM VALUATION REDUCTION"
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

    failures += test_frame_regression()

    failures += test_state_consistency(
        states
    )

    failures += test_U_equals_A_plus_B(
        states
    )

    failures += test_frame_specific_U(
        states
    )

    failures += test_symmetric_congruence(
        states
    )

    failures += test_gcd_reduction(
        states
    )

    failures += test_depth_min(
        states
    )

    failures += test_v2U_residual_structure(
        states
    )

    failures += test_U_vs_pq(
        states
    )

    test_possible_n_only_reduction(
        states
    )

    failures += test_exact_final_formula(
        states
    )

    print_examples(
        states
    )

    print_summary()

    print()
    print("=" * 90)
    print("EXPERIMENT 639 FINISHED")
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
