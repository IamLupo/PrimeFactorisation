#!/usr/bin/env python3

# ==============================================================================
# EXPERIMENT 643
# ==============================================================================
#
# 2-ADIC FACTOR-LIFT INFORMATION LOSS
#
# EXPERIMENT 642 established:
#
#   d = gcd(A,B)
#
#   depth =
#       v2(d)                  opposite parity
#       v2(d) + 1              both odd
#
# equivalently:
#
#   depth = v2(gcd(A,B)) + [v2(A)=v2(B)].
#
# This experiment studies ONLY the remaining quantity:
#
#   m = v2(gcd(A,B))
#
# and asks:
#
#   At each modulus 2^k,
#   how much of m is actually encoded by n mod 2^k?
#
# Instead of trying arbitrary n+c expressions, we reconstruct the
# complete factor-residue compatibility classes modulo 2^k.
#
# For each observed (frame, n mod 2^k), we determine the set of
# possible residual-minimum valuations:
#
#   m_k = min(v2(A), v2(B))
#
# truncated at k.
#
# This gives an exact information-theoretic description of the
# 2-adic projection:
#
#   factor pair (p,q)
#           |
#           v
#      (A,B)
#           |
#           v
#        m=v2(gcd(A,B))
#           |
#           X
#           |
#           v
#      n mod 2^k
#
# The experiment determines whether the n-projection preserves
# m at a given bit depth, and produces explicit collisions whenever
# it does not.
#
# IMPORTANT:
#
# We work with RESIDUE CLASSES, not the full integer factors.
#
# For a fixed modulus M=2^k:
#
#   q = n * p^{-1} mod M
#
# because p is odd and therefore invertible mod M.
#
# The factor-side m is then computed from the residues.
#
# ==============================================================================


from __future__ import annotations

from collections import defaultdict, Counter
from dataclasses import dataclass
from math import gcd


# ==============================================================================
# CONFIGURATION
# ==============================================================================

PRIME_LIMIT = 6000

# Keep this moderate: the residue-class analysis grows exponentially.
K_LEVELS = [
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    11,
    12,
    13,
    14,
    15,
    16,
]

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
    m: int


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

def ggcd(a: int, b: int) -> int:
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
# RESIDUALS
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
# X,Y
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
            f"frame={frame} "
            f"A={A} "
            f"B={B}"
        )

    return (
        nx // 2,
        ny // 2,
    )


# ==============================================================================
# PRIME SIEVE
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

        if not sieve[p]:
            continue

        start = p * p

        sieve[
            start:
            limit + 1:
            p
        ] = (
            b"\x00"
            * (
                (limit - start)
                // p
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

            d = ggcd(A, B)

            if d == 0:

                m = INF

            else:

                m = v2(d)

            depth = (
                1
                + v2(
                    ggcd(
                        X,
                        Y,
                    )
                )
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
                    m=m,
                )
            )

    return states


# ==============================================================================
# MODULAR RESIDUAL VALUATION
# ==============================================================================

def residue_v2(
    x: int,
    k: int,
) -> int:
    """
    v2(x) known from x modulo 2^k.

    Returns:
        exact value when < k
        k when x == 0 mod 2^k

    Thus this is a truncated valuation.
    """

    M = 1 << k

    r = x % M

    if r == 0:

        return k

    return v2(r)


# ==============================================================================
# RESIDUE M
# ==============================================================================

def residue_m(
    frame: str,
    p_mod: int,
    q_mod: int,
    k: int,
) -> int:

    if frame == "A":

        A = p_mod - 3
        B = q_mod + 3

    else:

        A = p_mod + 1
        B = q_mod - 3

    va = residue_v2(
        A,
        k,
    )

    vb = residue_v2(
        B,
        k,
    )

    return min(
        va,
        vb,
    )


# ==============================================================================
# BUILD POSSIBLE m SETS
# ==============================================================================

def build_possible_m_sets(
    k: int,
) -> dict[tuple[str, int], set[int]]:

    M = 1 << k

    print(
        f"    building residue classes "
        f"for 2^{k}={M}"
    )

    possible = defaultdict(set)

    # Every odd residue is invertible modulo 2^k.
    #
    # We enumerate p residues and determine q from:
    #
    #     q = n * p^{-1} mod M
    #
    # Thus every compatible factor residue pair is visited.
    #
    # Complexity:
    #
    #     O(2^(2k))
    #
    # for a single full table.
    #
    # Therefore we stop at moderate k.

    odd_residues = list(
        range(
            1,
            M,
            2,
        )
    )

    inverse = {
        p: pow(
            p,
            -1,
            M,
        )
        for p in odd_residues
    }

    for n_mod in range(
        M
    ):

        if (n_mod & 1) == 0:
            continue

        # FRAME A / B are determined by n mod 4.
        frame = (
            "A"
            if n_mod % 4 == 3
            else "B"
        )

        key = (
            frame,
            n_mod,
        )

        bucket = possible[key]

        for p_mod in odd_residues:

            q_mod = (
                n_mod
                * inverse[p_mod]
            ) % M

            m = residue_m(
                frame,
                p_mod,
                q_mod,
                k,
            )

            bucket.add(m)

    return possible


# ==============================================================================
# TEST 1
# ==============================================================================

def test_modular_information(
    states: list[State],
    k: int,
) -> dict:

    print("=" * 90)
    print(
        f"TEST 1: EXACT m INFORMATION AT "
        f"2^{k}"
    )
    print("=" * 90)

    possible = build_possible_m_sets(
        k
    )

    ambiguous = 0
    exact = 0

    max_width = 0

    for key, values in possible.items():

        width = len(values)

        max_width = max(
            max_width,
            width,
        )

        if width == 1:
            exact += 1
        else:
            ambiguous += 1

    print(
        f"    residue signatures="
        f"{len(possible)}"
    )

    print(
        f"    exact signatures="
        f"{exact}"
    )

    print(
        f"    ambiguous signatures="
        f"{ambiguous}"
    )

    print(
        f"    maximum m-set width="
        f"{max_width}"
    )

    # --------------------------------------------------------------
    # Compare every actual state
    # --------------------------------------------------------------

    state_ambiguous = 0
    state_incompatible = 0

    for s in states:

        key = (
            s.frame,
            s.n % (1 << k),
        )

        values = possible.get(
            key,
            set(),
        )

        if not values:

            state_incompatible += 1
            continue

        if len(values) > 1:

            state_ambiguous += 1

        # The actual truncated m must be present.
        expected = min(
            s.m,
            k,
        )

        if expected not in values:

            state_incompatible += 1

    print(
        f"    states with ambiguous m="
        f"{state_ambiguous}"
    )

    print(
        f"    incompatible states="
        f"{state_incompatible}"
    )

    print()

    return possible


# ==============================================================================
# TEST 2
# ==============================================================================

def test_m_ambiguity_frontier(
    states: list[State],
    possible: dict,
    k: int,
) -> None:

    print("=" * 90)
    print(
        "TEST 2: ACTUAL m VS MODULAR POSSIBILITY SET"
    )
    print("=" * 90)

    signature_examples = {}

    for s in states:

        key = (
            s.frame,
            s.n % (1 << k),
        )

        values = possible.get(
            key,
            set(),
        )

        if len(values) <= 1:
            continue

        signature_examples.setdefault(
            (
                key,
                tuple(
                    sorted(values)
                ),
            ),
            [],
        ).append(s)

    shown = 0

    for (
        key,
        values,
    ), bucket in sorted(
        signature_examples.items(),
        key=lambda item: (
            item[0][0],
            item[0][1],
        ),
    ):

        print(
            f"    signature={key}"
        )

        print(
            f"        possible_m="
            f"{list(values)}"
        )

        for s in bucket[:4]:

            print(
                f"        n={s.n:<10} "
                f"p={s.p:<6} "
                f"q={s.q:<6} "
                f"actual_m={s.m}"
            )

        print()

        shown += 1

        if shown >= 15:
            break


# ==============================================================================
# TEST 3
# ==============================================================================

def test_truncated_m_prediction(
    states: list[State],
    possible: dict,
    k: int,
) -> int:

    print("=" * 90)
    print(
        "TEST 3: CAN n mod 2^k DETERMINE "
        "THE TRUNCATED m?"
    )
    print("=" * 90)

    failures = 0

    for s in states:

        key = (
            s.frame,
            s.n % (1 << k),
        )

        values = possible.get(
            key,
            set(),
        )

        actual = min(
            s.m,
            k,
        )

        if len(values) == 1:

            predicted = next(
                iter(values)
            )

            if predicted != actual:
                failures += 1

    print(
        f"    checked={len(states)} "
        f"failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 4
# ==============================================================================

def test_residue_lift_monotonicity(
    states: list[State],
) -> int:

    print("=" * 90)
    print(
        "TEST 4: RESIDUE-LIFT MONOTONICITY"
    )
    print("=" * 90)

    failures = 0

    previous = None

    for k in K_LEVELS:

        if k > 13:
            # Avoid rebuilding enormous full tables in this
            # independent consistency check.
            break

        possible = build_possible_m_sets(
            k
        )

        ambiguous = sum(
            1
            for values in possible.values()
            if len(values) > 1
        )

        if previous is not None:

            if ambiguous > previous:
                failures += 1

                print(
                    f"    anomaly at k={k}: "
                    f"previous={previous} "
                    f"current={ambiguous}"
                )

        previous = ambiguous

        print(
            f"    k={k} "
            f"ambiguous={ambiguous}"
        )

    print(
        f"    failures={failures}"
    )
    print()

    return failures


# ==============================================================================
# TEST 5
# ==============================================================================

def test_m_vs_nplusc(
    states: list[State],
) -> int:

    print("=" * 90)
    print(
        "TEST 5: m VS v2(n+c)"
    )
    print("=" * 90)

    failures = 0

    distributions = Counter()

    for s in states:

        if s.frame == "A":
            c = 9
        else:
            c = 3

        nv = v2(
            s.n + c
        )

        if s.m >= INF:
            continue

        delta = nv - s.m

        distributions[
            (
                s.frame,
                delta,
            )
        ] += 1

        if delta < 0:
            failures += 1

    print(
        f"checked={len(states)} "
        f"failures={failures}"
    )

    print(
        "    delta=v2(n+c)-m:"
    )

    for key in sorted(
        distributions
    ):

        print(
            f"        {key}: "
            f"{distributions[key]}"
        )

    print()

    return failures


# ==============================================================================
# TEST 6
# ==============================================================================

def test_depth_reconstruction_from_m(
    states: list[State],
) -> int:

    print("=" * 90)
    print(
        "TEST 6: DEPTH = m + EQUALITY"
    )
    print("=" * 90)

    failures = 0

    for s in states:

        if s.A == 0:

            predicted = v2(
                abs(s.B)
            )

        elif s.B == 0:

            predicted = v2(
                abs(s.A)
            )

        else:

            d = ggcd(
                s.A,
                s.B,
            )

            m = v2(d)

            equality = (
                v2(s.A)
                ==
                v2(s.B)
            )

            predicted = (
                m
                + int(equality)
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
# TEST 7
# ==============================================================================

def test_low_bit_collision_examples(
    states: list[State],
    k: int,
) -> None:

    print("=" * 90)
    print(
        "TEST 7: EXPLICIT n-RESIDUE COLLISIONS"
    )
    print("=" * 90)

    buckets = defaultdict(list)

    M = 1 << k

    for s in states:

        key = (
            s.frame,
            s.n % M,
        )

        buckets[key].append(s)

    shown = 0

    for key, bucket in sorted(
        buckets.items(),
        key=lambda item: (
            item[0][0],
            item[0][1],
        ),
    ):

        ms = sorted(
            {
                min(
                    s.m,
                    k,
                )
                for s in bucket
            }
        )

        if len(ms) <= 1:
            continue

        print(
            f"    signature={key}"
        )

        print(
            f"        truncated_m={ms}"
        )

        for s in bucket:

            print(
                f"        n={s.n:<10} "
                f"p={s.p:<6} "
                f"q={s.q:<6} "
                f"m={s.m:<3} "
                f"depth={s.depth}"
            )

            if bucket.index(s) >= 5:
                break

        print()

        shown += 1

        if shown >= 10:
            break


# ==============================================================================
# EXAMPLES
# ==============================================================================

def print_examples(
    states: list[State],
) -> None:

    print("=" * 90)
    print(
        "EXAMPLES"
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

    for n in examples:

        if n not in lookup:
            continue

        s = lookup[n]

        c = (
            9
            if s.frame == "A"
            else 3
        )

        print(
            f"n={s.n} "
            f"p={s.p} "
            f"q={s.q} "
            f"frame={s.frame}"
        )

        print(
            f"    A={s.A} "
            f"B={s.B}"
        )

        print(
            f"    m=v2(gcd(A,B))="
            f"{s.m}"
        )

        print(
            f"    v2(n+c)="
            f"{v2(s.n+c)} "
            f"c={c}"
        )

        print(
            f"    depth={s.depth}"
        )

        print()


# ==============================================================================
# SUMMARY
# ==============================================================================

def print_summary():

    print("=" * 90)
    print(
        "FINAL STRUCTURAL SUMMARY"
    )
    print("=" * 90)

    print(
r"""
Experiment 642 established the exact normalized theorem:

    d = gcd(A,B)

    A=d*a
    B=d*b
    gcd(a,b)=1.

Only three normalized parity states exist:

    (a,b)=(0,1)
    (a,b)=(1,0)
    (a,b)=(1,1).

For both frames:

    opposite parity
        ->
    depth=v2(d).

    both odd
        ->
    depth=v2(d)+1.

Therefore:

    depth =
        v2(gcd(A,B))
        +
        [v2(A)=v2(B)].

The ONLY remaining factor-side quantity is:

    m=v2(gcd(A,B)).

Experiment 643 does not search arbitrary polynomial
expressions in n.

Instead it measures exactly how much information
about m survives in:

    n mod 2^k.

For fixed frame and residue n mod 2^k, every odd
factor residue p determines:

    q = n*p^{-1} mod 2^k.

Thus the full compatible 2-adic factor-pair space
can be enumerated exactly.

For each compatible pair we compute the truncated:

    m_k =
        min(
            v2(A),
            v2(B),
            k
        ).

This gives the exact set of m-values consistent with
the observable:

    (frame, n mod 2^k).

Interpretation:

    one possible m
        ->
    n's k low bits determine m_k exactly.

    multiple possible m
        ->
    the factor-side residual valuation is already
    non-unique at that modulus.

This is a structural information-loss test rather than
another arbitrary n+c signature search.
"""
    )


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print(
        "EXPERIMENT 643 START"
    )
    print("=" * 90)
    print()
    print(
        "2-ADIC FACTOR-LIFT INFORMATION LOSS"
    )
    print()

    print(
        "[1] PRIME SIEVE"
    )

    primes = prime_sieve(
        PRIME_LIMIT
    )

    print(
        f"    odd primes={len(primes)}"
    )
    print()

    print(
        "[2] SEMIPRIME GENERATION"
    )

    states = build_states(
        primes
    )

    print(
        f"    semiprimes={len(states)}"
    )
    print()

    total_failures = 0

    total_failures += test_m_vs_nplusc(
        states
    )

    total_failures += test_depth_reconstruction_from_m(
        states
    )

    # --------------------------------------------------------------
    # Main modular-information experiment.
    #
    # Full enumeration becomes expensive at high k.
    # Keep the default target moderate.
    # --------------------------------------------------------------

    target_k = 10

    possible = test_modular_information(
        states,
        target_k,
    )

    total_failures += test_truncated_m_prediction(
        states,
        possible,
        target_k,
    )

    test_m_ambiguity_frontier(
        states,
        possible,
        target_k,
    )

    test_low_bit_collision_examples(
        states,
        target_k,
    )

    # --------------------------------------------------------------
    # Explicitly verify the projection gets more informative
    # as k increases.
    # --------------------------------------------------------------

    total_failures += test_residue_lift_monotonicity(
        states
    )

    print_examples(
        states
    )

    print_summary()

    print()
    print("=" * 90)
    print(
        "EXPERIMENT 643 FINISHED"
    )
    print("=" * 90)
    print()

    print(
        f"TOTAL FAILURES={total_failures}"
    )

    if total_failures == 0:

        print(
            "STATUS=ALL CONSISTENCY TESTS PASSED"
        )

    else:

        print(
            "STATUS=COUNTEREXAMPLES FOUND"
        )


if __name__ == "__main__":
    main()
