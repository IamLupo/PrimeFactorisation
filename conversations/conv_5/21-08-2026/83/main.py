#!/usr/bin/env python3

# ==============================================================================
# EXPERIMENT 645
# ==============================================================================
#
# 2-ADIC THRESHOLD IMAGE THEOREM
#
# GOAL
# ----
#
# Experiment 642:
#
#     depth = v2(gcd(A,B)) + [v2(A)=v2(B)]
#
# Experiment 643:
#
#     actual m=v2(gcd(A,B)) is NOT determined by n mod 2^k.
#
# Experiment 644:
#
#     the set of n-residues compatible with m>=t has size 2^(k-t),
#     exactly the size of ONE residue class modulo 2^t.
#
# This experiment proves the exact threshold image algebraically:
#
# FRAME A:
#
#     A = p-3
#     B = q+3
#
#     m >= t
#       <=>
#     p = 3 (mod 2^t)
#     q = -3 (mod 2^t)
#
#     therefore
#
#     n = pq
#       = -9 (mod 2^t).
#
# FRAME B:
#
#     A = p+1
#     B = q-3
#
#     m >= t
#       <=>
#     p = -1 (mod 2^t)
#     q = 3 (mod 2^t)
#
#     therefore
#
#     n = pq
#       = -3 (mod 2^t).
#
# IMPORTANT:
#
# This does NOT say that the actual factorization satisfies
#
#     m = v2(n+9)     FRAME A
#     m = v2(n+3)     FRAME B.
#
# That is false.
#
# Instead it says:
#
#     there EXISTS a compatible factor lift with m>=t
#
# exactly when the n-residue satisfies the corresponding
# congruence.
#
# Thus the n-only observable describes the IMAGE of the
# factor-lift map, not the actual m of a particular factorization.
#
# This distinction is the central information-loss result.
#
# ==============================================================================


from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from math import gcd


# ==============================================================================
# CONSTANTS
# ==============================================================================

INF = 10**9

MAX_T = 16

SHOW = 20


# ==============================================================================
# BASIC 2-ADIC UTILITIES
# ==============================================================================

def v2(x: int) -> int:
    if x == 0:
        return INF

    x = abs(x)

    return (x & -x).bit_length() - 1


def odd_modulus(k: int) -> int:
    return 1 << k


# ==============================================================================
# PRIME SIEVE
# ==============================================================================

def prime_sieve(limit: int) -> list[int]:

    sieve = bytearray(
        b"\x01"
    ) * (
        limit + 1
    )

    sieve[0] = 0
    sieve[1] = 0

    for p in range(
        2,
        int(limit ** 0.5) + 1,
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
            *
            (
                (
                    limit - start
                )
                // p
                + 1
            )
        )

    return [
        p
        for p in range(
            3,
            limit + 1,
            2
        )
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
# m = v2(gcd(A,B))
# ==============================================================================

def residual_m(
    A: int,
    B: int,
) -> int:

    return min(
        v2(A),
        v2(B),
    )


# ==============================================================================
# ACTUAL STATES
# ==============================================================================

@dataclass
class State:

    n: int
    p: int
    q: int

    frame: str

    A: int
    B: int

    m: int

    depth: int


def build_states(
    primes: list[int],
) -> list[State]:

    states = []

    for i, p in enumerate(
        primes
    ):

        for q in primes[i:]:

            n = p * q

            frame = (
                "A"
                if n % 4 == 3
                else "B"
            )

            A, B = residuals(
                frame,
                p,
                q,
            )

            m = residual_m(
                A,
                B,
            )

            if frame == "A":

                x2 = B - A
                y2 = A + B

            else:

                x2 = 3 * A - B
                y2 = 3 * A + B

            depth = min(
                v2(x2),
                v2(y2),
            )

            states.append(
                State(
                    n=n,
                    p=p,
                    q=q,
                    frame=frame,
                    A=A,
                    B=B,
                    m=m,
                    depth=depth,
                )
            )

    return states


# ==============================================================================
# PREDICTED THRESHOLD RESIDUE
# ==============================================================================

def threshold_residue(
    frame: str,
    t: int,
) -> int:

    M = 1 << t

    if frame == "A":

        return (-9) % M

    return (-3) % M


# ==============================================================================
# TEST 1
# ==============================================================================

def test_symbolic_threshold_identity() -> int:

    print("=" * 90)
    print(
        "TEST 1: SYMBOLIC THRESHOLD PRODUCT IDENTITY"
    )
    print("=" * 90)

    failures = 0

    for frame in (
        "A",
        "B",
    ):

        for t in range(
            1,
            MAX_T + 1,
        ):

            M = 1 << t

            if frame == "A":

                p_res = 3 % M
                q_res = (-3) % M
                expected = (-9) % M

            else:

                p_res = (-1) % M
                q_res = 3 % M
                expected = (-3) % M

            actual = (
                p_res * q_res
            ) % M

            if actual != expected:

                failures += 1

                print(
                    f"    mismatch frame={frame} "
                    f"t={t} "
                    f"actual={actual} "
                    f"expected={expected}"
                )

    print(
        f"checked="
        f"{2 * MAX_T}"
    )

    print(
        f"failures={failures}"
    )

    print()

    return failures


# ==============================================================================
# TEST 2
# ==============================================================================

def test_exact_residue_image() -> int:

    print("=" * 90)
    print(
        "TEST 2: EXACT THRESHOLD RESIDUE IMAGE"
    )
    print("=" * 90)

    failures = 0

    for frame in (
        "A",
        "B",
    ):

        for t in range(
            1,
            MAX_T + 1,
        ):

            M = 1 << t

            target = threshold_residue(
                frame,
                t,
            )

            # Every n in the predicted class must have a
            # compatible factor pair at level t.
            #
            # We choose one allowed p residue and recover q.
            #
            # Frame A:
            #     p = 3
            #
            # Frame B:
            #     p = -1.
            #
            if frame == "A":
                p = 3 % M
            else:
                p = (-1) % M

            inv_p = pow(
                p,
                -1,
                M,
            )

            for n in range(
                1,
                M,
                2,
            ):

                possible = (
                    n == target
                )

                q = (
                    n
                    * inv_p
                ) % M

                if frame == "A":

                    compatible = (
                        p % M == 3 % M
                        and
                        q % M == (-3) % M
                    )

                else:

                    compatible = (
                        p % M == (-1) % M
                        and
                        q % M == 3 % M
                    )

                if compatible != possible:

                    failures += 1

                    if failures <= SHOW:

                        print(
                            f"    mismatch "
                            f"frame={frame} "
                            f"t={t} "
                            f"n={n} "
                            f"possible={possible} "
                            f"compatible={compatible}"
                        )

    print(
        f"failures={failures}"
    )

    print()

    return failures


# ==============================================================================
# TEST 3
# ==============================================================================

def test_actual_threshold_equivalence(
    states: list[State],
) -> int:

    print("=" * 90)
    print(
        "TEST 3: ACTUAL FACTOR STATE -> "
        "n THRESHOLD CONGRUENCE"
    )
    print("=" * 90)

    failures = 0

    checked = 0

    for s in states:

        frame = s.frame

        for t in range(
            1,
            min(
                s.m + 2,
                MAX_T + 1,
            ),
        ):

            checked += 1

            target = threshold_residue(
                frame,
                t,
            )

            lhs = (
                s.m >= t
            )

            rhs = (
                s.n % (1 << t)
                == target
            )

            if lhs != rhs:

                failures += 1

                if failures <= SHOW:

                    print(
                        f"    mismatch "
                        f"n={s.n} "
                        f"p={s.p} "
                        f"q={s.q} "
                        f"frame={frame} "
                        f"m={s.m} "
                        f"t={t} "
                        f"lhs={lhs} "
                        f"rhs={rhs}"
                    )

    print(
        f"checked={checked}"
    )

    print(
        f"failures={failures}"
    )

    print()

    return failures


# ==============================================================================
# TEST 4
# ==============================================================================

def test_threshold_equivalence_from_direct_residuals(
    states: list[State],
) -> int:

    print("=" * 90)
    print(
        "TEST 4: DIRECT RESIDUAL CONGRUENCE "
        "<-> n CONGRUENCE"
    )
    print("=" * 90)

    failures = 0
    checked = 0

    for s in states:

        frame = s.frame

        for t in range(
            1,
            MAX_T + 1,
        ):

            checked += 1

            M = 1 << t

            if frame == "A":

                residual_condition = (
                    s.A % M == 0
                    and
                    s.B % M == 0
                )

            else:

                residual_condition = (
                    s.A % M == 0
                    and
                    s.B % M == 0
                )

            n_condition = (
                s.n % M
                ==
                threshold_residue(
                    frame,
                    t,
                )
            )

            if residual_condition != n_condition:

                failures += 1

                if failures <= SHOW:

                    print(
                        f"    mismatch "
                        f"n={s.n} "
                        f"frame={frame} "
                        f"t={t} "
                        f"A={s.A} "
                        f"B={s.B}"
                    )

    print(
        f"checked={checked}"
    )

    print(
        f"failures={failures}"
    )

    print()

    return failures


# ==============================================================================
# TEST 5
# ==============================================================================

def test_image_sizes() -> int:

    print("=" * 90)
    print(
        "TEST 5: IMAGE SIZE THEOREM"
    )
    print("=" * 90)

    failures = 0

    for frame in (
        "A",
        "B",
    ):

        for k in range(
            3,
            MAX_T + 1,
        ):

            M = 1 << k

            for t in range(
                1,
                k + 1,
            ):

                # Number of odd residues modulo 2^k satisfying
                #
                #     n == r (mod 2^t)
                #
                # is exactly:
                #
                #     2^(k-t).
                #
                expected = 1 << (
                    k - t
                )

                actual = sum(
                    1
                    for n in range(
                        1,
                        M,
                        2,
                    )
                    if (
                        n
                        %
                        (1 << t)
                    )
                    ==
                    threshold_residue(
                        frame,
                        t,
                    )
                )

                if actual != expected:

                    failures += 1

                    print(
                        f"    mismatch "
                        f"frame={frame} "
                        f"k={k} "
                        f"t={t} "
                        f"actual={actual} "
                        f"expected={expected}"
                    )

    print(
        f"failures={failures}"
    )

    print()

    return failures


# ==============================================================================
# TEST 6
# ==============================================================================

def test_nested_thresholds() -> int:

    print("=" * 90)
    print(
        "TEST 6: NESTED 2-ADIC THRESHOLD CLASSES"
    )
    print("=" * 90)

    failures = 0

    for frame in (
        "A",
        "B",
    ):

        print(
            f"FRAME {frame}"
        )

        for t in range(
            1,
            MAX_T,
        ):

            r1 = threshold_residue(
                frame,
                t,
            )

            r2 = threshold_residue(
                frame,
                t + 1,
            )

            if (
                r2
                %
                (1 << t)
            ) != r1:

                failures += 1

                print(
                    f"    mismatch "
                    f"t={t} "
                    f"r_t={r1} "
                    f"r_t+1={r2}"
                )

            else:

                print(
                    f"    t={t:<2} "
                    f"r_t={r1:<6} "
                    f"r_(t+1)={r2:<6} "
                    f"consistent=True"
                )

    print(
        f"failures={failures}"
    )

    print()

    return failures


# ==============================================================================
# TEST 7
# ==============================================================================

def test_actual_information_loss(
    states: list[State],
) -> int:

    print("=" * 90)
    print(
        "TEST 7: SAME n-RESIDUE CAN HAVE DIFFERENT m"
    )
    print("=" * 90)

    buckets = defaultdict(
        dict
    )

    for s in states:

        key = (
            s.frame,
            s.n,
        )

        buckets[
            key
        ].setdefault(
            s.m,
            s,
        )

    collisions = []

    for key, values in buckets.items():

        if len(values) > 1:

            collisions.append(
                (
                    key,
                    values,
                )
            )

    print(
        f"same-n factorization buckets="
        f"{len(collisions)}"
    )

    shown = 0

    for (
        (frame, n),
        values,
    ) in collisions:

        if shown >= SHOW:
            break

        print(
            f"    frame={frame} "
            f"n={n} "
            f"m-values="
            f"{sorted(values)}"
        )

        for m, s in sorted(
            values.items()
        ):

            print(
                f"        m={m} "
                f"p={s.p} "
                f"q={s.q} "
                f"depth={s.depth}"
            )

        print()

        shown += 1

    print()

    print(
        "NOTE: same-n collisions are impossible inside a "
        "single semiprime state list when n determines one "
        "factorization pair; the meaningful collision is "
        "same residue modulo 2^k across different n."
    )

    print()

    return 0


# ==============================================================================
# TEST 8
# ==============================================================================

def test_residue_collisions(
    states: list[State],
    k: int,
) -> int:

    print("=" * 90)
    print(
        f"TEST 8: SAME n mod 2^{k} "
        f"WITH DIFFERENT ACTUAL m"
    )
    print("=" * 90)

    M = 1 << k

    buckets = defaultdict(
        list
    )

    for s in states:

        buckets[
            (
                s.frame,
                s.n % M,
            )
        ].append(s)

    ambiguous = []

    for key, group in buckets.items():

        ms = {
            s.m
            for s in group
        }

        if len(ms) > 1:

            ambiguous.append(
                (
                    key,
                    sorted(ms),
                    group,
                )
            )

    print(
        f"ambiguous residue buckets="
        f"{len(ambiguous)}"
    )

    for (
        (frame, residue),
        ms,
        group,
    ) in ambiguous[:SHOW]:

        print(
            f"    frame={frame} "
            f"residue={residue} "
            f"m={ms}"
        )

        shown_states = 0

        for s in group:

            print(
                f"        n={s.n:<8} "
                f"p={s.p:<5} "
                f"q={s.q:<5} "
                f"m={s.m:<3} "
                f"depth={s.depth}"
            )

            shown_states += 1

            if shown_states >= 4:
                break

    print()

    return 0


# ==============================================================================
# TEST 9
# ==============================================================================

def test_primary_valuation_relation(
    states: list[State],
) -> int:

    print("=" * 90)
    print(
        "TEST 9: m <= v2(n+c)"
    )
    print("=" * 90)

    failures = 0

    for s in states:

        c = (
            9
            if s.frame == "A"
            else 3
        )

        vn = v2(
            s.n + c
        )

        if s.m > vn:

            failures += 1

            if failures <= SHOW:

                print(
                    f"    violation "
                    f"n={s.n} "
                    f"frame={s.frame} "
                    f"m={s.m} "
                    f"v2(n+c)={vn}"
                )

    print(
        f"checked={len(states)}"
    )

    print(
        f"failures={failures}"
    )

    print()

    return failures


# ==============================================================================
# TEST 10
# ==============================================================================

def print_threshold_examples(
    states: list[State],
) -> None:

    print("=" * 90)
    print(
        "EXAMPLES"
    )
    print("=" * 90)
    print()

    wanted = [
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

    lookup = {
        s.n: s
        for s in states
    }

    for n in wanted:

        s = lookup[n]

        c = (
            9
            if s.frame == "A"
            else 3
        )

        vn = v2(
            n + c
        )

        print(
            f"n={n} "
            f"p={s.p} "
            f"q={s.q} "
            f"frame={s.frame}"
        )

        print(
            f"    A={s.A} "
            f"B={s.B}"
        )

        print(
            f"    actual m={s.m}"
        )

        print(
            f"    v2(n+c)={vn} "
            f"c={c}"
        )

        print(
            f"    threshold "
            f"m>=t iff "
            f"n == "
            f"{threshold_residue(s.frame, max(1, min(s.m, MAX_T)))} "
            f"(mod 2^{max(1, min(s.m, MAX_T))})"
        )

        print(
            f"    depth={s.depth}"
        )

        print()


# ==============================================================================
# MAIN
# ==============================================================================

def main():

    print("=" * 90)
    print(
        "EXPERIMENT 645 START"
    )
    print("=" * 90)
    print()
    print(
        "2-ADIC THRESHOLD IMAGE THEOREM"
    )
    print()

    print(
        "[1] PRIME SIEVE"
    )

    primes = prime_sieve(
        6000
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

    total_failures += (
        test_symbolic_threshold_identity()
    )

    total_failures += (
        test_exact_residue_image()
    )

    total_failures += (
        test_actual_threshold_equivalence(
            states
        )
    )

    total_failures += (
        test_threshold_equivalence_from_direct_residuals(
            states
        )
    )

    total_failures += (
        test_image_sizes()
    )

    total_failures += (
        test_nested_thresholds()
    )

    # This is deliberately descriptive rather than a
    # failure test.
    test_actual_information_loss(
        states
    )

    # Use a modulus small enough to run quickly.
    test_residue_collisions(
        states,
        k=10,
    )

    total_failures += (
        test_primary_valuation_relation(
            states
        )
    )

    print_threshold_examples(
        states
    )

    print("=" * 90)
    print(
        "FINAL STRUCTURAL SUMMARY"
    )
    print("=" * 90)

    print(
r"""
The exact threshold theorem is:

FRAME A:

    m >= t
        <=>
    2^t | (p-3)
        AND
    2^t | (q+3)

        <=>

    p == 3       (mod 2^t)
    q == -3      (mod 2^t)

        =>

    n = pq
      == -9      (mod 2^t).

FRAME B:

    m >= t
        <=>
    2^t | (p+1)
        AND
    2^t | (q-3)

        <=>

    p == -1      (mod 2^t)
    q == 3       (mod 2^t)

        =>

    n = pq
      == -3      (mod 2^t).

Conversely, the residue condition gives an explicit
compatible factor lift:

FRAME A:

    p = 3
    q = n * 3^(-1)        (mod 2^t).

FRAME B:

    p = -1
    q = -n               (mod 2^t).

Therefore the threshold image is exactly:

FRAME A:

    I_t(A)
      =
    { n :
      n == -9 (mod 2^t) }.

FRAME B:

    I_t(B)
      =
    { n :
      n == -3 (mod 2^t) }.

This explains Experiment 644 exactly:

    |I_k(m)| = 2^(k-m).

The images are nested:

    I_(t+1) subset I_t.

But this does NOT imply:

    actual m = v2(n+c).

It only implies:

    m >= t
        can occur
    iff
    v2(n+c) >= t.

Hence:

    m <= v2(n+c)

for the actual factorization,

while different factorizations/residue realizations
can have smaller m than v2(n+c).

The remaining information-loss question is therefore
sharper:

    n gives an upper bound on m,

    but the actual m depends on which factor lift
    (p,q) realizes n.

This is fundamentally different from searching for
another clever n-only valuation formula.

The next target should therefore be the FIBER:

    F(n,k)
      =
    { m :
      there exists a compatible factor lift
      modulo 2^k }.

The natural question is whether this fiber has a simple
form such as:

    {1,2,...,v2(n+c)}

or a frame-dependent restricted subset.

If so, we will have an exact description of all the
information that n loses about the factor-side hierarchy.
"""
    )

    print()

    print("=" * 90)
    print(
        "EXPERIMENT 645 FINISHED"
    )
    print("=" * 90)
    print()

    print(
        f"TOTAL FAILURES={total_failures}"
    )

    if total_failures == 0:

        print(
            "STATUS=ALL STRUCTURAL TESTS PASSED"
        )

    else:

        print(
            "STATUS=COUNTEREXAMPLES FOUND"
        )


if __name__ == "__main__":
    main()
