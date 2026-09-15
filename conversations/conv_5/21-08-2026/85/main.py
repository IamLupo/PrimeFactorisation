#!/usr/bin/env python3

from __future__ import annotations

from collections import defaultdict, Counter
from math import isqrt


# =============================================================================
# EXPERIMENT 647
# PRIME-SELECTION INSIDE THE 2-ADIC FIBER
# =============================================================================

INF = 10**9


# =============================================================================
# 2-ADIC UTILITIES
# =============================================================================

def v2(x: int) -> int:
    if x == 0:
        return INF

    x = abs(x)
    return (x & -x).bit_length() - 1


def frame_from_n(n: int) -> str:
    r = n & 3

    if r == 3:
        return "A"

    if r == 1:
        return "B"

    raise ValueError(f"unexpected odd n mod 4={r}")


def frame_data(frame: str):
    if frame == "A":
        return 3, -3, 9

    return -1, 3, 3


# =============================================================================
# PRIME SIEVE
# =============================================================================

def sieve(limit: int) -> list[int]:
    a = bytearray(b"\x01") * (limit + 1)

    a[0] = 0
    a[1] = 0

    for p in range(2, isqrt(limit) + 1):
        if a[p]:
            start = p * p
            a[start::p] = b"\x00" * (
                (limit - start) // p + 1
            )

    return [
        p
        for p in range(3, limit + 1, 2)
        if a[p]
    ]


# =============================================================================
# SEMIPRIME STATES
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

            m = min(v2(A), v2(B))

            states.append(
                (
                    n,
                    p,
                    q,
                    frame,
                    A,
                    B,
                    m,
                )
            )

    return states


# =============================================================================
# THEORETICAL FIBER
# =============================================================================

def fiber_max(residue: int, frame: str, k: int) -> int:
    _, _, c = frame_data(frame)

    x = (residue + c) % (1 << k)

    if x == 0:
        return k

    return min(v2(x), k)


def fiber(residue: int, frame: str, k: int) -> set[int]:
    T = fiber_max(
        residue,
        frame,
        k,
    )

    return set(range(1, T + 1))


# =============================================================================
# ACTUAL FACTOR-LIFT RESIDUE
# =============================================================================

def normalized_factor_residue(
    p: int,
    frame: str,
    k: int,
):
    modulus = 1 << k

    p_res = p % modulus

    if frame == "A":
        # A = p-3
        return (p_res - 3) % modulus

    # B-side uses p+1
    return (p_res + 1) % modulus


def q_residue_from_n(
    n: int,
    p: int,
    k: int,
) -> int:
    modulus = 1 << k

    p_res = p % modulus

    return (
        n
        * pow(
            p_res,
            -1,
            modulus,
        )
    ) % modulus


# =============================================================================
# TEST 1
# FACTOR-LIFT RESIDUE -> m_k
# =============================================================================

def test_factor_lift_signature(
    states,
    k: int,
) -> int:

    print("=" * 90)
    print(
        f"TEST 1: FACTOR-LIFT SIGNATURE k={k}"
    )
    print("=" * 90)

    failures = 0

    modulus = 1 << k

    # signature:
    #
    #   (frame, n mod 2^k, p mod 2^k)
    #
    # uniquely determines q mod 2^k.

    for (
        n,
        p,
        q,
        frame,
        A,
        B,
        m,
    ) in states:

        r = n % modulus
        p_res = p % modulus

        q_res = q_residue_from_n(
            n,
            p,
            k,
        )

        if (p_res * q_res) % modulus != r:
            failures += 1

            if failures <= 20:
                print(
                    "    mismatch"
                    f" n={n}"
                    f" frame={frame}"
                    f" k={k}"
                    f" p_res={p_res}"
                    f" q_res={q_res}"
                    f" residue={r}"
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
# ACTUAL m_k FROM FACTOR RESIDUES
# =============================================================================

def test_m_from_factor_residue(
    states,
    k: int,
) -> int:

    print("=" * 90)
    print(
        f"TEST 2: m_k FROM FACTOR RESIDUES k={k}"
    )
    print("=" * 90)

    failures = 0
    modulus = 1 << k

    for (
        n,
        p,
        q,
        frame,
        A,
        B,
        m,
    ) in states:

        p_res = p % modulus
        q_res = q % modulus

        p0, q0, _ = frame_data(frame)

        A_res = (
            p_res - p0
        ) % modulus

        B_res = (
            q_res - q0
        ) % modulus

        ma = (
            k
            if A_res == 0
            else v2(A_res)
        )

        mb = (
            k
            if B_res == 0
            else v2(B_res)
        )

        mk = min(
            m,
            k,
        )

        predicted_mk = min(
            ma,
            mb,
        )

        if predicted_mk != mk:

            failures += 1

            if failures <= 20:
                print(
                    "    mismatch"
                    f" n={n}"
                    f" frame={frame}"
                    f" k={k}"
                    f" actual_mk={mk}"
                    f" predicted={predicted_mk}"
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
# PRIME RESIDUE ALONE DETERMINES m_k
# =============================================================================

def test_prime_branch_map(
    states,
    k: int,
):

    print("=" * 90)
    print(
        f"TEST 3: PRIME BRANCH -> TRUNCATED m k={k}"
    )
    print("=" * 90)

    modulus = 1 << k

    buckets = defaultdict(set)

    examples = defaultdict(list)

    for (
        n,
        p,
        q,
        frame,
        A,
        B,
        m,
    ) in states:

        signature = (
            frame,
            n % modulus,
            p % modulus,
        )

        mk = min(
            m,
            k,
        )

        buckets[
            signature
        ].add(
            mk
        )

        if len(
            examples[signature]
        ) < 3:

            examples[signature].append(
                (
                    n,
                    p,
                    q,
                    m,
                )
            )

    ambiguous = {
        sig: vals
        for sig, vals in buckets.items()
        if len(vals) > 1
    }

    print(
        f"signatures={len(buckets)}"
    )
    print(
        f"ambiguous={len(ambiguous)}"
    )

    if ambiguous:

        print(
            "    first collisions:"
        )

        for sig, vals in list(
            ambiguous.items()
        )[:20]:

            print(
                f"        {sig}"
                f" -> {sorted(vals)}"
            )

            for x in examples[sig]:
                print(
                    f"            n={x[0]}"
                    f" p={x[1]}"
                    f" q={x[2]}"
                    f" m={x[3]}"
                )

    print()

    return len(ambiguous)


# =============================================================================
# TEST 4
# DOES n RESIDUE + PRIME BRANCH EXPLAIN THE ACTUAL m?
# =============================================================================

def test_branch_vs_fiber(
    states,
    k: int,
):

    print("=" * 90)
    print(
        f"TEST 4: ACTUAL BRANCH INSIDE FIBER k={k}"
    )
    print("=" * 90)

    modulus = 1 << k

    distribution = defaultdict(
        Counter
    )

    for (
        n,
        p,
        q,
        frame,
        A,
        B,
        m,
    ) in states:

        r = n % modulus

        mk = min(
            m,
            k,
        )

        T = fiber_max(
            r,
            frame,
            k,
        )

        # Relative position inside the fiber:
        #
        #   mk / T
        #
        # is not itself expected to be constant.
        #
        # Instead collect the complete distribution.

        distribution[
            (
                frame,
                T,
            )
        ][
            mk
        ] += 1

    print(
        "    branch distributions:"
    )

    for key in sorted(
        distribution
    )[:60]:

        frame, T = key

        print(
            f"        frame={frame}"
            f" T={T}"
            f" distribution="
            f"{dict(distribution[key])}"
        )

    print()

    return 0


# =============================================================================
# TEST 5
# PRIME SIZE / RESIDUE EFFECT
# =============================================================================

def test_prime_size_effect(
    states,
    k: int,
):

    print("=" * 90)
    print(
        f"TEST 5: PRIME-SIZE EFFECT k={k}"
    )
    print("=" * 90)

    modulus = 1 << k

    buckets = defaultdict(
        Counter
    )

    for (
        n,
        p,
        q,
        frame,
        A,
        B,
        m,
    ) in states:

        r = n % modulus
        T = fiber_max(
            r,
            frame,
            k,
        )

        mk = min(
            m,
            k,
        )

        p_class = min(
            p.bit_length(),
            16,
        )

        buckets[
            (
                frame,
                T,
                p_class,
            )
        ][
            mk
        ] += 1

    print(
        "    selected branch statistics:"
    )

    count = 0

    for key in sorted(
        buckets
    ):

        if count >= 80:
            break

        print(
            f"        {key}"
            f" -> {dict(buckets[key])}"
        )

        count += 1

    print()

    return 0


# =============================================================================
# TEST 6
# CONSTRUCTIVE BRANCH REPRESENTATIVES
# =============================================================================

def construct_branch(
    n: int,
    frame: str,
    m: int,
    k: int,
):

    modulus = 1 << k
    c = frame_data(frame)[2]

    if m > k:
        return None

    # We need p and q such that:
    #
    # FRAME A:
    #   p == 3 mod 2^m
    #   q == -3 mod 2^m
    #
    # FRAME B:
    #   p == -1 mod 2^m
    #   q == 3 mod 2^m
    #
    # while pq == n mod 2^k.

    p0, q0, _ = frame_data(frame)

    base = 1 << m

    possible_p = p0 % base

    if possible_p % 2 == 0:
        possible_p = (
            possible_p
            + base
        ) % (
            2 * base
        )

    # Lift p through odd residues modulo 2^k.
    for p_res in range(
        possible_p,
        modulus,
        base,
    ):

        if p_res % 2 == 0:
            continue

        q_res = q_residue_from_n(
            n,
            p_res,
            k,
        )

        A_res = (
            p_res - p0
        ) % modulus

        B_res = (
            q_res - q0
        ) % modulus

        ma = (
            k
            if A_res == 0
            else v2(A_res)
        )

        mb = (
            k
            if B_res == 0
            else v2(B_res)
        )

        if min(
            ma,
            mb,
        ) == m:

            return (
                p_res,
                q_res,
            )

    return None


def test_constructive_branches(
    max_k: int = 10,
):

    print("=" * 90)
    print(
        "TEST 6: CONSTRUCTIVE PRIME-RESIDUE BRANCHES"
    )
    print("=" * 90)

    failures = 0

    for k in (
        6,
        8,
        max_k,
    ):

        modulus = 1 << k

        for frame in (
            "A",
            "B",
        ):

            for residue in range(
                1,
                modulus,
                2,
            ):

                T = fiber_max(
                    residue,
                    frame,
                    k,
                )

                for m in range(
                    1,
                    T + 1,
                ):

                    lift = construct_branch(
                        residue,
                        frame,
                        m,
                        k,
                    )

                    if lift is None:

                        failures += 1

                        if failures <= 20:
                            print(
                                "    no lift"
                                f" frame={frame}"
                                f" k={k}"
                                f" residue={residue}"
                                f" m={m}"
                            )

                        continue

                    p_res, q_res = lift

                    if (
                        p_res
                        * q_res
                    ) % modulus != residue:

                        failures += 1

    print(
        f"failures={failures}"
    )
    print()

    return failures


# =============================================================================
# EXAMPLES
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
        ) = state

        if n not in wanted:
            continue

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
            f" actual_m={m}"
        )

        for k in (
            4,
            6,
            8,
            10,
        ):

            r = n % (
                1 << k
            )

            T = fiber_max(
                r,
                frame,
                k,
            )

            mk = min(
                m,
                k,
            )

            print(
                f"    k={k}"
                f" residue={r}"
                f" m_k={mk}"
                f" fiber=1..{T}"
            )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 90)
    print(
        "EXPERIMENT 647 START"
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

    # -------------------------------------------------------------------------
    # Exact factor-residue theorem
    # -------------------------------------------------------------------------

    for k in (
        6,
        8,
        10,
    ):

        failures += (
            test_factor_lift_signature(
                states,
                k,
            )
        )

        failures += (
            test_m_from_factor_residue(
                states,
                k,
            )
        )

    # -------------------------------------------------------------------------
    # Core branch analysis
    # -------------------------------------------------------------------------

    test_prime_branch_map(
        states,
        10,
    )

    test_branch_vs_fiber(
        states,
        10,
    )

    test_prime_size_effect(
        states,
        10,
    )

    failures += (
        test_constructive_branches(
            max_k=10,
        )
    )

    print_examples(
        states
    )

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------

    print()
    print("=" * 90)
    print(
        "FINAL STRUCTURAL SUMMARY"
    )
    print("=" * 90)

    print(
r"""
Experiment 646 established the existence fiber:

    F_k(n)
      =
    {1,...,T_k}

with:

    FRAME A:
        T_k = min(v2(n+9), k)

    FRAME B:
        T_k = min(v2(n+3), k).

The remaining question is not whether the fiber exists.
It does.

The remaining question is:

    WHICH ELEMENT OF THE FIBER
    is selected by the actual factorization?

For a known factor branch:

    p mod 2^k

the other factor is forced:

    q = n * p^(-1) mod 2^k.

Therefore the signature

    (frame, n mod 2^k, p mod 2^k)

contains exactly the missing branch information.

The experiment tests whether that signature determines:

    m_k = min(m,k)

without ambiguity.

If it does, then the factor-side valuation is not an
independent mysterious quantity. It is exactly a branch
coordinate in the 2-adic factorization fiber.

The most important possible outcome is therefore:

    n residue
        +
    one factor branch
        ->
    m.

That would identify the missing information lost by
projecting the factor pair (p,q) onto n=pq.

A stronger result would show that the branch can be
compressed even further, for example to only:

    v2(p - p0)

or:

    the first differing 2-adic bit of p.

That would turn the entire factor-depth hierarchy into
a finite-state 2-adic branch automaton.
"""
    )

    print()

    print("=" * 90)
    print(
        "EXPERIMENT 647 FINISHED"
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
