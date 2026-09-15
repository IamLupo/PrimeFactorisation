#!/usr/bin/env python3

"""
======================================================================
2^z DIRECT DIFFERENCE-OF-SQUARES BRANCH EXPERIMENT
======================================================================

For

    k = 2^(z-1)

    p = k(y-x) + a
    q = k(y+x) - a

we have the exact identity

    N = (k*y)^2 - (k*x-a)^2.

Define

    A = k*y
    B = k*x-a.

Then

    N = A^2-B^2.

For a fixed (z,a), therefore, we can search for A,B directly.

Additional congruence:

    B ≡ -a (mod k)

and

    A ≡ 0 (mod k).

Once A,B are found:

    y = A/k
    x = (B+a)/k

and then

    p = A-B
    q = A+B.

This avoids the previous 5-million-x brute-force loop.

======================================================================
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isqrt
import time

import sympy as sp


# ======================================================================
# CONFIGURATION
# ======================================================================

N_VALUES = [
    5959,
    10403,
    100127,
    10025616383,
]

Z_VALUES = [
    2, 3, 4, 5, 6, 7, 8,
    10, 12, 14, 16, 20, 24,
]

MAX_FERMAT_STEPS = 10_000_000

TIME_LIMIT_SECONDS = 3.0


# ======================================================================
# DATA
# ======================================================================

@dataclass
class Candidate:
    z: int
    a: int
    A: int
    B: int
    x: int
    y: int
    p: int
    q: int
    p_prime: bool
    q_prime: bool


# ======================================================================
# ADMISSIBLE a
# ======================================================================

def admissible_a(
    n: int,
    z: int,
) -> list[int]:

    modulus = 1 << z
    a_modulus = 1 << (z - 1)

    # We only test moderate z here.
    if a_modulus > 2_000_000:
        return []

    target = n % modulus

    return [
        a
        for a in range(a_modulus)
        if (-a * a) % modulus == target
    ]


# ======================================================================
# DIRECT BRANCH SEARCH
# ======================================================================

def search_branch(
    n: int,
    z: int,
    a: int,
) -> list[Candidate]:

    k = 1 << (z - 1)

    start = time.perf_counter()

    # --------------------------------------------------------------
    # We need
    #
    #     N = A² - B²
    #
    # with
    #
    #     A ≡ 0  (mod k)
    #     B ≡ -a (mod k)
    #
    # Let
    #
    #     A = A0 + step*k
    #
    # but rather than scanning every integer A, search Fermat-style
    # and filter the congruence.
    # --------------------------------------------------------------

    A = isqrt(n)

    if A * A < n:
        A += 1

    found: list[Candidate] = []

    steps = 0

    while steps < MAX_FERMAT_STEPS:

        # ----------------------------------------------------------
        # Time guard.
        # ----------------------------------------------------------

        if (
            time.perf_counter() - start
            > TIME_LIMIT_SECONDS
        ):
            break

        B2 = A * A - n

        if B2 >= 0:

            B = isqrt(B2)

            if B * B == B2:

                # --------------------------------------------------
                # Congruence restrictions.
                # --------------------------------------------------

                if A % k == 0 and B % k == (-a) % k:

                    y = A // k
                    x = (B + a) // k

                    if (B + a) % k == 0:

                        p = A - B
                        q = A + B

                        if (
                            p > 1
                            and q > 1
                            and p * q == n
                        ):

                            found.append(
                                Candidate(
                                    z=z,
                                    a=a,
                                    A=A,
                                    B=B,
                                    x=x,
                                    y=y,
                                    p=p,
                                    q=q,
                                    p_prime=bool(
                                        sp.isprime(p)
                                    ),
                                    q_prime=bool(
                                        sp.isprime(q)
                                    ),
                                )
                            )

                            break

        A += 1
        steps += 1

    return found


# ======================================================================
# BETTER CONGRUENCE-ALIGNED SEARCH
# ======================================================================

def search_branch_aligned(
    n: int,
    z: int,
    a: int,
) -> list[Candidate]:

    """
    Instead of scanning every A, write:

        A = k*y

    and use the fact that

        B = k*x-a.

    We search over A in increments of k.

    This is still a Fermat-style search, but aligned to the
    actual 2^z lattice.
    """

    k = 1 << (z - 1)

    # Smallest y for which A >= sqrt(N).
    sqrt_n = isqrt(n)

    y = (
        sqrt_n + k - 1
    ) // k

    start = time.perf_counter()

    for _ in range(MAX_FERMAT_STEPS):

        if (
            time.perf_counter() - start
            > TIME_LIMIT_SECONDS
        ):
            return []

        A = k * y

        if A * A < n:
            y += 1
            continue

        B2 = A * A - n

        B = isqrt(B2)

        if B * B == B2:

            # Required:
            #
            # B = k*x-a
            #
            if (
                B + a
            ) % k == 0:

                x = (
                    B + a
                ) // k

                p = A - B
                q = A + B

                if (
                    p > 1
                    and q > 1
                    and p * q == n
                ):

                    return [
                        Candidate(
                            z=z,
                            a=a,
                            A=A,
                            B=B,
                            x=x,
                            y=y,
                            p=p,
                            q=q,
                            p_prime=bool(
                                sp.isprime(p)
                            ),
                            q_prime=bool(
                                sp.isprime(q)
                            ),
                        )
                    ]

        y += 1

    return []


# ======================================================================
# RUN ONE N
# ======================================================================

def experiment_n(
    n: int,
):

    print()
    print("=" * 110)
    print(f"N = {n}")
    print("=" * 110)

    try:
        print(
            "factorint =",
            sp.factorint(n),
        )
    except Exception:
        pass

    print()
    print(
        "ALIGNED 2^z BRANCH SEARCH"
    )
    print("-" * 110)

    for z in Z_VALUES:

        values = admissible_a(
            n,
            z,
        )

        if not values:
            continue

        print()
        print(
            f"z={z:2d} "
            f"2^z={1 << z:<10d} "
            f"a={values}"
        )

        for a in values:

            start = time.perf_counter()

            results = search_branch_aligned(
                n,
                z,
                a,
            )

            elapsed = (
                time.perf_counter()
                - start
            )

            if not results:

                print(
                    f"  a={a:<8d} "
                    f"NO SOLUTION "
                    f"time={elapsed:.6f}s"
                )

                continue

            for candidate in results:

                print(
                    f"  a={a:<8d} "
                    f"A={candidate.A:<12d} "
                    f"B={candidate.B:<12d} "
                    f"x={candidate.x:<10d} "
                    f"y={candidate.y:<10d} "
                    f"p={candidate.p:<12d} "
                    f"q={candidate.q:<12d} "
                    f"prime=("
                    f"{candidate.p_prime},"
                    f"{candidate.q_prime}) "
                    f"time={elapsed:.6f}s"
                )


# ======================================================================
# VERIFY THE IDENTITY DIRECTLY
# ======================================================================

def verify_identity(
    p: int,
    q: int,
):

    print()
    print("=" * 110)
    print(
        f"IDENTITY CHECK: {p} × {q}"
    )
    print("=" * 110)

    n = p * q

    for z in range(
        2,
        10,
    ):

        k = 1 << (z - 1)

        a = p % k

        # Direct coordinates.
        y_num = p + q
        x_num = q - p + 2 * a

        if (
            y_num % (1 << z)
            != 0
        ):
            continue

        if (
            x_num % (1 << z)
            != 0
        ):
            continue

        x = x_num // (1 << z)
        y = y_num // (1 << z)

        A = k * y
        B = k * x - a

        print()
        print(
            f"z={z}"
        )

        print(
            f"  k       = {k}"
        )

        print(
            f"  a       = {a}"
        )

        print(
            f"  x       = {x}"
        )

        print(
            f"  y       = {y}"
        )

        print(
            f"  A=ky    = {A}"
        )

        print(
            f"  B=kx-a  = {B}"
        )

        print(
            f"  A²-B²   = "
            f"{A*A-B*B}"
        )

        print(
            f"  N       = {n}"
        )

        print(
            f"  PASS    = "
            f"{A*A-B*B == n}"
        )


# ======================================================================
# MAIN
# ======================================================================

def main():

    print("=" * 110)
    print(
        "2^z CONGRUENCE-ALIGNED "
        "DIFFERENCE-OF-SQUARES SEARCH"
    )
    print("=" * 110)

    verify_identity(
        59,
        101,
    )

    verify_identity(
        223,
        449,
    )

    for n in N_VALUES:
        experiment_n(n)


if __name__ == "__main__":
    main()
