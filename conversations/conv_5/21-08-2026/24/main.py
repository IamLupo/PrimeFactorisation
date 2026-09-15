#!/usr/bin/env python3

"""
======================================================================
N-ONLY 2^z PERFECT-SQUARE FACTOR SEARCH
======================================================================

We use

    p = k(y-x) + a
    q = k(y+x) - a

with

    k = 2^(z-1)

and

    N = p*q

giving

    N = k^2(y^2-x^2) + 2akx - a^2

so

    y^2
      = x^2
        + (N + a^2 - 2akx) / k^2

For a fixed (N,z,a), a valid factor pair exists iff

    T(x)
      = x^2
        + (N + a^2 - 2akx) / k^2

is an integer perfect square.

This experiment:

    1. Finds admissible a from N mod 2^z.
    2. Searches x using the exact square condition.
    3. Never uses known factors to construct x/y.
    4. Tests several z values.
    5. Compares candidate counts.
    6. Checks primality of recovered factors.

It also tests the "strongest z" strategy:

    choose the largest z for which N ≡ -a² (mod 2^z)

and solve there first.

======================================================================
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from math import isqrt

import sympy as sp


# ======================================================================
# CONFIGURATION
# ======================================================================

N_VALUES = [
    5959,
    3233,
    10403,
    18721,
    100127,
    10025616383,
]

# z values tested explicitly
Z_VALUES = [
    2, 3, 4, 5, 6, 7, 8,
    10, 12, 14, 16, 20, 24
]

# Maximum x values checked for a single branch.
# Increase only after small tests work.
MAX_X = 5_000_000

# Maximum number of solutions printed per branch.
MAX_RESULTS = 20

# Optional time protection.
MAX_SECONDS_PER_BRANCH = 5.0


# ======================================================================
# DATA
# ======================================================================

@dataclass
class Candidate:
    z: int
    a: int
    x: int
    y: int
    p: int
    q: int
    p_prime: bool
    q_prime: bool


@dataclass
class BranchResult:
    z: int
    a: int
    x_tested: int
    candidates: list[Candidate]
    elapsed: float
    stopped_by_limit: bool


# ======================================================================
# MODULAR HELPERS
# ======================================================================

def admissible_a(
    n: int,
    z: int,
) -> list[int]:
    """
    Find a modulo 2^(z-1) satisfying

        n ≡ -a² (mod 2^z)
    """

    modulus = 1 << z
    a_modulus = 1 << (z - 1)

    target = n % modulus

    # Keep this bounded for the experiment.
    if a_modulus > 2_000_000:
        return []

    return [
        a
        for a in range(a_modulus)
        if (-a * a) % modulus == target
    ]


# ======================================================================
# FACTORS FROM STATE
# ======================================================================

def factors_from_state(
    z: int,
    a: int,
    x: int,
    y: int,
):
    k = 1 << (z - 1)

    p = (
        k * (y - x)
        + a
    )

    q = (
        k * (y + x)
        - a
    )

    return p, q


# ======================================================================
# EXACT BRANCH SEARCH
# ======================================================================

def search_branch(
    n: int,
    z: int,
    a: int,
) -> BranchResult:
    """
    Search for integer x such that

        y² =
            x²
            + (N + a² - 2akx)/k²

    is a perfect square.
    """

    k = 1 << (z - 1)
    k2 = k * k

    start = time.perf_counter()

    candidates: list[Candidate] = []

    x_tested = 0
    stopped_by_limit = False

    for x in range(
        0,
        MAX_X + 1,
    ):

        x_tested += 1

        # ----------------------------------------------------------
        # Time protection
        # ----------------------------------------------------------

        if (
            time.perf_counter() - start
            >= MAX_SECONDS_PER_BRANCH
        ):

            stopped_by_limit = True
            break

        # ----------------------------------------------------------
        # Numerator for y²
        # ----------------------------------------------------------

        numerator = (
            n
            + a * a
            - 2 * a * k * x
        )

        if numerator <= 0:
            continue

        if numerator % k2 != 0:
            continue

        y2 = (
            x * x
            + numerator // k2
        )

        if y2 <= 0:
            continue

        y = isqrt(y2)

        if y * y != y2:
            continue

        # We normally want positive distinct factors.
        if y <= x:
            continue

        p, q = factors_from_state(
            z,
            a,
            x,
            y,
        )

        if p <= 1 or q <= 1:
            continue

        if p * q != n:
            continue

        candidate = Candidate(
            z=z,
            a=a,
            x=x,
            y=y,
            p=p,
            q=q,
            p_prime=bool(sp.isprime(p)),
            q_prime=bool(sp.isprime(q)),
        )

        candidates.append(candidate)

        if len(candidates) >= MAX_RESULTS:
            break

    elapsed = (
        time.perf_counter() - start
    )

    return BranchResult(
        z=z,
        a=a,
        x_tested=x_tested,
        candidates=candidates,
        elapsed=elapsed,
        stopped_by_limit=stopped_by_limit,
    )


# ======================================================================
# STRONGEST Z
# ======================================================================

def strongest_z(
    n: int,
    max_z: int = 32,
):
    """
    Find the largest z<=max_z having at least one admissible a.
    """

    best = None

    for z in range(
        2,
        max_z + 1,
    ):

        values = admissible_a(
            n,
            z,
        )

        if values:
            best = (
                z,
                values,
            )

    return best


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

    print()
    print(
        f"bits = {n.bit_length()}"
    )

    try:
        print(
            "factorint =",
            sp.factorint(n),
        )
    except Exception:
        pass

    # --------------------------------------------------------------
    # All tested z values
    # --------------------------------------------------------------

    print()
    print(
        "ADMISSIBLE z/a"
    )
    print("-" * 110)

    branches = []

    for z in Z_VALUES:

        values = admissible_a(
            n,
            z,
        )

        if not values:
            continue

        print(
            f"z={z:2d} "
            f"n mod 2^z={n % (1 << z):<10d} "
            f"a={values}"
        )

        for a in values:
            branches.append(
                (z, a)
            )

    # --------------------------------------------------------------
    # Search every branch
    # --------------------------------------------------------------

    print()
    print(
        "BRANCH SEARCH"
    )
    print("-" * 110)

    total_candidates = []

    for z, a in branches:

        result = search_branch(
            n,
            z,
            a,
        )

        print(
            f"z={z:2d} "
            f"a={a:<8d} "
            f"x_tested={result.x_tested:<10d} "
            f"time={result.elapsed:.4f}s "
            f"stopped={result.stopped_by_limit} "
            f"solutions={len(result.candidates)}"
        )

        for candidate in result.candidates:

            total_candidates.append(
                candidate
            )

            print(
                f"    x={candidate.x:<10d} "
                f"y={candidate.y:<10d} "
                f"p={candidate.p:<12d} "
                f"q={candidate.q:<12d} "
                f"prime=("
                f"{candidate.p_prime},"
                f"{candidate.q_prime})"
            )

    # --------------------------------------------------------------
    # Strongest z
    # --------------------------------------------------------------

    print()
    print(
        "STRONGEST-z TEST"
    )
    print("-" * 110)

    strongest = strongest_z(
        n,
        max_z=24,
    )

    if strongest is None:

        print(
            "No admissible z."
        )

        return

    z, values = strongest

    print(
        f"strongest z = {z}"
    )

    print(
        f"a values     = {values}"
    )

    for a in values:

        result = search_branch(
            n,
            z,
            a,
        )

        print(
            f"  z={z} "
            f"a={a} "
            f"x_tested={result.x_tested} "
            f"time={result.elapsed:.4f}s "
            f"solutions={len(result.candidates)}"
        )

        for candidate in result.candidates:

            print(
                f"    p={candidate.p} "
                f"q={candidate.q} "
                f"x={candidate.x} "
                f"y={candidate.y} "
                f"prime=("
                f"{candidate.p_prime},"
                f"{candidate.q_prime})"
            )


# ======================================================================
# DIRECT EXAMPLE
# ======================================================================

def verify_known_example():

    print()
    print("=" * 110)
    print(
        "KNOWN 223 × 449 EXAMPLE"
    )
    print("=" * 110)

    n = 223 * 449

    for z in [2, 3, 4, 5]:

        values = admissible_a(
            n,
            z,
        )

        print()
        print(
            f"z={z} "
            f"n mod 2^z={n % (1 << z)} "
            f"a={values}"
        )

        for a in values:

            result = search_branch(
                n,
                z,
                a,
            )

            for candidate in result.candidates:

                print(
                    f"  a={a} "
                    f"x={candidate.x} "
                    f"y={candidate.y} "
                    f"p={candidate.p} "
                    f"q={candidate.q}"
                )


# ======================================================================
# MAIN
# ======================================================================

def main():

    print("=" * 110)
    print(
        "N-ONLY 2^z PERFECT-SQUARE SEARCH"
    )
    print("=" * 110)

    print()
    print(
        "Equation:"
    )

    print(
        "  y² = x² + "
        "(N + a² - 2*a*k*x)/k²"
    )

    print()
    print(
        "where k = 2^(z-1)"
    )

    verify_known_example()

    for n in N_VALUES:

        experiment_n(n)


if __name__ == "__main__":
    main()
