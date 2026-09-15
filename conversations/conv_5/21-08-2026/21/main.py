#!/usr/bin/env python3

"""
======================================================================
FERMAT vs 2^z RECURSIVE STRUCTURE EXPERIMENT
======================================================================

For odd factors p < q:

    N = p*q

Fermat representation:

    N = A^2 - B^2

where

    A = (p+q)/2
    B = (q-p)/2

Our z=2 coordinates satisfy:

    A = 2y
    B = 2x - 1

and therefore:

    N = (2y)^2 - (2x-1)^2

The recursive 2^z system is:

    a' = a + (x mod 2)*2^(z-1)

    x' = ceil(x/2)

    y' = y/2

This experiment compares:

    1. Fermat's search distance.
    2. x_2.
    3. Number of valid 2^z levels.
    4. v2(p+q).
    5. Closed-form recursive coordinates.

The key question:

    Is the 2^z hierarchy computationally stronger than
    ordinary Fermat difference-of-squares, or is it a
    structured re-expression of the same factorization?

======================================================================
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isqrt, gcd

import sympy as sp


# ======================================================================
# CONFIGURATION
# ======================================================================

PRIME_LIMIT = 500

MAX_FERMAT_STEPS = 10_000_000

SHOW_CASES = 30


# ======================================================================
# DATA
# ======================================================================

@dataclass
class FermatResult:
    n: int
    A: int
    B: int
    p: int
    q: int
    iterations: int


@dataclass
class RecursiveResult:
    x2: int
    y2: int
    valid_levels: int
    transitions: int
    v2_sum: int
    last_z: int


# ======================================================================
# v2
# ======================================================================

def v2(n: int) -> int:

    n = abs(n)

    if n == 0:
        return 10**9

    result = 0

    while n % 2 == 0:

        n //= 2
        result += 1

    return result


# ======================================================================
# FERMAT
# ======================================================================

def fermat_factor(
    n: int,
) -> FermatResult | None:
    """
    Standard Fermat difference-of-squares search.

        N = A^2 - B^2
          = (A-B)(A+B)

    """

    A = isqrt(n)

    if A * A < n:
        A += 1

    for iterations in range(
        MAX_FERMAT_STEPS + 1
    ):

        B2 = A * A - n

        if B2 >= 0:

            B = isqrt(B2)

            if B * B == B2:

                p = A - B
                q = A + B

                if (
                    p > 1
                    and q > 1
                    and p * q == n
                ):

                    return FermatResult(
                        n=n,
                        A=A,
                        B=B,
                        p=min(p, q),
                        q=max(p, q),
                        iterations=iterations,
                    )

        A += 1

    return None


# ======================================================================
# z=2 COORDINATES
# ======================================================================

def z2_coordinates(
    p: int,
    q: int,
):
    """
    Our representation:

        p = 2y - 2x + 1
        q = 2y + 2x - 1

    """

    # Orientation-independent:
    # choose p < q.

    if p > q:
        p, q = q, p

    y = (p + q) // 4
    x = (q - p + 2) // 4

    if (
        (p + q) % 4 != 0
        or (q - p + 2) % 4 != 0
    ):
        return None

    return x, y


# ======================================================================
# RECURSIVE LEVELS
# ======================================================================

def recursive_levels(
    p: int,
    q: int,
):
    """
    Compute actual valid z-levels directly from p,q.
    """

    if p > q:
        p, q = q, p

    levels = []

    for z in range(
        2,
        64,
    ):

        modulus = 1 << z
        k = 1 << (z - 1)

        a = p % k

        y_num = p + q
        x_num = q - p + 2 * a

        if y_num % modulus != 0:
            break

        if x_num % modulus != 0:
            break

        x = x_num // modulus
        y = y_num // modulus

        p_check = (
            k * (y - x) + a
        )

        q_check = (
            k * (y + x) - a
        )

        if p_check != p or q_check != q:
            raise AssertionError(
                "parameterization check failed"
            )

        levels.append(
            (z, a, x, y)
        )

    return levels


# ======================================================================
# CASE ANALYSIS
# ======================================================================

def analyze(
    p: int,
    q: int,
):
    """

    Analyze one prime pair.
    """

    if p > q:
        p, q = q, p

    n = p * q

    fermat = fermat_factor(n)

    coords = z2_coordinates(
        p,
        q,
    )

    levels = recursive_levels(
        p,
        q,
    )

    if coords is None:

        return {
            "n": n,
            "p": p,
            "q": q,
            "fermat": fermat,
            "coords": None,
            "levels": levels,
            "enters_z2": False,
        }

    x2, y2 = coords

    sum_v2 = v2(
        p + q
    )

    transitions = max(
        0,
        len(levels) - 1
    )

    return {
        "n": n,
        "p": p,
        "q": q,
        "fermat": fermat,
        "coords": coords,
        "levels": levels,
        "enters_z2": True,
        "x2": x2,
        "y2": y2,
        "v2_sum": sum_v2,
        "transitions": transitions,
        "last_z": (
            levels[-1][0]
            if levels
            else None
        ),
    }


# ======================================================================
# PRINT CASE
# ======================================================================

def print_case(
    result,
):

    print()
    print("=" * 100)
    print(
        f"N={result['n']} "
        f"= {result['p']} × {result['q']}"
    )
    print("=" * 100)

    print()

    if result["fermat"] is not None:

        f = result["fermat"]

        print(
            f"Fermat A       = {f.A}"
        )

        print(
            f"Fermat B       = {f.B}"
        )

        print(
            f"Fermat steps   = {f.iterations}"
        )

    else:

        print(
            "Fermat search did not finish."
        )

    print()

    if not result["enters_z2"]:

        print(
            "z=2 parameterization: NO"
        )

        print(
            f"N mod 4 = "
            f"{result['n'] % 4}"
        )

        return

    x2 = result["x2"]
    y2 = result["y2"]

    print(
        f"x2             = {x2}"
    )

    print(
        f"y2             = {y2}"
    )

    print(
        f"v2(p+q)        = "
        f"{result['v2_sum']}"
    )

    print(
        f"valid levels   = "
        f"{len(result['levels'])}"
    )

    print(
        f"transitions    = "
        f"{result['transitions']}"
    )

    print(
        f"last z         = "
        f"{result['last_z']}"
    )

    # --------------------------------------------------------------
    # Fermat correspondence.
    # --------------------------------------------------------------

    fermat = result["fermat"]

    if fermat is not None:

        print()

        print(
            f"2*y2           = {2*y2}"
        )

        print(
            f"2*x2-1         = {2*x2-1}"
        )

        print(
            f"Fermat A       = {fermat.A}"
        )

        print(
            f"Fermat B       = {fermat.B}"
        )

        print(
            f"A match        = "
            f"{2*y2 == fermat.A}"
        )

        print(
            f"B match        = "
            f"{2*x2-1 == fermat.B}"
        )

    # --------------------------------------------------------------
    # Recursive levels.
    # --------------------------------------------------------------

    print()
    print(
        "RECURSIVE LEVELS"
    )

    print("-" * 100)

    for z, a, x, y in result["levels"]:

        print(
            f"z={z:2d} "
            f"a={a:<8d} "
            f"x={x:<8d} "
            f"y={y:<8d}"
        )


# ======================================================================
# BULK EXPERIMENT
# ======================================================================

def bulk():

    primes = list(
        sp.primerange(
            3,
            PRIME_LIMIT + 1
        )
    )

    total = 0

    enters = 0

    fermat_matches = 0

    x_matches = 0

    y_matches = 0

    v2_matches = 0

    last_z_matches = 0

    failures = []

    for i, p in enumerate(primes):

        for q in primes[i + 1:]:

            total += 1

            result = analyze(
                p,
                q,
            )

            if not result[
                "enters_z2"
            ]:
                continue

            enters += 1

            fermat = result[
                "fermat"
            ]

            if fermat is None:
                continue

            # --------------------------------------------------
            # Check Fermat correspondence.
            # --------------------------------------------------

            if (
                2 * result["y2"]
                == fermat.A
                and
                2 * result["x2"] - 1
                == fermat.B
            ):
                fermat_matches += 1
            else:
                failures.append(
                    (
                        p,
                        q,
                        "FERMAT",
                    )
                )

            # --------------------------------------------------
            # Closed forms.
            # --------------------------------------------------

            levels = result[
                "levels"
            ]

            x2 = result["x2"]
            y2 = result["y2"]

            ok_x = True
            ok_y = True

            for z, a, x, y in levels:

                r = z - 2
                divisor = 1 << r

                predicted_x = (
                    x2
                    + divisor
                    - 1
                ) // divisor

                if (
                    predicted_x
                    != x
                ):
                    ok_x = False

                if (
                    y2 % divisor
                    != 0
                ):
                    ok_y = False
                else:

                    predicted_y = (
                        y2 // divisor
                    )

                    if predicted_y != y:
                        ok_y = False

            if ok_x:
                x_matches += 1
            else:
                failures.append(
                    (
                        p,
                        q,
                        "X",
                    )
                )

            if ok_y:
                y_matches += 1
            else:
                failures.append(
                    (
                        p,
                        q,
                        "Y",
                    )
                )

            # --------------------------------------------------
            # v2 theorem.
            # --------------------------------------------------

            if (
                result["transitions"]
                == result["v2_sum"]
            ):
                v2_matches += 1
            else:
                failures.append(
                    (
                        p,
                        q,
                        "V2",
                    )
                )

            if (
                result["last_z"]
                == result["v2_sum"]
            ):
                last_z_matches += 1
            else:
                failures.append(
                    (
                        p,
                        q,
                        "LAST_Z",
                    )
                )

    # --------------------------------------------------------------
    # Summary
    # --------------------------------------------------------------

    print()
    print("=" * 100)
    print(
        "BULK RESULTS"
    )
    print("=" * 100)

    print(
        f"prime pairs tested       = {total}"
    )

    print(
        f"pairs entering z=2       = {enters}"
    )

    print(
        f"Fermat equivalence       = "
        f"{fermat_matches}/{enters}"
    )

    print(
        f"x closed form            = "
        f"{x_matches}/{enters}"
    )

    print(
        f"y closed form            = "
        f"{y_matches}/{enters}"
    )

    print(
        f"transition = v2(p+q)     = "
        f"{v2_matches}/{enters}"
    )

    print(
        f"last z = v2(p+q)         = "
        f"{last_z_matches}/{enters}"
    )

    print(
        f"failures                 = "
        f"{len(failures)}"
    )

    print(
        f"PASS                     = "
        f"{len(failures) == 0}"
    )

    if failures:

        print()
        print(
            "FIRST FAILURES"
        )

        for failure in failures[:20]:

            print(
                failure
            )


# ======================================================================
# MAIN
# ======================================================================

def main():

    print("=" * 100)
    print(
        "FERMAT vs 2^z RECURSIVE STRUCTURE"
    )
    print("=" * 100)

    # --------------------------------------------------------------
    # Detailed examples.
    # --------------------------------------------------------------

    for p, q in [
        (59, 101),
        (101, 103),
        (223, 449),
    ]:

        print_case(
            analyze(
                p,
                q,
            )
        )

    # --------------------------------------------------------------
    # Bulk test.
    # --------------------------------------------------------------

    bulk()


if __name__ == "__main__":
    main()
