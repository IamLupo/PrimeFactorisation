#!/usr/bin/env python3

"""
======================================================================
EXACT CEILING-IDENTITY / REMAINDER EXPERIMENT
======================================================================

For odd p < q and a valid level z:

    m = 2^z
    a = p mod 2^(z-1)

The parameterization is

    p = 2^(z-1)(y-x) + a
    q = 2^(z-1)(y+x) - a

Therefore:

    q - p + 2
        = 2^z*x - 2(a-1)

Define

    r = 2(a-1)

Then:

    q-p+2 = 2^z*x - r

and for odd p:

    1 <= a < 2^(z-1)

so

    0 <= r < 2^z.

Therefore:

    x = ceil((q-p+2)/2^z)

This experiment independently verifies:

    1. exact remainder identity
    2. remainder bound
    3. divisibility identity
    4. ceiling identity
    5. reconstruction of p and q
    6. recursive x transition
    7. a_z residue identity

No recursive formula is used to calculate x_z.
======================================================================
"""

from __future__ import annotations

from dataclasses import dataclass
from math import gcd


# ======================================================================
# CONFIGURATION
# ======================================================================

MAX_VALUE = 1000
MAX_Z = 32

TEST_COPRIME_ONLY = False
TEST_PRIMES_ONLY = False

SHOW_FIRST_FAILURES = 20


# ======================================================================
# DATA
# ======================================================================

@dataclass
class Level:
    z: int
    modulus: int
    a: int
    x: int
    y: int

    remainder: int
    quotient: int
    ceiling_x: int

    remainder_ok: bool
    exact_identity_ok: bool
    divisibility_ok: bool
    ceiling_ok: bool

    passed: bool


# ======================================================================
# BASIC HELPERS
# ======================================================================

def v2(n: int) -> int:
    if n == 0:
        return 10**9

    n = abs(n)
    result = 0

    while n % 2 == 0:
        n //= 2
        result += 1

    return result


def is_prime(n: int) -> bool:

    if n < 2:
        return False

    if n % 2 == 0:
        return n == 2

    d = 3

    while d * d <= n:

        if n % d == 0:
            return False

        d += 2

    return True


def ceil_div(a: int, b: int) -> int:

    return (a + b - 1) // b


# ======================================================================
# DIRECT LEVEL
# ======================================================================

def direct_level(
    p: int,
    q: int,
    z: int,
) -> Level | None:

    modulus = 1 << z
    half_modulus = 1 << (z - 1)

    # A valid level requires:
    #
    #   p+q ≡ 0 (mod 2^z)

    if (p + q) % modulus != 0:
        return None

    # --------------------------------------------------------------
    # a_z
    # --------------------------------------------------------------

    a = p % half_modulus

    # --------------------------------------------------------------
    # y_z
    # --------------------------------------------------------------

    y = (p + q) // modulus

    # --------------------------------------------------------------
    # x_z
    # --------------------------------------------------------------

    numerator = (
        q - p + 2 * a
    )

    if numerator % modulus != 0:
        raise AssertionError(
            "x numerator not divisible"
        )

    x = numerator // modulus

    # --------------------------------------------------------------
    # Exact remainder expression.
    #
    # q-p+2 = 2^z*x - 2(a-1)
    # --------------------------------------------------------------

    lhs = q - p + 2

    remainder = 2 * (a - 1)

    quotient = lhs // modulus

    ceiling_x = ceil_div(
        lhs,
        modulus,
    )

    # --------------------------------------------------------------
    # Checks
    # --------------------------------------------------------------

    remainder_ok = (
        0 <= remainder < modulus
    )

    exact_identity_ok = (
        lhs
        == modulus * x - remainder
    )

    divisibility_ok = (
        (
            lhs + remainder
        ) % modulus
        == 0
    )

    ceiling_ok = (
        ceiling_x == x
    )

    passed = (
        remainder_ok
        and exact_identity_ok
        and divisibility_ok
        and ceiling_ok
    )

    return Level(
        z=z,
        modulus=modulus,
        a=a,
        x=x,
        y=y,

        remainder=remainder,
        quotient=quotient,
        ceiling_x=ceiling_x,

        remainder_ok=remainder_ok,
        exact_identity_ok=exact_identity_ok,
        divisibility_ok=divisibility_ok,
        ceiling_ok=ceiling_ok,

        passed=passed,
    )


# ======================================================================
# TEST ONE FACTOR PAIR
# ======================================================================

def test_pair(
    p: int,
    q: int,
) -> list[Level]:

    if p > q:
        p, q = q, p

    levels = []

    for z in range(
        2,
        MAX_Z + 1,
    ):

        level = direct_level(
            p,
            q,
            z,
        )

        if level is None:
            break

        levels.append(level)

    return levels


# ======================================================================
# DETAILED CASE
# ======================================================================

def detailed_case(
    p: int,
    q: int,
):

    if p > q:
        p, q = q, p

    levels = test_pair(
        p,
        q,
    )

    print()
    print("=" * 110)
    print(
        f"N = {p*q} = {p} × {q}"
    )
    print("=" * 110)

    print()
    print(
        f"p+q = {p+q}"
    )

    print(
        f"v2(p+q) = {v2(p+q)}"
    )

    print()

    print(
        f"{'z':>3} "
        f"{'a':>10} "
        f"{'x':>10} "
        f"{'2(a-1)':>12} "
        f"{'2^z*x':>14} "
        f"{'q-p+2':>12} "
        f"{'ceil':>8} "
        f"{'PASS':>7}"
    )

    print("-" * 100)

    for L in levels:

        lhs = q - p + 2
        big = L.modulus * L.x

        print(
            f"{L.z:3d} "
            f"{L.a:10d} "
            f"{L.x:10d} "
            f"{L.remainder:12d} "
            f"{big:14d} "
            f"{lhs:12d} "
            f"{L.ceiling_x:8d} "
            f"{str(L.passed):>7}"
        )

    print()

    failures = [
        L
        for L in levels
        if not L.passed
    ]

    print(
        f"levels = {len(levels)}"
    )

    print(
        f"failures = {len(failures)}"
    )

    print(
        f"PASS = {len(failures) == 0}"
    )

    if levels:

        print()
        print(
            "IDENTITY DETAILS"
        )
        print("-" * 100)

        for L in levels:

            print(
                f"z={L.z:2d}: "
                f"{q-p+2} = "
                f"{L.modulus}*{L.x} "
                f"- {L.remainder}"
            )

            print(
                f"       "
                f"0 <= {L.remainder} "
                f"< {L.modulus} : "
                f"{L.remainder_ok}"
            )

            print(
                f"       "
                f"ceil({q-p+2}/{L.modulus}) "
                f"= {L.ceiling_x} "
                f"= x : "
                f"{L.ceiling_ok}"
            )


# ======================================================================
# TEST a-BOUND DIRECTLY
# ======================================================================

def test_a_bound(
    p: int,
    q: int,
):

    levels = test_pair(
        p,
        q,
    )

    failures = []

    for L in levels:

        half = 1 << (L.z - 1)

        condition = (
            1 <= L.a < half
        )

        remainder_condition = (
            0 <= 2 * (L.a - 1) < L.modulus
        )

        if not (
            condition
            and remainder_condition
        ):

            failures.append(
                (
                    L.z,
                    L.a,
                    half,
                    L.modulus,
                )
            )

    return failures


# ======================================================================
# BULK EXPERIMENT
# ======================================================================

def bulk_test():

    print()
    print("=" * 110)
    print("BULK EXACT REMAINDER TEST")
    print("=" * 110)

    total_pairs = 0
    entered_z2 = 0

    total_levels = 0

    identity_failures = 0
    remainder_failures = 0
    divisibility_failures = 0
    ceiling_failures = 0
    a_bound_failures = 0

    first_failures = []

    for p in range(
        3,
        MAX_VALUE + 1,
        2,
    ):

        if TEST_PRIMES_ONLY and not is_prime(p):
            continue

        for q in range(
            p,
            MAX_VALUE + 1,
            2,
        ):

            if TEST_PRIMES_ONLY and not is_prime(q):
                continue

            if (
                TEST_COPRIME_ONLY
                and gcd(p, q) != 1
            ):
                continue

            total_pairs += 1

            levels = test_pair(
                p,
                q,
            )

            if not levels:
                continue

            entered_z2 += 1

            for L in levels:

                total_levels += 1

                if not L.exact_identity_ok:

                    identity_failures += 1

                    if len(first_failures) < SHOW_FIRST_FAILURES:
                        first_failures.append(
                            (
                                p,
                                q,
                                L.z,
                                "identity",
                            )
                        )

                if not L.remainder_ok:

                    remainder_failures += 1

                    if len(first_failures) < SHOW_FIRST_FAILURES:
                        first_failures.append(
                            (
                                p,
                                q,
                                L.z,
                                "remainder",
                            )
                        )

                if not L.divisibility_ok:

                    divisibility_failures += 1

                    if len(first_failures) < SHOW_FIRST_FAILURES:
                        first_failures.append(
                            (
                                p,
                                q,
                                L.z,
                                "divisibility",
                            )
                        )

                if not L.ceiling_ok:

                    ceiling_failures += 1

                    if len(first_failures) < SHOW_FIRST_FAILURES:
                        first_failures.append(
                            (
                                p,
                                q,
                                L.z,
                                "ceiling",
                            )
                        )

                half = 1 << (L.z - 1)

                if not (
                    1 <= L.a < half
                ):

                    a_bound_failures += 1

                    if len(first_failures) < SHOW_FIRST_FAILURES:
                        first_failures.append(
                            (
                                p,
                                q,
                                L.z,
                                "a-bound",
                            )
                        )

    # --------------------------------------------------------------
    # Summary
    # --------------------------------------------------------------

    print()

    print(
        f"pairs tested          = {total_pairs}"
    )

    print(
        f"pairs entering z=2    = {entered_z2}"
    )

    print(
        f"levels tested         = {total_levels}"
    )

    print()

    print(
        f"identity failures     = "
        f"{identity_failures}"
    )

    print(
        f"remainder failures    = "
        f"{remainder_failures}"
    )

    print(
        f"divisibility failures = "
        f"{divisibility_failures}"
    )

    print(
        f"ceiling failures      = "
        f"{ceiling_failures}"
    )

    print(
        f"a-bound failures      = "
        f"{a_bound_failures}"
    )

    total_failures = (
        identity_failures
        + remainder_failures
        + divisibility_failures
        + ceiling_failures
        + a_bound_failures
    )

    print()

    print(
        f"TOTAL FAILURES        = "
        f"{total_failures}"
    )

    print(
        f"PASS                  = "
        f"{total_failures == 0}"
    )

    if first_failures:

        print()
        print(
            "FIRST FAILURES"
        )

        print("-" * 100)

        for failure in first_failures:

            print(
                failure
            )


# ======================================================================
# BIT-LEVEL VIEW OF a
# ======================================================================

def show_a_bits(
    p: int,
    q: int,
):

    levels = test_pair(
        p,
        q,
    )

    print()
    print("=" * 110)
    print(
        f"a_z BIT ACCUMULATION: {p} × {q}"
    )
    print("=" * 110)

    print()

    for L in levels:

        bits = format(
            L.a,
            f"0{L.z - 1}b"
        )

        p_low = (
            p & (L.modulus // 2 - 1)
        )

        print(
            f"z={L.z:2d} "
            f"a={L.a:<8d} "
            f"binary={bits} "
            f"p low bits={p_low:<8d}"
        )


# ======================================================================
# MAIN
# ======================================================================

def main():

    print("=" * 110)
    print(
        "EXACT q-p+2 REMAINDER / CEILING EXPERIMENT"
    )
    print("=" * 110)

    print()
    print(
        "Testing:"
    )

    print(
        "  q-p+2 = 2^z*x_z - 2(a_z-1)"
    )

    print(
        "  0 <= 2(a_z-1) < 2^z"
    )

    print(
        "  x_z = ceil((q-p+2)/2^z)"
    )

    # --------------------------------------------------------------
    # Detailed examples.
    # --------------------------------------------------------------

    detailed_case(
        59,
        101,
    )

    detailed_case(
        223,
        449,
    )

    detailed_case(
        100127,
        100129,
    )

    # --------------------------------------------------------------
    # Show a bits.
    # --------------------------------------------------------------

    show_a_bits(
        223,
        449,
    )

    # --------------------------------------------------------------
    # Bulk.
    # --------------------------------------------------------------

    bulk_test()


if __name__ == "__main__":
    main()
