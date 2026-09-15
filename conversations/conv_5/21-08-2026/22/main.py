#!/usr/bin/env python3

"""
======================================================================
FULL CLOSED-FORM 2^z THEOREM EXPERIMENT
======================================================================

For an odd factor pair:

    p < q
    N = p*q

Define:

    k_z = 2^(z-1)

    a_z = p mod 2^(z-1)

    y_z = (p+q) / 2^z

    x_z = (q-p+2*a_z) / 2^z

The parameterization is:

    p = k_z (y_z-x_z) + a_z
    q = k_z (y_z+x_z) - a_z

Test proposed closed forms:

    a_z = p mod 2^(z-1)

    y_z = (p+q)/2^z

    x_z = ceil((q-p+2)/2^z)

Validity:

    z is valid iff 2^z divides p+q

Therefore:

    z_max = v2(p+q)

This experiment DOES NOT use recursion to generate x,y.
Every level is calculated independently from p,q.

It then checks whether all formulas agree.

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

TEST_PRIMES_ONLY = False
TEST_COPRIME_ONLY = False

SHOW_FIRST_FAILURES = 20


# ======================================================================
# DATA
# ======================================================================

@dataclass
class LevelResult:
    z: int
    modulus: int

    actual_a: int
    predicted_a: int

    actual_x: int
    predicted_x: int

    actual_y: int
    predicted_y: int

    valid: bool
    passed: bool


@dataclass
class CaseResult:
    p: int
    q: int

    n: int
    sum_pq: int

    v2_sum: int
    predicted_last_z: int
    actual_last_z: int

    levels_tested: int
    failures: int

    passed: bool


# ======================================================================
# v2
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


# ======================================================================
# CEILING DIVISION
# ======================================================================

def ceil_div(a: int, b: int) -> int:
    return (a + b - 1) // b


# ======================================================================
# DIRECT LEVEL
# ======================================================================

def direct_level(
    p: int,
    q: int,
    z: int,
) -> LevelResult | None:

    modulus = 1 << z
    half_modulus = 1 << (z - 1)

    # --------------------------------------------------------------
    # A z-level exists exactly when:
    #
    # p+q ≡ 0 (mod 2^z)
    #
    # --------------------------------------------------------------

    sum_pq = p + q

    if sum_pq % modulus != 0:
        return None

    # --------------------------------------------------------------
    # Actual a_z
    # --------------------------------------------------------------

    actual_a = p % half_modulus

    # This is the definition itself.
    predicted_a = actual_a

    # --------------------------------------------------------------
    # Actual y_z
    # --------------------------------------------------------------

    actual_y = sum_pq // modulus

    # Closed form.
    predicted_y = sum_pq // modulus

    # --------------------------------------------------------------
    # Actual x_z from the parameterization.
    # --------------------------------------------------------------

    numerator_x = (
        q - p + 2 * actual_a
    )

    if numerator_x % modulus != 0:
        return LevelResult(
            z=z,
            modulus=modulus,
            actual_a=actual_a,
            predicted_a=predicted_a,
            actual_x=0,
            predicted_x=0,
            actual_y=actual_y,
            predicted_y=predicted_y,
            valid=False,
            passed=False,
        )

    actual_x = (
        numerator_x // modulus
    )

    # --------------------------------------------------------------
    # Proposed closed form:
    #
    # x_z = ceil((q-p+2)/2^z)
    # --------------------------------------------------------------

    predicted_x = ceil_div(
        q - p + 2,
        modulus,
    )

    # --------------------------------------------------------------
    # Check original factor equations.
    # --------------------------------------------------------------

    k = 1 << (z - 1)

    p_check = (
        k * (actual_y - actual_x)
        + actual_a
    )

    q_check = (
        k * (actual_y + actual_x)
        - actual_a
    )

    reconstruction_ok = (
        p_check == p
        and q_check == q
    )

    # --------------------------------------------------------------
    # Closed-form agreement.
    # --------------------------------------------------------------

    passed = (
        actual_a == predicted_a
        and actual_x == predicted_x
        and actual_y == predicted_y
        and reconstruction_ok
    )

    return LevelResult(
        z=z,
        modulus=modulus,
        actual_a=actual_a,
        predicted_a=predicted_a,
        actual_x=actual_x,
        predicted_x=predicted_x,
        actual_y=actual_y,
        predicted_y=predicted_y,
        valid=True,
        passed=passed,
    )


# ======================================================================
# TEST ONE CASE
# ======================================================================

def test_case(
    p: int,
    q: int,
) -> CaseResult:

    if p > q:
        p, q = q, p

    n = p * q
    sum_pq = p + q

    predicted_last_z = v2(
        sum_pq
    )

    actual_last_z = 1
    levels_tested = 0
    failures = 0

    # --------------------------------------------------------------
    # Test every z up to the theoretical last level.
    # --------------------------------------------------------------

    for z in range(
        2,
        min(MAX_Z, predicted_last_z) + 1,
    ):

        result = direct_level(
            p,
            q,
            z,
        )

        if result is None:
            continue

        levels_tested += 1
        actual_last_z = z

        if not result.passed:
            failures += 1

    # --------------------------------------------------------------
    # Check z_max theorem.
    # --------------------------------------------------------------

    last_z_ok = (
        actual_last_z == predicted_last_z
        if levels_tested > 0
        else True
    )

    if not last_z_ok:
        failures += 1

    return CaseResult(
        p=p,
        q=q,
        n=n,
        sum_pq=sum_pq,
        v2_sum=predicted_last_z,
        predicted_last_z=predicted_last_z,
        actual_last_z=actual_last_z,
        levels_tested=levels_tested,
        failures=failures,
        passed=(failures == 0),
    )


# ======================================================================
# DETAILED CASE
# ======================================================================

def detailed_case(
    p: int,
    q: int,
):

    if p > q:
        p, q = q, p

    result = test_case(
        p,
        q,
    )

    print()
    print("=" * 110)
    print(
        f"N = {result.n} = "
        f"{p} × {q}"
    )
    print("=" * 110)

    print()
    print(
        f"p+q             = {result.sum_pq}"
    )

    print(
        f"v2(p+q)         = {result.v2_sum}"
    )

    print(
        f"predicted z_max = {result.predicted_last_z}"
    )

    print(
        f"actual z_max    = {result.actual_last_z}"
    )

    print(
        f"levels tested    = {result.levels_tested}"
    )

    print()

    print(
        f"{'z':>3} "
        f"{'2^z':>10} "
        f"{'a':>10} "
        f"{'x':>10} "
        f"{'x_closed':>12} "
        f"{'y':>10} "
        f"{'PASS':>7}"
    )

    print("-" * 80)

    for z in range(
        2,
        min(
            MAX_Z,
            result.predicted_last_z
        ) + 1,
    ):

        level = direct_level(
            p,
            q,
            z,
        )

        if level is None:
            continue

        print(
            f"{z:3d} "
            f"{level.modulus:10d} "
            f"{level.actual_a:10d} "
            f"{level.actual_x:10d} "
            f"{level.predicted_x:12d} "
            f"{level.actual_y:10d} "
            f"{str(level.passed):>7}"
        )

    print()
    print(
        f"CASE PASS = {result.passed}"
    )


# ======================================================================
# GENERATE ODD NUMBERS
# ======================================================================

def odd_numbers(
    limit: int,
):
    return range(
        3,
        limit + 1,
        2,
    )


# ======================================================================
# SIMPLE PRIME TEST
# ======================================================================

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


# ======================================================================
# BULK TEST
# ======================================================================

def bulk_test():

    print()
    print("=" * 110)
    print("BULK CLOSED-FORM TEST")
    print("=" * 110)

    total_pairs = 0
    entered_z2 = 0
    failures = []

    for p in odd_numbers(MAX_VALUE):

        if TEST_PRIMES_ONLY and not is_prime(p):
            continue

        for q in odd_numbers(MAX_VALUE):

            if TEST_PRIMES_ONLY and not is_prime(q):
                continue

            if TEST_COPRIME_ONLY and gcd(p, q) != 1:
                continue

            # Avoid duplicate orientations.
            if p > q:
                continue

            total_pairs += 1

            result = test_case(
                p,
                q,
            )

            # N enters z=2 exactly when p+q divisible by 4.
            if result.levels_tested > 0:
                entered_z2 += 1

            if not result.passed:

                failures.append(
                    result
                )

    print()
    print(
        f"odd pairs tested     = {total_pairs}"
    )

    print(
        f"pairs entering z=2   = {entered_z2}"
    )

    print(
        f"failures             = {len(failures)}"
    )

    print(
        f"PASS                 = "
        f"{len(failures) == 0}"
    )

    if failures:

        print()
        print("FIRST FAILURES")
        print("-" * 110)

        for result in failures[
            :SHOW_FIRST_FAILURES
        ]:

            print(
                f"p={result.p} "
                f"q={result.q} "
                f"v2={result.v2_sum} "
                f"last={result.actual_last_z} "
                f"failures={result.failures}"
            )


# ======================================================================
# SPECIAL FACTOR CASES
# ======================================================================

def compare_known_cases():

    print()
    print("=" * 110)
    print("KNOWN CASES")
    print("=" * 110)

    cases = [
        (59, 101),
        (101, 103),
        (223, 449),
        (100127, 100129),
    ]

    for p, q in cases:
        detailed_case(
            p,
            q,
        )


# ======================================================================
# MAIN
# ======================================================================

def main():

    print("=" * 110)
    print(
        "FULL CLOSED-FORM 2^z THEOREM EXPERIMENT"
    )
    print("=" * 110)

    print()
    print(
        "Testing independently:"
    )

    print(
        "  a_z = p mod 2^(z-1)"
    )

    print(
        "  y_z = (p+q) / 2^z"
    )

    print(
        "  x_z = ceil((q-p+2) / 2^z)"
    )

    print(
        "  z_max = v2(p+q)"
    )

    compare_known_cases()

    bulk_test()


if __name__ == "__main__":
    main()
