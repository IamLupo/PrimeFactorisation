#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 247 — FAST ORDER-GROWTH EXCEPTION / NONEXCEPTIONAL CLASSIFICATION
==============================================================================

Experiment 244 found one exceptional prime <= 200:

    ell = 11,

with

    ord_11(3)=5,
    ord_121(3)=5,
    ord_1331(3)=55.

Experiment 246 confirmed that this delayed order growth genuinely changes
the exponent coordinate at the exceptional prime.

Experiment 247 therefore isolates the order theory itself.

For primes ell <= PRIME_LIMIT, compute

    d = ord_ell(3)

and

    D_ell = (3^d - 1)/ell mod ell.

Classify:

    NONEXCEPTIONAL:
        D_ell != 0 mod ell

    EXCEPTIONAL:
        D_ell == 0 mod ell.

For each exceptional prime, verify the actual order tower

    ord_{ell^r}(3)

through a small number of levels and record the first level at which
the order grows by ell.

For each nonexceptional prime, verify the expected law

    ord_{ell^r}(3)
      = d*ell^(r-1).

This experiment is independent of q1,q3. It identifies the arithmetic
exception set that must be excluded or treated separately in the final
local theorem.

Only main.py is used.
"""

from __future__ import annotations

import sys


PRIME_LIMIT = 2000
MAX_LEVEL = 4


# ============================================================================
# HELPERS
# ============================================================================

def is_prime(n: int) -> bool:

    if n < 2:
        return False

    if n == 2:
        return True

    if n % 2 == 0:
        return False

    d = 3

    while d * d <= n:
        if n % d == 0:
            return False
        d += 2

    return True


def factorize_small(
    n: int,
) -> list[int]:

    """
    Return prime divisors of n, enough for order computation.
    """

    factors = []

    d = 2
    x = n

    while d * d <= x:

        if x % d == 0:

            factors.append(d)

            while x % d == 0:
                x //= d

        d = 3 if d == 2 else d + 2

    if x > 1:
        factors.append(x)

    return factors


def multiplicative_order(
    base: int,
    modulus: int,
    phi: int | None = None,
) -> int:

    """
    Exact multiplicative order.

    If phi is supplied, factor phi and reduce the order efficiently.
    """

    if phi is None:
        # For prime modulus p.
        phi = modulus - 1

    order = int(phi)

    for prime_factor in factorize_small(order):

        while (
            order % prime_factor == 0
            and
            pow(
                base,
                order // prime_factor,
                modulus,
            ) == 1
        ):
            order //= prime_factor

    return int(order)


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 247 — FAST ORDER-GROWTH EXCEPTION / "
        "NONEXCEPTIONAL CLASSIFICATION"
    )
    print("=" * 78)

    exceptional = []
    nonexceptional = []
    failures = []

    primes = [
        p
        for p in range(
            2,
            PRIME_LIMIT + 1,
        )
        if is_prime(p)
        and p != 3
    ]

    print()
    print("=" * 78)
    print("1. PRIME CLASSIFICATION")
    print("=" * 78)

    for ell in primes:

        d = multiplicative_order(
            3,
            ell,
            ell - 1,
        )

        modulus2 = ell * ell

        increment = (
            pow(
                3,
                d,
                modulus2,
            )
            - 1
        )

        if increment % ell != 0:

            failures.append(
                (
                    ell,
                    "increment_not_divisible",
                )
            )

            continue

        D = int(
            (
                increment
                // ell
            )
            % ell
        )

        is_exceptional = (
            D == 0
        )

        if is_exceptional:
            exceptional.append(
                (
                    ell,
                    d,
                    D,
                )
            )
        else:
            nonexceptional.append(
                (
                    ell,
                    d,
                    D,
                )
            )

        if is_exceptional:
            print(
                f"  ell={ell}: "
                f"EXCEPTIONAL "
                f"ord1={d} "
                f"D={D}"
            )

    # ------------------------------------------------------------------
    # 2. NONEXCEPTIONAL CHECK
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. NONEXCEPTIONAL ORDER TOWERS")
    print("=" * 78)

    for ell, d, D in nonexceptional:

        actual_orders = []
        expected_orders = []

        previous_order = d
        actual_orders.append(d)
        expected_orders.append(d)

        ok = True

        for level in range(
            2,
            MAX_LEVEL + 1,
        ):

            modulus = ell ** level

            phi = (
                ell - 1
            ) * ell ** (
                level - 1
            )

            actual = multiplicative_order(
                3,
                modulus,
                phi,
            )

            expected = (
                d
                * ell ** (
                    level - 1
                )
            )

            actual_orders.append(actual)
            expected_orders.append(expected)

            if actual != expected:
                ok = False

        if not ok:
            failures.append(
                (
                    ell,
                    "nonexceptional_order_growth_failure",
                    actual_orders,
                    expected_orders,
                )
            )

        # Only print failures in the large nonexceptional set.
        if not ok:
            print(
                f"  ell={ell}: "
                f"actual={actual_orders} "
                f"expected={expected_orders}"
            )

    # ------------------------------------------------------------------
    # 3. EXCEPTIONAL CHECK
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. EXCEPTIONAL ORDER TOWERS")
    print("=" * 78)

    for ell, d, D in exceptional:

        actual_orders = []
        growth_ratios = []

        for level in range(
            1,
            MAX_LEVEL + 1,
        ):

            modulus = ell ** level

            if level == 1:
                order = d
            else:
                phi = (
                    ell - 1
                ) * ell ** (
                    level - 1
                )

                order = multiplicative_order(
                    3,
                    modulus,
                    phi,
                )

            actual_orders.append(order)

            if level == 1:
                growth_ratios.append(
                    None
                )
            else:
                previous = actual_orders[
                    level - 2
                ]

                growth_ratios.append(
                    order // previous
                )

        # First level at which the order grows.
        first_growth_level = None

        for level in range(
            2,
            MAX_LEVEL + 1,
        ):

            if (
                actual_orders[level - 1]
                >
                actual_orders[level - 2]
            ):

                first_growth_level = level
                break

        print(
            f"  ell={ell}: "
            f"ord={actual_orders} "
            f"growth_ratios={growth_ratios} "
            f"first_growth_level="
            f"{first_growth_level}"
        )

    # ------------------------------------------------------------------
    # 4. EXCEPTION SUMMARY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. EXCEPTION SUMMARY")
    print("=" * 78)

    print(
        f"  prime_limit={PRIME_LIMIT}"
    )

    print(
        f"  prime_count={len(primes)}"
    )

    print(
        f"  nonexceptional_count="
        f"{len(nonexceptional)}"
    )

    print(
        f"  exceptional_count="
        f"{len(exceptional)}"
    )

    print(
        f"  exceptional_primes="
        f"{[e[0] for e in exceptional]}"
    )

    print(
        f"  failures={len(failures)}"
    )

    # ------------------------------------------------------------------
    # 5. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The local principal-coordinate theorem depends on the growth law

    ord_{ell^r}(3)
      = ord_ell(3)*ell^(r-1).

The exact obstruction is

    D_ell
      = (3^ord_ell(3)-1)/ell mod ell.

If D_ell != 0, the order grows immediately and the exponent coordinate
has the ordinary base-ell form.

If D_ell = 0, the prime is exceptional. The first lift may fail to
split into ell distinct exponent classes, exactly as observed for
ell=11.

Experiment 247 determines the exceptional set in a larger finite
range independently of the source ratio.

That separates the universal local arithmetic of 3 from the particular
projective ratio q1/q3.
"""
    )

    # ------------------------------------------------------------------
    # 6. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        len(failures) == 0
    )

    print()
    print("=" * 78)
    print("6. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  nonexceptional_order_law_exact="
        f"{final_ok}"
    )

    print(
        f"  exceptional_set="
        f"{[e[0] for e in exceptional]}"
    )

    print(
        f"  failures={len(failures)}"
    )

    print(
        "  source_ratio_independent=True"
    )

    print(
        "  general_n_pq_theorem_established=False"
    )

    print(
        "  reason=single_known_n6_instance"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 247 COMPLETE")


if __name__ == "__main__":

    try:
        main()

    except KeyboardInterrupt:

        print(
            "\nInterrupted."
        )

        sys.exit(130)

    except Exception as exc:

        print(
            "\nFATAL ERROR: "
            + type(exc).__name__
            + ": "
            + str(exc)
        )

        raise

