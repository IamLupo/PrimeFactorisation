#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 248 — EXACT ODD-PRIME ORDER-GROWTH / EXCEPTIONAL-SET AUDIT
==============================================================================

Correction of Experiment 247:

    ell=2 is excluded because the local principal-unit/Hensel theory
    under study is for odd primes ell.

Experiment 248:

    * classifies all odd primes <= PRIME_LIMIT;
    * computes d = ord_ell(3);
    * computes
          D_ell = (3^d - 1)/ell mod ell;
    * classifies ell as nonexceptional or exceptional;
    * verifies the order tower through ell^4;
    * verifies that every nonexceptional prime satisfies

          ord_{ell^r}(3)
            = ord_ell(3) * ell^(r-1);

    * for every exceptional prime, records the delayed-growth tower.

This experiment is source-ratio independent.

Only main.py is used.
"""

from __future__ import annotations

import sys


PRIME_LIMIT = 20000
MAX_LEVEL = 4


# ============================================================================
# PRIME TEST
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


# ============================================================================
# FACTORIZATION OF SMALL PHI VALUES
# ============================================================================

def prime_factors(n: int) -> list[int]:

    factors = []

    if n % 2 == 0:

        factors.append(2)

        while n % 2 == 0:
            n //= 2

    d = 3

    while d * d <= n:

        if n % d == 0:

            factors.append(d)

            while n % d == 0:
                n //= d

        d += 2

    if n > 1:
        factors.append(n)

    return factors


# ============================================================================
# EXACT MULTIPLICATIVE ORDER
# ============================================================================

def multiplicative_order(
    base: int,
    modulus: int,
    phi: int,
) -> int:

    order = int(phi)

    for p in prime_factors(order):

        while (
            order % p == 0
            and
            pow(
                base,
                order // p,
                modulus,
            ) == 1
        ):

            order //= p

    return int(order)


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 248 — EXACT ODD-PRIME ORDER-GROWTH / "
        "EXCEPTIONAL-SET AUDIT"
    )
    print("=" * 78)

    failures = []

    nonexceptional = []
    exceptional = []

    primes = [
        p
        for p in range(
            3,
            PRIME_LIMIT + 1,
            2,
        )
        if is_prime(p)
    ]

    # ------------------------------------------------------------------
    # 1. CLASSIFICATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. ODD-PRIME CLASSIFICATION")
    print("=" * 78)

    for ell in primes:

        phi1 = ell - 1

        d = multiplicative_order(
            3,
            ell,
            phi1,
        )

        modulus2 = ell * ell

        delta = (
            pow(
                3,
                d,
                modulus2,
            )
            - 1
        )

        if delta % ell != 0:

            failures.append(
                (
                    ell,
                    "D_not_integral",
                )
            )

            continue

        D = int(
            (
                delta // ell
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

    # ------------------------------------------------------------------
    # 2. NONEXCEPTIONAL ORDER LAW
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. NONEXCEPTIONAL ORDER-LAW VERIFICATION")
    print("=" * 78)

    nonexceptional_failure_count = 0

    for ell, d, D in nonexceptional:

        actual = []
        expected = []

        for level in range(
            1,
            MAX_LEVEL + 1,
        ):

            modulus = ell ** level

            phi = (
                (ell - 1)
                * ell ** (
                    level - 1
                )
            )

            order = multiplicative_order(
                3,
                modulus,
                phi,
            )

            actual.append(order)

            expected.append(
                d
                * ell ** (
                    level - 1
                )
            )

        ok = (
            actual == expected
        )

        if not ok:

            nonexceptional_failure_count += 1

            failures.append(
                (
                    ell,
                    "nonexceptional_order_law",
                    actual,
                    expected,
                )
            )

            print(
                f"  ell={ell}: "
                f"actual={actual} "
                f"expected={expected}"
            )

    print(
        f"  tested={len(nonexceptional)}"
    )

    print(
        f"  failures={nonexceptional_failure_count}"
    )

    # ------------------------------------------------------------------
    # 3. EXCEPTIONAL TOWERS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. EXCEPTIONAL PRIME TOWERS")
    print("=" * 78)

    for ell, d, D in exceptional:

        orders = []

        growth = []

        for level in range(
            1,
            MAX_LEVEL + 1,
        ):

            modulus = ell ** level

            phi = (
                (ell - 1)
                * ell ** (
                    level - 1
                )
            )

            order = multiplicative_order(
                3,
                modulus,
                phi,
            )

            orders.append(order)

            if level == 1:

                growth.append(None)

            else:

                previous = orders[
                    level - 2
                ]

                if (
                    order
                    % previous
                    != 0
                ):

                    failures.append(
                        (
                            ell,
                            "exceptional_order_not_dividing_previous",
                        )
                    )

                    growth.append(
                        None
                    )

                else:

                    growth.append(
                        order // previous
                    )

        first_growth = None

        for level in range(
            2,
            MAX_LEVEL + 1,
        ):

            if (
                orders[level - 1]
                != orders[level - 2]
            ):

                first_growth = level
                break

        print(
            f"  ell={ell}: "
            f"ord={orders} "
            f"growth={growth} "
            f"first_growth_level="
            f"{first_growth}"
        )

    # ------------------------------------------------------------------
    # 4. SUMMARY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. SUMMARY")
    print("=" * 78)

    print(
        f"  prime_limit={PRIME_LIMIT}"
    )

    print(
        f"  odd_prime_count={len(primes)}"
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
        f"{[x[0] for x in exceptional]}"
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
Experiment 247 was contaminated by ell=2, which belongs to a different
local setting.

Experiment 248 restricts the analysis to odd primes.

For an odd prime ell let

    d = ord_ell(3),

and

    D_ell
      = (3^d - 1)/ell mod ell.

The generic branch is

    D_ell != 0.

Then the order grows maximally:

    ord_{ell^r}(3)
      = d*ell^(r-1).

The exceptional branch is

    D_ell = 0,

where the first order lift is delayed.

This distinction is exactly what controls whether the ordinary
base-ell principal-coordinate expansion is available.

The experiment therefore isolates the purely local arithmetic
exception set before returning to the projective source ratio.
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
        f"  odd_prime_order_law_exact="
        f"{final_ok}"
    )

    print(
        f"  exceptional_primes="
        f"{[x[0] for x in exceptional]}"
    )

    print(
        f"  excluded_ell_2=True"
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
    print("EXPERIMENT 248 COMPLETE")


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

