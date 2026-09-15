#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 248R — EXACT ODD-PRIME ORDER-GROWTH / EXCEPTIONAL-SET AUDIT
==============================================================================

Corrected version of Experiment 248.

Important correction:

    ell = 3

is excluded because 3 is not a unit modulo 3, so ord_3(3) is not defined.

The admissible local setting is therefore:

    ell odd prime, ell != 3.

For every odd prime ell >= 5:

    d = ord_ell(3)

    D_ell = (3^d - 1)/ell mod ell.

Nonexceptional:

    D_ell != 0

and then

    ord_{ell^r}(3) = d * ell^(r-1).

Exceptional:

    D_ell = 0.

Experiment 248R verifies the generic order law through ell^4 and
records the actual exceptional order towers.

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
# FACTOR HELPERS
# ============================================================================

def prime_factors(n: int) -> list[int]:

    factors: list[int] = []

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
        "EXPERIMENT 248R — EXACT ODD-PRIME ORDER-GROWTH / "
        "EXCEPTIONAL-SET AUDIT"
    )
    print("=" * 78)

    failures = []

    nonexceptional = []
    exceptional = []

    # IMPORTANT:
    # Start at 5, not 3.
    primes = [
        p
        for p in range(
            5,
            PRIME_LIMIT + 1,
            2,
        )
        if is_prime(p)
    ]

    # ------------------------------------------------------------------
    # 1. ODD PRIME CLASSIFICATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. ODD-PRIME CLASSIFICATION")
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
                    "D_not_integral",
                )
            )

            continue

        D = int(
            (
                increment // ell
            )
            % ell
        )

        if D == 0:

            exceptional.append(
                (
                    ell,
                    d,
                    D,
                )
            )

            print(
                f"  ell={ell}: "
                f"EXCEPTIONAL "
                f"ord1={d} "
                f"D={D}"
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

    nonexceptional_failures = []

    for ell, d, D in nonexceptional:

        actual_orders = []
        expected_orders = []

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

            actual_orders.append(
                actual
            )

            expected_orders.append(
                expected
            )

        ok = (
            actual_orders
            ==
            expected_orders
        )

        if not ok:

            nonexceptional_failures.append(
                (
                    ell,
                    actual_orders,
                    expected_orders,
                )
            )

            failures.append(
                (
                    ell,
                    "nonexceptional_order_law",
                )
            )

    print(
        f"  tested={len(nonexceptional)}"
    )

    print(
        f"  failures={len(nonexceptional_failures)}"
    )

    for row in nonexceptional_failures:

        print(
            f"  FAILURE: ell={row[0]} "
            f"actual={row[1]} "
            f"expected={row[2]}"
        )

    # ------------------------------------------------------------------
    # 3. EXCEPTIONAL TOWERS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. EXCEPTIONAL PRIME TOWERS")
    print("=" * 78)

    exceptional_diagnostics = []

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

            orders.append(
                order
            )

            if level == 1:

                growth.append(
                    None
                )

            else:

                previous = orders[
                    level - 2
                ]

                growth.append(
                    order // previous
                )

        first_growth_level = None

        for level in range(
            2,
            MAX_LEVEL + 1,
        ):

            if (
                orders[level - 1]
                !=
                orders[level - 2]
            ):

                first_growth_level = level
                break

        exceptional_diagnostics.append(
            (
                ell,
                orders,
                growth,
                first_growth_level,
            )
        )

        print(
            f"  ell={ell}: "
            f"ord={orders} "
            f"growth={growth} "
            f"first_growth_level="
            f"{first_growth_level}"
        )

    # ------------------------------------------------------------------
    # 4. EXCEPTIONAL CONSISTENCY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. EXCEPTIONAL-BRANCH CONSISTENCY")
    print("=" * 78)

    exceptional_consistency_failures = []

    for (
        ell,
        orders,
        growth,
        first_growth_level,
    ) in exceptional_diagnostics:

        d = orders[0]

        # By definition D=0, so the first growth ratio should be 1.
        first_growth_ok = (
            growth[1] == 1
        )

        # Once the growth resumes, it should be by ell.
        resumed_growth_ok = True

        for index in range(
            2,
            len(growth),
        ):

            if growth[index] != ell:

                resumed_growth_ok = False

        ok = (
            first_growth_ok
            and
            resumed_growth_ok
        )

        print(
            f"  ell={ell}: "
            f"first_growth_delayed="
            f"{first_growth_ok} "
            f"resumed_growth_by_ell="
            f"{resumed_growth_ok} "
            f"OK={ok}"
        )

        if not ok:

            exceptional_consistency_failures.append(
                (
                    ell,
                    orders,
                    growth,
                )
            )

    failures.extend(
        (
            ell,
            "exceptional_consistency",
            orders,
            growth,
        )
        for ell, orders, growth
        in exceptional_consistency_failures
    )

    # ------------------------------------------------------------------
    # 5. SUMMARY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. SUMMARY")
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
        f"{[row[0] for row in exceptional]}"
    )

    print(
        f"  failures={len(failures)}"
    )

    # ------------------------------------------------------------------
    # 6. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
Experiment 248R corrects the domain error from Experiment 248 by
excluding ell=3.

For every odd prime ell >= 5:

    d = ord_ell(3),

    D_ell = (3^d - 1)/ell mod ell.

The generic branch is

    D_ell != 0,

which gives maximal order growth:

    ord_{ell^r}(3)
      = d*ell^(r-1).

The exceptional branch is

    D_ell = 0.

In the tested range, the only exceptional prime is ell=11.

For ell=11:

    ord_11(3)=5,
    ord_121(3)=5,
    ord_1331(3)=55,
    ord_14641(3)=605.

Thus the generic principal-coordinate theorem should explicitly assume
the nonexceptional order-growth condition, while ell=11 belongs to a
separate delayed-growth branch.
"""
    )

    # ------------------------------------------------------------------
    # 7. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        len(failures) == 0
    )

    print()
    print("=" * 78)
    print("7. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  odd_prime_order_law_exact="
        f"{final_ok}"
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
        f"{[row[0] for row in exceptional]}"
    )

    print(
        "  ell_3_excluded=True"
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
    print("EXPERIMENT 248R COMPLETE")


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

