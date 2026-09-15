#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 244 — EXACT PRIME-LOCAL ORDER-GROWTH / NON-WIEFERICH AUDIT
==============================================================================

Experiment 243 established, for all 28 admissible primes <= 200,

    ord_{ell^r}(3) = ord_ell(3) * ell^(r-1)

through r=4, and consequently the Hensel digits became the ordinary
base-ell digits of the principal exponent coordinate.

Experiment 244 isolates the arithmetic condition responsible for that
order growth.

Let

    d = ord_ell(3).

Define

    D1 = (3^d - 1)/ell mod ell.

If

    D1 != 0 mod ell,

then the order grows by ell at every subsequent level:

    ord_{ell^r}(3)
      = d * ell^(r-1).

If

    D1 = 0 mod ell,

then ell is exceptional for the lifting of 3 and the simple principal
coordinate scaling requires special treatment.

This experiment therefore checks, for all primes <= 200:

    1. admissible/inadmissible projective orbit chart;
    2. D1;
    3. order growth through ell^4;
    4. whether the simple formula holds;
    5. whether any exceptional prime appears.

It does NOT assume that every admissible prime is non-exceptional.

Only main.py is used.
"""

from __future__ import annotations

import sys
from math import gcd


Q1 = 29144191
Q3 = 24794967

PRIME_LIMIT = 200


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


def multiplicative_order_prime(
    base: int,
    prime: int,
) -> int:

    current = 1

    for e in range(
        1,
        prime + 1,
    ):
        current = (
            current * base
        ) % prime

        if current == 1:
            return int(e)

    raise ArithmeticError(
        f"ord_{prime}({base}) not found."
    )


def lift_order_once(
    base: int,
    prime: int,
    previous_order: int,
    current_power: int,
) -> int:

    next_modulus = prime ** (
        current_power + 1
    )

    if (
        pow(
            base,
            previous_order,
            next_modulus,
        )
        == 1
    ):
        return int(previous_order)

    return int(
        previous_order * prime
    )


def ratio_mod(
    modulus: int,
) -> int | None:

    if gcd(
        Q3,
        modulus,
    ) != 1:
        return None

    return int(
        (
            (Q1 % modulus)
            * pow(
                Q3 % modulus,
                -1,
                modulus,
            )
        )
        % modulus
    )


def orbit_mod(
    k: int,
    modulus: int,
) -> int:

    return int(
        (
            2
            * pow(
                3,
                int(k),
                modulus,
            )
        )
        % modulus
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 244 — EXACT PRIME-LOCAL ORDER-GROWTH / "
        "NON-WIEFERICH AUDIT"
    )
    print("=" * 78)

    failures = []

    all_primes = []
    admissible = []
    exceptional = []
    nonexceptional = []

    # ------------------------------------------------------------------
    # 1. SOURCE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. SOURCE DATA")
    print("=" * 78)

    print(
        f"  q1={Q1}"
    )

    print(
        f"  q3={Q3}"
    )

    print(
        f"  gcd(q1,q3)={gcd(Q1,Q3)}"
    )

    print(
        f"  prime_limit={PRIME_LIMIT}"
    )

    # ------------------------------------------------------------------
    # 2. PRIME CLASSIFICATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. PRIME ORDER-GROWTH CLASSIFICATION")
    print("=" * 78)

    for ell in range(
        2,
        PRIME_LIMIT + 1,
    ):

        if not is_prime(ell):
            continue

        all_primes.append(ell)

        if ell in (2, 3):
            continue

        if Q3 % ell == 0:
            continue

        d = multiplicative_order_prime(
            3,
            ell,
        )

        modulus2 = ell * ell

        pow3d_mod_ell2 = pow(
            3,
            d,
            modulus2,
        )

        delta = (
            pow3d_mod_ell2 - 1
        )

        if delta % ell != 0:

            failures.append(
                (
                    ell,
                    "3^d-1_not_divisible_by_ell",
                )
            )

            continue

        D1 = int(
            (
                delta
                // ell
            )
            % ell
        )

        is_exceptional = (
            D1 == 0
        )

        if is_exceptional:
            exceptional.append(
                (
                    ell,
                    d,
                    D1,
                )
            )
        else:
            nonexceptional.append(
                (
                    ell,
                    d,
                    D1,
                )
            )

        rho = ratio_mod(
            ell
        )

        base_matches = [
            k
            for k in range(d)
            if orbit_mod(
                k,
                ell,
            ) == rho
        ]

        is_admissible = (
            len(base_matches) == 1
        )

        if is_admissible:
            admissible.append(
                ell
            )

        # --------------------------------------------------------------
        # Actual order growth.
        # --------------------------------------------------------------

        order1 = d

        order2 = lift_order_once(
            3,
            ell,
            order1,
            1,
        )

        order3 = lift_order_once(
            3,
            ell,
            order2,
            2,
        )

        order4 = lift_order_once(
            3,
            ell,
            order3,
            3,
        )

        expected_simple = [
            d,
            d * ell,
            d * ell * ell,
            d * ell * ell * ell,
        ]

        actual = [
            order1,
            order2,
            order3,
            order4,
        ]

        simple_growth = (
            actual == expected_simple
        )

        print(
            f"  ell={ell}: "
            f"ord1={order1} "
            f"D1={D1} "
            f"exceptional={is_exceptional} "
            f"admissible={is_admissible} "
            f"orders={actual} "
            f"simple_growth={simple_growth}"
        )

        # For the sample up to 200, an exceptional prime should force
        # failure of the simple order-growth identity.
        expected_exceptional_relation = (
            is_exceptional
            == (not simple_growth)
        )

        if not expected_exceptional_relation:
            failures.append(
                (
                    ell,
                    "exceptional_vs_growth_mismatch",
                    is_exceptional,
                    actual,
                    expected_simple,
                )
            )

    # ------------------------------------------------------------------
    # 3. SUMMARY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. SUMMARY")
    print("=" * 78)

    print(
        f"  primes_seen={len(all_primes)}"
    )

    print(
        f"  admissible_count={len(admissible)}"
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
        f"{exceptional}"
    )

    # ------------------------------------------------------------------
    # 4. ADMISSIBLE ORDER PROFILE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. ADMISSIBLE PRIME ORDER PROFILE")
    print("=" * 78)

    admissible_exceptions = []

    for ell in admissible:

        d = multiplicative_order_prime(
            3,
            ell,
        )

        D1 = int(
            (
                (
                    pow(
                        3,
                        d,
                        ell * ell,
                    )
                    - 1
                )
                // ell
            )
            % ell
        )

        order2 = lift_order_once(
            3,
            ell,
            d,
            1,
        )

        order3 = lift_order_once(
            3,
            ell,
            order2,
            2,
        )

        order4 = lift_order_once(
            3,
            ell,
            order3,
            3,
        )

        expected = [
            d,
            d * ell,
            d * ell ** 2,
            d * ell ** 3,
        ]

        actual = [
            d,
            order2,
            order3,
            order4,
        ]

        ok = (
            actual == expected
            and
            D1 != 0
        )

        if not ok:
            admissible_exceptions.append(
                (
                    ell,
                    D1,
                    actual,
                    expected,
                )
            )

        print(
            f"  ell={ell}: "
            f"D1={D1} "
            f"orders={actual} "
            f"expected={expected} "
            f"OK={ok}"
        )

    # ------------------------------------------------------------------
    # 5. CONSEQUENCE FOR PRINCIPAL COORDINATES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. PRINCIPAL-COORDINATE CONSEQUENCE")
    print("=" * 78)

    for ell in admissible:

        d = multiplicative_order_prime(
            3,
            ell,
        )

        D1 = int(
            (
                (
                    pow(
                        3,
                        d,
                        ell * ell,
                    )
                    - 1
                )
                // ell
            )
            % ell
        )

        simple = (
            D1 != 0
        )

        print(
            f"  ell={ell}: "
            f"canonical_base_ell_coordinate="
            f"{simple}"
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
Experiment 243 showed that, for all 28 admissible primes in the sample,
the Hensel digits are literally the base-ell digits of the principal
coordinate.

Experiment 244 identifies the arithmetic condition behind that fact.

Let

    d = ord_ell(3).

If

    D1 = (3^d - 1)/ell mod ell

is nonzero, then the order grows normally:

    ord_{ell^r}(3)
      = d*ell^(r-1).

Consequently the exponent recurrence

    k_{r+1}
      = k_r + t_r*ord_{ell^r}(3)

becomes

    m_{r+1}
      = m_r + t_r*ell^(r-1),

so the Hensel digits are the ordinary base-ell digits of m.

If D1 vanishes, the simple coordinate scaling breaks and a separate
exceptional analysis is required.

Thus the principal-coordinate theorem has a natural local hypothesis:
nonexceptional order growth of 3 at ell.
"""
    )

    # ------------------------------------------------------------------
    # 7. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        len(failures) == 0
        and
        len(admissible_exceptions) == 0
    )

    print()
    print("=" * 78)
    print("7. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  admissible_count="
        f"{len(admissible)}"
    )

    print(
        f"  exceptional_count="
        f"{len(exceptional)}"
    )

    print(
        f"  admissible_exceptional_conflicts="
        f"{admissible_exceptions}"
    )

    print(
        f"  order_growth_hypothesis_exact="
        f"{final_ok}"
    )

    print(
        f"  failures="
        f"{len(failures) + len(admissible_exceptions)}"
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
    print("EXPERIMENT 244 COMPLETE")


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

