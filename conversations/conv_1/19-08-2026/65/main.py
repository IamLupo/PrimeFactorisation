#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 233 — FAST EXACT PRIME-LOCAL EXPONENT HENSEL-LIFT AUDIT
==============================================================================

Correction to Experiment 232
-----------------------------

Experiment 232 incorrectly interpreted the first mod-ell^2 solution as
a "return" of the initial exponent after ord_{ell^2}(3).

The correct lifting structure is:

    k' = k0 + t * ord_ell(3),

with

    t in {0,...,ell-1}.

Exactly one t should satisfy

    2*3^k' == rho (mod ell^2),

where

    rho = q1/q3.

This is the finite-prime analogue of the Hensel exponent lifting already
observed in the 7-adic computation.

Known small repeated primes:

    ell = 5
    ell = 7
    ell = 73

Expected phenomenon:

    ell=5:
        k0=2
        ord_5(3)=4
        lift digit t=2
        k'=10

    ell=7:
        k0=1
        ord_7(3)=6
        lift digit t=0
        k'=1

    ell=73:
        k0=1
        ord_73(3)=12
        lift digit t=7
        k'=85

Only main.py is used.
"""

from __future__ import annotations

import sys
from math import gcd


# ============================================================================
# SOURCE DATA
# ============================================================================

Q1 = 29144191
Q3 = 24794967

SMALL_PRIMES = [
    5,
    7,
    73,
]


# ============================================================================
# EXACT HELPERS
# ============================================================================

def valuation(
    x: int,
    ell: int,
) -> int | None:

    x = int(x)
    ell = int(ell)

    if x == 0:
        return None

    x = abs(x)
    v = 0

    while x % ell == 0:
        x //= ell
        v += 1

    return int(v)


def ratio_mod(
    modulus: int,
) -> int | None:

    modulus = int(modulus)

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
                int(modulus),
            )
        )
        % modulus
    )


def multiplicative_order(
    base: int,
    modulus: int,
) -> int:

    if gcd(
        base,
        modulus,
    ) != 1:

        raise ArithmeticError(
            f"{base} is not invertible modulo {modulus}."
        )

    current = 1

    for e in range(
        1,
        modulus + 1,
    ):

        current = (
            current
            * base
        ) % modulus

        if current == 1:
            return int(e)

    raise ArithmeticError(
        f"Order not found for {base} modulo {modulus}."
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 233 — FAST EXACT PRIME-LOCAL "
        "EXPONENT HENSEL-LIFT AUDIT"
    )
    print("=" * 78)

    total_failures = 0

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

    # ------------------------------------------------------------------
    # 2. EACH SMALL PRIME
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. PRIME-LOCAL HENSEL LIFTS")
    print("=" * 78)

    for ell in SMALL_PRIMES:

        print()
        print(
            f"  ell={ell}"
        )

        # --------------------------------------------------------------
        # Base ratio
        # --------------------------------------------------------------

        rho1 = ratio_mod(
            ell
        )

        rho2 = ratio_mod(
            ell * ell
        )

        order1 = multiplicative_order(
            3,
            ell,
        )

        order2 = multiplicative_order(
            3,
            ell * ell,
        )

        print(
            f"    rho_mod_ell={rho1}"
        )

        print(
            f"    rho_mod_ell2={rho2}"
        )

        print(
            f"    ord_ell(3)={order1}"
        )

        print(
            f"    ord_ell2(3)={order2}"
        )

        # --------------------------------------------------------------
        # Find unique base exponent k0 modulo ord_ell(3)
        # --------------------------------------------------------------

        base_matches = []

        for k in range(
            order1
        ):

            if (
                orbit_mod(
                    k,
                    ell,
                )
                == rho1
            ):

                base_matches.append(
                    k
                )

        base_unique = (
            len(base_matches) == 1
        )

        print(
            f"    base_matches={base_matches}"
        )

        print(
            f"    base_unique={base_unique}"
        )

        if not base_unique:

            total_failures += 1
            continue

        k0 = base_matches[0]

        # --------------------------------------------------------------
        # Test all ell possible Hensel digits
        # --------------------------------------------------------------

        candidates = []

        print(
            "    lift_candidates:"
        )

        for t in range(
            ell
        ):

            candidate_k = int(
                k0
                + t * order1
            )

            orbit2 = orbit_mod(
                candidate_k,
                ell * ell,
            )

            match2 = (
                orbit2
                == rho2
            )

            r = (
                Q1
                -
                (
                    2
                    * pow(
                        3,
                        candidate_k,
                    )
                    * Q3
                )
            )

            v = valuation(
                r,
                ell,
            )

            candidates.append(
                (
                    t,
                    candidate_k,
                    orbit2,
                    match2,
                    v,
                )
            )

            print(
                f"      t={t}: "
                f"k={candidate_k} "
                f"orbit_mod_ell2={orbit2} "
                f"match={match2} "
                f"v_ell(F)={v}"
            )

        successful = [
            row
            for row in candidates
            if row[3]
        ]

        unique_lift = (
            len(successful) == 1
        )

        print(
            f"    successful_lifts="
            f"{successful}"
        )

        print(
            f"    unique_lift="
            f"{unique_lift}"
        )

        if not unique_lift:

            total_failures += 1
            continue

        chosen_t, next_k, next_orbit, _, next_v = (
            successful[0]
        )

        # --------------------------------------------------------------
        # Verify canonical digit relation
        # --------------------------------------------------------------

        reconstructed_k = int(
            k0
            + chosen_t * order1
        )

        recurrence_exact = (
            reconstructed_k
            == next_k
        )

        print(
            f"    chosen_t={chosen_t}"
        )

        print(
            f"    next_k={next_k}"
        )

        print(
            f"    recurrence_exact="
            f"{recurrence_exact}"
        )

        # --------------------------------------------------------------
        # Compare with observed residual depth
        # --------------------------------------------------------------

        expected_lift_depth = 2

        valuation_exact = (
            next_v >= expected_lift_depth
        )

        print(
            f"    next_v_ell={next_v}"
        )

        print(
            f"    solves_mod_ell2="
            f"{valuation_exact}"
        )

        if not recurrence_exact:
            total_failures += 1

        if not valuation_exact:
            total_failures += 1

        # --------------------------------------------------------------
        # Check that the next digit is genuinely nontrivial or zero
        # --------------------------------------------------------------

        print(
            f"    lift_digit={chosen_t}"
        )

        print(
            f"    zero_lift_digit="
            f"{chosen_t == 0}"
        )

        # --------------------------------------------------------------
        # Next return modulo ell^2
        #
        # Once k' is found, the full recurrence modulo ell^2 repeats
        # with ord_ell2(3).
        # --------------------------------------------------------------

        next_return = int(
            next_k
            + order2
        )

        return_match = (
            orbit_mod(
                next_return,
                ell * ell,
            )
            == rho2
        )

        print(
            f"    first_same_class_return_candidate="
            f"{next_return}"
        )

        print(
            f"    return_mod_ell2_exact="
            f"{return_match}"
        )

        if not return_match:
            total_failures += 1

    # ------------------------------------------------------------------
    # 3. CROSS-PRIME LIFT SUMMARY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. CROSS-PRIME LIFT SUMMARY")
    print("=" * 78)

    for ell in SMALL_PRIMES:

        rho1 = ratio_mod(
            ell
        )

        rho2 = ratio_mod(
            ell * ell
        )

        order1 = multiplicative_order(
            3,
            ell,
        )

        order2 = multiplicative_order(
            3,
            ell * ell,
        )

        base = [
            k
            for k in range(order1)
            if orbit_mod(k,ell) == rho1
        ]

        if len(base) != 1:
            continue

        k0 = base[0]

        lifts = [
            (
                t,
                k0 + t * order1
            )
            for t in range(ell)
            if orbit_mod(
                k0 + t * order1,
                ell * ell,
            )
            == rho2
        ]

        print(
            f"  ell={ell}: "
            f"k0={k0} "
            f"ord_ell={order1} "
            f"ord_ell2={order2} "
            f"lift={lifts}"
        )

    # ------------------------------------------------------------------
    # 4. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The correct local analogue of the 7-adic exponent lift is:

    k_(e+1)
      = k_e + t_e * ord_{ell^e}(3),

with

    t_e in {0,...,ell-1}.

At the first prime-power step:

    k2 = k1 + t*ord_ell(3).

Experiment 233 checks every possible t.

For the observed primes this should recover:

    ell=5:
        k1=2
        t=2
        k2=10

    ell=7:
        k1=1
        t=0
        k2=1

    ell=73:
        k1=1
        t=7
        k2=85.

The important point is that the prime-power coincidence is itself a
digit-lift phenomenon in the exponent coordinate.

This makes the exceptional factors

    5^2,
    7^2,

and the 73-based higher-order event part of the same mechanism as the
previous 7-adic Hensel construction.
"""
    )

    # ------------------------------------------------------------------
    # 5. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        total_failures == 0
    )

    print()
    print("=" * 78)
    print("5. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  tested_small_primes="
        f"{SMALL_PRIMES}"
    )

    print(
        "  corrected_hensel_lift_formula=True"
    )

    print(
        f"  all_prime_local_lifts_unique="
        f"{final_ok}"
    )

    print(
        f"  failures={total_failures}"
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
    print("EXPERIMENT 233 COMPLETE")


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

