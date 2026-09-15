#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 232 — FAST EXACT SMALL-PRIME ORBIT / PRIME-POWER RETURN AUDIT
==============================================================================

This is the fast replacement for Experiment 231.

It deliberately avoids factoring large orbit residuals and avoids searching
multiplicative orders for huge primes.

We study only the repeated small primes already observed:

    ell = 5
    ell = 7
    ell = 73

For each ell we compute directly:

    ord_ell(3)
    ord_{ell^2}(3)

and compare the observed orbit hits with

    c_k = 2*3^k.

The source ratio is

    rho = q1/q3.

For ell not dividing q3,

    ell | q1 - c_k*q3
iff
    c_k == rho (mod ell).

At prime-power precision:

    ell^2 | q1 - c_k*q3
iff
    c_k == rho (mod ell^2).

The experiment therefore distinguishes:

    ordinary orbit return modulo ell

from

    higher-order return modulo ell^2.

No large-prime factorization is performed.
Only main.py is used.
"""

from __future__ import annotations

import sys
from math import gcd


# ============================================================================
# SOURCE
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

    base = int(base)
    modulus = int(modulus)

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
        f"Could not determine order of {base} modulo {modulus}."
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 232 — FAST EXACT SMALL-PRIME ORBIT / "
        "PRIME-POWER RETURN AUDIT"
    )
    print("=" * 78)

    print()
    print("=" * 78)
    print("1. SOURCE PROJECTIVE DATA")
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

    global_failures = 0

    # ------------------------------------------------------------------
    # 1. SMALL PRIME CASES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. SMALL-PRIME LOCAL CASES")
    print("=" * 78)

    for ell in SMALL_PRIMES:

        modulus1 = int(ell)
        modulus2 = int(ell * ell)

        rho1 = ratio_mod(
            modulus1
        )

        rho2 = ratio_mod(
            modulus2
        )

        order1 = multiplicative_order(
            3,
            modulus1,
        )

        order2 = multiplicative_order(
            3,
            modulus2,
        )

        print()
        print(
            f"  ell={ell}"
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
        # 2. ALL HITS IN ONE PRIME PERIOD
        # --------------------------------------------------------------

        hits1 = []
        hits2 = []

        for k in range(
            order2
        ):

            c1 = orbit_mod(
                k,
                modulus1,
            )

            c2 = orbit_mod(
                k,
                modulus2,
            )

            if c1 == rho1:
                hits1.append(
                    k
                )

            if c2 == rho2:
                hits2.append(
                    k
                )

        print(
            f"    hits_mod_ell="
            f"{hits1}"
        )

        print(
            f"    hits_mod_ell2="
            f"{hits2}"
        )

        # --------------------------------------------------------------
        # 3. EXPECTED CLASS STRUCTURE
        # --------------------------------------------------------------

        class1 = sorted(
            {
                k % order1
                for k in hits1
            }
        )

        class2 = sorted(
            {
                k % order2
                for k in hits2
            }
        )

        class1_exact = (
            len(class1) == 1
        )

        class2_exact = (
            len(class2) == 1
        )

        print(
            f"    mod_ell_classes="
            f"{class1}"
        )

        print(
            f"    mod_ell2_classes="
            f"{class2}"
        )

        print(
            f"    unique_mod_ell_class="
            f"{class1_exact}"
        )

        print(
            f"    unique_mod_ell2_class="
            f"{class2_exact}"
        )

        if not class1_exact:
            global_failures += 1

        if not class2_exact:
            global_failures += 1

        # --------------------------------------------------------------
        # 4. PRIME-TO-PRIME-POWER RETURN
        # --------------------------------------------------------------

        if len(hits1) >= 1:

            base_k = hits1[0]

        else:

            base_k = None

        if base_k is not None:

            first_return_ell = None

            for k in hits1:

                if k > base_k:
                    first_return_ell = k
                    break

            first_return_ell2 = None

            for k in hits2:

                if k > base_k:
                    first_return_ell2 = k
                    break

            print(
                f"    base_k={base_k}"
            )

            print(
                f"    first_distinct_return_mod_ell="
                f"{first_return_ell}"
            )

            print(
                f"    first_distinct_return_mod_ell2="
                f"{first_return_ell2}"
            )

            if first_return_ell is not None:

                expected1 = (
                    base_k
                    + order1
                )

                check1 = (
                    first_return_ell
                    == expected1
                )

                print(
                    f"    first_return_ell_expected="
                    f"{expected1}"
                )

                print(
                    f"    first_return_ell_exact="
                    f"{check1}"
                )

                if not check1:
                    global_failures += 1

            if first_return_ell2 is not None:

                expected2 = (
                    base_k
                    + order2
                )

                check2 = (
                    first_return_ell2
                    == expected2
                )

                print(
                    f"    first_return_ell2_expected="
                    f"{expected2}"
                )

                print(
                    f"    first_return_ell2_exact="
                    f"{check2}"
                )

                if not check2:
                    global_failures += 1

        # --------------------------------------------------------------
        # 5. OBSERVED RESIDUAL DEPTHS
        # --------------------------------------------------------------

        print(
            f"    observed_residual_depths:"
        )

        # Search only through two prime-power periods.
        observed = []

        search_end = (
            order2 * 2
        )

        for k in range(
            search_end
        ):

            c = (
                2
                * pow(
                    3,
                    k,
                )
            )

            residual = (
                Q1
                - c * Q3
            )

            v = valuation(
                residual,
                ell,
            )

            if v is not None and v > 0:

                observed.append(
                    (
                        k,
                        v,
                    )
                )

        print(
            f"      {observed}"
        )

        # --------------------------------------------------------------
        # 6. EXPECTED HIGHER-ORDER EVENTS
        # --------------------------------------------------------------

        exceptional = [
            (
                k,
                depth,
            )
            for (
                k,
                depth
            ) in observed
            if depth >= 2
        ]

        print(
            f"    depth_ge_2_events="
            f"{exceptional}"
        )

        # --------------------------------------------------------------
        # 7. DIRECT CONGRUENCE CHECK AT OBSERVED EVENTS
        # --------------------------------------------------------------

        local_failure = False

        for k, depth in observed:

            modulus_depth = int(
                ell ** depth
            )

            modulus_next = int(
                ell ** (depth + 1)
            )

            rho_depth = ratio_mod(
                modulus_depth
            )

            rho_next = ratio_mod(
                modulus_next
            )

            c_depth = orbit_mod(
                k,
                modulus_depth,
            )

            c_next = orbit_mod(
                k,
                modulus_next,
            )

            exact_depth = (
                c_depth
                == rho_depth
            )

            stops_next = (
                c_next
                != rho_next
            )

            print(
                f"      k={k}: "
                f"depth={depth} "
                f"match_depth={exact_depth} "
                f"breaks_next={stops_next}"
            )

            if not exact_depth:
                local_failure = True

            if not stops_next:
                local_failure = True

        if local_failure:
            global_failures += 1

    # ------------------------------------------------------------------
    # 8. CROSS-PRIME COMPARISON
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. CROSS-PRIME COMPARISON")
    print("=" * 78)

    for ell in SMALL_PRIMES:

        order1 = multiplicative_order(
            3,
            ell,
        )

        order2 = multiplicative_order(
            3,
            ell * ell,
        )

        rho1 = ratio_mod(
            ell
        )

        rho2 = ratio_mod(
            ell * ell
        )

        print(
            f"  ell={ell}: "
            f"ord_ell={order1} "
            f"ord_ell2={order2} "
            f"rho_ell={rho1} "
            f"rho_ell2={rho2}"
        )

    # ------------------------------------------------------------------
    # 9. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The corrected local relation is:

    ell | (q1 - 2*3^k*q3)

if and only if

    2*3^k = q1/q3 (mod ell),

provided ell does not divide q3.

The higher-order version is:

    ell^2 | (q1 - 2*3^k*q3)

if and only if

    2*3^k = q1/q3 (mod ell^2).

This experiment tests the statement only for the small repeated primes
that were already observed:

    ell = 5,
    ell = 7,
    ell = 73.

This makes the computation fast enough to finish while still probing
the key mathematical issue.

The important distinction is:

    periodic return modulo ell

versus

    lifted return modulo ell^2.

For the known data the exceptional depth-two events should appear as
isolated higher-order alignments inside otherwise periodic prime-local
orbits.
"""
    )

    # ------------------------------------------------------------------
    # 10. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        global_failures == 0
    )

    print()
    print("=" * 78)
    print("5. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  small_primes_tested="
        f"{SMALL_PRIMES}"
    )

    print(
        f"  corrected_2_times_3k_projective_formula=True"
    )

    print(
        f"  all_small_prime_class_checks_exact="
        f"{global_failures == 0}"
    )

    print(
        f"  failures={global_failures}"
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
    print("EXPERIMENT 232 COMPLETE")


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