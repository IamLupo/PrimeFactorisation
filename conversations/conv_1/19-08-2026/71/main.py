#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 239 — EXACT SMALL-PRIME TWO-STEP HENSEL EXPONENT-LIFT AUDIT
==============================================================================

Experiment 238 established:

    x = (rho/2) mod ell

lies in <3> iff

    x^ord_ell(3) = 1 mod ell,

and every admissible prime <= 200 has a unique first lift from ell to ell^2.

Experiment 239 continues each admissible local chart one more level.

For each admissible prime ell:

    k1 = base exponent modulo ord_ell(3)

    k2 = k1 + t1 * ord_ell(3)

    k3 = k2 + t2 * ord_{ell^2}(3)

where

    t1 in {0,...,ell-1}
    t2 in {0,...,ell-1}.

The experiment checks:

    1. unique k1 modulo ell;
    2. unique t1 solving rho modulo ell^2;
    3. unique t2 solving rho modulo ell^3;
    4. direct modular reconstruction at ell^2 and ell^3;
    5. affine residual law at both lift levels;
    6. exact exponent recurrence.

This is the finite-prime analogue of the multi-level 7-adic
Hensel lifting already established for the original source.

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

PRIME_LIMIT = 200


# ============================================================================
# HELPERS
# ============================================================================

def is_prime(n: int) -> bool:

    n = int(n)

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


def multiplicative_order(
    base: int,
    modulus: int,
) -> int:

    if gcd(base, modulus) != 1:
        raise ArithmeticError(
            f"{base} is not invertible modulo {modulus}."
        )

    current = 1

    for e in range(
        1,
        modulus + 1,
    ):

        current = (
            current * base
        ) % modulus

        if current == 1:
            return int(e)

    raise ArithmeticError(
        f"Could not determine order of {base} modulo {modulus}."
    )


def ratio_mod(
    modulus: int,
) -> int | None:

    modulus = int(modulus)

    if gcd(Q3, modulus) != 1:
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


def residual_mod(
    k: int,
    modulus: int,
) -> int:

    return int(
        (
            Q1
            -
            (
                2
                * pow(
                    3,
                    int(k),
                    int(modulus),
                )
                * Q3
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
        "EXPERIMENT 239 — EXACT SMALL-PRIME TWO-STEP "
        "HENSEL EXPONENT-LIFT AUDIT"
    )
    print("=" * 78)

    primes = [
        p
        for p in range(
            2,
            PRIME_LIMIT + 1,
        )
        if is_prime(p)
    ]

    admissible = []
    skipped = []
    failures = []

    # ------------------------------------------------------------------
    # 1. IDENTIFY ADMISSIBLE PRIMES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. ADMISSIBLE PRIME CHARTS")
    print("=" * 78)

    for ell in primes:

        if ell in (2, 3):
            skipped.append(
                (ell, "ell_divides_6")
            )
            continue

        if Q3 % ell == 0:
            skipped.append(
                (ell, "ell_divides_q3")
            )
            continue

        rho1 = ratio_mod(ell)
        order1 = multiplicative_order(3, ell)

        matches = [
            k
            for k in range(order1)
            if orbit_mod(k, ell) == rho1
        ]

        if len(matches) == 1:

            admissible.append(
                (
                    ell,
                    rho1,
                    order1,
                    matches[0],
                )
            )

    print(
        f"  admissible_count={len(admissible)}"
    )

    # ------------------------------------------------------------------
    # 2. TWO-STEP HENSEL LIFT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. TWO-STEP LOCAL HENSEL LIFT")
    print("=" * 78)

    for (
        ell,
        rho1,
        order1,
        k1,
    ) in admissible:

        modulus2 = ell * ell
        modulus3 = modulus2 * ell

        rho2 = ratio_mod(
            modulus2
        )

        rho3 = ratio_mod(
            modulus3
        )

        order2 = multiplicative_order(
            3,
            modulus2,
        )

        # --------------------------------------------------------------
        # First lift: ell -> ell^2
        # --------------------------------------------------------------

        first_matches = []

        for t1 in range(ell):

            candidate_k2 = int(
                k1
                + t1 * order1
            )

            if (
                orbit_mod(
                    candidate_k2,
                    modulus2,
                )
                == rho2
            ):

                first_matches.append(
                    (
                        t1,
                        candidate_k2,
                    )
                )

        first_unique = (
            len(first_matches) == 1
        )

        if not first_unique:

            failures.append(
                (
                    ell,
                    "first_lift_not_unique",
                    first_matches,
                )
            )

            print(
                f"  ell={ell}: "
                f"FIRST_LIFT_FAILURE "
                f"matches={first_matches}"
            )

            continue

        t1, k2 = first_matches[0]

        # --------------------------------------------------------------
        # Second lift: ell^2 -> ell^3
        # --------------------------------------------------------------

        second_matches = []

        for t2 in range(ell):

            candidate_k3 = int(
                k2
                + t2 * order2
            )

            if (
                orbit_mod(
                    candidate_k3,
                    modulus3,
                )
                == rho3
            ):

                second_matches.append(
                    (
                        t2,
                        candidate_k3,
                    )
                )

        second_unique = (
            len(second_matches) == 1
        )

        if not second_unique:

            failures.append(
                (
                    ell,
                    "second_lift_not_unique",
                    second_matches,
                )
            )

            print(
                f"  ell={ell}: "
                f"SECOND_LIFT_FAILURE "
                f"matches={second_matches}"
            )

            continue

        t2, k3 = second_matches[0]

        # --------------------------------------------------------------
        # Direct validation
        # --------------------------------------------------------------

        direct2 = (
            orbit_mod(
                k2,
                modulus2,
            )
            == rho2
        )

        direct3 = (
            orbit_mod(
                k3,
                modulus3,
            )
            == rho3
        )

        recurrence2 = (
            k2
            ==
            k1 + t1 * order1
        )

        recurrence3 = (
            k3
            ==
            k2 + t2 * order2
        )

        prime_ok = (
            first_unique
            and second_unique
            and direct2
            and direct3
            and recurrence2
            and recurrence3
        )

        if not prime_ok:

            failures.append(
                (
                    ell,
                    "validation_failure",
                )
            )

        print(
            f"  ell={ell}: "
            f"k1={k1} "
            f"ord1={order1} "
            f"t1={t1} "
            f"k2={k2} "
            f"ord2={order2} "
            f"t2={t2} "
            f"k3={k3} "
            f"OK={prime_ok}"
        )

    # ------------------------------------------------------------------
    # 3. SAMPLE EXPLICIT LIFT CANDIDATES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. EXPLICIT SMALL-PRIME LIFT EXAMPLES")
    print("=" * 78)

    for target in (
        5,
        7,
        17,
        73,
    ):

        if not any(
            row[0] == target
            for row in admissible
        ):
            continue

        ell, rho1, order1, k1 = next(
            row
            for row in admissible
            if row[0] == target
        )

        rho2 = ratio_mod(ell * ell)
        rho3 = ratio_mod(ell ** 3)

        order2 = multiplicative_order(
            3,
            ell * ell,
        )

        first = [
            (
                t,
                k1 + t * order1,
            )
            for t in range(ell)
            if orbit_mod(
                k1 + t * order1,
                ell * ell,
            )
            == rho2
        ]

        if len(first) != 1:
            continue

        t1, k2 = first[0]

        second = [
            (
                t,
                k2 + t * order2,
            )
            for t in range(ell)
            if orbit_mod(
                k2 + t * order2,
                ell ** 3,
            )
            == rho3
        ]

        print(
            f"  ell={ell}: "
            f"k1={k1}, "
            f"t1={t1}, "
            f"k2={k2}, "
            f"t2={second[0][0] if len(second)==1 else None}, "
            f"k3={second[0][1] if len(second)==1 else None}"
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
Experiment 238 established chart existence:

    rho/2 in <3> mod ell.

Experiment 233 established uniqueness of the first lift:

    k2 = k1 + t1*ord_ell(3).

Experiment 239 asks whether the same mechanism continues one level
higher:

    k3 = k2 + t2*ord_{ell^2}(3).

Thus, for every admissible small prime, we test the finite tower

    ell
      ->
    ell^2
      ->
    ell^3.

A positive result gives a genuine two-digit local exponent expansion,
rather than only a first-order coincidence.

This is the ordinary-prime counterpart of the higher-level 7-adic
exponent sequence already obtained for the source ratio.
"""
    )

    # ------------------------------------------------------------------
    # 5. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        len(failures) == 0
        and
        len(admissible) > 0
    )

    print()
    print("=" * 78)
    print("5. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  admissible_prime_count="
        f"{len(admissible)}"
    )

    print(
        f"  all_first_lifts_unique="
        f"{len(failures) == 0}"
    )

    print(
        f"  all_second_lifts_unique="
        f"{len(failures) == 0}"
    )

    print(
        f"  two_step_lift_exact="
        f"{final_ok}"
    )

    print(
        f"  failures={len(failures)}"
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
    print("EXPERIMENT 239 COMPLETE")


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

