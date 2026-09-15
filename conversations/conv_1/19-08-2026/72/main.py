#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 240 — EXACT SECOND-LEVEL PRIME-LOCAL DERIVATIVE AUDIT
==============================================================================

Experiment 239 established, for every admissible prime ell <= 200:

    k1 -> k2 -> k3

with unique Hensel digits through ell^3.

Experiment 240 tests the local derivative law at the SECOND lift level.

Let

    L2 = ord_{ell^2}(3),

and

    k(t) = k2 + t*L2.

Since k2 already solves the projective equation modulo ell^2,

    F(k) = q1 - 2*3^k*q3

is divisible by ell^2 for every candidate k(t).

Define

    R2(t)
      = F(k(t))/ell^2
      (mod ell).

The predicted slope is

    B2
      = -q3*rho*D2
      (mod ell),

where

    rho = q1/q3 (mod ell),

and

    D2
      = (3^L2 - 1)/ell^2
      (mod ell).

The experiment checks:

    * complete affine law at the second lift;
    * predicted slope;
    * unique second lift digit;
    * reconstruction of k3;
    * agreement with direct modular validation.

This is the exact prime-local analogue of continuing the 7-adic
linearization one more level.

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

def is_prime(
    n: int,
) -> bool:

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
                modulus,
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
                    modulus,
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
        "EXPERIMENT 240 — EXACT SECOND-LEVEL PRIME-LOCAL "
        "DERIVATIVE AUDIT"
    )
    print("=" * 78)

    failures = []

    admissible = []

    # ------------------------------------------------------------------
    # 1. FIND ADMISSIBLE PRIME CHARTS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. ADMISSIBLE PRIME CHARTS")
    print("=" * 78)

    for ell in range(
        2,
        PRIME_LIMIT + 1,
    ):

        if not is_prime(ell):
            continue

        if ell in (
            2,
            3,
        ):
            continue

        if Q3 % ell == 0:
            continue

        rho = ratio_mod(
            ell
        )

        order1 = multiplicative_order(
            3,
            ell,
        )

        base_matches = [
            k
            for k in range(
                order1
            )
            if orbit_mod(
                k,
                ell,
            ) == rho
        ]

        if len(base_matches) != 1:
            continue

        k1 = base_matches[0]

        # First lift.
        modulus2 = ell * ell
        rho2 = ratio_mod(
            modulus2
        )

        lift1 = [
            (
                t,
                k1 + t * order1,
            )
            for t in range(
                ell
            )
            if orbit_mod(
                k1 + t * order1,
                modulus2,
            )
            == rho2
        ]

        if len(lift1) != 1:

            failures.append(
                (
                    ell,
                    "first_lift_not_unique",
                )
            )

            continue

        t1, k2 = lift1[0]

        admissible.append(
            (
                ell,
                rho,
                k1,
                order1,
                t1,
                k2,
            )
        )

        print(
            f"  ell={ell}: "
            f"k1={k1} "
            f"t1={t1} "
            f"k2={k2}"
        )

    print(
        f"  admissible_count="
        f"{len(admissible)}"
    )

    # ------------------------------------------------------------------
    # 2. SECOND-LEVEL LINEARIZATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. SECOND-LEVEL AFFINE RESIDUAL")
    print("=" * 78)

    for (
        ell,
        rho,
        k1,
        order1,
        t1,
        k2,
    ) in admissible:

        modulus2 = ell * ell
        modulus3 = ell ** 3

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

        pow3L2 = pow(
            3,
            order2,
            modulus3,
        )

        increment = (
            pow3L2 - 1
        )

        if (
            increment % modulus2
            != 0
        ):

            failures.append(
                (
                    ell,
                    "second_derivative_not_integral",
                )
            )

            print(
                f"  ell={ell}: "
                f"SECOND_DERIVATIVE_FAILURE"
            )

            continue

        D2 = int(
            (
                increment
                // modulus2
            )
            % ell
        )

        q3_mod = Q3 % ell

        predicted_B2 = int(
            (
                -q3_mod
                * rho
                * D2
            )
            % ell
        )

        # --------------------------------------------------------------
        # Evaluate all ell possible second-level digits.
        # --------------------------------------------------------------

        normalized = []

        division_ok = True

        for t2 in range(
            ell
        ):

            candidate_k3 = int(
                k2
                + t2 * order2
            )

            F_mod = residual_mod(
                candidate_k3,
                modulus3,
            )

            if F_mod % modulus2 != 0:

                division_ok = False
                break

            R2 = int(
                (
                    F_mod
                    // modulus2
                )
                % ell
            )

            normalized.append(
                (
                    t2,
                    R2,
                )
            )

        if not division_ok:

            failures.append(
                (
                    ell,
                    "second_normalization_failure",
                )
            )

            continue

        A2 = normalized[0][1]

        slope_values = [
            int(
                (
                    normalized[i][1]
                    - normalized[i - 1][1]
                )
                % ell
            )
            for i in range(
                1,
                len(normalized)
            )
        ]

        slope_set = sorted(
            set(
                slope_values
            )
        )

        affine_exact = (
            len(slope_set) == 1
        )

        measured_B2 = (
            slope_set[0]
            if affine_exact
            else None
        )

        slope_exact = (
            affine_exact
            and
            measured_B2 == predicted_B2
        )

        # --------------------------------------------------------------
        # Full affine identity.
        # --------------------------------------------------------------

        full_affine_exact = False

        if measured_B2 is not None:

            full_affine_exact = all(
                normalized[t][1]
                ==
                (
                    A2
                    + t * measured_B2
                ) % ell
                for t in range(
                    ell
                )
            )

        # --------------------------------------------------------------
        # Recover t2 algebraically.
        # --------------------------------------------------------------

        if (
            measured_B2 is not None
            and
            gcd(
                measured_B2,
                ell,
            ) == 1
        ):

            recovered_t2 = int(
                (
                    (-A2)
                    * pow(
                        measured_B2,
                        -1,
                        ell,
                    )
                )
                % ell
            )

        else:

            recovered_t2 = None

        # --------------------------------------------------------------
        # Direct lift validation.
        # --------------------------------------------------------------

        if recovered_t2 is not None:

            recovered_k3 = int(
                k2
                + recovered_t2 * order2
            )

            reconstructed = (
                orbit_mod(
                    recovered_k3,
                    modulus3,
                )
                == rho3
            )

            direct_matches = [
                (
                    t,
                    k2 + t * order2,
                )
                for t in range(
                    ell
                )
                if orbit_mod(
                    k2 + t * order2,
                    modulus3,
                )
                == rho3
            ]

            unique_direct = (
                len(direct_matches) == 1
            )

        else:

            recovered_k3 = None
            reconstructed = False
            direct_matches = []
            unique_direct = False

        prime_ok = (
            division_ok
            and
            affine_exact
            and
            slope_exact
            and
            full_affine_exact
            and
            recovered_t2 is not None
            and
            reconstructed
            and
            unique_direct
            and
            direct_matches[0][0]
            == recovered_t2
            if unique_direct
            and recovered_t2 is not None
            else False
        )

        if not prime_ok:

            failures.append(
                (
                    ell,
                    {
                        "affine_exact": affine_exact,
                        "slope_exact": slope_exact,
                        "full_affine_exact": full_affine_exact,
                        "recovered_t2": recovered_t2,
                        "reconstructed": reconstructed,
                        "unique_direct": unique_direct,
                    },
                )
            )

        print(
            f"  ell={ell}: "
            f"k2={k2} "
            f"ord_ell2={order2} "
            f"D2={D2} "
            f"A2={A2} "
            f"B_pred={predicted_B2} "
            f"B_meas={measured_B2} "
            f"t2={recovered_t2} "
            f"k3={recovered_k3} "
            f"OK={prime_ok}"
        )

    # ------------------------------------------------------------------
    # 3. KEY EXAMPLES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. EXPLICIT SECOND-LEVEL EXAMPLES")
    print("=" * 78)

    for target in (
        5,
        7,
        17,
        73,
    ):

        matching = [
            row
            for row in admissible
            if row[0] == target
        ]

        if not matching:
            continue

        ell, rho, k1, order1, t1, k2 = matching[0]

        modulus2 = ell * ell
        modulus3 = ell ** 3

        rho3 = ratio_mod(
            modulus3
        )

        order2 = multiplicative_order(
            3,
            modulus2,
        )

        second = [
            (
                t,
                k2 + t * order2,
            )
            for t in range(
                ell
            )
            if orbit_mod(
                k2 + t * order2,
                modulus3,
            )
            == rho3
        ]

        print(
            f"  ell={ell}: "
            f"k1={k1} "
            f"t1={t1} "
            f"k2={k2} "
            f"t2={second[0][0] if len(second)==1 else None} "
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
Experiment 239 established unique exponent lifting

    ell -> ell^2 -> ell^3.

Experiment 240 asks whether the linearized residual mechanism itself
persists at the second step.

At level ell^2, set

    L2 = ord_{ell^2}(3),

    k(t) = k2 + t*L2.

Since

    3^L2 = 1 + ell^2*D2 (mod ell^3),

the normalized residual

    R2(t)
      = F(k(t))/ell^2 mod ell

must be affine:

    R2(t) = A2 + B2*t mod ell.

The predicted derivative is

    B2 = -q3*rho*D2 mod ell.

If the experiment passes for all admissible primes, then the local
Hensel mechanism is not merely a first-order phenomenon. The same
derivative law continues one full prime-power level higher.
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
        f"  all_second_level_affine_laws_exact="
        f"{final_ok}"
    )

    print(
        f"  all_second_level_derivative_formulas_exact="
        f"{final_ok}"
    )

    print(
        f"  all_second_level_lift_digits_recovered="
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
    print("EXPERIMENT 240 COMPLETE")


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

