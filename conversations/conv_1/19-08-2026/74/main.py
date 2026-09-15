#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 242 — EXACT PRIME-LOCAL EXPONENT DIGIT RECONSTRUCTION AUDIT
==============================================================================

Experiment 239:
    unique lifts through ell^3.

Experiment 240:
    exact derivative law at the ell^2 -> ell^3 step.

Experiment 241:
    exact derivative law at the ell^3 -> ell^4 step.

Experiment 242 asks whether the complete exponent sequence can be
reconstructed solely from the local residual digits.

For each admissible prime ell:

    k1 = base exponent modulo ord_ell(3)

and recursively

    k_{r+1}
      = k_r
        + t_r * ord_{ell^r}(3),

where

    R_r(t)
      = A_r + B_r*t (mod ell),

and therefore

    t_r
      = -A_r / B_r (mod ell).

The experiment computes the digits algebraically rather than choosing
them by direct search.

It then independently validates the reconstructed exponent at every
tested modulus:

    ell,
    ell^2,
    ell^3,
    ell^4.

This is the finite-prime analogue of the canonical base-7 exponent
reconstruction already established for the original source.

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
        f"Could not determine order of {base} modulo {prime}."
    )


def lift_order_once(
    base: int,
    prime: int,
    previous_order: int,
    current_power: int,
) -> int:

    """
    Given ord_{p^r}(base)=previous_order,
    the next order is either previous_order or p*previous_order.
    """

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
# FIND AFFINE DIGIT ALGEBRAICALLY
# ============================================================================

def recover_digit(
    ell: int,
    k_current: int,
    order_current: int,
    level: int,
    rho: int,
) -> tuple[
    int,
    int,
    int,
    bool,
]:

    """
    Current level is ell^level.

    Candidate exponents:

        k_current + t*order_current.

    Compute

        R(t)=F/ell^level mod ell.

    Return:

        A,
        B,
        recovered_t,
        affine_exact.
    """

    modulus_current = ell ** level
    modulus_next = ell ** (
        level + 1
    )

    normalized = []

    for t in range(ell):

        candidate = int(
            k_current
            + t * order_current
        )

        F = residual_mod(
            candidate,
            modulus_next,
        )

        if F % modulus_current != 0:
            return (
                0,
                0,
                0,
                False,
            )

        R = int(
            (
                F // modulus_current
            )
            % ell
        )

        normalized.append(R)

    A = normalized[0]

    slopes = [
        int(
            (
                normalized[i]
                - normalized[i - 1]
            )
            % ell
        )
        for i in range(
            1,
            ell
        )
    ]

    slope_set = sorted(
        set(slopes)
    )

    if len(slope_set) != 1:

        return (
            A,
            0,
            0,
            False,
        )

    B = slope_set[0]

    if gcd(
        B,
        ell,
    ) != 1:

        return (
            A,
            B,
            0,
            False,
        )

    affine_exact = all(
        normalized[t]
        ==
        (
            A
            + t * B
        ) % ell
        for t in range(ell)
    )

    if not affine_exact:
        return (
            A,
            B,
            0,
            False,
        )

    recovered_t = int(
        (
            (-A)
            * pow(
                B,
                -1,
                ell,
            )
        )
        % ell
    )

    return (
        A,
        B,
        recovered_t,
        True,
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 242 — EXACT PRIME-LOCAL EXPONENT DIGIT "
        "RECONSTRUCTION AUDIT"
    )
    print("=" * 78)

    failures = []
    admissible_count = 0

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
    # PRIME LOOP
    # ------------------------------------------------------------------

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

        order1 = multiplicative_order_prime(
            3,
            ell,
        )

        # --------------------------------------------------------------
        # Base chart
        # --------------------------------------------------------------

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

        admissible_count += 1

        k = int(
            base_matches[0]
        )

        current_order = int(
            order1
        )

        exponent_sequence = [
            k
        ]

        digit_sequence = []

        intercept_sequence = []
        slope_sequence = []

        all_levels_exact = True

        print()
        print("=" * 78)
        print(
            f"2. LOCAL EXPONENT RECONSTRUCTION ell={ell}"
        )
        print("=" * 78)

        print(
            f"  level=1 "
            f"k={k} "
            f"order={current_order} "
            f"alignment="
            f"{orbit_mod(k,ell) == ratio_mod(ell)}"
        )

        # --------------------------------------------------------------
        # Levels 1 -> 2 -> 3 -> 4
        # --------------------------------------------------------------

        for level in (
            1,
            2,
            3,
        ):

            A, B, t, affine_exact = recover_digit(
                ell,
                k,
                current_order,
                level,
                rho,
            )

            if not affine_exact:

                failures.append(
                    (
                        ell,
                        level,
                        "affine_recovery_failure",
                    )
                )

                all_levels_exact = False

                print(
                    f"  level={level}: "
                    f"A={A} "
                    f"B={B} "
                    f"affine_exact=False"
                )

                break

            next_k = int(
                k
                + t * current_order
            )

            next_modulus = ell ** (
                level + 1
            )

            next_rho = ratio_mod(
                next_modulus
            )

            next_orbit = orbit_mod(
                next_k,
                next_modulus,
            )

            alignment = (
                next_orbit
                == next_rho
            )

            intercept_sequence.append(A)
            slope_sequence.append(B)
            digit_sequence.append(t)
            exponent_sequence.append(next_k)

            print(
                f"  level={level}: "
                f"k={k} "
                f"A={A} "
                f"B={B} "
                f"t={t} "
                f"next_k={next_k} "
                f"alignment={alignment}"
            )

            if not alignment:

                failures.append(
                    (
                        ell,
                        level,
                        "modular_alignment_failure",
                    )
                )

                all_levels_exact = False
                break

            # ----------------------------------------------------------
            # Advance order to ord_{ell^(level+1)}(3).
            # ----------------------------------------------------------

            current_order = lift_order_once(
                3,
                ell,
                current_order,
                level,
            )

            k = next_k

        # --------------------------------------------------------------
        # Independent direct validation of final sequence
        # --------------------------------------------------------------

        direct_checks = []

        for E, exponent in enumerate(
            exponent_sequence,
            start=1,
        ):

            modulus = ell ** E

            rho_E = ratio_mod(
                modulus
            )

            orbit_E = orbit_mod(
                exponent,
                modulus,
            )

            exact = (
                orbit_E
                == rho_E
            )

            direct_checks.append(
                exact
            )

        direct_exact = all(
            direct_checks
        )

        if not direct_exact:

            failures.append(
                (
                    ell,
                    "direct_sequence_validation",
                )
            )

            all_levels_exact = False

        print(
            f"  exponent_sequence="
            f"{exponent_sequence}"
        )

        print(
            f"  lift_digits="
            f"{digit_sequence}"
        )

        print(
            f"  intercepts="
            f"{intercept_sequence}"
        )

        print(
            f"  slopes="
            f"{slope_sequence}"
        )

        print(
            f"  direct_checks="
            f"{direct_checks}"
        )

        print(
            f"  local_reconstruction_exact="
            f"{all_levels_exact}"
        )

    # ------------------------------------------------------------------
    # SUMMARY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. GLOBAL SUMMARY")
    print("=" * 78)

    print(
        f"  admissible_prime_count="
        f"{admissible_count}"
    )

    print(
        f"  local_failures="
        f"{len(failures)}"
    )

    print(
        f"  failures_detail="
        f"{failures}"
    )

    # ------------------------------------------------------------------
    # IMPORTANT SPECIAL CASES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. SPECIAL PRIME COMPARISON")
    print("=" * 78)

    for target in (
        5,
        7,
        17,
        73,
    ):

        if target == 3:
            continue

        rho = ratio_mod(target)

        if rho is None:
            continue

        order1 = multiplicative_order_prime(
            3,
            target,
        )

        base = [
            k
            for k in range(order1)
            if orbit_mod(k,target) == rho
        ]

        if len(base) != 1:
            continue

        k = base[0]
        order_current = order1
        digits = [ ]

        for level in (
            1,
            2,
            3,
        ):

            A, B, t, exact = recover_digit(
                target,
                k,
                order_current,
                level,
                rho,
            )

            if not exact:
                break

            digits.append(t)

            k = int(
                k + t * order_current
            )

            order_current = lift_order_once(
                3,
                target,
                order_current,
                level,
            )

        print(
            f"  ell={target}: "
            f"digits={digits} "
            f"final_k={k}"
        )

    # ------------------------------------------------------------------
    # STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
Experiment 242 reverses the perspective of Experiments 239-241.

Instead of merely observing the exponents

    k1, k2, k3, k4,

it reconstructs them entirely from the local residual digits.

At each level:

    R_r(t)
      = A_r + B_r*t (mod ell),

and therefore

    t_r
      = -A_r/B_r (mod ell).

The exponent recurrence is then

    k_{r+1}
      = k_r
        + t_r*ord_{ell^r}(3).

Thus the local exponent is encoded by a finite digit sequence

    (t1,t2,t3,...).

For the tested primes this is the exact ordinary-prime analogue of
the canonical 7-adic digit expansion obtained earlier.

The important distinction is that the digits are reconstructed from
the residuals, while the resulting exponents are validated independently
against the original source ratio.
"""
    )

    # ------------------------------------------------------------------
    # FINAL
    # ------------------------------------------------------------------

    final_ok = (
        admissible_count > 0
        and len(failures) == 0
    )

    print()
    print("=" * 78)
    print("6. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  admissible_prime_count="
        f"{admissible_count}"
    )

    print(
        f"  all_prime_local_digit_reconstructions_exact="
        f"{final_ok}"
    )

    print(
        f"  all_independent_modular_validations_exact="
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
    print("EXPERIMENT 242 COMPLETE")


if __name__ == "__main__":

    try:
        main()

    except KeyboardInterrupt:

        print("\nInterrupted.")
        sys.exit(130)

    except Exception as exc:

        print(
            "\nFATAL ERROR: "
            + type(exc).__name__
            + ": "
            + str(exc)
        )

        raise

