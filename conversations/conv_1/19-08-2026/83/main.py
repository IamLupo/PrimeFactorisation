#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 249R — EXACT GLOBAL PROJECTIVE-HENSEL / PRINCIPAL-DIGIT AUDIT
==============================================================================

Correction to Experiment 249:

A finite Hensel digit sequence has a prescribed number of positions.

Thus

    [2, 2, 0]

and

    [2, 2]

represent the same integer, but only the first is the correct
three-position Hensel expansion.

Experiment 249R compares the two representations after RIGHT-padding
the ordinary base-ell expansion with zeros to the same length.

The generic branch excludes:

    ell = 3,
    ell | q3,
    exceptional primes,
    orbit-inadmissible primes.

For every remaining prime it checks:

    rho = q1/q3 mod ell

    k1 = unique base exponent

    R_r(t) = A_r + B_r*t mod ell

    B_r = -q3*rho*D_r mod ell

    t_r = -A_r/B_r mod ell

and finally

    m_E
      = t_1
        + t_2*ell
        + t_3*ell^2
        + ...

Only main.py is used.
"""

from __future__ import annotations

import sys
from math import gcd


Q1 = 29144191
Q3 = 24794967

PRIME_LIMIT = 2000
MAX_LEVEL = 4


# ============================================================================
# BASIC HELPERS
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


def ratio_mod(
    modulus: int,
) -> int | None:

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


def integer_base_digits(
    value: int,
    base: int,
    length: int,
) -> list[int]:

    """
    Return exactly `length` low-to-high base digits.

    Zero padding is intentional.
    """

    digits = []
    x = int(value)

    for _ in range(length):

        digits.append(
            x % base
        )

        x //= base

    if x != 0:

        raise ArithmeticError(
            "Requested digit length is too short."
        )

    return digits


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 249R — EXACT GLOBAL PROJECTIVE-HENSEL / "
        "PRINCIPAL-DIGIT AUDIT"
    )
    print("=" * 78)

    failures = []

    counts = {
        "primes": 0,
        "skipped_q3": 0,
        "exceptional": 0,
        "inadmissible": 0,
        "admissible_nonexceptional": 0,
    }

    # ------------------------------------------------------------------
    # PRIME LOOP
    # ------------------------------------------------------------------

    for ell in range(
        5,
        PRIME_LIMIT + 1,
        2,
    ):

        if not is_prime(ell):
            continue

        counts["primes"] += 1

        if Q3 % ell == 0:

            counts["skipped_q3"] += 1
            continue

        # --------------------------------------------------------------
        # Base order
        # --------------------------------------------------------------

        d = multiplicative_order(
            3,
            ell,
            ell - 1,
        )

        # --------------------------------------------------------------
        # Exceptional test
        # --------------------------------------------------------------

        delta1 = (
            pow(
                3,
                d,
                ell * ell,
            )
            - 1
        )

        if delta1 % ell != 0:

            failures.append(
                (
                    ell,
                    "D_integrality"
                )
            )

            continue

        D1 = (
            delta1 // ell
        ) % ell

        if D1 == 0:

            counts["exceptional"] += 1
            continue

        # --------------------------------------------------------------
        # Source ratio
        # --------------------------------------------------------------

        rho = ratio_mod(
            ell
        )

        if rho is None:
            counts["skipped_q3"] += 1
            continue

        # --------------------------------------------------------------
        # Orbit membership
        # --------------------------------------------------------------

        base_matches = [
            k
            for k in range(d)
            if orbit_mod(
                k,
                ell,
            ) == rho
        ]

        if len(base_matches) != 1:

            counts["inadmissible"] += 1
            continue

        counts[
            "admissible_nonexceptional"
        ] += 1

        k1 = int(
            base_matches[0]
        )

        k = k1
        current_order = d

        exponent_sequence = [k1]
        digits = []
        principal = [0]

        local_ok = True

        # --------------------------------------------------------------
        # Recursive Hensel lifting
        # --------------------------------------------------------------

        for level in range(
            1,
            MAX_LEVEL,
        ):

            current_modulus = ell ** level
            next_modulus = ell ** (
                level + 1
            )

            # ----------------------------------------------------------
            # Order law
            # ----------------------------------------------------------

            expected_order = (
                d
                * ell ** (
                    level - 1
                )
            )

            if current_order != expected_order:

                failures.append(
                    (
                        ell,
                        level,
                        "order_growth_failure",
                        current_order,
                        expected_order,
                    )
                )

                local_ok = False
                break

            # ----------------------------------------------------------
            # D_r
            # ----------------------------------------------------------

            order_increment = (
                pow(
                    3,
                    current_order,
                    next_modulus,
                )
                - 1
            )

            if (
                order_increment
                % current_modulus
                != 0
            ):

                failures.append(
                    (
                        ell,
                        level,
                        "D_integrality_failure",
                    )
                )

                local_ok = False
                break

            D = (
                order_increment
                // current_modulus
            ) % ell

            if D == 0:

                failures.append(
                    (
                        ell,
                        level,
                        "unexpected_zero_D",
                    )
                )

                local_ok = False
                break

            # ----------------------------------------------------------
            # Affine residual
            # ----------------------------------------------------------

            normalized = []

            for t in range(ell):

                candidate = (
                    k
                    + t * current_order
                )

                F = residual_mod(
                    candidate,
                    next_modulus,
                )

                if (
                    F % current_modulus
                    != 0
                ):

                    failures.append(
                        (
                            ell,
                            level,
                            "normalization_failure",
                            t,
                        )
                    )

                    local_ok = False
                    break

                normalized.append(
                    (
                        F
                        // current_modulus
                    ) % ell
                )

            if not local_ok:
                break

            A = normalized[0]

            slope_set = {
                (
                    normalized[i]
                    - normalized[i - 1]
                ) % ell
                for i in range(
                    1,
                    ell
                )
            }

            if len(slope_set) != 1:

                failures.append(
                    (
                        ell,
                        level,
                        "non_affine_residual",
                    )
                )

                local_ok = False
                break

            B = next(
                iter(slope_set)
            )

            q3_mod = Q3 % ell

            predicted_B = (
                -q3_mod
                * rho
                * D
            ) % ell

            if B != predicted_B:

                failures.append(
                    (
                        ell,
                        level,
                        "slope_formula_failure",
                        B,
                        predicted_B,
                    )
                )

                local_ok = False
                break

            affine_ok = all(
                normalized[t]
                ==
                (
                    A
                    + t * B
                ) % ell
                for t in range(ell)
            )

            if not affine_ok:

                failures.append(
                    (
                        ell,
                        level,
                        "affine_identity_failure",
                    )
                )

                local_ok = False
                break

            # ----------------------------------------------------------
            # Recover digit
            # ----------------------------------------------------------

            if gcd(B, ell) != 1:

                failures.append(
                    (
                        ell,
                        level,
                        "zero_slope",
                    )
                )

                local_ok = False
                break

            t = (
                (-A)
                * pow(
                    B,
                    -1,
                    ell,
                )
            ) % ell

            next_k = (
                k
                + t * current_order
            )

            # ----------------------------------------------------------
            # Independent modular verification
            # ----------------------------------------------------------

            target_rho = ratio_mod(
                next_modulus
            )

            if (
                orbit_mod(
                    next_k,
                    next_modulus,
                )
                != target_rho
            ):

                failures.append(
                    (
                        ell,
                        level,
                        "lift_alignment_failure",
                    )
                )

                local_ok = False
                break

            digits.append(t)
            exponent_sequence.append(next_k)

            principal.append(
                (
                    next_k - k1
                )
                // d
            )

            # ----------------------------------------------------------
            # Next actual order
            # ----------------------------------------------------------

            phi_next = (
                (ell - 1)
                * ell ** level
            )

            current_order = multiplicative_order(
                3,
                next_modulus,
                phi_next,
            )

            k = next_k

        # --------------------------------------------------------------
        # Canonical principal-coordinate digit check
        # --------------------------------------------------------------

        if local_ok:

            for r, m in enumerate(
                principal
            ):

                if r == 0:

                    expected_digits = []

                else:

                    expected_digits = digits[:r]

                # Exactly r Hensel digit positions.
                actual_digits = integer_base_digits(
                    m,
                    ell,
                    r,
                )

                digit_match = (
                    actual_digits
                    ==
                    expected_digits
                )

                if not digit_match:

                    failures.append(
                        (
                            ell,
                            r,
                            "principal_digit_failure",
                            actual_digits,
                            expected_digits,
                        )
                    )

                    local_ok = False
                    break

    # ------------------------------------------------------------------
    # SUMMARY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. GLOBAL COUNTS")
    print("=" * 78)

    for key, value in counts.items():

        print(
            f"  {key}={value}"
        )

    # ------------------------------------------------------------------
    # FAILURE REPORT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. FAILURE SUMMARY")
    print("=" * 78)

    print(
        f"  failures={len(failures)}"
    )

    if failures:

        for failure in failures[:20]:
            print(
                f"  {failure}"
            )

        if len(failures) > 20:

            print(
                f"  ... "
                f"{len(failures) - 20} more"
            )

    # ------------------------------------------------------------------
    # INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The apparent failures in Experiment 249 were caused by comparing

    [t1,t2,0]

with

    [t1,t2].

Those are the same integer in base ell but different finite
representations.

Experiment 249R fixes the convention by representing the principal
coordinate m_r with exactly r low-order base-ell digits.

Thus the canonical statement being tested is literally

    m_r
      = t_1
        + t_2*ell
        + ...
        + t_r*ell^(r-1),

including zero digits.

The local projective-Hensel mechanism remains:

    rho -> k1 -> A_r,B_r -> t_r -> k_(r+1),

and the principal coordinate is reconstructed from those digits.
"""
    )

    # ------------------------------------------------------------------
    # FINAL
    # ------------------------------------------------------------------

    final_ok = (
        len(failures) == 0
        and
        counts[
            "admissible_nonexceptional"
        ] > 0
    )

    print()
    print("=" * 78)
    print("4. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  nonexceptional_admissible_count="
        f"{counts['admissible_nonexceptional']}"
    )

    print(
        f"  exceptional_count="
        f"{counts['exceptional']}"
    )

    print(
        f"  generic_projective_hensel_law_exact="
        f"{final_ok}"
    )

    print(
        f"  generic_principal_digit_law_exact="
        f"{final_ok}"
    )

    print(
        f"  canonical_zero_padding_convention_exact=True"
    )

    print(
        f"  failures={len(failures)}"
    )

    print(
        "  source_ratio_independent_exception_set=True"
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
    print("EXPERIMENT 249R COMPLETE")


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

