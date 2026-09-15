#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 249 — FAST GLOBAL PROJECTIVE-HENSEL / PRINCIPAL-COORDINATE AUDIT
==============================================================================

This experiment combines the main local results proved experimentally so far.

For every odd prime ell >= 5 up to PRIME_LIMIT:

    rho = q1/q3 (mod ell)

If rho lies in the orbit 2*<3> modulo ell, there is a unique base exponent

    k1 mod ord_ell(3).

For nonexceptional ell:

    D_ell != 0,

and

    ord_{ell^r}(3)
        = ord_ell(3) * ell^(r-1).

At each level define

    R_r(t)
      = F(k_r + t*ord_{ell^r}(3)) / ell^r
      (mod ell),

where

    F(k) = q1 - 2*3^k*q3.

The local affine law should be

    R_r(t)
      = A_r + B*t
      (mod ell),

with constant slope

    B = -q3*rho*D_r
      (mod ell),

where

    D_r
      = (3^ord_{ell^r}(3)-1)/ell^r
      (mod ell).

The lift digit is therefore

    t_r = -A_r/B mod ell.

Finally, writing

    k_r = k1 + d*m_r,

    d = ord_ell(3),

the nonexceptional order growth implies

    m_E
      = t_1 + t_2*ell + t_3*ell^2 + ...

This experiment checks all of these statements together.

It deliberately skips:

    ell = 3,
    ell | q3,
    orbit-inadmissible primes,
    exceptional primes.

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
        "EXPERIMENT 249 — FAST GLOBAL PROJECTIVE-HENSEL / "
        "PRINCIPAL-COORDINATE AUDIT"
    )
    print("=" * 78)

    failures = []

    counts = {
        "primes": 0,
        "skipped_3": 0,
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

        # --------------------------------------------------------------
        # q3-unit condition
        # --------------------------------------------------------------

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
        # Base orbit membership
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

        k = base_matches[0]

        orders = [d]
        digits = []
        principal = [0]

        local_ok = True

        # --------------------------------------------------------------
        # Hensel levels
        # --------------------------------------------------------------

        for level in range(
            1,
            MAX_LEVEL,
        ):

            current_order = orders[
                level - 1
            ]

            current_modulus = ell ** level
            next_modulus = ell ** (
                level + 1
            )

            # ----------------------------------------------------------
            # Order growth
            # ----------------------------------------------------------

            expected_order = (
                d
                * ell ** (
                    level - 1
                )
            )

            if current_order != expected_order:

                local_ok = False

                failures.append(
                    (
                        ell,
                        level,
                        "order_growth_failure",
                        current_order,
                        expected_order,
                    )
                )

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

                local_ok = False

                failures.append(
                    (
                        ell,
                        level,
                        "D_integrality_failure",
                    )
                )

                break

            D = int(
                (
                    order_increment
                    // current_modulus
                )
                % ell
            )

            if D == 0:

                local_ok = False

                failures.append(
                    (
                        ell,
                        level,
                        "unexpected_zero_D",
                    )
                )

                break

            # ----------------------------------------------------------
            # Local affine residual
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
                    F
                    % current_modulus
                    != 0
                ):

                    local_ok = False

                    failures.append(
                        (
                            ell,
                            level,
                            "normalization_failure",
                            t,
                        )
                    )

                    break

                normalized.append(
                    (
                        F
                        // current_modulus
                    )
                    % ell
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

                local_ok = False

                failures.append(
                    (
                        ell,
                        level,
                        "non_affine_residual",
                    )
                )

                break

            measured_B = next(
                iter(slope_set)
            )

            rho_mod = rho
            q3_mod = Q3 % ell

            predicted_B = (
                -q3_mod
                * rho_mod
                * D
            ) % ell

            if (
                measured_B
                != predicted_B
            ):

                local_ok = False

                failures.append(
                    (
                        ell,
                        level,
                        "slope_formula_failure",
                        measured_B,
                        predicted_B,
                    )
                )

                break

            # ----------------------------------------------------------
            # Full affine identity
            # ----------------------------------------------------------

            affine_ok = all(
                normalized[t]
                ==
                (
                    A
                    + t * measured_B
                ) % ell
                for t in range(ell)
            )

            if not affine_ok:

                local_ok = False

                failures.append(
                    (
                        ell,
                        level,
                        "full_affine_failure",
                    )
                )

                break

            # ----------------------------------------------------------
            # Recover lift digit
            # ----------------------------------------------------------

            if gcd(
                measured_B,
                ell,
            ) != 1:

                local_ok = False

                failures.append(
                    (
                        ell,
                        level,
                        "zero_slope",
                    )
                )

                break

            t = (
                (-A)
                * pow(
                    measured_B,
                    -1,
                    ell,
                )
            ) % ell

            next_k = (
                k
                + t * current_order
            )

            # ----------------------------------------------------------
            # Independent modular check
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

                local_ok = False

                failures.append(
                    (
                        ell,
                        level,
                        "lift_alignment_failure",
                    )
                )

                break

            digits.append(
                t
            )

            k = next_k

            principal.append(
                (
                    k
                    - base_matches[0]
                )
                // d
            )

            # ----------------------------------------------------------
            # Actual next order
            # ----------------------------------------------------------

            phi_next = (
                (ell - 1)
                * ell ** level
            )

            next_order = multiplicative_order(
                3,
                next_modulus,
                phi_next,
            )

            orders.append(
                next_order
            )

        # --------------------------------------------------------------
        # Principal coordinate digit test
        # --------------------------------------------------------------

        if local_ok:

            for r, m in enumerate(
                principal
            ):

                value = m
                base_digits = []

                if value == 0:
                    base_digits = [0]
                else:
                    while value > 0:
                        base_digits.append(
                            value % ell
                        )
                        value //= ell

                expected_digits = digits[:r]

                # Empty representation is equivalent to zero.
                if r == 0:
                    digit_ok = (
                        m == 0
                    )
                else:
                    while (
                        len(base_digits) > 1
                        and
                        base_digits[-1] == 0
                    ):
                        base_digits.pop()

                    digit_ok = (
                        base_digits
                        ==
                        expected_digits
                    )

                if not digit_ok:

                    local_ok = False

                    failures.append(
                        (
                            ell,
                            r,
                            "principal_digit_failure",
                            base_digits,
                            expected_digits,
                        )
                    )

                    break

        if not local_ok:
            continue

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
    # FAILURE DETAILS
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
                f"{len(failures)-20} more"
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
Experiment 249 combines the local results into one finite audit.

For every odd prime ell in the nonexceptional, orbit-admissible branch:

    rho = q1/q3 mod ell

determines a unique exponent class

    k1 mod ord_ell(3).

At each successive level:

    R_r(t)
      = A_r + B*t mod ell,

with

    B
      = -q3*rho*D_r mod ell,

and

    t
      = -A_r/B mod ell.

The exponent then lifts recursively.

Because the order grows as

    ord_{ell^r}(3)
      = ord_ell(3)*ell^(r-1),

the principal coordinate

    k_r = k1 + d*m_r

has the ordinary base-ell expansion

    m_r
      = t_1
        + t_2*ell
        + t_3*ell^2
        + ...

The exceptional prime 11 is deliberately excluded from this generic
branch.

This experiment therefore tests the complete local theorem template,
rather than an isolated component.
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
        "  generic_projective_hensel_law_exact="
        f"{final_ok}"
    )

    print(
        "  generic_principal_digit_law_exact="
        f"{final_ok}"
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
    print("EXPERIMENT 249 COMPLETE")


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

