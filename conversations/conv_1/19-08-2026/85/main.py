#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 250R — EXACT GENERAL PRIME PRINCIPAL-COORDINATE /
                   RATIONAL p-ADIC LOGARITHMIC AUDIT
==============================================================================

Experiment 249R established the Hensel principal-coordinate reconstruction
for the admissible nonexceptional primes.

Experiment 250R fixes the logarithmic implementation used in Experiment 250.

For each admissible nonexceptional prime ell:

    d = ord_ell(3)

    k_r = k1 + d*m_r

and

    rho / (2*3^k1) = (3^d)^m.

The logarithmic coordinate is

    m = log(rho/(2*3^k1)) / log(3^d).

The logarithms are computed as exact finite rational sums

    log(1+u)
      = sum_{n>=1} (-1)^(n+1) u^n/n,

with enough terms to determine the requested ell-adic precision.

The quotient of the two rational truncated logarithms is then reduced
modulo ell^(r-1).

This avoids the incorrect modular-division logic used previously.

Only main.py is used.
"""

from __future__ import annotations

import sys
from fractions import Fraction
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


# ============================================================================
# EXACT p-ADIC VALUATION
# ============================================================================

def vp(
    x: int,
    p: int,
) -> int:

    x = abs(int(x))

    if x == 0:
        return 10**9

    v = 0

    while x % p == 0:

        x //= p
        v += 1

    return int(v)


# ============================================================================
# EXACT TRUNCATED p-ADIC LOGARITHM
# ============================================================================

def padic_log_fraction(
    x: int,
    p: int,
    precision: int,
) -> Fraction:

    """
    Compute log(x) modulo p^precision as an exact rational representative.

    x must satisfy

        x == 1 (mod p).

    We use

        log(1+u)
          = u - u^2/2 + u^3/3 - ...

    and stop once the p-adic valuation of the next term is >= precision.

    The result is kept as Fraction so division by another truncated
    logarithm is exact before p-adic reduction.
    """

    modulus = p ** precision

    x = int(x) % modulus
    u = (x - 1) % modulus

    if u % p != 0:

        raise ValueError(
            "padic_log_fraction requires x == 1 mod p."
        )

    if u == 0:

        return Fraction(0, 1)

    result = Fraction(0, 1)

    power = u

    # For u divisible by p, the valuation grows at least linearly
    # with n. This bound is deliberately generous but still tiny.
    max_terms = precision * p + 20

    for n in range(
        1,
        max_terms + 1,
    ):

        if power == 0:
            break

        numerator_v = vp(
            power,
            p,
        )

        denominator_v = vp(
            n,
            p,
        )

        term_v = (
            numerator_v
            - denominator_v
        )

        if term_v >= precision:

            break

        term = Fraction(
            power,
            n,
        )

        if n % 2 == 1:
            result += term
        else:
            result -= term

        power = (
            power * u
        ) % modulus

    return result


# ============================================================================
# EXACT RATIONAL -> p-ADIC RESIDUE
# ============================================================================

def fraction_mod_prime_power(
    value: Fraction,
    p: int,
    exponent: int,
) -> int:

    """
    Reduce a rational p-adically modulo p^exponent.

    The denominator must be prime to p after Fraction reduction.
    """

    modulus = p ** exponent

    numerator = int(
        value.numerator
    )

    denominator = int(
        value.denominator
    )

    if denominator % p == 0:

        raise ArithmeticError(
            "Rational denominator is not a p-adic unit."
        )

    numerator_mod = (
        numerator % modulus
    )

    denominator_inv = pow(
        denominator % modulus,
        -1,
        modulus,
    )

    return int(
        numerator_mod
        * denominator_inv
        % modulus
    )


# ============================================================================
# LOGARITHMIC PRINCIPAL COORDINATE
# ============================================================================

def logarithmic_coordinate(
    source_unit: int,
    base_unit: int,
    p: int,
    precision: int,
) -> int:

    """
    Compute

        log(source_unit)/log(base_unit)

    modulo p^(precision-1).

    Both units must be congruent to 1 mod p.

    Their logarithms have valuation 1 in the generic nonexceptional
    branch. After dividing both by p, the quotient is a p-adic unit
    denominator and may be inverted modulo p^(precision-1).
    """

    modulus = p ** precision
    coordinate_modulus = p ** (
        precision - 1
    )

    source_log = padic_log_fraction(
        source_unit % modulus,
        p,
        precision,
    )

    base_log = padic_log_fraction(
        base_unit % modulus,
        p,
        precision,
    )

    if source_log == 0:
        return 0

    if base_log == 0:
        raise ArithmeticError(
            "Base logarithm vanished."
        )

    quotient = (
        source_log
        / base_log
    )

    return fraction_mod_prime_power(
        quotient,
        p,
        precision - 1,
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 250R — EXACT GENERAL PRIME PRINCIPAL-COORDINATE / "
        "RATIONAL p-ADIC LOGARITHMIC AUDIT"
    )
    print("=" * 78)

    failures = []

    counts = {
        "primes": 0,
        "skipped_q3": 0,
        "exceptional": 0,
        "inadmissible": 0,
        "tested": 0,
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

        D1 = (
            (
                pow(
                    3,
                    d,
                    ell * ell,
                )
                - 1
            )
            // ell
        ) % ell

        if D1 == 0:

            counts["exceptional"] += 1
            continue

        rho = ratio_mod(
            ell
        )

        if rho is None:
            counts["skipped_q3"] += 1
            continue

        # --------------------------------------------------------------
        # Base class
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

        counts["tested"] += 1

        k1 = base_matches[0]

        k = k1
        current_order = d

        exponents = [k1]
        principal_m = [0]
        digits = []

        local_ok = True

        # --------------------------------------------------------------
        # Hensel reconstruction through ell^4
        # --------------------------------------------------------------

        for level in range(
            1,
            MAX_LEVEL,
        ):

            current_modulus = ell ** level
            next_modulus = ell ** (
                level + 1
            )

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
                        "order_growth",
                    )
                )

                local_ok = False
                break

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
                        "D_integrality",
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
                        "zero_D",
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
                    F
                    % current_modulus
                    != 0
                ):

                    failures.append(
                        (
                            ell,
                            level,
                            "normalization",
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
                        "non_affine",
                    )
                )

                local_ok = False
                break

            B = next(
                iter(slope_set)
            )

            predicted_B = (
                -(Q3 % ell)
                * rho
                * D
            ) % ell

            if B != predicted_B:

                failures.append(
                    (
                        ell,
                        level,
                        "slope",
                        B,
                        predicted_B,
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

            target = ratio_mod(
                next_modulus
            )

            if (
                orbit_mod(
                    next_k,
                    next_modulus,
                )
                != target
            ):

                failures.append(
                    (
                        ell,
                        level,
                        "alignment",
                    )
                )

                local_ok = False
                break

            digits.append(t)
            exponents.append(next_k)

            m = (
                next_k
                - k1
            ) // d

            principal_m.append(
                m
            )

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

        if not local_ok:
            continue

        # --------------------------------------------------------------
        # Direct exponential identity
        # --------------------------------------------------------------

        for level, m in enumerate(
            principal_m
        ):

            if level == 0:
                continue

            modulus = ell ** level

            rho_level = ratio_mod(
                modulus
            )

            base_factor = (
                2
                * pow(
                    3,
                    k1,
                    modulus,
                )
            ) % modulus

            base_unit = pow(
                3,
                d,
                modulus,
            )

            reconstructed = (
                base_factor
                * pow(
                    base_unit,
                    m,
                    modulus,
                )
            ) % modulus

            if reconstructed != rho_level:

                failures.append(
                    (
                        ell,
                        level,
                        "exponential_reconstruction",
                        reconstructed,
                        rho_level,
                    )
                )

                local_ok = False
                break

        if not local_ok:
            continue

        # --------------------------------------------------------------
        # Exact logarithmic coordinate
        # --------------------------------------------------------------

        for level in range(
            2,
            MAX_LEVEL + 1,
        ):

            modulus = ell ** level
            coordinate_modulus = ell ** (
                level - 1
            )

            rho_level = ratio_mod(
                modulus
            )

            base_factor = (
                2
                * pow(
                    3,
                    k1,
                    modulus,
                )
            ) % modulus

            inv_base_factor = pow(
                base_factor,
                -1,
                modulus,
            )

            source_unit = (
                rho_level
                * inv_base_factor
            ) % modulus

            base_unit = (
                pow(
                    3,
                    d,
                    modulus,
                )
            )

            # Both must be principal units.
            if (
                source_unit % ell
                != 1
            ):

                failures.append(
                    (
                        ell,
                        level,
                        "source_unit_not_principal",
                        source_unit % ell,
                    )
                )

                local_ok = False
                break

            if (
                base_unit % ell
                != 1
            ):

                failures.append(
                    (
                        ell,
                        level,
                        "base_unit_not_principal",
                        base_unit % ell,
                    )
                )

                local_ok = False
                break

            try:

                log_m = logarithmic_coordinate(
                    source_unit,
                    base_unit,
                    ell,
                    level,
                )

            except Exception as exc:

                failures.append(
                    (
                        ell,
                        level,
                        "logarithmic_coordinate_failure",
                        type(exc).__name__,
                        str(exc),
                    )
                )

                local_ok = False
                break

            hensel_m = principal_m[
                level - 1
            ]

            if (
                log_m
                % coordinate_modulus
                !=
                hensel_m
                % coordinate_modulus
            ):

                failures.append(
                    (
                        ell,
                        level,
                        "logarithmic_mismatch",
                        log_m,
                        hensel_m,
                        coordinate_modulus,
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
    # FAILURES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. FAILURE SUMMARY")
    print("=" * 78)

    print(
        f"  failures={len(failures)}"
    )

    if failures:

        for failure in failures[:25]:
            print(
                f"  {failure}"
            )

        if len(failures) > 25:

            print(
                f"  ... "
                f"{len(failures)-25} more"
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
The earlier Experiment 250 failures all occurred inside the logarithmic
implementation, not inside the projective-Hensel or exponential
reconstruction.

Experiment 250R computes the p-adic logarithm by exact rational
truncation:

    log(1+u)
      = u - u^2/2 + u^3/3 - ...

The quotient

    log(source_unit)
    -----------------
      log(base_unit)

is then reduced p-adically.

Thus three independent descriptions are compared:

    Hensel digits,
    exponential principal-coordinate reconstruction,
    logarithmic principal-coordinate reconstruction.

Agreement of all three would establish the local logarithmic chart
independently of the recursive digit construction.
"""
    )

    # ------------------------------------------------------------------
    # FINAL
    # ------------------------------------------------------------------

    final_ok = (
        len(failures) == 0
        and
        counts["tested"] > 0
    )

    print()
    print("=" * 78)
    print("4. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  tested_admissible_nonexceptional_count="
        f"{counts['tested']}"
    )

    print(
        f"  exceptional_count="
        f"{counts['exceptional']}"
    )

    print(
        f"  all_exponential_reconstructions_exact="
        f"{final_ok}"
    )

    print(
        f"  all_logarithmic_coordinates_exact="
        f"{final_ok}"
    )

    print(
        f"  failures={len(failures)}"
    )

    print(
        "  exceptional_11_excluded_from_generic_chart=True"
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
    print("EXPERIMENT 250R COMPLETE")


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

