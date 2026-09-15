#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 252 — EXACT SYNTHETIC UNIT-SCALE / PROJECTIVE-HENSEL UNIVERSALITY
==============================================================================

Experiment 251 proved the complete local projective/Hensel/logarithmic
mechanism for synthetic sources normalized by

    q3 = 1.

Experiment 252 removes that normalization.

For each selected nonexceptional prime ell and several exponent classes k1,
choose several unit values U and define

    q3 = U,

    rho = 2*3^k1,

    q1 = rho*q3

modulo ell^MAX_LEVEL.

Thus

    q1/q3 = rho

exactly at the tested ell-adic precision, but q3 is no longer 1.

The experiment checks:

    1. base orbit identification;
    2. affine residual law;
    3. full derivative

         B = -q3*rho*D mod ell;

    4. unique Hensel digits;
    5. principal-coordinate base-ell digits;
    6. exponential reconstruction;
    7. logarithmic reconstruction.

This isolates the projective invariance and verifies that the q3 factor
in the derivative is genuinely structural.

Only main.py is used.
"""

from __future__ import annotations

import sys
from fractions import Fraction
from math import gcd


PRIME_LIMIT = 200
MAX_LEVEL = 4

# First few nonexceptional primes.
MAX_PRIMES = 12

# Number of unit scales q3 tested per (ell,k1).
UNIT_SCALES = [1, 2, 3, 5, 7]


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
# EXACT p-ADIC LOGARITHMS
# ============================================================================

def padic_log_fraction(
    x: int,
    p: int,
    precision: int,
) -> Fraction:

    modulus = p ** precision

    x %= modulus

    u = (
        x - 1
    ) % modulus

    if u % p != 0:

        raise ValueError(
            "p-adic logarithm requires x == 1 mod p."
        )

    if u == 0:
        return Fraction(0, 1)

    result = Fraction(0, 1)

    power = u

    max_terms = (
        precision * p
        + 20
    )

    for n in range(
        1,
        max_terms + 1,
    ):

        if power == 0:
            break

        term_v = (
            vp(power, p)
            -
            vp(n, p)
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


def fraction_mod_prime_power(
    value: Fraction,
    p: int,
    exponent: int,
) -> int:

    modulus = p ** exponent

    numerator = (
        value.numerator
        % modulus
    )

    denominator = value.denominator

    if denominator % p == 0:

        raise ArithmeticError(
            "Denominator is not a p-adic unit."
        )

    return int(
        numerator
        * pow(
            denominator % modulus,
            -1,
            modulus,
        )
        % modulus
    )


def logarithmic_coordinate(
    source_unit: int,
    base_unit: int,
    p: int,
    precision: int,
) -> int:

    source_log = padic_log_fraction(
        source_unit,
        p,
        precision,
    )

    base_log = padic_log_fraction(
        base_unit,
        p,
        precision,
    )

    if base_log == 0:

        raise ArithmeticError(
            "Zero base logarithm."
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
# SYNTHETIC SOURCE
# ============================================================================

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
    q1: int,
    q3: int,
    k: int,
    modulus: int,
) -> int:

    return int(
        (
            q1
            -
            (
                2
                * pow(
                    3,
                    k,
                    modulus,
                )
                * q3
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
        "EXPERIMENT 252 — EXACT SYNTHETIC UNIT-SCALE / "
        "PROJECTIVE-HENSEL UNIVERSALITY"
    )
    print("=" * 78)

    failures = []

    selected_primes = []

    # --------------------------------------------------------------
    # Find first MAX_PRIMES nonexceptional primes.
    # --------------------------------------------------------------

    for ell in range(
        5,
        PRIME_LIMIT + 1,
        2,
    ):

        if not is_prime(ell):
            continue

        d = multiplicative_order(
            3,
            ell,
            ell - 1,
        )

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
            continue

        selected_primes.append(
            (
                ell,
                d,
            )
        )

        if len(selected_primes) >= MAX_PRIMES:
            break

    case_count = 0
    fully_exact = 0

    print()
    print("=" * 78)
    print("1. SELECTED NONEXCEPTIONAL PRIMES")
    print("=" * 78)

    print(
        f"  primes="
        f"{[ell for ell, _ in selected_primes]}"
    )

    # --------------------------------------------------------------
    # Synthetic unit-scaled cases.
    # --------------------------------------------------------------

    for ell, d in selected_primes:

        modulus_full = ell ** MAX_LEVEL

        for k1 in range(
            min(d, 4)
        ):

            rho = orbit_mod(
                k1,
                modulus_full,
            )

            for scale in UNIT_SCALES:

                q3 = scale

                if gcd(
                    q3,
                    ell,
                ) != 1:

                    continue

                # Exact source relation modulo ell^MAX_LEVEL.
                q1 = (
                    rho
                    * q3
                ) % modulus_full

                case_count += 1

                local_ok = True

                # --------------------------------------------------
                # Base class
                # --------------------------------------------------

                base_matches = [
                    k
                    for k in range(d)
                    if (
                        orbit_mod(
                            k,
                            ell,
                        )
                        ==
                        (
                            q1
                            * pow(
                                q3,
                                -1,
                                ell,
                            )
                        )
                        % ell
                    )
                ]

                if base_matches != [k1]:

                    failures.append(
                        (
                            ell,
                            k1,
                            scale,
                            "base_class",
                            base_matches,
                        )
                    )

                    continue

                # --------------------------------------------------
                # Recursive Hensel lifting
                # --------------------------------------------------

                k = k1
                current_order = d

                exponents = [k1]
                digits = []
                principal = [0]

                for level in range(
                    1,
                    MAX_LEVEL,
                ):

                    current_modulus = (
                        ell ** level
                    )

                    next_modulus = (
                        ell ** (
                            level + 1
                        )
                    )

                    expected_order = (
                        d
                        * ell ** (
                            level - 1
                        )
                    )

                    if (
                        current_order
                        != expected_order
                    ):

                        failures.append(
                            (
                                ell,
                                k1,
                                scale,
                                level,
                                "order_growth",
                            )
                        )

                        local_ok = False
                        break

                    # --------------------------------------------------
                    # D_r
                    # --------------------------------------------------

                    increment = (
                        pow(
                            3,
                            current_order,
                            next_modulus,
                        )
                        - 1
                    )

                    if (
                        increment
                        % current_modulus
                        != 0
                    ):

                        failures.append(
                            (
                                ell,
                                k1,
                                scale,
                                level,
                                "D_integrality",
                            )
                        )

                        local_ok = False
                        break

                    D = (
                        increment
                        // current_modulus
                    ) % ell

                    if D == 0:

                        failures.append(
                            (
                                ell,
                                k1,
                                scale,
                                level,
                                "zero_D",
                            )
                        )

                        local_ok = False
                        break

                    # --------------------------------------------------
                    # Affine residual
                    # --------------------------------------------------

                    normalized = []

                    for t in range(ell):

                        candidate = (
                            k
                            + t * current_order
                        )

                        F = residual_mod(
                            q1,
                            q3,
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
                                    k1,
                                    scale,
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
                                k1,
                                scale,
                                level,
                                "non_affine",
                            )
                        )

                        local_ok = False
                        break

                    measured_B = next(
                        iter(slope_set)
                    )

                    rho_mod = (
                        (
                            q1
                            * pow(
                                q3,
                                -1,
                                ell,
                            )
                        )
                        % ell
                    )

                    predicted_B = (
                        -(
                            q3
                            % ell
                        )
                        * rho_mod
                        * D
                    ) % ell

                    if (
                        measured_B
                        != predicted_B
                    ):

                        failures.append(
                            (
                                ell,
                                k1,
                                scale,
                                level,
                                "slope",
                                measured_B,
                                predicted_B,
                            )
                        )

                        local_ok = False
                        break

                    # --------------------------------------------------
                    # Recover digit
                    # --------------------------------------------------

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

                    target = (
                        q1
                        * pow(
                            q3,
                            -1,
                            next_modulus,
                        )
                    ) % next_modulus

                    recovered = (
                        orbit_mod(
                            next_k,
                            next_modulus,
                        )
                    )

                    if recovered != target:

                        failures.append(
                            (
                                ell,
                                k1,
                                scale,
                                level,
                                "lift_alignment",
                                recovered,
                                target,
                            )
                        )

                        local_ok = False
                        break

                    digits.append(t)
                    exponents.append(
                        next_k
                    )

                    principal.append(
                        (
                            next_k
                            - k1
                        )
                        // d
                    )

                    phi_next = (
                        (ell - 1)
                        * ell ** level
                    )

                    current_order = (
                        multiplicative_order(
                            3,
                            next_modulus,
                            phi_next,
                        )
                    )

                    k = next_k

                if not local_ok:
                    continue

                # --------------------------------------------------
                # Principal digit identity
                # --------------------------------------------------

                for r, m in enumerate(
                    principal
                ):

                    expected_digits = []

                    value = m

                    for _ in range(r):

                        expected_digits.append(
                            value % ell
                        )

                        value //= ell

                    if (
                        expected_digits
                        !=
                        digits[:r]
                    ):

                        failures.append(
                            (
                                ell,
                                k1,
                                scale,
                                r,
                                "principal_digits",
                                expected_digits,
                                digits[:r],
                            )
                        )

                        local_ok = False
                        break

                if not local_ok:
                    continue

                # --------------------------------------------------
                # Exponential reconstruction
                # --------------------------------------------------

                for level, m in enumerate(
                    principal
                ):

                    if level == 0:
                        continue

                    modulus = ell ** level

                    ratio = (
                        q1
                        * pow(
                            q3,
                            -1,
                            modulus,
                        )
                    ) % modulus

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

                    if reconstructed != ratio:

                        failures.append(
                            (
                                ell,
                                k1,
                                scale,
                                level,
                                "exponential",
                            )
                        )

                        local_ok = False
                        break

                if not local_ok:
                    continue

                # --------------------------------------------------
                # Logarithmic reconstruction
                # --------------------------------------------------

                for level in range(
                    2,
                    MAX_LEVEL + 1,
                ):

                    modulus = ell ** level

                    ratio = (
                        q1
                        * pow(
                            q3,
                            -1,
                            modulus,
                        )
                    ) % modulus

                    base_factor = (
                        2
                        * pow(
                            3,
                            k1,
                            modulus,
                        )
                    ) % modulus

                    source_unit = (
                        ratio
                        * pow(
                            base_factor,
                            -1,
                            modulus,
                        )
                    ) % modulus

                    base_unit = (
                        pow(
                            3,
                            d,
                            modulus,
                        )
                    )

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
                                k1,
                                scale,
                                level,
                                "logarithm",
                                type(exc).__name__,
                                str(exc),
                            )
                        )

                        local_ok = False
                        break

                    expected_m = principal[
                        level - 1
                    ]

                    coordinate_modulus = (
                        ell ** (
                            level - 1
                        )
                    )

                    if (
                        log_m
                        % coordinate_modulus
                        !=
                        expected_m
                        % coordinate_modulus
                    ):

                        failures.append(
                            (
                                ell,
                                k1,
                                scale,
                                level,
                                "log_coordinate",
                                log_m,
                                expected_m,
                            )
                        )

                        local_ok = False
                        break

                if local_ok:
                    fully_exact += 1

    # ------------------------------------------------------------------
    # SUMMARY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. UNIT-SCALE UNIVERSALITY SUMMARY")
    print("=" * 78)

    print(
        f"  cases_tested={case_count}"
    )

    print(
        f"  fully_exact_cases={fully_exact}"
    )

    print(
        f"  failures={len(failures)}"
    )

    # ------------------------------------------------------------------
    # FAILURE SUMMARY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. FAILURE SUMMARY")
    print("=" * 78)

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

    else:

        print(
            "  none"
        )

    # ------------------------------------------------------------------
    # STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
Experiment 251 tested the normalized chart

    q3 = 1.

Experiment 252 restores the general source scale

    q3 = unit,

while preserving the same projective ratio

    rho = q1/q3.

The local residual is

    F(k)
      = q1 - 2*3^k*q3,

so the derivative must contain q3:

    B
      = -q3*rho*D mod ell.

The experiment therefore separates:

    projective information:
        q1/q3,

from

    source scaling:
        q3.

If all unit-scaled cases pass, the local Hensel/logarithmic theorem is
genuinely projective and is not an artifact of choosing q3=1.
"""
    )

    # ------------------------------------------------------------------
    # FINAL
    # ------------------------------------------------------------------

    final_ok = (
        len(failures) == 0
        and
        case_count > 0
    )

    print()
    print("=" * 78)
    print("5. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  cases_tested={case_count}"
    )

    print(
        f"  fully_exact_cases={fully_exact}"
    )

    print(
        f"  projective_unit_scale_invariance_exact="
        f"{final_ok}"
    )

    print(
        f"  q3_derivative_factor_exact="
        f"{final_ok}"
    )

    print(
        f"  exponential_reconstruction_exact="
        f"{final_ok}"
    )

    print(
        f"  logarithmic_reconstruction_exact="
        f"{final_ok}"
    )

    print(
        f"  failures={len(failures)}"
    )

    print(
        "  exceptional_prime_11_excluded=True"
    )

    print(
        "  general_n_pq_theorem_established=False"
    )

    print(
        "  reason=single_actual_n6_source_family"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 252 COMPLETE")


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

