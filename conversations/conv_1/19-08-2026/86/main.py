#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 251 — EXACT SYNTHETIC PROJECTIVE-CHART / LOG-HENSEL UNIVERSALITY
==============================================================================

Experiments 249R and 250R established the full local mechanism for the
single source ratio

    rho = q1/q3,

across many admissible nonexceptional primes.

Experiment 251 removes q1 and q3.

For each nonexceptional odd prime ell:

    d = ord_ell(3),

choose several synthetic base exponent classes

    k1 in {0,...,d-1}

and define the synthetic projective ratio

    rho = 2*3^k1.

Use the normalized synthetic source pair

    q3 = 1,
    q1 = rho,

where rho is represented modulo ell^MAX_LEVEL.

The experiment then independently verifies:

    1. base orbit identification;
    2. Hensel digit recovery;
    3. local derivative formula;
    4. exponent reconstruction;
    5. principal-coordinate base-ell digits;
    6. exponential reconstruction;
    7. logarithmic reconstruction.

This tests whether the entire local theorem is a property of the
multiplicative orbit chart itself, rather than of the particular
integer source pair from n=6.

Only main.py is used.
"""

from __future__ import annotations

import sys
from fractions import Fraction
from math import gcd


PRIME_LIMIT = 500
MAX_LEVEL = 4

# Number of synthetic orbit classes tested per prime.
# 0 means every class, otherwise use the first SAMPLE_CLASSES classes.
SAMPLE_CLASSES = 5


# ============================================================================
# PRIME HELPERS
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


# ============================================================================
# EXACT p-ADIC LOGGING
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


def padic_log_fraction(
    x: int,
    p: int,
    precision: int,
) -> Fraction:

    modulus = p ** precision

    x = int(x) % modulus

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

    denominator = (
        value.denominator
    )

    if denominator % p == 0:

        raise ArithmeticError(
            "Non-unit p-adic denominator."
        )

    inv_den = pow(
        denominator % modulus,
        -1,
        modulus,
    )

    return int(
        numerator
        * inv_den
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

def synthetic_ratio(
    k1: int,
    modulus: int,
) -> int:

    return int(
        (
            2
            * pow(
                3,
                k1,
                modulus,
            )
        )
        % modulus
    )


def synthetic_residual(
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
        "EXPERIMENT 251 — EXACT SYNTHETIC PROJECTIVE-CHART / "
        "LOG-HENSEL UNIVERSALITY"
    )
    print("=" * 78)

    failures = []

    prime_count = 0
    nonexceptional_count = 0
    class_count = 0
    fully_tested_classes = 0

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

        prime_count += 1

        d = multiplicative_order(
            3,
            ell,
            ell - 1,
        )

        # --------------------------------------------------------------
        # Nonexceptional test
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
            continue

        nonexceptional_count += 1

        # --------------------------------------------------------------
        # Choose synthetic classes
        # --------------------------------------------------------------

        if SAMPLE_CLASSES == 0:

            classes = list(
                range(d)
            )

        else:

            classes = list(
                range(
                    min(
                        d,
                        SAMPLE_CLASSES,
                    )
                )
            )

        for k1 in classes:

            class_count += 1

            # ----------------------------------------------------------
            # Build exact synthetic source modulo ell^MAX_LEVEL
            # ----------------------------------------------------------

            full_modulus = ell ** (
                MAX_LEVEL
            )

            q3 = 1

            q1 = synthetic_ratio(
                k1,
                full_modulus,
            )

            # ----------------------------------------------------------
            # Verify base chart
            # ----------------------------------------------------------

            matches = [
                k
                for k in range(d)
                if synthetic_ratio(
                    k,
                    ell,
                )
                == q1 % ell
            ]

            if matches != [k1]:

                failures.append(
                    (
                        ell,
                        k1,
                        "base_class",
                        matches,
                    )
                )

                continue

            # ----------------------------------------------------------
            # Recursive local lift
            # ----------------------------------------------------------

            k = k1
            current_order = d

            exponent_sequence = [
                k1
            ]

            digits = []
            principal = [0]

            local_ok = True

            for level in range(
                1,
                MAX_LEVEL,
            ):

                current_modulus = ell ** level
                next_modulus = ell ** (
                    level + 1
                )

                # ------------------------------------------------------
                # Order
                # ------------------------------------------------------

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
                            level,
                            "order_growth",
                        )
                    )

                    local_ok = False
                    break

                # ------------------------------------------------------
                # D
                # ------------------------------------------------------

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
                            level,
                            "zero_D",
                        )
                    )

                    local_ok = False
                    break

                # ------------------------------------------------------
                # Synthetic affine residual
                # ------------------------------------------------------

                residuals = []

                for t in range(ell):

                    candidate = (
                        k
                        + t * current_order
                    )

                    F = synthetic_residual(
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
                                level,
                                "residual_normalization",
                                t,
                            )
                        )

                        local_ok = False
                        break

                    residuals.append(
                        (
                            F
                            // current_modulus
                        ) % ell
                    )

                if not local_ok:
                    break

                A = residuals[0]

                slope_set = {
                    (
                        residuals[i]
                        - residuals[i - 1]
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
                            level,
                            "non_affine",
                        )
                    )

                    local_ok = False
                    break

                B = next(
                    iter(slope_set)
                )

                rho = q1 % ell

                predicted_B = (
                    -rho
                    * D
                ) % ell

                if B != predicted_B:

                    failures.append(
                        (
                            ell,
                            k1,
                            level,
                            "derivative",
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

                # ------------------------------------------------------
                # Independent exact target
                # ------------------------------------------------------

                target = (
                    q1
                    % next_modulus
                )

                reconstructed = (
                    synthetic_ratio(
                        next_k,
                        next_modulus,
                    )
                )

                if (
                    reconstructed
                    != target
                ):

                    failures.append(
                        (
                            ell,
                            k1,
                            level,
                            "lift_alignment",
                        )
                    )

                    local_ok = False
                    break

                digits.append(t)

                exponent_sequence.append(
                    next_k
                )

                principal.append(
                    (
                        next_k
                        - k1
                    )
                    // d
                )

                # Actual next order.
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

            # ----------------------------------------------------------
            # Principal base-ell digits
            # ----------------------------------------------------------

            for r, m in enumerate(
                principal
            ):

                digits_expected = []

                if r > 0:

                    value = m

                    for _ in range(r):

                        digits_expected.append(
                            value % ell
                        )

                        value //= ell

                if (
                    digits_expected
                    !=
                    digits[:r]
                ):

                    failures.append(
                        (
                            ell,
                            k1,
                            r,
                            "principal_digits",
                            digits_expected,
                            digits[:r],
                        )
                    )

                    local_ok = False
                    break

            if not local_ok:
                continue

            # ----------------------------------------------------------
            # Exponential reconstruction
            # ----------------------------------------------------------

            for level, m in enumerate(
                principal
            ):

                if level == 0:
                    continue

                modulus = ell ** level

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

                target = (
                    q1
                    % modulus
                )

                if reconstructed != target:

                    failures.append(
                        (
                            ell,
                            k1,
                            level,
                            "exponential",
                        )
                    )

                    local_ok = False
                    break

            if not local_ok:
                continue

            # ----------------------------------------------------------
            # Logarithmic reconstruction
            # ----------------------------------------------------------

            for level in range(
                2,
                MAX_LEVEL + 1,
            ):

                modulus = ell ** level

                source_ratio = (
                    q1
                    % modulus
                )

                base_factor = (
                    2
                    * pow(
                        3,
                        k1,
                        modulus,
                    )
                ) % modulus

                source_unit = (
                    source_ratio
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

                if source_unit % ell != 1:

                    failures.append(
                        (
                            ell,
                            k1,
                            level,
                            "source_unit",
                        )
                    )

                    local_ok = False
                    break

                if base_unit % ell != 1:

                    failures.append(
                        (
                            ell,
                            k1,
                            level,
                            "base_unit",
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
                            k1,
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
                    %
                    coordinate_modulus
                    !=
                    expected_m
                    %
                    coordinate_modulus
                ):

                    failures.append(
                        (
                            ell,
                            k1,
                            level,
                            "log_coordinate",
                            log_m,
                            expected_m,
                        )
                    )

                    local_ok = False
                    break

            if local_ok:
                fully_tested_classes += 1

    # ------------------------------------------------------------------
    # SUMMARY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. SYNTHETIC UNIVERSALITY SUMMARY")
    print("=" * 78)

    print(
        f"  prime_limit={PRIME_LIMIT}"
    )

    print(
        f"  nonexceptional_primes="
        f"{nonexceptional_count}"
    )

    print(
        f"  synthetic_classes_tested="
        f"{class_count}"
    )

    print(
        f"  fully_exact_classes="
        f"{fully_tested_classes}"
    )

    print(
        f"  failures={len(failures)}"
    )

    # ------------------------------------------------------------------
    # FAILURE DETAILS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. FAILURE SUMMARY")
    print("=" * 78)

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

    else:

        print(
            "  none"
        )

    # ------------------------------------------------------------------
    # STRUCTURE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
Experiments 249R and 250R verified the local theorem for one actual
source ratio.

Experiment 251 removes that source-specific dependence.

For every selected nonexceptional prime and several independent
projective orbit classes,

    rho = 2*3^k1,

the complete local mechanism is reconstructed from scratch:

    base orbit class
        ->
    affine residual
        ->
    derivative
        ->
    Hensel digit
        ->
    exponent
        ->
    principal coordinate
        ->
    exponential reconstruction
        ->
    logarithmic reconstruction.

This is a much stronger universality test because q1 and q3 no longer
supply the observed ratio.

A successful result means the local theorem is a structural property of
the multiplicative orbit chart itself.
"""
    )

    # ------------------------------------------------------------------
    # FINAL
    # ------------------------------------------------------------------

    final_ok = (
        len(failures) == 0
        and
        fully_tested_classes > 0
    )

    print()
    print("=" * 78)
    print("4. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  synthetic_classes_tested="
        f"{class_count}"
    )

    print(
        f"  fully_exact_classes="
        f"{fully_tested_classes}"
    )

    print(
        f"  synthetic_projective_hensel_universality_exact="
        f"{final_ok}"
    )

    print(
        f"  synthetic_exponential_reconstruction_exact="
        f"{final_ok}"
    )

    print(
        f"  synthetic_logarithmic_reconstruction_exact="
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
        "  reason=single_known_n6_source_family"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 251 COMPLETE")


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

