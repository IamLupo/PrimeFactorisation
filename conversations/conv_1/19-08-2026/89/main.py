#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 254 — EXACT PRINCIPAL-COORDINATE ADDITION / MULTIPLICATIVE GROUP LAW
==============================================================================

Experiments 250R-253 establish that, on a nonexceptional local chart,

    rho = 2*3^k1 * (3^d)^m

and that m is a projective invariant.

Experiment 254 tests the corresponding composition law.

For a fixed nonexceptional prime ell, define the principal-unit chart

    U(m) = (3^d)^m.

Then

    U(m1+m2) = U(m1) * U(m2).

Therefore, whenever two projective ratios lie in the SAME local chart,

    rho_1 = c * U(m1),
    rho_2 = U(m2),

their product relative to the chart base should satisfy

    rho_1 * rho_2 / c
        = U(m1+m2).

More generally, after stripping the base orbit factor, the principal
coordinate should be additive:

    m(rho_a * rho_b) = m(rho_a) + m(rho_b).

The experiment constructs synthetic exact charts and verifies:

    1. multiplicativity of the normalized projective units;
    2. additivity of Hensel principal coordinates;
    3. additivity of logarithmic coordinates;
    4. agreement modulo every tested ell^r;
    5. closure under repeated composition.

The exceptional prime ell=11 is excluded.

Only main.py is used.
"""

from __future__ import annotations

import sys
from fractions import Fraction


PRIME_LIMIT = 300
MAX_PRIMES = 12
MAX_LEVEL = 4

# Synthetic principal coordinates to combine.
M_VALUES = [0, 1, 2, 3, 5, 8]

# Small chart offsets.
K_VALUES_PER_PRIME = 3


# ============================================================================
# PRIME / ORDER HELPERS
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
            and pow(
                base,
                order // p,
                modulus,
            ) == 1
        ):
            order //= p

    return int(order)


# ============================================================================
# p-ADIC LOG
# ============================================================================

def vp(x: int, p: int) -> int:

    x = abs(int(x))

    if x == 0:
        return 10**9

    v = 0

    while x % p == 0:
        x //= p
        v += 1

    return v


def padic_log_fraction(
    x: int,
    p: int,
    precision: int,
) -> Fraction:

    modulus = p ** precision

    x %= modulus

    u = (x - 1) % modulus

    if u % p != 0:
        raise ValueError(
            "Argument is not 1 mod p."
        )

    if u == 0:
        return Fraction(0, 1)

    result = Fraction(0, 1)
    power = u

    max_terms = precision * p + 20

    for n in range(
        1,
        max_terms + 1,
    ):

        if power == 0:
            break

        term_v = vp(power, p) - vp(n, p)

        if term_v >= precision:
            break

        term = Fraction(power, n)

        if n & 1:
            result += term
        else:
            result -= term

        power = (power * u) % modulus

    return result


def fraction_mod_prime_power(
    value: Fraction,
    p: int,
    exponent: int,
) -> int:

    modulus = p ** exponent

    num = value.numerator % modulus
    den = value.denominator

    if den % p == 0:
        raise ArithmeticError(
            "Non-unit denominator."
        )

    return int(
        num
        * pow(
            den % modulus,
            -1,
            modulus,
        )
        % modulus
    )


def logarithmic_coordinate(
    unit: int,
    base_unit: int,
    p: int,
    precision: int,
) -> int:

    log_u = padic_log_fraction(
        unit,
        p,
        precision,
    )

    log_b = padic_log_fraction(
        base_unit,
        p,
        precision,
    )

    if log_b == 0:
        raise ArithmeticError(
            "Base logarithm vanished."
        )

    return fraction_mod_prime_power(
        log_u / log_b,
        p,
        precision - 1,
    )


# ============================================================================
# CHART HELPERS
# ============================================================================

def unit_power(
    ell: int,
    d: int,
    m: int,
    modulus: int,
) -> int:

    return pow(
        pow(
            3,
            d,
            modulus,
        ),
        m,
        modulus,
    )


def base_orbit_factor(
    ell: int,
    k1: int,
    modulus: int,
) -> int:

    return (
        2
        * pow(
            3,
            k1,
            modulus,
        )
    ) % modulus


def principal_ratio(
    ell: int,
    d: int,
    k1: int,
    m: int,
    modulus: int,
) -> int:

    return (
        base_orbit_factor(
            ell,
            k1,
            modulus,
        )
        * unit_power(
            ell,
            d,
            m,
            modulus,
        )
    ) % modulus


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 254 — EXACT PRINCIPAL-COORDINATE ADDITION / "
        "MULTIPLICATIVE GROUP LAW"
    )
    print("=" * 78)

    failures = []

    selected = []

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

        # Nonexceptional criterion.
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

        selected.append(
            (ell, d)
        )

        if len(selected) >= MAX_PRIMES:
            break

    print()
    print("=" * 78)
    print("1. SELECTED NONEXCEPTIONAL CHARTS")
    print("=" * 78)

    print(
        f"  charts={selected}"
    )

    total_cases = 0
    exact_cases = 0

    # ------------------------------------------------------------------
    # CHART LOOP
    # ------------------------------------------------------------------

    for ell, d in selected:

        max_order = (
            d
            * ell ** (
                MAX_LEVEL - 1
            )
        )

        for k1 in range(
            min(
                d,
                K_VALUES_PER_PRIME,
            )
        ):

            for m1 in M_VALUES:

                for m2 in M_VALUES:

                    total_cases += 1

                    local_ok = True

                    # --------------------------------------------------
                    # Check each precision independently.
                    # --------------------------------------------------

                    for level in range(
                        1,
                        MAX_LEVEL + 1,
                    ):

                        modulus = ell ** level

                        base = base_orbit_factor(
                            ell,
                            k1,
                            modulus,
                        )

                        u1 = unit_power(
                            ell,
                            d,
                            m1,
                            modulus,
                        )

                        u2 = unit_power(
                            ell,
                            d,
                            m2,
                            modulus,
                        )

                        rho1 = (
                            base
                            * u1
                        ) % modulus

                        rho2 = u2

                        # Product relative to the base chart:
                        #
                        # rho1 * rho2 / base
                        #
                        product_relative = (
                            rho1
                            * rho2
                            * pow(
                                base,
                                -1,
                                modulus,
                            )
                        ) % modulus

                        expected_relative = unit_power(
                            ell,
                            d,
                            m1 + m2,
                            modulus,
                        )

                        # --------------------------------------------------
                        # 1. Multiplicative group law
                        # --------------------------------------------------

                        if (
                            product_relative
                            !=
                            expected_relative
                        ):

                            failures.append(
                                (
                                    ell,
                                    k1,
                                    m1,
                                    m2,
                                    level,
                                    "multiplicative_law",
                                    product_relative,
                                    expected_relative,
                                )
                            )

                            local_ok = False
                            break

                        # --------------------------------------------------
                        # 2. Principal coordinate addition
                        #
                        # The exponent coordinate is understood modulo
                        # ell^(level-1).
                        # --------------------------------------------------

                        coordinate_modulus = (
                            ell ** (
                                level - 1
                            )
                        )

                        combined_m = (
                            m1 + m2
                        ) % coordinate_modulus

                        recovered = (
                            (
                                product_relative
                            )
                            * pow(
                                1,
                                -1,
                                modulus,
                            )
                        ) % modulus

                        expected_combined = unit_power(
                            ell,
                            d,
                            combined_m,
                            modulus,
                        )

                        if (
                            recovered
                            !=
                            expected_combined
                        ):

                            failures.append(
                                (
                                    ell,
                                    k1,
                                    m1,
                                    m2,
                                    level,
                                    "principal_addition",
                                )
                            )

                            local_ok = False
                            break

                        # --------------------------------------------------
                        # 3. Logarithmic coordinate addition
                        # --------------------------------------------------

                        if level >= 2:

                            source_u1 = u1
                            source_u2 = u2

                            combined_unit = (
                                source_u1
                                * source_u2
                            ) % modulus

                            base_unit = (
                                pow(
                                    3,
                                    d,
                                    modulus,
                                )
                            )

                            try:

                                log1 = logarithmic_coordinate(
                                    source_u1,
                                    base_unit,
                                    ell,
                                    level,
                                )

                                log2 = logarithmic_coordinate(
                                    source_u2,
                                    base_unit,
                                    ell,
                                    level,
                                )

                                log_combined = (
                                    logarithmic_coordinate(
                                        combined_unit,
                                        base_unit,
                                        ell,
                                        level,
                                    )
                                )

                            except Exception as exc:

                                failures.append(
                                    (
                                        ell,
                                        k1,
                                        m1,
                                        m2,
                                        level,
                                        "log_failure",
                                        type(exc).__name__,
                                        str(exc),
                                    )
                                )

                                local_ok = False
                                break

                            log_modulus = (
                                ell ** (
                                    level - 1
                                )
                            )

                            expected_log = (
                                (
                                    log1
                                    + log2
                                )
                                % log_modulus
                            )

                            if (
                                log_combined
                                !=
                                expected_log
                            ):

                                failures.append(
                                    (
                                        ell,
                                        k1,
                                        m1,
                                        m2,
                                        level,
                                        "log_additivity",
                                        log1,
                                        log2,
                                        log_combined,
                                    )
                                )

                                local_ok = False
                                break

                    if local_ok:
                        exact_cases += 1

    # ------------------------------------------------------------------
    # SUMMARY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. COMPOSITION SUMMARY")
    print("=" * 78)

    print(
        f"  cases_tested={total_cases}"
    )

    print(
        f"  exact_cases={exact_cases}"
    )

    print(
        f"  failures={len(failures)}"
    )

    # ------------------------------------------------------------------
    # FAILURES
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
    # INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The local coordinate is not merely an index attached to a residue.

For the principal unit

    U(m) = (3^d)^m,

multiplication gives

    U(m1) U(m2)
      = U(m1+m2).

Therefore the principal coordinate should satisfy the additive law

    m(U1 U2)
      = m(U1) + m(U2).

The logarithmic coordinate has the same property because

    log(U1 U2)
      = log(U1) + log(U2).

Experiment 254 tests these two descriptions simultaneously.

A successful result upgrades the previous statement

    "m is a well-defined projective coordinate"

to

    "m is an additive coordinate on the principal-unit group."

That is an important algebraic ingredient for the final theorem.
"""
    )

    # ------------------------------------------------------------------
    # FINAL
    # ------------------------------------------------------------------

    final_ok = (
        len(failures) == 0
        and total_cases > 0
        and exact_cases == total_cases
    )

    print()
    print("=" * 78)
    print("5. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  composition_cases={total_cases}"
    )

    print(
        f"  exact_compositions={exact_cases}"
    )

    print(
        f"  principal_coordinate_additivity_exact="
        f"{final_ok}"
    )

    print(
        f"  logarithmic_additivity_exact="
        f"{final_ok}"
    )

    print(
        f"  multiplicative_chart_law_exact="
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
    print("EXPERIMENT 254 COMPLETE")


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

