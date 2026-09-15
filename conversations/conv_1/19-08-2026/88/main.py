#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 253 — EXACT PROJECTIVE REPRESENTATIVE / HIGHER-ORDER PERTURBATION
AUDIT
==============================================================================

Experiment 251:
    synthetic projective charts with q3 = 1.

Experiment 252:
    arbitrary unit q3 scales with the same projective ratio.

Experiment 253:
    changes the actual integer representatives q1 and q3 while preserving
    the same projective ratio modulo ell^MAX_LEVEL.

For a synthetic chart with

    rho = 2*3^k1,

choose a unit q3 and define

    q1 = rho*q3 + c*ell^MAX_LEVEL.

Because q3 is an ell-adic unit,

    q1/q3 = rho (mod ell^MAX_LEVEL).

Thus the perturbed pair represents the same projective class through the
entire tested precision.

The experiment compares the perturbed source against the canonical source:

    q3_base = q3
    q1_base = rho*q3

versus

    q3_pert = q3 + a*ell^MAX_LEVEL
    q1_pert = rho*q3_pert + b*ell^MAX_LEVEL.

It checks whether all of the following remain identical:

    * base exponent k1;
    * Hensel digits;
    * exponent sequence;
    * principal coordinates;
    * local slopes;
    * logarithmic coordinate;
    * reconstructed projective ratio.

Only main.py is used.
"""

from __future__ import annotations

import sys
from fractions import Fraction
from math import gcd


PRIME_LIMIT = 200
MAX_LEVEL = 4

MAX_PRIMES = 12

K_CLASSES_PER_PRIME = 4

UNIT_SCALES = [1, 2, 5]

PERTURB_A = [0, 1, 2, 3]
PERTURB_B = [0, 1, 2, 5]


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
# EXACT p-ADIC LOG
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
# ORBIT / RESIDUAL HELPERS
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


def ratio_mod(
    q1: int,
    q3: int,
    modulus: int,
) -> int:

    if gcd(
        q3,
        modulus,
    ) != 1:

        raise ArithmeticError(
            "q3 is not invertible."
        )

    return int(
        (
            (q1 % modulus)
            * pow(
                q3 % modulus,
                -1,
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
                    int(k),
                    modulus,
                )
                * q3
            )
        )
        % modulus
    )


# ============================================================================
# SOURCE RECONSTRUCTION
# ============================================================================

def reconstruct_chart(
    q1: int,
    q3: int,
    ell: int,
    max_level: int,
) -> dict:

    d = multiplicative_order(
        3,
        ell,
        ell - 1,
    )

    rho0 = ratio_mod(
        q1,
        q3,
        ell,
    )

    base_matches = [
        k
        for k in range(d)
        if orbit_mod(
            k,
            ell,
        ) == rho0
    ]

    if len(base_matches) != 1:

        raise ArithmeticError(
            "Base projective class is not uniquely represented."
        )

    k1 = base_matches[0]

    k = k1
    current_order = d

    exponents = [k1]
    digits = []
    slopes = []
    principal = [0]

    for level in range(
        1,
        max_level,
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

            raise ArithmeticError(
                "Unexpected order growth."
            )

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

            raise ArithmeticError(
                "Order increment is not divisible by the expected power."
            )

        D = (
            increment
            // current_modulus
        ) % ell

        if D == 0:

            raise ArithmeticError(
                "Unexpected exceptional level."
            )

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

                raise ArithmeticError(
                    "Residual is not divisible at the expected level."
                )

            normalized.append(
                (
                    F
                    // current_modulus
                ) % ell
            )

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

            raise ArithmeticError(
                "Residual is not affine."
            )

        B = next(
            iter(slope_set)
        )

        rho = ratio_mod(
            q1,
            q3,
            ell,
        )

        predicted_B = (
            -(q3 % ell)
            * rho
            * D
        ) % ell

        if B != predicted_B:

            raise ArithmeticError(
                "Derivative formula failed."
            )

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
            q1,
            q3,
            next_modulus,
        )

        if (
            orbit_mod(
                next_k,
                next_modulus,
            )
            != target
        ):

            raise ArithmeticError(
                "Lift alignment failed."
            )

        digits.append(t)
        slopes.append(B)

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

        current_order = multiplicative_order(
            3,
            next_modulus,
            phi_next,
        )

        k = next_k

    # --------------------------------------------------------------
    # Exponential reconstruction
    # --------------------------------------------------------------

    for level, m in enumerate(
        principal
    ):

        if level == 0:
            continue

        modulus = ell ** level

        rho = ratio_mod(
            q1,
            q3,
            modulus,
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

        if reconstructed != rho:

            raise ArithmeticError(
                "Exponential reconstruction failed."
            )

    # --------------------------------------------------------------
    # Logarithmic coordinate
    # --------------------------------------------------------------

    logarithmic = []

    for level in range(
        2,
        max_level + 1,
    ):

        modulus = ell ** level

        rho = ratio_mod(
            q1,
            q3,
            modulus,
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
            rho
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

        log_m = logarithmic_coordinate(
            source_unit,
            base_unit,
            ell,
            level,
        )

        logarithmic.append(
            log_m
        )

        expected = principal[
            level - 1
        ]

        coord_modulus = ell ** (
            level - 1
        )

        if (
            log_m
            % coord_modulus
            !=
            expected
            % coord_modulus
        ):

            raise ArithmeticError(
                "Logarithmic coordinate mismatch."
            )

    return {
        "k1": k1,
        "digits": digits,
        "slopes": slopes,
        "exponents": exponents,
        "principal": principal,
        "logarithmic": logarithmic,
    }


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 253 — EXACT PROJECTIVE REPRESENTATIVE / "
        "HIGHER-ORDER PERTURBATION AUDIT"
    )
    print("=" * 78)

    failures = []
    cases = 0
    exact = 0

    selected_primes = []

    # --------------------------------------------------------------
    # Select nonexceptional primes.
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
            ell
        )

        if len(selected_primes) >= MAX_PRIMES:
            break

    print()
    print("=" * 78)
    print("1. SELECTED NONEXCEPTIONAL PRIMES")
    print("=" * 78)

    print(
        f"  primes={selected_primes}"
    )

    # --------------------------------------------------------------
    # Test same projective class with different representatives.
    # --------------------------------------------------------------

    for ell in selected_primes:

        d = multiplicative_order(
            3,
            ell,
            ell - 1,
        )

        full_modulus = ell ** MAX_LEVEL

        for k1 in range(
            min(
                d,
                K_CLASSES_PER_PRIME,
            )
        ):

            rho = orbit_mod(
                k1,
                full_modulus,
            )

            for scale in UNIT_SCALES:

                if gcd(
                    scale,
                    ell,
                ) != 1:
                    continue

                # Canonical representative.
                q3_base = scale
                q1_base = (
                    rho
                    * q3_base
                ) % full_modulus

                try:

                    base_result = reconstruct_chart(
                        q1_base,
                        q3_base,
                        ell,
                        MAX_LEVEL,
                    )

                except Exception as exc:

                    failures.append(
                        (
                            ell,
                            k1,
                            scale,
                            "base_reconstruction",
                            type(exc).__name__,
                            str(exc),
                        )
                    )

                    continue

                for a in PERTURB_A:

                    q3_pert = (
                        q3_base
                        + a
                        * full_modulus
                    )

                    if gcd(
                        q3_pert,
                        ell,
                    ) != 1:

                        continue

                    for b in PERTURB_B:

                        q1_pert = (
                            rho
                            * q3_pert
                            + b
                            * full_modulus
                        )

                        cases += 1

                        # --------------------------------------------------
                        # Verify projective equivalence directly.
                        # --------------------------------------------------

                        for level in range(
                            1,
                            MAX_LEVEL + 1,
                        ):

                            modulus = ell ** level

                            ratio_base = (
                                ratio_mod(
                                    q1_base,
                                    q3_base,
                                    modulus,
                                )
                            )

                            ratio_pert = (
                                ratio_mod(
                                    q1_pert,
                                    q3_pert,
                                    modulus,
                                )
                            )

                            if (
                                ratio_base
                                !=
                                ratio_pert
                            ):

                                failures.append(
                                    (
                                        ell,
                                        k1,
                                        scale,
                                        a,
                                        b,
                                        level,
                                        "projective_ratio_changed",
                                        ratio_base,
                                        ratio_pert,
                                    )
                                )

                                break

                        else:

                            try:

                                pert_result = reconstruct_chart(
                                    q1_pert,
                                    q3_pert,
                                    ell,
                                    MAX_LEVEL,
                                )

                            except Exception as exc:

                                failures.append(
                                    (
                                        ell,
                                        k1,
                                        scale,
                                        a,
                                        b,
                                        "pert_reconstruction",
                                        type(exc).__name__,
                                        str(exc),
                                    )
                                )

                                continue

                            # ----------------------------------------------
                            # Compare the complete exponent data.
                            # ----------------------------------------------

                            same_k1 = (
                                pert_result["k1"]
                                ==
                                base_result["k1"]
                            )

                            same_digits = (
                                pert_result["digits"]
                                ==
                                base_result["digits"]
                            )

                            same_slopes = (
                                pert_result["slopes"]
                                ==
                                base_result["slopes"]
                            )

                            same_exponents = (
                                pert_result["exponents"]
                                ==
                                base_result["exponents"]
                            )

                            same_principal = (
                                pert_result["principal"]
                                ==
                                base_result["principal"]
                            )

                            same_log = (
                                pert_result["logarithmic"]
                                ==
                                base_result["logarithmic"]
                            )

                            if not all(
                                (
                                    same_k1,
                                    same_digits,
                                    same_slopes,
                                    same_exponents,
                                    same_principal,
                                    same_log,
                                )
                            ):

                                failures.append(
                                    (
                                        ell,
                                        k1,
                                        scale,
                                        a,
                                        b,
                                        "projective_data_changed",
                                        {
                                            "k1": same_k1,
                                            "digits": same_digits,
                                            "slopes": same_slopes,
                                            "exponents": same_exponents,
                                            "principal": same_principal,
                                            "logarithmic": same_log,
                                        },
                                    )
                                )

                                continue

                            exact += 1

    # ------------------------------------------------------------------
    # SUMMARY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. PROJECTIVE-REPRESENTATIVE SUMMARY")
    print("=" * 78)

    print(
        f"  cases_tested={cases}"
    )

    print(
        f"  fully_exact_cases={exact}"
    )

    print(
        f"  failures={len(failures)}"
    )

    # ------------------------------------------------------------------
    # FAILURE DETAILS
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
Experiment 252 showed invariance under unit scaling:

    (q1,q3) -> (u*q1,u*q3).

Experiment 253 tests a stronger statement.

Two source pairs may be different integers but represent the same
projective ell-adic class through the tested precision:

    q1/q3 = rho (mod ell^E).

The question is whether the complete exponent coordinate depends only
on that projective class.

If the experiment passes, then the following quantities are projective
invariants of the local chart:

    k1,
    Hensel digits,
    exponent sequence,
    principal coordinate,
    logarithmic coordinate.

In particular, higher-order changes in the integer representatives
cannot alter the recovered exponent as long as the projective ratio is
unchanged to the tested precision.

That is the precise mathematical content behind the phrase
"projective 7-adic coordinate."
"""
    )

    # ------------------------------------------------------------------
    # FINAL
    # ------------------------------------------------------------------

    final_ok = (
        len(failures) == 0
        and cases > 0
        and exact == cases
    )

    print()
    print("=" * 78)
    print("5. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  cases_tested={cases}"
    )

    print(
        f"  exact_projective_representative_invariance="
        f"{final_ok}"
    )

    print(
        f"  exponent_data_invariant="
        f"{final_ok}"
    )

    print(
        f"  logarithmic_coordinate_invariant="
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
    print("EXPERIMENT 253 COMPLETE")


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

