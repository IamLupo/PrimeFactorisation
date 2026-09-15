#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 250 — EXACT GENERAL PRIME PRINCIPAL-COORDINATE / LOGARITHMIC AUDIT
==============================================================================

Experiment 249R established, for 185 admissible nonexceptional primes <= 2000:

    rho -> k1 -> Hensel digits -> k_r

and

    m_r = t1 + t2*ell + ...

where

    k_r = k1 + d*m_r,

    d = ord_ell(3).

Experiment 250 tests the multiplicative principal-coordinate identity
directly.

For a fixed admissible nonexceptional prime ell, define

    d = ord_ell(3),

    k = k1 + d*m.

Then

    3^k
      = 3^k1 * (3^d)^m.

Therefore the normalized source ratio satisfies

    rho / (2*3^k1)
      = (3^d)^m.

The experiment computes the principal coordinate m from the Hensel
digits and independently checks this multiplicative identity modulo
ell^r.

Because 3^d = 1 mod ell, the finite logarithmic coordinate can also
be computed by truncated formal logarithms on the principal unit:

    log(rho/(2*3^k1))
    ------------------
          log(3^d)

modulo ell^(r-1),

provided the denominator is a principal unit and the standard
nonexceptional condition holds.

The main purpose is not to rely on floating point or an external
p-adic library. Everything is performed with exact modular arithmetic.

The experiment checks:

    1. Hensel principal coordinate;
    2. direct exponential reconstruction;
    3. formal-logarithmic coordinate reconstruction;
    4. equality of the two coordinates;
    5. precision-by-precision agreement.

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


def ratio_mod(modulus: int) -> int | None:

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


def orbit_mod(k: int, modulus: int) -> int:

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


def residual_mod(k: int, modulus: int) -> int:

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


def valuation(
    x: int,
    p: int,
    max_v: int,
) -> int:

    x = int(x)

    if x == 0:
        return max_v

    v = 0

    while (
        v < max_v
        and
        x % p == 0
    ):
        x //= p
        v += 1

    return int(v)


def formal_log(
    x: int,
    p: int,
    precision: int,
) -> int:
    """
    Exact truncated p-adic logarithm for x = 1 + u,
    with u divisible by p.

        log(1+u)
          = u - u^2/2 + u^3/3 - ...

    Terms are evaluated modulo p^precision.

    This routine uses exact rational arithmetic term-by-term:
    u^n / n is divisible enough by p that modular inversion is used
    only after removing the p-adic valuation of n.
    """

    modulus = p ** precision

    u = (x - 1) % modulus

    if u % p != 0:
        raise ValueError(
            "formal_log requires x == 1 mod p."
        )

    result = 0
    power = u

    # Terms with n >= precision+2 have sufficiently large p-adic
    # valuation for the present finite precision.
    for n in range(
        1,
        precision * p + 3,
    ):

        if power % modulus == 0:
            break

        denominator = n

        # Remove powers of p from denominator.
        vden = 0

        while denominator % p == 0:
            denominator //= p
            vden += 1

        # The p-adic valuation of power/n determines whether the term
        # contributes at the requested precision.
        vp_power = valuation(
            power,
            p,
            precision + vden + 2,
        )

        if vp_power - vden >= precision:
            break

        # power = p^vp_power * unit.
        unit_power = power // (p ** vp_power)

        exponent = precision - (
            vp_power - vden
        )

        mod_reduced = p ** exponent

        unit_den = denominator % mod_reduced

        inv_den = pow(
            unit_den,
            -1,
            mod_reduced,
        )

        term = (
            unit_power
            * inv_den
            * (
                p ** (
                    vp_power - vden
                )
            )
        ) % modulus

        if n % 2 == 0:
            result -= term
        else:
            result += term

        result %= modulus

        power *= u
        power %= modulus

    return int(result)


def formal_log_coordinate(
    source_unit: int,
    base_unit: int,
    p: int,
    precision: int,
) -> int:
    """
    Compute

        log(source_unit) / log(base_unit)

    modulo p^(precision-1).

    The denominator has valuation 1 in the nonexceptional branch.
    """

    modulus_log = p ** precision
    modulus_coordinate = p ** (
        precision - 1
    )

    log_source = formal_log(
        source_unit % modulus_log,
        p,
        precision,
    )

    log_base = formal_log(
        base_unit % modulus_log,
        p,
        precision,
    )

    if log_base % p == 0:
        raise ArithmeticError(
            "Unexpected zero principal logarithmic denominator."
        )

    # Divide p-adically by removing one p.
    source_unit_part = log_source // p
    base_unit_part = log_base // p

    numerator = (
        source_unit_part
        % modulus_coordinate
    )

    denominator = (
        base_unit_part
        % modulus_coordinate
    )

    if gcd(
        denominator,
        p,
    ) != 1:
        raise ArithmeticError(
            "Reduced logarithmic denominator is not a unit."
        )

    return int(
        numerator
        * pow(
            denominator,
            -1,
            modulus_coordinate,
        )
        % modulus_coordinate
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 250 — EXACT GENERAL PRIME PRINCIPAL-COORDINATE / "
        "LOGARITHMIC AUDIT"
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

        rho = ratio_mod(ell)

        if rho is None:
            counts["skipped_q3"] += 1
            continue

        # --------------------------------------------------------------
        # Base orbit class
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

        # --------------------------------------------------------------
        # Hensel reconstruction
        # --------------------------------------------------------------

        k = k1
        current_order = d

        exponent_sequence = [k1]
        digit_sequence = []
        principal_sequence = [0]

        local_ok = True

        for level in range(
            1,
            MAX_LEVEL,
        ):

            current_modulus = ell ** level
            next_modulus = ell ** (
                level + 1
            )

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
                        "order_increment",
                    )
                )
                break

            D = (
                order_increment
                // current_modulus
            ) % ell

            if D == 0:
                local_ok = False
                failures.append(
                    (
                        ell,
                        level,
                        "zero_D",
                    )
                )
                break

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
                    break

                normalized.append(
                    (
                        F
                        // current_modulus
                    ) % ell
                )

            if not local_ok:
                failures.append(
                    (
                        ell,
                        level,
                        "normalization",
                    )
                )
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
                        "non_affine",
                    )
                )
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
                local_ok = False
                failures.append(
                    (
                        ell,
                        level,
                        "slope",
                        B,
                        predicted_B,
                    )
                )
                break

            t = (
                (-A)
                * pow(
                    B,
                    -1,
                    ell,
                )
            ) % ell

            k_next = (
                k
                + t * current_order
            )

            target = ratio_mod(
                next_modulus
            )

            if (
                orbit_mod(
                    k_next,
                    next_modulus,
                )
                != target
            ):
                local_ok = False
                failures.append(
                    (
                        ell,
                        level,
                        "alignment",
                    )
                )
                break

            digit_sequence.append(t)
            exponent_sequence.append(
                k_next
            )

            principal_sequence.append(
                (
                    k_next - k1
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

            k = k_next

        if not local_ok:
            continue

        # --------------------------------------------------------------
        # Direct exponential reconstruction
        # --------------------------------------------------------------

        for level, m in enumerate(
            principal_sequence
        ):

            modulus = ell ** level

            if level == 0:
                continue

            source_unit = (
                rho
                * pow(
                    pow(
                        3,
                        k1,
                        modulus * ell,
                    ),
                    -1,
                    modulus * ell,
                )
            ) % (
                modulus * ell
            )

            # Simpler exact identity:
            # rho = 2*3^k1*(3^d)^m
            source_target = ratio_mod(
                ell ** level
            )

            base_factor = (
                2
                * pow(
                    3,
                    k1,
                    ell ** level,
                )
            ) % (
                ell ** level
            )

            base_unit = pow(
                3,
                d,
                ell ** level,
            )

            reconstructed = (
                base_factor
                * pow(
                    base_unit,
                    m,
                    ell ** level,
                )
            ) % (
                ell ** level
            )

            if reconstructed != source_target:

                failures.append(
                    (
                        ell,
                        level,
                        "exponential_reconstruction",
                        reconstructed,
                        source_target,
                    )
                )

                local_ok = False
                break

        if not local_ok:
            continue

        # --------------------------------------------------------------
        # Formal logarithmic comparison
        # --------------------------------------------------------------

        for level in range(
            2,
            MAX_LEVEL + 1,
        ):

            precision = level

            modulus = ell ** precision

            source_ratio = ratio_mod(
                modulus
            )

            base_class_factor = (
                2
                * pow(
                    3,
                    k1,
                    modulus,
                )
            ) % modulus

            inv_base_class = pow(
                base_class_factor,
                -1,
                modulus,
            )

            source_unit = (
                source_ratio
                * inv_base_class
            ) % modulus

            base_unit = pow(
                3,
                d,
                modulus,
            )

            # Both are principal units.
            if source_unit % ell != 1:

                failures.append(
                    (
                        ell,
                        level,
                        "source_not_principal_unit",
                    )
                )

                local_ok = False
                break

            if base_unit % ell != 1:

                failures.append(
                    (
                        ell,
                        level,
                        "base_not_principal_unit",
                    )
                )

                local_ok = False
                break

            try:

                log_m = formal_log_coordinate(
                    source_unit,
                    base_unit,
                    ell,
                    precision,
                )

            except Exception:

                failures.append(
                    (
                        ell,
                        level,
                        "formal_log_failure",
                    )
                )

                local_ok = False
                break

            hensel_m = principal_sequence[
                level - 1
            ]

            expected_modulus = ell ** (
                level - 1
            )

            if (
                log_m
                % expected_modulus
            ) != (
                hensel_m
                % expected_modulus
            ):

                failures.append(
                    (
                        ell,
                        level,
                        "log_coordinate_mismatch",
                        log_m,
                        hensel_m,
                    )
                )

                local_ok = False
                break

        # No per-prime output unless something fails; this keeps the run fast.

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

    # ------------------------------------------------------------------
    # INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
Experiment 249R established the recursive digit coordinate

    k_r = k1 + d*m_r,

    m_r = t1 + t2*ell + t3*ell^2 + ...

for a large admissible nonexceptional sample.

Experiment 250 tests whether this coordinate also has an independent
multiplicative/logarithmic characterization.

After removing the base orbit factor,

    rho / (2*3^k1)
      = (3^d)^m.

Hence the principal coordinate is characterized by

    m = log(rho/(2*3^k1)) / log(3^d),

whenever the ordinary nonexceptional principal-unit chart is valid.

The important point is independence:

    Hensel reconstruction
        versus
    direct exponential reconstruction
        versus
    formal logarithmic reconstruction.

Agreement of all three gives a much stronger local theorem template.
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
    print("EXPERIMENT 250 COMPLETE")


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

