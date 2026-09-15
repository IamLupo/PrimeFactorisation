#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 235 — EXACT PRIME-LOCAL RESIDUAL SLOPE / DERIVATIVE FORMULA AUDIT
==============================================================================

Experiment 234 found exact affine residual laws, but its predicted slope
formula was incomplete.

For

    F(k) = q1 - 2*3^k*q3,

fix a prime ell and let

    L = ord_ell(3),

    k(t) = k0 + t*L.

Because

    3^L = 1 + ell*D_ell (mod ell^2),

we have

    2*3^(k(t))
      = rho * (1 + ell*D_ell*t)
      (mod ell^2),

where

    rho = q1/q3 (mod ell).

Therefore

    F(k(t))/ell
      =
    F(k0)/ell
      -
    q3*rho*D_ell*t
      (mod ell).

Hence the correct affine slope is

    B = -q3*rho*D_ell (mod ell).

Experiment 235 tests this formula for:

    ell = 5,
    ell = 7,
    ell = 73.

It also verifies:

    * complete affine residual law;
    * unique lift digit;
    * reconstruction of the successful exponent;
    * agreement with direct modular validation.

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

SMALL_PRIMES = [
    5,
    7,
    73,
]


# ============================================================================
# HELPERS
# ============================================================================

def valuation(
    x: int,
    p: int,
) -> int | None:

    x = int(x)
    p = int(p)

    if x == 0:
        return None

    x = abs(x)
    v = 0

    while x % p == 0:
        x //= p
        v += 1

    return v


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
                int(modulus),
            )
        )
        % modulus
    )


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
            current * base
        ) % modulus

        if current == 1:
            return e

    raise ArithmeticError(
        f"Order not found for {base} modulo {modulus}."
    )


def residual(
    k: int,
) -> int:

    return int(
        Q1
        -
        (
            2
            * pow(
                3,
                int(k),
            )
            * Q3
        )
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 235 — EXACT PRIME-LOCAL RESIDUAL SLOPE / "
        "DERIVATIVE FORMULA AUDIT"
    )
    print("=" * 78)

    failures = 0

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

    for ell in SMALL_PRIMES:

        print()
        print("=" * 78)
        print(
            f"2. LOCAL ANALYSIS ell={ell}"
        )
        print("=" * 78)

        modulus2 = ell * ell

        rho = ratio_mod(
            ell
        )

        rho2 = ratio_mod(
            modulus2
        )

        q3_mod_ell = Q3 % ell

        order = multiplicative_order(
            3,
            ell,
        )

        pow3L_mod_ell2 = pow(
            3,
            order,
            modulus2,
        )

        if (
            (pow3L_mod_ell2 - 1)
            % ell
            != 0
        ):

            print(
                "  ERROR: local increment is not divisible by ell."
            )

            failures += 1
            continue

        D = int(
            (
                pow3L_mod_ell2 - 1
            )
            // ell
            % ell
        )

        # Correct derivative.
        predicted_B = int(
            (
                -q3_mod_ell
                * rho
                * D
            )
            % ell
        )

        print(
            f"  rho_mod_ell={rho}"
        )

        print(
            f"  rho_mod_ell2={rho2}"
        )

        print(
            f"  q3_mod_ell={q3_mod_ell}"
        )

        print(
            f"  ord_ell(3)={order}"
        )

        print(
            f"  3^L_mod_ell2={pow3L_mod_ell2}"
        )

        print(
            f"  D=((3^L-1)/ell)_mod_ell={D}"
        )

        print(
            f"  predicted_B=-q3*rho*D mod ell="
            f"{predicted_B}"
        )

        # --------------------------------------------------------------
        # Base exponent.
        # --------------------------------------------------------------

        base_matches = [
            k
            for k in range(order)
            if orbit_mod(
                k,
                ell,
            ) == rho
        ]

        print(
            f"  base_matches={base_matches}"
        )

        if len(base_matches) != 1:

            failures += 1
            continue

        k0 = base_matches[0]

        print(
            f"  k0={k0}"
        )

        # --------------------------------------------------------------
        # Normalized residual table.
        # --------------------------------------------------------------

        normalized = []

        print(
            "  normalized_residuals:"
        )

        local_division_ok = True

        for t in range(ell):

            k = int(
                k0 + t * order
            )

            F = residual(
                k
            )

            if F % ell != 0:

                local_division_ok = False

                print(
                    f"    t={t}: "
                    f"k={k} "
                    f"ERROR:not_divisible"
                )

                continue

            R = int(
                (F // ell) % ell
            )

            normalized.append(
                (
                    t,
                    R,
                )
            )

            print(
                f"    t={t}: "
                f"k={k} "
                f"R(t)={R}"
            )

        if not local_division_ok:
            failures += 1
            continue

        # --------------------------------------------------------------
        # Measured slope.
        # --------------------------------------------------------------

        A = normalized[0][1]

        measured_slopes = []

        for i in range(
            1,
            len(normalized),
        ):

            previous = normalized[i - 1][1]
            current = normalized[i][1]

            measured_slopes.append(
                (
                    current - previous
                ) % ell
            )

        measured_slope_set = sorted(
            set(
                measured_slopes
            )
        )

        affine_exact = (
            len(measured_slope_set) == 1
        )

        measured_B = (
            measured_slope_set[0]
            if affine_exact
            else None
        )

        print(
            f"  intercept_A={A}"
        )

        print(
            f"  measured_slope_set="
            f"{measured_slope_set}"
        )

        print(
            f"  affine_exact="
            f"{affine_exact}"
        )

        if not affine_exact:

            failures += 1
            continue

        slope_exact = (
            measured_B
            == predicted_B
        )

        print(
            f"  measured_B={measured_B}"
        )

        print(
            f"  slope_formula_exact="
            f"{slope_exact}"
        )

        if not slope_exact:
            failures += 1

        # --------------------------------------------------------------
        # Full affine law.
        # --------------------------------------------------------------

        affine_failures = []

        for t, actual in normalized:

            expected = int(
                (
                    A
                    + t * measured_B
                )
                % ell
            )

            if actual != expected:

                affine_failures.append(
                    (
                        t,
                        actual,
                        expected,
                    )
                )

        full_affine_exact = (
            len(affine_failures) == 0
        )

        print(
            f"  full_affine_exact="
            f"{full_affine_exact}"
        )

        if not full_affine_exact:

            print(
                f"  affine_failures="
                f"{affine_failures}"
            )

            failures += 1

        # --------------------------------------------------------------
        # Recover lift digit.
        # --------------------------------------------------------------

        if gcd(
            measured_B,
            ell,
        ) != 1:

            print(
                "  ERROR: measured slope is not invertible."
            )

            failures += 1
            continue

        t_recovered = int(
            (
                (-A)
                * pow(
                    measured_B,
                    -1,
                    ell,
                )
            )
            % ell
        )

        recovered_k = int(
            k0
            + t_recovered * order
        )

        recovered_orbit = orbit_mod(
            recovered_k,
            modulus2
        )

        lift_exact = (
            recovered_orbit
            == rho2
        )

        print(
            f"  recovered_t={t_recovered}"
        )

        print(
            f"  recovered_k={recovered_k}"
        )

        print(
            f"  recovered_orbit_mod_ell2="
            f"{recovered_orbit}"
        )

        print(
            f"  recovered_lift_exact="
            f"{lift_exact}"
        )

        if not lift_exact:
            failures += 1

        # --------------------------------------------------------------
        # Direct validation.
        # --------------------------------------------------------------

        direct = []

        for t in range(ell):

            candidate = int(
                k0
                + t * order
            )

            if (
                orbit_mod(
                    candidate,
                    modulus2,
                )
                == rho2
            ):

                direct.append(
                    (
                        t,
                        candidate,
                    )
                )

        direct_unique = (
            len(direct) == 1
        )

        affine_matches_direct = (
            direct_unique
            and
            direct[0][0]
            == t_recovered
        )

        print(
            f"  direct_lifts={direct}"
        )

        print(
            f"  direct_unique="
            f"{direct_unique}"
        )

        print(
            f"  affine_matches_direct="
            f"{affine_matches_direct}"
        )

        if not direct_unique:
            failures += 1

        if not affine_matches_direct:
            failures += 1

        # --------------------------------------------------------------
        # Compact local theorem line.
        # --------------------------------------------------------------

        print(
            "  LOCAL AFFINE LAW:"
        )

        print(
            f"    R(t)="
            f"{A}+{measured_B}*t mod {ell}"
        )

        print(
            f"    predicted slope="
            f"-({q3_mod_ell})*({rho})*({D})"
            f" mod {ell}"
        )

        print(
            f"    recovered t="
            f"{t_recovered}"
        )

    # ------------------------------------------------------------------
    # 3. UNIVERSAL LOCAL SLOPE STATEMENT
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. CROSS-PRIME DERIVATIVE FORMULA")
    print("=" * 78)

    print(
        """
For each tested prime ell:

    B_ell
      = -q3*rho*D_ell (mod ell),

where

    rho = q1/q3 (mod ell),

and

    D_ell
      = (3^ord_ell(3)-1)/ell (mod ell).

The factor q3 is essential because the residual is

    q1 - c_k*q3,

not merely

    rho - c_k.
"""
    )

    for ell in SMALL_PRIMES:

        rho = ratio_mod(
            ell
        )

        order = multiplicative_order(
            3,
            ell,
        )

        D = int(
            (
                (
                    pow(
                        3,
                        order,
                        ell * ell,
                    )
                    - 1
                )
                // ell
            )
            % ell
        )

        q3_mod = Q3 % ell

        B = int(
            (
                -q3_mod
                * rho
                * D
            )
            % ell
        )

        print(
            f"  ell={ell}: "
            f"q3={q3_mod} "
            f"rho={rho} "
            f"D={D} "
            f"B={B}"
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
Experiment 234 already showed that the local normalized residual is
affine in the Hensel digit.

The corrected derivative calculation explains the slope.

With

    k(t)=k0+tL,

    L=ord_ell(3),

and

    3^L = 1 + ell*D (mod ell^2),

we obtain

    2*3^(k(t))
      = rho + ell*rho*D*t (mod ell^2).

Multiplying by q3 and subtracting from q1 gives

    F(k(t))/ell
      =
    A - q3*rho*D*t (mod ell).

Hence

    B = -q3*rho*D mod ell.

The successful digit is therefore determined by

    t = -A/B mod ell.

This is the exact finite-prime analogue of the slope-one 7-adic
linearization. The slope need not equal 1; the general local derivative
is the unit B above.
"""
    )

    # ------------------------------------------------------------------
    # 5. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        failures == 0
    )

    print()
    print("=" * 78)
    print("5. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  primes_tested="
        f"{SMALL_PRIMES}"
    )

    print(
        "  corrected_derivative_formula="
        "B=-q3*rho*D"
    )

    print(
        f"  all_local_affine_laws_exact="
        f"{final_ok}"
    )

    print(
        f"  all_slope_formulas_exact="
        f"{final_ok}"
    )

    print(
        f"  all_affine_digit_recoveries_exact="
        f"{final_ok}"
    )

    print(
        f"  failures={failures}"
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
    print("EXPERIMENT 235 COMPLETE")


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

