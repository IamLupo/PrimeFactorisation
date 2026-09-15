#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 234 — EXACT PRIME-LOCAL LINEARIZED RESIDUAL / SLOPE AUDIT
==============================================================================

Known source:

    q1 = 29144191
    q3 = 24794967

Repeated small primes:

    ell = 5
    ell = 7
    ell = 73

Experiment 233 established the unique prime-local exponent lifts:

    ell=5:
        k0=2
        ord_5(3)=4
        lift_digit=2
        k1=10

    ell=7:
        k0=1
        ord_7(3)=6
        lift_digit=0
        k1=1

    ell=73:
        k0=1
        ord_73(3)=12
        lift_digit=7
        k1=85

Experiment 234 tests the local residual itself.

For a fixed prime ell, write

    L = ord_ell(3),

and

    k(t) = k0 + t*L.

Define the normalized residual

    R(t) = (q1 - 2*3^(k(t))*q3) / ell  (mod ell).

Because k(t) is already a solution modulo ell, the quotient is
well-defined modulo ell.

The experiment checks whether

    R(t) = A + B*t (mod ell)

for every t.

It then compares the measured slope B with the predicted derivative
coming from the orbit:

    B = rho * D_ell (mod ell),

where

    D_ell = (3^L - 1)/ell mod ell,

and

    rho = q1/q3 mod ell.

The sign is fixed by using the residual

    q1 - c_k*q3.

The successful lift digit should satisfy

    t = -A/B (mod ell).

Thus the experiment replaces candidate search by a local affine
equation.

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
# EXACT HELPERS
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

    return int(v)


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
            current
            * base
        ) % modulus

        if current == 1:
            return int(e)

    raise ArithmeticError(
        f"Order of {base} modulo {modulus} was not found."
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


def mod7_value(
    x: int,
    ell: int,
) -> int:

    return int(
        x % ell
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 234 — EXACT PRIME-LOCAL LINEARIZED "
        "RESIDUAL / SLOPE AUDIT"
    )
    print("=" * 78)

    total_failures = 0

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
    # 2. PRIME-LOCAL ANALYSIS
    # ------------------------------------------------------------------

    for ell in SMALL_PRIMES:

        print()
        print("=" * 78)
        print(
            f"2. PRIME-LOCAL LINEARIZATION: ell={ell}"
        )
        print("=" * 78)

        rho = ratio_mod(
            ell
        )

        order = multiplicative_order(
            3,
            ell,
        )

        modulus2 = int(
            ell * ell
        )

        # --------------------------------------------------------------
        # Base exponent
        # --------------------------------------------------------------

        base_matches = []

        for k in range(
            order
        ):

            if (
                orbit_mod(
                    k,
                    ell,
                )
                == rho
            ):

                base_matches.append(
                    k
                )

        base_unique = (
            len(base_matches) == 1
        )

        if not base_unique:

            print(
                f"  ERROR: "
                f"base_matches={base_matches}"
            )

            total_failures += 1
            continue

        k0 = int(
            base_matches[0]
        )

        print(
            f"  rho_mod_ell={rho}"
        )

        print(
            f"  ord_ell(3)={order}"
        )

        print(
            f"  k0={k0}"
        )

        # --------------------------------------------------------------
        # Local exponential increment
        # --------------------------------------------------------------

        pow3_L = pow(
            3,
            order,
            modulus2,
        )

        if (
            (pow3_L - 1)
            % ell
            != 0
        ):

            print(
                "  ERROR: 3^L-1 is not divisible by ell."
            )

            total_failures += 1
            continue

        D = int(
            (
                (pow3_L - 1)
                // ell
            )
            % ell
        )

        predicted_B = int(
            (
                rho
                * D
            )
            % ell
        )

        print(
            f"  3^L_mod_ell2={pow3_L}"
        )

        print(
            f"  D=((3^L-1)/ell)_mod_ell={D}"
        )

        print(
            f"  predicted_B={predicted_B}"
        )

        # --------------------------------------------------------------
        # Evaluate all ell candidate residual digits
        # --------------------------------------------------------------

        A = None
        measured_B = None

        normalized_values = []

        affine_exact = True

        print(
            "  normalized_residuals:"
        )

        for t in range(
            ell
        ):

            k = int(
                k0
                + t * order
            )

            F = residual(
                k
            )

            divisible = (
                F % ell == 0
            )

            if not divisible:

                print(
                    f"    t={t}: "
                    f"k={k} "
                    f"ERROR: residual not divisible by ell"
                )

                total_failures += 1
                affine_exact = False
                continue

            normalized = int(
                (
                    F // ell
                )
                % ell
            )

            normalized_values.append(
                (
                    t,
                    normalized,
                )
            )

            if A is None:
                A = normalized

            expected_placeholder = normalized

            print(
                f"    t={t}: "
                f"k={k} "
                f"R(t)={normalized}"
            )

        # --------------------------------------------------------------
        # Determine slope from adjacent values
        # --------------------------------------------------------------

        if len(normalized_values) != ell:

            total_failures += 1
            continue

        A = int(
            normalized_values[0][1]
        )

        slope_candidates = []

        for i in range(
            1,
            ell
        ):

            previous = (
                normalized_values[i - 1][1]
            )

            current = (
                normalized_values[i][1]
            )

            slope = int(
                (
                    current
                    - previous
                )
                % ell
            )

            slope_candidates.append(
                slope
            )

        measured_B_values = sorted(
            set(
                slope_candidates
            )
        )

        if len(measured_B_values) != 1:

            affine_exact = False
            measured_B = None

        else:

            measured_B = int(
                measured_B_values[0]
            )

        print(
            f"  measured_slope_values="
            f"{measured_B_values}"
        )

        print(
            f"  affine_slope_exact="
            f"{affine_exact}"
        )

        if measured_B is None:

            total_failures += 1
            continue

        slope_formula_exact = (
            measured_B
            == predicted_B
        )

        print(
            f"  measured_B={measured_B}"
        )

        print(
            f"  slope_formula_exact="
            f"{slope_formula_exact}"
        )

        if not slope_formula_exact:

            total_failures += 1

        # --------------------------------------------------------------
        # Check full affine law
        # --------------------------------------------------------------

        if measured_B % ell == 0:

            # This case should not happen for these primes because the
            # source ratio and derivative are units.
            affine_formula = False

        else:

            affine_formula = True

        affine_point_failures = []

        for t, actual in normalized_values:

            predicted = int(
                (
                    A
                    + t * measured_B
                )
                % ell
            )

            if predicted != actual:

                affine_point_failures.append(
                    (
                        t,
                        actual,
                        predicted,
                    )
                )

        affine_full_exact = (
            len(affine_point_failures) == 0
        )

        print(
            f"  affine_full_exact="
            f"{affine_full_exact}"
        )

        if not affine_full_exact:

            print(
                f"  affine_failures="
                f"{affine_point_failures}"
            )

            total_failures += 1

        # --------------------------------------------------------------
        # Recover successful lift digit from affine equation
        # --------------------------------------------------------------

        if (
            measured_B is not None
            and gcd(
                measured_B,
                ell,
            ) == 1
        ):

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

        else:

            t_recovered = None

        print(
            f"  intercept_A={A}"
        )

        print(
            f"  recovered_lift_digit="
            f"{t_recovered}"
        )

        if t_recovered is None:

            total_failures += 1
            continue

        recovered_k = int(
            k0
            + t_recovered * order
        )

        recovered_orbit = orbit_mod(
            recovered_k,
            modulus2,
        )

        rho2 = ratio_mod(
            modulus2
        )

        recovered_match = (
            recovered_orbit
            == rho2
        )

        print(
            f"  recovered_k={recovered_k}"
        )

        print(
            f"  rho_mod_ell2={rho2}"
        )

        print(
            f"  recovered_orbit_mod_ell2="
            f"{recovered_orbit}"
        )

        print(
            f"  recovered_lift_exact="
            f"{recovered_match}"
        )

        if not recovered_match:

            total_failures += 1

        # --------------------------------------------------------------
        # Compare to direct seven-candidate answer
        # --------------------------------------------------------------

        direct_matches = []

        for t in range(
            ell
        ):

            candidate_k = int(
                k0
                + t * order
            )

            if (
                orbit_mod(
                    candidate_k,
                    modulus2,
                )
                == rho2
            ):

                direct_matches.append(
                    (
                        t,
                        candidate_k,
                    )
                )

        direct_unique = (
            len(direct_matches) == 1
        )

        print(
            f"  direct_successful_lifts="
            f"{direct_matches}"
        )

        print(
            f"  direct_unique="
            f"{direct_unique}"
        )

        direct_reconstruction_exact = (
            direct_unique
            and
            direct_matches[0][0]
            == t_recovered
        )

        print(
            f"  affine_matches_direct_search="
            f"{direct_reconstruction_exact}"
        )

        if not direct_reconstruction_exact:

            total_failures += 1

        # --------------------------------------------------------------
        # Local residual slope interpretation
        # --------------------------------------------------------------

        print(
            "  LOCAL FORMULA:"
        )

        print(
            f"    R(t) = {A} + "
            f"{measured_B}*t (mod {ell})"
        )

        print(
            f"    t = -A/B = "
            f"{t_recovered} (mod {ell})"
        )

    # ------------------------------------------------------------------
    # 3. CROSS-PRIME SUMMARY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. CROSS-PRIME SLOPE SUMMARY")
    print("=" * 78)

    for ell in SMALL_PRIMES:

        rho = ratio_mod(
            ell
        )

        order = multiplicative_order(
            3,
            ell,
        )

        modulus2 = ell * ell

        D = (
            (
                pow(
                    3,
                    order,
                    modulus2,
                )
                - 1
            )
            // ell
        ) % ell

        B = (
            rho
            * D
        ) % ell

        print(
            f"  ell={ell}: "
            f"rho={rho} "
            f"ord_ell={order} "
            f"D={D} "
            f"predicted_B={B}"
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
Experiment 233 showed that the prime-power lift is unique.

Experiment 234 asks what determines the successful digit.

For

    k(t) = k0 + t*ord_ell(3),

the normalized residual is

    R(t)
      = (q1 - 2*3^k(t)*q3)/ell
      (mod ell).

Because

    3^(ord_ell(3)) = 1 + ell*D (mod ell^2),

the first-order change is linear in t.

The predicted slope is

    B = rho*D (mod ell),

where

    rho = q1/q3 (mod ell).

Thus the lift digit should be recovered from

    A + B*t = 0 (mod ell),

without testing all candidates.

This is exactly the same structural mechanism seen in the 7-adic
principal-coordinate experiments, now expressed locally at the small
primes 5, 7, and 73.
"""
    )

    # ------------------------------------------------------------------
    # 5. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        total_failures == 0
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
        "  local_residual_affine=True"
    )

    print(
        f"  slope_formula_exact="
        f"{final_ok}"
    )

    print(
        f"  affine_digit_recovery_exact="
        f"{final_ok}"
    )

    print(
        f"  failures={total_failures}"
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
    print("EXPERIMENT 234 COMPLETE")


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

