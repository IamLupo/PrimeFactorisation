#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 236 — FAST SMALL-PRIME UNIVERSAL LOCAL DERIVATIVE AUDIT
==============================================================================

Experiment 235 established the local derivative formula for the three
repeated primes 5, 7, and 73:

    B_ell = -q3 * rho * D_ell (mod ell),

where

    rho = q1/q3 (mod ell),

and

    D_ell =
        (3^ord_ell(3)-1)/ell (mod ell).

Experiment 236 tests the same formula across a broader finite set of
small primes, without factoring large orbit residuals.

For every prime ell <= PRIME_LIMIT with ell not dividing q3:

    1. compute ord_ell(3);
    2. find the unique base exponent k0 satisfying
           2*3^k0 = rho (mod ell);
    3. form
           k(t) = k0 + t*ord_ell(3),
           t=0,...,ell-1;
    4. compute the normalized residual
           R(t) =
             (q1 - 2*3^k(t)*q3)/ell (mod ell);
    5. verify that R(t) is affine;
    6. verify the predicted slope
           B = -q3*rho*D (mod ell);
    7. recover the unique lift digit
           t = -A/B (mod ell);
    8. independently verify the lifted exponent modulo ell^2.

No large integer factorization is performed.

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

PRIME_LIMIT = 200


# ============================================================================
# HELPERS
# ============================================================================

def is_prime(
    n: int,
) -> bool:

    n = int(n)

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

    value = 1

    for e in range(
        1,
        modulus + 1,
    ):

        value = (
            value
            * base
        ) % modulus

        if value == 1:
            return e

    raise ArithmeticError(
        f"Could not determine order of {base} modulo {modulus}."
    )


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
                    int(modulus),
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
        "EXPERIMENT 236 — FAST SMALL-PRIME UNIVERSAL "
        "LOCAL DERIVATIVE AUDIT"
    )
    print("=" * 78)

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

    primes = [
        p
        for p in range(
            2,
            PRIME_LIMIT + 1,
        )
        if is_prime(p)
    ]

    print(
        f"  prime_limit={PRIME_LIMIT}"
    )

    print(
        f"  primes_tested={len(primes)}"
    )

    # ------------------------------------------------------------------
    # Per-prime results
    # ------------------------------------------------------------------

    tested = []
    skipped = []

    total_failures = 0

    for ell in primes:

        # Skip primes dividing q3.
        if Q3 % ell == 0:

            skipped.append(
                (
                    ell,
                    "ell_divides_q3",
                )
            )

            continue

        # We need 3 to be a unit.
        if ell == 3:

            skipped.append(
                (
                    ell,
                    "ell_equals_3",
                )
            )

            continue

        rho = ratio_mod(
            ell
        )

        order = multiplicative_order(
            3,
            ell,
        )

        # --------------------------------------------------------------
        # Base class
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
                f"ell={ell}: "
                f"BASE_CLASS_FAILURE "
                f"matches={base_matches}"
            )

            total_failures += 1

            continue

        k0 = base_matches[0]

        # --------------------------------------------------------------
        # Local derivative increment
        # --------------------------------------------------------------

        modulus2 = ell * ell

        pow3L = pow(
            3,
            order,
            modulus2,
        )

        if (
            pow3L - 1
        ) % ell != 0:

            print(
                f"ell={ell}: "
                f"DERIVATIVE_DIVISIBILITY_FAILURE"
            )

            total_failures += 1

            continue

        D = int(
            (
                (pow3L - 1)
                // ell
            )
            % ell
        )

        q3_mod = Q3 % ell

        predicted_B = int(
            (
                -q3_mod
                * rho
                * D
            )
            % ell
        )

        # --------------------------------------------------------------
        # Affine residual test
        # --------------------------------------------------------------

        normalized = []

        local_division_ok = True

        for t in range(
            ell
        ):

            k = int(
                k0
                + t * order
            )

            # F(k) is known to be divisible by ell because k is in
            # the base projective class.
            F_mod_ell2 = residual_mod(
                k,
                modulus2,
            )

            if (
                F_mod_ell2 % ell
                != 0
            ):

                local_division_ok = False
                break

            R = int(
                (
                    F_mod_ell2
                    // ell
                )
                % ell
            )

            normalized.append(
                R
            )

        if not local_division_ok:

            print(
                f"ell={ell}: "
                f"NORMALIZATION_FAILURE"
            )

            total_failures += 1

            continue

        A = normalized[0]

        slopes = []

        for i in range(
            1,
            ell
        ):

            slopes.append(
                int(
                    (
                        normalized[i]
                        - normalized[i - 1]
                    )
                    % ell
                )
            )

        slope_set = sorted(
            set(
                slopes
            )
        )

        affine_exact = (
            len(slope_set) == 1
        )

        measured_B = (
            slope_set[0]
            if affine_exact
            else None
        )

        slope_formula_exact = (
            affine_exact
            and
            measured_B == predicted_B
        )

        # --------------------------------------------------------------
        # Full affine identity
        # --------------------------------------------------------------

        affine_identity_exact = True

        if measured_B is None:

            affine_identity_exact = False

        else:

            for t in range(
                ell
            ):

                expected = int(
                    (
                        A
                        + t * measured_B
                    )
                    % ell
                )

                if (
                    expected
                    != normalized[t]
                ):

                    affine_identity_exact = False
                    break

        # --------------------------------------------------------------
        # Recover lift digit
        # --------------------------------------------------------------

        invertible_slope = (
            measured_B is not None
            and
            gcd(
                measured_B,
                ell,
            ) == 1
        )

        if invertible_slope:

            recovered_t = int(
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

            recovered_t = None

        # --------------------------------------------------------------
        # Independent lift validation
        # --------------------------------------------------------------

        if recovered_t is not None:

            lifted_k = int(
                k0
                + recovered_t * order
            )

            rho2 = ratio_mod(
                modulus2
            )

            lifted_orbit = orbit_mod(
                lifted_k,
                modulus2,
            )

            lift_exact = (
                lifted_orbit
                == rho2
            )

            direct_matches = [
                t
                for t in range(
                    ell
                )
                if orbit_mod(
                    k0 + t * order,
                    modulus2,
                )
                == rho2
            ]

            unique_lift = (
                len(direct_matches) == 1
            )

            reconstruction_exact = (
                unique_lift
                and
                direct_matches[0]
                == recovered_t
            )

        else:

            lifted_k = None
            lift_exact = False
            direct_matches = []
            unique_lift = False
            reconstruction_exact = False

        prime_ok = (
            base_unique
            and
            local_division_ok
            and
            affine_exact
            and
            slope_formula_exact
            and
            affine_identity_exact
            and
            invertible_slope
            and
            lift_exact
            and
            unique_lift
            and
            reconstruction_exact
        )

        if not prime_ok:
            total_failures += 1

        tested.append(
            {
                "ell": ell,
                "rho": rho,
                "order": order,
                "k0": k0,
                "D": D,
                "predicted_B": predicted_B,
                "measured_B": measured_B,
                "A": A,
                "t": recovered_t,
                "lifted_k": lifted_k,
                "prime_ok": prime_ok,
            }
        )

        print(
            f"  ell={ell}: "
            f"rho={rho} "
            f"ord={order} "
            f"k0={k0} "
            f"D={D} "
            f"B_pred={predicted_B} "
            f"B_meas={measured_B} "
            f"A={A} "
            f"t={recovered_t} "
            f"k_lift={lifted_k} "
            f"OK={prime_ok}"
        )

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. SUMMARY")
    print("=" * 78)

    successful = [
        row
        for row in tested
        if row["prime_ok"]
    ]

    failed = [
        row
        for row in tested
        if not row["prime_ok"]
    ]

    print(
        f"  usable_primes={len(tested)}"
    )

    print(
        f"  successful_primes={len(successful)}"
    )

    print(
        f"  failed_primes={len(failed)}"
    )

    print(
        f"  skipped={skipped}"
    )

    print(
        f"  failures="
        f"{failed}"
    )

    # ------------------------------------------------------------------
    # Formula comparison examples
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. DERIVATIVE FORMULA CROSS-CHECK")
    print("=" * 78)

    for row in tested:

        print(
            f"  ell={row['ell']}: "
            f"B_pred={row['predicted_B']} "
            f"B_meas={row['measured_B']} "
            f"exact="
            f"{row['predicted_B'] == row['measured_B']}"
        )

    # ------------------------------------------------------------------
    # Structural interpretation
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
Experiment 235 established the corrected local derivative formula

    B_ell
      = -q3 * rho * D_ell mod ell,

with

    rho = q1/q3 mod ell,

and

    D_ell =
      (3^ord_ell(3)-1)/ell mod ell.

Experiment 236 tests that formula over a broader finite set of small
primes rather than only the three repeated primes 5, 7, and 73.

For every usable prime, the complete local statement is:

    R(t) = A + B_ell*t mod ell,

and the unique next exponent digit is

    t = -A/B_ell mod ell.

Thus the prime-local Hensel step has three independently testable
components:

    projective base class,
    affine local residual,
    unique lifted exponent.

A positive result across the full finite prime sample would strongly
support that this is the general local mechanism behind the repeated
prime-factor phenomena observed in the orbit residuals.
"""
    )

    # ------------------------------------------------------------------
    # Final exactness
    # ------------------------------------------------------------------

    final_ok = (
        total_failures == 0
        and
        len(tested) > 0
    )

    print()
    print("=" * 78)
    print("5. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  prime_limit={PRIME_LIMIT}"
    )

    print(
        f"  usable_prime_count={len(tested)}"
    )

    print(
        f"  all_base_classes_unique="
        f"{all(
            row['k0'] is not None
            for row in tested
        )}"
    )

    print(
        f"  all_local_affine_laws_exact="
        f"{len(failed) == 0}"
    )

    print(
        f"  universal_local_slope_formula_exact="
        f"{len(failed) == 0}"
    )

    print(
        f"  all_lift_digits_reconstructed="
        f"{len(failed) == 0}"
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
    print("EXPERIMENT 236 COMPLETE")


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

