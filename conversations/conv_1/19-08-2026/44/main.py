#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 211 — EXACT 7-ADIC PRINCIPAL-COORDINATE LINEARIZATION AUDIT
==============================================================================

Goal
----

The experiments have established the persistent normalized residual law

    R_e(t) = A_e + t  (mod 7).

Experiment 211 derives the coefficient "1" directly from the principal
coordinate

    k = 1 + 6m.

Then

    2*3^k
      = 6 * (3^6)^m
      = 6 * 729^m.

Since

    729 - 1 = 728 = 7*104,

we have

    v7(729-1) = 1,

and modulo 7,

    (729-1)/7 = 104 = 6 (mod 7).

Because

    rho = q1/q3 = 6 (mod 7),

the normalized one-digit response is predicted to have slope

    rho * ((729-1)/7)
      = 6*6
      = 36
      = 1 (mod 7).

This experiment verifies that identity directly at every tested level.

The key point is that the slope-one phenomenon is no longer merely
empirical: it is checked against the exact multiplicative linearization
of 729^m.

No giant powers of 3 are constructed.
No exponent sequence is hard-coded.
"""

from __future__ import annotations

import sys


# ============================================================================
# DATA
# ============================================================================

Q1 = 29144191
Q3 = 24794967

P = 7

# Principal-coordinate base:
#
#     k = 1 + 6m
#     2*3^k = 6*729^m
#
PRINCIPAL_BASE = 3 ** 6       # 729
PRINCIPAL_SCALE = 6           # 2*3


# ============================================================================
# BASIC HELPERS
# ============================================================================

def valuation_7(x: int):
    if x == 0:
        return None

    x = abs(x)
    e = 0

    while x % P == 0:
        x //= P
        e += 1

    return e


def inverse_mod(a: int, m: int) -> int:

    if m <= 1:
        return 0

    a %= m

    if a == 0:
        raise ArithmeticError(
            f"Cannot invert 0 modulo {m}."
        )

    old_r, r = a, m
    old_s, s = 1, 0

    while r:

        q = old_r // r

        old_r, r = (
            r,
            old_r - q * r,
        )

        old_s, s = (
            s,
            old_s - q * s,
        )

    if old_r != 1:
        raise ArithmeticError(
            f"{a} is not invertible modulo {m}."
        )

    return old_s % m


def ratio_mod(modulus: int) -> int:

    return (
        Q1
        * inverse_mod(
            Q3,
            modulus,
        )
    ) % modulus


def principal_orbit_mod(
    m: int,
    modulus: int,
) -> int:

    """
    2*3^(1+6m) mod modulus
      = 6*729^m mod modulus.
    """

    return (
        PRINCIPAL_SCALE
        * pow(
            PRINCIPAL_BASE,
            m,
            modulus,
        )
    ) % modulus


def order_3_mod_7e(e: int) -> int:

    return (
        6
        * (P ** (e - 1))
    )


def order_729_mod_7e(e: int) -> int:

    """
    Since 729 = 3^6,

        ord_{7^e}(729) = 7^(e-1).
    """

    return (
        P ** (e - 1)
    )


def defect_digit_m(
    e: int,
    m: int,
) -> int:

    """
    Principal-coordinate defect:

        delta_e =
          (6*729^m - rho) / 7^e mod 7.

    Only modulus 7^(e+1) is needed.
    """

    modulus_e = P ** e
    modulus_next = P ** (
        e + 1
    )

    rho_next = ratio_mod(
        modulus_next
    )

    orbit_next = principal_orbit_mod(
        m,
        modulus_next
    )

    diff = (
        orbit_next
        - rho_next
    ) % modulus_next

    if diff % modulus_e != 0:
        raise ArithmeticError(
            f"m={m} is not a solution modulo 7^{e}."
        )

    return (
        diff // modulus_e
    ) % P


def candidate_residual_digit(
    e: int,
    m: int,
    t: int,
) -> int:

    """
    Evaluate the normalized residual digit for

        m' = m + t*7^(e-1).
    """

    m_candidate = (
        m
        + t * order_729_mod_7e(e)
    )

    return defect_digit_m(
        e,
        m_candidate,
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 211 — EXACT 7-ADIC "
        "PRINCIPAL-COORDINATE LINEARIZATION AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. SOURCE / PRINCIPAL BASE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. SOURCE / PRINCIPAL-COORDINATE DATA")
    print("=" * 78)

    q1_v7 = valuation_7(Q1)
    q3_v7 = valuation_7(Q3)

    rho7 = ratio_mod(7)

    principal_base_minus_one = (
        PRINCIPAL_BASE - 1
    )

    base_v7 = valuation_7(
        principal_base_minus_one
    )

    principal_unit = (
        principal_base_minus_one // P
    )

    print(
        f"  q1={Q1}"
    )

    print(
        f"  q3={Q3}"
    )

    print(
        f"  v7(q1)={q1_v7}"
    )

    print(
        f"  v7(q3)={q3_v7}"
    )

    print(
        f"  rho_mod7={rho7}"
    )

    print(
        f"  principal_base={PRINCIPAL_BASE}"
    )

    print(
        f"  principal_base_minus_1="
        f"{principal_base_minus_one}"
    )

    print(
        f"  v7(729-1)={base_v7}"
    )

    print(
        f"  (729-1)/7={principal_unit}"
    )

    print(
        f"  ((729-1)/7) mod7="
        f"{principal_unit % P}"
    )

    # ------------------------------------------------------------------
    # 2. PREDICTED LINEARIZATION SLOPE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. PREDICTED LINEARIZATION SLOPE")
    print("=" * 78)

    rho_unit = rho7

    normalized_729_increment = (
        principal_unit
        % P
    )

    predicted_slope = (
        rho_unit
        * normalized_729_increment
    ) % P

    print(
        f"  rho_mod7={rho_unit}"
    )

    print(
        f"  normalized_729_increment_mod7="
        f"{normalized_729_increment}"
    )

    print(
        f"  predicted_slope="
        f"{predicted_slope}"
    )

    print(
        f"  predicted_slope_is_1="
        f"{predicted_slope == 1}"
    )

    # ------------------------------------------------------------------
    # 3. BASE PRINCIPAL EXPONENT
    # ------------------------------------------------------------------

    current_m = 0

    print()
    print("=" * 78)
    print("3. INITIAL PRINCIPAL COORDINATE")
    print("=" * 78)

    print(
        "  k=1+6m"
    )

    print(
        f"  m={current_m}"
    )

    print(
        f"  orbit_mod7="
        f"{principal_orbit_mod(current_m, 7)}"
    )

    print(
        f"  rho_mod7={rho7}"
    )

    # ------------------------------------------------------------------
    # 4. LINEARIZED SEVEN-CANDIDATE AUDIT
    # ------------------------------------------------------------------

    transition_rows = []

    print()
    print("=" * 78)
    print("4. EXACT SEVEN-CANDIDATE LINEARIZATION")
    print("=" * 78)

    for e in range(
        1,
        12,
    ):

        current_m = int(
            current_m
        )

        base_digit = defect_digit_m(
            e,
            current_m,
        )

        candidate_rows = []

        for t in range(7):

            digit = candidate_residual_digit(
                e,
                current_m,
                t,
            )

            candidate_rows.append(
                {
                    "t": t,
                    "digit": digit,
                }
            )

        A = candidate_rows[0]["digit"]

        B = (
            candidate_rows[1]["digit"]
            - candidate_rows[0]["digit"]
        ) % P

        affine_exact = all(
            row["digit"]
            ==
            (
                A
                + B * row["t"]
            ) % P
            for row in candidate_rows
        )

        predicted_by_unit = (
            A
            + predicted_slope * 0
        ) % P

        chosen_t = (
            (-A)
            % P
        )

        chosen_m = (
            current_m
            + chosen_t
            * order_729_mod_7e(e)
        )

        alignment = (
            principal_orbit_mod(
                chosen_m,
                P ** (e + 1),
            )
            ==
            ratio_mod(
                P ** (e + 1)
            )
        )

        transition_rows.append(
            {
                "e": e,
                "m": current_m,
                "A": A,
                "B": B,
                "affine": affine_exact,
                "chosen_t": chosen_t,
                "chosen_m": chosen_m,
                "alignment": alignment,
                "candidate_rows": candidate_rows,
            }
        )

        print(
            f"  e={e}: "
            f"m={current_m} "
            f"A={A} "
            f"B={B} "
            f"predicted_B={predicted_slope} "
            f"affine={affine_exact} "
            f"chosen_t={chosen_t} "
            f"next_m={chosen_m} "
            f"alignment={alignment}"
        )

        print(
            "    candidates="
            + str(
                [
                    row["digit"]
                    for row in candidate_rows
                ]
            )
        )

        current_m = (
            chosen_m
        )

    # ------------------------------------------------------------------
    # 5. SLOPE COMPARISON
    # ------------------------------------------------------------------

    actual_slopes = [
        row["B"]
        for row in transition_rows
    ]

    affine_flags = [
        row["affine"]
        for row in transition_rows
    ]

    slope_matches = [
        B == predicted_slope
        for B in actual_slopes
    ]

    print()
    print("=" * 78)
    print("5. ACTUAL VS PREDICTED SLOPE")
    print("=" * 78)

    print(
        f"  predicted_slope={predicted_slope}"
    )

    print(
        f"  actual_slopes={actual_slopes}"
    )

    print(
        f"  slope_matches={slope_matches}"
    )

    # ------------------------------------------------------------------
    # 6. LOCAL MULTIPLICATIVE LINEARIZATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. LOCAL MULTIPLICATIVE LINEARIZATION")
    print("=" * 78)

    exact_local_flags = []

    for e in range(
        1,
        8,
    ):

        modulus = P ** (
            e + 1
        )

        h = P ** (
            e - 1
        )

        # Direct local unit:
        #
        # 729^(t*7^(e-1))
        #
        # compared with
        #
        # 1 + t*(729-1)*7^(e-1).
        #
        local_base = pow(
            PRINCIPAL_BASE,
            h,
            modulus,
        )

        linear_unit = (
            1
            + (
                PRINCIPAL_BASE - 1
            )
            * h
        ) % modulus

        exact = (
            local_base
            == linear_unit
        )

        exact_local_flags.append(
            exact
        )

        print(
            f"  e={e}: "
            f"7^(e-1)={h} "
            f"729^(7^(e-1)) mod7^(e+1)="
            f"{local_base} "
            f"linearized={linear_unit} "
            f"exact={exact}"
        )

    # ------------------------------------------------------------------
    # 7. THEORETICAL SLOPE IDENTITY
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. SLOPE IDENTITY")
    print("=" * 78)

    print(
        """
For

    k = 1 + 6m,

the orbit is

    c(m) = 6*729^m.

A one-digit change

    m -> m + t*7^(e-1)

gives, modulo 7^(e+1),

    729^(t*7^(e-1))
      = 1
        + t*(729-1)*7^(e-1).

Since

    729-1 = 7*104,

the normalized first-order change is controlled by

    104 mod 7 = 6.

The orbit residue itself is

    rho = 6 mod 7.

Therefore the normalized residual slope is

    6*6 = 36 = 1 mod 7.

Experiment 211 checks every component of this derivation exactly.
"""
    )

    # ------------------------------------------------------------------
    # 8. FINAL
    # ------------------------------------------------------------------

    all_slope_match = all(
        slope_matches
    )

    all_affine = all(
        affine_flags
    )

    all_alignment = all(
        row["alignment"]
        for row in transition_rows
    )

    all_local_linearization = all(
        exact_local_flags
    )

    final_ok = (
        q1_v7 == 0
        and q3_v7 == 0
        and predicted_slope == 1
        and all_slope_match
        and all_affine
        and all_alignment
        and all_local_linearization
    )

    print()
    print("=" * 78)
    print("8. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  v7(q1)=0="
        f"{q1_v7 == 0}"
    )

    print(
        f"  v7(q3)=0="
        f"{q3_v7 == 0}"
    )

    print(
        f"  v7(729-1)=1="
        f"{base_v7 == 1}"
    )

    print(
        f"  predicted_slope={predicted_slope}"
    )

    print(
        f"  predicted_slope_is_1="
        f"{predicted_slope == 1}"
    )

    print(
        f"  actual_slopes={actual_slopes}"
    )

    print(
        f"  all_actual_slopes_match_prediction="
        f"{all_slope_match}"
    )

    print(
        f"  all_affine_exact="
        f"{all_affine}"
    )

    print(
        f"  all_lift_alignments_exact="
        f"{all_alignment}"
    )

    print(
        f"  local_multiplicative_linearization_exact="
        f"{all_local_linearization}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 211 COMPLETE")


if __name__ == "__main__":

    try:
        main()

    except KeyboardInterrupt:
        print("\nInterrupted.")
        sys.exit(130)

    except Exception as exc:
        print(
            "\nFATAL ERROR: "
            + type(exc).__name__
            + ": "
            + str(exc)
        )
        raise

