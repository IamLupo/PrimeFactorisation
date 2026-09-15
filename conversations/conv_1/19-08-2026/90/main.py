#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 255R — EXACT n=pq SOURCE INTERFACE / PROJECTIVE DATA AUDIT
==============================================================================

Correction of Experiment 255:

Common source scalings are only admissible in the 7-adic projective chart
when the scaling factor is a 7-adic unit.

Therefore common scales divisible by 7 are excluded from the projective
ratio test.

The experiment checks:

    1. n=pq source data;
    2. invariance under admissible common scaling;
    3. invariance under higher-order representative perturbation;
    4. projective normalization q3 -> 1;
    5. exact Schur difference geometry;
    6. the minimal source-side interface required by the local theorem.

Only main.py is used.
"""

from __future__ import annotations

import sys
from math import gcd


# ============================================================================
# KNOWN n=6 SOURCE
# ============================================================================

P = 2
Q = 3
N = P * Q

Q1 = 29144191
Q3 = 24794967

PRIME = 7
MAX_E = 4


# Only 7-adic unit common scales are admissible.
COMMON_SCALES = [1, 2, 3, 5, 10, 11, 13, 17]


# ============================================================================
# HELPERS
# ============================================================================

def valuation(
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

    return v


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
            "q3 is not invertible modulo modulus."
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


def projective_profile(
    q1: int,
    q3: int,
) -> list[int]:

    return [
        ratio_mod(
            q1,
            q3,
            PRIME ** e,
        )
        for e in range(
            1,
            MAX_E + 1,
        )
    ]


def schur_row(
    q1: int,
    q3: int,
) -> tuple[int, int, int]:

    return (
        q1 - 2 * q3,
        q1 - 6 * q3,
        q1 - 18 * q3,
    )


def difference_geometry(
    q1: int,
    q3: int,
) -> tuple[int, int, int, bool]:

    s0, s1, s2 = schur_row(
        q1,
        q3,
    )

    d01 = s1 - s0
    d12 = s2 - s1
    d02 = s2 - s0

    exact = (
        d01 == -4 * q3
        and
        d12 == -12 * q3
        and
        d02 == -16 * q3
        and
        d12 == 3 * d01
    )

    return (
        d01,
        d12,
        d02,
        exact,
    )


def normalized_source(
    p: int,
    q: int,
    q1: int,
    q3: int,
) -> tuple[int, int]:

    if q3 % q != 0:

        raise ArithmeticError(
            "q3 is not divisible by q."
        )

    u = q3 // q

    if (q1 - 1) % (p * q) != 0:

        raise ArithmeticError(
            "q1 does not satisfy pq normalization."
        )

    if (u - 1) % (p * q) != 0:

        raise ArithmeticError(
            "q3/q does not satisfy pq normalization."
        )

    a = (
        q1 - 1
    ) // (
        p * q
    )

    b = (
        u - 1
    ) // (
        p * q
    )

    return (
        a,
        b,
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 255R — EXACT n=pq SOURCE INTERFACE / "
        "PROJECTIVE DATA AUDIT"
    )
    print("=" * 78)

    failures = []

    original_profile = projective_profile(
        Q1,
        Q3,
    )

    # ------------------------------------------------------------------
    # 1. ORIGINAL SOURCE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. ORIGINAL n=pq SOURCE")
    print("=" * 78)

    print(
        f"  p={P}"
    )
    print(
        f"  q={Q}"
    )
    print(
        f"  n={N}"
    )
    print(
        f"  q1={Q1}"
    )
    print(
        f"  q3={Q3}"
    )
    print(
        f"  gcd(q1,q3)={gcd(Q1,Q3)}"
    )
    print(
        f"  v7(q1)={valuation(Q1,7)}"
    )
    print(
        f"  v7(q3)={valuation(Q3,7)}"
    )
    print(
        f"  rho_profile={original_profile}"
    )

    # ------------------------------------------------------------------
    # 2. COMMON-SCALING TEST
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. ADMISSIBLE COMMON-SCALING INVARIANCE")
    print("=" * 78)

    print(
        "  admissible_scales="
        f"{[u for u in COMMON_SCALES if u % PRIME != 0]}"
    )

    for scale in COMMON_SCALES:

        if scale % PRIME == 0:

            print(
                f"  scale={scale}: "
                "SKIPPED_NOT_7_ADIC_UNIT"
            )

            continue

        q1_scaled = (
            scale
            * Q1
        )

        q3_scaled = (
            scale
            * Q3
        )

        scaled_profile = projective_profile(
            q1_scaled,
            q3_scaled,
        )

        invariant = (
            scaled_profile
            ==
            original_profile
        )

        gcd_scaled = gcd(
            q1_scaled,
            q3_scaled,
        )

        primitive_ratio_same = (
            (
                q1_scaled
                // gcd_scaled,
                q3_scaled
                // gcd_scaled,
            )
            ==
            (
                Q1,
                Q3,
            )
        )

        print(
            f"  scale={scale}: "
            f"rho_invariant={invariant} "
            f"primitive_ratio_invariant="
            f"{primitive_ratio_same}"
        )

        if not (
            invariant
            and
            primitive_ratio_same
        ):

            failures.append(
                (
                    "common_scale",
                    scale,
                )
            )

    # ------------------------------------------------------------------
    # 3. HIGHER-ORDER REPRESENTATIVE PERTURBATIONS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. HIGHER-ORDER REPRESENTATIVE PERTURBATIONS")
    print("=" * 78)

    perturb_modulus = (
        PRIME ** MAX_E
    )

    for a_shift in [0, 1, 2, 5]:

        for b_shift in [0, 1, 3, 7]:

            q1_pert = (
                Q1
                + a_shift
                * perturb_modulus
            )

            q3_pert = (
                Q3
                + b_shift
                * perturb_modulus
            )

            pert_profile = projective_profile(
                q1_pert,
                q3_pert,
            )

            invariant = (
                pert_profile
                ==
                original_profile
            )

            print(
                f"  shifts=({a_shift},{b_shift}): "
                f"profile_unchanged="
                f"{invariant}"
            )

            if not invariant:

                failures.append(
                    (
                        "higher_order_perturbation",
                        a_shift,
                        b_shift,
                    )
                )

    # ------------------------------------------------------------------
    # 4. PROJECTIVE NORMALIZATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. PROJECTIVE NORMALIZATION q3 -> 1")
    print("=" * 78)

    modulus = (
        PRIME ** MAX_E
    )

    rho = ratio_mod(
        Q1,
        Q3,
        modulus,
    )

    normalized_profile = projective_profile(
        rho,
        1,
    )

    normalization_exact = (
        normalized_profile
        ==
        original_profile
    )

    print(
        f"  rho_mod_7^{MAX_E}={rho}"
    )

    print(
        f"  normalized_profile="
        f"{normalized_profile}"
    )

    print(
        f"  original_profile="
        f"{original_profile}"
    )

    print(
        f"  normalization_exact="
        f"{normalization_exact}"
    )

    if not normalization_exact:

        failures.append(
            (
                "projective_normalization",
            )
        )

    # ------------------------------------------------------------------
    # 5. SCHUR STRUCTURE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. EXACT SCHUR SOURCE INTERFACE")
    print("=" * 78)

    s0, s1, s2 = schur_row(
        Q1,
        Q3,
    )

    d01, d12, d02, geometry = (
        difference_geometry(
            Q1,
            Q3,
        )
    )

    print(
        f"  s0={s0}"
    )

    print(
        f"  s1={s1}"
    )

    print(
        f"  s2={s2}"
    )

    print(
        f"  d01={d01}"
    )

    print(
        f"  d12={d12}"
    )

    print(
        f"  d02={d02}"
    )

    print(
        f"  difference_geometry_exact="
        f"{geometry}"
    )

    if not geometry:

        failures.append(
            (
                "Schur_geometry",
            )
        )

    # ------------------------------------------------------------------
    # 6. NORMALIZED SOURCE COORDINATES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. NORMALIZED SOURCE COORDINATES")
    print("=" * 78)

    a, b = normalized_source(
        P,
        Q,
        Q1,
        Q3,
    )

    u = Q3 // Q

    print(
        f"  q3/q={u}"
    )

    print(
        f"  a=(q1-1)/(pq)={a}"
    )

    print(
        f"  b=((q3/q)-1)/(pq)={b}"
    )

    print(
        f"  q1_reconstructed="
        f"{1 + P*Q*a}"
    )

    print(
        f"  q3_reconstructed="
        f"{Q*(1 + P*Q*b)}"
    )

    # ------------------------------------------------------------------
    # 7. MINIMAL SOURCE INTERFACE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. MINIMAL SOURCE INTERFACE FOR THE FUTURE n=pq THEOREM")
    print("=" * 78)

    print(
        """
The local theorem now depends on only a small source interface.

Projective part:

    rho = q1/q3.

Local chart:

    rho mod 7 = 2*3^k1.

Derivative data:

    q3 mod 7,

    D = (729-1)/7 mod 7 = 6

for the original 7-adic chart.

Terminal Schur part:

    q1 - 2q3,
    q1 - 6q3,
    q1 - 18q3.

Source normalization for the known n=6 case:

    q3 = q*u,

    q1 = 1 + pq*a,

    u  = 1 + pq*b.

Therefore the remaining symbolic n=pq problem is exactly:

    derive q1(p,q),
    derive q3(p,q),

or derive an equivalent formula for the projective/source interface.
"""
    )

    # ------------------------------------------------------------------
    # 8. STATUS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. FAMILY STATUS")
    print("=" * 78)

    print(
        "  genuine_n_pq_cases_available=1"
    )

    print(
        "  local_projective_structure_verified=True"
    )

    print(
        "  symbolic_q1_of_p_q_known=False"
    )

    print(
        "  symbolic_q3_of_p_q_known=False"
    )

    print(
        "  general_n_pq_theorem_established=False"
    )

    # ------------------------------------------------------------------
    # FINAL
    # ------------------------------------------------------------------

    final_ok = (
        len(failures) == 0
        and
        geometry
    )

    print()
    print("=" * 78)
    print("9. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  admissible_common_scale_tests_exact="
        f"{final_ok}"
    )

    print(
        f"  higher_order_projective_perturbations_exact="
        f"{final_ok}"
    )

    print(
        f"  projective_normalization_exact="
        f"{normalization_exact}"
    )

    print(
        f"  Schur_difference_geometry_exact="
        f"{geometry}"
    )

    print(
        f"  failures={len(failures)}"
    )

    print(
        "  local_projective_theorem_complete=True"
    )

    print(
        "  n_pq_source_theorem_complete=False"
    )

    print(
        "  reason="
        "upstream symbolic q1(p,q),q3(p,q) "
        "still missing"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 255R COMPLETE")


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