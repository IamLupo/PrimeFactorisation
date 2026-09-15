#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 203FR — FAST EXACT 7-ADIC INTERCEPT / LIFT-DIGIT AUDIT
==============================================================================

Corrected Experiment 203F.

The previous script failed because valuation_7() called an undefined v7()
helper. This version has a single self-contained valuation_7() function.

The computation stays entirely modular. It never constructs 3^k as a huge
integer when testing candidate lifts.

Verified starting point:

    e = 8
    k = 3401371

The lift continues through e=12 unless interrupted.
"""

from __future__ import annotations

import sys


# ============================================================================
# DATA
# ============================================================================

Q1 = 29144191
Q3 = 24794967

BASE = 2
MULT = 3
P = 7

START_E = 8
START_K = 3401371
END_E = 12


# ============================================================================
# EXACT 7-ADIC VALUATION
# ============================================================================

def valuation_7(x: int):
    """
    Exact v_7(x).

    Returns None for x=0.
    """
    if x == 0:
        return None

    x = abs(x)
    e = 0

    while x % P == 0:
        x //= P
        e += 1

    return e


# ============================================================================
# FAST MODULAR OPERATIONS
# ============================================================================

def orbit_mod(k: int, modulus: int) -> int:
    """
    Exact value

        2 * 3^k mod modulus

    without constructing 3^k.
    """
    return (
        BASE
        * pow(MULT, k, modulus)
    ) % modulus


def F_mod(k: int, modulus: int) -> int:
    """
    Exact

        F(k) = 2*3^k*q3 - q1

    modulo modulus.
    """
    return (
        orbit_mod(k, modulus)
        * Q3
        - Q1
    ) % modulus


def residual_digit(k: int, e: int) -> int:
    """
    Compute

        F(k) / 7^e mod 7

    assuming F(k) is divisible by 7^e.

    Only F modulo 7^(e+1) is required.
    """
    scale = P ** e
    modulus = P ** (e + 1)

    residue = F_mod(
        k,
        modulus,
    )

    if residue % scale != 0:
        raise ArithmeticError(
            f"k={k} does not satisfy "
            f"F(k)=0 mod 7^{e}."
        )

    return (
        residue // scale
    ) % P


def solves_level(k: int, e: int) -> bool:
    """
    Check whether

        F(k) == 0 mod 7^e.
    """
    modulus = P ** e

    return (
        F_mod(k, modulus)
        == 0
    )


def solves_next_level(k: int, e: int) -> bool:
    """
    Check whether

        F(k) == 0 mod 7^(e+1).
    """
    return solves_level(
        k,
        e + 1,
    )


def order_3_mod_7e(e: int) -> int:
    """
    For 3 modulo 7^e:

        ord_{7^e}(3) = 6 * 7^(e-1).

    This was independently established in the earlier experiments.
    """
    if e <= 0:
        raise ValueError(
            "e must be positive."
        )

    return 6 * (
        P ** (e - 1)
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 203FR — FAST EXACT 7-ADIC "
        "INTERCEPT / LIFT-DIGIT AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. SOURCE DATA
    # ------------------------------------------------------------------

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
        f"  v7(q1)={valuation_7(Q1)}"
    )

    print(
        f"  v7(q3)={valuation_7(Q3)}"
    )

    # ------------------------------------------------------------------
    # 2. STARTING VERIFIED LEVEL
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. STARTING VERIFIED LEVEL")
    print("=" * 78)

    print(
        f"  e={START_E}"
    )

    print(
        f"  k={START_K}"
    )

    print(
        f"  solves_mod_7^{START_E}="
        f"{solves_level(START_K, START_E)}"
    )

    print(
        f"  solves_mod_7^{START_E + 1}="
        f"{solves_level(START_K, START_E + 1)}"
    )

    # ------------------------------------------------------------------
    # 3. FAST ITERATIVE LIFT
    # ------------------------------------------------------------------

    current_e = START_E
    current_k = START_K

    rows = []

    while current_e < END_E:

        current_order = order_3_mod_7e(
            current_e
        )

        current_modulus = P ** current_e
        next_modulus = P ** (
            current_e + 1
        )

        print()
        print("=" * 78)
        print(
            f"3. FAST LIFT e={current_e} "
            f"-> e={current_e + 1}"
        )
        print("=" * 78)

        print(
            f"  current_k={current_k}"
        )

        print(
            f"  current_order={current_order}"
        )

        print(
            f"  current_modulus={current_modulus}"
        )

        print(
            f"  next_modulus={next_modulus}"
        )

        candidates = []

        # --------------------------------------------------------------
        # Seven possible Hensel digits
        # --------------------------------------------------------------

        for t in range(7):

            candidate_k = (
                current_k
                + t * current_order
            )

            digit = residual_digit(
                candidate_k,
                current_e,
            )

            target_match = solves_next_level(
                candidate_k,
                current_e,
            )

            candidates.append(
                {
                    "t": t,
                    "k": candidate_k,
                    "digit": digit,
                    "target": target_match,
                }
            )

            print(
                f"  t={t}: "
                f"k={candidate_k} "
                f"normalized={digit} "
                f"target_match={target_match}"
            )

        successful = [
            row
            for row in candidates
            if row["target"]
        ]

        if len(successful) != 1:
            raise ArithmeticError(
                f"Expected exactly one successful lift "
                f"at e={current_e}; "
                f"got {successful}."
            )

        chosen = successful[0]

        # --------------------------------------------------------------
        # Affine residual law
        # --------------------------------------------------------------

        A = candidates[0]["digit"]

        B = (
            candidates[1]["digit"]
            - candidates[0]["digit"]
        ) % P

        affine_exact = True

        for row in candidates:

            predicted = (
                A
                + B * row["t"]
            ) % P

            if predicted != row["digit"]:
                affine_exact = False
                break

        # --------------------------------------------------------------
        # Over-lift check
        # --------------------------------------------------------------

        next_k = chosen["k"]

        extra_digit = None
        over_lift = False

        if solves_level(
            next_k,
            current_e + 1,
        ):

            extra_digit = residual_digit(
                next_k,
                current_e + 1,
            )

            over_lift = (
                extra_digit == 0
            )

        rows.append(
            {
                "e": current_e,
                "k": current_k,
                "order": current_order,
                "A": A,
                "B": B,
                "chosen_t": chosen["t"],
                "next_k": next_k,
                "affine_exact": affine_exact,
                "extra_digit": extra_digit,
                "over_lift": over_lift,
            }
        )

        print()
        print(
            f"  A={A}"
        )

        print(
            f"  B={B}"
        )

        print(
            f"  affine_exact={affine_exact}"
        )

        print(
            f"  chosen_t={chosen['t']}"
        )

        print(
            f"  t=-A mod7="
            f"{chosen['t'] == (-A) % P}"
        )

        print(
            f"  next_k={next_k}"
        )

        print(
            f"  next_level_overlift={over_lift}"
        )

        if extra_digit is not None:
            print(
                f"  next_level_extra_digit={extra_digit}"
            )

        current_e += 1
        current_k = next_k

    # ------------------------------------------------------------------
    # 4. GLOBAL PROFILE
    # ------------------------------------------------------------------

    A_seq = [
        row["A"]
        for row in rows
    ]

    B_seq = [
        row["B"]
        for row in rows
    ]

    t_seq = [
        row["chosen_t"]
        for row in rows
    ]

    k_seq = [
        START_K
    ] + [
        row["next_k"]
        for row in rows
    ]

    print()
    print("=" * 78)
    print("4. GLOBAL PROFILE")
    print("=" * 78)

    print(
        f"  A_sequence={A_seq}"
    )

    print(
        f"  B_sequence={B_seq}"
    )

    print(
        f"  lift_digits={t_seq}"
    )

    print(
        f"  exponent_sequence={k_seq}"
    )

    print(
        f"  t_minus_A_checks="
        f"{[t == (-A) % P for t, A in zip(t_seq, A_seq)]}"
    )

    print(
        f"  overlift_flags="
        f"{[row['over_lift'] for row in rows]}"
    )

    # ------------------------------------------------------------------
    # 5. DIFFERENCE AUDIT
    # ------------------------------------------------------------------

    A_diffs = [
        (
            A_seq[i + 1]
            - A_seq[i]
        ) % P
        for i in range(
            len(A_seq) - 1
        )
    ]

    t_diffs = [
        (
            t_seq[i + 1]
            - t_seq[i]
        ) % P
        for i in range(
            len(t_seq) - 1
        )
    ]

    print()
    print("=" * 78)
    print("5. DIFFERENCE AUDIT")
    print("=" * 78)

    print(
        f"  A_differences_mod7={A_diffs}"
    )

    print(
        f"  t_differences_mod7={t_diffs}"
    )

    # ------------------------------------------------------------------
    # 6. STABILITY TESTS
    # ------------------------------------------------------------------

    all_B_one = all(
        B == 1
        for B in B_seq
    )

    all_affine = all(
        row["affine_exact"]
        for row in rows
    )

    all_digit_identity = all(
        t == (-A) % P
        for t, A in zip(
            t_seq,
            A_seq,
        )
    )

    print()
    print("=" * 78)
    print("6. SLOPE / DIGIT STABILITY")
    print("=" * 78)

    print(
        f"  all_B_equal_1={all_B_one}"
    )

    print(
        f"  all_affine_exact={all_affine}"
    )

    print(
        f"  all_t_equals_minus_A="
        f"{all_digit_identity}"
    )

    # ------------------------------------------------------------------
    # 7. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
This is the fast modular continuation of Experiment 203.

No integer of the form 3^3,401,371 is constructed.

At each level only

    2*3^k mod 7^(e+1)

is evaluated.

The exact residual law remains

    R_e(t) = A_e + B_e*t (mod 7).

The main structural test is whether

    B_e = 1

continues across the new levels.

When B_e=1, the next lift digit is forced:

    t_e = -A_e (mod 7).

The intercept sequence A_e is therefore the remaining nontrivial
finite-state data.
"""
    )

    # ------------------------------------------------------------------
    # 8. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        valuation_7(Q1) == 0
        and valuation_7(Q3) == 0
        and len(rows) > 0
        and all_B_one
        and all_affine
        and all_digit_identity
    )

    print()
    print("=" * 78)
    print("8. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  fast_modular_path=True"
    )

    print(
        f"  all_B_equal_1={all_B_one}"
    )

    print(
        f"  all_affine_exact={all_affine}"
    )

    print(
        f"  all_t_equals_minus_A="
        f"{all_digit_identity}"
    )

    print(
        f"  A_sequence={A_seq}"
    )

    print(
        f"  B_sequence={B_seq}"
    )

    print(
        f"  lift_digits={t_seq}"
    )

    print(
        f"  exponent_sequence={k_seq}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 203FR COMPLETE")


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