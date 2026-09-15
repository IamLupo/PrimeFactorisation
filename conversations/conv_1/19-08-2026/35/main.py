#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 203 — EXACT 7-ADIC INTERCEPT / LIFT-DIGIT RECURRENCE AUDIT
==============================================================================

Experiment 202 established, through three additional levels:

    slopes B = [1, 1, 1]

and

    intercepts A = [0, 2, 3]

with lift digits

    t = [0, 5, 4].

Because

    R(t) = A + t  (mod 7),

the successful digit always satisfies

    t == -A (mod 7).

Experiment 203 asks the next structural question:

    Does the intercept sequence A itself obey a simple recurrence?

The experiment continues the exact lift from the verified exponent

    k = 3401371

and records at each new level:

    * current exponent;
    * order ord_{7^e}(3);
    * normalized residual intercept A_e;
    * slope B_e;
    * chosen lift digit t_e;
    * the identity t_e == -A_e mod 7;
    * first differences of A;
    * first differences of t;
    * A_{e+1} - A_e;
    * whether A is itself affine, periodic, or a simple transform of t.

No fitting beyond exact finite identities.
No SymPy.
No floating point.
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

END_E = 11


# ============================================================================
# HELPERS
# ============================================================================

def v7(x: int):
    if x == 0:
        return None

    x = abs(x)
    e = 0

    while x % 7 == 0:
        x //= 7
        e += 1

    return e


def orbit_value(k: int) -> int:
    return BASE * (MULT ** k)


def orbit_value_mod(
    k: int,
    modulus: int,
) -> int:
    return (
        BASE
        * pow(MULT, k, modulus)
    ) % modulus


def inverse_mod(
    a: int,
    m: int,
) -> int:

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


def ratio_mod(
    a: int,
    b: int,
    modulus: int,
) -> int:

    if modulus == 1:
        return 0

    return (
        a
        * inverse_mod(
            b,
            modulus,
        )
    ) % modulus


def multiplicative_order(
    a: int,
    modulus: int,
) -> int:

    if modulus == 1:
        return 1

    a %= modulus

    if a == 0:
        raise ArithmeticError(
            f"{a} is not a unit modulo {modulus}."
        )

    x = a

    for n in range(1, 100_000_000):

        if x == 1:
            return n

        x = (
            x * a
        ) % modulus

    raise ArithmeticError(
        f"Could not determine order of {a} modulo {modulus}."
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 203 — EXACT 7-ADIC INTERCEPT / "
        "LIFT-DIGIT RECURRENCE AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. SOURCE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. SOURCE DATA")
    print("=" * 78)

    print(f"  q1={Q1}")
    print(f"  q3={Q3}")
    print(f"  v7(q1)={v7(Q1)}")
    print(f"  v7(q3)={v7(Q3)}")

    # ------------------------------------------------------------------
    # 2. STARTING LEVEL
    # ------------------------------------------------------------------

    start_F = (
        orbit_value(START_K)
        * Q3
        - Q1
    )

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
        f"  v7(F(k))={v7(start_F)}"
    )

    print(
        f"  solves_mod_7^{START_E}="
        f"{start_F % (P ** START_E) == 0}"
    )

    # ------------------------------------------------------------------
    # 3. ITERATIVE LIFT
    # ------------------------------------------------------------------

    current_e = START_E
    current_k = START_K

    rows = []

    while current_e < END_E:

        current_modulus = P ** current_e
        next_modulus = P ** (current_e + 1)

        order = multiplicative_order(
            MULT,
            current_modulus,
        )

        rho = ratio_mod(
            Q1,
            Q3,
            next_modulus,
        )

        candidates = []

        for t in range(7):

            candidate_k = (
                current_k
                + t * order
            )

            F = (
                orbit_value(candidate_k)
                * Q3
                - Q1
            )

            if (
                F % current_modulus
                != 0
            ):
                raise ArithmeticError(
                    f"Candidate k={candidate_k} "
                    f"does not satisfy current level."
                )

            A_t = (
                F // current_modulus
            ) % P

            target_match = (
                orbit_value_mod(
                    candidate_k,
                    next_modulus,
                )
                == rho
            )

            candidates.append(
                {
                    "t": t,
                    "k": candidate_k,
                    "F": F,
                    "A": A_t,
                    "target": target_match,
                    "v7": v7(F),
                }
            )

        successful = [
            x
            for x in candidates
            if x["target"]
        ]

        if len(successful) != 1:
            raise ArithmeticError(
                f"Expected one successful lift at e={current_e}, "
                f"got {successful}."
            )

        chosen = successful[0]

        A = candidates[0]["A"]

        B = (
            candidates[1]["A"]
            - candidates[0]["A"]
        ) % P

        affine_exact = all(
            (
                row["A"]
                == (
                    A
                    + B * row["t"]
                ) % P
            )
            for row in candidates
        )

        digit_relation = (
            chosen["t"]
            == (-A) % P
        )

        extra_depth = (
            chosen["v7"]
            - (current_e + 1)
        )

        rows.append(
            {
                "e": current_e,
                "k": current_k,
                "order": order,
                "A": A,
                "B": B,
                "affine": affine_exact,
                "chosen_t": chosen["t"],
                "digit_relation": digit_relation,
                "next_k": chosen["k"],
                "next_v7": chosen["v7"],
                "extra_depth": extra_depth,
            }
        )

        current_k = chosen["k"]
        current_e += 1

    # ------------------------------------------------------------------
    # 4. LEVEL TABLE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. EXACT INTERCEPT / DIGIT TABLE")
    print("=" * 78)

    for row in rows:

        print(
            f"  e={row['e']}: "
            f"k={row['k']} "
            f"order={row['order']} "
            f"A={row['A']} "
            f"B={row['B']} "
            f"t={row['chosen_t']} "
            f"t=-A mod7={row['digit_relation']} "
            f"next_k={row['next_k']} "
            f"next_v7={row['next_v7']} "
            f"extra_depth={row['extra_depth']}"
        )

    # ------------------------------------------------------------------
    # 5. SEQUENCE EXTRACTION
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

    v_seq = [
        row["next_v7"]
        for row in rows
    ]

    print()
    print("=" * 78)
    print("4. GLOBAL SEQUENCES")
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
        f"  residual_v7_sequence={v_seq}"
    )

    # ------------------------------------------------------------------
    # 6. DIRECT DIGIT IDENTITY
    # ------------------------------------------------------------------

    digit_identity = all(
        row["digit_relation"]
        for row in rows
    )

    print()
    print("=" * 78)
    print("5. DIRECT INTERCEPT -> DIGIT IDENTITY")
    print("=" * 78)

    print(
        f"  t_e == -A_e (mod7) "
        f"for every level = "
        f"{digit_identity}"
    )

    for row in rows:

        rhs = (
            -row["A"]
        ) % P

        print(
            f"  e={row['e']}: "
            f"A={row['A']} "
            f"-A mod7={rhs} "
            f"t={row['chosen_t']} "
            f"exact={rhs == row['chosen_t']}"
        )

    # ------------------------------------------------------------------
    # 7. A-DIFFERENCE SEQUENCE
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
    print("6. FIRST DIFFERENCE AUDIT")
    print("=" * 78)

    print(
        f"  A_differences_mod7={A_diffs}"
    )

    print(
        f"  t_differences_mod7={t_diffs}"
    )

    # ------------------------------------------------------------------
    # 8. TEST SIMPLE RELATIONS
    # ------------------------------------------------------------------

    A_equals_neg_t = all(
        (
            A_seq[i]
            + t_seq[i]
        ) % P == 0
        for i in range(
            len(A_seq)
        )
    )

    constant_B = (
        len(set(B_seq)) == 1
    )

    B_one = all(
        B == 1
        for B in B_seq
    )

    affine_all = all(
        row["affine"]
        for row in rows
    )

    print()
    print("=" * 78)
    print("7. SIMPLE RECURRENCE TESTS")
    print("=" * 78)

    print(
        f"  A+t == 0 mod7 everywhere="
        f"{A_equals_neg_t}"
    )

    print(
        f"  B_constant={constant_B}"
    )

    print(
        f"  B_all_one={B_one}"
    )

    print(
        f"  affine_every_level={affine_all}"
    )

    # ------------------------------------------------------------------
    # 9. OVERLIFT PROFILE
    # ------------------------------------------------------------------

    overlifts = [
        row
        for row in rows
        if row["extra_depth"] > 0
    ]

    print()
    print("=" * 78)
    print("8. OVERLIFT PROFILE")
    print("=" * 78)

    print(
        f"  overlift_levels="
        f"{[row['e'] + 1 for row in overlifts]}"
    )

    print(
        f"  overlift_count="
        f"{len(overlifts)}"
    )

    # ------------------------------------------------------------------
    # 10. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The experiments now show a persistent slope-one residual law:

    R_e(t) = A_e + t (mod 7).

Therefore the lift digit is forced by

    t_e = -A_e (mod 7).

Experiment 203 asks whether the intercept sequence itself has a simpler
recurrence.

The data are kept as exact residues in F_7.

A constant or periodic relation in A would be meaningful because it
would reduce the exponent lifting to a small finite-state rule.

A failure of the simple tests would mean that the stable slope-one
phenomenon is real, but the intercepts still carry nontrivial
level-dependent information.

No general recurrence is asserted from the finite sequence alone.
"""
    )

    # ------------------------------------------------------------------
    # 11. FINAL
    # ------------------------------------------------------------------

    expected_minimum_levels = (
        len(rows) > 0
    )

    final_ok = (
        v7(Q1) == 0
        and v7(Q3) == 0
        and expected_minimum_levels
        and digit_identity
        and B_one
        and affine_all
    )

    print()
    print("=" * 78)
    print("10. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  direct_digit_identity_exact="
        f"{digit_identity}"
    )

    print(
        f"  all_slopes_one="
        f"{B_one}"
    )

    print(
        f"  affine_law_exact="
        f"{affine_all}"
    )

    print(
        f"  A_sequence="
        f"{A_seq}"
    )

    print(
        f"  lift_digits="
        f"{t_seq}"
    )

    print(
        f"  exponent_sequence="
        f"{k_seq}"
    )

    print(
        f"  A_differences="
        f"{A_diffs}"
    )

    print(
        f"  t_differences="
        f"{t_diffs}"
    )

    print(
        f"  overlift_count="
        f"{len(overlifts)}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 203 COMPLETE")


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

