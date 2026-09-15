#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 202 — EXACT LONG 7-ADIC LIFT / SLOPE-STABILITY AUDIT
==============================================================================

Established through Experiment 201:

    k1 = 1
    k2 = 1
    k3 = 85
    k4 = 1555
    k5 = 73585

with lift digits

    t = [0, 2, 5, 5].

At each successful lift, the normalized residual over the seven
candidate exponents has been exactly affine:

    R(t) = A + B*t (mod 7),

and so far

    B = 1

at every tested level.

Experiment 202 continues the Hensel lift for several more powers of 7
and records:

    * the unique exponent k_e;
    * the lift digit t_e;
    * the affine coefficients A_e, B_e;
    * whether B_e remains 1;
    * the residual valuation of the successful exponent;
    * whether any new over-lifts occur.

The experiment deliberately derives every next exponent from the
previous level. It does not hard-code future exponents.

No floating point.
No SymPy.
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

# Start from the last verified exponent.
START_E = 5
START_K = 73585

# Continue several additional levels.
END_E = 8


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
        "EXPERIMENT 202 — EXACT LONG 7-ADIC LIFT / "
        "SLOPE-STABILITY AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. SOURCE
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
        f"  v7(q1)={v7(Q1)}"
    )

    print(
        f"  v7(q3)={v7(Q3)}"
    )

    # ------------------------------------------------------------------
    # 2. STARTING LEVEL
    # ------------------------------------------------------------------

    start_modulus = P ** START_E

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
        f"  modulus={start_modulus}"
    )

    print(
        f"  k={START_K}"
    )

    print(
        f"  v7(F(k))={v7(start_F)}"
    )

    print(
        f"  solves_mod_7^5="
        f"{start_F % (P ** 5) == 0}"
    )

    print(
        f"  solves_mod_7^6="
        f"{start_F % (P ** 6) == 0}"
    )

    print(
        f"  solves_mod_7^7="
        f"{start_F % (P ** 7) == 0}"
    )

    # ------------------------------------------------------------------
    # 3. ITERATIVE HENSEL LIFT
    # ------------------------------------------------------------------

    current_e = START_E
    current_k = START_K

    lift_rows = []

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

        print()
        print("=" * 78)
        print(
            f"3. LIFT e={current_e} -> e={current_e + 1}"
        )
        print("=" * 78)

        print(
            f"  current_k={current_k}"
        )

        print(
            f"  current_order={order}"
        )

        print(
            f"  current_modulus={current_modulus}"
        )

        print(
            f"  next_modulus={next_modulus}"
        )

        candidate_rows = []

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
                    f"does not solve current modulus."
                )

            normalized = (
                (F // current_modulus)
                % P
            )

            target_match = (
                orbit_value_mod(
                    candidate_k,
                    next_modulus,
                )
                == rho
            )

            candidate_rows.append(
                {
                    "t": t,
                    "k": candidate_k,
                    "F": F,
                    "normalized": normalized,
                    "target_match": target_match,
                    "v7": v7(F),
                }
            )

            print(
                f"  t={t}: "
                f"k={candidate_k} "
                f"normalized={normalized} "
                f"v7(F)={v7(F)} "
                f"target_match={target_match}"
            )

        successful = [
            row
            for row in candidate_rows
            if row["target_match"]
        ]

        if len(successful) != 1:
            raise ArithmeticError(
                f"Expected exactly one successful lift, "
                f"got {successful}."
            )

        chosen = successful[0]

        # --------------------------------------------------------------
        # 4. AFFINE LAW
        # --------------------------------------------------------------

        A = candidate_rows[0]["normalized"]

        B = (
            candidate_rows[1]["normalized"]
            - candidate_rows[0]["normalized"]
        ) % P

        affine_exact = True

        for row in candidate_rows:

            predicted = (
                A
                + B * row["t"]
            ) % P

            if predicted != row["normalized"]:
                affine_exact = False

        # --------------------------------------------------------------
        # 5. SUCCESSFUL LIFT DATA
        # --------------------------------------------------------------

        next_k = chosen["k"]

        next_F = chosen["F"]
        next_v7 = chosen["v7"]

        extra_depth = (
            next_v7
            - (current_e + 1)
            if next_v7 is not None
            else None
        )

        lift_rows.append(
            {
                "from_e": current_e,
                "to_e": current_e + 1,
                "current_k": current_k,
                "order": order,
                "A": A,
                "B": B,
                "affine_exact": affine_exact,
                "chosen_t": chosen["t"],
                "next_k": next_k,
                "next_v7": next_v7,
                "extra_depth": extra_depth,
            }
        )

        print()
        print(
            f"  affine_A={A}"
        )

        print(
            f"  affine_B={B}"
        )

        print(
            f"  affine_exact={affine_exact}"
        )

        print(
            f"  chosen_t={chosen['t']}"
        )

        print(
            f"  next_k={next_k}"
        )

        print(
            f"  next_v7(F)={next_v7}"
        )

        print(
            f"  extra_depth={extra_depth}"
        )

        current_e += 1
        current_k = next_k

    # ------------------------------------------------------------------
    # 6. GLOBAL SLOPE PROFILE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. GLOBAL LIFT PROFILE")
    print("=" * 78)

    print(
        f"  exponent_start={START_K}"
    )

    print(
        f"  final_exponent={current_k}"
    )

    print(
        f"  lift_digits="
        f"{[row['chosen_t'] for row in lift_rows]}"
    )

    print(
        f"  slopes="
        f"{[row['B'] for row in lift_rows]}"
    )

    print(
        f"  intercepts="
        f"{[row['A'] for row in lift_rows]}"
    )

    print(
        f"  affine_flags="
        f"{[row['affine_exact'] for row in lift_rows]}"
    )

    print(
        f"  residual_v7="
        f"{[row['next_v7'] for row in lift_rows]}"
    )

    print(
        f"  extra_depth="
        f"{[row['extra_depth'] for row in lift_rows]}"
    )

    # ------------------------------------------------------------------
    # 7. SLOPE-STABILITY
    # ------------------------------------------------------------------

    slopes = [
        row["B"]
        for row in lift_rows
    ]

    affine_flags = [
        row["affine_exact"]
        for row in lift_rows
    ]

    slope_one = all(
        B == 1
        for B in slopes
    )

    affine_all = all(
        affine_flags
    )

    print()
    print("=" * 78)
    print("5. SLOPE-STABILITY")
    print("=" * 78)

    print(
        f"  slopes={slopes}"
    )

    print(
        f"  all_slopes_equal_1={slope_one}"
    )

    print(
        f"  all_affine_exact={affine_all}"
    )

    # ------------------------------------------------------------------
    # 8. NEW OVERLIFT AUDIT
    # ------------------------------------------------------------------

    overlifts = [
        row
        for row in lift_rows
        if row["extra_depth"] is not None
        and row["extra_depth"] > 0
    ]

    print()
    print("=" * 78)
    print("6. NEW OVERLIFT AUDIT")
    print("=" * 78)

    print(
        f"  overlift_levels="
        f"{[row['to_e'] for row in overlifts]}"
    )

    print(
        f"  overlift_count={len(overlifts)}"
    )

    # ------------------------------------------------------------------
    # 9. FINAL STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The lift is now continued beyond the exceptional e=3 -> e=4
over-lift.

At every new level the seven candidates are tested and the normalized
residual is fitted exactly by

    R(t) = A + B*t mod 7.

The central question is whether

    B = 1

persists.

If so, the exponent lift has a stable unit slope through all tested
levels, rather than merely producing slope one in the first few
examples.

The experiment also records any further over-lifts separately. A
future over-lift would be a zero next residual digit, not a failure of
the lifting mechanism.

This remains a finite exact observation about the specific source ratio.
"""
    )

    # ------------------------------------------------------------------
    # 10. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        v7(Q1) == 0
        and v7(Q3) == 0
        and affine_all
        and slope_one
        and len(overlifts) >= 0
    )

    print()
    print("=" * 78)
    print("8. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  all_affine_exact={affine_all}"
    )

    print(
        f"  all_slopes_equal_1={slope_one}"
    )

    print(
        f"  lift_digits="
        f"{[row['chosen_t'] for row in lift_rows]}"
    )

    print(
        f"  final_exponent={current_k}"
    )

    print(
        f"  overlift_count={len(overlifts)}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 202 COMPLETE")


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

