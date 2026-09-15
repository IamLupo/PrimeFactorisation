#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 180 — EXACT TERMINAL SCHUR INVERSE / CONTENT-MAP AUDIT
==============================================================================

The terminal Schur coefficients are

    s0 = q1 - 2 q3
    s1 = q1 - 6 q3
    s2 = q1 - 18 q3.

Correct inverse identities are

    q3 = (s0 - s1) / 4
       = (s1 - s2) / 12,

and

    q1 = s0 + 2 q3
       = (3 s1 - s2) / 2.

Experiment 179 used incorrect signs/factors in part of the inverse
calculation. This experiment corrects that and then studies the exact
integer-content map.

For

    (q1,q3) -> (s0,s1,s2),

the coefficient matrix is

    T = [[1,-2],
         [1,-6],
         [1,-18]].

Its 2x2 minors are

    -4, -16, -12,

whose gcd is 4.

Therefore the only possible new common content introduced by the
linear map is constrained by this determinantal factor.

Experiment 180 checks, for the actual source pair:

    * exact inverse recovery;
    * primitive inverse recovery;
    * gcd(q1,q3);
    * gcd(s0,s1,s2);
    * whether the quotient gcd(s)/gcd(q) is 1 or a power of 2;
    * exact parity criterion;
    * all 2x2 minors of T;
    * exact reconstruction after dividing the source by its gcd.

No floating point.
No SymPy.
No extrapolation.
"""

from __future__ import annotations

import sys
from math import gcd


# ============================================================================
# DATA
# ============================================================================

Q = {
    1: [
        -126258696,
        -11600759760,
        2668721436,
        1764373740,
        -1338089411,
        495451247,
    ],
    3: [
        9955176,
        -1263551016,
        -152369292,
        -128667196,
        421514439,
    ],
    5: [
        -62398,
        4771718,
        16027881,
    ],
    7: [
        1,
    ],
}


# ============================================================================
# HELPERS
# ============================================================================

def gcd_list(values):
    g = 0
    for x in values:
        g = gcd(g, abs(x))
    return g


def v2(x):
    if x == 0:
        return None

    x = abs(x)
    e = 0

    while x % 2 == 0:
        x //= 2
        e += 1

    return e


def q_terminal(p):
    return Q[p][-1]


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 180 — EXACT TERMINAL SCHUR INVERSE / "
        "CONTENT-MAP AUDIT"
    )
    print("=" * 78)

    q1 = q_terminal(1)
    q3 = q_terminal(3)

    # ------------------------------------------------------------------
    # 1. SOURCE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. EXACT SOURCE PAIR")
    print("=" * 78)

    print(f"  q1={q1}")
    print(f"  q3={q3}")

    source_gcd = gcd(
        abs(q1),
        abs(q3),
    )

    print(f"  gcd(q1,q3)={source_gcd}")
    print(f"  v2(gcd(q1,q3))={v2(source_gcd)}")

    # ------------------------------------------------------------------
    # 2. FORWARD MAP
    # ------------------------------------------------------------------

    s0 = q1 - 2 * q3
    s1 = q1 - 6 * q3
    s2 = q1 - 18 * q3

    s = [s0, s1, s2]

    print()
    print("=" * 78)
    print("2. FORWARD SCHUR MAP")
    print("=" * 78)

    print(f"  s0=q1-2*q3={s0}")
    print(f"  s1=q1-6*q3={s1}")
    print(f"  s2=q1-18*q3={s2}")

    schur_gcd = gcd_list(s)

    print(f"  gcd(s0,s1,s2)={schur_gcd}")
    print(f"  v2(gcd(s))={v2(schur_gcd)}")

    # ------------------------------------------------------------------
    # 3. CORRECT q3 INVERSION
    # ------------------------------------------------------------------

    num_q3_a = s0 - s1
    num_q3_b = s1 - s2

    q3_a_exact = (
        num_q3_a % 4 == 0
    )

    q3_b_exact = (
        num_q3_b % 12 == 0
    )

    q3_a = (
        num_q3_a // 4
        if q3_a_exact
        else None
    )

    q3_b = (
        num_q3_b // 12
        if q3_b_exact
        else None
    )

    q3_recovery_exact = (
        q3_a == q3
        and q3_b == q3
    )

    print()
    print("=" * 78)
    print("3. CORRECT q3 RECOVERY")
    print("=" * 78)

    print(
        f"  (s0-s1)/4={q3_a}"
    )

    print(
        f"  (s1-s2)/12={q3_b}"
    )

    print(
        f"  first_division_exact={q3_a_exact}"
    )

    print(
        f"  second_division_exact={q3_b_exact}"
    )

    print(
        f"  exact_q3_recovery={q3_recovery_exact}"
    )

    # ------------------------------------------------------------------
    # 4. CORRECT q1 INVERSION
    # ------------------------------------------------------------------

    q1_from_s0 = s0 + 2 * q3_a if q3_a is not None else None

    num_q1_b = 3 * s1 - s2

    q1_b_exact = (
        num_q1_b % 2 == 0
    )

    q1_from_s1 = (
        num_q1_b // 2
        if q1_b_exact
        else None
    )

    q1_recovery_exact = (
        q1_from_s0 == q1
        and q1_from_s1 == q1
    )

    print()
    print("=" * 78)
    print("4. CORRECT q1 RECOVERY")
    print("=" * 78)

    print(
        f"  s0+2*q3={q1_from_s0}"
    )

    print(
        f"  (3*s1-s2)/2={q1_from_s1}"
    )

    print(
        f"  second_division_exact={q1_b_exact}"
    )

    print(
        f"  exact_q1_recovery={q1_recovery_exact}"
    )

    # ------------------------------------------------------------------
    # 5. PRIMITIVE SOURCE INVERSION
    # ------------------------------------------------------------------

    q1_primitive = q1 // source_gcd
    q3_primitive = q3 // source_gcd

    s_primitive = [
        x // source_gcd
        for x in s
    ]

    recovered_primitive_q3_a = (
        (s_primitive[0] - s_primitive[1]) // 4
    )

    recovered_primitive_q3_b = (
        (s_primitive[1] - s_primitive[2]) // 12
    )

    recovered_primitive_q1_a = (
        s_primitive[0]
        + 2 * recovered_primitive_q3_a
    )

    recovered_primitive_q1_b = (
        (3 * s_primitive[1] - s_primitive[2]) // 2
    )

    primitive_recovery_exact = (
        recovered_primitive_q1_a == q1_primitive
        and recovered_primitive_q1_b == q1_primitive
        and recovered_primitive_q3_a == q3_primitive
        and recovered_primitive_q3_b == q3_primitive
    )

    print()
    print("=" * 78)
    print("5. PRIMITIVE SOURCE RECOVERY")
    print("=" * 78)

    print(
        f"  q1_primitive={q1_primitive}"
    )

    print(
        f"  q3_primitive={q3_primitive}"
    )

    print(
        f"  s_primitive={s_primitive}"
    )

    print(
        f"  recovered_q1_a={recovered_primitive_q1_a}"
    )

    print(
        f"  recovered_q1_b={recovered_primitive_q1_b}"
    )

    print(
        f"  recovered_q3_a={recovered_primitive_q3_a}"
    )

    print(
        f"  recovered_q3_b={recovered_primitive_q3_b}"
    )

    print(
        f"  primitive_recovery_exact="
        f"{primitive_recovery_exact}"
    )

    # ------------------------------------------------------------------
    # 6. TRANSFORMATION MATRIX MINORS
    # ------------------------------------------------------------------

    T = [
        [1, -2],
        [1, -6],
        [1, -18],
    ]

    minors = []

    for i in range(3):
        for j in range(i + 1, 3):

            determinant = (
                T[i][0] * T[j][1]
                - T[i][1] * T[j][0]
            )

            minors.append(
                determinant
            )

    minors_gcd = gcd_list(minors)

    print()
    print("=" * 78)
    print("6. EXACT TRANSFORMATION-MATRIX MINORS")
    print("=" * 78)

    print(
        f"  T={T}"
    )

    print(
        f"  2x2_minors={minors}"
    )

    print(
        f"  gcd_of_2x2_minors={minors_gcd}"
    )

    # ------------------------------------------------------------------
    # 7. CONTENT RATIO
    # ------------------------------------------------------------------

    if source_gcd != 0:

        content_ratio_exact = (
            schur_gcd % source_gcd == 0
        )

        content_ratio = (
            schur_gcd // source_gcd
            if content_ratio_exact
            else None
        )

    else:

        content_ratio_exact = False
        content_ratio = None

    print()
    print("=" * 78)
    print("7. CONTENT RATIO")
    print("=" * 78)

    print(
        f"  gcd_source={source_gcd}"
    )

    print(
        f"  gcd_schur={schur_gcd}"
    )

    print(
        f"  gcd_schur/gcd_source={content_ratio}"
    )

    if content_ratio is not None:

        print(
            f"  v2(content_ratio)="
            f"{v2(content_ratio)}"
        )

    # ------------------------------------------------------------------
    # 8. PARITY CRITERION
    # ------------------------------------------------------------------

    source_primitive_both_odd = (
        q1_primitive % 2 == 1
        and q3_primitive % 2 == 1
    )

    schur_primitive_all_odd = all(
        x % 2 == 1
        for x in s_primitive
    )

    print()
    print("=" * 78)
    print("8. PARITY / CONTENT CRITERION")
    print("=" * 78)

    print(
        f"  primitive_q1_odd="
        f"{q1_primitive % 2 == 1}"
    )

    print(
        f"  primitive_q3_odd="
        f"{q3_primitive % 2 == 1}"
    )

    print(
        f"  primitive_source_both_odd="
        f"{source_primitive_both_odd}"
    )

    print(
        f"  primitive_schur_all_odd="
        f"{schur_primitive_all_odd}"
    )

    # ------------------------------------------------------------------
    # 9. EXACT RECONSTRUCTION
    # ------------------------------------------------------------------

    reconstructed_s = [
        q1 - 2 * q3,
        q1 - 6 * q3,
        q1 - 18 * q3,
    ]

    forward_reconstruction_exact = (
        reconstructed_s == s
    )

    print()
    print("=" * 78)
    print("9. FORWARD / INVERSE CONSISTENCY")
    print("=" * 78)

    print(
        f"  forward_exact="
        f"{forward_reconstruction_exact}"
    )

    print(
        f"  inverse_q1_exact="
        f"{q1_recovery_exact}"
    )

    print(
        f"  inverse_q3_exact="
        f"{q3_recovery_exact}"
    )

    # ------------------------------------------------------------------
    # 10. INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
Experiment 179's numerical data were correct; only its inverse
formulas were coded incorrectly.

The exact terminal map is

    (q1,q3)
        ->
    (q1-2q3, q1-6q3, q1-18q3).

Its inverse can be recovered using any two adjacent Schur coefficients:

    q3 = (s0-s1)/4
       = (s1-s2)/12,

and

    q1 = s0 + 2q3
       = (3s1-s2)/2.

This experiment verifies those identities directly.

The transformation matrix has 2x2-minor gcd 4. That is the exact
determinantal content controlling possible divisibility introduced by
the linear transformation.

For the actual primitive source pair here, both source entries are odd,
and the three primitive Schur coefficients are also odd. Therefore no
additional factor of 2 appears in the observed content.

The exact common factor 17 is therefore inherited from the source pair
for this dataset.

This is an arithmetic statement about the observed terminal map, not
a general theorem about other p-values.
"""
    )

    # ------------------------------------------------------------------
    # 11. FINAL
    # ------------------------------------------------------------------

    content_match = (
        source_gcd == schur_gcd
    )

    final_ok = (
        q3_recovery_exact
        and q1_recovery_exact
        and primitive_recovery_exact
        and forward_reconstruction_exact
        and minors_gcd == 4
        and content_match
        and source_primitive_both_odd
        and schur_primitive_all_odd
    )

    print()
    print("=" * 78)
    print("11. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  exact_q3_recovery={q3_recovery_exact}"
    )

    print(
        f"  exact_q1_recovery={q1_recovery_exact}"
    )

    print(
        f"  primitive_recovery_exact="
        f"{primitive_recovery_exact}"
    )

    print(
        f"  forward_reconstruction_exact="
        f"{forward_reconstruction_exact}"
    )

    print(
        f"  transformation_minor_gcd_is_4="
        f"{minors_gcd == 4}"
    )

    print(
        f"  source_schur_content_match="
        f"{content_match}"
    )

    print(
        f"  primitive_source_both_odd="
        f"{source_primitive_both_odd}"
    )

    print(
        f"  primitive_schur_all_odd="
        f"{schur_primitive_all_odd}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 180 COMPLETE")


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

