#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 179 — EXACT TERMINAL SCHUR SOURCE-RECOVERY /
                GEOMETRIC p-COEFFICIENT + CONTENT PROVENANCE AUDIT
==============================================================================

Experiment 178 found:

    s0 = -347577631
    s1 = -2033635387
    s2 = -7091808655

with

    gcd(s0,s1,s2) = 17.

The exact construction is

    s0 = q1 - 2 q3
    s1 = q1 - 6 q3
    s2 = q1 - 18 q3.

Experiment 179 audits the inverse direction.

It asks whether q1 and q3 can be recovered exactly from the three
terminal Schur coefficients alone, and whether the complete factor 17
comes directly from the common content of q1 and q3.

The exact identities being tested include:

    s0 - s1 = 4 q3
    s1 - s2 = 12 q3

hence

    (s1-s2) = 3(s0-s1).

Also:

    q3 = (s0-s1)/4

and

    q1 = (3 s0 - s1)/2.

Independent reconstructions are performed using different pairs.

The experiment additionally computes:

    gcd(q1,q3)
    gcd(s0,s1,s2)

and checks whether the latter is exactly the former.

It also strips the common source content and verifies that the
normalized terminal coefficients retain the same exact geometric law.

No extrapolation is performed beyond the three existing terminal
coefficients.

No SymPy.
No floating point.
No recurrence fitting.
"""

from __future__ import annotations

import sys
from math import gcd


# ============================================================================
# EXACT SOURCE DATA
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

def v2(x):
    if x == 0:
        return None

    x = abs(x)
    e = 0

    while x % 2 == 0:
        x //= 2
        e += 1

    return e


def vp(x, p):
    if x == 0:
        return None

    x = abs(x)
    e = 0

    while x % p == 0:
        x //= p
        e += 1

    return e


def gcd_list(values):
    g = 0

    for x in values:
        g = gcd(g, abs(x))

    return g


def q_terminal(p):
    return Q[p][-1]


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 179 — EXACT TERMINAL SCHUR SOURCE-RECOVERY / "
        "GEOMETRIC p-COEFFICIENT + CONTENT PROVENANCE AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. SOURCE VALUES
    # ------------------------------------------------------------------

    q1 = q_terminal(1)
    q3 = q_terminal(3)

    print()
    print("=" * 78)
    print("1. EXACT TERMINAL SOURCE VALUES")
    print("=" * 78)

    print(
        f"  q1={q1}"
    )

    print(
        f"  q3={q3}"
    )

    print(
        f"  gcd(q1,q3)={gcd(q1, q3)}"
    )

    print(
        f"  v2(q1)={v2(q1)}"
    )

    print(
        f"  v2(q3)={v2(q3)}"
    )

    print(
        f"  v17(q1)={vp(q1, 17)}"
    )

    print(
        f"  v17(q3)={vp(q3, 17)}"
    )

    # ------------------------------------------------------------------
    # 2. TERMINAL SCHUR COEFFICIENTS
    # ------------------------------------------------------------------

    s0 = q1 - 2 * q3
    s1 = q1 - 6 * q3
    s2 = q1 - 18 * q3

    s = [s0, s1, s2]

    print()
    print("=" * 78)
    print("2. EXACT TERMINAL SCHUR COEFFICIENTS")
    print("=" * 78)

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
        f"  gcd(s0,s1,s2)={gcd_list(s)}"
    )

    # ------------------------------------------------------------------
    # 3. EXACT FIRST-DIFFERENCE GEOMETRY
    # ------------------------------------------------------------------

    d0 = s1 - s0
    d1 = s2 - s1

    print()
    print("=" * 78)
    print("3. EXACT DIFFERENCE GEOMETRY")
    print("=" * 78)

    print(
        f"  d0=s1-s0={d0}"
    )

    print(
        f"  d1=s2-s1={d1}"
    )

    print(
        f"  d0/(q3)={d0 // q3 if d0 % q3 == 0 else 'nonintegral'}"
    )

    print(
        f"  d1/(q3)={d1 // q3 if d1 % q3 == 0 else 'nonintegral'}"
    )

    ratio_3_exact = (
        d1 == 3 * d0
    )

    print(
        f"  d1=3*d0={ratio_3_exact}"
    )

    # ------------------------------------------------------------------
    # 4. INVERSE RECOVERY OF q3
    # ------------------------------------------------------------------

    q3_from_01_num = s1 - s0
    q3_from_01_exact = (
        q3_from_01_num % 4 == 0
    )

    q3_from_01 = (
        q3_from_01_num // 4
        if q3_from_01_exact
        else None
    )

    q3_from_12_num = s2 - s1
    q3_from_12_exact = (
        q3_from_12_num % 12 == 0
    )

    q3_from_12 = (
        q3_from_12_num // 12
        if q3_from_12_exact
        else None
    )

    print()
    print("=" * 78)
    print("4. EXACT q3 RECOVERY")
    print("=" * 78)

    print(
        f"  q3_from_s0_s1={q3_from_01}"
    )

    print(
        f"  q3_from_s1_s2={q3_from_12}"
    )

    print(
        f"  recovery_01_exact="
        f"{q3_from_01_exact}"
    )

    print(
        f"  recovery_12_exact="
        f"{q3_from_12_exact}"
    )

    print(
        f"  both_recover_q3="
        f"{q3_from_01 == q3 and q3_from_12 == q3}"
    )

    # ------------------------------------------------------------------
    # 5. INVERSE RECOVERY OF q1
    # ------------------------------------------------------------------

    q1_from_s0s1_num = (
        3 * s0 - s1
    )

    q1_from_s0s1_exact = (
        q1_from_s0s1_num % 2 == 0
    )

    q1_from_s0s1 = (
        q1_from_s0s1_num // 2
        if q1_from_s0s1_exact
        else None
    )

    q1_from_s1s2_num = (
        3 * s1 - s2
    )

    q1_from_s1s2_exact = (
        q1_from_s1s2_num % 4 == 0
    )

    q1_from_s1s2 = (
        q1_from_s1s2_num // 4
        if q1_from_s1s2_exact
        else None
    )

    print()
    print("=" * 78)
    print("5. EXACT q1 RECOVERY")
    print("=" * 78)

    print(
        f"  q1_from_s0_s1={q1_from_s0s1}"
    )

    print(
        f"  q1_from_s1_s2={q1_from_s1s2}"
    )

    print(
        f"  recovery_01_exact="
        f"{q1_from_s0s1_exact}"
    )

    print(
        f"  recovery_12_exact="
        f"{q1_from_s1s2_exact}"
    )

    print(
        f"  both_recover_q1="
        f"{q1_from_s0s1 == q1 and q1_from_s1s2 == q1}"
    )

    # ------------------------------------------------------------------
    # 6. SOURCE-CONTENT PROVENANCE
    # ------------------------------------------------------------------

    source_gcd = gcd(q1, q3)
    schur_gcd = gcd_list(s)

    content_matches = (
        source_gcd == schur_gcd
    )

    print()
    print("=" * 78)
    print("6. COMMON-CONTENT PROVENANCE")
    print("=" * 78)

    print(
        f"  gcd(q1,q3)={source_gcd}"
    )

    print(
        f"  gcd(s0,s1,s2)={schur_gcd}"
    )

    print(
        f"  source_content_equals_schur_content="
        f"{content_matches}"
    )

    if source_gcd != 0:

        print(
            f"  source_gcd_v17={vp(source_gcd,17)}"
        )

        print(
            f"  source_gcd_v2={vp(source_gcd,2)}"
        )

    # ------------------------------------------------------------------
    # 7. REMOVE SOURCE CONTENT
    # ------------------------------------------------------------------

    if source_gcd == 0:
        raise ArithmeticError(
            "Unexpected zero source gcd."
        )

    q1_prim = q1 // source_gcd
    q3_prim = q3 // source_gcd

    s_prim = [
        x // source_gcd
        for x in s
    ]

    primitive_source_gcd = gcd(
        abs(q1_prim),
        abs(q3_prim)
    )

    primitive_schur_gcd = gcd_list(
        s_prim
    )

    print()
    print("=" * 78)
    print("7. SOURCE-NORMALIZED SCHUR ROW")
    print("=" * 78)

    print(
        f"  q1_primitive={q1_prim}"
    )

    print(
        f"  q3_primitive={q3_prim}"
    )

    print(
        f"  primitive_source_gcd="
        f"{primitive_source_gcd}"
    )

    print(
        f"  primitive_schur={s_prim}"
    )

    print(
        f"  primitive_schur_gcd="
        f"{primitive_schur_gcd}"
    )

    # ------------------------------------------------------------------
    # 8. NORMALIZED DIFFERENCE GEOMETRY
    # ------------------------------------------------------------------

    pd0 = s_prim[1] - s_prim[0]
    pd1 = s_prim[2] - s_prim[1]

    normalized_ratio_3 = (
        pd1 == 3 * pd0
    )

    print()
    print("=" * 78)
    print("8. NORMALIZED DIFFERENCE GEOMETRY")
    print("=" * 78)

    print(
        f"  primitive_d0={pd0}"
    )

    print(
        f"  primitive_d1={pd1}"
    )

    print(
        f"  primitive_d1=3*primitive_d0="
        f"{normalized_ratio_3}"
    )

    print(
        f"  primitive_d0/(-4*q3_primitive)="
        f"{pd0 == -4 * q3_prim}"
    )

    print(
        f"  primitive_d1/(-12*q3_primitive)="
        f"{pd1 == -12 * q3_prim}"
    )

    # ------------------------------------------------------------------
    # 9. EXACT RECONSTRUCTION FROM SOURCE
    # ------------------------------------------------------------------

    reconstructed_s = [
        source_gcd * (
            q1_prim - 2 * q3_prim
        ),
        source_gcd * (
            q1_prim - 6 * q3_prim
        ),
        source_gcd * (
            q1_prim - 18 * q3_prim
        ),
    ]

    reconstruction_exact = (
        reconstructed_s == s
    )

    print()
    print("=" * 78)
    print("9. EXACT SOURCE-TO-SCHUR RECONSTRUCTION")
    print("=" * 78)

    print(
        f"  reconstructed={reconstructed_s}"
    )

    print(
        f"  original={s}"
    )

    print(
        f"  reconstruction_exact="
        f"{reconstruction_exact}"
    )

    # ------------------------------------------------------------------
    # 10. PRIME CONTENT OF NORMALIZED SOURCE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. NORMALIZED SOURCE PRIME SPECTRUM")
    print("=" * 78)

    for label, value in [
        ("q1_primitive", q1_prim),
        ("q3_primitive", q3_prim),
    ]:

        print(
            f"  {label}: "
            f"v2={v2(value)} "
            f"v3={vp(value,3)} "
            f"v5={vp(value,5)} "
            f"v7={vp(value,7)} "
            f"v17={vp(value,17)}"
        )

    # ------------------------------------------------------------------
    # 11. INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("11. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The terminal Schur row has the exact form

    s_k = q1 - 2*3^k*q3,

for k = 0,1,2.

Experiment 179 checks the reverse implication:

    the three Schur coefficients are sufficient to reconstruct q1 and
    q3 exactly.

The first differences satisfy

    s1-s0 = -4 q3,
    s2-s1 = -12 q3,

so

    s2-s1 = 3(s1-s0).

This is an exact geometric difference law coming directly from the
powers 1,3,9 in the terminal p^2 basis.

The common factor 17 is then traced back to the source pair

    (q1,q3).

If

    gcd(q1,q3) = gcd(s0,s1,s2),

then the entire common terminal-row content is inherited from the
source q-values rather than being generated by a new Schur-specific
factor.

No statement about other p-values is inferred.
"""
    )

    # ------------------------------------------------------------------
    # 12. FINAL
    # ------------------------------------------------------------------

    source_gcd_ok = (
        source_gcd == 17
    )

    schur_gcd_ok = (
        schur_gcd == 17
    )

    recovery_ok = (
        q3_from_01 == q3
        and q3_from_12 == q3
        and q1_from_s0s1 == q1
        and q1_from_s1s2 == q1
    )

    geometry_ok = (
        ratio_3_exact
        and normalized_ratio_3
    )

    primitive_ok = (
        primitive_source_gcd == 1
        and primitive_schur_gcd == 1
    )

    final_ok = (
        source_gcd_ok
        and schur_gcd_ok
        and content_matches
        and recovery_ok
        and geometry_ok
        and primitive_ok
        and reconstruction_exact
    )

    print()
    print("=" * 78)
    print("12. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  source_gcd_is_17={source_gcd_ok}"
    )

    print(
        f"  schur_gcd_is_17={schur_gcd_ok}"
    )

    print(
        f"  source_schur_content_match="
        f"{content_matches}"
    )

    print(
        f"  exact_source_recovery="
        f"{recovery_ok}"
    )

    print(
        f"  geometric_difference_law="
        f"{geometry_ok}"
    )

    print(
        f"  normalized_source_primitive="
        f"{primitive_ok}"
    )

    print(
        f"  source_to_schur_reconstruction="
        f"{reconstruction_exact}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 179 COMPLETE")


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

