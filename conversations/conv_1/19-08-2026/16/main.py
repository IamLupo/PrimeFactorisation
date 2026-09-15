#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 185 — EXACT 3-ADIC TERMINAL SOURCE / SCHUR-TRANSFER AUDIT
==============================================================================

Experiment 184 established an exact terminal counterfactual:

    C_original = 17 C_modified
    K_original = 17 K_modified,

and after removing 17:

    gcd(C_modified) = 1
    gcd(K_modified) = 1.

Thus the factor 17 is completely localized to the terminal source pair.

The primitive terminal sources are

    q1' = q1 / 17 = 29144191
    q3' = q3 / 17 = 24794967.

Experiment 185 studies the next visible odd prime, 3.

The exact questions are:

    * What are v3(q1') and v3(q3')?
    * What are the v3-valuations of the primitive Schur coefficients?
    * Is any common 3-content created by the terminal Schur map?
    * Is the 3-adic profile inherited entrywise from the source pair?
    * Does the primitive Schur row retain the asymmetric 3-adic structure?
    * What happens to the 3-adic profile under the same terminal
      source normalization used in Experiment 184?

No determinant computation is required.
No floating point.
No SymPy.
No extrapolation.
No general theorem is inferred from this single terminal pair.
"""

from __future__ import annotations

import sys
from math import gcd


# ============================================================================
# DATA
# ============================================================================

Q1_ORIGINAL = 495451247
Q3_ORIGINAL = 421514439

CONTENT_17 = 17

Q1 = Q1_ORIGINAL // CONTENT_17
Q3 = Q3_ORIGINAL // CONTENT_17

# Primitive terminal Schur row from the established exact construction.
S0 = Q1 - 2 * Q3
S1 = Q1 - 6 * Q3
S2 = Q1 - 18 * Q3

PRIMITIVE_S = [S0, S1, S2]


# ============================================================================
# HELPERS
# ============================================================================

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


def odd_part_at_prime(x, p):
    if x == 0:
        return 0

    e = vp(x, p)

    return x // (p ** e)


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 185 — EXACT 3-ADIC TERMINAL SOURCE / "
        "SCHUR-TRANSFER AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. PRIMITIVE SOURCE PAIR
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. PRIMITIVE TERMINAL SOURCE")
    print("=" * 78)

    print(
        f"  q1={Q1}"
    )

    print(
        f"  q3={Q3}"
    )

    print(
        f"  gcd(q1,q3)={gcd(Q1, Q3)}"
    )

    print(
        f"  v3(q1)={vp(Q1, 3)}"
    )

    print(
        f"  v3(q3)={vp(Q3, 3)}"
    )

    print(
        f"  q1_mod3={Q1 % 3}"
    )

    print(
        f"  q3_mod3={Q3 % 3}"
    )

    # ------------------------------------------------------------------
    # 2. PRIMITIVE SCHUR COEFFICIENTS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. PRIMITIVE TERMINAL SCHUR COEFFICIENTS")
    print("=" * 78)

    for i, value in enumerate(PRIMITIVE_S):

        print(
            f"  s{i}={value} "
            f"v3={vp(value,3)} "
            f"mod3={value % 3}"
        )

    primitive_s_gcd = gcd_list(
        PRIMITIVE_S
    )

    print(
        f"  gcd(s0,s1,s2)={primitive_s_gcd}"
    )

    print(
        f"  v3(gcd(s))="
        f"{vp(primitive_s_gcd,3)}"
    )

    # ------------------------------------------------------------------
    # 3. EXACT 3-ADIC SOURCE/OUTPUT CONTENT
    # ------------------------------------------------------------------

    source_v3 = min(
        vp(Q1, 3),
        vp(Q3, 3)
    )

    schur_v3 = min(
        vp(x, 3)
        for x in PRIMITIVE_S
        if x != 0
    )

    print()
    print("=" * 78)
    print("3. 3-ADIC CONTENT TRANSFER")
    print("=" * 78)

    print(
        f"  source_min_v3={source_v3}"
    )

    print(
        f"  schur_min_v3={schur_v3}"
    )

    print(
        f"  source_common_3_content="
        f"{3 ** source_v3}"
    )

    print(
        f"  schur_common_3_content="
        f"{3 ** schur_v3}"
    )

    print(
        f"  content_preserved="
        f"{source_v3 == schur_v3}"
    )

    # ------------------------------------------------------------------
    # 4. SOURCE-TO-SCHUR FORMULAS
    # ------------------------------------------------------------------

    reconstructed = [
        Q1 - 2 * Q3,
        Q1 - 6 * Q3,
        Q1 - 18 * Q3,
    ]

    reconstruction_exact = (
        reconstructed == PRIMITIVE_S
    )

    print()
    print("=" * 78)
    print("4. EXACT SOURCE-TO-SCHUR RECONSTRUCTION")
    print("=" * 78)

    print(
        f"  reconstructed={reconstructed}"
    )

    print(
        f"  original={PRIMITIVE_S}"
    )

    print(
        f"  reconstruction_exact="
        f"{reconstruction_exact}"
    )

    # ------------------------------------------------------------------
    # 5. 3-ADIC DIFFERENCE GEOMETRY
    # ------------------------------------------------------------------

    d0 = S1 - S0
    d1 = S2 - S1

    print()
    print("=" * 78)
    print("5. 3-ADIC DIFFERENCE GEOMETRY")
    print("=" * 78)

    print(
        f"  d0=s1-s0={d0} "
        f"v3={vp(d0,3)}"
    )

    print(
        f"  d1=s2-s1={d1} "
        f"v3={vp(d1,3)}"
    )

    print(
        f"  d1=3*d0="
        f"{d1 == 3 * d0}"
    )

    print(
        f"  d0/(-4*q3)="
        f"{d0 == -4 * Q3}"
    )

    print(
        f"  d1/(-12*q3)="
        f"{d1 == -12 * Q3}"
    )

    # ------------------------------------------------------------------
    # 6. 3-ADIC RESIDUE LAYERS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. 3-ADIC RESIDUE PROFILE")
    print("=" * 78)

    for i, value in enumerate(PRIMITIVE_S):

        residues = [
            (3 ** e, value % (3 ** e))
            for e in range(1, 7)
        ]

        print(
            f"  s{i}: {residues}"
        )

    # ------------------------------------------------------------------
    # 7. REMOVE THE SOURCE 3-FACTOR FROM q3
    # ------------------------------------------------------------------

    q3_v3 = vp(Q3, 3)

    q3_unit = (
        Q3 // (3 ** q3_v3)
    )

    q1_unit = Q1

    unit_s = [
        q1_unit - 2 * q3_unit,
        q1_unit - 6 * q3_unit,
        q1_unit - 18 * q3_unit,
    ]

    print()
    print("=" * 78)
    print("7. q3 3-FACTOR REMOVAL")
    print("=" * 78)

    print(
        f"  q3_v3={q3_v3}"
    )

    print(
        f"  q3_unit={q3_unit}"
    )

    print(
        f"  q1_unit={q1_unit}"
    )

    for i, value in enumerate(unit_s):

        print(
            f"  unit_s{i}={value} "
            f"v3={vp(value,3)}"
        )

    # ------------------------------------------------------------------
    # 8. COMPARISON OF 3-ADIC STRUCTURE
    # ------------------------------------------------------------------

    source_asymmetric = (
        vp(Q1, 3) != vp(Q3, 3)
    )

    schur_asymmetric = len({
        vp(x, 3)
        for x in PRIMITIVE_S
    }) > 1

    print()
    print("=" * 78)
    print("8. 3-ADIC ASYMMETRY TEST")
    print("=" * 78)

    print(
        f"  source_v3_pair="
        f"[{vp(Q1,3)}, {vp(Q3,3)}]"
    )

    print(
        f"  source_asymmetric="
        f"{source_asymmetric}"
    )

    print(
        f"  schur_v3_triple="
        f"{[vp(x,3) for x in PRIMITIVE_S]}"
    )

    print(
        f"  schur_asymmetric="
        f"{schur_asymmetric}"
    )

    # ------------------------------------------------------------------
    # 9. COMMON 3 CONTENT
    # ------------------------------------------------------------------

    source_gcd = gcd(
        abs(Q1),
        abs(Q3)
    )

    schur_gcd = gcd_list(
        PRIMITIVE_S
    )

    print()
    print("=" * 78)
    print("9. COMMON 3-CONTENT")
    print("=" * 78)

    print(
        f"  gcd_source={source_gcd}"
    )

    print(
        f"  gcd_schur={schur_gcd}"
    )

    print(
        f"  source_has_common_3="
        f"{source_gcd % 3 == 0}"
    )

    print(
        f"  schur_has_common_3="
        f"{schur_gcd % 3 == 0}"
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
After removing the localized factor 17, the terminal source pair is

    q1 = 29144191
    q3 = 24794967.

The two source entries have different 3-adic orders:

    v3(q1) = 0,
    v3(q3) = 1.

Experiment 185 checks whether that asymmetry survives the exact terminal
Schur map

    s_k = q1 - 2*3^k*q3.

Because q3 carries one factor of 3 while q1 does not, the three Schur
coefficients can have different 3-adic behavior even though their
common gcd is 1.

A positive transfer result would show that the 3-adic profile is
directly inherited from the source pair.

A negative result would indicate cancellation between q1 and the
3^k-weighted q3 terms.

No claim is made beyond this exact three-coefficient calculation.
"""
    )

    # ------------------------------------------------------------------
    # 11. FINAL
    # ------------------------------------------------------------------

    final_ok = (
        reconstruction_exact
        and vp(Q1,3) == 0
        and vp(Q3,3) == 1
        and gcd(Q1,Q3) == 1
        and gcd_list(PRIMITIVE_S) == 1
        and d1 == 3 * d0
    )

    print()
    print("=" * 78)
    print("11. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  source_reconstruction_exact="
        f"{reconstruction_exact}"
    )

    print(
        f"  v3(q1)_is_0="
        f"{vp(Q1,3) == 0}"
    )

    print(
        f"  v3(q3)_is_1="
        f"{vp(Q3,3) == 1}"
    )

    print(
        f"  source_pair_primitive="
        f"{gcd(Q1,Q3) == 1}"
    )

    print(
        f"  Schur_triplet_primitive="
        f"{gcd_list(PRIMITIVE_S) == 1}"
    )

    print(
        f"  exact_geometric_difference="
        f"{d1 == 3*d0}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 185 COMPLETE")


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
