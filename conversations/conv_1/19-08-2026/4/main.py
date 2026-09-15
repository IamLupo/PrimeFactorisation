#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 173 — EXACT TERMINAL-BOUNDARY DIFFERENCE /
                BASIS-COORDINATE PROVENANCE AUDIT
==============================================================================

Experiment 172 established:

    C0 == C1 (mod 2)

but

    C0 != C1 (mod 4),

with

    v2(C0-C1) = 1.

Moreover,

    (C0-C1)/2 mod 2

has support only at coordinate 3 of the six-dimensional operator basis.

This experiment traces that first difference back to the exact basis
coordinates.

The terminal equations have d=0 and therefore each C row has the form

    q_p(0) * [1, 0, 0, p, 0, p^2]

inside the twelve F coordinates, because the q(r+1) term vanishes.

Thus the difference between the two terminal rows can be decomposed
coordinate-by-coordinate.

The six basis coordinates are

    0 : 1
    1 : d
    2 : d^2
    3 : p
    4 : p*d
    5 : p^2

The experiment asks:

    * Which basis coordinates actually contribute to C0-C1?
    * What are their exact 2-adic valuations?
    * Which coordinate gives the minimum valuation?
    * Does the first nonzero mod-2 layer really come from one basis
      coordinate?
    * Is that behavior intrinsic to the exact terminal rows?

No SymPy.
No floating point.
No extrapolation.
No recurrence fitting.
"""

from __future__ import annotations

import sys


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

D = {
    1: 5,
    3: 4,
    5: 2,
    7: 0,
}

TRANSITIONS = [
    (1, 3),
    (3, 5),
    (5, 7),
]

FULL_N = 14
CORE_N = 12

BOUNDARY_ROWS = (5, 10)
BOUNDARY_COLS = (12, 13)


# ============================================================================
# BASIC HELPERS
# ============================================================================

def q_value(p, r):
    if r < 0 or r >= len(Q[p]):
        return 0
    return Q[p][r]


def basis(p, d):
    return [
        1,
        d,
        d * d,
        p,
        p * d,
        p * p,
    ]


def v2(x):
    if x == 0:
        return None

    x = abs(x)
    e = 0

    while x % 2 == 0:
        x //= 2
        e += 1

    return e


def support(v):
    return [
        i
        for i, x in enumerate(v)
        if x != 0
    ]


# ============================================================================
# FULL SYSTEM
# ============================================================================

def build_system():

    M = []
    metadata = []

    for p, p_next in TRANSITIONS:

        for r in range(D[p] + 1):

            d = D[p] - r

            q0 = q_value(p, r)
            q1 = q_value(p, r + 1)

            b = basis(p, d)

            row = [0] * FULL_N

            for j in range(6):
                row[j] = b[j] * q0

            for j in range(6):
                row[6 + j] = b[j] * q1

            if d == 0:
                row[12] = 1
                row[13] = p

            M.append(row)

            metadata.append(
                {
                    "full_row": len(M) - 1,
                    "p": p,
                    "p_next": p_next,
                    "r": r,
                    "d": d,
                    "q0": q0,
                    "q1": q1,
                    "basis": b,
                }
            )

    return M, metadata


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 173 — EXACT TERMINAL-BOUNDARY DIFFERENCE / "
        "BASIS-COORDINATE PROVENANCE AUDIT"
    )
    print("=" * 78)

    M, metadata = build_system()

    m0 = metadata[BOUNDARY_ROWS[0]]
    m1 = metadata[BOUNDARY_ROWS[1]]

    # ------------------------------------------------------------------
    # 1. TERMINAL DATA
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. EXACT TERMINAL DATA")
    print("=" * 78)

    print(
        f"  row0: full_row={m0['full_row']} "
        f"p={m0['p']} r={m0['r']} d={m0['d']} "
        f"q0={m0['q0']}"
    )

    print(
        f"  row1: full_row={m1['full_row']} "
        f"p={m1['p']} r={m1['r']} d={m1['d']} "
        f"q0={m1['q0']}"
    )

    print(
        f"  basis0={m0['basis']}"
    )

    print(
        f"  basis1={m1['basis']}"
    )

    # ------------------------------------------------------------------
    # 2. ACTUAL C ROWS
    # ------------------------------------------------------------------

    C0 = M[BOUNDARY_ROWS[0]][:CORE_N]
    C1 = M[BOUNDARY_ROWS[1]][:CORE_N]

    delta = [
        C0[j] - C1[j]
        for j in range(CORE_N)
    ]

    print()
    print("=" * 78)
    print("2. EXACT C-ROW DIFFERENCE")
    print("=" * 78)

    print(
        f"  C0={C0}"
    )

    print(
        f"  C1={C1}"
    )

    print(
        f"  DeltaC={delta}"
    )

    print(
        f"  support={support(delta)}"
    )

    # ------------------------------------------------------------------
    # 3. SIX BASIS COORDINATES
    # ------------------------------------------------------------------

    labels = [
        "1",
        "d",
        "d^2",
        "p",
        "p*d",
        "p^2",
    ]

    coordinate_differences = []

    for k, label in enumerate(labels):

        x0 = m0["q0"] * m0["basis"][k]
        x1 = m1["q0"] * m1["basis"][k]

        diff = x0 - x1

        coordinate_differences.append(diff)

        print(
            f"  coordinate={label}: "
            f"row0={x0} "
            f"row1={x1} "
            f"difference={diff} "
            f"v2={v2(diff)}"
        )

    # ------------------------------------------------------------------
    # 4. VERIFY RECONSTRUCTION
    # ------------------------------------------------------------------

    reconstructed = [0] * CORE_N

    for k in range(6):
        reconstructed[k] = coordinate_differences[k]

    reconstruction_exact = (
        reconstructed == delta[:6]
    )

    print()
    print("=" * 78)
    print("4. BASIS-COORDINATE RECONSTRUCTION")
    print("=" * 78)

    print(
        f"  reconstructed_first_six="
        f"{reconstructed}"
    )

    print(
        f"  actual_first_six="
        f"{delta[:6]}"
    )

    print(
        f"  reconstruction_exact="
        f"{reconstruction_exact}"
    )

    # ------------------------------------------------------------------
    # 5. COORDINATE VALUATION PROFILE
    # ------------------------------------------------------------------

    coordinate_v2 = [
        v2(x)
        for x in coordinate_differences
    ]

    finite_v2 = [
        e
        for e in coordinate_v2
        if e is not None
    ]

    min_v2 = min(finite_v2)

    minimizers = [
        k
        for k, e in enumerate(coordinate_v2)
        if e == min_v2
    ]

    print()
    print("=" * 78)
    print("5. COORDINATE 2-ADIC PROFILE")
    print("=" * 78)

    print(
        f"  coordinate_v2={coordinate_v2}"
    )

    print(
        f"  minimum_v2={min_v2}"
    )

    print(
        f"  minimum_coordinates={minimizers}"
    )

    # ------------------------------------------------------------------
    # 6. FIRST CORRECTION LAYER
    # ------------------------------------------------------------------

    first_layer = [
        0 if x % (2 ** min_v2) != 0
        else (x // (2 ** min_v2)) % 2
        for x in coordinate_differences
    ]

    print()
    print("=" * 78)
    print("6. FIRST NONZERO 2-ADIC LAYER")
    print("=" * 78)

    print(
        f"  DeltaC / 2^{min_v2} mod 2="
        f"{first_layer}"
    )

    print(
        f"  first_layer_support="
        f"{support(first_layer)}"
    )

    # ------------------------------------------------------------------
    # 7. MOD-2 / MOD-4 / MOD-8 COMPARISON
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. TERMINAL ROW CONGRUENCE")
    print("=" * 78)

    for e in range(1, 8):

        modulus = 2 ** e

        equal = all(
            (C0[j] - C1[j]) % modulus == 0
            for j in range(CORE_N)
        )

        print(
            f"  modulus={modulus} "
            f"C0_equals_C1={equal}"
        )

    # ------------------------------------------------------------------
    # 8. BASIS EXPLANATION TEST
    # ------------------------------------------------------------------

    only_first_layer_coordinate = (
        len(minimizers) == 1
    )

    print()
    print("=" * 78)
    print("8. BASIS-LEVEL FIRST-LAYER TEST")
    print("=" * 78)

    print(
        f"  unique_minimum_coordinate="
        f"{only_first_layer_coordinate}"
    )

    if only_first_layer_coordinate:

        print(
            f"  selected_coordinate="
            f"{labels[minimizers[0]]}"
        )

    # ------------------------------------------------------------------
    # 9. SECOND LAYER
    # ------------------------------------------------------------------

    second_layer = []

    for x in coordinate_differences:

        if x == 0:
            second_layer.append(None)
            continue

        if x % (2 ** min_v2) != 0:
            second_layer.append(None)
            continue

        y = x // (2 ** min_v2)

        second_layer.append(
            y % 4
        )

    print()
    print("=" * 78)
    print("9. FIRST-LAYER MOD-4 REFINEMENT")
    print("=" * 78)

    print(
        f"  DeltaC / 2^{min_v2} mod 4="
        f"{second_layer}"
    )

    # ------------------------------------------------------------------
    # 10. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
Experiment 172 showed:

    C0 = C1 (mod 2),

but

    C0 != C1 (mod 4).

The present experiment asks where that first distinction lives in the
six-dimensional operator basis.

Because the terminal rows have d=0, their basis evaluations simplify
to

    [1, 0, 0, p, 0, p^2].

Thus only

    1,
    p,
    p^2

can contribute directly to the terminal C row.

The coordinate valuation profile determines which of those basis
coordinates actually produces the first nonzero 2-adic difference.

A unique minimum would localize the first boundary-row distinction to
one basis coordinate.

A non-unique minimum would show that several coordinates contribute at
the same first 2-adic order.

The purpose is descriptive: isolate the arithmetic source of the
mod-4 separation rather than infer a general law.
"""
    )

    # ------------------------------------------------------------------
    # 11. FINAL
    # ------------------------------------------------------------------

    data_exact = (
        len(C0) == CORE_N
        and len(C1) == CORE_N
    )

    mod2_equal = all(
        x % 2 == 0
        for x in delta
    )

    mod4_different = any(
        x % 4 != 0
        for x in delta
    )

    final_ok = (
        data_exact
        and reconstruction_exact
        and mod2_equal
        and mod4_different
        and min_v2 == 1
    )

    print()
    print("=" * 78)
    print("11. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  terminal_rows_exact={data_exact}"
    )

    print(
        f"  basis_reconstruction_exact="
        f"{reconstruction_exact}"
    )

    print(
        f"  C0_equals_C1_mod2={mod2_equal}"
    )

    print(
        f"  C0_differs_mod4={mod4_different}"
    )

    print(
        f"  minimum_v2_equals_1="
        f"{min_v2 == 1}"
    )

    print(
        f"  unique_first_layer_coordinate="
        f"{only_first_layer_coordinate}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 173 COMPLETE")


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

