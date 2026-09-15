#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 174 — EXACT TERMINAL q0 / p-BASIS / 2-ADIC FIRST-LAYER AUDIT
==============================================================================

Experiment 173 established:

    C0 = q_1(0) [1,0,0,1,0,1,0,...]
    C1 = q_3(0) [1,0,0,3,0,9,0,...]

and

    C0-C1 != 0 mod 4,

with coordinate valuations

    1-coordinate : v2 = 3
    p-coordinate : v2 = 1
    p^2-coordinate: v2 = 4.

Thus the first nonzero 2-adic layer is uniquely the p-coordinate.

Experiment 174 isolates the scalar source of those three coordinates.

For terminal rows:

    C(p) = q_p(0) * [1, p, p^2]

in the active coordinates.

Therefore:

    Delta_1  = q_1(0) - q_3(0)
    Delta_p  = q_1(0) - 3 q_3(0)
    Delta_p2 = q_1(0) - 9 q_3(0).

The experiment computes these three scalars exactly and studies:

    * their exact values;
    * their 2-adic valuations;
    * their odd parts;
    * their residues modulo powers of 2;
    * the exact reconstruction of C0-C1;
    * the first 2-adic layer;
    * whether the dominance of the p-coordinate is forced by the
      scalar q-values alone.

No SymPy.
No floating point.
No extrapolation.
No recurrence fitting.
"""

from __future__ import annotations

import sys


# ============================================================================
# EXACT DATA
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

def q0(p):
    return Q[p][
        len(Q[p]) - 1
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


def odd_part(x):
    if x == 0:
        return 0

    e = v2(x)

    return x // (2 ** e)


def mod_profile(x, max_e=8):
    return [
        (2 ** e, x % (2 ** e))
        for e in range(1, max_e + 1)
    ]


def support(v):
    return [
        i
        for i, x in enumerate(v)
        if x != 0
    ]


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 174 — EXACT TERMINAL q0 / p-BASIS / "
        "2-ADIC FIRST-LAYER AUDIT"
    )
    print("=" * 78)

    p0 = 1
    p1 = 3

    q0_1 = q0(p0)
    q0_3 = q0(p1)

    # ------------------------------------------------------------------
    # 1. TERMINAL SCALARS
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. EXACT TERMINAL q0 DATA")
    print("=" * 78)

    print(
        f"  q_1(0)={q0_1}"
    )

    print(
        f"  q_3(0)={q0_3}"
    )

    print(
        f"  v2(q_1(0))={v2(q0_1)}"
    )

    print(
        f"  v2(q_3(0))={v2(q0_3)}"
    )

    print(
        f"  odd_part(q_1(0))={odd_part(q0_1)}"
    )

    print(
        f"  odd_part(q_3(0))={odd_part(q0_3)}"
    )

    # ------------------------------------------------------------------
    # 2. THREE ACTIVE COORDINATE DIFFERENCES
    # ------------------------------------------------------------------

    delta_1 = (
        q0_1
        - q0_3
    )

    delta_p = (
        q0_1
        - 3 * q0_3
    )

    delta_p2 = (
        q0_1
        - 9 * q0_3
    )

    labels = [
        "1",
        "p",
        "p^2",
    ]

    values = [
        delta_1,
        delta_p,
        delta_p2,
    ]

    print()
    print("=" * 78)
    print("2. EXACT ACTIVE-COORDINATE DIFFERENCES")
    print("=" * 78)

    for label, value in zip(
        labels,
        values,
    ):

        print(
            f"  coordinate={label}: "
            f"value={value} "
            f"v2={v2(value)} "
            f"odd_part={odd_part(value)}"
        )

    # ------------------------------------------------------------------
    # 3. EXACT TERMINAL ROW RECONSTRUCTION
    # ------------------------------------------------------------------

    C0 = [
        q0_1,
        0,
        0,
        q0_1,
        0,
        q0_1,
        0,
        0,
        0,
        0,
        0,
        0,
    ]

    C1 = [
        q0_3,
        0,
        0,
        3 * q0_3,
        0,
        9 * q0_3,
        0,
        0,
        0,
        0,
        0,
        0,
    ]

    delta = [
        C0[j] - C1[j]
        for j in range(12)
    ]

    reconstructed = [
        delta_1,
        0,
        0,
        delta_p,
        0,
        delta_p2,
        0,
        0,
        0,
        0,
        0,
        0,
    ]

    reconstruction_exact = (
        reconstructed == delta
    )

    print()
    print("=" * 78)
    print("3. EXACT FULL 12-DIMENSIONAL RECONSTRUCTION")
    print("=" * 78)

    print(
        f"  DeltaC={delta}"
    )

    print(
        f"  reconstructed={reconstructed}"
    )

    print(
        f"  reconstruction_exact="
        f"{reconstruction_exact}"
    )

    print(
        f"  DeltaC_support="
        f"{support(delta)}"
    )

    # ------------------------------------------------------------------
    # 4. 2-ADIC ORDERING
    # ------------------------------------------------------------------

    coordinate_v2 = [
        v2(delta_1),
        v2(delta_p),
        v2(delta_p2),
    ]

    minimum_v2 = min(
        x
        for x in coordinate_v2
        if x is not None
    )

    minimizers = [
        labels[i]
        for i, e in enumerate(coordinate_v2)
        if e == minimum_v2
    ]

    print()
    print("=" * 78)
    print("4. EXACT 2-ADIC ORDERING")
    print("=" * 78)

    print(
        f"  coordinate_v2={coordinate_v2}"
    )

    print(
        f"  minimum_v2={minimum_v2}"
    )

    print(
        f"  minimum_coordinates={minimizers}"
    )

    # ------------------------------------------------------------------
    # 5. FIRST NONZERO LAYER
    # ------------------------------------------------------------------

    first_layer = [
        (
            value // (2 ** minimum_v2)
        ) % 2
        for value in values
    ]

    print()
    print("=" * 78)
    print("5. FIRST NONZERO 2-ADIC LAYER")
    print("=" * 78)

    print(
        f"  Delta_active / 2^{minimum_v2} mod 2="
        f"{first_layer}"
    )

    print(
        f"  first_layer_support="
        f"{support(first_layer)}"
    )

    # ------------------------------------------------------------------
    # 6. RESIDUE PROFILES
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. SCALAR 2-ADIC RESIDUE PROFILES")
    print("=" * 78)

    for label, value in zip(
        labels,
        values,
    ):

        print(
            f"  {label}: "
            f"{mod_profile(value, 8)}"
        )

    # ------------------------------------------------------------------
    # 7. DOMINANCE TEST
    # ------------------------------------------------------------------

    p_coordinate_unique_minimum = (
        v2(delta_p) == minimum_v2
        and v2(delta_1) > minimum_v2
        and v2(delta_p2) > minimum_v2
    )

    print()
    print("=" * 78)
    print("7. p-COORDINATE DOMINANCE TEST")
    print("=" * 78)

    print(
        f"  p_coordinate_unique_minimum="
        f"{p_coordinate_unique_minimum}"
    )

    # ------------------------------------------------------------------
    # 8. MOD-4 EXPLANATION
    # ------------------------------------------------------------------

    mod2_equal = all(
        x % 2 == 0
        for x in delta
    )

    mod4_different = any(
        x % 4 != 0
        for x in delta
    )

    print()
    print("=" * 78)
    print("8. MOD-2 / MOD-4 CONSEQUENCE")
    print("=" * 78)

    print(
        f"  C0_equals_C1_mod2={mod2_equal}"
    )

    print(
        f"  C0_differs_mod4={mod4_different}"
    )

    print(
        f"  first_layer_explains_mod4_split="
        f"{p_coordinate_unique_minimum}"
    )

    # ------------------------------------------------------------------
    # 9. SECOND LAYER
    # ------------------------------------------------------------------

    stripped_mod4 = [
        (
            value // (2 ** minimum_v2)
        ) % 4
        for value in values
    ]

    print()
    print("=" * 78)
    print("9. FIRST-LAYER MOD-4 REFINEMENT")
    print("=" * 78)

    print(
        f"  Delta_active / 2^{minimum_v2} mod 4="
        f"{stripped_mod4}"
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
For terminal rows d=0, the active basis vector is

    [1, p, p^2].

Thus the exact difference between the p=1 and p=3 terminal rows is
controlled by three scalar combinations:

    q_1(0) - q_3(0),
    q_1(0) - 3 q_3(0),
    q_1(0) - 9 q_3(0).

The previous experiment found valuations

    3, 1, 4.

Experiment 174 verifies those valuations directly from the scalar
terminal q-values.

Therefore the first nonzero 2-adic distinction is caused specifically
by the p-coordinate combination

    q_1(0) - 3 q_3(0),

rather than by a cancellation among several basis coordinates.

This is still only a finite exact observation. It does not establish
a general rule for other p-values.
"""
    )

    # ------------------------------------------------------------------
    # 11. FINAL
    # ------------------------------------------------------------------

    data_exact = (
        q0_1 == 495451247
        and q0_3 == 421514439
    )

    final_ok = (
        data_exact
        and reconstruction_exact
        and mod2_equal
        and mod4_different
        and minimum_v2 == 1
        and p_coordinate_unique_minimum
        and first_layer == [0, 1, 0]
    )

    print()
    print("=" * 78)
    print("11. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  terminal_q_data_exact={data_exact}"
    )

    print(
        f"  reconstruction_exact="
        f"{reconstruction_exact}"
    )

    print(
        f"  C0_equals_C1_mod2="
        f"{mod2_equal}"
    )

    print(
        f"  C0_differs_mod4="
        f"{mod4_different}"
    )

    print(
        f"  minimum_v2_equals_1="
        f"{minimum_v2 == 1}"
    )

    print(
        f"  p_coordinate_unique_minimum="
        f"{p_coordinate_unique_minimum}"
    )

    print(
        f"  first_layer_is_0_1_0="
        f"{first_layer == [0, 1, 0]}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 174 COMPLETE")


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
