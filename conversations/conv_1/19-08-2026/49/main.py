#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 216RRR — EXACT PROJECTIVE 7-ADIC COORDINATE /
                    COMMON-CONTENT DECOUPLING AUDIT
==============================================================================

Corrected Experiment 216RR.

Why the previous control failed
--------------------------------

For the original source pair,

    rho = q1/q3 = 6 mod 7,

and the principal-coordinate lift has slope 1, so

    t = -delta mod 7.

After one-sided scaling,

    (q1,q3) -> (17*q1,q3),

the ratio changes. The lift slope is no longer guaranteed to be 1.

Therefore the general reconstruction must test all seven candidates

    k_candidate = k_current + 6*t*7^(e-1),

    t = 0,...,6,

and select the unique candidate satisfying the next congruence.

This script does exactly that.

The experiment tests:

    1. common source scaling preserves q1/q3;
    2. common source scaling preserves the complete 7-adic exponent data;
    3. the terminal Schur source row scales linearly;
    4. one-sided scaling changes the projective 7-adic coordinate;
    5. the general seven-candidate lift remains exact even when slope != 1.

No hard-coded exponent sequence.
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

P = 7

END_E = 16

COMMON_SCALES = [
    2,
    3,
    5,
    17,
    17 ** 2,
]


# ============================================================================
# GENERIC EXACT HELPERS
# ============================================================================

def valuation_p(
    x: int,
    p: int,
) -> int | None:

    x = int(x)
    p = int(p)

    if p < 2:
        raise ValueError(
            "p must be >= 2."
        )

    if x == 0:
        return None

    x = abs(x)
    v = 0

    while x % p == 0:
        x //= p
        v += 1

    return v


def valuation_7(x: int) -> int | None:
    return valuation_p(
        x,
        7,
    )


def valuation_17(x: int) -> int | None:
    return valuation_p(
        x,
        17,
    )


def gcd_int(
    a: int,
    b: int,
) -> int:

    a = abs(int(a))
    b = abs(int(b))

    while b:
        a, b = (
            b,
            a % b,
        )

    return a


def inverse_mod(
    a: int,
    m: int,
) -> int:

    a = int(a)
    m = int(m)

    if m <= 1:
        return 0

    a %= m

    if a == 0:
        raise ArithmeticError(
            f"Cannot invert 0 modulo {m}."
        )

    old_r = a
    r = m

    old_s = 1
    s = 0

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

    return int(
        old_s % m
    )


def ratio_mod(
    q1: int,
    q3: int,
    modulus: int,
) -> int:

    q1 = int(q1)
    q3 = int(q3)
    modulus = int(modulus)

    return int(
        (
            q1
            * inverse_mod(
                q3,
                modulus,
            )
        )
        % modulus
    )


def full_orbit_mod(
    k: int,
    modulus: int,
) -> int:

    k = int(k)
    modulus = int(modulus)

    return int(
        (
            2
            * pow(
                3,
                k,
                modulus,
            )
        )
        % modulus
    )


# ============================================================================
# GENERAL EXPONENT RECONSTRUCTION
# ============================================================================

def find_base_k1(
    q1: int,
    q3: int,
) -> int:

    rho = ratio_mod(
        q1,
        q3,
        7,
    )

    matches = []

    for k in range(6):

        value = int(
            (
                2
                * pow(
                    3,
                    k,
                    7,
                )
            )
            % 7
        )

        if value == rho:
            matches.append(
                int(k)
            )

    if len(matches) != 1:
        raise ArithmeticError(
            f"Expected one mod-7 exponent class; got {matches}."
        )

    return int(
        matches[0]
    )


def candidate_lifts(
    current_k: int,
    e: int,
) -> list[int]:

    """
    At level e,

        ord_{7^e}(3) = 6*7^(e-1).

    The seven candidate lifts are

        k + t*ord_{7^e}(3).
    """

    current_k = int(current_k)
    e = int(e)

    order = int(
        6
        * (
            P ** (
                e - 1
            )
        )
    )

    return [
        int(
            current_k
            + t * order
        )
        for t in range(7)
    ]


def reconstruct_general(
    q1: int,
    q3: int,
    end_e: int,
):

    q1 = int(q1)
    q3 = int(q3)

    k1 = find_base_k1(
        q1,
        q3,
    )

    current_k = int(
        k1
    )

    k_sequence = [
        current_k
    ]

    t_digits = []

    delta_sequence = []

    slope_sequence = []

    transition_rows = []

    for e in range(
        1,
        end_e,
    ):

        modulus_current = int(
            P ** e
        )

        modulus_next = int(
            P ** (
                e + 1
            )
        )

        target = ratio_mod(
            q1,
            q3,
            modulus_next,
        )

        current_orbit = full_orbit_mod(
            current_k,
            modulus_next,
        )

        current_residual = int(
            (
                current_orbit
                - target
            )
            % modulus_next
        )

        if current_residual % modulus_current != 0:
            raise ArithmeticError(
                f"Current k={current_k} "
                f"does not solve modulo 7^{e}."
            )

        delta = int(
            (
                current_residual
                // modulus_current
            )
            % P
        )

        matches = []

        candidate_values = []

        for t in range(7):

            candidate_k = int(
                candidate_lifts(
                    current_k,
                    e,
                )[t]
            )

            candidate_orbit = full_orbit_mod(
                candidate_k,
                modulus_next,
            )

            candidate_match = (
                candidate_orbit
                == target
            )

            candidate_values.append(
                {
                    "t": int(t),
                    "k": candidate_k,
                    "orbit": candidate_orbit,
                    "match": candidate_match,
                }
            )

            if candidate_match:
                matches.append(
                    int(t)
                )

        if len(matches) != 1:
            raise ArithmeticError(
                f"At e={e}, expected exactly one lift "
                f"but found matches={matches}."
            )

        chosen_t = int(
            matches[0]
        )

        next_k = int(
            candidate_values[
                chosen_t
            ]["k"]
        )

        # Local residual slope diagnostic:
        #
        # R(t) = R(0) + B*t mod 7.
        #
        r0 = int(
            (
                candidate_values[0]["orbit"]
                - target
            )
            // modulus_current
        ) % P

        r1 = int(
            (
                candidate_values[1]["orbit"]
                - target
            )
            // modulus_current
        ) % P

        slope = int(
            (
                r1
                - r0
            )
            % P
        )

        delta_sequence.append(
            delta
        )

        t_digits.append(
            chosen_t
        )

        slope_sequence.append(
            slope
        )

        transition_rows.append(
            {
                "e": int(e),
                "current_k": int(current_k),
                "delta": delta,
                "t": chosen_t,
                "next_k": next_k,
                "slope": slope,
                "matches": matches,
            }
        )

        k_sequence.append(
            next_k
        )

        current_k = int(
            next_k
        )

    return {
        "k1": int(k1),
        "k": k_sequence,
        "t": t_digits,
        "delta": delta_sequence,
        "slope": slope_sequence,
        "rows": transition_rows,
    }


# ============================================================================
# TERMINAL SCHUR SOURCE ROW
# ============================================================================

def terminal_schur_triplet(
    q1: int,
    q3: int,
):

    return [
        int(
            q1
            - 2 * q3
        ),
        int(
            q1
            - 6 * q3
        ),
        int(
            q1
            - 18 * q3
        ),
    ]


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 216RRR — EXACT PROJECTIVE 7-ADIC "
        "COORDINATE / COMMON-CONTENT DECOUPLING AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. ORIGINAL SOURCE
    # ------------------------------------------------------------------

    q1_v7 = valuation_7(Q1)
    q3_v7 = valuation_7(Q3)

    q1_v17 = valuation_17(Q1)
    q3_v17 = valuation_17(Q3)

    original_gcd = gcd_int(
        Q1,
        Q3,
    )

    print()
    print("=" * 78)
    print("1. ORIGINAL SOURCE")
    print("=" * 78)

    print(
        f"  q1={Q1}"
    )

    print(
        f"  q3={Q3}"
    )

    print(
        f"  gcd(q1,q3)={original_gcd}"
    )

    print(
        f"  v7(q1)={q1_v7}"
    )

    print(
        f"  v7(q3)={q3_v7}"
    )

    print(
        f"  v17(q1)={q1_v17}"
    )

    print(
        f"  v17(q3)={q3_v17}"
    )

    # ------------------------------------------------------------------
    # 2. ORIGINAL GENERAL 7-ADIC LIFT
    # ------------------------------------------------------------------

    original = reconstruct_general(
        Q1,
        Q3,
        END_E,
    )

    print()
    print("=" * 78)
    print("2. ORIGINAL GENERAL 7-ADIC LIFT")
    print("=" * 78)

    print(
        f"  k1={original['k1']}"
    )

    print(
        f"  k_final={original['k'][-1]}"
    )

    print(
        f"  t_digits={original['t']}"
    )

    print(
        f"  delta_digits={original['delta']}"
    )

    print(
        f"  slopes={original['slope']}"
    )

    print(
        f"  all_original_lifts_unique="
        f"{all(len(row['matches']) == 1 for row in original['rows'])}"
    )

    # ------------------------------------------------------------------
    # 3. COMMON SOURCE SCALES
    # ------------------------------------------------------------------

    scale_rows = []

    print()
    print("=" * 78)
    print("3. COMMON SOURCE-SCALE TESTS")
    print("=" * 78)

    for scale in COMMON_SCALES:

        scale = int(scale)

        q1_scaled = int(
            scale * Q1
        )

        q3_scaled = int(
            scale * Q3
        )

        scaled = reconstruct_general(
            q1_scaled,
            q3_scaled,
            END_E,
        )

        ratio_invariant = True

        for e in range(
            1,
            END_E + 1,
        ):

            modulus = int(
                P ** e
            )

            original_ratio = ratio_mod(
                Q1,
                Q3,
                modulus,
            )

            scaled_ratio = ratio_mod(
                q1_scaled,
                q3_scaled,
                modulus,
            )

            if original_ratio != scaled_ratio:
                ratio_invariant = False
                break

        k1_invariant = (
            scaled["k1"]
            == original["k1"]
        )

        k_invariant = (
            scaled["k"]
            == original["k"]
        )

        t_invariant = (
            scaled["t"]
            == original["t"]
        )

        delta_invariant = (
            scaled["delta"]
            == original["delta"]
        )

        scale_rows.append(
            {
                "scale": scale,
                "ratio": ratio_invariant,
                "k1": k1_invariant,
                "k": k_invariant,
                "t": t_invariant,
                "delta": delta_invariant,
            }
        )

        print(
            f"  scale={scale}: "
            f"gcd={gcd_int(q1_scaled,q3_scaled)} "
            f"v17(q1)={valuation_17(q1_scaled)} "
            f"v17(q3)={valuation_17(q3_scaled)} "
            f"k1={scaled['k1']} "
            f"k_final={scaled['k'][-1]} "
            f"ratio_invariant={ratio_invariant} "
            f"k1_invariant={k1_invariant} "
            f"k_invariant={k_invariant} "
            f"t_invariant={t_invariant} "
            f"delta_invariant={delta_invariant}"
        )

    # ------------------------------------------------------------------
    # 4. EXPLICIT 17 SCALING
    # ------------------------------------------------------------------

    q1_17 = int(
        17 * Q1
    )

    q3_17 = int(
        17 * Q3
    )

    scaled17 = reconstruct_general(
        q1_17,
        q3_17,
        END_E,
    )

    print()
    print("=" * 78)
    print("4. EXPLICIT 17-ADIC CONTENT / 7-ADIC INVARIANCE")
    print("=" * 78)

    print(
        f"  scaled_gcd="
        f"{gcd_int(q1_17,q3_17)}"
    )

    print(
        f"  scaled_v17(q1)="
        f"{valuation_17(q1_17)}"
    )

    print(
        f"  scaled_v17(q3)="
        f"{valuation_17(q3_17)}"
    )

    print(
        f"  original_k1={original['k1']}"
    )

    print(
        f"  scaled17_k1={scaled17['k1']}"
    )

    print(
        f"  original_k_final={original['k'][-1]}"
    )

    print(
        f"  scaled17_k_final={scaled17['k'][-1]}"
    )

    print(
        f"  k_sequence_unchanged="
        f"{scaled17['k'] == original['k']}"
    )

    print(
        f"  t_digits_unchanged="
        f"{scaled17['t'] == original['t']}"
    )

    print(
        f"  delta_digits_unchanged="
        f"{scaled17['delta'] == original['delta']}"
    )

    # ------------------------------------------------------------------
    # 5. TERMINAL SCHUR SOURCE SCALING
    # ------------------------------------------------------------------

    original_s = terminal_schur_triplet(
        Q1,
        Q3,
    )

    scaled_s = terminal_schur_triplet(
        q1_17,
        q3_17,
    )

    expected_s = [
        int(
            17 * x
        )
        for x in original_s
    ]

    schur_scaling_exact = (
        scaled_s
        == expected_s
    )

    print()
    print("=" * 78)
    print("5. TERMINAL SCHUR SOURCE SCALING")
    print("=" * 78)

    print(
        f"  original_terminal_row={original_s}"
    )

    print(
        f"  scaled_terminal_row={scaled_s}"
    )

    print(
        f"  expected_17x_row={expected_s}"
    )

    print(
        f"  exact_17_scaling="
        f"{schur_scaling_exact}"
    )

    # ------------------------------------------------------------------
    # 6. ONE-SIDED CONTROL
    # ------------------------------------------------------------------

    one_q1 = int(
        17 * Q1
    )

    one_q3 = int(
        Q3
    )

    one = reconstruct_general(
        one_q1,
        one_q3,
        END_E,
    )

    original_rho7 = ratio_mod(
        Q1,
        Q3,
        7,
    )

    one_rho7 = ratio_mod(
        one_q1,
        one_q3,
        7,
    )

    ratio_changed = (
        one_rho7
        != original_rho7
    )

    k1_changed = (
        one["k1"]
        != original["k1"]
    )

    k_changed = (
        one["k"]
        != original["k"]
    )

    t_changed = (
        one["t"]
        != original["t"]
    )

    print()
    print("=" * 78)
    print("6. ONE-SIDED SOURCE-SCALING CONTROL")
    print("=" * 78)

    print(
        f"  one_sided_q1={one_q1}"
    )

    print(
        f"  one_sided_q3={one_q3}"
    )

    print(
        f"  original_rho_mod7={original_rho7}"
    )

    print(
        f"  one_sided_rho_mod7={one_rho7}"
    )

    print(
        f"  ratio_changed={ratio_changed}"
    )

    print(
        f"  original_k1={original['k1']}"
    )

    print(
        f"  one_sided_k1={one['k1']}"
    )

    print(
        f"  k1_changed={k1_changed}"
    )

    print(
        f"  k_sequence_changed={k_changed}"
    )

    print(
        f"  t_digits_changed={t_changed}"
    )

    print(
        f"  one_sided_t_digits={one['t']}"
    )

    print(
        f"  one_sided_slopes={one['slope']}"
    )

    # ------------------------------------------------------------------
    # 7. SLOPE COMPARISON
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. GENERAL LIFT-SLOPE CONTROL")
    print("=" * 78)

    original_slope_one = all(
        slope == 1
        for slope in original["slope"]
    )

    one_sided_nontrivial_slope = any(
        slope != 1
        for slope in one["slope"]
    )

    print(
        f"  original_all_slopes_1="
        f"{original_slope_one}"
    )

    print(
        f"  one_sided_any_slope_not_1="
        f"{one_sided_nontrivial_slope}"
    )

    print(
        f"  original_slopes={original['slope']}"
    )

    print(
        f"  one_sided_slopes={one['slope']}"
    )

    # ------------------------------------------------------------------
    # 8. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
        """
The failure of the previous control was caused by assuming the
original slope-one formula

    t = -delta mod 7

for a changed projective ratio.

That shortcut is special to the original source ratio.

The correct general statement is:

    At every level there are seven possible exponent lifts,
    and exactly one satisfies the next congruence.

For the original pair the normalized derivative happens to be 1 mod 7,
so the lift digit has the simple form

    t = -delta.

For the one-sided control the projective ratio changes, so the local
slope may change. The seven-candidate reconstruction handles that case
without assumptions.

The common-scaling experiment is therefore a clean projective/content
control:

    common source content
        -> unchanged 7-adic exponent orbit,

while

    one-sided source modification
        -> changed projective orbit and potentially changed lift slope.

The 17-scaling simultaneously tests linear scaling of the terminal
Schur source row.
"""
    )

    # ------------------------------------------------------------------
    # 9. FINAL
    # ------------------------------------------------------------------

    all_common_invariant = all(
        row["ratio"]
        and row["k1"]
        and row["k"]
        and row["t"]
        and row["delta"]
        for row in scale_rows
    )

    explicit_17 = (
        scaled17["k1"] == original["k1"]
        and scaled17["k"] == original["k"]
        and scaled17["t"] == original["t"]
        and scaled17["delta"] == original["delta"]
    )

    control_pass = (
        ratio_changed
        and k_changed
        and t_changed
    )

    final_ok = (
        q1_v7 == 0
        and q3_v7 == 0
        and all_common_invariant
        and explicit_17
        and schur_scaling_exact
        and control_pass
        and original_slope_one
    )

    print()
    print("=" * 78)
    print("9. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  all_common_scales_preserve_7adic_orbit="
        f"{all_common_invariant}"
    )

    print(
        f"  common_17_scaling_preserves_exponent_data="
        f"{explicit_17}"
    )

    print(
        f"  terminal_schur_row_scales_by_17="
        f"{schur_scaling_exact}"
    )

    print(
        f"  one_sided_ratio_changes="
        f"{ratio_changed}"
    )

    print(
        f"  one_sided_exponent_changes="
        f"{k_changed}"
    )

    print(
        f"  one_sided_digits_change="
        f"{t_changed}"
    )

    print(
        f"  original_slope_one="
        f"{original_slope_one}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 216RRR COMPLETE")


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