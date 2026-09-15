#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 217 — EXACT GENERAL MOD-7 HENSEL-SLOPE FORMULA AUDIT
==============================================================================

Candidate general formula:

    B(rho)
      = rho * ((729-1)/7)
      = 6*rho
      (mod 7).

Experiments 211 and 216RRR observed:

    original rho = 6  -> slope 1,
    one-sided rho = 4 -> slope 3.

Experiment 217 removes dependence on the original dataset.

For every nonzero residue

    rho in F_7^* = {1,2,3,4,5,6},

we construct the exact source pair

    q1 = rho,
    q3 = 1,

so that

    q1/q3 = rho.

For each source ratio we:

    1. find the unique exponent class k1 modulo 6;
    2. perform exact Hensel lifting;
    3. enumerate all seven possible lifts;
    4. measure the normalized residual law
           R(t) = A + B*t mod 7;
    5. verify affine exactness;
    6. verify uniqueness of the successful lift;
    7. compare measured B with

           B_pred = rho * ((729-1)/7) mod 7
                 = 6*rho mod 7.

The test is repeated through several 7-adic precision levels.

No floating point.
No SymPy.
No hard-coded exponent sequence.
"""


from __future__ import annotations

import sys


# ============================================================================
# CONSTANTS
# ============================================================================

P = 7

PRINCIPAL_BASE = 729

# Since k = k1 + 6m,
# the order of 3 modulo 7^e is 6*7^(e-1).
END_E = 9

NONZERO_RESIDUES = [1, 2, 3, 4, 5, 6]


# ============================================================================
# BASIC EXACT HELPERS
# ============================================================================

def inverse_mod(a: int, m: int) -> int:
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

    return int(old_s % m)


def ratio_mod(
    q1: int,
    q3: int,
    modulus: int,
) -> int:

    return int(
        (
            int(q1)
            * inverse_mod(
                int(q3),
                int(modulus),
            )
        )
        % int(modulus)
    )


def orbit_mod(
    k: int,
    modulus: int,
) -> int:

    return int(
        (
            2
            * pow(
                3,
                int(k),
                int(modulus),
            )
        )
        % int(modulus)
    )


def order_3_mod_7e(
    e: int,
) -> int:

    return int(
        6
        * (
            P ** (
                int(e) - 1
            )
        )
    )


# ============================================================================
# BASE EXPONENT
# ============================================================================

def find_k1(
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

        if (
            orbit_mod(
                k,
                7,
            )
            == rho
        ):
            matches.append(
                int(k)
            )

    if len(matches) != 1:
        raise ArithmeticError(
            f"Expected one base exponent for rho={rho}, "
            f"got {matches}."
        )

    return int(matches[0])


# ============================================================================
# GENERAL SEVEN-CANDIDATE LIFT
# ============================================================================

def lift_candidates(
    current_k: int,
    e: int,
) -> list[int]:

    order = order_3_mod_7e(
        e
    )

    return [
        int(
            current_k
            + t * order
        )
        for t in range(7)
    ]


def normalized_residual_digit(
    q1: int,
    q3: int,
    e: int,
    k: int,
) -> int:

    modulus_e = int(
        P ** int(e)
    )

    modulus_next = int(
        P ** (
            int(e) + 1
        )
    )

    target = ratio_mod(
        q1,
        q3,
        modulus_next,
    )

    orbit = orbit_mod(
        k,
        modulus_next,
    )

    residual = int(
        (
            orbit
            - target
        )
        % modulus_next
    )

    if residual % modulus_e != 0:
        raise ArithmeticError(
            f"k={k} does not solve the equation "
            f"modulo 7^{e}."
        )

    return int(
        (
            residual
            // modulus_e
        )
        % P
    )


# ============================================================================
# SINGLE SOURCE-RATIO AUDIT
# ============================================================================

def audit_ratio(
    rho: int,
    end_e: int,
):

    q1 = int(rho)
    q3 = 1

    k1 = find_k1(
        q1,
        q3,
    )

    current_k = int(
        k1
    )

    predicted_slope = int(
        (
            rho
            * (
                (PRINCIPAL_BASE - 1)
                // P
            )
        )
        % P
    )

    slopes = []
    intercepts = []
    digits = []
    exponents = [
        current_k
    ]

    all_affine = True
    all_unique = True
    all_alignment = True

    rows = []

    for e in range(
        1,
        end_e,
    ):

        candidates = lift_candidates(
            current_k,
            e,
        )

        residuals = []

        for candidate in candidates:

            value = normalized_residual_digit(
                q1,
                q3,
                e,
                candidate,
            )

            residuals.append(
                value
            )

        A = int(
            residuals[0]
        )

        B = int(
            (
                residuals[1]
                - residuals[0]
            )
            % P
        )

        affine = True

        for t in range(7):

            expected = int(
                (
                    A
                    + B * t
                )
                % P
            )

            if residuals[t] != expected:
                affine = False
                break

        matching_t = []

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

        for t in range(7):

            candidate = int(
                candidates[t]
            )

            orbit = orbit_mod(
                candidate,
                modulus_next,
            )

            if orbit == target:
                matching_t.append(
                    int(t)
                )

        unique = (
            len(matching_t) == 1
        )

        if unique:
            chosen_t = int(
                matching_t[0]
            )
        else:
            chosen_t = -1

        if chosen_t >= 0:

            next_k = int(
                candidates[
                    chosen_t
                ]
            )

            exact_alignment = (
                orbit_mod(
                    next_k,
                    modulus_next,
                )
                == target
            )

        else:
            next_k = int(
                current_k
            )

            exact_alignment = False

        slopes.append(
            B
        )

        intercepts.append(
            A
        )

        digits.append(
            chosen_t
        )

        exponents.append(
            next_k
        )

        all_affine = (
            all_affine
            and affine
        )

        all_unique = (
            all_unique
            and unique
        )

        all_alignment = (
            all_alignment
            and exact_alignment
        )

        rows.append(
            {
                "e": int(e),
                "current_k": int(current_k),
                "A": A,
                "B": B,
                "predicted_B": predicted_slope,
                "residuals": residuals,
                "matching_t": matching_t,
                "chosen_t": chosen_t,
                "next_k": next_k,
                "affine": affine,
                "unique": unique,
                "alignment": exact_alignment,
            }
        )

        current_k = int(
            next_k
        )

    slopes_match = all(
        B == predicted_slope
        for B in slopes
    )

    return {
        "rho": int(rho),
        "k1": int(k1),
        "predicted_slope": int(
            predicted_slope
        ),
        "slopes": slopes,
        "intercepts": intercepts,
        "digits": digits,
        "exponents": exponents,
        "rows": rows,
        "all_affine": all_affine,
        "all_unique": all_unique,
        "all_alignment": all_alignment,
        "slopes_match": slopes_match,
    }


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 217 — EXACT GENERAL MOD-7 "
        "HENSEL-SLOPE FORMULA AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------
    # 1. THEORETICAL CONSTANT
    # ------------------------------------------------------------------

    normalized_increment = int(
        (
            PRINCIPAL_BASE - 1
        )
        // P
    )

    increment_mod7 = int(
        normalized_increment % P
    )

    print()
    print("=" * 78)
    print("1. PRINCIPAL-LINEARIZATION CONSTANT")
    print("=" * 78)

    print(
        f"  principal_base={PRINCIPAL_BASE}"
    )

    print(
        f"  principal_base_minus_1="
        f"{PRINCIPAL_BASE - 1}"
    )

    print(
        f"  (729-1)/7={normalized_increment}"
    )

    print(
        f"  ((729-1)/7)_mod7="
        f"{increment_mod7}"
    )

    print(
        "  predicted_general_formula="
        "B = rho * 6 mod 7"
    )

    # ------------------------------------------------------------------
    # 2. ALL NONZERO RESIDUES
    # ------------------------------------------------------------------

    results = []

    print()
    print("=" * 78)
    print("2. COMPLETE F_7^* SLOPE AUDIT")
    print("=" * 78)

    for rho in NONZERO_RESIDUES:

        result = audit_ratio(
            rho,
            END_E,
        )

        results.append(
            result
        )

        print()
        print(
            f"  rho={rho}"
        )

        print(
            f"    k1={result['k1']}"
        )

        print(
            f"    predicted_B="
            f"{result['predicted_slope']}"
        )

        print(
            f"    measured_slopes="
            f"{result['slopes']}"
        )

        print(
            f"    intercepts="
            f"{result['intercepts']}"
        )

        print(
            f"    lift_digits="
            f"{result['digits']}"
        )

        print(
            f"    all_affine="
            f"{result['all_affine']}"
        )

        print(
            f"    all_unique="
            f"{result['all_unique']}"
        )

        print(
            f"    all_alignment="
            f"{result['all_alignment']}"
        )

        print(
            f"    slope_formula_exact="
            f"{result['slopes_match']}"
        )

    # ------------------------------------------------------------------
    # 3. EXPLICIT FORMULA TABLE
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. EXPLICIT SLOPE FORMULA TABLE")
    print("=" * 78)

    formula_rows = []

    for result in results:

        rho = result["rho"]

        expected = int(
            (
                rho
                * increment_mod7
            )
            % P
        )

        measured_values = sorted(
            set(
                result["slopes"]
            )
        )

        exact = (
            measured_values
            == [expected]
        )

        formula_rows.append(
            exact
        )

        print(
            f"  rho={rho}: "
            f"rho*6 mod7={expected} "
            f"measured={measured_values} "
            f"exact={exact}"
        )

    # ------------------------------------------------------------------
    # 4. AFFINE LAW
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. GENERAL AFFINE RESIDUAL LAW")
    print("=" * 78)

    all_affine = all(
        result["all_affine"]
        for result in results
    )

    all_unique = all(
        result["all_unique"]
        for result in results
    )

    all_alignment = all(
        result["all_alignment"]
        for result in results
    )

    all_slope_formula = all(
        result["slopes_match"]
        for result in results
    )

    print(
        f"  all_ratios_affine="
        f"{all_affine}"
    )

    print(
        f"  all_ratios_unique_lifts="
        f"{all_unique}"
    )

    print(
        f"  all_ratios_exactly_aligned="
        f"{all_alignment}"
    )

    print(
        f"  all_ratios_match_general_slope_formula="
        f"{all_slope_formula}"
    )

    # ------------------------------------------------------------------
    # 5. SLOPE MAP ON F_7^*
    # ------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. SLOPE MAP rho -> B")
    print("=" * 78)

    slope_map = {}

    for result in results:

        rho = result["rho"]

        unique_slopes = sorted(
            set(
                result["slopes"]
            )
        )

        slope_map[
            rho
        ] = unique_slopes

        print(
            f"  rho={rho} "
            f"-> B={unique_slopes}"
        )

    # ------------------------------------------------------------------
    # 6. NONZERO-SLOPE PROPERTY
    # ------------------------------------------------------------------

    nonzero_slope_property = all(
        all(
            B != 0
            for B in result["slopes"]
        )
        for result in results
    )

    print()
    print("=" * 78)
    print("6. NONZERO-DERIVATIVE PROPERTY")
    print("=" * 78)

    print(
        f"  all_measured_slopes_nonzero="
        f"{nonzero_slope_property}"
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
For any nonzero source ratio rho modulo 7, let k1 be the unique
exponent satisfying

    2*3^k1 = rho (mod 7).

Then write

    k = k1 + 6m.

The orbit becomes

    2*3^k
      = (2*3^k1) * 729^m.

For one Hensel digit at level e,

    m -> m + t*7^(e-1).

The first-order change is controlled by

    729 - 1 = 7*104,

with

    104 = 6 (mod 7).

Therefore the normalized residual slope is

    B
      = rho * 6
      (mod 7).

Experiment 217 tests this formula for every nonzero residue rho in F_7.

This separates the universal local mechanism from the special source
ratio rho=6 of the original experiment.
"""
    )

    # ------------------------------------------------------------------
    # 8. FINAL
    # ------------------------------------------------------------------

    all_formula_rows = all(
        formula_rows
    )

    final_ok = (
        all_formula_rows
        and all_affine
        and all_unique
        and all_alignment
        and all_slope_formula
        and nonzero_slope_property
    )

    print()
    print("=" * 78)
    print("8. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  residues_tested="
        f"{NONZERO_RESIDUES}"
    )

    print(
        f"  normalized_729_increment_mod7="
        f"{increment_mod7}"
    )

    print(
        f"  all_affine_exact="
        f"{all_affine}"
    )

    print(
        f"  all_lift_digits_unique="
        f"{all_unique}"
    )

    print(
        f"  all_modular_alignments_exact="
        f"{all_alignment}"
    )

    print(
        f"  general_B_formula_exact="
        f"{all_slope_formula}"
    )

    print(
        f"  all_slopes_nonzero="
        f"{nonzero_slope_property}"
    )

    print(
        f"  failures={0 if final_ok else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS={final_ok}"
    )

    print()
    print("EXPERIMENT 217 COMPLETE")


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

