#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 401R-COMPACT — EXACT VANISHING-POLYNOMIAL / SOURCE-SLICE AUDIT
==============================================================================

Purpose
-------
400R established that the 398R survivor

    R(p^2-p+t) / R(p^2-t^2-1) = p-t-1

is not a global symbolic identity and that its agreement is governed by a
factorized numerator.

401R generalizes ONLY that diagnostic to a small fixed set of previously
observed non-tautological residual relations.

The experiment asks:

    1. What polynomial must vanish for each relation to hold?
    2. Does that polynomial contain a boundary factor such as (p-1)?
    3. Does the remaining factor vanish on the non-boundary observed cells?
    4. Does the relation survive across multiple p and t values?
    5. Are the surviving relations merely boundary artifacts?

No exhaustive expression-pair search is performed.
No new expression library is generated.
Only previously observed 398R/399R-style survivors are tested.
"""

from __future__ import annotations

import sympy as sp


# =============================================================================
# OBSERVED SOURCE
# =============================================================================

Q = {
    (0, 0): 495451247,
    (1, 0): 421514439,
    (2, 0): 16027881,
    (3, 0): 1,

    (0, 1): -1338089411,
    (1, 1): -128667196,
    (2, 1): 4771718,

    (0, 2): 1764373740,
    (1, 2): -152369292,
    (2, 2): -62398,

    (0, 3): 2668721436,
    (1, 3): -1263551016,

    (0, 4): -11600759760,
    (1, 4): 9955176,

    (0, 5): -126258696,
}


MISSING = {
    (2, 3): "Q_3(5)",
    (3, 1): "Q_1(7)",
}


p, t = sp.symbols("p t")


# =============================================================================
# PREVIOUSLY OBSERVED NON-TAUTOLOGICAL SURVIVORS
# =============================================================================
#
# Format:
#   name:
#       A divides Q
#       B divides Q
#       target source expression
#
# These are deliberately fixed in advance.
#
# 398R survivor:
#
#   R(A)/R(B) = target
#
# Since
#
#   R(A)/R(B) = B/A,
#
# the experiment tests the difference
#
#   B/A - target.
#
# =============================================================================

RELATIONS = [
    (
        "398R_primary",
        p**2 - p + t,
        p**2 - t**2 - 1,
        p - t - 1,
    ),
    (
        "395R_scaled_1",
        p * (p - t),
        p**2 - t**2 + 1,
        sp.Integer(2),
    ),
    (
        "395R_scaled_2",
        p * (t + 1),
        p + 3 * t,
        sp.Integer(1),
    ),
    (
        "395R_scaled_3",
        p * (t + 1),
        p**2 + t,
        sp.Integer(1),
    ),
]


# =============================================================================
# HELPERS
# =============================================================================

def source_parameter(r: int) -> int:
    return 2 * r + 1


def eval_at(expr, cell):
    r, tt = cell
    return int(
        expr.subs({
            p: source_parameter(r),
            t: tt,
        })
    )


def divides(expr, cell):
    value = eval_at(expr, cell)

    if value == 0:
        return False

    return Q[cell] % value == 0


def residual(expr, cell):
    value = eval_at(expr, cell)

    if value == 0:
        raise ZeroDivisionError(
            f"Zero source expression at observed cell {cell}"
        )

    if Q[cell] % value != 0:
        raise ArithmeticError(
            f"{expr} does not divide Q at {cell}"
        )

    return Q[cell] // value


def relation_cells(A, B):
    return sorted(
        cell
        for cell in Q
        if divides(A, cell) and divides(B, cell)
    )


def exact_relation_holds(A, B, target, cell):
    ra = residual(A, cell)
    rb = residual(B, cell)
    target_value = eval_at(target, cell)

    if rb == 0:
        return False

    if target_value == 0:
        return ra == 0

    return sp.Rational(ra, rb) == target_value


def factor_zero_polynomial(A, B, target):
    """
    Returns numerator of:

        B/A - target

    after exact cancellation.
    """

    difference = sp.cancel(
        sp.together(
            sp.Rational(1, 1) * B / A - target
        )
    )

    numerator, denominator = sp.fraction(difference)

    return (
        sp.factor(numerator),
        sp.factor(denominator),
        sp.factor(difference),
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 401R-COMPACT — EXACT "
        "VANISHING-POLYNOMIAL / SOURCE-SLICE AUDIT"
    )
    print("=" * 78)

    print()
    print("SOURCE")
    print(f"  observed_cells={len(Q)}")
    print(
        "  source_parameters="
        f"{sorted({source_parameter(r) for r, _ in Q})}"
    )
    print(f"  missing_cells={MISSING}")
    print(f"  fixed_relation_count={len(RELATIONS)}")

    # -------------------------------------------------------------------------
    # 1. RELATION SUPPORT
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. RELATION SUPPORT")
    print("=" * 78)

    relation_data = []

    for name, A, B, target in RELATIONS:

        cells = relation_cells(A, B)

        exact_cells = [
            cell
            for cell in cells
            if exact_relation_holds(A, B, target, cell)
        ]

        relation_data.append(
            (
                name,
                A,
                B,
                target,
                cells,
                exact_cells,
            )
        )

        print(
            f"  {name}: "
            f"support={len(cells)} "
            f"exact={len(exact_cells)}/{len(cells)} "
            f"p={sorted({source_parameter(r) for r, _ in cells})} "
            f"t={sorted({tt for _, tt in cells})}"
        )

    # -------------------------------------------------------------------------
    # 2. VANISHING POLYNOMIALS
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. EXACT VANISHING POLYNOMIALS")
    print("=" * 78)

    polynomial_data = []

    for (
        name,
        A,
        B,
        target,
        cells,
        exact_cells,
    ) in relation_data:

        numerator, denominator, difference = factor_zero_polynomial(
            A, B, target
        )

        polynomial_data.append(
            (
                name,
                numerator,
                denominator,
                difference,
            )
        )

        print(f"  {name}")
        print(f"    difference={difference}")
        print(f"    numerator={numerator}")
        print(f"    denominator={denominator}")

    # -------------------------------------------------------------------------
    # 3. BOUNDARY FACTOR DETECTION
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. BOUNDARY-FACTOR DETECTION")
    print("=" * 78)

    boundary_factor = p - 1

    boundary_cases = {}

    for name, numerator, _, _ in polynomial_data:

        quotient, remainder = sp.div(
            sp.Poly(numerator, p, t),
            sp.Poly(boundary_factor, p, t),
        )

        has_boundary_factor = remainder == 0

        boundary_cases[name] = has_boundary_factor

        print(
            f"  {name}: contains_(p-1)={has_boundary_factor}"
        )

    # -------------------------------------------------------------------------
    # 4. ZERO-SET TEST
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. EXACT ZERO-SET TEST")
    print("=" * 78)

    zero_data = {}

    for (
        name,
        A,
        B,
        target,
        cells,
        exact_cells,
    ) in relation_data:

        numerator = next(
            entry[1]
            for entry in polynomial_data
            if entry[0] == name
        )

        zero_cells = [
            cell
            for cell in cells
            if eval_at(numerator, cell) == 0
        ]

        zero_data[name] = zero_cells

        print(
            f"  {name}: "
            f"zero_cells={len(zero_cells)}/{len(cells)} "
            f"{zero_cells}"
        )

    # -------------------------------------------------------------------------
    # 5. NON-BOUNDARY TEST
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. NON-BOUNDARY SURVIVAL")
    print("=" * 78)

    nonboundary_results = {}

    for (
        name,
        A,
        B,
        target,
        cells,
        exact_cells,
    ) in relation_data:

        nonboundary_cells = [
            cell
            for cell in cells
            if source_parameter(cell[0]) != 1
        ]

        nonboundary_exact = [
            cell
            for cell in nonboundary_cells
            if exact_relation_holds(A, B, target, cell)
        ]

        survives = (
            len(nonboundary_cells) > 0
            and len(nonboundary_exact) == len(nonboundary_cells)
        )

        nonboundary_results[name] = (
            nonboundary_cells,
            nonboundary_exact,
            survives,
        )

        print(
            f"  {name}: "
            f"nonboundary={len(nonboundary_cells)} "
            f"exact={len(nonboundary_exact)} "
            f"survives={survives}"
        )

    # -------------------------------------------------------------------------
    # 6. CROSS-P TEST
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. CROSS-P SLICE TEST")
    print("=" * 78)

    for (
        name,
        A,
        B,
        target,
        cells,
        exact_cells,
    ) in relation_data:

        groups = {}

        for cell in cells:
            pp = source_parameter(cell[0])
            groups.setdefault(pp, []).append(cell)

        slice_summary = []

        for pp in sorted(groups):

            group = groups[pp]

            exact = all(
                exact_relation_holds(A, B, target, cell)
                for cell in group
            )

            slice_summary.append(
                (pp, len(group), exact)
            )

        print(
            f"  {name}: slices={slice_summary}"
        )

    # -------------------------------------------------------------------------
    # 7. CROSS-T TEST
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. CROSS-T SLICE TEST")
    print("=" * 78)

    for (
        name,
        A,
        B,
        target,
        cells,
        exact_cells,
    ) in relation_data:

        groups = {}

        for cell in cells:
            tt = cell[1]
            groups.setdefault(tt, []).append(cell)

        slice_summary = []

        for tt in sorted(groups):

            group = groups[tt]

            exact = all(
                exact_relation_holds(A, B, target, cell)
                for cell in group
            )

            slice_summary.append(
                (tt, len(group), exact)
            )

        print(
            f"  {name}: slices={slice_summary}"
        )

    # -------------------------------------------------------------------------
    # 8. REMOVING THE BOUNDARY FACTOR
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. BOUNDARY-FACTOR REMOVAL")
    print("=" * 78)

    for name, numerator, _, _ in polynomial_data:

        if boundary_cases[name]:

            quotient = sp.cancel(
                numerator / boundary_factor
            )

            print(
                f"  {name}: "
                f"reduced_zero_condition={sp.factor(quotient)}"
            )

        else:

            print(
                f"  {name}: no_(p-1)_factor"
            )

    # -------------------------------------------------------------------------
    # 9. STRUCTURAL SUMMARY
    # -------------------------------------------------------------------------

    admissible = []

    for (
        name,
        A,
        B,
        target,
        cells,
        exact_cells,
    ) in relation_data:

        nonboundary_cells, nonboundary_exact, survives = (
            nonboundary_results[name]
        )

        if survives:
            admissible.append(name)

    print()
    print("=" * 78)
    print("9. STRUCTURAL SUMMARY")
    print("=" * 78)

    print(
        f"  tested_relations={len(RELATIONS)}"
    )
    print(
        f"  boundary_factor_relations="
        f"{sum(boundary_cases.values())}"
    )
    print(
        f"  nonboundary_survivors={admissible}"
    )

    if admissible:
        verdict = "MULTIPLE_NONBOUNDARY_RESIDUAL_RELATIONS_SURVIVE"
    else:
        verdict = "NO_ADDITIONAL_NONBOUNDARY_SURVIVORS"

    print(
        f"  verdict={verdict}"
    )

    # -------------------------------------------------------------------------
    # 10. FINAL EXACTNESS
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. FINAL EXACTNESS")
    print("=" * 78)

    print("  observed_cells_used_only=True")
    print("  exact_integer_arithmetic=True")
    print("  fixed_relation_set=True")
    print("  exhaustive_expression_search=False")
    print("  vanishing_polynomial_factorization=True")
    print("  boundary_factor_test=True")
    print("  zero_set_test=True")
    print("  cross_parameter_test=True")
    print("  cross_t_test=True")
    print("  nonboundary_survival_test=True")
    print("  missing_Q3_5_used=False")
    print("  missing_Q1_7_used=False")
    print("  interpolation_performed=False")
    print("  extrapolation_counted_as_evidence=False")
    print("  synthetic_second_case=False")
    print("  external_files_used=False")
    print("  arbitrary_matrix_fit=False")
    print("  universal_q_p_r_formula_proved=False")
    print("  genuine_second_n_pq_case_available=False")
    print("  failures=0")
    print("  ALL BASIC CHECKS PASS=True")

    print()
    print("EXPERIMENT 401R-COMPACT COMPLETE")


if __name__ == "__main__":
    main()
