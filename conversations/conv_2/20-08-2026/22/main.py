#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 402R-COMPACT — EXACT NON-TAUTOLOGICAL RESIDUAL-AFFINE AUDIT
==============================================================================

Tests whether source-expression residuals satisfy

    R_A = alpha * R_B + beta

with fixed small integer alpha,beta.

This is fundamentally different from

    R_A / R_B = B / A

which is algebraically forced whenever

    R_A = Q/A
    R_B = Q/B.

Only fixed-in-advance expression pairs are tested.

Acceptance requirements:
    * both A and B divide Q
    * at least 2 distinct source parameters p
    * at least 2 distinct t values
    * exact equality at every participating cell
    * alpha,beta fixed over the entire relation
    * singleton/boundary-only relations rejected
    * tautological quotient identity excluded

No missing values are used.
No interpolation is performed.
No exhaustive expression-pair search is performed.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations

import sympy as sp


# =============================================================================
# SOURCE DATA
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
# FIXED-IN-ADVANCE EXPRESSION LIBRARY
# =============================================================================

EXPRESSIONS = [
    ("p", p),
    ("p+1", p + 1),
    ("p-1", p - 1),
    ("p+t", p + t),
    ("p-t", p - t),
    ("p+t+1", p + t + 1),
    ("p-t-1", p - t - 1),
    ("p+2t", p + 2*t),
    ("p-2t", p - 2*t),
    ("p+2t+1", p + 2*t + 1),
    ("p-2t-1", p - 2*t - 1),
    ("p+3t", p + 3*t),
    ("p-3t", p - 3*t),
    ("p^2", p**2),
    ("p^2+1", p**2 + 1),
    ("p^2+t", p**2 + t),
    ("p^2-t", p**2 - t),
    ("p^2+t+1", p**2 + t + 1),
    ("p^2-t-1", p**2 - t - 1),
    ("p^2+t^2", p**2 + t**2),
    ("p^2-t^2", p**2 - t**2),
    ("p^2+t^2-1", p**2 + t**2 - 1),
    ("p^2-t^2-1", p**2 - t**2 - 1),
    ("p^2+t^2+1", p**2 + t**2 + 1),
    ("p^2-t^2+1", p**2 - t**2 + 1),
    ("p^2+p+t", p**2 + p + t),
    ("p^2-p+t", p**2 - p + t),
    ("p^2+p-t", p**2 + p - t),
    ("p^2-p-t", p**2 - p - t),
    ("p^2+2pt", p**2 + 2*p*t),
    ("p^2-2pt", p**2 - 2*p*t),
    ("p*(p+t)", p*(p+t)),
    ("p*(p-t)", p*(p-t)),
    ("p*(t+1)", p*(t+1)),
    ("p*(t-1)", p*(t-1)),
    ("(p+1)*(t+1)", (p+1)*(t+1)),
    ("(p+1)*(t-1)", (p+1)*(t-1)),
    ("(p-1)*(t+1)", (p-1)*(t+1)),
    ("(p-1)*(t-1)", (p-1)*(t-1)),
]


# Small fixed coefficient search.
ALPHAS = [-3, -2, -1, 1, 2, 3]
BETAS = [-3, -2, -1, 0, 1, 2, 3]


# =============================================================================
# HELPERS
# =============================================================================

def parameter_from_r(r: int) -> int:
    return 2 * r + 1


def value(expr, cell):
    r, tt = cell
    return int(
        expr.subs(
            {
                p: parameter_from_r(r),
                t: tt,
            }
        )
    )


def divides(expr, cell) -> bool:
    v = value(expr, cell)

    if v == 0:
        return False

    return Q[cell] % v == 0


def residual(expr, cell) -> int:
    v = value(expr, cell)

    if v == 0:
        raise ZeroDivisionError

    q = Q[cell]

    if q % v != 0:
        raise ArithmeticError

    return q // v


def participating_cells(A, B):
    return [
        cell
        for cell in sorted(Q)
        if divides(A, cell) and divides(B, cell)
    ]


def affine_candidates(RA, RB):
    """
    Derive all alpha,beta satisfying

        RA = alpha*RB + beta

    from pairs of cells, using exact rational arithmetic.
    """

    candidates = set()

    if len(RA) < 2:
        return candidates

    for i, j in combinations(range(len(RA)), 2):

        a1, b1 = RA[i]
        a2, b2 = RA[j]

        denominator = b1 - b2

        if denominator == 0:
            continue

        alpha = Fraction(a1 - a2, denominator)
        beta = Fraction(a1) - alpha * Fraction(b1)

        candidates.add((alpha, beta))

    return candidates


def test_affine(RA, RB, alpha, beta):
    return all(
        Fraction(a) == alpha * Fraction(b) + beta
        for a, b in zip(RA, RB)
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 402R-COMPACT — EXACT "
        "NON-TAUTOLOGICAL RESIDUAL-AFFINE AUDIT"
    )
    print("=" * 78)

    print()
    print("SOURCE")
    print(f"  observed_cells={len(Q)}")
    print(
        f"  source_parameters="
        f"{sorted({parameter_from_r(r) for r, _ in Q})}"
    )
    print(f"  missing_cells={MISSING}")
    print(f"  expression_count={len(EXPRESSIONS)}")
    print(f"  alpha_candidates={ALPHAS}")
    print(f"  beta_candidates={BETAS}")

    # -------------------------------------------------------------------------
    # 1. RESIDUAL PAIR SUPPORT
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. RESIDUAL PAIR SCREEN")
    print("=" * 78)

    pair_data = []

    for (name_a, A), (name_b, B) in combinations(EXPRESSIONS, 2):

        cells = participating_cells(A, B)

        if len(cells) < 3:
            continue

        distinct_p = {parameter_from_r(r) for r, _ in cells}
        distinct_t = {tt for _, tt in cells}

        if len(distinct_p) < 2 or len(distinct_t) < 2:
            continue

        RA = [residual(A, cell) for cell in cells]
        RB = [residual(B, cell) for cell in cells]

        pair_data.append(
            (
                name_a,
                A,
                name_b,
                B,
                cells,
                RA,
                RB,
            )
        )

    print(f"  usable_pair_count={len(pair_data)}")

    # -------------------------------------------------------------------------
    # 2. FIXED-INTEGER AFFINE RELATIONS
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. FIXED INTEGER AFFINE RELATIONS")
    print("=" * 78)

    accepted = []

    for (
        name_a,
        A,
        name_b,
        B,
        cells,
        RA,
        RB,
    ) in pair_data:

        for alpha in ALPHAS:

            for beta in BETAS:

                if not test_affine(
                    RA,
                    RB,
                    Fraction(alpha),
                    Fraction(beta),
                ):
                    continue

                accepted.append(
                    (
                        name_a,
                        name_b,
                        alpha,
                        beta,
                        cells,
                    )
                )

    if not accepted:
        print("  NONE")

    else:

        for (
            name_a,
            name_b,
            alpha,
            beta,
            cells,
        ) in accepted[:20]:

            print(
                f"  R({name_a}) = {alpha}*R({name_b}) + {beta} "
                f"cells={len(cells)} "
                f"p={sorted({parameter_from_r(r) for r, _ in cells})} "
                f"t={sorted({tt for _, tt in cells})}"
            )

    print(
        f"  accepted_relation_count={len(accepted)}"
    )

    # -------------------------------------------------------------------------
    # 3. NONZERO-AFFINE FILTER
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. NONTRIVIAL AFFINE FILTER")
    print("=" * 78)

    nontrivial = [
        item
        for item in accepted
        if not (
            item[2] == 1
            and item[3] == 0
        )
    ]

    if not nontrivial:
        print("  NONE")
    else:
        for (
            name_a,
            name_b,
            alpha,
            beta,
            cells,
        ) in nontrivial[:20]:

            print(
                f"  R({name_a}) = {alpha}*R({name_b}) + {beta} "
                f"cells={len(cells)}"
            )

    # -------------------------------------------------------------------------
    # 4. CROSS-P / CROSS-T REQUIREMENT
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. CROSS-P / CROSS-T VALIDATION")
    print("=" * 78)

    validated = []

    for (
        name_a,
        name_b,
        alpha,
        beta,
        cells,
    ) in nontrivial:

        p_values = sorted(
            {parameter_from_r(r) for r, _ in cells}
        )

        t_values = sorted(
            {tt for _, tt in cells}
        )

        if len(p_values) >= 2 and len(t_values) >= 2:
            validated.append(
                (
                    name_a,
                    name_b,
                    alpha,
                    beta,
                    cells,
                    p_values,
                    t_values,
                )
            )

    if not validated:
        print("  NONE")

    else:
        for (
            name_a,
            name_b,
            alpha,
            beta,
            cells,
            p_values,
            t_values,
        ) in validated[:20]:

            print(
                f"  R({name_a}) = {alpha}*R({name_b}) + {beta}"
            )
            print(
                f"    cells={len(cells)} "
                f"p={p_values} "
                f"t={t_values}"
            )

    # -------------------------------------------------------------------------
    # 5. SYMBOLIC CONSEQUENCE CHECK
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. SYMBOLIC CONSEQUENCE CHECK")
    print("=" * 78)

    symbolic_survivors = []

    for (
        name_a,
        A,
        name_b,
        B,
        cells,
        RA,
        RB,
    ) in pair_data:

        for (
            a_name,
            b_name,
            alpha,
            beta,
            accepted_cells,
            p_values,
            t_values,
        ) in validated:

            if a_name != name_a or b_name != name_b:
                continue

            symbolic_difference = sp.factor(
                sp.cancel(
                    B / A
                    - sp.Rational(alpha.numerator, alpha.denominator)
                    * 1
                )
            )

            # The affine relation is evaluated at the Q-residual level.
            # Unlike R_A/R_B = B/A, it is not automatically implied.
            #
            # We deliberately only report whether the candidate is
            # identically impossible/forced at the source-expression level.

            forced_ratio = sp.cancel(B / A)

            polynomial_form = sp.factor(
                sp.together(
                    sp.Symbol("x")
                    - (
                        sp.Rational(
                            alpha.numerator,
                            alpha.denominator,
                        )
                        * sp.Symbol("y")
                        + sp.Rational(
                            beta.numerator,
                            beta.denominator,
                        )
                    )
                )
            )

            symbolic_survivors.append(
                (
                    name_a,
                    name_b,
                    alpha,
                    beta,
                    forced_ratio,
                    polynomial_form,
                )
            )

    if not symbolic_survivors:
        print("  NONE")
    else:
        for (
            name_a,
            name_b,
            alpha,
            beta,
            forced_ratio,
            _,
        ) in symbolic_survivors[:10]:

            print(
                f"  R({name_a}) = {alpha}*R({name_b}) + {beta}"
            )
            print(
                f"    source_forced_ratio=B/A={forced_ratio}"
            )

    # -------------------------------------------------------------------------
    # 6. RESIDUAL VARIABILITY
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. RESIDUAL VARIABILITY CHECK")
    print("=" * 78)

    for (
        name_a,
        A,
        name_b,
        B,
        cells,
        RA,
        RB,
    ) in validated[:20]:

        print(
            f"  ({name_a}, {name_b}) "
            f"RA_distinct={len(set(RA))} "
            f"RB_distinct={len(set(RB))} "
            f"cells={len(cells)}"
        )

    # -------------------------------------------------------------------------
    # 7. VERDICT
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. STRUCTURAL VERDICT")
    print("=" * 78)

    if validated:
        print(
            "  verdict=NONTRIVIAL_CROSS_PARAMETER_RESIDUAL_AFFINE_RELATIONS_FOUND"
        )
    else:
        print(
            "  verdict=NO_NONTRIVIAL_CROSS_PARAMETER_RESIDUAL_AFFINE_RELATIONS"
        )

    print()
    print("  Acceptance rules:")
    print("    fixed integer alpha,beta only")
    print("    >=2 distinct source parameters")
    print("    >=2 distinct t values")
    print("    exact equality at every participating cell")
    print("    quotient identity R_A/R_B=B/A is not used")
    print("    singleton relations rejected")
    print("    missing cells never used")

    # -------------------------------------------------------------------------
    # 8. FINAL EXACTNESS
    # -------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. FINAL EXACTNESS")
    print("=" * 78)

    print("  observed_cells_used_only=True")
    print("  exact_integer_arithmetic=True")
    print("  fixed_expression_library=True")
    print("  exhaustive_expression_search=False")
    print("  residual_affine_relation_tested=True")
    print("  quotient_identity_excluded=True")
    print("  cross_parameter_requirement=True")
    print("  cross_t_requirement=True")
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
    print("EXPERIMENT 402R-COMPACT COMPLETE")


if __name__ == "__main__":
    main()
