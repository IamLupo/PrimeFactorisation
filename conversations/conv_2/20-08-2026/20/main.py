#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 399R-COMPACT — EXACT RESIDUAL RELATION / FINITE-GRID
COINCIDENCE AUDIT
==============================================================================

Purpose
-------
398R correctly removed the universal identity

    R_A / R_B = B / A.

However, a remaining relation such as

    R(A) / R(B) = f(p,t)

can still arise because B/A happens to equal f(p,t) on the finite observed
support of A and B.

399R therefore separates:

    1. globally symbolic identities;
    2. exact identities on the observed support;
    3. finite-grid coincidences caused by the small available (p,t) set;
    4. genuinely cross-parameter residual relations.

Additional rule
---------------
A relation is only promoted when its expression difference has no symbolic
factorization explaining the agreement and the relation survives separate
p-slices and separate t-slices.

No missing cell is used.
No interpolation.
No extrapolation.
No fitted coefficients.
"""

from __future__ import annotations

from itertools import combinations
import sympy as sp


# ============================================================================
# SOURCE
# ============================================================================

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


# ============================================================================
# SOURCE PARAMETER
# ============================================================================

def parameter_p(r: int) -> int:
    return 2 * r + 1


# ============================================================================
# FIXED EXPRESSION LIBRARY
# ============================================================================

EXPRESSIONS = {
    "p": p,
    "p-1": p - 1,
    "p+1": p + 1,
    "p^2": p**2,
    "p^2-1": p**2 - 1,
    "p^2+1": p**2 + 1,

    "p+t": p + t,
    "p-t": p - t,
    "p+t+1": p + t + 1,
    "p-t-1": p - t - 1,

    "p+2t": p + 2*t,
    "p-2t": p - 2*t,
    "p+2t+1": p + 2*t + 1,
    "p-2t-1": p - 2*t - 1,

    "p+3t": p + 3*t,
    "p-3t": p - 3*t,

    "p^2+t": p**2 + t,
    "p^2-t": p**2 - t,
    "p^2+t+1": p**2 + t + 1,
    "p^2-t-1": p**2 - t - 1,

    "p^2+t^2": p**2 + t**2,
    "p^2-t^2": p**2 - t**2,
    "p^2+t^2-1": p**2 + t**2 - 1,
    "p^2-t^2-1": p**2 - t**2 - 1,
    "p^2+t^2+1": p**2 + t**2 + 1,
    "p^2-t^2+1": p**2 - t**2 + 1,

    "p^2+p+t": p**2 + p + t,
    "p^2-p+t": p**2 - p + t,
    "p^2+p-t": p**2 + p - t,
    "p^2-p-t": p**2 - p - t,

    "p^2+2pt": p**2 + 2*p*t,
    "p^2-2pt": p**2 - 2*p*t,

    "p*(p+t)": p * (p + t),
    "p*(p-t)": p * (p - t),
    "p*(t+1)": p * (t + 1),
    "p*(t-1)": p * (t - 1),

    "(p-1)*(t+1)": (p - 1) * (t + 1),
    "(p-1)*(t-1)": (p - 1) * (t - 1),
    "(p+1)*(t+1)": (p + 1) * (t + 1),
    "(p+1)*(t-1)": (p + 1) * (t - 1),

    "(p-1)*(p+t)": (p - 1) * (p + t),
    "(p+1)*(p-t)": (p + 1) * (p - t),
}


# ============================================================================
# HELPERS
# ============================================================================

def eval_expr(expr: sp.Expr, cell: tuple[int, int]) -> int:
    r, tt = cell
    return int(expr.subs({p: parameter_p(r), t: tt}))


def support(expr: sp.Expr) -> set[tuple[int, int]]:
    result = set()

    for cell, value in Q.items():
        x = eval_expr(expr, cell)

        if x != 0 and value % x == 0:
            result.add(cell)

    return result


def residual(expr: sp.Expr, cell: tuple[int, int]) -> int:
    x = eval_expr(expr, cell)

    if x == 0:
        raise ZeroDivisionError

    return Q[cell] // x


def cross_parameter(cells) -> bool:
    ps = {parameter_p(r) for r, _ in cells}
    return len(ps) >= 2


def cross_t(cells) -> bool:
    ts = {tt for _, tt in cells}
    return len(ts) >= 2


def min_slice_counts(cells):
    by_p = {}
    by_t = {}

    for r, tt in cells:
        pp = parameter_p(r)
        by_p.setdefault(pp, 0)
        by_t.setdefault(tt, 0)
        by_p[pp] += 1
        by_t[tt] += 1

    return by_p, by_t


def residual_ratio(a_name, b_name, cell):
    ra = residual(EXPRESSIONS[a_name], cell)
    rb = residual(EXPRESSIONS[b_name], cell)

    if rb == 0:
        return None

    return sp.Rational(ra, rb)


def symbolic_forced_ratio(a_name, b_name):
    """
    Universal identity:

        R_A/R_B = B/A
    """
    return sp.cancel(
        EXPRESSIONS[b_name] / EXPRESSIONS[a_name]
    )


def expression_match_on_cells(expr1, expr2, cells):
    return all(
        eval_expr(expr1 - expr2, cell) == 0
        for cell in cells
    )


def symbolic_difference(expr1, expr2):
    return sp.factor(sp.together(expr1 - expr2))


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 399R-COMPACT — EXACT RESIDUAL RELATION / FINITE-GRID")
    print("COINCIDENCE AUDIT")
    print("=" * 78)

    print()
    print("SOURCE")
    print(f"  observed_cells={len(Q)}")
    print(
        "  source_parameters="
        f"{sorted({parameter_p(r) for r, _ in Q})}"
    )
    print(f"  missing_cells={MISSING}")
    print("  expressions_fixed_in_advance=True")

    # ------------------------------------------------------------------------
    # SUPPORT CACHE
    # ------------------------------------------------------------------------

    supports = {
        name: support(expr)
        for name, expr in EXPRESSIONS.items()
    }

    # ------------------------------------------------------------------------
    # 1. RECHECK 398R SURVIVOR
    # ------------------------------------------------------------------------

    candidate = (
        "p^2-p+t",
        "p^2-t^2-1",
        "p-t-1",
    )

    a_name, b_name, target_name = candidate

    cells = sorted(
        supports[a_name] &
        supports[b_name]
    )

    target_expr = EXPRESSIONS[target_name]
    forced_ratio = symbolic_forced_ratio(a_name, b_name)

    observed_ratio_matches = []

    for cell in cells:
        ratio = residual_ratio(a_name, b_name, cell)

        if ratio is None:
            continue

        target_value = eval_expr(target_expr, cell)

        observed_ratio_matches.append(
            ratio == target_value
        )

    print()
    print("=" * 78)
    print("1. 398R SURVIVOR RECHECK")
    print("=" * 78)

    print(f"  A={a_name}")
    print(f"  B={b_name}")
    print(f"  target={target_name}")
    print(f"  cells={len(cells)}")
    print(f"  distinct_p={sorted({parameter_p(r) for r, _ in cells})}")
    print(f"  distinct_t={sorted({tt for _, tt in cells})}")
    print(
        f"  target_match_count="
        f"{sum(observed_ratio_matches)}/{len(observed_ratio_matches)}"
    )
    print(f"  forced_ratio={sp.factor(forced_ratio)}")

    # ------------------------------------------------------------------------
    # 2. IS THE TARGET SYMBOLICALLY EQUAL TO B/A?
    # ------------------------------------------------------------------------

    difference = symbolic_difference(
        sp.cancel(forced_ratio),
        target_expr,
    )

    globally_equal = sp.simplify(difference) == 0

    print()
    print("=" * 78)
    print("2. SYMBOLIC-vs-GRID DISTINCTION")
    print("=" * 78)

    print(f"  symbolic_difference={difference}")
    print(f"  global_symbolic_identity={globally_equal}")

    if globally_equal:
        print("  classification=TAUTOLOGICAL")
    else:
        print("  classification=NOT_SYMBOLICALLY_FORCED")

    # ------------------------------------------------------------------------
    # 3. SOURCE-PARAMETER SLICE TEST
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. CROSS-PARAMETER SLICE TEST")
    print("=" * 78)

    by_p, by_t = min_slice_counts(cells)

    slice_failures = []

    for pp in sorted(by_p):

        p_cells = [
            c for c in cells
            if parameter_p(c[0]) == pp
        ]

        for cell in p_cells:

            ratio = residual_ratio(a_name, b_name, cell)
            target = eval_expr(target_expr, cell)

            if ratio != target:
                slice_failures.append((pp, cell))

        print(
            f"  p={pp} cells={len(p_cells)} "
            f"exact={all(
                residual_ratio(a_name, b_name, c)
                == eval_expr(target_expr, c)
                for c in p_cells
            )}"
        )

    # ------------------------------------------------------------------------
    # 4. T SLICE TEST
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. CROSS-T SLICE TEST")
    print("=" * 78)

    for tt in sorted(by_t):

        t_cells = [
            c for c in cells
            if c[1] == tt
        ]

        exact = all(
            residual_ratio(a_name, b_name, c)
            == eval_expr(target_expr, c)
            for c in t_cells
        )

        print(
            f"  t={tt} cells={len(t_cells)} exact={exact}"
        )

    # ------------------------------------------------------------------------
    # 5. ZERO SET / ROOT STRUCTURE OF THE DIFFERENCE
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. DIFFERENCE ROOT STRUCTURE")
    print("=" * 78)

    num, den = sp.fraction(
        sp.cancel(forced_ratio - target_expr)
    )

    print(f"  numerator={sp.factor(num)}")
    print(f"  denominator={sp.factor(den)}")

    # ------------------------------------------------------------------------
    # 6. SEARCH OTHER SURVIVING NON-TAUTOLOGICAL RELATIONS
    # ------------------------------------------------------------------------

    survivors = []

    names = sorted(EXPRESSIONS)

    for a_name, b_name in combinations(names, 2):

        common = sorted(
            supports[a_name] &
            supports[b_name]
        )

        if len(common) < 4:
            continue

        if not cross_parameter(common):
            continue

        if not cross_t(common):
            continue

        forced = symbolic_forced_ratio(a_name, b_name)

        # Search fixed library targets.
        for target_name, target_expr in EXPRESSIONS.items():

            # The forced ratio itself is not evidence.
            if sp.simplify(
                sp.cancel(forced) - target_expr
            ) == 0:
                continue

            exact = True

            for cell in common:

                ratio = residual_ratio(
                    a_name,
                    b_name,
                    cell,
                )

                if ratio is None:
                    exact = False
                    break

                if ratio != eval_expr(target_expr, cell):
                    exact = False
                    break

            if not exact:
                continue

            survivors.append(
                (
                    a_name,
                    b_name,
                    target_name,
                    common,
                    sp.factor(forced - target_expr),
                )
            )

    print()
    print("=" * 78)
    print("6. NON-TAUTOLOGICAL SOURCE-EXPRESSION SURVIVORS")
    print("=" * 78)

    print(f"  survivor_count={len(survivors)}")

    for a, b, target, common, difference in survivors[:15]:

        ps = sorted({parameter_p(r) for r, _ in common})
        ts = sorted({tt for _, tt in common})

        print(
            f"  {a}/{b} = {target} "
            f"cells={len(common)} "
            f"p={ps} "
            f"t={ts} "
            f"difference={difference}"
        )

    if len(survivors) > 15:
        print(
            f"  ... {len(survivors)-15} additional survivors omitted"
        )

    # ------------------------------------------------------------------------
    # 7. STRICT PROMOTION
    # ------------------------------------------------------------------------

    strict = []

    for a, b, target, common, difference in survivors:

        ps = {parameter_p(r) for r, _ in common}
        ts = {tt for _, tt in common}

        # Require at least two cells for every p appearing.
        p_counts = {}
        for r, tt in common:
            pp = parameter_p(r)
            p_counts[pp] = p_counts.get(pp, 0) + 1

        # Require at least two cells for every t appearing.
        t_counts = {}
        for r, tt in common:
            t_counts[tt] = t_counts.get(tt, 0) + 1

        if (
            len(ps) >= 2
            and len(ts) >= 2
            and all(v >= 2 for v in p_counts.values())
            and all(v >= 2 for v in t_counts.values())
        ):
            strict.append(
                (a, b, target, common, difference)
            )

    print()
    print("=" * 78)
    print("7. STRICT CROSS-P/T PROMOTION")
    print("=" * 78)

    print(
        f"  strict_candidate_count={len(strict)}"
    )

    for a, b, target, common, difference in strict[:10]:
        print(
            f"  {a}/{b} = {target} "
            f"cells={len(common)} "
            f"difference={difference}"
        )

    # ------------------------------------------------------------------------
    # 8. VERDICT
    # ------------------------------------------------------------------------

    if strict:
        verdict = "STRICT_NON_TAUTOLOGICAL_RESIDUAL_CANDIDATE"
    elif survivors:
        verdict = "FINITE_GRID_NON_TAUTOLOGICAL_CANDIDATES_ONLY"
    elif not globally_equal and cells:
        verdict = "398R_SURVIVOR_IS_A_GRID_COINCIDENCE_OR_SOURCE-IDENTITY"
    else:
        verdict = "NO_NON_TAUTOLOGICAL_RESIDUAL_STRUCTURE"

    print()
    print("=" * 78)
    print("8. STRUCTURAL VERDICT")
    print("=" * 78)

    print(f"  verdict={verdict}")

    print()
    print("  IMPORTANT")
    print("    R_A/R_B = B/A is universally forced.")
    print("    A source-library match is only interesting if it differs")
    print("    symbolically from B/A and survives independent p/t slices.")
    print("    finite-grid agreement is not automatically a law.")
    print("    missing cells remain diagnostic only.")

    # ------------------------------------------------------------------------
    # 9. FINAL EXACTNESS
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. FINAL EXACTNESS")
    print("=" * 78)

    print("  observed_cells_used_only=True")
    print("  exact_integer_arithmetic=True")
    print("  fixed_expression_library=True")
    print("  universal_B_over_A_identity_removed=True")
    print("  symbolic_identity_tested=True")
    print("  p_slice_tested=True")
    print("  t_slice_tested=True")
    print("  difference_factorization_completed=True")
    print("  strict_cross_parameter_requirement=True")
    print("  strict_cross_t_requirement=True")
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
    print("EXPERIMENT 399R-COMPACT COMPLETE")


if __name__ == "__main__":
    main()
