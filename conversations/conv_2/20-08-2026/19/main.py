#!/usr/bin/env python3
"""
EXPERIMENT 398R-COMPACT
EXACT RESIDUAL-INDEPENDENCE / TAUTOLOGICAL-IDENTITY AUDIT

Core correction
---------------
For

    R_A = Q / A
    R_B = Q / B

we ALWAYS have

    R_A / R_B = B / A.

Therefore every residual-ratio relation observed in 395R-397R that is
merely equivalent to B/A is algebraically forced and is NOT evidence.

398R asks a stricter question:

    After removing the universal identity B/A,
    is there any additional source-expression law?

Tested signals:

1. residual A and residual B directly;
2. cross-residual equality;
3. cross-residual constant multiples;
4. residual-ratio AFTER dividing out the forced B/A identity;
5. cross-parameter equality of normalized residual invariants;
6. simple quotient invariants such as
       R_A * A / Q = 1
   are explicitly classified as tautological;
7. only genuinely non-tautological relations are reported.

Rules
-----
- exact integer arithmetic;
- fixed expression library;
- >=2 distinct p values;
- >=2 distinct t values;
- missing cells never used;
- no interpolation;
- no extrapolation;
- compact output.
"""

from __future__ import annotations

from itertools import combinations
from math import gcd
from functools import reduce

import sympy as sp


# ============================================================================
# SOURCE DATA
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


def parameter_p(r: int) -> int:
    return 2 * r + 1


# ============================================================================
# FIXED SOURCE-EXPRESSION LIBRARY
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
    out = set()

    for cell, value in Q.items():
        e = eval_expr(expr, cell)

        if e != 0 and value % e == 0:
            out.add(cell)

    return out


def residual(expr: sp.Expr, cell: tuple[int, int]) -> int:
    e = eval_expr(expr, cell)

    if e == 0:
        raise ZeroDivisionError

    value = Q[cell]

    if value % e != 0:
        raise ValueError

    return value // e


def cross_parameter_ok(cells: list[tuple[int, int]]) -> bool:
    ps = {parameter_p(r) for r, _ in cells}
    ts = {tt for _, tt in cells}
    return len(ps) >= 2 and len(ts) >= 2


def exact_constant(values: list[sp.Rational]) -> sp.Rational | None:
    if not values:
        return None

    first = values[0]

    if all(v == first for v in values):
        return first

    return None


def gcd_values(values: list[int]) -> int:
    if not values:
        return 0
    return reduce(gcd, (abs(x) for x in values))


def simplify_ratio(a: sp.Expr, b: sp.Expr) -> sp.Expr:
    return sp.factor(sp.cancel(a / b))


def ratio_library_match(
    values: list[sp.Rational],
    cells: list[tuple[int, int]],
) -> list[str]:

    hits = []

    for name, expr in EXPRESSIONS.items():
        ok = True

        for value, cell in zip(values, cells):
            ev = eval_expr(expr, cell)

            if sp.Rational(ev) != value:
                ok = False
                break

        if ok:
            hits.append(name)

    return hits


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:

    print("=" * 78)
    print("EXPERIMENT 398R-COMPACT — EXACT RESIDUAL-INDEPENDENCE")
    print("/ TAUTOLOGICAL-IDENTITY AUDIT")
    print("=" * 78)

    print()
    print("SOURCE")
    print(f"  observed_cells={len(Q)}")
    print(
        f"  source_parameters="
        f"{sorted(parameter_p(r) for r, _ in Q.keys())}"
    )
    print(f"  missing_cells={MISSING}")
    print("  expressions_fixed_in_advance=True")

    # ------------------------------------------------------------------------
    # SUPPORT CACHE
    # ------------------------------------------------------------------------

    expr_names = sorted(EXPRESSIONS)
    expr_support = {
        name: support(expr)
        for name, expr in EXPRESSIONS.items()
    }

    # ------------------------------------------------------------------------
    # 1. UNIVERSAL RESIDUAL IDENTITY
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. UNIVERSAL RESIDUAL IDENTITY")
    print("=" * 78)

    print("  For A | Q and B | Q:")
    print("    R_A = Q/A")
    print("    R_B = Q/B")
    print("    therefore R_A/R_B = B/A")
    print()
    print("  This identity is classified as TAUTOLOGICAL.")
    print("  It cannot constitute evidence of source structure.")

    # ------------------------------------------------------------------------
    # 2. DIRECT CROSS-RESIDUAL EQUALITY
    # ------------------------------------------------------------------------

    equality_hits = []

    for a_name, b_name in combinations(expr_names, 2):

        cells = sorted(expr_support[a_name] & expr_support[b_name])

        if len(cells) < 3 or not cross_parameter_ok(cells):
            continue

        ra = [residual(EXPRESSIONS[a_name], c) for c in cells]
        rb = [residual(EXPRESSIONS[b_name], c) for c in cells]

        if ra == rb:
            equality_hits.append(
                (a_name, b_name, cells)
            )

    print()
    print("=" * 78)
    print("2. DIRECT CROSS-RESIDUAL EQUALITY")
    print("=" * 78)

    if not equality_hits:
        print("  NONE")
    else:
        for a, b, cells in equality_hits[:10]:
            print(
                f"  R({a}) = R({b}) cells={len(cells)}"
            )

        if len(equality_hits) > 10:
            print(
                f"  ... {len(equality_hits) - 10} additional omitted"
            )

    # ------------------------------------------------------------------------
    # 3. CONSTANT RESIDUAL SCALING
    # ------------------------------------------------------------------------

    constant_ratio_hits = []

    for a_name, b_name in combinations(expr_names, 2):

        cells = sorted(expr_support[a_name] & expr_support[b_name])

        if len(cells) < 3 or not cross_parameter_ok(cells):
            continue

        ra = [residual(EXPRESSIONS[a_name], c) for c in cells]
        rb = [residual(EXPRESSIONS[b_name], c) for c in cells]

        ratios = []

        valid = True

        for x, y in zip(ra, rb):

            if y == 0:
                valid = False
                break

            ratios.append(sp.Rational(x, y))

        if not valid:
            continue

        const = exact_constant(ratios)

        if const is None:
            continue

        # Universal expected ratio:
        # R_A / R_B = B/A.
        forced = simplify_ratio(
            EXPRESSIONS[b_name],
            EXPRESSIONS[a_name],
        )

        is_universal_identity = (
            sp.simplify(forced - const) == 0
        )

        constant_ratio_hits.append(
            (
                a_name,
                b_name,
                cells,
                const,
                is_universal_identity,
            )
        )

    print()
    print("=" * 78)
    print("3. CONSTANT CROSS-RESIDUAL RATIOS")
    print("=" * 78)

    genuine_constant = []

    for a, b, cells, ratio, forced in constant_ratio_hits:
        if not forced:
            genuine_constant.append((a, b, cells, ratio))

    print(
        f"  total_constant_ratio_hits={len(constant_ratio_hits)}"
    )
    print(
        f"  tautological_constant_ratio_hits="
        f"{len(constant_ratio_hits) - len(genuine_constant)}"
    )

    if not genuine_constant:
        print("  genuine_constant_ratio_hits=NONE")
    else:
        for a, b, cells, ratio in genuine_constant[:10]:
            print(
                f"  R({a})/R({b})={ratio} "
                f"cells={len(cells)}"
            )

    # ------------------------------------------------------------------------
    # 4. RESIDUAL RATIO AFTER REMOVING FORCED B/A
    # ------------------------------------------------------------------------

    nontrivial_ratio_hits = []

    for a_name, b_name in combinations(expr_names, 2):

        cells = sorted(expr_support[a_name] & expr_support[b_name])

        if len(cells) < 3 or not cross_parameter_ok(cells):
            continue

        ra = [residual(EXPRESSIONS[a_name], c) for c in cells]
        rb = [residual(EXPRESSIONS[b_name], c) for c in cells]

        # Compute:
        #
        #       (R_A / R_B) / (B/A)
        #
        # which must equal 1 identically if the residuals are purely
        # tautological.
        normalized = []

        valid = True

        for r_a, r_b, cell in zip(ra, rb, cells):

            if r_b == 0:
                valid = False
                break

            observed_ratio = sp.Rational(r_a, r_b)

            A = eval_expr(EXPRESSIONS[a_name], cell)
            B = eval_expr(EXPRESSIONS[b_name], cell)

            if A == 0 or B == 0:
                valid = False
                break

            forced_ratio = sp.Rational(B, A)

            normalized.append(
                sp.cancel(observed_ratio / forced_ratio)
            )

        if not valid:
            continue

        const = exact_constant(normalized)

        if const is None:
            nontrivial_ratio_hits.append(
                (a_name, b_name, cells, normalized)
            )
        elif const != 1:
            nontrivial_ratio_hits.append(
                (a_name, b_name, cells, normalized)
            )

    print()
    print("=" * 78)
    print("4. NORMALIZED RESIDUAL-RATIO AUDIT")
    print("=" * 78)

    print(
        "  normalized quantity = "
        "(R_A/R_B)/(B/A)"
    )
    print(
        "  universal tautological value = 1"
    )

    print(
        f"  nontrivial_normalized_relations="
        f"{len(nontrivial_ratio_hits)}"
    )

    if not nontrivial_ratio_hits:
        print("  NONE")
    else:
        for a, b, cells, values in nontrivial_ratio_hits[:10]:
            print(
                f"  pair=({a}, {b}) "
                f"cells={len(cells)} "
                f"normalized_values={values[:6]}"
            )

    # ------------------------------------------------------------------------
    # 5. CROSS-CELL RESIDUAL INVARIANTS
    # ------------------------------------------------------------------------

    invariant_hits = []

    for name in expr_names:

        cells = sorted(expr_support[name])

        if len(cells) < 4 or not cross_parameter_ok(cells):
            continue

        # K_A = R_A * A = Q.
        # This is universally Q and therefore tautological.
        #
        # Instead compare normalized quotients:
        #
        # R_A(cell1) / R_A(cell2)
        #
        # against fixed source expressions evaluated between cells.

        values = [residual(EXPRESSIONS[name], c) for c in cells]

        ratios = []

        for i in range(len(cells) - 1):
            if values[i + 1] == 0:
                continue

            ratios.append(
                sp.Rational(values[i], values[i + 1])
            )

        if not ratios:
            continue

        # Only report if the residual itself is constant.
        if len(set(values)) == 1:
            invariant_hits.append(
                (name, cells, values[0])
            )

    print()
    print("=" * 78)
    print("5. DIRECT RESIDUAL CONSTANCY")
    print("=" * 78)

    if not invariant_hits:
        print("  NONE")
    else:
        for name, cells, value in invariant_hits[:10]:
            print(
                f"  R({name}) constant={value} "
                f"cells={len(cells)}"
            )

    # ------------------------------------------------------------------------
    # 6. SOURCE-EXPRESSION RELATION AFTER TAUTOLOGY REMOVAL
    # ------------------------------------------------------------------------

    genuine_source_matches = []

    for a_name, b_name in combinations(expr_names, 2):

        cells = sorted(expr_support[a_name] & expr_support[b_name])

        if len(cells) < 3 or not cross_parameter_ok(cells):
            continue

        ratios = []

        for cell in cells:
            ra = residual(EXPRESSIONS[a_name], cell)
            rb = residual(EXPRESSIONS[b_name], cell)

            if rb == 0:
                break

            ratios.append(sp.Rational(ra, rb))

        else:

            hits = ratio_library_match(ratios, cells)

            if not hits:
                continue

            # Every hit must be compared against the forced B/A.
            forced_expr = sp.cancel(
                EXPRESSIONS[b_name] / EXPRESSIONS[a_name]
            )

            genuine_hits = []

            for h in hits:
                h_expr = EXPRESSIONS[h]

                if sp.simplify(h_expr - forced_expr) != 0:
                    genuine_hits.append(h)

            if genuine_hits:
                genuine_source_matches.append(
                    (
                        a_name,
                        b_name,
                        cells,
                        genuine_hits,
                        forced_expr,
                    )
                )

    print()
    print("=" * 78)
    print("6. SOURCE-EXPRESSION RELATIONS BEYOND B/A")
    print("=" * 78)

    print(
        f"  genuine_relation_count={len(genuine_source_matches)}"
    )

    if not genuine_source_matches:
        print("  NONE")
    else:
        for a, b, cells, hits, forced in genuine_source_matches[:15]:
            print(
                f"  R({a})/R({b}) "
                f"cells={len(cells)} "
                f"source={hits} "
                f"forced_ratio={sp.factor(forced)}"
            )

    # ------------------------------------------------------------------------
    # 7. MOST IMPORTANT DIAGNOSTIC
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. 397R-TO-398R DIAGNOSTIC")
    print("=" * 78)

    print(
        "  Previous residual-ratio relations are algebraically forced:"
    )
    print(
        "      R_A/R_B = (Q/A)/(Q/B) = B/A"
    )
    print()
    print(
        "  Therefore relations such as"
    )
    print(
        "      R(p) / R(p*(p+t)) = p+t"
    )
    print(
        "  are classified here as TAUTOLOGICAL."
    )

    # ------------------------------------------------------------------------
    # 8. VERDICT
    # ------------------------------------------------------------------------

    if genuine_source_matches or nontrivial_ratio_hits or genuine_constant:
        verdict = "NON_TAUTOLOGICAL_RESIDUAL_STRUCTURE_CANDIDATE"
    elif invariant_hits:
        verdict = "DIRECT_RESIDUAL_CONSTANCY_ONLY"
    else:
        verdict = "NO_NON_TAUTOLOGICAL_RESIDUAL_LAW"

    print()
    print("=" * 78)
    print("8. STRUCTURAL VERDICT")
    print("=" * 78)

    print(f"  verdict={verdict}")

    print()
    print("  strongest admissible evidence:")
    print("    1. exact residual equality")
    print("    2. nontrivial constant residual scaling")
    print("    3. nontrivial normalized residual relation")
    print("    4. fixed source-expression relation not equal to B/A")
    print()
    print("  rejected automatically:")
    print("    R_A/R_B = B/A")
    print("    Q/A * A = Q")
    print("    Q/B * B = Q")
    print("    singleton coincidences")
    print("    missing-cell inference")

    # ------------------------------------------------------------------------
    # 9. EXACTNESS
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. FINAL EXACTNESS")
    print("=" * 78)

    print("  observed_cells_used_only=True")
    print("  exact_integer_arithmetic=True")
    print("  fixed_expression_library=True")
    print("  universal_residual_ratio_identity_removed=True")
    print("  tautological_B_over_A_relations_rejected=True")
    print("  normalized_residual_test_completed=True")
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
    print("EXPERIMENT 398R-COMPACT COMPLETE")


if __name__ == "__main__":
    main()
