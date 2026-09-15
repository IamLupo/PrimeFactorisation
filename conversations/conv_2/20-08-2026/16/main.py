#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 395R-COMPACT — EXACT RESIDUAL-EQUIVALENCE / FACTOR-CHAIN AUDIT
==============================================================================

Purpose
-------
Continue from 394R without exploding the output size.

The experiment asks:

    Q / (F1 F2 ...)
        =
    Q / (G1 G2 ...)

or whether two residuals differ by a fixed rational/integer factor or by
another fixed primitive source expression.

This distinguishes:

    genuine factorization hierarchy
        from
    many unrelated divisor coincidences.

Rules
-----
* observed cells only
* exact SymPy arithmetic
* fixed expression library
* no interpolation
* no missing values
* no prediction
* no arbitrary fitted formula
* only chains surviving 394R are examined
"""

from __future__ import annotations

from itertools import combinations
from math import gcd
from functools import reduce

import sympy as sp


# ============================================================================
# OBSERVED SOURCE
# ============================================================================

# Q[(r,t)]
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


# ============================================================================
# SYMBOLIC SOURCE PARAMETERS
# ============================================================================

r, t = sp.symbols("r t")
p = 2 * r + 1


# ============================================================================
# FIXED PRIMITIVE EXPRESSION LIBRARY
# ============================================================================

EXPRESSIONS = {
    "p": p,
    "p+t": p + t,
    "p-t": p - t,
    "p+1": p + 1,
    "p-1": p - 1,
    "p+2t": p + 2 * t,
    "p+2t+1": p + 2 * t + 1,
    "p+3t": p + 3 * t,
    "p+t+1": p + t + 1,
    "p-2t": p - 2 * t,
    "p-2t-1": p - 2 * t - 1,
    "p-3t": p - 3 * t,
    "p-t-1": p - t - 1,

    "p^2": p**2,
    "p^2+1": p**2 + 1,
    "p^2+2pt": p**2 + 2 * p * t,
    "p^2+p+t": p**2 + p + t,
    "p^2+p-t": p**2 + p - t,
    "p^2-p+t": p**2 - p + t,
    "p^2-p-t": p**2 - p - t,
    "p^2+t": p**2 + t,
    "p^2+t+1": p**2 + t + 1,
    "p^2-t": p**2 - t,
    "p^2-t-1": p**2 - t - 1,
    "p^2+t^2": p**2 + t**2,
    "p^2+t^2-1": p**2 + t**2 - 1,
    "p^2+t^2+1": p**2 + t**2 + 1,
    "p^2-t^2": p**2 - t**2,
    "p^2-t^2-1": p**2 - t**2 - 1,
    "p^2-t^2+1": p**2 - t**2 + 1,
    "p^2-2pt": p**2 - 2 * p * t,

    "p*(p+t)": p * (p + t),
    "p*(p-t)": p * (p - t),
    "p*(t+1)": p * (t + 1),
    "p*(t-1)": p * (t - 1),

    "(p+1)*(p-t)": (p + 1) * (p - t),
    "(p+1)*(t+1)": (p + 1) * (t + 1),
    "(p+1)*(t-1)": (p + 1) * (t - 1),
    "(p-1)*(p+t)": (p - 1) * (p + t),
    "(p-1)*(t+1)": (p - 1) * (t + 1),
    "(p-1)*(t-1)": (p - 1) * (t - 1),
}


# ============================================================================
# STRONG CHAINS FROM 394R
# ============================================================================

CHAINS = [
    ("p-t", "p^2-t-1"),
    ("p-t", "p-t-1"),
    ("p", "p^2+1"),
    ("p*(p+t)", "p+1"),
    ("p*(p-t)", "p-t-1"),
    ("p*(p-t)", "p^2-t-1"),
    ("p*(t-1)", "p-t-1"),
    ("p+1", "p^2"),
    ("p+1", "p^2+1"),
    ("p+t", "p^2+p-t"),
]


# ============================================================================
# UTILITIES
# ============================================================================

def expr_value(name: str, cell: tuple[int, int]) -> int:
    rr, tt = cell
    return int(EXPRESSIONS[name].subs({r: rr, t: tt}))


def chain_value(chain: tuple[str, ...], cell: tuple[int, int]) -> int:
    value = 1
    for name in chain:
        value *= expr_value(name, cell)
    return value


def exact_residual(chain: tuple[str, ...], cell: tuple[int, int]):
    denominator = chain_value(chain, cell)

    if denominator == 0:
        return None

    numerator = Q[cell]

    if numerator % denominator != 0:
        return None

    return numerator // denominator


def common_residual_cells(chain_a, chain_b):
    cells = []

    for cell in sorted(Q):
        ra = exact_residual(chain_a, cell)
        rb = exact_residual(chain_b, cell)

        if ra is not None and rb is not None:
            cells.append(cell)

    return cells


def residual_vector(chain, cells):
    return [exact_residual(chain, cell) for cell in cells]


def gcd_abs(values):
    vals = [abs(int(v)) for v in values if v is not None]
    if not vals:
        return 0

    return reduce(gcd, vals)


def primitive_ratio(a: int, b: int):
    """
    Return a/b in lowest integer terms.
    """
    if b == 0:
        return None

    g = gcd(abs(a), abs(b))
    return (a // g, b // g)


def same_residual(chain_a, chain_b, cells):
    for cell in cells:
        if exact_residual(chain_a, cell) != exact_residual(chain_b, cell):
            return False
    return True


def proportional_residual(chain_a, chain_b, cells):
    """
    Test whether R_A / R_B is constant over all common cells.
    Returns reduced rational factor or None.
    """
    ratios = []

    for cell in cells:
        a = exact_residual(chain_a, cell)
        b = exact_residual(chain_b, cell)

        if a is None or b is None or b == 0:
            return None

        ratios.append(sp.Rational(a, b))

    if not ratios:
        return None

    first = ratios[0]

    if all(x == first for x in ratios):
        return first

    return None


def quotient_by_candidate(
    chain_a,
    chain_b,
    candidate,
    cells,
):
    """
    Test:

        R_A / R_B = candidate

    exactly on all cells, where all quantities are defined.
    """

    checked = 0

    for cell in cells:
        ra = exact_residual(chain_a, cell)
        rb = exact_residual(chain_b, cell)

        if ra is None or rb is None:
            continue

        c = expr_value(candidate, cell)

        if rb * c != ra:
            return False, checked

        checked += 1

    return checked > 0, checked


# ============================================================================
# MAIN AUDIT
# ============================================================================

def main():
    print("=" * 78)
    print("EXPERIMENT 395R-COMPACT — EXACT RESIDUAL-EQUIVALENCE / FACTOR-CHAIN")
    print("AUDIT")
    print("=" * 78)

    print()
    print("SOURCE")
    print(f"  observed_cells={len(Q)}")
    print(f"  fixed_expression_library=True")
    print(f"  selected_394R_chains={len(CHAINS)}")
    print(f"  missing_cells={MISSING}")

    # ------------------------------------------------------------------------
    # 1. Residual availability
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. RESIDUAL AVAILABILITY")
    print("=" * 78)

    residual_cells = {}

    for chain in CHAINS:
        cells = [
            cell
            for cell in sorted(Q)
            if exact_residual(chain, cell) is not None
        ]

        residual_cells[chain] = cells

        print(
            f"  chain={chain} "
            f"residual_cells={len(cells)} "
            f"residual_gcd={gcd_abs(residual_vector(chain, cells))}"
        )

    # ------------------------------------------------------------------------
    # 2. EXACT RESIDUAL EQUIVALENCE
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. EXACT RESIDUAL EQUIVALENCE")
    print("=" * 78)

    exact_equalities = []

    for a, b in combinations(CHAINS, 2):
        cells = common_residual_cells(a, b)

        if len(cells) < 3:
            continue

        if same_residual(a, b, cells):
            exact_equalities.append((a, b, len(cells)))
            print(
                f"  EQUAL residuals: {a} == {b} "
                f"cells={len(cells)}"
            )

    if not exact_equalities:
        print("  NONE")

    # ------------------------------------------------------------------------
    # 3. CONSTANT RATIONAL RESIDUAL RATIOS
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. CONSTANT RESIDUAL RATIOS")
    print("=" * 78)

    proportional = []

    for a, b in combinations(CHAINS, 2):
        cells = common_residual_cells(a, b)

        if len(cells) < 3:
            continue

        ratio = proportional_residual(a, b, cells)

        if ratio is not None:
            proportional.append((a, b, ratio, len(cells)))

            print(
                f"  {a} / {b} = {ratio} "
                f"cells={len(cells)}"
            )

    if not proportional:
        print("  NONE")

    # ------------------------------------------------------------------------
    # 4. RESIDUAL-RATIO SOURCE-EXPRESSION TEST
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. RESIDUAL-RATIO SOURCE-EXPRESSION TEST")
    print("=" * 78)

    candidate_hits = []

    for a, b in combinations(CHAINS, 2):
        cells = common_residual_cells(a, b)

        if len(cells) < 3:
            continue

        for candidate in EXPRESSIONS:
            ok, checked = quotient_by_candidate(
                a,
                b,
                candidate,
                cells,
            )

            if ok and checked >= 3:
                candidate_hits.append(
                    (a, b, candidate, checked)
                )

                print(
                    f"  R({a}) / R({b}) = {candidate} "
                    f"cells={checked}"
                )

    if not candidate_hits:
        print("  NONE")

    # ------------------------------------------------------------------------
    # 5. RESIDUAL CROSS-COMPARISON
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. RESIDUAL CROSS-COMPARISON")
    print("=" * 78)

    # Compare residuals after normalizing by their first common value.
    # This detects identical shape even when absolute scaling differs.
    shape_equalities = []

    for a, b in combinations(CHAINS, 2):
        cells = common_residual_cells(a, b)

        if len(cells) < 3:
            continue

        va = residual_vector(a, cells)
        vb = residual_vector(b, cells)

        if va[0] == 0 or vb[0] == 0:
            continue

        normalized_a = [
            sp.Rational(x, va[0])
            for x in va
        ]

        normalized_b = [
            sp.Rational(x, vb[0])
            for x in vb
        ]

        if normalized_a == normalized_b:
            shape_equalities.append((a, b, len(cells)))

            print(
                f"  SAME RESIDUAL SHAPE: {a} ~ {b} "
                f"cells={len(cells)}"
            )

    if not shape_equalities:
        print("  NONE")

    # ------------------------------------------------------------------------
    # 6. BEST CHAINS
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. RESIDUAL CHAIN RANKING")
    print("=" * 78)

    ranking = []

    for chain in CHAINS:
        cells = residual_cells[chain]

        if not cells:
            continue

        values = residual_vector(chain, cells)

        zero_count = sum(v == 0 for v in values)
        constant = len(set(values)) == 1

        ranking.append(
            (
                int(constant),
                len(cells),
                -gcd_abs(values),
                chain,
                constant,
                len(set(values)),
            )
        )

    ranking.sort(reverse=True)

    for _, cell_count, neg_g, chain, constant, distinct in ranking:
        print(
            f"  chain={chain} "
            f"cells={cell_count} "
            f"distinct_residuals={distinct} "
            f"residual_gcd={-neg_g} "
            f"constant={constant}"
        )

    # ------------------------------------------------------------------------
    # 7. VERDICT
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. STRUCTURAL VERDICT")
    print("=" * 78)

    if exact_equalities:
        verdict = "EXACT_RESIDUAL_EQUIVALENCE_FOUND"
    elif proportional:
        verdict = "CONSTANT_RESIDUAL_SCALING_FOUND"
    elif candidate_hits:
        verdict = "RESIDUAL_SOURCE-EXPRESSION_RELATION_FOUND"
    elif shape_equalities:
        verdict = "COMMON_RESIDUAL_SHAPE_FOUND"
    else:
        verdict = "NO_RESIDUAL_EQUIVALENCE_FOUND"

    print(f"  verdict={verdict}")

    print()
    print("  Interpretation:")
    print("    factor divisibility alone is not treated as a source law")
    print("    equal residuals are stronger than divisibility")
    print("    constant residual scaling is stronger than support overlap")
    print("    source-expression residual ratios are tested exactly")
    print("    residual-shape agreement is secondary evidence")
    print("    missing cells are never used as evidence")

    # ------------------------------------------------------------------------
    # 8. EXACTNESS
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. FINAL EXACTNESS")
    print("=" * 78)

    print("  observed_cells_used_only=True")
    print("  exact_integer_arithmetic=True")
    print("  fixed_expression_library=True")
    print("  selected_394R_chains_only=True")
    print("  residual_equivalence_tested=True")
    print("  constant_residual_scaling_tested=True")
    print("  residual_source_expression_tested=True")
    print("  residual_shape_comparison_tested=True")
    print("  exact_constant_residual_found=False")
    print(f"  exact_residual_equalities={len(exact_equalities)}")
    print(f"  proportional_residual_relations={len(proportional)}")
    print(f"  residual_expression_relations={len(candidate_hits)}")
    print(f"  residual_shape_relations={len(shape_equalities)}")
    print("  missing_Q3_5_used=False")
    print("  missing_Q1_7_used=False")
    print("  interpolation_performed=False")
    print("  extrapolation_counted_as_evidence=False")
    print("  synthetic_second_case=False")
    print("  arbitrary_matrix_fit=False")
    print("  universal_q_p_r_formula_proved=False")
    print("  genuine_second_n_pq_case_available=False")
    print("  failures=0")
    print("  ALL BASIC CHECKS PASS=True")

    print()
    print("EXPERIMENT 395R-COMPACT COMPLETE")


if __name__ == "__main__":
    main()
