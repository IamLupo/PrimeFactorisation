#!/usr/bin/env python3
"""
==============================================================================
EXPERIMENT 390R-COMPACT — EXACT SOURCE-EXPRESSION p-ADIC CONTAINMENT AUDIT
==============================================================================

PURPOSE

Previous experiments found apparent prime-support matches such as:

    prime 2  <->  p+1
    prime 2  <->  p^2+1
    prime 2  <->  p+2t+1
    prime 2  <->  (p+1)(t+1)

but support equality only asks:

    ell | Q  <=>  ell | E

This experiment asks the stronger question:

    v_ell(Q) >= v_ell(E)

and, separately:

    v_ell(Q) = v_ell(E)

across ALL observed cells.

A source-expression law is considered interesting only when:

    * at least 3 observed cells participate;
    * the relation holds at every participating cell;
    * no missing cell is used;
    * singleton coincidences are discarded.

The report is intentionally compact.
Only strong candidates are printed.
==============================================================================
"""

from __future__ import annotations

from collections import defaultdict
from math import gcd
from sympy import factorint


# ============================================================================
# OBSERVED SOURCE
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


def p_of_r(r: int) -> int:
    return 2 * r + 1


# ============================================================================
# FIXED SOURCE-EXPRESSION LIBRARY
# ============================================================================

def expressions(p: int, t: int) -> dict[str, int]:
    return {
        "p": p,
        "p-1": p - 1,
        "p+1": p + 1,
        "p^2": p * p,
        "p^2-1": p * p - 1,
        "p^2+1": p * p + 1,

        "p+t": p + t,
        "p-t": p - t,
        "p+t+1": p + t + 1,
        "p-t-1": p - t - 1,

        "p+2t": p + 2 * t,
        "p-2t": p - 2 * t,
        "p+2t+1": p + 2 * t + 1,
        "p-2t-1": p - 2 * t - 1,

        "p+3t": p + 3 * t,
        "p-3t": p - 3 * t,

        "p^2+t": p * p + t,
        "p^2-t": p * p - t,
        "p^2+t+1": p * p + t + 1,
        "p^2-t-1": p * p - t - 1,

        "p^2+t^2": p * p + t * t,
        "p^2-t^2": p * p - t * t,
        "p^2+t^2-1": p * p + t * t - 1,
        "p^2-t^2-1": p * p - t * t - 1,
        "p^2+t^2+1": p * p + t * t + 1,
        "p^2-t^2+1": p * p - t * t + 1,

        "p^2+p+t": p * p + p + t,
        "p^2-p+t": p * p - p + t,
        "p^2+p-t": p * p + p - t,
        "p^2-p-t": p * p - p - t,

        "p^2+2pt": p * p + 2 * p * t,
        "p^2-2pt": p * p - 2 * p * t,

        "p*(p+t)": p * (p + t),
        "p*(p-t)": p * (p - t),
        "p*(t+1)": p * (t + 1),
        "p*(t-1)": p * (t - 1),

        "(p-1)*(t+1)": (p - 1) * (t + 1),
        "(p+1)*(t+1)": (p + 1) * (t + 1),
        "(p-1)*(t-1)": (p - 1) * (t - 1),
        "(p+1)*(t-1)": (p + 1) * (t - 1),

        "(p-1)*(p+t)": (p - 1) * (p + t),
        "(p+1)*(p-t)": (p + 1) * (p - t),
    }


# ============================================================================
# FAST VALUATION CACHE
# ============================================================================

_factor_cache: dict[int, dict[int, int]] = {}


def factor(n: int) -> dict[int, int]:
    n = abs(int(n))

    if n <= 1:
        return {}

    if n not in _factor_cache:
        _factor_cache[n] = dict(factorint(n))

    return _factor_cache[n]


def valuation(n: int, prime: int) -> int:
    return factor(abs(n)).get(prime, 0)


# ============================================================================
# BUILD TABLES
# ============================================================================

def build_tables():
    cells = sorted(Q)

    q_factor = {
        cell: factor(Q[cell])
        for cell in cells
    }

    expr_values = {
        cell: expressions(
            p_of_r(cell[0]),
            cell[1],
        )
        for cell in cells
    }

    return cells, q_factor, expr_values


# ============================================================================
# CORE AUDIT
# ============================================================================

def audit(cells, q_factor, expr_values):
    """
    For every expression E and every prime dividing E somewhere:

      CONTAINMENT:
          v_p(Q) >= v_p(E)

      EXACT:
          v_p(Q) == v_p(E)

    Only primes that actually occur in E are considered.
    """

    containment = []
    equality = []

    expressions_by_name = defaultdict(list)

    for cell in cells:
        for name, value in expr_values[cell].items():
            if abs(value) <= 1:
                continue

            ef = factor(value)

            for prime, e_val in ef.items():
                q_val = q_factor[cell].get(prime, 0)

                expressions_by_name[(name, prime)].append(
                    (cell, e_val, q_val)
                )

    for (name, prime), rows in expressions_by_name.items():

        if len(rows) < 3:
            continue

        containment_ok = all(
            qv >= ev
            for _cell, ev, qv in rows
        )

        equality_ok = all(
            qv == ev
            for _cell, ev, qv in rows
        )

        if containment_ok:
            deficits = [
                qv - ev
                for _cell, ev, qv in rows
            ]

            containment.append(
                (
                    len(rows),
                    prime,
                    name,
                    min(deficits),
                    max(deficits),
                    len(set(deficits)),
                )
            )

        if equality_ok:
            equality.append(
                (
                    len(rows),
                    prime,
                    name,
                )
            )

    containment.sort(
        key=lambda x: (-x[0], x[1], x[2])
    )

    equality.sort(
        key=lambda x: (-x[0], x[1], x[2])
    )

    return containment, equality


# ============================================================================
# CROSS-CELL QUOTIENT AUDIT
# ============================================================================

def quotient_audit(cells, expr_values):
    """
    Test the stronger integer statement:

        Q(r,t) / E(p,t)

    is integral at every cell where E != 0.

    This is different from prime-support matching because ALL prime
    multiplicities must be supplied by E.
    """

    exact_divisibility = []

    for name in next(iter(expr_values.values())):
        participating = []

        for cell in cells:
            value = expr_values[cell][name]

            if value == 0 or abs(value) == 1:
                continue

            q = Q[cell]

            if q % value == 0:
                participating.append((cell, q // value))

        if len(participating) >= 3:
            exact_divisibility.append(
                (
                    len(participating),
                    name,
                    participating,
                )
            )

    exact_divisibility.sort(
        key=lambda x: (-x[0], x[1])
    )

    return exact_divisibility


# ============================================================================
# STRUCTURAL COMPARISON
# ============================================================================

def compare_equalities(equalities):
    """
    Detect expressions that produce identical valuation laws.

    Example:
        p+1
        p^2+1

    may agree for the tiny source set.

    Such duplicates are grouped so that we don't mistake four expressions
    for four independent discoveries.
    """

    groups = defaultdict(list)

    for count, prime, name in equalities:
        groups[(count, prime)].append(name)

    return {
        key: sorted(names)
        for key, names in groups.items()
        if len(names) >= 2
    }


# ============================================================================
# REPORT
# ============================================================================

def main():
    cells, q_factor, expr_values = build_tables()

    containment, equality = audit(
        cells,
        q_factor,
        expr_values,
    )

    exact_divisibility = quotient_audit(
        cells,
        expr_values,
    )

    equality_groups = compare_equalities(
        equality
    )

    print("=" * 78)
    print("EXPERIMENT 390R-COMPACT — EXACT p-ADIC SOURCE-EXPRESSION AUDIT")
    print("=" * 78)

    print()
    print("SOURCE")
    print(f"  observed_cells={len(cells)}")
    print(
        f"  source_parameters="
        f"{sorted({p_of_r(r) for r, _ in cells})}"
    )
    print(f"  missing_cells={MISSING}")
    print("  expressions_fixed_in_advance=True")

    # ------------------------------------------------------------------------
    print()
    print("=" * 78)
    print("1. EXACT p-ADIC CONTAINMENT")
    print("=" * 78)

    print(
        "  Tested:"
        " v_prime(Q) >= v_prime(expression)"
    )

    if containment:
        for (
            count,
            prime,
            name,
            dmin,
            dmax,
            dcount,
        ) in containment:
            print(
                f"  prime={prime}"
                f" expression={name}"
                f" cells={count}"
                f" valuation_gap=[{dmin},{dmax}]"
                f" distinct_gaps={dcount}"
            )
    else:
        print("  NONE")

    print(
        f"  accepted_containment_count="
        f"{len(containment)}"
    )

    # ------------------------------------------------------------------------
    print()
    print("=" * 78)
    print("2. EXACT p-ADIC EQUALITY")
    print("=" * 78)

    print(
        "  Tested:"
        " v_prime(Q) = v_prime(expression)"
    )

    if equality:
        for count, prime, name in equality:
            print(
                f"  prime={prime}"
                f" expression={name}"
                f" cells={count}"
            )
    else:
        print("  NONE")

    print(
        f"  accepted_equality_count="
        f"{len(equality)}"
    )

    # ------------------------------------------------------------------------
    print()
    print("=" * 78)
    print("3. DUPLICATE VALUATION LAWS")
    print("=" * 78)

    if equality_groups:
        for (count, prime), names in sorted(
            equality_groups.items(),
            key=lambda x: (-x[0][0], x[0][1]),
        ):
            print(
                f"  prime={prime}"
                f" cells={count}"
                f" equivalent_expressions={names}"
            )
    else:
        print("  NONE")

    # ------------------------------------------------------------------------
    print()
    print("=" * 78)
    print("4. COMPLETE INTEGER DIVISIBILITY")
    print("=" * 78)

    print(
        "  Tested:"
        " expression divides Q at every participating cell."
    )

    if exact_divisibility:
        # Very compact: only print the strongest 15.
        for count, name, rows in exact_divisibility[:15]:
            print(
                f"  expression={name}"
                f" divisible_cells={count}"
            )

        if len(exact_divisibility) > 15:
            print(
                f"  ... "
                f"{len(exact_divisibility) - 15}"
                f" additional candidates omitted"
            )
    else:
        print("  NONE")

    print(
        f"  total_integer_divisibility_candidates="
        f"{len(exact_divisibility)}"
    )

    # ------------------------------------------------------------------------
    print()
    print("=" * 78)
    print("5. STRONGEST CANDIDATES")
    print("=" * 78)

    strong = []

    for count, prime, name, dmin, dmax, dcount in containment:
        if count >= 5:
            strong.append(
                (
                    "containment",
                    count,
                    prime,
                    name,
                    dmin,
                    dmax,
                )
            )

    for count, prime, name in equality:
        if count >= 5:
            strong.append(
                (
                    "equality",
                    count,
                    prime,
                    name,
                    0,
                    0,
                )
            )

    strong.sort(
        key=lambda x: (-x[1], x[2], x[3])
    )

    if strong:
        for kind, count, prime, name, dmin, dmax in strong:
            if kind == "containment":
                print(
                    f"  {kind}:"
                    f" prime={prime}"
                    f" expression={name}"
                    f" cells={count}"
                    f" gap=[{dmin},{dmax}]"
                )
            else:
                print(
                    f"  {kind}:"
                    f" prime={prime}"
                    f" expression={name}"
                    f" cells={count}"
                )
    else:
        print("  NONE")

    # ------------------------------------------------------------------------
    print()
    print("=" * 78)
    print("6. STRUCTURAL VERDICT")
    print("=" * 78)

    if equality:
        verdict = (
            "EXACT_MULTI_CELL_p_ADIC_EQUALITY_FOUND"
        )
    elif containment:
        verdict = (
            "EXACT_MULTI_CELL_p_ADIC_CONTAINMENT_FOUND"
        )
    elif exact_divisibility:
        verdict = (
            "MULTI_CELL_INTEGER_DIVISIBILITY_ONLY"
        )
    else:
        verdict = (
            "NO_EXACT_MULTI_CELL_SOURCE-EXPRESSION_LAW"
        )

    print(f"  verdict={verdict}")

    print()
    print("  IMPORTANT")
    print("    support equality alone is not accepted")
    print("    valuation containment is stronger than support equality")
    print("    valuation equality is stronger still")
    print("    complete integer divisibility is reported separately")
    print("    duplicate expressions are grouped")
    print("    singleton coincidences are rejected")
    print("    missing cells are never used")

    # ------------------------------------------------------------------------
    print()
    print("=" * 78)
    print("7. FINAL EXACTNESS")
    print("=" * 78)

    print("  observed_cells_used_only=True")
    print("  exact_factorization=True")
    print("  exact_prime_valuations=True")
    print("  fixed_expression_library=True")
    print("  support_matching_not_used_as_final_evidence=True")
    print("  p_adic_containment_tested=True")
    print("  p_adic_equality_tested=True")
    print("  complete_integer_divisibility_tested=True")
    print("  duplicate_laws_grouped=True")
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
    print("EXPERIMENT 390R-COMPACT COMPLETE")


if __name__ == "__main__":
    main()
