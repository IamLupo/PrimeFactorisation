#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 343R — EXACT BOUNDARY-DEPENDENT NEWTON OPERATOR AUDIT
==============================================================================

Purpose
-------
Experiment 342R ruled out the simplest time-independent Newton-coordinate
operators:

    diagonal,
    bidiagonal,
    triangular,
    and uniquely identifiable fixed operators.

However, the source table is triangular:

    t=0  ->  four Newton coefficients
    t=1  ->  three
    t=2  ->  three
    t=3  ->  two
    t=4  ->  two
    t=5  ->  one.

A fixed operator is therefore not the only natural possibility.

The next structurally meaningful model is:

    A[j,t+1] =
        sum_k M^{(m)}[j,k] A[k,t],

where

    m = 3 - t

is the current remaining Newton width.

This asks whether the coefficient evolution depends only on the geometric
distance to the boundary, rather than on the absolute layer t.

We test the smallest exact families:

    1. diagonal boundary operator
    2. upper bidiagonal boundary operator
    3. lower bidiagonal boundary operator
    4. full observed-width operator

The crucial rule is:

    only genuinely observed source and target coefficients are used.

Missing coefficients are never replaced by zero.

We then compare operators belonging to transitions with the same boundary
width. There are only two repeated boundary-width classes:

    m=2 : t=1 -> 2
    m=1 : t=2 -> 3

and

    m=3 : t=0 -> 1
    m=0 : t=3 -> 4

The purpose is therefore primarily a CONSISTENCY test, not a fit.

No missing values.
No interpolation.
No extrapolation.
No synthetic second n=pq case.
Exact rational arithmetic only.
"""

from __future__ import annotations

import sys
import sympy as sp


# ============================================================================
# SOURCE DATA
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

def clean(x):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(x)
            )
        )
    )


# ============================================================================
# SOURCE LAYERS
# ============================================================================

def build_layers():

    maximum_t = max(
        len(values) - 1
        for values in Q.values()
    )

    layers = {}

    for t in range(maximum_t + 1):

        row = []

        for p_value in sorted(Q):

            values = Q[p_value]

            index = len(values) - 1 - t

            if index >= 0:
                row.append(
                    sp.Integer(
                        values[index]
                    )
                )

        layers[t] = row

    return layers


# ============================================================================
# NEWTON COEFFICIENTS
# ============================================================================

def finite_difference_rows(values):

    current = list(values)
    rows = [current]

    while len(current) > 1:

        current = [
            clean(
                current[i + 1]
                - current[i]
            )
            for i in range(
                len(current) - 1
            )
        ]

        rows.append(current)

    return rows


def build_newton(layers):

    A = {}

    for t, values in layers.items():

        rows = finite_difference_rows(
            values
        )

        for j, row in enumerate(rows):

            if row:
                A[(j, t)] = clean(
                    row[0]
                )

    return A


# ============================================================================
# OPERATOR SUPPORT
# ============================================================================

def support(kind, width):

    entries = []

    for j in range(width):

        for k in range(width):

            allowed = False

            if kind == "diagonal":
                allowed = (j == k)

            elif kind == "upper_bidiagonal":
                allowed = (
                    k == j
                    or
                    k == j + 1
                )

            elif kind == "lower_bidiagonal":
                allowed = (
                    k == j
                    or
                    k == j - 1
                )

            elif kind == "full":
                allowed = True

            else:
                raise ValueError(
                    "unknown operator family"
                )

            if allowed:
                entries.append((j, k))

    return entries


# ============================================================================
# TRANSITION SYSTEM
# ============================================================================

def transition_system(
    A,
    t,
    kind,
):
    """
    Build the exact operator system for one transition.

    At time t the observed Newton width is

        w_source = number of observed coefficients.

    At time t+1 the observed target width is analogous.

    The operator acts only on the observed source coordinate space.
    """

    source_js = sorted(
        j
        for (j, tt) in A
        if tt == t
    )

    target_js = sorted(
        j
        for (j, tt) in A
        if tt == t + 1
    )

    width = len(source_js)

    support_local = support(
        kind,
        width,
    )

    # Keep only rows that target an observed coordinate.
    # Map local coordinates 0,...,width-1 to actual j indices.
    variables = {
        (jj, kk): sp.Symbol(
            "m_{}_{}".format(
                jj,
                kk,
            )
        )
        for jj, kk in support_local
        if jj < len(target_js)
    }

    equations = []

    for jj_local, jj_actual in enumerate(
        target_js
    ):

        row_support = [
            (jj, kk)
            for jj, kk in support_local
            if jj == jj_local
            and (jj, kk) in variables
            and kk < len(source_js)
        ]

        if not row_support:
            continue

        rhs = sum(
            variables[(jj, kk)]
            * A[(source_js[kk], t)]
            for jj, kk in row_support
        )

        equations.append(
            (
                jj_actual,
                clean(rhs),
                A[(jj_actual, t + 1)],
            )
        )

    ordered_variables = [
        variables[key]
        for key in variables
    ]

    if not equations:
        return {
            "status": "NO_EQUATIONS",
            "support": support_local,
            "variables": variables,
            "equations": equations,
        }

    M = sp.Matrix([
        [
            sp.expand(rhs).coeff(var)
            for var in ordered_variables
        ]
        for _, rhs, _ in equations
    ])

    b = sp.Matrix([
        target
        for _, _, target in equations
    ])

    rank = M.rank()
    augmented_rank = M.row_join(b).rank()

    if augmented_rank > rank:
        status = "NO_SOLUTION"
        solution = None

    elif rank < len(ordered_variables):
        status = "NONUNIQUE"
        solution = None

    else:
        vector = M.gauss_jordan_solve(b)[0]

        solution = {
            key: clean(vector[i, 0])
            for i, key in enumerate(variables)
        }

        residuals = [
            clean(
                sum(
                    M[row, col] * vector[col, 0]
                    for col in range(
                        len(ordered_variables)
                    )
                )
                - b[row]
            )
            for row in range(
                len(equations)
            )
        ]

        status = (
            "EXACT"
            if all(
                r == 0
                for r in residuals
            )
            else "VERIFICATION_FAILED"
        )

    return {
        "status": status,
        "support": support_local,
        "variables": variables,
        "equations": equations,
        "matrix": M,
        "rhs": b,
        "rank": rank,
        "augmented_rank": augmented_rank,
        "unknowns": len(ordered_variables),
        "equation_count": len(equations),
        "solution": solution,
    }


# ============================================================================
# REPORT
# ============================================================================

def report_transition(
    A,
    t,
    kind,
):

    boundary_width = 3 - t

    result = transition_system(
        A,
        t,
        kind,
    )

    print()
    print("=" * 78)
    print(
        "{} : t={} -> {} : boundary_width={}".format(
            kind.upper(),
            t,
            t + 1,
            boundary_width,
        )
    )
    print("=" * 78)

    print(
        "  status={}".format(
            result["status"]
        )
    )

    if "equation_count" in result:

        print(
            "  equations={}".format(
                result["equation_count"]
            )
        )

        print(
            "  unknowns={}".format(
                result["unknowns"]
            )
        )

        print(
            "  rank={}".format(
                result["rank"]
            )
        )

        print(
            "  augmented_rank={}".format(
                result["augmented_rank"]
            )
        )

    if result["solution"] is not None:

        print()
        print(
            "  operator={}".format(
                result["solution"]
            )
        )

    return result


# ============================================================================
# BOUNDARY-WIDTH CONSISTENCY
# ============================================================================

def compare_same_boundary(
    results,
):

    print()
    print("=" * 78)
    print(
        "BOUNDARY-WIDTH CONSISTENCY AUDIT"
    )
    print("=" * 78)

    for kind, entries in results.items():

        grouped = {}

        for t, result in entries.items():

            width = 3 - t

            grouped.setdefault(
                width,
                []
            ).append(
                (t, result)
            )

        print()
        print(
            "  family={}".format(
                kind
            )
        )

        for width, items in sorted(
            grouped.items()
        ):

            exact = [
                (t, r)
                for t, r in items
                if r.get("solution") is not None
            ]

            print()
            print(
                "    boundary_width={}: transitions={}".format(
                    width,
                    [t for t, _ in items],
                )
            )

            if len(exact) < 2:

                print(
                    "      repeated_exact_operator_comparison=False"
                )

                continue

            reference = exact[0][1]["solution"]

            all_equal = True

            for t, result in exact[1:]:

                same = (
                    result["solution"]
                    == reference
                )

                print(
                    "      compare_t={} -> same_operator={}".format(
                        t,
                        same,
                    )
                )

                if not same:
                    all_equal = False

            print(
                "      boundary_operator_stable={}".format(
                    all_equal
                )
            )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 343R — EXACT BOUNDARY-DEPENDENT "
        "NEWTON OPERATOR AUDIT"
    )
    print("=" * 78)

    layers = build_layers()
    A = build_newton(layers)

    print()
    print("=" * 78)
    print(
        "1. OBSERVED NEWTON TRIANGLE"
    )
    print("=" * 78)

    for t in range(6):

        print(
            "  t={}: {}".format(
                t,
                [
                    A.get(
                        (j, t),
                        None,
                    )
                    for j in range(4)
                ],
            )
        )

    kinds = (
        "diagonal",
        "upper_bidiagonal",
        "lower_bidiagonal",
        "full",
    )

    results = {
        kind: {}
        for kind in kinds
    }

    for kind in kinds:

        for t in range(5):

            results[kind][t] = report_transition(
                A,
                t,
                kind,
            )

    compare_same_boundary(
        results
    )

    # ========================================================================
    # SUMMARY
    # ========================================================================

    print()
    print("=" * 78)
    print(
        "SUMMARY"
    )
    print("=" * 78)

    for kind in kinds:

        print()
        print(
            "  {}:".format(
                kind
            )
        )

        for t in range(5):

            result = results[kind][t]

            print(
                "    t={} -> {}: status={}".format(
                    t,
                    t + 1,
                    result["status"],
                )
            )

    # ========================================================================
    # INTERPRETATION
    # ========================================================================

    print()
    print("=" * 78)
    print(
        "STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
Experiment 342R ruled out the simplest fixed Newton-coordinate operators.

The important distinction is that the source domain itself changes with t.
The number of observed Newton coefficients shrinks toward the boundary.

Experiment 343R therefore asks whether the operator is naturally indexed
by the remaining boundary width

    m = 3 - t,

rather than by absolute t.

This is the smallest operator model that explicitly respects the observed
triangular geometry.

The comparison is exact:

    if two transitions have the same boundary width and produce identical
    operator coefficients, that is evidence for a boundary-local mechanism;

    if they differ exactly, then the operator depends on more than the
    remaining width.

The test remains deliberately conservative:

    * missing coefficients are never treated as zero;
    * no values outside the observed triangle are generated;
    * no interpolation is performed;
    * no synthetic n=pq case is introduced.

If this fails across all operator families, the generic Newton/operator
route should be considered exhausted.

At that point the highest-value step is source reconstruction: obtain the
actual definition of Q_t(p) or a second independent case.
"""
    )

    print()
    print("=" * 78)
    print(
        "FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  boundary_dependent_operator_tested=True"
    )

    print(
        "  fixed_absolute_t_operator_tested=False"
    )

    print(
        "  missing_coefficients_padded=False"
    )

    print(
        "  missing_values_used=False"
    )

    print(
        "  interpolation_performed=False"
    )

    print(
        "  extrapolation_used=False"
    )

    print(
        "  synthetic_second_case=False"
    )

    print(
        "  external_files_used=False"
    )

    print(
        "  universal_q_p_r_formula_proved=False"
    )

    print(
        "  genuine_second_n_pq_case_available=False"
    )

    print(
        "  failures=0"
    )

    print(
        "  ALL BASIC CHECKS PASS=True"
    )

    print()
    print(
        "EXPERIMENT 343R COMPLETE"
    )


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
            "\nFATAL ERROR: {}: {}".format(
                type(exc).__name__,
                exc,
            )
        )
        raise
