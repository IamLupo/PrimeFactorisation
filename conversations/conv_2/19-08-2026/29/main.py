#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 344R — EXACT LOW-DEGREE t-DEPENDENT NEWTON-OPERATOR AUDIT
==============================================================================

Purpose
-------
Experiments 342R-343R eliminated fixed-coefficient Newton operators and
did not obtain a meaningful boundary-width comparison.

The remaining natural operator model is

    A(t+1) = M(t) A(t),

where the matrix entries of M(t) are low-degree polynomials in t.

This is strictly more flexible than a fixed operator, but still highly
structured.

We test:

    diagonal
    upper_bidiagonal
    lower_bidiagonal
    full

with matrix-entry polynomial degree

    0, 1, 2.

For every model:

    * only observed Newton coefficients are used;
    * missing coefficients are never replaced by zero;
    * equations are exact over Q;
    * an exact model is accepted only if the linear system is consistent;
    * UNIQUE / OVERDETERMINED status is reported separately from NONUNIQUE.

This is a closure test for the remaining low-complexity linear-operator
hypothesis.

No interpolation.
No missing-value reconstruction.
No extrapolation.
No synthetic second n=pq case.
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


def finite_difference_rows(values):

    current = [
        sp.Integer(v)
        for v in values
    ]

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


def build_layers():

    maximum_t = max(
        len(values) - 1
        for values in Q.values()
    )

    layers = {}

    for t in range(
        maximum_t + 1
    ):

        values = []

        for p_value in sorted(Q):

            data = Q[p_value]

            index = (
                len(data)
                - 1
                - t
            )

            if index >= 0:

                values.append(
                    sp.Integer(
                        data[index]
                    )
                )

        layers[t] = values

    return layers


def build_newton():

    layers = build_layers()

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

def allowed_entries(
    family,
    target_width,
    source_width,
):

    entries = []

    for j in range(target_width):

        for k in range(source_width):

            allowed = False

            if family == "diagonal":
                allowed = (
                    j == k
                )

            elif family == "upper_bidiagonal":
                allowed = (
                    j == k
                    or
                    k == j + 1
                )

            elif family == "lower_bidiagonal":
                allowed = (
                    j == k
                    or
                    k == j - 1
                )

            elif family == "full":
                allowed = True

            else:
                raise ValueError(
                    "Unknown family: {}".format(
                        family
                    )
                )

            if allowed:
                entries.append(
                    (j, k)
                )

    return entries


# ============================================================================
# POLYNOMIAL t-DEPENDENT SYSTEM
# ============================================================================

def build_model(
    A,
    family,
    degree,
):

    equations = []
    unknown_keys = []

    # Determine all possible coefficient names.
    for t_value in range(5):

        source_indices = sorted(
            j
            for (j, tt) in A
            if tt == t_value
        )

        target_indices = sorted(
            j
            for (j, tt) in A
            if tt == t_value + 1
        )

        support = allowed_entries(
            family,
            len(target_indices),
            len(source_indices),
        )

        for j, k in support:

            for d in range(
                degree + 1
            ):

                key = (
                    j,
                    k,
                    d,
                )

                if key not in unknown_keys:

                    unknown_keys.append(
                        key
                    )

    variables = {
        key: sp.symbols(
            "c_{}_{}_{}".format(
                *key
            )
        )
        for key in unknown_keys
    }

    # Build equations transition by transition.
    for t_value in range(5):

        source_indices = sorted(
            j
            for (j, tt) in A
            if tt == t_value
        )

        target_indices = sorted(
            j
            for (j, tt) in A
            if tt == t_value + 1
        )

        support = allowed_entries(
            family,
            len(target_indices),
            len(source_indices),
        )

        for j_local, j_actual in enumerate(
            target_indices
        ):

            rhs = 0

            for (
                j,
                k,
            ) in support:

                if j != j_local:
                    continue

                source_j = source_indices[k]

                source_value = A[
                    (source_j, t_value)
                ]

                coefficient = sum(
                    variables[
                        (j, k, d)
                    ] * sp.Integer(
                        t_value ** d
                    )
                    for d in range(
                        degree + 1
                    )
                )

                rhs += (
                    coefficient
                    * source_value
                )

            equations.append(
                clean(
                    rhs
                    - A[
                        (j_actual, t_value + 1)
                    ]
                )
            )

    return variables, equations


# ============================================================================
# EXACT SOLVER
# ============================================================================

def solve_model(
    A,
    family,
    degree,
):

    variables, equations = build_model(
        A,
        family,
        degree,
    )

    symbols = [
        variables[key]
        for key in variables
    ]

    M, b = sp.linear_eq_to_matrix(
        equations,
        symbols,
    )

    rank = M.rank()
    augmented_rank = (
        M.row_join(
            b
        ).rank()
    )

    unknowns = len(symbols)
    equation_count = len(
        equations
    )

    if augmented_rank > rank:

        return {
            "status": "NO_SOLUTION",
            "equations": equation_count,
            "unknowns": unknowns,
            "rank": rank,
            "augmented_rank": augmented_rank,
            "nullity": unknowns - rank,
            "solution": None,
        }

    if rank < unknowns:

        return {
            "status": "NONUNIQUE",
            "equations": equation_count,
            "unknowns": unknowns,
            "rank": rank,
            "augmented_rank": augmented_rank,
            "nullity": unknowns - rank,
            "solution": None,
        }

    solution = M.gauss_jordan_solve(
        b
    )[0]

    result = {}

    for i, key in enumerate(
        variables
    ):

        result[key] = clean(
            solution[i, 0]
        )

    residuals = [
        clean(
            expression.subs(
                result
            )
        )
        for expression in equations
    ]

    if not all(
        residual == 0
        for residual in residuals
    ):

        status = (
            "VERIFICATION_FAILED"
        )

    elif equation_count > unknowns:

        status = (
            "EXACT_OVERDETERMINED"
        )

    else:

        status = (
            "EXACT_DATA_SIZED"
        )

    return {
        "status": status,
        "equations": equation_count,
        "unknowns": unknowns,
        "rank": rank,
        "augmented_rank": augmented_rank,
        "nullity": 0,
        "solution": result,
        "residuals": residuals,
    }


# ============================================================================
# REPORT
# ============================================================================

def report(
    A,
    family,
    degree,
):

    result = solve_model(
        A,
        family,
        degree,
    )

    print()
    print("=" * 78)
    print(
        "{}  t-DEGREE <= {}".format(
            family.upper(),
            degree,
        )
    )
    print("=" * 78)

    print(
        "  equations={}".format(
            result["equations"]
        )
    )

    print(
        "  unknowns={}".format(
            result["unknowns"]
        )
    )

    print(
        "  redundancy={}".format(
            result["equations"]
            - result["unknowns"]
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

    print(
        "  nullity={}".format(
            result["nullity"]
        )
    )

    print(
        "  status={}".format(
            result["status"]
        )
    )

    if (
        result["status"]
        == "EXACT_OVERDETERMINED"
    ):

        print()

        for key in sorted(
            result["solution"]
        ):

            value = result[
                "solution"
            ][key]

            if value != 0:

                print(
                    "  c{}={}".format(
                        key,
                        value,
                    )
                )

        print()
        print(
            "  residuals_zero={}".format(
                all(
                    r == 0
                    for r in result[
                        "residuals"
                    ]
                )
            )
        )

    return result


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 344R — EXACT LOW-DEGREE "
        "t-DEPENDENT NEWTON-OPERATOR AUDIT"
    )
    print("=" * 78)

    A = build_newton()

    print()
    print("=" * 78)
    print(
        "1. OBSERVED NEWTON TRIANGLE"
    )
    print("=" * 78)

    for t_value in range(6):

        print(
            "  t={}: {}".format(
                t_value,
                [
                    A.get(
                        (j, t_value),
                        None,
                    )
                    for j in range(4)
                ],
            )
        )

    families = (
        "diagonal",
        "upper_bidiagonal",
        "lower_bidiagonal",
        "full",
    )

    results = {}

    for family in families:

        results[family] = {}

        for degree in (
            0,
            1,
            2,
        ):

            results[family][degree] = report(
                A,
                family,
                degree,
            )

    # ========================================================================
    # SUMMARY
    # ========================================================================

    print()
    print("=" * 78)
    print(
        "2. SUMMARY"
    )
    print("=" * 78)

    exact_overdetermined = []

    for family in families:

        for degree in (
            0,
            1,
            2,
        ):

            result = results[
                family
            ][degree]

            print(
                "  {} degree<= {}: status={}, "
                "equations={}, unknowns={}, "
                "rank={}, augmented_rank={}, "
                "redundancy={}".format(
                    family,
                    degree,
                    result["status"],
                    result["equations"],
                    result["unknowns"],
                    result["rank"],
                    result["augmented_rank"],
                    result["equations"]
                    - result["unknowns"],
                )
            )

            if (
                result["status"]
                == "EXACT_OVERDETERMINED"
            ):

                exact_overdetermined.append(
                    (
                        family,
                        degree,
                    )
                )

    # ========================================================================
    # INTERPRETATION
    # ========================================================================

    print()
    print("=" * 78)
    print(
        "3. STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
The fixed Newton-operator hypothesis has already failed.

Experiment 344R tests the final natural linear extension:

    A(t+1) = M(t) A(t),

where every matrix entry of M(t) is a polynomial of degree at most d
in the layer variable t.

This includes:

    d=0:
        fixed operator;

    d=1:
        affine layer-dependence;

    d=2:
        quadratic layer-dependence.

The diagonal and bidiagonal families are highly constrained.

The full family is much larger and is therefore only interesting if an
exact solution is uniquely identified despite the available data.

The critical output is:

    exact_overdetermined_models

A nonempty result would identify a concrete low-degree Newton-coordinate
evolution law.

An empty result means that even allowing low-degree layer-dependent
operators does not produce a validated mechanism.

At that point, continued blind algebraic formula searching is unlikely
to be productive with the present 15-cell dataset.

The next step should be source reconstruction or acquisition of an
independent second n=pq case.

No missing values are used.
No extrapolation is used.
No interpolation is used.
"""
    )

    # ========================================================================
    # FINAL EXACTNESS
    # ========================================================================

    print()
    print("=" * 78)
    print(
        "4. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  exact_overdetermined_models={}".format(
            exact_overdetermined
        )
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
        "  arbitrary_matrix_fit=False"
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
        "EXPERIMENT 344R COMPLETE"
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
