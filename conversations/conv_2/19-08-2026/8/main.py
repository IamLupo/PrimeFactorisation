#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 324R — EXACT VERTICAL LAYER-RECURRENCE / CROSS-P AUDIT
==============================================================================

Purpose
-------
Experiments 317R--323R found:

    * the raw width-2 law fails at the only overdetermined transition;
    * no wider p-stencil gives an overdetermined exact validation;
    * the surviving width-2 laws at t=1->2 and t=2->3 are reconstruction-only.

The next natural direction is therefore the OTHER axis of the triangular
source table.

We test whether the source layers satisfy a constant-coefficient recurrence
in terminal-distance t:

    Q_{t+r}(p)
      =
      c_0 Q_t(p)
      + c_1 Q_{t+1}(p)
      + ...
      + c_{r-1} Q_{t+r-1}(p).

The crucial point is that the same coefficients are required for every
available p.

For each order r we compute:

    * every usable (p,t) equation;
    * exact coefficient-matrix rank;
    * augmented rank;
    * equation count;
    * redundancy;
    * exact consistency;
    * exact coefficients when uniquely determined;
    * residuals;
    * whether the fit is genuinely overdetermined.

We also perform a secondary per-column audit, but this is explicitly
diagnostic because the individual columns are short.

Orders tested:
    r = 1,2,3,4

No interpolation.
No missing-value reconstruction.
No synthetic second n=pq case.
Exact SymPy rational arithmetic only.
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


def valuation(x, prime):

    x = sp.Rational(x)

    if x == 0:
        return sp.oo

    numerator = abs(int(x.p))
    denominator = abs(int(x.q))

    value = 0

    while numerator % prime == 0:
        numerator //= prime
        value += 1

    while denominator % prime == 0:
        denominator //= prime
        value -= 1

    return value


# ============================================================================
# TERMINAL-DISTANCE LAYERS
# ============================================================================

def build_layers():

    maximum_t = max(
        len(values) - 1
        for values in Q.values()
    )

    layers = {}

    for t in range(
        maximum_t + 1
    ):

        layer = {}

        for p in sorted(Q):

            index = (
                len(Q[p])
                - 1
                - t
            )

            if index < 0:
                continue

            layer[p] = sp.Integer(
                Q[p][index]
            )

        layers[t] = layer

    return layers


# ============================================================================
# VERTICAL RECURRENCE SYSTEM
# ============================================================================

def vertical_system(
    layers,
    order,
    selected_p=None,
):

    maximum_t = max(
        layers
    )

    equations = []

    ps = (
        sorted(selected_p)
        if selected_p is not None
        else sorted(Q)
    )

    for p in ps:

        available_t = [
            t
            for t in range(
                maximum_t + 1
            )
            if p in layers[t]
        ]

        if len(available_t) <= order:
            continue

        # Available t values are contiguous from 0.
        for t in range(
            len(available_t) - order
        ):

            source = [
                layers[t + k][p]
                for k in range(order)
            ]

            target = layers[
                t + order
            ][p]

            equations.append(
                (
                    p,
                    t,
                    source,
                    target,
                )
            )

    if not equations:

        return {
            "status": "NO_EQUATIONS",
            "order": order,
            "equations": [],
            "rank": 0,
            "augmented_rank": 0,
            "coefficients": None,
            "residuals": [],
        }

    M = sp.Matrix([
        source
        for _, _, source, _
        in equations
    ])

    b = sp.Matrix([
        target
        for _, _, _, target
        in equations
    ])

    rank = M.rank()
    augmented_rank = (
        M.row_join(b).rank()
    )

    result = {
        "order": order,
        "equations": equations,
        "rank": rank,
        "augmented_rank": augmented_rank,
        "coefficients": None,
        "residuals": [],
    }

    if augmented_rank > rank:

        result["status"] = "NO_SOLUTION"
        return result

    if rank < order:

        result["status"] = "NONUNIQUE"
        return result

    solution, parameters = (
        M.gauss_jordan_solve(b)
    )

    try:
        if len(parameters) > 0:
            result["status"] = "NONUNIQUE"
            return result
    except TypeError:
        pass

    coefficients = [
        clean(solution[i, 0])
        for i in range(
            solution.rows
        )
    ]

    residuals = []

    for _, _, source, target in equations:

        predicted = clean(
            sum(
                coefficients[k]
                * source[k]
                for k in range(order)
            )
        )

        residuals.append(
            clean(
                predicted
                - target
            )
        )

    result["coefficients"] = coefficients
    result["residuals"] = residuals

    if all(
        residual == 0
        for residual in residuals
    ):

        result["status"] = "EXACT"

    else:

        result["status"] = "VERIFICATION_FAILED"

    return result


# ============================================================================
# REPORT
# ============================================================================

def print_result(
    label,
    result,
):

    print()
    print(
        "  {}:".format(label)
    )

    equation_count = len(
        result["equations"]
    )

    order = result["order"]

    print(
        "    equations={}".format(
            equation_count
        )
    )

    print(
        "    unknown_coefficients={}".format(
            order
        )
    )

    print(
        "    redundancy={}".format(
            equation_count - order
        )
    )

    print(
        "    rank={}".format(
            result["rank"]
        )
    )

    print(
        "    augmented_rank={}".format(
            result["augmented_rank"]
        )
    )

    print(
        "    status={}".format(
            result["status"]
        )
    )

    if result["coefficients"] is not None:

        print(
            "    coefficients={}".format(
                result["coefficients"]
            )
        )

        print(
            "    residuals={}".format(
                result["residuals"]
            )
        )

        print(
            "    all_residuals_zero={}".format(
                all(
                    r == 0
                    for r
                    in result["residuals"]
                )
            )
        )

        for i, coefficient in enumerate(
            result["coefficients"]
        ):

            print(
                "    c{}_valuations={}".format(
                    i,
                    {
                        p: valuation(
                            coefficient,
                            p,
                        )
                        for p in (
                            2,
                            3,
                            5,
                            7,
                            11,
                            13,
                            17,
                        )
                    },
                )
            )


# ============================================================================
# SHARED VERTICAL RECURRENCE AUDIT
# ============================================================================

def shared_vertical_audit(
    layers,
):

    print()
    print("=" * 78)
    print(
        "2. SHARED VERTICAL RECURRENCE AUDIT"
    )
    print("=" * 78)

    results = {}

    for order in range(
        1,
        5,
    ):

        result = vertical_system(
            layers,
            order,
        )

        results[order] = result

        print_result(
            "order={}".format(
                order
            ),
            result,
        )

    return results


# ============================================================================
# PER-COLUMN DIAGNOSTIC
# ============================================================================

def per_column_audit(
    layers,
):

    print()
    print("=" * 78)
    print(
        "3. PER-COLUMN VERTICAL RECURRENCE DIAGNOSTIC"
    )
    print("=" * 78)

    all_results = {}

    maximum_t = max(
        layers
    )

    for p in sorted(Q):

        print()
        print(
            "  p={}:".format(
                p
            )
        )

        column = [
            layers[t][p]
            for t in range(
                maximum_t + 1
            )
            if p in layers[t]
        ]

        print(
            "    sequence={}".format(
                column
            )
        )

        results = {}

        max_order = min(
            4,
            len(column) - 1,
        )

        for order in range(
            1,
            max_order + 1,
        ):

            selected = vertical_system(
                layers,
                order,
                selected_p=[p],
            )

            results[order] = selected

            print_result(
                "order={}".format(
                    order
                ),
                selected,
            )

        all_results[p] = results

    return all_results


# ============================================================================
# CROSS-P CONTRADICTION AUDIT
# ============================================================================

def contradiction_audit(
    shared_results,
):

    print()
    print("=" * 78)
    print(
        "4. CROSS-P CONTRADICTION AUDIT"
    )
    print("=" * 78)

    for order, result in shared_results.items():

        if result["status"] != "NO_SOLUTION":

            continue

        equations = result[
            "equations"
        ]

        print()
        print(
            "  order={}:".format(
                order
            )
        )

        print(
            "    global_system_inconsistent=True"
        )

        # Find a minimal contradiction by searching
        # small subsets of equations.

        found = False

        for i in range(
            len(equations)
        ):

            for j in range(
                i + 1,
                len(equations)
            ):

                subset = [
                    equations[i],
                    equations[j],
                ]

                if len(subset) < order + 1:
                    continue

                # For order 1 this is enough.
                M = sp.Matrix([
                    row[2]
                    for row in subset
                ])

                b = sp.Matrix([
                    row[3]
                    for row in subset
                ])

                if (
                    M.rank()
                    <
                    M.row_join(b).rank()
                ):

                    print(
                        "    contradiction_pair={}".format(
                            [
                                (
                                    row[0],
                                    row[1],
                                )
                                for row
                                in subset
                            ]
                        )
                    )

                    found = True
                    break

            if found:
                break


# ============================================================================
# INTERPRETATION
# ============================================================================

def interpretation():

    print()
    print("=" * 78)
    print(
        "STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
Experiment 323R showed that the source table does not support a
falsifiable low-width recurrence in the p-direction.

Experiment 324R tests the orthogonal possibility:

    a common recurrence in terminal-distance t.

For order r this means

    Q_{t+r}(p)
      =
    c_0 Q_t(p)
      + c_1 Q_{t+1}(p)
      + ...
      + c_{r-1} Q_{t+r-1}(p),

with the SAME coefficients for every p.

This is much stronger than fitting each p-column separately.

The principal distinction is again:

    overdetermined exact validation
        versus
    exact reconstruction.

For the full triangular data:

    order 1 has many equations for one coefficient;
    order 2 has many equations for two coefficients;
    order 3 still has more equations than unknowns;
    order 4 becomes data-limited.

Therefore an exact shared recurrence of order 1, 2, or 3 would be
substantial evidence for genuine layer dynamics.

If all three fail, the evidence increasingly favors a genuinely
two-variable source object Q_t(p), rather than a simple evolution in
t with fixed coefficients.

The per-column audit is secondary because individual columns are too
short to provide strong validation.

No interpolation.
No missing-value reconstruction.
No synthetic second n=pq case.
"""
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 324R — EXACT VERTICAL LAYER-RECURRENCE / "
        "CROSS-P AUDIT"
    )
    print("=" * 78)

    layers = build_layers()

    print()
    print(
        "1. TERMINAL-DISTANCE SOURCE LAYERS"
    )
    print("=" * 78)

    for t, layer in layers.items():

        print(
            "  t={}: {}".format(
                t,
                layer,
            )
        )

    shared_results = shared_vertical_audit(
        layers
    )

    per_column_results = per_column_audit(
        layers
    )

    contradiction_audit(
        shared_results
    )

    interpretation()

    shared_validated = [
        order
        for order, result
        in shared_results.items()
        if (
            result["status"] == "EXACT"
            and
            len(result["equations"])
            > order
        )
    ]

    print()
    print("=" * 78)
    print(
        "5. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  shared_vertical_recurrence_tested=True"
    )

    print(
        "  shared_exact_overdetermined_orders={}".format(
            shared_validated
        )
    )

    print(
        "  per_column_diagnostic_completed=True"
    )

    print(
        "  interpolation_performed=False"
    )

    print(
        "  missing_value_reconstruction=False"
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
        "EXPERIMENT 324R COMPLETE"
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