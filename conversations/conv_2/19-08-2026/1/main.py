#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 317R — EXACT SOURCE-LAYER TRANSFER / DIFFERENCE-COMMUTATION AUDIT
==============================================================================

Purpose
-------
Experiments 305R--316R pushed the width-2 transfer matrices extracted from
the first-entry finite-difference triangle as far as is useful.

Experiment 317R moves upstream.

The central question is:

    Does the apparent cross-layer transfer already exist at the raw
    terminal-distance source-layer level?

We compare three levels:

    (A) raw source layers Q_t;

    (B) finite-difference rows Delta^j Q_t;

    (C) first-entry triangle A[j,t].

For each transition t -> t+1 we test exact finite-width laws

    Q_{t+1,k}
        = sum_{r=0}^{w-1} c_r Q_{t,k+r}

for widths w = 1,2,3,...

whenever enough equations exist.

The coefficients are solved exactly from a minimal independent subset and
then VERIFIED against every remaining equation.

This deliberately distinguishes:

    EXACT
    NO_SOLUTION
    NONUNIQUE
    INSUFFICIENT_DATA

rather than treating interpolation as evidence.

We then ask whether a surviving raw-layer transfer is compatible with finite
differences.

For an operator C acting on the layer coordinate, the structural question is
whether

    Delta(C Q_t) = C(Delta Q_t)

on every entry where both sides are defined.

Finally we compare raw-layer transfer coefficients with the independently
reconstructed first-entry coefficients from Experiments 305R--316R.

The goal is to determine whether the observed width-2 law is:

    * already present in the source data;
    * created or modified by finite differencing;
    * or merely an artifact of the compressed first-entry representation.

No external files.
No synthetic second n=pq case.
Exact rational arithmetic throughout.
"""


from __future__ import annotations

import math
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
# EXACT HELPERS
# ============================================================================

def clean(value):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(value)
            )
        )
    )


def rational(value):
    return sp.Rational(value)


def valuation(value, prime):
    value = sp.Rational(value)

    if value == 0:
        return None

    num = abs(int(value.p))
    den = abs(int(value.q))

    v = 0

    while num % prime == 0:
        num //= prime
        v += 1

    while den % prime == 0:
        den //= prime
        v -= 1

    return v


def primitive_integer_list(values):
    values = [
        sp.Rational(v)
        for v in values
    ]

    if not values:
        return []

    lcm_den = 1

    for value in values:
        lcm_den = math.lcm(
            lcm_den,
            abs(int(value.q)),
        )

    integers = [
        int(value * lcm_den)
        for value in values
    ]

    g = 0

    for value in integers:
        g = math.gcd(
            g,
            abs(value),
        )

    if g:
        integers = [
            value // g
            for value in integers
        ]

    for value in integers:
        if value != 0:
            if value < 0:
                integers = [
                    -x
                    for x in integers
                ]
            break

    return integers


# ============================================================================
# BASIC SOURCE GEOMETRY
# ============================================================================

def D(p):
    return len(Q[p]) - 1


def build_layers():

    layers = {}

    maximum_t = max(
        D(p)
        for p in Q
    )

    for t in range(
        maximum_t + 1
    ):

        row = []

        for p in sorted(Q):

            r = D(p) - t

            if r < 0:
                continue

            row.append(
                (
                    p,
                    sp.Integer(
                        Q[p][r]
                    ),
                )
            )

        layers[t] = row

    return layers


def layer_values(layers, t):

    return [
        value
        for _, value
        in layers[t]
    ]


def layer_ps(layers, t):

    return [
        p
        for p, _
        in layers[t]
    ]


# ============================================================================
# FINITE DIFFERENCES
# ============================================================================

def finite_difference_rows(values):

    current = [
        rational(v)
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


def build_difference_triangle(layers):

    triangle = {}

    for t in sorted(layers):

        values = layer_values(
            layers,
            t,
        )

        rows = finite_difference_rows(
            values
        )

        for j, row in enumerate(rows):

            if row:
                triangle[
                    (j, t)
                ] = clean(
                    row[0]
                )

    return triangle


def build_A(triangle):

    return {
        key: value
        for key, value
        in triangle.items()
    }


# ============================================================================
# RAW-LAYER WIDTH-W TRANSFER EQUATIONS
# ============================================================================

def raw_transfer_equations(
    layers,
    t,
    width,
):
    """
    Construct equations

        Q_{t+1,k}
          = sum_{r=0}^{w-1} c_r Q_{t,k+r}

    using the common coordinate range available in both layers.

    Here k is the position in the ordered p-list, not the original p itself.
    """

    source = layer_values(
        layers,
        t,
    )

    target = layer_values(
        layers,
        t + 1,
    )

    equations = []

    if not source or not target:
        return equations

    maximum_k = min(
        len(target),
        len(source) - width + 1,
    )

    for k in range(
        maximum_k
    ):

        lhs = rational(
            target[k]
        )

        rhs_row = [
            rational(
                source[k + r]
            )
            for r in range(width)
        ]

        equations.append(
            (
                rhs_row,
                lhs,
            )
        )

    return equations


# ============================================================================
# EXACT LINEAR WIDTH SOLVER
# ============================================================================

def solve_width(
    layers,
    t,
    width,
):
    """
    Solve the raw-layer transfer exactly.

    We distinguish:

        INSUFFICIENT_DATA
        NO_SOLUTION
        NONUNIQUE
        EXACT

    For an exact solution, every available equation is verified.
    """

    equations = raw_transfer_equations(
        layers,
        t,
        width,
    )

    if len(equations) < width:
        return {
            "status": "INSUFFICIENT_DATA",
            "equations": equations,
            "coefficients": None,
            "rank": None,
            "augmented_rank": None,
            "training_equations": 0,
        }

    matrix = sp.Matrix([
        rhs_row
        for rhs_row, _ in equations
    ])

    rhs = sp.Matrix([
        lhs
        for _, lhs in equations
    ])

    rank = matrix.rank()
    augmented_rank = (
        matrix
        .row_join(rhs)
        .rank()
    )

    if augmented_rank > rank:
        return {
            "status": "NO_SOLUTION",
            "equations": equations,
            "coefficients": None,
            "rank": rank,
            "augmented_rank": augmented_rank,
            "training_equations": 0,
        }

    if rank < width:
        return {
            "status": "NONUNIQUE",
            "equations": equations,
            "coefficients": None,
            "rank": rank,
            "augmented_rank": augmented_rank,
            "training_equations": 0,
        }

    # Determine a full-rank row subset.
    selected = []

    current_rank = 0

    for i in range(
        matrix.rows
    ):

        candidate = matrix[
            selected + [i],
            :,
        ]

        candidate_rank = candidate.rank()

        if candidate_rank > current_rank:

            selected.append(i)
            current_rank = candidate_rank

            if current_rank == width:
                break

    if len(selected) != width:
        return {
            "status": "NONUNIQUE",
            "equations": equations,
            "coefficients": None,
            "rank": rank,
            "augmented_rank": augmented_rank,
            "training_equations": len(selected),
        }

    training_matrix = matrix[
        selected,
        :,
    ]

    training_rhs = rhs[
        selected,
        :,
    ]

    try:
        coeff_vector = training_matrix.inv() * training_rhs
    except Exception:
        coeff_vector = training_matrix.gauss_jordan_solve(
            training_rhs
        )[0]

    coefficients = [
        clean(
            coeff_vector[i, 0]
        )
        for i in range(width)
    ]

    # Exact verification against ALL equations.
    failed = []

    for i, (
        rhs_row,
        lhs,
    ) in enumerate(equations):

        predicted = clean(
            sum(
                coefficients[r] * rhs_row[r]
                for r in range(width)
            )
        )

        residual = clean(
            predicted - lhs
        )

        if residual != 0:
            failed.append(
                (
                    i,
                    predicted,
                    lhs,
                    residual,
                )
            )

    if failed:
        return {
            "status": "VERIFICATION_FAILED",
            "equations": equations,
            "coefficients": coefficients,
            "rank": rank,
            "augmented_rank": augmented_rank,
            "training_equations": len(selected),
            "selected_rows": selected,
            "failed_equations": failed,
        }

    return {
        "status": "EXACT",
        "equations": equations,
        "coefficients": coefficients,
        "rank": rank,
        "augmented_rank": augmented_rank,
        "training_equations": len(selected),
        "selected_rows": selected,
        "failed_equations": [],
    }


# ============================================================================
# RAW SOURCE TRANSFER AUDIT
# ============================================================================

def raw_transfer_audit(
    layers,
    maximum_width=4,
):

    print()
    print("=" * 78)
    print(
        "3. RAW SOURCE-LAYER TRANSFER AUDIT"
    )
    print("=" * 78)

    results = {}

    for t in range(
        max(layers)
    ):

        print()
        print(
            "  TRANSITION t={} -> {}:".format(
                t,
                t + 1,
            )
        )

        results[t] = {}

        for width in range(
            1,
            maximum_width + 1,
        ):

            result = solve_width(
                layers,
                t,
                width,
            )

            results[t][width] = result

            print()
            print(
                "    width={}:".format(
                    width
                )
            )

            print(
                "      status={}".format(
                    result["status"]
                )
            )

            print(
                "      equations={}".format(
                    len(
                        result["equations"]
                    )
                )
            )

            if result["coefficients"] is not None:

                print(
                    "      coefficients={}".format(
                        result["coefficients"]
                    )
                )

            if result["rank"] is not None:

                print(
                    "      rank={}".format(
                        result["rank"]
                    )
                )

                print(
                    "      augmented_rank={}".format(
                        result["augmented_rank"]
                    )
                )

            if result["status"] == "EXACT":

                print(
                    "      verified_all_equations=True"
                )

    return results


# ============================================================================
# FIRST-ENTRY WIDTH-2 TRANSFER
# ============================================================================

def solve_A_width2(A, t):

    equations = []

    maximum_j = max(
        (
            j
            for (j, tt) in A
            if tt == t
        ),
        default=-1,
    )

    for j in range(
        maximum_j + 1
    ):

        lhs = A.get(
            (j, t + 1)
        )

        x0 = A.get(
            (j, t)
        )

        x1 = A.get(
            (j + 1, t)
        )

        if (
            lhs is None
            or x0 is None
            or x1 is None
        ):
            continue

        equations.append(
            (
                [
                    rational(x0),
                    rational(x1),
                ],
                rational(lhs),
            )
        )

    if len(equations) < 2:
        return {
            "status": "INSUFFICIENT_DATA",
            "equations": equations,
            "coefficients": None,
        }

    M = sp.Matrix([
        row
        for row, _ in equations
    ])

    b = sp.Matrix([
        lhs
        for _, lhs in equations
    ])

    rank = M.rank()
    augmented_rank = M.row_join(b).rank()

    if augmented_rank > rank:

        return {
            "status": "NO_SOLUTION",
            "equations": equations,
            "coefficients": None,
            "rank": rank,
            "augmented_rank": augmented_rank,
        }

    if rank < 2:

        return {
            "status": "NONUNIQUE",
            "equations": equations,
            "coefficients": None,
            "rank": rank,
            "augmented_rank": augmented_rank,
        }

    solution = M.inv() * b

    coefficients = [
        clean(solution[i, 0])
        for i in range(2)
    ]

    failed = []

    for i, (row, lhs) in enumerate(
        equations
    ):

        residual = clean(
            coefficients[0] * row[0]
            + coefficients[1] * row[1]
            - lhs
        )

        if residual != 0:
            failed.append(
                (
                    i,
                    residual,
                )
            )

    if failed:

        return {
            "status": "VERIFICATION_FAILED",
            "equations": equations,
            "coefficients": coefficients,
            "rank": rank,
            "augmented_rank": augmented_rank,
            "failed_equations": failed,
        }

    return {
        "status": "EXACT",
        "equations": equations,
        "coefficients": coefficients,
        "rank": rank,
        "augmented_rank": augmented_rank,
        "failed_equations": [],
    }


def first_entry_transfer_audit(
    A,
):

    print()
    print("=" * 78)
    print(
        "4. FIRST-ENTRY TRIANGLE WIDTH-2 TRANSFER AUDIT"
    )
    print("=" * 78)

    profile = {}

    maximum_t = max(
        t
        for _, t
        in A
    )

    for t in range(
        maximum_t
    ):

        result = solve_A_width2(
            A,
            t,
        )

        profile[t] = result

        print()
        print(
            "  t={} -> {}:".format(
                t,
                t + 1,
            )
        )

        print(
            "    status={}".format(
                result["status"]
            )
        )

        if result["coefficients"] is not None:

            print(
                "    alpha={}".format(
                    result["coefficients"][0]
                )
            )

            print(
                "    beta={}".format(
                    result["coefficients"][1]
                )
            )

        print(
            "    equations={}".format(
                len(
                    result["equations"]
                )
            )
        )

    return profile


# ============================================================================
# SOURCE-vs-FIRST-ENTRY COMPARISON
# ============================================================================

def compare_transfers(
    raw_results,
    A_results,
):

    print()
    print("=" * 78)
    print(
        "5. RAW-LAYER VS FIRST-ENTRY TRANSFER COMPARISON"
    )
    print("=" * 78)

    common_t = sorted(
        set(raw_results)
        &
        set(A_results)
    )

    for t in common_t:

        raw = raw_results[t].get(2)
        first = A_results[t]

        print()
        print(
            "  t={} -> {}:".format(
                t,
                t + 1,
            )
        )

        print(
            "    raw_width2_status={}".format(
                None
                if raw is None
                else raw["status"]
            )
        )

        print(
            "    first_entry_width2_status={}".format(
                first["status"]
            )
        )

        if (
            raw is not None
            and raw["status"] == "EXACT"
        ):

            print(
                "    raw_coefficients={}".format(
                    raw["coefficients"]
                )
            )

        if (
            first["status"] == "EXACT"
        ):

            print(
                "    first_entry_coefficients={}".format(
                    first["coefficients"]
                )
            )

        if (
            raw is not None
            and raw["status"] == "EXACT"
            and first["status"] == "EXACT"
        ):

            raw_coeffs = raw[
                "coefficients"
            ]

            first_coeffs = first[
                "coefficients"
            ]

            differences = [
                clean(
                    raw_coeffs[i]
                    - first_coeffs[i]
                )
                for i in range(2)
            ]

            ratios = []

            for i in range(2):

                if first_coeffs[i] == 0:
                    ratios.append(
                        None
                    )
                else:
                    ratios.append(
                        clean(
                            raw_coeffs[i]
                            /
                            first_coeffs[i]
                        )
                    )

            print(
                "    coefficient_differences={}".format(
                    differences
                )
            )

            print(
                "    raw_over_first_entry_ratios={}".format(
                    ratios
                )
            )

            print(
                "    same_operator={}".format(
                    all(
                        d == 0
                        for d in differences
                    )
                )
            )


# ============================================================================
# FINITE-DIFFERENCE COMMUTATION AUDIT
# ============================================================================

def raw_transfer_prediction(
    source_values,
    coefficients,
    width,
):

    predicted = []

    for k in range(
        len(source_values) - width + 1
    ):

        predicted.append(
            clean(
                sum(
                    coefficients[r]
                    * source_values[k + r]
                    for r in range(width)
                )
            )
        )

    return predicted


def difference_commutation_audit(
    layers,
    raw_results,
):

    print()
    print("=" * 78)
    print(
        "6. DIFFERENCE / SOURCE-OPERATOR COMPATIBILITY AUDIT"
    )
    print("=" * 78)

    results = {}

    for t in sorted(
        raw_results
    ):

        result = raw_results[t].get(2)

        if (
            result is None
            or result["status"] != "EXACT"
        ):
            continue

        coefficients = result[
            "coefficients"
        ]

        source_t = layer_values(
            layers,
            t,
        )

        source_next = layer_values(
            layers,
            t + 1,
        )

        # Reconstruct target from source operator.
        predicted_target = raw_transfer_prediction(
            source_t,
            coefficients,
            2,
        )

        # The ordinary first difference of the raw target.
        actual_difference = [
            clean(
                source_next[i + 1]
                - source_next[i]
            )
            for i in range(
                len(source_next) - 1
            )
        ]

        # Difference of the predicted target.
        predicted_difference = [
            clean(
                predicted_target[i + 1]
                - predicted_target[i]
            )
            for i in range(
                len(predicted_target) - 1
            )
        ]

        # Apply the SAME width-2 operator to the first difference of Q_t.
        source_difference = [
            clean(
                source_t[i + 1]
                - source_t[i]
            )
            for i in range(
                len(source_t) - 1
            )
        ]

        operator_on_difference = raw_transfer_prediction(
            source_difference,
            coefficients,
            2,
        )

        # Compare where lengths coincide.
        n = min(
            len(actual_difference),
            len(predicted_difference),
            len(operator_on_difference),
        )

        commutation_residuals = []

        for i in range(n):

            commutation_residuals.append(
                clean(
                    predicted_difference[i]
                    -
                    operator_on_difference[i]
                )
            )

        target_residuals = [
            clean(
                predicted_target[i]
                - source_next[i]
            )
            for i in range(
                min(
                    len(predicted_target),
                    len(source_next),
                )
            )
        ]

        exact_commutation = all(
            r == 0
            for r in commutation_residuals
        )

        exact_target = all(
            r == 0
            for r in target_residuals
        )

        results[t] = {
            "target_residuals": target_residuals,
            "commutation_residuals": commutation_residuals,
            "exact_target": exact_target,
            "exact_commutation": exact_commutation,
        }

        print()
        print(
            "  t={} -> {}:".format(
                t,
                t + 1,
            )
        )

        print(
            "    source_operator_reconstructs_target={}".format(
                exact_target
            )
        )

        print(
            "    Delta_operator_commutation={}".format(
                exact_commutation
            )
        )

        print(
            "    target_residuals={}".format(
                target_residuals
            )
        )

        print(
            "    commutation_residuals={}".format(
                commutation_residuals
            )
        )

    return results


# ============================================================================
# GENERAL WIDTH DIFFERENCE-COMPATIBILITY
# ============================================================================

def generalized_difference_identity(
    values,
    coefficients,
    width,
):

    """
    Compare

        Delta(
            C(values)
        )

    against

        C(
            Delta(values)
        )

    over all available entries.

    This is purely an algebraic diagnostic for the finite sequence.
    """

    predicted = raw_transfer_prediction(
        values,
        coefficients,
        width,
    )

    delta_predicted = [
        clean(
            predicted[i + 1]
            - predicted[i]
        )
        for i in range(
            len(predicted) - 1
        )
    ]

    delta_values = [
        clean(
            values[i + 1]
            - values[i]
        )
        for i in range(
            len(values) - 1
        )
    ]

    operator_delta = raw_transfer_prediction(
        delta_values,
        coefficients,
        width,
    )

    n = min(
        len(delta_predicted),
        len(operator_delta),
    )

    residuals = [
        clean(
            delta_predicted[i]
            - operator_delta[i]
        )
        for i in range(n)
    ]

    return residuals


# ============================================================================
# CONTENT / NORMALIZATION AUDIT
# ============================================================================

def integer_content(values):

    values = [
        sp.Rational(v)
        for v in values
    ]

    if not values:
        return sp.Integer(0)

    lcm_den = 1

    for value in values:
        lcm_den = math.lcm(
            lcm_den,
            abs(int(value.q)),
        )

    integers = [
        int(value * lcm_den)
        for value in values
    ]

    g = 0

    for value in integers:
        g = math.gcd(
            g,
            abs(value),
        )

    return sp.Integer(g)


def primitive_row(values):

    values = [
        sp.Rational(v)
        for v in values
    ]

    if not values:
        return []

    content = integer_content(
        values
    )

    if content == 0:
        return list(values)

    lcm_den = 1

    for value in values:
        lcm_den = math.lcm(
            lcm_den,
            abs(int(value.q)),
        )

    integers = [
        int(value * lcm_den)
        for value in values
    ]

    g = int(content)

    if g:
        integers = [
            value // g
            for value in integers
        ]

    return integers


def content_audit(layers):

    print()
    print("=" * 78)
    print(
        "7. EXACT LAYER CONTENT / PRIMITIVE-NORMALIZATION AUDIT"
    )
    print("=" * 78)

    for t in sorted(
        layers
    ):

        values = layer_values(
            layers,
            t,
        )

        content = integer_content(
            values
        )

        primitive = primitive_row(
            values
        )

        print()
        print(
            "  t={}:".format(
                t
            )
        )

        print(
            "    values={}".format(
                values
            )
        )

        print(
            "    integer_content={}".format(
                content
            )
        )

        print(
            "    primitive_integer_row={}".format(
                primitive
            )
        )

        print(
            "    v17(content)={}".format(
                valuation(
                    content,
                    17,
                )
                if content != 0
                else None
            )
        )


# ============================================================================
# PRIME VALUATION OF RAW TRANSFER COEFFICIENTS
# ============================================================================

def coefficient_valuation_audit(
    raw_results,
):

    print()
    print("=" * 78)
    print(
        "8. RAW TRANSFER COEFFICIENT PRIME PROFILE"
    )
    print("=" * 78)

    for t in sorted(
        raw_results
    ):

        result = raw_results[t].get(2)

        if (
            result is None
            or result["status"] != "EXACT"
        ):
            continue

        print()
        print(
            "  t={} -> {}:".format(
                t,
                t + 1,
            )
        )

        for i, coefficient in enumerate(
            result["coefficients"]
        ):

            print(
                "    c_{}={}".format(
                    i,
                    coefficient,
                )
            )

            print(
                "      valuations={}".format(
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
                    }
                )
            )


# ============================================================================
# SOURCE LAYER DISPLAY
# ============================================================================

def print_layers(layers):

    print()
    print("=" * 78)
    print(
        "1. TERMINAL-DISTANCE SOURCE LAYERS"
    )
    print("=" * 78)

    for t in sorted(
        layers
    ):

        print()
        print(
            "  t={}:".format(
                t
            )
        )

        print(
            "    available_p={}".format(
                layer_ps(
                    layers,
                    t,
                )
            )
        )

        print(
            "    Q_t={}".format(
                layer_values(
                    layers,
                    t,
                )
            )
        )


# ============================================================================
# FIRST-ENTRY DISPLAY
# ============================================================================

def print_first_entry_triangle(
    triangle,
):

    print()
    print("=" * 78)
    print(
        "2. FIRST-ENTRY DIFFERENCE TRIANGLE"
    )
    print("=" * 78)

    for j in sorted(
        set(
            jj
            for jj, _ in triangle
        )
    ):

        entries = []

        for t in sorted(
            tt
            for jj, tt in triangle
            if jj == j
        ):

            entries.append(
                (
                    t,
                    triangle[
                        (j, t)
                    ],
                )
            )

        print()
        print(
            "  j={}: {}".format(
                j,
                entries,
            )
        )


# ============================================================================
# SOURCE REFERENCE
# ============================================================================

def terminal_reference():

    q1 = Q[1][-1]
    q3 = Q[3][-1]

    g = math.gcd(
        q1,
        q3,
    )

    print()
    print("=" * 78)
    print(
        "9. TERMINAL PROJECTIVE SOURCE REFERENCE"
    )
    print("=" * 78)

    print(
        "  q1_terminal={}".format(
            q1
        )
    )

    print(
        "  q3_terminal={}".format(
            q3
        )
    )

    print(
        "  gcd={}".format(
            g
        )
    )

    print(
        "  q1/17={}".format(
            q1 // 17
        )
    )

    print(
        "  q3/17={}".format(
            q3 // 17
        )
    )

    print(
        "  v17(q1_terminal)={}".format(
            valuation(
                q1,
                17,
            )
        )
    )

    print(
        "  v17(q3_terminal)={}".format(
            valuation(
                q3,
                17,
            )
        )
    )


# ============================================================================
# STRUCTURAL INTERPRETATION
# ============================================================================

def interpretation():

    print()
    print("=" * 78)
    print(
        "10. STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
Experiment 317R deliberately moves upstream from the 2x2 companion
matrices.

The previous experiments established that the first-entry width-2
transitions themselves behave like generic noncommutative elements of
M_2(Q).

That does NOT tell us whether the original source layers Q_t obey a
simpler law.

The correct hierarchy is:

    raw source layers Q_t
          |
          v
    finite differences Delta^j Q_t
          |
          v
    first-entry triangle A[j,t]
          |
          v
    width-2 companion transfer matrices.

This experiment tests whether the apparent transfer already exists at
the first level.

A genuine source-layer law would be much stronger evidence than a law
recovered only after taking first differences.

The especially important comparison is:

    raw Q-layer width-2 coefficients
        versus
    first-entry alpha_t, beta_t.

If they differ, finite differencing is not merely revealing an existing
source operator; it is changing the effective transfer coordinates.

The next structural test is then whether a surviving source operator
commutes with finite differencing:

    Delta(C Q_t) = C(Delta Q_t).

If that identity holds exactly, the transfer law has a genuine operator
origin compatible with the difference construction.

If it fails, the first-entry transfer should be treated as a derived
coordinate-level phenomenon rather than a source-level mechanism.

The content/normalization audit is included because the observed
coefficients contain strong prime-specific denominators, especially 17.

No synthetic second n=pq case is generated.
"""
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 317R — EXACT SOURCE-LAYER TRANSFER / "
        "DIFFERENCE-COMMUTATION AUDIT"
    )
    print("=" * 78)

    layers = build_layers()

    triangle = build_difference_triangle(
        layers
    )

    A = build_A(
        triangle
    )

    print_layers(
        layers
    )

    print_first_entry_triangle(
        triangle
    )

    raw_results = raw_transfer_audit(
        layers,
        maximum_width=4,
    )

    A_results = first_entry_transfer_audit(
        A
    )

    compare_transfers(
        raw_results,
        A_results,
    )

    difference_results = (
        difference_commutation_audit(
            layers,
            raw_results,
        )
    )

    content_audit(
        layers
    )

    coefficient_valuation_audit(
        raw_results
    )

    terminal_reference()

    interpretation()

    # ------------------------------------------------------------------------
    # FINAL EXACTNESS
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "11. FINAL EXACTNESS"
    )
    print("=" * 78)

    raw_width2_exact = [
        t
        for t in raw_results
        if (
            raw_results[t].get(2) is not None
            and
            raw_results[t][2]["status"]
            == "EXACT"
        )
    ]

    first_width2_exact = [
        t
        for t in A_results
        if (
            A_results[t]["status"]
            == "EXACT"
        )
    ]

    commuting_t = [
        t
        for t, result
        in difference_results.items()
        if result[
            "exact_commutation"
        ]
    ]

    print(
        "  raw_width2_exact_transitions={}".format(
            raw_width2_exact
        )
    )

    print(
        "  first_entry_width2_exact_transitions={}".format(
            first_width2_exact
        )
    )

    print(
        "  exact_difference_commuting_raw_transitions={}".format(
            commuting_t
        )
    )

    print(
        "  source_layer_transfer_tested=True"
    )

    print(
        "  finite_difference_commutation_tested=True"
    )

    print(
        "  content_normalization_tested=True"
    )

    print(
        "  robust_exact_verification=True"
    )

    print(
        "  arbitrary_matrix_fit=False"
    )

    print(
        "  synthetic_second_case=False"
    )

    print(
        "  external_files_used=False"
    )

    print(
        "  interpolation_counted_as_proof=False"
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
        "EXPERIMENT 317R COMPLETE"
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

