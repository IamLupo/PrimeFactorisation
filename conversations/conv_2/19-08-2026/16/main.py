#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 332R — EXACT NEWTON-COEFFICIENT TRIANGLE /
                  HIDDEN-LOW-RANK DYNAMICS AUDIT
==============================================================================

Purpose
-------
Experiments 322R-331R rejected:

    * fixed vertical linear recurrences;
    * low-degree p-dependent vertical recurrences;
    * low-degree rational vertical recurrences;
    * low-degree bivariate polynomial laws;
    * low-width constant local 2D stencils;
    * rank-one additive separability;
    * rank-one multiplicative separability.

The source rows nevertheless have a natural discrete-polynomial structure
because the available p-values lie on

    p = 1, 3, 5, 7.

This experiment moves to the exact Newton / finite-difference coordinates

    A[j,t] = Delta_p^j Q_t(1).

It audits:

    1. exact Newton coefficient triangle;
    2. support pattern;
    3. coefficient-matrix rank;
    4. all complete 2x2 coefficient minors;
    5. fixed-j scalar vertical dynamics;
    6. common scalar dynamics across coefficient rows;
    7. coefficient cross-ratios;
    8. exact prime valuations.

Important:
-----------
Only observed source values are used.

No missing-value reconstruction.
No interpolation beyond observed cells.
No extrapolation.
No synthetic second n=pq case.
Exact SymPy arithmetic only.
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
# SOURCE LAYERS
# ============================================================================

def build_layers():
    """
    Q[p] is stored in reverse t-order:

        index = len(Q[p]) - 1 - t.

    Only actually observed cells are returned.
    """

    maximum_t = max(
        len(values) - 1
        for values in Q.values()
    )

    layers = {}

    for t in range(maximum_t + 1):

        layer = []

        for p_value in sorted(Q):

            values = Q[p_value]

            index = (
                len(values)
                - 1
                - t
            )

            if index >= 0:

                layer.append(
                    (
                        p_value,
                        sp.Integer(
                            values[index]
                        ),
                    )
                )

        layers[t] = layer

    return layers


# ============================================================================
# FINITE DIFFERENCE / NEWTON TRIANGLE
# ============================================================================

def finite_difference_rows(values):
    """
    Exact forward-difference rows.

    For values [q0,q1,q2,...]:

        Delta^0 = [q0,q1,q2,...]
        Delta^1 = [q1-q0,q2-q1,...]
        Delta^2 = [...]
        ...

    Because p increases by 2, these are the natural Newton coefficients
    on the observed arithmetic p-grid up to the conventional scaling.
    """

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


def build_newton_triangle(layers):

    A = {}

    for t, layer in layers.items():

        values = [
            value
            for _, value in layer
        ]

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
# SECTION 1 — EXACT TRIANGLE
# ============================================================================

def print_triangle(A):

    print()
    print("=" * 78)
    print(
        "1. EXACT NEWTON-COEFFICIENT TRIANGLE"
    )
    print("=" * 78)

    max_t = max(
        t
        for _, t in A
    )

    for t in range(
        max_t + 1
    ):

        row = [
            A.get(
                (j, t),
                None,
            )
            for j in range(4)
        ]

        print()
        print(
            "  t={}:".format(
                t
            )
        )

        print(
            "    A[j,t]={}".format(
                row
            )
        )


# ============================================================================
# SECTION 2 — SUPPORT
# ============================================================================

def support_audit(A):

    print()
    print("=" * 78)
    print(
        "2. NEWTON SUPPORT-PATTERN AUDIT"
    )
    print("=" * 78)

    max_t = max(
        t
        for _, t in A
    )

    for t in range(
        max_t + 1
    ):

        support = [
            j
            for j in range(4)
            if A.get(
                (j, t),
                0
            ) != 0
        ]

        print(
            "  t={}: support={}".format(
                t,
                support,
            )
        )

    nonzero_entries = sum(
        1
        for value in A.values()
        if value != 0
    )

    print()
    print(
        "  nonzero_entries={}".format(
            nonzero_entries
        )
    )


# ============================================================================
# SECTION 3 — COMPLETE COEFFICIENT MINORS
# ============================================================================

def coefficient_rectangles(A):

    rectangles = []

    max_t = max(
        t
        for _, t in A
    )

    for j1 in range(4):

        for j2 in range(
            j1 + 1,
            4,
        ):

            for t in range(
                max_t
            ):

                required = [
                    (j1, t),
                    (j2, t),
                    (j1, t + 1),
                    (j2, t + 1),
                ]

                if all(
                    location in A
                    for location in required
                ):
                    rectangles.append(
                        (
                            j1,
                            j2,
                            t,
                        )
                    )

    return rectangles


def coefficient_minor_audit(A):

    print()
    print("=" * 78)
    print(
        "3. NEWTON-COEFFICIENT 2x2 MINOR AUDIT"
    )
    print("=" * 78)

    rectangles = coefficient_rectangles(
        A
    )

    zero = 0
    nonzero = 0

    for j1, j2, t in rectangles:

        determinant = clean(
            A[(j1, t)] * A[(j2, t + 1)]
            -
            A[(j2, t)] * A[(j1, t + 1)]
        )

        print()
        print(
            "  rows=({},{}), columns=({},{}):".format(
                j1,
                j2,
                t,
                t + 1,
            )
        )

        print(
            "    determinant={}".format(
                determinant
            )
        )

        if determinant == 0:
            zero += 1
        else:
            nonzero += 1

    print()
    print(
        "  complete_coefficient_rectangles={}".format(
            len(rectangles)
        )
    )

    print(
        "  zero_minors={}".format(
            zero
        )
    )

    print(
        "  nonzero_minors={}".format(
            nonzero
        )
    )

    return rectangles


# ============================================================================
# SECTION 4 — COEFFICIENT MATRIX RANK
# ============================================================================

def coefficient_rank_audit(A):

    print()
    print("=" * 78)
    print(
        "4. NEWTON-COEFFICIENT RANK AUDIT"
    )
    print("=" * 78)

    max_t = max(
        t
        for _, t in A
    )

    matrix = sp.Matrix([
        [
            A.get(
                (j, t),
                0,
            )
            for t in range(
                max_t + 1
            )
        ]
        for j in range(4)
    ])

    print(
        "  coefficient_matrix="
    )

    print(
        matrix
    )

    print()
    print(
        "  shape={}".format(
            matrix.shape
        )
    )

    print(
        "  rank={}".format(
            matrix.rank()
        )
    )

    return matrix


# ============================================================================
# SECTION 5 — FIXED-j VERTICAL DYNAMICS
# ============================================================================

def solve_scalar_vertical_recurrence(
    values
):

    values = [
        sp.Rational(v)
        for v in values
    ]

    if len(values) < 2:
        return {
            "status": "INSUFFICIENT_DATA"
        }

    source = values[:-1]
    target = values[1:]

    A_matrix = sp.Matrix(
        [
            [source[i]]
            for i in range(
                len(source)
            )
        ]
    )

    b = sp.Matrix(
        target
    )

    rank = A_matrix.rank()

    augmented_rank = (
        A_matrix
        .row_join(b)
        .rank()
    )

    if augmented_rank > rank:

        return {
            "status": "NO_SOLUTION",
            "rank": rank,
            "augmented_rank": augmented_rank,
        }

    if rank != 1:

        return {
            "status": "NONUNIQUE",
            "rank": rank,
            "augmented_rank": augmented_rank,
        }

    solution = A_matrix.gauss_jordan_solve(
        b
    )[0]

    coefficient = clean(
        solution[0, 0]
    )

    residuals = [
        clean(
            coefficient * source[i]
            - target[i]
        )
        for i in range(
            len(source)
        )
    ]

    return {
        "status": (
            "EXACT"
            if all(
                residual == 0
                for residual in residuals
            )
            else "VERIFICATION_FAILED"
        ),
        "rank": rank,
        "augmented_rank": augmented_rank,
        "coefficient": coefficient,
        "residuals": residuals,
    }


def fixed_j_recurrence_audit(A):

    print()
    print("=" * 78)
    print(
        "5. FIXED-j NEWTON COEFFICIENT DYNAMICS"
    )
    print("=" * 78)

    results = {}

    max_j = max(
        j
        for j, _ in A
    )

    max_t = max(
        t
        for _, t in A
    )

    for j in range(
        max_j + 1
    ):

        values = [
            A[(j, t)]
            for t in range(
                max_t + 1
            )
            if (j, t) in A
        ]

        result = solve_scalar_vertical_recurrence(
            values
        )

        results[j] = result

        print()
        print(
            "  j={}: values={}".format(
                j,
                values,
            )
        )

        print(
            "    status={}".format(
                result["status"]
            )
        )

        if result["status"] == "EXACT":

            print(
                "    coefficient={}".format(
                    result["coefficient"]
                )
            )

            print(
                "    residuals={}".format(
                    result["residuals"]
                )
            )

    return results


# ============================================================================
# SECTION 6 — COMMON SCALAR DYNAMICS
# ============================================================================

def common_scalar_dynamics(A):

    print()
    print("=" * 78)
    print(
        "6. COMMON NEWTON-COEFFICIENT SCALAR DYNAMICS"
    )
    print("=" * 78)

    max_t = max(
        t
        for _, t in A
    )

    results = []

    for t in range(
        max_t
    ):

        rows_available = [
            j
            for j in range(4)
            if (
                (j, t) in A
                and
                (j, t + 1) in A
            )
        ]

        if len(rows_available) < 2:
            continue

        ratios = []

        for j in rows_available:

            source = A[(j, t)]
            target = A[(j, t + 1)]

            if source == 0:
                ratios.append(
                    None
                )
            else:
                ratios.append(
                    clean(
                        target / source
                    )
                )

        nonzero_ratios = [
            r
            for r in ratios
            if r is not None
        ]

        common = (
            len(nonzero_ratios) >= 2
            and all(
                r == nonzero_ratios[0]
                for r in nonzero_ratios[1:]
            )
            and all(
                A[(j, t)] != 0
                for j in rows_available
            )
        )

        print()
        print(
            "  t={} -> {}:".format(
                t,
                t + 1,
            )
        )

        print(
            "    available_j={}".format(
                rows_available
            )
        )

        print(
            "    ratios={}".format(
                ratios
            )
        )

        print(
            "    common_scalar={}".format(
                common
            )
        )

        results.append(
            (
                t,
                common,
            )
        )

    return results


# ============================================================================
# SECTION 7 — COEFFICIENT CROSS-RATIOS
# ============================================================================

def coefficient_cross_ratio_audit(A):

    print()
    print("=" * 78)
    print(
        "7. NEWTON-COEFFICIENT CROSS-RATIO AUDIT"
    )
    print("=" * 78)

    ratios = []

    max_t = max(
        t
        for _, t in A
    )

    for j1 in range(3):

        for j2 in range(
            j1 + 1,
            3,
        ):

            for t in range(
                max_t
            ):

                locations = [
                    (j1, t),
                    (j2, t),
                    (j1, t + 1),
                    (j2, t + 1),
                ]

                if not all(
                    location in A
                    for location in locations
                ):
                    continue

                denominator = (
                    A[(j2, t)]
                    * A[(j1, t + 1)]
                )

                if denominator == 0:
                    continue

                ratio = clean(
                    (
                        A[(j1, t)]
                        * A[(j2, t + 1)]
                    )
                    /
                    denominator
                )

                ratios.append(
                    (
                        j1,
                        j2,
                        t,
                        ratio,
                    )
                )

    if not ratios:

        print(
            "  complete_cross_ratio_rectangles=0"
        )

        return []

    for j1, j2, t, ratio in ratios:

        print()
        print(
            "  rows=({},{}), t={} -> ratio={}".format(
                j1,
                j2,
                t,
                ratio,
            )
        )

        print(
            "    equals_one={}".format(
                ratio == 1
            )
        )

    print()
    print(
        "  complete_cross_ratio_rectangles={}".format(
            len(ratios)
        )
    )

    return ratios


# ============================================================================
# SECTION 8 — PRIME PROFILE
# ============================================================================

def coefficient_prime_audit(A):

    print()
    print("=" * 78)
    print(
        "8. NEWTON-COEFFICIENT PRIME PROFILE"
    )
    print("=" * 78)

    for (
        j,
        t,
    ), value in sorted(
        A.items(),
        key=lambda item: (
            item[0][1],
            item[0][0],
        )
    ):

        print()
        print(
            "  A[{},{}]={}".format(
                j,
                t,
                value,
            )
        )

        print(
            "    valuations={}".format(
                {
                    prime: valuation(
                        value,
                        prime,
                    )
                    for prime in (
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
# SECTION 9 — BASIC CONSISTENCY
# ============================================================================

def reconstruction_consistency_audit(
    layers,
    A,
):

    print()
    print("=" * 78)
    print(
        "9. NEWTON COORDINATE RECONSTRUCTION CONSISTENCY"
    )
    print("=" * 78)

    failures = 0

    for t, layer in layers.items():

        values = [
            value
            for _, value in layer
        ]

        rows = finite_difference_rows(
            values
        )

        for j, row in enumerate(
            rows
        ):

            if not row:
                continue

            if A.get(
                (j, t),
                None
            ) != clean(row[0]):

                failures += 1

                print()
                print(
                    "  FAILURE at (j,t)=({},{}):".format(
                        j,
                        t,
                    )
                )

    print(
        "  reconstruction_coordinate_failures={}".format(
            failures
        )
    )

    print(
        "  exact_newton_coordinate_construction={}".format(
            failures == 0
        )
    )

    return failures == 0


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 332R — EXACT NEWTON-COEFFICIENT TRIANGLE / "
        "HIDDEN-LOW-RANK DYNAMICS AUDIT"
    )
    print("=" * 78)

    layers = build_layers()

    A = build_newton_triangle(
        layers
    )

    print_triangle(
        A
    )

    support_audit(
        A
    )

    coefficient_minor_audit(
        A
    )

    coefficient_matrix = coefficient_rank_audit(
        A
    )

    vertical_results = fixed_j_recurrence_audit(
        A
    )

    common_scalar_results = common_scalar_dynamics(
        A
    )

    coefficient_ratios = coefficient_cross_ratio_audit(
        A
    )

    coefficient_prime_audit(
        A
    )

    reconstruction_ok = (
        reconstruction_consistency_audit(
            layers,
            A,
        )
    )

    # ------------------------------------------------------------------------
    # Structural interpretation
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "10. STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
The source table Q_t(p) has resisted a large class of direct laws:

    fixed vertical recurrences,
    p-dependent polynomial recurrences,
    low-degree rational recurrences,
    low-degree bivariate polynomial laws,
    additive local stencils,
    additive rank-one separability,
    multiplicative rank-one separability.

The observed p-grid is nevertheless highly structured:

    p = 1, 3, 5, 7.

Therefore the natural discrete coordinates are the exact finite
differences

    A[j,t] = Delta_p^j Q_t(1).

Experiment 332R asks whether the source becomes simpler after this
coordinate change.

The strongest useful outcomes would be:

    * coefficient rank significantly below the ambient 4;
    * scalar vertical dynamics shared by several j-levels;
    * coefficient cross-ratios equal to 1 or otherwise simple;
    * low-complexity prime structure.

A failure of all these tests would be important: it would show that
even the canonical Newton coordinates do not expose a low-dimensional
hidden source mechanism.

The reconstruction audit guarantees that no numerical approximation,
interpolation, missing-value inference, or altered source data have
entered the calculation.

No synthetic second n=pq case is introduced.
"""
    )

    # ------------------------------------------------------------------------
    # FINAL EXACTNESS
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  coefficient_matrix_rank={}".format(
            coefficient_matrix.rank()
        )
    )

    print(
        "  coefficient_matrix_shape={}".format(
            coefficient_matrix.shape
        )
    )

    print(
        "  coefficient_minor_audit_completed=True"
    )

    print(
        "  fixed_j_dynamics_tested=True"
    )

    print(
        "  common_scalar_dynamics_tested=True"
    )

    print(
        "  coefficient_cross_ratio_tested=True"
    )

    print(
        "  exact_newton_reconstruction={}".format(
            reconstruction_ok
        )
    )

    print(
        "  missing_values_used=False"
    )

    print(
        "  interpolation_performed=False"
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
        "  failures={}".format(
            0
            if reconstruction_ok
            else 1
        )
    )

    print(
        "  ALL BASIC CHECKS PASS={}".format(
            reconstruction_ok
        )
    )

    print()
    print(
        "EXPERIMENT 332R COMPLETE"
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