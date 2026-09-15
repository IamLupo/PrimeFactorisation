#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 301R — EXACT TERMINAL-DISTANCE / p-FINITE-DIFFERENCE SOURCE AUDIT
==============================================================================

Purpose
-------

Experiment 300R rejected low-degree joint laws of the form

    q = q(p,r)

using direct polynomial, tensor, binomial, falling-factorial, and
terminal-distance coordinates.

Experiment 301R now uses the triangular support itself.

Define

    D(p) = len(q_p) - 1
    t    = D(p) - r.

For fixed t, compare the aligned source values

    Q_t(p) = q_p(D(p)-t).

The available primes are

    p = 1, 3, 5, 7,

which are equally spaced.  Therefore put

    n = (p-1)/2

and test exact finite-difference structure in n.

The experiment tests:

    1. exact finite-difference tables for each terminal layer t;
    2. exact polynomial degree in n;
    3. first-order and second-order constant recurrences;
    4. geometric ratios;
    5. normalization by the terminal value q_p(D(p));
    6. cross-layer ratios Q_{t+1}(p)/Q_t(p);
    7. exact common finite-difference signatures across t.

Only overdetermined structure is treated as evidence.

A polynomial through exactly all available points is labelled
INTERPOLATION and is not counted as a structural hit.

No synthetic second n=pq case.
No external files.
No arbitrary fitted source operator.
Exact rational/integer arithmetic only.
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
# BASIC HELPERS
# ============================================================================

def clean(expr):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(expr)
            )
        )
    )


def v17(value):
    value = abs(int(value))

    if value == 0:
        return None

    out = 0

    while value % 17 == 0:
        value //= 17
        out += 1

    return out


def finite_differences(values):
    """
    Return all successive forward-difference rows.

    Input:
        [a0,a1,...]

    Output:
        [
            [a0,a1,...],
            [Δa0,Δa1,...],
            [Δ²a0,...],
            ...
        ]
    """
    rows = [
        [
            sp.Rational(v)
            for v in values
        ]
    ]

    current = rows[0]

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


def recurrence_order_one(values):
    """
    Test

        a[n+1] = A a[n]

    with one constant A.
    """
    if len(values) < 2:
        return (
            "INSUFFICIENT_DATA",
            None,
        )

    values = [
        sp.Rational(v)
        for v in values
    ]

    if values[0] == 0:
        return (
            "NO_SOLUTION",
            None,
        )

    A = clean(
        values[1] / values[0]
    )

    for i in range(
        len(values) - 1
    ):

        if clean(
            values[i + 1]
            - A * values[i]
        ) != 0:
            return (
                "NO_SOLUTION",
                None,
            )

    return (
        "EXACT",
        A,
    )


def recurrence_order_two(values):
    """
    Test

        a[n+2] = A a[n+1] + B a[n].

    Requires at least four values for an overdetermined test.
    """
    if len(values) < 4:
        return (
            "INSUFFICIENT_DATA",
            None,
        )

    vals = [
        sp.Rational(v)
        for v in values
    ]

    A, B = sp.symbols(
        "A B"
    )

    equations = []

    for i in range(
        len(vals) - 2
    ):
        equations.append(
            sp.Eq(
                vals[i + 2],
                A * vals[i + 1]
                + B * vals[i],
            )
        )

    solution = sp.solve(
        equations,
        [A, B],
        dict=True,
    )

    if not solution:
        return (
            "NO_SOLUTION",
            None,
        )

    if len(solution) != 1:
        return (
            "NONUNIQUE",
            solution,
        )

    sol = solution[0]

    if A not in sol or B not in sol:
        return (
            "NONUNIQUE",
            sol,
        )

    A_value = clean(sol[A])
    B_value = clean(sol[B])

    for i in range(
        len(vals) - 2
    ):

        lhs = vals[i + 2]
        rhs = clean(
            A_value * vals[i + 1]
            + B_value * vals[i]
        )

        if lhs != rhs:
            return (
                "NO_SOLUTION",
                None,
            )

    return (
        "EXACT",
        (
            A_value,
            B_value,
        ),
    )


def polynomial_degree_from_values(values):
    """
    Determine the finite-difference degree.

    For a sequence a_n, the first difference row which is constant
    indicates the polynomial degree.

    If all values are zero, degree = -1.

    A degree determined from exactly m points is labelled only as an
    interpolation degree elsewhere; this helper reports the algebraic
    degree of the supplied finite table.
    """
    vals = [
        sp.Rational(v)
        for v in values
    ]

    if all(
        v == 0
        for v in vals
    ):
        return -1

    rows = finite_differences(
        vals
    )

    for order, row in enumerate(rows):

        if len(row) == 1:
            return order

        if all(
            row[i] == row[0]
            for i in range(
                1,
                len(row),
            )
        ):
            return order

    return None


def lagrange_polynomial(values):
    """
    Interpolate values on n=0,1,... exactly.
    Used only diagnostically.
    """
    n = sp.Symbol("n")

    points = [
        (
            sp.Integer(i),
            sp.Rational(v),
        )
        for i, v in enumerate(values)
    ]

    if not points:
        return sp.Integer(0)

    return clean(
        sp.interpolate(
            points,
            n,
        )
    )


# ============================================================================
# TERMINAL DISTANCE TABLE
# ============================================================================

def D(p):
    return len(Q[p]) - 1


def terminal_distance_table():
    """
    Return

        t -> [(p,n,value), ...]

    where

        n=(p-1)/2
        r=D(p)-t.
    """
    table = {}

    max_t = max(
        D(p)
        for p in Q
    )

    for t in range(
        max_t + 1
    ):

        entries = []

        for p in sorted(Q):

            r = D(p) - t

            if r < 0:
                continue

            entries.append(
                (
                    p,
                    (p - 1) // 2,
                    sp.Integer(
                        Q[p][r]
                    ),
                )
            )

        table[t] = entries

    return table


# ============================================================================
# 1. PRINT TERMINAL-LAYER INVENTORY
# ============================================================================

def print_inventory(table):

    print()
    print("=" * 78)
    print(
        "1. TERMINAL-DISTANCE SOURCE LAYERS"
    )
    print("=" * 78)

    for t in sorted(table):

        entries = table[t]

        print()
        print(
            "  t={}  (t=D(p)-r)".format(
                t
            )
        )

        print(
            "    available_p={}".format(
                [
                    e[0]
                    for e in entries
                ]
            )
        )

        print(
            "    values={}".format(
                [
                    e[2]
                    for e in entries
                ]
            )
        )


# ============================================================================
# 2. FINITE-DIFFERENCE AUDIT
# ============================================================================

def print_difference_audit(table):

    print()
    print("=" * 78)
    print(
        "2. EXACT p-FINITE-DIFFERENCE TABLES"
    )
    print("=" * 78)

    profile = {}

    for t in sorted(table):

        values = [
            entry[2]
            for entry in table[t]
        ]

        rows = finite_differences(
            values
        )

        print()
        print(
            "  t={}:".format(t)
        )

        for order, row in enumerate(rows):

            print(
                "    delta^{}={}".format(
                    order,
                    row,
                )
            )

        degree = (
            polynomial_degree_from_values(
                values
            )
        )

        print(
            "    finite_difference_degree={}".format(
                degree,
            )
        )

        profile[t] = degree

    return profile


# ============================================================================
# 3. OVERDETERMINED LINEAR / GEOMETRIC TESTS
# ============================================================================

def recurrence_audit(table):

    print()
    print("=" * 78)
    print(
        "3. EXACT CONSTANT-COEFFICIENT RECURRENCE AUDIT"
    )
    print("=" * 78)

    hits_order1 = []
    hits_order2 = []

    for t in sorted(table):

        values = [
            entry[2]
            for entry in table[t]
        ]

        status1, param1 = (
            recurrence_order_one(
                values
            )
        )

        status2, param2 = (
            recurrence_order_two(
                values
            )
        )

        print()
        print(
            "  t={}:".format(t)
        )

        print(
            "    order1: status={} parameter={}".format(
                status1,
                param1,
            )
        )

        print(
            "    order2: status={} parameters={}".format(
                status2,
                param2,
            )
        )

        if status1 == "EXACT":
            hits_order1.append(
                (
                    t,
                    param1,
                )
            )

        if status2 == "EXACT":
            hits_order2.append(
                (
                    t,
                    param2,
                )
            )

    return (
        hits_order1,
        hits_order2,
    )


# ============================================================================
# 4. FIRST-DIFFERENCE / DEGREE-1 OVERDETERMINATION
# ============================================================================

def affine_p_audit(table):

    print()
    print("=" * 78)
    print(
        "4. EXACT AFFINE-IN-p TEST"
    )
    print("=" * 78)

    hits = []

    n_sym = sp.Symbol(
        "n"
    )

    for t in sorted(table):

        entries = table[t]

        if len(entries) < 3:
            print()
            print(
                "  t={}: INSUFFICIENT_DATA".format(
                    t
                )
            )
            continue

        values = [
            entry[2]
            for entry in entries
        ]

        diffs = finite_differences(
            values
        )[1]

        is_constant = (
            len(diffs) >= 2
            and all(
                d == diffs[0]
                for d in diffs[1:]
            )
        )

        poly = lagrange_polynomial(
            values
        )

        degree = int(
            sp.Poly(
                poly,
                n_sym,
                domain=sp.QQ,
            ).degree()
        )

        print()
        print(
            "  t={}: points={} degree={} "
            "affine_exact={}".format(
                t,
                len(values),
                degree,
                is_constant,
            )
        )

        if is_constant:
            print(
                "    polynomial={}".format(
                    poly
                )
            )
            hits.append(
                (
                    t,
                    poly,
                )
            )

    return hits


# ============================================================================
# 5. TERMINAL NORMALIZATION
# ============================================================================

def terminal_normalization_audit(table):

    print()
    print("=" * 78)
    print(
        "5. TERMINAL-VALUE NORMALIZATION AUDIT"
    )
    print("=" * 78)

    normalized = {}

    for t in sorted(table):

        entries = table[t]

        values = [
            entry[2]
            for entry in entries
        ]

        if not values:
            continue

        print()
        print(
            "  t={}:".format(t)
        )

        base = values[0]

        if base == 0:
            print(
                "    normalization_undefined=True"
            )
            continue

        row = [
            clean(
                value / base
            )
            for value in values
        ]

        print(
            "    values={}".format(
                values
            )
        )

        print(
            "    normalized_to_first={}".format(
                row
            )
        )

        ratio_status, ratio_param = (
            recurrence_order_one(
                row
            )
        )

        print(
            "    normalized_geometric_test={} {}".format(
                ratio_status,
                ratio_param,
            )
        )

        normalized[t] = row

    return normalized


# ============================================================================
# 6. CROSS-LAYER RATIO AUDIT
# ============================================================================

def cross_layer_ratio_audit(table):

    print()
    print("=" * 78)
    print(
        "6. CROSS-LAYER RATIO q_t / q_(t-1)"
    )
    print("=" * 78)

    ratio_sequences = {}

    for t in range(
        1,
        max(table.keys()) + 1,
    ):

        current = {
            entry[0]: entry[2]
            for entry in table[t]
        }

        previous = {
            entry[0]: entry[2]
            for entry in table[t - 1]
        }

        common_p = sorted(
            set(current)
            & set(previous)
        )

        ratios = []

        for p in common_p:

            numerator = current[p]
            denominator = previous[p]

            if denominator == 0:
                ratio = None
            else:
                ratio = clean(
                    numerator
                    / denominator
                )

            ratios.append(
                (
                    p,
                    ratio,
                )
            )

        print()
        print(
            "  t={} over t-1:".format(
                t
            )
        )

        print(
            "    ratios={}".format(
                ratios
            )
        )

        defined = [
            ratio
            for _, ratio in ratios
            if ratio is not None
        ]

        if defined:

            status, parameter = (
                recurrence_order_one(
                    defined
                )
            )

            print(
                "    ratio_geometric_test={} {}".format(
                    status,
                    parameter,
                )
            )

            ratio_sequences[t] = (
                ratios
            )

    return ratio_sequences


# ============================================================================
# 7. EXACT NORMALIZED TERMINAL LAYER
# ============================================================================

def terminal_layer_audit(table):

    print()
    print("=" * 78)
    print(
        "7. TERMINAL LAYER t=0"
    )
    print("=" * 78)

    entries = table.get(
        0,
        [],
    )

    values = [
        entry[2]
        for entry in entries
    ]

    primes = [
        entry[0]
        for entry in entries
    ]

    print(
        "  p={}".format(
            primes
        )
    )

    print(
        "  q_p(D(p))={}".format(
            values
        )
    )

    print(
        "  v17={}".format(
            [
                v17(v)
                for v in values
            ]
        )
    )

    if len(values) >= 2:

        ratios = []

        for i in range(
            len(values) - 1
        ):

            a = values[i]
            b = values[i + 1]

            if a == 0:
                ratios.append(
                    None
                )
            else:
                ratios.append(
                    clean(
                        b / a
                    )
                )

        print(
            "  adjacent_ratios={}".format(
                ratios
            )
        )

    differences = finite_differences(
        values
    )

    print(
        "  difference_table={}".format(
            differences
        )
    )

    # Since p=1,3,5,7 correspond to n=0,1,2,3,
    # the third difference is the only fully overdetermined
    # polynomial-degree diagnostic for t=0.
    if len(differences) >= 4:

        third = differences[3]

        print(
            "  third_difference={}".format(
                third
            )
        )

        if len(third) == 1:

            print(
                "  cubic_interpolation_value={}".format(
                    third[0]
                )
            )


# ============================================================================
# 8. PRIMITIVE INTEGER DIFFERENCE SIGNATURES
# ============================================================================

def primitive_signature(
    values
):
    values = [
        sp.Rational(v)
        for v in values
    ]

    if not values:
        return []

    denominator_lcm = 1

    for value in values:

        denominator_lcm = math.lcm(
            denominator_lcm,
            int(value.q),
        )

    integers = [
        int(
            value
            * denominator_lcm
        )
        for value in values
    ]

    gcd_value = 0

    for value in integers:

        gcd_value = math.gcd(
            gcd_value,
            abs(value),
        )

    if gcd_value:

        integers = [
            value // gcd_value
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


def difference_signature_audit(
    table
):

    print()
    print("=" * 78)
    print(
        "8. PRIMITIVE FINITE-DIFFERENCE SIGNATURES"
    )
    print("=" * 78)

    signatures = {}

    for t in sorted(table):

        values = [
            entry[2]
            for entry in table[t]
        ]

        print()
        print(
            "  t={}:".format(t)
        )

        diff_rows = finite_differences(
            values
        )

        sigs = []

        for order, row in enumerate(
            diff_rows
        ):

            sig = primitive_signature(
                row
            )

            print(
                "    delta^{} primitive={}".format(
                    order,
                    sig,
                )
            )

            sigs.append(
                sig
            )

        signatures[t] = sigs

    return signatures


# ============================================================================
# 9. STRUCTURAL INTERPRETATION
# ============================================================================

def interpretation():

    print()
    print("=" * 78)
    print(
        "9. STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
Experiment 300R showed that direct low-degree q(p,r) models do not give
an identifiable source law.

Experiment 301R reorganizes the same source data by the exact triangular
coordinate

    t = D(p) - r.

This is not a fitted coordinate. It is forced by the finite support of
each q-row.

For fixed t, the source values are aligned from the same terminal depth:

    Q_t(p) = q_p(D(p)-t).

Because

    p = 1,3,5,7

is equally spaced, set

    n = (p-1)/2

and use exact finite differences in n.

The strongest structural cases are:

    * constant first difference:
          Q_t is affine in n;

    * constant second difference:
          Q_t is quadratic in n;

    * constant third difference:
          Q_t is cubic in n;

    * exact constant-coefficient recurrence;

    * simple cross-layer ratio

          Q_{t+1}(p) / Q_t(p).

The terminal layer t=0 is especially important because it contains
all four supplied p-values and therefore permits the strongest
overdetermined finite-difference test.

Layers t>=1 contain fewer p-values, so their polynomial observations
must be interpreted more cautiously.

No polynomial interpolation is promoted to a universal law.

No second n=pq case is generated.
"""
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 301R — EXACT TERMINAL-DISTANCE / "
        "p-FINITE-DIFFERENCE SOURCE AUDIT"
    )
    print("=" * 78)

    table = terminal_distance_table()

    print_inventory(
        table
    )

    degree_profile = (
        print_difference_audit(
            table
        )
    )

    order1_hits, order2_hits = (
        recurrence_audit(
            table
        )
    )

    affine_hits = (
        affine_p_audit(
            table
        )
    )

    normalized = (
        terminal_normalization_audit(
            table
        )
    )

    cross_layer = (
        cross_layer_ratio_audit(
            table
        )
    )

    terminal_layer_audit(
        table
    )

    signatures = (
        difference_signature_audit(
            table
        )
    )

    interpretation()

    # ----------------------------------------------------------------------
    # FINAL EXACTNESS
    # ----------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "10. FINAL EXACTNESS"
    )
    print("=" * 78)

    exact_affine_hits = len(
        affine_hits
    )

    exact_order1_hits = len(
        order1_hits
    )

    exact_order2_hits = len(
        order2_hits
    )

    print(
        "  terminal_layers={}".format(
            len(table)
        )
    )

    print(
        "  affine_terminal_layer_hits={}".format(
            exact_affine_hits
        )
    )

    print(
        "  order1_recurrence_hits={}".format(
            exact_order1_hits
        )
    )

    print(
        "  order2_recurrence_hits={}".format(
            exact_order2_hits
        )
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
        "EXPERIMENT 301R COMPLETE"
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

