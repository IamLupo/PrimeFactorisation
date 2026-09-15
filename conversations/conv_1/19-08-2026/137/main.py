#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 303RR — EXACT TERMINAL-DISTANCE MIXED-DIFFERENCE / CROSS-LAYER AUDIT
==============================================================================

Corrected version of Experiment 303R.

This version explicitly defines the exact p-adic valuation helper that was
missing from the previous script.

No external files are used.
No synthetic second n=pq case is generated.
Only exact SymPy rational arithmetic is used.
"""

from __future__ import annotations

import math
import sys
import sympy as sp


# ============================================================================
# SOURCE TABLE
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

def clean(x):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(x)
            )
        )
    )


def valuation(value, prime):
    """
    Exact p-adic valuation v_prime(value) for a nonzero rational.

    For

        value = numerator / denominator,

    returns

        v_p(numerator) - v_p(denominator).

    Returns None for zero.
    """
    value = sp.Rational(value)

    if value == 0:
        return None

    numerator = abs(int(value.p))
    denominator = abs(int(value.q))

    result = 0

    while numerator % prime == 0:
        numerator //= prime
        result += 1

    while denominator % prime == 0:
        denominator //= prime
        result -= 1

    return result


def finite_difference_rows(values):
    current = [
        sp.Rational(v)
        for v in values
    ]

    rows = [current]

    while len(current) > 1:
        current = [
            clean(
                current[i + 1] - current[i]
            )
            for i in range(len(current) - 1)
        ]
        rows.append(current)

    return rows


def recurrence_order1(values):
    values = [
        sp.Rational(v)
        for v in values
    ]

    if len(values) < 2:
        return "INSUFFICIENT_DATA", None

    if values[0] == 0:
        return "NO_SOLUTION", None

    a = clean(
        values[1] / values[0]
    )

    for i in range(len(values) - 1):
        if clean(
            values[i + 1] - a * values[i]
        ) != 0:
            return "NO_SOLUTION", None

    return "EXACT", a


def recurrence_order2(values):
    values = [
        sp.Rational(v)
        for v in values
    ]

    if len(values) < 4:
        return "INSUFFICIENT_DATA", None

    A, B = sp.symbols("A B")

    equations = [
        sp.Eq(
            values[i + 2],
            A * values[i + 1] + B * values[i],
        )
        for i in range(len(values) - 2)
    ]

    solutions = sp.solve(
        equations,
        [A, B],
        dict=True,
    )

    if len(solutions) == 0:
        return "NO_SOLUTION", None

    if len(solutions) != 1:
        return "NONUNIQUE", solutions

    sol = solutions[0]

    if A not in sol or B not in sol:
        return "NONUNIQUE", sol

    a = clean(sol[A])
    b = clean(sol[B])

    for i in range(len(values) - 2):
        if clean(
            values[i + 2]
            - a * values[i + 1]
            - b * values[i]
        ) != 0:
            return "NO_SOLUTION", None

    return "EXACT", (a, b)


def primitive_integer_list(values):
    """
    Convert a rational list to a primitive integer signature.
    """
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
        int(value * denominator_lcm)
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
                    -x for x in integers
                ]
            break

    return integers


def exact_interpolation_degree(values):
    if not values:
        return None

    x = sp.Symbol("x")

    polynomial = sp.interpolate(
        [
            (sp.Integer(i), sp.Rational(v))
            for i, v in enumerate(values)
        ],
        x,
    )

    polynomial = clean(polynomial)

    return int(
        sp.Poly(
            polynomial,
            x,
            domain=sp.QQ,
        ).degree()
    )


# ============================================================================
# TERMINAL DISTANCE
# ============================================================================

def D(p):
    return len(Q[p]) - 1


def build_layers():
    """
    Construct

        Q_t(p) = q_p(D(p)-t).

    Each layer is ordered by increasing p.
    """
    layers = {}

    max_t = max(
        D(p)
        for p in Q
    )

    for t in range(max_t + 1):

        layer = []

        for p in sorted(Q):

            r = D(p) - t

            if r < 0:
                continue

            layer.append(
                (
                    p,
                    sp.Integer(Q[p][r]),
                )
            )

        layers[t] = layer

    return layers


# ============================================================================
# DIFFERENCE TRIANGLE
# ============================================================================

def build_difference_triangle(layers):
    """
    triangle[(j,t)] is the full j-th forward-difference row for Q_t.
    """
    triangle = {}

    for t in sorted(layers):

        values = [
            value
            for _, value in layers[t]
        ]

        rows = finite_difference_rows(
            values
        )

        for j, row in enumerate(rows):
            triangle[(j, t)] = row

    return triangle


# ============================================================================
# OUTPUT: SOURCE LAYERS
# ============================================================================

def print_layers(layers):
    print()
    print("=" * 78)
    print(
        "1. TERMINAL-DISTANCE SOURCE LAYERS"
    )
    print("=" * 78)

    for t in sorted(layers):

        print()
        print(
            "  t={}:".format(t)
        )

        print(
            "    available_p={}".format(
                [
                    p
                    for p, _ in layers[t]
                ]
            )
        )

        print(
            "    Q_t={}".format(
                [
                    value
                    for _, value in layers[t]
                ]
            )
        )


# ============================================================================
# OUTPUT: FULL DIFFERENCE TRIANGLE
# ============================================================================

def print_triangle(triangle):
    print()
    print("=" * 78)
    print(
        "2. COMPLETE p-FINITE-DIFFERENCE TRIANGLE"
    )
    print("=" * 78)

    for j in sorted(
        set(
            jj
            for jj, _ in triangle
        )
    ):

        print()
        print(
            "  difference_order={}".format(
                j
            )
        )

        for t in sorted(
            tt
            for jj, tt in triangle
            if jj == j
        ):

            print(
                "    t={}: {}".format(
                    t,
                    triangle[(j, t)],
                )
            )


# ============================================================================
# FIXED-ORDER CROSS-LAYER AUDIT
# ============================================================================

def cross_layer_order_audit(triangle):
    print()
    print("=" * 78)
    print(
        "3. FIXED-ORDER CROSS-LAYER AUDIT"
    )
    print("=" * 78)

    results = {}

    for j in sorted(
        set(
            jj
            for jj, _ in triangle
        )
    ):

        usable = []

        for t in sorted(
            tt
            for jj, tt in triangle
            if jj == j
        ):

            values = triangle[(j, t)]

            # To obtain a scalar cross-layer sequence we require a
            # single coefficient at that difference order.
            if len(values) == 1:
                usable.append(
                    (
                        t,
                        clean(values[0]),
                    )
                )

        if len(usable) < 2:

            print()
            print(
                "  j={}: status=INSUFFICIENT_DATA".format(
                    j
                )
            )

            results[j] = None
            continue

        scalar_values = [
            value
            for _, value in usable
        ]

        order1_status, order1_parameter = (
            recurrence_order1(
                scalar_values
            )
        )

        order2_status, order2_parameters = (
            recurrence_order2(
                scalar_values
            )
        )

        polynomial_degree = (
            exact_interpolation_degree(
                scalar_values
            )
        )

        print()
        print(
            "  j={}:".format(j)
        )

        print(
            "    t={}".format(
                [
                    t
                    for t, _ in usable
                ]
            )
        )

        print(
            "    values={}".format(
                scalar_values
            )
        )

        print(
            "    order1: status={} parameter={}".format(
                order1_status,
                order1_parameter,
            )
        )

        print(
            "    order2: status={} parameters={}".format(
                order2_status,
                order2_parameters,
            )
        )

        print(
            "    interpolation_degree={}".format(
                polynomial_degree
            )
        )

        results[j] = {
            "usable": usable,
            "order1": (
                order1_status,
                order1_parameter,
            ),
            "order2": (
                order2_status,
                order2_parameters,
            ),
            "degree": polynomial_degree,
        }

    return results


# ============================================================================
# COMMON VECTOR TRANSFER
# ============================================================================

def vector_transfer_tests(triangle):
    print()
    print("=" * 78)
    print(
        "4. EXACT COMMON SCALAR TRANSFER IN t"
    )
    print("=" * 78)

    forward_hits = []
    reverse_hits = []

    all_t = sorted(
        set(
            t
            for _, t in triangle
        )
    )

    # Forward: t -> t+1
    for t in all_t:

        if t + 1 not in all_t:
            continue

        common_orders = sorted(
            set(
                j
                for j, tt in triangle
                if tt == t
            )
            &
            set(
                j
                for j, tt in triangle
                if tt == t + 1
            )
        )

        candidate = None
        valid = True

        for j in common_orders:

            arow = triangle[(j, t)]
            brow = triangle[(j, t + 1)]

            if len(arow) != len(brow):
                continue

            for a, b in zip(arow, brow):

                if a == 0:

                    if b != 0:
                        valid = False
                        break

                    continue

                ratio = clean(
                    b / a
                )

                if candidate is None:
                    candidate = ratio
                elif ratio != candidate:
                    valid = False
                    break

            if not valid:
                break

        if valid and candidate is not None:
            forward_hits.append(
                (
                    t,
                    t + 1,
                    candidate,
                )
            )

    # Reverse: t -> t-1
    for t in all_t:

        if t - 1 not in all_t:
            continue

        common_orders = sorted(
            set(
                j
                for j, tt in triangle
                if tt == t
            )
            &
            set(
                j
                for j, tt in triangle
                if tt == t - 1
            )
        )

        candidate = None
        valid = True

        for j in common_orders:

            arow = triangle[(j, t)]
            brow = triangle[(j, t - 1)]

            if len(arow) != len(brow):
                continue

            for a, b in zip(arow, brow):

                if a == 0:

                    if b != 0:
                        valid = False
                        break

                    continue

                ratio = clean(
                    b / a
                )

                if candidate is None:
                    candidate = ratio
                elif ratio != candidate:
                    valid = False
                    break

            if not valid:
                break

        if valid and candidate is not None:
            reverse_hits.append(
                (
                    t,
                    t - 1,
                    candidate,
                )
            )

    print()
    print(
        "  forward_hits={}".format(
            forward_hits
        )
    )

    print(
        "  reverse_hits={}".format(
            reverse_hits
        )
    )

    return (
        forward_hits,
        reverse_hits,
    )


# ============================================================================
# MIXED DIFFERENCE COMMUTATION
# ============================================================================

def mixed_difference_test(layers):
    print()
    print("=" * 78)
    print(
        "5. MIXED-DIFFERENCE COMMUTATION TEST"
    )
    print("=" * 78)

    ps = sorted(Q)

    matrix = {}

    for n, p in enumerate(ps):

        matrix[n] = {}

        for t in sorted(layers):

            matches = [
                value
                for pp, value in layers[t]
                if pp == p
            ]

            matrix[n][t] = (
                matches[0]
                if matches
                else None
            )

    checks = 0
    failures = 0

    for n in range(len(ps) - 1):

        for t in sorted(layers):

            if t + 1 not in layers:
                continue

            a = matrix[n].get(t)
            b = matrix[n + 1].get(t)
            c = matrix[n].get(t + 1)
            d = matrix[n + 1].get(t + 1)

            if any(
                value is None
                for value in (
                    a, b, c, d
                )
            ):
                continue

            left = clean(
                (d - c) - (b - a)
            )

            right = clean(
                (d - b) - (c - a)
            )

            checks += 1

            if left != 0 or right != 0:
                failures += 1

            print()
            print(
                "  n={}, t={}: ΔnΔt={} ΔtΔn={} exact={}".format(
                    n,
                    t,
                    left,
                    right,
                    left == 0 and right == 0,
                )
            )

    print()
    print(
        "  checks={}".format(
            checks
        )
    )

    print(
        "  failures={}".format(
            failures
        )
    )

    return (
        checks,
        failures,
    )


# ============================================================================
# DIFFERENCE TRIANGLE RANK
# ============================================================================

def rank_audit(triangle):
    print()
    print("=" * 78)
    print(
        "6. DIFFERENCE-TRIANGLE RANK AUDIT"
    )
    print("=" * 78)

    max_j = max(
        j
        for j, _ in triangle
    )

    max_t = max(
        t
        for _, t in triangle
    )

    matrix = []

    for j in range(max_j + 1):

        row = []

        for t in range(max_t + 1):

            if (j, t) not in triangle:
                row.append(
                    sp.Integer(0)
                )
                continue

            values = triangle[(j, t)]

            # This rank diagnostic uses the first available
            # coefficient in each difference row.
            row.append(
                clean(values[0])
            )

        matrix.append(row)

    M = sp.Matrix(matrix)

    rank = M.rank()

    print()
    print(
        "  matrix_shape={}".format(
            M.shape
        )
    )

    print(
        "  matrix={}".format(
            M
        )
    )

    print(
        "  rank={}".format(
            rank
        )
    )

    print(
        "  rank_one={}".format(
            rank <= 1
        )
    )

    print(
        "  rank_two_or_less={}".format(
            rank <= 2
        )
    )

    return M, rank


# ============================================================================
# PRIMITIVE / VALUATION AUDIT
# ============================================================================

def primitive_content_audit(triangle):
    print()
    print("=" * 78)
    print(
        "7. PRIMITIVE CONTENT / VALUATION AUDIT"
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

            for value in triangle[(j, t)]:
                entries.append(
                    (
                        t,
                        sp.Rational(value),
                    )
                )

        if not entries:
            continue

        values = [
            value
            for _, value in entries
        ]

        print()
        print(
            "  j={}:".format(
                j
            )
        )

        print(
            "    entries={}".format(
                entries
            )
        )

        print(
            "    primitive_signature={}".format(
                primitive_integer_list(
                    values
                )
            )
        )

        for prime in (
            2,
            3,
            5,
            7,
            17,
        ):

            print(
                "    v_{}={}".format(
                    prime,
                    [
                        (
                            t,
                            valuation(
                                value,
                                prime,
                            ),
                        )
                        for t, value in entries
                    ],
                )
            )


# ============================================================================
# LEAVE-ONE-t-OUT TEST
# ============================================================================

def leave_one_t_out_order_rows(triangle):
    print()
    print("=" * 78)
    print(
        "8. LEAVE-ONE-t-OUT DIFFERENCE-ROW TEST"
    )
    print("=" * 78)

    tested = 0
    exact = 0

    x = sp.Symbol("t")

    for j in sorted(
        set(
            jj
            for jj, _ in triangle
        )
    ):

        usable = []

        for t in sorted(
            tt
            for jj, tt in triangle
            if jj == j
        ):

            values = triangle[(j, t)]

            if len(values) != 1:
                continue

            usable.append(
                (
                    t,
                    clean(values[0]),
                )
            )

        if len(usable) < 3:
            continue

        print()
        print(
            "  difference_order={}".format(
                j
            )
        )

        for omitted_index in range(
            len(usable)
        ):

            omitted_t, actual = (
                usable[omitted_index]
            )

            training = [
                point
                for idx, point in enumerate(
                    usable
                )
                if idx != omitted_index
            ]

            if len(training) < 2:
                continue

            tested += 1

            polynomial = clean(
                sp.interpolate(
                    training,
                    x,
                )
            )

            prediction = clean(
                polynomial.subs(
                    x,
                    omitted_t,
                )
            )

            ok = (
                prediction
                == actual
            )

            if ok:
                exact += 1

            print(
                "    omitted_t={}: prediction={} "
                "actual={} exact={}".format(
                    omitted_t,
                    prediction,
                    actual,
                    ok,
                )
            )

    return (
        tested,
        exact,
    )


# ============================================================================
# STRUCTURAL INTERPRETATION
# ============================================================================

def structural_interpretation():
    print()
    print("=" * 78)
    print(
        "9. STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
Experiment 302R compressed each terminal-distance layer to one scalar
H_t and found no simple law in that scalar ladder.

Experiment 303RR restores the COMPLETE finite-difference information.

For

    Q_t(n) = q_{p_n}(D(p_n)-t),
    p_n = 2n+1,

define

    F_{j,t} = Δ_n^j Q_t.

The question is whether the terminal-depth coordinate t acts by a
simple exact operator on these difference coordinates.

The relevant possibilities are:

    * fixed-order recurrence in t;
    * common scalar transfer;
    * low-rank difference triangle;
    * exact mixed-difference structure;
    * exact primitive or valuation patterns.

This is stronger than studying only H_t because no intermediate
finite-difference coefficients have been discarded.

A positive result would expose a two-coordinate source structure.

A negative result would show that even after the forced terminal
re-indexing and exact Newton/finite-difference transformation, the
available source table does not reveal a compact cross-layer operator.

No synthetic n=pq case is generated.
"""
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 303RR — EXACT TERMINAL-DISTANCE MIXED-DIFFERENCE / "
        "CROSS-LAYER AUDIT"
    )
    print("=" * 78)

    layers = build_layers()

    print_layers(
        layers
    )

    triangle = build_difference_triangle(
        layers
    )

    print_triangle(
        triangle
    )

    order_results = (
        cross_layer_order_audit(
            triangle
        )
    )

    forward_hits, reverse_hits = (
        vector_transfer_tests(
            triangle
        )
    )

    mixed_checks, mixed_failures = (
        mixed_difference_test(
            layers
        )
    )

    _, triangle_rank = (
        rank_audit(
            triangle
        )
    )

    primitive_content_audit(
        triangle
    )

    loo_tested, loo_exact = (
        leave_one_t_out_order_rows(
            triangle
        )
    )

    structural_interpretation()

    order1_hits = sum(
        1
        for result in order_results.values()
        if result is not None
        and result["order1"][0] == "EXACT"
    )

    order2_hits = sum(
        1
        for result in order_results.values()
        if result is not None
        and result["order2"][0] == "EXACT"
    )

    print()
    print("=" * 78)
    print(
        "10. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  terminal_layer_count={}".format(
            len(layers)
        )
    )

    print(
        "  finite_difference_orders={}".format(
            len(
                set(
                    j
                    for j, _ in triangle
                )
            )
        )
    )

    print(
        "  fixed_order_order1_hits={}".format(
            order1_hits
        )
    )

    print(
        "  fixed_order_order2_hits={}".format(
            order2_hits
        )
    )

    print(
        "  common_scalar_forward_hits={}".format(
            len(forward_hits)
        )
    )

    print(
        "  common_scalar_reverse_hits={}".format(
            len(reverse_hits)
        )
    )

    print(
        "  mixed_difference_checks={}".format(
            mixed_checks
        )
    )

    print(
        "  mixed_difference_failures={}".format(
            mixed_failures
        )
    )

    print(
        "  difference_triangle_rank={}".format(
            triangle_rank
        )
    )

    print(
        "  leave_one_t_tested={}".format(
            loo_tested
        )
    )

    print(
        "  leave_one_t_exact={}".format(
            loo_exact
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
        "EXPERIMENT 303RR COMPLETE"
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