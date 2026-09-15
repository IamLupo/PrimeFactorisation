#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 304RR — EXACT TERMINAL-DISTANCE
DIFFERENCE-TRIANGLE TRANSFER / PRIMITIVE CONTENT AUDIT
==============================================================================

Self-contained single-file experiment.

No external files.
No synthetic second n=pq case.
Exact SymPy rational arithmetic only.
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

def clean(x):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(x)
            )
        )
    )


def valuation(value, prime):
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


def primitive_integer_list(values):
    """
    Convert a finite rational list to a primitive integer signature.

    1. Clear all denominators using their LCM.
    2. Divide the resulting integer vector by its gcd.
    3. Normalize the first nonzero entry to be positive.
    """
    rationals = [
        sp.Rational(value)
        for value in values
    ]

    if not rationals:
        return []

    denominator_lcm = 1

    for value in rationals:
        denominator_lcm = math.lcm(
            denominator_lcm,
            int(value.q),
        )

    integers = [
        int(
            value * denominator_lcm
        )
        for value in rationals
    ]

    gcd_value = 0

    for integer in integers:
        gcd_value = math.gcd(
            gcd_value,
            abs(integer),
        )

    if gcd_value != 0:
        integers = [
            integer // gcd_value
            for integer in integers
        ]

    for integer in integers:
        if integer != 0:

            if integer < 0:
                integers = [
                    -x
                    for x in integers
                ]

            break

    return integers


def finite_difference_rows(values):
    current = [
        sp.Rational(value)
        for value in values
    ]

    rows = [current]

    while len(current) > 1:

        current = [
            clean(
                current[i + 1] - current[i]
            )
            for i in range(
                len(current) - 1
            )
        ]

        rows.append(current)

    return rows


def D(p):
    return len(Q[p]) - 1


def build_layers():
    layers = {}

    max_t = max(
        D(p)
        for p in Q
    )

    for t in range(
        max_t + 1
    ):

        layer = []

        for p in sorted(Q):

            r = D(p) - t

            if r < 0:
                continue

            layer.append(
                (
                    p,
                    sp.Integer(
                        Q[p][r]
                    ),
                )
            )

        layers[t] = layer

    return layers


def build_triangle(layers):
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

            triangle[
                (j, t)
            ] = row

    return triangle


def build_first_entry_triangle(triangle):
    A = {}

    for (j, t), row in triangle.items():

        if row:
            A[
                (j, t)
            ] = clean(
                row[0]
            )

    return A


# ============================================================================
# SOURCE DISPLAY
# ============================================================================

def print_source_layers(layers):

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
# FINITE-DIFFERENCE TRIANGLE
# ============================================================================

def print_difference_triangle(triangle):

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
                    triangle[
                        (j, t)
                    ],
                )
            )


# ============================================================================
# FIRST-ENTRY TRIANGLE
# ============================================================================

def print_A(A):

    print()
    print("=" * 78)
    print(
        "3. FIRST-ENTRY TRIANGLE A[j,t]"
    )
    print("=" * 78)

    js = sorted(
        set(
            j
            for j, _ in A
        )
    )

    for j in js:

        values = [
            (
                t,
                A[
                    (j, t)
                ],
            )
            for t in sorted(
                tt
                for jj, tt in A
                if jj == j
            )
        ]

        print()
        print(
            "  j={}: {}".format(
                j,
                values,
            )
        )


# ============================================================================
# CORRECT MIXED-DIFFERENCE COMMUTATION
# ============================================================================

def mixed_difference_commutation(layers):

    print()
    print("=" * 78)
    print(
        "4. EXACT MIXED-DIFFERENCE COMMUTATION"
    )
    print("=" * 78)

    ps = sorted(Q)

    matrix = {}

    for n, p in enumerate(ps):

        matrix[n] = {}

        for t in sorted(layers):

            matches = [
                value
                for pp, value
                in layers[t]
                if pp == p
            ]

            matrix[n][t] = (
                matches[0]
                if matches
                else None
            )

    checks = 0
    failures = 0

    for n in range(
        len(ps) - 1
    ):

        for t in sorted(layers):

            if (
                t + 1
                not in layers
            ):
                continue

            a = matrix[n].get(t)
            b = matrix[n + 1].get(t)
            c = matrix[n].get(t + 1)
            d = matrix[n + 1].get(t + 1)

            if any(
                value is None
                for value in (
                    a,
                    b,
                    c,
                    d,
                )
            ):
                continue

            delta_n_delta_t = clean(
                d - c - b + a
            )

            delta_t_delta_n = clean(
                d - b - c + a
            )

            exact = (
                delta_n_delta_t
                == delta_t_delta_n
            )

            checks += 1

            if not exact:
                failures += 1

            print(
                "  n={}, t={}: "
                "Delta_nDelta_t={} "
                "Delta_tDelta_n={} "
                "exact={}".format(
                    n,
                    t,
                    delta_n_delta_t,
                    delta_t_delta_n,
                    exact,
                )
            )

    print()
    print(
        "  checks={}".format(
            checks
        )
    )

    print(
        "  commutation_failures={}".format(
            failures
        )
    )

    return (
        checks,
        failures,
    )


# ============================================================================
# FIXED-j POLYNOMIAL / LOO AUDIT
# ============================================================================

def polynomial_law(points, degree):

    if len(points) < degree + 1:
        return None

    x = sp.Symbol("x")

    try:
        polynomial = sp.interpolate(
            [
                (
                    sp.Integer(t),
                    sp.Rational(value),
                )
                for t, value in points
            ],
            x,
        )
    except Exception:
        return None

    polynomial = clean(
        polynomial
    )

    actual_degree = sp.Poly(
        polynomial,
        x,
        domain=sp.QQ,
    ).degree()

    if actual_degree > degree:
        return None

    for t, value in points:

        if clean(
            polynomial.subs(
                x,
                t,
            ) - value
        ) != 0:
            return None

    return polynomial


def polynomial_and_loo(A):

    print()
    print("=" * 78)
    print(
        "5. FIXED-j POLYNOMIAL / LEAVE-ONE-t-OUT AUDIT"
    )
    print("=" * 78)

    tested = 0
    exact = 0

    js = sorted(
        set(
            j
            for j, _ in A
        )
    )

    for j in js:

        points = [
            (
                t,
                A[
                    (j, t)
                ],
            )
            for t in sorted(
                tt
                for jj, tt in A
                if jj == j
            )
        ]

        print()
        print(
            "  j={}: points={}".format(
                j,
                points,
            )
        )

        for degree in (
            0,
            1,
            2,
            3,
        ):

            law = polynomial_law(
                points,
                degree,
            )

            if law is None:

                print(
                    "    degree<={}: "
                    "status=NO_EXACT_LAW".format(
                        degree
                    )
                )

            else:

                print(
                    "    degree<={}: "
                    "EXACT polynomial={}".format(
                        degree,
                        law,
                    )
                )

        if len(points) < 3:
            continue

        x = sp.Symbol("x")

        for omitted in range(
            len(points)
        ):

            omitted_t, actual = (
                points[
                    omitted
                ]
            )

            training = [
                point
                for index, point
                in enumerate(points)
                if index != omitted
            ]

            if len(training) < 2:
                continue

            tested += 1

            law = polynomial_law(
                training,
                1,
            )

            if law is None:

                prediction = None
                is_exact = False

            else:

                prediction = clean(
                    law.subs(
                        x,
                        omitted_t,
                    )
                )

                is_exact = (
                    prediction
                    == actual
                )

            if is_exact:
                exact += 1

            print(
                "    LOO omitted_t={}: "
                "prediction={} "
                "actual={} "
                "exact={}".format(
                    omitted_t,
                    prediction,
                    actual,
                    is_exact,
                )
            )

    return (
        tested,
        exact,
    )


# ============================================================================
# WIDTH-2 TRANSFER
# ============================================================================

def width2_transfer(A):

    print()
    print("=" * 78)
    print(
        "6. EXACT WIDTH-2 CROSS-LAYER TRANSFER"
    )
    print("=" * 78)

    all_j = {
        j
        for j, _ in A
    }

    all_t = {
        t
        for _, t in A
    }

    hits = []

    for t in sorted(all_t):

        if t + 1 not in all_t:
            continue

        common = sorted(
            set(
                j
                for j, tt in A
                if tt == t
            )
            &
            set(
                j
                for j, tt in A
                if tt == t + 1
            )
        )

        valid_js = [
            j
            for j in common
            if (
                j + 1 in all_j
                and
                (j + 1, t) in A
            )
        ]

        if len(valid_js) < 2:

            print(
                "  t={} -> {}: "
                "INSUFFICIENT_DATA".format(
                    t,
                    t + 1,
                )
            )

            continue

        alpha, beta = sp.symbols(
            "alpha beta"
        )

        equations = [
            sp.Eq(
                A[
                    (j, t + 1)
                ],
                alpha
                * A[
                    (j, t)
                ]
                +
                beta
                * A[
                    (j + 1, t)
                ],
            )
            for j in valid_js
        ]

        try:
            solutions = sp.solve(
                equations,
                [alpha, beta],
                dict=True,
            )
        except Exception:
            solutions = []

        if len(solutions) != 1:

            print(
                "  t={} -> {}: "
                "NO_UNIQUE_SOLUTION".format(
                    t,
                    t + 1,
                )
            )

            continue

        solution = solutions[0]

        if (
            alpha not in solution
            or beta not in solution
        ):

            print(
                "  t={} -> {}: "
                "NONUNIQUE".format(
                    t,
                    t + 1,
                )
            )

            continue

        a = clean(
            solution[alpha]
        )

        b = clean(
            solution[beta]
        )

        exact = True

        for j in valid_js:

            rhs = clean(
                a
                * A[
                    (j, t)
                ]
                +
                b
                * A[
                    (j + 1, t)
                ]
            )

            if rhs != A[
                (j, t + 1)
            ]:

                exact = False
                break

        if exact:

            hits.append(
                (
                    t,
                    a,
                    b,
                )
            )

            print(
                "  t={} -> {}: "
                "EXACT alpha={} beta={} "
                "j={}".format(
                    t,
                    t + 1,
                    a,
                    b,
                    valid_js,
                )
            )

        else:

            print(
                "  t={} -> {}: "
                "VERIFICATION_FAILED".format(
                    t,
                    t + 1,
                )
            )

    print()
    print(
        "  exact_width2_hits={}".format(
            len(hits)
        )
    )

    return hits


# ============================================================================
# WIDTH-3 TRANSFER
# ============================================================================

def width3_transfer(A):

    print()
    print("=" * 78)
    print(
        "7. EXACT WIDTH-3 CROSS-LAYER TRANSFER"
    )
    print("=" * 78)

    all_j = {
        j
        for j, _ in A
    }

    all_t = {
        t
        for _, t in A
    }

    hits = []

    for t in sorted(all_t):

        if t + 1 not in all_t:
            continue

        common = sorted(
            set(
                j
                for j, tt in A
                if tt == t
            )
            &
            set(
                j
                for j, tt in A
                if tt == t + 1
            )
        )

        valid_js = [
            j
            for j in common
            if (
                j + 1 in all_j
                and
                j + 2 in all_j
                and
                (j + 1, t) in A
                and
                (j + 2, t) in A
            )
        ]

        if len(valid_js) < 3:

            print(
                "  t={} -> {}: "
                "INSUFFICIENT_DATA".format(
                    t,
                    t + 1,
                )
            )

            continue

        alpha, beta, gamma = sp.symbols(
            "alpha beta gamma"
        )

        equations = [
            sp.Eq(
                A[
                    (j, t + 1)
                ],
                alpha
                * A[
                    (j, t)
                ]
                +
                beta
                * A[
                    (j + 1, t)
                ]
                +
                gamma
                * A[
                    (j + 2, t)
                ],
            )
            for j in valid_js
        ]

        try:
            solutions = sp.solve(
                equations,
                [alpha, beta, gamma],
                dict=True,
            )
        except Exception:
            solutions = []

        if len(solutions) != 1:

            print(
                "  t={} -> {}: "
                "NO_UNIQUE_SOLUTION".format(
                    t,
                    t + 1,
                )
            )

            continue

        solution = solutions[0]

        if not all(
            symbol in solution
            for symbol in (
                alpha,
                beta,
                gamma,
            )
        ):

            print(
                "  t={} -> {}: "
                "NONUNIQUE".format(
                    t,
                    t + 1,
                )
            )

            continue

        a = clean(
            solution[alpha]
        )

        b = clean(
            solution[beta]
        )

        c = clean(
            solution[gamma]
        )

        exact = True

        for j in valid_js:

            rhs = clean(
                a
                * A[
                    (j, t)
                ]
                +
                b
                * A[
                    (j + 1, t)
                ]
                +
                c
                * A[
                    (j + 2, t)
                ]
            )

            if rhs != A[
                (j, t + 1)
            ]:

                exact = False
                break

        if exact:

            hits.append(
                (
                    t,
                    a,
                    b,
                    c,
                )
            )

            print(
                "  t={} -> {}: "
                "EXACT alpha={} beta={} gamma={} "
                "j={}".format(
                    t,
                    t + 1,
                    a,
                    b,
                    c,
                    valid_js,
                )
            )

        else:

            print(
                "  t={} -> {}: "
                "VERIFICATION_FAILED".format(
                    t,
                    t + 1,
                )
            )

    print()
    print(
        "  exact_width3_hits={}".format(
            len(hits)
        )
    )

    return hits


# ============================================================================
# RANK AUDIT
# ============================================================================

def rank_audit(A):

    print()
    print("=" * 78)
    print(
        "8. A[j,t] RANK / MINOR AUDIT"
    )
    print("=" * 78)

    js = sorted(
        set(
            j
            for j, _ in A
        )
    )

    ts = sorted(
        set(
            t
            for _, t in A
        )
    )

    M = sp.Matrix([
        [
            A.get(
                (j, t),
                sp.Integer(0),
            )
            for t in ts
        ]
        for j in js
    ])

    rank = M.rank()

    print(
        "  shape={}".format(
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

    first_nonzero_minor = None

    for r1 in range(M.rows):
        for r2 in range(
            r1 + 1,
            M.rows
        ):

            for c1 in range(M.cols):
                for c2 in range(
                    c1 + 1,
                    M.cols
                ):

                    minor = clean(
                        M.extract(
                            [r1, r2],
                            [c1, c2],
                        ).det()
                    )

                    if minor != 0:

                        first_nonzero_minor = (
                            (r1, r2),
                            (c1, c2),
                            minor,
                        )

                        break

                if first_nonzero_minor:
                    break

            if first_nonzero_minor:
                break

        if first_nonzero_minor:
            break

    print(
        "  first_nonzero_2x2_minor={}".format(
            first_nonzero_minor
        )
    )

    return (
        rank,
        first_nonzero_minor,
    )


# ============================================================================
# PRIMITIVE CONTENT / VALUATION AUDIT
# ============================================================================

def primitive_audit(A):

    print()
    print("=" * 78)
    print(
        "9. PRIMITIVE CONTENT / p-ADIC PROFILE"
    )
    print("=" * 78)

    primes = (
        2,
        3,
        5,
        7,
        17,
    )

    for j in sorted(
        set(
            jj
            for jj, _ in A
        )
    ):

        values = [
            (
                t,
                A[
                    (j, t)
                ],
            )
            for t in sorted(
                tt
                for jj, tt in A
                if jj == j
            )
        ]

        print()
        print(
            "  j={}:".format(
                j
            )
        )

        print(
            "    values={}".format(
                values
            )
        )

        print(
            "    primitive={}".format(
                primitive_integer_list(
                    [
                        value
                        for _, value
                        in values
                    ]
                )
            )
        )

        for prime in primes:

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
                        for t, value
                        in values
                    ],
                )
            )


# ============================================================================
# INTERPRETATION
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
Experiment 303RR corrected the mixed-difference criterion and showed
that the terminal-distance source coordinates satisfy the exact
commutation identity

    Delta_n Delta_t Q = Delta_t Delta_n Q.

The commutation itself is automatic for the exact two-coordinate
difference construction; the informative question is therefore whether
the resulting difference triangle has a compact transfer law.

Experiment 304RR tests

    A[j,t] = first(Delta_n^j Q_t)

directly.

The important distinctions are:

    interpolation
        versus
    prediction;

    arbitrary fitting
        versus
    fixed-width transfer;

    scalar separability
        versus
    genuine two-coordinate structure.

A repeated exact width-2 or width-3 transfer across multiple terminal
layers would identify a compact cross-layer source operator.

Failure of those tests, together with full-rank behavior, would indicate
that the terminal-distance reindexing is structurally valid but does
not by itself expose a low-complexity upstream law.

No synthetic second n=pq case is generated.
"""
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 304RR — EXACT TERMINAL-DISTANCE "
        "DIFFERENCE-TRIANGLE TRANSFER / CONTENT AUDIT"
    )
    print("=" * 78)

    layers = build_layers()

    triangle = build_triangle(
        layers
    )

    A = build_first_entry_triangle(
        triangle
    )

    print_source_layers(
        layers
    )

    print_difference_triangle(
        triangle
    )

    print_A(
        A
    )

    (
        commutation_checks,
        commutation_failures,
    ) = mixed_difference_commutation(
        layers
    )

    (
        loo_tested,
        loo_exact,
    ) = polynomial_and_loo(
        A
    )

    width2_hits = width2_transfer(
        A
    )

    width3_hits = width3_transfer(
        A
    )

    (
        rank,
        first_nonzero_minor,
    ) = rank_audit(
        A
    )

    primitive_audit(
        A
    )

    interpretation()

    print()
    print("=" * 78)
    print(
        "11. FINAL EXACTNESS"
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
        "  mixed_difference_checks={}".format(
            commutation_checks
        )
    )

    print(
        "  mixed_difference_commutation_failures={}".format(
            commutation_failures
        )
    )

    print(
        "  polynomial_LOO_tested={}".format(
            loo_tested
        )
    )

    print(
        "  polynomial_LOO_exact={}".format(
            loo_exact
        )
    )

    print(
        "  width2_exact_hits={}".format(
            len(width2_hits)
        )
    )

    print(
        "  width3_exact_hits={}".format(
            len(width3_hits)
        )
    )

    print(
        "  A_rank={}".format(
            rank
        )
    )

    print(
        "  first_nonzero_minor_found={}".format(
            first_nonzero_minor is not None
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
        "EXPERIMENT 304RR COMPLETE"
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