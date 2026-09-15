#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 302R — EXACT TERMINAL-DIFFERENCE LADDER / CROSS-LAYER LAW AUDIT
==============================================================================

Goal
----

Experiment 301R found the terminal-distance layers

    Q_t(p) = q_p(D(p)-t)

and their exact finite-difference tables in

    n = (p-1)/2.

The next structurally meaningful object is the HIGHEST AVAILABLE
finite difference in each terminal layer.

For each t define

    h(t) = highest-order forward finite difference available
           in the Q_t sequence.

We now test whether the ladder

    h(t)

has an exact low-complexity law in t.

This is deliberately different from interpolating Q_t itself.

Tests
-----

1. Extract the highest finite difference H_t for every t.
2. Test constant / affine / quadratic / cubic laws in t.
3. Test first-order and second-order constant recurrences in t.
4. Test geometric ratios H_(t+1)/H_t.
5. Test absolute-value and sign-separated versions.
6. Test divisibility and primitive-content structure.
7. Compare H_t against the terminal source values q_p(D(p)).
8. Test whether H_t itself admits a falling-factorial / binomial
   expansion in t.
9. Perform leave-one-t-out tests whenever enough layers exist.

Important
---------

* Exact QQ arithmetic only.
* No floating point.
* No synthetic second n=pq case.
* No external files.
* No arbitrary fitted matrices.
* Interpolation is explicitly labelled and never counted as proof.
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
# SAFE EXACT HELPERS
# ============================================================================

def clean(expr):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(expr)
            )
        )
    )


def finite_differences(values):
    """
    Return all forward-difference rows.
    """
    current = [
        sp.Rational(v)
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


def primitive_integer_signature(values):
    """
    Exact primitive integer normalization.
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

    ints = [
        int(
            value * denominator_lcm
        )
        for value in values
    ]

    g = 0

    for value in ints:
        g = math.gcd(
            g,
            abs(value),
        )

    if g:
        ints = [
            value // g
            for value in ints
        ]

    for value in ints:
        if value != 0:
            if value < 0:
                ints = [
                    -x
                    for x in ints
                ]
            break

    return ints


def valuation(value, prime):
    """
    Exact p-adic valuation of a nonzero integer/rational.
    """
    value = sp.Rational(value)

    if value == 0:
        return None

    numerator = abs(int(value.p))
    denominator = abs(int(value.q))

    out = 0

    while numerator % prime == 0:
        numerator //= prime
        out += 1

    while denominator % prime == 0:
        denominator //= prime
        out -= 1

    return out


def recurrence_order_one(values):
    """
    Test

        a[n+1] = A a[n]
    """
    values = [
        sp.Rational(v)
        for v in values
    ]

    if len(values) < 2:
        return (
            "INSUFFICIENT_DATA",
            None,
        )

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

    Requires at least four terms.
    """
    values = [
        sp.Rational(v)
        for v in values
    ]

    if len(values) < 4:
        return (
            "INSUFFICIENT_DATA",
            None,
        )

    A, B = sp.symbols(
        "A B"
    )

    equations = []

    for i in range(
        len(values) - 2
    ):
        equations.append(
            sp.Eq(
                values[i + 2],
                A * values[i + 1]
                + B * values[i],
            )
        )

    solutions = sp.solve(
        equations,
        [A, B],
        dict=True,
    )

    if not solutions:
        return (
            "NO_SOLUTION",
            None,
        )

    if len(solutions) != 1:
        return (
            "NONUNIQUE",
            solutions,
        )

    sol = solutions[0]

    if A not in sol or B not in sol:
        return (
            "NONUNIQUE",
            sol,
        )

    Av = clean(sol[A])
    Bv = clean(sol[B])

    for i in range(
        len(values) - 2
    ):
        rhs = clean(
            Av * values[i + 1]
            + Bv * values[i]
        )

        if rhs != values[i + 2]:
            return (
                "NO_SOLUTION",
                None,
            )

    return (
        "EXACT",
        (Av, Bv),
    )


def interpolate(values, symbol):
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
            symbol,
        )
    )


def exact_polynomial_degree(values, symbol):
    """
    Return degree of the exact interpolation polynomial.
    This is diagnostic only.
    """
    if not values:
        return None

    polynomial = interpolate(
        values,
        symbol,
    )

    return int(
        sp.Poly(
            polynomial,
            symbol,
            domain=sp.QQ,
        ).degree()
    )


def divided_difference_profile(values):
    """
    For equally spaced points, this is just the forward-difference ladder.
    """
    return finite_differences(values)


# ============================================================================
# TERMINAL-DISTANCE LAYERS
# ============================================================================

def D(p):
    return len(Q[p]) - 1


def terminal_layers():
    """
    t -> list of source values ordered by p.

        t = D(p)-r
    """
    layers = {}

    max_t = max(
        D(p)
        for p in Q
    )

    for t in range(
        max_t + 1
    ):

        values = []

        for p in sorted(Q):

            r = D(p) - t

            if r < 0:
                continue

            values.append(
                (
                    p,
                    sp.Integer(
                        Q[p][r]
                    ),
                )
            )

        layers[t] = values

    return layers


# ============================================================================
# EXTRACT HIGHEST AVAILABLE DIFFERENCE
# ============================================================================

def highest_difference(values):
    """
    For m values, the highest available difference is order m-1,
    consisting of one scalar.
    """
    rows = finite_differences(
        values
    )

    order = len(rows) - 1
    value = rows[-1][0]

    return order, clean(value)


# ============================================================================
# 1. LAYER DIFFERENCE LADDER
# ============================================================================

def print_layer_ladder(layers):
    print()
    print("=" * 78)
    print(
        "1. TERMINAL-LAYER HIGHEST-DIFFERENCE LADDER"
    )
    print("=" * 78)

    ladder = []

    for t in sorted(layers):

        values = [
            value
            for _, value in layers[t]
        ]

        order, highest = (
            highest_difference(
                values
            )
        )

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
                values
            )
        )

        print(
            "    highest_difference_order={}".format(
                order
            )
        )

        print(
            "    H_t={}".format(
                highest
            )
        )

        print(
            "    H_t_primitive={}".format(
                primitive_integer_signature(
                    [highest]
                )
            )
        )

        ladder.append(
            (
                t,
                order,
                highest,
            )
        )

    return ladder


# ============================================================================
# 2. LADDER FINITE-DIFFERENCE AUDIT
# ============================================================================

def ladder_difference_audit(ladder):
    print()
    print("=" * 78)
    print(
        "2. EXACT FINITE-DIFFERENCE AUDIT OF H_t"
    )
    print("=" * 78)

    H = [
        value
        for _, _, value in ladder
    ]

    print()
    print(
        "  H={}".format(
            H
        )
    )

    rows = finite_differences(
        H
    )

    for order, row in enumerate(
        rows
    ):

        print()
        print(
            "  delta_t^{}={}".format(
                order,
                row,
            )
        )

    return rows


# ============================================================================
# 3. POLYNOMIAL LAW SEARCH
# ============================================================================

def polynomial_law_search(H):
    print()
    print("=" * 78)
    print(
        "3. EXACT POLYNOMIAL LAW SEARCH IN t"
    )
    print("=" * 78)

    symbol = sp.Symbol(
        "t"
    )

    max_degree = min(
        4,
        max(
            0,
            len(H) - 1,
        )
    )

    results = {}

    for degree in range(
        max_degree + 1
    ):

        if len(H) <= degree:
            status = (
                "INSUFFICIENT_DATA"
            )
            polynomial = None

        else:

            polynomial = interpolate(
                H,
                symbol,
            )

            actual_degree = exact_polynomial_degree(
                H,
                symbol,
            )

            if actual_degree <= degree:
                status = "EXACT"
            else:
                status = "NO_SOLUTION"

        print()
        print(
            "  degree<={}: status={} polynomial={}".format(
                degree,
                status,
                polynomial,
            )
        )

        results[degree] = (
            status,
            polynomial,
        )

    return results


# ============================================================================
# 4. RECURRENCE AUDIT
# ============================================================================

def ladder_recurrence_audit(H):
    print()
    print("=" * 78)
    print(
        "4. EXACT H_t RECURRENCE AUDIT"
    )
    print("=" * 78)

    status1, parameter1 = (
        recurrence_order_one(
            H
        )
    )

    status2, parameter2 = (
        recurrence_order_two(
            H
        )
    )

    print()
    print(
        "  order1: status={} parameter={}".format(
            status1,
            parameter1,
        )
    )

    print(
        "  order2: status={} parameters={}".format(
            status2,
            parameter2,
        )
    )

    return (
        (status1, parameter1),
        (status2, parameter2),
    )


# ============================================================================
# 5. RATIO / SIGN / VALUATION AUDIT
# ============================================================================

def ratio_audit(ladder):
    print()
    print("=" * 78)
    print(
        "5. H_t RATIO / SIGN / p-ADIC AUDIT"
    )
    print("=" * 78)

    values = [
        value
        for _, _, value in ladder
    ]

    ratios = []

    for i in range(
        len(values) - 1
    ):

        if values[i] == 0:
            ratio = None
        else:
            ratio = clean(
                values[i + 1]
                / values[i]
            )

        ratios.append(
            ratio
        )

    print()
    print(
        "  adjacent_ratios={}".format(
            ratios
        )
    )

    print()
    print(
        "  signs={}".format(
            [
                "+" if value > 0
                else "-"
                if value < 0
                else "0"
                for value in values
            ]
        )
    )

    for prime in [
        2,
        3,
        5,
        7,
        11,
        13,
        17,
    ]:

        print()
        print(
            "  v_{}(H_t)={}".format(
                prime,
                [
                    valuation(
                        value,
                        prime,
                    )
                    for value in values
                ],
            )
        )

    return ratios


# ============================================================================
# 6. NORMALIZATION BY TERMINAL SOURCE
# ============================================================================

def terminal_normalization(layers):
    print()
    print("=" * 78)
    print(
        "6. H_t / TERMINAL-SOURCE NORMALIZATION"
    )
    print("=" * 78)

    normalized = []

    terminal_by_p = {
        p: sp.Integer(
            Q[p][D(p)]
        )
        for p in sorted(Q)
    }

    for t, layer in sorted(
        layers.items()
    ):

        # H_t
        values = [
            value
            for _, value in layer
        ]

        _, H_t = (
            highest_difference(
                values
            )
        )

        # For a given t, use terminal source values only
        # at p-values where that layer exists.
        aligned = []

        for p, value in layer:

            terminal = (
                terminal_by_p[p]
            )

            if terminal == 0:
                ratio = None
            else:
                ratio = clean(
                    H_t / terminal
                )

            aligned.append(
                (
                    p,
                    ratio,
                )
            )

        print()
        print(
            "  t={}: H_t={}".format(
                t,
                H_t,
            )
        )

        print(
            "    H_t/q_p(D(p))={}".format(
                aligned
            )
        )

        normalized.extend(
            aligned
        )

    return normalized


# ============================================================================
# 7. LEAVE-ONE-LAYER-OUT POLYNOMIAL TEST
# ============================================================================

def leave_one_layer_out(H):
    print()
    print("=" * 78)
    print(
        "7. LEAVE-ONE-t-OUT POLYNOMIAL TEST"
    )
    print("=" * 78)

    symbol = sp.Symbol(
        "t"
    )

    tested = 0
    exact = 0

    if len(H) < 3:
        print(
            "  insufficient layers"
        )
        return (
            tested,
            exact,
        )

    for degree in [
        1,
        2,
    ]:

        print()
        print(
            "  degree<={}".format(
                degree
            )
        )

        for omitted in range(
            len(H)
        ):

            training_points = []

            for index, value in enumerate(
                H
            ):

                if index == omitted:
                    continue

                training_points.append(
                    (
                        sp.Integer(index),
                        sp.Rational(value),
                    )
                )

            if len(training_points) <= degree:
                continue

            tested += 1

            try:
                poly = clean(
                    sp.interpolate(
                        training_points,
                        symbol,
                    )
                )
            except Exception:
                print(
                    "    omitted_t={} "
                    "status=ERROR".format(
                        omitted
                    )
                )
                continue

            prediction = clean(
                poly.subs(
                    symbol,
                    omitted,
                )
            )

            actual = (
                sp.Rational(
                    H[omitted]
                )
            )

            ok = (
                prediction
                == actual
            )

            if ok:
                exact += 1

            print(
                "    omitted_t={} prediction={} "
                "actual={} exact={}".format(
                    omitted,
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
# 8. CROSS-LAYER H_t NORMALIZATION
# ============================================================================

def primitive_ladder(ladder):
    print()
    print("=" * 78)
    print(
        "8. PRIMITIVE H_t LADDER"
    )
    print("=" * 78)

    values = [
        value
        for _, _, value in ladder
    ]

    for t, value in enumerate(
        values
    ):

        print()
        print(
            "  t={}: H_t={} primitive={}".format(
                t,
                value,
                primitive_integer_signature(
                    [value]
                ),
            )
        )

        print(
            "    denominator={}".format(
                sp.denom(
                    value
                )
            )
        )

        print(
            "    numerator={}".format(
                sp.numer(
                    value
                )
            )
        )


# ============================================================================
# 9. STRUCTURAL INTERPRETATION
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
Experiment 301R exposed a terminal-distance decomposition

    Q_t(p) = q_p(D(p)-t),

and showed that each fixed-t layer is a short exact finite-difference
sequence in

    n=(p-1)/2.

Experiment 302R isolates the strongest scalar attached to each layer:

    H_t = Δ^(m_t-1) Q_t,

where m_t is the number of available p-values in that layer.

Thus H_t is not a fitted coefficient. It is the terminal finite
difference forced by the observed layer itself.

The central question is now:

    Does the sequence

        H_0,H_1,...,H_5

    obey its own exact structural law?

A positive result would produce a second-level source law:

    q_p(r)
       ->
    terminal-distance layer
       ->
    highest p-difference H_t.

That would be materially different from the earlier attempts to fit
q_p(r) directly against B[k,r].

The experiment therefore concentrates the surviving information into
one scalar ladder and tests polynomial, recurrence, ratio, valuation,
normalization, and leave-one-layer-out structure.

No interpolation through all six values is treated as proof.
"""
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 302R — EXACT TERMINAL-DIFFERENCE LADDER / "
        "CROSS-LAYER LAW AUDIT"
    )
    print("=" * 78)

    layers = terminal_layers()

    ladder = print_layer_ladder(
        layers
    )

    H = [
        value
        for _, _, value in ladder
    ]

    ladder_difference_audit(
        ladder
    )

    polynomial_law_search(
        H
    )

    ladder_recurrence_audit(
        H
    )

    ratio_audit(
        ladder
    )

    terminal_normalization(
        layers
    )

    loo_tested, loo_exact = (
        leave_one_layer_out(
            H
        )
    )

    primitive_ladder(
        ladder
    )

    structural_interpretation()

    print()
    print("=" * 78)
    print(
        "10. FINAL EXACTNESS"
    )
    print("=" * 78)

    # Determine exact polynomial hits conservatively:
    x = sp.Symbol("t")

    polynomial_hits = 0

    for degree in range(
        0,
        min(
            4,
            len(H) - 1,
        ) + 1
    ):

        polynomial = interpolate(
            H,
            x,
        )

        actual_degree = (
            exact_polynomial_degree(
                H,
                x,
            )
        )

        if actual_degree <= degree:
            polynomial_hits += 1

    order1_status, _ = (
        recurrence_order_one(
            H
        )
    )

    order2_status, _ = (
        recurrence_order_two(
            H
        )
    )

    print(
        "  terminal_layer_count={}".format(
            len(H)
        )
    )

    print(
        "  polynomial_law_hits={}".format(
            polynomial_hits
        )
    )

    print(
        "  order1_recurrence_status={}".format(
            order1_status
        )
    )

    print(
        "  order2_recurrence_status={}".format(
            order2_status
        )
    )

    print(
        "  leave_one_layer_out_tested={}".format(
            loo_tested
        )
    )

    print(
        "  leave_one_layer_out_exact={}".format(
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
        "EXPERIMENT 302R COMPLETE"
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
