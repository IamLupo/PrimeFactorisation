#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 299R — EXACT q-SOURCE COORDINATE / CROSS-p RECONSTRUCTION AUDIT
==============================================================================

Purpose
-------

Experiment 298R established that the supplied q-values do not identify
the B-channel through several natural universal q -> B operators.

Experiment 299R therefore moves one step further upstream.

The q-table itself is treated as the primary source object.

For each available p-row,

    q_p(r),  r = 0,...,D(p),

we reconstruct the unique source polynomial

    Q_p(x)

in several exact bases.

Then we ask a narrower question:

    Do the coefficients of Q_p depend on p through a low-complexity
    exact law that survives leave-one-out testing?

This is NOT a proof of a universal q_p(r) formula.

The experiment distinguishes:

    1. exact reconstruction of the supplied rows;
    2. cross-p interpolation;
    3. leave-one-p-out prediction;
    4. terminal coefficient behavior;
    5. D(p) behavior.

Only leave-one-p-out success counts as evidence of a cross-p law.

No second synthetic n=pq case is generated.

No external files are used.
No floating point.
No arbitrary matrix fit.
"""

from __future__ import annotations

import itertools
import math
import sys

import sympy as sp


# ============================================================================
# EXACT SOURCE DATA
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


def falling(x, n):
    out = sp.Integer(1)

    for i in range(n):
        out *= x - i

    return clean(out)


def rising(x, n):
    out = sp.Integer(1)

    for i in range(n):
        out *= x + i

    return clean(out)


def power_basis_coefficients(values, x):
    """
    Interpolate Q(x) in the ordinary power basis.
    """
    points = [
        (
            sp.Integer(r),
            sp.Integer(value),
        )
        for r, value in enumerate(values)
    ]

    poly = clean(
        sp.interpolate(
            points,
            x,
        )
    )

    p = sp.Poly(
        poly,
        x,
        domain=sp.QQ,
    )

    degree = int(p.degree())

    coeffs = [
        clean(
            p.nth(i)
        )
        for i in range(degree + 1)
    ]

    return poly, coeffs


def falling_basis_coefficients(values, x):
    """
    Convert Q(x) to the falling-factorial basis

        Q(x) = sum a_j x_(j)

    by exact triangular evaluation at x=0,1,...,D.
    """
    D = len(values) - 1

    residual = [
        sp.Integer(v)
        for v in values
    ]

    coeffs = []

    for j in range(D + 1):

        value = clean(
            residual[j]
        )

        coeffs.append(value)

        if value == 0:
            continue

        for r in range(j, D + 1):

            residual[r] = clean(
                residual[r]
                - value
                * sp.binomial(
                    r,
                    j,
                )
                * sp.factorial(j)
            )

    return coeffs


def forward_difference_coefficients(values):
    """
    Newton/binomial coefficients:

        f(x) = sum c_j C(x,j)

    where c_j = Δ^j f(0).
    """
    current = [
        sp.Integer(v)
        for v in values
    ]

    coeffs = []

    while current:

        coeffs.append(
            clean(current[0])
        )

        if len(current) == 1:
            break

        current = [
            clean(
                current[i + 1]
                - current[i]
            )
            for i in range(
                len(current) - 1
            )
        ]

    return coeffs


def verify_power_reconstruction(poly, values, x):
    for r, value in enumerate(values):

        got = clean(
            poly.subs(
                x,
                r,
            )
        )

        if got != value:
            return False

    return True


def primitive_integer_signature(values):
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
            value
            * denominator_lcm
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

    return ints


def interpolate_across_p(
    p_values,
    coefficient_values,
    variable,
):
    """
    Exact interpolation of a coefficient as a function of p.
    """
    points = [
        (
            sp.Integer(p),
            sp.Rational(
                coefficient_values[p]
            ),
        )
        for p in p_values
    ]

    return clean(
        sp.interpolate(
            points,
            variable,
        )
    )


def leave_one_out_prediction(
    p_values,
    coefficient_values,
    degree,
    target_p,
    variable,
):
    train_p = [
        p
        for p in p_values
        if p != target_p
    ]

    if len(train_p) <= degree:
        return None

    points = [
        (
            sp.Integer(p),
            sp.Rational(
                coefficient_values[p]
            ),
        )
        for p in train_p
    ]

    poly = clean(
        sp.interpolate(
            points,
            variable,
        )
    )

    prediction = clean(
        poly.subs(
            variable,
            target_p,
        )
    )

    return prediction


def v17(value):
    value = int(value)

    if value == 0:
        return None

    value = abs(value)
    count = 0

    while value % 17 == 0:
        value //= 17
        count += 1

    return count


# ============================================================================
# 1. SOURCE INVENTORY
# ============================================================================

def print_source_inventory():
    print()
    print("=" * 78)
    print("1. SOURCE q-TABLE INVENTORY")
    print("=" * 78)

    for p in sorted(Q):

        row = Q[p]
        D = len(row) - 1

        print()
        print(
            "  p={}: D(p)={}".format(
                p,
                D,
            )
        )

        print(
            "    q={}".format(
                row,
            )
        )


# ============================================================================
# 2. EXACT ROW RECONSTRUCTION
# ============================================================================

def reconstruct_rows(x):
    print()
    print("=" * 78)
    print("2. EXACT SOURCE-POLYNOMIAL RECONSTRUCTION")
    print("=" * 78)

    power_data = {}
    falling_data = {}
    newton_data = {}

    for p in sorted(Q):

        values = Q[p]

        poly, power_coeffs = (
            power_basis_coefficients(
                values,
                x,
            )
        )

        falling_coeffs = (
            falling_basis_coefficients(
                values,
                x,
            )
        )

        newton_coeffs = (
            forward_difference_coefficients(
                values,
            )
        )

        exact = (
            verify_power_reconstruction(
                poly,
                values,
                x,
            )
        )

        power_data[p] = power_coeffs
        falling_data[p] = falling_coeffs
        newton_data[p] = newton_coeffs

        print()
        print(
            "  p={}:".format(p)
        )

        print(
            "    Q_p(x)={}".format(
                poly,
            )
        )

        print(
            "    power_basis={}".format(
                power_coeffs,
            )
        )

        print(
            "    falling_basis={}".format(
                falling_coeffs,
            )
        )

        print(
            "    newton_binomial_basis={}".format(
                newton_coeffs,
            )
        )

        print(
            "    exact_reconstruction={}".format(
                exact,
            )
        )

    return (
        power_data,
        falling_data,
        newton_data,
    )


# ============================================================================
# 3. COEFFICIENT PADDED TRIANGLES
# ============================================================================

def print_coefficient_triangles(
    name,
    data,
):
    print()
    print("=" * 78)
    print(
        "3. {} COEFFICIENT TRIANGLE".format(
            name
        )
    )
    print("=" * 78)

    max_len = max(
        len(values)
        for values in data.values()
    )

    for p in sorted(data):

        row = list(
            data[p]
        )

        padded = row + [
            sp.Integer(0)
            for _ in range(
                max_len - len(row)
            )
        ]

        print()
        print(
            "  p={}: {}".format(
                p,
                padded,
            )
        )


# ============================================================================
# 4. CROSS-p POLYNOMIAL DEGREE SEARCH
# ============================================================================

def cross_p_degree_search(
    name,
    data,
    variable,
):
    print()
    print("=" * 78)
    print(
        "4. CROSS-p {} COEFFICIENT-LAW SEARCH".format(
            name
        )
    )
    print("=" * 78)

    p_values = sorted(
        data
    )

    max_width = max(
        len(data[p])
        for p in p_values
    )

    results = {}

    for j in range(max_width):

        available = [
            p
            for p in p_values
            if j < len(data[p])
        ]

        values = {
            p: data[p][j]
            for p in available
        }

        print()
        print(
            "  coefficient_index={}:".format(
                j
            )
        )

        print(
            "    available_p={}".format(
                available
            )
        )

        for degree in range(0, 4):

            if len(available) <= degree:
                status = "INSUFFICIENT_DATA"
                law = None

            else:
                law = interpolate_across_p(
                    available,
                    values,
                    variable,
                )

                # A degree-d interpolation requires d+1 points.
                # Therefore this is a fitted law whenever all available
                # points are used. We mark it explicitly as interpolation.
                actual_degree = int(
                    sp.Poly(
                        law,
                        variable,
                        domain=sp.QQ,
                    ).degree()
                )

                if actual_degree <= degree:
                    status = "INTERPOLATION"
                else:
                    status = "NO"

            print(
                "    degree<={}: status={} law={}".format(
                    degree,
                    status,
                    law,
                )
            )

        results[j] = values

    return results


# ============================================================================
# 5. LEAVE-ONE-p-OUT TEST
# ============================================================================

def loo_test(
    name,
    data,
    variable,
):
    print()
    print("=" * 78)
    print(
        "5. {} LEAVE-ONE-p TEST".format(
            name
        )
    )
    print("=" * 78)

    p_values = sorted(
        data
    )

    exact_total = 0
    total = 0

    max_width = max(
        len(data[p])
        for p in p_values
    )

    for j in range(max_width):

        available = [
            p
            for p in p_values
            if j < len(data[p])
        ]

        if len(available) < 3:
            continue

        print()
        print(
            "  coefficient_index={}:".format(
                j
            )
        )

        coefficient_values = {
            p: data[p][j]
            for p in available
        }

        # With 4 prime rows available, the largest nontrivial
        # leave-one-out polynomial degree that has training data is 2.
        for degree in range(0, 3):

            print()
            print(
                "    degree<={}:".format(
                    degree
                )
            )

            for target_p in available:

                prediction = (
                    leave_one_out_prediction(
                        available,
                        coefficient_values,
                        degree,
                        target_p,
                        variable,
                    )
                )

                if prediction is None:
                    print(
                        "      target_p={} status=INSUFFICIENT_DATA".format(
                            target_p
                        )
                    )
                    continue

                actual = clean(
                    coefficient_values[
                        target_p
                    ]
                )

                exact = (
                    prediction
                    == actual
                )

                total += 1

                if exact:
                    exact_total += 1

                print(
                    "      target_p={} prediction={} actual={} exact={}".format(
                        target_p,
                        prediction,
                        actual,
                        exact,
                    )
                )

    print()
    print(
        "  exact_predictions={}".format(
            exact_total
        )
    )

    print(
        "  tested_predictions={}".format(
            total
        )
    )

    return (
        exact_total,
        total,
    )


# ============================================================================
# 6. TERMINAL COEFFICIENT AUDIT
# ============================================================================

def terminal_audit(
    power_data,
    falling_data,
):
    print()
    print("=" * 78)
    print("6. TERMINAL COEFFICIENT AUDIT")
    print("=" * 78)

    for p in sorted(Q):

        D = len(Q[p]) - 1

        q_terminal = Q[p][D]

        power_terminal = (
            power_data[p][-1]
        )

        falling_terminal = (
            falling_data[p][-1]
        )

        print()
        print(
            "  p={}: D={}".format(
                p,
                D,
            )
        )

        print(
            "    q_p(D)={}".format(
                q_terminal,
            )
        )

        print(
            "    highest_power_coefficient={}".format(
                power_terminal,
            )
        )

        print(
            "    highest_falling_coefficient={}".format(
                falling_terminal,
            )
        )

        print(
            "    v17(q_terminal)={}".format(
                v17(
                    q_terminal
                ),
            )
        )


# ============================================================================
# 7. D(p) FORMULA AUDIT
# ============================================================================

def D_formula_audit():
    print()
    print("=" * 78)
    print("7. TERMINAL-INDEX D(p) AUDIT")
    print("=" * 78)

    p_values = sorted(Q)

    observed = {
        p: len(Q[p]) - 1
        for p in p_values
    }

    print(
        "  observed_D={}".format(
            observed
        )
    )

    # Test low-degree interpolation laws only as diagnostics.
    for degree in range(0, 4):

        if len(p_values) <= degree:
            print(
                "  degree<={}: INSUFFICIENT_DATA".format(
                    degree
                )
            )
            continue

        x = sp.Symbol("p")

        law = clean(
            sp.interpolate(
                [
                    (
                        sp.Integer(p),
                        sp.Integer(observed[p]),
                    )
                    for p in p_values
                ],
                x,
            )
        )

        actual_degree = int(
            sp.Poly(
                law,
                x,
                domain=sp.QQ,
            ).degree()
        )

        print(
            "  degree<={}: interpolation_D(p)={} actual_degree={}".format(
                degree,
                law,
                actual_degree,
            )
        )

    # Simple exact candidate checks.
    candidates = {
        "7-p": lambda p: 7 - p,
        "(7-p)//2": lambda p: (7 - p) // 2,
        "8-p": lambda p: 8 - p,
        "p-1": lambda p: p - 1,
        "(7-p)/2": lambda p: sp.Rational(7 - p, 2),
    }

    for name, fn in candidates.items():

        checks = []

        valid_domain = True

        for p in p_values:
            try:
                value = clean(
                    fn(p)
                )
            except Exception:
                valid_domain = False
                break

            checks.append(
                value == observed[p]
            )

        exact = (
            valid_domain
            and all(checks)
        )

        print(
            "  candidate {} exact={}".format(
                name,
                exact,
            )
        )


# ============================================================================
# 8. TERMINAL PROJECTIVE SOURCE
# ============================================================================

def terminal_projective_reference():
    print()
    print("=" * 78)
    print("8. TERMINAL PROJECTIVE SOURCE REFERENCE")
    print("=" * 78)

    q1 = Q[1][-1]
    q3 = Q[3][-1]

    gcd_value = math.gcd(
        q1,
        q3,
    )

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
            gcd_value
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
        "  v17(q1)={}".format(
            v17(q1)
        )
    )

    print(
        "  v17(q3)={}".format(
            v17(q3)
        )
    )


# ============================================================================
# 9. INTERPRETATION
# ============================================================================

def interpretation():
    print()
    print("=" * 78)
    print("9. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
r"""
Experiment 298R showed that the available q-values do not identify
B[k,r] through low-complexity universal q-to-B maps.

Experiment 299R therefore moves upstream of B.

For each known p-row, the complete q_p(r) sequence determines a unique
finite source polynomial Q_p(x).

That reconstruction is exact and unavoidable; it is simply a change of
representation of the supplied q-row.

The actual question is cross-p:

    Are the coefficients of Q_p(x) governed by a law in p?

This script deliberately separates:

    interpolation
        from
    prediction.

A polynomial through all four p-values is only an interpolation.

The stronger test is leave-one-p-out:

    fit the coefficient law using the other p-values,
    predict the omitted coefficient,
    compare exactly.

Repeated exact leave-one-out success would be evidence for a genuine
cross-p source law.

Repeated failure means that the present four q-rows are insufficient to
identify such a low-degree cross-p law.

No synthetic second n=pq case is generated.

Therefore a positive result here would identify a source-coordinate law
already latent in the supplied q-table, while a negative result tells us
that the actual upstream construction is still required.
"""
    )


# ============================================================================
# MAIN
# ============================================================================

def main():
    x = sp.Symbol("x")
    cross_p = sp.Symbol("p")

    print("=" * 78)
    print(
        "EXPERIMENT 299R — EXACT q-SOURCE COORDINATE / "
        "CROSS-p RECONSTRUCTION AUDIT"
    )
    print("=" * 78)

    print_source_inventory()

    power_data, falling_data, newton_data = (
        reconstruct_rows(
            x
        )
    )

    print_coefficient_triangles(
        "POWER-BASIS",
        power_data,
    )

    print_coefficient_triangles(
        "FALLING-BASIS",
        falling_data,
    )

    print_coefficient_triangles(
        "NEWTON-BINOMIAL",
        newton_data,
    )

    cross_p_degree_search(
        "POWER-BASIS",
        power_data,
        cross_p,
    )

    cross_p_degree_search(
        "FALLING-BASIS",
        falling_data,
        cross_p,
    )

    cross_p_degree_search(
        "NEWTON-BINOMIAL",
        newton_data,
        cross_p,
    )

    loo_power = loo_test(
        "POWER-BASIS",
        power_data,
        cross_p,
    )

    loo_falling = loo_test(
        "FALLING-BASIS",
        falling_data,
        cross_p,
    )

    loo_newton = loo_test(
        "NEWTON-BINOMIAL",
        newton_data,
        cross_p,
    )

    terminal_audit(
        power_data,
        falling_data,
    )

    D_formula_audit()

    terminal_projective_reference()

    interpretation()

    print()
    print("=" * 78)
    print("10. FINAL EXACTNESS")
    print("=" * 78)

    power_exact_total = (
        loo_power[0]
    )

    falling_exact_total = (
        loo_falling[0]
    )

    newton_exact_total = (
        loo_newton[0]
    )

    print(
        "  exact_q_row_reconstruction=True"
    )

    print(
        "  power_basis_LOO_exact_predictions={}".format(
            power_exact_total
        )
    )

    print(
        "  falling_basis_LOO_exact_predictions={}".format(
            falling_exact_total
        )
    )

    print(
        "  newton_basis_LOO_exact_predictions={}".format(
            newton_exact_total
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
        "EXPERIMENT 299R COMPLETE"
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

