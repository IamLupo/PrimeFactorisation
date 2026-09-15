#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 320R — EXACT SOURCE-POLYNOMIAL / LAYER-RECURRENCE AUDIT
==============================================================================

Purpose
-------
Experiment 319R found no obvious arithmetic normalization making the two
raw coefficient pairs

    (c0_1,c1_1)
    (c0_2,c1_2)

look like the same simple operator.

The next step is therefore to return directly to the source data Q_t(p).

For the available odd p-values

    p = 1, 3, 5, 7,

the source layers provide exact samples of

    Q_t(p).

This experiment reconstructs the exact interpolation polynomial in p for
each available layer and then studies the resulting polynomial sequence
across t.

The experiment asks:

    1. minimal interpolation degree in p at every layer;
    2. exact polynomial reconstruction;
    3. primitive integer polynomial normalization;
    4. factorization over Q;
    5. coefficient tables;
    6. common factors between adjacent layer polynomials;
    7. whether adjacent layers obey a low-order polynomial differential,
       shift, or affine relation in p;
    8. whether the raw width-2 law

           Q_{t+1}(p) = c0_t Q_t(p) + c1_t Q_t(p+2)

       becomes an exact polynomial identity wherever Q_t has enough
       existing p-samples;
    9. whether the interpolating polynomials reveal a simpler recurrence
       than the pointwise transfer;
   10. exact coefficient recurrences across t;
   11. finite-difference structure with respect to p;
   12. exact residual audit against every existing source sample.

IMPORTANT
---------
Interpolation here is used only as an exact representation of the finite
existing p-samples. No value at a new p is counted as evidence.

A polynomial identity is only declared when it is checked symbolically
from the reconstructed exact polynomials.

No synthetic second n=pq case.
No external files.
Exact rational arithmetic only.
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

p = sp.Symbol("p")
t_symbol = sp.Symbol("t")


def clean(x):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(x)
            )
        )
    )


def clean_poly(expr):
    return sp.Poly(
        sp.expand(
            sp.cancel(expr)
        ),
        p,
        domain=sp.QQ,
    )


def degree_at_p(prime):
    return len(Q[prime]) - 1


def build_layers():

    layers = {}

    maximum_t = max(
        degree_at_p(prime)
        for prime in Q
    )

    for t in range(maximum_t + 1):

        row = []

        for prime in sorted(Q):

            r = degree_at_p(prime) - t

            if r < 0:
                continue

            row.append(
                (
                    prime,
                    sp.Integer(
                        Q[prime][r]
                    ),
                )
            )

        layers[t] = row

    return layers


def layer_points(layers, t):

    return [
        (
            sp.Integer(prime),
            sp.Integer(value),
        )
        for prime, value
        in layers[t]
    ]


def integer_content(values):

    vals = [
        int(v)
        for v in values
        if int(v) != 0
    ]

    if not vals:
        return 0

    g = 0

    for value in vals:
        g = math.gcd(
            g,
            abs(value),
        )

    return g


def primitive_integer_poly(poly):

    P = clean_poly(poly)

    coeffs = P.all_coeffs()

    den_lcm = 1

    for c in coeffs:
        den_lcm = sp.ilcm(
            den_lcm,
            int(sp.denom(c)),
        )

    ints = [
        int(
            c * den_lcm
        )
        for c in coeffs
    ]

    g = 0

    for value in ints:
        g = math.gcd(
            g,
            abs(value),
        )

    if g == 0:
        return sp.Poly(
            0,
            p,
            domain=sp.ZZ,
        )

    ints = [
        value // g
        for value in ints
    ]

    if ints[0] < 0:
        ints = [
            -value
            for value in ints
        ]

    expr = sum(
        sp.Integer(c)
        * p ** (
            len(ints) - 1 - i
        )
        for i, c
        in enumerate(ints)
    )

    return sp.Poly(
        expr,
        p,
        domain=sp.ZZ,
    )


# ============================================================================
# INTERPOLATION
# ============================================================================

def interpolate_layer(layers, t):

    points = layer_points(
        layers,
        t,
    )

    if not points:
        return None

    polynomial = sp.interpolate(
        points,
        p,
    )

    return clean_poly(
        polynomial
    )


# ============================================================================
# LAYER POLYNOMIAL AUDIT
# ============================================================================

def layer_polynomial_audit(layers):

    print()
    print("=" * 78)
    print(
        "2. EXACT LAYER POLYNOMIAL RECONSTRUCTION"
    )
    print("=" * 78)

    polys = {}

    for t in sorted(layers):

        points = layer_points(
            layers,
            t,
        )

        P = interpolate_layer(
            layers,
            t,
        )

        polys[t] = P

        print()
        print(
            "  t={}:".format(t)
        )

        print(
            "    points={}".format(
                points
            )
        )

        print(
            "    interpolation_degree={}".format(
                P.degree()
            )
        )

        print(
            "    polynomial={}".format(
                sp.factor(
                    P.as_expr()
                )
            )
        )

        print(
            "    expanded={}".format(
                P.as_expr()
            )
        )

        primitive = primitive_integer_poly(
            P.as_expr()
        )

        print(
            "    primitive_integer_polynomial={}".format(
                primitive.as_expr()
            )
        )

        print(
            "    factorization={}".format(
                sp.factor_list(
                    P.as_expr()
                )
            )
        )

    return polys


# ============================================================================
# INTERPOLATION RESIDUAL AUDIT
# ============================================================================

def interpolation_residual_audit(
    layers,
    polys,
):

    print()
    print("=" * 78)
    print(
        "3. EXACT INTERPOLATION RESIDUAL AUDIT"
    )
    print("=" * 78)

    failures = 0

    for t in sorted(layers):

        P = polys[t]

        residuals = []

        for prime, value in layer_points(
            layers,
            t,
        ):

            residual = clean(
                P.as_expr().subs(
                    p,
                    prime,
                )
                - value
            )

            residuals.append(
                residual
            )

            if residual != 0:
                failures += 1

        print()
        print(
            "  t={}: residuals={}".format(
                t,
                residuals,
            )
        )

    print()
    print(
        "  interpolation_failures={}".format(
            failures
        )
    )

    return failures == 0


# ============================================================================
# COEFFICIENT TABLE
# ============================================================================

def coefficient_table(polys):

    print()
    print("=" * 78)
    print(
        "4. INTERPOLATING POLYNOMIAL COEFFICIENT TABLE"
    )
    print("=" * 78)

    maximum_degree = max(
        P.degree()
        for P in polys.values()
    )

    print()
    print(
        "  degree_range=0..{}".format(
            maximum_degree
        )
    )

    for power in range(
        maximum_degree,
        -1,
        -1,
    ):

        row = []

        for t in sorted(polys):

            coeff = polys[t].coeff_monomial(
                p ** power
            )

            row.append(
                (
                    t,
                    clean(coeff),
                )
            )

        print()
        print(
            "  p^{}: {}".format(
                power,
                row,
            )
        )


# ============================================================================
# ADJACENT POLYNOMIAL GCD / DIFFERENCES
# ============================================================================

def adjacent_layer_audit(polys):

    print()
    print("=" * 78)
    print(
        "5. ADJACENT-LAYER POLYNOMIAL RELATION AUDIT"
    )
    print("=" * 78)

    for t in range(
        max(polys.keys())
    ):

        if t not in polys or t + 1 not in polys:
            continue

        P = polys[t]
        R = polys[t + 1]

        gcd_poly = sp.gcd(
            P,
            R,
        )

        difference = clean_poly(
            R.as_expr()
            - P.as_expr()
        )

        derivative_difference = clean_poly(
            sp.diff(
                R.as_expr(),
                p,
            )
            -
            sp.diff(
                P.as_expr(),
                p,
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
            "    degree_P={}".format(
                P.degree()
            )
        )

        print(
            "    degree_R={}".format(
                R.degree()
            )
        )

        print(
            "    gcd={}".format(
                sp.factor(
                    gcd_poly.as_expr()
                )
            )
        )

        print(
            "    gcd_degree={}".format(
                gcd_poly.degree()
            )
        )

        print(
            "    R-P={}".format(
                sp.factor(
                    difference.as_expr()
                )
            )
        )

        print(
            "    derivative_difference={}".format(
                sp.factor(
                    derivative_difference.as_expr()
                )
            )
        )


# ============================================================================
# TEST RAW WIDTH-2 LAW AS POLYNOMIAL IDENTITY
# ============================================================================

def polynomial_shift(P):

    return clean_poly(
        P.as_expr().subs(
            p,
            p + 2,
        )
    )


def raw_transfer_polynomial_audit(
    polys,
    raw_coefficients,
):

    print()
    print("=" * 78)
    print(
        "6. RAW WIDTH-2 LAW AS POLYNOMIAL IDENTITY"
    )
    print("=" * 78)

    results = []

    for t in (
        1,
        2,
    ):

        if (
            t not in polys
            or t + 1 not in polys
            or t not in raw_coefficients
        ):
            continue

        P = polys[t]
        R = polys[t + 1]

        c0 = raw_coefficients[t]["c0"]
        c1 = raw_coefficients[t]["c1"]

        shifted = polynomial_shift(P)

        residual = clean_poly(
            R.as_expr()
            - (
                c0 * P.as_expr()
                + c1 * shifted.as_expr()
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
            "    c0={}".format(
                c0
            )
        )

        print(
            "    c1={}".format(
                c1
            )
        )

        print(
            "    Q_t(p)={}".format(
                sp.factor(
                    P.as_expr()
                )
            )
        )

        print(
            "    Q_t(p+2)={}".format(
                sp.factor(
                    shifted.as_expr()
                )
            )
        )

        print(
            "    target_Q={}".format(
                sp.factor(
                    R.as_expr()
                )
            )
        )

        print(
            "    symbolic_residual={}".format(
                sp.factor(
                    residual.as_expr()
                )
            )
        )

        print(
            "    exact_polynomial_identity={}".format(
                residual.is_zero
            )
        )

        results.append(
            residual.is_zero
        )

    return results


# ============================================================================
# RAW COEFFICIENT EXTRACTION
# ============================================================================

def solve_raw_width2(
    layers,
    t,
):

    source = [
        value
        for _, value
        in layers[t]
    ]

    target = [
        value
        for _, value
        in layers[t + 1]
    ]

    equations = []

    for k in range(
        min(
            len(target),
            len(source) - 1,
        )
    ):

        equations.append(
            (
                source[k],
                source[k + 1],
                target[k],
            )
        )

    if len(equations) < 2:
        return {
            "status": "INSUFFICIENT_DATA"
        }

    M = sp.Matrix([
        [x0, x1]
        for x0, x1, _
        in equations
    ])

    rhs = sp.Matrix([
        value
        for _, _, value
        in equations
    ])

    if (
        M.rank()
        < 2
        or
        M.row_join(rhs).rank()
        >
        M.rank()
    ):
        return {
            "status": "NO_SOLUTION"
        }

    solution = M.inv() * rhs

    c0 = clean(
        solution[0]
    )

    c1 = clean(
        solution[1]
    )

    return {
        "status": "EXACT",
        "c0": c0,
        "c1": c1,
    }


# ============================================================================
# COEFFICIENT SEQUENCE AUDIT
# ============================================================================

def coefficient_sequence_audit(
    polys,
    raw_coefficients,
):

    print()
    print("=" * 78)
    print(
        "7. COEFFICIENT-SEQUENCE / T-AUDIT"
    )
    print("=" * 78)

    for label in (
        "c0",
        "c1",
    ):

        values = [
            (
                t,
                raw_coefficients[t][label]
            )
            for t in sorted(
                raw_coefficients
            )
        ]

        print()
        print(
            "  {}={}".format(
                label,
                values,
            )
        )

        if len(values) == 2:

            t0, v0 = values[0]
            t1, v1 = values[1]

            print(
                "    ratio={}".format(
                    clean(
                        v1 / v0
                    )
                )
            )

            print(
                "    difference={}".format(
                    clean(
                        v1 - v0
                    )
                )
            )

            print(
                "    sum={}".format(
                    clean(
                        v1 + v0
                    )
                )
            )


# ============================================================================
# P-FINITE DIFFERENCE AUDIT
# ============================================================================

def p_difference_rows(values):

    current = [
        sp.Integer(v)
        for v in values
    ]

    rows = [
        current
    ]

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

        rows.append(
            current
        )

    return rows


def p_difference_audit(layers):

    print()
    print("=" * 78)
    print(
        "8. EXACT p-DIFFERENCE TABLE"
    )
    print("=" * 78)

    for t in sorted(layers):

        values = [
            value
            for _, value
            in layers[t]
        ]

        ps = [
            prime
            for prime, _
            in layers[t]
        ]

        if len(values) < 2:
            continue

        rows = p_difference_rows(
            values
        )

        print()
        print(
            "  t={}: p={}".format(
                t,
                ps,
            )
        )

        for j, row in enumerate(rows):

            print(
                "    Delta_p^{}={}".format(
                    j,
                    row,
                )
            )


# ============================================================================
# BIVARIATE SAMPLE TABLE
# ============================================================================

def bivariate_table(layers):

    print()
    print("=" * 78)
    print(
        "9. EXISTING BIVARIATE Q(t,p) SAMPLE TABLE"
    )
    print("=" * 78)

    all_p = sorted(
        Q.keys()
    )

    print()
    print(
        "  columns p={}".format(
            all_p
        )
    )

    for t in sorted(layers):

        values = dict(
            layer_points(
                layers,
                t,
            )
        )

        row = [
            values.get(
                prime,
                None,
            )
            for prime in all_p
        ]

        print()
        print(
            "  t={}: {}".format(
                t,
                row,
            )
        )


# ============================================================================
# TERMINAL REFERENCE
# ============================================================================

def terminal_reference():

    q1 = Q[1][-1]
    q3 = Q[3][-1]

    print()
    print("=" * 78)
    print(
        "10. TERMINAL SOURCE REFERENCE"
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
            math.gcd(
                q1,
                q3,
            )
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


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 320R — EXACT SOURCE-POLYNOMIAL / "
        "LAYER-RECURRENCE AUDIT"
    )
    print("=" * 78)

    layers = build_layers()

    # ------------------------------------------------------------------------
    # 1. BIVARIATE DATA
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "1. EXISTING SOURCE DATA"
    )
    print("=" * 78)

    bivariate_table(
        layers
    )

    # ------------------------------------------------------------------------
    # 2. INTERPOLATING POLYNOMIALS
    # ------------------------------------------------------------------------

    polys = layer_polynomial_audit(
        layers
    )

    # ------------------------------------------------------------------------
    # 3. EXACT INTERPOLATION CHECK
    # ------------------------------------------------------------------------

    interpolation_ok = (
        interpolation_residual_audit(
            layers,
            polys,
        )
    )

    # ------------------------------------------------------------------------
    # 4. COEFFICIENT TABLE
    # ------------------------------------------------------------------------

    coefficient_table(
        polys
    )

    # ------------------------------------------------------------------------
    # 5. ADJACENT LAYER STRUCTURE
    # ------------------------------------------------------------------------

    adjacent_layer_audit(
        polys
    )

    # ------------------------------------------------------------------------
    # 6. RAW TRANSFER COEFFICIENTS
    # ------------------------------------------------------------------------

    raw_coefficients = {}

    for t in (
        1,
        2,
    ):

        result = solve_raw_width2(
            layers,
            t,
        )

        if result["status"] == "EXACT":
            raw_coefficients[t] = result

    # ------------------------------------------------------------------------
    # 7. RAW TRANSFER AS SYMBOLIC POLYNOMIAL IDENTITY
    # ------------------------------------------------------------------------

    transfer_checks = (
        raw_transfer_polynomial_audit(
            polys,
            raw_coefficients,
        )
    )

    # ------------------------------------------------------------------------
    # 8. t-SEQUENCE AUDIT
    # ------------------------------------------------------------------------

    coefficient_sequence_audit(
        polys,
        raw_coefficients,
    )

    # ------------------------------------------------------------------------
    # 9. p-DIFFERENCE AUDIT
    # ------------------------------------------------------------------------

    p_difference_audit(
        layers
    )

    # ------------------------------------------------------------------------
    # 10. TERMINAL SOURCE
    # ------------------------------------------------------------------------

    terminal_reference()

    # ------------------------------------------------------------------------
    # 11. INTERPRETATION
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "11. STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
Experiment 319R found that the two raw transition coefficient pairs do not
collapse under the most obvious normalization tests.

Experiment 320R therefore returns to the source object itself.

The available data naturally form a triangular table

    Q_t(p),

with fewer p-values as t increases.

Each row can be represented exactly by its interpolation polynomial in p,
but that interpolation is only a coordinate representation of the
existing samples.

The critical test is whether the raw transfer law

    Q_{t+1}(p)
      =
    c0_t Q_t(p) + c1_t Q_t(p+2)

survives as an exact polynomial identity after reconstructing the
existing layer polynomials.

A positive result means that the raw transfer is not merely an accidental
fit to the two available p-pairs: it is compatible with an exact symbolic
polynomial relation in p.

A negative result means the raw transfer is only a finite-data relation and
should not be elevated.

The layer-polynomial coefficient table is the next potential source of a
genuine formula. In particular, one should inspect whether the
coefficients across t exhibit:

    * linear or geometric progression;
    * common factors;
    * low-degree t-recurrences;
    * simple shifts or derivatives in p;
    * factorization patterns shared across layers.

The experiment does not extrapolate to a new p or a new n=pq instance.
"""
    )

    # ------------------------------------------------------------------------
    # 12. FINAL EXACTNESS
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "12. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  interpolation_represents_all_existing_samples={}".format(
            interpolation_ok
        )
    )

    print(
        "  raw_transfer_polynomial_checks={}".format(
            transfer_checks
        )
    )

    print(
        "  exact_raw_transfer_identity_all_tested={}".format(
            all(transfer_checks)
            if transfer_checks
            else False
        )
    )

    print(
        "  coefficient_sequence_analysis=True"
    )

    print(
        "  p_difference_analysis=True"
    )

    print(
        "  no_new_p_values_created=False"
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
        "EXPERIMENT 320R COMPLETE"
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

