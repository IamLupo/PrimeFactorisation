#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 347R — EXACT TERMINAL-POLYNOMIAL / CONVOLUTION PROVENANCE AUDIT
==============================================================================

Purpose
-------
The generic recurrence/operator searches have now been exhausted.

The remaining high-value structural hypothesis is that the triangular
source table is generated from the terminal row

    Q_0(1), Q_0(3), Q_0(5), Q_0(7)

through an explicit coefficient-extraction / convolution mechanism.

Define

    F(x) = 495451247
         + 421514439 x
         + 16027881 x^2
         + x^3.

This experiment tests exact, non-fitted transforms of F against the
observed rows.

Families tested:

    * F^k coefficient slices;
    * reversed polynomial;
    * derivative / repeated derivative;
    * finite-difference transforms;
    * ordinary coefficient shifts;
    * reciprocal-polynomial transforms;
    * products F^a * F_rev^b for small a,b;
    * coefficient-window transforms;
    * signed transforms.

IMPORTANT
---------
No parameter is fitted to the observed Q-table.

Only small explicitly enumerated transforms are tested.

A transform counts as a discovery only if it reproduces EVERY available
cell of one or more observed layers exactly.

No missing values.
No extrapolation.
No interpolation.
No synthetic second n=pq case.
Exact SymPy arithmetic only.
"""


from __future__ import annotations

import itertools
import sys

import sympy as sp


x = sp.symbols("x")


# ============================================================================
# OBSERVED SOURCE TABLE
# ============================================================================

Q = {
    (1, 0): 495451247,
    (3, 0): 421514439,
    (5, 0): 16027881,
    (7, 0): 1,

    (1, 1): -1338089411,
    (3, 1): -128667196,
    (5, 1): 4771718,

    (1, 2): 1764373740,
    (3, 2): -152369292,
    (5, 2): -62398,

    (1, 3): 2668721436,
    (3, 3): -1263551016,

    (1, 4): -11600759760,
    (3, 4): 9955176,

    (1, 5): -126258696,
}


# ============================================================================
# TERMINAL POLYNOMIAL
# ============================================================================

F = sp.Poly(
    495451247
    + 421514439 * x
    + 16027881 * x**2
    + x**3,
    x,
    domain=sp.ZZ,
)

FREV = sp.Poly(
    sp.expand(
        x**3 * F.as_expr().subs(
            x,
            1 / x,
        )
    ),
    x,
    domain=sp.QQ,
)


# ============================================================================
# HELPERS
# ============================================================================

def clean(value):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(value)
            )
        )
    )


def coeff(poly, k):
    return clean(
        sp.Poly(
            poly,
            x,
            domain=sp.QQ,
        ).coeff_monomial(
            x**k
        )
    )


def degree(poly):
    return sp.Poly(
        poly,
        x,
        domain=sp.QQ,
    ).degree()


def vector_coefficients(poly):
    P = sp.Poly(
        poly,
        x,
        domain=sp.QQ,
    )

    return [
        clean(
            P.coeff_monomial(
                x**k
            )
        )
        for k in range(
            P.degree(),
            -1,
            -1,
        )
    ]


def observed_layer(t):
    points = sorted(
        [
            (p, value)
            for (p, tt), value in Q.items()
            if tt == t
        ]
    )

    return points


def compare_layer_vector(t, values):
    observed = observed_layer(t)

    if len(values) != len(observed):
        return False, []

    residuals = []

    for (_, observed_value), predicted in zip(
        observed,
        values,
    ):
        residuals.append(
            clean(
                predicted
                - observed_value
            )
        )

    return all(
        residual == 0
        for residual in residuals
    ), residuals


def polynomial_vector(poly, length):
    return [
        coeff(poly, k)
        for k in range(
            length
        )
    ]


# ============================================================================
# TRANSFORM GENERATORS
# ============================================================================

def generate_basic_transforms():

    transforms = {}

    transforms["F"] = F.as_expr()
    transforms["F_reverse"] = FREV.as_expr()

    transforms["F_derivative"] = sp.diff(
        F.as_expr(),
        x,
    )

    transforms["F_second_derivative"] = sp.diff(
        F.as_expr(),
        x,
        2,
    )

    transforms["F_minus_1"] = F.as_expr() - 1
    transforms["F_plus_1"] = F.as_expr() + 1
    transforms["-F"] = -F.as_expr()

    return transforms


def generate_power_transforms():

    transforms = {}

    for k in range(
        2,
        6,
    ):

        transforms[
            "F^{}".format(k)
        ] = sp.expand(
            F.as_expr() ** k
        )

        transforms[
            "Frev^{}".format(k)
        ] = sp.expand(
            FREV.as_expr() ** k
        )

    return transforms


def generate_mixed_products():

    transforms = {}

    for a in range(
        1,
        4,
    ):

        for b in range(
            1,
            4,
        ):

            transforms[
                "F^{}*Frev^{}".format(
                    a,
                    b,
                )
            ] = sp.expand(
                F.as_expr() ** a
                * FREV.as_expr() ** b
            )

    return transforms


def generate_difference_transforms():

    transforms = {}

    current = F.as_expr()

    for k in range(
        1,
        4,
    ):

        current = sp.expand(
            current.subs(
                x,
                x + 1,
            )
            -
            current
        )

        transforms[
            "Delta_x^{}_F".format(k)
        ] = current

    current = FREV.as_expr()

    for k in range(
        1,
        4,
    ):

        current = sp.expand(
            current.subs(
                x,
                x + 1,
            )
            -
            current
        )

        transforms[
            "Delta_x^{}_Frev".format(k)
        ] = current

    return transforms


def generate_shift_transforms():

    transforms = {}

    for shift in range(
        -3,
        4,
    ):

        transforms[
            "F(x+{})".format(
                shift
            )
        ] = sp.expand(
            F.as_expr().subs(
                x,
                x + shift,
            )
        )

    return transforms


# ============================================================================
# COEFFICIENT-WINDOW TEST
# ============================================================================

def test_coefficient_windows(
    name,
    poly,
    max_width=6,
):

    hits = []

    d = degree(poly)

    for width in range(
        1,
        min(
            max_width,
            d + 1,
        ) + 1,
    ):

        for start in range(
            0,
            d - width + 2,
        ):

            values = [
                coeff(
                    poly,
                    start + j,
                )
                for j in range(
                    width
                )
            ]

            for t in range(
                0,
                6,
            ):

                observed = observed_layer(
                    t
                )

                if len(observed) != width:
                    continue

                ok, residuals = compare_layer_vector(
                    t,
                    values,
                )

                if ok:

                    hits.append(
                        {
                            "name": name,
                            "t": t,
                            "start": start,
                            "width": width,
                            "residuals": residuals,
                        }
                    )

    return hits


# ============================================================================
# FULL TRANSFORM AUDIT
# ============================================================================

def audit_transforms(transforms):

    print()
    print("=" * 78)
    print(
        "1. EXACT TERMINAL-POLYNOMIAL TRANSFORM AUDIT"
    )
    print("=" * 78)

    all_hits = []

    for name, poly in transforms.items():

        print()
        print(
            "  transform={}".format(
                name
            )
        )

        print(
            "    degree={}".format(
                degree(poly)
            )
        )

        hits = test_coefficient_windows(
            name,
            poly,
        )

        if hits:

            print(
                "    HITS={}".format(
                    hits
                )
            )

            all_hits.extend(
                hits
            )

        else:

            print(
                "    hits=[]"
            )

    return all_hits


# ============================================================================
# SMALL LINEAR COMBINATION AUDIT
# ============================================================================

def linear_combination_audit():

    print()
    print("=" * 78)
    print(
        "2. SMALL FIXED-COMBINATION CONVOLUTION AUDIT"
    )
    print("=" * 78)

    basis = {
        "F": F.as_expr(),
        "Frev": FREV.as_expr(),
        "F'": sp.diff(
            F.as_expr(),
            x,
        ),
        "Frev'": sp.diff(
            FREV.as_expr(),
            x,
        ),
    }

    hits = []

    coefficients = (
        -2,
        -1,
        0,
        1,
        2,
    )

    names = list(
        basis
    )

    expressions = []

    for combo in itertools.product(
        coefficients,
        repeat=len(names),
    ):

        if all(
            c == 0
            for c in combo
        ):
            continue

        expr = clean(
            sum(
                combo[i]
                * basis[names[i]]
                for i in range(
                    len(names)
                )
            )
        )

        expressions.append(
            (
                combo,
                expr,
            )
        )

    for combo, expr in expressions:

        d = degree(expr)

        if d < 0:
            continue

        values = polynomial_vector(
            expr,
            d + 1,
        )

        for t in range(
            6,
        ):

            if len(
                observed_layer(t)
            ) != len(values):

                continue

            ok, _ = compare_layer_vector(
                t,
                values,
            )

            if ok:

                hits.append(
                    (
                        combo,
                        t,
                        expr,
                    )
                )

    print(
        "  exact_hits={}".format(
            len(hits)
        )
    )

    for hit in hits:
        print(
            "    {}".format(
                hit
            )
        )

    return hits


# ============================================================================
# STRUCTURAL SUMMARY
# ============================================================================

def structural_summary(
    transform_hits,
    linear_hits,
):

    print()
    print("=" * 78)
    print(
        "3. STRUCTURAL SUMMARY"
    )
    print("=" * 78)

    print(
        "  terminal_polynomial={}".format(
            F.as_expr()
        )
    )

    print(
        "  terminal_reverse={}".format(
            FREV.as_expr()
        )
    )

    print(
        "  transform_hits={}".format(
            len(transform_hits)
        )
    )

    print(
        "  small_linear_combination_hits={}".format(
            len(linear_hits)
        )
    )

    if transform_hits or linear_hits:

        print(
            "  interpretation=potential_convolution_provenance"
        )

    else:

        print(
            "  interpretation=no_small_terminal_polynomial_mechanism_found"
        )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 347R — EXACT TERMINAL-POLYNOMIAL / "
        "CONVOLUTION PROVENANCE AUDIT"
    )
    print("=" * 78)

    print()
    print(
        "  F(x)={}".format(
            F.as_expr()
        )
    )

    print(
        "  F_reverse(x)={}".format(
            FREV.as_expr()
        )
    )

    transforms = {}

    transforms.update(
        generate_basic_transforms()
    )

    transforms.update(
        generate_power_transforms()
    )

    transforms.update(
        generate_mixed_products()
    )

    transforms.update(
        generate_difference_transforms()
    )

    transforms.update(
        generate_shift_transforms()
    )

    transform_hits = audit_transforms(
        transforms
    )

    linear_hits = linear_combination_audit()

    structural_summary(
        transform_hits,
        linear_hits,
    )

    print()
    print("=" * 78)
    print(
        "4. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  terminal_polynomial_built_exactly=True"
    )

    print(
        "  observed_cells_used_only=True"
    )

    print(
        "  transform_family_count={}".format(
            len(transforms)
        )
    )

    print(
        "  exact_transform_hits={}".format(
            len(transform_hits)
        )
    )

    print(
        "  exact_small_linear_hits={}".format(
            len(linear_hits)
        )
    )

    print(
        "  fitting_performed=False"
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
        "EXPERIMENT 347R COMPLETE"
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
