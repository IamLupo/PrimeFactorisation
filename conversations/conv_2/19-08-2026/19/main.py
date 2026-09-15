#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 335R — EXACT ADJACENT ROW-POLYNOMIAL GCD / SHIFT /
                  DIVISIBILITY / BOUNDARY-FACTOR AUDIT
==============================================================================

Purpose
-------
Experiment 334R found no unique fixed low-order differential operator

    P_{t+1} = sum a_k(p) D^k P_t

for any tested derivative/order degree <= 2 model.

The most conspicuous remaining feature of the source table is not a
recurrence but the shrinking p-support:

    t=0 : p=1,3,5,7
    t=1 : p=1,3,5
    t=2 : p=1,3,5
    t=3 : p=1,3
    t=4 : p=1,3
    t=5 : p=1

Experiment 335R asks whether this boundary structure has an exact
algebraic signature in the row polynomials.

For each observed row polynomial P_t(p), test:

    1. factorization over Q;
    2. squarefreeness;
    3. gcd(P_t, P_{t+1});
    4. divisibility P_{t+1} | P_t and P_t | P_{t+1};
    5. polynomial remainders;
    6. resultants;
    7. shifted gcds:
           gcd(P_t(p), P_{t+1}(p+s))
       for small integer shifts s;
    8. shifted divisibility;
    9. whether a small linear boundary factor divides any row;
   10. exact degree-drop and leading-coefficient signatures.

This is deliberately NOT interpreted as a universal law.

The row polynomials are reconstructed from the observed cells only and
are treated as exact coordinate representatives of those observed rows.

No missing value is inserted.
No new p-value is used as data.
No synthetic second n=pq case.
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


p = sp.symbols("p")


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


def primitive_integer_poly(expr):
    poly = sp.Poly(
        sp.together(expr),
        p,
        domain=sp.QQ,
    )

    if poly.is_zero:
        return sp.Poly(
            0,
            p,
            domain=sp.ZZ,
        )

    coeffs = poly.all_coeffs()

    denominator_lcm = 1

    for c in coeffs:
        denominator_lcm = sp.ilcm(
            denominator_lcm,
            int(sp.denom(c)),
        )

    ints = [
        int(c * denominator_lcm)
        for c in coeffs
    ]

    g = 0

    for x in ints:
        g = math.gcd(
            g,
            abs(x),
        )

    ints = [
        x // g
        for x in ints
    ]

    if ints[0] < 0:
        ints = [
            -x
            for x in ints
        ]

    expr_int = sum(
        sp.Integer(c)
        * p ** (
            len(ints) - 1 - i
        )
        for i, c in enumerate(ints)
    )

    return sp.Poly(
        sp.expand(expr_int),
        p,
        domain=sp.ZZ,
    )


def valuation(x, prime):
    x = sp.Rational(x)

    if x == 0:
        return sp.oo

    n = abs(int(x.p))
    d = abs(int(x.q))

    vn = 0
    vd = 0

    while n % prime == 0:
        n //= prime
        vn += 1

    while d % prime == 0:
        d //= prime
        vd += 1

    return vn - vd


# ============================================================================
# OBSERVED LAYERS
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
                        sp.Integer(p_value),
                        sp.Integer(values[index]),
                    )
                )

        layers[t] = layer

    return layers


# ============================================================================
# ROW POLYNOMIALS
# ============================================================================

def row_polynomial(layer):

    return clean(
        sp.interpolate(
            layer,
            p,
        )
    )


def build_row_polynomials(layers):

    return {
        t: row_polynomial(layer)
        for t, layer in layers.items()
    }


# ============================================================================
# BASIC ROW STRUCTURE
# ============================================================================

def factorization_audit(P):

    poly = sp.Poly(
        P,
        p,
        domain=sp.QQ,
    )

    content, factors = sp.factor_list(
        P,
        p,
    )

    return {
        "degree": sp.degree(P, p),
        "content": clean(content),
        "factors": factors,
        "squarefree": (
            sp.gcd(
                poly,
                poly.diff(),
            ).degree()
            == 0
        ),
    }


# ============================================================================
# ADJACENT GCD / DIVISIBILITY
# ============================================================================

def adjacent_pair_audit(
    P,
    t,
):

    A = P[t]
    B = P[t + 1]

    poly_A = sp.Poly(
        A,
        p,
        domain=sp.QQ,
    )

    poly_B = sp.Poly(
        B,
        p,
        domain=sp.QQ,
    )

    gcd_poly = sp.gcd(
        poly_A,
        poly_B,
    )

    rem_A_mod_B = sp.rem(
        poly_A,
        poly_B,
    )

    rem_B_mod_A = sp.rem(
        poly_B,
        poly_A,
    )

    resultant = sp.resultant(
        A,
        B,
        p,
    )

    return {
        "gcd": clean(
            gcd_poly.as_expr()
        ),
        "gcd_degree": gcd_poly.degree(),
        "A_div_B": (
            rem_A_mod_B.is_zero
        ),
        "B_div_A": (
            rem_B_mod_A.is_zero
        ),
        "rem_A_mod_B": clean(
            rem_A_mod_B.as_expr()
        ),
        "rem_B_mod_A": clean(
            rem_B_mod_A.as_expr()
        ),
        "resultant": clean(
            resultant
        ),
    }


# ============================================================================
# SHIFTED GCD / DIVISIBILITY
# ============================================================================

def shifted_pair_audit(
    A,
    B,
    shift,
):

    shifted_B = clean(
        B.subs(
            p,
            p + shift,
        )
    )

    poly_A = sp.Poly(
        A,
        p,
        domain=sp.QQ,
    )

    poly_B = sp.Poly(
        shifted_B,
        p,
        domain=sp.QQ,
    )

    gcd_poly = sp.gcd(
        poly_A,
        poly_B,
    )

    rem_A_mod_B = sp.rem(
        poly_A,
        poly_B,
    )

    rem_B_mod_A = sp.rem(
        poly_B,
        poly_A,
    )

    resultant = sp.resultant(
        A,
        shifted_B,
        p,
    )

    return {
        "shifted_B": shifted_B,
        "gcd": clean(
            gcd_poly.as_expr()
        ),
        "gcd_degree": gcd_poly.degree(),
        "A_div_shifted_B": (
            rem_A_mod_B.is_zero
        ),
        "shifted_B_div_A": (
            rem_B_mod_A.is_zero
        ),
        "resultant": clean(
            resultant
        ),
    }


# ============================================================================
# BOUNDARY LINEAR-FACTOR AUDIT
# ============================================================================

def linear_factor_audit(
    P,
    t,
):

    candidates = {
        "p-(7-2t)": p - (
            7 - 2*t
        ),
        "p-(5-2t)": p - (
            5 - 2*t
        ),
        "p-(3-2t)": p - (
            3 - 2*t
        ),
        "p-1": p - 1,
        "p-3": p - 3,
        "p-5": p - 5,
        "p-7": p - 7,
    }

    result = {}

    poly = sp.Poly(
        P,
        p,
        domain=sp.QQ,
    )

    for name, factor in candidates.items():

        f = sp.Poly(
            factor,
            p,
            domain=sp.QQ,
        )

        result[name] = (
            sp.rem(
                poly,
                f,
            ).is_zero
        )

    return result


# ============================================================================
# DEGREE-DROP SIGNATURE
# ============================================================================

def degree_signature(P):

    rows = []

    for t in sorted(P):

        poly = sp.Poly(
            P[t],
            p,
            domain=sp.QQ,
        )

        degree = poly.degree()
        leading = clean(
            poly.LC()
        )
        constant = clean(
            poly.eval(0)
        )

        rows.append(
            (
                t,
                degree,
                leading,
                constant,
            )
        )

    return rows


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 335R — EXACT ADJACENT ROW-POLYNOMIAL "
        "GCD / SHIFT / DIVISIBILITY / BOUNDARY-FACTOR AUDIT"
    )
    print("=" * 78)

    layers = build_layers()

    P = build_row_polynomials(
        layers
    )

    # ------------------------------------------------------------------------
    # 1. Row polynomial inventory
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. ROW POLYNOMIAL INVENTORY")
    print("=" * 78)

    for t in sorted(P):

        info = factorization_audit(
            P[t]
        )

        print()
        print(
            "  t={}: P_t(p)={}".format(
                t,
                P[t],
            )
        )

        print(
            "    degree={}".format(
                info["degree"]
            )
        )

        print(
            "    content={}".format(
                info["content"]
            )
        )

        print(
            "    factors={}".format(
                info["factors"]
            )
        )

        print(
            "    squarefree={}".format(
                info["squarefree"]
            )
        )

        print(
            "    primitive_integer_polynomial={}".format(
                primitive_integer_poly(
                    P[t]
                ).as_expr()
            )
        )

    # ------------------------------------------------------------------------
    # 2. Adjacent gcd/divisibility
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. ADJACENT GCD / DIVISIBILITY AUDIT")
    print("=" * 78)

    adjacent = {}

    for t in range(
        max(P)
    ):

        result = adjacent_pair_audit(
            P,
            t,
        )

        adjacent[t] = result

        print()
        print(
            "  t={} -> {}:".format(
                t,
                t + 1,
            )
        )

        print(
            "    gcd={}".format(
                result["gcd"]
            )
        )

        print(
            "    gcd_degree={}".format(
                result["gcd_degree"]
            )
        )

        print(
            "    P_t_divides_P_next={}".format(
                result["B_div_A"]
            )
        )

        print(
            "    P_next_divides_P_t={}".format(
                result["A_div_B"]
            )
        )

        print(
            "    remainder_P_t_mod_P_next={}".format(
                result["rem_A_mod_B"]
            )
        )

        print(
            "    remainder_P_next_mod_P_t={}".format(
                result["rem_B_mod_A"]
            )
        )

        print(
            "    resultant={}".format(
                result["resultant"]
            )
        )

    # ------------------------------------------------------------------------
    # 3. Small integer shift search
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "3. SMALL INTEGER SHIFT GCD / DIVISIBILITY SEARCH"
    )
    print("=" * 78)

    shift_hits = []

    for t in range(
        max(P)
    ):

        A = P[t]
        B = P[t + 1]

        print()
        print(
            "  t={} -> {}:".format(
                t,
                t + 1,
            )
        )

        for shift in range(
            -6,
            7,
        ):

            result = shifted_pair_audit(
                A,
                B,
                shift,
            )

            if (
                result["gcd_degree"] > 0
                or result["A_div_shifted_B"]
                or result["shifted_B_div_A"]
            ):

                shift_hits.append(
                    (
                        t,
                        shift,
                        result,
                    )
                )

                print()
                print(
                    "    SHIFT={}:".format(
                        shift
                    )
                )

                print(
                    "      shifted_P_next={}".format(
                        result["shifted_B"]
                    )
                )

                print(
                    "      gcd={}".format(
                        result["gcd"]
                    )
                )

                print(
                    "      gcd_degree={}".format(
                        result["gcd_degree"]
                    )
                )

                print(
                    "      P_t_divides_shifted_P_next={}".format(
                        result["shifted_B_div_A"]
                    )
                )

                print(
                    "      shifted_P_next_divides_P_t={}".format(
                        result["A_div_shifted_B"]
                    )
                )

    print()
    print(
        "  total_shift_hits={}".format(
            len(shift_hits)
        )
    )

    # ------------------------------------------------------------------------
    # 4. Boundary factor audit
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "4. BOUNDARY LINEAR-FACTOR AUDIT"
    )
    print("=" * 78)

    boundary_hits = []

    for t in sorted(P):

        result = linear_factor_audit(
            P[t],
            t,
        )

        print()
        print(
            "  t={}:".format(
                t
            )
        )

        for name, divides in result.items():

            if divides:

                boundary_hits.append(
                    (
                        t,
                        name,
                    )
                )

            print(
                "    {} divides={}".format(
                    name,
                    divides,
                )
            )

    print()
    print(
        "  boundary_factor_hits={}".format(
            boundary_hits
        )
    )

    # ------------------------------------------------------------------------
    # 5. Degree-drop signature
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "5. DEGREE-DROP / LEADING-COEFFICIENT SIGNATURE"
    )
    print("=" * 78)

    signature = degree_signature(
        P
    )

    for t, degree, leading, constant in signature:

        print()
        print(
            "  t={}:".format(
                t
            )
        )

        print(
            "    degree={}".format(
                degree
            )
        )

        print(
            "    leading_coefficient={}".format(
                leading
            )
        )

        print(
            "    constant_term={}".format(
                constant
            )
        )

        print(
            "    leading_valuations={}".format(
                {
                    prime: valuation(
                        leading,
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

    # ------------------------------------------------------------------------
    # 6. Structural interpretation
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "6. STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
Experiment 334R ruled out a broad class of fixed differential operators
acting on the reconstructed p-polynomial rows.

Experiment 335R asks a different question:

    Are consecutive row-polynomials algebraically related at all?

The three strongest simple possibilities are:

    common factors,
    shifts,
    divisibility.

A nontrivial gcd between P_t and P_{t+1} would indicate a persistent
algebraic factor in the source profile.

A nontrivial shifted gcd

    gcd(P_t(p), P_{t+1}(p+s))

would indicate translation structure in the p-coordinate.

Divisibility would be even stronger and could reveal a boundary
factorization mechanism.

The linear-factor audit specifically tests whether the shrinking support
boundary corresponds to factors such as

    p-(7-2t).

A complete absence of these signatures would argue that the changing
support is a data-availability boundary rather than an algebraic factor
of the underlying source law.

This is still diagnostic because each P_t is reconstructed only from
the observed row.

No unobserved value is promoted to data.
No synthetic second n=pq instance is introduced.
"""
    )

    # ------------------------------------------------------------------------
    # 7. Final exactness
    # ------------------------------------------------------------------------

    nontrivial_adjacent_gcd = [
        t
        for t, result
        in adjacent.items()
        if result["gcd_degree"] > 0
    ]

    any_divisibility = any(
        result["A_div_B"]
        or result["B_div_A"]
        for result in adjacent.values()
    )

    print()
    print("=" * 78)
    print(
        "FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  adjacent_pairs_tested={}".format(
            len(adjacent)
        )
    )

    print(
        "  nontrivial_adjacent_gcd_pairs={}".format(
            nontrivial_adjacent_gcd
        )
    )

    print(
        "  any_adjacent_divisibility={}".format(
            any_divisibility
        )
    )

    print(
        "  shift_hits={}".format(
            len(shift_hits)
        )
    )

    print(
        "  boundary_factor_hits={}".format(
            len(boundary_hits)
        )
    )

    print(
        "  observed_row_polynomials_only=True"
    )

    print(
        "  missing_values_used=False"
    )

    print(
        "  extrapolation_used=False"
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
        "EXPERIMENT 335R COMPLETE"
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
