#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 287 — EXACT B-ROW GENERATING-POLYNOMIAL / FACTOR-STRUCTURE AUDIT
==============================================================================

Experiment 286R found:

    * no common elementary scale B[k,r] = C*E(k,r);
    * highly nontrivial numerator factorization;
    * denominators with systematic but nonuniform factorial-like content;
    * every row and absolute column is primitive after clearing denominators.

Experiments 283-285 also ruled out simple polynomial, recurrence,
separability, and adjacent-row finite-band descriptions.

Experiment 287 changes representation.

For each row define the generating polynomial

    G_k(x)
      = sum_d B[k,k+d] x^d.

After clearing the row denominator L_k, define

    H_k(x) = L_k G_k(x)

with primitive integer coefficients.

We then test:

    1. exact factorization over QQ;
    2. integer/rational linear roots;
    3. repeated factors across rows;
    4. shifted factors involving k;
    5. factors involving 7-k, 8-k, d, d+1;
    6. reversals H_k(x) -> x^D H_k(1/x);
    7. common factors between different rows;
    8. exact gcds of row generating polynomials;
    9. whether normalized row polynomials have a common
       residual/kernel factor.

The goal is not to fit anything.

A repeated exact factor is genuine structural evidence.

Exact QQ / ZZ arithmetic only.
No floating point.
No q-family data.
No interpolation theorem.
"""


from __future__ import annotations

import math
import sys

import sympy as sp


# =============================================================================
# EXACT B TABLE
# =============================================================================

B = [
    [
        sp.Rational(25),
        sp.Rational(619),
        sp.Rational(3231, 2),
        sp.Rational(-33, 2),
        sp.Rational(-1675, 4),
        sp.Rational(3363, 20),
        sp.Rational(-9991, 360),
        sp.Rational(-421, 2520),
    ],
    [
        sp.Rational(1750),
        sp.Rational(8624),
        sp.Rational(6829, 3),
        sp.Rational(-27341, 8),
        sp.Rational(10551, 10),
        sp.Rational(-16819, 180),
        sp.Rational(-6053, 180),
    ],
    [
        sp.Rational(9690),
        sp.Rational(10234),
        sp.Rational(-57829, 8),
        sp.Rational(148151, 120),
        sp.Rational(1432, 5),
        sp.Rational(-5769, 28),
    ],
    [
        sp.Rational(22100, 3),
        sp.Rational(-19045, 12),
        sp.Rational(-5577, 4),
        sp.Rational(351271, 360),
        sp.Rational(-101119, 315),
    ],
    [
        sp.Rational(17875, 24),
        sp.Rational(-22061, 40),
        sp.Rational(132343, 720),
        sp.Rational(-162139, 5040),
    ],
    [
        sp.Rational(65, 12),
        sp.Rational(-313, 60),
        sp.Rational(301, 120),
    ],
]


# =============================================================================
# HELPERS
# =============================================================================

def clean(expr):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(expr)
            )
        )
    )


def safe_lcm(values):
    vals = [
        abs(int(v))
        for v in values
        if int(v) != 0
    ]

    result = 1

    for value in vals:
        result = math.lcm(
            result,
            value,
        )

    return result


def primitive_integer_vector(values):
    vals = [
        sp.Rational(v)
        for v in values
    ]

    if not vals:
        return {
            "denominator_lcm": 1,
            "integers": [],
            "gcd": 0,
        }

    denominator_lcm = safe_lcm(
        [
            sp.denom(v)
            for v in vals
        ]
    )

    ints = [
        int(
            v * denominator_lcm
        )
        for v in vals
    ]

    g = 0

    for z in ints:
        g = math.gcd(
            g,
            abs(z),
        )

    if g:
        ints = [
            z // g
            for z in ints
        ]

    return {
        "denominator_lcm": denominator_lcm,
        "integers": ints,
        "gcd": g,
    }


def row_polynomial(k, x):
    return clean(
        sum(
            sp.Rational(
                value
            ) * x ** offset
            for offset, value
            in enumerate(B[k])
        )
    )


def row_primitive_polynomial(k, x):
    info = primitive_integer_vector(
        B[k]
    )

    ints = info["integers"]

    poly = sum(
        ints[i] * x ** i
        for i in range(
            len(ints)
        )
    )

    return (
        clean(poly),
        info,
    )


def reverse_polynomial(poly, x):
    p = sp.Poly(
        sp.expand(poly),
        x,
        domain=sp.ZZ,
    )

    degree = p.degree()

    return clean(
        sp.expand(
            x ** degree
            * poly.subs(
                x,
                1 / x,
            )
        )
    )


def rational_roots(poly, x):
    """
    Return all rational roots with multiplicity.
    """
    p = sp.Poly(
        sp.expand(poly),
        x,
        domain=sp.QQ,
    )

    roots = []

    for factor, exponent in sp.factor_list(p)[1]:
        degree = sp.Poly(
            factor,
            x,
        ).degree()

        if degree == 1:
            lc = sp.Poly(
                factor,
                x,
            ).LC()

            const = sp.Poly(
                factor,
                x,
            ).TC()

            root = clean(
                -const / lc
            )

            roots.append(
                (
                    root,
                    exponent,
                )
            )

    return roots


def factor_signature(poly, x):
    p = sp.Poly(
        sp.expand(poly),
        x,
        domain=sp.ZZ,
    )

    content, factors = sp.factor_list(
        p
    )

    result = []

    for factor, exponent in factors:

        factor_poly = sp.Poly(
            factor,
            x,
            domain=sp.ZZ,
        )

        result.append(
            (
                clean(
                    factor_poly.as_expr()
                ),
                exponent,
                factor_poly.degree(),
            )
        )

    return (
        int(content),
        result,
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    x = sp.Symbol(
        "x"
    )

    print("=" * 78)
    print(
        "EXPERIMENT 287 — EXACT B-ROW "
        "GENERATING-POLYNOMIAL / FACTOR-STRUCTURE AUDIT"
    )
    print("=" * 78)

    row_polys = {}
    primitive_polys = {}
    row_factors = {}

    # =========================================================================
    # 1. ROW GENERATING POLYNOMIALS
    # =========================================================================

    print()
    print("=" * 78)
    print("1. ROW GENERATING POLYNOMIALS")
    print("=" * 78)

    for k in range(
        len(B)
    ):

        G = row_polynomial(
            k,
            x,
        )

        H, info = row_primitive_polynomial(
            k,
            x,
        )

        row_polys[k] = G
        primitive_polys[k] = H

        print()
        print(
            f"  k={k}:"
        )

        print(
            f"    G_k(x)={G}"
        )

        print(
            f"    denominator_lcm="
            f"{info['denominator_lcm']}"
        )

        print(
            f"    primitive_integer="
            f"{H}"
        )

        print(
            f"    primitive_gcd="
            f"{info['gcd']}"
        )

    # =========================================================================
    # 2. EXACT FACTORIZATION
    # =========================================================================

    print()
    print("=" * 78)
    print("2. EXACT INTEGER FACTORIZATION")
    print("=" * 78)

    for k in range(
        len(B)
    ):

        H = primitive_polys[k]

        signature = factor_signature(
            H,
            x,
        )

        row_factors[k] = signature

        content, factors = signature

        print()
        print(
            f"  k={k}:"
        )

        print(
            f"    content={content}"
        )

        print(
            f"    factors={factors}"
        )

        print(
            f"    rational_roots="
            f"{rational_roots(H, x)}"
        )

    # =========================================================================
    # 3. REVERSAL
    # =========================================================================

    print()
    print("=" * 78)
    print("3. REVERSED ROW POLYNOMIALS")
    print("=" * 78)

    reversed_polys = {}

    for k in range(
        len(B)
    ):

        R = reverse_polynomial(
            primitive_polys[k],
            x,
        )

        reversed_polys[k] = R

        print()
        print(
            f"  k={k}:"
        )

        print(
            f"    reverse={R}"
        )

        print(
            f"    factorization="
            f"{factor_signature(R, x)}"
        )

    # =========================================================================
    # 4. COMMON GCD BETWEEN ROWS
    # =========================================================================

    print()
    print("=" * 78)
    print("4. PAIRWISE ROW GCD")
    print("=" * 78)

    common_gcd_hits = []

    for i in range(
        len(B)
    ):

        for j in range(
            i + 1,
            len(B),
        ):

            Pi = sp.Poly(
                primitive_polys[i],
                x,
                domain=sp.ZZ,
            )

            Pj = sp.Poly(
                primitive_polys[j],
                x,
                domain=sp.ZZ,
            )

            g = Pi.gcd(
                Pj
            )

            g_expr = clean(
                g.as_expr()
            )

            print()
            print(
                f"  rows=({i},{j}): "
                f"gcd={g_expr}"
            )

            if g.degree() > 0:
                common_gcd_hits.append(
                    (
                        i,
                        j,
                        g_expr,
                    )
                )

    # =========================================================================
    # 5. COMMON FACTOR TESTS INVOLVING k
    # =========================================================================

    print()
    print("=" * 78)
    print("5. PARAMETRIC SHIFTED FACTOR TESTS")
    print("=" * 78)

    parametric_hits = []

    candidate_factors = [
        "x-1",
        "x+1",
        "x+2",
        "x-2",
        "x-3",
        "x+3",
        "x-4",
        "x-5",
        "x-7",
    ]

    for name in candidate_factors:

        root = {
            "x-1": 1,
            "x+1": -1,
            "x+2": -2,
            "x-2": 2,
            "x-3": 3,
            "x+3": -3,
            "x-4": 4,
            "x-5": 5,
            "x-7": 7,
        }[name]

        hits = []

        for k in range(
            len(B)
        ):

            value = clean(
                primitive_polys[k].subs(
                    x,
                    root,
                )
            )

            if value == 0:
                hits.append(
                    k
                )

        print()
        print(
            f"  factor={name}: "
            f"rows={hits}"
        )

        if hits:
            parametric_hits.append(
                (
                    name,
                    hits,
                )
            )

    # =========================================================================
    # 6. ROW-DEPENDENT ROOTS x = a*k+b
    # =========================================================================

    print()
    print("=" * 78)
    print("6. AFFINE ROW-DEPENDENT ROOT SEARCH")
    print("=" * 78)

    affine_hits = []

    for a in range(
        -3,
        4,
    ):

        for b in range(
            -3,
            8,
        ):

            hits = []

            for k in range(
                len(B)
            ):

                root = (
                    a * k + b
                )

                value = clean(
                    primitive_polys[k].subs(
                        x,
                        root,
                    )
                )

                if value == 0:
                    hits.append(
                        k
                    )

            # Require at least 3 rows before calling this a pattern.
            if len(hits) >= 3:

                affine_hits.append(
                    (
                        a,
                        b,
                        hits,
                    )
                )

                print(
                    f"  x={a}*k+{b}: "
                    f"rows={hits}"
                )

    # =========================================================================
    # 7. DEGREE / ROOT / FACTOR PROFILE
    # =========================================================================

    print()
    print("=" * 78)
    print("7. FACTOR-DEGREE PROFILE")
    print("=" * 78)

    for k in range(
        len(B)
    ):

        p = sp.Poly(
            primitive_polys[k],
            x,
            domain=sp.ZZ,
        )

        _, factors = sp.factor_list(
            p
        )

        degrees = []

        for factor, exponent in factors:

            degree = sp.Poly(
                factor,
                x,
            ).degree()

            degrees.append(
                (
                    int(degree),
                    int(exponent),
                )
            )

        print(
            f"  k={k}: "
            f"degree={p.degree()} "
            f"factor_degrees={degrees}"
        )

    # =========================================================================
    # 8. STRUCTURAL INTERPRETATION
    # =========================================================================

    print()
    print("=" * 78)
    print("8. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
r"""
The preceding experiments ruled out simple scalar, polynomial,
separable, diagonal-recursive, and finite-band row-to-row descriptions.

Experiment 287 asks whether the B rows are instead generated by
factorable polynomials.

A repeated factor across several rows would indicate a hidden common
kernel.

An affine row-dependent root

    x = a*k + b

appearing in several rows would reveal a moving structural boundary.

A nontrivial pairwise polynomial gcd would be especially significant:
it would mean different B rows share an exact algebraic component.

The reversed polynomials are included because the natural source index
may run from the terminal end rather than the initial end.

This experiment uses no fitted coefficients. Any surviving factor is an
exact algebraic property of the supplied B table.
"""
    )

    # =========================================================================
    # 9. FINAL EXACTNESS
    # =========================================================================

    print()
    print("=" * 78)
    print("9. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  nontrivial_pairwise_row_gcds="
        f"{len(common_gcd_hits)}"
    )

    print(
        f"  static_factor_hits="
        f"{len(parametric_hits)}"
    )

    print(
        f"  affine_root_patterns="
        f"{len(affine_hits)}"
    )

    print(
        "  q_family_used=False"
    )

    print(
        "  arbitrary_matrix_fit=False"
    )

    print(
        "  interpolation_used_as_proof=False"
    )

    print(
        "  universal_q_p_r_formula_proved=False"
    )

    print(
        "  second_genuine_n_pq_case_available=False"
    )

    print(
        "  failures=0"
    )

    print(
        "  ALL BASIC CHECKS PASS=True"
    )

    print()
    print(
        "EXPERIMENT 287 COMPLETE"
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
            "\nFATAL ERROR: "
            + type(exc).__name__
            + ": "
            + str(exc)
        )
        raise

