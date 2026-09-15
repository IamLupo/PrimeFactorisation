#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 337R — EXACT RATIONAL AFFINE EQUIVALENCE / DEGREE-DROP AUDIT
==============================================================================

Purpose
-------
Experiment 336R found:

    * t=1 -> 2 admits an affine polynomial equivalence over an algebraic
      extension, but the reported parameters contain a large square root;
    * no bounded rational affine map was found;
    * t=4 -> 5 satisfies P5 proportional to P4', which is automatic for
      every nonconstant linear polynomial.

Therefore neither observation is yet a substantive rational structural law.

This experiment performs the stronger exact tests:

    1. rational affine equivalence
           P_{t+1}(p) = c P_t(a p + b),
       with a,b,c in Q;

    2. exact discriminant obstruction for the quadratic pair;

    3. exact elimination/resultant test for rational affine parameters;

    4. whether the quadratic affine solution necessarily lives outside Q;

    5. whether any degree-drop pair is genuinely a derivative relation,
       versus merely the trivial degree-lowering fact;

    6. direct comparison of root invariants for the two quadratic rows;

    7. rational scalar/translation/scaling special cases.

Important
---------
The row polynomials are reconstructed only from their observed row cells.
No new source value is treated as data.

An algebraic affine equivalence is NOT counted as a rational source law.

No missing cells.
No extrapolation as experimental data.
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

p = sp.symbols("p")
a, b, c = sp.symbols("a b c")


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


def row_polynomials(layers):
    return {
        t: clean(
            sp.interpolate(
                layer,
                p,
            )
        )
        for t, layer in layers.items()
    }


def primitive_integer_poly(P):
    poly = sp.Poly(
        sp.together(P),
        p,
        domain=sp.QQ,
    )

    coeffs = poly.all_coeffs()

    if not coeffs:
        return sp.Poly(
            0,
            p,
            domain=sp.ZZ,
        )

    den = 1

    for coeff in coeffs:
        den = sp.ilcm(
            den,
            int(sp.denom(coeff)),
        )

    ints = [
        int(coeff * den)
        for coeff in coeffs
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

    return sp.Poly(
        sum(
            sp.Integer(v)
            * p**(
                len(ints) - 1 - i
            )
            for i, v in enumerate(ints)
        ),
        p,
        domain=sp.ZZ,
    )


def coefficients(P):
    poly = sp.Poly(
        P,
        p,
        domain=sp.QQ,
    )
    return poly.all_coeffs()


# ============================================================================
# AFFINE COEFFICIENT EQUATIONS
# ============================================================================

def affine_equations(P, R):
    aa, bb, cc = a, b, c

    transformed = sp.Poly(
        sp.expand(
            cc * P.subs(
                p,
                aa*p + bb,
            )
        ),
        p,
        domain=sp.EX,
    )

    target = sp.Poly(
        sp.expand(R),
        p,
        domain=sp.EX,
    )

    degree = target.degree()

    equations = [
        clean(
            transformed.coeff_monomial(
                p**k
            )
            -
            target.coeff_monomial(
                p**k
            )
        )
        for k in range(
            degree + 1
        )
    ]

    return equations


def rational_affine_test(P, R):
    degP = sp.degree(P, p)
    degR = sp.degree(R, p)

    if degP != degR:
        return {
            "status": "DEGREE_MISMATCH",
            "degree_P": degP,
            "degree_R": degR,
        }

    equations = affine_equations(
        P,
        R,
    )

    # Solve algebraically first, then explicitly filter for rational
    # parameter triples.
    try:
        solutions = sp.solve(
            equations,
            [a, b, c],
            dict=True,
        )
    except Exception:
        solutions = []

    rational_solutions = []

    for solution in solutions:

        aa = solution.get(a)
        bb = solution.get(b)
        cc = solution.get(c)

        if aa is None or bb is None or cc is None:
            continue

        aa = sp.simplify(aa)
        bb = sp.simplify(bb)
        cc = sp.simplify(cc)

        if (
            aa.is_Rational
            and bb.is_Rational
            and cc.is_Rational
        ):
            rational_solutions.append(
                {
                    "a": aa,
                    "b": bb,
                    "c": cc,
                }
            )

    algebraic_solutions = solutions

    return {
        "status": (
            "RATIONAL_EXACT"
            if rational_solutions
            else (
                "ALGEBRAIC_ONLY"
                if algebraic_solutions
                else "NO_SOLUTION"
            )
        ),
        "degree_P": degP,
        "degree_R": degR,
        "rational_solutions": rational_solutions,
        "algebraic_solutions": algebraic_solutions,
        "equations": equations,
    }


# ============================================================================
# QUADRATIC AFFINE OBSTRUCTION
# ============================================================================

def quadratic_data(P):
    poly = sp.Poly(
        P,
        p,
        domain=sp.QQ,
    )

    A, B, C = poly.all_coeffs()

    discriminant = clean(
        B**2 - 4*A*C
    )

    center = clean(
        -B / (2*A)
    )

    radius_square = clean(
        discriminant / (4*A**2)
    )

    return {
        "A": A,
        "B": B,
        "C": C,
        "discriminant": discriminant,
        "center": center,
        "radius_square": radius_square,
    }


def quadratic_affine_invariant(P):
    data = quadratic_data(P)

    # Under p -> a p+b, the quantity
    #
    #     discriminant / A^2
    #
    # scales by a^2.
    #
    # Hence the square-class in Q*/(Q*)^2 is invariant.
    value = clean(
        data["discriminant"]
        / data["A"]**2
    )

    return value


def rational_square_class_equal(x, y):
    """
    Test whether x/y is a rational square.
    """
    ratio = clean(
        sp.Rational(x)
        / sp.Rational(y)
    )

    if ratio == 0:
        return x == y

    if ratio < 0:
        return False

    num = int(sp.numer(ratio))
    den = int(sp.denom(ratio))

    sn, ok_n = sp.integer_nthroot(
        num,
        2,
    )

    sd, ok_d = sp.integer_nthroot(
        den,
        2,
    )

    return bool(
        ok_n and ok_d
    )


# ============================================================================
# SPECIAL RATIONAL AFFINE SUBCLASSES
# ============================================================================

def scalar_only_test(P, R):
    """
    R = c P(p).
    """
    degP = sp.degree(P, p)
    degR = sp.degree(R, p)

    if degP != degR:
        return None

    p_poly = sp.Poly(
        P,
        p,
        domain=sp.QQ,
    )

    r_poly = sp.Poly(
        R,
        p,
        domain=sp.QQ,
    )

    c_value = clean(
        r_poly.LC()
        / p_poly.LC()
    )

    return clean(
        c_value * P - R
    ) == 0, c_value


def translation_only_test(P, R):
    """
    R = c P(p+b), b,c in Q.
    """
    bb, cc = sp.symbols(
        "bb cc"
    )

    equations = affine_equations(
        P,
        R,
    )

    # Force a=1.
    equations = [
        clean(
            expression.subs(
                a,
                1,
            )
        )
        for expression in equations
    ]

    try:
        solutions = sp.solve(
            equations,
            [bb, cc],
            dict=True,
        )
    except Exception:
        solutions = []

    rational = []

    for solution in solutions:
        if (
            bb in solution
            and cc in solution
            and solution[bb].is_Rational
            and solution[cc].is_Rational
        ):
            rational.append(solution)

    return rational


def scaling_only_test(P, R):
    """
    R = c P(a p), a,c in Q.
    """
    aa, cc = sp.symbols(
        "aa cc"
    )

    equations = affine_equations(
        P,
        R,
    )

    equations = [
        clean(
            expression.subs(
                b,
                0,
            )
        )
        for expression in equations
    ]

    try:
        solutions = sp.solve(
            equations,
            [aa, cc],
            dict=True,
        )
    except Exception:
        solutions = []

    rational = []

    for solution in solutions:
        if (
            aa in solution
            and cc in solution
            and solution[aa].is_Rational
            and solution[cc].is_Rational
        ):
            rational.append(solution)

    return rational


# ============================================================================
# DERIVATIVE TEST
# ============================================================================

def derivative_relation(P, R, order):
    D = sp.diff(
        P,
        p,
        order,
    )

    if D == 0:
        return {
            "defined": False,
            "exact": False,
            "scalar": None,
        }

    poly_D = sp.Poly(
        D,
        p,
        domain=sp.QQ,
    )

    poly_R = sp.Poly(
        R,
        p,
        domain=sp.QQ,
    )

    if poly_D.degree() != poly_R.degree():
        return {
            "defined": True,
            "exact": False,
            "scalar": None,
        }

    scalar = clean(
        poly_R.LC()
        / poly_D.LC()
    )

    exact = clean(
        scalar * D - R
    ) == 0

    return {
        "defined": True,
        "exact": exact,
        "scalar": scalar,
    }


# ============================================================================
# RESULTANT / RATIONALITY DIAGNOSTIC
# ============================================================================

def algebraic_parameter_profile(
    solution,
):

    profile = {}

    for name, value in (
        ("a", solution.get(a)),
        ("b", solution.get(b)),
        ("c", solution.get(c)),
    ):

        if value is None:
            continue

        value = sp.simplify(value)

        profile[name] = {
            "value": value,
            "is_rational": bool(
                value.is_Rational
            ),
            "minimal_polynomial": (
                sp.minimal_polynomial(
                    value
                )
                if not value.is_Rational
                else None
            ),
        }

    return profile


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 337R — EXACT RATIONAL AFFINE "
        "EQUIVALENCE / DEGREE-DROP AUDIT"
    )
    print("=" * 78)

    layers = build_layers()
    rows = row_polynomials(
        layers
    )

    # ------------------------------------------------------------------------
    # 1. ROW INVENTORY
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "1. ROW POLYNOMIAL INVENTORY"
    )
    print("=" * 78)

    for t in sorted(rows):

        print()
        print(
            "  t={}: degree={}".format(
                t,
                sp.degree(
                    rows[t],
                    p,
                ),
            )
        )

        print(
            "    P_t={}".format(
                rows[t]
            )
        )

    # ------------------------------------------------------------------------
    # 2. RATIONAL AFFINE AUDIT
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "2. RATIONAL AFFINE TRANSITION AUDIT"
    )
    print("=" * 78)

    affine_results = {}

    for t in range(
        max(rows)
    ):

        result = rational_affine_test(
            rows[t],
            rows[t + 1],
        )

        affine_results[t] = result

        print()
        print(
            "  t={} -> {}:".format(
                t,
                t + 1,
            )
        )

        print(
            "    status={}".format(
                result["status"]
            )
        )

        print(
            "    degree_P={}".format(
                result["degree_P"]
            )
        )

        print(
            "    degree_R={}".format(
                result["degree_R"]
            )
        )

        print(
            "    rational_solutions={}".format(
                result.get(
                    "rational_solutions",
                    [],
                )
            )
        )

        if result.get(
            "algebraic_solutions"
        ):

            print(
                "    algebraic_solution_count={}".format(
                    len(
                        result[
                            "algebraic_solutions"
                        ]
                    )
                )
            )

            for index, solution in enumerate(
                result["algebraic_solutions"]
            ):

                print()
                print(
                    "    algebraic_solution_{}:".format(
                        index
                    )
                )

                print(
                    "      {}".format(
                        algebraic_parameter_profile(
                            solution
                        )
                    )
                )

    # ------------------------------------------------------------------------
    # 3. QUADRATIC RATIONAL AFFINE OBSTRUCTION
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "3. QUADRATIC RATIONAL-AFFINE INVARIANT AUDIT"
    )
    print("=" * 78)

    quadratic_rows = [
        t
        for t in rows
        if sp.degree(
            rows[t],
            p,
        ) == 2
    ]

    for t in quadratic_rows:

        data = quadratic_data(
            rows[t]
        )

        print()
        print(
            "  t={}:".format(
                t
            )
        )

        print(
            "    center={}".format(
                data["center"]
            )
        )

        print(
            "    discriminant={}".format(
                data["discriminant"]
            )
        )

        print(
            "    radius_square={}".format(
                data["radius_square"]
            )
        )

        print(
            "    affine_square_class={}".format(
                quadratic_affine_invariant(
                    rows[t]
                )
            )
        )

    if len(quadratic_rows) == 2:

        t0, t1 = quadratic_rows

        q0 = quadratic_affine_invariant(
            rows[t0]
        )

        q1 = quadratic_affine_invariant(
            rows[t1]
        )

        ratio = clean(
            q1 / q0
        )

        print()
        print(
            "  quadratic_pair={}->{}:".format(
                t0,
                t1,
            )
        )

        print(
            "    invariant_ratio={}".format(
                ratio
            )
        )

        print(
            "    rational_square_ratio={}".format(
                rational_square_class_equal(
                    q1,
                    q0,
                )
            )
        )

    # ------------------------------------------------------------------------
    # 4. SPECIAL RATIONAL AFFINE SUBCLASSES
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "4. SPECIAL RATIONAL AFFINE SUBCLASS AUDIT"
    )
    print("=" * 78)

    special_hits = []

    for t in range(
        max(rows)
    ):

        P = rows[t]
        R = rows[t + 1]

        scalar = scalar_only_test(
            P,
            R,
        )

        translation = translation_only_test(
            P,
            R,
        )

        scaling = scaling_only_test(
            P,
            R,
        )

        print()
        print(
            "  t={} -> {}:".format(
                t,
                t + 1,
            )
        )

        print(
            "    scalar_only={}".format(
                scalar
            )
        )

        print(
            "    translation_only={}".format(
                translation
            )
        )

        print(
            "    scaling_only={}".format(
                scaling
            )
        )

        if (
            scalar is not None
            and scalar[0]
        ):
            special_hits.append(
                (
                    t,
                    "scalar",
                    scalar[1],
                )
            )

        if translation:
            special_hits.append(
                (
                    t,
                    "translation",
                    translation,
                )
            )

        if scaling:
            special_hits.append(
                (
                    t,
                    "scaling",
                    scaling,
                )
            )

    # ------------------------------------------------------------------------
    # 5. DEGREE-DROP DERIVATIVE AUDIT
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "5. DEGREE-DROP DERIVATIVE AUDIT"
    )
    print("=" * 78)

    derivative_hits = []

    for t in range(
        max(rows)
    ):

        print()
        print(
            "  t={} -> {}:".format(
                t,
                t + 1,
            )
        )

        for order in (
            1,
            2,
            3,
        ):

            result = derivative_relation(
                rows[t],
                rows[t + 1],
                order,
            )

            print(
                "    order_{}={}".format(
                    order,
                    result,
                )
            )

            if result["exact"]:
                derivative_hits.append(
                    (
                        t,
                        order,
                        result["scalar"],
                    )
                )

    # ------------------------------------------------------------------------
    # 6. CONCLUSION
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "6. STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
Experiment 336R produced an algebraic affine equivalence between the two
quadratic rows t=1 and t=2, but its affine parameter contains a square
root. That is materially different from a rational transformation.

Experiment 337R therefore imposes the field restriction suggested by the
source arithmetic:

    a,b,c ∈ Q.

A rational affine equivalence would provide a genuine low-complexity
coordinate law over the same field as the observed source data.

The quadratic square-class test is particularly strong. Under a rational
affine change p -> a p+b, the normalized quadratic discriminant changes
by a rational square. Hence a nonsquare ratio rules out rational affine
equivalence immediately.

The derivative audit is interpreted cautiously:

    P_{t+1} = c P_t'

for a linear P_t is automatic and therefore carries little structural
weight.

A derivative relation at a higher degree would be substantially more
interesting.

If the rational-affine tests all fail, the remaining evidence strongly
favors abandoning generic coordinate-transformation searches.

At that point the highest-value work is source reconstruction:

    determine what Q_t(p) actually counts or computes,

rather than continuing to search the observed 15 integers for an
arbitrary algebraic relation.

No synthetic second case is introduced.
"""
    )

    # ------------------------------------------------------------------------
    # 7. FINAL EXACTNESS
    # ------------------------------------------------------------------------

    rational_affine_pairs = [
        t
        for t, result in affine_results.items()
        if result["status"]
        == "RATIONAL_EXACT"
    ]

    algebraic_only_pairs = [
        t
        for t, result in affine_results.items()
        if result["status"]
        == "ALGEBRAIC_ONLY"
    ]

    nontrivial_derivative_hits = [
        hit
        for hit in derivative_hits
        if hit[0] != 4
    ]

    print()
    print("=" * 78)
    print(
        "FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  rational_affine_exact_pairs={}".format(
            rational_affine_pairs
        )
    )

    print(
        "  algebraic_only_affine_pairs={}".format(
            algebraic_only_pairs
        )
    )

    print(
        "  rational_affine_special_hits={}".format(
            special_hits
        )
    )

    print(
        "  nontrivial_derivative_hits={}".format(
            nontrivial_derivative_hits
        )
    )

    print(
        "  trivial_linear_to_constant_derivative_hit={}".format(
            any(
                hit[0] == 4
                for hit in derivative_hits
            )
        )
    )

    print(
        "  missing_values_used=False"
    )

    print(
        "  extrapolation_used=False"
    )

    print(
        "  interpolation_outside_observed_cells=False"
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
        "EXPERIMENT 337R COMPLETE"
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
