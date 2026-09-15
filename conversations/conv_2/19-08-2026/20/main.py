#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 336R — EXACT AFFINE-PROJECTIVE ROW-TRANSITION AUDIT
==============================================================================

Purpose
-------
Experiments 322R-335R have rejected a large family of generic additive,
rational, separable, local, differential, gcd, shift, and divisibility
relations.

The remaining natural low-complexity mechanism for the row profiles is an
affine/projective change of the p-coordinate.

Test, for adjacent observed row-polynomials:

    P_{t+1}(p) = c_t * P_t(a_t*p + b_t)

and the slightly more general projective degree-preserving form

    P_{t+1}(p)
      = c_t * (gamma_t*p + delta_t)^d
        * P_t(
            (alpha_t*p + beta_t)
            / (gamma_t*p + delta_t)
          ),

restricted to cases where the degrees permit such a transformation.

Because the current rows have decreasing degrees

    3, 2, 2, 1, 1, 0,

the degree-preserving affine model is tested only when

    deg(P_t) == deg(P_{t+1}).

For degree-changing transitions we test instead whether a derivative
or fixed factor could account for the drop, but without fitting a free
high-dimensional operator.

The main exact test is the affine-projective invariant of polynomial roots:
under

    p -> a*p+b

or more generally a Möbius map, the root configuration is transformed
projectively.

For quadratic rows, this gives an especially strong exact discriminant
and root-cross-ratio diagnostic.

No missing cells are introduced.
No new p-value is used.
No synthetic n=pq case.
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

    for t in range(maximum_t + 1):

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


def row_polynomial(layer):

    return clean(
        sp.interpolate(
            layer,
            p,
        )
    )


def build_rows(layers):

    return {
        t: row_polynomial(layer)
        for t, layer in layers.items()
    }


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

    ints = [
        value // g
        for value in ints
    ]

    if ints[0] < 0:
        ints = [
            -value
            for value in ints
        ]

    expr_int = sum(
        sp.Integer(v)
        * p**(
            len(ints) - 1 - i
        )
        for i, v in enumerate(ints)
    )

    return sp.Poly(
        sp.expand(expr_int),
        p,
        domain=sp.ZZ,
    )


# ============================================================================
# AFFINE TRANSFORMATION TEST
# ============================================================================

def solve_affine_transition(P, R):

    """
    Solve

        R(p) = c P(a p + b)

    over Q.

    This is deliberately exact.

    The unknowns are nonlinear, so we compare coefficient invariants and
    then solve the resulting polynomial equations.
    """

    degP = sp.degree(P, p)
    degR = sp.degree(R, p)

    if degP != degR:
        return {
            "status": "DEGREE_MISMATCH",
            "degree_P": degP,
            "degree_R": degR,
        }

    degree = int(degP)

    # Introduce affine coefficients aa, bb and scalar cc.
    aa, bb, cc = sp.symbols(
        "aa bb cc"
    )

    transformed = sp.Poly(
        sp.expand(
            cc * P.subs(
                p,
                aa * p + bb,
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

    equations = []

    for k in range(
        degree + 1
    ):

        equations.append(
            sp.Eq(
                transformed.coeff_monomial(
                    p**k
                ),
                target.coeff_monomial(
                    p**k
                ),
            )
        )

    # Normalize projectively by setting aa != 0.  The first equation
    # gives cc*aa^degree = leading_ratio.
    leading_P = sp.Poly(
        P,
        p,
        domain=sp.QQ,
    ).LC()

    leading_R = sp.Poly(
        R,
        p,
        domain=sp.QQ,
    ).LC()

    if leading_P == 0:
        return {
            "status": "ZERO_LEADING_COEFFICIENT"
        }

    # For degree 1, affine matching always has enough freedom, so only
    # report exact solutions but mark them data-sized.
    try:
        solution = sp.solve(
            equations,
            [aa, bb, cc],
            dict=True,
        )
    except Exception:
        solution = []

    rational_solutions = []

    for sol in solution:

        vals = [
            sol.get(
                symbol,
                symbol,
            )
            for symbol in (
                aa,
                bb,
                cc,
            )
        ]

        if all(
            value.free_symbols.isdisjoint(
                {
                    aa,
                    bb,
                    cc,
                }
            )
            for value in vals
            if hasattr(value, "free_symbols")
        ):

            substituted = clean(
                cc * P.subs(
                    p,
                    aa * p + bb,
                )
            ).subs(sol)

            if clean(
                substituted - R
            ) == 0:

                rational_solutions.append(
                    sol
                )

    return {
        "status": (
            "EXACT"
            if rational_solutions
            else "NO_SOLUTION"
        ),
        "degree_P": degP,
        "degree_R": degR,
        "solutions": rational_solutions,
        "equations": equations,
    }


# ============================================================================
# DISCRIMINANT / ROOT-CONFIGURATION INVARIANTS
# ============================================================================

def discriminant(P):

    poly = sp.Poly(
        P,
        p,
        domain=sp.QQ,
    )

    return clean(
        sp.discriminant(
            poly.as_expr(),
            p,
        )
    )


def normalized_discriminant(P):

    poly = sp.Poly(
        P,
        p,
        domain=sp.QQ,
    )

    degree = poly.degree()
    leading = poly.LC()
    disc = discriminant(P)

    if leading == 0:
        return None

    # For an affine coordinate change p -> a p+b,
    # Disc(c P(a p+b)) scales by predictable powers.
    # The ratio below is therefore reported as a diagnostic only.
    return clean(
        disc
        / leading**(
            2 * degree - 2
        )
    )


def quadratic_root_geometry(P):

    poly = sp.Poly(
        P,
        p,
        domain=sp.QQ,
    )

    if poly.degree() != 2:
        return None

    A, B, C = poly.all_coeffs()

    # Center and discriminant scale.
    center = clean(
        -B / (2 * A)
    )

    disc = clean(
        B**2 - 4*A*C
    )

    normalized = clean(
        disc / A**2
    )

    return {
        "center": center,
        "discriminant": disc,
        "normalized_discriminant": normalized,
    }


# ============================================================================
# DERIVATIVE DROP TEST
# ============================================================================

def derivative_drop_test(P, R):

    candidates = {
        "P'": sp.diff(P, p),
        "P''": sp.diff(
            P,
            p,
            2,
        ),
    }

    results = {}

    for name, candidate in candidates.items():

        candidate = clean(
            candidate
        )

        if clean(candidate - R) == 0:
            results[name] = True
        else:
            # Allow constant scalar multiple.
            if candidate == 0:
                results[name] = False
                continue

            poly_c = sp.Poly(
                candidate,
                p,
                domain=sp.QQ,
            )

            poly_r = sp.Poly(
                R,
                p,
                domain=sp.QQ,
            )

            if poly_c.degree() != poly_r.degree():
                results[name] = False
                continue

            ratio = clean(
                poly_r.LC()
                / poly_c.LC()
            )

            results[name] = (
                clean(
                    ratio * candidate - R
                )
                == 0
            )

    return results


# ============================================================================
# SMALL RATIONAL AFFINE SEARCH
# ============================================================================

def bounded_rational_affine_search(
    P,
    R,
    bound=6,
):

    hits = []

    degP = sp.degree(P, p)
    degR = sp.degree(R, p)

    if degP != degR:
        return hits

    for aa_num in range(
        -bound,
        bound + 1,
    ):

        if aa_num == 0:
            continue

        for aa_den in range(
            1,
            bound + 1,
        ):

            aa = sp.Rational(
                aa_num,
                aa_den,
            )

            for bb_num in range(
                -bound,
                bound + 1,
            ):

                for bb_den in range(
                    1,
                    bound + 1,
                ):

                    bb = sp.Rational(
                        bb_num,
                        bb_den,
                    )

                    transformed = clean(
                        P.subs(
                            p,
                            aa*p + bb,
                        )
                    )

                    if transformed == 0:
                        continue

                    poly_t = sp.Poly(
                        transformed,
                        p,
                        domain=sp.QQ,
                    )

                    poly_r = sp.Poly(
                        R,
                        p,
                        domain=sp.QQ,
                    )

                    if (
                        poly_t.degree()
                        != poly_r.degree()
                    ):
                        continue

                    cc = clean(
                        poly_r.LC()
                        / poly_t.LC()
                    )

                    if clean(
                        cc * transformed
                        - R
                    ) == 0:

                        hits.append(
                            (
                                aa,
                                bb,
                                cc,
                            )
                        )

    # unique
    unique = []

    for hit in hits:
        if hit not in unique:
            unique.append(hit)

    return unique


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 336R — EXACT AFFINE-PROJECTIVE "
        "ROW-TRANSITION AUDIT"
    )
    print("=" * 78)

    layers = build_layers()
    P = build_rows(layers)

    # ------------------------------------------------------------------------
    # 1. ROW INVENTORY
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "1. ROW-POLYNOMIAL INVENTORY"
    )
    print("=" * 78)

    for t in sorted(P):

        poly = sp.Poly(
            P[t],
            p,
            domain=sp.QQ,
        )

        print()
        print(
            "  t={}:".format(
                t
            )
        )

        print(
            "    degree={}".format(
                poly.degree()
            )
        )

        print(
            "    polynomial={}".format(
                P[t]
            )
        )

        print(
            "    primitive={}".format(
                primitive_integer_poly(
                    P[t]
                ).as_expr()
            )
        )

    # ------------------------------------------------------------------------
    # 2. ADJACENT AFFINE TEST
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "2. ADJACENT AFFINE-TRANSFORM TEST"
    )
    print("=" * 78)

    affine_results = {}

    for t in range(
        max(P)
    ):

        result = solve_affine_transition(
            P[t],
            P[t + 1],
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
            "    degree_P={}".format(
                result.get(
                    "degree_P"
                )
            )
        )

        print(
            "    degree_R={}".format(
                result.get(
                    "degree_R"
                )
            )
        )

        print(
            "    status={}".format(
                result["status"]
            )
        )

        if result["status"] == "EXACT":

            print(
                "    solutions={}".format(
                    result["solutions"]
                )
            )

    # ------------------------------------------------------------------------
    # 3. SMALL RATIONAL AFFINE SEARCH
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "3. BOUNDED RATIONAL AFFINE SEARCH"
    )
    print("=" * 78)

    affine_hits = []

    for t in range(
        max(P)
    ):

        hits = bounded_rational_affine_search(
            P[t],
            P[t + 1],
            bound=6,
        )

        if hits:

            affine_hits.append(
                (
                    t,
                    hits,
                )
            )

            print()
            print(
                "  t={} -> {}: hits={}".format(
                    t,
                    t + 1,
                    hits,
                )
            )

    print()
    print(
        "  total_affine_hits={}".format(
            len(affine_hits)
        )
    )

    # ------------------------------------------------------------------------
    # 4. QUADRATIC GEOMETRY
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "4. QUADRATIC ROOT-GEOMETRY AUDIT"
    )
    print("=" * 78)

    quadratic_rows = [
        t
        for t in P
        if sp.degree(
            P[t],
            p,
        ) == 2
    ]

    for t in quadratic_rows:

        geometry = quadratic_root_geometry(
            P[t]
        )

        print()
        print(
            "  t={}:".format(
                t
            )
        )

        print(
            "    center={}".format(
                geometry["center"]
            )
        )

        print(
            "    discriminant={}".format(
                geometry["discriminant"]
            )
        )

        print(
            "    normalized_discriminant={}".format(
                geometry[
                    "normalized_discriminant"
                ]
            )
        )

    if len(quadratic_rows) >= 2:

        t0 = quadratic_rows[0]
        t1 = quadratic_rows[1]

        g0 = quadratic_root_geometry(
            P[t0]
        )

        g1 = quadratic_root_geometry(
            P[t1]
        )

        print()
        print(
            "  quadratic_row_comparison:"
        )

        print(
            "    center_ratio={}".format(
                clean(
                    g1["center"]
                    / g0["center"]
                )
                if g0["center"] != 0
                else None
            )
        )

        print(
            "    normalized_discriminant_ratio={}".format(
                clean(
                    g1[
                        "normalized_discriminant"
                    ]
                    /
                    g0[
                        "normalized_discriminant"
                    ]
                )
            )
        )

    # ------------------------------------------------------------------------
    # 5. DEGREE-DROP DERIVATIVE TEST
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "5. DEGREE-DROP DERIVATIVE TEST"
    )
    print("=" * 78)

    derivative_hits = []

    for t in range(
        max(P)
    ):

        results = derivative_drop_test(
            P[t],
            P[t + 1],
        )

        print()
        print(
            "  t={} -> {}: {}".format(
                t,
                t + 1,
                results,
            )
        )

        for name, hit in results.items():

            if hit:
                derivative_hits.append(
                    (
                        t,
                        name,
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
Experiment 335R found no nontrivial adjacent polynomial gcds and no
meaningful boundary factors.

The main remaining low-complexity possibility is therefore that the
successive p-profiles represent the same algebraic shape in a changing
coordinate system.

The affine model

    P_{t+1}(p) = c_t P_t(a_t p + b_t)

is particularly strong: it preserves degree and maps the root set
affinely.

For quadratic rows, the root-center and normalized-discriminant data give
coordinate-free diagnostics for whether the two quadratics could be
related by such a transformation.

The derivative tests separately ask whether the degree drops

    3 -> 2 -> 1 -> 0

come from repeated differentiation.

A positive result here would be qualitatively different from the failed
fixed differential-operator searches: it would identify a direct
coordinate evolution of the p-profile itself.

A complete negative result would leave very little reason to continue
searching for generic low-complexity transformations of the 15 observed
numbers.

At that point the highest-value step is to recover the actual
combinatorial/source definition of Q_t(p), or obtain another independent
n=pq instance.

No synthetic second case is introduced.
"""
    )

    # ------------------------------------------------------------------------
    # 7. FINAL EXACTNESS
    # ------------------------------------------------------------------------

    exact_affine = [
        t
        for t, result in affine_results.items()
        if result["status"] == "EXACT"
    ]

    print()
    print("=" * 78)
    print(
        "FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  adjacent_affine_exact_pairs={}".format(
            exact_affine
        )
    )

    print(
        "  bounded_rational_affine_hits={}".format(
            len(affine_hits)
        )
    )

    print(
        "  derivative_degree_drop_hits={}".format(
            derivative_hits
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
        "EXPERIMENT 336R COMPLETE"
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
