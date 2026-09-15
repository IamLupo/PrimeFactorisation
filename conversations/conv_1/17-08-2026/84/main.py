#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 148
EXACT QUADRATIC DETECTOR-INVARIANT SEARCH

PURPOSE
-------
Experiment 147R ruled out low-degree LINEAR detector syzygies:

    sum_i A_i(N) Q_i(N,S) = R(N)

through degree_N <= 4.

The next question is whether the first S-independent invariant is
NONLINEAR in the detector coordinates.

We therefore search for an exact identity of total detector degree <= 2:

    P(N,Q13,Q15,Q17) = 0

where

    P = sum_m A_m(N) M_m(Q13,Q15,Q17),

and M_m runs over

    1,
    Qi,
    Qi^2,
    Qi*Qj.

Each A_m(N) is a polynomial in N of bounded degree.

The condition is:

    all coefficients of S^r, r >= 1, vanish identically.

The surviving S^0 part is the N-only polynomial relation.

This is exact linear algebra over QQ.

IMPORTANT
---------
This experiment asks whether a low-degree algebraic relation exists
among detector values and N.

It does NOT claim that such a relation gives a factorisation algorithm.

A successful relation is only the next algebraic bridge.

NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
NO FULL NEWTON TENSOR
NO FULL C/D TENSOR
SAFE SYMPY DOMAINS
==============================================================================
"""

import sys
import sympy as sp


# ============================================================================
# Symbols
# ============================================================================

p, q = sp.symbols("p q")
N, S = sp.symbols("N S")


# ============================================================================
# Exact kernel / quotient
# ============================================================================

def kernel_F(k, ell):
    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


def quotient_Q_pq(k, ell):
    F = kernel_F(k, ell)

    PF = sp.Poly(
        F,
        p,
        domain=sp.QQ.frac_field(q),
    )

    PD = sp.Poly(
        p + q + 1,
        p,
        domain=sp.QQ.frac_field(q),
    )

    Q, R = sp.div(PF, PD)

    rem = sp.expand(R.as_expr())

    if rem != 0:
        raise ArithmeticError(
            f"quotient remainder ({k},{ell}) = "
            f"{sp.factor(rem)}"
        )

    return sp.expand(Q.as_expr())


# ============================================================================
# Symmetric reduction
# ============================================================================

def symmetric_to_NS(expr):
    expr = sp.expand(
        expr.subs(q, S - p)
    )

    modulus = sp.Poly(
        p**2 - S*p + N,
        p,
        domain=sp.QQ.frac_field(S, N),
    )

    poly = sp.Poly(
        expr,
        p,
        domain=sp.QQ.frac_field(S, N),
    )

    rem = sp.expand(
        sp.rem(poly, modulus).as_expr()
    )

    if sp.expand(rem.coeff(p, 1)) != 0:
        raise ArithmeticError(
            "symmetric reduction retained p-dependence"
        )

    return sp.expand(
        rem.coeff(p, 0)
    )


def detector(k, ell):
    return symmetric_to_NS(
        quotient_Q_pq(k, ell)
    )


# ============================================================================
# Exact detector family
# ============================================================================

def build_detectors():
    """
    Keep the first three detectors initially.

    Their S-degrees are:

        Q13 : 2
        Q15 : 4
        Q17 : 6

    This is enough to test whether a genuinely quadratic detector
    invariant appears before the algebra becomes unnecessarily large.
    """

    specs = [
        (1, 3),
        (1, 5),
        (1, 7),
    ]

    result = []

    for k, ell in specs:

        print(
            f"building Q_({k},{ell}) ..."
        )

        Q = detector(k, ell)

        result.append(
            {
                "k": k,
                "ell": ell,
                "Q": sp.expand(Q),
            }
        )

    return result


# ============================================================================
# Detector monomials of total degree <= 2
# ============================================================================

def detector_monomials(Qs):
    """
    Return:

        1
        Q1,Q2,Q3
        Q1^2,Q1Q2,Q1Q3,Q2^2,Q2Q3,Q3^2
    """

    q1, q2, q3 = Qs

    return [
        ("1", sp.Integer(1)),

        ("Q1", q1),
        ("Q2", q2),
        ("Q3", q3),

        ("Q1^2", sp.expand(q1*q1)),
        ("Q1Q2", sp.expand(q1*q2)),
        ("Q1Q3", sp.expand(q1*q3)),
        ("Q2^2", sp.expand(q2*q2)),
        ("Q2Q3", sp.expand(q2*q3)),
        ("Q3^2", sp.expand(q3*q3)),
    ]


# ============================================================================
# Build coefficient matrix
# ============================================================================

def build_matrix(
    monomials,
    degree_N
):
    """
    Unknown coefficient attached to every pair

        (detector monomial, N^d).

    For every S^r, r>=1, require exact cancellation.

    The system is built from rational coefficients only.
    """

    n_monomials = len(monomials)

    block = degree_N + 1
    n_unknowns = n_monomials * block

    rows = []

    # Maximum S degree across all detector monomials.
    max_S_degree = 0

    coeff_vectors = []

    for _, expr in monomials:

        poly = sp.Poly(
            sp.expand(expr),
            S,
            domain=sp.QQ.frac_field(N),
        )

        if poly.is_zero:
            vec = [sp.Integer(0)]
        else:
            deg = poly.degree()

            vec = [
                sp.expand(
                    poly.coeff_monomial(
                        S**r
                    )
                )
                for r in range(deg + 1)
            ]

            max_S_degree = max(
                max_S_degree,
                deg
            )

        coeff_vectors.append(vec)

    # For each nonconstant S power, construct equations in N.
    for s_power in range(
        1,
        max_S_degree + 1
    ):

        basis = []

        for _, _expr in monomials:

            q_coeff = (
                coeff_vectors[len(basis)][s_power]
                if False
                else None
            )

        for i, (_, _expr) in enumerate(
            monomials
        ):

            vec = coeff_vectors[i]

            q_coeff = (
                vec[s_power]
                if s_power < len(vec)
                else sp.Integer(0)
            )

            for d in range(
                block
            ):

                basis.append(
                    sp.expand(
                        N**d * q_coeff
                    )
                )

        max_N_degree = -1

        for expr in basis:

            poly_N = sp.Poly(
                sp.expand(expr),
                N,
                domain=sp.QQ,
            )

            if not poly_N.is_zero:
                max_N_degree = max(
                    max_N_degree,
                    poly_N.degree()
                )

        if max_N_degree < 0:
            continue

        for n_power in range(
            max_N_degree + 1
        ):

            row = []

            for expr in basis:

                poly_N = sp.Poly(
                    sp.expand(expr),
                    N,
                    domain=sp.QQ,
                )

                row.append(
                    sp.Rational(
                        poly_N.coeff_monomial(
                            N**n_power
                        )
                    )
                )

            if any(
                value != 0
                for value in row
            ):
                rows.append(row)

    if not rows:
        return sp.zeros(
            0,
            n_unknowns
        )

    return sp.Matrix(rows)


# ============================================================================
# Convert vector -> N-polynomial coefficient functions
# ============================================================================

def vector_to_polys(
    vector,
    n_monomials,
    degree_N
):
    block = degree_N + 1

    result = []

    for i in range(
        n_monomials
    ):

        coeff_block = vector[
            i*block:(i+1)*block
        ]

        expr = sp.expand(
            sum(
                coeff_block[d] * N**d
                for d in range(block)
            )
        )

        result.append(
            sp.factor(expr)
        )

    return result


# ============================================================================
# Primitive normalization
# ============================================================================

def normalize_vector(vector):

    rational_values = [
        sp.Rational(v)
        for v in vector
    ]

    denominators = [
        int(v.q)
        for v in rational_values
    ]

    lcm_den = 1

    for den in denominators:
        lcm_den = sp.ilcm(
            lcm_den,
            den
        )

    ints = [
        int(v*lcm_den)
        for v in rational_values
    ]

    nonzero = [
        abs(v)
        for v in ints
        if v != 0
    ]

    if not nonzero:
        return sp.Matrix(ints)

    g = nonzero[0]

    for value in nonzero[1:]:
        g = sp.igcd(
            g,
            value
        )

    if g > 1:
        ints = [
            v // g
            for v in ints
        ]

    first = next(
        v for v in ints
        if v != 0
    )

    if first < 0:
        ints = [
            -v
            for v in ints
        ]

    return sp.Matrix(
        ints
    )


# ============================================================================
# Construct relation
# ============================================================================

def construct_relation(
    monomials,
    coeff_polys
):
    expr = sp.Integer(0)

    for (_, monomial), A in zip(
        monomials,
        coeff_polys
    ):
        expr += A * monomial

    return sp.expand(expr)


# ============================================================================
# Verify relation
# ============================================================================

def verify_relation(expr):

    poly = sp.Poly(
        sp.expand(expr),
        S,
        domain=sp.QQ.frac_field(N),
    )

    if poly.is_zero:
        return {
            "ok": True,
            "degree_S": -sp.oo,
            "remainder": sp.Integer(0),
        }

    degree_S = poly.degree()

    if degree_S == 0:
        return {
            "ok": True,
            "degree_S": 0,
            "remainder": sp.factor(
                poly.coeff_monomial(1)
            ),
        }

    return {
        "ok": False,
        "degree_S": degree_S,
        "remainder": sp.factor(
            poly.as_expr()
        ),
    }


# ============================================================================
# Remove relations that are trivial because they contain only the constant
# monomial 1.
# ============================================================================

def is_nontrivial_detector_relation(
    coeff_polys
):
    return any(
        sp.expand(A) != 0
        for A in coeff_polys[1:]
    )


# ============================================================================
# Search one degree bound
# ============================================================================

def search_degree(
    monomials,
    degree_N
):

    matrix = build_matrix(
        monomials,
        degree_N
    )

    print(
        f"    matrix = "
        f"{matrix.rows} x {matrix.cols}"
    )

    if matrix.cols == 0:
        return []

    nullspace = matrix.nullspace()

    print(
        f"    nullspace dimension = "
        f"{len(nullspace)}"
    )

    winners = []

    for vector in nullspace:

        normalized = normalize_vector(
            vector
        )

        coeff_polys = vector_to_polys(
            normalized,
            len(monomials),
            degree_N
        )

        if not is_nontrivial_detector_relation(
            coeff_polys
        ):
            continue

        relation = construct_relation(
            monomials,
            coeff_polys
        )

        verification = verify_relation(
            relation
        )

        if verification["ok"]:

            winners.append(
                {
                    "coeffs": coeff_polys,
                    "relation": relation,
                    "R": verification["remainder"],
                }
            )

    return winners


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 148")
    print("EXACT QUADRATIC DETECTOR-INVARIANT SEARCH")
    print("=" * 78)
    print()

    # ------------------------------------------------------------------------
    # Build detector family.
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("1. BUILDING DETECTOR FAMILY")
    print("=" * 78)

    detectors = build_detectors()

    Qs = [
        item["Q"]
        for item in detectors
    ]

    print()

    for item in detectors:
        print(
            f"Q_({item['k']},{item['ell']}) = "
            f"{sp.factor(item['Q'])}"
        )

    # ------------------------------------------------------------------------
    # Monomials.
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. DETECTOR MONOMIAL BASIS")
    print("=" * 78)

    monomials = detector_monomials(
        Qs
    )

    for name, _ in monomials:
        print(
            " ",
            name
        )

    print(
        "monomial count =",
        len(monomials)
    )

    # ------------------------------------------------------------------------
    # Search.
    # ------------------------------------------------------------------------

    MAX_DEGREE_N = 4

    print()
    print("=" * 78)
    print("3. NONLINEAR SYZYGY SEARCH")
    print("=" * 78)

    winners_all = []

    for degree_N in range(
        MAX_DEGREE_N + 1
    ):

        print()
        print(
            f"degree_N <= {degree_N}"
        )

        winners = search_degree(
            monomials,
            degree_N
        )

        print(
            f"    S-independent nonlinear relations = "
            f"{len(winners)}"
        )

        for idx, winner in enumerate(
            winners[:10],
            start=1
        ):

            print(
                f"  RELATION {idx}"
            )

            for (
                (name, _),
                A
            ) in zip(
                monomials,
                winner["coeffs"]
            ):

                if sp.expand(A) != 0:
                    print(
                        f"    {name}: "
                        f"{sp.factor(A)}"
                    )

            print(
                "    N-only remainder =",
                sp.factor(
                    winner["R"]
                )
            )

            print(
                "    relation ="
            )

            print(
                sp.factor(
                    winner["relation"]
                )
            )

            winners_all.append(
                (
                    degree_N,
                    winner
                )
            )

    # ------------------------------------------------------------------------
    # Final.
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    if winners_all:

        degree_N, winner = winners_all[0]

        print(
            "STATUS = PASS"
        )

        print()
        print(
            "A genuine quadratic S-independent detector relation"
        )

        print(
            "was found."
        )

        print()
        print(
            "lowest searched N-degree =",
            degree_N
        )

        print()
        print(
            "The next task is to determine whether this relation"
        )

        print(
            "can be inverted using N alone."
        )

    else:

        print(
            "STATUS = NO-QUADRATIC-RELATION"
        )

        print()
        print(
            "No nontrivial S-independent detector relation"
        )

        print(
            "of detector degree <= 2 was found through"
        )

        print(
            f"N-degree <= {MAX_DEGREE_N}."
        )

        print()
        print(
            "This rules out the simplest nonlinear invariant class."
        )

        print()
        print(
            "NEXT TARGET:"
        )

        print(
            "either search cubic detector relations or"
        )

        print(
            "switch to a recurrence/generating-function route"
        )

        print(
            "for detector evaluation from N."
        )

    print("=" * 78)


if __name__ == "__main__":

    try:
        main()

    except Exception as exc:

        print()
        print(
            "FATAL:",
            type(exc).__name__,
            str(exc)
        )

        sys.exit(1)

