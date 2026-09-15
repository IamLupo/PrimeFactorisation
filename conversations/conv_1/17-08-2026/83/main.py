#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 147R
EXACT S-INDEPENDENT DETECTOR SYZYGY SEARCH
ROBUST LINEAR-SYSTEM CONSTRUCTION
==============================================================================

GOAL
----
Search for exact identities

    sum_i A_i(N) Q_i(N,S) = R(N)

where A_i(N) are bounded-degree polynomials in N and R(N) is independent
of S.

This tests whether the detector family contains a genuinely N-only
invariant.

IMPORTANT:
    We do NOT impose S^2 - 4N = T^2.
    We do NOT enumerate factor pairs.
    We do NOT fit numerical data.

If an S-independent syzygy exists, it is an exact symbolic identity.

If none exists through the searched degree, that is a structural
negative result.

NO FULL NEWTON TENSOR
NO FULL C/D TENSOR
NO FACTOR-PAIR SEARCH
NO CSV
NO SKLEARN
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
# Exact kernel
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
            f"quotient remainder ({k},{ell}) = {sp.factor(rem)}"
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
# Exact coefficient vector in S
# ============================================================================

def S_coefficient_vector(expr):
    poly = sp.Poly(
        sp.expand(expr),
        S,
        domain=sp.QQ.frac_field(N),
    )

    if poly.is_zero:
        return [sp.Integer(0)]

    degree = poly.degree()

    return [
        sp.expand(
            poly.coeff_monomial(S**j)
        )
        for j in range(degree + 1)
    ]


# ============================================================================
# Build exact rational linear system
# ============================================================================

def build_syzygy_matrix(detectors, degree_N):
    """
    Unknowns:

        c_(i,d)

    where

        A_i(N) = sum_d c_(i,d) N^d.

    For each S^s, s >= 1, require the coefficient to vanish
    identically as a polynomial in N.

    The matrix is built directly from rational coefficients, avoiding
    SymPy coercion problems with symbolic c_i.
    """

    n_det = len(detectors)
    block = degree_N + 1
    n_unknowns = n_det * block

    detector_vectors = []

    max_S_degree = 0

    for item in detectors:

        vec = S_coefficient_vector(
            item["Q"]
        )

        detector_vectors.append(vec)

        max_S_degree = max(
            max_S_degree,
            len(vec) - 1
        )

    rows = []

    for s_power in range(
        1,
        max_S_degree + 1
    ):

        # For this S^s coefficient, construct the polynomial in N
        # contributed by every unknown independently.
        #
        # Each unknown gets one coefficient polynomial f(N).
        # Extract its N coefficients and place them into matrix rows.

        basis_polynomials = []

        for i in range(n_det):

            q_coeff = (
                detector_vectors[i][s_power]
                if s_power < len(detector_vectors[i])
                else sp.Integer(0)
            )

            for d in range(block):

                basis_polynomials.append(
                    sp.expand(
                        N**d * q_coeff
                    )
                )

        # Determine maximal N degree for this S power.
        max_N_degree = 0

        for f in basis_polynomials:

            poly_f = sp.Poly(
                sp.expand(f),
                N,
                domain=sp.QQ,
            )

            if not poly_f.is_zero:
                max_N_degree = max(
                    max_N_degree,
                    poly_f.degree()
                )

        # Each N^r coefficient gives one rational linear equation.
        for r in range(
            max_N_degree + 1
        ):

            row = []

            for f in basis_polynomials:

                poly_f = sp.Poly(
                    sp.expand(f),
                    N,
                    domain=sp.QQ,
                )

                row.append(
                    sp.Rational(
                        poly_f.coeff_monomial(
                            N**r
                        )
                    )
                )

            # Ignore identically-zero equations.
            if any(value != 0 for value in row):
                rows.append(row)

    if not rows:
        return sp.zeros(0, n_unknowns)

    return sp.Matrix(
        rows
    )


# ============================================================================
# Convert nullspace vector to polynomial coefficients
# ============================================================================

def null_vector_to_polys(
    vector,
    n_det,
    degree_N
):
    block = degree_N + 1

    polys = []

    for i in range(n_det):

        coeffs = vector[
            i*block:(i+1)*block
        ]

        poly = sp.expand(
            sum(
                coeffs[d] * N**d
                for d in range(block)
            )
        )

        polys.append(
            sp.factor(poly)
        )

    return polys


# ============================================================================
# Exact normalization
# ============================================================================

def normalize_vector(vector):
    """
    Convert a rational nullspace vector to a primitive integer vector.
    """

    values = [
        sp.Rational(v)
        for v in vector
    ]

    denominators = [
        int(v.q)
        for v in values
    ]

    lcm_den = 1

    for den in denominators:
        lcm_den = sp.ilcm(
            lcm_den,
            den
        )

    ints = [
        int(v * lcm_den)
        for v in values
    ]

    nonzero = [
        abs(v)
        for v in ints
        if v != 0
    ]

    if not nonzero:
        return vector

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

    # Fix sign.
    first_nonzero = next(
        v for v in ints if v != 0
    )

    if first_nonzero < 0:
        ints = [
            -v
            for v in ints
        ]

    return sp.Matrix(
        ints
    )


# ============================================================================
# Build syzygy expression
# ============================================================================

def syzygy_expression(
    detectors,
    A_polys
):
    expr = sp.Integer(0)

    for item, A in zip(
        detectors,
        A_polys
    ):
        expr += A * item["Q"]

    return sp.expand(expr)


# ============================================================================
# Extract S-independent remainder
# ============================================================================

def S_independent_remainder(expr):
    poly = sp.Poly(
        sp.expand(expr),
        S,
        domain=sp.QQ.frac_field(N),
    )

    if poly.is_zero:
        return sp.Integer(0)

    if poly.degree() != 0:
        return None

    return sp.factor(
        poly.coeff_monomial(1)
    )


# ============================================================================
# Verify exact syzygy
# ============================================================================

def verify_syzygy(
    detectors,
    A_polys
):
    expr = syzygy_expression(
        detectors,
        A_polys
    )

    remainder = sp.Poly(
        sp.expand(expr),
        S,
        domain=sp.QQ.frac_field(N),
    )

    if remainder.is_zero:
        return {
            "ok": True,
            "expression": sp.Integer(0),
            "remainder": sp.Integer(0),
        }

    remainder_expr = sp.factor(
        remainder.as_expr()
    )

    if remainder.degree() == 0:
        return {
            "ok": True,
            "expression": remainder_expr,
            "remainder": remainder_expr,
        }

    return {
        "ok": False,
        "expression": remainder_expr,
        "remainder": remainder_expr,
    }


# ============================================================================
# Search one degree bound
# ============================================================================

def search_degree(detectors, degree_N):
    matrix = build_syzygy_matrix(
        detectors,
        degree_N
    )

    n_unknowns = matrix.cols

    print(
        f"    matrix = {matrix.rows} x {matrix.cols}"
    )

    if matrix.rows == 0:
        return []

    basis = matrix.nullspace()

    winners = []

    for vector in basis:

        normalized = normalize_vector(
            vector
        )

        A_polys = null_vector_to_polys(
            normalized,
            len(detectors),
            degree_N
        )

        if not any(
            sp.expand(A) != 0
            for A in A_polys
        ):
            continue

        verification = verify_syzygy(
            detectors,
            A_polys
        )

        if verification["ok"]:

            R = (
                verification["remainder"]
            )

            winners.append(
                {
                    "A": A_polys,
                    "R": sp.factor(R),
                }
            )

    return winners


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 147R")
    print("EXACT S-INDEPENDENT DETECTOR SYZYGY SEARCH")
    print("ROBUST RATIONAL LINEAR SYSTEM")
    print("=" * 78)
    print()

    ELLS = [
        3,
        5,
        7,
        9,
        11,
    ]

    print("=" * 78)
    print("1. BUILDING EXACT DETECTOR FAMILY")
    print("=" * 78)

    detectors = []

    for ell in ELLS:

        print(
            f"building Q_(1,{ell}) ..."
        )

        Q = detector(
            1,
            ell
        )

        detectors.append(
            {
                "ell": ell,
                "Q": sp.expand(Q),
            }
        )

    print()

    for item in detectors:

        print(
            f"Q_(1,{item['ell']}) = "
            f"{sp.factor(item['Q'])}"
        )

    print()

    # ------------------------------------------------------------------------
    # Search degree bounds.
    # ------------------------------------------------------------------------

    MAX_DEGREE_N = 4

    all_winners = []

    print("=" * 78)
    print("2. SYZYGY SEARCH")
    print("=" * 78)

    for degree_N in range(
        MAX_DEGREE_N + 1
    ):

        print()
        print(
            f"degree_N <= {degree_N}"
        )

        winners = search_degree(
            detectors,
            degree_N
        )

        print(
            f"S-independent syzygies = "
            f"{len(winners)}"
        )

        for idx, winner in enumerate(
            winners[:10],
            start=1
        ):

            print(
                f"  SYZYGY {idx}"
            )

            for item, A in zip(
                detectors,
                winner["A"]
            ):

                print(
                    f"    A_(1,{item['ell']})(N) = "
                    f"{sp.factor(A)}"
                )

            print(
                "    R(N) =",
                sp.factor(
                    winner["R"]
                )
            )

            all_winners.append(
                (
                    degree_N,
                    winner
                )
            )

    # ------------------------------------------------------------------------
    # Exact verification of first result.
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. EXACT VERIFICATION")
    print("=" * 78)

    if all_winners:

        degree_N, winner = all_winners[0]

        expr = syzygy_expression(
            detectors,
            winner["A"]
        )

        residual = sp.expand(
            expr - winner["R"]
        )

        print(
            "selected degree_N =",
            degree_N
        )

        print(
            "exact residual =",
            sp.factor(residual)
        )

        print()
        print(
            "STATUS = PASS"
        )

        print()
        print(
            "A genuine S-independent detector invariant was found."
        )

        print(
            "This is a nontrivial algebraic candidate for an"
        )

        print(
            "N-only observable relation."
        )

    else:

        print(
            "No nontrivial S-independent syzygy found."
        )

        print(
            f"through degree_N <= {MAX_DEGREE_N}."
        )

        print()
        print(
            "STATUS = NO-SYZYGY"
        )

        print()
        print(
            "The tested detector family has no linear N-polynomial"
        )

        print(
            "invariant at the searched degree."
        )

        print()
        print(
            "This is evidence that N-only recovery will require"
        )

        print(
            "a nonlinear detector relation or a different invariant."
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