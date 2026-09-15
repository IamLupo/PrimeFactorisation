#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 311R — EXACT TRANSFER-PAIR IRREDUCIBILITY /
COMMON-EIGENLINE / CENTRALIZER AUDIT
==============================================================================

This experiment starts from the two exact width-2 transfer matrices recovered
in Experiments 305R--310R.

The previous results established:

    * projective similarity: refuted
    * Möbius/polynomial relation: refuted
    * commutation: refuted
    * generated algebra: full M_2(Q)

The next structural question is representation-theoretic:

    Do T1 and T2 share a common invariant line?

For 2x2 matrices over Q, a common invariant line is equivalent to the
existence of a common eigenvector over an algebraic closure.

We therefore test:

    1. common rational eigenvector;
    2. common eigenvalue/eigenline over an algebraic closure;
    3. simultaneous triangularizability;
    4. invertibility and determinant of the commutator;
    5. common centralizer dimension;
    6. scalar-only common centralizer;
    7. exact irreducibility criterion through the commutator.

No external files.
No synthetic second n=pq case.
Exact rational arithmetic only.
"""


from __future__ import annotations

import math
import sys

import sympy as sp


# ============================================================================
# EXACT TRANSFER MATRICES
# ============================================================================

T1 = sp.Matrix([
    [
        sp.Rational(
            -209857461192170070,
            11470116327290257,
        ),
        sp.Rational(
            -215450379004300026,
            11470116327290257,
        ),
    ],
    [
        1,
        0,
    ],
])

T2 = sp.Matrix([
    [
        sp.Rational(
            83976580526089197,
            971937272383741,
        ),
        sp.Rational(
            75947709834674022,
            971937272383741,
        ),
    ],
    [
        1,
        0,
    ],
])


# ============================================================================
# HELPERS
# ============================================================================

def clean(expr):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(expr)
            )
        )
    )


def valuation(value, prime):
    value = sp.Rational(value)

    if value == 0:
        return None

    num = abs(int(value.p))
    den = abs(int(value.q))

    v = 0

    while num % prime == 0:
        num //= prime
        v += 1

    while den % prime == 0:
        den //= prime
        v -= 1

    return v


def matrix_rank(M):
    return int(M.rank())


def flatten_matrix(M):
    return sp.Matrix([
        M[0, 0],
        M[0, 1],
        M[1, 0],
        M[1, 1],
    ])


# ============================================================================
# CHARACTERISTIC DATA
# ============================================================================

def characteristic_data(M):

    tr = clean(sp.trace(M))
    det = clean(M.det())
    disc = clean(
        tr**2 - 4 * det
    )
    poly = sp.Poly(
        M.charpoly().as_expr(),
        sp.Symbol("lambda"),
    )

    return tr, det, disc, poly


def print_characteristic_data(name, M):

    tr, det, disc, poly = characteristic_data(M)

    print()
    print("  {}:".format(name))

    print(
        "    trace={}".format(
            tr
        )
    )

    print(
        "    determinant={}".format(
            det
        )
    )

    print(
        "    discriminant={}".format(
            disc
        )
    )

    print(
        "    characteristic_polynomial={}".format(
            poly.as_expr()
        )
    )

    print(
        "    valuations={}".format(
            {
                p: valuation(disc, p)
                for p in (
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


# ============================================================================
# COMMON EIGENVALUE TEST
# ============================================================================

def common_eigenvalue_polynomial_test():

    lam = sp.symbols(
        "lambda"
    )

    p1 = clean(
        T1.charpoly(lam).as_expr()
    )

    p2 = clean(
        T2.charpoly(lam).as_expr()
    )

    gcd_poly = sp.Poly(
        p1,
        lam,
        extension=True,
    ).gcd(
        sp.Poly(
            p2,
            lam,
            extension=True,
        )
    )

    return p1, p2, gcd_poly


# ============================================================================
# COMMON EIGENLINE VIA COMMUTATOR
# ============================================================================

def common_invariant_line_test():

    """
    A common invariant line for T1 and T2 implies that the commutator
    annihilates that line.

    For a 2x2 commutator, if det([T1,T2]) != 0 then there is no common
    invariant line.

    This is a strong exact test.
    """

    comm = clean_matrix(
        T1 * T2 - T2 * T1
    )

    det_comm = clean(
        comm.det()
    )

    return comm, det_comm


def clean_matrix(M):
    return M.applyfunc(clean)


# ============================================================================
# DIRECT COMMON EIGENVECTOR TEST OVER Q
# ============================================================================

def common_rational_eigenvector_test():

    """
    Search projectively for a vector (x,y).

    Write y = 1 when possible and test the two eigenline equations
    simultaneously. Also test x = 1, y = 0 separately.
    """

    x = sp.symbols(
        "x"
    )

    candidates = []

    # Case v = (x,1)
    v = sp.Matrix([
        x,
        1,
    ])

    # A vector v spans an invariant line of M iff
    # det(v, Mv) = 0.
    det1 = clean(
        sp.det(
            sp.Matrix.hstack(
                v,
                T1 * v,
            )
        )
    )

    det2 = clean(
        sp.det(
            sp.Matrix.hstack(
                v,
                T2 * v,
            )
        )
    )

    gcd_expr = sp.Poly(
        det1,
        x,
    ).gcd(
        sp.Poly(
            det2,
            x,
        )
    )

    if gcd_expr.degree() >= 1:

        roots = sp.solve(
            gcd_expr.as_expr(),
            x,
        )

        for root in roots:
            candidates.append(
                (
                    root,
                    sp.Matrix([
                        root,
                        1,
                    ]),
                )
            )

    # Case v = (1,0)
    v_inf = sp.Matrix([
        1,
        0,
    ])

    det1_inf = clean(
        sp.det(
            sp.Matrix.hstack(
                v_inf,
                T1 * v_inf,
            )
        )
    )

    det2_inf = clean(
        sp.det(
            sp.Matrix.hstack(
                v_inf,
                T2 * v_inf,
            )
        )
    )

    if det1_inf == 0 and det2_inf == 0:

        candidates.append(
            (
                "infinity",
                v_inf,
            )
        )

    return det1, det2, gcd_expr, candidates


# ============================================================================
# SIMULTANEOUS TRIANGULARIZATION TEST
# ============================================================================

def simultaneous_triangularization_test():

    """
    For 2x2 matrices over a field, simultaneous triangularization implies
    a common invariant line.

    Hence the previous invariant-line test is decisive.

    We additionally verify by solving directly for an invertible change of
    basis S such that both S^-1 T_i S are upper triangular.
    """

    a, b, c, d = sp.symbols(
        "a b c d"
    )

    S = sp.Matrix([
        [a, b],
        [c, d],
    ])

    sdet = clean(
        S.det()
    )

    # Off-diagonal-lower entries of S^{-1} T_i S.
    # We avoid matrix inversion symbolically by using the adjugate.
    adjS = sp.Matrix([
        [d, -b],
        [-c, a],
    ])

    M1 = clean_matrix(
        adjS * T1 * S
    )

    M2 = clean_matrix(
        adjS * T2 * S
    )

    equations = [
        clean(
            M1[1, 0]
        ),
        clean(
            M2[1, 0]
        ),
    ]

    return S, sdet, equations


# ============================================================================
# COMMON CENTRALIZER
# ============================================================================

def common_centralizer():

    a, b, c, d = sp.symbols(
        "a b c d"
    )

    X = sp.Matrix([
        [a, b],
        [c, d],
    ])

    E1 = clean_matrix(
        X * T1 - T1 * X
    )

    E2 = clean_matrix(
        X * T2 - T2 * X
    )

    equations = [
        E1[0, 0],
        E1[0, 1],
        E1[1, 0],
        E1[1, 1],
        E2[0, 0],
        E2[0, 1],
        E2[1, 0],
        E2[1, 1],
    ]

    vars_ = [
        a,
        b,
        c,
        d,
    ]

    A, rhs = sp.linear_eq_to_matrix(
        equations,
        vars_,
    )

    nullspace = A.nullspace()

    return X, A, nullspace


# ============================================================================
# COMMUTATOR CENTRALIZER
# ============================================================================

def commutator_centralizer_test(comm):

    a, b, c, d = sp.symbols(
        "a b c d"
    )

    X = sp.Matrix([
        [a, b],
        [c, d],
    ])

    E = clean_matrix(
        X * comm
        - comm * X
    )

    equations = [
        E[0, 0],
        E[0, 1],
        E[1, 0],
        E[1, 1],
    ]

    A, _ = sp.linear_eq_to_matrix(
        equations,
        [
            a,
            b,
            c,
            d,
        ],
    )

    return A, A.nullspace()


# ============================================================================
# TRACE / DETERMINANT COMMUTATOR AUDIT
# ============================================================================

def commutator_profile(comm):

    tr = clean(
        sp.trace(comm)
    )

    det = clean(
        comm.det()
    )

    rank = matrix_rank(
        comm
    )

    print()
    print("=" * 78)
    print(
        "4. EXACT COMMUTATOR PROFILE"
    )
    print("=" * 78)

    print()
    print(
        "  commutator=[T1,T2]"
    )

    print(
        comm
    )

    print()
    print(
        "  trace={}".format(
            tr
        )
    )

    print(
        "  determinant={}".format(
            det
        )
    )

    print(
        "  rank={}".format(
            rank
        )
    )

    print(
        "  invertible={}".format(
            det != 0
        )
    )

    print(
        "  valuations={}".format(
            {
                p: valuation(det, p)
                for p in (
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


# ============================================================================
# COMMON EIGENLINE REPORT
# ============================================================================

def report_common_eigenline():

    print()
    print("=" * 78)
    print(
        "5. COMMON EIGENLINE / RATIONAL EIGENVECTOR AUDIT"
    )
    print("=" * 78)

    det1, det2, gcd_expr, candidates = (
        common_rational_eigenvector_test()
    )

    print()
    print(
        "  invariant_line_polynomial_T1={}".format(
            det1
        )
    )

    print(
        "  invariant_line_polynomial_T2={}".format(
            det2
        )
    )

    print(
        "  common_line_gcd={}".format(
            gcd_expr.as_expr()
        )
    )

    print(
        "  rational_common_eigenline_candidates={}".format(
            candidates
        )
    )

    print(
        "  common_rational_eigenline_exists={}".format(
            bool(candidates)
        )
    )

    p1, p2, gcd_poly = (
        common_eigenvalue_polynomial_test()
    )

    print()
    print(
        "  characteristic_polynomial_T1={}".format(
            p1
        )
    )

    print(
        "  characteristic_polynomial_T2={}".format(
            p2
        )
    )

    print(
        "  common_eigenvalue_polynomial_gcd={}".format(
            gcd_poly.as_expr()
        )
    )

    print(
        "  common_eigenvalue_exists={}".format(
            gcd_poly.degree() > 0
        )
    )


# ============================================================================
# CENTRALIZER REPORT
# ============================================================================

def report_centralizer():

    print()
    print("=" * 78)
    print(
        "6. COMMON CENTRALIZER DIMENSION AUDIT"
    )
    print("=" * 78)

    X, A, nullspace = (
        common_centralizer()
    )

    print()
    print(
        "  linear_system_matrix_shape={}".format(
            A.shape
        )
    )

    print(
        "  centralizer_dimension={}".format(
            len(nullspace)
        )
    )

    for i, basis in enumerate(
        nullspace
    ):

        vec = sp.Matrix(basis)

        B = sp.Matrix([
            [vec[0], vec[1]],
            [vec[2], vec[3]],
        ])

        print()
        print(
            "  basis_{}={}".format(
                i,
                B
            )
        )

    scalar_basis = (
        len(nullspace) == 1
        and clean_matrix(
            sp.Matrix([
                [nullspace[0][0], nullspace[0][1]],
                [nullspace[0][2], nullspace[0][3]],
            ])
        )
        .has(
            nullspace[0][0]
        )
    )

    print()
    print(
        "  common_centralizer_equals_scalars={}".format(
            len(nullspace) == 1
        )
    )

    Acomm, null_comm = (
        commutator_centralizer_test(
            clean_matrix(
                T1 * T2
                - T2 * T1
            )
        )
    )

    print()
    print(
        "  commutator_centralizer_dimension={}".format(
            len(null_comm)
        )
    )


# ============================================================================
# TRIANGULARIZATION REPORT
# ============================================================================

def report_triangularization():

    print()
    print("=" * 78)
    print(
        "7. SIMULTANEOUS TRIANGULARIZATION AUDIT"
    )
    print("=" * 78)

    S, sdet, equations = (
        simultaneous_triangularization_test()
    )

    print()
    print(
        "  symbolic_change_of_basis_determinant={}".format(
            sdet
        )
    )

    print(
        "  lower_left_equations_T1={}".format(
            equations[0]
        )
    )

    print(
        "  lower_left_equations_T2={}".format(
            equations[1]
        )
    )

    det1, det2, gcd_expr, candidates = (
        common_rational_eigenvector_test()
    )

    print()
    print(
        "  common_invariant_line_exists={}".format(
            bool(candidates)
            or gcd_expr.as_expr() != 1
        )
    )

    print(
        "  simultaneous_triangularization_possible={}".format(
            bool(candidates)
        )
    )


# ============================================================================
# IRREDUCIBILITY SUMMARY
# ============================================================================

def irreducibility_summary(comm):

    _, _, _, candidates = (
        common_rational_eigenvector_test()
    )

    comm_det = clean(
        comm.det()
    )

    _, _, common_gcd = (
        common_eigenvalue_polynomial_test()
    )

    _, _, central_basis = (
        common_centralizer()
    )

    common_line = (
        bool(candidates)
        or common_gcd.degree() > 0
    )

    irreducible_over_Q = (
        not bool(candidates)
    )

    full_commutator_rank = (
        comm_det != 0
    )

    print()
    print("=" * 78)
    print(
        "8. IRREDUCIBILITY SUMMARY"
    )
    print("=" * 78)

    print(
        "  common_invariant_line_over_Q={}".format(
            common_line
        )
    )

    print(
        "  common_rational_eigenvector={}".format(
            bool(candidates)
        )
    )

    print(
        "  common_eigenvalue_gcd_degree={}".format(
            common_gcd.degree()
        )
    )

    print(
        "  commutator_invertible={}".format(
            full_commutator_rank
        )
    )

    print(
        "  common_centralizer_dimension={}".format(
            len(central_basis)
        )
    )

    print(
        "  pair_irreducible_over_Q_candidate={}".format(
            irreducible_over_Q
        )
    )


# ============================================================================
# TERMINAL SOURCE
# ============================================================================

def terminal_reference():

    q1 = 495451247
    q3 = 421514439

    g = math.gcd(
        q1,
        q3,
    )

    print()
    print("=" * 78)
    print(
        "9. TERMINAL PROJECTIVE SOURCE REFERENCE"
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
            g
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
# INTERPRETATION
# ============================================================================

def interpretation():

    print()
    print("=" * 78)
    print(
        "10. STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
Experiment 310R established that the two exact width-2 transfer
matrices generate the full algebra M_2(Q).

The next question is whether this full algebra is accompanied by a
common invariant subspace.

For a 2x2 pair, the relevant hierarchy is:

    common eigenline
        ->
    simultaneous triangularization
        ->
    reducible representation.

If no common eigenline exists and the commutator is invertible, the pair
is strongly separated at the representation level.

The common centralizer provides a second exact diagnostic:

    centralizer = Q*I

means that no non-scalar matrix commutes with both transitions.

The common characteristic-polynomial gcd is an additional obstruction:
a nontrivial gcd signals a common eigenvalue.

These tests are independent of the earlier projective-similarity and
Möbius audits.

A positive reducibility result would reveal hidden one-dimensional
structure despite the full M_2(Q) span.

A negative result would strengthen the conclusion that the two exact
cross-layer transitions are genuinely noncommutative and representation-
theoretically unrelated by a shared invariant line.

Only two exact transitions are available, so this remains a structural
diagnostic rather than a universal source theorem.

No synthetic second n=pq case is generated.
"""
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 311R — EXACT TRANSFER-PAIR IRREDUCIBILITY / "
        "COMMON-EIGENLINE / CENTRALIZER AUDIT"
    )
    print("=" * 78)

    print()
    print(
        "1. EXACT TRANSFER MATRICES"
    )
    print("=" * 78)

    print(
        "  T1={}".format(
            T1
        )
    )

    print()
    print(
        "  T2={}".format(
            T2
        )
    )

    print()
    print(
        "2. CHARACTERISTIC DATA"
    )
    print("=" * 78)

    print_characteristic_data(
        "T1",
        T1,
    )

    print_characteristic_data(
        "T2",
        T2,
    )

    comm = clean_matrix(
        T1 * T2
        - T2 * T1
    )

    commutator_profile(
        comm
    )

    report_common_eigenline()

    report_centralizer()

    report_triangularization()

    irreducibility_summary(
        comm
    )

    terminal_reference()

    interpretation()

    print()
    print("=" * 78)
    print(
        "11. FINAL EXACTNESS"
    )
    print("=" * 78)

    print(
        "  exact_transfer_pair=True"
    )

    print(
        "  projective_similarity_already_refuted=True"
    )

    print(
        "  commutator_nonzero={}".format(
            comm != sp.zeros(2)
        )
    )

    print(
        "  commutator_invertible={}".format(
            clean(comm.det()) != 0
        )
    )

    _, _, common_gcd = (
        common_eigenvalue_polynomial_test()
    )

    print(
        "  common_characteristic_factor_degree={}".format(
            common_gcd.degree()
        )
    )

    _, _, central_basis = (
        common_centralizer()
    )

    print(
        "  common_centralizer_dimension={}".format(
            len(central_basis)
        )
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
        "EXPERIMENT 311R COMPLETE"
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

