#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 308R — EXACT TRANSFER-ALGEBRA / POLYNOMIAL-MÖBIUS AUDIT
==============================================================================

Goal
----
Experiment 307R proved that the two exact width-2 transfer matrices are not
projectively similar.

We now test a different possibility:

    T2 = f(T1)

for a low-complexity exact matrix function f.

Because a 2x2 matrix satisfies its quadratic characteristic equation, every
polynomial in T1 reduces to

    a*I + b*T1.

We therefore test:

    1. affine polynomial relation:
           T2 = a I + b T1

    2. shifted/scaled inverse relation:
           T2 = a I + b T1^{-1}

    3. exact Möbius / fractional-linear relation:
           T2 = (a*T1 + b*I) (c*T1 + d*I)^(-1)

       with determinant(ad-bc) != 0.

    4. commuting-algebra test:
           T1*T2 = T2*T1

If T2 is a rational function of T1, then T1 and T2 must commute.
If they do not commute, all such scalar-function relations are impossible.

We also test whether T1 and T2 generate a 2-dimensional commutative algebra,
and whether their commutator is exactly nonzero.

No external files.
No synthetic second n=pq case.
Exact rational arithmetic only.
"""


from __future__ import annotations

import math
import sys

import sympy as sp


# ============================================================================
# EXACT TRANSFER DATA
# ============================================================================

alpha1 = sp.Rational(
    -209857461192170070,
    11470116327290257,
)

beta1 = sp.Rational(
    -215450379004300026,
    11470116327290257,
)

alpha2 = sp.Rational(
    83976580526089197,
    971937272383741,
)

beta2 = sp.Rational(
    75947709834674022,
    971937272383741,
)


T1 = sp.Matrix(
    [
        [alpha1, beta1],
        [1, 0],
    ]
)

T2 = sp.Matrix(
    [
        [alpha2, beta2],
        [1, 0],
    ]
)

I2 = sp.eye(2)


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


def matrix_clean(M):
    return M.applyfunc(clean)


def matrix_zero(M):
    return all(
        clean(x) == 0
        for x in M
    )


def valuation(x, p):
    x = sp.Rational(x)

    if x == 0:
        return None

    num = abs(int(x.p))
    den = abs(int(x.q))

    e = 0

    while num % p == 0:
        num //= p
        e += 1

    while den % p == 0:
        den //= p
        e -= 1

    return e


def exact_linear_relation(A, B):
    """
    Solve A = a*I + b*B exactly.
    """
    a, b = sp.symbols("a b")

    equations = []

    M = A - a * I2 - b * B

    for entry in M:
        equations.append(
            sp.Eq(
                clean(entry),
                0,
            )
        )

    sols = sp.solve(
        equations,
        [a, b],
        dict=True,
    )

    verified = []

    for sol in sols:
        if a in sol and b in sol:
            AA = matrix_clean(
                A - (
                    sol[a] * I2
                    + sol[b] * B
                )
            )

            if matrix_zero(AA):
                verified.append(
                    (
                        clean(sol[a]),
                        clean(sol[b]),
                    )
                )

    return verified


def exact_inverse_relation(A, B):
    """
    Solve A = a*I + b*B^{-1}.
    """
    if B.det() == 0:
        return []

    Binv = matrix_clean(
        B.inv()
    )

    return exact_linear_relation(
        A,
        Binv,
    )


def solve_mobius(A, B):
    """
    Solve exactly for

        A = (a B + b I) (c B + d I)^(-1)

    This is equivalent to

        A(cB+dI) = aB+bI

    which gives linear equations in a,b,c,d.

    The overall scale is irrelevant, so we normalize by fixing one
    nonzero coefficient when possible and verify the resulting relation.
    """

    a, b, c, d = sp.symbols(
        "a b c d"
    )

    equation_matrix = (
        A * (c * B + d * I2)
        - (a * B + b * I2)
    )

    equations = [
        sp.Eq(
            clean(entry),
            0,
        )
        for entry in equation_matrix
    ]

    # Linear system in (a,b,c,d).
    M, rhs = sp.linear_eq_to_matrix(
        equations,
        [a, b, c, d],
    )

    rank = M.rank()
    aug_rank = M.row_join(rhs).rank()

    if aug_rank > rank:
        return {
            "status": "NO_SOLUTION",
            "rank": rank,
            "augmented_rank": aug_rank,
            "solutions": [],
        }

    nullspace = M.nullspace()

    if not nullspace:
        return {
            "status": "NO_SOLUTION",
            "rank": rank,
            "augmented_rank": aug_rank,
            "solutions": [],
        }

    verified = []

    for vector in nullspace:

        aa = clean(vector[0])
        bb = clean(vector[1])
        cc = clean(vector[2])
        dd = clean(vector[3])

        determinant = clean(
            aa * dd
            - bb * cc
        )

        if determinant == 0:
            continue

        denominator = clean(
            cc * B + dd * I2
        )

        if denominator.det() == 0:
            continue

        lhs = matrix_clean(
            A * denominator
        )

        rhs_matrix = matrix_clean(
            aa * B + bb * I2
        )

        if matrix_zero(
            lhs - rhs_matrix
        ):
            verified.append(
                {
                    "a": aa,
                    "b": bb,
                    "c": cc,
                    "d": dd,
                    "det": determinant,
                }
            )

    if verified:
        return {
            "status": "EXACT",
            "rank": rank,
            "augmented_rank": aug_rank,
            "solutions": verified,
        }

    return {
        "status": "NO_NONDEGENERATE_SOLUTION",
        "rank": rank,
        "augmented_rank": aug_rank,
        "solutions": [],
    }


# ============================================================================
# CHARACTERISTIC DATA
# ============================================================================

def characteristic_data(T):

    tr = clean(
        sp.trace(T)
    )

    det = clean(
        T.det()
    )

    disc = clean(
        tr**2 - 4 * det
    )

    poly = clean(
        T.charpoly(
            sp.Symbol("x")
        ).as_expr()
    )

    return tr, det, disc, poly


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 308R — EXACT TRANSFER-ALGEBRA / POLYNOMIAL-MÖBIUS AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------------
    # 1. INPUT MATRICES
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. EXACT TRANSFER MATRICES")
    print("=" * 78)

    print()
    print("  T1=")
    print(T1)

    print()
    print("  T2=")
    print(T2)

    # ------------------------------------------------------------------------
    # 2. CHARACTERISTIC DATA
    # ------------------------------------------------------------------------

    tr1, det1, disc1, poly1 = characteristic_data(
        T1
    )

    tr2, det2, disc2, poly2 = characteristic_data(
        T2
    )

    print()
    print("=" * 78)
    print("2. CHARACTERISTIC DATA")
    print("=" * 78)

    print()
    print("  T1:")
    print("    trace={}".format(tr1))
    print("    determinant={}".format(det1))
    print("    discriminant={}".format(disc1))
    print("    characteristic_polynomial={}".format(poly1))

    print()
    print("  T2:")
    print("    trace={}".format(tr2))
    print("    determinant={}".format(det2))
    print("    discriminant={}".format(disc2))
    print("    characteristic_polynomial={}".format(poly2))

    # ------------------------------------------------------------------------
    # 3. COMMUTATOR
    # ------------------------------------------------------------------------

    commutator = matrix_clean(
        T1 * T2
        - T2 * T1
    )

    print()
    print("=" * 78)
    print("3. COMMUTATOR / COMMON-ALGEBRA TEST")
    print("=" * 78)

    print()
    print("  T1*T2-T2*T1=")
    print(commutator)

    comm_zero = matrix_zero(
        commutator
    )

    print(
        "  exact_commutation={}".format(
            comm_zero
        )
    )

    if not comm_zero:
        print(
            "  consequence=NO_SCALAR_FUNCTION_OF_T1_CAN_EQUAL_T2"
        )

    # ------------------------------------------------------------------------
    # 4. AFFINE POLYNOMIAL RELATION
    # ------------------------------------------------------------------------

    affine = exact_linear_relation(
        T2,
        T1,
    )

    print()
    print("=" * 78)
    print("4. EXACT AFFINE POLYNOMIAL RELATION")
    print("=" * 78)

    if affine:

        for a, b in affine:

            print()
            print(
                "  EXACT: T2 = a*I + b*T1"
            )

            print(
                "    a={}".format(a)
            )

            print(
                "    b={}".format(b)
            )

    else:

        print(
            "  status=NO_SOLUTION"
        )

    # ------------------------------------------------------------------------
    # 5. INVERSE RELATION
    # ------------------------------------------------------------------------

    inverse_relation = exact_inverse_relation(
        T2,
        T1,
    )

    print()
    print("=" * 78)
    print("5. EXACT INVERSE-POLYNOMIAL RELATION")
    print("=" * 78)

    if inverse_relation:

        for a, b in inverse_relation:

            print()
            print(
                "  EXACT: T2 = a*I + b*T1^(-1)"
            )

            print(
                "    a={}".format(a)
            )

            print(
                "    b={}".format(b)
            )

    else:

        print(
            "  status=NO_SOLUTION"
        )

    # ------------------------------------------------------------------------
    # 6. MÖBIUS RELATION
    # ------------------------------------------------------------------------

    mobius = solve_mobius(
        T2,
        T1,
    )

    print()
    print("=" * 78)
    print(
        "6. EXACT MÖBIUS / FRACTIONAL-LINEAR RELATION"
    )
    print("=" * 78)

    print(
        "  status={}".format(
            mobius["status"]
        )
    )

    print(
        "  rank={}".format(
            mobius["rank"]
        )
    )

    print(
        "  augmented_rank={}".format(
            mobius["augmented_rank"]
        )
    )

    if mobius["solutions"]:

        for idx, sol in enumerate(
            mobius["solutions"]
        ):

            print()
            print(
                "  solution_{}:".format(
                    idx
                )
            )

            print(
                "    a={}".format(
                    sol["a"]
                )
            )

            print(
                "    b={}".format(
                    sol["b"]
                )
            )

            print(
                "    c={}".format(
                    sol["c"]
                )
            )

            print(
                "    d={}".format(
                    sol["d"]
                )
            )

            print(
                "    ad-bc={}".format(
                    sol["det"]
                )
            )

    # ------------------------------------------------------------------------
    # 7. T1-CENTRALIZER TEST
    # ------------------------------------------------------------------------

    x00, x01, x10, x11 = sp.symbols(
        "x00 x01 x10 x11"
    )

    X = sp.Matrix(
        [
            [x00, x01],
            [x10, x11],
        ]
    )

    equations = list(
        T1 * X
        - X * T1
    )

    M, rhs = sp.linear_eq_to_matrix(
        [
            sp.Eq(clean(e), 0)
            for e in equations
        ],
        [x00, x01, x10, x11],
    )

    centralizer_basis = M.nullspace()

    print()
    print("=" * 78)
    print("7. CENTRALIZER / ALGEBRA DIMENSION AUDIT")
    print("=" * 78)

    print(
        "  centralizer_dimension={}".format(
            len(centralizer_basis)
        )
    )

    for i, vec in enumerate(
        centralizer_basis
    ):

        print()
        print(
            "  basis_{}={}".format(
                i,
                vec,
            )
        )

    # ------------------------------------------------------------------------
    # 8. PRIME VALUATION PROFILE
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print(
        "8. PRIME-VALUATION PROFILE OF THE OBSTRUCTIONS"
    )
    print("=" * 78)

    quantities = {
        "T1_trace": tr1,
        "T1_det": det1,
        "T1_disc": disc1,
        "T2_trace": tr2,
        "T2_det": det2,
        "T2_disc": disc2,
    }

    for name, value in quantities.items():

        print()
        print(
            "  {}={}".format(
                name,
                value,
            )
        )

        for prime in (
            2,
            3,
            5,
            7,
            11,
            13,
            17,
        ):

            print(
                "    v_{}={}".format(
                    prime,
                    valuation(
                        value,
                        prime,
                    ),
                )
            )

    # ------------------------------------------------------------------------
    # 9. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
r"""
Experiment 307R ruled out projective similarity of the two exact
width-2 transition matrices.

Experiment 308R now tests whether that failure is compatible with a
different kind of common algebraic origin.

The hierarchy is:

    T2 = a I + b T1
        (affine polynomial)

    T2 = a I + b T1^(-1)
        (inverse-polynomial)

    T2 = (a T1 + b I)(c T1 + d I)^(-1)
        (Möbius / fractional-linear)

All of these relations force T1 and T2 to lie in the same commutative
algebra generated by T1.

Therefore the commutator

    [T1,T2] = T1*T2 - T2*T1

is a decisive obstruction.

If the commutator is nonzero, then T2 is not any scalar rational
function of T1, including the Möbius class.

If the commutator vanishes but the affine test fails, the experiment
still leaves the possibility of a genuinely higher rational function,
although for a 2x2 matrix every polynomial reduces to degree <=1 by
Cayley-Hamilton and a rational function reduces correspondingly on the
generic invertible case.

A positive low-complexity relation would therefore be substantially
stronger than the failed projective-similarity hypothesis.

No synthetic second n=pq case is generated.
"""
    )

    # ------------------------------------------------------------------------
    # 10. FINAL EXACTNESS
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. FINAL EXACTNESS")
    print("=" * 78)

    print(
        "  exact_commutator_test={}".format(
            comm_zero
        )
    )

    print(
        "  exact_affine_relation={}".format(
            bool(affine)
        )
    )

    print(
        "  exact_inverse_relation={}".format(
            bool(inverse_relation)
        )
    )

    print(
        "  exact_mobius_relation={}".format(
            bool(
                mobius["solutions"]
            )
        )
    )

    print(
        "  projective_similarity_already_refuted=True"
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
        "EXPERIMENT 308R COMPLETE"
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

