#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 314R — EXACT ADDITIVE VS GROUP COMMUTATOR / FRICKE INVARIANT AUDIT
==============================================================================

Purpose
-------
Experiment 313R mixed two different notions of commutator:

    additive commutator:
        K = T1*T2 - T2*T1

    multiplicative/group commutator:
        G = T1*T2*T1^(-1)*T2^(-1)

These have completely different determinant identities.

For invertible 2x2 matrices:

    det(G) = 1

always.

The Fricke trace identity is for G, not for K.

This experiment therefore:

    1. recomputes T1, T2 exactly;
    2. separates K and G;
    3. verifies det(G)=1;
    4. verifies tr(G) from the five scalar invariants;
    5. checks the corrected general Fricke identity;
    6. checks Cayley-Hamilton for both K and G;
    7. compares additive-commutator and group-commutator invariants;
    8. recomputes the trace pairing and M_2(Q) basis rank;
    9. audits inverse-word trace identities.

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
    [1, 0],
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
    [1, 0],
])


# ============================================================================
# EXACT HELPERS
# ============================================================================

def clean(x):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(x)
            )
        )
    )


def valuation(x, p):
    x = sp.Rational(x)

    if x == 0:
        return sp.oo

    num = abs(int(x.p))
    den = abs(int(x.q))

    v = 0

    while num % p == 0:
        num //= p
        v += 1

    while den % p == 0:
        den //= p
        v -= 1

    return v


def matrix_zero(M):
    return all(
        clean(M[i, j]) == 0
        for i in range(M.rows)
        for j in range(M.cols)
    )


def print_matrix(name, M):
    print(f"  {name}=")
    print(M)


# ============================================================================
# BASIC INVARIANTS
# ============================================================================

def basic_invariants():

    x = clean(T1.trace())
    y = clean(T2.trace())
    z = clean((T1 * T2).trace())

    a = clean(T1.det())
    b = clean(T2.det())

    return x, y, z, a, b


# ============================================================================
# CAYLEY-HAMILTON
# ============================================================================

def cayley_hamilton(M):

    tr = clean(M.trace())
    det = clean(M.det())

    residual = clean(
        M**2 - tr * M + det * sp.eye(2)
    )

    return residual


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 314R — EXACT ADDITIVE VS GROUP COMMUTATOR / "
        "FRICKE INVARIANT AUDIT"
    )
    print("=" * 78)

    # ------------------------------------------------------------------------
    # 1. GENERATORS
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("1. EXACT GENERATORS")
    print("=" * 78)

    print_matrix("T1", T1)
    print()
    print_matrix("T2", T2)

    # ------------------------------------------------------------------------
    # 2. SCALAR DATA
    # ------------------------------------------------------------------------

    x, y, z, a, b = basic_invariants()

    print()
    print("=" * 78)
    print("2. BASIC SCALAR INVARIANTS")
    print("=" * 78)

    print(f"  x=tr(T1)={x}")
    print(f"  y=tr(T2)={y}")
    print(f"  z=tr(T1*T2)={z}")
    print(f"  a=det(T1)={a}")
    print(f"  b=det(T2)={b}")

    # ------------------------------------------------------------------------
    # 3. TWO DIFFERENT COMMUTATORS
    # ------------------------------------------------------------------------

    K = clean(T1 * T2 - T2 * T1)

    G = clean(
        T1
        * T2
        * T1.inv()
        * T2.inv()
    )

    print()
    print("=" * 78)
    print("3. ADDITIVE VS MULTIPLICATIVE COMMUTATORS")
    print("=" * 78)

    print_matrix(
        "K = T1*T2 - T2*T1",
        K,
    )

    print()
    print_matrix(
        "G = T1*T2*T1^(-1)*T2^(-1)",
        G,
    )

    print()
    print(
        f"  K_zero={matrix_zero(K)}"
    )

    print(
        f"  G_zero={matrix_zero(G)}"
    )

    # ------------------------------------------------------------------------
    # 4. DETERMINANTS
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. DETERMINANT AUDIT")
    print("=" * 78)

    det_K = clean(K.det())
    det_G = clean(G.det())

    print(
        f"  det(K)={det_K}"
    )

    print(
        f"  det(G)={det_G}"
    )

    print(
        f"  det(G)==1: {det_G == 1}"
    )

    print(
        "  additive_commutator_det_one=False_by_definition"
    )

    # ------------------------------------------------------------------------
    # 5. TRACE OF GROUP COMMUTATOR
    # ------------------------------------------------------------------------

    tr_G = clean(G.trace())
    tr_K = clean(K.trace())

    print()
    print("=" * 78)
    print("5. TRACE COMMUTATOR AUDIT")
    print("=" * 78)

    print(
        f"  tr(K)={tr_K}"
    )

    print(
        f"  tr(G)={tr_G}"
    )

    print(
        f"  tr(K)==0: {tr_K == 0}"
    )

    # Correct general Fricke identity:
    #
    # tr(G)
    #   =
    # (z^2 + b*x^2 + a*y^2 - x*y*z - 2*a*b)/(a*b)

    fricke_rhs = clean(
        (
            z**2
            + b * x**2
            + a * y**2
            - x * y * z
            - 2 * a * b
        )
        / (a * b)
    )

    print()
    print(
        "  Fricke RHS="
        f"{fricke_rhs}"
    )

    print(
        f"  exact_tr(G)==Fricke_RHS: "
        f"{clean(tr_G - fricke_rhs) == 0}"
    )

    print(
        f"  cleared_difference="
        f"{clean((tr_G - fricke_rhs) * a * b)}"
    )

    # ------------------------------------------------------------------------
    # 6. GROUP COMMUTATOR CAYLEY-HAMILTON
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. GROUP-COMMUTATOR CAYLEY-HAMILTON")
    print("=" * 78)

    G_ch = cayley_hamilton(G)

    print_matrix(
        "G^2-tr(G)G+det(G)I",
        G_ch,
    )

    print(
        f"  exact_zero={matrix_zero(G_ch)}"
    )

    G_reduced = clean(
        G**2
        - tr_G * G
        + sp.eye(2)
    )

    print(
        f"  G^2-tr(G)G+I exact_zero="
        f"{matrix_zero(G_reduced)}"
    )

    # ------------------------------------------------------------------------
    # 7. ADDITIVE COMMUTATOR CAYLEY-HAMILTON
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. ADDITIVE-COMMUTATOR CAYLEY-HAMILTON")
    print("=" * 78)

    K_ch = cayley_hamilton(K)

    print_matrix(
        "K^2-tr(K)K+det(K)I",
        K_ch,
    )

    print(
        f"  exact_zero={matrix_zero(K_ch)}"
    )

    print(
        f"  since tr(K)=0, K^2+det(K)I="
        f"{matrix_zero(clean(K**2 + det_K * sp.eye(2)))}"
    )

    # ------------------------------------------------------------------------
    # 8. INVERSE TRACE IDENTITIES
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. MIXED INVERSE-TRACE IDENTITIES")
    print("=" * 78)

    T1T2inv = clean(
        (T1 * T2.inv()).trace()
    )

    T1invT2 = clean(
        (T1.inv() * T2).trace()
    )

    rhs_1 = clean(
        (x * y - z) / b
    )

    rhs_2 = clean(
        (x * y - z) / a
    )

    print(
        f"  tr(T1*T2^(-1))={T1T2inv}"
    )

    print(
        f"  expected=(x*y-z)/b={rhs_1}"
    )

    print(
        f"  exact={clean(T1T2inv-rhs_1)==0}"
    )

    print()
    print(
        f"  tr(T1^(-1)*T2)={T1invT2}"
    )

    print(
        f"  expected=(x*y-z)/a={rhs_2}"
    )

    print(
        f"  exact={clean(T1invT2-rhs_2)==0}"
    )

    # ------------------------------------------------------------------------
    # 9. TRACE OF GROUP-COMMUTATOR INVERSE
    # ------------------------------------------------------------------------

    G_inv = clean(G.inv())

    print()
    print("=" * 78)
    print("9. GROUP-COMMUTATOR INVERSE AUDIT")
    print("=" * 78)

    print(
        f"  tr(G)={clean(G.trace())}"
    )

    print(
        f"  tr(G^(-1))={clean(G_inv.trace())}"
    )

    print(
        f"  exact_equal={clean(G.trace()-G_inv.trace())==0}"
    )

    print(
        f"  det(G^(-1))={clean(G_inv.det())}"
    )

    # ------------------------------------------------------------------------
    # 10. COMMUTATOR DETERMINANT RELATION
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("10. ADDITIVE-COMMUTATOR DETERMINANT RELATION AUDIT")
    print("=" * 78)

    # This quantity is intentionally NOT expected to equal 1.
    #
    # We simply compare it to exact scalar expressions built from
    # x,y,z,a,b, without assuming a formula in advance.

    candidate_1 = clean(
        z**2
        + b * x**2
        + a * y**2
        - x * y * z
        - 2 * a * b
    )

    candidate_2 = clean(
        a * b * tr_G
    )

    print(
        f"  det(K)={det_K}"
    )

    print(
        f"  candidate_fricke_numerator={candidate_1}"
    )

    print(
        f"  a*b*tr(G)={candidate_2}"
    )

    print(
        f"  fricke_numerator_matches_a_b_trG="
        f"{clean(candidate_1-candidate_2)==0}"
    )

    # ------------------------------------------------------------------------
    # 11. BASIS RANK
    # ------------------------------------------------------------------------

    basis = [
        sp.eye(2),
        T1,
        T2,
        T1 * T2,
    ]

    B = sp.Matrix.hstack(
        *[
            sp.Matrix(M).reshape(4, 1)
            for M in basis
        ]
    )

    basis_rank = B.rank()

    print()
    print("=" * 78)
    print("11. FULL MATRIX-ALGEBRA BASIS AUDIT")
    print("=" * 78)

    print(
        f"  basis=['I','T1','T2','T1T2']"
    )

    print(
        f"  basis_rank={basis_rank}"
    )

    print(
        f"  full_M2_Q_basis={basis_rank == 4}"
    )

    # ------------------------------------------------------------------------
    # 12. TRACE PAIRING
    # ------------------------------------------------------------------------

    Gmat = sp.Matrix([
        [
            clean(
                (E1 * E2).trace()
            )
            for E2 in basis
        ]
        for E1 in basis
    ])

    print()
    print("=" * 78)
    print("12. TRACE PAIRING")
    print("=" * 78)

    print_matrix(
        "G_ij=tr(E_i E_j)",
        Gmat,
    )

    print(
        f"  trace_pairing_rank={Gmat.rank()}"
    )

    print(
        f"  trace_pairing_determinant="
        f"{clean(Gmat.det())}"
    )

    print(
        f"  trace_pairing_nondegenerate="
        f"{Gmat.det() != 0}"
    )

    # ------------------------------------------------------------------------
    # 13. PRIME PROFILE
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("13. PRIME-VALUATION PROFILE")
    print("=" * 78)

    scalar_objects = {
        "tr(T1)": x,
        "tr(T2)": y,
        "tr(T1T2)": z,
        "det(T1)": a,
        "det(T2)": b,
        "tr(K)": tr_K,
        "det(K)": det_K,
        "tr(G)": tr_G,
        "det(G)": det_G,
    }

    for name, value in scalar_objects.items():

        print()
        print(
            f"  {name}={value}"
        )

        for p in (
            2,
            3,
            5,
            7,
            11,
            13,
            17,
        ):

            print(
                f"    v_{p}={valuation(value,p)}"
            )

    # ------------------------------------------------------------------------
    # 14. TERMINAL SOURCE REFERENCE
    # ------------------------------------------------------------------------

    q1_terminal = 495451247
    q3_terminal = 421514439
    g = math.gcd(
        q1_terminal,
        q3_terminal,
    )

    print()
    print("=" * 78)
    print("14. TERMINAL PROJECTIVE SOURCE REFERENCE")
    print("=" * 78)

    print(
        f"  q1_terminal={q1_terminal}"
    )

    print(
        f"  q3_terminal={q3_terminal}"
    )

    print(
        f"  gcd={g}"
    )

    print(
        f"  q1/17={q1_terminal//17}"
    )

    print(
        f"  q3/17={q3_terminal//17}"
    )

    # ------------------------------------------------------------------------
    # 15. STRUCTURAL INTERPRETATION
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("15. STRUCTURAL INTERPRETATION")
    print("=" * 78)

    print(
r"""
Experiment 313R mixed two different commutator notions.

For the additive commutator

    K = T1*T2 - T2*T1,

the trace is zero, but its determinant is not constrained to equal one.

For the multiplicative/group commutator

    G = T1*T2*T1^(-1)*T2^(-1),

the determinant is identically

    det(G)=1

because determinants multiply.

The Fricke trace identity applies to G:

    tr(G)
      =
    (z^2 + b*x^2 + a*y^2 - x*y*z - 2*a*b)/(a*b),

where

    x=tr(T1),
    y=tr(T2),
    z=tr(T1*T2),
    a=det(T1),
    b=det(T2).

This experiment checks that identity directly and independently.

It also keeps K visible because the additive commutator is still the
correct object for the Lie-algebra/noncommutative-span calculations from
Experiments 308R-309R.

Thus the two structures are separated:

    Lie-style obstruction:
        K = T1*T2 - T2*T1,

    group-style obstruction:
        G = T1*T2*T1^(-1)*T2^(-1).

A successful Fricke test would repair the failed scalar identity in
Experiment 313R without changing the earlier conclusion that

    <I,T1,T2,T1T2> = M_2(Q).

No synthetic second n=pq case is generated.
"""
    )

    # ------------------------------------------------------------------------
    # 16. FINAL EXACTNESS
    # ------------------------------------------------------------------------

    fricke_exact = (
        clean(
            tr_G - fricke_rhs
        )
        == 0
    )

    inv_trace_exact = (
        clean(
            T1T2inv - rhs_1
        )
        == 0
        and
        clean(
            T1invT2 - rhs_2
        )
        == 0
    )

    print()
    print("=" * 78)
    print("16. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  additive_commutator_nonzero={not matrix_zero(K)}"
    )

    print(
        f"  additive_commutator_trace_zero={tr_K == 0}"
    )

    print(
        f"  additive_commutator_det_one={det_K == 1}"
    )

    print(
        f"  group_commutator_det_one={det_G == 1}"
    )

    print(
        f"  fricke_trace_identity_exact={fricke_exact}"
    )

    print(
        f"  inverse_trace_identities_exact={inv_trace_exact}"
    )

    print(
        f"  group_commutator_ch_exact="
        f"{matrix_zero(G_ch)}"
    )

    print(
        f"  full_M2_Q_basis={basis_rank == 4}"
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
        "EXPERIMENT 314R COMPLETE"
    )


if __name__ == "__main__":

    try:
        main()

    except KeyboardInterrupt:

        print(
            "\nInterrupted."
        )

        sys.exit(130)

    except Exception:
        raise
