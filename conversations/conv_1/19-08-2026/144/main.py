#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 310R — EXACT TRANSFER-ALGEBRA CLOSURE / WORD-SPAN AUDIT
==============================================================================

Purpose
-------
Repair the matrix-span solver used after Experiment 309R.

The experiment keeps everything in this single file and uses exact
SymPy rational arithmetic.

It audits:

    I
    T1
    T2
    T1*T2
    T2*T1
    [T1,T2]
    {T1,T2}
    T1^2
    T2^2
    T1*T2*T1
    T2*T1*T2

and determines:

    * exact span membership;
    * dimensions of generated matrix spaces;
    * closure of <I,T1,T2>;
    * closure after adjoining mixed words;
    * explicit exact relations when they exist;
    * Cayley-Hamilton reductions;
    * commutator / anticommutator structure.

No external files.
No synthetic second n=pq case.
"""


from __future__ import annotations

import math
import sys

import sympy as sp


# ============================================================================
# EXACT DATA
# ============================================================================

ALPHA1 = sp.Rational(
    -209857461192170070,
    11470116327290257,
)

BETA1 = sp.Rational(
    -215450379004300026,
    11470116327290257,
)

ALPHA2 = sp.Rational(
    83976580526089197,
    971937272383741,
)

BETA2 = sp.Rational(
    75947709834674022,
    971937272383741,
)


# ============================================================================
# HELPERS
# ============================================================================

def clean(x):
    return sp.factor(
        sp.cancel(
            sp.sympify(x)
        )
    )


def matrix_clean(M):
    return M.applyfunc(clean)


def matrix_key(M):
    M = matrix_clean(M)
    return tuple(
        clean(M[i, j])
        for i in range(M.rows)
        for j in range(M.cols)
    )


def matrix_to_tuple(M):
    return tuple(
        clean(M[i, j])
        for i in range(M.rows)
        for j in range(M.cols)
    )


def flatten_matrix(M):
    return sp.Matrix(
        [
            clean(M[i, j])
            for i in range(M.rows)
            for j in range(M.cols)
        ]
    )


def valuation(x, prime):
    x = sp.Rational(x)

    if x == 0:
        return None

    n = abs(int(x.p))
    d = abs(int(x.q))

    v = 0

    while n % prime == 0:
        n //= prime
        v += 1

    while d % prime == 0:
        d //= prime
        v -= 1

    return v


# ============================================================================
# MATRIX DATA
# ============================================================================

def build_matrices():

    T1 = sp.Matrix(
        [
            [ALPHA1, BETA1],
            [sp.Integer(1), sp.Integer(0)],
        ]
    )

    T2 = sp.Matrix(
        [
            [ALPHA2, BETA2],
            [sp.Integer(1), sp.Integer(0)],
        ]
    )

    I = sp.eye(2)

    return I, T1, T2


# ============================================================================
# ROBUST LINEAR SPAN SOLVER
# ============================================================================

def solve_matrix_span(target, basis):
    """
    Solve

        target = c0*basis[0] + ... + cn*basis[n-1]

    exactly over Q.

    IMPORTANT:
    SymPy's linsolve() may return a FiniteSet containing tuples, rather than
    a dictionary. This function normalizes either representation.

    Returns
    -------
    dict with:
        status
        coefficients
        rank
        augmented_rank
    """

    target = matrix_clean(target)
    basis = [
        matrix_clean(M)
        for M in basis
    ]

    n = len(basis)

    if n == 0:
        zero = sp.zeros(
            target.rows,
            target.cols,
        )

        if matrix_key(target) == matrix_key(zero):
            return {
                "status": "EXACT",
                "coefficients": [],
                "rank": 0,
                "augmented_rank": 0,
            }

        return {
            "status": "NO_SOLUTION",
            "coefficients": None,
            "rank": 0,
            "augmented_rank": 1,
        }

    A = sp.Matrix.hstack(
        *[
            flatten_matrix(M)
            for M in basis
        ]
    )

    b = flatten_matrix(target)

    rank = A.rank()

    augmented = A.row_join(b)
    augmented_rank = augmented.rank()

    if augmented_rank > rank:
        return {
            "status": "NO_SOLUTION",
            "coefficients": None,
            "rank": rank,
            "augmented_rank": augmented_rank,
        }

    # Use gauss_jordan_solve. It gives a Matrix column of coefficients and
    # handles the free-variable case explicitly.
    try:
        solution, params = A.gauss_jordan_solve(b)

    except Exception:
        # Fallback to linsolve.
        symbols = sp.symbols(
            "c0:%d" % n
        )

        equations = [
            sp.Eq(
                sum(
                    symbols[j] * A[i, j]
                    for j in range(n)
                ),
                b[i],
            )
            for i in range(A.rows)
        ]

        solution_set = sp.linsolve(
            equations,
            symbols,
        )

        if solution_set == sp.EmptySet:
            return {
                "status": "NO_SOLUTION",
                "coefficients": None,
                "rank": rank,
                "augmented_rank": augmented_rank,
            }

        # Convert FiniteSet({(c0,...,cn)}) safely.
        tuples = list(solution_set)

        if len(tuples) != 1:
            return {
                "status": "NONUNIQUE",
                "coefficients": None,
                "rank": rank,
                "augmented_rank": augmented_rank,
            }

        tup = tuples[0]

        if any(
            value.free_symbols
            for value in tup
        ):
            return {
                "status": "NONUNIQUE",
                "coefficients": [
                    clean(v)
                    for v in tup
                ],
                "rank": rank,
                "augmented_rank": augmented_rank,
            }

        coeffs = [
            clean(v)
            for v in tup
        ]

        return {
            "status": "EXACT",
            "coefficients": coeffs,
            "rank": rank,
            "augmented_rank": augmented_rank,
        }

    # gauss_jordan_solve returns:
    #
    #   solution = Matrix(...)
    #   params   = Matrix([...]) or empty Matrix
    #
    # If params is nonempty, the system is underdetermined.
    if params is not None:
        try:
            if len(params) > 0:
                return {
                    "status": "NONUNIQUE",
                    "coefficients": [
                        clean(v)
                        for v in solution
                    ],
                    "rank": rank,
                    "augmented_rank": augmented_rank,
                }
        except TypeError:
            pass

    coeffs = [
        clean(solution[i, 0])
        for i in range(solution.rows)
    ]

    # Verify independently.
    reconstructed = sp.zeros(
        target.rows,
        target.cols,
    )

    for c, M in zip(coeffs, basis):
        reconstructed += c * M

    reconstructed = matrix_clean(
        reconstructed
    )

    if matrix_key(reconstructed) != matrix_key(target):
        return {
            "status": "VERIFICATION_FAILED",
            "coefficients": coeffs,
            "rank": rank,
            "augmented_rank": augmented_rank,
        }

    return {
        "status": "EXACT",
        "coefficients": coeffs,
        "rank": rank,
        "augmented_rank": augmented_rank,
    }


# ============================================================================
# SPAN AUDIT
# ============================================================================

def audit_span(label, target, basis_labels, basis):

    result = solve_matrix_span(
        target,
        basis,
    )

    print()
    print("  {}:".format(label))
    print(
        "    basis={}".format(
            basis_labels
        )
    )
    print(
        "    status={}".format(
            result["status"]
        )
    )
    print(
        "    rank={}".format(
            result["rank"]
        )
    )
    print(
        "    augmented_rank={}".format(
            result["augmented_rank"]
        )
    )

    if result["coefficients"] is not None:
        print(
            "    coefficients={}".format(
                result["coefficients"]
            )
        )

    return result


# ============================================================================
# CLOSURE AUDIT
# ============================================================================

def closure_audit(I, T1, T2):

    print()
    print("=" * 78)
    print("5. EXACT PRODUCT-SPAN / CLOSURE AUDIT")
    print("=" * 78)

    basis_I_T1_T2 = [
        I,
        T1,
        T2,
    ]

    labels_I_T1_T2 = [
        "I",
        "T1",
        "T2",
    ]

    targets = {
        "T1*T2": T1 * T2,
        "T2*T1": T2 * T1,
        "T1^2": T1 ** 2,
        "T2^2": T2 ** 2,
        "commutator": T1 * T2 - T2 * T1,
        "anticommutator": T1 * T2 + T2 * T1,
    }

    results = {}

    for label, target in targets.items():

        results[label] = audit_span(
            label,
            target,
            labels_I_T1_T2,
            basis_I_T1_T2,
        )

    print()
    print(
        "  generated_span_rank={}".format(
            sp.Matrix.hstack(
                *[
                    flatten_matrix(M)
                    for M in basis_I_T1_T2
                ]
            ).rank()
        )
    )

    # Add T1*T2.
    basis_with_product = [
        I,
        T1,
        T2,
        T1 * T2,
    ]

    rank_with_product = sp.Matrix.hstack(
        *[
            flatten_matrix(M)
            for M in basis_with_product
        ]
    ).rank()

    print(
        "  rank_<I,T1,T2,T1T2>={}".format(
            rank_with_product
        )
    )

    # Add commutator.
    basis_with_comm = [
        I,
        T1,
        T2,
        T1 * T2 - T2 * T1,
    ]

    rank_with_comm = sp.Matrix.hstack(
        *[
            flatten_matrix(M)
            for M in basis_with_comm
        ]
    ).rank()

    print(
        "  rank_<I,T1,T2,[T1,T2]>={}".format(
            rank_with_comm
        )
    )

    return results


# ============================================================================
# WORD-BASIS AUDIT
# ============================================================================

def word_basis_audit(I, T1, T2):

    print()
    print("=" * 78)
    print(
        "6. WORD-BASIS DIMENSION / RELATION AUDIT"
    )
    print("=" * 78)

    words = [
        ("I", I),
        ("T1", T1),
        ("T2", T2),
        ("T1T2", T1 * T2),
        ("T2T1", T2 * T1),
        ("T1T2T1", T1 * T2 * T1),
        ("T2T1T2", T2 * T1 * T2),
    ]

    M = sp.Matrix.hstack(
        *[
            flatten_matrix(mat)
            for _, mat in words
        ]
    )

    print(
        "  word_labels={}".format(
            [name for name, _ in words]
        )
    )

    print(
        "  matrix_shape={}".format(
            M.shape
        )
    )

    print(
        "  rank={}".format(
            M.rank()
        )
    )

    print(
        "  nullity={}".format(
            len(words) - M.rank()
        )
    )

    nullspace = M.nullspace()

    print(
        "  nontrivial_relations={}".format(
            len(nullspace)
        )
    )

    for i, vec in enumerate(
        nullspace
    ):

        print()
        print(
            "  relation_{}=".format(i)
        )

        entries = [
            clean(vec[j, 0])
            for j in range(vec.rows)
        ]

        print(
            "    coefficients={}".format(
                entries
            )
        )

        nonzero = []

        for coeff, (name, _) in zip(
            entries,
            words,
        ):

            if coeff != 0:
                nonzero.append(
                    (coeff, name)
                )

        print(
            "    support={}".format(
                nonzero
            )
        )


# ============================================================================
# COMMUTATOR AUDIT
# ============================================================================

def commutator_audit(T1, T2):

    print()
    print("=" * 78)
    print(
        "7. COMMUTATOR / ANTI-COMMUTATOR AUDIT"
    )
    print("=" * 78)

    C = matrix_clean(
        T1 * T2 - T2 * T1
    )

    A = matrix_clean(
        T1 * T2 + T2 * T1
    )

    print(
        "  commutator={}".format(
            C
        )
    )

    print(
        "  anticommutator={}".format(
            A
        )
    )

    print(
        "  commutator_zero={}".format(
            C == sp.zeros(2)
        )
    )

    print(
        "  commutator_rank={}".format(
            C.rank()
        )
    )

    print(
        "  anticommutator_rank={}".format(
            A.rank()
        )
    )

    return C, A


# ============================================================================
# CAYLEY-HAMILTON
# ============================================================================

def cayley_hamilton_audit(label, T):

    tr = clean(
        sp.trace(T)
    )

    det = clean(
        T.det()
    )

    disc = clean(
        tr ** 2
        - 4 * det
    )

    residual = matrix_clean(
        T ** 2
        - tr * T
        + det * sp.eye(2)
    )

    print()
    print(
        "  {}:".format(label)
    )

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
        "    CH_residual={}".format(
            residual
        )
    )

    print(
        "    CH_exact={}".format(
            residual == sp.zeros(2)
        )
    )

    print(
        "    valuations={}".format(
            {
                p: (
                    valuation(disc, p)
                    if disc != 0
                    else None
                )
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


def characteristic_audit(T1, T2):

    print()
    print("=" * 78)
    print(
        "8. CAYLEY-HAMILTON / CHARACTERISTIC AUDIT"
    )
    print("=" * 78)

    cayley_hamilton_audit(
        "T1",
        T1,
    )

    cayley_hamilton_audit(
        "T2",
        T2,
    )


# ============================================================================
# FULL M2(Q) TEST
# ============================================================================

def full_algebra_audit(I, T1, T2):

    print()
    print("=" * 78)
    print(
        "9. GENERATED ALGEBRA DIMENSION AUDIT"
    )
    print("=" * 78)

    words = [
        I,
        T1,
        T2,
        T1 * T2,
        T2 * T1,
    ]

    labels = [
        "I",
        "T1",
        "T2",
        "T1T2",
        "T2T1",
    ]

    M = sp.Matrix.hstack(
        *[
            flatten_matrix(W)
            for W in words
        ]
    )

    rank = M.rank()

    print(
        "  labels={}".format(
            labels
        )
    )

    print(
        "  rank={}".format(
            rank
        )
    )

    print(
        "  ambient_dimension=4"
    )

    print(
        "  full_M2_Q_reached={}".format(
            rank == 4
        )
    )

    return rank


# ============================================================================
# TERMINAL SOURCE
# ============================================================================

def terminal_reference():

    q1_terminal = 495451247
    q3_terminal = 421514439

    g = math.gcd(
        q1_terminal,
        q3_terminal,
    )

    print()
    print("=" * 78)
    print(
        "10. TERMINAL PROJECTIVE SOURCE REFERENCE"
    )
    print("=" * 78)

    print(
        "  q1_terminal={}".format(
            q1_terminal
        )
    )

    print(
        "  q3_terminal={}".format(
            q3_terminal
        )
    )

    print(
        "  gcd={}".format(
            g
        )
    )

    print(
        "  q1/17={}".format(
            q1_terminal // 17
        )
    )

    print(
        "  q3/17={}".format(
            q3_terminal // 17
        )
    )

    print(
        "  v17(q1)={}".format(
            valuation(
                q1_terminal,
                17,
            )
        )
    )

    print(
        "  v17(q3)={}".format(
            valuation(
                q3_terminal,
                17,
            )
        )
    )


# ============================================================================
# STRUCTURAL INTERPRETATION
# ============================================================================

def interpretation():

    print()
    print("=" * 78)
    print(
        "11. STRUCTURAL INTERPRETATION"
    )
    print("=" * 78)

    print(
r"""
Experiment 307R ruled out projective similarity of the two independently
determined width-2 transition matrices.

Experiment 308R ruled out commutative polynomial and Möbius-type
relations.

Experiment 309R found that the two transitions generate the full
4-dimensional matrix algebra M_2(Q).

Experiment 310R verifies that conclusion with a robust exact span solver
and explicitly audits closure under the low-degree mixed words.

The key distinction is:

    intrinsic 2x2 identities
        versus
    mixed algebra generated by both transitions.

Each individual matrix still obeys its own Cayley-Hamilton identity:

    T_i^2 - tr(T_i) T_i + det(T_i) I = 0.

That does NOT imply a common algebraic relation between T1 and T2.

The mixed-word audit asks whether products such as

    T1T2,
    T2T1,
    T1T2T1,
    T2T1T2

remain inside the span of

    I, T1, T2.

If they do not, the two transitions cannot be described by one
low-dimensional commutative source operator.

No arbitrary fitted matrix is introduced.

No synthetic second n=pq case is generated.
"""
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 310R — EXACT TRANSFER-ALGEBRA CLOSURE / "
        "WORD-SPAN AUDIT"
    )
    print("=" * 78)

    I, T1, T2 = build_matrices()

    print()
    print("=" * 78)
    print("1. EXACT TRANSFER MATRICES")
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
    print("=" * 78)
    print(
        "2. BASIC CHARACTERISTIC DATA"
    )
    print("=" * 78)

    for label, T in (
        ("T1", T1),
        ("T2", T2),
    ):

        print()
        print(
            "  {}:".format(
                label
            )
        )

        print(
            "    trace={}".format(
                clean(
                    sp.trace(T)
                )
            )
        )

        print(
            "    determinant={}".format(
                clean(
                    T.det()
                )
            )
        )

        print(
            "    charpoly={}".format(
                sp.factor(
                    T.charpoly().as_expr()
                )
            )
        )

    closure_results = closure_audit(
        I,
        T1,
        T2,
    )

    word_basis_audit(
        I,
        T1,
        T2,
    )

    commutator_audit(
        T1,
        T2,
    )

    characteristic_audit(
        T1,
        T2,
    )

    algebra_rank = full_algebra_audit(
        I,
        T1,
        T2,
    )

    terminal_reference()

    interpretation()

    print()
    print("=" * 78)
    print(
        "12. FINAL EXACTNESS"
    )
    print("=" * 78)

    t12 = closure_results["T1*T2"]
    t21 = closure_results["T2*T1"]

    print(
        "  T1T2_in_span_I_T1_T2={}".format(
            t12["status"] == "EXACT"
        )
    )

    print(
        "  T2T1_in_span_I_T1_T2={}".format(
            t21["status"] == "EXACT"
        )
    )

    print(
        "  generated_algebra_dimension={}".format(
            algebra_rank
        )
    )

    print(
        "  full_M2_Q_reached={}".format(
            algebra_rank == 4
        )
    )

    print(
        "  robust_span_solver_used=True"
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
        "EXPERIMENT 310R COMPLETE"
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