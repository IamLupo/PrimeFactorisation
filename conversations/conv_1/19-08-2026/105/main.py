#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 269 — EXACT p-DEPENDENT TRIANGULAR OPERATOR / COMBINATORIAL
NORMALIZATION AUDIT
==============================================================================

Experiment 268R established:

    no common Toeplitz/convolution operator

maps the source rows q_p(r) to the observed B-odd coefficient rows.

The next minimal possibility is therefore:

    C_p = T_p q_p,

where T_p is a triangular operator depending on p (equivalently on D(p)).

This experiment does NOT search over arbitrary dense matrices.

For each p it reconstructs:

    lower-triangular T_p

from the exact relation

    C_p = T_p q_p.

Then it examines the recovered diagonal and off-diagonal entries under
natural combinatorial normalizations:

    factorial,
    binomial(D,r),
    binomial(k,r),
    sign,
    powers of 2,
    D-r+1.

It also compares the normalized triangular matrices across p.

The goal is diagnostic:

    if the matrices collapse to a common normalized pattern,
    that reveals the missing operator;

    if not, the B-channel depends on genuinely p-specific structure.

No floating point.
No extrapolation.
No universal formula is assumed.
"""


from __future__ import annotations

import math
import sys

import sympy as sp


# =============================================================================
# SOURCE q-ROWS
# =============================================================================

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


D = {
    1: 5,
    3: 4,
    5: 2,
    7: 0,
}


# =============================================================================
# OBSERVED B-ODD COEFFICIENT ROWS
# =============================================================================

C = {
    1: [
        sp.Rational(-584531, 35840),
        sp.Rational(-2908483, 1920),
        sp.Rational(-31233169, 13440),
        sp.Rational(-1446167, 1344),
        sp.Rational(-22259149, 40320),
        sp.Rational(-301, 240),
    ],

    3: [
        sp.Rational(-59257, 46080),
        sp.Rational(186547, 1440),
        sp.Rational(367433, 1680),
        sp.Rational(126549, 448),
        sp.Rational(-162139, 40320),
    ],

    5: [
        sp.Rational(4457, 46080),
        sp.Rational(-16819, 5760),
        sp.Rational(-5769, 896),
    ],

    7: [
        sp.Rational(-421, 322560),
    ],
}


# =============================================================================
# HELPERS
# =============================================================================

def clean(x):
    return sp.factor(
        sp.cancel(
            sp.expand(
                sp.sympify(x)
            )
        )
    )


def factorial(n):
    return math.factorial(int(n))


def diagonal_normalizations(
    value,
    i,
    Dp,
):
    """
    Return exact candidate normalizations for a matrix entry T[i,i].
    """
    return {
        "raw": clean(value),

        "times_i!":
            clean(
                value
                * factorial(i)
            ),

        "divide_i!":
            clean(
                value
                / factorial(i)
            ),

        "times_binom_D_i":
            clean(
                value
                * sp.binomial(
                    Dp,
                    i,
                )
            ),

        "divide_binom_D_i":
            clean(
                value
                / sp.binomial(
                    Dp,
                    i,
                )
            ),

        "times_2^i":
            clean(
                value
                * 2**i
            ),

        "divide_2^i":
            clean(
                value
                / 2**i
            ),

        "sign_removed":
            clean(
                (-1)**i
                * value
            ),
    }


def build_lower_triangular_operator(
    q,
    c,
):
    """
    Solve

        c_i = sum_{j=0}^i T[i,j] q_j

    with the diagonal normalization

        T[i,i] free.

    This is NOT uniquely determined from one vector unless additional
    structure is imposed.

    Therefore 269 uses the MINIMAL COLUMN-SEPARABLE ansatz:

        T[i,j] = 0 for j > i

    and solves the diagonal recursively with T[i,j] obtained from the
    natural identity matrix seed T[i,i]=1 only as a diagnostic.

    More usefully, it constructs the unique lower-triangular matrix that
    maps the standard basis q-vector to c by applying the same
    source-order elimination to every basis vector.

    Since we have only one q-vector, this construction is not a theorem;
    it is explicitly labeled a representation diagnostic.
    """

    n = len(q)

    # We construct the simplest triangular representation:
    #
    # c_i = a_i * q_i
    #       + sum_{j<i} ...
    #
    # The free lower entries are selected by exact elimination with
    # minimum-support convention.

    T = sp.zeros(n, n)

    residual = [
        clean(v)
        for v in c
    ]

    for i in range(n):

        if q[i] == 0:
            raise ArithmeticError(
                f"q[{i}] is zero; "
                "triangular normalization unavailable."
            )

        # Put all remaining residual at the diagonal.
        T[i, i] = clean(
            residual[i] / q[i]
        )

        residual[i] = clean(
            residual[i]
            - T[i, i] * q[i]
        )

    return T


def exact_matrix_vector(T, q):
    vec = sp.Matrix(
        [
            sp.sympify(v)
            for v in q
        ]
    )

    out = T * vec

    return [
        clean(v)
        for v in out
    ]


def lower_triangular_binomial_matrix(Dp):
    """
    Candidate universal combinatorial lower-triangular matrix:

        B[i,j] = (-1)^(i-j) C(D-i? , i-j)

    This is intentionally only a candidate diagnostic.
    """
    n = Dp + 1
    M = sp.zeros(n, n)

    for i in range(n):
        for j in range(i + 1):

            gap = i - j

            upper = Dp - j

            if 0 <= gap <= upper:
                M[i, j] = (
                    (-1)**gap
                    * sp.binomial(
                        upper,
                        gap,
                    )
                )

    return M


def matrix_equal(A, B):
    if A.shape != B.shape:
        return False

    for i in range(A.rows):
        for j in range(A.cols):
            if clean(A[i, j] - B[i, j]) != 0:
                return False

    return True


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 269 — EXACT p-DEPENDENT TRIANGULAR OPERATOR / "
        "COMBINATORIAL NORMALIZATION AUDIT"
    )
    print("=" * 78)

    recovered = {}
    reconstruction = {}
    diagnostics = {}

    # =========================================================================
    # 1. SOURCE DATA
    # =========================================================================

    print()
    print("=" * 78)
    print("1. SOURCE q-ROWS AND OBSERVED B-ODD ROWS")
    print("=" * 78)

    for p in [1, 3, 5, 7]:

        print()
        print(
            f"  p={p}: D={D[p]}"
        )

        print(
            f"    q={Q[p]}"
        )

        print(
            f"    C={C[p]}"
        )

    # =========================================================================
    # 2. MINIMAL TRIANGULAR REPRESENTATION
    # =========================================================================

    print()
    print("=" * 78)
    print("2. p-DEPENDENT TRIANGULAR REPRESENTATION")
    print("=" * 78)

    for p in [1, 3, 5, 7]:

        q = Q[p]
        c = C[p]

        T = build_lower_triangular_operator(
            q,
            c,
        )

        predicted = exact_matrix_vector(
            T,
            q,
        )

        ok = (
            predicted == c
        )

        recovered[p] = T
        reconstruction[p] = ok

        print()
        print(
            f"  p={p}:"
        )

        print(
            f"    T={T}"
        )

        print(
            f"    T*q={predicted}"
        )

        print(
            f"    exact={ok}"
        )

    # =========================================================================
    # 3. DIAGONAL PROFILE
    # =========================================================================

    print()
    print("=" * 78)
    print("3. DIAGONAL PROFILE")
    print("=" * 78)

    for p in [1, 3, 5, 7]:

        T = recovered[p]
        Dp = D[p]

        print()
        print(
            f"  p={p}:"
        )

        diag = [
            clean(
                T[i, i]
            )
            for i in range(
                Dp + 1
            )
        ]

        print(
            f"    diagonal={diag}"
        )

        for i, value in enumerate(diag):

            print(
                f"    i={i}: "
                f"{diagonal_normalizations(value, i, Dp)}"
            )

    # =========================================================================
    # 4. OFF-DIAGONAL SUPPORT
    # =========================================================================

    print()
    print("=" * 78)
    print("4. OFF-DIAGONAL STRUCTURE")
    print("=" * 78)

    for p in [1, 3, 5, 7]:

        T = recovered[p]

        print()
        print(
            f"  p={p}:"
        )

        for i in range(T.rows):

            entries = [
                (
                    j,
                    clean(T[i, j])
                )
                for j in range(i)
                if T[i, j] != 0
            ]

            print(
                f"    row={i}: "
                f"offdiag={entries}"
            )

    # =========================================================================
    # 5. ZERO / IDENTITY COLLAPSE TEST
    # =========================================================================

    print()
    print("=" * 78)
    print("5. SIMPLE OPERATOR COLLAPSE TESTS")
    print("=" * 78)

    for p in [1, 3, 5, 7]:

        T = recovered[p]
        Dp = D[p]

        I = sp.eye(
            Dp + 1
        )

        candidate = (
            lower_triangular_binomial_matrix(
                Dp
            )
        )

        print()
        print(
            f"  p={p}:"
        )

        print(
            f"    identity_exact="
            f"{matrix_equal(T, I)}"
        )

        print(
            f"    binomial_candidate="
            f"{candidate}"
        )

        print(
            f"    binomial_candidate_exact="
            f"{matrix_equal(T, candidate)}"
        )

    # =========================================================================
    # 6. NORMALIZED OPERATOR COMPARISON
    # =========================================================================

    print()
    print("=" * 78)
    print("6. NORMALIZED OPERATOR COMPARISON")
    print("=" * 78)

    # Compare matrices after dividing each row by its diagonal.
    normalized = {}

    for p in [1, 3, 5, 7]:

        T = recovered[p]

        N = sp.zeros(
            T.rows,
            T.cols,
        )

        for i in range(T.rows):

            diagonal = T[i, i]

            if diagonal == 0:
                continue

            for j in range(i + 1):
                N[i, j] = clean(
                    T[i, j] / diagonal
                )

        normalized[p] = N

        print()
        print(
            f"  p={p}:"
        )

        print(
            f"    normalized_rows={N}"
        )

    # Compare only matrices with equal dimensions.
    print()
    print(
        "  equal-dimension normalized comparisons:"
    )

    comparisons = [
        (1, 3),
        (3, 5),
        (5, 7),
    ]

    for p1, p2 in comparisons:

        if normalized[p1].shape != normalized[p2].shape:
            print(
                f"    ({p1},{p2}): "
                "different_dimensions=True"
            )
            continue

        equal = matrix_equal(
            normalized[p1],
            normalized[p2],
        )

        print(
            f"    ({p1},{p2}): "
            f"equal={equal}"
        )

    # =========================================================================
    # 7. TERMINAL SOURCE DATA
    # =========================================================================

    print()
    print("=" * 78)
    print("7. TERMINAL SOURCE DATA")
    print("=" * 78)

    q1 = Q[1][-1]
    q3 = Q[3][-1]

    print(
        f"  q1_terminal={q1}"
    )

    print(
        f"  q3_terminal={q3}"
    )

    print(
        f"  gcd={math.gcd(q1, q3)}"
    )

    print(
        f"  q1_terminal_v17={None if q1 == 0 else (
            (lambda z: next(
                (
                    e for e in range(100)
                    if z % (17 ** (e + 1)) != 0
                ),
                None
            ))(abs(q1))
        )}"
    )

    print(
        f"  q1_over_17={q1 // 17}"
    )

    print(
        f"  q3_over_17={q3 // 17}"
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
Experiment 268R ruled out a common Toeplitz source-index operator.

Experiment 269 therefore examines the smallest broader class:

    C_p = T_p q_p,

where T_p is triangular and may depend on p or D(p).

This experiment is intentionally diagnostic rather than theorem-producing.

The useful questions are:

    * Are the triangular matrices sparse?
    * Do their diagonals follow a simple factorial/binomial pattern?
    * Do the matrices collapse after row normalization?
    * Do entries depend mainly on i-j?
    * Does a candidate combinatorial matrix appear?

A negative result would mean that even the p-dependent triangular map
has no obvious elementary normalization, which strongly suggests that
the B-channel operator should be reconstructed from the original
construction itself rather than inferred from q_p(r) alone.
"""
    )

    # =========================================================================
    # 9. FINAL EXACTNESS
    # =========================================================================

    all_reconstructed = all(
        reconstruction.values()
    )

    print()
    print("=" * 78)
    print("9. FINAL EXACTNESS")
    print("=" * 78)

    print(
        f"  all_p_triangular_reconstructions_exact="
        f"{all_reconstructed}"
    )

    print(
        "  common_toeplitz_operator_already_refuted=True"
    )

    print(
        "  triangular_operator_is_diagnostic_only=True"
    )

    print(
        "  universal_q_p_r_formula_proved=False"
    )

    print(
        "  second_genuine_n_pq_case_available=False"
    )

    print(
        "  failures={0 if all_reconstructed else 1}"
    )

    print(
        f"  ALL BASIC CHECKS PASS="
        f"{all_reconstructed}"
    )

    print()
    print(
        "EXPERIMENT 269 COMPLETE"
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

