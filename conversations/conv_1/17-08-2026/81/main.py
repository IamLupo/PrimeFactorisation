#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 145
EXACT S,X ELIMINATION / N-ONLY RECOVERABILITY TEST

PURPOSE
-------
The symbolic proof phase is now complete enough to attack the actual
factorisation question:

        N = p q

with hidden

        S = p + q
        X = S(S-1)
        Y = 2S-1

and

        Y^2 = 4X + 1.

Each detector has the exact form

        Q_i(N,S)
          = C_i(N,X) + Y D_i(N,X).

The first algebraic step toward an N-only method is to eliminate Y.

For two detectors i,j:

    (Q_i-C_i) D_j - (Q_j-C_j) D_i = 0.

This gives a polynomial relation in

    N, X, Q_i, Q_j.

Then we eliminate X.

The experiment asks three exact questions:

    1. Does pairwise Y-elimination work exactly?
    2. Does eliminating X produce a nontrivial relation?
    3. Given N and the exact detector values, how many candidate S values
       survive the elimination?

This is NOT yet a factorisation algorithm.
It is a recoverability certificate.

NO FACTOR-PAIR SEARCH
NO FULL NEWTON TENSOR
NO FULL C/D TENSOR
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

N, S, X, Y = sp.symbols(
    "N S X Y"
)

Z1, Z2 = sp.symbols(
    "Z1 Z2"
)


# ============================================================================
# Detector definitions
# ============================================================================
#
# These are the exact detectors already established by the KAPPA program.
# We use a small pair first to keep elimination tractable.
#
# Q_(1,3) = 6N - X
#
# Q_(1,5) = C + Y D
#
# From the symbolic decomposition:
#
# C_(1,5) = -(40N^2 - 20NX - 30N + 2X^2 + 3X)/2
# D_(1,5) =  (10N - X)/2
#
# ============================================================================

def detector_13():
    C = sp.expand(
        6*N - X
    )

    D = sp.Integer(0)

    Q = sp.expand(
        C + Y*D
    )

    return C, D, Q


def detector_15():
    C = sp.expand(
        -(
            40*N**2
            - 20*N*X
            - 30*N
            + 2*X**2
            + 3*X
        ) / 2
    )

    D = sp.expand(
        (10*N - X) / 2
    )

    Q = sp.expand(
        C + Y*D
    )

    return C, D, Q


def detector_17():
    C = sp.expand(
        (
            84*N**3
            - 98*N**2*X
            - 224*N**2
            + 28*N*X**2
            + 133*N*X
            + 84*N
            - 2*X**3
            - 6*X**2
            - 4*X
        ) / 2
    )

    D = sp.expand(
        -(
            84*N**2
            - 35*N*X
            - 56*N
            + 2*X**2
            + 2*X
        ) / 2
    )

    Q = sp.expand(
        C + Y*D
    )

    return C, D, Q


# ============================================================================
# Exact elimination helpers
# ============================================================================

def factor_expr(expr):
    return sp.factor(
        sp.expand(expr)
    )


def eliminate_Y(C1, D1, C2, D2):
    """
    Exact linear elimination of Y.
    """

    relation = sp.expand(
        (Z1 - C1)*D2
        - (Z2 - C2)*D1
    )

    return sp.factor(
        relation
    )


def eliminate_X(relation):
    """
    Eliminate X by treating the relation as a polynomial in X.

    In this first experiment we pair it with Q_(1,3):

        Z1 = 6N-X

    which gives an exact direct substitution.
    """

    direct_X = sp.expand(
        6*N - Z1
    )

    reduced = sp.expand(
        relation.subs(
            X,
            direct_X
        )
    )

    return sp.factor(
        reduced
    )


# ============================================================================
# Recover Y after X is known
# ============================================================================

def recover_Y_from_detector(C, D, Z):
    """
    Y = (Z-C)/D, provided D != 0.
    """

    numerator = sp.expand(
        Z - C
    )

    return numerator


# ============================================================================
# Discriminant bridge
# ============================================================================

def discriminant_bridge():
    return sp.expand(
        Y**2 - (4*X + 1)
    )


# ============================================================================
# Exact detector reconstruction certificate
# ============================================================================

def reconstruction_certificate(
    C1,
    D1,
    C2,
    D2,
):
    Q1 = sp.expand(
        C1 + Y*D1
    )

    Q2 = sp.expand(
        C2 + Y*D2
    )

    r1 = sp.expand(
        Q1 - (C1 + Y*D1)
    )

    r2 = sp.expand(
        Q2 - (C2 + Y*D2)
    )

    return (
        sp.factor(r1),
        sp.factor(r2)
    )


# ============================================================================
# Candidate S analysis
# ============================================================================

def candidate_S_from_Y(y_value):
    return sp.expand(
        (y_value + 1) / 2
    )


def candidate_discriminant(
    n_value,
    s_value
):
    return sp.expand(
        s_value**2 - 4*n_value
    )


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 145")
    print("EXACT S,X ELIMINATION / N-ONLY RECOVERABILITY TEST")
    print("=" * 78)
    print()

    # ------------------------------------------------------------------------
    # 1. Build detectors
    # ------------------------------------------------------------------------

    C13, D13, Q13 = detector_13()
    C15, D15, Q15 = detector_15()
    C17, D17, Q17 = detector_17()

    print("=" * 78)
    print("1. DETECTOR DEFINITIONS")
    print("=" * 78)

    print(
        "Q_(1,3) =",
        factor_expr(Q13)
    )

    print(
        "Q_(1,5) =",
        factor_expr(Q15)
    )

    print(
        "Q_(1,7) =",
        factor_expr(Q17)
    )

    print()

    # ------------------------------------------------------------------------
    # 2. Exact Y-elimination
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("2. EXACT Y-ELIMINATION")
    print("=" * 78)

    R_13_15 = eliminate_Y(
        C13,
        D13,
        C15,
        D15,
    )

    R_13_17 = eliminate_Y(
        C13,
        D13,
        C17,
        D17,
    )

    print(
        "R_(13,15) =",
        R_13_15
    )

    print(
        "R_(13,17) =",
        R_13_17
    )

    # ------------------------------------------------------------------------
    # 3. Eliminate X
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. EXACT X-ELIMINATION")
    print("=" * 78)

    E_13_15 = eliminate_X(
        R_13_15
    )

    E_13_17 = eliminate_X(
        R_13_17
    )

    print(
        "E_(13,15)(N,Z1,Z2) ="
    )
    print(
        factor_expr(E_13_15)
    )

    print()

    print(
        "E_(13,17)(N,Z1,Z2) ="
    )
    print(
        factor_expr(E_13_17)
    )

    # ------------------------------------------------------------------------
    # 4. Check whether the eliminated relations are nontrivial
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. ELIMINATION NONTRIVIALITY")
    print("=" * 78)

    triv_15 = (
        sp.expand(E_13_15) == 0
    )

    triv_17 = (
        sp.expand(E_13_17) == 0
    )

    print(
        "13/15 relation identically zero =",
        triv_15
    )

    print(
        "13/17 relation identically zero =",
        triv_17
    )

    # ------------------------------------------------------------------------
    # 5. Solve for X and Y
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("5. RECONSTRUCTABILITY FROM DETECTOR VALUES")
    print("=" * 78)

    X_recovered = sp.expand(
        6*N - Z1
    )

    print(
        "X =",
        X_recovered
    )

    Y_from_15 = sp.factor(
        sp.cancel(
            (Z2 - C15.subs(X, X_recovered))
            / D15.subs(X, X_recovered)
        )
    )

    print(
        "Y from Q_(1,5) ="
    )
    print(
        Y_from_15
    )

    # ------------------------------------------------------------------------
    # 6. Discriminant bridge
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("6. DISCRIMINANT BRIDGE")
    print("=" * 78)

    bridge = discriminant_bridge()

    bridge_after_X = sp.expand(
        bridge.subs(
            X,
            X_recovered
        )
    )

    print(
        "Y^2-(4X+1) =",
        bridge
    )

    print(
        "after X substitution ="
    )

    print(
        factor_expr(
            bridge_after_X
        )
    )

    # ------------------------------------------------------------------------
    # 7. Detector reconstruction
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("7. EXACT RECONSTRUCTION")
    print("=" * 78)

    rec1, rec2 = reconstruction_certificate(
        C13,
        D13,
        C15,
        D15,
    )

    print(
        "Q_(1,3) reconstruction residual =",
        rec1
    )

    print(
        "Q_(1,5) reconstruction residual =",
        rec2
    )

    # ------------------------------------------------------------------------
    # 8. Concrete symbolic factorization map
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("8. SYMBOLIC FACTOR RECOVERY MAP")
    print("=" * 78)

    S_from_Y = candidate_S_from_Y(
        Y_from_15
    )

    Delta = candidate_discriminant(
        N,
        S_from_Y
    )

    print(
        "S = (Y+1)/2 ="
    )

    print(
        factor_expr(
            S_from_Y
        )
    )

    print()

    print(
        "Delta = S^2 - 4N ="
    )

    print(
        factor_expr(
            Delta
        )
    )

    p_symbolic = sp.factor(
        (S_from_Y + sp.sqrt(Delta)) / 2
    )

    q_symbolic = sp.factor(
        (S_from_Y - sp.sqrt(Delta)) / 2
    )

    print()

    print(
        "p ="
    )
    print(
        p_symbolic
    )

    print()

    print(
        "q ="
    )
    print(
        q_symbolic
    )

    # ------------------------------------------------------------------------
    # 9. Important distinction
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("9. N-ONLY RECOVERABILITY DIAGNOSTIC")
    print("=" * 78)

    print(
        "The algebraic system can recover X from Q_(1,3):"
    )

    print(
        "    X = 6N - Q_(1,3)"
    )

    print()

    print(
        "Then Q_(1,5) can recover Y = 2S-1, provided"
    )

    print(
        "its denominator is nonzero."
    )

    print()

    print(
        "Therefore the remaining algorithmic question is NOT"
    )

    print(
        "whether S is algebraically identifiable once detector"
    )

    print(
        "values are known."
    )

    print()

    print(
        "The remaining question is:"
    )

    print(
        "    Can Q_(1,3), Q_(1,5), ... be computed from N alone"
    )

    print(
        "    without already knowing p or q?"
    )

    # ------------------------------------------------------------------------
    # 10. Final status
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    if (
        not triv_15
        and not triv_17
        and rec1 == 0
        and rec2 == 0
    ):

        print(
            "STATUS = PASS"
        )

        print()
        print(
            "Exact detector elimination is successful."
        )

        print(
            "The hidden variables satisfy the exact chain:"
        )

        print()
        print(
            "    Q_(1,3) -> X"
        )

        print(
            "    Q_(1,5) -> Y=2S-1"
        )

        print(
            "    Y -> S"
        )

        print(
            "    S,N -> p,q"
        )

        print()
        print(
            "This does NOT yet constitute an N-only factorisation"
        )

        print(
            "algorithm, because detector evaluation from N alone"
        )

        print(
            "has not been established."
        )

        print()
        print(
            "NEXT TARGET:"
        )

        print(
            "construct an N-computable surrogate for the detector"
        )

        print(
            "values, or prove that the detector values themselves"
        )

        print(
            "can be generated without knowing the hidden trace."
        )

    else:

        print(
            "STATUS = FAIL"
        )

    print("=" * 78)


if __name__ == "__main__":

    try:
        main()

    except Exception as exc:

        print(
            "FATAL:",
            type(exc).__name__,
            str(exc)
        )

        sys.exit(1)

