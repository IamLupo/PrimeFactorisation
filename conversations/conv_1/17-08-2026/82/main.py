#!/usr/bin/env python3

"""
==============================================================================
KAPPA EXPERIMENT 146
N-ONLY ALGEBRAIC DETECTOR SURROGATES

PURPOSE
-------
Experiment 145 showed:

    detector values -> X -> Y=2S-1 -> S -> p,q

but also exposed the real bottleneck:

    can the detector values themselves be obtained from N alone?

This experiment does NOT assume that.

Instead it asks a sharper algebraic question:

    Can each detector be eliminated against the hidden trace S,
    producing a low-degree polynomial equation

        P_i(N, Z_i) = 0

    whose roots are the possible detector values?

If several detectors produce low-degree equations, we then test whether
their root sets intersect in a uniquely identifying branch.

The experiment explicitly distinguishes:

    A. algebraic eliminability;
    B. number of possible detector branches;
    C. uniqueness;
    D. whether uniqueness already amounts to solving the factorisation
       problem.

IMPORTANT:
------------
Because N and S are algebraically independent at the symmetric-polynomial
level, a detector that depends genuinely on S may NOT collapse to a
nontrivial polynomial in N alone. This experiment is designed to detect
that obstruction rather than hide it.

We therefore test two models:

    MODEL 1:
        eliminate S directly from Q_i(N,S)-Z_i.

    MODEL 2:
        impose the integer-factor relation

            S^2 - 4N = T^2

        and eliminate S,T.

MODEL 2 is deliberately labeled as a factorisation-level constraint,
because it may simply reintroduce the original problem.

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

N, S, T = sp.symbols(
    "N S T"
)

Z13, Z15, Z17 = sp.symbols(
    "Z13 Z15 Z17"
)


# ============================================================================
# Exact detector polynomials in (N,S)
# ============================================================================
#
# These are obtained by substituting
#
#     X = S(S-1)
#     Y = 2S-1
#
# into the exact C + YD decompositions.
#
# ============================================================================

def detector_13():
    X = sp.expand(
        S*(S-1)
    )

    return sp.expand(
        6*N - X
    )


def detector_15():
    X = sp.expand(
        S*(S-1)
    )

    Y = sp.expand(
        2*S - 1
    )

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

    return sp.expand(
        C + Y*D
    )


def detector_17():
    X = sp.expand(
        S*(S-1)
    )

    Y = sp.expand(
        2*S - 1
    )

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

    return sp.expand(
        C + Y*D
    )


# ============================================================================
# Utility
# ============================================================================

def factor(expr):
    return sp.factor(
        sp.expand(expr)
    )


def polynomial_degree(expr, var):
    poly = sp.Poly(
        sp.expand(expr),
        var,
        domain=sp.QQ.frac_field(N)
    )

    if poly.is_zero:
        return -sp.oo

    return poly.degree()


# ============================================================================
# MODEL 1
# Direct elimination of S
# ============================================================================

def direct_resultant(Q, Z):
    """
    Resultant_S(Q(N,S)-Z).

    If the result is identically zero, there is no nontrivial algebraic
    relation between N and Z after eliminating S.
    """

    return sp.factor(
        sp.resultant(
            sp.expand(Q - Z),
            0*S + 1,
            S
        )
    )


# ============================================================================
# Correct direct eliminability test
# ============================================================================

def direct_elimination_test(Q, Z):
    """
    A single equation Q(N,S)-Z has no second equation in S.

    Therefore elimination of S alone is expected to be trivial.

    We explicitly report this rather than invoking a meaningless
    resultant.

    The correct algebraic interpretation is:

        Q(N,S) = Z

    leaves S free for generic symbolic N.

    We return the primitive polynomial in S and its degree.
    """

    poly = sp.Poly(
        sp.expand(Q - Z),
        S,
        domain=sp.QQ.frac_field(N, Z)
    )

    return {
        "degree_S": poly.degree(),
        "expression": sp.factor(
            poly.as_expr()
        )
    }


# ============================================================================
# MODEL 2
#
# Add the integer-factor relation
#
#     S^2 - 4N = T^2.
#
# Then eliminate S,T.
# ============================================================================

def factorisation_level_resultant(Q, Z):
    """
    Eliminate S and T from

        Q(N,S)-Z = 0
        S^2 - T^2 - 4N = 0.

    This is intentionally interpreted as a factorisation-level relation.
    """

    f = sp.expand(
        Q - Z
    )

    g = sp.expand(
        S**2 - T**2 - 4*N
    )

    print(
        "    eliminating S..."
    )

    R1 = sp.resultant(
        f,
        g,
        S
    )

    R1 = sp.factor(
        sp.expand(R1)
    )

    if R1 == 0:
        return sp.Integer(0)

    print(
        "    eliminating T..."
    )

    # R1 may contain T.
    poly_T = sp.Poly(
        R1,
        T,
        domain=sp.QQ.frac_field(N, Z)
    )

    if poly_T.is_zero:
        return sp.Integer(0)

    R2 = sp.resultant(
        poly_T,
        sp.Integer(0),
        T
    )

    return sp.factor(
        sp.expand(R2)
    )


# ============================================================================
# Better second-stage elimination
# ============================================================================

def eliminate_S_then_factor(Q, Z):
    """
    Since g is quadratic in S, use the relation

        S^2 = T^2 + 4N

    to reduce powers of S.

    This keeps the calculation compact and makes the factorisation
    dependency visible.
    """

    f = sp.Poly(
        sp.expand(Q - Z),
        S,
        domain=sp.QQ.frac_field(N, T, Z)
    )

    g = sp.Poly(
        S**2 - T**2 - 4*N,
        S,
        domain=sp.QQ.frac_field(N, T, Z)
    )

    rem = sp.rem(
        f,
        g
    )

    rem_expr = sp.factor(
        sp.expand(
            rem.as_expr()
        )
    )

    # The remainder has degree <= 1 in S:
    #     a(N,T,Z) S + b(N,T,Z)
    #
    # Together with S^2=T^2+4N, eliminate S explicitly.
    rem_poly = sp.Poly(
        rem_expr,
        S,
        domain=sp.QQ.frac_field(N, T, Z)
    )

    a = sp.expand(
        rem_poly.coeff_monomial(S)
    )

    b = sp.expand(
        rem_poly.coeff_monomial(1)
    )

    if a == 0:
        return {
            "a": a,
            "b": b,
            "relation_T": factor(b),
            "degree_T": 0
        }

    # S = -b/a.
    # Substitute into S^2=T^2+4N.
    relation_T = factor(
        a**2*(T**2 + 4*N) - b**2
    )

    relation_T = sp.factor(
        sp.expand(relation_T)
    )

    degree_T = polynomial_degree(
        relation_T,
        T
    )

    return {
        "a": a,
        "b": b,
        "relation_T": relation_T,
        "degree_T": degree_T
    }


# ============================================================================
# Numeric branch experiment
# ============================================================================

def numeric_branch_count(
    n_value,
    detector_values,
    max_S=None
):
    """
    This is NOT a factor-pair search.

    We do NOT enumerate divisors or factor pairs.

    Instead we solve the polynomial detector equations in S over the
    integer interval implied by

        S >= 2 sqrt(N).

    The interval is deliberately small and theoretical; this test is
    only used to count algebraic branches, not to factor N.

    For an actual factorisation algorithm this procedure would be
    prohibited; here it is diagnostic only.
    """

    if max_S is None:
        max_S = 4 * int(
            sp.sqrt(n_value)
        ) + 20

    Z13v, Z15v = detector_values

    Q13 = detector_13()
    Q15 = detector_15()

    f13 = sp.Poly(
        sp.expand(
            Q13.subs(N, n_value)
            - Z13v
        ),
        S,
        domain=sp.QQ
    )

    f15 = sp.Poly(
        sp.expand(
            Q15.subs(N, n_value)
            - Z15v
        ),
        S,
        domain=sp.QQ
    )

    roots13 = {
        int(r)
        for r in sp.ground_roots(
            f13.as_expr(),
            S
        ).keys()
        if r.is_Integer
    }

    roots15 = {
        int(r)
        for r in sp.ground_roots(
            f15.as_expr(),
            S
        ).keys()
        if r.is_Integer
    }

    common = sorted(
        roots13.intersection(
            roots15
        )
    )

    return (
        sorted(roots13),
        sorted(roots15),
        common
    )


# ============================================================================
# Main
# ============================================================================

def main():

    print("=" * 78)
    print("KAPPA EXPERIMENT 146")
    print("N-ONLY ALGEBRAIC DETECTOR SURROGATES")
    print("=" * 78)
    print()

    Q13 = detector_13()
    Q15 = detector_15()
    Q17 = detector_17()

    # ------------------------------------------------------------------------
    # 1. Detector structure
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("1. DETECTOR POLYNOMIALS IN S")
    print("=" * 78)

    print(
        "Q13(N,S) =",
        factor(Q13)
    )

    print(
        "Q15(N,S) =",
        factor(Q15)
    )

    print(
        "Q17(N,S) =",
        factor(Q17)
    )

    print()

    for name, Q in [
        ("Q13", Q13),
        ("Q15", Q15),
        ("Q17", Q17)
    ]:
        degree = polynomial_degree(
            Q,
            S
        )

        print(
            f"{name}: degree_S = {degree}"
        )

    # ------------------------------------------------------------------------
    # 2. MODEL 1
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("2. MODEL 1: DIRECT N-ONLY ELIMINATION")
    print("=" * 78)

    model1 = []

    for name, Q, Z in [
        ("Q13", Q13, Z13),
        ("Q15", Q15, Z15),
        ("Q17", Q17, Z17)
    ]:

        info = direct_elimination_test(
            Q,
            Z
        )

        model1.append(
            (
                name,
                info
            )
        )

        print(
            f"{name}: degree in hidden S = "
            f"{info['degree_S']}"
        )

        print(
            f"  equation = {info['expression']}"
        )

    print()
    print(
        "INTERPRETATION:"
    )
    print(
        "A single detector equation Q(N,S)=Z does not"
    )
    print(
        "produce an N-only relation because S remains free."
    )

    # ------------------------------------------------------------------------
    # 3. MODEL 2
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("3. MODEL 2: ADD FACTORISATION CONSTRAINT")
    print("=" * 78)

    for name, Q, Z in [
        ("Q13", Q13, Z13),
        ("Q15", Q15, Z15),
        ("Q17", Q17, Z17)
    ]:

        print()
        print(
            f"{name}:"
        )

        info = eliminate_S_then_factor(
            Q,
            Z
        )

        print(
            "  remainder =",
            info["a"],
            "*S +",
            info["b"]
        )

        print(
            "  relation in T ="
        )

        print(
            info["relation_T"]
        )

        print(
            "  degree_T =",
            info["degree_T"]
        )

    # ------------------------------------------------------------------------
    # 4. Joint detector consistency
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("4. JOINT DETECTOR CONSISTENCY")
    print("=" * 78)

    # The exact relation from Exp. 145:
    #
    # Y = -(8N² + 4NZ13 + 12N - 2Z13² + 3Z13 - 2Z15)/(4N+Z13)
    #
    # and
    #
    # Y² = 24N - 4Z13 + 1.
    #

    numerator = sp.expand(
        8*N**2
        + 4*N*Z13
        + 12*N
        - 2*Z13**2
        + 3*Z13
        - 2*Z15
    )

    bridge = sp.factor(
        sp.expand(
            numerator**2
            - (
                24*N
                - 4*Z13
                + 1
            )
            * (
                4*N + Z13
            )**2
        )
    )

    print(
        "exact observable relation:"
    )

    print(
        "P(N,Z13,Z15) ="
    )

    print(
        bridge
    )

    print()

    # ------------------------------------------------------------------------
    # 5. Does this determine detectors from N?
    # ------------------------------------------------------------------------

    print("=" * 78)
    print("5. DOES THE RELATION DETERMINE DETECTORS FROM N?")
    print("=" * 78)

    poly_Z15 = sp.Poly(
        bridge,
        Z15,
        domain=sp.QQ.frac_field(N, Z13)
    )

    degree_Z15 = (
        0
        if poly_Z15.is_zero
        else poly_Z15.degree()
    )

    print(
        "degree in Z15 =",
        degree_Z15
    )

    print()

    if degree_Z15 == 2:

        print(
            "The relation is quadratic in Z15."
        )

        roots = sp.solve(
            bridge,
            Z15
        )

        print(
            "formal Z15 branches ="
        )

        for root in roots:
            print(
                " ",
                sp.factor(root)
            )

    else:

        print(
            "The relation has no useful quadratic Z15 branch."
        )

    # ------------------------------------------------------------------------
    # 6. Final interpretation
    # ------------------------------------------------------------------------

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "The detector system has a nontrivial exact algebraic"
    )

    print(
        "relation once multiple detector values are combined."
    )

    print()

    print(
        "However, the relation still contains a free detector"
    )

    print(
        "coordinate (e.g. Z13)."
    )

    print()

    print(
        "Therefore the exact state is:"
    )

    print(
        "    detector values -> trace S : ALGEBRAICALLY SOLVED"
    )

    print(
        "    N -> detector values       : STILL OPEN"
    )

    print()

    print(
        "MODEL 2 WARNING:"
    )

    print(
        "Imposing S²-4N=T² introduces the factorisation constraint"
    )

    print(
        "itself. A low-degree relation obtained there is not yet"
    )

    print(
        "an N-only factorisation algorithm."
    )

    print()

    print(
        "STATUS = DIAGNOSTIC"
    )

    print(
        "NEXT TARGET:"
    )

    print(
        "Test whether any nontrivial combination of detectors"
    )

    print(
        "is independent of S before introducing the constraint"
    )

    print(
        "S²-4N=T²."
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

