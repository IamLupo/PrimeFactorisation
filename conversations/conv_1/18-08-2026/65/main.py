#!/usr/bin/env python3

"""
==============================================================================
EXPERIMENT 68 — DISCRIMINANT / HOMOGENEOUS-LAYER REDUCTION
==============================================================================

Exact symbolic experiment for the k=9, ell=16 kernel.

Main questions
--------------
1. Convert F(p,q) exactly to G(N,X), where

       N = p*q
       X = p+q+1

2. Write

       G(N,S+1) = A(N,Delta) + S*B(N,Delta)

   with

       S     = p+q
       Delta = S^2 - 4N = (p-q)^2.

3. Reduce every homogeneous layer independently.

4. Test whether the kernel is Delta-only, S-dependent,
   or has useful layer structure.

5. Verify everything on a fresh random semiprime.

No fitting.
No interpolation.
No regression.
No floating-point arithmetic.
==============================================================================

"""

import random
import sympy as sp


# ============================================================================
# SYMBOLS
# ============================================================================

p, q = sp.symbols("p q")
N, X = sp.symbols("N X")
S, Delta = sp.symbols("S Delta")


# ============================================================================
# EXACT k=9, ell=16 KERNEL
# ============================================================================

def exact_F_9_16(pv, qv):
    return sp.expand(
        -9*pv**16*qv**8
        -36*pv**16*qv**7
        -84*pv**16*qv**6
        -126*pv**16*qv**5
        -126*pv**16*qv**4
        -84*pv**16*qv**3
        -36*pv**16*qv**2
        -9*pv**16*qv
        -pv**16

        +16*pv**15*qv**9
        +120*pv**14*qv**9
        +560*pv**13*qv**9
        +1820*pv**12*qv**9
        +4368*pv**11*qv**9
        +8008*pv**10*qv**9

        +16*pv**9*qv**15
        +120*pv**9*qv**14
        +560*pv**9*qv**13
        +1820*pv**9*qv**12
        +4368*pv**9*qv**11
        +8008*pv**9*qv**10
        +22880*pv**9*qv**9
        +12870*pv**9*qv**8
        +11440*pv**9*qv**7
        +8008*pv**9*qv**6
        +4368*pv**9*qv**5
        +1820*pv**9*qv**4
        +560*pv**9*qv**3
        +120*pv**9*qv**2
        +16*pv**9*qv
        +pv**9

        -9*pv**8*qv**16
        +12870*pv**8*qv**9

        -36*pv**7*qv**16
        +11440*pv**7*qv**9

        -84*pv**6*qv**16
        +8008*pv**6*qv**9

        -126*pv**5*qv**16
        +4368*pv**5*qv**9

        -126*pv**4*qv**16
        +1820*pv**4*qv**9

        -84*pv**3*qv**16
        +560*pv**3*qv**9

        -36*pv**2*qv**16
        +120*pv**2*qv**9

        -9*pv*qv**16
        +16*pv*qv**9

        -qv**16
        +qv**9
    )


# ============================================================================
# SYMMETRIC REDUCTION
# ============================================================================

def construct_G():
    """
    Exact conversion

        F(p,q) -> G(N,X)

    using elementary symmetric variables.
    """

    F = exact_F_9_16(p, q)

    symmetric_part, remainder, mapping = sp.symmetrize(
        F,
        [p, q],
        formal=True,
    )

    if sp.expand(remainder) != 0:
        raise RuntimeError(
            f"Symmetrization produced nonzero remainder: {remainder}"
        )

    mapping_dict = dict(mapping)

    sum_symbol = None
    product_symbol = None

    for generated, original in mapping_dict.items():
        original = sp.expand(original)

        if original == p + q:
            sum_symbol = generated

        elif original == p*q:
            product_symbol = generated

    if sum_symbol is None or product_symbol is None:
        raise RuntimeError(
            f"Could not identify elementary symmetric variables: "
            f"{mapping_dict}"
        )

    G_NS = sp.expand(
        symmetric_part.subs({
            sum_symbol: S,
            product_symbol: N,
        })
    )

    G_NX = sp.expand(
        G_NS.subs(S, X - 1)
    )

    return G_NX


# ============================================================================
# HOMOGENEOUS LAYERS
# ============================================================================

def homogeneous_layers(poly, variables):
    """
    Split a polynomial into ordinary total-degree layers.
    """

    P = sp.Poly(
        sp.expand(poly),
        *variables,
    )

    layers = {}

    for exponents, coeff in P.terms():

        total_degree = sum(exponents)

        monomial = coeff

        for variable, exponent in zip(
            variables,
            exponents,
        ):
            monomial *= variable**exponent

        layers[total_degree] = sp.expand(
            layers.get(
                total_degree,
                sp.Integer(0),
            ) + monomial
        )

    return {
        degree: sp.expand(expr)
        for degree, expr in sorted(
            layers.items(),
            reverse=True,
        )
    }


# ============================================================================
# EXACT S^r REDUCTION
# ============================================================================

def reduce_S_power(power):
    """
    Under

        S^2 = 4N + Delta,

    write

        S^power = A(N,Delta) + S*B(N,Delta).
    """

    power = int(power)

    if power < 0:
        raise ValueError("Negative S-power.")

    base = sp.expand(4*N + Delta)

    if power % 2 == 0:
        return (
            sp.expand(base**(power // 2)),
            sp.Integer(0),
        )

    return (
        sp.Integer(0),
        sp.expand(base**((power - 1) // 2)),
    )


# ============================================================================
# REDUCE P(N,S) -> A(N,Delta) + S B(N,Delta)
# ============================================================================

def reduce_to_N_Delta(poly_NS):

    P = sp.Poly(
        sp.expand(poly_NS),
        N,
        S,
    )

    A = sp.Integer(0)
    B = sp.Integer(0)

    for (n_exp, s_exp), coeff in P.terms():

        A_part, B_part = reduce_S_power(
            s_exp
        )

        A += (
            coeff
            * N**n_exp
            * A_part
        )

        B += (
            coeff
            * N**n_exp
            * B_part
        )

    return (
        sp.expand(A),
        sp.expand(B),
    )


# ============================================================================
# RECONSTRUCTION
# ============================================================================

def reconstruct_from_N_Delta(A, B):

    Delta_sub = S**2 - 4*N

    return sp.expand(
        A.subs(Delta, Delta_sub)
        +
        S * B.subs(Delta, Delta_sub)
    )


# ============================================================================
# POLYNOMIAL INFORMATION
# ============================================================================

def total_degree(poly, *variables):

    poly = sp.expand(poly)

    if poly == 0:
        return sp.S.NegativeInfinity

    return sp.Poly(
        poly,
        *variables,
    ).total_degree()


def degree_string(poly, *variables):

    degree = total_degree(
        poly,
        *variables,
    )

    if degree == sp.S.NegativeInfinity:
        return "-oo"

    return str(int(degree))


def monomial_count(poly, *variables):

    poly = sp.expand(poly)

    if poly == 0:
        return 0

    return len(
        sp.Poly(
            poly,
            *variables,
        ).terms()
    )


# ============================================================================
# ROBUST EXACT INTEGER HELPERS
# ============================================================================

def lcm_all(values):
    """
    Robust LCM.

    Important:
    SymPy's sp.ilcm(*args) raises TypeError when exactly
    one argument is supplied in the installed environment.

    This implementation never makes that call.
    """

    values = list(values)

    if not values:
        return sp.Integer(1)

    result = sp.Integer(1)

    for value in values:
        value = sp.Integer(abs(int(value)))

        if value == 0:
            continue

        result = sp.ilcm(
            result,
            value,
        )

    return sp.Integer(result)


def gcd_all(values):

    values = [
        sp.Integer(abs(int(v)))
        for v in values
        if int(v) != 0
    ]

    if not values:
        return sp.Integer(0)

    result = values[0]

    for value in values[1:]:
        result = sp.igcd(
            result,
            value,
        )

    return sp.Integer(result)


def primitive_integer_coefficients(
    poly,
    *variables,
):
    """
    Return primitive integer coefficients
    of a rational polynomial.

    Handles:
      zero polynomial
      constants
      one-coefficient polynomials
      arbitrary multivariate polynomials
    """

    poly = sp.expand(poly)

    if poly == 0:
        return [0]

    P = sp.Poly(
        poly,
        *variables,
    )

    coeffs = [
        sp.Rational(c)
        for c in P.coeffs()
    ]

    denominators = [
        sp.Integer(c.q)
        for c in coeffs
    ]

    denominator_lcm = lcm_all(
        denominators
    )

    integer_coeffs = [
        int(c * denominator_lcm)
        for c in coeffs
    ]

    content = gcd_all(
        integer_coeffs
    )

    if content == 0:
        return [0]

    primitive = [
        value // int(content)
        for value in integer_coeffs
    ]

    # Canonical sign convention:
    # first nonzero coefficient positive.
    for value in primitive:
        if value != 0:

            if value < 0:
                primitive = [
                    -v
                    for v in primitive
                ]

            break

    return primitive


# ============================================================================
# RANDOM PRIME PAIR
# ============================================================================

def random_prime_pair():

    rng = random.SystemRandom()

    while True:

        p_value = int(
            sp.randprime(
                500_000,
                1_000_000,
            )
        )

        q_value = int(
            sp.randprime(
                500_000,
                1_000_000,
            )
        )

        if p_value != q_value:
            return (
                p_value,
                q_value,
            )


# ============================================================================
# EXACT INTEGER EVALUATION
# ============================================================================

def exact_eval(
    poly,
    substitutions,
):

    value = sp.expand(
        poly.subs(substitutions)
    )

    if value.has(sp.Float):
        raise RuntimeError(
            f"Floating-point value detected: {value}"
        )

    return sp.Integer(value)


# ============================================================================
# MAIN EXPERIMENT
# ============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 68 — "
        "DISCRIMINANT / HOMOGENEOUS-LAYER REDUCTION"
    )
    print("=" * 78)

    # ----------------------------------------------------------------------
    # 0. CONSTRUCT EXACT G
    # ----------------------------------------------------------------------

    print("\n0. EXACT SYMBOLIC SETUP")

    F_symbolic = exact_F_9_16(
        p,
        q,
    )

    G = construct_G()

    symbolic_G_check = (
        sp.expand(
            F_symbolic
            -
            G.subs({
                N: p*q,
                X: p + q + 1,
            })
        )
        == 0
    )

    print(
        "  kernel = exact k=9, ell=16"
    )

    print(
        f"  F(p,q) == G(pq,p+q+1): "
        f"{symbolic_G_check}"
    )

    if not symbolic_G_check:
        raise RuntimeError(
            "Symbolic F/G identity failed."
        )

    # ----------------------------------------------------------------------
    # 1. HOMOGENEOUS LAYERS
    # ----------------------------------------------------------------------

    print("\n1. HOMOGENEOUS LAYERS")

    layers = homogeneous_layers(
        G,
        (N, X),
    )

    degrees = sorted(
        layers.keys(),
        reverse=True,
    )

    print(
        f"  highest total degree = {degrees[0]}"
    )

    print(
        f"  lowest total degree = {degrees[-1]}"
    )

    print(
        f"  number of layers = {len(degrees)}"
    )

    for degree in degrees:
        print(
            f"  degree={degree}: PRESENT"
        )

    # ----------------------------------------------------------------------
    # 2. FULL DISCRIMINANT REDUCTION
    # ----------------------------------------------------------------------

    print(
        "\n2. FULL DISCRIMINANT REDUCTION"
    )

    G_NS = sp.expand(
        G.subs(X, S + 1)
    )

    A_full, B_full = reduce_to_N_Delta(
        G_NS
    )

    print(
        "  canonical form:"
    )

    print(
        "    G(N,S+1) = "
        "A(N,Delta) + S*B(N,Delta)"
    )

    print(
        f"\n  degree(A) = "
        f"{degree_string(A_full, N, Delta)}"
    )

    print(
        f"  degree(B) = "
        f"{degree_string(B_full, N, Delta)}"
    )

    print(
        f"  monomial_count(A) = "
        f"{monomial_count(A_full, N, Delta)}"
    )

    print(
        f"  monomial_count(B) = "
        f"{monomial_count(B_full, N, Delta)}"
    )

    B_is_zero = (
        sp.expand(B_full) == 0
    )

    print(
        f"  B == 0: {B_is_zero}"
    )

    print(
        "\n  A(N,Delta) ="
    )
    print(
        f"    {A_full}"
    )

    print(
        "\n  B(N,Delta) ="
    )
    print(
        f"    {B_full}"
    )

    full_reconstruction = (
        sp.expand(
            reconstruct_from_N_Delta(
                A_full,
                B_full,
            )
            -
            G_NS
        )
        == 0
    )

    print(
        "\n  exact reduction reconstruction = "
        f"{full_reconstruction}"
    )

    if not full_reconstruction:
        raise RuntimeError(
            "Full discriminant reduction reconstruction failed."
        )

    # ----------------------------------------------------------------------
    # 3. LAYER-BY-LAYER REDUCTION
    # ----------------------------------------------------------------------

    print(
        "\n3. LAYER-BY-LAYER DISCRIMINANT REDUCTION"
    )

    reduced_layers = {}

    for degree in degrees:

        layer_NS = sp.expand(
            layers[degree].subs(
                X,
                S + 1,
            )
        )

        A_d, B_d = reduce_to_N_Delta(
            layer_NS
        )

        reduced_layers[degree] = (
            A_d,
            B_d,
        )

        print(
            f"\n  degree={degree}"
        )

        print(
            f"    A degree="
            f"{degree_string(A_d, N, Delta)}"
        )

        print(
            f"    B degree="
            f"{degree_string(B_d, N, Delta)}"
        )

        print(
            f"    A monomials="
            f"{monomial_count(A_d, N, Delta)}"
        )

        print(
            f"    B monomials="
            f"{monomial_count(B_d, N, Delta)}"
        )

        print(
            f"    B_zero="
            f"{sp.expand(B_d) == 0}"
        )

    # ----------------------------------------------------------------------
    # 4. SPARSITY SUMMARY
    # ----------------------------------------------------------------------

    print(
        "\n4. LAYER SPARSITY SUMMARY"
    )

    print(
        "  degree | A_terms | B_terms | B_zero"
    )

    print(
        "  " + "-" * 52
    )

    for degree in degrees:

        A_d, B_d = reduced_layers[degree]

        print(
            f"  {degree:6d} | "
            f"{monomial_count(A_d, N, Delta):7d} | "
            f"{monomial_count(B_d, N, Delta):7d} | "
            f"{sp.expand(B_d) == 0}"
        )

    # ----------------------------------------------------------------------
    # 5. PRIMITIVE SIGNATURES
    # ----------------------------------------------------------------------

    print(
        "\n5. PRIMITIVE INTEGER SIGNATURES"
    )

    for degree in degrees:

        A_d, B_d = reduced_layers[degree]

        print(
            f"\n  degree={degree}"
        )

        print(
            "    A primitive="
            f"{primitive_integer_coefficients(A_d, N, Delta)}"
        )

        print(
            "    B primitive="
            f"{primitive_integer_coefficients(B_d, N, Delta)}"
        )

    # ----------------------------------------------------------------------
    # 6. FRESH RANDOM SEMIPRIME
    # ----------------------------------------------------------------------

    p_value, q_value = (
        random_prime_pair()
    )

    n_value = (
        p_value
        *
        q_value
    )

    s_value = (
        p_value
        +
        q_value
    )

    x_value = (
        s_value
        +
        1
    )

    phi_value = (
        (p_value - 1)
        *
        (q_value - 1)
    )

    delta_value = (
        s_value**2
        -
        4*n_value
    )

    delta_direct = (
        p_value
        -
        q_value
    )**2

    print(
        "\n6. FRESH RANDOM SEMIPRIME"
    )

    print(
        f"  p = {p_value}"
    )

    print(
        f"  q = {q_value}"
    )

    print(
        f"  p != q = "
        f"{p_value != q_value}"
    )

    print(
        f"  n = p*q = {n_value}"
    )

    print(
        f"  S = p+q = {s_value}"
    )

    print(
        f"  X = p+q+1 = {x_value}"
    )

    print(
        f"  phi(n) = {phi_value}"
    )

    print(
        "\n  discriminant"
    )

    print(
        f"    Delta = S^2 - 4N = "
        f"{delta_value}"
    )

    print(
        f"    (p-q)^2 = "
        f"{delta_direct}"
    )

    delta_identity = (
        delta_value
        ==
        delta_direct
    )

    print(
        f"    exact equality = "
        f"{delta_identity}"
    )

    print(
        "\n  totient bridge"
    )

    phi_from_N_X = (
        n_value
        -
        x_value
        +
        2
    )

    phi_identity = (
        phi_from_N_X
        ==
        phi_value
    )

    print(
        f"    N-X+2 = {phi_from_N_X}"
    )

    print(
        f"    phi(n) = {phi_value}"
    )

    print(
        f"    exact equality = "
        f"{phi_identity}"
    )

    # ----------------------------------------------------------------------
    # 7. FRESH EVALUATION
    # ----------------------------------------------------------------------

    print(
        "\n7. FRESH SEMIPRIME EVALUATION"
    )

    F_value = exact_eval(
        F_symbolic,
        {
            p: p_value,
            q: q_value,
        },
    )

    G_value = exact_eval(
        G,
        {
            N: n_value,
            X: x_value,
        },
    )

    A_value = exact_eval(
        A_full,
        {
            N: n_value,
            Delta: delta_value,
        },
    )

    B_value = exact_eval(
        B_full,
        {
            N: n_value,
            Delta: delta_value,
        },
    )

    SB_value = sp.Integer(
        s_value
        *
        B_value
    )

    reduced_value = sp.Integer(
        A_value
        +
        SB_value
    )

    print(
        f"  F(p,q) = {F_value}"
    )

    print(
        f"  G(N,X) = {G_value}"
    )

    print(
        f"  A(N,Delta) = {A_value}"
    )

    print(
        f"  B(N,Delta) = {B_value}"
    )

    print(
        f"  S*B = {SB_value}"
    )

    print(
        f"  A + S*B = {reduced_value}"
    )

    fresh_FG = (
        F_value
        ==
        G_value
    )

    fresh_reduced = (
        reduced_value
        ==
        G_value
    )

    print(
        f"\n  F == G: {fresh_FG}"
    )

    print(
        f"  reduced == G: "
        f"{fresh_reduced}"
    )

    # ----------------------------------------------------------------------
    # 8. DELTA-ONLY TEST
    # ----------------------------------------------------------------------

    print(
        "\n8. DELTA-ONLY TEST"
    )

    if B_is_zero:
        print(
            "  RESULT: EXACT DELTA-ONLY FORM"
        )

        print(
            "  G is a polynomial in N and Delta alone."
        )

    else:
        print(
            "  RESULT: NONZERO S-COMPONENT"
        )

        print(
            "  G is NOT a polynomial in N and Delta alone."
        )

    # ----------------------------------------------------------------------
    # 9. LAYER DELTA-ONLY PROFILE
    # ----------------------------------------------------------------------

    print(
        "\n9. LAYER DELTA-ONLY PROFILE"
    )

    delta_only_degrees = []

    for degree in degrees:

        _, B_d = reduced_layers[degree]

        if sp.expand(B_d) == 0:
            delta_only_degrees.append(
                degree
            )

    print(
        f"  B_d=0 layers = "
        f"{delta_only_degrees}"
    )

    print(
        f"  count = "
        f"{len(delta_only_degrees)} / {len(degrees)}"
    )

    # ----------------------------------------------------------------------
    # 10. FACTORIZATION
    # ----------------------------------------------------------------------

    print(
        "\n10. FACTORIZATION AUDIT"
    )

    print(
        "\n  factor(A_full) ="
    )

    print(
        f"    {sp.factor(A_full)}"
    )

    print(
        "\n  factor(B_full) ="
    )

    print(
        f"    {sp.factor(B_full)}"
    )

    # ----------------------------------------------------------------------
    # 11. EXACT n=pq BRIDGE
    # ----------------------------------------------------------------------

    print(
        "\n11. EXACT n=pq BRIDGE"
    )

    print(
        "  N = p*q"
    )

    print(
        "  S = p+q"
    )

    print(
        "  X = S+1"
    )

    print(
        "  Delta = S^2 - 4N = (p-q)^2"
    )

    print(
        "  phi(n) = N-S+1 = N-X+2"
    )

    print(
        "\n  factor reconstruction equation:"
    )

    print(
        "    t^2 - S*t + N = 0"
    )

    print(
        "    discriminant = Delta"
    )

    # ----------------------------------------------------------------------
    # 12. FIXED-N X CONTROL
    # ----------------------------------------------------------------------

    print(
        "\n12. FIXED-N X CONTROL"
    )

    offsets = [-2, -1, 0, 1, 2]

    control_results = []

    for offset in offsets:

        X_control = (
            x_value
            +
            offset
        )

        value = exact_eval(
            G,
            {
                N: n_value,
                X: X_control,
            },
        )

        control_results.append(
            (
                offset,
                X_control,
                value,
            )
        )

    for offset, X_control, value in control_results:

        print(
            f"  offset={offset:+d} "
            f"X={X_control} "
            f"G={value}"
        )

    independent_of_X = all(
        value == control_results[0][2]
        for _, _, value in control_results
    )

    print(
        f"\n  G independent of X: "
        f"{independent_of_X}"
    )

    # This is NOT an error condition.
    # In fact, for the present kernel we expect False.
    x_dependence_detected = (
        not independent_of_X
    )

    # ----------------------------------------------------------------------
    # 13. EXACTNESS CONTROLS
    # ----------------------------------------------------------------------

    print(
        "\n13. EXACTNESS CONTROLS"
    )

    checks = {
        "G contains Float": (
            G.has(sp.Float)
        ),

        "A contains Float": (
            A_full.has(sp.Float)
        ),

        "B contains Float": (
            B_full.has(sp.Float)
        ),

        "F fresh contains Float": (
            F_value.has(sp.Float)
        ),

        "G fresh contains Float": (
            G_value.has(sp.Float)
        ),

        "Delta identity FAILED": (
            not delta_identity
        ),

        "Totient identity FAILED": (
            not phi_identity
        ),

        "Symbolic F/G identity FAILED": (
            not symbolic_G_check
        ),

        "Full A/B reconstruction FAILED": (
            not full_reconstruction
        ),

        "Fresh F/G equality FAILED": (
            not fresh_FG
        ),

        "Fresh reduced equality FAILED": (
            not fresh_reduced
        ),
    }

    failures = [
        name
        for name, failed in checks.items()
        if failed
    ]

    for name, failed in checks.items():

        print(
            f"  {name}: {failed}"
        )

    # ----------------------------------------------------------------------
    # FINAL STATUS
    # ----------------------------------------------------------------------

    print(
        "\n14. FINAL STATUS"
    )

    if failures:

        print(
            "  FAILED CHECKS:"
        )

        for failure in failures:
            print(
                f"    - {failure}"
            )

        raise RuntimeError(
            "One or more exactness checks failed."
        )

    print(
        "  All exactness checks passed."
    )

    print(
        f"  X-dependence detected: "
        f"{x_dependence_detected}"
    )

    print(
        "  NOTE: X-dependence is an observed structural property, "
        "not an error."
    )

    print(
        "\n" + "=" * 78
    )

    print(
        "EXPERIMENT 68 COMPLETE — "
        "ALL EXACT CHECKS PASSED"
    )

    print(
        "=" * 78
    )


if __name__ == "__main__":
    main()