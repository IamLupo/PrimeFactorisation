import sympy as sp


# =============================================================================
# EXPERIMENT 215
# TRUE L2 BOUNDARY: a=k AND a=k+1
# =============================================================================
#
# No previous experiment output is read.
# No unrestricted candidate enumeration is performed.
#
# Established facts:
#
#   G_top = -X^(ell-k) * ((X+N)^k - N^k)
#
#   L1_a = ell*C(k+1,a) - k*C(k,a-1)
#
#   L2 has degree ell-2 and support
#
#       a = 0,...,k+1
#
# The previous experiment incorrectly treated a=k as the final endpoint.
#
# This experiment:
#
#   1. verifies the complete L2 support;
#   2. derives D_k(ell);
#   3. derives the true boundary D_(k+1)(ell);
#   4. checks whether both are quadratic in ell;
#   5. compares D_k with the interior formula;
#   6. extracts the exact correction at a=k;
#   7. tests whether D_(k+1) is constant, linear, or quadratic in ell;
#   8. performs exact reconstruction of L2 using the two boundary rows.
#
# The computation is deliberately bounded.
# =============================================================================


p, q = sp.symbols("p q")
N, S, X = sp.symbols("N S X")
e = sp.symbols("e")


# =============================================================================
# BASIC POLYNOMIAL UTILITIES
# =============================================================================

def poly_NX(expr):
    expr = sp.cancel(sp.together(sp.expand(expr)))

    den = sp.expand(sp.denom(expr))

    if den != 1:
        raise ValueError(
            f"Non-polynomial expression in N,X: denominator={den}"
        )

    return sp.Poly(
        sp.expand(expr),
        N,
        X,
        domain=sp.QQ,
    )


def coefficient(expr, a, b):
    if a < 0 or b < 0:
        return sp.Integer(0)

    return sp.expand(
        poly_NX(expr).coeff_monomial(
            N**a * X**b
        )
    )


def homogeneous(expr, degree):
    if degree < 0:
        return sp.Integer(0)

    poly = poly_NX(expr)

    result = sp.Integer(0)

    for (a, b), c in poly.terms():

        if a + b == degree:
            result += (
                c
                * N**a
                * X**b
            )

    return sp.expand(result)


# =============================================================================
# EXACT pq KERNEL
# =============================================================================

def exact_F(k, ell):

    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


# =============================================================================
# DIRECT SYMMETRIC REDUCTION
# =============================================================================

def exact_G(k, ell):

    F = exact_F(k, ell)

    result, remainder, mapping = sp.symmetrize(
        F,
        [p, q],
        formal=True,
    )

    if sp.expand(remainder) != 0:
        raise ValueError(
            f"Symmetrization remainder: {remainder}"
        )

    sum_symbol = None
    prod_symbol = None

    for symbol, expr in mapping:

        expr = sp.expand(expr)

        if sp.expand(expr - (p + q)) == 0:
            sum_symbol = symbol

        elif sp.expand(expr - p*q) == 0:
            prod_symbol = symbol

    if sum_symbol is None or prod_symbol is None:
        raise ValueError(
            f"Could not identify symmetric variables: {mapping}"
        )

    return sp.expand(
        result.subs(
            {
                sum_symbol: S,
                prod_symbol: N,
            }
        )
    )


def exact_GX(k, ell):

    return sp.expand(
        exact_G(k, ell).subs(
            S,
            X - 1,
        )
    )


# =============================================================================
# TOP LAYER
# =============================================================================

def top_layer(k, ell):

    return sp.expand(
        -X**(ell - k)
        * (
            (X + N)**k
            - N**k
        )
    )


# =============================================================================
# L1
# =============================================================================

def L1_coefficient(k, ell, a):

    if a < 0 or a > k:
        return sp.Integer(0)

    return sp.expand(
        ell * sp.binomial(k + 1, a)
        - k * sp.binomial(k, a - 1)
    )


def L1(k, ell):

    degree = ell - 1

    result = sp.Integer(0)

    for a in range(k + 1):

        b = degree - a

        if b < 0:
            continue

        result += (
            L1_coefficient(k, ell, a)
            * N**a
            * X**b
        )

    return sp.expand(result)


# =============================================================================
# EXACT L2
# =============================================================================

L2_CACHE = {}


def exact_L2(k, ell):

    key = (k, ell)

    if key in L2_CACHE:
        return L2_CACHE[key]

    G = exact_GX(k, ell)

    residual = sp.expand(
        G
        - top_layer(k, ell)
        - L1(k, ell)
    )

    result = homogeneous(
        residual,
        ell - 2,
    )

    L2_CACHE[key] = result

    return result


def D(k, ell, a):

    b = ell - 2 - a

    if b < 0:
        return sp.Integer(0)

    return coefficient(
        exact_L2(k, ell),
        a,
        b,
    )


# =============================================================================
# QUADRATIC EXTRACTION
# =============================================================================

def fit_quadratic(k, a, ell0):

    values = [
        (
            ell0,
            D(k, ell0, a),
        ),
        (
            ell0 + 2,
            D(k, ell0 + 2, a),
        ),
        (
            ell0 + 4,
            D(k, ell0 + 4, a),
        ),
    ]

    polynomial = sp.interpolate(
        values,
        e,
    )

    polynomial = sp.factor(
        sp.expand(polynomial)
    )

    P = sp.Poly(
        polynomial,
        e,
    )

    A = sp.expand(
        P.coeff_monomial(e**2)
    )

    B = sp.expand(
        P.coeff_monomial(e)
    )

    C = sp.expand(
        P.coeff_monomial(1)
    )

    return (
        polynomial,
        A,
        B,
        C,
    )


# =============================================================================
# INTERIOR FORMULAS FROM 213
# =============================================================================

def A_interior(k, a):

    return sp.expand(
        -sp.binomial(k + 2, a) / 2
    )


def B_interior(k, a):

    return sp.expand(
        sp.binomial(k + 2, a)
        * (
            2 * (k + 1) * a
            + k + 2
        )
        / (
            2 * (k + 2)
        )
    )


def C_interior(k, a):

    return sp.expand(
        -sp.Rational(1, 2)
        * k
        * a
        * (a + 1)
        * sp.binomial(k + 2, a)
        / (k + 2)
    )


# =============================================================================
# 1. SUPPORT AUDIT
# =============================================================================

def support_audit():

    print()
    print("=" * 78)
    print("1. COMPLETE L2 SUPPORT AUDIT")
    print("=" * 78)

    cases = [
        (3, 7),
        (3, 9),
        (5, 11),
        (5, 19),
        (7, 15),
        (9, 21),
        (11, 23),
    ]

    failures = 0

    for k, ell in cases:

        L = exact_L2(k, ell)
        P = poly_NX(L)

        support = sorted(
            {
                a
                for (a, b), c in P.terms()
                if c != 0
            }
        )

        expected = list(
            range(
                0,
                min(k + 1, ell - 2) + 1,
            )
        )

        # For the chosen cases ell is comfortably larger
        # than k, so the expected support is 0,...,k+1.
        ok = support == expected

        if not ok:
            failures += 1

        print(
            f"k={k:2d} ell={ell:2d} "
            f"support={support} "
            f"expected={expected} "
            f"{'PASS' if ok else 'FAIL'}"
        )

    print()
    print(
        f"support failures = {failures}"
    )


# =============================================================================
# 2. a=k-1, k, k+1 TABLE
# =============================================================================

def boundary_rows():

    print()
    print("=" * 78)
    print("2. LAST THREE L2 ROWS")
    print("=" * 78)

    cases = {
        3: [7, 9, 11, 13],
        5: [11, 13, 15, 17],
        7: [15, 17],
        9: [21, 23],
        11: [23, 25],
    }

    for k, ells in cases.items():

        print()
        print(f"k={k}")

        for ell in ells:

            vals = []

            for a in [
                k - 1,
                k,
                k + 1,
            ]:

                vals.append(
                    (
                        a,
                        D(k, ell, a),
                    )
                )

            print(
                f"  ell={ell:2d}: "
                + "  ".join(
                    f"a={a}: {value}"
                    for a, value in vals
                )
            )


# =============================================================================
# 3. EXACT a=k QUADRATIC LAW
# =============================================================================

def a_k_audit():

    print()
    print("=" * 78)
    print("3. EXACT a=k QUADRATIC LAW")
    print("=" * 78)

    cases = {
        3: 7,
        5: 11,
        7: 15,
        9: 21,
        11: 23,
    }

    failures = 0

    for k, ell0 in cases.items():

        polynomial, A, B, C = fit_quadratic(
            k,
            k,
            ell0,
        )

        Aint = A_interior(
            k,
            k,
        )

        Bint = B_interior(
            k,
            k,
        )

        Cint = C_interior(
            k,
            k,
        )

        correction_B = sp.expand(
            B - Bint
        )

        correction_C = sp.expand(
            C - Cint
        )

        print()
        print(
            f"k={k}"
        )

        print(
            f"  D_k(ell) = {polynomial}"
        )

        print(
            f"  A = {A}"
        )

        print(
            f"  A_interior = {Aint}"
        )

        print(
            f"  B = {B}"
        )

        print(
            f"  B_interior = {Bint}"
        )

        print(
            f"  B correction = {correction_B}"
        )

        print(
            f"  C = {C}"
        )

        print(
            f"  C_interior = {Cint}"
        )

        print(
            f"  C correction = {correction_C}"
        )


# =============================================================================
# 4. TRUE ENDPOINT a=k+1
# =============================================================================

def a_k_plus_one_audit():

    print()
    print("=" * 78)
    print("4. TRUE ENDPOINT a=k+1")
    print("=" * 78)

    cases = {
        3: [7, 9, 11, 13, 15],
        5: [11, 13, 15, 17],
        7: [15, 17, 19],
        9: [21, 23],
        11: [23, 25],
    }

    for k, ells in cases.items():

        print()
        print(f"k={k}")

        for ell in ells:

            value = D(
                k,
                ell,
                k + 1,
            )

            print(
                f"  ell={ell:2d} "
                f"D_(k+1)={value}"
            )


# =============================================================================
# 5. TRUE ENDPOINT ELL-DEPENDENCE
# =============================================================================

def endpoint_ell_dependence():

    print()
    print("=" * 78)
    print("5. TRUE ENDPOINT ELL-DEPENDENCE")
    print("=" * 78)

    cases = {
        3: [7, 9, 11, 13, 15],
        5: [11, 13, 15, 17],
        7: [15, 17, 19],
        9: [21, 23],
        11: [23, 25],
    }

    for k, ells in cases.items():

        values = [
            D(
                k,
                ell,
                k + 1,
            )
            for ell in ells
        ]

        print()
        print(
            f"k={k}"
        )

        print(
            f"  ell={ells}"
        )

        print(
            f"  D_(k+1)={values}"
        )

        if len(values) >= 2:

            first = [
                sp.expand(
                    values[i + 1]
                    - values[i]
                )
                for i in range(
                    len(values) - 1
                )
            ]

            print(
                f"  Delta={first}"
            )

            second = [
                sp.expand(
                    first[i + 1]
                    - first[i]
                )
                for i in range(
                    len(first) - 1
                )
            ]

            if second:
                print(
                    f"  Delta^2={second}"
                )


# =============================================================================
# 6. QUADRATIC FIT FOR TRUE ENDPOINT
# =============================================================================

def endpoint_quadratic_fit():

    print()
    print("=" * 78)
    print("6. TRUE ENDPOINT QUADRATIC FIT")
    print("=" * 78)

    cases = {
        3: 7,
        5: 11,
        7: 15,
        9: 21,
        11: 23,
    }

    for k, ell0 in cases.items():

        polynomial, A, B, C = fit_quadratic(
            k,
            k + 1,
            ell0,
        )

        print()
        print(
            f"k={k}"
        )

        print(
            f"  D_(k+1)(ell) = {polynomial}"
        )

        print(
            f"  A={A}"
        )

        print(
            f"  B={B}"
        )

        print(
            f"  C={C}"
        )


# =============================================================================
# 7. TRUE ENDPOINT FACTORIZATION
# =============================================================================

def endpoint_factorization():

    print()
    print("=" * 78)
    print("7. TRUE ENDPOINT FACTORIZATION")
    print("=" * 78)

    cases = {
        3: 7,
        5: 11,
        7: 15,
        9: 21,
        11: 23,
    }

    for k, ell0 in cases.items():

        polynomial, A, B, C = fit_quadratic(
            k,
            k + 1,
            ell0,
        )

        print()
        print(
            f"k={k}"
        )

        print(
            "  expanded =",
            polynomial,
        )

        print(
            "  factored =",
            sp.factor(polynomial),
        )

        print(
            "  roots =",
            sp.solve(
                sp.Eq(
                    polynomial,
                    0,
                ),
                e,
            ),
        )


# =============================================================================
# 8. SMALL NATURAL ENDPOINT NORMALIZATIONS
# =============================================================================

def endpoint_normalizations():

    print()
    print("=" * 78)
    print("8. TRUE ENDPOINT NORMALIZATIONS")
    print("=" * 78)

    cases = {
        3: 11,
        5: 15,
        7: 19,
        9: 23,
    }

    for k, ell in cases.items():

        value = D(
            k,
            ell,
            k + 1,
        )

        print()
        print(
            f"k={k} ell={ell} D={value}"
        )

        candidates = {
            "D": value,

            "D / k": sp.factor(
                value / k
            ),

            "D / (k+1)": sp.factor(
                value / (k + 1)
            ),

            "D / C(k+2,2)": sp.factor(
                value
                / sp.binomial(
                    k + 2,
                    2,
                )
            ),

            "D / C(k+2,3)": sp.factor(
                value
                / sp.binomial(
                    k + 2,
                    3,
                )
            ),

            "D / C(k+2,k)": sp.factor(
                value
                / sp.binomial(
                    k + 2,
                    k,
                )
            ),
        }

        for name, normalized in candidates.items():

            print(
                f"  {name:22s} = {normalized}"
            )


# =============================================================================
# 9. COMPLETE L2 RECONSTRUCTION
# =============================================================================

def complete_reconstruction():

    print()
    print("=" * 78)
    print("9. COMPLETE L2 RECONSTRUCTION")
    print("=" * 78)

    cases = [
        (3, 7),
        (3, 9),
        (3, 11),

        (5, 11),
        (5, 13),
        (5, 15),

        (7, 15),
        (7, 17),

        (9, 21),

        (11, 23),
    ]

    failures = 0

    for k, ell in cases:

        exact = exact_L2(
            k,
            ell,
        )

        predicted = sp.Integer(0)

        # Interior rows a=0,...,k-1
        for a in range(k):

            b = ell - 2 - a

            if b < 0:
                continue

            A = A_interior(
                k,
                a,
            )

            B = B_interior(
                k,
                a,
            )

            C = C_interior(
                k,
                a,
            )

            value = sp.expand(
                A * ell**2
                + B * ell
                + C
            )

            predicted += (
                value
                * N**a
                * X**b
            )

        # a=k
        _, A, B, C = fit_quadratic(
            k,
            k,
            ell,
        )

        value = sp.expand(
            A * ell**2
            + B * ell
            + C
        )

        predicted += (
            value
            * N**k
            * X**(ell - 2 - k)
        )

        # a=k+1: use the exact fitted quadratic
        # temporarily, only to determine whether
        # the support decomposition itself is correct.
        _, A, B, C = fit_quadratic(
            k,
            k + 1,
            ell,
        )

        value = sp.expand(
            A * ell**2
            + B * ell
            + C
        )

        predicted += (
            value
            * N**(k + 1)
            * X**(ell - 3 - k)
        )

        ok = (
            sp.expand(
                exact - predicted
            )
            == 0
        )

        if not ok:
            failures += 1

        print(
            f"k={k:2d} ell={ell:2d} "
            f"reconstruction="
            f"{'PASS' if ok else 'FAIL'}"
        )

    print()
    print(
        f"complete L2 reconstruction failures = {failures}"
    )


# =============================================================================
# 10. FINAL DIAGNOSTIC
# =============================================================================

def final_diagnostic():

    print()
    print("=" * 78)
    print("10. FINAL DIAGNOSTIC")
    print("=" * 78)

    print()
    print(
        "The correct L2 support is:"
    )

    print(
        "  a = 0,...,k+1"
    )

    print()
    print(
        "Therefore the previous 'endpoint a=k' analysis"
    )

    print(
        "was one step short of the actual boundary."
    )

    print()
    print(
        "The next decisive questions are:"
    )

    print(
        "  1. Does a=k satisfy the same interior formula"
    )

    print(
        "     or require a correction?"
    )

    print(
        "  2. What is the exact quadratic law for a=k+1?"
    )

    print(
        "  3. Can the a=k+1 law be derived directly from"
    )

    print(
        "     the pq kernel rather than interpolated?"
    )

    print()
    print(
        "If the true endpoint simplifies, we will have"
    )

    print(
        "the complete L2 layer."
    )

    print()
    print(
        "Only then do we peel L3."
    )

    print()
    print(
        "Methodological rule:"
    )

    print(
        "  pq kernel -> homogeneous layer -> boundary term"
    )

    print(
        "not"
    )

    print(
        "  anchors -> fitted formula"
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print(
        "EXPERIMENT 215"
    )
    print(
        "TRUE L2 BOUNDARY: a=k AND a=k+1"
    )
    print("=" * 78)

    print()
    print(
        "No previous experiment output is read."
    )

    print(
        "No unrestricted candidate enumeration is performed."
    )

    support_audit()
    boundary_rows()
    a_k_audit()
    a_k_plus_one_audit()
    endpoint_ell_dependence()
    endpoint_quadratic_fit()
    endpoint_factorization()
    endpoint_normalizations()
    complete_reconstruction()
    final_diagnostic()

    print()
    print("=" * 78)
    print(
        "END EXPERIMENT 215"
    )
    print("=" * 78)


if __name__ == "__main__":
    main()

