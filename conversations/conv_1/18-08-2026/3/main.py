import sympy as sp


# =============================================================================
# EXPERIMENT 214
# EXACT L2 ENDPOINT / BOUNDARY-CORRECTION AUDIT
# =============================================================================
#
# Purpose:
#
#   Experiments 209R-213 established:
#
#       G = G_top + L1 + L2 + ...
#
#   with
#
#       G_top = -X^(ell-k) ((X+N)^k - N^k)
#
#   and
#
#       [N^a X^(ell-1-a)] L1
#           = ell*C(k+1,a) - k*C(k,a-1).
#
#   Experiment 213 established:
#
#       D_a(k,ell)
#         = [N^a X^(ell-2-a)] L2
#         = A(k,a) ell^2 + B(k,a) ell + C(k,a).
#
#   The exact interior laws are:
#
#       A(k,a) = -1/2*C(k+2,a)
#
#       B(k,a)
#         = C(k+2,a) * (2(k+1)a + k+2)
#             / (2(k+2))
#
#   for 0 <= a < k.
#
#   Experiment 213 then showed:
#
#       B(k,k) differs from that interior formula
#       by exactly -(k+1).
#
#   The same experiment showed that
#
#       C(k,a)
#         = - k*a*(a+1)/(2*(k+2))*C(k+2,a)
#
#   for all tested 0 <= a < k,
#
#   but did NOT explain C(k,k).
#
# This experiment therefore:
#
#   1. verifies the complete interior L2 formula;
#   2. isolates a=k-1 and a=k;
#   3. verifies that the endpoint coefficients are independent of ell;
#   4. measures the B endpoint correction;
#   5. extracts C(k,k) for many k;
#   6. factors and finite-differences the endpoint sequence;
#   7. tests a small number of structurally motivated endpoint formulas;
#   8. reconstructs the complete L2 layer using the resulting piecewise law.
#
# No unrestricted candidate search is used.
# No previous result files are read.
# =============================================================================


p, q = sp.symbols("p q")
N, S, X = sp.symbols("N S X")
e = sp.symbols("e")


# =============================================================================
# POLYNOMIAL UTILITIES
# =============================================================================

def poly_NX(expr):
    expr = sp.cancel(sp.together(sp.expand(expr)))

    if sp.expand(sp.denom(expr) - 1) != 0:
        raise ValueError(
            f"Non-polynomial expression in N,X: denominator={sp.denom(expr)}"
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

    out = sp.Integer(0)

    for (a, b), c in poly.terms():
        if a + b == degree:
            out += c * N**a * X**b

    return sp.expand(out)


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

    out = sp.Integer(0)

    for a in range(k + 1):

        b = degree - a

        if b < 0:
            continue

        out += (
            L1_coefficient(k, ell, a)
            * N**a
            * X**b
        )

    return sp.expand(out)


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

    L2_expr = homogeneous(
        residual,
        ell - 2,
    )

    L2_CACHE[key] = L2_expr

    return L2_expr


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

def quadratic_coefficients(k, a, ell0):

    points = [
        (ell0, D(k, ell0, a)),
        (ell0 + 2, D(k, ell0 + 2, a)),
        (ell0 + 4, D(k, ell0 + 4, a)),
    ]

    polynomial = sp.interpolate(
        points,
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

    return polynomial, A, B, C


# =============================================================================
# INTERIOR FORMULAS
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
# 1. INTERIOR FORMULA AUDIT
# =============================================================================

def interior_audit():

    print()
    print("=" * 78)
    print("1. COMPLETE INTERIOR L2 FORMULA")
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

        print()
        print(f"k={k}")

        for a in range(k):

            _, A, B, C = quadratic_coefficients(
                k,
                a,
                ell0,
            )

            Aexp = A_interior(k, a)
            Bexp = B_interior(k, a)
            Cexp = C_interior(k, a)

            okA = sp.expand(A - Aexp) == 0
            okB = sp.expand(B - Bexp) == 0
            okC = sp.expand(C - Cexp) == 0

            if not (okA and okB and okC):
                failures += 1

            print(
                f"  a={a:2d} "
                f"A={'PASS' if okA else 'FAIL'} "
                f"B={'PASS' if okB else 'FAIL'} "
                f"C={'PASS' if okC else 'FAIL'}"
            )

    print()
    print(
        f"interior formula failures = {failures}"
    )


# =============================================================================
# 2. ENDPOINT B CORRECTION
# =============================================================================

def endpoint_B_audit():

    print()
    print("=" * 78)
    print("2. ENDPOINT B CORRECTION")
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

        _, A, B, C = quadratic_coefficients(
            k,
            k,
            ell0,
        )

        interior = B_interior(
            k,
            k,
        )

        correction = sp.expand(
            B - interior
        )

        expected = -sp.Integer(
            k + 1
        )

        ok = (
            sp.expand(
                correction - expected
            )
            == 0
        )

        if not ok:
            failures += 1

        print(
            f"k={k:2d} "
            f"B_endpoint={B} "
            f"B_interior={interior} "
            f"correction={correction} "
            f"expected={expected} "
            f"{'PASS' if ok else 'FAIL'}"
        )

    print()
    print(
        f"B endpoint failures = {failures}"
    )


# =============================================================================
# 3. ENDPOINT C EXTRACTION
# =============================================================================

def endpoint_C_table():

    print()
    print("=" * 78)
    print("3. EXACT ENDPOINT C(k,k)")
    print("=" * 78)

    cases = {
        1: 5,
        3: 7,
        5: 11,
        7: 15,
        9: 21,
        11: 23,
        13: 27,
    }

    values = []

    for k, ell0 in cases.items():

        _, A, B, C = quadratic_coefficients(
            k,
            k,
            ell0,
        )

        values.append(
            (k, sp.expand(C))
        )

        print(
            f"k={k:2d} "
            f"A={A} "
            f"B={B} "
            f"C(k,k)={C} "
            f"factor={sp.factor(C)}"
        )

    return values


# =============================================================================
# 4. ENDPOINT C INDEPENDENCE FROM ell
# =============================================================================

def endpoint_C_ell_audit():

    print()
    print("=" * 78)
    print("4. ENDPOINT C INDEPENDENCE FROM ell")
    print("=" * 78)

    cases = {
        3: [7, 9, 11, 13],
        5: [11, 13, 15, 17],
        7: [15, 17, 19],
        9: [21, 23],
    }

    failures = 0

    for k, ells in cases.items():

        observed = []

        for ell0 in ells:

            _, _, _, C = quadratic_coefficients(
                k,
                k,
                ell0,
            )

            observed.append(
                sp.expand(C)
            )

        ok = all(
            value == observed[0]
            for value in observed
        )

        if not ok:
            failures += 1

        print(
            f"k={k:2d} "
            f"C-values={observed} "
            f"{'PASS' if ok else 'FAIL'}"
        )

    print()
    print(
        f"endpoint ell-independence failures = {failures}"
    )


# =============================================================================
# 5. NATURAL ENDPOINT C NORMALIZATIONS
# =============================================================================

def endpoint_C_normalization():

    print()
    print("=" * 78)
    print("5. ENDPOINT C NORMALIZATION")
    print("=" * 78)

    cases = {
        1: 5,
        3: 7,
        5: 11,
        7: 15,
        9: 21,
        11: 23,
        13: 27,
    }

    for k, ell0 in cases.items():

        _, _, _, C = quadratic_coefficients(
            k,
            k,
            ell0,
        )

        candidates = {
            "C / k": sp.factor(C / k),

            "C / C(k+2,k)": sp.factor(
                C / sp.binomial(k + 2, k)
            ),

            "C / C(k+1,k)": sp.factor(
                C / sp.binomial(k + 1, k)
            ),

            "C / C(k+1,2)": sp.factor(
                C / sp.binomial(k + 1, 2)
            ),

            "C / C(k+2,3)": sp.factor(
                C / sp.binomial(k + 2, 3)
            ),

            "C / (k*C(k+1,2))": sp.factor(
                C
                / (
                    k
                    * sp.binomial(k + 1, 2)
                )
            ),

            "C / (k*C(k+2,3))": sp.factor(
                C
                / (
                    k
                    * sp.binomial(k + 2, 3)
                )
            ),
        }

        print()
        print(
            f"k={k:2d} C={C}"
        )

        for name, value in candidates.items():
            print(
                f"  {name:25s} = {value}"
            )


# =============================================================================
# 6. FINITE-DIFFERENCE TEST OF C(k,k)
# =============================================================================

def endpoint_C_difference():

    print()
    print("=" * 78)
    print("6. FINITE-DIFFERENCE TEST OF C(k,k)")
    print("=" * 78)

    cases = {
        1: 5,
        3: 7,
        5: 11,
        7: 15,
        9: 21,
        11: 23,
        13: 27,
    }

    values = []

    for k, ell0 in cases.items():

        _, _, _, C = quadratic_coefficients(
            k,
            k,
            ell0,
        )

        values.append(
            (k, sp.Integer(C))
        )

    print()
    print(
        "k-values:",
        [k for k, _ in values],
    )

    print(
        "C(k,k):",
        [c for _, c in values],
    )

    current = [
        c for _, c in values
    ]

    order = 1

    while len(current) >= 2:

        current = [
            sp.expand(
                current[i + 1]
                - current[i]
            )
            for i in range(
                len(current) - 1
            )
        ]

        print(
            f"Delta^{order}:",
            current,
        )

        order += 1


# =============================================================================
# 7. ENDPOINT C SYMBOLIC INTERPOLATION
# =============================================================================
#
# This is deliberately used only as a diagnostic:
# if the endpoint sequence is polynomial in k on the odd branch,
# factor the resulting polynomial.
#
# It is NOT used as evidence by itself.
# =============================================================================

def endpoint_C_interpolation():

    print()
    print("=" * 78)
    print("7. ENDPOINT C ODD-k INTERPOLATION DIAGNOSTIC")
    print("=" * 78)

    data = []

    for k in [1, 3, 5, 7, 9, 11, 13]:

        ell0 = {
            1: 5,
            3: 7,
            5: 11,
            7: 15,
            9: 21,
            11: 23,
            13: 27,
        }[k]

        _, _, _, C = quadratic_coefficients(
            k,
            k,
            ell0,
        )

        data.append(
            (sp.Integer(k), sp.Integer(C))
        )

    polynomial = sp.interpolate(
        data,
        sp.Symbol("K"),
    )

    polynomial = sp.factor(
        sp.expand(polynomial)
    )

    print()
    print(
        "Interpolating polynomial through the"
    )
    print(
        "computed odd-k endpoint values:"
    )
    print()
    print(
        polynomial
    )

    print()
    print(
        "IMPORTANT:"
    )
    print(
        "This interpolation is diagnostic only."
    )
    print(
        "It must be independently derived before"
    )
    print(
        "being accepted as the endpoint formula."
    )


# =============================================================================
# 8. PIECEWISE L2 RECONSTRUCTION
# =============================================================================

def complete_L2_candidate():

    print()
    print("=" * 78)
    print("8. PIECEWISE L2 RECONSTRUCTION DIAGNOSTIC")
    print("=" * 78)

    failures = 0

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

    for k, ell in cases:

        exact = exact_L2(
            k,
            ell,
        )

        predicted = sp.Integer(0)

        for a in range(k + 1):

            b = ell - 2 - a

            if b < 0:
                continue

            if a < k:

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

            else:

                # At this point the endpoint C(k,k)
                # is intentionally taken from the exact
                # computation rather than guessed.
                #
                # This makes this section a structural
                # reconstruction check, not a formula claim.
                _, A, B, C = quadratic_coefficients(
                    k,
                    k,
                    ell,
                )

            D_pred = sp.expand(
                A * ell**2
                + B * ell
                + C
            )

            predicted += (
                D_pred
                * N**a
                * X**b
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
        f"piecewise reconstruction failures = {failures}"
    )


# =============================================================================
# 9. FINAL DIAGNOSTIC
# =============================================================================

def final_diagnostic():

    print()
    print("=" * 78)
    print("9. FINAL DIAGNOSTIC")
    print("=" * 78)

    print()
    print(
        "The exact L2 coefficient has now separated into:"
    )

    print()
    print(
        "  INTERIOR: 0 <= a < k"
    )

    print(
        "  ENDPOINT: a = k"
    )

    print()
    print(
        "Interior laws from Experiment 213:"
    )

    print(
        "  A = -1/2*C(k+2,a)"
    )

    print(
        "  B = C(k+2,a)"
    )

    print(
        "      * (2(k+1)a + k+2)"
    )

    print(
        "      / (2(k+2))"
    )

    print(
        "  C = -k*a*(a+1)"
    )

    print(
        "      * C(k+2,a)"
    )

    print(
        "      / (2(k+2))"
    )

    print()
    print(
        "Experiment 214 must determine whether the endpoint"
    )

    print(
        "a=k has a clean independently derivable correction."
    )

    print()
    print(
        "Only after that endpoint is understood should we"
    )

    print(
        "write the complete closed form for L2 and proceed"
    )

    print(
        "to the L3 peel."
    )

    print()
    print(
        "The crucial distinction remains:"
    )

    print(
        "  exact pq kernel -> homogeneous layers -> formula"
    )

    print(
        "not"
    )

    print(
        "  theorem anchors -> fitted coefficient formula"
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 214")
    print("EXACT L2 ENDPOINT / BOUNDARY-CORRECTION AUDIT")
    print("=" * 78)

    print()
    print(
        "No previous experiment output is read."
    )

    print(
        "No unrestricted affine candidate enumeration is performed."
    )

    interior_audit()
    endpoint_B_audit()
    endpoint_C_table()
    endpoint_C_ell_audit()
    endpoint_C_normalization()
    endpoint_C_difference()
    endpoint_C_interpolation()
    complete_L2_candidate()
    final_diagnostic()

    print()
    print("=" * 78)
    print("END EXPERIMENT 214")
    print("=" * 78)


if __name__ == "__main__":
    main()

