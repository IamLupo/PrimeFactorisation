import sympy as sp

# =============================================================================
# EXPERIMENT 219
# EXACT INTERIOR L3 CLOSED-FORM VALIDATION
# =============================================================================
#
# Goal:
#
#   Validate the complete interior formula discovered from Experiment 218:
#
#   D_a(k,ell) = [N^a X^(ell-3-a)] L3
#
#   for 0 <= a < k.
#
# No arbitrary formula search.
# No theorem-anchor fitting.
# No repeated reconstruction of the same G.
#
# One exact G is constructed per (k,ell), then the homogeneous L3 layer
# is extracted once and all coefficients are checked against the formula.
#
# The endpoint a >= k is deliberately excluded.
# That is a separate boundary problem.
# =============================================================================


p, q = sp.symbols("p q")
N, S, X = sp.symbols("N S X")


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
# EXACT F -> G(N,S)
# =============================================================================

def exact_G(k, ell):

    F = exact_F(k, ell)

    sym_expr, remainder, mapping = sp.symmetrize(
        F,
        [p, q],
        formal=True,
    )

    if sp.expand(remainder) != 0:
        raise ValueError(
            f"Nonzero symmetric remainder for k={k}, ell={ell}: "
            f"{remainder}"
        )

    sum_symbol = None
    prod_symbol = None

    for symbol, replacement in mapping:

        replacement = sp.expand(replacement)

        if sp.expand(
            replacement - (p + q)
        ) == 0:
            sum_symbol = symbol

        elif sp.expand(
            replacement - p*q
        ) == 0:
            prod_symbol = symbol

    if sum_symbol is None or prod_symbol is None:
        raise ValueError(
            f"Could not identify symmetric coordinates: {mapping}"
        )

    return sp.expand(
        sym_expr.subs(
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
# POLYNOMIAL UTILITIES
# =============================================================================

def poly_NX(expr):

    expr = sp.cancel(
        sp.together(
            sp.expand(expr)
        )
    )

    den = sp.denom(expr)

    if den != 1:
        raise ValueError(
            f"Non-polynomial expression: denominator={den}"
        )

    return sp.Poly(
        sp.expand(expr),
        N,
        X,
        domain=sp.QQ,
    )


def homogeneous(expr, degree):

    P = poly_NX(expr)

    out = sp.Integer(0)

    for (a, b), coeff in P.terms():

        if a + b == degree:
            out += (
                coeff
                * N**a
                * X**b
            )

    return sp.expand(out)


def coefficient(expr, a, b):

    if a < 0 or b < 0:
        return sp.Integer(0)

    return poly_NX(expr).coeff_monomial(
        N**a * X**b
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
# FIRST LAYER
# =============================================================================

def L1(k, ell):

    degree = ell - 1

    result = sp.Integer(0)

    for a in range(k + 1):

        b = degree - a

        if b < 0:
            continue

        coeff = (
            ell * sp.binomial(k + 1, a)
            - k * sp.binomial(k, a - 1)
        )

        result += (
            coeff
            * N**a
            * X**b
        )

    return sp.expand(result)


# =============================================================================
# SECOND LAYER
# =============================================================================

def L2_interior(k, ell, a):

    A = (
        -sp.Rational(1, 2)
        * sp.binomial(k + 2, a)
    )

    B = (
        sp.binomial(k + 2, a)
        * (
            2 * (k + 1) * a
            + k + 2
        )
        / (
            2 * (k + 2)
        )
    )

    C = (
        -k
        * a
        * (a + 1)
        * sp.binomial(k + 2, a)
        / (
            2 * (k + 2)
        )
    )

    return sp.expand(
        A * ell**2
        + B * ell
        + C
    )


def L2_endpoint_k(k, ell):

    # Interior value plus the proven boundary correction.
    return sp.expand(
        L2_interior(
            k,
            ell,
            k,
        )
        + sp.Rational(
            k * (k + 1),
            2,
        )
    )


def L2_endpoint_kplus1(k, ell):

    return sp.expand(
        sp.Rational(
            2,
            k + 1,
        )
        * L2_endpoint_k(
            k,
            ell,
        )
    )


def L2(k, ell):

    degree = ell - 2

    result = sp.Integer(0)

    # -------------------------------------------------------------------------
    # Interior a = 0,...,k-1
    # -------------------------------------------------------------------------

    for a in range(k):

        b = degree - a

        if b < 0:
            continue

        result += (
            L2_interior(
                k,
                ell,
                a,
            )
            * N**a
            * X**b
        )

    # -------------------------------------------------------------------------
    # Boundary a=k
    # -------------------------------------------------------------------------

    a = k
    b = degree - a

    if b >= 0:

        result += (
            L2_endpoint_k(
                k,
                ell,
            )
            * N**a
            * X**b
        )

    # -------------------------------------------------------------------------
    # True final boundary a=k+1
    # Only include when the exact layer has this support.
    # -------------------------------------------------------------------------

    a = k + 1
    b = degree - a

    if b >= 0:

        result += (
            L2_endpoint_kplus1(
                k,
                ell,
            )
            * N**a
            * X**b
        )

    return sp.expand(result)


# =============================================================================
# EXACT L3 LAYER
# =============================================================================

def exact_L3(k, ell):

    G = exact_GX(
        k,
        ell,
    )

    residual = sp.expand(
        G
        - top_layer(k, ell)
        - L1(k, ell)
        - L2(k, ell)
    )

    return homogeneous(
        residual,
        ell - 3,
    )


# =============================================================================
# DISCOVERED INTERIOR FORMULA
# =============================================================================

def predicted_L3_coefficient(k, ell, a):

    choose = sp.binomial(
        k + 3,
        a,
    )

    A = sp.Rational(1, 6) * choose

    B = (
        -choose
        * (
            a * (k + 2)
            + k + 3
        )
        / (
            2 * (k + 3)
        )
    )

    C = (
        choose
        * (
            sp.Rational(1, 3)
            + (
                a
                * (
                    (k + 1) * a
                    + 2 * k
                    + 3
                )
                / (
                    2 * (k + 3)
                )
            )
        )
    )

    D = (
        -choose
        * k
        * a
        * (a + 1)
        * (a + 2)
        / (
            6 * (k + 3)
        )
    )

    return sp.expand(
        A * ell**3
        + B * ell**2
        + C * ell
        + D
    )


# =============================================================================
# 1. DIRECT COEFFICIENT VALIDATION
# =============================================================================

def direct_validation():

    print()
    print("=" * 78)
    print("1. DIRECT INTERIOR L3 FORMULA VALIDATION")
    print("=" * 78)

    cases = [
        # Fresh validation cases rather than only the previous table.
        (3, 7),
        (3, 9),
        (3, 11),
        (3, 13),
        (3, 15),

        (5, 11),
        (5, 13),
        (5, 15),
        (5, 17),

        (7, 15),
        (7, 17),
        (7, 19),

        (9, 21),
        (9, 23),

        (11, 23),
    ]

    failures = 0
    tested = 0

    for k, ell in cases:

        exact = exact_L3(
            k,
            ell,
        )

        degree = ell - 3

        local_failures = 0

        for a in range(k):

            b = degree - a

            if b < 0:
                continue

            actual = coefficient(
                exact,
                a,
                b,
            )

            predicted = predicted_L3_coefficient(
                k,
                ell,
                a,
            )

            tested += 1

            if sp.expand(
                actual - predicted
            ) != 0:

                failures += 1
                local_failures += 1

                print(
                    f"FAIL "
                    f"k={k} ell={ell} a={a}"
                )

                print(
                    f"  actual   = {actual}"
                )

                print(
                    f"  predicted= {predicted}"
                )

        print(
            f"k={k:2d} ell={ell:2d} "
            f"interior coefficients={k} "
            f"failures={local_failures} "
            f"{'PASS' if local_failures == 0 else 'FAIL'}"
        )

    print()
    print(
        f"tested coefficients = {tested}"
    )

    print(
        f"formula failures = {failures}"
    )


# =============================================================================
# 2. COMPONENT-BY-COMPONENT AUDIT
# =============================================================================

def component_audit():

    print()
    print("=" * 78)
    print("2. A/B/C/D COMPONENT AUDIT")
    print("=" * 78)

    # These cases are deliberately separated from the direct formula audit.
    cases = [
        (3, 13),
        (5, 17),
        (7, 19),
        (9, 23),
    ]

    e = sp.symbols(
        "e",
    )

    for k, ell in cases:

        exact = exact_L3(
            k,
            ell,
        )

        print()
        print(
            f"k={k} ell={ell}"
        )

        for a in range(k):

            # Recover the cubic polynomial using four exact ell values.
            ell_values = [
                ell,
                ell + 2,
                ell + 4,
                ell + 6,
            ]

            values = []

            for ee in ell_values:

                L3 = exact_L3(
                    k,
                    ee,
                )

                values.append(
                    coefficient(
                        L3,
                        a,
                        ee - 3 - a,
                    )
                )

            cubic = sp.interpolate(
                list(
                    zip(
                        ell_values,
                        values,
                    )
                ),
                e,
            )

            cubic = sp.expand(
                cubic
            )

            P = sp.Poly(
                cubic,
                e,
            )

            A = P.coeff_monomial(
                e**3
            )

            B = P.coeff_monomial(
                e**2
            )

            C = P.coeff_monomial(
                e
            )

            D = P.coeff_monomial(
                1
            )

            pred = sp.expand(
                predicted_L3_coefficient(
                    k,
                    e,
                    a,
                )
            )

            PP = sp.Poly(
                pred,
                e,
            )

            pA = PP.coeff_monomial(
                e**3
            )

            pB = PP.coeff_monomial(
                e**2
            )

            pC = PP.coeff_monomial(
                e
            )

            pD = PP.coeff_monomial(
                1
            )

            print(
                f"  a={a:2d} "
                f"A={A} "
                f"B={B} "
                f"C={C} "
                f"D={D}"
            )

            print(
                f"       expected "
                f"A={pA} "
                f"B={pB} "
                f"C={pC} "
                f"D={pD}"
            )

            ok = (
                sp.expand(A - pA) == 0
                and sp.expand(B - pB) == 0
                and sp.expand(C - pC) == 0
                and sp.expand(D - pD) == 0
            )

            print(
                f"       "
                f"{'PASS' if ok else 'FAIL'}"
            )


# =============================================================================
# 3. A-COEFFICIENT BINOMIAL LAW
# =============================================================================

def A_law_audit():

    print()
    print("=" * 78)
    print("3. A-COEFFICIENT BINOMIAL LAW")
    print("=" * 78)

    failures = 0

    for k in [
        3,
        5,
        7,
        9,
        11,
        13,
    ]:

        ell = k + 6

        if ell % 2 == 0:
            ell += 1

        exact = exact_L3(
            k,
            ell,
        )

        degree = ell - 3

        for a in range(k):

            b = degree - a

            actual = coefficient(
                exact,
                a,
                b,
            )

            # Extract the cubic leading coefficient
            # from three/four exact ell values through
            # the known finite-difference normalization.
            #
            # Since ell advances by 2:
            #
            #   Delta_2^3 P = 48 A
            #
            values = []

            for j in range(4):

                ee = ell + 2 * j

                L3 = exact_L3(
                    k,
                    ee,
                )

                values.append(
                    coefficient(
                        L3,
                        a,
                        ee - 3 - a,
                    )
                )

            delta1 = [
                values[i + 1] - values[i]
                for i in range(3)
            ]

            delta2 = [
                delta1[i + 1] - delta1[i]
                for i in range(2)
            ]

            delta3 = (
                delta2[1]
                - delta2[0]
            )

            A = sp.expand(
                delta3 / 48
            )

            expected = sp.Rational(
                1,
                6,
            ) * sp.binomial(
                k + 3,
                a,
            )

            ok = (
                sp.expand(
                    A - expected
                )
                == 0
            )

            if not ok:

                failures += 1

            print(
                f"k={k:2d} a={a:2d} "
                f"A={A} "
                f"expected={expected} "
                f"{'PASS' if ok else 'FAIL'}"
            )

    print()
    print(
        f"A-law failures = {failures}"
    )


# =============================================================================
# 4. ENDPOINT SEPARATION
# =============================================================================

def endpoint_separation():

    print()
    print("=" * 78)
    print("4. L3 ENDPOINT SEPARATION")
    print("=" * 78)

    cases = [
        (3, 9),
        (5, 13),
        (7, 17),
        (9, 23),
        (11, 25),
    ]

    for k, ell in cases:

        exact = exact_L3(
            k,
            ell,
        )

        degree = ell - 3

        print()
        print(
            f"k={k} ell={ell}"
        )

        for a in [
            k - 1,
            k,
            k + 1,
            k + 2,
            k + 3,
        ]:

            if a < 0:
                continue

            b = degree - a

            if b < 0:
                continue

            actual = coefficient(
                exact,
                a,
                b,
            )

            interior = predicted_L3_coefficient(
                k,
                ell,
                a,
            ) if a < k else None

            print(
                f"  a={a:2d} "
                f"[N^{a}X^{b}]={actual}"
                + (
                    f" interior-formula={interior}"
                    if interior is not None
                    else ""
                )
            )


# =============================================================================
# 5. FINAL DIAGNOSTIC
# =============================================================================

def final_diagnostic():

    print()
    print("=" * 78)
    print("5. FINAL DIAGNOSTIC")
    print("=" * 78)

    print()
    print(
        "Experiment 218 supplied enough data to identify "
        "the interior L3 cubic."
    )

    print()
    print(
        "The candidate interior law is:"
    )

    print()
    print(
        "A = C(k+3,a)/6"
    )

    print(
        "B = -C(k+3,a)"
        " * (a(k+2)+k+3)"
        " / (2(k+3))"
    )

    print(
        "C = C(k+3,a)"
        " * [1/3"
        " + a((k+1)a+2k+3)/(2(k+3))]"
    )

    print(
        "D = -k*a(a+1)(a+2)"
        " C(k+3,a)"
        " / (6(k+3))"
    )

    print()
    print(
        "This experiment tests that structure directly "
        "against fresh exact pq-derived L3 coefficients."
    )

    print()
    print(
        "The endpoint a >= k is deliberately isolated."
    )

    print(
        "It must not be mixed into the interior formula."
    )

    print()
    print(
        "If the interior audit passes, the next mathematical "
        "task is endpoint derivation, followed by L4."
    )

    print()
    print(
        "Method:"
    )

    print(
        "  exact pq kernel"
    )

    print(
        "      -> G(N,X)"
    )

    print(
        "      -> top"
    )

    print(
        "      -> L1"
    )

    print(
        "      -> L2"
    )

    print(
        "      -> L3"
    )

    print(
        "      -> exact closed interior formula"
    )

    print(
        "      -> boundary analysis"
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 219")
    print("EXACT INTERIOR L3 CLOSED-FORM VALIDATION")
    print("=" * 78)

    print()
    print(
        "No previous experiment output is read."
    )

    print(
        "No unrestricted candidate enumeration is performed."
    )

    direct_validation()
    component_audit()
    A_law_audit()
    endpoint_separation()
    final_diagnostic()

    print()
    print("=" * 78)
    print("END EXPERIMENT 219")
    print("=" * 78)


if __name__ == "__main__":
    main()

