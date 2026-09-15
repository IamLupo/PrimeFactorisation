import sympy as sp

# =============================================================================
# EXPERIMENT 218
# EXACT L3 CUBIC-LAW / BINOMIAL STRUCTURE AUDIT
# =============================================================================
#
# No previous experiment output is read.
# No unrestricted affine candidate enumeration is performed.
#
# Goal:
#
#   G -> G_top -> L1 -> L2 -> L3
#
# and determine the exact dependence of
#
#   D_a(k, ell) = [N^a X^(ell-3-a)] L3
#
# on ell.
#
# The exact pq kernel remains the source of every coefficient.
# =============================================================================


p, q = sp.symbols("p q")
N, S, X = sp.symbols("N S X")


# =============================================================================
# POLYNOMIAL HELPERS
# =============================================================================

def poly_NX(expr):
    expr = sp.cancel(sp.together(sp.expand(expr)))

    den = sp.denom(expr)

    if den != 1:
        raise ValueError(
            f"Expression is not polynomial in N,X: denominator={den}"
        )

    return sp.Poly(
        sp.expand(expr),
        N,
        X,
        domain=sp.QQ,
    )


def coeff_NX(expr, a, b):
    if a < 0 or b < 0:
        return sp.Integer(0)

    return sp.expand(
        poly_NX(expr).coeff_monomial(
            N**a * X**b
        )
    )


def homogeneous(expr, degree):
    P = poly_NX(expr)

    out = sp.Integer(0)

    for (a, b), c in P.terms():
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
# EXACT SYMMETRIC REDUCTION F -> G(N,S)
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
            f"Unexpected symmetrization remainder: {remainder}"
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
# TOP LAYER
# =============================================================================

def G_top(k, ell):

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

def L1(k, ell):

    degree = ell - 1

    out = sp.Integer(0)

    for a in range(k + 1):

        b = degree - a

        if b < 0:
            continue

        c = (
            ell * sp.binomial(k + 1, a)
            - k * sp.binomial(k, a - 1)
        )

        out += (
            c
            * N**a
            * X**b
        )

    return sp.expand(out)


# =============================================================================
# CORRECTED COMPLETE L2
# =============================================================================

def L2_interior_coeff(k, ell, a):

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


def L2_k_coeff(k, ell):

    return sp.expand(
        L2_interior_coeff(
            k,
            ell,
            k,
        )
        + sp.Rational(
            k * (k + 1),
            2,
        )
    )


def L2_k_plus_1_coeff(k, ell):

    return sp.expand(
        sp.Rational(
            2,
            k + 1,
        )
        * L2_k_coeff(
            k,
            ell,
        )
    )


def L2(k, ell):

    degree = ell - 2

    out = sp.Integer(0)

    # Interior: a = 0,...,k-1
    for a in range(k):

        b = degree - a

        if b < 0:
            continue

        out += (
            L2_interior_coeff(
                k,
                ell,
                a,
            )
            * N**a
            * X**b
        )

    # a = k
    a = k
    b = degree - a

    if b >= 0:

        out += (
            L2_k_coeff(
                k,
                ell,
            )
            * N**a
            * X**b
        )

    # a = k+1
    a = k + 1
    b = degree - a

    if b >= 0:

        # This exists only when the exact homogeneous layer
        # actually has this monomial.
        #
        # We include it as a formula, but the exact support
        # audit below decides whether the monomial is present.
        out += (
            L2_k_plus_1_coeff(
                k,
                ell,
            )
            * N**a
            * X**b
        )

    return sp.expand(out)


# =============================================================================
# EXACT L3
# =============================================================================

def exact_L3(k, ell):

    G = exact_GX(
        k,
        ell,
    )

    residual = sp.expand(
        G
        - G_top(k, ell)
        - L1(k, ell)
        - L2(k, ell)
    )

    return homogeneous(
        residual,
        ell - 3,
    )


# =============================================================================
# L3 SUPPORT
# =============================================================================

def L3_support(k, ell):

    expr = exact_L3(
        k,
        ell,
    )

    degree = ell - 3

    support = []

    for a in range(degree + 1):

        b = degree - a

        c = coeff_NX(
            expr,
            a,
            b,
        )

        if c != 0:
            support.append(a)

    return support


# =============================================================================
# L3 COEFFICIENT
# =============================================================================

def L3_coeff(k, ell, a):

    degree = ell - 3
    b = degree - a

    if b < 0:
        return sp.Integer(0)

    return coeff_NX(
        exact_L3(
            k,
            ell,
        ),
        a,
        b,
    )


# =============================================================================
# THIRD FINITE DIFFERENCE
# =============================================================================

def third_difference(values):

    d1 = [
        values[i + 1] - values[i]
        for i in range(
            len(values) - 1
        )
    ]

    d2 = [
        d1[i + 1] - d1[i]
        for i in range(
            len(d1) - 1
        )
    ]

    d3 = [
        d2[i + 1] - d2[i]
        for i in range(
            len(d2) - 1
        )
    ]

    return d3


# =============================================================================
# EXACT CUBIC FROM FOUR VALUES
# =============================================================================

def cubic_from_four_points(ell_values, values):

    if len(ell_values) != 4:
        raise ValueError(
            "Exactly four ell values are required."
        )

    e = sp.symbols(
        "e",
        integer=True,
    )

    poly = sp.interpolate(
        list(
            zip(
                ell_values,
                values,
            )
        ),
        e,
    )

    return sp.factor(
        sp.expand(poly)
    )


# =============================================================================
# EXTRACT A,B,C,D
# =============================================================================

def cubic_coefficients(poly):

    e = sp.symbols(
        "e",
        integer=True,
    )

    P = sp.Poly(
        sp.expand(poly),
        e,
    )

    return (
        sp.expand(P.coeff_monomial(e**3)),
        sp.expand(P.coeff_monomial(e**2)),
        sp.expand(P.coeff_monomial(e)),
        sp.expand(P.coeff_monomial(1)),
    )


# =============================================================================
# 1. SUPPORT AUDIT
# =============================================================================

def support_audit():

    print()
    print("=" * 78)
    print("1. EXACT L3 SUPPORT AUDIT")
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
        (9, 23),

        (11, 23),
    ]

    for k, ell in cases:

        support = L3_support(
            k,
            ell,
        )

        print(
            f"k={k:2d} ell={ell:2d} "
            f"degree={ell - 3:2d} "
            f"support={support}"
        )


# =============================================================================
# 2. EXACT L3 TABLES
# =============================================================================

def coefficient_table():

    print()
    print("=" * 78)
    print("2. EXACT L3 COEFFICIENT TABLES")
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

    for k, ell in cases:

        expr = exact_L3(
            k,
            ell,
        )

        print()
        print(
            f"k={k} ell={ell}"
        )

        degree = ell - 3

        terms = []

        for a in range(
            degree + 1
        ):

            b = degree - a

            c = coeff_NX(
                expr,
                a,
                b,
            )

            if c != 0:
                terms.append(
                    (a, b, c)
                )

        for a, b, c in terms:

            print(
                f"  [N^{a}X^{b}] = {c}"
            )


# =============================================================================
# 3. CUBIC-IN-ELL AUDIT
# =============================================================================

def cubic_audit():

    print()
    print("=" * 78)
    print("3. EXACT CUBIC-IN-ELL AUDIT")
    print("=" * 78)

    cases = [
        (3, 7, 15),
        (5, 11, 19),
        (7, 15, 21),
        (9, 21, 25),
        (11, 23, 27),
    ]

    total_failures = 0

    for k, ell_start, ell_extra in cases:

        print()
        print(
            f"k={k}"
        )

        # Use five odd ell values.
        ells = [
            ell_start,
            ell_start + 2,
            ell_start + 4,
            ell_start + 6,
            ell_start + 8,
        ]

        # Use only the first few supported a values.
        support = L3_support(
            k,
            ell_start,
        )

        max_a = min(
            max(support),
            k + 2,
        )

        for a in range(
            max_a + 1
        ):

            values = [
                L3_coeff(
                    k,
                    ell,
                    a,
                )
                for ell in ells
            ]

            # Four points determine the candidate cubic.
            cubic = cubic_from_four_points(
                ells[:4],
                values[:4],
            )

            e = sp.symbols(
                "e",
                integer=True,
            )

            predicted_extra = sp.expand(
                cubic.subs(
                    e,
                    ells[4],
                )
            )

            actual_extra = values[4]

            d3 = third_difference(
                values
            )

            # A cubic in ell has constant third
            # finite difference when ell advances by 2.
            d3_constant = (
                len(set(d3)) == 1
            )

            ok_extra = (
                sp.expand(
                    predicted_extra
                    - actual_extra
                )
                == 0
            )

            if not ok_extra:
                total_failures += 1

            print(
                f"  a={a:2d} "
                f"values={values} "
                f"Delta3={d3} "
                f"cubic-extra="
                f"{'PASS' if ok_extra else 'FAIL'} "
                f"constant-Delta3="
                f"{'PASS' if d3_constant else 'FAIL'}"
            )

    print()
    print(
        f"cubic validation failures = "
        f"{total_failures}"
    )


# =============================================================================
# 4. EXTRACT A,B,C,D
# =============================================================================

def coefficient_formula_extraction():

    print()
    print("=" * 78)
    print("4. EXTRACTED L3 CUBIC COEFFICIENTS")
    print("=" * 78)

    cases = [
        (3, 7),
        (5, 11),
        (7, 15),
        (9, 21),
        (11, 23),
    ]

    e = sp.symbols(
        "e",
        integer=True,
    )

    for k, ell0 in cases:

        print()
        print(
            f"k={k}"
        )

        ells = [
            ell0,
            ell0 + 2,
            ell0 + 4,
            ell0 + 6,
        ]

        # Only print coefficients through the actual
        # observed support of the base case.
        support = L3_support(
            k,
            ell0,
        )

        for a in support:

            if a > k + 2:
                continue

            values = [
                L3_coeff(
                    k,
                    ell,
                    a,
                )
                for ell in ells
            ]

            cubic = cubic_from_four_points(
                ells,
                values,
            )

            A, B, C, D = cubic_coefficients(
                cubic
            )

            print(
                f"  a={a:2d}: "
                f"A={A}, "
                f"B={B}, "
                f"C={C}, "
                f"D={D}"
            )


# =============================================================================
# 5. LEADING-COEFFICIENT BINOMIAL LAW
# =============================================================================

def leading_binomial_audit():

    print()
    print("=" * 78)
    print("5. L3 LEADING-COEFFICIENT BINOMIAL AUDIT")
    print("=" * 78)

    cases = [
        (3, 7),
        (5, 11),
        (7, 15),
        (9, 21),
        (11, 23),
    ]

    failures = 0

    for k, ell0 in cases:

        e = sp.symbols(
            "e",
            integer=True,
        )

        ells = [
            ell0,
            ell0 + 2,
            ell0 + 4,
            ell0 + 6,
        ]

        support = L3_support(
            k,
            ell0,
        )

        for a in support:

            if a > k + 2:
                continue

            values = [
                L3_coeff(
                    k,
                    ell,
                    a,
                )
                for ell in ells
            ]

            cubic = cubic_from_four_points(
                ells,
                values,
            )

            A, B, C, D = cubic_coefficients(
                cubic
            )

            expected = (
                sp.binomial(
                    k + 3,
                    a,
                )
                / 6
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
        f"leading-coefficient failures = "
        f"{failures}"
    )


# =============================================================================
# 6. ENDPOINT / SUPPORT AUDIT
# =============================================================================

def endpoint_audit():

    print()
    print("=" * 78)
    print("6. L3 ENDPOINT / SUPPORT AUDIT")
    print("=" * 78)

    cases = [
        (3, 7),
        (5, 11),
        (7, 15),
        (9, 21),
        (11, 23),
    ]

    for k, ell in cases:

        support = L3_support(
            k,
            ell,
        )

        print()
        print(
            f"k={k} ell={ell}"
        )

        print(
            f"  support={support}"
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

            b = ell - 3 - a

            if b < 0:
                continue

            value = coeff_NX(
                exact_L3(
                    k,
                    ell,
                ),
                a,
                b,
            )

            print(
                f"  a={a:2d} "
                f"[N^{a}X^{b}]={value}"
            )


# =============================================================================
# 7. FINAL DIAGNOSTIC
# =============================================================================

def final_diagnostic():

    print()
    print("=" * 78)
    print("7. FINAL DIAGNOSTIC")
    print("=" * 78)

    print()
    print(
        "The corrected L2 endpoint relation is:"
    )

    print(
        "  D_(k+1) = 2 D_k / (k+1)"
    )

    print()
    print(
        "The previous L2 audit reported artificial"
    )

    print(
        "coefficient failures by checking monomials"
    )

    print(
        "outside the actual homogeneous support."
    )

    print()
    print(
        "Experiment 218 now leaves L2 fixed and"
    )

    print(
        "studies only the exact next homogeneous layer:"
    )

    print()
    print(
        "  L3 = homogeneous_(ell-3)"
    )

    print(
        "       (G - G_top - L1 - L2)."
    )

    print()
    print(
        "The key structural test is whether"
    )

    print(
        "  D_a(k,ell)"
    )

    print(
        "is exactly cubic in ell."
    )

    print()
    print(
        "If that passes, the next task is to derive"
    )

    print(
        "B(k,a), C(k,a), D(k,a) from the pq kernel."
    )

    print()
    print(
        "Do not return to theorem anchors yet."
    )

    print(
        "The correct order remains:"
    )

    print(
        "  pq kernel"
    )

    print(
        "    -> G(N,X)"
    )

    print(
        "    -> top layer"
    )

    print(
        "    -> L1"
    )

    print(
        "    -> L2"
    )

    print(
        "    -> L3"
    )

    print(
        "    -> closed structural formulas"
    )

    print(
        "    -> theorem-index identification"
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 218")
    print("EXACT L3 CUBIC-LAW / BINOMIAL STRUCTURE AUDIT")
    print("=" * 78)

    print()
    print(
        "No previous experiment output is read."
    )

    print(
        "No unrestricted candidate enumeration is performed."
    )

    support_audit()
    coefficient_table()
    cubic_audit()
    coefficient_formula_extraction()
    leading_binomial_audit()
    endpoint_audit()
    final_diagnostic()

    print()
    print("=" * 78)
    print("END EXPERIMENT 218")
    print("=" * 78)


if __name__ == "__main__":
    main()

