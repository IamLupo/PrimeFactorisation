import sympy as sp


# =============================================================================
# EXPERIMENT 212
# EXACT L2 QUADRATIC-LAW RECONSTRUCTION
# =============================================================================
#
# No previous experiment output is read.
# No unrestricted affine search is performed.
#
# Target:
#
#   D_a(k,ell) = [N^a X^(ell-2-a)] L2
#
# Experiment 211 showed that, for fixed k,a,
#
#   D_a(k,ell)
#
# is quadratic in ell, with step-2 second difference
#
#   -4*C(k+2,a).
#
# This experiment:
#
#   1. reconstructs the exact quadratic from three ell values;
#   2. validates it on additional ell values;
#   3. checks the leading coefficient;
#   4. extracts the linear coefficient;
#   5. extracts the constant term;
#   6. looks for exact binomial normalization of those coefficients.
#
# The point is derivation of L2, not theorem-anchor fitting.
# =============================================================================


p, q = sp.symbols("p q")
N, S, X = sp.symbols("N S X")


# =============================================================================
# SAFE POLYNOMIAL HELPERS
# =============================================================================

def poly_NX(expr):
    expr = sp.cancel(sp.together(sp.expand(expr)))

    if sp.expand(sp.denom(expr) - 1) != 0:
        raise ValueError(
            f"Expression is not polynomial in N,X: denominator={sp.denom(expr)}"
        )

    return sp.Poly(
        sp.expand(expr),
        N,
        X,
        domain=sp.ZZ,
    )


def coefficient(expr, a, b):
    if a < 0 or b < 0:
        return sp.Integer(0)

    poly = poly_NX(expr)

    return sp.expand(
        poly.coeff_monomial(
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
# EXACT SYMMETRIC REDUCTION
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
            f"Symmetrization remainder for k={k}, ell={ell}: {remainder}"
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
# PROVED TOP LAYER
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
# PROVED L1
#
# From Experiment 210R:
#
# [N^a X^(ell-1-a)] L1
#     = ell*C(k+1,a) - k*C(k,a-1)
#
# =============================================================================

def L1_coefficient(k, ell, a):

    if a < 0 or a > k:
        return sp.Integer(0)

    return sp.expand(
        ell * sp.binomial(k + 1, a)
        - k * sp.binomial(k, a - 1)
    )


def L1(k, ell):

    out = sp.Integer(0)

    degree = ell - 1

    for a in range(k + 1):

        b = degree - a

        if b < 0:
            continue

        out += (
            L1_coefficient(
                k,
                ell,
                a,
            )
            * N**a
            * X**b
        )

    return sp.expand(out)


# =============================================================================
# EXACT L2
# =============================================================================

def L2(k, ell):

    G = exact_GX(
        k,
        ell,
    )

    residual = sp.expand(
        G
        - top_layer(k, ell)
        - L1(k, ell)
    )

    return homogeneous(
        residual,
        ell - 2,
    )


# =============================================================================
# CACHE
# =============================================================================

L2_CACHE = {}


def get_L2(k, ell):

    key = (k, ell)

    if key not in L2_CACHE:
        L2_CACHE[key] = L2(
            k,
            ell,
        )

    return L2_CACHE[key]


def D(k, ell, a):
    """
    Exact coefficient

        D_a(k,ell)
        = [N^a X^(ell-2-a)] L2.
    """

    b = ell - 2 - a

    if b < 0:
        return sp.Integer(0)

    return coefficient(
        get_L2(k, ell),
        a,
        b,
    )


# =============================================================================
# QUADRATIC RECONSTRUCTION
# =============================================================================

def quadratic_from_three_points(
    k,
    a,
    ell0,
):
    """
    Fit D_a(k,ell) using three consecutive odd ell values

        ell0, ell0+2, ell0+4.
    """

    e = sp.symbols("e")

    points = [
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

    poly = sp.interpolate(
        points,
        e,
    )

    return sp.factor(
        sp.expand(poly)
    )


# =============================================================================
# 1. QUADRATIC LAW AUDIT
# =============================================================================

def quadratic_audit():

    print()
    print("=" * 78)
    print("1. EXACT QUADRATIC RECONSTRUCTION")
    print("=" * 78)

    cases = {
        3: [0, 1, 2, 3, 4],
        5: [0, 1, 2, 3, 4, 5, 6],
        7: [0, 1, 2, 3, 4, 5, 6, 7, 8],
        9: [0, 1, 2, 3, 4, 5],
    }

    # These are chosen so that several validation ell values remain modest.
    starts = {
        3: 7,
        5: 11,
        7: 15,
        9: 21,
    }

    validation_failures = 0

    for k, a_values in cases.items():

        ell0 = starts[k]

        print()
        print(f"k={k}")

        for a in a_values:

            if ell0 + 6 > 24:
                continue

            poly = quadratic_from_three_points(
                k,
                a,
                ell0,
            )

            print(
                f"  a={a}: "
                f"D(ell)={poly}"
            )

            tests = [
                ell0 + 6,
                ell0 + 8,
            ]

            for ell in tests:

                actual = D(
                    k,
                    ell,
                    a,
                )

                predicted = sp.expand(
                    poly.subs(
                        sp.Symbol("e"),
                        ell,
                    )
                )

                # The interpolation variable created in
                # quadratic_from_three_points may have been
                # reconstructed with a fresh symbol, so use
                # direct substitution by extracting free symbol.
                symbols = list(poly.free_symbols)

                if symbols:

                    var = symbols[0]

                    predicted = sp.expand(
                        poly.subs(
                            var,
                            ell,
                        )
                    )

                else:
                    predicted = poly

                ok = (
                    sp.expand(
                        actual
                        - predicted
                    )
                    == 0
                )

                if not ok:
                    validation_failures += 1

                print(
                    f"      ell={ell:2d}: "
                    f"actual={str(actual):>8} "
                    f"predicted={str(predicted):>8} "
                    f"{'PASS' if ok else 'FAIL'}"
                )

    print()
    print(
        f"quadratic validation failures = "
        f"{validation_failures}"
    )


# =============================================================================
# 2. LEADING COEFFICIENT AUDIT
# =============================================================================

def leading_coefficient_audit():

    print()
    print("=" * 78)
    print("2. LEADING COEFFICIENT AUDIT")
    print("=" * 78)

    cases = {
        3: [0, 1, 2, 3, 4],
        5: [0, 1, 2, 3, 4, 5, 6],
        7: [0, 1, 2, 3, 4, 5, 6, 7, 8],
    }

    failures = 0

    e = sp.symbols("e")

    for k, a_values in cases.items():

        ell0 = {
            3: 7,
            5: 11,
            7: 15,
        }[k]

        print()
        print(f"k={k}")

        for a in a_values:

            poly = quadratic_from_three_points(
                k,
                a,
                ell0,
            )

            var = list(poly.free_symbols)[0]

            lc = sp.Poly(
                poly,
                var,
            ).LC()

            expected = -sp.binomial(
                k + 2,
                a,
            )

            ok = (
                sp.expand(
                    lc - expected
                )
                == 0
            )

            if not ok:
                failures += 1

            print(
                f"  a={a:2d} "
                f"leading={str(lc):>8} "
                f"expected={str(expected):>8} "
                f"{'PASS' if ok else 'FAIL'}"
            )

    print()
    print(
        f"leading coefficient failures = {failures}"
    )


# =============================================================================
# 3. EXPLICIT FORM
# =============================================================================
#
# Write
#
#   D_a(k,ell)
#      = -C(k+2,a) * ell^2
#        + B(k,a) * ell
#        + C(k,a)
#
# and extract B,C exactly.
# =============================================================================

def extract_quadratic_coefficients():

    print()
    print("=" * 78)
    print("3. EXTRACTED QUADRATIC COEFFICIENTS")
    print("=" * 78)

    cases = {
        3: [0, 1, 2, 3, 4],
        5: [0, 1, 2, 3, 4, 5, 6],
        7: [0, 1, 2, 3, 4, 5, 6, 7, 8],
    }

    for k, a_values in cases.items():

        ell0 = {
            3: 7,
            5: 11,
            7: 15,
        }[k]

        print()
        print(f"k={k}")

        for a in a_values:

            poly = quadratic_from_three_points(
                k,
                a,
                ell0,
            )

            var = list(poly.free_symbols)[0]

            P = sp.Poly(
                poly,
                var,
            )

            A = P.coeff_monomial(
                var**2
            )

            B = P.coeff_monomial(
                var
            )

            C = P.coeff_monomial(
                1
            )

            print(
                f"  a={a:2d}: "
                f"A={sp.factor(A)} "
                f"B={sp.factor(B)} "
                f"C={sp.factor(C)}"
            )


# =============================================================================
# 4. SECOND-DIFFERENCE FORMULA AUDIT
# =============================================================================

def second_difference_audit():

    print()
    print("=" * 78)
    print("4. SECOND-DIFFERENCE FORMULA")
    print("=" * 78)

    cases = {
        3: [0, 1, 2, 3, 4],
        5: [0, 1, 2, 3, 4, 5, 6],
        7: [0, 1, 2, 3, 4, 5, 6, 7, 8],
        9: [0, 1, 2, 3, 4, 5],
    }

    starts = {
        3: 7,
        5: 11,
        7: 15,
        9: 21,
    }

    failures = 0

    for k, a_values in cases.items():

        ell0 = starts[k]

        print()
        print(f"k={k}")

        for a in a_values:

            v0 = D(
                k,
                ell0,
                a,
            )

            v1 = D(
                k,
                ell0 + 2,
                a,
            )

            v2 = D(
                k,
                ell0 + 4,
                a,
            )

            second = sp.expand(
                v2
                - 2 * v1
                + v0
            )

            expected = sp.expand(
                -4
                * sp.binomial(
                    k + 2,
                    a,
                )
            )

            ok = (
                sp.expand(
                    second
                    - expected
                )
                == 0
            )

            if not ok:
                failures += 1

            print(
                f"  a={a:2d}: "
                f"second={str(second):>8} "
                f"expected={str(expected):>8} "
                f"{'PASS' if ok else 'FAIL'}"
            )

    print()
    print(
        f"second-difference failures = {failures}"
    )


# =============================================================================
# 5. NATURAL BINOMIAL NORMALIZATIONS OF THE LINEAR TERM
# =============================================================================
#
# We do NOT search arbitrary formulas.
#
# We compare B(k,a) against a small set of structurally natural
# binomial quantities forced by the previous layers.
# =============================================================================

def linear_normalization_audit():

    print()
    print("=" * 78)
    print("5. LINEAR-TERM NORMALIZATION AUDIT")
    print("=" * 78)

    cases = {
        3: [0, 1, 2, 3, 4],
        5: [0, 1, 2, 3, 4, 5, 6],
        7: [0, 1, 2, 3, 4, 5, 6, 7, 8],
    }

    for k, a_values in cases.items():

        ell0 = {
            3: 7,
            5: 11,
            7: 15,
        }[k]

        print()
        print(f"k={k}")

        for a in a_values:

            poly = quadratic_from_three_points(
                k,
                a,
                ell0,
            )

            var = list(poly.free_symbols)[0]

            P = sp.Poly(
                poly,
                var,
            )

            B = sp.expand(
                P.coeff_monomial(
                    var
                )
            )

            candidates = []

            cb = sp.binomial(k + 2, a)

            if cb != 0:
                candidates.append(
                    (
                        "B/C(k+2,a)",
                        sp.factor(B / cb),
                    )
                )

            if a >= 1:
                c1 = sp.binomial(
                    k + 1,
                    a - 1,
                )

                if c1 != 0:
                    candidates.append(
                        (
                            "B/C(k+1,a-1)",
                            sp.factor(B / c1),
                        )
                    )

            c2 = sp.binomial(
                k + 1,
                a,
            )

            if c2 != 0:
                candidates.append(
                    (
                        "B/C(k+1,a)",
                        sp.factor(B / c2),
                    )
                )

            print(
                f"  a={a:2d}: B={B} "
                + " ".join(
                    f"{name}={value}"
                    for name, value in candidates
                )
            )


# =============================================================================
# 6. RECONSTRUCTION CHECK
# =============================================================================

def reconstruction_audit():

    print()
    print("=" * 78)
    print("6. THREE-LAYER RECONSTRUCTION")
    print("=" * 78)

    cases = [
        (3, 7),
        (3, 9),
        (3, 11),
        (3, 13),

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

        G = exact_GX(
            k,
            ell,
        )

        exact_top = homogeneous(
            G,
            ell,
        )

        exact_L1 = homogeneous(
            G - top_layer(k, ell),
            ell - 1,
        )

        exact_L2 = homogeneous(
            G
            - top_layer(k, ell)
            - exact_L1,
            ell - 2,
        )

        rebuilt = sp.expand(
            exact_top
            + exact_L1
            + exact_L2
        )

        exact_three = sp.expand(
            homogeneous(G, ell)
            + homogeneous(G, ell - 1)
            + homogeneous(G, ell - 2)
        )

        ok = (
            sp.expand(
                rebuilt
                - exact_three
            )
            == 0
        )

        print(
            f"k={k:2d} ell={ell:2d} "
            f"{'PASS' if ok else 'FAIL'}"
        )

        if not ok:
            failures += 1

    print()
    print(
        f"reconstruction failures = {failures}"
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
        "Experiment 209R established:"
    )

    print(
        "  G_top = -X^(ell-k)((X+N)^k - N^k)"
    )

    print()
    print(
        "Experiment 210R established the complete L1 layer."
    )

    print()
    print(
        "Experiment 211 showed that L2 coefficients are quadratic"
    )

    print(
        "in ell for fixed k,a."
    )

    print()
    print(
        "The most important structural observation is:"
    )

    print(
        "  Delta_(ell,2)^2 D_a"
    )

    print(
        "      = -4*C(k+2,a)"
    )

    print()
    print(
        "Experiment 212 therefore extracts the exact quadratic"
    )

    print(
        "coefficients and asks whether the remaining linear and"
    )

    print(
        "constant terms collapse to natural binomial expressions."
    )

    print()
    print(
        "If the linear/constant terms stabilize, the next step is"
    )

    print(
        "to write an exact closed form for L2."
    )

    print()
    print(
        "Only after that should we peel L3."
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 212")
    print("EXACT L2 QUADRATIC-LAW RECONSTRUCTION")
    print("=" * 78)

    print()
    print(
        "No previous experiment output is read."
    )

    print(
        "No unrestricted affine candidate enumeration is performed."
    )

    quadratic_audit()
    leading_coefficient_audit()
    extract_quadratic_coefficients()
    second_difference_audit()
    linear_normalization_audit()
    reconstruction_audit()
    final_diagnostic()

    print()
    print("=" * 78)
    print("END EXPERIMENT 212")
    print("=" * 78)


if __name__ == "__main__":
    main()

