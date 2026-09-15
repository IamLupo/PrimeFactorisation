import sympy as sp


# =============================================================================
# EXPERIMENT 217
# CORRECTED COMPLETE L2 + FIRST L3 AUDIT
# =============================================================================
#
# Purpose:
#   1. Correct the endpoint a=k+1 law exposed by Experiment 216.
#   2. Verify the complete L2 layer exactly.
#   3. Only after L2 passes, peel L2.
#   4. Record the first L3 homogeneous layer.
#
# No previous experiment output is read.
# No unrestricted candidate enumeration is performed.
# Everything is derived from the exact pq kernel.
# =============================================================================


p, q = sp.symbols("p q")
N, S, X = sp.symbols("N S X")


# =============================================================================
# BASIC POLYNOMIAL HELPERS
# =============================================================================

def poly_NX(expr):
    expr = sp.cancel(sp.together(sp.expand(expr)))

    if sp.denom(expr) != 1:
        raise ValueError(
            f"Non-polynomial expression in N,X: {sp.denom(expr)}"
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
# EXACT SYMMETRIC REDUCTION
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

        if sp.expand(replacement - (p + q)) == 0:
            sum_symbol = symbol

        elif sp.expand(replacement - p*q) == 0:
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

def L1_coeff(k, ell, a):

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
            L1_coeff(k, ell, a)
            * N**a
            * X**b
        )

    return sp.expand(out)


# =============================================================================
# EXACT L2
# =============================================================================

def exact_L2(k, ell):

    G = exact_GX(k, ell)

    residual = sp.expand(
        G
        - G_top(k, ell)
        - L1(k, ell)
    )

    return homogeneous(
        residual,
        ell - 2,
    )


# =============================================================================
# CORRECTED L2 COEFFICIENT FORMULAS
# =============================================================================

def L2_interior_coeff(k, ell, a):

    # 0 <= a < k

    A = sp.Rational(
        -1,
        2,
    ) * sp.binomial(
        k + 2,
        a,
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

    # Correct a=k boundary row.
    #
    # This is the interior formula plus
    # k(k+1)/2.

    A = (
        -sp.Rational(1, 2)
        * sp.binomial(k + 2, k)
    )

    B = (
        sp.binomial(k + 2, k)
        * (
            2 * (k + 1) * k
            + k + 2
        )
        / (
            2 * (k + 2)
        )
    )

    C = (
        -k
        * k
        * (k + 1)
        * sp.binomial(k + 2, k)
        / (
            2 * (k + 2)
        )
        + sp.Rational(
            k * (k + 1),
            2,
        )
    )

    return sp.expand(
        A * ell**2
        + B * ell
        + C
    )


def L2_k_plus_1_coeff(k, ell):

    # Correct endpoint law:
    #
    #   D_(k+1) = 2/(k+1) * D_k

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


def predicted_L2(k, ell):

    out = sp.Integer(0)

    degree = ell - 2

    # -------------------------------------------------------------------------
    # Interior: a = 0,...,k-1
    # -------------------------------------------------------------------------

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

    # -------------------------------------------------------------------------
    # Boundary a=k
    # -------------------------------------------------------------------------

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

    # -------------------------------------------------------------------------
    # True endpoint a=k+1
    # -------------------------------------------------------------------------

    a = k + 1
    b = degree - a

    if b >= 0:

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
# 1. ENDPOINT RATIO AUDIT
# =============================================================================

def endpoint_ratio_audit():

    print()
    print("=" * 78)
    print("1. CORRECTED ENDPOINT RATIO AUDIT")
    print("=" * 78)

    failures = 0

    cases = [
        (3, [7, 9, 11, 13]),
        (5, [11, 13, 15, 17, 19]),
        (7, [15, 17, 19]),
        (9, [21, 23]),
        (11, [23, 25]),
    ]

    for k, ells in cases:

        for ell in ells:

            exact = exact_L2(k, ell)

            Dk = coeff_NX(
                exact,
                k,
                ell - 2 - k,
            )

            Dkp1 = coeff_NX(
                exact,
                k + 1,
                ell - 3 - k,
            )

            predicted = sp.expand(
                sp.Rational(
                    2,
                    k + 1,
                )
                * Dk
            )

            ok = (
                sp.expand(
                    Dkp1 - predicted
                )
                == 0
            )

            if not ok:
                failures += 1

            print(
                f"k={k:2d} ell={ell:2d} "
                f"Dk={str(Dk):>8} "
                f"Dk+1={str(Dkp1):>8} "
                f"expected={str(predicted):>8} "
                f"{'PASS' if ok else 'FAIL'}"
            )

    print()
    print(
        f"endpoint ratio failures = {failures}"
    )


# =============================================================================
# 2. DIRECT COMPLETE L2 AUDIT
# =============================================================================

def complete_L2_audit():

    print()
    print("=" * 78)
    print("2. COMPLETE CORRECTED L2 AUDIT")
    print("=" * 78)

    cases = [
        (3, 7),
        (3, 9),
        (3, 11),
        (3, 13),

        (5, 11),
        (5, 13),
        (5, 15),
        (5, 17),
        (5, 19),

        (7, 15),
        (7, 17),
        (7, 19),

        (9, 21),
        (9, 23),

        (11, 23),
        (11, 25),
    ]

    coefficient_failures = 0
    reconstruction_failures = 0

    for k, ell in cases:

        exact = exact_L2(k, ell)
        predicted = predicted_L2(k, ell)

        reconstruction_ok = (
            sp.expand(
                exact - predicted
            )
            == 0
        )

        if not reconstruction_ok:
            reconstruction_failures += 1

        print()
        print(
            f"k={k} ell={ell} "
            f"reconstruction="
            f"{'PASS' if reconstruction_ok else 'FAIL'}"
        )

        degree = ell - 2

        for a in range(
            min(k + 2, degree) + 1
        ):

            b = degree - a

            actual = coeff_NX(
                exact,
                a,
                b,
            )

            if a < k:

                expected = L2_interior_coeff(
                    k,
                    ell,
                    a,
                )

            elif a == k:

                expected = L2_k_coeff(
                    k,
                    ell,
                )

            else:

                expected = L2_k_plus_1_coeff(
                    k,
                    ell,
                )

            ok = (
                sp.expand(
                    actual - expected
                )
                == 0
            )

            if not ok:
                coefficient_failures += 1

            if not ok:
                print(
                    f"  FAIL a={a}: "
                    f"actual={actual} "
                    f"expected={expected}"
                )

    print()
    print(
        f"coefficient failures = "
        f"{coefficient_failures}"
    )

    print(
        f"reconstruction failures = "
        f"{reconstruction_failures}"
    )


# =============================================================================
# 3. SYMBOLIC CORRECTED ENDPOINT FORM
# =============================================================================

def symbolic_endpoint_form():

    print()
    print("=" * 78)
    print("3. SYMBOLIC CORRECTED ENDPOINT FORM")
    print("=" * 78)

    k, ell = sp.symbols(
        "k ell",
        integer=True,
        positive=True,
    )

    Dk = sp.expand(
        -(
            (k + 1)
            * (k + 2)
            / 4
        )
        * ell**2
        + (
            (k + 1)
            * (
                2 * k**2
                + 3 * k
                + 2
            )
            / 4
        )
        * ell
        - (
            k
            * (k + 1)
            * (
                k * (k + 1)
                - 2
            )
            / 4
        )
    )

    Dkp1 = sp.expand(
        sp.Rational(
            2,
            1,
        )
        / (k + 1)
        * Dk
    )

    print()
    print("D_k(ell) =")
    print(
        sp.factor(Dk)
    )

    print()
    print("D_(k+1)(ell) =")
    print(
        sp.factor(Dkp1)
    )

    print()
    print(
        "Check:"
    )

    print(
        sp.factor(
            Dkp1
            - sp.Rational(
                2,
                1,
            )
            / (k + 1)
            * Dk
        )
    )


# =============================================================================
# 4. PEEL COMPLETE L2 AND EXTRACT L3
# =============================================================================

def L3(k, ell):

    G = exact_GX(
        k,
        ell,
    )

    residual = sp.expand(
        G
        - G_top(k, ell)
        - L1(k, ell)
        - predicted_L2(k, ell)
    )

    return homogeneous(
        residual,
        ell - 3,
    )


def L3_sequence(k, ell):

    expr = L3(
        k,
        ell,
    )

    degree = ell - 3

    seq = []

    for a in range(
        0,
        degree + 1,
    ):

        b = degree - a

        c = coeff_NX(
            expr,
            a,
            b,
        )

        if c != 0:
            seq.append(
                (a, b, c)
            )

    return seq


# =============================================================================
# 5. FIRST L3 AUDIT
# =============================================================================

def first_L3_audit():

    print()
    print("=" * 78)
    print("5. FIRST EXACT L3 AUDIT")
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

        print()
        print(
            f"k={k} ell={ell}"
        )

        seq = L3_sequence(
            k,
            ell,
        )

        print(
            "  L3 nonzero coefficients:"
        )

        for a, b, c in seq:

            print(
                f"    [N^{a}X^{b}] = {c}"
            )


# =============================================================================
# 6. L3 SUPPORT AUDIT
# =============================================================================

def L3_support_audit():

    print()
    print("=" * 78)
    print("6. L3 SUPPORT AUDIT")
    print("=" * 78)

    cases = [
        (3, 7),
        (3, 9),
        (5, 11),
        (5, 13),
        (7, 15),
        (9, 21),
        (11, 23),
    ]

    for k, ell in cases:

        expr = L3(
            k,
            ell,
        )

        degree = ell - 3

        support = []

        for a in range(
            degree + 1
        ):

            b = degree - a

            if coeff_NX(
                expr,
                a,
                b,
            ) != 0:

                support.append(a)

        print(
            f"k={k:2d} ell={ell:2d} "
            f"degree={degree:2d} "
            f"support={support}"
        )


# =============================================================================
# 7. L3 ELL-FINITE-DIFFERENCE AUDIT
# =============================================================================

def L3_ell_difference_audit():

    print()
    print("=" * 78)
    print("7. L3 ELL-FINITE-DIFFERENCE AUDIT")
    print("=" * 78)

    ks = [3, 5, 7]

    for k in ks:

        print()
        print(
            f"k={k}"
        )

        # Four odd ell values.
        ell_values = [
            k + 4,
            k + 6,
            k + 8,
            k + 10,
        ]

        # Keep only modest values.
        rows = {}

        for ell in ell_values:

            expr = L3(
                k,
                ell,
            )

            degree = ell - 3

            for a in range(
                min(k + 2, degree) + 1
            ):

                b = degree - a

                rows.setdefault(
                    a,
                    [],
                ).append(
                    coeff_NX(
                        expr,
                        a,
                        b,
                    )
                )

        for a, values in rows.items():

            if len(values) >= 4:

                d1 = [
                    values[i + 1]
                    - values[i]
                    for i in range(
                        len(values) - 1
                    )
                ]

                d2 = [
                    d1[i + 1]
                    - d1[i]
                    for i in range(
                        len(d1) - 1
                    )
                ]

                d3 = [
                    d2[i + 1]
                    - d2[i]
                    for i in range(
                        len(d2) - 1
                    )
                ]

                print(
                    f"  a={a}:"
                )

                print(
                    f"    values={values}"
                )

                print(
                    f"    Delta2={d2}"
                )

                print(
                    f"    Delta3={d3}"
                )


# =============================================================================
# 8. FINAL DIAGNOSTIC
# =============================================================================

def final_diagnostic():

    print()
    print("=" * 78)
    print("8. FINAL DIAGNOSTIC")
    print("=" * 78)

    print()
    print(
        "Experiment 216 found that the proposed"
    )

    print(
        "D_(k+1) = D_k / 2"
    )

    print(
        "was false for k >= 5."
    )

    print()
    print(
        "The exact observed law is:"
    )

    print(
        "  D_(k+1) = 2 D_k / (k+1)"
    )

    print()
    print(
        "The corrected a=k boundary coefficient is:"
    )

    print(
        "  D_k(ell)"
    )

    print(
        "    = -(k+1)(k+2)/4 * ell^2"
    )

    print(
        "      + (k+1)(2k^2+3k+2)/4 * ell"
    )

    print(
        "      - k(k+1)(k(k+1)-2)/4"
    )

    print()
    print(
        "Therefore the corrected a=k+1 coefficient is:"
    )

    print(
        "  D_(k+1)(ell)"
    )

    print(
        "    = -(k+2)/2 * ell^2"
    )

    print(
        "      + (2k^2+3k+2)/2 * ell"
    )

    print(
        "      - k(k(k+1)-2)/2"
    )

    print()
    print(
        "The purpose of this run is to establish"
    )

    print(
        "the complete L2 layer before trusting any"
    )

    print(
        "L3 pattern."
    )

    print()
    print(
        "If L2 reconstruction is exact, the"
    )

    print(
        "L3 output becomes the next genuine"
    )

    print(
        "mathematical object to analyze."
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 217")
    print("CORRECTED COMPLETE L2 + FIRST L3 AUDIT")
    print("=" * 78)

    print()
    print(
        "No previous experiment output is read."
    )

    print(
        "No unrestricted candidate enumeration is performed."
    )

    endpoint_ratio_audit()
    complete_L2_audit()
    symbolic_endpoint_form()
    first_L3_audit()
    L3_support_audit()
    L3_ell_difference_audit()
    final_diagnostic()

    print()
    print("=" * 78)
    print("END EXPERIMENT 217")
    print("=" * 78)


if __name__ == "__main__":
    main()

