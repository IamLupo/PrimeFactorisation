import sympy as sp


# =============================================================================
# EXPERIMENT 216
# COMPLETE EXACT L2 CLOSED-FORM AUDIT
# =============================================================================
#
# No previous experiment output is read.
# No unrestricted candidate enumeration is performed.
#
# Established:
#
#   G_top = -X^(ell-k) ((X+N)^k - N^k)
#
#   L1_a = ell*C(k+1,a) - k*C(k,a-1)
#
# L2 support:
#
#   a = 0,...,k+1
#
# Experiment 215 found:
#
#   0 <= a < k:
#
#      D_a =
#        -1/2*C(k+2,a)*ell^2
#        + C(k+2,a)*(2(k+1)a+k+2)/(2(k+2))*ell
#        - k*a*(a+1)*C(k+2,a)/(2(k+2))
#
#   a = k:
#
#      same A and B,
#      but C receives +k(k+1)/2.
#
#   a = k+1:
#
#      D_(k+1) = D_k / 2
#
# This experiment verifies those laws directly from the exact pq kernel.
# It does NOT infer the formulas from the theorem anchors.
# =============================================================================


p, q = sp.symbols("p q")
N, S, X = sp.symbols("N S X")
ell_symbol = sp.symbols("ell")


# =============================================================================
# POLYNOMIAL HELPERS
# =============================================================================

def poly_NX(expr):
    expr = sp.cancel(sp.together(sp.expand(expr)))

    if sp.denom(expr) != 1:
        raise ValueError(
            f"Expression is not polynomial in N,X: {sp.denom(expr)}"
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

    result = sp.Integer(0)

    for (a, b), c in P.terms():
        if a + b == degree:
            result += c * N**a * X**b

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
# EXACT SYMMETRIC REDUCTION F -> G(N,S)
# =============================================================================

def exact_G(k, ell):

    F = exact_F(k, ell)

    sym_result, remainder, mapping = sp.symmetrize(
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
            f"Could not identify symmetric variables: {mapping}"
        )

    return sp.expand(
        sym_result.subs(
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


def exact_L1(k, ell):

    result = sp.Integer(0)
    degree = ell - 1

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

L2_cache = {}


def exact_L2(k, ell):

    key = (k, ell)

    if key in L2_cache:
        return L2_cache[key]

    G = exact_GX(k, ell)

    residual = sp.expand(
        G
        - top_layer(k, ell)
        - exact_L1(k, ell)
    )

    result = homogeneous(
        residual,
        ell - 2,
    )

    L2_cache[key] = result

    return result


# =============================================================================
# PROPOSED CLOSED FORM
# =============================================================================

def D_interior(k, ell, a):

    A = sp.expand(
        -sp.binomial(k + 2, a) / 2
    )

    B = sp.expand(
        sp.binomial(k + 2, a)
        * (
            2 * (k + 1) * a
            + k + 2
        )
        / (
            2 * (k + 2)
        )
    )

    C = sp.expand(
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


def D_k_closed(k, ell):

    a = k

    A = sp.expand(
        -sp.binomial(k + 2, k) / 2
    )

    B = sp.expand(
        sp.binomial(k + 2, k)
        * (
            2 * (k + 1) * k
            + k + 2
        )
        / (
            2 * (k + 2)
        )
    )

    C = sp.expand(
        -k * (k + 1)
        * (
            k * (k + 1) - 2
        )
        / 4
    )

    return sp.expand(
        A * ell**2
        + B * ell
        + C
    )


def D_k_plus_1_closed(k, ell):

    return sp.expand(
        D_k_closed(k, ell) / 2
    )


def predicted_L2(k, ell):

    result = sp.Integer(0)

    # ---------------------------------------------------------
    # a = 0,...,k-1
    # ---------------------------------------------------------

    for a in range(k):

        b = ell - 2 - a

        if b < 0:
            continue

        result += (
            D_interior(
                k,
                ell,
                a,
            )
            * N**a
            * X**b
        )

    # ---------------------------------------------------------
    # a = k
    # ---------------------------------------------------------

    if ell - 2 - k >= 0:

        result += (
            D_k_closed(
                k,
                ell,
            )
            * N**k
            * X**(ell - 2 - k)
        )

    # ---------------------------------------------------------
    # a = k+1
    # ---------------------------------------------------------

    if ell - 3 - k >= 0:

        result += (
            D_k_plus_1_closed(
                k,
                ell,
            )
            * N**(k + 1)
            * X**(ell - 3 - k)
        )

    return sp.expand(result)


# =============================================================================
# 1. DIRECT COEFFICIENT AUDIT
# =============================================================================

def direct_coefficient_audit():

    print()
    print("=" * 78)
    print("1. DIRECT CLOSED-FORM COEFFICIENT AUDIT")
    print("=" * 78)

    cases = [
        (3, 7),
        (3, 9),
        (3, 11),
        (3, 13),
        (3, 15),

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

    failures = 0
    tested = 0

    for k, ell in cases:

        exact = exact_L2(k, ell)

        print()
        print(f"k={k} ell={ell}")

        for a in range(
            0,
            min(k + 1, ell - 2) + 1,
        ):

            b = ell - 2 - a

            actual = coeff_NX(
                exact,
                a,
                b,
            )

            if a < k:

                predicted = D_interior(
                    k,
                    ell,
                    a,
                )

            elif a == k:

                predicted = D_k_closed(
                    k,
                    ell,
                )

            else:

                predicted = D_k_plus_1_closed(
                    k,
                    ell,
                )

            tested += 1

            ok = sp.expand(
                actual - predicted
            ) == 0

            if not ok:
                failures += 1

            print(
                f"  a={a:2d} "
                f"actual={str(actual):>12} "
                f"predicted={str(predicted):>12} "
                f"{'PASS' if ok else 'FAIL'}"
            )

    print()
    print(
        f"tested coefficients = {tested}"
    )

    print(
        f"coefficient failures = {failures}"
    )


# =============================================================================
# 2. a=k CORRECTION AUDIT
# =============================================================================

def endpoint_correction_audit():

    print()
    print("=" * 78)
    print("2. a=k BOUNDARY CORRECTION AUDIT")
    print("=" * 78)

    for k in [3, 5, 7, 9, 11, 13]:

        ell = k + 8

        exact = coeff_NX(
            exact_L2(k, ell),
            k,
            ell - 2 - k,
        )

        naive = D_interior(
            k,
            ell,
            k,
        )

        correction = sp.expand(
            exact - naive
        )

        predicted_correction = sp.expand(
            sp.Rational(1, 2)
            * k
            * (k + 1)
        )

        print()
        print(
            f"k={k} ell={ell}"
        )

        print(
            f"  exact       = {exact}"
        )

        print(
            f"  naive       = {naive}"
        )

        print(
            f"  correction   = {correction}"
        )

        print(
            f"  predicted    = {predicted_correction}"
        )

        print(
            "  "
            + (
                "PASS"
                if correction == predicted_correction
                else "FAIL"
            )
        )


# =============================================================================
# 3. TRUE ENDPOINT HALF-LAW AUDIT
# =============================================================================

def half_endpoint_audit():

    print()
    print("=" * 78)
    print("3. TRUE ENDPOINT D_(k+1) = D_k / 2 AUDIT")
    print("=" * 78)

    failures = 0

    cases = [
        (3, [7, 9, 11, 13, 15]),
        (5, [11, 13, 15, 17, 19]),
        (7, [15, 17, 19]),
        (9, [21, 23]),
        (11, [23, 25]),
    ]

    for k, ells in cases:

        for ell in ells:

            actual_k = coeff_NX(
                exact_L2(k, ell),
                k,
                ell - 2 - k,
            )

            actual_k1 = coeff_NX(
                exact_L2(k, ell),
                k + 1,
                ell - 3 - k,
            )

            predicted = sp.expand(
                actual_k / 2
            )

            ok = (
                sp.expand(
                    actual_k1 - predicted
                )
                == 0
            )

            if not ok:
                failures += 1

            print(
                f"k={k:2d} ell={ell:2d} "
                f"Dk={str(actual_k):>8} "
                f"Dk+1={str(actual_k1):>8} "
                f"expected={str(predicted):>8} "
                f"{'PASS' if ok else 'FAIL'}"
            )

    print()
    print(
        f"half-law failures = {failures}"
    )


# =============================================================================
# 4. SYMBOLIC RELATION BETWEEN LAST TWO ROWS
# =============================================================================

def symbolic_boundary_relation():

    print()
    print("=" * 78)
    print("4. SYMBOLIC BOUNDARY RELATION")
    print("=" * 78)

    k = sp.symbols(
        "k",
        integer=True,
        positive=True,
    )

    L = sp.symbols(
        "L",
        integer=True,
        positive=True,
    )

    Dk = sp.expand(
        -(
            k * (k + 1) / 4
        ) * L**2
        + (
            (k + 1)
            * (
                2 * (k + 1) * k
                + k + 2
            )
            / 4
        ) * L
        - (
            k
            * (k + 1)
            * (
                k * (k + 1) - 2
            )
            / 4
        )
    )

    Dkp1 = sp.expand(
        Dk / 2
    )

    print()
    print(
        "D_k(L) ="
    )

    print(
        sp.factor(Dk)
    )

    print()
    print(
        "D_(k+1)(L) ="
    )

    print(
        sp.factor(Dkp1)
    )

    print()
    print(
        "D_(k+1) - D_k/2 ="
    )

    print(
        sp.factor(
            Dkp1 - Dk / 2
        )
    )


# =============================================================================
# 5. COMPLETE SYMBOLIC L2 RECONSTRUCTION
# =============================================================================

def reconstruction_audit():

    print()
    print("=" * 78)
    print("5. COMPLETE L2 RECONSTRUCTION")
    print("=" * 78)

    cases = [
        (3, 7),
        (3, 9),
        (3, 11),
        (3, 13),

        (5, 11),
        (5, 13),
        (5, 15),
        (5, 19),

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

        predicted = predicted_L2(
            k,
            ell,
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
        f"reconstruction failures = {failures}"
    )


# =============================================================================
# 6. MANY-k GENERALIZATION AUDIT
# =============================================================================

def many_k_audit():

    print()
    print("=" * 78)
    print("6. EXTENDED k AUDIT")
    print("=" * 78)

    failures = 0
    tested = 0

    # Keep ell modest so the experiment stays fast.
    for k in [1, 3, 5, 7, 9, 11, 13, 15]:

        ell_values = [
            k + 4,
            k + 6,
            k + 8,
        ]

        for ell in ell_values:

            if ell < k + 2:
                continue

            exact = exact_L2(
                k,
                ell,
            )

            predicted = predicted_L2(
                k,
                ell,
            )

            tested += 1

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
                f"{'PASS' if ok else 'FAIL'}"
            )

    print()
    print(
        f"extended cases = {tested}"
    )

    print(
        f"extended failures = {failures}"
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
        "The L2 support is:"
    )

    print(
        "  a = 0,...,k+1"
    )

    print()
    print(
        "The proposed complete structure is:"
    )

    print(
        "  0 <= a < k:"
    )

    print(
        "    D_a ="
    )

    print(
        "      -1/2*C(k+2,a)*ell^2"
    )

    print(
        "      + C(k+2,a)"
    )

    print(
        "        * (2(k+1)a+k+2)"
    )

    print(
        "        / (2(k+2)) * ell"
    )

    print(
        "      - k*a*(a+1)*C(k+2,a)"
    )

    print(
        "        / (2(k+2))"
    )

    print()
    print(
        "  a = k:"
    )

    print(
        "    same A and B,"
    )

    print(
        "    C receives the boundary correction"
    )

    print(
        "      + k(k+1)/2"
    )

    print()
    print(
        "  a = k+1:"
    )

    print(
        "      D_(k+1) = D_k / 2"
    )

    print()
    print(
        "The decisive outcome is no longer an anchor match."
    )

    print(
        "It is whether these formulas reproduce the exact"
    )

    print(
        "pq-derived homogeneous layer identically."
    )

    print()
    print(
        "If this experiment passes, L2 is closed and"
    )

    print(
        "Experiment 217 should peel L2 and derive L3."
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
        "      -> homogeneous layers"
    )

    print(
        "      -> exact closed L2"
    )

    print(
        "      -> L3"
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 216")
    print("COMPLETE EXACT L2 CLOSED-FORM AUDIT")
    print("=" * 78)

    print()
    print(
        "No previous experiment output is read."
    )

    print(
        "No unrestricted affine candidate enumeration is performed."
    )

    direct_coefficient_audit()
    endpoint_correction_audit()
    half_endpoint_audit()
    symbolic_boundary_relation()
    reconstruction_audit()
    many_k_audit()
    final_diagnostic()

    print()
    print("=" * 78)
    print("END EXPERIMENT 216")
    print("=" * 78)


if __name__ == "__main__":
    main()
