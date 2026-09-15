import sympy as sp


# =============================================================================
# EXPERIMENT 210R
# CORRECTED EXACT NEXT-LAYER BINOMIAL DERIVATION AUDIT
# =============================================================================
#
# Corrections from Experiment 210:
#
#   * The degenerate ell == k boundary is excluded.
#   * Only ell > k is tested.
#   * The coefficient extractor refuses negative exponents.
#   * Every homogeneous layer is explicitly checked to be polynomial.
#   * No arbitrary candidate search is performed.
#
# Main question:
#
#   After removing
#
#       G_top = -X^(ell-k) ((X+N)^k - N^k),
#
#   does the degree-(ell-1) layer have the observed exact formula
#
#       [N^a X^(ell-1-a)] L1
#
#       = ell*C(k+1,a) - k*C(k,a-1)
#
#   for 0 <= a < k,
#
#   with the observed endpoint a=k:
#
#       ell*(k+1) - k*(k-1) ?
#
# If this survives, peel it and inspect L2.
# =============================================================================


p, q = sp.symbols("p q")
N, S, X = sp.symbols("N S X")


# =============================================================================
# THEOREM ANCHORS
# =============================================================================

ANCHORS = [
    ("A1", 1, 10, 2, 27),
    ("A2", 3, 20, 6, 935),
    ("A3", 7, 40, 15, -1797818),
    ("B1", 3, 7, 3, -3),
    ("B2", 3, 9, 4, 9),
    ("B3", 3, 11, 5, -16),
    ("B4", 5, 11, 5, -5),
    ("B5", 5, 19, 9, -196),
    ("B6", 9, 21, 11, -84),
]


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
            f"Nonzero symmetrization remainder "
            f"for k={k}, ell={ell}: {remainder}"
        )

    sym_sum = None
    sym_prod = None

    for symbol, expr in mapping:
        expr = sp.expand(expr)

        if sp.expand(expr - (p + q)) == 0:
            sym_sum = symbol

        elif sp.expand(expr - p*q) == 0:
            sym_prod = symbol

    if sym_sum is None or sym_prod is None:
        raise ValueError(
            f"Could not identify symmetric variables "
            f"for k={k}, ell={ell}: {mapping}"
        )

    return sp.expand(
        result.subs(
            {
                sym_sum: S,
                sym_prod: N,
            }
        )
    )


# =============================================================================
# SHIFT TO X = S+1
# =============================================================================

def exact_GX(k, ell):
    expr = sp.expand(
        exact_G(k, ell).subs(
            S,
            X - 1,
        )
    )

    expr = sp.cancel(
        sp.together(expr)
    )

    # We expect an ordinary polynomial in N,X.
    if expr.has(sp.Pow):

        denom = sp.denom(
            sp.together(expr)
        )

        if sp.expand(denom - 1) != 0:
            raise ValueError(
                f"Non-polynomial G(N,X-1) "
                f"for k={k}, ell={ell}: denominator={denom}"
            )

    return sp.expand(expr)


# =============================================================================
# SAFE POLYNOMIAL CONVERSION
# =============================================================================

def as_polynomial(expr):
    expr = sp.cancel(
        sp.together(
            sp.expand(expr)
        )
    )

    denom = sp.denom(expr)

    if sp.expand(denom - 1) != 0:
        raise ValueError(
            f"Expression is not polynomial in N,X: "
            f"denominator={denom}"
        )

    return sp.Poly(
        sp.expand(expr),
        N,
        X,
        domain=sp.ZZ,
    )


# =============================================================================
# SAFE COEFFICIENT EXTRACTION
# =============================================================================

def coefficient(expr, a, b):

    # Negative exponents are not polynomial coefficients.
    if a < 0 or b < 0:
        return sp.Integer(0)

    poly = as_polynomial(expr)

    return sp.expand(
        poly.coeff_monomial(
            N**a * X**b
        )
    )


# =============================================================================
# HOMOGENEOUS COMPONENT
# =============================================================================

def homogeneous_layer(expr, total_degree):

    if total_degree < 0:
        return sp.Integer(0)

    poly = as_polynomial(expr)

    result = sp.Integer(0)

    for (a, b), coeff in poly.terms():

        if a + b == total_degree:

            result += (
                coeff
                * N**a
                * X**b
            )

    return sp.expand(result)


# =============================================================================
# PROVED TOP LAYER FROM 209R
# =============================================================================

def top_layer(k, ell):

    # For the valid regime ell > k this is an ordinary polynomial.
    return sp.expand(
        -X**(ell - k)
        * (
            (X + N)**k
            - N**k
        )
    )


# =============================================================================
# RESIDUAL AFTER TOP LAYER
# =============================================================================

def residual_after_top(k, ell):

    G = exact_GX(
        k,
        ell,
    )

    top = top_layer(
        k,
        ell,
    )

    residual = sp.expand(
        G - top
    )

    # Explicit polynomial validation.
    as_polynomial(residual)

    return residual


# =============================================================================
# EXACT FIRST RESIDUAL LAYER
# =============================================================================

def exact_L1(k, ell):

    residual = residual_after_top(
        k,
        ell,
    )

    # The intended non-degenerate regime is ell > k.
    if ell <= k:
        raise ValueError(
            f"L1 experiment requires ell > k, "
            f"got k={k}, ell={ell}"
        )

    return homogeneous_layer(
        residual,
        ell - 1,
    )


# =============================================================================
# PROPOSED L1 COEFFICIENT
# =============================================================================

def proposed_L1_coefficient(k, ell, a):

    if a < 0 or a > k:
        return sp.Integer(0)

    if a < k:

        previous = (
            sp.binomial(k, a - 1)
            if a >= 1
            else 0
        )

        return sp.expand(
            ell * sp.binomial(
                k + 1,
                a,
            )
            - k * previous
        )

    # Endpoint a = k observed in the exact data.
    return sp.expand(
        ell * (k + 1)
        - k * (k - 1)
    )


# =============================================================================
# PROPOSED L1 POLYNOMIAL
# =============================================================================

def proposed_L1(k, ell):

    expr = sp.Integer(0)

    degree = ell - 1

    for a in range(k + 1):

        b = degree - a

        if b < 0:
            continue

        value = proposed_L1_coefficient(
            k,
            ell,
            a,
        )

        if value == 0:
            continue

        expr += (
            value
            * N**a
            * X**b
        )

    return sp.expand(expr)


# =============================================================================
# RECONSTRUCTION CHECK
# =============================================================================

def l1_reconstruction_check(k, ell):

    exact = exact_L1(
        k,
        ell,
    )

    candidate = proposed_L1(
        k,
        ell,
    )

    return (
        sp.expand(
            exact - candidate
        ) == 0
    )


# =============================================================================
# 1. COMPLETE COEFFICIENT AUDIT
# =============================================================================

def complete_L1_audit():

    print()
    print("=" * 78)
    print("1. COMPLETE NEXT-LAYER COEFFICIENT AUDIT")
    print("=" * 78)

    cases = []

    # IMPORTANT:
    # Start at ell = k+2.
    # The ell == k case is degenerate and is not part of this layer test.
    odd_ks = [
        1,
        3,
        5,
        7,
        9,
        11,
        13,
    ]

    for k in odd_ks:

        for ell in range(
            k + 2,
            k + 12,
        ):

            cases.append(
                (k, ell)
            )

    coefficient_failures = 0
    reconstruction_failures = 0

    for k, ell in cases:

        exact = exact_L1(
            k,
            ell,
        )

        for a in range(k + 1):

            b = ell - 1 - a

            actual = coefficient(
                exact,
                a,
                b,
            )

            expected = proposed_L1_coefficient(
                k,
                ell,
                a,
            )

            if actual != expected:

                coefficient_failures += 1

                print()
                print(
                    f"FAIL k={k} ell={ell} a={a}"
                )

                print(
                    f"  actual   = {actual}"
                )

                print(
                    f"  expected = {expected}"
                )

        ok = l1_reconstruction_check(
            k,
            ell,
        )

        if not ok:

            reconstruction_failures += 1

            print()
            print(
                f"RECONSTRUCTION FAIL "
                f"k={k} ell={ell}"
            )

            print(
                "  exact = "
                f"{exact}"
            )

            print(
                "  candidate = "
                f"{proposed_L1(k, ell)}"
            )

    print()
    print(
        f"tested (k,ell) pairs = {len(cases)}"
    )

    print(
        f"coefficient failures = "
        f"{coefficient_failures}"
    )

    print(
        f"reconstruction failures = "
        f"{reconstruction_failures}"
    )


# =============================================================================
# 2. EXACT L1 ROWS
# =============================================================================

def l1_row_audit():

    print()
    print("=" * 78)
    print("2. EXACT L1 ROWS")
    print("=" * 78)

    cases = [
        (3, 7),
        (3, 9),
        (3, 11),
        (5, 11),
        (5, 19),
        (7, 15),
        (9, 21),
        (11, 23),
    ]

    for k, ell in cases:

        L1 = exact_L1(
            k,
            ell,
        )

        degree = ell - 1

        print()
        print(
            f"k={k} ell={ell}"
        )

        for a in range(k + 1):

            b = degree - a

            if b < 0:
                continue

            value = coefficient(
                L1,
                a,
                b,
            )

            expected = proposed_L1_coefficient(
                k,
                ell,
                a,
            )

            print(
                f"  a={a:2d} "
                f"[N^{a}X^{b}]={str(value):>8} "
                f"expected={str(expected):>8}"
            )


# =============================================================================
# 3. PARITY AUDIT
# =============================================================================

def parity_audit():

    print()
    print("=" * 78)
    print("3. PARITY-UNIFIED L1 AUDIT")
    print("=" * 78)

    cases = [
        (3, 7),
        (3, 8),
        (3, 9),
        (3, 10),
        (3, 11),
        (5, 11),
        (5, 12),
        (5, 13),
        (5, 19),
        (5, 20),
        (7, 15),
        (7, 16),
        (7, 17),
        (9, 21),
        (9, 22),
        (11, 23),
        (11, 24),
    ]

    failures = 0

    for k, ell in cases:

        if ell <= k:
            continue

        parity = (
            "same"
            if (k - ell) % 2 == 0
            else "different"
        )

        ok = l1_reconstruction_check(
            k,
            ell,
        )

        print(
            f"k={k:2d} ell={ell:2d} "
            f"parity={parity:9s} "
            f"L1="
            f"{'PASS' if ok else 'FAIL'}"
        )

        if not ok:
            failures += 1

    print()
    print(
        f"parity failures = {failures}"
    )


# =============================================================================
# 4. PEEL L1
# =============================================================================

def residual_after_L1(k, ell):

    residual = residual_after_top(
        k,
        ell,
    )

    L1 = exact_L1(
        k,
        ell,
    )

    result = sp.expand(
        residual - L1
    )

    as_polynomial(result)

    return result


# =============================================================================
# 5. EXACT SECOND RESIDUAL LAYER
# =============================================================================

def exact_L2(k, ell):

    residual = residual_after_L1(
        k,
        ell,
    )

    return homogeneous_layer(
        residual,
        ell - 2,
    )


# =============================================================================
# 6. NEXT LAYER AUDIT
# =============================================================================

def next_layer_audit():

    print()
    print("=" * 78)
    print("4. NEXT DEGREE-(ell-2) LAYER")
    print("=" * 78)

    cases = [
        (3, 7),
        (3, 9),
        (3, 11),
        (5, 11),
        (5, 19),
        (7, 15),
        (9, 21),
        (11, 23),
    ]

    for k, ell in cases:

        L2 = exact_L2(
            k,
            ell,
        )

        degree = ell - 2

        print()
        print(
            f"k={k} ell={ell}"
        )

        terms = []

        for a in range(
            degree + 1
        ):

            b = degree - a

            value = coefficient(
                L2,
                a,
                b,
            )

            if value != 0:
                terms.append(
                    (a, b, value)
                )

        if not terms:

            print(
                "  L2 = ZERO"
            )

            continue

        for a, b, value in terms:

            print(
                f"  [N^{a}X^{b}] = {value}"
            )


# =============================================================================
# 7. L2 SUPPORT AUDIT
# =============================================================================

def L2_support_audit():

    print()
    print("=" * 78)
    print("5. L2 SUPPORT AUDIT")
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

    for k, ell in cases:

        L2 = exact_L2(
            k,
            ell,
        )

        degree = ell - 2

        support = []

        for a in range(
            degree + 1
        ):

            b = degree - a

            value = coefficient(
                L2,
                a,
                b,
            )

            if value != 0:
                support.append(a)

        print()
        print(
            f"k={k} ell={ell}"
        )

        print(
            f"  degree = {degree}"
        )

        print(
            f"  N-exponent support = "
            f"{support}"
        )

        if support:

            print(
                f"  lowest a = {min(support)}"
            )

            print(
                f"  highest a = {max(support)}"
            )


# =============================================================================
# 8. ANCHOR DIAGNOSTIC
# =============================================================================

def anchor_audit():

    print()
    print("=" * 78)
    print("6. THEOREM-ANCHOR DIAGNOSTIC AFTER TWO PEELS")
    print("=" * 78)

    for name, k, ell, s, expected in ANCHORS:

        # Skip any invalid layer regime.
        if ell <= k:
            print(
                f"{name}: skipped because ell <= k"
            )
            continue

        h = (k - 1) // 2
        r = s - h

        L1 = exact_L1(
            k,
            ell,
        )

        L2 = exact_L2(
            k,
            ell,
        )

        print()
        print(
            f"{name}: "
            f"k={k} ell={ell} s={s} "
            f"h={h} r={r} "
            f"expected={expected}"
        )

        # Natural post-boundary coordinate in L1.
        a1 = r - 1
        b1 = ell - 1 - a1

        v1 = coefficient(
            L1,
            a1,
            b1,
        )

        # Same N-index in L2.
        a2 = r - 1
        b2 = ell - 2 - a2

        v2 = coefficient(
            L2,
            a2,
            b2,
        )

        print(
            f"  L1 [N^{a1}X^{b1}] = {v1}"
        )

        print(
            f"  L2 [N^{a2}X^{b2}] = {v2}"
        )

        print(
            f"  expected          = {expected}"
        )


# =============================================================================
# 9. TOP-LAYER + L1 RECONSTRUCTION
# =============================================================================

def two_layer_reconstruction_audit():

    print()
    print("=" * 78)
    print("7. TOP + L1 TWO-LAYER RECONSTRUCTION")
    print("=" * 78)

    cases = [
        (1, 4),
        (1, 10),
        (3, 8),
        (3, 10),
        (3, 20),
        (5, 10),
        (5, 12),
        (5, 19),
        (7, 14),
        (7, 16),
        (9, 20),
        (11, 24),
    ]

    failures = 0

    for k, ell in cases:

        G = exact_GX(
            k,
            ell,
        )

        top = top_layer(
            k,
            ell,
        )

        L1 = exact_L1(
            k,
            ell,
        )

        degree = ell - 2

        L2 = homogeneous_layer(
            G - top - L1,
            degree,
        )

        rebuilt_partial = sp.expand(
            top
            + L1
            + L2
        )

        exact_partial = (
            homogeneous_layer(
                G,
                ell,
            )
            + homogeneous_layer(
                G,
                ell - 1,
            )
            + homogeneous_layer(
                G,
                ell - 2,
            )
        )

        ok = (
            sp.expand(
                rebuilt_partial
                - exact_partial
            )
            == 0
        )

        print(
            f"k={k:2d} ell={ell:2d} "
            f"2-layer peel="
            f"{'PASS' if ok else 'FAIL'}"
        )

        if not ok:
            failures += 1

    print()
    print(
        f"two-layer reconstruction failures = "
        f"{failures}"
    )


# =============================================================================
# 10. FINAL DIAGNOSTIC
# =============================================================================

def final_diagnostic():

    print()
    print("=" * 78)
    print("8. FINAL DIAGNOSTIC")
    print("=" * 78)

    print()
    print(
        "Experiment 209R established the exact top layer:"
    )

    print(
        "  G_top = -X^(ell-k) ((X+N)^k - N^k)"
    )

    print()
    print(
        "Experiment 210R now tests the next homogeneous layer"
    )

    print(
        "only in the non-degenerate regime ell > k."
    )

    print()
    print(
        "The proposed exact law is:"
    )

    print(
        "  [N^a X^(ell-1-a)] L1"
    )

    print(
        "    = ell*C(k+1,a) - k*C(k,a-1)"
    )

    print(
        "      for 0 <= a < k,"
    )

    print(
        "  with endpoint"
    )

    print(
        "    a=k:"
    )

    print(
        "      ell*(k+1) - k*(k-1)."
    )

    print()
    print(
        "The ell == k cases from the failed run are deliberately"
    )

    print(
        "excluded because they are degenerate for this peeling."
    )

    print()
    print(
        "If the complete L1 audit passes, the next object is:"
    )

    print()
    print(
        "  L2 = homogeneous_(ell-2)"
    )

    print(
        "       (G - G_top - L1)."
    )

    print()
    print(
        "The purpose of L2 inspection is derivation, not fitting."
    )

    print(
        "The intended sequence remains:"
    )

    print()
    print(
        "  exact pq kernel"
    )

    print(
        "      -> top homogeneous layer"
    )

    print(
        "      -> next homogeneous layer"
    )

    print(
        "      -> further peeling"
    )

    print(
        "      -> theorem coefficient."
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 210R")
    print("CORRECTED EXACT NEXT-LAYER BINOMIAL DERIVATION AUDIT")
    print("=" * 78)

    print()
    print(
        "No previous experiment output is read."
    )

    print(
        "No unrestricted affine candidate enumeration is performed."
    )

    complete_L1_audit()
    l1_row_audit()
    parity_audit()
    next_layer_audit()
    L2_support_audit()
    anchor_audit()
    two_layer_reconstruction_audit()
    final_diagnostic()

    print()
    print("=" * 78)
    print("END EXPERIMENT 210R")
    print("=" * 78)


if __name__ == "__main__":
    main()

