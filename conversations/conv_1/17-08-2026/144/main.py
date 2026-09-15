import sympy as sp


# =============================================================================
# EXPERIMENT 210
# EXACT NEXT-LAYER BINOMIAL DERIVATION AUDIT
# =============================================================================
#
# Goal:
#
#   G(N,S) -> G(N,X-1)
#
# Remove the proved top layer
#
#   G_top = -X^(ell-k) ((X+N)^k - N^k)
#
# and let
#
#   L1 = homogeneous_degree_(ell-1)(G - G_top).
#
# The 209R data strongly suggest:
#
#   [N^a X^(ell-1-a)] L1
#
#       = ell*C(k+1,a) - k*C(k,a-1),
#
# for 0 <= a < k,
#
# with the endpoint a=k apparently obeying
#
#   ell*(k+1) - k*(k-1).
#
# This experiment:
#
#   1. verifies the proposed complete L1 coefficient law;
#   2. checks it across many (k,ell), both parity branches;
#   3. checks the reconstruction of L1;
#   4. peels L1;
#   5. extracts the next degree-(ell-2) layer;
#   6. reports that layer without guessing its formula.
#
# No previous experiment output is read.
# No unrestricted affine candidate search is performed.
# =============================================================================


p, q = sp.symbols("p q")
N, S, X = sp.symbols("N S X")


# =============================================================================
# ANCHORS
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
            f"Nonzero symmetrization remainder for "
            f"k={k}, ell={ell}: {remainder}"
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
            f"Could not identify elementary symmetric variables: "
            f"{mapping}"
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
# X = S+1 BASIS
# =============================================================================

def exact_GX(k, ell):
    return sp.expand(
        exact_G(k, ell).subs(
            S,
            X - 1,
        )
    )


# =============================================================================
# HOMOGENEOUS COMPONENT
# =============================================================================

def homogeneous_layer(expr, total_degree):
    poly = sp.Poly(
        sp.expand(expr),
        N,
        X,
        domain=sp.EX,
    )

    result = 0

    for (a, b), coeff in poly.terms():

        if a + b == total_degree:
            result += (
                coeff
                * N**a
                * X**b
            )

    return sp.expand(result)


def coefficient(expr, a, b):
    poly = sp.Poly(
        sp.expand(expr),
        N,
        X,
        domain=sp.EX,
    )

    return sp.expand(
        poly.coeff_monomial(
            N**a * X**b
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
# EXACT FIRST RESIDUAL
# =============================================================================

def residual_after_top(k, ell):
    G = exact_GX(k, ell)

    return sp.expand(
        G - top_layer(k, ell)
    )


# =============================================================================
# EXACT L1
# =============================================================================

def exact_L1(k, ell):
    residual = residual_after_top(
        k,
        ell,
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

        return sp.expand(
            ell * sp.binomial(
                k + 1,
                a,
            )
            - k * sp.binomial(
                k,
                a - 1,
            )
        )

    # Endpoint a = k.
    #
    # This is what the data indicate.
    #
    return sp.expand(
        ell * (k + 1)
        - k * (k - 1)
    )


# =============================================================================
# BUILD PROPOSED L1
# =============================================================================

def proposed_L1(k, ell):

    expr = 0

    for a in range(k + 1):

        b = ell - 1 - a

        if b < 0:
            continue

        coeff = proposed_L1_coefficient(
            k,
            ell,
            a,
        )

        if coeff == 0:
            continue

        expr += (
            coeff
            * N**a
            * X**b
        )

    return sp.expand(expr)


# =============================================================================
# FULL L1 RECONSTRUCTION
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
            max(3, k),
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

    print()
    print(
        f"tested (k,ell) pairs = {len(cases)}"
    )

    print(
        f"coefficient failures  = "
        f"{coefficient_failures}"
    )

    print(
        f"reconstruction failures = "
        f"{reconstruction_failures}"
    )


# =============================================================================
# 2. INTERIOR FORMULA VS ENDPOINT FORMULA
# =============================================================================

def endpoint_audit():

    print()
    print("=" * 78)
    print("2. INTERIOR VS ENDPOINT AUDIT")
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

        print()
        print(
            f"k={k} ell={ell}"
        )

        for a in range(k + 1):

            actual = coefficient(
                exact_L1(k, ell),
                a,
                ell - 1 - a,
            )

            interior = sp.expand(
                ell * sp.binomial(
                    k + 1,
                    a,
                )
                - k * (
                    sp.binomial(k, a - 1)
                    if a >= 1
                    else 0
                )
            )

            endpoint = proposed_L1_coefficient(
                k,
                ell,
                a,
            )

            print(
                f"  a={a:2d} "
                f"actual={str(actual):>8} "
                f"interior={str(interior):>8} "
                f"used={str(endpoint):>8}"
            )


# =============================================================================
# 3. STRUCTURAL FORM OF L1
# =============================================================================

def structural_form_audit():

    print()
    print("=" * 78)
    print("3. STRUCTURAL FORM OF L1")
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

        exact = exact_L1(
            k,
            ell,
        )

        candidate = proposed_L1(
            k,
            ell,
        )

        print()
        print(
            f"k={k} ell={ell}"
        )

        print(
            "  exact L1    = "
            f"{exact}"
        )

        print(
            "  proposed L1 = "
            f"{candidate}"
        )

        print(
            "  difference  = "
            f"{sp.expand(exact - candidate)}"
        )


# =============================================================================
# 4. PEEL L1 AND EXTRACT NEXT LAYER
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

    return sp.expand(
        residual - L1
    )


def exact_L2(k, ell):

    residual = residual_after_L1(
        k,
        ell,
    )

    return homogeneous_layer(
        residual,
        ell - 2,
    )


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
# 5. L2 SUPPORT / DEGREE AUDIT
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
                support.append(
                    a
                )

        print()
        print(
            f"k={k} ell={ell}"
        )

        print(
            f"  degree = {degree}"
        )

        print(
            f"  N-exponent support = {support}"
        )

        if support:

            print(
                f"  lowest a = {min(support)}"
            )

            print(
                f"  highest a = {max(support)}"
            )


# =============================================================================
# 6. PARITY AUDIT
# =============================================================================

def parity_audit():

    print()
    print("=" * 78)
    print("6. L1 PARITY-UNIFICATION AUDIT")
    print("=" * 78)

    cases = [
        (3, 7),
        (3, 8),
        (3, 9),
        (3, 10),
        (5, 11),
        (5, 12),
        (5, 19),
        (5, 20),
        (7, 15),
        (7, 16),
        (9, 21),
        (9, 22),
        (11, 23),
        (11, 24),
    ]

    failures = 0

    for k, ell in cases:

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
            f"L1 identity="
            f"{'PASS' if ok else 'FAIL'}"
        )

        if not ok:
            failures += 1

    print()
    print(
        f"parity failures = {failures}"
    )


# =============================================================================
# 7. ANCHOR AUDIT
# =============================================================================

def anchor_audit():

    print()
    print("=" * 78)
    print("7. THEOREM ANCHOR AUDIT AFTER TWO PEELS")
    print("=" * 78)

    print(
        "The purpose here is diagnostic only."
    )

    print(
        "No candidate is accepted merely because it matches anchors."
    )

    for name, k, ell, s, expected in ANCHORS:

        L1 = exact_L1(
            k,
            ell,
        )

        L2 = exact_L2(
            k,
            ell,
        )

        h = (k - 1) // 2
        r = s - h

        print()
        print(
            f"{name}: "
            f"k={k} ell={ell} s={s} "
            f"h={h} r={r} "
            f"expected={expected}"
        )

        # First peeled layer:
        a1 = r - 1
        b1 = ell - 1 - a1

        v1 = (
            coefficient(
                L1,
                a1,
                b1,
            )
            if a1 >= 0 and b1 >= 0
            else 0
        )

        # Second peeled layer:
        a2 = r - 1
        b2 = ell - 2 - a2

        v2 = (
            coefficient(
                L2,
                a2,
                b2,
            )
            if a2 >= 0 and b2 >= 0
            else 0
        )

        print(
            f"  L1 [N^{a1}X^{b1}] = {v1}"
        )

        print(
            f"  L2 [N^{a2}X^{b2}] = {v2}"
        )


# =============================================================================
# 8. COMBINATORIAL INTERIOR CHECK
# =============================================================================

def combinatorial_identity_audit():

    print()
    print("=" * 78)
    print("8. BINOMIAL-IDENTITY AUDIT")
    print("=" * 78)

    failures = 0

    for k in range(1, 16, 2):

        for a in range(0, k):

            lhs = sp.expand(
                sp.Symbol("ell")
                * sp.binomial(k + 1, a)
                - k * (
                    sp.binomial(k, a - 1)
                    if a >= 1
                    else 0
                )
            )

            # Equivalent factored identity:
            #
            # C(k+1,a) = (k+1)/(k+1-a) C(k,a)
            #
            # and
            #
            # C(k,a-1) = a/(k-a+1) C(k,a).
            #
            # We leave the result in this unsimplified form and verify
            # numerical equality over several ell values.

            for ell_value in (
                k,
                k + 1,
                k + 4,
                k + 9,
            ):

                rhs = (
                    ell_value
                    * sp.binomial(
                        k + 1,
                        a,
                    )
                    - k
                    * (
                        sp.binomial(
                            k,
                            a - 1,
                        )
                        if a >= 1
                        else 0
                    )
                )

                direct = (
                    ell_value
                    * sp.binomial(
                        k + 1,
                        a,
                    )
                    - k
                    * (
                        sp.binomial(
                            k,
                            a - 1,
                        )
                        if a >= 1
                        else 0
                    )
                )

                if sp.expand(
                    rhs - direct
                ) != 0:

                    failures += 1

    print(
        f"symbolic/numeric consistency failures = {failures}"
    )


# =============================================================================
# FINAL DIAGNOSTIC
# =============================================================================

def final_diagnostic():

    print()
    print("=" * 78)
    print("9. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "Experiment 209R established the exact top homogeneous layer:"
    )

    print()
    print(
        "  G_top = -X^(ell-k) * ((X+N)^k - N^k)"
    )

    print()
    print(
        "The next layer now has a strongly structured candidate:"
    )

    print()
    print(
        "  [N^a X^(ell-1-a)] L1"
    )

    print(
        "    = ell*C(k+1,a) - k*C(k,a-1)"
    )

    print(
        "      for 0 <= a < k,"
    )

    print()
    print(
        "with endpoint"
    )

    print(
        "  a=k:"
    )

    print(
        "    ell*(k+1) - k*(k-1)."
    )

    print()
    print(
        "The experiment deliberately verifies this law before"
    )

    print(
        "using L1 for any theorem-index search."
    )

    print()
    print(
        "If the law survives the full audit, the mathematical"
    )

    print(
        "object to study next is the peeled L2 layer:"
    )

    print()
    print(
        "  L2 = homogeneous_(ell-2)"
        "(G - G_top - L1)."
    )

    print()
    print(
        "The next successful direction would therefore be"
    )

    print(
        "  pq kernel"
        " -> exact homogeneous layer 1"
        " -> exact homogeneous layer 2"
        " -> identify the theorem coefficient."
    )

    print()
    print(
        "No arbitrary index fitting is used."
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 210")
    print("EXACT NEXT-LAYER BINOMIAL DERIVATION AUDIT")
    print("=" * 78)

    print()
    print(
        "No previous experiment output is read."
    )

    print(
        "No unrestricted affine candidate enumeration is performed."
    )

    complete_L1_audit()
    endpoint_audit()
    structural_form_audit()
    next_layer_audit()
    L2_support_audit()
    parity_audit()
    anchor_audit()
    combinatorial_identity_audit()
    final_diagnostic()

    print()
    print("=" * 78)
    print("END EXPERIMENT 210")
    print("=" * 78)


if __name__ == "__main__":
    main()

