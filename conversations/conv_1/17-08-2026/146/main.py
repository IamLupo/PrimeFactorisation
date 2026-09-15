import sympy as sp

p, q = sp.symbols("p q")
N, S, X = sp.symbols("N S X")


# =============================================================================
# EXPERIMENT 211
# EXACT SECOND-LAYER STRUCTURE AUDIT
# =============================================================================
#
# Goal:
#
#   1. Compute L2 = hom_(ell-2)(G - G_top - L1)
#      only for bounded cases.
#
#   2. Record the exact coefficient sequence
#
#        [N^a X^(ell-2-a)] L2.
#
#   3. Compare the sequence across different ell for fixed k.
#
#   4. Test whether each coefficient has a simple binomial dependence
#      on ell and k.
#
# No unrestricted search.
# No large ell.
# No theorem-anchor fitting.
#
# The purpose is DERIVATION:
#
#     top layer
#        ->
#     L1 exact closed form
#        ->
#     identify L2 exact closed form.
# =============================================================================


# =============================================================================
# SAFE POLYNOMIAL UTILITIES
# =============================================================================

def poly_NX(expr):
    expr = sp.cancel(sp.together(sp.expand(expr)))

    denom = sp.denom(expr)

    if sp.expand(denom - 1) != 0:
        raise ValueError(
            f"Non-polynomial expression: denominator={denom}"
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
            f"Symmetrization remainder nonzero "
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
            f"Could not identify symmetric variables: {mapping}"
        )

    return sp.expand(
        result.subs(
            {
                sym_sum: S,
                sym_prod: N,
            }
        )
    )


def exact_GX(k, ell):
    G = exact_G(k, ell)

    expr = sp.expand(
        G.subs(
            S,
            X - 1,
        )
    )

    return sp.expand(expr)


# =============================================================================
# PROVED TOP LAYER FROM EXPERIMENT 209R
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
# PROVED L1 FORMULA FROM EXPERIMENT 210R
# =============================================================================

def L1_coefficient(k, ell, a):

    if a < 0 or a > k:
        return sp.Integer(0)

    if a < k:
        return sp.expand(
            ell * sp.binomial(k + 1, a)
            - k * sp.binomial(k, a - 1)
        )

    return sp.expand(
        ell * (k + 1)
        - k * (k - 1)
    )


def L1(k, ell):

    degree = ell - 1

    out = sp.Integer(0)

    for a in range(k + 1):

        b = degree - a

        if b < 0:
            continue

        c = L1_coefficient(
            k,
            ell,
            a,
        )

        out += c * N**a * X**b

    return sp.expand(out)


# =============================================================================
# EXACT L2
# =============================================================================

def L2(k, ell):

    # IMPORTANT:
    # Keep ell modest so SymPy never performs the expensive large
    # expansions that caused the previous run to become slow.
    if ell > 24:
        raise ValueError(
            f"Experiment 211 intentionally limits ell <= 24; "
            f"got ell={ell}"
        )

    G = exact_GX(
        k,
        ell,
    )

    top = top_layer(
        k,
        ell,
    )

    first = L1(
        k,
        ell,
    )

    residual = sp.expand(
        G - top - first
    )

    return homogeneous(
        residual,
        ell - 2,
    )


# =============================================================================
# 1. EXACT L2 SEQUENCES
# =============================================================================

def sequence_audit():

    print()
    print("=" * 78)
    print("1. EXACT L2 COEFFICIENT SEQUENCES")
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
        (9, 23),

        (11, 23),
    ]

    for k, ell in cases:

        layer = L2(
            k,
            ell,
        )

        degree = ell - 2

        print()
        print(
            f"k={k} ell={ell}"
        )

        values = []

        for a in range(
            k + 2
        ):

            b = degree - a

            if b < 0:
                continue

            value = coefficient(
                layer,
                a,
                b,
            )

            if value != 0:
                values.append(
                    (a, value)
                )

        print(
            "  "
            + " ".join(
                f"a={a}:{v}"
                for a, v in values
            )
        )


# =============================================================================
# 2. SUPPORT AUDIT
# =============================================================================

def support_audit():

    print()
    print("=" * 78)
    print("2. L2 SUPPORT AUDIT")
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

        layer = L2(
            k,
            ell,
        )

        degree = ell - 2

        support = []

        for a in range(
            degree + 1
        ):

            b = degree - a

            if coefficient(
                layer,
                a,
                b,
            ) != 0:

                support.append(a)

        print(
            f"k={k:2d} ell={ell:2d} "
            f"support={support}"
        )


# =============================================================================
# 3. FIXED-k / VARY-ell COMPARISON
# =============================================================================

def fixed_k_comparison():

    print()
    print("=" * 78)
    print("3. FIXED-k / VARY-ell COMPARISON")
    print("=" * 78)

    groups = {
        3: [7, 9, 11, 13],
        5: [11, 13, 15, 19],
        7: [15, 17],
    }

    for k, ells in groups.items():

        print()
        print(
            f"k={k}"
        )

        max_a = k + 1

        for a in range(
            max_a + 1
        ):

            row = []

            for ell in ells:

                layer = L2(
                    k,
                    ell,
                )

                b = ell - 2 - a

                value = coefficient(
                    layer,
                    a,
                    b,
                )

                row.append(
                    value
                )

            print(
                f"  a={a}: {row}"
            )


# =============================================================================
# 4. NORMALIZED BINOMIAL AUDIT
# =============================================================================
#
# The raw sequences strongly suggest that L2 is not arbitrary.
#
# We inspect ratios against natural binomial factors without searching
# an unrestricted family.
#
# =============================================================================

def normalized_audit():

    print()
    print("=" * 78)
    print("4. NORMALIZED L2 AUDIT")
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

        layer = L2(
            k,
            ell,
        )

        degree = ell - 2

        print()
        print(
            f"k={k} ell={ell}"
        )

        for a in range(
            k + 1
        ):

            b = degree - a

            value = coefficient(
                layer,
                a,
                b,
            )

            if value == 0:
                continue

            candidates = []

            cka = sp.binomial(k, a)
            ell_choose_2 = sp.binomial(ell, 2)
            ell_choose_a = sp.binomial(ell, a + 2)

            if cka != 0:
                candidates.append(
                    (
                        "D/C(k,a)",
                        sp.simplify(
                            value / cka
                        ),
                    )
                )

            if ell_choose_2 != 0:
                candidates.append(
                    (
                        "D/C(ell,2)",
                        sp.simplify(
                            value / ell_choose_2
                        ),
                    )
                )

            if ell_choose_a != 0:
                candidates.append(
                    (
                        "D/C(ell,a+2)",
                        sp.simplify(
                            value / ell_choose_a
                        ),
                    )
                )

            print(
                f"  a={a}: D={value} "
                + " ".join(
                    f"{name}={ratio}"
                    for name, ratio in candidates
                )
            )


# =============================================================================
# 5. SECOND-DIFFERENCE AUDIT IN ell
# =============================================================================
#
# For fixed k and a we inspect whether the dependence on ell has low
# polynomial degree. This is a structural check, not a fitting search.
#
# =============================================================================

def finite_difference(values):

    rows = [values]

    while len(rows[-1]) > 1:

        current = rows[-1]

        nxt = [
            sp.expand(
                current[i + 1]
                - current[i]
            )
            for i in range(
                len(current) - 1
            )
        ]

        rows.append(nxt)

    return rows


def ell_difference_audit():

    print()
    print("=" * 78)
    print("5. FINITE-DIFFERENCE AUDIT IN ell")
    print("=" * 78)

    groups = {
        3: [7, 9, 11, 13],
        5: [11, 13, 15, 17, 19],
        7: [15, 17, 19, 21],
    }

    for k, ells in groups.items():

        print()
        print(
            f"k={k}"
        )

        for a in range(
            k + 1
        ):

            values = []

            for ell in ells:

                layer = L2(
                    k,
                    ell,
                )

                b = ell - 2 - a

                values.append(
                    coefficient(
                        layer,
                        a,
                        b,
                    )
                )

            rows = finite_difference(
                values
            )

            print(
                f"  a={a}: values={values}"
            )

            for order, row in enumerate(
                rows[1:],
                start=1,
            ):

                if len(set(row)) == 1:

                    print(
                        f"       constant "
                        f"difference order={order}: "
                        f"{row[0]}"
                    )

                    break


# =============================================================================
# 6. DIRECT RECONSTRUCTION AUDIT
# =============================================================================

def reconstruction_audit():

    print()
    print("=" * 78)
    print("6. TOP + L1 + L2 RECONSTRUCTION")
    print("=" * 78)

    cases = [
        (3, 7),
        (3, 9),
        (3, 11),

        (5, 11),
        (5, 13),
        (5, 19),

        (7, 15),
        (9, 21),
        (11, 23),
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

        first = L1(
            k,
            ell,
        )

        second = L2(
            k,
            ell,
        )

        exact_three_layers = sp.expand(
            homogeneous(G, ell)
            + homogeneous(G, ell - 1)
            + homogeneous(G, ell - 2)
        )

        rebuilt = sp.expand(
            top
            + first
            + second
        )

        ok = (
            sp.expand(
                rebuilt
                - exact_three_layers
            )
            == 0
        )

        print(
            f"k={k:2d} ell={ell:2d} "
            f"reconstruction="
            f"{'PASS' if ok else 'FAIL'}"
        )

        if not ok:
            failures += 1

    print()
    print(
        f"three-layer reconstruction failures = "
        f"{failures}"
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
        "Experiment 209R proved the top layer:"
    )

    print(
        "  G_top = -X^(ell-k) ((X+N)^k - N^k)"
    )

    print()
    print(
        "Experiment 210R proved the complete first residual layer:"
    )

    print(
        "  [N^a X^(ell-1-a)] L1"
    )

    print(
        "    = ell*C(k+1,a) - k*C(k,a-1)"
    )

    print()
    print(
        "Experiment 211 now asks whether L2 has an equally"
    )

    print(
        "simple binomial structure."
    )

    print()
    print(
        "The bounded computation intentionally avoids the"
    )

    print(
        "large ell cases that caused the previous run to stall."
    )

    print()
    print(
        "The desired outcome is an exact closed form for"
    )

    print(
        "  L2 = homogeneous_(ell-2)"
    )

    print(
        "       (G - G_top - L1)."
    )

    print()
    print(
        "Only after L2 is structurally derived should we"
    )

    print(
        "return to the theorem-index anchors."
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 211")
    print("EXACT SECOND-LAYER STRUCTURE AUDIT")
    print("=" * 78)

    print()
    print(
        "No previous experiment output is read."
    )

    print(
        "No unrestricted affine candidate enumeration is performed."
    )

    sequence_audit()
    support_audit()
    fixed_k_comparison()
    normalized_audit()
    ell_difference_audit()
    reconstruction_audit()
    final_diagnostic()

    print()
    print("=" * 78)
    print("END EXPERIMENT 211")
    print("=" * 78)


if __name__ == "__main__":
    main()

