import sympy as sp


# =============================================================================
# EXPERIMENT 209R
# CORRECTED BOUNDARY PEELING AND NEXT-HOMOGENEOUS-LAYER AUDIT
# =============================================================================
#
# Correction from Experiment 209:
#
# The observed outer sequence is
#
#   -1, -C(k,1), -C(k,2), ..., -C(k,k-1), 0
#
# when indexed by
#
#   [N^(r-1) X^(ell-r+1)].
#
# Therefore
#
#   [N^a X^(ell-a)] G = -C(k,a)
#
# for a = 0,...,k-1.
#
# The correct top homogeneous layer is
#
#   G_top(N,X)
#       = -X^(ell-k) ((X+N)^k - N^k).
#
# This experiment:
#
#   1. verifies that identity;
#   2. subtracts the exact top layer;
#   3. extracts the degree-(ell-1) layer;
#   4. examines that next layer structurally;
#   5. compares the theorem anchors only at structurally motivated
#      near-boundary coordinates.
#
# No unrestricted affine search is performed.
# No previous experiment output is read.
# =============================================================================


p, q = sp.symbols("p q")
N, S, X = sp.symbols("N S X")


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
            f"Nonzero symmetrization remainder for "
            f"k={k}, ell={ell}: {remainder}"
        )

    first_sym = None
    second_sym = None

    for sym, expr in mapping:
        expr = sp.expand(expr)

        if sp.expand(expr - (p + q)) == 0:
            first_sym = sym

        elif sp.expand(expr - p*q) == 0:
            second_sym = sym

    if first_sym is None or second_sym is None:
        raise ValueError(
            f"Could not identify elementary symmetric variables: "
            f"{mapping}"
        )

    return sp.expand(
        result.subs(
            {
                first_sym: S,
                second_sym: N,
            }
        )
    )


# =============================================================================
# SHIFT TO X = S+1
# =============================================================================

def exact_GX(k, ell):
    G = exact_G(k, ell)
    return sp.expand(
        G.subs(S, X - 1)
    )


# =============================================================================
# HOMOGENEOUS LAYER
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
            result += coeff * N**a * X**b

    return sp.expand(result)


def coeff_NX(expr, a, b):
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
# CORRECT TOP LAYER
# =============================================================================

def proved_top_layer(k, ell):
    #
    # Exact observed law:
    #
    #   coefficient of N^a X^(ell-a)
    #       = -C(k,a)
    #
    # for a=0,...,k-1.
    #
    # Therefore
    #
    #   G_top
    #      = -sum_{a=0}^{k-1} C(k,a) N^a X^(ell-a)
    #
    #      = -X^(ell-k) ((X+N)^k - N^k).
    #

    return sp.expand(
        -X**(ell - k)
        * (
            (X + N)**k
            - N**k
        )
    )


# =============================================================================
# TOP-LAYER IDENTITY
# =============================================================================

def verify_top_layer(k, ell):
    GX = exact_GX(k, ell)

    top_exact = homogeneous_layer(
        GX,
        ell,
    )

    top_formula = proved_top_layer(
        k,
        ell,
    )

    return (
        sp.expand(top_exact - top_formula) == 0,
        top_exact,
        top_formula,
    )


# =============================================================================
# RESIDUAL AFTER PEELING
# =============================================================================

def residual_after_top(k, ell):
    GX = exact_GX(k, ell)
    top = proved_top_layer(k, ell)

    return sp.expand(
        GX - top
    )


# =============================================================================
# NEXT LAYER
# =============================================================================

def next_layer(k, ell):
    residual = residual_after_top(k, ell)

    return homogeneous_layer(
        residual,
        ell - 1,
    )


# =============================================================================
# RECONSTRUCTION CHECK
# =============================================================================

def reconstruction_check(k, ell):
    GX = exact_GX(k, ell)
    top = proved_top_layer(k, ell)
    residual = residual_after_top(k, ell)

    return (
        sp.expand(
            top + residual - GX
        ) == 0
    )


# =============================================================================
# TOP-LAYER BINOMIAL AUDIT
# =============================================================================

def top_binomial_audit():

    print()
    print("=" * 78)
    print("1. CORRECTED TOP-LAYER BINOMIAL AUDIT")
    print("=" * 78)

    cases = [
        (1, 3),
        (1, 4),
        (1, 10),
        (3, 7),
        (3, 8),
        (3, 9),
        (3, 20),
        (5, 10),
        (5, 11),
        (5, 19),
        (7, 15),
        (7, 16),
        (9, 21),
        (11, 23),
    ]

    failures = 0

    for k, ell in cases:

        ok, top_exact, top_formula = verify_top_layer(
            k,
            ell,
        )

        print(
            f"k={k:2d} ell={ell:2d} "
            f"top-layer identity = "
            f"{'PASS' if ok else 'FAIL'}"
        )

        if not ok:
            failures += 1

            print(
                f"  exact    = {top_exact}"
            )

            print(
                f"  formula  = {top_formula}"
            )

    print()
    print(
        f"top-layer identity failures = {failures}"
    )


# =============================================================================
# BOUNDARY COEFFICIENT SEQUENCE
# =============================================================================

def boundary_sequence_audit():

    print()
    print("=" * 78)
    print("2. CORRECTED OUTER-DIAGONAL SEQUENCE")
    print("=" * 78)

    cases = [
        (3, 9),
        (5, 11),
        (7, 15),
        (9, 21),
        (11, 23),
    ]

    failures = 0

    for k, ell in cases:

        print()
        print(
            f"k={k} ell={ell}"
        )

        top = proved_top_layer(k, ell)

        for r in range(1, k + 2):

            a = r - 1
            b = ell - a

            coeff = coeff_NX(
                top,
                a,
                b,
            )

            expected = (
                -sp.binomial(k, a)
                if 0 <= a <= k - 1
                else sp.Integer(0)
            )

            ok = coeff == expected

            print(
                f"  r={r:2d} "
                f"[N^{a}X^{b}]={str(coeff):>6} "
                f"-C(k,{a})={str(expected):>6} "
                f"{'PASS' if ok else 'FAIL'}"
            )

            if not ok:
                failures += 1

    print()
    print(
        f"boundary-sequence failures = {failures}"
    )


# =============================================================================
# PEELING AUDIT
# =============================================================================

def peeling_audit():

    print()
    print("=" * 78)
    print("3. EXACT BOUNDARY PEELING AUDIT")
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

        ok = reconstruction_check(
            k,
            ell,
        )

        residual = residual_after_top(
            k,
            ell,
        )

        residual_top = homogeneous_layer(
            residual,
            ell,
        )

        print()
        print(
            f"k={k:2d} ell={ell:2d} "
            f"reconstruction="
            f"{'PASS' if ok else 'FAIL'}"
        )

        print(
            "  residual degree-ell layer = "
            f"{residual_top}"
        )

        if not ok or residual_top != 0:
            failures += 1

    print()
    print(
        f"peeling failures = {failures}"
    )


# =============================================================================
# NEXT HOMOGENEOUS LAYER
# =============================================================================

def print_next_layer(k, ell):

    layer = next_layer(
        k,
        ell,
    )

    degree = ell - 1

    terms = []

    for a in range(degree + 1):

        b = degree - a

        coeff = coeff_NX(
            layer,
            a,
            b,
        )

        if coeff != 0:
            terms.append(
                (a, b, coeff)
            )

    return terms


def next_layer_audit():

    print()
    print("=" * 78)
    print("4. NEXT HOMOGENEOUS LAYER")
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

        print()
        print(
            f"k={k} ell={ell}"
        )

        terms = print_next_layer(
            k,
            ell,
        )

        if not terms:
            print(
                "  NEXT LAYER = ZERO"
            )
            continue

        for a, b, coeff in terms:

            print(
                f"  [N^{a}X^{b}] = {coeff}"
            )


# =============================================================================
# NEXT-LAYER NORMALIZATION AUDIT
# =============================================================================

def next_layer_structural_audit():

    print()
    print("=" * 78)
    print("5. STRUCTURAL NEXT-LAYER AUDIT")
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

        layer = next_layer(
            k,
            ell,
        )

        degree = ell - 1

        print()
        print(
            f"k={k} ell={ell}"
        )

        for a in range(
            0,
            min(k + 1, degree + 1),
        ):

            b = degree - a

            coeff = coeff_NX(
                layer,
                a,
                b,
            )

            if coeff == 0:
                continue

            c_a = sp.binomial(
                k,
                a,
            )

            c_prev = (
                sp.binomial(k, a - 1)
                if a >= 1
                else sp.Integer(0)
            )

            candidates = {
                "-a*C(k,a)": -a * c_a,
                "-(k-a)*C(k,a)": -(k - a) * c_a,
                "a*C(k,a)": a * c_a,
                "(k-a)*C(k,a)": (k - a) * c_a,
                "-k*C(k,a)": -k * c_a,
                "-ell*C(k,a)": -ell * c_a,
                "-C(k,a-1)": -c_prev,
                "C(k,a-1)": c_prev,
            }

            print(
                f"  a={a:2d} "
                f"[N^{a}X^{b}]={coeff}"
            )

            matches = []

            for name, value in candidates.items():
                if coeff == value:
                    matches.append(
                        f"{name}={value}"
                    )

            if matches:
                print(
                    "      "
                    + " | ".join(matches)
                )


# =============================================================================
# THEOREM ANCHOR AUDIT
# =============================================================================

def anchor_next_layer_audit():

    print()
    print("=" * 78)
    print("6. ANCHOR AUDIT AGAINST PEELED NEXT LAYER")
    print("=" * 78)

    exact_hits = 0
    abs_hits = 0

    for name, k, ell, s, expected in ANCHORS:

        h = (k - 1) // 2
        r = s - h

        degree = ell - 1

        if r < 1:
            print()
            print(
                f"{name}: invalid r={r}"
            )
            continue

        a = r - 1
        b = degree - a

        value = coeff_NX(
            next_layer(k, ell),
            a,
            b,
        )

        exact = (
            value == expected
        )

        absolute = (
            abs(value) == abs(expected)
        )

        if exact:
            exact_hits += 1

        if absolute:
            abs_hits += 1

        print()
        print(
            f"{name}: "
            f"k={k} ell={ell} s={s} "
            f"h={h} r={r}"
        )

        print(
            f"  [N^{a}X^{b}] next-layer = {value}"
        )

        print(
            f"  expected = {expected}"
        )

        print(
            f"  exact = {exact}"
        )

        print(
            f"  abs   = {absolute}"
        )

    print()
    print(
        f"next-layer exact matches = {exact_hits}/9"
    )

    print(
        f"next-layer absolute matches = {abs_hits}/9"
    )


# =============================================================================
# LEADING NEXT-LAYER COEFFICIENT AUDIT
# =============================================================================

def leading_next_layer_audit():

    print()
    print("=" * 78)
    print("7. LEADING COEFFICIENT OF THE NEXT LAYER")
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

        layer = next_layer(
            k,
            ell,
        )

        degree = ell - 1

        nonzero = []

        for a in range(
            degree + 1
        ):

            b = degree - a

            coeff = coeff_NX(
                layer,
                a,
                b,
            )

            if coeff != 0:
                nonzero.append(
                    (a, b, coeff)
                )

        print()
        print(
            f"k={k} ell={ell}"
        )

        if not nonzero:
            print(
                "  next layer = ZERO"
            )
            continue

        print(
            f"  first nonzero = {nonzero[0]}"
        )

        print(
            f"  last nonzero  = {nonzero[-1]}"
        )

        max_abs = max(
            abs(c)
            for _, _, c in nonzero
        )

        max_terms = [
            x for x in nonzero
            if abs(x[2]) == max_abs
        ]

        print(
            f"  max absolute coefficient = {max_abs}"
        )

        print(
            f"  max positions = {max_terms}"
        )


# =============================================================================
# CLOSED FORM TOP-LAYER CHECK
# =============================================================================

def closed_form_audit():

    print()
    print("=" * 78)
    print("8. CLOSED-FORM TOP-LAYER CHECK")
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

        top = proved_top_layer(
            k,
            ell,
        )

        direct = 0

        for a in range(
            0,
            k,
        ):

            direct += (
                -sp.binomial(k, a)
                * N**a
                * X**(ell - a)
            )

        direct = sp.expand(
            direct
        )

        ok = (
            sp.expand(
                top - direct
            ) == 0
        )

        print(
            f"k={k:2d} ell={ell:2d} "
            f"closed-form equality="
            f"{'PASS' if ok else 'FAIL'}"
        )

        if not ok:
            failures += 1

    print()
    print(
        f"closed-form failures = {failures}"
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
        "The previous Experiment 209 failed because the binomial "
        "index on the outer diagonal was shifted."
    )

    print()
    print(
        "The corrected exact law is:"
    )

    print(
        "    [N^a X^(ell-a)] G = -C(k,a),"
    )

    print(
        "    for a = 0,...,k-1."
    )

    print()
    print(
        "Therefore:"
    )

    print(
        "    G_top = -X^(ell-k) * ((X+N)^k - N^k)."
    )

    print()
    print(
        "This is now a proved structural layer, not a candidate."
    )

    print()
    print(
        "After removing G_top, the new mathematical object is:"
    )

    print(
        "    L1 = homogeneous_degree_(ell-1)(G - G_top)."
    )

    print()
    print(
        "The next experiment question is:"
    )

    print(
        "    Does L1 have an equally simple closed form?"
    )

    print()
    print(
        "If yes, continue peeling the homogeneous expansion."
    )

    print(
        "If not, inspect the first degree-(ell-1) coefficients "
        "directly from the pq kernel."
    )

    print()
    print(
        "The objective remains:"
    )

    print(
        "    derive -> peel -> identify,"
    )

    print(
        "not fit arbitrary coordinate transformations."
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 209R")
    print("CORRECTED BOUNDARY PEELING AND NEXT-HOMOGENEOUS-LAYER AUDIT")
    print("=" * 78)

    print()
    print(
        "No previous experiment output is read."
    )

    print(
        "No unrestricted affine candidate enumeration is performed."
    )

    top_binomial_audit()
    boundary_sequence_audit()
    peeling_audit()
    next_layer_audit()
    next_layer_structural_audit()
    anchor_next_layer_audit()
    leading_next_layer_audit()
    closed_form_audit()
    final_diagnostic()

    print()
    print("=" * 78)
    print("END EXPERIMENT 209R")
    print("=" * 78)


if __name__ == "__main__":
    main()