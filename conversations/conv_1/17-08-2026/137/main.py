# ============================================================================
# EXPERIMENT 203
# DIRECT SYMMETRIC REDUCTION F(p,q) -> G(N,S)
# ============================================================================
#
# Purpose:
#
#   Experiment 202 established:
#
#       F(k,l;p,q) is symmetric in p,q,
#
#   but
#
#       (p+q+1) | F
#
#   only when k and l have the same parity.
#
# Therefore we now abandon the Q = F/(p+q+1) normalization completely.
#
# We construct the exact symmetric polynomial
#
#       G_{k,l}(N,S)
#
# satisfying
#
#       F(k,l;p,q) = G_{k,l}(pq,p+q).
#
# No guessed pq -> (k,l,s) map.
# No division by p+q+1.
# No unrestricted affine search.
#
# The experiment then:
#
#   1. verifies the symmetric reduction exactly;
#   2. prints compact coefficient tables;
#   3. checks coefficient symmetry/degree structure;
#   4. tests natural theorem-index coordinates against the nine anchors;
#   5. identifies which coefficient locations reproduce anchor values.
#
# ============================================================================

import sympy as sp


# ----------------------------------------------------------------------------
# SYMBOLS
# ----------------------------------------------------------------------------

p, q = sp.symbols("p q")
S, N = sp.symbols("S N")


# ----------------------------------------------------------------------------
# ANCHORS
# ----------------------------------------------------------------------------

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


# ----------------------------------------------------------------------------
# EXACT pq KERNEL
# ----------------------------------------------------------------------------

def exact_F(k, ell):
    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


# ----------------------------------------------------------------------------
# DIRECT SYMMETRIC REDUCTION
# ----------------------------------------------------------------------------

def exact_G_NS(k, ell):
    """
    Convert the exact symmetric polynomial F(p,q) into
    elementary symmetric coordinates:

        S = p + q
        N = p q

    SymPy's symmetrize performs the exact algebraic reduction.
    """

    F = exact_F(k, ell)

    result = sp.symmetrize(
        F,
        [p, q],
        formal=True
    )

    sym_poly, remainder, mapping = result

    if remainder != 0:
        raise ValueError(
            f"Unexpected non-symmetric remainder for k={k}, ell={ell}: "
            f"{remainder}"
        )

    # For two variables the formal mapping is typically:
    #
    #   s1 -> p+q
    #   s2 -> p*q
    #
    # We do not trust symbol names; use the mapping returned by SymPy.

    formal_symbols = list(sym_poly.free_symbols)

    # Determine which formal symbol corresponds to p+q and pq.
    s1 = None
    s2 = None

    for formal_symbol, original_expr in mapping:
        if sp.expand(original_expr - (p + q)) == 0:
            s1 = formal_symbol
        elif sp.expand(original_expr - p*q) == 0:
            s2 = formal_symbol

    if s1 is None or s2 is None:
        raise ValueError(
            f"Could not identify symmetric coordinates: {mapping}"
        )

    G = sp.expand(
        sym_poly.subs({
            s1: S,
            s2: N,
        })
    )

    return G


# ----------------------------------------------------------------------------
# EXACT RECONSTRUCTION CHECK
# ----------------------------------------------------------------------------

def reconstruction_check(k, ell):

    G = exact_G_NS(k, ell)

    reconstructed = sp.expand(
        G.subs({
            N: p*q,
            S: p + q,
        })
    )

    F = exact_F(k, ell)

    return sp.expand(F - reconstructed) == 0


# ----------------------------------------------------------------------------
# COEFFICIENT TABLE
# ----------------------------------------------------------------------------

def coefficient_table(G):

    poly = sp.Poly(G, N, S)

    terms = []

    for (a, b), coeff in sorted(
        poly.terms(),
        key=lambda item: (-item[0][0], -item[0][1])
    ):
        coeff = sp.expand(coeff)

        if coeff != 0:
            terms.append((a, b, coeff))

    return terms


# ----------------------------------------------------------------------------
# COMPACT TABLE PRINTER
# ----------------------------------------------------------------------------

def print_compact_table(k, ell):

    G = exact_G_NS(k, ell)

    print()
    print(f"k={k} ell={ell}")
    print(f"  total degree = {sp.Poly(G, N, S).total_degree()}")
    print(f"  degree_N     = {sp.Poly(G, N, S).degree(N)}")
    print(f"  degree_S     = {sp.Poly(G, N, S).degree(S)}")

    terms = coefficient_table(G)

    for a, b, coeff in terms:
        print(
            f"  [N^{a} S^{b}] = {coeff}"
        )


# ----------------------------------------------------------------------------
# COEFFICIENT LOOKUP
# ----------------------------------------------------------------------------

def coeff(G, a, b):

    poly = sp.Poly(G, N, S)

    return sp.expand(
        poly.coeff_monomial(N**a * S**b)
    )


# ----------------------------------------------------------------------------
# NATURAL COORDINATE FAMILIES
# ----------------------------------------------------------------------------
#
# These are deliberately small and transparent.
#
# The experiment does NOT search arbitrary affine expressions.
#
# We ask whether the anchor value occurs at obvious coefficient locations.
#
# ----------------------------------------------------------------------------

def coordinate_candidates(k, ell, s):

    return [
        ("N^s", s, 0),

        ("N^s S^1", s, 1),

        ("N^s S^2", s, 2),

        ("N^(s-1) S^1", s - 1, 1),

        ("N^(s-1) S^2", s - 1, 2),

        ("N^(s-1) S^3", s - 1, 3),

        ("N^(s-2) S^2", s - 2, 2),

        ("N^(s-2) S^3", s - 2, 3),

        ("N^(s-2) S^4", s - 2, 4),

        ("N^(s+1) S^0", s + 1, 0),

        ("N^(s+1) S^1", s + 1, 1),
    ]


# ----------------------------------------------------------------------------
# ANCHOR COORDINATE AUDIT
# ----------------------------------------------------------------------------

def anchor_coordinate_audit():

    print()
    print("=" * 78)
    print("3. NATURAL COEFFICIENT COORDINATE AUDIT")
    print("=" * 78)

    for name, k, ell, s, expected in ANCHORS:

        G = exact_G_NS(k, ell)

        print()
        print(
            f"{name}: k={k} ell={ell} s={s} "
            f"expected={expected}"
        )

        found = []

        for label, a, b in coordinate_candidates(k, ell, s):

            if a < 0 or b < 0:
                continue

            value = coeff(G, a, b)

            if value == expected:
                found.append(
                    (label, a, b, value)
                )

        if found:

            for item in found:
                label, a, b, value = item

                print(
                    f"  MATCH: {label:18s} "
                    f"[N^{a} S^{b}] = {value}"
                )

        else:

            print("  no natural-coordinate match")


# ----------------------------------------------------------------------------
# ALL NONZERO COEFFICIENTS NEAR THE TARGET REGION
# ----------------------------------------------------------------------------

def local_anchor_neighborhood():

    print()
    print("=" * 78)
    print("4. LOCAL COEFFICIENT NEIGHBORHOODS")
    print("=" * 78)

    for name, k, ell, s, expected in ANCHORS:

        G = exact_G_NS(k, ell)

        print()
        print(
            f"{name}: k={k} ell={ell} s={s} "
            f"expected={expected}"
        )

        # Inspect a small neighborhood around N-degree s.
        for a in range(max(0, s - 3), s + 3):

            row = []

            poly = sp.Poly(G, N, S)

            degree_s = poly.degree(S)

            if degree_s is None:
                continue

            for b in range(0, min(6, degree_s + 1)):

                value = coeff(G, a, b)

                if value != 0:
                    row.append(
                        f"[N^{a}S^{b}]={value}"
                    )

            if row:
                print(
                    "  " + "  ".join(row)
                )


# ----------------------------------------------------------------------------
# PARITY COMPARISON
# ----------------------------------------------------------------------------

def parity_comparison():

    print()
    print("=" * 78)
    print("5. SAME-PARITY VS OPPOSITE-PARITY G STRUCTURE")
    print("=" * 78)

    cases = [
        (1, 9),
        (1, 10),
        (3, 9),
        (3, 10),
        (5, 11),
        (5, 12),
        (7, 15),
        (7, 16),
        (9, 21),
        (9, 22),
    ]

    for k, ell in cases:

        G = exact_G_NS(k, ell)

        poly = sp.Poly(G, N, S)

        print()
        print(
            f"k={k:2d} ell={ell:2d} "
            f"parity={'same' if (k-ell)%2 == 0 else 'different'}"
        )

        print(
            f"  degree_N = {poly.degree(N)}"
        )

        print(
            f"  degree_S = {poly.degree(S)}"
        )

        # Coefficients on the lowest N-degree layer.
        min_N = min(
            a for (a, b), c in poly.terms()
        )

        lowest = []

        for (a, b), c in poly.terms():

            if a == min_N:
                lowest.append(
                    (b, c)
                )

        lowest.sort()

        print(
            f"  lowest N-degree = {min_N}"
        )

        print(
            f"  lowest-layer S coefficients = {lowest}"
        )


# ----------------------------------------------------------------------------
# TOP N-DEGREE LAYERS
# ----------------------------------------------------------------------------

def top_layers():

    print()
    print("=" * 78)
    print("6. TOP N-DEGREE LAYERS")
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

        G = exact_G_NS(k, ell)

        poly = sp.Poly(G, N, S)

        max_N = poly.degree(N)

        print()
        print(
            f"k={k} ell={ell} "
            f"max N-degree={max_N}"
        )

        for a in range(
            max_N,
            max(-1, max_N - 4),
            -1
        ):

            row = []

            for b in range(
                0,
                poly.degree(S) + 1
            ):

                value = coeff(G, a, b)

                if value != 0:
                    row.append(
                        (b, value)
                    )

            if row:
                print(
                    f"  N^{a}: {row}"
                )


# ----------------------------------------------------------------------------
# EVALUATION TESTS
# ----------------------------------------------------------------------------

def coefficient_evaluation_checks():

    print()
    print("=" * 78)
    print("7. EXACT COEFFICIENT EVALUATION CHECK")
    print("=" * 78)

    # This verifies that extracting a coefficient from G(N,S)
    # and evaluating the full polynomial are distinct operations.
    #
    # It is deliberately included to prevent accidentally comparing
    # the complete F value with a partial coefficient.

    tests = [
        (3, 7, 3),
        (3, 9, 4),
        (5, 11, 5),
        (5, 19, 9),
    ]

    for k, ell, s in tests:

        G = exact_G_NS(k, ell)

        print()
        print(
            f"k={k} ell={ell} s={s}"
        )

        print(
            f"  G(2,3) = "
            f"{G.subs({N: 2, S: 3})}"
        )

        print(
            f"  [N^{s}]G = "
            f"{coeff(G, s, 0)}"
        )

        print(
            f"  [N^{s}S^1]G = "
            f"{coeff(G, s, 1)}"
        )

        print(
            f"  [N^{s-1}S^1]G = "
            f"{coeff(G, s-1, 1)}"
        )


# ----------------------------------------------------------------------------
# FINAL DIAGNOSTIC
# ----------------------------------------------------------------------------

def final_diagnostic():

    print()
    print("=" * 78)
    print("8. FINAL DIAGNOSTIC")
    print("=" * 78)
    print()

    failures = 0

    for k in range(1, 12, 2):
        for ell in range(3, 18):

            try:
                ok = reconstruction_check(k, ell)
            except Exception:
                ok = False

            if not ok:
                failures += 1

    print(
        f"direct F -> G(N,S) reconstruction failures = {failures}"
    )

    print()

    if failures == 0:

        print(
            "The exact kernel has now been reduced directly to"
        )

        print(
            "elementary symmetric coordinates without dividing by"
        )

        print(
            "p+q+1."
        )

        print()

        print(
            "Therefore the next theorem-identification problem is:"
        )

        print(
            "identify which exact coefficient of G_{k,l}(N,S)"
        )

        print(
            "corresponds to the theorem index s."
        )

        print()

        print(
            "This keeps the full p,q dependence and does not assume"
        )

        print(
            "the parity-restricted quotient normalization."
        )

    else:

        print(
            "The symmetric reduction itself failed."
        )

        print(
            "Do not proceed to coefficient identification."
        )


# ----------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------

def main():

    print("=" * 78)
    print("EXPERIMENT 203")
    print("DIRECT SYMMETRIC REDUCTION F(p,q) -> G(N,S)")
    print("=" * 78)

    print()
    print(
        "No previous experiment output is read."
    )

    print(
        "No division by p+q+1 is used."
    )

    print(
        "No unrestricted affine candidate search is performed."
    )

    print()
    print("=" * 78)
    print("1. RECONSTRUCTION SANITY CHECK")
    print("=" * 78)

    for k, ell in [
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
        (7, 40),
        (9, 21),
        (11, 23),
    ]:

        result = reconstruction_check(k, ell)

        print(
            f"k={k:2d} ell={ell:2d} "
            f"{'PASS' if result else 'FAIL'}"
        )

    for k, ell in [
        (3, 7),
        (3, 9),
        (5, 11),
        (5, 19),
        (7, 15),
        (9, 21),
        (11, 23),
    ]:

        print_compact_table(k, ell)

    anchor_coordinate_audit()
    local_anchor_neighborhood()
    parity_comparison()
    top_layers()
    coefficient_evaluation_checks()
    final_diagnostic()

    print()
    print("=" * 78)
    print("END EXPERIMENT 203")
    print("=" * 78)


if __name__ == "__main__":
    main()

