import sympy as sp


# ============================================================================
# EXPERIMENT 204
# PARITY-UNIFIED S+1 DECOMPOSITION OF THE EXACT SYMMETRIC KERNEL
# ============================================================================
#
# Exact object:
#
#     F(k,l;p,q)
#       = p^k(1+q)^l + q^k(1+p)^l
#         - p^l(1+q)^k - q^l(1+p)^k
#
# Exact symmetric coordinates:
#
#     N = p*q
#     S = p+q
#
# Therefore
#
#     F = G(N,S)
#
# Experiment 202 showed:
#
#     S+1 divides F  <=>  k,l have the same parity.
#
# So instead of defining Q only on that branch, perform the exact
# polynomial division in S for every case:
#
#     G(N,S) = (S+1) H(N,S) + R(N)
#
# with
#
#     R(N) = G(N,-1).
#
# This gives one parity-unified representation.
#
# The experiment then:
#
#   1. verifies F -> G;
#   2. verifies G = (S+1)H + R;
#   3. checks that R=0 exactly on the same-parity branch;
#   4. inspects H and R near theorem-index locations;
#   5. tests natural coefficient coordinates against the nine anchors;
#   6. searches only a bounded neighborhood, not arbitrary affine formulas.
#
# ============================================================================


p, q = sp.symbols("p q")
N, S = sp.symbols("N S")


# ============================================================================
# KNOWN ANCHORS
# ============================================================================

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


# ============================================================================
# EXACT p,q KERNEL
# ============================================================================

def exact_F(k, ell):
    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


# ============================================================================
# EXACT SYMMETRIC REDUCTION F -> G(N,S)
# ============================================================================

def exact_G_NS(k, ell):

    F = exact_F(k, ell)

    sym_poly, remainder, mapping = sp.symmetrize(
        F,
        [p, q],
        formal=True,
    )

    if sp.expand(remainder) != 0:
        raise ValueError(
            f"Non-symmetric remainder for k={k}, ell={ell}: {remainder}"
        )

    s1 = None
    s2 = None

    for formal_symbol, original_expr in mapping:
        if sp.expand(original_expr - (p + q)) == 0:
            s1 = formal_symbol
        elif sp.expand(original_expr - p*q) == 0:
            s2 = formal_symbol

    if s1 is None or s2 is None:
        raise ValueError(
            f"Could not identify symmetric coordinates: mapping={mapping}"
        )

    G = sp.expand(
        sym_poly.subs({
            s1: S,
            s2: N,
        })
    )

    return G


# ============================================================================
# RECONSTRUCTION CHECK
# ============================================================================

def reconstruction_check(k, ell):

    G = exact_G_NS(k, ell)

    reconstructed = sp.expand(
        G.subs({
            N: p*q,
            S: p + q,
        })
    )

    return sp.expand(
        exact_F(k, ell) - reconstructed
    ) == 0


# ============================================================================
# DIVISION BY S+1
# ============================================================================

def split_S_plus_1(G):

    poly = sp.Poly(G, S, domain=sp.EX)

    quotient, remainder = sp.div(
        poly,
        sp.Poly(S + 1, S, domain=sp.EX),
    )

    H = sp.expand(quotient.as_expr())
    R = sp.expand(remainder.as_expr())

    return H, R


# ============================================================================
# VERIFY SPLIT
# ============================================================================

def split_check(k, ell):

    G = exact_G_NS(k, ell)

    H, R = split_S_plus_1(G)

    reconstructed = sp.expand(
        (S + 1) * H + R
    )

    return (
        sp.expand(G - reconstructed) == 0
        and sp.expand(R - G.subs(S, -1)) == 0
    )


# ============================================================================
# COEFFICIENT LOOKUP
# ============================================================================

def coefficient(poly_expr, n_deg, s_deg):

    poly = sp.Poly(poly_expr, N, S)

    return sp.expand(
        poly.coeff_monomial(
            N**n_deg * S**s_deg
        )
    )


# ============================================================================
# R(N) COEFFICIENT
# ============================================================================

def remainder_coefficient(R, n_deg):

    poly = sp.Poly(R, N)

    return sp.expand(
        poly.coeff_monomial(
            N**n_deg
        )
    )


# ============================================================================
# PRINT DIVISION SUMMARY
# ============================================================================

def print_division_summary(k, ell):

    G = exact_G_NS(k, ell)

    H, R = split_S_plus_1(G)

    print()
    print(f"k={k} ell={ell}")

    print(
        f"  parity = "
        f"{'same' if (k-ell) % 2 == 0 else 'different'}"
    )

    print(
        f"  H degree_N = {sp.Poly(H, N, S).degree(N)}"
    )

    print(
        f"  H degree_S = {sp.Poly(H, N, S).degree(S)}"
    )

    print(
        f"  R(N) = {sp.expand(R)}"
    )


# ============================================================================
# PARITY AUDIT
# ============================================================================

def parity_audit():

    print()
    print("=" * 78)
    print("2. PARITY-UNIFIED DECOMPOSITION AUDIT")
    print("=" * 78)

    cases = [
        (1, 3),
        (1, 4),
        (1, 9),
        (1, 10),
        (3, 7),
        (3, 8),
        (3, 9),
        (3, 10),
        (5, 11),
        (5, 12),
        (5, 19),
        (7, 15),
        (7, 16),
        (9, 21),
        (9, 22),
    ]

    failures = 0

    for k, ell in cases:

        G = exact_G_NS(k, ell)
        H, R = split_S_plus_1(G)

        expected_zero = ((k - ell) % 2 == 0)

        actual_zero = sp.expand(R) == 0

        ok = (expected_zero == actual_zero)

        if not ok:
            failures += 1

        print(
            f"k={k:2d} ell={ell:2d} "
            f"same_parity={expected_zero} "
            f"R_zero={actual_zero} "
            f"{'PASS' if ok else 'FAIL'}"
        )

    print()
    print(
        f"parity/remainder classification failures = {failures}"
    )


# ============================================================================
# PRINT R(N)
# ============================================================================

def print_remainders():

    print()
    print("=" * 78)
    print("3. REMAINDER POLYNOMIALS R(N)=G(N,-1)")
    print("=" * 78)

    cases = [
        (1, 4),
        (1, 10),
        (3, 8),
        (3, 10),
        (5, 10),
        (5, 12),
        (7, 14),
        (9, 20),
        (11, 24),
    ]

    for k, ell in cases:

        G = exact_G_NS(k, ell)
        H, R = split_S_plus_1(G)

        print()
        print(
            f"k={k} ell={ell}"
        )
        print(
            f"  R(N) = {sp.factor(R)}"
        )

        poly = sp.Poly(R, N)

        terms = []

        for (a,), c in poly.terms():

            terms.append(
                f"[N^{a}]={c}"
            )

        print(
            "  coefficients: "
            + "  ".join(terms)
        )


# ============================================================================
# H TOP N-LAYERS
# ============================================================================

def print_H_layers():

    print()
    print("=" * 78)
    print("4. TOP N-LAYERS OF H")
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
        H, R = split_S_plus_1(G)

        poly = sp.Poly(H, N, S)

        max_N = poly.degree(N)

        print()
        print(
            f"k={k} ell={ell} max_N_degree={max_N}"
        )

        # Print a bounded number of upper N-layers.
        for a in range(
            max_N,
            max(-1, max_N - 4),
            -1
        ):

            row = []

            for b in range(
                0,
                int(poly.degree(S)) + 1
            ):

                c = coefficient(H, a, b)

                if c != 0:
                    row.append(
                        (b, c)
                    )

            if row:
                print(
                    f"  N^{a}: {row}"
                )


# ============================================================================
# NATURAL ANCHOR COORDINATES
# ============================================================================
#
# We test only local coordinates suggested by the theorem index s.
#
# Both H and R are tested.
#
# ============================================================================

def natural_H_candidates(s):

    return [
        ("H[N^s S^0]", s, 0),
        ("H[N^s S^1]", s, 1),
        ("H[N^s S^2]", s, 2),

        ("H[N^(s-1) S^0]", s - 1, 0),
        ("H[N^(s-1) S^1]", s - 1, 1),
        ("H[N^(s-1) S^2]", s - 1, 2),
        ("H[N^(s-1) S^3]", s - 1, 3),

        ("H[N^(s-2) S^0]", s - 2, 0),
        ("H[N^(s-2) S^1]", s - 2, 1),
        ("H[N^(s-2) S^2]", s - 2, 2),
        ("H[N^(s-2) S^3]", s - 2, 3),
        ("H[N^(s-2) S^4]", s - 2, 4),

        ("H[N^(s+1) S^0]", s + 1, 0),
        ("H[N^(s+1) S^1]", s + 1, 1),

        ("H[N^(s+2) S^0]", s + 2, 0),
        ("H[N^(s+2) S^1]", s + 2, 1),
    ]


def natural_R_candidates(s):

    return [
        ("R[N^s]", s),
        ("R[N^(s-1)]", s - 1),
        ("R[N^(s-2)]", s - 2),
        ("R[N^(s+1)]", s + 1),
        ("R[N^(s+2)]", s + 2),
    ]


# ============================================================================
# ANCHOR AUDIT
# ============================================================================

def anchor_audit():

    print()
    print("=" * 78)
    print("5. THEOREM-ANCHOR AUDIT AGAINST H AND R")
    print("=" * 78)

    exact_hits_H = 0
    exact_hits_R = 0

    for name, k, ell, s, expected in ANCHORS:

        G = exact_G_NS(k, ell)
        H, R = split_S_plus_1(G)

        print()
        print(
            f"{name}: "
            f"k={k} ell={ell} s={s} "
            f"expected={expected}"
        )

        H_hits = []

        for label, a, b in natural_H_candidates(s):

            if a < 0 or b < 0:
                continue

            value = coefficient(H, a, b)

            if value == expected:
                H_hits.append(
                    (label, value)
                )

        R_hits = []

        for label, a in natural_R_candidates(s):

            if a < 0:
                continue

            value = remainder_coefficient(R, a)

            if value == expected:
                R_hits.append(
                    (label, value)
                )

        if H_hits:

            exact_hits_H += 1

            for label, value in H_hits:
                print(
                    f"  H MATCH: {label} = {value}"
                )
        else:

            print(
                "  H: no natural-coordinate match"
            )

        if R_hits:

            exact_hits_R += 1

            for label, value in R_hits:
                print(
                    f"  R MATCH: {label} = {value}"
                )
        else:

            print(
                "  R: no natural-coordinate match"
            )

    print()
    print(
        f"H natural exact matches = "
        f"{exact_hits_H}/9"
    )

    print(
        f"R natural exact matches = "
        f"{exact_hits_R}/9"
    )


# ============================================================================
# LOCAL H/R NEIGHBORHOOD AROUND EACH ANCHOR
# ============================================================================

def local_anchor_neighborhood():

    print()
    print("=" * 78)
    print("6. LOCAL H/R COEFFICIENT NEIGHBORHOODS")
    print("=" * 78)

    for name, k, ell, s, expected in ANCHORS:

        G = exact_G_NS(k, ell)
        H, R = split_S_plus_1(G)

        print()
        print(
            f"{name}: "
            f"k={k} ell={ell} s={s} "
            f"expected={expected}"
        )

        print("  H:")

        for a in range(
            max(0, s - 2),
            s + 3
        ):

            row = []

            for b in range(0, 5):

                value = coefficient(H, a, b)

                if value != 0:

                    row.append(
                        f"[N^{a}S^{b}]={value}"
                    )

            if row:

                print(
                    "    " + "  ".join(row)
                )

        print("  R:")

        for a in range(
            max(0, s - 2),
            s + 3
        ):

            value = remainder_coefficient(R, a)

            if value != 0:

                print(
                    f"    [N^{a}]={value}"
                )


# ============================================================================
# SMALL LINEAR-COMBINATION TEST
# ============================================================================
#
# Only very simple combinations are considered:
#
#   H[N^a S^0] + H[N^a S^1]
#   H[N^a S^0] - H[N^a S^1]
#   H[N^a S^1] + H[N^(a-1) S^0]
#   etc.
#
# This is intentionally tiny and interpretable.
#
# ============================================================================

def bounded_combination_audit():

    print()
    print("=" * 78)
    print("7. BOUNDED LOCAL H-COMBINATION AUDIT")
    print("=" * 78)

    patterns = [
        ("c00", [(0, 0, 1)]),

        ("c00+c01", [(0, 0, 1), (0, 1, 1)]),

        ("c00-c01", [(0, 0, 1), (0, 1, -1)]),

        ("c01+c10", [(0, 1, 1), (-1, 0, 1)]),

        ("c01-c10", [(0, 1, 1), (-1, 0, -1)]),

        ("c00+c10", [(0, 0, 1), (-1, 0, 1)]),

        ("c00-c10", [(0, 0, 1), (-1, 0, -1)]),

        ("c00+c01+c10",
         [(0, 0, 1), (0, 1, 1), (-1, 0, 1)]),

        ("c00-c01+c10",
         [(0, 0, 1), (0, 1, -1), (-1, 0, 1)]),

        ("c00+c01-c10",
         [(0, 0, 1), (0, 1, 1), (-1, 0, -1)]),
    ]

    global_hits = {}

    for pattern_name, pattern in patterns:

        hits = 0

        for name, k, ell, s, expected in ANCHORS:

            G = exact_G_NS(k, ell)
            H, R = split_S_plus_1(G)

            total = 0

            valid = True

            for da, db, sign in pattern:

                a = s + da
                b = db

                if a < 0:
                    valid = False
                    break

                total += (
                    sign
                    * coefficient(H, a, b)
                )

            if valid and total == expected:

                hits += 1

                print(
                    f"{pattern_name:20s} "
                    f"{name}: MATCH"
                )

        global_hits[pattern_name] = hits

    print()
    print("SUMMARY")

    for pattern_name, hits in global_hits.items():

        print(
            f"  {pattern_name:20s} "
            f"{hits}/9"
        )


# ============================================================================
# FINAL DIAGNOSTIC
# ============================================================================

def final_diagnostic():

    print()
    print("=" * 78)
    print("8. FINAL DIAGNOSTIC")
    print("=" * 78)

    reconstruction_failures = 0
    split_failures = 0

    for k in range(1, 12, 2):

        for ell in range(3, 19):

            if not reconstruction_check(k, ell):
                reconstruction_failures += 1

            if not split_check(k, ell):
                split_failures += 1

    print(
        f"F -> G reconstruction failures = "
        f"{reconstruction_failures}"
    )

    print(
        f"G -> (S+1)H + R failures = "
        f"{split_failures}"
    )

    print()

    if (
        reconstruction_failures == 0
        and split_failures == 0
    ):

        print(
            "The exact kernel now has a parity-unified decomposition:"
        )

        print()
        print(
            "    F = G(N,S) = (S+1) H(N,S) + R(N)"
        )

        print()

        print(
            "Same-parity cases have R(N)=0."
        )

        print(
            "Opposite-parity cases retain the exact obstruction R(N)."
        )

        print()

        print(
            "The next question is therefore no longer:"
        )

        print(
            "    'what is Q?'"
        )

        print(
            "but:"
        )

        print(
            "    'does the theorem coefficient live in H, in R,"
        )

        print(
            "     or in a specific boundary combination of H?'"
        )

    else:

        print(
            "The decomposition audit failed."
        )

        print(
            "Do not use coefficient comparisons yet."
        )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 204")
    print("PARITY-UNIFIED S+1 DECOMPOSITION OF G(N,S)")
    print("=" * 78)

    print()
    print(
        "No previous experiment output is read."
    )

    print(
        "No unrestricted affine candidate search is performed."
    )

    print(
        "The exact symmetric kernel is decomposed for both parity branches."
    )

    print()
    print("=" * 78)
    print("1. BASIC RECONSTRUCTION")
    print("=" * 78)

    for k, ell in [
        (1, 3),
        (1, 4),
        (1, 10),
        (3, 7),
        (3, 8),
        (3, 9),
        (3, 10),
        (3, 20),
        (5, 10),
        (5, 11),
        (5, 12),
        (5, 19),
        (7, 15),
        (7, 16),
        (7, 40),
        (9, 21),
        (11, 23),
    ]:

        print(
            f"k={k:2d} ell={ell:2d} "
            f"{'PASS' if reconstruction_check(k, ell) else 'FAIL'}"
        )

    parity_audit()

    print_remainders()

    for k, ell in [
        (3, 7),
        (3, 9),
        (5, 11),
        (5, 19),
        (7, 15),
        (9, 21),
        (11, 23),
    ]:

        print_division_summary(k, ell)

    print_H_layers()

    anchor_audit()

    local_anchor_neighborhood()

    bounded_combination_audit()

    final_diagnostic()

    print()
    print("=" * 78)
    print("END EXPERIMENT 204")
    print("=" * 78)


if __name__ == "__main__":
    main()

