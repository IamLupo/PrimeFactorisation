# ============================================================================
# EXPERIMENT 200
# EXACT N,S-COEFFICIENT IDENTIFICATION
# ============================================================================
#
# Purpose:
#   Continue directly from Experiment 199.
#
#   We start from the exact pq-kernel
#
#       F(k,l;p,q)
#         = p^k(1+q)^l + q^k(1+p)^l
#           - p^l(1+q)^k - q^l(1+p)^k
#
#   and the exact factorization
#
#       F = (p+q+1) Q(pq,p+q).
#
#   We then extract the exact polynomial
#
#       Q(N,S),   N=pq, S=p+q,
#
#   and inspect its individual coefficients
#
#       [N^a S^b] Q.
#
#   The experiment does NOT:
#       - use results.txt
#       - use previous experiment output
#       - search arbitrary affine families
#       - guess a pq -> (k,l,s) formula
#
#   The only question is:
#
#       Which natural coefficient coordinate in Q could correspond to
#       the theorem's s-index?
#
#   Tested coordinate relations are deliberately small and structural:
#
#       b = s
#       a = s
#       a+b = s
#       2a+b = s
#       b-a = s
#       a-b = s
#
#   For each relation we inspect every exact coefficient satisfying it.
#   A relation is considered interesting only if it repeatedly reproduces
#   the known anchor values.
#
# ============================================================================

import sympy as sp


# ----------------------------------------------------------------------------
# SYMBOLS
# ----------------------------------------------------------------------------

p, q = sp.symbols("p q")
N, S = sp.symbols("N S")


# ----------------------------------------------------------------------------
# KNOWN EXACT ANCHORS
# ----------------------------------------------------------------------------
#
# These are the nine nonzero exact values accumulated in Experiments 191-194.
#
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
    """
    Exact pq-level kernel from Experiment 199.
    """
    return (
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


# ----------------------------------------------------------------------------
# EXACT QUOTIENT IN p,q
# ----------------------------------------------------------------------------

def exact_Q_pq(k, ell):
    """
    Compute the exact quotient Q(pq,p+q) from

        F = (p+q+1) Q.

    No interpolation is used.
    """
    F = sp.Poly(sp.expand(exact_F(k, ell)), p, q)
    D = sp.Poly(p + q + 1, p, q)

    Q, R = sp.div(F, D)

    if not R.is_zero:
        raise ValueError(
            f"Nonzero quotient remainder for k={k}, ell={ell}: {R.as_expr()}"
        )

    return sp.expand(Q.as_expr())


# ----------------------------------------------------------------------------
# CONVERT EXACT SYMMETRIC POLYNOMIAL TO Q(N,S)
# ----------------------------------------------------------------------------

def exact_Q_NS(k, ell):
    """
    Convert the exact symmetric Q(p,q) into Q(N,S), where

        N = p q
        S = p + q.

    SymPy's exact symmetric reduction is used only as algebraic
    coordinate conversion.
    """
    Qpq = exact_Q_pq(k, ell)

    sym_expr, remainder, substitutions = sp.symmetrize(
        Qpq,
        [p, q],
        formal=True,
    )

    if sp.expand(remainder) != 0:
        raise ValueError(
            f"Symmetric remainder is nonzero for k={k}, ell={ell}: "
            f"{remainder}"
        )

    if len(substitutions) != 2:
        raise ValueError(
            f"Unexpected symmetric substitution data: {substitutions}"
        )

    s1_symbol = substitutions[0][0]  # corresponds to p+q
    s2_symbol = substitutions[1][0]  # corresponds to p*q

    expr = sp.expand(sym_expr)
    expr = expr.subs(s1_symbol, S)
    expr = expr.subs(s2_symbol, N)

    return sp.Poly(sp.expand(expr), N, S)


# ----------------------------------------------------------------------------
# EXACT RECONSTRUCTION CHECK
# ----------------------------------------------------------------------------

def reconstruction_check(k, ell):
    """
    Verify that Q(N,S), after substituting N=pq and S=p+q,
    reproduces the exact quotient Q(p,q).
    """
    Qpq = exact_Q_pq(k, ell)
    QNS = exact_Q_NS(k, ell)

    reconstructed = sp.expand(
        QNS.as_expr().subs({
            N: p * q,
            S: p + q,
        })
    )

    return sp.expand(Qpq - reconstructed) == 0


# ----------------------------------------------------------------------------
# COEFFICIENT DICTIONARY
# ----------------------------------------------------------------------------

def coefficient_dict(QNS):
    """
    Return {(a,b): coefficient} for

        Q(N,S) = sum c[a,b] N^a S^b.
    """
    result = {}

    for (a, b), coeff in QNS.terms():
        result[(int(a), int(b))] = int(coeff)

    return result


# ----------------------------------------------------------------------------
# NATURAL COORDINATE RELATIONS
# ----------------------------------------------------------------------------

def relation_matches(relation_name, a, b, s):
    """
    Deliberately small structural relations.

    Returns True exactly when (a,b,s) satisfies the chosen relation.
    """
    if relation_name == "b=s":
        return b == s

    if relation_name == "a=s":
        return a == s

    if relation_name == "a+b=s":
        return a + b == s

    if relation_name == "2a+b=s":
        return 2 * a + b == s

    if relation_name == "b-a=s":
        return b - a == s

    if relation_name == "a-b=s":
        return a - b == s

    raise ValueError(f"Unknown relation: {relation_name}")


RELATIONS = [
    "b=s",
    "a=s",
    "a+b=s",
    "2a+b=s",
    "b-a=s",
    "a-b=s",
]


# ----------------------------------------------------------------------------
# ANCHOR INSPECTION
# ----------------------------------------------------------------------------

def inspect_anchor(name, k, ell, s, expected, relation_name):
    """
    For one anchor and one coordinate relation, list all Q coefficients
    satisfying that relation.

    We deliberately do NOT select one coefficient by another heuristic.
    Every qualifying coefficient is reported.
    """
    QNS = exact_Q_NS(k, ell)
    coeffs = coefficient_dict(QNS)

    candidates = []

    for (a, b), coeff in sorted(coeffs.items()):
        if relation_matches(relation_name, a, b, s):
            candidates.append((a, b, coeff))

    exact_hits = [
        (a, b, coeff)
        for (a, b, coeff) in candidates
        if coeff == expected
    ]

    abs_hits = [
        (a, b, coeff)
        for (a, b, coeff) in candidates
        if abs(coeff) == abs(expected)
    ]

    return candidates, exact_hits, abs_hits


# ----------------------------------------------------------------------------
# DIRECT ANCHOR SUMMARY
# ----------------------------------------------------------------------------

def anchor_summary():
    print("=" * 78)
    print("EXPERIMENT 200")
    print("EXACT N,S-COEFFICIENT IDENTIFICATION")
    print("=" * 78)
    print()
    print("No previous experiment output is read.")
    print("All pq, Q(N,S), and coefficient data are generated locally.")
    print()
    print("Target:")
    print("    F(p,q) -> Q(N,S),  N=pq, S=p+q")
    print("    identify natural coefficient coordinates for theorem index s")
    print()

    print("=" * 78)
    print("1. EXACT Q(N,S) RECONSTRUCTION SANITY CHECK")
    print("=" * 78)

    test_cases = [
        (3, 7),
        (3, 9),
        (5, 11),
        (5, 19),
        (7, 15),
        (9, 21),
        (11, 23),
    ]

    failures = 0

    for k, ell in test_cases:
        ok = reconstruction_check(k, ell)
        print(
            f"k={k:2d} ell={ell:2d} "
            f"reconstruction={'PASS' if ok else 'FAIL'}"
        )
        if not ok:
            failures += 1

    print()
    print(f"reconstruction failures = {failures}")

    print()
    print("=" * 78)
    print("2. EXACT Q(N,S) COEFFICIENT TABLES")
    print("=" * 78)

    table_cases = [
        (3, 7),
        (3, 9),
        (5, 11),
        (5, 19),
    ]

    for k, ell in table_cases:
        print()
        print(f"k={k} ell={ell}")

        QNS = exact_Q_NS(k, ell)
        coeffs = coefficient_dict(QNS)

        for (a, b), coeff in sorted(coeffs.items(), reverse=True):
            print(
                f"  [N^{a} S^{b}] = {coeff}"
            )

    print()
    print("=" * 78)
    print("3. NATURAL COORDINATE SEARCH AGAINST NINE ANCHORS")
    print("=" * 78)

    relation_scores = {}

    for relation_name in RELATIONS:
        exact_match_count = 0
        abs_match_count = 0
        anchors_with_candidates = 0
        total_candidate_coefficients = 0

        print()
        print("-" * 78)
        print(relation_name)
        print("-" * 78)

        for name, k, ell, s, expected in ANCHORS:
            candidates, exact_hits, abs_hits = inspect_anchor(
                name,
                k,
                ell,
                s,
                expected,
                relation_name,
            )

            if candidates:
                anchors_with_candidates += 1

            total_candidate_coefficients += len(candidates)

            if exact_hits:
                exact_match_count += 1

            if abs_hits:
                abs_match_count += 1

            print(
                f"{name}: k={k:2d} ell={ell:2d} s={s:2d} "
                f"expected={expected:10d}"
            )

            if not candidates:
                print("    no coefficients satisfy the coordinate relation")
                continue

            print("    coefficients:")
            for a, b, coeff in candidates:
                marker = ""
                if coeff == expected:
                    marker = "  <-- EXACT"
                elif abs(coeff) == abs(expected):
                    marker = "  <-- ABS"

                print(
                    f"      (a={a:2d}, b={b:2d}) -> {coeff:12d}{marker}"
                )

        relation_scores[relation_name] = (
            exact_match_count,
            abs_match_count,
            anchors_with_candidates,
            total_candidate_coefficients,
        )

        print()
        print(
            f"SUMMARY {relation_name}: "
            f"exact={exact_match_count}/9, "
            f"absolute={abs_match_count}/9, "
            f"anchors_with_candidates={anchors_with_candidates}/9, "
            f"candidate_coefficients={total_candidate_coefficients}"
        )

    print()
    print("=" * 78)
    print("4. RELATION RANKING")
    print("=" * 78)

    ranking = sorted(
        relation_scores.items(),
        key=lambda item: (
            item[1][0],   # exact matches
            item[1][1],   # absolute matches
            item[1][2],   # anchors covered
            -item[1][3],  # fewer candidate coefficients is preferable
        ),
        reverse=True,
    )

    for rank, (relation_name, score) in enumerate(ranking, 1):
        exact_count, abs_count, covered, candidates = score

        print(
            f"{rank:2d}. {relation_name:10s} "
            f"exact={exact_count}/9 "
            f"absolute={abs_count}/9 "
            f"covered={covered}/9 "
            f"candidate_coefficients={candidates}"
        )

    print()
    print("=" * 78)
    print("5. DIRECT ANCHOR COEFFICIENT LOCATIONS")
    print("=" * 78)

    print()
    print("For every exact anchor, report every coefficient in Q(N,S)")
    print("that equals the anchor value, regardless of the s-relation.")
    print()

    for name, k, ell, s, expected in ANCHORS:
        QNS = exact_Q_NS(k, ell)
        coeffs = coefficient_dict(QNS)

        exact_locations = [
            (a, b, coeff)
            for (a, b), coeff in sorted(coeffs.items())
            if coeff == expected
        ]

        absolute_locations = [
            (a, b, coeff)
            for (a, b), coeff in sorted(coeffs.items())
            if abs(coeff) == abs(expected)
        ]

        print(
            f"{name}: k={k:2d} ell={ell:2d} s={s:2d} "
            f"expected={expected}"
        )

        if exact_locations:
            print("  EXACT coefficient locations:")
            for a, b, coeff in exact_locations:
                print(
                    f"    (a={a}, b={b}) -> {coeff}"
                )
        else:
            print("  EXACT coefficient locations: none")

        if absolute_locations:
            print("  ABSOLUTE coefficient locations:")
            for a, b, coeff in absolute_locations:
                print(
                    f"    (a={a}, b={b}) -> {coeff}"
                )
        else:
            print("  ABSOLUTE coefficient locations: none")

    print()
    print("=" * 78)
    print("6. TOP-LAYER COEFFICIENT DIAGNOSTIC")
    print("=" * 78)

    for k, ell in test_cases:
        QNS = exact_Q_NS(k, ell)
        coeffs = coefficient_dict(QNS)

        max_b = max(b for (a, b) in coeffs)
        top_terms = [
            (a, b, coeff)
            for (a, b), coeff in sorted(
                coeffs.items(),
                key=lambda item: (item[0][1], item[0][0]),
                reverse=True,
            )
            if b >= max_b - 2
        ]

        print()
        print(f"k={k} ell={ell}  top S-degree={max_b}")

        for a, b, coeff in top_terms:
            print(
                f"  [N^{a} S^{b}] = {coeff}"
            )

    print()
    print("=" * 78)
    print("7. FINAL DIAGNOSTIC")
    print("=" * 78)

    best_relation, best_score = ranking[0]
    best_exact, best_abs, best_covered, best_candidates = best_score

    print()
    print(
        f"best natural coordinate relation = {best_relation}"
    )
    print(
        f"exact anchor agreement              = "
        f"{best_exact}/9"
    )
    print(
        f"absolute-value agreement             = "
        f"{best_abs}/9"
    )
    print(
        f"anchors having at least one candidate = "
        f"{best_covered}/9"
    )

    print()
    if best_exact >= 5:
        print("A potentially meaningful coefficient coordinate has emerged.")
        print()
        print("NEXT QUESTION:")
        print("Determine whether the same (a,b) coordinate is selected")
        print("uniformly across k, ell, and the full support.")
    elif best_abs >= 5:
        print("Only a possible global sign issue is suggested.")
        print()
        print("NEXT QUESTION:")
        print("Test the sign using the exact pq normalization.")
    else:
        print("No natural coefficient coordinate explains the anchors.")
        print()
        print("The theorem's C-value is therefore probably not a single")
        print("raw coefficient [N^a S^b]Q.")
        print()
        print("NEXT QUESTION:")
        print("Test short linear functionals of one exact Q-layer,")
        print("starting from the structural boundary layers in Experiment 199.")
        print()
        print("Do NOT return to unrestricted affine candidate searches.")

    print()
    print("=" * 78)
    print("END EXPERIMENT 200")
    print("=" * 78)


# ----------------------------------------------------------------------------
# MAIN
# ----------------------------------------------------------------------------

if __name__ == "__main__":
    anchor_summary()

