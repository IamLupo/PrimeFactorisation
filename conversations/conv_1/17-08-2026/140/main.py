import sympy as sp


# =============================================================================
# EXPERIMENT 206
# THEOREM-INDEX DIAGONAL AUDIT IN G(N,S+1)
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
# EXACT SYMMETRIC REDUCTION
# =============================================================================

def exact_G(k, ell):
    F = exact_F(k, ell)

    sym_part, remainder, mapping = sp.symmetrize(
        F,
        [p, q],
        formal=True,
    )

    remainder = sp.expand(remainder)

    if remainder != 0:
        raise ValueError(
            f"Symmetric reduction failed for "
            f"k={k}, ell={ell}: remainder={remainder}"
        )

    formal_sum = None
    formal_prod = None

    for sym, expr in mapping:
        expr = sp.expand(expr)

        if sp.expand(expr - (p + q)) == 0:
            formal_sum = sym

        elif sp.expand(expr - p*q) == 0:
            formal_prod = sym

    if formal_sum is None or formal_prod is None:
        raise ValueError(
            f"Could not identify symmetric variables: "
            f"mapping={mapping}"
        )

    return sp.expand(
        sym_part.subs(
            {
                formal_sum: S,
                formal_prod: N,
            }
        )
    )


# =============================================================================
# SHIFTED X BASIS
# =============================================================================

def shifted_coefficients(k, ell):
    """
    G(N,S) = sum_j C_j(N) (S+1)^j.
    """

    G = exact_G(k, ell)

    shifted = sp.expand(
        G.subs(S, X - 1)
    )

    if shifted == 0:
        return {0: sp.Integer(0)}

    poly = sp.Poly(
        shifted,
        X,
        domain=sp.EX,
    )

    degree = poly.degree()

    if degree is None:
        return {0: sp.Integer(0)}

    degree = int(degree)

    return {
        j: sp.expand(
            poly.coeff_monomial(X**j)
        )
        for j in range(degree + 1)
    }


# =============================================================================
# N COEFFICIENT
# =============================================================================

def N_coeff(poly_expr, d):
    """
    Return [N^d] poly_expr safely.

    Important:
    Do not use sp.NegativeInfinity. Some SymPy versions do not expose
    that name at the module level.
    """

    if d < 0:
        return sp.Integer(0)

    expr = sp.expand(poly_expr)

    if expr == 0:
        return sp.Integer(0)

    poly = sp.Poly(
        expr,
        N,
        domain=sp.EX,
    )

    return sp.expand(
        poly.coeff_monomial(N**d)
    )


# =============================================================================
# SAFE COEFFICIENT LOOKUP
# =============================================================================

def get_coeff(coeffs, j, d):
    if j < 0 or d < 0:
        return sp.Integer(0)

    if j not in coeffs:
        return sp.Integer(0)

    return N_coeff(
        coeffs[j],
        d,
    )


# =============================================================================
# MOTIVATED COORDINATES
# =============================================================================

def motivated_coordinates(k, ell, s):
    """
    Bounded local neighborhood around

        j = ell-s+1
        d = s-1
    """

    coords = [
        ("D0",          ell - s + 1, s - 1),

        ("D_left",      ell - s,     s - 1),
        ("D_right",     ell - s + 2, s - 1),

        ("D_down",      ell - s + 1, s - 2),
        ("D_up",        ell - s + 1, s),

        ("D_diag_left", ell - s,     s),
        ("D_diag_right",ell - s + 2, s - 2),

        ("D_SE",        ell - s + 2, s),
        ("D_NW",        ell - s,     s - 2),
    ]

    return coords


# =============================================================================
# FINITE DIFFERENCE CANDIDATES
# =============================================================================

def finite_difference_candidates(k, ell, s, coeffs):

    j = ell - s + 1
    d = s - 1

    a = get_coeff(coeffs, j, d)
    b = get_coeff(coeffs, j + 1, d)
    c = get_coeff(coeffs, j - 1, d)

    a2 = get_coeff(coeffs, j, d - 1)
    b2 = get_coeff(coeffs, j + 1, d - 1)

    return [
        ("FD_X_plus",    b - a),
        ("FD_X_minus",   a - c),
        ("FD_X_sum",     a + b),
        ("FD_X_sym",     b + c),
        ("FD_N_plus",    a - a2),
        ("FD_N_reverse", a2 - a),
        ("FD_XN",        (b - a) - (b2 - a2)),
    ]


# =============================================================================
# RECONSTRUCTION
# =============================================================================

def reconstruction_audit():

    print()
    print("=" * 78)
    print("1. SHIFTED-BASIS RECONSTRUCTION")
    print("=" * 78)

    cases = [
        (1, 3),
        (1, 4),
        (1, 10),
        (3, 7),
        (3, 8),
        (3, 9),
        (3, 10),
        (5, 10),
        (5, 11),
        (5, 12),
        (5, 19),
        (7, 14),
        (7, 15),
        (7, 16),
        (9, 20),
        (9, 21),
        (11, 23),
    ]

    failures = 0

    for k, ell in cases:

        G = exact_G(k, ell)
        coeffs = shifted_coefficients(k, ell)

        rebuilt = sp.Integer(0)

        for j, Cj in coeffs.items():
            rebuilt += Cj * (S + 1)**j

        ok = (
            sp.expand(G - rebuilt) == 0
        )

        if not ok:
            failures += 1

        print(
            f"k={k:2d} ell={ell:2d} "
            f"{'PASS' if ok else 'FAIL'}"
        )

    print()
    print(
        f"reconstruction failures = {failures}"
    )


# =============================================================================
# CENTRAL DIAGONAL AUDIT
# =============================================================================

def diagonal_audit():

    print()
    print("=" * 78)
    print("2. CENTRAL THEOREM-INDEX DIAGONAL AUDIT")
    print("=" * 78)

    labels = [
        "D0",
        "D_left",
        "D_right",
        "D_down",
        "D_up",
        "D_diag_left",
        "D_diag_right",
        "D_SE",
        "D_NW",
    ]

    exact_counts = {
        x: 0
        for x in labels
    }

    abs_counts = {
        x: 0
        for x in labels
    }

    for name, k, ell, s, expected in ANCHORS:

        coeffs = shifted_coefficients(
            k,
            ell,
        )

        print()
        print(
            f"{name}: k={k} ell={ell} "
            f"s={s} expected={expected}"
        )

        found = False

        for label, j, d in motivated_coordinates(
            k,
            ell,
            s,
        ):

            value = get_coeff(
                coeffs,
                j,
                d,
            )

            print(
                f"  {label:13s} "
                f"[N^{d} X^{j}] = {value}"
            )

            if value == expected:
                exact_counts[label] += 1
                found = True

            if abs(value) == abs(expected):
                abs_counts[label] += 1

        if not found:
            print(
                "  no exact central-diagonal match"
            )

    print()
    print("EXACT MATCH COUNTS")

    for label in labels:
        print(
            f"  {label:13s}: "
            f"{exact_counts[label]}/9"
        )

    print()
    print("ABSOLUTE-VALUE MATCH COUNTS")

    for label in labels:
        print(
            f"  {label:13s}: "
            f"{abs_counts[label]}/9"
        )


# =============================================================================
# FINITE DIFFERENCE AUDIT
# =============================================================================

def finite_difference_audit():

    print()
    print("=" * 78)
    print("3. LOCAL FINITE-DIFFERENCE AUDIT")
    print("=" * 78)

    labels = [
        "FD_X_plus",
        "FD_X_minus",
        "FD_X_sum",
        "FD_X_sym",
        "FD_N_plus",
        "FD_N_reverse",
        "FD_XN",
    ]

    exact_counts = {
        x: 0
        for x in labels
    }

    abs_counts = {
        x: 0
        for x in labels
    }

    for name, k, ell, s, expected in ANCHORS:

        coeffs = shifted_coefficients(
            k,
            ell,
        )

        print()
        print(
            f"{name}: expected={expected}"
        )

        for label, value in finite_difference_candidates(
            k,
            ell,
            s,
            coeffs,
        ):

            print(
                f"  {label:15s} = {value}"
            )

            if value == expected:
                exact_counts[label] += 1

            if abs(value) == abs(expected):
                abs_counts[label] += 1

    print()
    print("FINITE-DIFFERENCE EXACT COUNTS")

    for label in labels:
        print(
            f"  {label:15s}: "
            f"{exact_counts[label]}/9"
        )

    print()
    print("FINITE-DIFFERENCE ABS COUNTS")

    for label in labels:
        print(
            f"  {label:15s}: "
            f"{abs_counts[label]}/9"
        )


# =============================================================================
# RATIO / CONSISTENCY AUDIT
# =============================================================================

def consistency_audit():

    print()
    print("=" * 78)
    print("4. CROSS-ANCHOR CONSISTENCY")
    print("=" * 78)

    for name, k, ell, s, expected in ANCHORS:

        coeffs = shifted_coefficients(
            k,
            ell,
        )

        print()
        print(
            f"{name}: expected={expected}"
        )

        for label, j, d in motivated_coordinates(
            k,
            ell,
            s,
        ):

            value = get_coeff(
                coeffs,
                j,
                d,
            )

            if expected != 0:
                ratio = sp.factor(
                    sp.Rational(
                        value,
                        expected,
                    )
                )
            else:
                ratio = "undefined"

            print(
                f"  {label:13s} "
                f"value={value} "
                f"ratio={ratio}"
            )


# =============================================================================
# MATCHING ANCHOR DIAGNOSTIC
# =============================================================================

def matching_anchor_diagnostics():

    print()
    print("=" * 78)
    print("5. KNOWN MATCHING-ANCHOR DIAGNOSTICS")
    print("=" * 78)

    for name, k, ell, s, expected in ANCHORS:

        coeffs = shifted_coefficients(
            k,
            ell,
        )

        matches = []

        for label, j, d in motivated_coordinates(
            k,
            ell,
            s,
        ):

            value = get_coeff(
                coeffs,
                j,
                d,
            )

            if value == expected:
                matches.append(
                    (
                        label,
                        j,
                        d,
                        value,
                    )
                )

        if matches:

            print()
            print(
                f"{name}: expected={expected}"
            )

            for match in matches:
                print(
                    f"  {match}"
                )


# =============================================================================
# PARITY COMPARISON
# =============================================================================

def parity_comparison():

    print()
    print("=" * 78)
    print("6. PARITY COMPARISON ON THE CENTRAL DIAGONAL")
    print("=" * 78)

    cases = [
        (1, 9, 4),
        (1, 10, 4),
        (3, 9, 4),
        (3, 10, 4),
        (5, 11, 5),
        (5, 12, 5),
        (7, 15, 7),
        (7, 16, 7),
        (9, 21, 10),
        (9, 22, 10),
    ]

    for k, ell, s in cases:

        coeffs = shifted_coefficients(
            k,
            ell,
        )

        j = ell - s + 1
        d = s - 1

        value = get_coeff(
            coeffs,
            j,
            d,
        )

        parity = (
            "same"
            if (k - ell) % 2 == 0
            else "different"
        )

        print(
            f"k={k:2d} ell={ell:2d} "
            f"s={s:2d} "
            f"parity={parity:9s} "
            f"[N^{d}X^{j}]={value}"
        )


# =============================================================================
# FINAL DIAGNOSTIC
# =============================================================================

def final_diagnostic():

    print()
    print("=" * 78)
    print("7. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "The shifted basis is"
    )

    print(
        "    G(N,S) = sum_j C_j(N) (S+1)^j."
    )

    print()
    print(
        "The central coordinate being tested is"
    )

    print(
        "    [N^(s-1) (S+1)^(ell-s+1)] G."
    )

    print()
    print(
        "The bounded neighboring coordinates and first finite"
    )

    print(
        "differences are being tested because Experiment 205"
    )

    print(
        "produced only isolated exact matches."
    )

    print()
    print(
        "No finite-anchor match is accepted as the theorem formula"
    )

    print(
        "unless the same coordinate relation becomes structurally"
    )

    print(
        "stable across many anchors and both parity branches."
    )

    print()
    print(
        "The next derivation should use the strongest surviving"
    )

    print(
        "coordinate relation and trace it back to the exact pq kernel."
    )


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 206")
    print("THEOREM-INDEX DIAGONAL AUDIT IN G(N,S+1)")
    print("=" * 78)

    print()
    print(
        "No previous experiment output is read."
    )

    print(
        "No unrestricted affine candidate search is performed."
    )

    reconstruction_audit()
    diagonal_audit()
    finite_difference_audit()
    consistency_audit()
    matching_anchor_diagnostics()
    parity_comparison()
    final_diagnostic()

    print()
    print("=" * 78)
    print("END EXPERIMENT 206")
    print("=" * 78)


if __name__ == "__main__":
    main()