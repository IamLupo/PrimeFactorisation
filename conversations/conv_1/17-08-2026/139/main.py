import sympy as sp


# ============================================================================
# EXPERIMENT 205
# BOUNDARY-ADAPTED SHIFTED-S COEFFICIENT AUDIT
# ============================================================================
#
# Exact kernel:
#
#   F(p,q)
#     = p^k (1+q)^ell + q^k (1+p)^ell
#       - p^ell (1+q)^k - q^ell (1+p)^k
#
# Exact symmetric reduction:
#
#   N = p q
#   S = p + q
#
#   F(p,q) = G(N,S)
#
# New basis:
#
#   X = S + 1
#
#   G(N,S) = sum_j C_j(N) X^j
#
# The purpose is to isolate the parity obstruction in
#
#   C_0(N) = G(N,-1)
#
# while studying the higher post-boundary layers
#
#   C_1(N), C_2(N), ...
#
# No previous experiment output is read.
# No unrestricted affine candidate search is performed.
#
# ============================================================================


p, q = sp.symbols("p q")
N, S, X = sp.symbols("N S X")


# ============================================================================
# THEOREM ANCHORS
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
# EXACT pq KERNEL
# ============================================================================

def exact_F(k, ell):
    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


# ============================================================================
# EXACT SYMMETRIC REDUCTION
# ============================================================================

def exact_G(k, ell):
    """
    Convert the exact symmetric polynomial F(p,q) into
    G(N,S), where

        N = p*q
        S = p+q.

    Symmetrize is used directly; no division by p+q+1 occurs.
    """

    F = exact_F(k, ell)

    symmetric_part, remainder, mapping = sp.symmetrize(
        F,
        [p, q],
        formal=True,
    )

    remainder = sp.expand(remainder)

    if remainder != 0:
        raise ValueError(
            "Symmetric reduction failed for "
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
            "Could not identify elementary symmetric variables. "
            f"mapping={mapping}"
        )

    G = sp.expand(
        symmetric_part.subs(
            {
                formal_sum: S,
                formal_prod: N,
            }
        )
    )

    return G


# ============================================================================
# SHIFTED-S COEFFICIENTS
# ============================================================================

def shifted_S_coefficients(k, ell):
    """
    Return the exact coefficients C_j(N) in

        G(N,S) = sum_j C_j(N) (S+1)^j.

    We use X=S+1, hence S=X-1.
    """

    G = exact_G(k, ell)

    shifted = sp.expand(
        G.subs(S, X - 1)
    )

    # Explicit zero-polynomial handling.
    if shifted == 0:
        return {0: sp.Integer(0)}

    poly = sp.Poly(
        shifted,
        X,
        domain=sp.EX,
    )

    degree = poly.degree(X)

    # SymPy returns -oo for the zero polynomial.
    if degree == sp.S.NegativeInfinity:
        return {0: sp.Integer(0)}

    degree_int = int(degree)

    coeffs = {}

    for j in range(degree_int + 1):
        coeffs[j] = sp.expand(
            poly.coeff_monomial(X**j)
        )

    return coeffs


# ============================================================================
# EXACT SHIFTED-BASIS RECONSTRUCTION
# ============================================================================

def shifted_reconstruction_check(k, ell):
    """
    Verify

        G(N,S) = sum_j C_j(N) (S+1)^j.
    """

    G = exact_G(k, ell)

    coeffs = shifted_S_coefficients(k, ell)

    rebuilt = sp.Integer(0)

    for j, Cj in coeffs.items():
        rebuilt += Cj * (S + 1)**j

    rebuilt = sp.expand(rebuilt)

    return sp.expand(G - rebuilt) == 0


# ============================================================================
# RECONSTRUCTION AUDIT
# ============================================================================

def reconstruction_audit():

    print()
    print("=" * 78)
    print("1. SHIFTED-S RECONSTRUCTION AUDIT")
    print("=" * 78)

    failures = 0

    for k in range(1, 12, 2):
        for ell in range(3, 19):

            ok = shifted_reconstruction_check(
                k,
                ell,
            )

            if not ok:
                failures += 1

            print(
                f"k={k:2d} ell={ell:2d} "
                f"{'PASS' if ok else 'FAIL'}"
            )

    print()
    print(
        f"shifted-basis reconstruction failures = {failures}"
    )


# ============================================================================
# PARITY AUDIT
# ============================================================================

def parity_audit():

    print()
    print("=" * 78)
    print("2. SHIFTED CONSTANT TERM / PARITY AUDIT")
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

        G = exact_G(k, ell)
        coeffs = shifted_S_coefficients(k, ell)

        C0 = sp.expand(
            coeffs.get(0, 0)
        )

        direct = sp.expand(
            G.subs(S, -1)
        )

        direct_ok = (
            sp.expand(C0 - direct) == 0
        )

        same_parity = (
            (k - ell) % 2 == 0
        )

        zero_expected = same_parity
        zero_observed = (
            C0 == 0
        )

        ok = (
            direct_ok
            and zero_observed == zero_expected
        )

        if not ok:
            failures += 1

        print(
            f"k={k:2d} ell={ell:2d} "
            f"same_parity={same_parity} "
            f"C0_zero={zero_observed} "
            f"{'PASS' if ok else 'FAIL'}"
        )

        if not ok:
            print(
                f"  C0={C0}"
            )
            print(
                f"  G(N,-1)={direct}"
            )

    print()
    print(
        f"shifted parity failures = {failures}"
    )


# ============================================================================
# MOTIVATED X-LAYER INDICES
# ============================================================================

def motivated_j_indices(k, ell):

    candidates = set()

    def add(value):
        if value >= 0:
            candidates.add(int(value))

    # Constant / first few layers.
    add(0)
    add(1)
    add(2)
    add(3)

    # Boundary locations suggested by the previous experiments.
    add(ell - k)
    add(ell - k - 1)
    add(ell - k + 1)

    # Half-difference location.
    if (ell - k) % 2 == 0:
        h = (ell - k) // 2
        add(h - 1)
        add(h)
        add(h + 1)

    # Half-sum location.
    if (ell + k) % 2 == 0:
        h = (ell + k) // 2
        add(h - 1)
        add(h)
        add(h + 1)

    return sorted(candidates)


# ============================================================================
# MOTIVATED N-DEGREES
# ============================================================================

def motivated_N_indices(s, k, ell, Cj):

    if Cj == 0:
        return []

    poly = sp.Poly(
        Cj,
        N,
        domain=sp.EX,
    )

    degree = poly.degree(N)

    if degree == sp.S.NegativeInfinity:
        return []

    max_degree = int(degree)

    candidates = set()

    # Direct neighborhood around s.
    for d in range(s - 3, s + 4):
        if 0 <= d <= max_degree:
            candidates.add(d)

    # Natural boundary-related transforms.
    values = [
        ell - s,
        ell - s - 1,
        ell - s + 1,
        s - k,
        s - k - 1,
        s - k + 1,
        s + k,
        s + k - 1,
        s + k + 1,
    ]

    for d in values:
        if 0 <= d <= max_degree:
            candidates.add(d)

    return sorted(candidates)


# ============================================================================
# COEFFICIENT OF N^d
# ============================================================================

def coeff_N(poly_expr, degree):

    if poly_expr == 0:
        return sp.Integer(0)

    if degree < 0:
        return sp.Integer(0)

    poly = sp.Poly(
        poly_expr,
        N,
        domain=sp.EX,
    )

    return sp.expand(
        poly.coeff_monomial(
            N**degree
        )
    )


# ============================================================================
# PRINT SELECTED SHIFTED LAYERS
# ============================================================================

def print_shifted_rows():

    print()
    print("=" * 78)
    print("3. SELECTED SHIFTED-S COEFFICIENT ROWS")
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

        coeffs = shifted_S_coefficients(
            k,
            ell,
        )

        for j in motivated_j_indices(k, ell):

            if j not in coeffs:
                continue

            Cj = sp.expand(
                coeffs[j]
            )

            if Cj == 0:
                print(
                    f"  X^{j}: ZERO"
                )
                continue

            poly = sp.Poly(
                Cj,
                N,
                domain=sp.EX,
            )

            degree = poly.degree(N)

            if degree == sp.S.NegativeInfinity:
                print(
                    f"  X^{j}: ZERO"
                )
                continue

            degree = int(degree)

            # Show a bounded low-N window.
            terms = []

            for d in range(
                min(8, degree + 1)
            ):
                c = coeff_N(
                    Cj,
                    d,
                )
                if c != 0:
                    terms.append(
                        (d, c)
                    )

            print(
                f"  X^{j}: degree_N={degree}"
            )

            if terms:
                print(
                    f"      low-N coefficients={terms}"
                )


# ============================================================================
# THEOREM-ANCHOR SEARCH
# ============================================================================

def anchor_shifted_search():

    print()
    print("=" * 78)
    print("4. THEOREM-ANCHOR SEARCH IN SHIFTED-S BASIS")
    print("=" * 78)

    exact_hits = 0

    for name, k, ell, s, expected in ANCHORS:

        coeffs = shifted_S_coefficients(
            k,
            ell,
        )

        print()
        print(
            f"{name}: "
            f"k={k} ell={ell} s={s} "
            f"expected={expected}"
        )

        hits = []

        for j in motivated_j_indices(
            k,
            ell,
        ):

            if j not in coeffs:
                continue

            Cj = coeffs[j]

            for d in motivated_N_indices(
                s,
                k,
                ell,
                Cj,
            ):

                value = coeff_N(
                    Cj,
                    d,
                )

                if value == expected:
                    hits.append(
                        (j, d, value)
                    )

        if hits:

            exact_hits += 1

            for j, d, value in hits:
                print(
                    f"  MATCH: "
                    f"[X^{j} N^{d}]={value}"
                )

        else:

            print(
                "  no motivated shifted-basis match"
            )

    print()
    print(
        f"shifted-basis exact anchor matches = "
        f"{exact_hits}/9"
    )


# ============================================================================
# ABSOLUTE-VALUE DIAGNOSTIC
# ============================================================================

def anchor_absolute_search():

    print()
    print("=" * 78)
    print("5. ABSOLUTE-VALUE / SIGN DIAGNOSTIC")
    print("=" * 78)

    for name, k, ell, s, expected in ANCHORS:

        coeffs = shifted_S_coefficients(
            k,
            ell,
        )

        matches = []

        for j in motivated_j_indices(
            k,
            ell,
        ):

            if j not in coeffs:
                continue

            Cj = coeffs[j]

            for d in motivated_N_indices(
                s,
                k,
                ell,
                Cj,
            ):

                value = coeff_N(
                    Cj,
                    d,
                )

                if abs(value) == abs(expected):

                    matches.append(
                        (
                            j,
                            d,
                            value,
                            (
                                "same-sign"
                                if value == expected
                                else "opposite-sign"
                            ),
                        )
                    )

        print()
        print(
            f"{name}: expected={expected}"
        )

        if matches:

            for item in matches:
                print(
                    f"  {item}"
                )

        else:

            print(
                "  no absolute-value matches"
            )


# ============================================================================
# POST-BOUNDARY X-LAYER AUDIT
# ============================================================================

def boundary_layer_audit():

    print()
    print("=" * 78)
    print("6. POST-BOUNDARY X-LAYER AUDIT")
    print("=" * 78)

    for name, k, ell, s, expected in ANCHORS:

        coeffs = shifted_S_coefficients(
            k,
            ell,
        )

        base_j = ell - k

        print()
        print(
            f"{name}: "
            f"k={k} ell={ell} s={s} "
            f"expected={expected}"
        )

        for j in [
            base_j - 2,
            base_j - 1,
            base_j,
            base_j + 1,
            base_j + 2,
        ]:

            if j < 0 or j not in coeffs:
                continue

            Cj = coeffs[j]

            if Cj == 0:
                print(
                    f"  X^{j}: ZERO"
                )
                continue

            poly = sp.Poly(
                Cj,
                N,
                domain=sp.EX,
            )

            degree = poly.degree(N)

            if degree == sp.S.NegativeInfinity:
                print(
                    f"  X^{j}: ZERO"
                )
                continue

            degree = int(degree)

            interesting = []

            for d in range(
                max(0, s - 3),
                min(degree, s + 3) + 1,
            ):

                value = coeff_N(
                    Cj,
                    d,
                )

                if value != 0:
                    interesting.append(
                        (d, value)
                    )

            print(
                f"  X^{j}: {interesting}"
            )


# ============================================================================
# SUCCESSIVE X-LAYER RECURRENCE AUDIT
# ============================================================================

def layer_recurrence_audit():

    print()
    print("=" * 78)
    print("7. SUCCESSIVE X-LAYER RECURRENCE AUDIT")
    print("=" * 78)

    for name, k, ell, s, expected in ANCHORS:

        coeffs = shifted_S_coefficients(
            k,
            ell,
        )

        print()
        print(
            f"{name}: expected={expected}"
        )

        found = False

        for j in motivated_j_indices(
            k,
            ell,
        ):

            if j not in coeffs:
                continue

            if j + 1 not in coeffs:
                continue

            A = coeffs[j]
            B = coeffs[j + 1]

            if A == 0 and B == 0:
                continue

            for d in range(
                max(0, s - 3),
                s + 4,
            ):

                a = coeff_N(A, d)
                b = coeff_N(B, d)

                tests = [
                    ("A", a),
                    ("B", b),
                    ("A+B", a + b),
                    ("A-B", a - b),
                    ("B-A", b - a),
                ]

                for label, value in tests:

                    if value == expected:

                        print(
                            f"  MATCH: "
                            f"j={j}, d={d}, "
                            f"{label}={value}"
                        )

                        found = True

        if not found:

            print(
                "  no local recurrence match"
            )


# ============================================================================
# DIRECT C_0 / HIGHER-LAYER SPLIT
# ============================================================================

def c0_higher_layer_audit():

    print()
    print("=" * 78)
    print("8. C0 VS HIGHER-LAYER STRUCTURE")
    print("=" * 78)

    cases = [
        (1, 4),
        (1, 10),
        (3, 8),
        (3, 10),
        (5, 12),
        (7, 14),
        (9, 20),
        (11, 24),
    ]

    for k, ell in cases:

        coeffs = shifted_S_coefficients(
            k,
            ell,
        )

        C0 = sp.expand(
            coeffs.get(0, 0)
        )

        higher = [
            j
            for j, Cj in coeffs.items()
            if j > 0 and Cj != 0
        ]

        print()
        print(
            f"k={k} ell={ell} "
            f"same_parity={((k-ell) % 2 == 0)}"
        )

        print(
            f"  C0(N) = {sp.factor(C0)}"
        )

        print(
            f"  nonzero higher X-layers = "
            f"{higher[:12]}"
            + (" ..." if len(higher) > 12 else "")
        )


# ============================================================================
# FINAL DIAGNOSTIC
# ============================================================================

def final_diagnostic():

    print()
    print("=" * 78)
    print("9. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "The exact representation being tested is"
    )

    print()
    print(
        "    G(N,S) = sum_j C_j(N) (S+1)^j"
    )

    print()
    print(
        "where"
    )

    print(
        "    C_0(N) = G(N,-1)."
    )

    print()
    print(
        "This isolates the parity obstruction entirely in C_0."
    )

    print()
    print(
        "The higher C_j(N) are uniform post-boundary layers"
    )

    print(
        "available in both parity branches."
    )

    print()
    print(
        "The decisive diagnostics are:"
    )

    print(
        "  1. shifted-basis anchor matches,"
    )

    print(
        "  2. absolute-value matches,"
    )

    print(
        "  3. post-boundary X-layer structure,"
    )

    print(
        "  4. local recurrence between successive X-layers."
    )

    print()
    print(
        "A successful pattern must still be interpreted from the exact"
    )

    print(
        "pq construction rather than accepted merely because it matches"
    )

    print(
        "a finite anchor set."
    )


# ============================================================================
# MAIN
# ============================================================================

def main():

    print("=" * 78)
    print("EXPERIMENT 205")
    print("BOUNDARY-ADAPTED SHIFTED-S COEFFICIENT AUDIT")
    print("=" * 78)

    print()
    print(
        "No previous experiment output is read."
    )

    print(
        "No unrestricted affine candidate enumeration is performed."
    )

    reconstruction_audit()
    parity_audit()
    print_shifted_rows()
    anchor_shifted_search()
    anchor_absolute_search()
    boundary_layer_audit()
    layer_recurrence_audit()
    c0_higher_layer_audit()
    final_diagnostic()

    print()
    print("=" * 78)
    print("END EXPERIMENT 205")
    print("=" * 78)


if __name__ == "__main__":
    main()