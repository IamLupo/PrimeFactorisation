# ==============================================================================
# EXPERIMENT 240
# EXACT KERNEL -> COMPLETE L2 BOUNDARY + L3 CONSISTENCY AUDIT
# ==============================================================================
#
# Standalone main.py
# Exact arithmetic over QQ
# No previous experiment imported
# No filesystem access
#
# Mathematical source:
#
#   F_{k,ell}(p,q) =
#       p^k (1+q)^ell
#     + q^k (1+p)^ell
#     - p^ell (1+q)^k
#     - q^ell (1+p)^k
#
# Transform:
#
#   N = p q
#   X = p + q + 1
#
# The experiment:
#
#   1. reconstruct exact G(N,X) directly from the pq kernel;
#   2. extract homogeneous layers L0,L1,L2,L3;
#   3. verify the complete interior product law for r=1,2,3;
#   4. derive/test the corrected L2 boundary rows:
#
#        a = k:
#          correction = -C(k+1,3)
#
#        a = k+1:
#          correction = C(k+1,2) (ell-k-2)
#
#        a >= k+2:
#          coefficient = 0
#
#   5. verify the known L3 boundary laws;
#   6. test k=13 explicitly, removing the previous false anomaly;
#   7. reconstruct L2 exactly from interior + boundary pieces.
#
# No L4 analysis is performed.
# ==============================================================================

import sympy as sp


# ==============================================================================
# SYMBOLS
# ==============================================================================

p, q = sp.symbols("p q")
N, X = sp.symbols("N X")

K, A, L, R = sp.symbols("K A L R", integer=True, nonnegative=True)


# ==============================================================================
# BASIC EXACT UTILITIES
# ==============================================================================

def simp(expr):
    return sp.factor(sp.cancel(sp.expand(expr)))


def exact_int(expr):
    """
    Convert an exact SymPy integer/rational to a clean SymPy value.
    """
    return sp.Integer(expr) if isinstance(expr, int) else sp.simplify(expr)


# ==============================================================================
# EXACT pq KERNEL
# ==============================================================================

def exact_F(k, ell):
    """
    Exact pq kernel supplied by the established construction.
    """
    return sp.expand(
        p**k * (1 + q)**ell
        + q**k * (1 + p)**ell
        - p**ell * (1 + q)**k
        - q**ell * (1 + p)**k
    )


# ==============================================================================
# EXACT SYMMETRIC TRANSFORMATION
# ==============================================================================

def exact_G(k, ell):
    """
    Convert the exact symmetric pq kernel to G(N,X), where

        N = pq
        X = p+q+1.

    SymPy's formal symmetric reduction is used only as algebra,
    not as imported project functionality.
    """
    F = exact_F(k, ell)

    sym_expr, remainder, mapping = sp.symmetrize(
        F,
        [p, q],
        formal=True,
    )

    if sp.expand(remainder) != 0:
        raise ValueError(
            f"Kernel is not symmetric after expansion: remainder={remainder}"
        )

    # mapping contains:
    #   (s1, p+q)
    #   (s2, p*q)
    s1 = mapping[0][0]
    s2 = mapping[1][0]

    G = sym_expr.subs(
        {
            s1: X - 1,
            s2: N,
        }
    )

    return sp.expand(G)


# ==============================================================================
# HOMOGENEOUS LAYER EXTRACTION
# ==============================================================================

def homogeneous_part(expr, degree):
    """
    Extract all terms N^a X^b with a+b=degree.
    """
    poly = sp.Poly(sp.expand(expr), N, X)

    result = sp.Integer(0)

    for (a, b), coeff in poly.terms():
        if a + b == degree:
            result += coeff * N**a * X**b

    return sp.expand(result)


def exact_layers(k, ell):
    """
    Return G and the first four homogeneous layers:

      L0 = degree ell
      L1 = degree ell-1
      L2 = degree ell-2
      L3 = degree ell-3
    """
    G = exact_G(k, ell)

    L0 = homogeneous_part(G, ell)
    L1 = homogeneous_part(G, ell - 1)
    L2 = homogeneous_part(G, ell - 2)
    L3 = homogeneous_part(G, ell - 3)

    return G, L0, L1, L2, L3


# ==============================================================================
# COEFFICIENT EXTRACTION
# ==============================================================================

def coeff_layer(layer, a, degree):
    """
    Coefficient of N^a X^(degree-a).
    """
    b = degree - a

    if a < 0 or b < 0:
        return sp.Integer(0)

    return sp.expand(layer).coeff(N, a).coeff(X, b)


# ==============================================================================
# ESTABLISHED INTERIOR PRODUCT LAW
# ==============================================================================

def P_interior(r, k, a, ell):
    """
    Established exact interior law:

      P_r(k,a,ell)
        = (-1)^(r+1)/r!
          * C(k+r,a)
          * product_{j=1}^{r-1}(ell-a-j)
          * ((k+r)ell-ka)/(k+r)

    Validity domain used here:
        0 <= a < k
    """
    r = sp.Integer(r)
    k = sp.sympify(k)
    a = sp.sympify(a)
    ell = sp.sympify(ell)

    product_part = sp.Integer(1)

    for j in range(1, int(r)):
        product_part *= ell - a - j

    return simp(
        (-1) ** (int(r) + 1)
        / sp.factorial(r)
        * sp.binomial(k + r, a)
        * product_part
        * ((k + r) * ell - k * a)
        / (k + r)
    )


# ==============================================================================
# CORRECTED L2 BOUNDARY FORMULAS
# ==============================================================================

def L2_corrected_coefficient(k, a, ell):
    """
    Complete proposed L2 coefficient:

      a < k:
        interior P2

      a = k:
        interior + (-C(k+1,3))

      a = k+1:
        interior + C(k+1,2)*(ell-k-2)

      a >= k+2:
        0
    """
    k = sp.sympify(k)
    a = sp.sympify(a)
    ell = sp.sympify(ell)

    if a < k:
        return simp(P_interior(2, k, a, ell))

    interior = P_interior(2, k, a, ell)

    if a == k:
        correction = -sp.binomial(k + 1, 3)
        return simp(interior + correction)

    if a == k + 1:
        correction = sp.binomial(k + 1, 2) * (ell - k - 2)
        return simp(interior + correction)

    return sp.Integer(0)


# ==============================================================================
# KNOWN L3 BOUNDARY FORMULAS
# ==============================================================================

def L3_known_boundary(k, a, ell):
    """
    Known exact L3 boundary structure.

      a = k:
        correction = C(k+2,3)

      a = k+1:
        correction =
          -C(k+2,2)*ell
          + k(k+2)(k+3)/2

      a = k+2:
        correction =
          (ell-k-4) *
          ((2k-1)ell-k(k+3))/2

    a >= k+3 is not asserted here.
    """
    k = sp.sympify(k)
    a = sp.sympify(a)
    ell = sp.sympify(ell)

    interior = P_interior(3, k, a, ell)

    if a == k:
        return simp(
            interior
            + sp.binomial(k + 2, 3)
        )

    if a == k + 1:
        correction = (
            -sp.binomial(k + 2, 2) * ell
            + k * (k + 2) * (k + 3) / 2
        )
        return simp(interior + correction)

    if a == k + 2:
        correction = (
            (ell - k - 4)
            * ((2 * k - 1) * ell - k * (k + 3))
            / 2
        )
        return simp(interior + correction)

    return None


# ==============================================================================
# 1. KERNEL SANITY AUDIT
# ==============================================================================

def kernel_sanity_audit():
    print("=" * 78)
    print("1. EXACT pq KERNEL SANITY AUDIT")
    print("=" * 78)

    cases = [
        (1, 3),
        (1, 4),
        (3, 7),
        (3, 9),
        (5, 11),
        (7, 15),
        (9, 21),
        (11, 23),
        (13, 25),
    ]

    failures = 0

    for k, ell in cases:
        F = exact_F(k, ell)

        symmetric = sp.expand(
            F - F.xreplace({p: q, q: p})
        ) == 0

        # Correct sanity condition:
        # F(p,p) is NOT expected to vanish when k != ell.
        diagonal = sp.factor(F.subs(q, p))

        expected_diagonal = sp.expand(
            2 * p**k * (1 + p)**ell
            - 2 * p**ell * (1 + p)**k
        )

        diagonal_ok = sp.expand(
            diagonal - expected_diagonal
        ) == 0

        ok = symmetric and diagonal_ok

        print(
            f"k={k:2d} ell={ell:2d} "
            f"symmetric={symmetric} "
            f"diagonal-form={diagonal_ok} "
            f"{'PASS' if ok else 'FAIL'}"
        )

        if not ok:
            failures += 1

    print()
    print(f"kernel sanity failures = {failures}")
    print()


# ==============================================================================
# 2. TOP-LAYER AUDIT
# ==============================================================================

def top_layer_audit():
    print("=" * 78)
    print("2. TOP-LAYER AUDIT")
    print("=" * 78)

    cases = [
        (1, 3),
        (1, 10),
        (3, 7),
        (3, 20),
        (5, 11),
        (5, 19),
        (7, 15),
        (9, 21),
        (11, 23),
        (13, 25),
    ]

    failures = 0

    for k, ell in cases:
        _, L0, _, _, _ = exact_layers(k, ell)

        expected = (
            -X**(ell - k)
            * ((X + N)**k - N**k)
        )

        residual = simp(L0 - expected)

        ok = residual == 0

        print(
            f"k={k:2d} ell={ell:2d} "
            f"residual={str(residual):>4} "
            f"{'PASS' if ok else 'FAIL'}"
        )

        if not ok:
            failures += 1

    print()
    print(f"top-layer failures = {failures}")
    print()


# ==============================================================================
# 3. COMPLETE INTERIOR L1/L2/L3 AUDIT
# ==============================================================================

def interior_product_audit():
    print("=" * 78)
    print("3. EXACT INTERIOR L1/L2/L3 PRODUCT-LAW AUDIT")
    print("=" * 78)

    cases = [
        (3, 7),
        (3, 9),
        (3, 11),
        (3, 13),
        (5, 11),
        (5, 13),
        (5, 15),
        (5, 17),
        (7, 15),
        (7, 17),
        (9, 21),
        (9, 23),
        (11, 23),
        (13, 25),
    ]

    total = 0
    failures = 0

    for k, ell in cases:
        _, _, L1, L2, L3 = exact_layers(k, ell)

        local_failures = 0

        for r, layer in [
            (1, L1),
            (2, L2),
            (3, L3),
        ]:
            degree = ell - r

            for a in range(k):
                actual = coeff_layer(layer, a, degree)
                expected = simp(P_interior(r, k, a, ell))

                total += 1

                if simp(actual - expected) != 0:
                    local_failures += 1
                    failures += 1

        print(
            f"k={k:2d} ell={ell:2d} "
            f"interior failures={local_failures} "
            f"{'PASS' if local_failures == 0 else 'FAIL'}"
        )

    print()
    print(f"tested interior coefficients = {total}")
    print(f"interior formula failures = {failures}")
    print()


# ==============================================================================
# 4. COMPLETE L2 BOUNDARY AUDIT
# ==============================================================================

def l2_boundary_audit():
    print("=" * 78)
    print("4. COMPLETE L2 BOUNDARY AUDIT")
    print("=" * 78)

    cases = [
        (3, 9),
        (3, 11),
        (3, 13),
        (5, 13),
        (5, 15),
        (5, 17),
        (7, 17),
        (7, 19),
        (9, 21),
        (9, 23),
        (11, 23),
        (11, 25),
        (13, 25),
        (13, 27),
        (13, 29),
    ]

    failures = 0

    for k, ell in cases:
        _, _, _, L2, _ = exact_layers(k, ell)

        print(f"k={k:2d} ell={ell:2d}")

        local = 0

        for a in [k, k + 1, k + 2]:
            degree = ell - 2

            actual = coeff_layer(L2, a, degree)
            expected = L2_corrected_coefficient(k, a, ell)
            residual = simp(actual - expected)

            ok = residual == 0

            print(
                f"  a={a:2d} "
                f"actual={str(actual):>8} "
                f"predicted={str(expected):>8} "
                f"residual={str(residual):>6} "
                f"{'PASS' if ok else 'FAIL'}"
            )

            if not ok:
                failures += 1
                local += 1

        print(
            f"  local failures = {local}"
        )
        print()

    print(f"L2 boundary failures = {failures}")
    print()


# ==============================================================================
# 5. L2 SUPPORT AUDIT
# ==============================================================================

def l2_support_audit():
    print("=" * 78)
    print("5. L2 SUPPORT AUDIT")
    print("=" * 78)

    cases = [
        (3, 9),
        (5, 13),
        (7, 17),
        (9, 21),
        (11, 23),
        (13, 25),
        (13, 27),
        (13, 29),
    ]

    failures = 0

    for k, ell in cases:
        _, _, _, L2, _ = exact_layers(k, ell)
        degree = ell - 2

        poly = sp.Poly(sp.expand(L2), N, X)

        support = sorted(
            a
            for (a, b), coeff in poly.terms()
            if a + b == degree and coeff != 0
        )

        expected_support = list(range(k + 2))

        ok = support == expected_support

        print(
            f"k={k:2d} ell={ell:2d} "
            f"support={support} "
            f"expected={expected_support} "
            f"{'PASS' if ok else 'FAIL'}"
        )

        if not ok:
            failures += 1

    print()
    print(f"L2 support failures = {failures}")
    print()


# ==============================================================================
# 6. EXPLICIT L2 BOUNDARY CORRECTION EXTRACTION
# ==============================================================================

def l2_correction_extraction():
    print("=" * 78)
    print("6. EXACT L2 BOUNDARY CORRECTION EXTRACTION")
    print("=" * 78)

    cases = [
        (3, 9),
        (5, 13),
        (7, 17),
        (9, 21),
        (11, 23),
        (13, 25),
        (13, 27),
        (13, 29),
    ]

    failures = 0

    for k, ell in cases:
        _, _, _, L2, _ = exact_layers(k, ell)

        # a=k
        actual_k = coeff_layer(L2, k, ell - 2)
        interior_k = P_interior(2, k, k, ell)
        delta_k = simp(actual_k - interior_k)

        expected_k = -sp.binomial(k + 1, 3)

        # a=k+1
        actual_k1 = coeff_layer(L2, k + 1, ell - 2)
        interior_k1 = P_interior(2, k, k + 1, ell)
        delta_k1 = simp(actual_k1 - interior_k1)

        expected_k1 = (
            sp.binomial(k + 1, 2)
            * (ell - k - 2)
        )

        ok1 = simp(delta_k - expected_k) == 0
        ok2 = simp(delta_k1 - expected_k1) == 0

        print(f"k={k:2d} ell={ell:2d}")
        print(
            f"  a=k   delta={str(delta_k):>8} "
            f"expected={str(expected_k):>8} "
            f"{'PASS' if ok1 else 'FAIL'}"
        )
        print(
            f"  a=k+1 delta={str(delta_k1):>8} "
            f"expected={str(expected_k1):>8} "
            f"{'PASS' if ok2 else 'FAIL'}"
        )

        if not ok1 or not ok2:
            failures += 1

        print()

    print(f"L2 correction failures = {failures}")
    print()


# ==============================================================================
# 7. L2 FULL RECONSTRUCTION
# ==============================================================================

def l2_full_reconstruction():
    print("=" * 78)
    print("7. COMPLETE L2 RECONSTRUCTION")
    print("=" * 78)

    cases = [
        (3, 7),
        (3, 9),
        (3, 11),
        (5, 11),
        (5, 13),
        (5, 15),
        (7, 15),
        (7, 17),
        (9, 21),
        (11, 23),
        (13, 25),
        (13, 27),
        (13, 29),
    ]

    failures = 0

    for k, ell in cases:
        _, _, _, L2, _ = exact_layers(k, ell)

        degree = ell - 2
        predicted = sp.Integer(0)

        # Entire degree-(ell-2) layer.
        for a in range(k + 2):
            coeff = L2_corrected_coefficient(k, a, ell)
            predicted += coeff * N**a * X**(degree - a)

        predicted = sp.expand(predicted)

        residual = simp(L2 - predicted)

        ok = residual == 0

        print(
            f"k={k:2d} ell={ell:2d} "
            f"reconstruction={'PASS' if ok else 'FAIL'}"
        )

        if not ok:
            print(f"  residual = {residual}")
            failures += 1

    print()
    print(f"L2 reconstruction failures = {failures}")
    print()


# ==============================================================================
# 8. L3 KNOWN BOUNDARY AUDIT
# ==============================================================================

def l3_boundary_audit():
    print("=" * 78)
    print("8. L3 KNOWN BOUNDARY AUDIT")
    print("=" * 78)

    cases = [
        (3, 9),
        (3, 11),
        (5, 13),
        (5, 15),
        (7, 17),
        (9, 21),
        (11, 23),
        (13, 25),
        (13, 27),
        (13, 29),
    ]

    failures = 0

    for k, ell in cases:
        _, _, _, _, L3 = exact_layers(k, ell)
        degree = ell - 3

        print(f"k={k:2d} ell={ell:2d}")

        for a in [k, k + 1, k + 2]:
            actual = coeff_layer(L3, a, degree)
            expected = L3_known_boundary(k, a, ell)

            residual = simp(actual - expected)

            ok = residual == 0

            print(
                f"  a={a:2d} "
                f"actual={str(actual):>8} "
                f"predicted={str(expected):>8} "
                f"residual={str(residual):>6} "
                f"{'PASS' if ok else 'FAIL'}"
            )

            if not ok:
                failures += 1

        print()

    print(f"L3 boundary failures = {failures}")
    print()


# ==============================================================================
# 9. K=13 LOCAL AUDIT
# ==============================================================================

def k13_local_audit():
    print("=" * 78)
    print("9. k=13 LOCAL BOUNDARY AUDIT")
    print("=" * 78)

    for ell in [25, 27, 29, 31]:
        _, _, _, L2, L3 = exact_layers(13, ell)

        print(f"k=13 ell={ell}")

        print("  L2:")
        for a in [13, 14, 15]:
            actual = coeff_layer(L2, a, ell - 2)
            expected = L2_corrected_coefficient(13, a, ell)
            residual = simp(actual - expected)

            print(
                f"    a={a:2d} "
                f"actual={str(actual):>8} "
                f"expected={str(expected):>8} "
                f"residual={str(residual):>6}"
            )

        print("  L3:")
        for a in [13, 14, 15]:
            actual = coeff_layer(L3, a, ell - 3)
            expected = L3_known_boundary(13, a, ell)
            residual = simp(actual - expected)

            print(
                f"    a={a:2d} "
                f"actual={str(actual):>8} "
                f"expected={str(expected):>8} "
                f"residual={str(residual):>6}"
            )

        print()


# ==============================================================================
# 10. SYMBOLIC BOUNDARY FORMULAS
# ==============================================================================

def symbolic_boundary_audit():
    print("=" * 78)
    print("10. SYMBOLIC L2 BOUNDARY FORMULAS")
    print("=" * 78)

    l2_k = simp(
        P_interior(2, K, K, L)
        - sp.binomial(K + 1, 3)
    )

    l2_k1 = simp(
        P_interior(2, K, K + 1, L)
        + sp.binomial(K + 1, 2) * (L - K - 2)
    )

    print("L2 a=k:")
    print(f"  {l2_k}")
    print()

    print("L2 a=k+1:")
    print(f"  {l2_k1}")
    print()

    print("L2 a>=k+2:")
    print("  coefficient = 0")
    print()

    print("Corrections:")
    print("  Delta_2(k)   =", -sp.binomial(K + 1, 3))
    print(
        "  Delta_2(k+1) =",
        sp.binomial(K + 1, 2) * (L - K - 2),
    )
    print()


# ==============================================================================
# 11. FINAL SUMMARY
# ==============================================================================

def final_summary():
    print("=" * 78)
    print("11. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        """
The exact pq kernel is now used as the sole mathematical source.

The experiment tests the corrected L2 boundary structure:

  a < k:
    P2 =
      -1/2 * C(k+2,a)
      * (L-a-1)
      * ((k+2)L-ka)/(k+2)

  a = k:
    interior + Delta_2(k)

    Delta_2(k) = -C(k+1,3)

  a = k+1:
    interior + Delta_2(k+1)

    Delta_2(k+1)
      = C(k+1,2)(L-k-2)

  a >= k+2:
    coefficient = 0

The experiment independently reconstructs the complete
degree-(ell-2) layer from these formulas.

It also checks that the established L3 boundary formulas
continue to hold at k=13 and fresh ell-values.

If all sections pass, the corrected L2 boundary is closed.

Only after that should the project move to L4.
"""
    )


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    print("=" * 78)
    print("EXPERIMENT 240")
    print("EXACT KERNEL -> COMPLETE L2 BOUNDARY + L3 CONSISTENCY AUDIT")
    print("=" * 78)
    print()
    print("Standalone main.py")
    print("Exact arithmetic over QQ")
    print("No previous experiment imported")
    print("No filesystem access")
    print("No L4 analysis")
    print()

    kernel_sanity_audit()
    top_layer_audit()
    interior_product_audit()
    l2_boundary_audit()
    l2_support_audit()
    l2_correction_extraction()
    l2_full_reconstruction()
    l3_boundary_audit()
    k13_local_audit()
    symbolic_boundary_audit()
    final_summary()


if __name__ == "__main__":
    main()

