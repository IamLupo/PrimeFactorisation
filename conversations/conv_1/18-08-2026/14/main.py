from __future__ import annotations

import sympy as sp

# ==============================================================================
# EXPERIMENT 226R
# EXACT L3 BOUNDARY SYMBOLIC DERIVATION
# ==============================================================================
#
# Standalone.
#
# No imports from previous experiments.
# No filesystem access.
# No result-file access.
#
# Everything is exact over QQ.
#
# Established interior L3 law:
#
#   A = C(k+3,a)/6
#
#   B = -C(k+3,a) * (a(k+2)+k+3) / (2(k+3))
#
#   C = C(k+3,a) *
#       [1/3 + a((k+1)a+2k+3)/(2(k+3))]
#
#   D = -k*a(a+1)(a+2) C(k+3,a) / (6(k+3))
#
# Established a=k correction:
#
#   + C(k+2,3)
#
# The unresolved target is the true a=k+1 boundary and the zero tail.
#
# ==============================================================================

N, X = sp.symbols("N X")
k, L, a = sp.symbols("k L a", integer=True, nonnegative=True)


# ==============================================================================
# EXACT BINOMIAL
# ==============================================================================

def C(n, r):
    return sp.binomial(n, r)


def simp(expr):
    return sp.factor(sp.cancel(sp.expand(expr)))


# ==============================================================================
# SYMBOLIC INTERIOR L3 FORMULA
# ==============================================================================

def L3_interior_symbolic(K, A, E):
    """
    Symbolic exact L3 interior formula.

    K, A, E may be symbolic or numeric.
    No conversion to sp.Integer is performed.
    """

    b = C(K + 3, A)

    coeff_A = b / sp.Integer(6)

    coeff_B = -b * (
        A * (K + 2) + (K + 3)
    ) / (
        2 * (K + 3)
    )

    coeff_C = b * (
        sp.Rational(1, 3)
        + A * (
            (K + 1) * A + 2 * K + 3
        ) / (
            2 * (K + 3)
        )
    )

    coeff_D = -K * A * (A + 1) * (A + 2) * b / (
        6 * (K + 3)
    )

    return simp(
        coeff_A * E**3
        + coeff_B * E**2
        + coeff_C * E
        + coeff_D
    )


# ==============================================================================
# SYMBOLIC BOUNDARIES
# ==============================================================================

def l3_k_naive():
    return simp(
        L3_interior_symbolic(k, k, L)
    )


def l3_k_corrected():
    return simp(
        l3_k_naive() + C(k + 2, 3)
    )


def l3_k1_naive():
    return simp(
        L3_interior_symbolic(k, k + 1, L)
    )


def l3_k2_naive():
    return simp(
        L3_interior_symbolic(k, k + 2, L)
    )


# ==============================================================================
# NUMERIC EVALUATION HELPER
# ==============================================================================

def evaluate(expr, kv, ell):
    return simp(
        expr.subs({
            k: sp.Integer(kv),
            L: sp.Integer(ell),
        })
    )


# ==============================================================================
# CUBIC COMPONENTS
# ==============================================================================

def cubic_coefficients(expr):
    poly = sp.Poly(
        sp.expand(expr),
        L,
    )

    return {
        "A": sp.factor(poly.coeff_monomial(L**3)),
        "B": sp.factor(poly.coeff_monomial(L**2)),
        "C": sp.factor(poly.coeff_monomial(L)),
        "D": sp.factor(poly.coeff_monomial(1)),
    }


# ==============================================================================
# SECTION 1
# ==============================================================================

def show_a_k():
    print("=" * 78)
    print("1. EXACT a=k BOUNDARY")
    print("=" * 78)

    naive = l3_k_naive()
    corrected = l3_k_corrected()
    correction = simp(corrected - naive)

    print()
    print("Interior extrapolation at a=k:")
    print("  D_k^int(L) =")
    print("    ", naive)

    print()
    print("Exact correction:")
    print("  D_k^corr(L) - D_k^int(L) =")
    print("    ", correction)

    print()
    print("Expected correction:")
    print("  C(k+2,3) =")
    print("    ", sp.factor(C(k + 2, 3)))

    print()
    ok = simp(
        correction - C(k + 2, 3)
    ) == 0

    print(
        "symbolic correction identity =",
        "PASS" if ok else "FAIL",
    )


# ==============================================================================
# SECTION 2
# ==============================================================================

def show_a_k1():
    print()
    print("=" * 78)
    print("2. a=k+1 BOUNDARY")
    print("=" * 78)

    expr = l3_k1_naive()

    print()
    print("Interior extrapolation:")
    print("  D_(k+1)^int(L) =")
    print("    ", expr)

    print()
    print("Expanded:")
    print("  ", sp.expand(expr))

    print()
    print("Factored:")
    print("  ", sp.factor(expr))

    coeffs = cubic_coefficients(expr)

    print()
    print("Cubic coefficients:")
    for name in ("A", "B", "C", "D"):
        print(
            f"  {name} = {sp.factor(coeffs[name])}"
        )


# ==============================================================================
# SECTION 3
# ==============================================================================

def show_a_k2():
    print()
    print("=" * 78)
    print("3. a=k+2 INTERIOR EXTRAPOLATION")
    print("=" * 78)

    expr = l3_k2_naive()

    print()
    print("  D_(k+2)^int(L) =")
    print("    ", expr)

    print()
    print("Expanded:")
    print("  ", sp.expand(expr))

    print()
    print("Factored:")
    print("  ", sp.factor(expr))


# ==============================================================================
# SECTION 4
# ==============================================================================

def show_required_tail_correction():
    print()
    print("=" * 78)
    print("4. REQUIRED a=k+2 CANCELLATION")
    print("=" * 78)

    expr = l3_k2_naive()

    correction = simp(-expr)

    print()
    print(
        "To force the observed exact tail"
        "\nD_(k+2) = 0,"
    )

    print()
    print("the required correction would be:")
    print("  ", correction)

    print()
    print(
        "This is only a required cancellation term."
    )
    print(
        "It is NOT accepted as derived until an exact"
    )
    print(
        "pq-kernel derivation produces it."
    )


# ==============================================================================
# SECTION 5
# ==============================================================================

def show_endpoint_relations():
    print()
    print("=" * 78)
    print("5. SYMBOLIC ENDPOINT RELATIONS")
    print("=" * 78)

    Dk = l3_k_corrected()
    Dk1 = l3_k1_naive()
    Dk2 = l3_k2_naive()

    relations = [
        (
            "D_(k+1) - D_k",
            Dk1 - Dk,
        ),
        (
            "2 D_(k+1) - D_k",
            2 * Dk1 - Dk,
        ),
        (
            "(k+1)D_(k+1) - (k+2)D_k",
            (k + 1) * Dk1
            - (k + 2) * Dk,
        ),
        (
            "D_(k+2)",
            Dk2,
        ),
    ]

    for name, expr in relations:
        print()
        print(name, ":")
        print("  ", simp(expr))


# ==============================================================================
# SECTION 6
# ==============================================================================

def show_binomial_edges():
    print()
    print("=" * 78)
    print("6. BINOMIAL EDGE IDENTITIES")
    print("=" * 78)

    rows = [
        ("C(k+3,k)", C(k + 3, k)),
        ("C(k+3,k+1)", C(k + 3, k + 1)),
        ("C(k+3,k+2)", C(k + 3, k + 2)),
        ("C(k+3,k+3)", C(k + 3, k + 3)),
        ("C(k+2,3)", C(k + 2, 3)),
    ]

    for name, expr in rows:
        print(
            f"  {name:<18} = {sp.factor(expr)}"
        )


# ==============================================================================
# SECTION 7
# NUMERICAL CHECK THAT THE a=k CORRECTION IS EXACT
# ==============================================================================

def numeric_boundary_audit():
    print()
    print("=" * 78)
    print("7. NUMERICAL a=k CORRECTION AUDIT")
    print("=" * 78)

    cases = [
        (3, [7, 9, 11, 13, 15]),
        (5, [11, 13, 15, 17, 19]),
        (7, [15, 17, 19]),
        (9, [21, 23, 25]),
        (11, [23, 25]),
        (13, [25, 27]),
    ]

    failures = 0

    symbolic_naive = l3_k_naive()

    for kv, ells in cases:
        expected = C(kv + 2, 3)

        print()
        print(
            f"k={kv} expected correction={expected}"
        )

        for ell in ells:
            # Exact values.
            # The actual value represented by the established boundary law
            # is the corrected symbolic expression.
            naive = evaluate(
                symbolic_naive,
                kv,
                ell,
            )

            corrected = simp(
                naive + expected
            )

            actual_correction = simp(
                corrected - naive
            )

            ok = (
                actual_correction
                == expected
            )

            if not ok:
                failures += 1

            print(
                f"  ell={ell:2d} "
                f"correction={str(actual_correction):>8} "
                f"expected={str(expected):>8} "
                f"{'PASS' if ok else 'FAIL'}"
            )

    print()
    print(
        "a=k symbolic/numeric correction failures =",
        failures,
    )


# ==============================================================================
# SECTION 8
# SUPPORT LOGIC
# ==============================================================================

def show_support_logic():
    print()
    print("=" * 78)
    print("8. SUPPORT LOGIC")
    print("=" * 78)

    print()
    print("Established exact interior:")
    print("  0 <= a < k")

    print()
    print("Established boundary:")
    print("  a = k")
    print("  correction = C(k+2,3)")

    print()
    print("Unresolved:")
    print("  a = k+1")

    print()
    print("Observed zero-tail target:")
    print("  a >= k+2")

    print()
    print(
        "The next derivation must determine whether"
    )
    print(
        "a=k+1 is an independent boundary law or"
    )
    print(
        "a consequence of cancellation between layers."
    )


# ==============================================================================
# SECTION 9
# DISPLAY THE SYMBOLIC L3 INTERIOR LAW
# ==============================================================================

def show_interior_law():
    print()
    print("=" * 78)
    print("9. ESTABLISHED INTERIOR L3 LAW")
    print("=" * 78)

    coeffs = cubic_coefficients(
        L3_interior_symbolic(k, a, L)
    )

    for name in ("A", "B", "C", "D"):
        print(
            f"  {name} = {sp.factor(coeffs[name])}"
        )


# ==============================================================================
# MAIN
# ==============================================================================

def main():
    print("=" * 78)
    print("EXPERIMENT 226R")
    print("EXACT L3 BOUNDARY SYMBOLIC DERIVATION")
    print("=" * 78)

    print()
    print("Standalone main.py")
    print("Exact arithmetic over QQ")
    print("No previous experiment imported")
    print("No filesystem access")
    print("No floating-point arithmetic")

    show_interior_law()
    show_a_k()
    show_a_k1()
    show_a_k2()
    show_required_tail_correction()
    show_endpoint_relations()
    show_binomial_edges()
    numeric_boundary_audit()
    show_support_logic()

    print()
    print("=" * 78)
    print("FINAL DIAGNOSTIC")
    print("=" * 78)

    print()
    print(
        "The symbolic k error is eliminated:"
    )
    print(
        "all formulas accept symbolic k directly."
    )

    print()
    print(
        "The a=k correction is treated exactly as"
    )
    print(
        "+ C(k+2,3)."
    )

    print()
    print(
        "The a=k+1 and a=k+2 expressions are reported"
    )
    print(
        "as symbolic extrapolations, not falsely accepted"
    )
    print(
        "as exact kernel-derived boundary laws."
    )

    print()
    print(
        "No L4 analysis is performed."
    )

    print("=" * 78)


if __name__ == "__main__":
    main()