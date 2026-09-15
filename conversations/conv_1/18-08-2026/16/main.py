import sympy as sp


# ============================================================================
# EXPERIMENT 228
# NORMALIZED LAYER POLYNOMIAL / GENERATING-MECHANISM AUDIT
# ============================================================================
#
# Standalone main.py
# Exact arithmetic over QQ
# No previous experiment imported
# No filesystem access
# No exact_F required
# No unrestricted candidate enumeration
#
# Established interior layers used as exact mathematical input:
#
# L0:
#   P0(k,a) = -C(k,a)
#
# L1:
#   P1(k,a,L)
#     = L*C(k+1,a) - k*C(k,a-1)
#
# L2:
#   P2(k,a,L)
#     = -1/2*C(k+2,a)*L^2
#       + C(k+2,a)*(2(k+1)a+k+2)/(2(k+2))*L
#       - k*a(a+1)*C(k+2,a)/(2(k+2))
#
# L3:
#   P3(k,a,L)
#     = C(k+3,a)/6 * L^3
#       - C(k+3,a)*(a(k+2)+k+3)/(2(k+3))*L^2
#       + C(k+3,a)*
#           (3a^2 k + 3a^2 + 6ak + 9a + 2k + 6)
#           /(6(k+3))*L
#       - k*a(a+1)(a+2)*C(k+3,a)/(6(k+3))
#
# The experiment asks:
#
#   Q1 = P1 / C(k+1,a)
#   Q2 = -2 P2 / C(k+2,a)
#   Q3 =  6 P3 / C(k+3,a)
#
# Do Q1,Q2,Q3 reveal a common multiplicative structure?
#
# ============================================================================


def simp(expr):
    return sp.factor(sp.cancel(sp.expand(expr)))


def C(n, r):
    """
    Exact binomial coefficient.
    Return 0 outside the ordinary combinatorial range.
    """
    if r < 0 or r > n:
        return sp.Integer(0)
    return sp.binomial(n, r)


def P0(k, a):
    return -C(k, a)


def P1(k, a, L):
    return simp(
        L * C(k + 1, a)
        - k * C(k, a - 1)
    )


def P2(k, a, L):
    K = sp.Integer(k)
    A = sp.Integer(a)

    return simp(
        -sp.Rational(1, 2) * C(k + 2, a) * L**2
        + C(k + 2, a)
        * (2 * (K + 1) * A + K + 2)
        / (2 * (K + 2))
        * L
        - K * A * (A + 1)
        * C(k + 2, a)
        / (2 * (K + 2))
    )


def P3(k, a, L):
    K = sp.Integer(k)
    A = sp.Integer(a)

    return simp(
        C(k + 3, a) * L**3 / 6
        - C(k + 3, a)
        * (A * (K + 2) + K + 3)
        / (2 * (K + 3))
        * L**2
        + C(k + 3, a)
        * (
            3 * A**2 * K
            + 3 * A**2
            + 6 * A * K
            + 9 * A
            + 2 * K
            + 6
        )
        / (6 * (K + 3))
        * L
        - K * A * (A + 1) * (A + 2)
        * C(k + 3, a)
        / (6 * (K + 3))
    )


# ============================================================================
# 1. LEADING BINOMIAL LAW
# ============================================================================

def leading_law_audit():
    print("=" * 78)
    print("1. LEADING BINOMIAL LAW")
    print("=" * 78)

    failures = 0

    samples = [
        (3, range(0, 4)),
        (5, range(0, 6)),
        (7, range(0, 8)),
        (9, range(0, 8)),
        (11, range(0, 8)),
        (13, range(0, 8)),
    ]

    for k, arange in samples:
        for a in arange:
            p1 = sp.Poly(P1(k, a, sp.Symbol("L")), sp.Symbol("L"))
            p2 = sp.Poly(P2(k, a, sp.Symbol("L")), sp.Symbol("L"))
            p3 = sp.Poly(P3(k, a, sp.Symbol("L")), sp.Symbol("L"))

            c1 = p1.coeff_monomial(sp.Symbol("L"))
            c2 = p2.coeff_monomial(sp.Symbol("L") ** 2)
            c3 = p3.coeff_monomial(sp.Symbol("L") ** 3)

            e1 = C(k + 1, a)
            e2 = -sp.Rational(1, 2) * C(k + 2, a)
            e3 = sp.Rational(1, 6) * C(k + 3, a)

            ok = (
                sp.simplify(c1 - e1) == 0
                and sp.simplify(c2 - e2) == 0
                and sp.simplify(c3 - e3) == 0
            )

            if not ok:
                failures += 1
                print(
                    f"FAIL k={k} a={a}"
                )

    print(f"leading-law failures = {failures}")
    print()


# ============================================================================
# 2. NORMALIZED POLYNOMIALS
# ============================================================================

def normalized_polynomials():
    print("=" * 78)
    print("2. NORMALIZED LAYER POLYNOMIALS")
    print("=" * 78)

    L = sp.Symbol("L")

    samples = [
        (3, 1),
        (3, 2),
        (5, 1),
        (5, 2),
        (5, 3),
        (7, 2),
        (7, 4),
        (9, 3),
        (9, 4),
        (11, 5),
    ]

    for k, a in samples:
        q1 = simp(P1(k, a, L) / C(k + 1, a))
        q2 = simp(-2 * P2(k, a, L) / C(k + 2, a))
        q3 = simp(6 * P3(k, a, L) / C(k + 3, a))

        print(f"k={k:2d} a={a:2d}")
        print(f"  Q1 = {q1}")
        print(f"  Q2 = {q2}")
        print(f"  Q3 = {q3}")
        print()


# ============================================================================
# 3. EXACT FACTORIZATION AUDIT
# ============================================================================

def factorization_audit():
    print("=" * 78)
    print("3. EXACT FACTORIZATION AUDIT")
    print("=" * 78)

    L = sp.Symbol("L")

    samples = [
        (3, 1),
        (3, 2),
        (5, 2),
        (5, 3),
        (7, 2),
        (7, 4),
        (9, 3),
        (9, 4),
        (11, 5),
    ]

    for k, a in samples:
        q1 = simp(P1(k, a, L) / C(k + 1, a))
        q2 = simp(-2 * P2(k, a, L) / C(k + 2, a))
        q3 = simp(6 * P3(k, a, L) / C(k + 3, a))

        print(f"k={k:2d} a={a:2d}")
        print(f"  factor(Q1) = {sp.factor(q1)}")
        print(f"  factor(Q2) = {sp.factor(q2)}")
        print(f"  factor(Q3) = {sp.factor(q3)}")

        roots1 = sp.solve(sp.Eq(q1, 0), L)
        roots2 = sp.solve(sp.Eq(q2, 0), L)
        roots3 = sp.solve(sp.Eq(q3, 0), L)

        print(f"  roots(Q1)  = {roots1}")
        print(f"  roots(Q2)  = {roots2}")
        print(f"  roots(Q3)  = {roots3}")
        print()


# ============================================================================
# 4. VIETA STRUCTURE
#
# For monic normalized polynomials:
#
# Q1 = L - R1
# Q2 = L^2 - E1 L + E2
# Q3 = L^3 - E1 L^2 + E2 L - E3
#
# The experiment checks whether the coefficients behave like elementary
# symmetric functions of a coherent collection of hidden shifts.
# ============================================================================

def vieta_audit():
    print("=" * 78)
    print("4. VIETA / HIDDEN-SHIFT STRUCTURE")
    print("=" * 78)

    L = sp.Symbol("L")

    samples = [
        (3, 1),
        (3, 2),
        (5, 2),
        (5, 3),
        (7, 3),
        (7, 4),
        (9, 4),
        (11, 5),
    ]

    for k, a in samples:
        q1 = sp.Poly(
            simp(P1(k, a, L) / C(k + 1, a)),
            L
        )

        q2 = sp.Poly(
            simp(-2 * P2(k, a, L) / C(k + 2, a)),
            L
        )

        q3 = sp.Poly(
            simp(6 * P3(k, a, L) / C(k + 3, a)),
            L
        )

        r1 = simp(-q1.all_coeffs()[1])

        coeff2 = q2.all_coeffs()
        e1_2 = simp(-coeff2[1])
        e2_2 = simp(coeff2[2])

        coeff3 = q3.all_coeffs()
        e1_3 = simp(-coeff3[1])
        e2_3 = simp(coeff3[2])
        e3_3 = simp(-coeff3[3])

        print(f"k={k:2d} a={a:2d}")
        print(f"  Q1 root-sum      = {r1}")
        print(f"  Q2 E1           = {e1_2}")
        print(f"  Q2 E2           = {e2_2}")
        print(f"  Q3 E1           = {e1_3}")
        print(f"  Q3 E2           = {e2_3}")
        print(f"  Q3 E3           = {e3_3}")
        print(
            "  E1(Q1,Q2,Q3) consistency = "
            + ("YES" if sp.simplify(r1 - e1_2) == 0
               and sp.simplify(r1 - e1_3) == 0 else "NO")
        )
        print()


# ============================================================================
# 5. GENERATING-BASIS TEST
#
# Test the natural candidate:
#
#   H_r = (-1)^(r+1) / r! * C(k+r,a) * Q_r
#
# with Q_r normalized above.
#
# Instead of fitting anything, we inspect whether Q_r itself admits
# a common product structure.
# ============================================================================

def generating_basis_audit():
    print("=" * 78)
    print("5. GENERATING-BASIS AUDIT")
    print("=" * 78)

    L = sp.Symbol("L")

    samples = [
        (3, 1),
        (3, 2),
        (5, 2),
        (5, 3),
        (7, 3),
        (7, 4),
        (9, 4),
        (11, 5),
    ]

    for k, a in samples:
        q1 = simp(P1(k, a, L) / C(k + 1, a))
        q2 = simp(-2 * P2(k, a, L) / C(k + 2, a))
        q3 = simp(6 * P3(k, a, L) / C(k + 3, a))

        print(f"k={k:2d} a={a:2d}")

        # Polynomial discriminant information is exact and may reveal
        # whether Q2/Q3 have unusually simple root structures.
        disc2 = simp(sp.discriminant(q2, L))
        disc3 = simp(sp.discriminant(q3, L))

        print(f"  Q1 = {sp.factor(q1)}")
        print(f"  Q2 = {sp.factor(q2)}")
        print(f"  disc(Q2) = {disc2}")
        print(f"  Q3 = {sp.factor(q3)}")
        print(f"  disc(Q3) = {disc3}")
        print()


# ============================================================================
# 6. SYMBOLIC COEFFICIENT IDENTITIES
# ============================================================================

def symbolic_identity_audit():
    print("=" * 78)
    print("6. SYMBOLIC NORMALIZED IDENTITIES")
    print("=" * 78)

    K, A, L = sp.symbols("K A L")

    # Work symbolically with binomial symbols represented as B1/B2/B3.
    B1 = sp.Symbol("B1")
    B2 = sp.Symbol("B2")
    B3 = sp.Symbol("B3")

    # Exact normalized expressions already established in the interior.
    Q1 = simp(
        L - K * A / (K + 1)
    )

    Q2 = simp(
        L**2
        - (2 * (K + 1) * A + K + 2) / (K + 2) * L
        + K * A * (A + 1) / (K + 2)
    )

    Q3 = simp(
        L**3
        - 3 * (A * (K + 2) + K + 3) / (K + 3) * L**2
        + (
            3 * A**2 * K
            + 3 * A**2
            + 6 * A * K
            + 9 * A
            + 2 * K
            + 6
        ) / (K + 3) * L
        - K * A * (A + 1) * (A + 2) / (K + 3)
    )

    print("Q1 =")
    print(sp.factor(Q1))
    print()

    print("Q2 =")
    print(sp.factor(Q2))
    print()

    print("Q3 =")
    print(sp.factor(Q3))
    print()

    # Compare the first coefficient at each order.
    expected_q1 = -K * A / (K + 1)
    expected_q2 = -(2 * (K + 1) * A + K + 2) / (K + 2)
    expected_q3 = -3 * (A * (K + 2) + K + 3) / (K + 3)

    print(
        "Q1 first-root coefficient identity = "
        + (
            "PASS"
            if sp.simplify(Q1.coeff(L, 0) - expected_q1) == 0
            else "FAIL"
        )
    )

    print(
        "Q2 first-root coefficient identity = "
        + (
            "PASS"
            if sp.simplify(Q2.coeff(L, 1) - expected_q2) == 0
            else "FAIL"
        )
    )

    print(
        "Q3 first-root coefficient identity = "
        + (
            "PASS"
            if sp.simplify(Q3.coeff(L, 2) - expected_q3) == 0
            else "FAIL"
        )
    )

    print()


# ============================================================================
# 7. DIRECT SIGNED LEADING-LAYER AUDIT
#
# This explicitly corrects the sign mistake detected in Experiment 227.
# ============================================================================

def corrected_leading_layer_audit():
    print("=" * 78)
    print("7. CORRECTED SIGNED LEADING-LAYER LAW")
    print("=" * 78)

    failures = 0
    L = sp.Symbol("L")

    samples = [
        (3, 0),
        (3, 1),
        (3, 2),
        (5, 0),
        (5, 1),
        (5, 2),
        (5, 3),
        (7, 0),
        (7, 1),
        (7, 2),
        (7, 3),
        (7, 4),
        (9, 0),
        (9, 1),
        (9, 2),
        (9, 3),
        (9, 4),
        (11, 0),
        (11, 1),
        (11, 2),
        (11, 3),
        (11, 4),
    ]

    polys = [None, P1, P2, P3]

    for k, a in samples:
        # r = 1
        p = sp.Poly(P1(k, a, L), L)
        actual = p.coeff_monomial(L)
        expected = C(k + 1, a)

        if sp.simplify(actual - expected) != 0:
            failures += 1
            print(f"FAIL r=1 k={k} a={a}")

        # r = 2
        p = sp.Poly(P2(k, a, L), L)
        actual = p.coeff_monomial(L**2)
        expected = -sp.Rational(1, 2) * C(k + 2, a)

        if sp.simplify(actual - expected) != 0:
            failures += 1
            print(f"FAIL r=2 k={k} a={a}")

        # r = 3
        p = sp.Poly(P3(k, a, L), L)
        actual = p.coeff_monomial(L**3)
        expected = sp.Rational(1, 6) * C(k + 3, a)

        if sp.simplify(actual - expected) != 0:
            failures += 1
            print(f"FAIL r=3 k={k} a={a}")

    print(f"corrected signed-law failures = {failures}")
    print()


# ============================================================================
# 8. FINAL DIAGNOSTIC
# ============================================================================

def final_diagnostic():
    print("=" * 78)
    print("8. FINAL DIAGNOSTIC")
    print("=" * 78)

    print(
        "The experiment does not attempt to reconstruct exact_F."
    )
    print()

    print(
        "It treats the already-derived interior L1/L2/L3 formulas "
        "as exact algebraic constraints."
    )
    print()

    print(
        "The key objects are the normalized polynomials:"
    )
    print(
        "  Q1 = P1 / C(k+1,a)"
    )
    print(
        "  Q2 = -2 P2 / C(k+2,a)"
    )
    print(
        "  Q3 =  6 P3 / C(k+3,a)"
    )
    print()

    print(
        "The central question is whether Q1,Q2,Q3 are successive "
        "members of one hidden product/generating mechanism."
    )
    print()

    print(
        "The corrected leading coefficient law is:"
    )
    print(
        "  coefficient of ell^r in P_r"
    )
    print(
        "    = (-1)^(r+1) / r! * C(k+r,a)"
    )
    print()

    print(
        "If the factorizations and Vieta data reveal a coherent "
        "shift sequence, the next experiment should derive that "
        "mechanism symbolically rather than fitting new layers."
    )
    print()

    print(
        "No boundary formula is accepted here."
    )
    print(
        "No L4 analysis is performed."
    )


def main():
    print("=" * 78)
    print("EXPERIMENT 228")
    print("NORMALIZED LAYER POLYNOMIAL / GENERATING-MECHANISM AUDIT")
    print("=" * 78)
    print()
    print("Standalone main.py")
    print("Exact arithmetic over QQ")
    print("No previous experiment imported")
    print("No filesystem access")
    print("No exact_F required")
    print()

    leading_law_audit()
    normalized_polynomials()
    factorization_audit()
    vieta_audit()
    generating_basis_audit()
    symbolic_identity_audit()
    corrected_leading_layer_audit()
    final_diagnostic()


if __name__ == "__main__":
    main()
