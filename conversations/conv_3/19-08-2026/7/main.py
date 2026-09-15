import sympy as sp


# =============================================================================
# EXPERIMENT 382 START
# RESULTANT COEFFICIENT / DISCRIMINANT STRUCTURE
# =============================================================================

print("=" * 110)
print("EXPERIMENT 382 START")
print("=" * 110)
print()
print("EXACT RESULTANT COEFFICIENT STRUCTURE")
print()
print("Goal:")
print("  1. Extract the resultant as A(K) R^2 + B(K) R + C(K).")
print("  2. Determine degree patterns of A, B, C.")
print("  3. Factor A, B, C as far as practical.")
print("  4. Analyze the quadratic discriminant B^2 - 4AC.")
print("  5. Compare the normalized structures across independent N.")
print()


# =============================================================================
# SYMBOLS
# =============================================================================

t = sp.symbols("t")
K = sp.symbols("K")
R = sp.symbols("R")


# =============================================================================
# EXACT LAYER POLYNOMIALS
# =============================================================================

def h16(t):
    return (
        -9*t**8
        -36*t**7
        -84*t**6
        -126*t**5
        -126*t**4
        -84*t**3
        -36*t**2
        -9*t
        -1
    )


def h15(t):
    return (
        88*t**9
        +396*t**8
        +1164*t**7
        +2226*t**6
        +2898*t**5
        +2604*t**4
        +1596*t**3
        +639*t**2
        +151*t
        +16
    )


def h14(t):
    return (
        -276*t**10
        -1380*t**9
        -5460*t**8
        -13560*t**7
        -23058*t**6
        -27510*t**5
        -23100*t**4
        -13410*t**3
        -5135*t**2
        -1169*t
        -120
    )


H16 = sp.expand(h16(t))
H15 = sp.expand(h15(t))
H14 = sp.expand(h14(t))

NUM_R = sp.expand(H15**2)
DEN_R = sp.expand(H16 * H14)

LAYER_POLY = sp.expand(
    R * DEN_R - NUM_R
)


# =============================================================================
# TEST INSTANCES
# =============================================================================

TEST_PAIRS = (
    (50387, 282589),
    (1009, 10007),
    (10007, 1000003),
    (10007, 10009),
    (50021, 50047),
    (100003, 100019),
    (200003, 200009),
    (300007, 900001),
    (500009, 700001),
    (1000003, 1000033),
    (2000003, 3000017),
)


# =============================================================================
# EXACT KAPPA
# =============================================================================

def kappa_from_x(N, X):
    S = X - 1

    return (
        -S**2
        + (N + 1)*S
        - N**2
        + N
    )


# =============================================================================
# COMPACT POLYNOMIAL SUMMARY
# =============================================================================

def summarize_poly(label, poly, variable):
    poly = sp.Poly(sp.expand(poly), variable)

    print(f"  {label}:")
    print(f"    degree     = {poly.degree()}")
    print(f"    terms      = {len(poly.terms())}")
    print(f"    LC         = {poly.LC()}")
    print(f"    TC         = {poly.TC()}")
    print(f"    content    = {poly.content()}")

    return poly


# =============================================================================
# PER-INSTANCE ANALYSIS
# =============================================================================

degree_signatures = []
discriminant_degree_signatures = []

for instance_id, (p, q) in enumerate(TEST_PAIRS, start=1):

    print()
    print("=" * 110)
    print(f"INSTANCE {instance_id}")
    print("=" * 110)

    N = p * q
    X_true = p + q + 1

    K_true = sp.Integer(
        kappa_from_x(N, X_true)
    )

    print()
    print("BASIC DATA")
    print("  N      =", N)
    print("  X_TRUE =", X_true)
    print("  K_TRUE =", K_true)

    # -------------------------------------------------------------------------
    # Correct Kappa polynomial
    # -------------------------------------------------------------------------

    kappa_poly = sp.expand(
        (K + N**2 + 2)*t**2
        - N*(N + 3)*t
        + N**2
    )

    # -------------------------------------------------------------------------
    # Exact resultant
    # -------------------------------------------------------------------------

    print()
    print("COMPUTING RESULTANT")

    resultant = sp.resultant(
        kappa_poly,
        LAYER_POLY,
        t
    )

    resultant = sp.Poly(
        sp.expand(resultant),
        R,
        K
    )

    # -------------------------------------------------------------------------
    # Extract coefficients in R
    # -------------------------------------------------------------------------

    A = sp.expand(resultant.coeff_monomial(R**2))
    B = sp.expand(resultant.coeff_monomial(R))
    C = sp.expand(resultant.coeff_monomial(1))

    print()
    print("RESULTANT AS:")
    print("  A(K) R^2 + B(K) R + C(K)")

    print()

    A_poly = summarize_poly("A(K)", A, K)
    B_poly = summarize_poly("B(K)", B, K)
    C_poly = summarize_poly("C(K)", C, K)

    degree_signature = (
        A_poly.degree(),
        B_poly.degree(),
        C_poly.degree()
    )

    degree_signatures.append(degree_signature)

    print()
    print("R-COEFFICIENT DEGREE SIGNATURE:")
    print(" ", degree_signature)

    # -------------------------------------------------------------------------
    # Factor coefficient polynomials
    # -------------------------------------------------------------------------

    print()
    print("FACTOR STRUCTURE")

    for label, expr in (
        ("A", A),
        ("B", B),
        ("C", C)
    ):

        factored = sp.factor(expr)

        term_count = len(
            sp.Poly(sp.expand(expr), K).terms()
        )

        print()
        print(f"  {label}:")

        if term_count <= 12:
            print("    ", factored)
        else:
            factor_list = sp.factor_list(expr)

            print(
                "    factor count =",
                len(factor_list[1])
            )

            for factor, multiplicity in factor_list[1]:
                fp = sp.Poly(factor, K)

                print(
                    f"      degree={fp.degree()} "
                    f"multiplicity={multiplicity} "
                    f"terms={len(fp.terms())}"
                )

    # -------------------------------------------------------------------------
    # Quadratic discriminant in R
    # -------------------------------------------------------------------------

    print()
    print("QUADRATIC DISCRIMINANT IN R")

    disc_R = sp.expand(
        B**2 - 4*A*C
    )

    disc_poly = sp.Poly(
        disc_R,
        K
    )

    discriminant_degree_signatures.append(
        (
            disc_poly.degree(),
            len(disc_poly.terms())
        )
    )

    print(
        "  degree =",
        disc_poly.degree()
    )

    print(
        "  terms  =",
        len(disc_poly.terms())
    )

    print(
        "  content =",
        disc_poly.content()
    )

    # -------------------------------------------------------------------------
    # Factor discriminant
    # -------------------------------------------------------------------------

    disc_factor = sp.factor(disc_R)

    if len(disc_poly.terms()) <= 20:
        print()
        print("  FACTORED DISCRIMINANT:")
        print("   ", disc_factor)
    else:
        factors = sp.factor_list(disc_R)

        print()
        print("  DISCRIMINANT FACTOR COUNT =",
              len(factors[1]))

        for factor, multiplicity in factors[1]:
            fp = sp.Poly(factor, K)

            print(
                "    degree=",
                fp.degree(),
                "multiplicity=",
                multiplicity,
                "terms=",
                len(fp.terms())
            )

    # -------------------------------------------------------------------------
    # TRUE R roots of resultant
    #
    # We know the true R exactly. Check whether the second R-root has a
    # simple interpretation.
    # -------------------------------------------------------------------------

    R_true = sp.cancel(
        NUM_R.subs(t, sp.Rational(N, X_true))
        /
        DEN_R.subs(t, sp.Rational(N, X_true))
    )

    quadratic_at_true = sp.cancel(
        A*R_true**2
        + B*R_true
        + C
    )

    print()
    print("TRUE POINT CHECK:")
    print("  F(K_TRUE,R_TRUE) =",
          sp.cancel(
              resultant.as_expr().subs({
                  K: K_true,
                  R: R_true
              })
          )
          )

    # -------------------------------------------------------------------------
    # Find the second R root exactly through Vieta.
    #
    # If R_true is one root:
    #
    #   R_2 = C / (A R_true)
    #
    # provided A and R_true are nonzero.
    # -------------------------------------------------------------------------

    if A != 0 and R_true != 0:

        R_second = sp.cancel(
            C / (A * R_true)
        )

        print()
        print("SECOND R-ROOT ANALYSIS")

        print("  second root term count =",
              len(
                  sp.Poly(
                      sp.together(R_second).as_numer_denom()[0],
                      K
                  ).terms()
              )
              if sp.together(R_second).as_numer_denom()[0].is_polynomial(K)
              else "non-polynomial")

        # Check simple relations against R_true.
        sum_roots = sp.cancel(
            R_true + R_second
        )

        product_roots = sp.cancel(
            R_true * R_second
        )

        print(
            "  sum-root expression size =",
            len(str(sum_roots))
        )

        print(
            "  product-root expression size =",
            len(str(product_roots))
        )

        # Check whether second root equals R at the second Kappa-root.
        other_t = sp.solve(
            kappa_poly.subs(K, K_true),
            t
        )

        if len(other_t) == 2:

            print()
            print("  KAPPA HAS", len(other_t),
                  "t-ROOTS AT K_TRUE")

            other_t_value = None

            for candidate in other_t:
                if sp.simplify(candidate - sp.Rational(N, X_true)) != 0:
                    other_t_value = candidate
                    break

            if other_t_value is not None:

                R_other_t = sp.cancel(
                    NUM_R.subs(t, other_t_value)
                    /
                    DEN_R.subs(t, other_t_value)
                )

                equality_check = sp.cancel(
                    R_other_t - R_second
                )

                print(
                    "  second resultant root == R(other t):",
                    equality_check == 0
                )

    # -------------------------------------------------------------------------
    # Sample leading coefficients
    # -------------------------------------------------------------------------

    print()
    print("LEADING COEFFICIENTS IN K")

    print("  LC(A) =", A_poly.LC())
    print("  LC(B) =", B_poly.LC())
    print("  LC(C) =", C_poly.LC())


# =============================================================================
# CROSS-INSTANCE SUMMARY
# =============================================================================

print()
print("=" * 110)
print("CROSS-INSTANCE SUMMARY")
print("=" * 110)

print()
print("R-coefficient degree signatures:")

for i, sig in enumerate(
    degree_signatures,
    start=1
):
    print(
        f"  instance {i:2d}: "
        f"A={sig[0]}, B={sig[1]}, C={sig[2]}"
    )

print()
print(
    "All R-coefficient degree signatures identical:",
    len(set(degree_signatures)) == 1
)

print()
print("Discriminant degree signatures:")

for i, sig in enumerate(
    discriminant_degree_signatures,
    start=1
):
    print(
        f"  instance {i:2d}: "
        f"degree={sig[0]}, terms={sig[1]}"
    )

print()
print(
    "All discriminant signatures identical:",
    len(set(discriminant_degree_signatures)) == 1
)

print()
print("=" * 110)
print("EXPERIMENT 382 FINAL STATUS")
print("=" * 110)

print()
print("The key object is now:")
print()
print("  F_N(K,R) = A_N(K) R^2 + B_N(K) R + C_N(K)")
print()
print("The next structural question is whether these three coefficient")
print("polynomials admit a much simpler common normalization in N and K.")
print()
print("In particular, inspect:")
print()
print("  1. A_N(K)")
print("  2. B_N(K)")
print("  3. C_N(K)")
print("  4. B_N(K)^2 - 4 A_N(K) C_N(K)")
print("  5. The second R-root")
print()
print("EXPERIMENT 382 FINISHED")
print("=" * 110)
