import sympy as sp


# =============================================================================
# EXPERIMENT 383 START
# =============================================================================

print("=" * 110)
print("EXPERIMENT 383 START")
print("=" * 110)
print()
print("CORRECTED RESULTANT COEFFICIENT STRUCTURE")
print()
print("Purpose:")
print("  1. Correctly extract A_N(K), B_N(K), C_N(K).")
print("  2. Determine their true K-degree signatures.")
print("  3. Determine gcd structure.")
print("  4. Analyze the R-quadratic discriminant.")
print("  5. Verify the second R-root exactly.")
print("  6. Compare the second R-root with the other Kappa t-root.")
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

h16 = (
    -9*t**8
    - 36*t**7
    - 84*t**6
    - 126*t**5
    - 126*t**4
    - 84*t**3
    - 36*t**2
    - 9*t
    - 1
)

h15 = (
    88*t**9
    + 396*t**8
    + 1164*t**7
    + 2226*t**6
    + 2898*t**5
    + 2604*t**4
    + 1596*t**3
    + 639*t**2
    + 151*t
    + 16
)

h14 = (
    -276*t**10
    - 1380*t**9
    - 5460*t**8
    - 13560*t**7
    - 23058*t**6
    - 27510*t**5
    - 23100*t**4
    - 13410*t**3
    - 5135*t**2
    - 1169*t
    - 120
)


# =============================================================================
# EXACT LAYER EQUATION
#
# R = h15(t)^2 / (h16(t)*h14(t))
#
# Therefore:
#
#     R*h16*h14 - h15^2 = 0
# =============================================================================

num_R = sp.expand(h15**2)
den_R = sp.expand(h16 * h14)

layer_poly = sp.expand(
    R * den_R - num_R
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
# KAPPA
# =============================================================================

def kappa_value(N, X):
    S = X - 1

    return (
        -S**2
        + (N + 1)*S
        - N**2
        + N
    )


def kappa_polynomial(N):
    return sp.expand(
        (K + N**2 + 2)*t**2
        - N*(N + 3)*t
        + N**2
    )


# =============================================================================
# SAFE INTEGER CONVERSION
# =============================================================================

def safe_int(value):
    """
    Convert a SymPy integer-like object to Python int.
    """
    value = sp.sympify(value)

    if value.is_Integer:
        return int(value)

    if value.is_Rational and value.q == 1:
        return int(value.p)

    return int(value)


def bit_length_of_integer(value):
    """
    Safe bit-length for SymPy/Python integer-like values.
    """
    value = abs(safe_int(value))
    return value.bit_length()


# =============================================================================
# POLYNOMIAL INFORMATION
# =============================================================================

def polynomial_info(expr, variable):
    poly = sp.Poly(
        sp.expand(expr),
        variable,
        domain=sp.QQ
    )

    return {
        "poly": poly,
        "degree": poly.degree(),
        "terms": len(poly.terms()),
        "content": poly.content(),
        "LC": poly.LC(),
        "TC": poly.TC(),
    }


# =============================================================================
# FACTOR INFORMATION
# =============================================================================

def factor_signature(expr, variable):
    poly = sp.Poly(
        sp.expand(expr),
        variable,
        domain=sp.QQ
    )

    coeff, factors = sp.factor_list(
        poly.as_expr()
    )

    return {
        "constant": coeff,
        "factor_count": len(factors),
        "factors": [
            {
                "degree": sp.Poly(factor, variable).degree(),
                "multiplicity": multiplicity,
                "terms": len(
                    sp.Poly(factor, variable).terms()
                ),
            }
            for factor, multiplicity in factors
        ],
    }


# =============================================================================
# CROSS-INSTANCE STORAGE
# =============================================================================

degree_signatures = []
gcd_degree_signatures = []
discriminant_signatures = []
true_root_failures = 0
second_root_failures = 0


# =============================================================================
# MAIN LOOP
# =============================================================================

for instance_no, (p, q) in enumerate(
    TEST_PAIRS,
    start=1
):

    print()
    print("=" * 110)
    print(f"INSTANCE {instance_no}")
    print("=" * 110)

    N = p * q
    X_true = p + q + 1
    K_true = sp.Integer(
        kappa_value(N, X_true)
    )

    t_true = sp.Rational(
        N,
        X_true
    )

    print()
    print("BASIC DATA")
    print("  p      =", p)
    print("  q      =", q)
    print("  N      =", N)
    print("  X_TRUE =", X_true)
    print("  K_TRUE =", K_true)
    print("  t_TRUE =", t_true)

    # =========================================================================
    # KAPPA EQUATION
    # =========================================================================

    kappa_poly = kappa_polynomial(N)

    kappa_true_check = sp.cancel(
        kappa_poly.subs(t, t_true)
    )

    print()
    print("KAPPA ROOT CHECK")
    print("  polynomial(t_TRUE) =", kappa_true_check)
    print(
        "  STATUS =",
        "PASS" if kappa_true_check == 0 else "FAIL"
    )

    if kappa_true_check != 0:
        true_root_failures += 1

    # =========================================================================
    # LAYER EQUATION
    # =========================================================================

    R_true = sp.cancel(
        num_R.subs(t, t_true)
        /
        den_R.subs(t, t_true)
    )

    layer_true_check = sp.cancel(
        layer_poly.subs({
            t: t_true,
            R: R_true
        })
    )

    # layer_poly contains t, so the above is the direct root test.
    print()
    print("LAYER ROOT CHECK")
    print("  polynomial(t_TRUE,R_TRUE) =", layer_true_check)
    print(
        "  STATUS =",
        "PASS" if layer_true_check == 0 else "FAIL"
    )

    if layer_true_check != 0:
        true_root_failures += 1

    # =========================================================================
    # RESULTANT
    # =========================================================================

    print()
    print("COMPUTING EXACT RESULTANT")

    resultant_expr = sp.resultant(
        kappa_poly,
        layer_poly,
        t
    )

    resultant_expr = sp.expand(
        resultant_expr
    )

    resultant_poly_KR = sp.Poly(
        resultant_expr,
        K,
        R,
        domain=sp.QQ
    )

    print()
    print("RESULTANT STRUCTURE")
    print("  total degree =", resultant_poly_KR.total_degree())
    print("  K degree     =", resultant_poly_KR.degree(K))
    print("  R degree     =", resultant_poly_KR.degree(R))
    print("  term count   =", len(resultant_poly_KR.terms()))

    # =========================================================================
    # CORRECT R-COEFFICIENT EXTRACTION
    #
    # We MUST first create a polynomial in R.
    #
    # Each coefficient is then a polynomial in K.
    # =========================================================================

    resultant_in_R = sp.Poly(
        resultant_expr,
        R,
        domain=sp.QQ[K]
    )

    A = sp.expand(
        resultant_in_R.coeff_monomial(
            R**2
        )
    )

    B = sp.expand(
        resultant_in_R.coeff_monomial(
            R
        )
    )

    C = sp.expand(
        resultant_in_R.coeff_monomial(
            1
        )
    )

    A_info = polynomial_info(A, K)
    B_info = polynomial_info(B, K)
    C_info = polynomial_info(C, K)

    A_poly = A_info["poly"]
    B_poly = B_info["poly"]
    C_poly = C_info["poly"]

    signature = (
        A_info["degree"],
        B_info["degree"],
        C_info["degree"]
    )

    degree_signatures.append(
        signature
    )

    print()
    print("CORRECT R-COEFFICIENT STRUCTURE")
    print("  Resultant = A(K) R^2 + B(K) R + C(K)")
    print()
    print("  A(K):")
    print("    degree =", A_info["degree"])
    print("    terms  =", A_info["terms"])
    print("    LC     =", A_info["LC"])
    print("    content=", A_info["content"])
    print()
    print("  B(K):")
    print("    degree =", B_info["degree"])
    print("    terms  =", B_info["terms"])
    print("    LC     =", B_info["LC"])
    print("    content=", B_info["content"])
    print()
    print("  C(K):")
    print("    degree =", C_info["degree"])
    print("    terms  =", C_info["terms"])
    print("    LC     =", C_info["LC"])
    print("    content=", C_info["content"])

    # =========================================================================
    # COMMON GCD OF A, B, C
    # =========================================================================

    gcd_AB = sp.gcd(
        A_poly,
        B_poly
    )

    gcd_ABC = sp.gcd(
        gcd_AB,
        C_poly
    )

    gcd_ABC = sp.Poly(
        gcd_ABC.as_expr(),
        K,
        domain=sp.QQ
    )

    gcd_degree_signatures.append(
        gcd_ABC.degree()
    )

    print()
    print("COMMON POLYNOMIAL GCD")
    print(
        "  degree =",
        gcd_ABC.degree()
    )
    print(
        "  content =",
        gcd_ABC.content()
    )

    # =========================================================================
    # PRIMITIVE PART OF THE WHOLE RESULTANT WITH RESPECT TO K
    #
    # Important:
    # Poly.primitive() returns:
    #
    #   (content, primitive_poly)
    #
    # not the reverse.
    # =========================================================================

    resultant_as_K_poly = sp.Poly(
        resultant_expr,
        K,
        domain=sp.QQ[R]
    )

    primitive_content, primitive_poly = (
        resultant_as_K_poly.primitive()
    )

    primitive_content = sp.sympify(
        primitive_content
    )

    primitive_poly = sp.Poly(
        primitive_poly,
        K,
        domain=sp.QQ[R]
    )

    print()
    print("WHOLE RESULTANT K-CONTENT")
    print(
        "  content bit-length =",
        bit_length_of_integer(
            primitive_content
        )
    )

    # =========================================================================
    # DISCRIMINANT WITH RESPECT TO R
    # =========================================================================

    disc_R = sp.expand(
        B**2 - 4*A*C
    )

    disc_info = polynomial_info(
        disc_R,
        K
    )

    disc_factor_info = factor_signature(
        disc_R,
        K
    )

    discriminant_signatures.append(
        (
            disc_info["degree"],
            disc_info["terms"],
        )
    )

    print()
    print("R-QUADRATIC DISCRIMINANT")
    print(
        "  degree =",
        disc_info["degree"]
    )
    print(
        "  terms  =",
        disc_info["terms"]
    )
    print(
        "  content bit-length =",
        bit_length_of_integer(
            disc_info["content"]
        )
    )
    print(
        "  nonconstant factor count =",
        disc_factor_info["factor_count"]
    )

    print(
        "  factor degree pattern =",
        [
            (
                f["degree"],
                f["multiplicity"],
                f["terms"]
            )
            for f in disc_factor_info["factors"]
        ]
    )

    # =========================================================================
    # TRUE RESULTANT COMPATIBILITY
    # =========================================================================

    true_resultant_check = sp.cancel(
        resultant_expr.subs({
            K: K_true,
            R: R_true
        })
    )

    print()
    print("RESULTANT TRUE-POINT CHECK")
    print(
        "  resultant(K_TRUE,R_TRUE) =",
        true_resultant_check
    )
    print(
        "  STATUS =",
        "PASS"
        if true_resultant_check == 0
        else "FAIL"
    )

    if true_resultant_check != 0:
        true_root_failures += 1

    # =========================================================================
    # SECOND R-ROOT
    # =========================================================================

    A_true = sp.cancel(
        A.subs(K, K_true)
    )

    B_true = sp.cancel(
        B.subs(K, K_true)
    )

    C_true = sp.cancel(
        C.subs(K, K_true)
    )

    print()
    print("SECOND R-ROOT ANALYSIS")

    if A_true == 0:

        print("  A(K_TRUE) = 0")
        print("  Quadratic degenerates: cannot use Vieta.")

    else:

        second_R_vieta = sp.cancel(
            C_true
            /
            (A_true * R_true)
        )

        r_quadratic_at_true = sp.Poly(
            sp.expand(
                resultant_expr.subs(
                    K,
                    K_true
                )
            ),
            R,
            domain=sp.QQ
        )

        direct_R_roots = sp.solve(
            r_quadratic_at_true.as_expr(),
            R
        )

        print(
            "  number of R roots =",
            len(direct_R_roots)
        )

        direct_second_R = None

        for root in direct_R_roots:
            if sp.cancel(root - R_true) != 0:
                direct_second_R = sp.cancel(
                    root
                )
                break

        if direct_second_R is None:

            print(
                "  second distinct R-root = NONE"
            )

        else:

            vieta_match = (
                sp.cancel(
                    direct_second_R
                    - second_R_vieta
                )
                == 0
            )

            print(
                "  Vieta second root "
                "matches direct root:",
                vieta_match
            )

            # ================================================================
            # KAPPA'S OTHER t ROOT
            # ================================================================

            kappa_true_poly = sp.Poly(
                kappa_poly.subs(
                    K,
                    K_true
                ),
                t,
                domain=sp.QQ
            )

            t_roots = sp.solve(
                kappa_true_poly.as_expr(),
                t
            )

            print(
                "  number of Kappa t-roots =",
                len(t_roots)
            )

            other_t = None

            for root in t_roots:
                if sp.cancel(
                    root - t_true
                ) != 0:
                    other_t = sp.cancel(
                        root
                    )
                    break

            if other_t is None:

                print(
                    "  other Kappa t-root = NONE"
                )

            else:

                R_other_t = sp.cancel(
                    num_R.subs(
                        t,
                        other_t
                    )
                    /
                    den_R.subs(
                        t,
                        other_t
                    )
                )

                branch_match = (
                    sp.cancel(
                        R_other_t
                        - direct_second_R
                    )
                    == 0
                )

                print(
                    "  second R-root = "
                    "R(other Kappa t-root):",
                    branch_match
                )

                if not branch_match:
                    second_root_failures += 1

                print(
                    "  other t numerator digits =",
                    len(
                        str(
                            sp.numer(
                                other_t
                            )
                        )
                    )
                )

                print(
                    "  other t denominator digits =",
                    len(
                        str(
                            sp.denom(
                                other_t
                            )
                        )
                    )
                )

    # =========================================================================
    # TRUE-K NORMALIZED STRUCTURE
    # =========================================================================

    print()
    print("TRUE-K NORMALIZED RATIOS")

    if A_true != 0 and C_true != 0:

        ratio_BA = sp.cancel(
            B_true / A_true
        )

        ratio_CA = sp.cancel(
            C_true / A_true
        )

        ratio_B2AC = sp.cancel(
            B_true**2
            /
            (A_true * C_true)
        )

        print(
            "  B/A:",
            "num_digits=",
            len(
                str(
                    abs(
                        safe_int(
                            sp.numer(
                                ratio_BA
                            )
                        )
                    )
                )
            ),
            "den_digits=",
            len(
                str(
                    abs(
                        safe_int(
                            sp.denom(
                                ratio_BA
                            )
                        )
                    )
                )
            )
        )

        print(
            "  C/A:",
            "num_digits=",
            len(
                str(
                    abs(
                        safe_int(
                            sp.numer(
                                ratio_CA
                            )
                        )
                    )
                )
            ),
            "den_digits=",
            len(
                str(
                    abs(
                        safe_int(
                            sp.denom(
                                ratio_CA
                            )
                        )
                    )
                )
            )
        )

        print(
            "  B^2/(A*C):",
            "num_digits=",
            len(
                str(
                    abs(
                        safe_int(
                            sp.numer(
                                ratio_B2AC
                            )
                        )
                    )
                )
            ),
            "den_digits=",
            len(
                str(
                    abs(
                        safe_int(
                            sp.denom(
                                ratio_B2AC
                            )
                        )
                    )
                )
            )
        )


# =============================================================================
# CROSS-INSTANCE SUMMARY
# =============================================================================

print()
print("=" * 110)
print("CROSS-INSTANCE SUMMARY")
print("=" * 110)

print()
print("R-COEFFICIENT DEGREE SIGNATURES")

for i, signature in enumerate(
    degree_signatures,
    start=1
):
    print(
        f"  instance {i:2d}: "
        f"A={signature[0]} "
        f"B={signature[1]} "
        f"C={signature[2]}"
    )

print()
print(
    "All A/B/C degree signatures identical:",
    len(set(degree_signatures)) == 1
)

print()
print("COMMON GCD DEGREE SIGNATURES")

for i, degree in enumerate(
    gcd_degree_signatures,
    start=1
):
    print(
        f"  instance {i:2d}: "
        f"gcd degree={degree}"
    )

print()
print(
    "All gcd degrees identical:",
    len(set(gcd_degree_signatures)) == 1
)

print()
print("DISCRIMINANT DEGREE SIGNATURES")

for i, signature in enumerate(
    discriminant_signatures,
    start=1
):
    print(
        f"  instance {i:2d}: "
        f"degree={signature[0]} "
        f"terms={signature[1]}"
    )

print()
print(
    "All discriminant degree/term signatures identical:",
    len(set(discriminant_signatures)) == 1
)

print()
print("GLOBAL CHECKS")
print(
    "  true-root failures =",
    true_root_failures
)
print(
    "  second-root branch mismatches =",
    second_root_failures
)

print()
print("=" * 110)
print("EXPERIMENT 383 FINAL STATUS")
print("=" * 110)
print()
print("Experiment 382 accidentally extracted only the K^0 monomial")
print("coefficients of R^2, R, and 1.")
print()
print("Experiment 383 correctly treats:")
print()
print("    F_N(K,R) = A_N(K) R^2 + B_N(K) R + C_N(K)")
print()
print("where A_N, B_N and C_N are genuine polynomials in K.")
print()
print("The decisive questions are now:")
print("  1. What are deg_K(A_N), deg_K(B_N), deg_K(C_N)?")
print("  2. Do A_N, B_N, C_N have a common polynomial factor?")
print("  3. Does the discriminant have a stable factor structure?")
print("  4. Does the second R-root correspond to the other Kappa t-root?")
print("  5. Does the normalized coefficient structure repeat across N?")
print()
print("EXPERIMENT 383 FINISHED")
print("=" * 110)