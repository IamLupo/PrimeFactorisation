# =============================================================================
# EXPERIMENT 385 START
# =============================================================================
#
# UNIVERSAL NORMALIZED RESULTANT / EXACT TWO-ROOT BRANCH
#
# Corrections from Experiment 384:
#
#   1. Do NOT use sympy.solve() to identify the Kappa roots.
#   2. Use the exact quadratic/Vieta identity:
#
#        (K + N^2 + 2)t^2 - N(N+3)t + N^2 = 0
#
#      Product of roots:
#
#        t1*t2 = N^2 / (K + N^2 + 2)
#
#      Therefore:
#
#        t_other =
#          N^2 / ((K_TRUE + N^2 + 2) * t_TRUE)
#
#   3. Verify the Kappa root after substituting BOTH t_TRUE and K_TRUE.
#
# Main goals:
#   - verify exact Kappa root
#   - extract A_N(K), B_N(K), C_N(K)
#   - verify exact resultant compatibility
#   - determine normalized coefficient structure
#   - determine discriminant factor structure
#   - verify Vieta second R-root
#   - verify second R-root == R(t_other)
#   - compare normalized structures across all N
#
# Exact arithmetic only.
# Huge polynomial expressions are never printed.
#
# =============================================================================

import sympy as sp


# -----------------------------------------------------------------------------
# SYMBOLS
# -----------------------------------------------------------------------------

K, R, t = sp.symbols("K R t")


# -----------------------------------------------------------------------------
# INSTANCES
# -----------------------------------------------------------------------------

INSTANCES = [
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
]


# -----------------------------------------------------------------------------
# KAPPA
# -----------------------------------------------------------------------------

def kappa_value(p, q):
    """
    Exact Kappa invariant:

        K = 1 - (p^2-p+1)(q^2-q+1)
    """
    return sp.Integer(1) - (
        p**2 - p + 1
    ) * (
        q**2 - q + 1
    )


def kappa_polynomial(N_value):
    """
    Exact quadratic:

        (K + N^2 + 2)t^2
        - N(N+3)t
        + N^2
        = 0
    """
    return sp.expand(
        (K + N_value**2 + 2) * t**2
        - N_value * (N_value + 3) * t
        + N_value**2
    )


# -----------------------------------------------------------------------------
# LAYER R(t)
# -----------------------------------------------------------------------------

def layer_function():

    h16 = (
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

    h15 = (
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

    h14 = (
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

    numerator = sp.expand(
        h15**2 - h16*h14
    )

    denominator = sp.expand(
        h16**2
    )

    return sp.cancel(
        numerator / denominator
    )


R_t = layer_function()


# -----------------------------------------------------------------------------
# RESULTANT
# -----------------------------------------------------------------------------

def make_resultant(N_value):

    kappa_eq = kappa_polynomial(N_value)

    numerator, denominator = sp.fraction(R_t)

    layer_eq = sp.expand(
        denominator * R - numerator
    )

    resultant = sp.resultant(
        kappa_eq,
        layer_eq,
        t,
    )

    return sp.Poly(
        sp.expand(resultant),
        K,
        R,
        domain=sp.QQ,
    )


# -----------------------------------------------------------------------------
# EXTRACT A(K), B(K), C(K)
# -----------------------------------------------------------------------------

def extract_ABC(resultant):

    poly_R = sp.Poly(
        resultant.as_expr(),
        R,
        domain=sp.QQ.frac_field(K),
    )

    A_expr = poly_R.coeff_monomial(R**2)
    B_expr = poly_R.coeff_monomial(R)
    C_expr = poly_R.coeff_monomial(1)

    A = sp.Poly(
        A_expr,
        K,
        domain=sp.QQ,
    )

    B = sp.Poly(
        B_expr,
        K,
        domain=sp.QQ,
    )

    C = sp.Poly(
        C_expr,
        K,
        domain=sp.QQ,
    )

    return A, B, C


# -----------------------------------------------------------------------------
# POLYNOMIAL HELPERS
# -----------------------------------------------------------------------------

def monic_poly(poly):
    poly = sp.Poly(poly, K, domain=sp.QQ)

    return sp.Poly(
        sp.cancel(poly.as_expr() / poly.LC()),
        K,
        domain=sp.QQ,
    )


def primitive_content(poly):
    """
    Safe extraction of integer polynomial content.
    """
    poly = sp.Poly(poly, K, domain=sp.QQ)

    integer_poly = poly.clear_denoms()[1]

    integer_poly = sp.Poly(
        integer_poly,
        K,
        domain=sp.ZZ,
    )

    content, primitive = sp.polys.polytools.primitive(
        integer_poly
    )

    return content, primitive


# -----------------------------------------------------------------------------
# DISCRIMINANT
# -----------------------------------------------------------------------------

def discriminant_R(A, B, C):

    return sp.Poly(
        sp.expand(
            B.as_expr()**2
            - 4*A.as_expr()*C.as_expr()
        ),
        K,
        domain=sp.QQ,
    )


# -----------------------------------------------------------------------------
# EXACT KAPPA ROOT CHECK
# -----------------------------------------------------------------------------

def check_kappa_root(N_value, K_true, t_true):

    expr = kappa_polynomial(N_value)

    value = sp.cancel(
        expr.subs({
            K: K_true,
            t: t_true,
        })
    )

    return value == 0, value


# -----------------------------------------------------------------------------
# EXACT OTHER KAPPA ROOT
# -----------------------------------------------------------------------------
#
# For:
#
#     a t^2 + b t + c = 0
#
# product of roots:
#
#     t1*t2 = c/a
#
# Here:
#
#     a = K + N^2 + 2
#     c = N^2
#
# Therefore:
#
#     t_other = N^2 / (a*t_true)
#
# -----------------------------------------------------------------------------

def exact_other_kappa_root(N_value, K_true, t_true):

    a = sp.Integer(K_true) + N_value**2 + 2

    t_other = sp.cancel(
        sp.Rational(N_value**2, 1)
        / (
            a * t_true
        )
    )

    return t_other


# -----------------------------------------------------------------------------
# VERIFY OTHER KAPPA ROOT
# -----------------------------------------------------------------------------

def verify_other_kappa_root(N_value, K_true, t_other):

    value = sp.cancel(
        kappa_polynomial(N_value).subs({
            K: K_true,
            t: t_other,
        })
    )

    return value == 0, value


# -----------------------------------------------------------------------------
# LAYER EVALUATION
# -----------------------------------------------------------------------------

def layer_at(value):

    return sp.cancel(
        R_t.subs(t, value)
    )


# -----------------------------------------------------------------------------
# SECOND R ROOT BY VIETA
# -----------------------------------------------------------------------------

def second_R_root(A, B, R_true):

    return sp.cancel(
        -B.as_expr() / A.as_expr()
        - R_true
    )


# -----------------------------------------------------------------------------
# NORMALIZED SHAPE TEST
# -----------------------------------------------------------------------------

def same_normalized_shape(poly1, poly2):

    return sp.expand(
        monic_poly(poly1).as_expr()
        - monic_poly(poly2).as_expr()
    ) == 0


# -----------------------------------------------------------------------------
# FACTOR DEGREE STRUCTURE
# -----------------------------------------------------------------------------

def factor_degree_signature(poly):

    monic = monic_poly(poly)

    _, factors = sp.factor_list(
        monic.as_expr(),
        K,
    )

    result = []

    for factor, multiplicity in factors:

        degree = sp.degree(
            factor,
            K,
        )

        if degree > 0:
            result.append(
                (
                    int(degree),
                    int(multiplicity),
                )
            )

    return result


# -----------------------------------------------------------------------------
# MAIN
# -----------------------------------------------------------------------------

def main():

    print("=" * 110)
    print("EXPERIMENT 385 START")
    print("=" * 110)
    print()
    print("UNIVERSAL NORMALIZED RESULTANT / EXACT TWO-ROOT BRANCH")
    print()
    print("Corrected Kappa root mechanism:")
    print()
    print("  (K + N^2 + 2)t^2 - N(N+3)t + N^2 = 0")
    print()
    print("  t_TRUE * t_OTHER = N^2 / (K + N^2 + 2)")
    print()
    print("No numerical approximations are used.")
    print()

    records = []

    reference_A = None
    reference_B = None
    reference_C = None
    reference_D = None

    # =====================================================================
    # INSTANCES
    # =====================================================================

    for index, (p, q) in enumerate(
        INSTANCES,
        start=1,
    ):

        print("=" * 110)
        print(f"INSTANCE {index}")
        print("=" * 110)
        print()

        N = p * q
        S = p + q
        X_true = S + 1

        t_true = sp.Rational(
            N,
            X_true,
        )

        K_true = kappa_value(
            p,
            q,
        )

        R_true = layer_at(
            t_true
        )

        print("BASIC DATA")
        print(f"  p       = {p}")
        print(f"  q       = {q}")
        print(f"  N       = {N}")
        print(f"  X_TRUE  = {X_true}")
        print(f"  t_TRUE  = {t_true}")
        print()

        # -----------------------------------------------------------------
        # KAPPA ROOT
        # -----------------------------------------------------------------

        kappa_pass, kappa_value_at_root = (
            check_kappa_root(
                N,
                K_true,
                t_true,
            )
        )

        print("KAPPA ROOT CHECK")
        print(
            f"  polynomial(t_TRUE,K_TRUE) = "
            f"{kappa_value_at_root}"
        )
        print(
            f"  STATUS = "
            f"{'PASS' if kappa_pass else 'FAIL'}"
        )
        print()

        # -----------------------------------------------------------------
        # RESULTANT
        # -----------------------------------------------------------------

        print("COMPUTING EXACT RESULTANT")

        resultant = make_resultant(
            N
        )

        print()

        print("RESULTANT STRUCTURE")
        print(
            f"  total degree = "
            f"{resultant.total_degree()}"
        )
        print(
            f"  K degree     = "
            f"{resultant.degree(K)}"
        )
        print(
            f"  R degree     = "
            f"{resultant.degree(R)}"
        )
        print(
            f"  term count   = "
            f"{len(resultant.terms())}"
        )
        print()

        # -----------------------------------------------------------------
        # A/B/C
        # -----------------------------------------------------------------

        A, B, C = extract_ABC(
            resultant
        )

        A_content, _ = primitive_content(A)
        B_content, _ = primitive_content(B)
        C_content, _ = primitive_content(C)

        print("A/B/C STRUCTURE")

        print(
            f"  A(K): degree={A.degree()} "
            f"terms={len(A.terms())} "
            f"LC={A.LC()} "
            f"content={A_content}"
        )

        print(
            f"  B(K): degree={B.degree()} "
            f"terms={len(B.terms())} "
            f"LC={B.LC()} "
            f"content={B_content}"
        )

        print(
            f"  C(K): degree={C.degree()} "
            f"terms={len(C.terms())} "
            f"LC={C.LC()} "
            f"content={C_content}"
        )

        print()

        # -----------------------------------------------------------------
        # COMMON GCD
        # -----------------------------------------------------------------

        gcd_abc = sp.gcd(
            sp.gcd(A, B),
            C,
        )

        gcd_content, gcd_primitive = (
            primitive_content(
                gcd_abc
            )
        )

        print("COMMON POLYNOMIAL GCD")
        print(
            f"  degree  = "
            f"{gcd_abc.degree()}"
        )
        print(
            f"  content = "
            f"{gcd_content}"
        )
        print()

        # -----------------------------------------------------------------
        # DISCRIMINANT
        # -----------------------------------------------------------------

        D = discriminant_R(
            A,
            B,
            C,
        )

        D_content, D_primitive = (
            primitive_content(
                D
            )
        )

        print("R-QUADRATIC DISCRIMINANT")
        print(
            f"  degree = {D.degree()}"
        )
        print(
            f"  terms  = {len(D.terms())}"
        )
        print(
            f"  content bit-length = "
            f"{abs(int(D_content)).bit_length()}"
        )

        D_signature = factor_degree_signature(
            D
        )

        print(
            f"  factor degree signature = "
            f"{D_signature}"
        )
        print()

        # -----------------------------------------------------------------
        # RESULTANT TRUE POINT
        # -----------------------------------------------------------------

        resultant_at_true = sp.cancel(
            resultant.as_expr().subs({
                K: K_true,
                R: R_true,
            })
        )

        resultant_pass = (
            resultant_at_true == 0
        )

        print("RESULTANT TRUE-POINT CHECK")
        print(
            f"  resultant(K_TRUE,R_TRUE) = "
            f"{resultant_at_true}"
        )
        print(
            f"  STATUS = "
            f"{'PASS' if resultant_pass else 'FAIL'}"
        )
        print()

        # -----------------------------------------------------------------
        # OTHER KAPPA ROOT
        # -----------------------------------------------------------------

        t_other = exact_other_kappa_root(
            N,
            K_true,
            t_true,
        )

        other_root_pass, other_root_value = (
            verify_other_kappa_root(
                N,
                K_true,
                t_other,
            )
        )

        print("SECOND KAPPA ROOT")
        print(
            f"  exact verification = "
            f"{'PASS' if other_root_pass else 'FAIL'}"
        )
        print(
            f"  polynomial(t_OTHER) = "
            f"{other_root_value}"
        )
        print(
            f"  numerator digits = "
            f"{len(str(abs(sp.numer(t_other))))}"
        )
        print(
            f"  denominator digits = "
            f"{len(str(abs(sp.denom(t_other))))}"
        )
        print()

        # -----------------------------------------------------------------
        # SECOND R ROOT
        # -----------------------------------------------------------------

        R_other_direct = layer_at(
            t_other
        )

        R_second = second_R_root(
            A,
            B,
            R_true,
        )

        second_R_matches_direct = (
            sp.cancel(
                R_second
                - R_other_direct
            )
            == 0
        )

        second_R_vieta_check = sp.cancel(
            A.as_expr() * R_second**2
            + B.as_expr() * R_second
            + C.as_expr()
        ) == 0

        print("SECOND R-ROOT ANALYSIS")
        print(
            f"  Vieta quadratic check = "
            f"{'PASS' if second_R_vieta_check else 'FAIL'}"
        )
        print(
            f"  second R-root == R(t_OTHER) = "
            f"{'PASS' if second_R_matches_direct else 'FAIL'}"
        )
        print()

        # -----------------------------------------------------------------
        # NORMALIZATION
        # -----------------------------------------------------------------

        A_monic = monic_poly(A)
        B_monic = monic_poly(B)
        C_monic = monic_poly(C)
        D_monic = monic_poly(D)

        if reference_A is None:

            reference_A = A_monic
            reference_B = B_monic
            reference_C = C_monic
            reference_D = D_monic

            print("NORMALIZED SHAPES")
            print("  reference instance initialized")

        else:

            print("NORMALIZED SHAPES")

            print(
                "  A shape == instance 1:",
                same_normalized_shape(
                    A,
                    reference_A,
                )
            )

            print(
                "  B shape == instance 1:",
                same_normalized_shape(
                    B,
                    reference_B,
                )
            )

            print(
                "  C shape == instance 1:",
                same_normalized_shape(
                    C,
                    reference_C,
                )
            )

            print(
                "  D shape == instance 1:",
                same_normalized_shape(
                    D,
                    reference_D,
                )
            )

        print()

        # -----------------------------------------------------------------
        # SAVE
        # -----------------------------------------------------------------

        records.append({
            "N": N,
            "K": K_true,
            "R": R_true,
            "A": A,
            "B": B,
            "C": C,
            "D": D,
            "kappa_pass": kappa_pass,
            "resultant_pass": resultant_pass,
            "other_root_pass": other_root_pass,
            "second_R_vieta_pass": second_R_vieta_check,
            "second_R_direct_pass": second_R_matches_direct,
            "D_signature": D_signature,
        })

    # =========================================================================
    # CROSS INSTANCE
    # =========================================================================

    print("=" * 110)
    print("CROSS-INSTANCE SUMMARY")
    print("=" * 110)
    print()

    print("A/B/C DEGREE SIGNATURES")

    for i, rec in enumerate(
        records,
        start=1,
    ):

        print(
            f"  instance {i:2d}: "
            f"A={rec['A'].degree()} "
            f"B={rec['B'].degree()} "
            f"C={rec['C'].degree()}"
        )

    print()

    print("DISCRIMINANT SIGNATURES")

    for i, rec in enumerate(
        records,
        start=1,
    ):

        print(
            f"  instance {i:2d}: "
            f"D={rec['D'].degree()} "
            f"terms={len(rec['D'].terms())} "
            f"factor={rec['D_signature']}"
        )

    print()

    all_abc_degree_equal = all(
        (
            rec["A"].degree(),
            rec["B"].degree(),
            rec["C"].degree(),
        )
        ==
        (
            records[0]["A"].degree(),
            records[0]["B"].degree(),
            records[0]["C"].degree(),
        )
        for rec in records
    )

    all_D_signature_equal = all(
        rec["D_signature"]
        == records[0]["D_signature"]
        for rec in records
    )

    all_kappa_pass = all(
        rec["kappa_pass"]
        for rec in records
    )

    all_resultant_pass = all(
        rec["resultant_pass"]
        for rec in records
    )

    all_other_root_pass = all(
        rec["other_root_pass"]
        for rec in records
    )

    all_second_R_vieta_pass = all(
        rec["second_R_vieta_pass"]
        for rec in records
    )

    all_second_R_direct_pass = all(
        rec["second_R_direct_pass"]
        for rec in records
    )

    print("GLOBAL CHECKS")
    print(
        f"  all Kappa roots pass        = "
        f"{all_kappa_pass}"
    )
    print(
        f"  all resultant checks pass   = "
        f"{all_resultant_pass}"
    )
    print(
        f"  all other-t roots pass      = "
        f"{all_other_root_pass}"
    )
    print(
        f"  all second-R Vieta pass     = "
        f"{all_second_R_vieta_pass}"
    )
    print(
        f"  all second-R direct pass    = "
        f"{all_second_R_direct_pass}"
    )
    print(
        f"  all A/B/C degrees identical = "
        f"{all_abc_degree_equal}"
    )
    print(
        f"  all D signatures identical  = "
        f"{all_D_signature_equal}"
    )

    print()
    print("=" * 110)
    print("EXPERIMENT 385 FINAL STATUS")
    print("=" * 110)
    print()

    print("The central exact correspondence being tested is:")
    print()
    print("  Kappa quadratic:")
    print("      (K + N^2 + 2)t^2 - N(N+3)t + N^2 = 0")
    print()
    print("  Vieta:")
    print("      t_TRUE * t_OTHER = N^2/(K + N^2 + 2)")
    print()
    print("  Resultant:")
    print("      F_N(K,R) = A_N(K)R^2 + B_N(K)R + C_N(K)")
    print()
    print("  Layer branches:")
    print("      R_TRUE   = R(t_TRUE)")
    print("      R_OTHER  = R(t_OTHER)")
    print()
    print("If every branch passes, then the apparent two-root structure")
    print("is exactly the image of the two Kappa t-roots under R(t).")
    print()
    print("The next structural target is then the normalized map")
    print()
    print("      K  -> {R(t_1), R(t_2)}")
    print()
    print("and whether its symmetric quantities")
    print()
    print("      R1 + R2")
    print("      R1 * R2")
    print()
    print("have substantially simpler closed forms than A_N, B_N, C_N.")
    print()
    print("=" * 110)
    print("EXPERIMENT 385 FINISHED")
    print("=" * 110)


if __name__ == "__main__":
    main()