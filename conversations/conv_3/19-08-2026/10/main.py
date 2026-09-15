# =============================================================================
# EXPERIMENT 387 START
# =============================================================================
#
# EXACT RESULTANT SPECIALIZATION AT K_TRUE
#
# Goal:
#
#   Experiment 385 showed:
#
#       Kappa roots                  PASS
#       resultant(K_TRUE,R_TRUE)     PASS
#       second-R Vieta/direct        FAIL
#
#   Experiment 386 introduced symbolic K-rational functions and therefore
#   failed during polynomial coercion.
#
#   This experiment avoids symbolic division entirely.
#
#   For each instance:
#
#       1. Build the exact resultant F_N(K,R).
#       2. Substitute K = K_TRUE immediately.
#       3. Obtain an ordinary quadratic:
#
#              a R^2 + b R + c
#
#       4. Compute its exact two R-roots.
#       5. Compute R(t_TRUE) directly.
#       6. Compute R(t_OTHER) directly.
#       7. Compare the two sets exactly.
#
#   This determines whether the resultant itself really has the two
#   expected layer branches.
#
# No floating point.
# No symbolic division by A(K).
# No candidate-prime search.
#
# =============================================================================

import sympy as sp


# =============================================================================
# SYMBOLS
# =============================================================================

K, R, t = sp.symbols("K R t")


# =============================================================================
# TEST INSTANCES
# =============================================================================

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


# =============================================================================
# KAPPA
# =============================================================================

def kappa_value(p, q):

    return sp.Integer(1) - (
        p**2 - p + 1
    ) * (
        q**2 - q + 1
    )


def kappa_polynomial(N):

    return sp.expand(
        (K + N**2 + 2) * t**2
        - N * (N + 3) * t
        + N**2
    )


# =============================================================================
# LAYER R(t)
# =============================================================================

def build_R():

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


R_t = build_R()


# =============================================================================
# RESULTANT
# =============================================================================

def build_resultant(N):

    kappa_eq = kappa_polynomial(N)

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


# =============================================================================
# EXACT OTHER KAPPA ROOT
# =============================================================================

def other_kappa_root(
    N,
    K_true,
    t_true,
):

    a = sp.Integer(K_true) + N**2 + 2

    return sp.cancel(
        sp.Rational(N**2, 1)
        / (a * t_true)
    )


# =============================================================================
# EXACT LAYER VALUE
# =============================================================================

def layer_at(t_value):

    return sp.cancel(
        R_t.subs(
            t,
            t_value,
        )
    )


# =============================================================================
# NORMALIZE A RATIONAL NUMBER
# =============================================================================

def rational_signature(x):

    x = sp.cancel(x)

    numerator = sp.numer(x)
    denominator = sp.denom(x)

    return (
        numerator,
        denominator,
    )


# =============================================================================
# EXACT ROOT SET COMPARISON
# =============================================================================

def root_set_match(
    root_a,
    root_b,
    target_a,
    target_b,
):

    direct_same = (
        root_a == target_a
        and root_b == target_b
    )

    swapped_same = (
        root_a == target_b
        and root_b == target_a
    )

    return direct_same or swapped_same


# =============================================================================
# MAIN
# =============================================================================

def main():

    print("=" * 110)
    print("EXPERIMENT 387 START")
    print("=" * 110)
    print()
    print("EXACT RESULTANT SPECIALIZATION AT K_TRUE")
    print()
    print("No symbolic division by A(K) is performed.")
    print("The resultant is specialized at exact K_TRUE first.")
    print()

    results = []

    for index, (p, q) in enumerate(
        INSTANCES,
        start=1,
    ):

        print("=" * 110)
        print(f"INSTANCE {index}")
        print("=" * 110)
        print()

        # ---------------------------------------------------------------------
        # BASIC VALUES
        # ---------------------------------------------------------------------

        N = p * q

        X_true = p + q + 1

        K_true = kappa_value(
            p,
            q,
        )

        t_true = sp.Rational(
            N,
            X_true,
        )

        t_other = other_kappa_root(
            N,
            K_true,
            t_true,
        )

        R_true = layer_at(
            t_true
        )

        R_other = layer_at(
            t_other
        )

        # ---------------------------------------------------------------------
        # ROOT CHECK
        # ---------------------------------------------------------------------

        kappa_eq = kappa_polynomial(
            N
        )

        kappa_true_check = sp.cancel(
            kappa_eq.subs(
                {
                    K: K_true,
                    t: t_true,
                }
            )
        )

        kappa_other_check = sp.cancel(
            kappa_eq.subs(
                {
                    K: K_true,
                    t: t_other,
                }
            )
        )

        # ---------------------------------------------------------------------
        # BUILD RESULTANT
        # ---------------------------------------------------------------------

        resultant = build_resultant(
            N
        )

        # ---------------------------------------------------------------------
        # SPECIALIZE K FIRST
        # ---------------------------------------------------------------------

        specialized_expr = sp.expand(
            resultant.as_expr().subs(
                K,
                K_true,
            )
        )

        specialized_poly = sp.Poly(
            specialized_expr,
            R,
            domain=sp.QQ,
        )

        degree_R = specialized_poly.degree()

        # ---------------------------------------------------------------------
        # COEFFICIENTS
        # ---------------------------------------------------------------------

        a = sp.cancel(
            specialized_poly.coeff_monomial(
                R**2
            )
        )

        b = sp.cancel(
            specialized_poly.coeff_monomial(
                R
            )
        )

        c = sp.cancel(
            specialized_poly.coeff_monomial(
                1
            )
        )

        # ---------------------------------------------------------------------
        # DIRECT ROOT CHECK
        # ---------------------------------------------------------------------

        resultant_at_R_true = sp.cancel(
            specialized_poly.as_expr().subs(
                R,
                R_true,
            )
        )

        resultant_at_R_other = sp.cancel(
            specialized_poly.as_expr().subs(
                R,
                R_other,
            )
        )

        # ---------------------------------------------------------------------
        # DISCRIMINANT
        # ---------------------------------------------------------------------

        discriminant = sp.cancel(
            b**2
            - 4*a*c
        )

        discriminant_square = sp.sqrt(
            discriminant
        )

        # ---------------------------------------------------------------------
        # QUADRATIC ROOTS
        # ---------------------------------------------------------------------

        resultant_roots = []

        if degree_R == 2:

            root_plus = sp.cancel(
                (-b + discriminant_square)
                / (2*a)
            )

            root_minus = sp.cancel(
                (-b - discriminant_square)
                / (2*a)
            )

            resultant_roots = [
                root_plus,
                root_minus,
            ]

        # ---------------------------------------------------------------------
        # EXACT ROOT SET MATCH
        # ---------------------------------------------------------------------

        resultant_root_set_match = False

        root_true_matches = False
        root_other_matches = False

        if degree_R == 2:

            root_true_matches = (
                root_plus == R_true
                or root_minus == R_true
            )

            root_other_matches = (
                root_plus == R_other
                or root_minus == R_other
            )

            resultant_root_set_match = root_set_match(
                root_plus,
                root_minus,
                R_true,
                R_other,
            )

        # ---------------------------------------------------------------------
        # Vieta CHECKS
        # ---------------------------------------------------------------------

        if degree_R == 2:

            direct_sum = sp.cancel(
                R_true + R_other
            )

            direct_product = sp.cancel(
                R_true * R_other
            )

            resultant_sum = sp.cancel(
                -b / a
            )

            resultant_product = sp.cancel(
                c / a
            )

            vieta_sum_pass = (
                direct_sum == resultant_sum
            )

            vieta_product_pass = (
                direct_product == resultant_product
            )

        else:

            direct_sum = None
            direct_product = None
            resultant_sum = None
            resultant_product = None
            vieta_sum_pass = False
            vieta_product_pass = False

        # ---------------------------------------------------------------------
        # RESULTANT FACTOR CHECK
        # ---------------------------------------------------------------------

        factor_list = sp.factor_list(
            specialized_poly.as_expr(),
            R,
        )

        nonconstant_factors = []

        for factor, multiplicity in factor_list[1]:

            if sp.Poly(
                factor,
                R,
                domain=sp.QQ,
            ).degree() > 0:

                nonconstant_factors.append(
                    (
                        sp.Poly(
                            factor,
                            R,
                            domain=sp.QQ,
                        ).degree(),
                        multiplicity,
                    )
                )

        # ---------------------------------------------------------------------
        # PRINT
        # ---------------------------------------------------------------------

        print("BASIC DATA")
        print(f"  p       = {p}")
        print(f"  q       = {q}")
        print(f"  N       = {N}")
        print(f"  X_TRUE  = {X_true}")
        print(f"  K_TRUE  = {K_true}")
        print()

        print("KAPPA ROOT CHECK")
        print(
            f"  t_TRUE = "
            f"{'PASS' if kappa_true_check == 0 else 'FAIL'}"
        )
        print(
            f"  t_OTHER = "
            f"{'PASS' if kappa_other_check == 0 else 'FAIL'}"
        )
        print()

        print("LAYER BRANCHES")
        print(
            f"  R_TRUE numerator digits = "
            f"{len(str(abs(sp.numer(R_true))))}"
        )
        print(
            f"  R_TRUE denominator digits = "
            f"{len(str(abs(sp.denom(R_true))))}"
        )
        print(
            f"  R_OTHER numerator digits = "
            f"{len(str(abs(sp.numer(R_other))))}"
        )
        print(
            f"  R_OTHER denominator digits = "
            f"{len(str(abs(sp.denom(R_other))))}"
        )
        print()

        print("SPECIALIZED RESULTANT")
        print(
            f"  degree in R = {degree_R}"
        )
        print(
            f"  A = leading coefficient"
        )
        print(
            f"  A bit-length = {abs(int(a)).bit_length()}"
        )
        print(
            f"  B bit-length = {abs(int(b)).bit_length()}"
        )
        print(
            f"  C bit-length = {abs(int(c)).bit_length()}"
        )
        print()

        print("RESULTANT DIRECT EVALUATION")
        print(
            f"  F(K_TRUE,R_TRUE) = "
            f"{'0' if resultant_at_R_true == 0 else 'NONZERO'}"
        )
        print(
            f"  F(K_TRUE,R_OTHER) = "
            f"{'0' if resultant_at_R_other == 0 else 'NONZERO'}"
        )
        print()

        print("RESULTANT ROOTS")
        print(
            f"  root count = {len(resultant_roots)}"
        )

        if degree_R == 2:

            print(
                f"  root 1 == R_TRUE  : "
                f"{root_plus == R_true}"
            )

            print(
                f"  root 1 == R_OTHER : "
                f"{root_plus == R_other}"
            )

            print(
                f"  root 2 == R_TRUE  : "
                f"{root_minus == R_true}"
            )

            print(
                f"  root 2 == R_OTHER : "
                f"{root_minus == R_other}"
            )

        print()

        print("SET COMPARISON")
        print(
            f"  resultant roots == "
            f"{{R_TRUE,R_OTHER}} : "
            f"{resultant_root_set_match}"
        )
        print()

        print("VIETA COMPARISON")
        print(
            f"  direct R sum == resultant sum     : "
            f"{vieta_sum_pass}"
        )
        print(
            f"  direct R product == resultant product : "
            f"{vieta_product_pass}"
        )
        print()

        print("DISCRIMINANT")
        print(
            f"  discriminant zero = "
            f"{discriminant == 0}"
        )
        print(
            f"  discriminant bit-length = "
            f"{abs(int(discriminant)).bit_length()}"
        )
        print()

        print("SPECIALIZED FACTOR STRUCTURE")
        print(
            f"  nonconstant factors = "
            f"{len(nonconstant_factors)}"
        )
        print(
            f"  factor signature = "
            f"{nonconstant_factors}"
        )
        print()

        results.append({
            "kappa_true": (
                kappa_true_check == 0
            ),
            "kappa_other": (
                kappa_other_check == 0
            ),
            "R_true_resultant": (
                resultant_at_R_true == 0
            ),
            "R_other_resultant": (
                resultant_at_R_other == 0
            ),
            "root_set_match": (
                resultant_root_set_match
            ),
            "vieta_sum": (
                vieta_sum_pass
            ),
            "vieta_product": (
                vieta_product_pass
            ),
            "degree": (
                degree_R
            ),
            "factor_signature": (
                nonconstant_factors
            ),
        })

    # =========================================================================
    # GLOBAL SUMMARY
    # =========================================================================

    print("=" * 110)
    print("GLOBAL SUMMARY")
    print("=" * 110)
    print()

    print(
        "  all TRUE Kappa roots pass = "
        f"{all(r['kappa_true'] for r in results)}"
    )

    print(
        "  all OTHER Kappa roots pass = "
        f"{all(r['kappa_other'] for r in results)}"
    )

    print(
        "  all R_TRUE resultant checks pass = "
        f"{all(r['R_true_resultant'] for r in results)}"
    )

    print(
        "  all R_OTHER resultant checks pass = "
        f"{all(r['R_other_resultant'] for r in results)}"
    )

    print(
        "  all resultant root sets match = "
        f"{all(r['root_set_match'] for r in results)}"
    )

    print(
        "  all Vieta sum checks pass = "
        f"{all(r['vieta_sum'] for r in results)}"
    )

    print(
        "  all Vieta product checks pass = "
        f"{all(r['vieta_product'] for r in results)}"
    )

    print(
        "  all specialized resultants quadratic = "
        f"{all(r['degree'] == 2 for r in results)}"
    )

    print()

    print("=" * 110)
    print("EXPERIMENT 387 FINAL STATUS")
    print("=" * 110)
    print()

    print("Interpretation:")
    print()
    print("  The crucial distinction is now:")
    print()
    print("    F_N(K_TRUE,R) = 0")
    print()
    print("  versus")
    print()
    print("    R = R(t_TRUE) or R(t_OTHER).")
    print()
    print("  Because K_TRUE is substituted before any normalization,")
    print("  no rational-function coefficient extraction is involved.")
    print()
    print("  CASE A:")
    print("    Both direct layer branches are roots.")
    print()
    print("    Then the previous Vieta failure was purely an artifact")
    print("    of the symbolic normalization procedure.")
    print()
    print("  CASE B:")
    print("    R_TRUE is a root but R_OTHER is not.")
    print()
    print("    Then the resultant contains only one layer branch at")
    print("    K_TRUE, despite the two Kappa t-roots.")
    print()
    print("  CASE C:")
    print("    Neither branch matches.")
    print()
    print("    Then the layer equation or resultant construction has")
    print("    a structural mismatch that must be resolved first.")
    print()
    print("=" * 110)
    print("EXPERIMENT 387 FINISHED")
    print("=" * 110)


if __name__ == "__main__":
    main()