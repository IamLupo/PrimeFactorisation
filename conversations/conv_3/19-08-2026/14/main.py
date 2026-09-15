import sympy as sp


# =============================================================================
# EXPERIMENT 391 START
# =============================================================================

print("=" * 110)
print("EXPERIMENT 391 START")
print("=" * 110)
print()
print("EXACT MOEBIUS REDUCTION / TWO-ROOT SECANT STRUCTURE")
print()
print("Goal:")
print("  1. Reduce R(t) modulo the exact Kappa quadratic.")
print("  2. Obtain:")
print("       R(t) = (a*t + b)/(c*t + d)")
print("  3. Verify the exact secant identity:")
print("       (R1-R2)/(t1-t2)")
print("         = (a*d-b*c)/((c*t1+d)(c*t2+d))")
print("  4. Express the denominator using:")
print("       t1+t2 = N(N+3)/(K+N^2+2)")
print("       t1*t2 = N^2/(K+N^2+2)")
print("  5. Test the symmetric R-sum and R-product.")
print()


# =============================================================================
# SYMBOL
# =============================================================================

t = sp.symbols("t")


# =============================================================================
# INSTANCES
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
# LAYER FUNCTION
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

    # Same rational layer used in the previous experiments.
    numerator = sp.expand(h15**2 - h16*h14)
    denominator = sp.expand(h16**2)

    return numerator, denominator


R_NUM, R_DEN = build_R()


# =============================================================================
# INSTANCE DATA
# =============================================================================

def instance_data(p, q):
    N = sp.Integer(p * q)

    X_TRUE = sp.Integer(p + q + 1)

    K_TRUE = sp.Integer(
        1
        - (p**2 - p + 1)
        * (q**2 - q + 1)
    )

    t_true = sp.cancel(N / X_TRUE)

    return N, X_TRUE, K_TRUE, t_true


# =============================================================================
# KAPPA QUADRATIC
# =============================================================================

def kappa_quadratic(N, K):
    return sp.Poly(
        sp.expand(
            (K + N**2 + 2) * t**2
            - N * (N + 3) * t
            + N**2
        ),
        t,
        domain=sp.QQ,
    )


# =============================================================================
# POLYNOMIAL REDUCTION
# =============================================================================

def reduce_mod_quadratic(expr, Q):
    poly = sp.Poly(
        sp.expand(expr),
        t,
        domain=sp.QQ,
    )

    remainder = sp.rem(poly, Q)

    return sp.expand(remainder.as_expr())


# =============================================================================
# MOBIUS REDUCTION
# =============================================================================

def reduce_R_mod_Q(N, K):
    Q = kappa_quadratic(N, K)

    num_red = reduce_mod_quadratic(
        R_NUM,
        Q,
    )

    den_red = reduce_mod_quadratic(
        R_DEN,
        Q,
    )

    num_poly = sp.Poly(
        num_red,
        t,
        domain=sp.QQ,
    )

    den_poly = sp.Poly(
        den_red,
        t,
        domain=sp.QQ,
    )

    a = sp.cancel(
        num_poly.coeff_monomial(t)
    )

    b = sp.cancel(
        num_poly.coeff_monomial(1)
    )

    c = sp.cancel(
        den_poly.coeff_monomial(t)
    )

    d = sp.cancel(
        den_poly.coeff_monomial(1)
    )

    return Q, a, b, c, d


# =============================================================================
# OTHER KAPPA ROOT
# =============================================================================

def other_kappa_root(N, K, t_true):
    coefficient = K + N**2 + 2

    t_other = sp.cancel(
        N**2 / (coefficient * t_true)
    )

    return t_other


# =============================================================================
# DIRECT R VALUE
# =============================================================================

def R_value(tv):
    numerator = sp.cancel(
        R_NUM.subs(t, tv)
    )

    denominator = sp.cancel(
        R_DEN.subs(t, tv)
    )

    return sp.cancel(
        numerator / denominator
    )


# =============================================================================
# DIGIT HELPERS
# =============================================================================

def num_den_digits(value):
    numerator = abs(int(sp.numer(value)))
    denominator = abs(int(sp.denom(value)))

    return (
        len(str(numerator)),
        len(str(denominator)),
    )


def numerator_digits(value):
    return num_den_digits(value)[0]


def denominator_digits(value):
    return num_den_digits(value)[1]


# =============================================================================
# MAIN
# =============================================================================

def main():

    all_roots_pass = True
    all_reduction_pass = True
    all_mobius_pass = True
    all_secant_pass = True
    all_symmetric_pass = True

    determinant_sizes = []
    denominator_sizes = []
    H_sizes = []

    for index, (p, q) in enumerate(
        INSTANCES,
        start=1,
    ):

        print()
        print("=" * 110)
        print(f"INSTANCE {index}")
        print("=" * 110)
        print()

        # ---------------------------------------------------------------------
        # Basic instance
        # ---------------------------------------------------------------------

        N, X_TRUE, K_TRUE, t_true = instance_data(
            p,
            q,
        )

        t_other = other_kappa_root(
            N,
            K_TRUE,
            t_true,
        )

        Q = kappa_quadratic(
            N,
            K_TRUE,
        )

        root_true_ok = (
            sp.expand(
                Q.as_expr().subs(t, t_true)
            )
            == 0
        )

        root_other_ok = (
            sp.expand(
                Q.as_expr().subs(t, t_other)
            )
            == 0
        )

        roots_ok = (
            root_true_ok
            and root_other_ok
        )

        # ---------------------------------------------------------------------
        # Möbius reduction
        # ---------------------------------------------------------------------

        Q, a, b, c, d = reduce_R_mod_Q(
            N,
            K_TRUE,
        )

        reduced_num = sp.expand(
            a * t + b
        )

        reduced_den = sp.expand(
            c * t + d
        )

        reduced_R = sp.cancel(
            reduced_num / reduced_den
        )

        # Verify modulo Q:
        #
        # R_NUM/R_DEN == reduced_R modulo Q
        #
        # Cross multiply and reduce the numerator modulo Q.

        cross_expr = sp.together(
            R_NUM
            - reduced_R * R_DEN
        )

        cross_num = cross_expr.as_numer_denom()[0]

        cross_poly = sp.Poly(
            sp.expand(cross_num),
            t,
            domain=sp.QQ,
        )

        cross_remainder = sp.rem(
            cross_poly,
            Q,
        )

        reduction_ok = (
            sp.expand(
                cross_remainder.as_expr()
            )
            == 0
        )

        # ---------------------------------------------------------------------
        # Direct branches
        # ---------------------------------------------------------------------

        R_true = R_value(t_true)
        R_other = R_value(t_other)

        R_true_reduced = sp.cancel(
            reduced_R.subs(
                t,
                t_true,
            )
        )

        R_other_reduced = sp.cancel(
            reduced_R.subs(
                t,
                t_other,
            )
        )

        branch_true_ok = (
            R_true_reduced == R_true
        )

        branch_other_ok = (
            R_other_reduced == R_other
        )

        mobius_ok = (
            branch_true_ok
            and branch_other_ok
        )

        # ---------------------------------------------------------------------
        # Symmetric Kappa root quantities
        # ---------------------------------------------------------------------

        gamma = sp.cancel(
            K_TRUE + N**2 + 2
        )

        t_sum = sp.cancel(
            N * (N + 3) / gamma
        )

        t_product = sp.cancel(
            N**2 / gamma
        )

        direct_t_sum = sp.cancel(
            t_true + t_other
        )

        direct_t_product = sp.cancel(
            t_true * t_other
        )

        t_sum_ok = (
            t_sum == direct_t_sum
        )

        t_product_ok = (
            t_product == direct_t_product
        )

        # ---------------------------------------------------------------------
        # Möbius determinant
        # ---------------------------------------------------------------------

        determinant = sp.cancel(
            a * d - b * c
        )

        # ---------------------------------------------------------------------
        # Symmetric denominator
        #
        # D = (c*t1+d)(c*t2+d)
        #   = c^2*t1*t2 + cd*(t1+t2) + d^2
        # ---------------------------------------------------------------------

        symmetric_denominator = sp.cancel(
            c**2 * t_product
            + c * d * t_sum
            + d**2
        )

        direct_denominator = sp.cancel(
            (
                c * t_true + d
            )
            * (
                c * t_other + d
            )
        )

        symmetric_denominator_ok = (
            symmetric_denominator
            == direct_denominator
        )

        # ---------------------------------------------------------------------
        # Secant slope
        # ---------------------------------------------------------------------

        delta_t = sp.cancel(
            t_true - t_other
        )

        delta_R = sp.cancel(
            R_true - R_other
        )

        H_direct = sp.cancel(
            delta_R / delta_t
        )

        H_mobius = sp.cancel(
            determinant
            / symmetric_denominator
        )

        H_ok = (
            H_direct == H_mobius
        )

        # ---------------------------------------------------------------------
        # Symmetric R sum
        #
        # R1+R2 =
        #
        # [2ac*t1*t2 + (ad+bc)(t1+t2) + 2bd] / D
        # ---------------------------------------------------------------------

        R_sum_predicted = sp.cancel(
            (
                2 * a * c * t_product
                + (a * d + b * c) * t_sum
                + 2 * b * d
            )
            / symmetric_denominator
        )

        R_sum_direct = sp.cancel(
            R_true + R_other
        )

        R_sum_ok = (
            R_sum_predicted
            == R_sum_direct
        )

        # ---------------------------------------------------------------------
        # Symmetric R product
        #
        # R1*R2 =
        #
        # [a^2*t1*t2 + ab(t1+t2) + b^2] / D
        # ---------------------------------------------------------------------

        R_product_predicted = sp.cancel(
            (
                a**2 * t_product
                + a * b * t_sum
                + b**2
            )
            / symmetric_denominator
        )

        R_product_direct = sp.cancel(
            R_true * R_other
        )

        R_product_ok = (
            R_product_predicted
            == R_product_direct
        )

        symmetric_ok = (
            t_sum_ok
            and t_product_ok
            and symmetric_denominator_ok
            and R_sum_ok
            and R_product_ok
        )

        # ---------------------------------------------------------------------
        # Size data
        # ---------------------------------------------------------------------

        determinant_size = num_den_digits(
            determinant
        )

        denominator_size = num_den_digits(
            symmetric_denominator
        )

        H_size = num_den_digits(
            H_direct
        )

        determinant_sizes.append(
            determinant_size
        )

        denominator_sizes.append(
            denominator_size
        )

        H_sizes.append(
            H_size
        )

        # ---------------------------------------------------------------------
        # Print basic data
        # ---------------------------------------------------------------------

        print("BASIC DATA")
        print(f"  p       = {p}")
        print(f"  q       = {q}")
        print(f"  N       = {N}")
        print(f"  X_TRUE  = {X_TRUE}")
        print(f"  K_TRUE  = {K_TRUE}")
        print()

        print("KAPPA ROOTS")
        print(f"  t_TRUE  = {t_true}")
        print(f"  t_OTHER = {t_other}")
        print(
            "  t_TRUE root verification  = "
            f"{root_true_ok}"
        )
        print(
            "  t_OTHER root verification = "
            f"{root_other_ok}"
        )
        print()

        print("INSTANCE-LOCAL MOBIUS REDUCTION")
        print(
            "  R(t) mod Q = (a*t+b)/(c*t+d)"
        )
        print(
            "  numerator degree   = "
            f"{sp.Poly(reduced_num, t).degree()}"
        )
        print(
            "  denominator degree = "
            f"{sp.Poly(reduced_den, t).degree()}"
        )
        print(
            "  exact reduction identity = "
            f"{reduction_ok}"
        )
        print()

        print("MOEBIUS COEFFICIENT SIZES")
        print(
            "  a numerator digits = "
            f"{numerator_digits(a)}"
        )
        print(
            "  b numerator digits = "
            f"{numerator_digits(b)}"
        )
        print(
            "  c numerator digits = "
            f"{numerator_digits(c)}"
        )
        print(
            "  d numerator digits = "
            f"{numerator_digits(d)}"
        )
        print()

        print("MOEBIUS BRANCH RECONSTRUCTION")
        print(
            "  R_TRUE reduced branch  = "
            f"{branch_true_ok}"
        )
        print(
            "  R_OTHER reduced branch = "
            f"{branch_other_ok}"
        )
        print(
            "  both branches exact     = "
            f"{mobius_ok}"
        )
        print()

        print("KAPPA SYMMETRIC ROOT DATA")
        print(
            "  t1+t2 identity = "
            f"{t_sum_ok}"
        )
        print(
            "  t1*t2 identity = "
            f"{t_product_ok}"
        )
        print()

        print("MOEBIUS DETERMINANT")
        print(
            "  determinant = a*d - b*c"
        )
        print(
            "  numerator digits   = "
            f"{determinant_size[0]}"
        )
        print(
            "  denominator digits = "
            f"{determinant_size[1]}"
        )
        print()

        print("SYMMETRIC MOEBIUS DENOMINATOR")
        print(
            "  D = (c*t1+d)(c*t2+d)"
        )
        print(
            "  symmetric reconstruction = "
            f"{symmetric_denominator_ok}"
        )
        print(
            "  numerator digits   = "
            f"{denominator_size[0]}"
        )
        print(
            "  denominator digits = "
            f"{denominator_size[1]}"
        )
        print()

        print("SECANT SLOPE")
        print(
            "  H_direct == H_Mobius = "
            f"{H_ok}"
        )
        print(
            "  numerator digits   = "
            f"{H_size[0]}"
        )
        print(
            "  denominator digits = "
            f"{H_size[1]}"
        )
        print()

        print("SYMMETRIC R-BRANCH TESTS")
        print(
            "  R1+R2 formula = "
            f"{R_sum_ok}"
        )
        print(
            "  R1*R2 formula = "
            f"{R_product_ok}"
        )
        print(
            "  all symmetric formulas = "
            f"{symmetric_ok}"
        )
        print()

        # ---------------------------------------------------------------------
        # Simple normalization sizes
        # ---------------------------------------------------------------------

        determinant_over_gamma = sp.cancel(
            determinant / gamma
        )

        H_times_gamma = sp.cancel(
            H_direct * gamma
        )

        H_over_N = sp.cancel(
            H_direct / N
        )

        determinant_over_N = sp.cancel(
            determinant / N
        )

        print("NORMALIZATION CANDIDATES")
        print(
            "  determinant/gamma numerator digits = "
            f"{numerator_digits(determinant_over_gamma)}"
        )
        print(
            "  determinant/gamma denominator digits = "
            f"{denominator_digits(determinant_over_gamma)}"
        )

        print(
            "  H*gamma numerator digits = "
            f"{numerator_digits(H_times_gamma)}"
        )
        print(
            "  H*gamma denominator digits = "
            f"{denominator_digits(H_times_gamma)}"
        )

        print(
            "  H/N numerator digits = "
            f"{numerator_digits(H_over_N)}"
        )
        print(
            "  H/N denominator digits = "
            f"{denominator_digits(H_over_N)}"
        )

        print(
            "  determinant/N numerator digits = "
            f"{numerator_digits(determinant_over_N)}"
        )
        print(
            "  determinant/N denominator digits = "
            f"{denominator_digits(determinant_over_N)}"
        )
        print()

        # ---------------------------------------------------------------------
        # Accumulate
        # ---------------------------------------------------------------------

        all_roots_pass &= roots_ok
        all_reduction_pass &= reduction_ok
        all_mobius_pass &= mobius_ok
        all_secant_pass &= H_ok
        all_symmetric_pass &= symmetric_ok

    # =========================================================================
    # GLOBAL SUMMARY
    # =========================================================================

    print("=" * 110)
    print("GLOBAL SUMMARY")
    print("=" * 110)
    print()

    print(
        "  all Kappa root checks pass       = "
        f"{all_roots_pass}"
    )

    print(
        "  all Möbius reductions pass       = "
        f"{all_reduction_pass}"
    )

    print(
        "  all branch reconstructions pass  = "
        f"{all_mobius_pass}"
    )

    print(
        "  all secant identities pass       = "
        f"{all_secant_pass}"
    )

    print(
        "  all symmetric R identities pass  = "
        f"{all_symmetric_pass}"
    )

    print()

    print("DETERMINANT SIZE RANGE")
    print(
        "  minimum numerator digits   = "
        f"{min(v[0] for v in determinant_sizes)}"
    )
    print(
        "  maximum numerator digits   = "
        f"{max(v[0] for v in determinant_sizes)}"
    )
    print(
        "  minimum denominator digits = "
        f"{min(v[1] for v in determinant_sizes)}"
    )
    print(
        "  maximum denominator digits = "
        f"{max(v[1] for v in determinant_sizes)}"
    )
    print()

    print("SYMMETRIC DENOMINATOR SIZE RANGE")
    print(
        "  minimum numerator digits   = "
        f"{min(v[0] for v in denominator_sizes)}"
    )
    print(
        "  maximum numerator digits   = "
        f"{max(v[0] for v in denominator_sizes)}"
    )
    print(
        "  minimum denominator digits = "
        f"{min(v[1] for v in denominator_sizes)}"
    )
    print(
        "  maximum denominator digits = "
        f"{max(v[1] for v in denominator_sizes)}"
    )
    print()

    print("SECANT H SIZE RANGE")
    print(
        "  minimum numerator digits   = "
        f"{min(v[0] for v in H_sizes)}"
    )
    print(
        "  maximum numerator digits   = "
        f"{max(v[0] for v in H_sizes)}"
    )
    print(
        "  minimum denominator digits = "
        f"{min(v[1] for v in H_sizes)}"
    )
    print(
        "  maximum denominator digits = "
        f"{max(v[1] for v in H_sizes)}"
    )
    print()

    print("=" * 110)
    print("EXPERIMENT 391 FINAL STATUS")
    print("=" * 110)
    print()
    print("The experiment tests the exact reduction:")
    print()
    print("  R(t) mod Q(t) = (a*t+b)/(c*t+d)")
    print()
    print("and the induced two-root identity:")
    print()
    print("  (R1-R2)/(t1-t2)")
    print("      = (a*d-b*c)")
    print("        ----------------")
    print("        (c*t1+d)(c*t2+d)")
    print()
    print("The two-root symmetric denominator is evaluated only through")
    print("the exact Kappa Vieta quantities.")
    print()
    print("No degree-20 resultant is constructed.")
    print()
    print("=" * 110)
    print("EXPERIMENT 391 FINISHED")
    print("=" * 110)


if __name__ == "__main__":
    main()